"""
Core crosslinking operations for connecting polymer chains
"""
import random
import numpy as np
import MDAnalysis as mda
from .chemistry import (
    unpack_linking_site, 
    calculate_reactive_direction, 
    get_reactive_group,
    create_chain_cylinder,
    scan_chain
)
from .geometry import align_molecule, apply_random_rotation, rotate_by_angle, r12, point_segment_distance
from .loops import coarse_circle_radius, find_possible_crosses, does_cyl_cross_plane, find_point_of_intersection
from .config import MAX_ROTATION_ATTEMPTS, RAND_CHAIN_ROTATION, MIN_RAND_ANGLE, MAX_RAND_ANGLE, MAX_BIS_ORIENTATION_ANGLE, ENABLE_LOOPING


def replace_with_bis(linking_site, chain_universe):
    """
    Replace a reactive site with a BIS crosslinker molecule
    
    Args:
        linking_site: Reactive site dictionary from scan_chain()
        chain_universe: MDAnalysis Universe containing the chain
        
    Returns:
        Tuple of (merged_universe, bis_universe)
    """
    # Unpack geometry
    site = unpack_linking_site(linking_site)
    reactive_dir = calculate_reactive_direction(site)
    
    # Remove reactive group
    reactive_group = get_reactive_group(site)
    atoms_to_keep = chain_universe.select_atoms('all').subtract(reactive_group)
    
    # Load and orient BIS
    bis_univ = mda.Universe("BIS.pdb", to_guess=['elements', 'bonds'])
    bis_attachment = bis_univ.select_atoms('index 5')[0]
    end_carbon = bis_univ.select_atoms('index 6')[0]
    bis_dir = end_carbon.position - bis_attachment.position
    
    # Align BIS with reactive site
    align_molecule(bis_univ, bis_dir, reactive_dir)
    
    # Position BIS at reactive site
    displacement = site['reactive_carbon'].position - bis_attachment.position
    bis_univ.atoms.translate(displacement)
    
    return mda.Merge(atoms_to_keep, bis_univ.atoms), bis_univ


def connect_PAAm(bis_univ, structure_univ, network, template,
                 max_rotation_attempts=MAX_ROTATION_ATTEMPTS, 
                 random_rotation=RAND_CHAIN_ROTATION):
    """
    Connect a new PAAm chain to an existing BIS crosslinker
    
    Args:
        bis_univ: MDAnalysis Universe containing the BIS crosslinker
        structure_univ: Current network structure
        network: CrosslinkNetwork instance
        template: PAAmTemplate instance
        max_rotation_attempts: Maximum number of rotation attempts per site
        random_rotation: Whether to use random rotation (vs systematic)
        
    Returns:
        Tuple of (merged_universe, new_chain_cylinder)
        
    Raises:
        RuntimeError: If unable to place chain without collision
    """
    max_site_attempts = 10
    
    # Cache BIS structure
    bis_attachment = bis_univ.select_atoms('index 6')[0]
    end_carbon = bis_univ.select_atoms('index 5')[0]
    bis_attachment_pos = bis_attachment.position.copy()
    bis_dir_template = (end_carbon.position - bis_attachment_pos).copy()

    # The cylinder the BIS is attached to is the most recently added one
    parent_cylinder = network.cylinders[-1] if network.cylinders else None

    for site_attempt in range(max_site_attempts):
        # Use pre-computed valid sites from template
        if not template.valid_site_indices:
            continue
            
        # Pick which site index to use (from pre-scanned list)
        rand_group_index = random.choice(template.valid_site_indices)
        
        # Create universe once per site attempt
        new_PAAm_univ = template.get_fresh_universe()
        
        # Scan once per site attempt  
        linking_groups_new = scan_chain(new_PAAm_univ)
        linking_site = linking_groups_new[rand_group_index]
        
        # Get initial site geometry
        site = unpack_linking_site(linking_site)
        initial_reactive_dir = calculate_reactive_direction(site)

        ANGLE = 0

        for rotation_attempt in range(max_rotation_attempts):            
            # Reset to original positions
            if rotation_attempt > 0:
                new_PAAm_univ.atoms.positions = template.original_positions.copy()
            
            # Recalculate reactive direction (positions reset, so recalc needed)
            site = unpack_linking_site(linking_site)
            reactive_dir = calculate_reactive_direction(site)
            
            # Use cached BIS direction
            bis_dir = bis_dir_template.copy()

            align_molecule(new_PAAm_univ, reactive_dir, bis_dir)

            if random_rotation:
                apply_random_rotation(new_PAAm_univ, bis_dir, MIN_RAND_ANGLE, MAX_RAND_ANGLE)
            else:
                rotate_by_angle(new_PAAm_univ.select_atoms('all'), bis_dir, angle=ANGLE)
                ANGLE += 360 / MAX_ROTATION_ATTEMPTS + 1

            reactive_group = get_reactive_group(site)
            atoms_to_keep = new_PAAm_univ.select_atoms('all').subtract(reactive_group)

            displacement = bis_attachment_pos - site['reactive_carbon'].position
            atoms_to_keep.translate(displacement)

            cts_atoms = atoms_to_keep.select_atoms('name CTS')
            new_cylinder = create_chain_cylinder(
                atoms_to_keep,
                start_index=cts_atoms[0].index,
                end_index=cts_atoms[1].index
            )
            monomer_length = new_cylinder.length / 25 #HARD CODED ==> NEED TO CHANGE FOR DIFF CHAIN LENGTHS
            # NAIVE WAY TO MEASURE MONOMER LENGTH, BETTER TO GET DISTANCES FROM CIRCLE CENTER

            # Use network's collision checking, skipping the parent cylinder
            if network.check_collision(new_cylinder, skip_cylinder=parent_cylinder):
                continue  # try again
            
# LOOP LOGIC ==============================================================================
            if ENABLE_LOOPING:
                circle_radius = coarse_circle_radius(end_carbon.position, new_cylinder)
                possible_hits = find_possible_crosses(end_carbon.position, circle_radius, network)
                has_looped = False
                
                for cylinder in possible_hits:
                    does_cross_plane, intersection_point = does_cyl_cross_plane(bis_dir, cylinder, end_carbon.position)
                    
                    if not has_looped and does_cross_plane:
                        # Find best reactive site on existing chain near intersection point
                        from .loops import find_best_existing_site, find_best_new_chain_site, calculate_loop_rotation_angle, check_bis_orientation
                        
                        existing_site, distance_to_pivot_existing = find_best_existing_site(
                            intersection_point, 
                            end_carbon.position, 
                            structure_univ,
                            max_search_radius=monomer_length * 2  # Search within ~2 monomers
                        )
                        
                        if not existing_site:
                            continue  # No valid site found on existing chain
                        
                        # Find best reactive site on new chain
                        new_chain_site, distance_to_pivot_new = find_best_new_chain_site(
                            new_PAAm_univ,
                            end_carbon.position,
                            max_distance_from_pivot=distance_to_pivot_existing
                        )
                        
                        if not new_chain_site:
                            continue  # No valid site found on new chain
                        
                        # Get positions for angle calculation
                        existing_site_pos = existing_site['central_carbon'].position
                        new_chain_site_pos = new_chain_site['central_carbon'].position
                        pivot_pos = end_carbon.position
                        
                        # Calculate rotation angle using law of cosines
                        rotation_angle = calculate_loop_rotation_angle(
                            pivot_pos,
                            new_chain_site_pos,
                            existing_site_pos,
                            BIS_LENGTH=7.0
                        )
                        
                        if rotation_angle is None:
                            continue  # Invalid triangle -> skip this loop opportunity
                        
                        # Determine the direction of rotation
                        vec_to_new_site = new_chain_site_pos - pivot_pos
                        vec_to_existing_site = existing_site_pos - pivot_pos
                        
                        # If cross(vec_to_new, vec_to_existing) points along bis_dir, rotate positive
                        # If it points opposite, rotate negative
                        cross_product = np.cross(vec_to_new_site, vec_to_existing_site)
                        bis_dir_normalized = bis_dir / np.linalg.norm(bis_dir)
                        
                        if np.dot(cross_product, bis_dir_normalized) < 0:
                            rotation_angle = -rotation_angle  # Reverse direction
                        
                        # Rotate the entire new chain around the pivot point
                        rotate_by_angle(
                            atoms_to_keep,
                            bis_dir_normalized,
                            rotation_angle,
                            point=pivot_pos
                        )

                        # Point roughly toward each other before accepting this loop
                        new_chain_site_unpacked      = unpack_linking_site(new_chain_site)
                        existing_chain_site_unpacked = unpack_linking_site(existing_site)

                        if not check_bis_orientation(
                            new_chain_site_unpacked,
                            existing_chain_site_unpacked,
                            MAX_BIS_ORIENTATION_ANGLE
                        ):
                            print(f'  Reactive sites not oriented toward each other, skipping.')
                            continue
                        
                        # TODO: Need to actually connect the chains chemically
                        # 1. Removing the reactive groups from both sites
                        # 2. Adding a BIS molecule between them
                        # 3. Merging everything together
                        
                        has_looped = True
                        print('  Chain rotated into loop position! ===========================================================')
                    
                    # Re-check collision after rotation, still skipping the parent cylinder
                    cts_atoms = atoms_to_keep.select_atoms('name CTS')
                    new_cylinder = create_chain_cylinder(
                        atoms_to_keep,
                        start_index=cts_atoms[0].index,
                        end_index=cts_atoms[1].index
                    )
                    if network.check_collision(new_cylinder, skip_cylinder=parent_cylinder):
                        has_looped = False
                        print('  Collision detected after rotation, continuing search...')
                        continue  # try again
# LOOP LOGIC ==============================================================================

            return mda.Merge(atoms_to_keep, structure_univ.atoms), new_cylinder

    raise RuntimeError("Failed to place chain without collision")