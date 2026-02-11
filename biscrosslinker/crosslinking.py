"""
Core crosslinking operations for connecting polymer chains
"""
import random
import MDAnalysis as mda
from .chemistry import (
    unpack_linking_site, 
    calculate_reactive_direction, 
    get_reactive_group,
    create_chain_cylinder,
    scan_chain
)
from .geometry import align_molecule, apply_random_rotation, rotate_by_angle
from .loops import coarse_circle_radius, find_possible_crosses, does_cyl_cross_plane, find_point_of_intersection
from .config import MAX_ROTATION_ATTEMPTS, RAND_CHAIN_ROTATION, MIN_RAND_ANGLE, MAX_RAND_ANGLE


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
    bis_attachment_pos = bis_attachment.position.copy()  # Cache position
    bis_dir_template = (end_carbon.position - bis_attachment_pos).copy()

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
            
            # Use cached bis direction
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

            new_cylinder = create_chain_cylinder(atoms_to_keep)

            # Use network's collision checking
            if network.check_collision(new_cylinder):
                continue  # try again

# LOOP LOGIC ==============================================================================
            circle_radius = coarse_circle_radius(end_carbon.position, new_cylinder)
            possible_hits = find_possible_crosses(end_carbon.position, circle_radius, network)
            good_hits, good_points = [], []
            has_looped = False
            for cylinder in possible_hits:
                does_cross_plane, t = does_cyl_cross_plane(bis_dir, cylinder, end_carbon.position)
                if not has_looped and does_cross_plane:
                    good_hits.append(cylinder)
                    good_points.append(find_point_of_intersection(cylinder, t))
                    rotate_by_angle(atoms_to_keep, bis_dir, 0., end_carbon.position)
                    has_looped = True
                new_cylinder = create_chain_cylinder(atoms_to_keep)
                if network.check_collision(new_cylinder):
                    has_looped = False
                    continue  # try again
# LOOP LOGIC ==============================================================================

            return mda.Merge(atoms_to_keep, structure_univ.atoms), new_cylinder

    raise RuntimeError("Failed to place chain without collision")