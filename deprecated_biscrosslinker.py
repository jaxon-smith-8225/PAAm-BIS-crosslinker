import MDAnalysis as mda
import numpy as np
import random
import warnings
warnings.filterwarnings('ignore')
import time

OUTPUT_STRUCTURE = "output.pdb"
MIN_RAND_ANGLE = 20
MAX_RAND_ANGLE = 200
RAND_CHAIN_ROTATION = True
CHAIN_RADIUS = 3.7
MAX_NUM_CHAINS = 20
MAX_ROATAION_ATTEMPTS = 25

class ChainCylinder:
    """Represents a cylindrical bounding volume for a polymer chain"""
    def __init__(self, start_point, end_point, radius=CHAIN_RADIUS):
        self.start = np.array(start_point)
        self.end = np.array(end_point)
        self.radius = radius

    @property
    def axis(self):
        """Vector along cylinder axis"""
        return self.end - self.start

    @property
    def length(self):
        return np.linalg.norm(self.axis)

# Math utils

def segment_segment_distance(a1, a2, b1, b2):
    a1, a2, b1, b2 = map(np.array, (a1, a2, b1, b2))
    da = a2 - a1
    db = b2 - b1
    r = b1 - a1

    cross = np.cross(da, db)
    if np.linalg.norm(cross) < 1e-8:
        distance = np.abs(np.cross(r, da)) / np.linalg.norm(da)
    else:
         distance = abs(np.dot(r, cross)) / np.linalg.norm(cross)
    return distance

def cylinders_collide(cyl1, cyl2):
    """
    Returns True if two ChainCylinder objects overlap
    """
    dist = segment_segment_distance(
        cyl1.start, cyl1.end,
        cyl2.start, cyl2.end
    )
    # print(dist)
    return dist < (cyl1.radius + cyl2.radius)

def check_collision(new_cylinder, cylinders):
    """
    Check new cylinder against all existing cylinders
    Return True if collision
    """
    for cyl in cylinders:
        if cylinders_collide(new_cylinder, cyl):
            return True
        else:
            continue
    return False

def angle_between_vectors(v1, v2): # theta = arccos( dot(v1, v2) / (norm(v1) * norm(v2)) )
    """Calculate angle (in degrees) between two vectors"""
    dot_prod = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    cosine_angle = dot_prod / (norm_v1 * norm_v2)
    clipped_cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    angle_radians = np.arccos(clipped_cosine_angle)
    angle_degrees = angle_radians * 180 / np.pi
    return angle_degrees

def rotate_by_angle(universe, axis, angle):
    """Rotate universe by specified angle around axis"""
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    universe.atoms.rotateby(angle, axis, point=[0, 0, 0])
    return

def align_molecule(universe, source_vector, target_vector):
    """
    Rotate universe to align source_vector with target_vector
    
    Args:
        universe: MDAnalysis Universe to rotate
        source_vector: Current orientation vector
        target_vector: Desired orientation vector
    """
    rotation_axis = np.cross(source_vector, target_vector)
    
    # Handle parallel/antiparallel vectors
    if np.linalg.norm(rotation_axis) < 1e-6:
        # Vectors are parallel or antiparallel
        if np.dot(source_vector, target_vector) < 0:
            # Find perpendicular axis for 180° rotation
            rotation_axis = np.array([1, 0, 0]) if abs(source_vector[0]) < 0.9 else np.array([0, 1, 0])
            rotation_angle = 180.0
        else:
            # Already aligned
            return
    else:
        rotation_angle = angle_between_vectors(target_vector, source_vector)
    
    rotate_by_angle(universe, rotation_axis, rotation_angle)

def apply_random_rotation(universe, axis, min_angle=MIN_RAND_ANGLE, max_angle=MAX_RAND_ANGLE):
    """
    Apply random rotation about an axis
    
    Args:
        universe: MDAnalysis Universe to rotate
        axis: Rotation axis vector
        min_angle: Minimum rotation angle (degrees)
        max_angle: Maximum rotation angle (degrees)
        
    Returns:
        The rotation angle that was applied
    """
    rand_rotation_angle = random.randint(min_angle, max_angle)
    # print(f'RANDOM ROTATION ANGLE: {rand_rotation_angle}')
    rotate_by_angle(universe, axis, rand_rotation_angle)
    return rand_rotation_angle

def is_site_too_close(candidate_site, crosslink_positions, min_distance=10.0):
    """
    Check if a candidate linking site is too close to existing crosslinks
    
    Args:
        candidate_site: linking site dictionary from scan_chain()
        crosslink_positions: list of numpy arrays with crosslink positions
        min_distance: minimum allowed distance in Angstroms
        
    Returns:
        True if too close (reject), False if okay (accept)
    """
    if not crosslink_positions:  # No existing crosslinks yet
        return False
    
    # Get position of this candidate's backbone carbon
    candidate_pos = candidate_site['backbone_carbon'][0].position
    
    # Check distance to all existing crosslinks
    for crosslink_pos in crosslink_positions:
        distance = np.linalg.norm(candidate_pos - crosslink_pos)
        if distance < min_distance:
            return True  # Too close -> reject this site
    
    return False  # All distances okay

# Structure analysis

def scan_chain(chain):
    """
    Scan a PAAm chain for reactive carbon sites
    
    Args:
        chain: MDAnalysis Universe containing the polymer chain
        
    Returns:
        List of dictionaries containing reactive group information:
        - central_carbon: reactive carbon atom
        - backbone_carbon: attached backbone carbon
        - oxygen: attached oxygen atom
        - nitrogen: attached nitrogen atom
        - hydrogens: hydrogen atoms bonded to nitrogen
    """
    carbons = chain.select_atoms('element C')
    reactive_groups = []

    for n_atom in carbons:

        if not hasattr(n_atom, 'bonds'):
            print('Warning: No bond information')
            return []
        
        bonded_atoms = n_atom.bonded_atoms
        if (bonded_atoms.select_atoms('element C').n_atoms == 1 and
            bonded_atoms.select_atoms('element O').n_atoms == 1 and
            bonded_atoms.select_atoms('element N').n_atoms == 1):

            backbone_carbon = bonded_atoms.select_atoms('element C')
            oxygen = bonded_atoms.select_atoms('element O')
            nitrogen = bonded_atoms.select_atoms('element N')
            hydrogens = nitrogen[0].bonded_atoms.select_atoms('element H')

            if hydrogens.n_atoms != 2:
                continue

            reactive_groups.append({
                'central_carbon': n_atom,
                'backbone_carbon': backbone_carbon,
                'oxygen': oxygen,
                'nitrogen': nitrogen,
                'hydrogens': hydrogens

            })

    return reactive_groups

def unpack_linking_site(linking_site):
    """Extract atoms from linking site dictionary"""
    return {
        'reactive_carbon': linking_site['central_carbon'],
        'backbone_carbon': linking_site['backbone_carbon'][0],
        'oxygen': linking_site['oxygen'][0],
        'nitrogen': linking_site['nitrogen'][0],
        'hydrogen1': linking_site['hydrogens'][0],
        'hydrogen2': linking_site['hydrogens'][1]
    }

def get_reactive_group(unpacked_site):
    """Create AtomGroup of atoms to remove"""
    return mda.AtomGroup([
        unpacked_site['reactive_carbon'],
        unpacked_site['oxygen'],
        unpacked_site['nitrogen'],
        unpacked_site['hydrogen1'],
        unpacked_site['hydrogen2']
    ])

def calculate_reactive_direction(unpacked_site):
    """Calculate direction vector from backbone to reactive carbon"""
    return (unpacked_site['reactive_carbon'].position - 
            unpacked_site['backbone_carbon'].position)

def create_chain_cylinder(chain_AtomGroup, radius=CHAIN_RADIUS):
    """
    chain_AtomGroup will almost always be the atoms_to_keep AtomGroup after it's been rotated and translated
    -> rotate + translate -> create cylinder around atoms_to_keep -> check collision -> continue
    """
    backbone_carbons = chain_AtomGroup.select_atoms('element C') 
    if len(backbone_carbons) < 2:
        # Fallback: use first and last atoms
        start_point = chain_AtomGroup[0].position
        end_point = chain_AtomGroup[-1].position
    else:
        # Use first and last backbone carbons
        start_point = backbone_carbons[0].position
        end_point = backbone_carbons[-1].position
    
    return ChainCylinder(start_point, end_point, radius)

# Cross linking operations

def replace_with_bis(linking_site, chain_universe):
    """Replace a reactive site with a BIS crosslinker molecule."""
    
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

def connect_PAAm(bis_univ, structure_univ, cylinders, max_rotation_attempts=MAX_ROATAION_ATTEMPTS, random_rotation=RAND_CHAIN_ROTATION):

    """
    Connect a new PAAm chain to an existing BIS crosslinker
    
    Args:
        bis_univ: MDAnalysis Universe containing the BIS crosslinker
        structure_univ: MDAnalysis Universe of the current structure
        random_rotation: If True, apply random rotation about linking axis
        
    Returns:
        Merged MDAnalysis Universe with new PAAm chain attached
    """
    max_site_attempts = 10
    # load new PAAm chain
    for site_attempt in range(max_site_attempts):
        temp_univ = mda.Universe("PAAm25mer.pdb", context='default', to_guess=['elements', 'bonds'])
        linking_groups = scan_chain(temp_univ)
        rand_group = random.randint(1, len(linking_groups) - 2)
        linking_site = linking_groups[rand_group]

        ANGLE = 0

        for rotation_attempt in range(max_rotation_attempts):            
            # Load fresh universe for each rotation attempt
            new_PAAm_univ = mda.Universe("PAAm25mer.pdb", context='default', to_guess=['elements', 'bonds'])
            
            linking_groups = scan_chain(new_PAAm_univ)
            linking_site = linking_groups[rand_group]
            site = unpack_linking_site(linking_site)

            reactive_dir = calculate_reactive_direction(site)

            bis_attachment = bis_univ.select_atoms('index 6')[0]
            end_carbon = bis_univ.select_atoms('index 5')[0]
            bis_dir = end_carbon.position - bis_attachment.position

            align_molecule(new_PAAm_univ, reactive_dir, bis_dir)

            if random_rotation:
                apply_random_rotation(new_PAAm_univ, bis_dir)
            else:
                # print(f'NON RANDOM ROTATION ANGLE: {ANGLE} degrees')
                rotate_by_angle(new_PAAm_univ, bis_dir, angle=ANGLE)
                ANGLE += 360 / MAX_ROATAION_ATTEMPTS + 1

            reactive_group = get_reactive_group(site)
            atoms_to_keep = new_PAAm_univ.select_atoms('all').subtract(reactive_group)

            displacement = bis_attachment.position - site['reactive_carbon'].position
            atoms_to_keep.translate(displacement)

            new_cylinder = create_chain_cylinder(atoms_to_keep)

            if check_collision(new_cylinder, cylinders):
                continue  # try again

            structure_atoms = structure_univ.select_atoms('all')
            return mda.Merge(atoms_to_keep, structure_atoms), new_cylinder

    raise RuntimeError("Failed to place PAAm without collision")


# Main ======================================================================================================

def main():
    cylinders = []
    crosslink_positions = []
    num_added_chains = 0
    while num_added_chains < MAX_NUM_CHAINS:
        if num_added_chains == 0:
            INPUT_STRUCTURE = "PAAm25mer.pdb"
        else:
            INPUT_STRUCTURE = OUTPUT_STRUCTURE

        # load polymer chain
        chain = mda.Universe(INPUT_STRUCTURE, context='default', to_guess=['elements', 'bonds'])

        if num_added_chains == 0:
            chain_cylinder = create_chain_cylinder(chain)
            cylinders.append(chain_cylinder)

        # find reactive sites and pick random middle residue
        linking_groups = scan_chain(chain)
        
        # filter out sites too close to existing crosslinks
        valid_sites = [site for site in linking_groups[1:-1]  # exclude first/last
                       if not is_site_too_close(site, crosslink_positions)]
        
        if not valid_sites:
            print("No valid sites found - all too close to existing crosslinks")
            break
        
        # pick from valid sites only
        linking_site = random.choice(valid_sites)
        backbone_carbon_pos = linking_site['backbone_carbon'][0].position.copy()
        crosslink_positions.append(backbone_carbon_pos)

        # replace reactive site w/ BIS and save

        first_modified_univ, bis_univ = replace_with_bis(linking_site, chain)
        try:
            final_modified_univ, new_cylinder = connect_PAAm(bis_univ, first_modified_univ, cylinders)
        except:
            continue
        # print(f'Added cylinder {new_cylinder} to list')
        cylinders.append(new_cylinder)
        # print(f'Cylinder lis: {cylinders}')
        final_modified_univ.select_atoms('all').write(OUTPUT_STRUCTURE)
        # print(f'Wrote to {OUTPUT_STRUCTURE}')
        num_added_chains += 1
        print(f'Hit! Number of chains: {num_added_chains + 1}')

if __name__ == '__main__':

    start_time = time.perf_counter()
    main()
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f'Execution time: {elapsed_time:.4f} seconds')