import numpy as np
from .geometry import point_segment_distance, r12

def find_point_of_intersection(cylinder, t):
    return cylinder.start + t*(cylinder.axis)

def does_cyl_cross_plane(bis_dir, cylinder, end_carbon_pos):
    bis_dir_normalized = bis_dir / np.linalg.norm(bis_dir)
    f_1 = np.dot(bis_dir_normalized, (cylinder.start - end_carbon_pos))
    f_2 = np.dot(bis_dir_normalized, (cylinder.end - end_carbon_pos))
    t = (f_1 * (-1)) / (f_2 - f_1)
    return f_1 * f_2 < 0, find_point_of_intersection(cylinder, t)

def coarse_circle_radius(circle_center, new_cylinder):
    dist_from_start = np.linalg.norm(circle_center - new_cylinder.start)
    dist_from_end = np.linalg.norm(circle_center - new_cylinder.end)
    dist_from_end > dist_from_start

    return max(dist_from_end, dist_from_start)

def find_possible_crosses(circle_center, coarse_circle_radius, cylinder_network):
    """
    In:
        all the cylinders to check = network.cylinders
        circle_center for point_segment_distance() = end_carbon.position
    """
    possible_hits = []
    for cylinder in cylinder_network.cylinders:
        min_dist = point_segment_distance(circle_center, cylinder)
        if min_dist > 20. and min_dist < coarse_circle_radius:
            possible_hits.append(cylinder)
    return possible_hits

def fine_circle_radius(cylinder, status, hit_distance):

    return

def find_best_existing_site(intersection_point, pivot_point, universe, max_search_radius=10.0):
    """
    Find the best reactive site on an existing chain near the plane intersection point
    
    Args:
        intersection_point: Where the cylinder crosses the rotation plane
        pivot_point: The BIS attachment point (end_carbon.position)
        universe: MDAnalysis universe (or atom selection) to search for reactive sites
        max_search_radius: Maximum distance from intersection point to consider
        
    Returns:
        Tuple of (best_site_dict, distance_to_pivot) or (None, None) if no valid site found
    """
    from .chemistry import scan_chain
    
    # Search for sites near the intersection point
    sites_near_intersection = scan_chain(universe.select_atoms(
        f'point {intersection_point[0]:.1f} {intersection_point[1]:.1f} {intersection_point[2]:.1f} {max_search_radius:.1f}'
    ))
    
    if not sites_near_intersection:
        return None, None
    
    # Find the site closest to the intersection point
    best_site = None
    best_distance_to_intersection = float('inf')
    best_distance_to_pivot = None
    
    for site in sites_near_intersection:
        site_pos = site['central_carbon'].position
        dist_to_intersection = r12(site_pos, intersection_point)
        
        if dist_to_intersection < best_distance_to_intersection:
            best_distance_to_intersection = dist_to_intersection
            best_site = site
            best_distance_to_pivot = r12(site_pos, pivot_point)
    
    return best_site, best_distance_to_pivot


def find_best_new_chain_site(new_chain_universe, pivot_point, max_distance_from_pivot, BIS_LENGTH=5.4):
    """
    Find the best reactive site on the new chain for forming a loop
    
    Want a site where:
    - Distance from pivot (side c) < max_distance_from_pivot (side b)
    - Distance from pivot > BIS_LENGTH (side a)
    - Among valid sites, pick the one with largest distance from pivot (tightest triangle)
    
    Args:
        new_chain_universe: MDAnalysis universe of the new chain
        pivot_point: The BIS attachment point (end_carbon.position)
        max_distance_from_pivot: Distance from pivot to existing chain's reactive site (side b)
        BIS_LENGTH: Length of BIS molecule (side a), default 5.4 Angstroms
        
    Returns:
        Tuple of (best_site_dict, distance_to_pivot) or (None, None) if no valid site found
    """
    from .chemistry import scan_chain
    
    # Scan all reactive sites on new chain
    all_sites = scan_chain(new_chain_universe)
    
    # Skip first and last sites
    valid_sites = all_sites[1:-1]
    
    if not valid_sites:
        return None, None
    
    # Find sites that satisfy triangle inequality
    best_site = None
    best_distance = 0  # Want the largest valid distance
    
    for site in valid_sites:
        site_pos = site['central_carbon'].position
        dist_to_pivot = r12(site_pos, pivot_point)
        
        # Triangle inequality constraints:
        # 1. side c < side b (this site closer to pivot than existing site)
        # 2. side c > side a (this site farther from pivot than BIS length)
        # 3. side b < side c + side a (automatically satisfied if 1 and 2 are true)
        
        if BIS_LENGTH < dist_to_pivot < max_distance_from_pivot:
            # Valid site - keep if it's farther from pivot than previous best
            if dist_to_pivot > best_distance:
                best_distance = dist_to_pivot
                best_site = site
    
    if best_site is None:
        return None, None
    
    return best_site, best_distance


def check_bis_orientation(new_chain_site, existing_chain_site, max_angle_degrees):
    """
    Check that the two reactive sites are pointing roughly toward each other,
    so a BIS molecule can connect them without straining against the chain geometry

    For each site, the reactive direction vector (backbone carbon -> reactive carbon)
    should point roughly toward the other site. Check this by measuring the angle
    between each site's reactive direction and the vector pointing to the other site.
    Both angles must be within max_angle_degrees

    Args:
        new_chain_site: Unpacked site dict from unpack_linking_site() on the new chain
        existing_chain_site: Unpacked site dict from unpack_linking_site() on the existing chain
        max_angle_degrees: Maximum allowed angle between reactive direction and
                           site-to-site vector (degrees)

    Returns:
        True if both sites are oriented acceptably, False otherwise
    """
    from .chemistry import calculate_reactive_direction
    from .geometry import angle_between_vectors

    new_site_pos      = new_chain_site['reactive_carbon'].position
    existing_site_pos = existing_chain_site['reactive_carbon'].position

    # Vector from new site pointing toward existing site, and vice versa
    vec_new_to_existing = existing_site_pos - new_site_pos
    vec_existing_to_new = new_site_pos - existing_site_pos

    # Reactive direction for each site (backbone -> reactive carbon)
    new_reactive_dir      = calculate_reactive_direction(new_chain_site)
    existing_reactive_dir = calculate_reactive_direction(existing_chain_site)

    # Angle between each site's reactive direction and the vector toward the other site
    angle_new      = angle_between_vectors(new_reactive_dir,      vec_new_to_existing)
    angle_existing = angle_between_vectors(existing_reactive_dir, vec_existing_to_new)

    return angle_new <= max_angle_degrees and angle_existing <= max_angle_degrees


def calculate_loop_rotation_angle(pivot_point, new_chain_site_pos, existing_chain_site_pos, BIS_LENGTH=5.4):
    """
    Calculate the angle to rotate the new chain to form a loop connection
    
    Uses law of cosines on the triangle formed by:
    - Pivot point (end_carbon, the BIS attachment)
    - New chain's reactive site
    - Existing chain's reactive site
    
    Triangle sides:
    - side a: BIS_LENGTH (distance between the two reactive carbons after connection)
    - side b: distance from pivot to existing chain's site (longest side)
    - side c: distance from pivot to new chain's site
    
    Want the angle at the pivot point (angle between sides b and c).
    
    Law of cosines: a**2 = b**2 + c**2 - 2bc·cos(angle)
    Rearranged: angle = arccos((b**2 + c**2 - a**2) / (2bc))
    
    Args:
        pivot_point: Position of the BIS attachment point (end_carbon.position)
        new_chain_site_pos: Position of reactive carbon on new chain
        existing_chain_site_pos: Position of reactive carbon on existing chain
        BIS_LENGTH: Length of BIS molecule, default 5.4 Angstroms
        
    Returns:
        Rotation angle in degrees, or None if invalid triangle
    """
    # Calculate the three sides of the triangle
    side_b = r12(pivot_point, existing_chain_site_pos)  # pivot to existing site (longest)
    side_c = r12(pivot_point, new_chain_site_pos)        # pivot to new site
    side_a = BIS_LENGTH                                   # BIS length (what connects the two sites)
    
    # Verify triangle inequality (should already be satisfied, but double-check)
    if not (side_a < side_b + side_c and side_b < side_a + side_c and side_c < side_a + side_b):
        # print(f"Warning: Invalid triangle - sides {side_a:.2f}, {side_b:.2f}, {side_c:.2f}")
        return None
    
    # Law of cosines: cos(angle) = (b**2 + c**2 - a**2) / (2bc)
    cos_angle = (side_b**2 + side_c**2 - side_a**2) / (2 * side_b * side_c)
    
    # Clamp to [-1, 1] to handle floating point errors
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    # Calculate angle in radians, then convert to degrees
    angle_radians = np.arccos(cos_angle)
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees