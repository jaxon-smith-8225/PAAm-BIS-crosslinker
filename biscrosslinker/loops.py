import numpy as np
from .geometry import point_segment_distance

def does_cyl_cross_plane(bis_dir, cylinder, end_carbon_pos):
    bis_dir_normalized = bis_dir / np.linalg.norm(bis_dir)
    f_1 = np.dot(bis_dir_normalized, (cylinder.start - end_carbon_pos))
    f_2 = np.dot(bis_dir_normalized, (cylinder.end - end_carbon_pos))
    t = (f_1 * (-1)) / (f_2 - f_1)
    return f_1 * f_2 < 0, t

def coarse_circle_radius(circle_center, new_cylinder):
    dist_from_start = np.linalg.norm(circle_center - new_cylinder.start)
    dist_from_end = np.linalg.norm(circle_center - new_cylinder.end)
    return max(dist_from_end, dist_from_start)

def find_possible_crosses(circle_center, coarse_circle_radius, cylinder_network):
    """
    In:
        all the cylinders to check = network.cylinders
        circle_center for point_segment_distance() = end_carbon.position
        new
    """
    possible_hits = []
    for cylinder in cylinder_network.cylinders:
        min_dist = point_segment_distance(circle_center, cylinder)
        if min_dist > 20. and min_dist < coarse_circle_radius:
            possible_hits.append(cylinder)
    return possible_hits