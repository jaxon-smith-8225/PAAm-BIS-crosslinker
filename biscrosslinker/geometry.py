"""
Geometric utilities for collision detection and molecular alignment
"""

import numpy as np
from .config import CHAIN_RADIUS


class ChainCylinder:
    """Represents a cylindrical bounding volume for a PAAm chain"""
    
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


def segment_segment_distance(a1, a2, b1, b2):
    """
    Calculate the minimum distance between two line segments in 3D space
    
    Credit: Arfana
    
    Args:
        a1, a2: Start and end points of first segment
        b1, b2: Start and end points of second segment
        
    Returns:
        Minimum distance between the segments
    """
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

def point_segment_distance(circle_center, cylinder_to_check):
    """
    In:
        center of circle = bis_univ.select_atoms('index 5)[0] = end_carbon
        all the cylinders to check = network
        the cylinder just added = new_cylinder

    FIND LATER:
        radius of circle = distance from pivot to longest end of chain 
        = site['backbone_carbon'].position - new_cylinder.start/end

    Out:
        set of cylinders that are close enough to worry about
    """

    """
    Equation for shortest distance btwn point and line segment

    |(B-A) cross (P - A)| / |B-A|
    || (new_cylinder.axis) x (end_carbon - new_cylinder.start) || / || (new_cylinder.length) ||
    top_vector = np.cross(new_cylinder.axis, (end_carbon - new_cylinder.start))
    numerator = np.linalg.norm(top_vector)
    denominator = new_cylinder.length
    min_dist = numerator / denominator
    """
    top_vector = np.linalg.cross(cylinder_to_check.axis, (circle_center - cylinder_to_check.start))
    numerator = np.linalg.norm(top_vector)
    denominator = cylinder_to_check.length
    min_dist_inclusive = numerator / denominator
    check_start = np.linalg.norm(circle_center - cylinder_to_check.start)
    check_end = np.linalg.norm(circle_center - cylinder_to_check.end)
    return min(min_dist_inclusive, check_start, check_end)


def cylinders_collide(cyl1, cyl2):
    """
    Check if two ChainCylinder objects overlap
    
    Args:
        cyl1, cyl2: ChainCylinder objects to check
        
    Returns:
        True if cylinders collide, False otherwise
    """
    dist = segment_segment_distance(
        cyl1.start, cyl1.end,
        cyl2.start, cyl2.end
    )
    return dist < (cyl1.radius + cyl2.radius)


def angle_between_vectors(v1, v2):
    """
    Calculate angle (in degrees) between two vectors
    
    Args:
        v1, v2: Input vectors
        
    Returns:
        Angle in degrees
    """
    dot_prod = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    cosine_angle = dot_prod / (norm_v1 * norm_v2)
    clipped_cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    angle_radians = np.arccos(clipped_cosine_angle)
    angle_degrees = angle_radians * 180 / np.pi
    return angle_degrees


def rotate_by_angle(universe, axis, angle):
    """
    Rotate universe by specified angle around axis
    
    Args:
        universe: MDAnalysis Universe to rotate
        axis: Rotation axis vector
        angle: Rotation angle in degrees
    """
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    universe.atoms.rotateby(angle, axis, point=[0, 0, 0])


def align_molecule(universe, source_vector, target_vector):
    """
    Rotate universe to align source_vector with target_vector
    
    Args:
        universe: mda.Universe to rotate
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


def apply_random_rotation(universe, axis, min_angle, max_angle):
    """
    Apply random rotation about an axis
    
    Args:
        universe: mda.Universe to rotate
        axis: Rotation axis vector
        min_angle: Minimum rotation angle (degrees)
        max_angle: Maximum rotation angle (degrees)
        
    Returns:
        The rotation angle that was applied
    """
    import random
    rand_rotation_angle = random.randint(min_angle, max_angle)
    rotate_by_angle(universe, axis, rand_rotation_angle)
    return rand_rotation_angle