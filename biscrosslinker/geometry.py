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
    Credit: Arfana

    Calculate the minimum distance between two finite line segments in 3D space.

    Uses Ericson's method (Real-Time Collision Detection, Ch. 5) which correctly
    clamps to segment extents, unlike the infinite-line formula.

    Args:
        a1, a2: Start and end points of first segment
        b1, b2: Start and end points of second segment

    Returns:
        Minimum distance between the segments
    """
    a1, a2, b1, b2 = map(np.array, (a1, a2, b1, b2))
    da = a2 - a1   # Direction of segment A
    db = b2 - b1   # Direction of segment B
    r  = a1 - b1

    len_a_sq = np.dot(da, da)  # Squared length of A
    len_b_sq = np.dot(db, db)  # Squared length of B
    f = np.dot(db, r)

    # Handle degenerate cases (zero-length segments)
    if len_a_sq < 1e-10 and len_b_sq < 1e-10:
        # Both segments are points
        return np.linalg.norm(r)

    if len_a_sq < 1e-10:
        # Segment A is a point - clamp s=0, find closest point on B
        s = 0.0
        t = np.clip(f / len_b_sq, 0.0, 1.0)
    else:
        c = np.dot(da, r)
        if len_b_sq < 1e-10:
            # Segment B is a point - clamp t=0, find closest point on A
            t = 0.0
            s = np.clip(-c / len_a_sq, 0.0, 1.0)
        else:
            # General non-degenerate case
            b_dot = np.dot(da, db)
            denom = len_a_sq * len_b_sq - b_dot * b_dot  # Always >= 0

            if denom > 1e-10:
                # Segments are not parallel: compute closest point and clamp to [0,1]
                s = np.clip((b_dot * f - c * len_b_sq) / denom, 0.0, 1.0)
            else:
                # Segments are parallel: pick s=0 arbitrarily
                s = 0.0

            # Compute t for the clamped s, then clamp t too
            t = (b_dot * s + f) / len_b_sq
            if t < 0.0:
                t = 0.0
                s = np.clip(-c / len_a_sq, 0.0, 1.0)
            elif t > 1.0:
                t = 1.0
                s = np.clip((b_dot - c) / len_a_sq, 0.0, 1.0)

    closest_a = a1 + s * da
    closest_b = b1 + t * db
    return np.linalg.norm(closest_a - closest_b)

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


def rotate_by_angle(universe_AtomGroup, axis, angle, point=[0, 0, 0]):
    """
    Rotate universe by specified angle around axis
    
    Args:
        universe: MDAnalysis Universe to rotate
        axis: Rotation axis vector
        angle: Rotation angle in degrees
    """
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    universe_AtomGroup.rotateby(angle, axis, point=point)


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
    
    universe = universe.select_atoms('all')
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
    universe = universe.select_atoms('all')
    rotate_by_angle(universe, axis, rand_rotation_angle)
    return rand_rotation_angle

def r12(point1, point2):
    dist = 0
    for i in range(len(point1)):
        dist += (point2[i] - point1[i]) ** 2
    return np.sqrt(dist)
