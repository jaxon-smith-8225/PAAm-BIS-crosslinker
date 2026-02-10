import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np
import warnings
warnings.filterwarnings('ignore')

u = mda.Universe("BIS.pdb", context='default', to_guess=['elements', 'bonds'])
atom1 = u.select_atoms('index 5')
atom2 = u.select_atoms('index 6')
current_vector = atom2.positions[0] - atom1.positions[0]
target_length = 6
reference_point = atom1.positions[0]

def stretch_universe(universe, vector_to_scale, reference_point, target_length):
    current_length = np.linalg.norm(vector_to_scale)
    stretch_factor = target_length / current_length
    universe.atoms.positions = reference_point + (universe.atoms.positions - reference_point) * stretch_factor
    return

stretch_universe(u, current_vector, reference_point, target_length)
u.select_atoms('all').write("BIS_stretched.pdb")
