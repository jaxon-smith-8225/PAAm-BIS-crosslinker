Overview: biscrosslinker.py generates crosslinked polyacrylamide (PAAm) hydrogel structures by connecting polymer chains with BIS crosslinker molecules while avoiding steric clashes

Features:
- Identifies reactive sites on PAAm chains
- Inserts BIS crosslinker molecules at selected sites
- Connects new PAAm chains to crosslinkers w/ proper orientation
- Uses cylindrical bounding volumes to prevent chain collisions
- Applies random rotations for more realistic chain conformations

Prereqs:
- Python 3.7+
- Required packages: MDAnalysis, numpy

Required Input Files:
- PAAm25mer.pdb - Template PAAm polymer chain (25 monomers)
- BIS.pdb - Modified BIS crosslinker molecule structure

Usage: python biscrosslinker.py
- The script will generate an output file containing the crosslinked structure

Configuration:
- OUTPUT_STRUCTURE: Output filename (defualt: "output.pdb")
- MAX_NUM_CHAINS: Maximum number of chains to add (default: 20)
- CHAIN_RADIUS: Cylinder radius for collision detection (default: 3.7 Angstroms)
- MIN_RAND_ANGLE/MAX_RAND_ANGLE: Random rotation range (degrees)
- RAND_CHAIN_ROTATION: Enable/disable random rotations
- MAX_ROTATION_ATTEMPTS: Max attempts to avoid collisions

Output:
- PDB file w/ crosslinked polymer network (TO DO: loops and misc. network topology)
- Console log showing chain addition progress
- Final chain count and execution time

Algorithm:
1. Load initial PAAm chain
2. Scan for reactive carbon sites (-CONH_2 groups)
3. Select valid site (not too close to existing groups)
4. Replace reactive group w/ BIS crosslinker
5. Attach new PAAm chain to BIS w/ collision checking
6. Repeat until reaching targe chain count or no valid sites remain

Notes:
- Reactive sites are identified by specific bonding patterns (C-C, C-O, C-N)
- Excludes first and last monomers to maintain chain integrity
- Maintains minimum distance (10 Angstroms) between crosslinking sites
- Uses geometric alignment for bond formation
