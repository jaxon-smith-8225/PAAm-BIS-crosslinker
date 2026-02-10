#!/usr/bin/env python3
"""
Main script for running BIS crosslinker simulation

This script builds a crosslinked polymer network by iteratively adding
PAAm chains connected via BIS crosslinkers.
"""

import random
import time
import warnings
import MDAnalysis as mda

from biscrosslinker import (
    PAAmTemplate,
    CrosslinkNetwork,
    replace_with_bis,
    connect_PAAm
)
from biscrosslinker.config import (
    OUTPUT_STRUCTURE,
    MAX_NUM_CHAINS,
    MIN_CROSSLINK_DISTANCE
)

# Suppress MDAnalysis warnings
warnings.filterwarnings('ignore')


def main():
    """
    Main simulation loop for building crosslinked network
    """
    # Create template class
    print("Loading PAAm template...")
    template = PAAmTemplate("PAAm25mer.pdb")
    
    # Initialize network with first chain
    print("Initializing network with first chain...")
    initial_structure = mda.Universe("PAAm25mer.pdb", context='default', 
                                     to_guess=['elements', 'bonds'])
    network = CrosslinkNetwork(initial_structure)
    print(f'Initial chain added. Total chains: {network.num_chains()}')
    
    # Main loop - keep adding chains
    while network.num_chains() < MAX_NUM_CHAINS:
        
        # Use network method to find valid sites
        valid_sites = network.find_valid_sites(min_distance=MIN_CROSSLINK_DISTANCE)
        
        if not valid_sites:
            print("No valid sites found - all too close to existing crosslinks")
            break
        
        # Pick a random valid site
        linking_site = random.choice(valid_sites)
        backbone_carbon_pos = linking_site['backbone_carbon'][0].position.copy()

        # Replace reactive site with BIS
        first_modified_univ, bis_univ = replace_with_bis(linking_site, network.structure)
        
        try:
            # Pass network and template instead of individual lists
            final_modified_univ, new_cylinder = connect_PAAm(
                bis_univ, first_modified_univ, network, template
            )
        except RuntimeError as e:
            print(f"{e}: Trying again")
            continue
        
        # Use network method to add chain w/ all metadata
        network.add_chain(final_modified_univ, new_cylinder, backbone_carbon_pos)
        
        print(f'Hit! Total chains: {network.num_chains()}')
    
    # Write final structure
    print(f"\nWriting final structure to {OUTPUT_STRUCTURE}...")
    network.structure.select_atoms('all').write(OUTPUT_STRUCTURE)
    print(f"Final structure saved with {network.num_chains()} total chains!")


if __name__ == '__main__':
    start_time = time.perf_counter()
    main()
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f'Execution time: {elapsed_time:.4f} seconds')