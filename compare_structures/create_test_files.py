#!/usr/bin/env python3
"""
Create sample PDB files for testing the comparison tool.
"""

import numpy as np
from protein_structure_comparison import create_sample_pdb

# Create test structures
create_sample_pdb("test_helix1.pdb", 20, 'helix')
create_sample_pdb("test_helix2.pdb", 25, 'helix')  # Different length
create_sample_pdb("test_sheet1.pdb", 20, 'sheet')
create_sample_pdb("test_random1.pdb", 20, 'random')

print("Created test PDB files:")
print("- test_helix1.pdb (20 residues, helix)")
print("- test_helix2.pdb (25 residues, helix)")
print("- test_sheet1.pdb (20 residues, sheet)")
print("- test_random1.pdb (20 residues, random)")