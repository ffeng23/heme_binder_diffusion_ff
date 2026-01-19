#!/usr/bin/env python3
"""
Simple demonstration of the protein structure comparison tool.
Creates sample PDB files and shows basic functionality without external dependencies.
"""

import numpy as np
import sys
import os

# Add current directory to path
sys.path.append('.')

from protein_structure_comparison import (
    parse_pdb_coords, validate_pdb_files, interpret_tm_score,
    create_sample_pdb
)


def simple_rmsd_comparison(coords1, coords2):
    """
    Simple RMSD calculation for demonstration.
    
    Args:
        coords1 (numpy array): First structure coordinates
        coords2 (numpy array): Second structure coordinates
        
    Returns:
        float: RMSD value
    """
    # Make coordinates the same length by truncating to shorter
    min_len = min(len(coords1), len(coords2))
    coords1_trim = coords1[:min_len]
    coords2_trim = coords2[:min_len]
    
    # Calculate RMSD
    diff = coords1_trim - coords2_trim
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd


def simple_tm_score_estimate(coords1, coords2, rmsd):
    """
    Simple TM-score estimation based on RMSD and length.
    
    Args:
        coords1 (numpy array): First structure coordinates
        coords2 (numpy array): Second structure coordinates
        rmsd (float): RMSD value
        
    Returns:
        float: Estimated TM-score
    """
    # Simple approximation: TM-score ≈ 1 / (1 + (RMSD/L)^2)
    # where L is the length of the shorter structure
    min_len = min(len(coords1), len(coords2))
    
    if min_len == 0:
        return 0.0
    
    tm_score = 1.0 / (1.0 + (rmsd / min_len)**2)
    return tm_score


def demo_comparison():
    """Demonstrate the comparison functionality."""
    print("Protein Structure Comparison Tool Demo")
    print("=" * 50)
    
    # Create test structures if they don't exist
    test_files = ['test_helix1.pdb', 'test_helix2.pdb', 'test_sheet1.pdb', 'test_random1.pdb']
    
    for filename in test_files:
        if not os.path.exists(filename):
            pattern = 'helix' if 'helix' in filename else ('sheet' if 'sheet' in filename else 'random')
            num_residues = 25 if 'helix2' in filename else 20
            create_sample_pdb(filename, num_residues, pattern)
    
    print("Test files created/verified:")
    for f in test_files:
        coords, seq = parse_pdb_coords(f)
        print(f"  {f}: {len(coords)} residues")
    
    print("\nComparison Results:")
    print("-" * 50)
    
    # Test comparisons
    comparisons = [
        ('test_helix1.pdb', 'test_helix2.pdb', 'Similar structures (helix vs helix)'),
        ('test_helix1.pdb', 'test_sheet1.pdb', 'Different structures (helix vs sheet)'),
        ('test_helix1.pdb', 'test_random1.pdb', 'Random vs structured'),
    ]
    
    for pdb1, pdb2, description in comparisons:
        print(f"\n{description}")
        print(f"Comparing {pdb1} vs {pdb2}")
        
        try:
            validate_pdb_files(pdb1, pdb2)
            
            # Parse coordinates
            coords1, seq1 = parse_pdb_coords(pdb1)
            coords2, seq2 = parse_pdb_coords(pdb2)
            
            # Calculate metrics
            rmsd = simple_rmsd_comparison(coords1, coords2)
            tm_score = simple_tm_score_estimate(coords1, coords2, rmsd)
            
            print(f"  Length 1: {len(coords1)} residues")
            print(f"  Length 2: {len(coords2)} residues")
            print(f"  Length difference: {abs(len(coords1) - len(coords2))} residues")
            print(f"  RMSD: {rmsd:.3f} Å")
            print(f"  Estimated TM-score: {tm_score:.3f}")
            print(f"  Interpretation: {interpret_tm_score(tm_score)}")
            
        except Exception as e:
            print(f"  Error: {e}")
    
    print("\n" + "=" * 50)
    print("Demo completed!")
    print("\nNote: This is a simplified demonstration.")
    print("For accurate TM-score calculations, install:")
    print("  - PyRosetta: https://www.pyrosetta.org/")
    print("  - tmtools: pip install tmtools")
    print("  - Mammoth: SBGrid installation or compile from source")


if __name__ == '__main__':
    demo_comparison()