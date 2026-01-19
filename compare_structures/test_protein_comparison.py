#!/usr/bin/env python3
"""
Test script for protein structure comparison tool.

Creates sample PDB files and tests the comparison functionality.
"""

import os
import sys
import numpy as np
from protein_structure_comparison import (
    parse_pdb_coords, validate_pdb_files, interpret_tm_score,
    compare_tmtools, TMTOOLS_AVAILABLE
)


def create_sample_pdb(filename, num_residues, pattern='helix'):
    """
    Create a sample PDB file for testing.
    
    Args:
        filename (str): Output PDB filename
        num_residues (int): Number of residues to create
        pattern (str): Structure pattern ('helix', 'sheet', 'random')
    """
    with open(filename, 'w') as f:
        for i in range(num_residues):
            res_num = i + 1
            res_name = 'ALA'  # Use alanine for simplicity
            
            if pattern == 'helix':
                # Alpha helix coordinates
                angle = i * 100  # degrees
                x = 10 * np.cos(np.radians(angle))
                y = 10 * np.sin(np.radians(angle))
                z = i * 1.5
            elif pattern == 'sheet':
                # Beta sheet coordinates
                x = i * 3.5
                y = 0
                z = 5 * np.sin(i * 0.5)
            else:  # random
                x = np.random.uniform(-10, 10)
                y = np.random.uniform(-10, 10)
                z = np.random.uniform(-10, 10)
            
            # Write ATOM record for C-alpha
            f.write(f"ATOM  {5*i+1:5d}  CA  {res_name} A{res_num:4d}    "
                   f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n")


def test_pdb_parsing():
    """Test PDB coordinate parsing."""
    print("Testing PDB parsing...")
    
    # Create test PDB
    test_pdb = "test_structure.pdb"
    create_sample_pdb(test_pdb, 10, 'helix')
    
    try:
        coords, sequence = parse_pdb_coords(test_pdb)
        assert len(coords) == 10, f"Expected 10 residues, got {len(coords)}"
        assert len(sequence) == 10, f"Expected sequence length 10, got {len(sequence)}"
        print("  ✓ PDB parsing works correctly")
    except Exception as e:
        print(f"  ✗ PDB parsing failed: {e}")
    finally:
        if os.path.exists(test_pdb):
            os.remove(test_pdb)


def test_file_validation():
    """Test file validation."""
    print("Testing file validation...")
    
    # Create test files
    create_sample_pdb("test1.pdb", 5, 'helix')
    create_sample_pdb("test2.pdb", 7, 'sheet')
    
    try:
        validate_pdb_files("test1.pdb", "test2.pdb")
        print("  ✓ File validation works correctly")
    except Exception as e:
        print(f"  ✗ File validation failed: {e}")
    finally:
        for f in ["test1.pdb", "test2.pdb"]:
            if os.path.exists(f):
                os.remove(f)


def test_tm_score_interpretation():
    """Test TM-score interpretation."""
    print("Testing TM-score interpretation...")
    
    test_cases = [
        (0.1, "Random/unrelated structures"),
        (0.3, "Different folds, some similarity"),
        (0.7, "Same structural fold")
    ]
    
    for tm_score, expected in test_cases:
        result = interpret_tm_score(tm_score)
        assert result == expected, f"Expected '{expected}', got '{result}'"
    
    print("  ✓ TM-score interpretation works correctly")


def test_tmtools_comparison():
    """Test tmtools comparison if available."""
    if not TMTOOLS_AVAILABLE:
        print("Skipping tmtools test (not available)")
        return
    
    print("Testing tmtools comparison...")
    
    # Create test structures
    create_sample_pdb("test_helix.pdb", 15, 'helix')
    create_sample_pdb("test_sheet.pdb", 15, 'sheet')
    
    try:
        results = compare_tmtools("test_helix.pdb", "test_sheet.pdb")
        
        # Check required fields
        required_fields = ['method', 'tm_score', 'length1', 'length2', 'interpretation']
        for field in required_fields:
            assert field in results, f"Missing required field: {field}"
        
        assert results['method'] == 'tmtools'
        assert results['length1'] == 15
        assert results['length2'] == 15
        assert 0 <= results['tm_score'] <= 1
        
        print("  ✓ tmtools comparison works correctly")
        
    except Exception as e:
        print(f"  ✗ tmtools comparison failed: {e}")
    finally:
        for f in ["test_helix.pdb", "test_sheet.pdb"]:
            if os.path.exists(f):
                os.remove(f)


def test_command_line_interface():
    """Test command line interface."""
    print("Testing command line interface...")
    
    # Create test structures
    create_sample_pdb("cli_test1.pdb", 12, 'helix')
    create_sample_pdb("cli_test2.pdb", 12, 'helix')
    
    try:
        # Test help
        result = os.system("python3 protein_structure_comparison.py --help > /dev/null 2>&1")
        if result == 0:
            print("  ✓ Command line interface works correctly")
        else:
            print("  ✗ Command line interface failed")
    except Exception as e:
        print(f"  ✗ Command line interface test failed: {e}")
    finally:
        for f in ["cli_test1.pdb", "cli_test2.pdb"]:
            if os.path.exists(f):
                os.remove(f)


def main():
    """Run all tests."""
    print("Running protein structure comparison tests...")
    print("=" * 50)
    
    test_pdb_parsing()
    test_file_validation()
    test_tm_score_interpretation()
    test_tmtools_comparison()
    test_command_line_interface()
    
    print("=" * 50)
    print("Testing completed!")


if __name__ == '__main__':
    main()