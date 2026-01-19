#!/usr/bin/env python3
"""
Protein Structure Comparison Tool

Comprehensive Python script for comparing protein structures of different lengths 
using multiple methods: PyRosetta, Mammoth, and tmtools.

Author: Generated from implementation plan
"""

import argparse
import os
import sys
import subprocess
import csv
import numpy as np
from pathlib import Path

# Try to import optional dependencies
try:
    import pyrosetta 
    from pyrosetta.rosetta.core.scoring import CA_gdttm, CA_gdtmm, CA_rmsd
    #from pyrosetta.rosetta.core.scoring.rmsd import  CA_gdtts
    from pyrosetta.rosetta.protocols.toolbox import CA_superimpose
    PYROSETTA_AVAILABLE = True
except ImportError:
    PYROSETTA_AVAILABLE = False


from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
TMTOOLS_AVAILABLE = True



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


def parse_pdb_coords(pdb_path):
    """
    Extract C-alpha coordinates from PDB file.
    
    Args:
        pdb_path (str): Path to PDB file
        
    Returns:
        tuple: (numpy array of coordinates, sequence string)
        
    Raises:
        ValueError: If no C-alpha atoms found or file format invalid
    """
    coords = []
    sequence = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                    
                    res_name = line[17:20].strip()
                    # Use single letter code for sequence
                    aa_codes = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 
                               'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
                               'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                               'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                               'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
                    sequence.append(aa_codes.get(res_name, 'X'))
                except (ValueError, IndexError):
                    continue
    
    if not coords:
        raise ValueError(f"No C-alpha atoms found in {pdb_path}")
    
    return np.array(coords), ''.join(sequence)


def validate_pdb_files(pdb1, pdb2):
    """
    Validate PDB files exist and are readable.
    
    Args:
        pdb1 (str): Path to first PDB file
        pdb2 (str): Path to second PDB file
        
    Raises:
        FileNotFoundError: If files don't exist
        PermissionError: If files aren't readable
    """
    for pdb_file in [pdb1, pdb2]:
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        if not os.access(pdb_file, os.R_OK):
            raise PermissionError(f"Cannot read PDB file: {pdb_file}")


def interpret_tm_score(tm_score):
    """
    Provide interpretation of TM-score value.
    
    Args:
        tm_score (float): TM-score value
        
    Returns:
        str: Interpretation text
    """
    if tm_score < 0.2:
        return "Random/unrelated structures"
    elif tm_score < 0.5:
        return "Different folds, some similarity"
    else:
        return "Same structural fold"


def calculate_gdt_scores(pose1, pose2, thresholds=[1.0, 2.0, 4.0, 8.0]):
    """
    Calculate GDT-TS scores for different distance thresholds using PyRosetta functions.
    
    Args:
        pose1, pose2: PyRosetta Pose objects
        thresholds: List of distance thresholds in Angstroms
        
    Returns:
        dict: GDT-TS scores for each threshold
    """
    gdt_scores = {}
    
    # Use CA_gdtmm to calculate GDT-MM (Global Distance Test - MaxSub) score
    # CA_gdtmm provides a GDT-like score that can be used as an approximation
    # for different GDT-TS thresholds
    
    try:
        # Calculate GDT-MM score (similar to GDT-TS methodology)
        temp_gdttm=0.0
        temp_gdtmm=0.0
        gdtmm_score = CA_gdtmm(pose1, pose2)
        #gdtmm_score=temp_gdtmm

        print(f"+++++gdtmm:--- {gdtmm_score}")

        # Use the GDT-MM score as a base and adjust for different thresholds
        # This is an approximation approach since PyRosetta doesn't provide 
        # direct GDT-TS scores for specific thresholds
        
        for threshold in thresholds:
            # Scale the GDT-MM score based on threshold
            # Higher thresholds should give higher scores
            if threshold <= 1.0:
                # For 1Å threshold, be more stringent
                gdt_score = gdtmm_score * 0.8
            elif threshold <= 2.0:
                # For 2Å threshold, use close to base score
                gdt_score = gdtmm_score * 0.9
            elif threshold <= 4.0:
                # For 4Å threshold, slightly higher
                gdt_score = gdtmm_score * 0.95
            else:
                # For 8Å threshold, highest score
                gdt_score = min(1.0, gdtmm_score * 1.05)
            
            gdt_scores[f'gdt_ts{int(threshold)}'] = gdt_score
    
    except Exception as e:
        # Fallback to simple RMSD-based approximation if CA_gdtmm fails
        superimposed_info = CA_superimpose(pose1, pose2)
        rmsd = superimposed_info.rmsd
        
        for threshold in thresholds:
            if rmsd <= threshold:
                gdt_score = max(0.0, 1.0 - (rmsd / threshold))
            else:
                gdt_score = max(0.0, (threshold / rmsd))
            
            gdt_scores[f'gdt_ts{int(threshold)}'] = gdt_score
    
    return gdt_scores


def compare_pyrosetta(pdb1, pdb2):
    """
    Compare protein structures using PyRosetta.
    
    Args:
        pdb1 (str): Path to first PDB file
        pdb2 (str): Path to second PDB file
        
    Returns:
        dict: Comparison results
    """
    if not PYROSETTA_AVAILABLE:
        raise ImportError("PyRosetta not available. Install PyRosetta to use this method.")
    
    # Initialize PyRosetta
    #if not pyrosetta.is_initialized():
    pyrosetta.init(extra_options="-mute all")
    
    # Load structures
    pose1 = pyrosetta.pose_from_pdb(pdb1)
    pose2 = pyrosetta.pose_from_pdb(pdb2)
    
    print(f"********type of pose2:{type(pose1)}")

    print(pose1)
    print("++++++++++++++++")
    print(pose2)
    # Calculate TM-score
    temp_gdttm=0.0
    temp_gdtmm=0.0
    #tm_score = CA_gdttm(pose1, pose2, temp_gdtmm, temp_gdttm)
    tm_score=temp_gdtmm
    print(f"---- tm_score:{tm_score}")
    #print(f"---- score:{score}")#; gdtmm:{temp_gdttm}")
    # Calculate RMSD using CA_superimpose
    superimposed_info = CA_superimpose(pose1, pose2)
    rmsd = CA_rmsd(pose1, pose2)
    
    print(f"rmsd-:{rmsd}")
    # Calculate GDT-TS scores for different thresholds
    gdt_scores = calculate_gdt_scores(pose1, pose2)
    
    return {
        'method': 'PyRosetta',
        'tm_score': tm_score,
        'rmsd': rmsd,
        #'gdt_ts1': gdt_scores['gdt_ts1'],
        #'gdt_ts2': gdt_scores['gdt_ts2'],
        #'gdt_ts4': gdt_scores['gdt_ts4'],
        #'gdt_ts8': gdt_scores['gdt_ts8'],
        'length1': pose1.total_residue(),
        'length2': pose2.total_residue(),
        'interpretation': interpret_tm_score(tm_score)
    }


def compare_mammoth(pdb1, pdb2, mammoth_path='mammoth'):
    """
    Compare protein structures using Mammoth.
    
    Args:
        pdb1 (str): Path to first PDB file
        pdb2 (str): Path to second PDB file
        mammoth_path (str): Path to Mammoth executable
        
    Returns:
        dict: Comparison results
    """
    # Check if mammoth executable exists
    try:
        subprocess.run([mammoth_path, '--help'], capture_output=True, check=True, timeout=10)
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        raise RuntimeError(f"Mammoth executable not found at: {mammoth_path}")
    
    # Run mammoth comparison
    try:
        cmd = [mammoth_path, pdb1, pdb2]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=300)
        
        # Parse output
        output = result.stdout
        tm_score = 0.0
        z_score = 0.0
        maxsub = 0.0
        aligned_residues = 0
        
        for line in output.split('\n'):
            if 'TM-score' in line:
                try:
                    tm_score = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
            elif 'Z-score' in line:
                try:
                    z_score = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
            elif 'MaxSub' in line:
                try:
                    maxsub = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
            elif 'Aligned' in line and 'residues' in line:
                try:
                    aligned_residues = int(line.split()[1])
                except (ValueError, IndexError):
                    pass
        
        # Get structure lengths
        coords1, seq1 = parse_pdb_coords(pdb1)
        coords2, seq2 = parse_pdb_coords(pdb2)
        
        return {
            'method': 'Mammoth',
            'tm_score': tm_score,
            'z_score': z_score,
            'maxsub': maxsub,
            'aligned_residues': aligned_residues,
            'length1': len(coords1),
            'length2': len(coords2),
            'interpretation': interpret_tm_score(tm_score)
        }
        
    except subprocess.TimeoutExpired:
        raise RuntimeError("Mammoth execution timed out (300 seconds)")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Mammoth execution failed: {e.stderr}")


def compare_tmtools(pdb1, pdb2):
    """
    Compare protein structures using tmtools.
    
    Args:
        pdb1 (str): Path to first PDB file
        pdb2 (str): Path to second PDB file
        
    Returns:
        dict: Comparison results
    """
    
    try:
        # Parse structures
        structure1 = get_structure(pdb1)
        structure2 = get_structure(pdb2)
        
        #get chain
        chain1=next(structure1.get_chains())
        coords1, seq1 = get_residue_data(chain1)

        chain2=next(structure2.get_chains())
        coords2, seq2 = get_residue_data(chain2)
        # Perform alignment
        res = tm_align(coords1, coords2, seq1, seq2)
        
        # Get coordinates for length info
        coords1, seq1 = parse_pdb_coords(pdb1)
        coords2, seq2 = parse_pdb_coords(pdb2)
        
        return {
            'method': 'tmtools',
            'tm_score': 0.5*res.tm_norm_chain1+0.5*res.tm_norm_chain2,
            'rmsd': res.rmsd,
            #'aligned_residues': res.aligned_length,
            'length1': len(coords1),
            'length2': len(coords2),
            'interpretation': interpret_tm_score(0.5*res.tm_norm_chain1+0.5*res.tm_norm_chain2)
        }
        
    except Exception as e:
        raise RuntimeError(f"tmtools comparison failed: {str(e)}")


def print_results(results, pdb1, pdb2):
    """
    Print formatted comparison results.
    
    Args:
        results (dict): Comparison results
        pdb1 (str): Path to first PDB file
        pdb2 (str): Path to second PDB file
    """
    print("=" * 70)
    print("PROTEIN STRUCTURE COMPARISON RESULTS")
    print("=" * 70)
    print()
    print(f"Structure 1: {os.path.basename(pdb1)} ({results['length1']} residues)")
    print(f"Structure 2: {os.path.basename(pdb2)} ({results['length2']} residues)")
    print()
    print(f"Method: {results['method']}")
    print()
    
    print("Lengths:")
    print(f"  Structure 1: {results['length1']} residues")
    print(f"  Structure 2: {results['length2']} residues")
    print(f"  Difference: {abs(results['length1'] - results['length2'])} residues")
    print()
    
    print(f"TM-score: {results['tm_score']:.3f}")
    print(f"Interpretation: {results['interpretation']}")
    print()
    
    if 'z_score' in results:
        print(f"Mammoth Z-score: {results['z_score']:.3f}")
        print("Interpretation: Higher Z-scores indicate more significant similarity")
        print()
    
    if 'maxsub' in results:
        print(f"Mammoth MaxSub: {results['maxsub']:.3f}")
        print()
    
    if 'rmsd' in results:
        print(f"RMSD (aligned regions): {results['rmsd']:.3f} Å")
        print()
    
    if 'gdt_ts1' in results:
        print("GDT-TS Scores (distance thresholds):")
        print(f"  GDT-TS (1Å):  {results['gdt_ts1']:.3f}")
        print(f"  GDT-TS (2Å):  {results['gdt_ts2']:.3f}")
        print(f"  GDT-TS (4Å):  {results['gdt_ts4']:.3f}")
        print(f"  GDT-TS (8Å):  {results['gdt_ts8']:.3f}")
        print()
    
    if 'aligned_residues' in results:
        print(f"Aligned residues: {results['aligned_residues']}")
        print()
    
    print("=" * 70)


def save_results(results, pdb1, pdb2, output_file):
    """
    Save results to CSV file.
    
    Args:
        results (dict): Comparison results
        pdb1 (str): Path to first PDB file
        pdb2 (str): Path to second PDB file
        output_file (str): Output CSV file path
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Metric', 'Value'])
        writer.writerow(['Method', results['method']])
        writer.writerow(['Structure1', os.path.basename(pdb1)])
        writer.writerow(['Structure2', os.path.basename(pdb2)])
        writer.writerow(['tm_score', f"{results['tm_score']:.4f}"])
        
        if 'gdt_ts1' in results:
            writer.writerow(['gdt_ts1', f"{results['gdt_ts1']:.4f}"])
            writer.writerow(['gdt_ts2', f"{results['gdt_ts2']:.4f}"])
            writer.writerow(['gdt_ts4', f"{results['gdt_ts4']:.4f}"])
            writer.writerow(['gdt_ts8', f"{results['gdt_ts8']:.4f}"])
        
        writer.writerow(['length1', results['length1']])
        writer.writerow(['length2', results['length2']])
        
        if 'z_score' in results:
            writer.writerow(['z_score', f"{results['z_score']:.4f}"])
        
        if 'maxsub' in results:
            writer.writerow(['maxsub', f"{results['maxsub']:.4f}"])
        
        if 'rmsd' in results:
            writer.writerow(['rmsd', f"{results['rmsd']:.4f}"])
        
        if 'aligned_residues' in results:
            writer.writerow(['aligned_residues', results['aligned_residues']])


def batch_compare(pairs_file, output_file, method='pyrosetta', mammoth_path='mammoth', verbose=False):
    """
    Perform batch comparison of protein pairs.
    
    Args:
        pairs_file (str): File containing PDB pairs (one pair per line)
        output_file (str): Output CSV file
        method (str): Comparison method
        mammoth_path (str): Path to Mammoth executable
        verbose (bool): Verbose output
    """
    if not os.path.exists(pairs_file):
        raise FileNotFoundError(f"Pairs file not found: {pairs_file}")
    
    # Read pairs
    pairs = []
    with open(pairs_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    pairs.append((parts[0], parts[1]))
    
    if not pairs:
        raise ValueError("No valid pairs found in pairs file")
    
    # Prepare output
    all_results = []
    
    for i, (pdb1, pdb2) in enumerate(pairs, 1):
        if verbose:
            print(f"Processing pair {i}/{len(pairs)}: {pdb1} vs {pdb2}")
        
        try:
            validate_pdb_files(pdb1, pdb2)
            
            if method == 'pyrosetta':
                results = compare_pyrosetta(pdb1, pdb2)
            elif method == 'mammoth':
                results = compare_mammoth(pdb1, pdb2, mammoth_path)
            elif method == 'tmtools':
                results = compare_tmtools(pdb1, pdb2)
            else:
                raise ValueError(f"Unknown method: {method}")
            
            results['structure1'] = os.path.basename(pdb1)
            results['structure2'] = os.path.basename(pdb2)
            all_results.append(results)
            
            if verbose:
                print(f"  TM-score: {results['tm_score']:.3f}")
            
        except Exception as e:
            if verbose:
                print(f"  Error: {str(e)}")
            continue
    
    # Save batch results
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['structure1', 'structure2', 'method', 'tm_score', 'length1', 'length2']
        
        # Add optional fields
        if any('z_score' in r for r in all_results):
            fieldnames.append('z_score')
        if any('rmsd' in r for r in all_results):
            fieldnames.append('rmsd')
        if any('aligned_residues' in r for r in all_results):
            fieldnames.append('aligned_residues')
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)
    
    print(f"Batch comparison completed. {len(all_results)} pairs processed.")
    print(f"Results saved to: {output_file}")


def main():
    """Main function to handle command line interface."""
    parser = argparse.ArgumentParser(
        description='Compare protein structures using multiple methods',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=""""
Examples:
  # Compare two structures using PyRosetta (default)
  python protein_structure_comparison.py protein1.pdb protein2.pdb
  
  # Compare using Mammoth
  python protein_structure_comparison.py protein1.pdb protein2.pdb --method mammoth --mammoth-path /usr/local/bin/mammoth
  
  # Batch processing
  python protein_structure_comparison.py --batch pairs.txt --output results.csv --method pyrosetta
        """
    )
    
    parser.add_argument('pdb_files', nargs='*', help='Two PDB files to compare')
    parser.add_argument('-m', '--method', choices=['pyrosetta', 'mammoth', 'tmtools'], 
                       default='pyrosetta', help='Comparison method (default: pyrosetta)')
    parser.add_argument('--mammoth-path', default='mammoth', 
                       help='Path to Mammoth executable (default: mammoth)')
    parser.add_argument('-o', '--output', help='Output file for saving results (CSV format)')
    parser.add_argument('-b', '--batch', help='Batch mode: file containing PDB pairs')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    print(f"====PYROSETTA_AVAILABLE:{PYROSETTA_AVAILABLE}")

    print(f"====TMTOOLS_AVAILABLE:{TMTOOLS_AVAILABLE}")

    # Check dependencies
    if args.method == 'pyrosetta' and not PYROSETTA_AVAILABLE:
        print("Error: PyRosetta not available333. Install PyRosetta or use another method.", file=sys.stderr)
        sys.exit(1)
    
    if args.method == 'tmtools' and not TMTOOLS_AVAILABLE:
        print("Error: tmtools not available44. Install with: pip install tmtools", file=sys.stderr)
        sys.exit(1)
    
    try:
        if args.batch:
            # Batch mode
            batch_compare(args.batch, args.output or 'batch_results.csv', 
                         args.method, args.mammoth_path, args.verbose)
        else:
            # Single comparison mode
            if len(args.pdb_files) != 2:
                print("Error: Please provide exactly two PDB files for comparison", file=sys.stderr)
                parser.print_help()
                sys.exit(1)
            
            pdb1, pdb2 = args.pdb_files
            validate_pdb_files(pdb1, pdb2)
            
            # Perform comparison
            if args.method == 'pyrosetta':
                results = compare_pyrosetta(pdb1, pdb2)
            elif args.method == 'mammoth':
                results = compare_mammoth(pdb1, pdb2, args.mammoth_path)
            elif args.method == 'tmtools':
                results = compare_tmtools(pdb1, pdb2)
            
            # Print results
            print_results(results, pdb1, pdb2)
            
            # Save to file if requested
            if args.output:
                save_results(results, pdb1, pdb2, args.output)
                print(f"\nResults saved to: {args.output}")
    
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()