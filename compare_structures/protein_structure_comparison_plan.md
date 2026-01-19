# Protein Structure Comparison Tool - Implementation Plan

## Overview
Comprehensive Python script for comparing protein structures of different lengths using multiple methods.

## Three Comparison Methods

### Method 1: PyRosetta Native TM-score (Primary)
**Advantages:** Pure Python, fast, similar methodology to Mammoth

**Key Functions from `rosetta.core.scoring`:**
- `CA_gdttm()` - Primary TM-score
- `CA_gdtts1/2/4/8()` - GDT-TS scores with distance thresholds
- `calculate_rmsd_superimposed()` - RMSD (if needed)

**TM-score Interpretation:**
- < 0.2: Random/unrelated structures
- 0.2-0.5: Different folds, some similarity
- > 0.5: Same structural fold
- 1.0: Perfect match

### Method 2: Mammoth Subprocess Function
**Implementation:** Call Mammoth executable via Python subprocess

**Features:**
- Extract TM-score (0-1 scale)
- Extract Z-score (statistical significance)
- Extract aligned residue counts
- Parse sequence alignment information

**Mammoth Installation:**
- Available via SBGrid: `sbgrid-cli install mammoth-mult`
- Or manual compilation from source (C++ code)
- No official Python wrapper exists

### Method 3: tmtools (TM-align) - Python Native
**Alternative if Mammoth installation fails**

**Features:**
- Python native implementation of TM-align algorithm
- Similar methodology to Mammoth
- Easy installation: `pip install tmtools`

## Complete Script Structure

### File Structure
```
protein_structure_comparison.py
├── Import statements (PyRosetta, tmtools, subprocess, numpy)
├── Helper functions
│   ├── parse_pdb_coords()
│   ├── validate_pdb_files()
│   └── interpret_tm_score()
├── Comparison methods
│   ├── compare_pyrosetta()
│   ├── compare_mammoth()
│   └── compare_tmtools()
├── Output functions
│   ├── print_results()
│   ├── save_results()
│   └── batch_compare()
└── Main execution
    └── main()
```

### Key Functions

**parse_pdb_coords(pdb_path)**
- Extract C-alpha coordinates from PDB file
- Extract amino acid sequence
- Return NumPy array and sequence string

**validate_pdb_files(pdb1, pdb2)**
- Check file existence
- Verify file permissions
- Raise descriptive errors

**compare_pyrosetta(pdb1, pdb2)**
- Initialize PyRosetta
- Load both structures as Pose objects
- Calculate TM-score and GDT-TS metrics
- Return results dictionary

**compare_mammoth(pdb1, pdb2, mammoth_path)**
- Validate Mammoth executable exists
- Run Mammoth via subprocess
- Parse stdout output
- Extract TM-score, Z-score, alignment info

**compare_tmtools(pdb1, pdb2)**
- Parse PDB coordinates
- Call tm_align function
- Extract TM-score, RMSD, alignment
- Return structured results

### Output Format

```
======================================================================
PROTEIN STRUCTURE COMPARISON RESULTS
======================================================================

Structure 1: proteinA.pdb (150 residues)
Structure 2: proteinB.pdb (200 residues)

Method: PyRosetta

Lengths:
  Structure 1: 150 residues
  Structure 2: 200 residues
  Difference: 50 residues

TM-score: 0.785
Interpretation: Same structural fold (TM-score >= 0.5)

Mammoth Z-score: 5.234
Interpretation: Higher Z-scores indicate more significant similarity

Mammoth MaxSub: 0.812

RMSD (aligned regions): 2.456 Å

GDT-TS Scores (distance thresholds):
  GDT-TS (1Å):  0.823
  GDT-TS (2Å):  0.891
  GDT-TS (4Å):  0.954
  GDT-TS (8Å):  0.978

Aligned residues: 142

Sequence alignment:
  A L V I N A - - - -
  A L V I N A - - - -
  ...
======================================================================
```

### CSV Output Format

```csv
Metric,Value
Method,PyRosetta
Structure1,proteinA.pdb
Structure2,proteinB.pdb
tm_score,0.7854
gdt_ts1,0.8231
gdt_ts2,0.8912
gdt_ts4,0.9538
gdt_ts8,0.9784
length1,150
length2,200
z_score,5.2341
```

## Usage Examples

### Single Comparison
```bash
# Using PyRosetta (default)
python protein_structure_comparison.py protein1.pdb protein2.pdb

# Using Mammoth
python protein_structure_comparison.py protein1.pdb protein2.pdb --method mammoth --mammoth-path /usr/local/bin/mammoth

# Using tmtools
python protein_structure_comparison.py protein1.pdb protein2.pdb --method tmtools
```

### Batch Processing
```bash
# Create pairs.txt file:
proteinA.pdb proteinB.pdb
proteinC.pdb proteinD.pdb
proteinE.pdb proteinF.pdb

# Run batch comparison
python protein_structure_comparison.py --batch pairs.txt --output results.csv --method pyrosetta
```

### Command Line Options

```
positional arguments:
  pdb_files           Two PDB files to compare

optional arguments:
  -h, --help          Show help message and exit
  -m, --method         Comparison method (default: pyrosetta)
                       Options: pyrosetta, mammoth, tmtools
  --mammoth-path       Path to Mammoth executable (default: mammoth)
  -o, --output         Output file for saving results (CSV format)
  -b, --batch          Batch mode: file containing PDB pairs
  -v, --verbose         Verbose output
```

## Dependencies

### Required
- Python 3.7+
- NumPy
- At least one of:
  - PyRosetta (for pyrosetta method)
  - Mammoth executable (for mammoth method)
  - tmtools package (for tmtools method)

### Installation

```bash
# Install NumPy
pip install numpy

# Install PyRosetta
# Download from: https://www.pyrosetta.org/
# Follow installation instructions

# Install tmtools (alternative)
pip install tmtools

# Install Mammoth (requires SBGrid or manual compilation)
sbgrid-cli install mammoth-mult
# Or compile from source code
```

## Error Handling

The script includes comprehensive error handling for:

1. **File Errors**
   - File not found
   - Permission denied
   - Invalid PDB format

2. **Import Errors**
   - PyRosetta not installed
   - tmtools not available
   - Missing dependencies

3. **Execution Errors**
   - Mammoth executable not found
   - Mammoth execution timeout
   - Coordinate extraction failures

4. **Data Errors**
   - No C-alpha atoms in PDB
   - Empty structures
   - Invalid residue information

## Features

### Core Features
- Three comparison methods (PyRosetta, Mammoth, tmtools)
- Length-independent TM-score calculation
- Multiple GDT-TS distance thresholds
- Batch processing capability
- CSV output for analysis
- Statistical interpretation (Z-score)
- Comprehensive error handling
- Detailed formatted output

### Advanced Features
- Sequence alignment display
- Residue mapping information
- Progress indicators for batch processing
- Verbose mode for debugging
- Interpretation of TM-score significance

## Advantages Over RMSD

### Why TM-score > RMSD for Different Length Proteins

1. **Length Normalization**
   - TM-score: Normalized by protein length
   - RMSD: Dominated by largest errors, length-dependent

2. **Sensitivity to Global Fold**
   - TM-score: Weights small distance errors more heavily
   - RMSD: Equal weighting, dominated by large deviations

3. **Interpretability**
   - TM-score: Standard 0-1 scale with clear thresholds
   - RMSD: No universal scale, context-dependent

4. **Statistical Significance**
   - TM-score: Can derive Z-score for significance
   - RMSD: No statistical framework

## Testing Recommendations

### Test with Sample Structures
1. Use PDB structures with known similarities
2. Compare proteins of same fold (TM-score > 0.5)
3. Compare proteins of different folds (TM-score < 0.2)
4. Verify batch processing with multiple pairs

### Validation
- Compare results from different methods
- Validate against known structural classifications
- Check TM-score interpretation accuracy
- Verify batch output integrity

## Citation Information

When using this tool, please cite:

### For PyRosetta Method
- PyRosetta: Chaudhury, S., Lyskov, S., Gray, J.J. (2010)
  PyRosetta: a Python-based interactive platform for developing protein modeling software.

### For Mammoth Method
- Ortiz, A.R., Strauss, C.E.M., Olmea, O. (2002)
  MAMMOTH: Matching molecular models obtained from theory.
  Protein Science, 11(11): 2606-2621.
  DOI: 10.1110/ps.0215902

### For TM-align Method
- Zhang, Y., Skolnick, J. (2005)
  TM-align: A protein structure alignment algorithm based on the TM-score.
  Nucleic Acids Research, 33(11): 2302-2309.
  DOI: 10.1093/nar/gki342

## Future Enhancements

Potential improvements:
1. Add support for multiple chain comparison
2. Implement alternative alignment algorithms (CE, DALI)
3. Add visualization support (PyMOL integration)
4. Include secondary structure comparison
5. Add domain-level comparison
6. Implement structure database searching
7. Add machine learning-based similarity prediction
8. Support for non-standard amino acids