# Protein Structure Comparison Tool

A comprehensive Python script for comparing protein structures of different lengths using multiple methods: PyRosetta, Mammoth, and tmtools.

## Features

- **Three Comparison Methods**: PyRosetta, Mammoth, and tmtools
- **Length-independent TM-score calculation**: Suitable for proteins of different lengths
- **Multiple GDT-TS distance thresholds**: 1Å, 2Å, 4Å, 8Å
- **Batch processing capability**: Compare multiple protein pairs
- **CSV output format**: For data analysis and record keeping
- **Comprehensive error handling**: Clear error messages and graceful degradation
- **Statistical interpretation**: TM-score significance and Z-score support

## Installation

### Required Dependencies
```bash
pip install numpy
```

### Optional Dependencies (choose at least one)

#### PyRosetta (Recommended)
Download from: https://www.pyrosetta.org/
Follow the installation instructions for your platform.

#### tmtools (Python-native)
```bash
pip install tmtools
```

#### Mammoth (C++ executable)
Install via SBGrid:
```bash
sbgrid-cli install mammoth-mult
```
Or compile from source code.

## Usage

### Single Comparison
```bash
# Using PyRosetta (default)
python3 protein_structure_comparison.py protein1.pdb protein2.pdb

# Using tmtools
python3 protein_structure_comparison.py protein1.pdb protein2.pdb --method tmtools

# Using Mammoth
python3 protein_structure_comparison.py protein1.pdb protein2.pdb --method mammoth --mammoth-path /usr/local/bin/mammoth
```

### Batch Processing
```bash
# Create pairs.txt file:
proteinA.pdb proteinB.pdb
proteinC.pdb proteinD.pdb
proteinE.pdb proteinF.pdb

# Run batch comparison
python3 protein_structure_comparison.py --batch pairs.txt --output results.csv --method pyrosetta
```

### Command Line Options
```
positional arguments:
  pdb_files             Two PDB files to compare

options:
  -h, --help            Show help message and exit
  -m, --method {pyrosetta,mammoth,tmtools}
                        Comparison method (default: pyrosetta)
  --mammoth-path MAMMOTH_PATH
                        Path to Mammoth executable (default: mammoth)
  -o, --output OUTPUT   Output file for saving results (CSV format)
  -b, --batch BATCH     Batch mode: file containing PDB pairs
  -v, --verbose         Verbose output
```

## Output Format

### Console Output
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

RMSD (aligned regions): 2.456 Å

GDT-TS Scores (distance thresholds):
  GDT-TS (1Å):  0.823
  GDT-TS (2Å):  0.891
  GDT-TS (4Å):  0.954
  GDT-TS (8Å):  0.978

======================================================================
```

### CSV Output
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
```

## TM-score Interpretation

- **< 0.2**: Random/unrelated structures
- **0.2-0.5**: Different folds, some similarity
- **> 0.5**: Same structural fold
- **1.0**: Perfect match

## Testing

### Run Tests
```bash
python3 test_protein_comparison.py
```

### Run Demo
```bash
python3 demo_comparison.py
```

The demo creates sample PDB files and shows basic functionality without requiring external dependencies.

## Methods Comparison

### PyRosetta (Recommended)
- **Advantages**: Pure Python, fast, comprehensive metrics
- **Metrics**: TM-score (CA_gdttm), GDT-TS calculated using CA_gdtmm with threshold scaling, RMSD
- **Installation**: Requires PyRosetta license and installation
- **Functions**: Uses CA_gdttm (from core.scoring) for TM-score, CA_gdtmm for GDT-based scores, and CA_superimpose (from protocols.toolbox) for structural alignment

### Mammoth
- **Advantages**: Established algorithm, statistical significance (Z-score)
- **Metrics**: TM-score, Z-score, MaxSub, alignment info
- **Installation**: Requires C++ executable

### tmtools
- **Advantages**: Python-native, easy installation
- **Metrics**: TM-score, RMSD, alignment length
- **Installation**: Simple pip install

## Why TM-score over RMSD?

1. **Length Normalization**: TM-score is normalized by protein length, making it suitable for comparing proteins of different sizes
2. **Global Fold Sensitivity**: TM-score weights small distance errors more heavily, making it more sensitive to overall fold similarity
3. **Interpretability**: Standard 0-1 scale with clear thresholds
4. **Statistical Framework**: Can derive statistical significance (Z-score)

## Error Handling

The script includes comprehensive error handling for:
- File not found or permission errors
- Invalid PDB format
- Missing dependencies
- Execution timeouts
- Coordinate extraction failures

## File Structure

```
protein_structure_comparison/
├── protein_structure_comparison.py  # Main script
├── test_protein_comparison.py       # Test suite
├── demo_comparison.py               # Demo with sample data
├── requirements.txt                 # Dependencies
└── README.md                       # This file
```

## Citation

When using this tool, please cite the appropriate method:

### PyRosetta
Chaudhury, S., Lyskov, S., Gray, J.J. (2010) PyRosetta: a Python-based interactive platform for developing protein modeling software.

### Mammoth
Ortiz, A.R., Strauss, C.E.M., Olmea, O. (2002) MAMMOTH: Matching molecular models obtained from theory. Protein Science, 11(11): 2606-2621.

### tmtools (TM-align)
Zhang, Y., Skolnick, J. (2005) TM-align: A protein structure alignment algorithm based on the TM-score. Nucleic Acids Research, 33(11): 2302-2309.

## License

This implementation follows the academic use license terms of the underlying comparison methods. Please ensure you have appropriate licenses for PyRosetta and other commercial tools.