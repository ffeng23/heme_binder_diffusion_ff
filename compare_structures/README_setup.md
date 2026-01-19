# README to set up to run compare protein structures

## NOTE:
In this folder, the python code were initially started by AI under opencode Big Pickle.

Of course, the original code were extensively modified by Feng to make sure it works as intended.

+ change the calling signature for pyrosetta GDT scoring

+ change the calling signature for pyrosetta CA_rmsd scoring.

+ change the calling signature for tmtool tm_align

+ change the calling signature and add averaging of tmtools tm_score

There are also other small things to make sure the project work.


## creat a new conda env by

  	conda env create -f environment.yml

## activate the env by

   	conda activate pyrosetta_env

## call the script to compare protein structures

There are example pdb files in the sub-directory, "data".

+ compare with pyrosetta to get GDT score


  	  python3 protein_compare_structure.py data/HEM_1.D12.pdb data/HEM_3.G7.pdb -m pyrosetta


+ compare with tmtools to get tm score

  	  python3 protein_compare_structure.py data/HEM_1.D12.pdb data/HEM_3.G7.pdb -m tmtools


+ compare with Mammoth to get z-score etc.

  	  python3 protein_compare_structure.py data/HEM_1.D12.pdb data/HEM_3.G7.pdb -m mammoth

Note: I haven't spent time working making mammoth work. This one need to install Mammoth (see README.md).


## the functions in "protein_compare_structure.py" can also be imported to other python script.
