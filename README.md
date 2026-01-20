#RDKit Molecular Modification Pipeline
Overview
This repository contains a Python script that systematically generates chemical analogs from an input SDF file using SMARTS-based transformations and 3D structure generation with RDKit.
It is designed for medicinal chemistry optimization, bioisosteric replacement, and structure–property exploration.

The script performs two rounds of molecular transformations, generating optimized 3D conformations and saving all results as individual SDF files with a detailed transformation log.

#Key Features

==> SMARTS-driven functional group replacements

==> Two-round iterative molecular expansion

==> Bioisosteric, polarity, and lipophilicity-driven transformations

==> Automatic hydrogen addition, 3D embedding, and force-field optimization

==> Structured output directories for each round

==> CSV log of all transformations and generated molecules

#Transformation Categories

The script includes a curated library of medicinal chemistry transformations, including:

1.  Targeted group:
-CH3
    Modified groups:
-NH2, -NH(CH3), -N(CH3)2, -CF3, -CHO, -COOH, 
-CONH2, -CH2OH, -C6H4CN, -SCN, -SC6H5,
-C6H4NH2
2.  Targeted group:
-CN
    Modified groups:
-COOH, -CONH2, -NH2, -NH(CH3), -N(CH3)2, 
-NHC6H5, -C6H4CN, -CHO
3.  Targeted group:
    Gamma-lactam ring (2-pyrrolidinone)
4.  Modified groups
2-piperazinone
    Targeted group:
-CONH
    Modified groups
-NH2,-CONH2, -O�C(=O), �SO2?�NH2, -NH2, 
-NH(CH3), -N(CH3)2, -CHO, -COOH, -CONH2, -NHC6H5,-C6H4CN
5.  Targeted group:
-CF3
    Modified groups
-SO2NH2, -CN, -COOH, -CH3, -OCH3, -NH2, -CHO, -CONH2, 
-NH(CH3), -N(CH3)2, -NHC6H5, -C6H4CN
6.  Targeted group:
-F
    Modified groups
-NH2, -CHO, -COOH, CONH2, -NHCH3, 
-N(CH3)2, -C6H4NH2, -C6H4CN, -SCN, -SC6H5

#Requirements
1. Python Version

Python 3.8+

2. Required Libraries
pip install rdkit pandas


RDKit is recommended to be installed via conda:

conda install -c conda-forge rdkit

#Input

initial_molecules.sdf or Multi-molecule SDF file

#Output Files
1. SDF Files

==> One file per generated molecule

==> Fully hydrogenated

==> 3D embedded using ETKDG

==> Optimized with UFF

2. Transformation Log (transformation_log.csv)

Includes:

==> Transformation round

==> Parent molecule

==> Applied transformation

==> Product name

==> Output file path

#How It Works

==> Load molecules from the input SDF

==> Apply all SMARTS reactions (Round 1)

==> Generate 3D coordinates and optimize structures

==> Apply all reactions again to Round 1 products (Round 2)

==> Log all transformations to CSV

#Running the Script
python Structure_modification.py


Note: Ensure initial_molecules.sdf is in the same directory or update the path in the script.


If you use this workflow in research or publications, please cite RDKit:

RDKit: Open-source cheminformatics; http://www.rdkit.org

