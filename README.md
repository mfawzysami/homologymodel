# Protein Homology Modeling Service

###Introduction
This python script will automate all homology modeling tasks for predicting tertiary structure for unknown protein sequence just by knowing the primary structure of the protein (Amino acids Sequence).

It starts with a fasta file containing the amino acids sequence for the query protein. Then, it performs the following steps:

1. Connect to NCBI Blast Server to perform Blastp Search using the query protein sequence.
2. It downloads the top matched homologs.
3. It builds an internal MODELLER database of top matched homologs.
4. It performs Multiple Sequence alignment for the retrieved homologs.
5. It performs structural alignment for the downloaded homologs.
6. It performs Template-based structure homology modelling by assigning secondary structure fragments to different regions of the query sequence.
7. It performs iterative refinement of the model using geometrical calculations and energy based calculations using CHARMM ("Chemistry at Harvard for Macromolecular Mechanics")
8. It generates dendrogram of the Top similar homologs.
9. You can specify the closer PDB file to use in order to generate the predicted 3D model for the given query protein.

## How to Install

1. Download and install MODELLER, from this link (https://salilab.org/modeller/download_installation.html)

2. Install BioPython package `pip install biopython`

3. Run this script `python pms.py --help`

```text
usage: pms.py [-h] [-i INPUT] [-o OUTPUT] [-r RESULTS] [-d [NODOWNLOADS]]
              [-s MIN] [-m MAX] [-c [CLEAN]] [-t ITERATIONS] [-k [CHECK]]
              [-e EVALUE] [-n COUNT] [-b [BUILD]] [-u TAKE] [-z MODELS]

Protein Homology Modelling Service.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Query protein Fasta file
  -o OUTPUT, --output OUTPUT
                        Output Directory where all downloaded files are stored
                        alongside predicted query protein model(s). Defaults
                        to current working directory
  -r RESULTS, --results RESULTS
                        Blast results XML file if exists for the query
                        sequence. if not provided, the software will perform
                        blast search against the query protein sequence
  -d [NODOWNLOADS], --nodownloads [NODOWNLOADS]
                        If specified PMS will not download PDB files but it
                        will search the output directory for these files
                        during downstream analysis
  -s MIN, --min MIN     Sequences below this length will be ignored during the
                        modelling process. Defaults : 30 residues
  -m MAX, --max MAX     Sequences above this length will be ignored during the
                        modeling process. Defaults: 4000 residues
  -c [CLEAN], --clean [CLEAN]
                        Clean sequences from all non-standard protein
                        residues. Defaults: True
  -t ITERATIONS, --iterations ITERATIONS
                        Number of Search iterations
  -k [CHECK], --check [CHECK]
                        Check the modeling profile for deviations . Defaults:
                        False
  -e EVALUE, --evalue EVALUE
                        Include sequences whose E-value is larger or equal to
                        this provided value. Defaults: 0.01
  -n COUNT, --count COUNT
                        Number of blast search results to use for model
                        building.
  -b [BUILD], --build [BUILD]
                        Instruct PMS to build a homology model for a specific
                        PDB
  -u TAKE, --take TAKE  The specific PDB ID to use to create a homology model
                        for the unknown query protein
  -z MODELS, --models MODELS
                        Number of predicted models to generate for the unknown
                        protein
```

All Parameters are self-explanatory.



## How to Use

- You should have your unknown protein sequence in a fasta file.

`python pms.py --input=yourfile.fa --output=./myfiles`

- If you have blast results in XML file, you can directly instruct the script to use this file instead of connecting to NCBI blast service.

`python pms.py --input=yourfile.fa --output=./myfiles  --results=results.xml`

- The script will perform NCBI blast Search and will also download the top matched Homologs' PDB files. If you have already these files in "myfiles" directory. you can instruct the script not to download these files again using `--nodownloads` switch


`python pms.py --input=yourfile.fa --output=./myfiles  --results=results.xml --nodownloads`

- you can specify an exact number of homologs to use by specifying `--count` switch

`python pms.py --input=yourfile.fa --output=./myfiles  --results=results.xml --count 30`


- The script will output a dendrogram of the best structurally aligned homologs. Choose the best one. and then instruct the script to align and generate the predicted model using the following command.

`python pms.py --input=6LU7.fa --results=results.xml --output=./myfiles --nodownloads --count 10 --build --take 5R7Y`

Here, we instruct the script to build the homology model using `5R7Y.pdb` file which was the best in my experiment.

###Note:

5R7Y.pdb file should be in the same directory as the running script for it to work properly.

