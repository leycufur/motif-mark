# Motif Marking Tool

## Description

The Motif Marking Tool is a Python-based utility for visualizing gene structures and motif occurrences in DNA sequences. It generates graphical representations highlighting exons, introns, and motif positions on a gene.

## Installation

To install the Motif Marking Tool, follow these steps:

1. Clone the repository to your local machine and navigate to the directory.
    ```bash
   git clone https://github.com/leycufur/motif-mark.git
   ```
   ```bash
   cd motif-mark
   ```
2. Install requirements
    ```bash
    conda create -n "motif-mark_env"
    conda activate motif-mark_env
    ```
    ```bash
    pip install pycairo==1.26.0
    ```
##

## Usage
### The following command will run the script using the input FASTA and motif files to render a .png image that graphically illustrates each gene along with its introns regions, exon regions, and motifs. 
    ```bash
    ./motif-mark.py -f input.fasta -m input_motif.txt
    ```