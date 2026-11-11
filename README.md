# DiceyTree

## Background
Create a tree for phage genomes based on homology of shared proteins across the full genome.

This is mainly to deal with the problem that phages don't have an universal marker gene and even when they do, multiple sequence alignments have more holes than swiss cheese. 

The solution presented here is to compute the similarity between two phage genomes based on the homology across all the shared genes. 

A very similar method was first used for phage genome trees in the Phage Proteomic Tree by Rohwer and Edwards.  
Rohwer F, Edwards R. The Phage Proteomic Tree: a genome-based taxonomy for phage. J Bacteriol. 2002 Aug;184(16):4529-35. doi: 10.1128/JB.184.16.4529-4535.2002

However, this is only based on the presence/absence of shared genes. 

The homology we use is represented by local alignment(blastp) and the similarity is calculated with a Dice coefficient (also known as the Sorensen-Dice Index) using a summation of the bitscores of all proteins.

$$
d(A,B) = 1 - \frac{2 \cdot \max\left(S_{AB},\, S_{BA}\right)}{S_{AA} + S_{BB}}
$$

where $S_{AB}$ is the sum of bitscores of all hits from proteins of genome $A$
against proteins of genome $B$ (after filtering), and $S_{AA}$, $S_{BB}$ are the
self-hit bitscore sums of each genome.

## Inputs
1. Call your genes (like with Pharokka)
2. [Diamond blastp]() of all your genes vs. all your genes  
    For this you need to make all your genes into a database and then run all-v-all
    ```bash
    diamond makedb --in all_genes.faa -d ./all_genes_diamond_db

    diamond blastp -d ./all_genes_diamond_db -q all_genes.faa -o all_v_all_diamond.out \
    -k 0 --evalue 0.0001 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp
    ```


## Outputs

## Command summary

```bash
usage: new_dice_tree.py [-h] -i INPUT [-out_csv OUT_CSV] [-out_tree OUT_TREE] [-e EVALUE] [-p PIDENT] [-tc TCOVER] [-qc QCOVER] [--sep SEP]

Input blastp output for all-v-all genes between any genomes --> output is a table of distances between genomes based on shared genes and a neighbour-joining (distance-based) tree

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to all-v-all blastp file
  -out_csv OUT_CSV, --out_csv OUT_CSV
                        tab-separated file of pairwise distances between genomes
  -out_tree OUT_TREE, --out_tree OUT_TREE
                        tree file (in newick format)
  -e EVALUE, --evalue EVALUE
                        evalue threshold to use (default: 1e-5)
  -p PIDENT, --pident PIDENT
                        minimum percent identity (default: 40)
  -tc TCOVER, --tcover TCOVER
                        minimum target cover in % (default: none; needs column 14 of the blast file)
  -qc QCOVER, --qcover QCOVER
                        minimum query cover in % (default: none; needs column 13 of the blast file)
  --sep SEP             text separating genome name from protein number in the IDs (default: _CDS_)
```