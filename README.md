# DiceyTree

## Background
Create a tree for phage genomes based on homology of shared proteins across the full genome.

**Problem**: This is mainly to deal with the problem that phages don't have an universal marker gene and even when they do, multiple sequence alignments have more holes than swiss cheese. 

**Solution**: The solution presented here is to compute the similarity between two phage genomes based on the homology across all the shared genes. 

A very similar method was first used for phage genome trees in the Phage Proteomic Tree by Rohwer and Edwards.  
Rohwer F, Edwards R. The Phage Proteomic Tree: a genome-based taxonomy for phage. J Bacteriol. 2002 Aug;184(16):4529-35. doi: 10.1128/JB.184.16.4529-4535.2002

However, this is only based on the presence/absence of shared genes. 

The homology we use is represented by local alignment(blastp) and the similarity is calculated with a Dice coefficient (also known as the Sorensen-Dice Index) using a summation of the bitscores of all proteins.

$$
d(A,B) = 1 - \frac{2 \cdot \max\left(S_{AB} \, S_{BA}\right)}{S_{AA} + S_{BB}}
$$

where $S_{AB}$ is the sum of bitscores of all hits from proteins of genome $A$
against proteins of genome $B$ (after filtering), and $S_{AA}$, $S_{BB}$ are the
self-hit bitscore sums of each genome.

**FYI** : This is NOT a phylogenetic tree. It is a relatedness tree and is therefore naturally unrooted. 

## Requirements

Make sure you the following packages installed in your python/conda environment:
- argparse, sys
- pandas
- numpy
- biopython

## Inputs
1. Call your genes (like with Pharokka)
2. [Diamond blastp]() of all your genes vs. all your genes  
    For this you need to make all your genes into a database and then run all-v-all
    ```bash
    diamond makedb --in all_genes.faa -d ./all_genes_diamond_db

    # I choose to prefilter my diamond blast with an evalue because otherwise it reports way too many hits. You can also filter by pident >=20 or mismatches etc

    diamond blastp -d ./all_genes_diamond_db -q all_genes.faa -o all_v_all_diamond.out \
    -k 0 --evalue 0.0001 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp
    ```

The input into DiceyTree is a blastp output with 12 - 14 columns (depending on if you want the query and subject coverage or not). The script accepts 

## Outputs

1. "dice_distances.tsv" : pairwise distances of genomes  
    genomeA | genomeB | 0.65  
    genomeA | genomeC | 0.3  
    genomeA | genomeD | 1  

A distance of 1 = no shared proteins

If you remove the ones, you can also represent this as a network in Cytoscape instead of a tree

2. "dice_tree.nwk" : a neighbour joining tree in newick format

## Command summary

```
usage: dice_tree.py [-h] -i INPUT [-out_csv OUT_CSV] [-out_tree OUT_TREE] [-e EVALUE] [-p PIDENT] [-tc TCOVER] [-qc QCOVER] [--sep SEP]

All-v-all blastp output -> Dice distances between genomes (from summed bitscores of shared genes) -> neighbour-joining tree.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to all-v-all blastp file
  -out_csv OUT_CSV, --out_csv OUT_CSV
                        tab-separated file of pairwise distances between genomes (default: dice_distances.tsv)
  -out_tree OUT_TREE, --out_tree OUT_TREE
                        tree file in newick format (default: dice_tree.nwk)
  -e EVALUE, --evalue EVALUE
                        evalue threshold to use (default: 1e-5)
  -p PIDENT, --pident PIDENT
                        minimum percent identity (default: 40)
  -tc TCOVER, --tcover TCOVER
                        minimum target cover in % (default: none; needs column 14 of the blast file)
  -qc QCOVER, --qcover QCOVER
                        minimum query cover in % (default: none; needs column 13 of the blast file)
  --sep SEP             text separating genome name from protein number in the IDs (default: _CDS_)

input format:
  Tab-separated BLAST output (-outfmt 6), one hit per line:

    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore [qcov tcov]

  Sequence IDs must look like <genome>_CDS_<number>, e.g. genome_A_icnght_CDS_0001
  (genome = genome_A_icnght, protein = CDS_0001). Use --sep if your IDs use a
  different separator. Columns 13 and 14 (query / target cover in %) are optional
  and only needed for -qc / -tc. For example, DIAMOND writes them with:

    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

  (NCBI blastp has qcovs for query cover but no target-cover field.)

method:
  1. keep hits passing the evalue, identity and (optional) cover thresholds
  2. sum bitscores for every genome pair -> S[a,b]
  3. Dice distance: d(a,b) = 1 - 2*max(S[a,b], S[b,a]) / (S[a,a] + S[b,b])
     (pairs with no hits get distance 1)
  4. neighbour-joining tree from the distance matrix

output:
  out_csv   tab-separated table with columns seq1, seq2, dice
  out_tree  neighbour-joining tree in newick format

examples:
  python3 dice_tree.py -i allvall.blastp
  python3 dice_tree.py -i allvall.blastp -out_csv dist.tsv -out_tree tree.nwk -e 1e-10 -p 50
  python3 dice_tree.py -i allvall.blastp -qc 70 -tc 70
```