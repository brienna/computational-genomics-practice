Implements the [Nei-Saitou neighbor-joining algorithm](https://en.wikipedia.org/wiki/Neighbor_joining) for phylogeny construction, with estimation of bootstrap support. 

To run the program: `python main.py -m [MSA.fna]` where MSA.fna is the FASTA file containing the multiple sequence alignment. (Bootstrapping takes a while.)

To visualize tree with bootstrap confidences: `Rscript hw3-plot-edges.r edges/edges.txt tip-labels.txt nj-solution-bootstrap.txt`. The resulting file `tree.pdf` overwrites the original file from Step 3, so be sure to manually rename the original file as soon as it's made, or change the R script to take an ID.

Files:

- `main.py` — Entry point. 
- `boot.py` Code for bootstrapping, given a neighbor joined tree 
- `hw3.fna` — MSA of 61 bacterial 16S subunit ribosomal RNA sequences
- `pairwise_dissimilarity.txt` — Pairwise dissimilarity matrix based on MSA
- `edges/edges.txt` — [ancestor, descendant, branch_length] representation of preordered tree
- `newick/tree.txt` — NEWICK representation of postordered tree
- `hw3-plot-edges.r` — External script to visualize tree with `edges.txt`
- `hw3-plot-newick.r` — External script to visualize tree with `tree.txt`
- `tree.pdf` — tree visualization
- `tip-labels.txt` — tab-delimited file containing [seqID, Phylum, color]
- `bootstrap_tree.pdf` — Visualization of tree with bootstrap confidences


