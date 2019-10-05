
To run the program: `python main.py -q query.fa -r ref.fa [-m matches.txt]`

Implements an anchored version of the Needleman-Wunsch algorithm to align two proteins from human and another species. The anchored global sequence alignment assumes known matched regions between two sequences and applies Needleman-Wunsch to align the unaligned regions between the matched regions. 

Also permutes the amino acids in the sequences and repeats the alignment 10,000 times. Reports the score distribution with a histogram

Input: 
- query sequence formatted as FASTA 
- reference sequence formatted as FASTA 
- matches text file [optional]

Output:
- `alignment_[sequence1]_[sequence2].txt` — An optimal alignment, along with score.
- `histogram_[sequence1]_[sequence2].pdf` — Histogram showing score distribution for 10,000 alignments of randomly permuted reference and query sequences, along with a dashed line marking the score of the original alignment. 

Notes: 
- If you have a large sequence, you may not want to permute it 10,000 times. Comment out `nw.permute()`.
- Any special characters are treated the same as the ones in the alphabet, i.e. use the same match and mismatch costs.
- If the matches file is not provided, the program will run the standard Needleman-Wunsch on the entire sequences. 

