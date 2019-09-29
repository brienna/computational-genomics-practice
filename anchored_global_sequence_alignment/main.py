import argparse
import os
from Bio import SeqIO
import numpy as np

def make_arg_parser():
	'''Parses arguments.'''

	parser = argparse.ArgumentParser(prog='main.py',
									 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-q', '--query',
						required=True,
						default=argparse.SUPPRESS, # no default value set, not even None
						help='Path to query sequence fasta [required]')
	parser.add_argument('-r', '--ref',
						required=True,
						default=argparse.SUPPRESS,
						help='Path to reference sequence fasta [required]')
	parser.add_argument('-m', '--matches',
						required=False,
						help='Path to matches text file')
	return parser


def run_nw(query, ref, match=1, gap=-1, mismatch=-1):
	r_seq = 'GCATGCU'
	q_seq = 'GATTACA'

	# Set up grid, ref on top, query on left side
	grid = np.empty((len(q_seq) + 1, len(r_seq) + 1)) # +1 for a gap at beginning of sequences
	grid[:] = np.nan

	# Populate first row and column based on gap penalty
	grid[0] = np.arange(0, -len(grid[0]), -1) # 1st row
	grid[:,0] = np.arange(0, -len(grid[:,0]), -1) # 1st column

	grid_iterator = np.nditer(grid, flags=['multi_index']) # flag allows to track individual cells, providing tuple
	for cell in grid_iterator:
		r = grid_iterator.multi_index[1]
		q = grid_iterator.multi_index[0]

		# Skip the first cell in the top left
		if q < 1 or r < 1: 
			continue
		
		# previous diagonal, top, left scores
		scores = np.array([grid[q - 1, r - 1], grid[q - 1, r], grid[q, r - 1]])
		
		# If pairing is a match, then
		if r_seq[r - 1] == q_seq[q - 1]:
			scores[0] += match # moving diagonally, only add +1
			scores[1] += gap # moving down represents an indel, so add score for indel
			scores[2] += gap # moving right also represents an indel
		else:
			scores[0] += mismatch # moving diagonally, only add +1
			scores[1] += gap # moving down represents an indel, so add score for indel
			scores[2] += gap # moving right also represents an indel
		
		grid[q, r] = scores.max()
	print(grid)
	

if __name__ == '__main__':
	# Parse command-line arguments
	parser = make_arg_parser()
	args = parser.parse_args()

	# Assign scores
	match = 1
	gap = -2
	mismatch = -3

	# Parse query and ref data
	query = list(SeqIO.parse(args.query, 'fasta'))[0]
	ref = list(SeqIO.parse(args.ref, 'fasta'))[0]

	# If matches file is given, 
	if args.matches is None:
		run_nw(query, ref, match, gap, mismatch)
	else: 
		print('Matches supplied')

