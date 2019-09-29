import argparse
import os
from Bio import SeqIO

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


def run_nw(query, ref):
	print(query.seq)
	print(ref.seq)
	

if __name__ == '__main__':
	# Parse command-line arguments
	parser = make_arg_parser()
	args = parser.parse_args()

	# Parse query and ref data
	query = list(SeqIO.parse(args.query, 'fasta'))[0]
	ref = list(SeqIO.parse(args.ref, 'fasta'))[0]

	# If matches file is given, 
	if args.matches is None:
		run_nw(query, ref)
	else: 
		print('Matches supplied')

