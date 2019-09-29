import argparse
import os

def make_arg_parser():
	'''Parses arguments.'''

	parser = argparse.ArgumentParser(prog='main.py',
									 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-q', '--query',
						required=True,
						default=argparse.SUPPRESS, # stop --help from displaying default: None
						help='Path to query sequence fasta [required]')
	parser.add_argument('-r', '--ref',
						required=True,
						default=argparse.SUPPRESS,
						help='Path to reference sequence fasta [required]')
	parser.add_argument('-m', '--matches',
						required=False,
						default=argparse.SUPPRESS,
						help='Path to matches text file')
	return parser


if __name__ == '__main__':
	# Parse command-line arguments
	parser = make_arg_parser()
	args = parser.parse_args()

