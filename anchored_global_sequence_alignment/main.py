import argparse
import os
from Bio import SeqIO
import numpy as np

class Needleman_Wunsch_Executable(object):
	'''Template for Needleman Wunsch executable. '''

	def __init__(self, match_award=1, gap_penalty=-1, mismatch_penalty=-1):
		'''
		Initializes Needleman Wunsch executable.
		Required arguments: ref, query
		'''
		self.score_grid = None
		self.traceback_grid = None
		self.candidates = []
		self.match_award = match_award
		self.gap_penalty = gap_penalty
		self.mismatch_penalty = mismatch_penalty

		self.args = self.make_arg_parser().parse_args()
		self.ref = self.args.ref
		self.query = self.args.query
		self.parse_sequences(self.ref, self.query)


	def parse_sequences(self, ref, query):
		# Parse query and ref data / can assign this as a function later to check for fasta instead of passed through executable init
		#query = list(SeqIO.parse(args.query, 'fasta'))[0]
		#ref = list(SeqIO.parse(args.ref, 'fasta'))[0]
		self.ref = 'GCATGCU'
		self.query = 'GATTACA'


	def make_arg_parser(self):
		'''
		Parses command line arguments.
		'''
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


	def trace_diagonally(self, a, b, candidate):
		a -= 1
		b -= 1
		candidate['ref'] = self.ref[b] + candidate['ref']
		candidate['query'] = self.query[a] + candidate['query']
		if self.ref[b] == self.query[a]:
			candidate['score'] += match_award
		else:
			candidate['score'] += mismatch_penalty
		return a, b


	def trace_up(self, a, candidate):
		a -= 1
		candidate['ref'] = '-' + candidate['ref']
		candidate['query'] = self.query[a] + candidate['query']
		candidate['score'] += gap_penalty
		return a
		

	def trace_horizontally(self, b, candidate):
		b -= 1
		candidate['ref'] = self.ref[b] + candidate['ref']
		candidate['query'] = '-' + candidate['query']
		candidate['score'] += gap_penalty
		return b


	def build_grids(self):
		# Set up score grid with reference on top, query on left
		score_grid = np.empty((len(self.query)+1, len(self.ref)+1)) # size grid to account for 1 gap at beginning of both sequences
		score_grid[:] = np.nan 
		score_grid[0] = np.arange(0, -len(score_grid[0]), self.gap_penalty) # Populate 1st row based on gap penalty
		score_grid[:,0] = np.arange(0, -len(score_grid[:,0]), self.gap_penalty) # Populate 1st column based on gap penalty

		# Set up traceback grid with same shape as score grid
		traceback_grid = np.empty((len(self.query)+1, len(self.ref)+1), dtype=object)
		traceback_grid[:1,1:] = 'H' # 1st row is just horizontal movements 
		traceback_grid[1:,:1] = 'U' # 1st column is just up movements 

		# Iterate through score grid
		score_grid_iterator = np.nditer(score_grid, flags=['multi_index']) # flag allows to track individual cells, providing tuple
		for cell in score_grid_iterator:
			r = score_grid_iterator.multi_index[1] 
			q = score_grid_iterator.multi_index[0]

			# Skip the first cell at (0, 0)
			if q < 1 or r < 1: 
				continue

			# Get scores of 3 cells neighboring current cell (q, r)
			score_d = score_grid[q-1, r-1] # diagonal
			score_u = score_grid[q-1, r] # up
			score_h = score_grid[q, r-1] # horizontal (left)

			# Update scores based on possible movements and sequence similarity
			score_u += self.gap_penalty # moving down represents an indel, so add score for indel
			score_h += self.gap_penalty # moving right also represents an indel
			if self.ref[r-1] == self.query[q-1]:
				score_d += self.match_award # diagonal score adds match score
			else:
				score_d += self.mismatch_penalty # diagonal score adds mismatch penalty

			# Enter highest score into cell
			score_grid[q, r] = max(score_d, score_u, score_h)

			# Assign arrows to traceback matrix
			arrow = ''
			if score_d > score_u and score_d > score_h:
				arrow = 'D'
			elif score_u > score_d and score_u > score_h: 
				arrow = 'U'
			elif score_h > score_d and score_h > score_u:
				arrow = 'H'
			elif score_h == score_d and score_h > score_u: 
				arrow = 'DH'
			elif score_u == score_h and score_u > score_d:
				arrow = 'UH'
			elif score_d == score_u and score_d > score_h:
				arrow = 'DU'
			else: 
				print('Arrow not accounted for.') # will happen if 3 scores are all equal...
			traceback_grid[q, r] = arrow
				
		print(score_grid)
		print(traceback_grid)
		self.score_grid = score_grid
		self.traceback_grid = traceback_grid


	def traceback(self, a, b, candidate={'ref': '', 'query': '', 'score': 0}):
		while a > 0 and b > 0:
			# Get current arrow(s)
			arrows = list(self.traceback_grid[a][b])
			two_choices = len(arrows) > 1 # True if two arrows
			
			# If there is a second arrow
			if two_choices:
				# Split path into two 
				candidate2 = {'ref': candidate['ref'], 'query': candidate['query'], 'score': candidate['score']}
				arrow2 = arrows[1]
				a2 = a
				b2 = b
				if arrow2 == 'D':
					a2, b2 = self.trace_diagonally(a2, b2, candidate2)
				elif arrow2 == 'U':
					a2 = self.trace_up(a2, candidate2)
				elif arrow2 == 'H':
					b2 = self.trace_horizontally(b2, candidate2)
				self.traceback(a2, b2, candidate2)
			
			# Deal with the first arrow 
			arrow = arrows[0]
			if arrow == 'U':
				a = self.trace_up(a, candidate)
			elif arrow == 'D':
				a, b = self.trace_diagonally(a, b, candidate)
			elif arrow == 'H':
				b = self.trace_horizontally(b, candidate) 
		self.candidates.append(candidate)


	def run(self):
		# If matches file is not given, 
		if self.args.matches is not None:
			print('Matches supplied')
		else: 
			self.build_grids()
			b = len(self.ref) 
			a = len(self.query)
			self.traceback(a, b)
			print(self.candidates)


if __name__ == '__main__':
	match_award = 1
	gap_penalty = -1
	mismatch_penalty = -1

	# Create and run Needleman Wunsch executable 
	nw = Needleman_Wunsch_Executable(match_award, gap_penalty, mismatch_penalty)
	nw.run()



