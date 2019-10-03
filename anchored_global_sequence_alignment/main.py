import argparse
import os
from Bio import SeqIO
import numpy as np
import re

class Needleman_Wunsch_Executable(object):
	'''Template for Needleman Wunsch executable. '''

	def __init__(self, max_alignments, match_award=1, gap_penalty=-1, mismatch_penalty=-1):
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
		self.max_alignments = max_alignments
		self.ref_matches = None
		self.query_matches = None

		self.args = self.make_arg_parser().parse_args()
		self.ref_obj, self.query_obj = self.parse_sequences(self.args.ref, self.args.query)
		self.ref = self.ref_obj.seq 
		self.query = self.query_obj.seq

	def parse_sequences(self, ref, query):
		# Parse query and ref data / can assign this as a function later to check for fasta instead of passed through executable init
		print('Parsing ' + query + '...')
		query_obj = list(SeqIO.parse(query, 'fasta'))[0]
		print('Parsing ' + ref + '...')
		ref_obj = list(SeqIO.parse(ref, 'fasta'))[0]
		#self.ref = 'GCATGCUE'
		#self.query = 'GATTACA'
		return ref_obj, query_obj

	def parse_matches(self, m_file):
		print('Parsing ' + m_file + '...')
		
		fruit_fly_matches = []
		human_matches = []

		with open(m_file, 'r') as infile:
			for line in infile:
				positions = re.split('\t|\n', line)
		
				h_region = {'start': int(positions[0]),
							'end': int(positions[1])}
				ff_region = {'start': int(positions[2]),
							'end': int(positions[3])}
		
				human_matches.append(h_region)
				fruit_fly_matches.append(ff_region)

		# Identify which sequence is human
		if 'human' in self.query_obj.description.lower():
			self.query_matches = human_matches
			self.ref_matches = fruit_fly_matches
		else:
			self.query_matches = fruit_fly_matches
			self.ref_matches = human_matches

		print(self.ref_matches)
		print(self.query_matches)


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
		print('Building score and traceback grids...')

		# Set up score grid with reference on top, query on left
		score_grid = np.empty((len(self.query)+1, len(self.ref)+1)) # size grid to account for 1 gap at beginning of both sequences
		score_grid[:] = np.nan 
		score_grid[0] = np.arange(0, -len(score_grid[0]) * -self.gap_penalty, self.gap_penalty) # Populate 1st row based on gap penalty
		score_grid[:,0] = np.arange(0, -len(score_grid[:,0]) * -self.gap_penalty, self.gap_penalty) # Populate 1st column based on gap penalty

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
				arrow = 'DUH' 
			traceback_grid[q, r] = arrow
				
		print(score_grid)
		print(traceback_grid)
		self.score_grid = score_grid
		self.traceback_grid = traceback_grid
		print(score_grid[-1,-1])

	def traceback(self, start_a, a, start_b, b, candidate, trace_id=1):
		'''
		Probably could write this section better, with regard to checking arrows.
		Been getting wrong alignments if I do only one loop over all arrows.
		'''

		print('Tracing...')

		if a == start_a and b == start_b:
			self.trace_diagonally(a+1, b+1, candidate)
			return candidate
		
		# need a condition if one sequence has only one letter and the other sequence has more than one letter
		
		while a > start_a and b > start_b:
			# Get current arrow(s)
			arrows = list(self.traceback_grid[a, b])
			num_arrows = len(arrows)
			
			if num_arrows == 3:
				# Handle third arrow
				new_candidate3 = {'ref': candidate['ref'], 'query': candidate['query'], 'score': candidate['score']}
				arrow = arrows[2]
				a3 = a
				b3 = b
				if arrow == 'D':
					a3, b3 = self.trace_diagonally(a3, b3, new_candidate3)
				elif arrow == 'U':
					a3 = self.trace_up(a3, new_candidate3)
				elif arrow == 'H':
					b3 = self.trace_horizontally(b3, new_candidate3)
				if trace_id < self.max_alignments:
					trace_id += 1
					self.traceback(start_a, a3, start_b, b3, new_candidate3, trace_id)
				
				# Handle second arrow
				new_candidate2 = {'ref': candidate['ref'], 'query': candidate['query'], 'score': candidate['score']}
				arrow = arrows[1]
				a2 = a
				b2 = b
				if arrow == 'D':
					a2, b2 = self.trace_diagonally(a2, b2, new_candidate2)
				elif arrow == 'U':
					a2 = self.trace_up(a2, new_candidate2)
				elif arrow == 'H':
					b2 = self.trace_horizontally(b2, new_candidate2)
				if trace_id < self.max_alignments:
					trace_id += 1
					self.traceback(start_a, a2, start_b, b2, new_candidate2, trace_id)
				
			elif num_arrows == 2:
				# Handle second arrow
				new_candidate2 = {'ref': candidate['ref'], 'query': candidate['query'], 'score': candidate['score']}
				arrow = arrows[1]
				a2 = a
				b2 = b
				if arrow == 'D':
					a2, b2 = self.trace_diagonally(a2, b2, new_candidate2)
				elif arrow == 'U':
					a2 = self.trace_up(a2, new_candidate2)
				elif arrow == 'H':
					b2 = self.trace_horizontally(b2, new_candidate2)
				if trace_id < self.max_alignments:
					trace_id += 1
					self.traceback(start_a, a2, start_b, b2, new_candidate2, trace_id)
				
			# Handle first arrow
			arrow = arrows[0]
			if arrow == 'U':
				a = self.trace_up(a, candidate)
			elif arrow == 'D':
				a, b = self.trace_diagonally(a, b, candidate)
			elif arrow == 'H':
				b = self.trace_horizontally(b, candidate)  
		
		while a > start_a:
			a = self.trace_up(a, candidate)     
		while b > start_b:
			b = self.trace_horizontally(b, candidate)
		
		self.candidates.append(candidate)
		return candidate # return candidates[0] # return first alignment
	

	def traceback_with_matches(self):
		if len(self.query_matches) != len(self.ref_matches):
			print('Please fix your matches file and rerun.')
			return

		global_ref = ''
		global_query = ''

		# Identify non matching region at start, if any
		if self.ref_matches[0]['start'] != 0:
			region_end_ref = self.ref_matches[0]['start']
			print('Aligning 0—' + str(region_end_ref))
		
			if self.query_matches[0]['start'] != 0:
				region_end_query = self.query_matches[0]['start']
				print('Aligning 0—' + str(region_end_query))
				# Submit NW request on this segment
				alignment = self.traceback(0, region_end_query, 0, region_end_ref, {'ref': '', 'query': '', 'score': 0})
				global_ref += alignment['ref'] 
				global_query += alignment['query'] 
				print(alignment)
		
		for i, match in enumerate(self.query_matches): # just pick one to start iterating over, since they are shaped the same
			ref_segment = self.ref[self.ref_matches[i]['start'] : self.ref_matches[i]['end']]
			query_segment = self.query[match['start'] : match['end']]
			global_ref += ref_segment 
			global_query += query_segment 
			print('Adding ' + str(self.ref_matches[i]['start']) + '—' + str(self.ref_matches[i]['end']))
			print('Adding ' + str(self.query_matches[i]['start']) + '—' + str(self.query_matches[i]['end']))
			print(ref_segment)
			print(query_segment)
			
			# Identify region between current match and next match or sequence end
			region_start_ref = self.ref_matches[i]['end']
			if i+1 < len(self.ref_matches): # If we have a next match 
				region_end_ref = self.ref_matches[i+1]['start']
			else:
				region_end_ref = len(self.ref) 
			region_start_query = match['end']
			if i+1 < len(self.query_matches):
				region_end_query = self.query_matches[i + 1]['start']
			else:
				region_end_query = len(self.query) 
				
			# Submit NW requests on these segments
			print('Aligning ' + str(region_start_ref) + '—' + str(region_end_ref))
			print('Aligning ' + str(region_start_query) + '—' + str(region_end_query))
			alignment2 = self.traceback(region_start_query, region_end_query, region_start_ref, region_end_ref, {'ref': '', 'query': '', 'score': 0})
			global_ref += alignment2['ref'] 
			global_query += alignment2['query'] 

		return {'ref': global_ref, 'query': global_query}
		

	def run(self):
		# Build score and traceback grids. (Don't know if it is faster to build a complete grid if there are matches or to build multiple smaller grids.)
		self.build_grids()

		# If matches file is not given, 
		if self.args.matches is None:
			a = len(self.query)
			b = len(self.ref) 
			self.traceback(0, a, 0, b, {'ref': '', 'query': '', 'score': 0})
			print('Number of alignments found: ' + str(len(self.candidates)))
			print(self.candidates)
		else: 
			self.parse_matches(self.args.matches)
			alignment = self.traceback_with_matches()
			print(alignment['ref'])
			print(alignment['query'])
			# Ensure all letters have been captured
			print(str(alignment['ref']).replace('-', '') == self.ref)
			print(str(alignment['query']).replace('-', '') == self.query)
			print(len(str(alignment['ref']).replace('-','')))

if __name__ == '__main__':
	match_award = 1
	gap_penalty = -2
	mismatch_penalty = -3
	max_alignments = 1

	# Can set max_alignments = None, but this will lead to memory blowout w/ larger sequences
	# Need to fix this if want to return more alignments. 
	# Plus we are only interested in 1 alignment for now. 

	# Assumes that one of the sequences is human and that the description contains 'human' somewhere.

	# Create and run Needleman Wunsch executable 
	nw = Needleman_Wunsch_Executable(max_alignments, match_award, gap_penalty, mismatch_penalty)
	nw.run()




