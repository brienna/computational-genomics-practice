# coding: utf-8

import argparse
import os
from Bio import SeqIO
import numpy as np
import re
import matplotlib.pyplot as plt


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
		self.alignment = None
		self.args = self.make_arg_parser().parse_args()
		self.ref_obj, self.query_obj = self.parse_sequences(self.args.ref, self.args.query) # contains description and other info in the FASTA file
		self.ref = self.ref_obj.seq 
		self.query = self.query_obj.seq


	def parse_sequences(self, ref, query):
		''''
		Parses given FASTA files for ref and query sequences.
		'''

		query_obj = list(SeqIO.parse(query, 'fasta'))[0]
		ref_obj = list(SeqIO.parse(ref, 'fasta'))[0]
		return ref_obj, query_obj


	def parse_matches(self, m_file):
		'''
		Parses a match file, in which the first 2 columns
		represent the start and end positions of matched regions
		for the human sequence, and the last 2 columns represent
		the start and end positions of the matched regions for 
		the other species.
		'''
		
		fruit_fly_matches = []
		human_matches = []

		# Iterate over each line in match file,
		# assigning positions to human & other species as specified 
		with open(m_file, 'r') as infile:
			for line in infile:
				positions = re.split('\t|\n', line) # split columns based on tab or newline
		
				h_region = {'start': int(positions[0]), # start of human matching region
							'end': int(positions[1])} # end of human matching region
				ff_region = {'start': int(positions[2]), # start of other species matching region
							'end': int(positions[3])} # end of other species matching region
		
				human_matches.append(h_region)
				fruit_fly_matches.append(ff_region)

		# Identify which sequence is human & assign matches to appropriate variables
		if 'human' in self.query_obj.description.lower():
			self.query_matches = human_matches
			self.ref_matches = fruit_fly_matches
		else:
			self.query_matches = fruit_fly_matches
			self.ref_matches = human_matches

		# print(self.ref_matches)
		# print(self.query_matches)


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
		'''
		Moves pointer to the diagonal neighboring cell during tracing.
		Updates candidate alignment based on diagonal movement & sequence match/mismatch.
		Returns updated part of pointer, a and b (grid column and row indexes).
		'''

		a -= 1
		b -= 1
		candidate['ref'] = self.ref[b] + candidate['ref'] 
		candidate['query'] = self.query[a] + candidate['query']
		if self.ref[b] == self.query[a]:
			candidate['score'] += match_award # if sequences match at this letter, add match award to alignment score
		else:
			candidate['score'] += mismatch_penalty # otherwise add mismatch penalty to alignment score
		return a, b


	def trace_up(self, a, candidate):
		'''
		Moves pointer to the top neighboring cell during tracing.
		Updates candidate alignment based on upward movement, indel. 
		Returns updated part of pointer, a (grid column index).
		'''

		a -= 1
		candidate['ref'] = '-' + candidate['ref'] # insert gap in reference sequence
		candidate['query'] = self.query[a] + candidate['query']
		candidate['score'] += gap_penalty # add gap penalty to alignment score
		return a
		

	def trace_horizontally(self, b, candidate):
		'''
		Moves pointer to the left neighboring cell during tracing.
		Updates candidate alignment based on leftward movement, indel. 
		Returns updated part of pointer, b (grid row index).
		'''

		b -= 1
		candidate['ref'] = self.ref[b] + candidate['ref']
		candidate['query'] = '-' + candidate['query'] # insert gap in query sequence
		candidate['score'] += gap_penalty # add gap penalty to alignment score
		return b


	def build_grids(self, ref, query):
		'''
		Builds score and traceback grids, with reference sequence
		on the top, query sequence on the side.
		'''

		# print('Building score and traceback grids...')

		# Set up score grid with reference on top, query on left
		score_grid = np.empty((len(query)+1, len(ref)+1)) # size grid to account for 1 gap at beginning of both sequences
		score_grid[:] = np.nan 
		score_grid[0] = np.arange(0, -len(score_grid[0]) * -self.gap_penalty, self.gap_penalty) # Populate 1st row based on gap penalty
		score_grid[:,0] = np.arange(0, -len(score_grid[:,0]) * -self.gap_penalty, self.gap_penalty) # Populate 1st column based on gap penalty

		# Set up traceback grid with same shape as score grid
		traceback_grid = np.empty((len(query)+1, len(ref)+1), dtype=object)
		traceback_grid[:1,1:] = 'H' # 1st row is just horizontal movements 
		traceback_grid[1:,:1] = 'U' # 1st column is just up movements 

		# Iterate through score grid
		score_grid_iterator = np.nditer(score_grid, flags=['multi_index']) # flag allows to track individual cells, providing tuple
		for cell in score_grid_iterator:
			r = score_grid_iterator.multi_index[1] # Column index
			q = score_grid_iterator.multi_index[0] # Row index

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
			if self.ref[r-1] == query[q-1]:
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
				
		# print(score_grid)
		# print(traceback_grid)
		self.score_grid = score_grid
		self.traceback_grid = traceback_grid


	def traceback(self, start_a, a, start_b, b, candidate, trace_id=1):
		'''
		Traces back through traceback grid to find optimal alignment.
		'''

		'''
		Probably could write this section better, with regard to checking arrows.
		Been getting wrong alignments if I do only one loop over all arrows.
		Something with scope probably.
		'''

		# print('Finding optimal alignment...')

		if a == start_a and b == start_b:
			self.trace_diagonally(a+1, b+1, candidate)
			return candidate
		
		while a > start_a and b > start_b:
			# Get current arrow(s)
			arrows = list(self.traceback_grid[a, b])
			num_arrows = len(arrows)
			
			# If all 3 neighbors have the same value, 
			if num_arrows == 3:
				# Handle third arrow
				new_candidate3 = {'ref': candidate['ref'], 'query': candidate['query'], 'score': candidate['score']}
				arrow = arrows[2]
				a3 = a
				b3 = b
				# Move traceback pointer based on direction given by arrow
				if arrow == 'D':
					a3, b3 = self.trace_diagonally(a3, b3, new_candidate3)
				elif arrow == 'U':
					a3 = self.trace_up(a3, new_candidate3)
				elif arrow == 'H':
					b3 = self.trace_horizontally(b3, new_candidate3)
				# Split into new 
				if trace_id < self.max_alignments:
					trace_id += 1
					self.traceback(start_a, a3, start_b, b3, new_candidate3, trace_id)
				
				# Handle second arrow
				new_candidate2 = {'ref': candidate['ref'], 'query': candidate['query'], 'score': candidate['score']}
				arrow = arrows[1]
				a2 = a
				b2 = b
				# Move traceback pointer based on direction given by arrow
				if arrow == 'D':
					a2, b2 = self.trace_diagonally(a2, b2, new_candidate2)
				elif arrow == 'U':
					a2 = self.trace_up(a2, new_candidate2)
				elif arrow == 'H':
					b2 = self.trace_horizontally(b2, new_candidate2)
				if trace_id < self.max_alignments:
					trace_id += 1
					self.traceback(start_a, a2, start_b, b2, new_candidate2, trace_id)
			
			# If two neighbors have the same value,
			elif num_arrows == 2:
				# Handle second arrow
				new_candidate2 = {'ref': candidate['ref'], 'query': candidate['query'], 'score': candidate['score']}
				arrow = arrows[1]
				a2 = a
				b2 = b
				# Move traceback pointer based on direction given by arrow
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
			# Move traceback pointer based on direction given by arrow
			if arrow == 'U':
				a = self.trace_up(a, candidate)
			elif arrow == 'D':
				a, b = self.trace_diagonally(a, b, candidate)
			elif arrow == 'H':
				b = self.trace_horizontally(b, candidate)  
		
		# Finish traceback if path hits first column or row and is not yet at the top left cell
		while a > start_a:
			a = self.trace_up(a, candidate)     
		while b > start_b:
			b = self.trace_horizontally(b, candidate)
		
		self.candidates.append(candidate)
		return candidate 


	def traceback_with_matches(self):
		'''
		Calls traceback() to find optimal alignment for segments that 
		are not already matched, to put together the overall alignment along
		with matched segments. 
		'''

		# If matches file is misshapen, quit
		if len(self.query_matches) != len(self.ref_matches):
			print('Please fix your matches file and rerun.')
			return

		global_ref = ''
		global_query = ''
		total_score = 0

		# Identify non matching region at start, if any
		if self.ref_matches[0]['start'] != 0:
			region_end_ref = self.ref_matches[0]['start']
			# print('Aligning 0—' + str(region_end_ref))
		
			if self.query_matches[0]['start'] != 0:
				region_end_query = self.query_matches[0]['start']
				# print('Aligning 0—' + str(region_end_query))
				# Submit NW request on this segment
				alignment = self.traceback(0, region_end_query, 0, region_end_ref, {'ref': '', 'query': '', 'score': 0})
				global_ref += alignment['ref'] 
				global_query += alignment['query'] 
				total_score += alignment['score']
				# print(alignment)
		
		# Iterate over matches
		for i, match in enumerate(self.query_matches): # just pick one species to start iterating over, since both species have the same amount of matches
			ref_segment = self.ref[self.ref_matches[i]['start'] : self.ref_matches[i]['end']]
			query_segment = self.query[match['start'] : match['end']]
			global_ref += ref_segment 
			global_query += query_segment 
			# print('Adding ' + str(self.ref_matches[i]['start']) + '—' + str(self.ref_matches[i]['end']))
			# print('Adding ' + str(self.query_matches[i]['start']) + '—' + str(self.query_matches[i]['end']))
			print(ref_segment)
			print(query_segment)
			total_score += (match['end'] - match['start']) * match_award # only need one pair since both segments are assumed to be the same length
			
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
			# print('Aligning ' + str(region_start_ref) + '—' + str(region_end_ref))
			# print('Aligning ' + str(region_start_query) + '—' + str(region_end_query))
			alignment2 = self.traceback(region_start_query, region_end_query, region_start_ref, region_end_ref, {'ref': '', 'query': '', 'score': 0})
			global_ref += alignment2['ref'] 
			global_query += alignment2['query'] 
			total_score += alignment2['score']

		return {'ref': str(global_ref), 'query': str(global_query), 'score': total_score}


	def run(self):
		'''
		Main method of executable.
		'''

		# Build score and traceback grids. (Don't know if it is faster to build a complete grid if there are matches or to build multiple smaller grids.)
		self.build_grids(self.ref, self.query)

		# If matches file is not given, 
		if self.args.matches is None:
			a = len(self.query)
			b = len(self.ref) 
			self.alignment = self.traceback(0, a, 0, b, {'ref': '', 'query': '', 'score': 0})
			print('Number of alignments found: ' + str(len(self.candidates)))
			print(self.candidates)
		else: 
			self.parse_matches(self.args.matches)
			self.alignment = self.traceback_with_matches()
			# print(self.alignment['ref'])
			# print(self.alignment['query'])
			# Ensure all letters have been captured
			# print(str(alignment['ref']).replace('-', '') == self.ref)
			# print(str(alignment['query']).replace('-', '') == self.query)
			# print(len(str(alignment['ref']).replace('-','')))

		# Print alignment to file
		filename = 'alignment_' + self.ref_obj.name + '_' + self.query_obj.name + '.txt'
		with open(filename, 'w') as output_file:
			ref_output = self.ref_obj.name + '\n' + self.alignment['ref'] + '\n\n' 
			query_output =  self.query_obj.name + '\n' + self.alignment['query'] + '\n\n'
			score_output = 'Alignment score: ' + str(self.alignment['score'])
			output = ref_output + query_output + score_output
			output_file.write(output)

	def permute(self):
		'''
		Permutes the amino acids in both sequences
		repeats the alignment 10,000 times, 
		and visualizes the score distribution.
		'''

		scores = []

		for i in range(0, 10000):
			# Print progress so we aren't wondering if the program froze
			if i % 100 == 0:
				print('... ' + str(i) + '/10000...')
			# Permute both sequences
			permuted_ref = ''.join(np.random.permutation(self.ref))
			permuted_query = ''.join(np.random.permutation(self.query))
			# print(permuted_ref)
			# print(permuted_query)
			
			# Find optimal alignment of permuted sequences 
			self.build_grids(permuted_ref, permuted_query)
			alignment = self.traceback(0, len(self.query), 0, len(self.ref), {'ref': '', 'query': '', 'score': 0})
			scores.append(alignment['score'])
		
		# Report score distribution with a histogram
		plot = plt.hist(x=scores, bins=10, facecolor='pink', alpha=0.5)
		plt.title('Score distribution for 10,000 alignments of randomly permuted ' + self.query_obj.name + ' and ' + self.ref_obj.name)
		plt.ylabel('Frequency')
		plt.xlabel('Alignment score')
		# Mark the alignment score of the original sequences
		plt.axvline(x=self.alignment['score'], color='gray', linestyle='dashed', linewidth=2)
		plt.legend(['Original sequences'])
		filename = 'histogram_' + self.ref_obj.name + '_' + self.query_obj.name + '.pdf'
		plt.savefig(filename)
		plt.show()
		# Save to PDF file 
		# https://stackoverflow.com/questions/38118510/append-page-to-existing-pdf-file-using-python-and-matplotlib


if __name__ == '__main__':
	match_award = 1
	gap_penalty = -2
	mismatch_penalty = -3
	max_alignments = 1 # How many optimal alignments we want to return --
	# Can set max_alignments = None, but this will lead to memory blowout w/ larger sequences
	# We are only interested in 1 optimal alignment for now. 

	# Create and run Needleman Wunsch executable 
	nw = Needleman_Wunsch_Executable(max_alignments, match_award, gap_penalty, mismatch_penalty)
	nw.run()
	nw.permute() 

	# This program assumes that one of the sequences is human and that its FASTA description contains 'human' somewhere.





