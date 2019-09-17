# Usage: python homework1.py -h
import sys
import os
import argparse
import collections
import re
import subprocess


def make_arg_parser():
	parser = argparse.ArgumentParser(prog='template-python3.py',
						  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', action='version', version='%prog 2.0')
	parser.add_argument("-q","--query",
					  default=argparse.SUPPRESS,
					  required=True,
					  help="Path to query fasta [required]") 
	parser.add_argument("-r","--ref",
					  default=argparse.SUPPRESS,
					  required=True,
					  help="Path to reference fasta [required]")
	parser.add_argument("-t","--taxonomy",
					  default=None,
					  required=True,
					  help="Path to taxonomy file [required]")
	parser.add_argument("-o","--output",
					  default=None,
					  required=True,
					  help="Path to output file [required]")
	parser.add_argument("-c","--command",
					  default='./burst',
					  help="Path to BURST command") 
	parser.add_argument("-V","--verbose",
					  action="store_true",
					  help="Verbose output")
	return parser


def run_burst(query, ref, taxonomy, output, burst_cmd='./burst', verbose=False):
	'''Thread worker function, runs BURST to search query sequences against reference sequences.'''
	# Concatenate BURST command based on expected order, which is
	# ./burst -r ref.fna -q query.fna --taxonomy taxonomy.txt -o output.txt
	cmd = ' '.join([burst_cmd, '-r', ref, '-q', query, '--taxonomy', taxonomy, '-o', output])
	return run_command(cmd, verbose=verbose)


def run_command(cmd, verbose=False):
	'''Runs the given command and returns return value and output.'''
	if verbose:
		print(cmd)
	proc = subprocess.Popen(cmd,shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	return_value = proc.returncode
	return return_value, stdout, stderr


def analyze(output_file):
	'''Analyzes output file.'''

	# Initialize variables
	percent_similarity_threshold = 97.0
	matches = 0  # that are above percent similarity threshold
	total = 0  # intermediate variable to help calculate average percent similarity 
	species_list = []
	species_regex = 's__(\\w+)'
	# s__					species prefix
	# (\\w+)				capture the species name
	percentage_regex = '^(?:[^\t]+\t){2}([\\d|\\\\.]*)'
	# ^                     beginning of line
	#  (?:                  start non capture group
	#    [^\t]+             1 or more character that is not tab 
	#      \t               a tabulation
	#        ){2}           2 tabs must appear 
	#          ([\d|\\.]*)  capture the percentage

	# Read output file line by line (without loading all of it into memory)
	with open(output_file, 'r') as infile:
		for line in infile:
			# Grab percentage similarity & species name
			percentage = float(re.search(percentage_regex, line)[1])
			species = re.search(species_regex, line)

			# Add species name to list, if specified
			# NOTE: Unclear if homework wants to know the most common bacterial species 
			# of matches at ≥97%, or the full query set. Either way, we get the same answer.
			if species:
				species_list.append(species[1])

			# If percentage similarity is above preassigned threshold, 
			is_above_threshold = percentage > percent_similarity_threshold 
			if is_above_threshold:
				# Update matches & total_percentage tallies
				matches += 1
				total += percentage

	# Print the number of sequences that had a match in the database at ≥97%
	# as a fraction of the original input query sequences
	total_original_query = subprocess.check_output('grep -c ">" query.fna', shell=True).decode('utf-8').strip()
	print('\na. Fraction of sequences that had a match in the database at ≥97%: ' + str(matches) + 
		'/' + str(total_original_query))

	# Find and print the most common species in the query set
	counter = collections.Counter(species_list)
	most_common_species = counter.most_common(1)[0][0]
	print('b. Most common bacterial species in the query set: ' + most_common_species)

	# Calculate and print the average percent similarity
	average = total/matches
	print(('c. Average percent similarity: {:.5f}').format(average))


if __name__ == '__main__':
	parser = make_arg_parser()
	args = parser.parse_args() 

	# Run BURST, logging result to console
	output_file = args.output
	result = run_burst(args.query, args.ref, args.taxonomy, output_file,
					   args.command, args.verbose)
	print(result)

	# Analyze BURST output file
	analyze(output_file)

	




