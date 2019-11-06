from njt import Neighbor_Joining_Tree
from boot import Bootstrapper
import argparse

def make_arg_parser():
	parser = argparse.ArgumentParser(prog='main.py',
						  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-m","--msa",
				  default=argparse.SUPPRESS,
				  required=True,
				  help="Path to multiple sequence alignment FASTA") 

	return parser

if __name__ == '__main__':
	parser = make_arg_parser()
	args = parser.parse_args() 

	# Create neighbor joining tree based on given alignment
	njt = Neighbor_Joining_Tree() 
	njt.upload_msa(args.msa)
	njt.parse_D()
	njt.print_D('pairwise_dissimilarity.txt') # HW3.1
	njt.run() # HW3.2
	njt.print_tree('') # HW3.3 and #HW3.4

	# Run bootstrapping 100 times (BONUS)
	b = Bootstrapper(njt)
	b.run_bootstrap(100)
	b.partition_all()
	b.calculate_bootstrap_confidences()
	b.print_bootstrap('nj-solution-bootstrap.txt')