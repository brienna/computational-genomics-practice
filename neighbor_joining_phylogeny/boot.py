from njt import Neighbor_Joining_Tree
import random
import numpy as np
import subprocess

class Bootstrapper(object):
    '''Performs inferences of an original neighbor joining tree using bootstrap samples.'''

    def __init__(self, njt):
        self.njt = njt # original tree
        self.partitions = []
        self.frequencies = None

    def run_bootstrap(self, num_trees):
        '''
        Runs bootstrap for given number of trees, randomly arranging the original MSA 
        columns before creating each new tree. Since this is repeated sampling with replacement, 
        some columns may appear more than once in a bootstrap sample. 
        '''

        random.seed(112)

        for i in range(0, num_trees):
            print('DOING BOOTSTRAP ID: ' + str(i))
            bootstrap = []
            # Create bootstramp sample of original MSA columns
            for j in range(0, len(self.njt.msa[0])):
                bootstrap.append(random.randint(0,len(self.njt.msa[0])-1))
            permuted_msa = []
            for seq in self.njt.msa:
                seq = ''.join([list(seq)[i] for i in bootstrap])
                permuted_msa.append(seq)

            # Create the bootstrapped tree 
            njt2 = Neighbor_Joining_Tree()
            njt2.set_msa(permuted_msa)
            print(njt2.msa[0])
            njt2.set_identifiers(self.njt.identifiers)
            njt2.parse_D()
            njt2.run()
            njt2.print_tree(i)
            
    def partition(self, nodes):
        '''
        Partitions one tree by internal node and all of its descendant tips.
        Returns partitions as dictionary {internal_node : [descendant_tips]}
        '''

        p = {}

        for node in nodes:
            node = node.split('\t')

            # If internal node, track its immediate children 
            root = int(node[0])
            if root > len(self.njt.msa):
                descendant = int(node[1])
                if root in p:
                    p[root].append(descendant)
                else:
                    p[root] = [descendant]

        # Rearrange partitions to track all descendants
        for node, children in p.items():    
            # Check each child
            for child in children:
                # If it is also a root itself
                if child in p:
                    # Track its children as descendants of original root
                    p[node] += p[child]

            # Filter out children that are internal nodes, keeping only tips
            p[node] = list(filter(lambda x: x < 61, p[node]))
            p[node].sort()

        return p

    def partition_all(self):
        '''
        Reads in edges files and calls partition() to
        partition original tree along with all bootstrapped trees.
        '''

        # Partition original tree
        with open('edges/edges.txt', 'r') as sample_file:
            nodes = [line.rstrip('\n') for line in sample_file.readlines()]
            self.partitions.append(self.partition(nodes))
            # Note: Important to partition this first to maintain original order of nodes in dict (>Python 3.6 maintains dict order)

        # Partition bootstrapped trees
        for i in range(0, 100): # should change 100 to num_trees 
            with open('edges/edges' + str(i) + '.txt', 'r') as sample_file:
                nodes = [line.rstrip('\n') for line in sample_file.readlines()]
                self.partitions.append(self.partition(nodes))   

    def calculate_bootstrap_confidences(self):
        '''
        Calculate bootstrap confidences for internal nodes of each tree.
        '''

        self.frequencies = np.zeros(shape=(len(self.partitions[0]),))

        # Get partitions for original tree
        original_p = self.partitions[0] 

        # Compare original partition with bootstrapped partitions, counting confidence 
        # for each internal node (the fraction of bootstrap trees in which it has the 
        # exact same partition as the original tree)
        for p in self.partitions[1:]:
            # Loop through the original partition
            for i, node in enumerate(original_p):
                # Compare its partition with the bootstrapped partition
                if original_p[node] == p[node]:
                    self.frequencies[i] += 1

        # Turn frequencies into fraction of bootstrap trees
        self.frequencies[:] = [x / 100 for x in self.frequencies]

    def print_bootstrap(self, filename):
        '''
        Print to given filename the bootstrap confidences for internal nodes.
        '''

        # Print frequencies to file
        with open(filename, 'w') as bootstrap_file:
            for freq in self.frequencies:
                bootstrap_file.write(str(freq) + '\n')


