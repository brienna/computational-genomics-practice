from Bio import SeqIO
import numpy as np
import itertools
import os
import subprocess
import pandas as pd

class Neighbor_Joining_Tree(object):
    '''
    Given a multiple sequence alignment, constructs a neighbor joined tree.
    All equations are from https://en.wikipedia.org/wiki/Neighbor_joining
    Stores tree as either NEWICK or as a tab-delimited file with the following format:
    [ancestor_node, descendant_node, branch_length]
    '''
    
    def __init__(self):
        np.set_printoptions(suppress=True)
        self.D = None
        self.identifiers = None
        self.tree = np.zeros(shape=(1,3))
        self.sums = None
        self.n = 0
        self.node_ids = None
        self.internal_node_id_tracker = 0
        self.Q = None
        self.msa = None
        self.root_id = 0
    
    def initialize_other_variables(self):
        '''
        Initializes remaining variables, 
        after D has been initialized or uploaded.
        '''

        self.sums = self.D.sum(axis=0)
        self.n = len(self.D[0])
        self.node_ids = np.arange(1, len(self.D)+1, 1) # 1 based counting, not 0 based
        self.internal_node_id_tracker = self.n + 1
        self.root_id = self.n * 2 - 2
        
    def upload_D(self, filename):
        '''
        Uploads given tab-delimited distance matrix file.
        '''

        with open(filename, 'r') as file:
            D = pd.read_csv(file, delimiter='\t', index_col=0)
            self.identifiers = D.columns
            self.msa = D.columns # for now (fix later)
            self.D = np.asarray(D.rename_axis().values)

    def print_D(self, filename):
        '''
        Saves pairwise dissimilarities to file.
        '''

        distances_df = pd.DataFrame(self.D, columns=self.identifiers) 
        distances_df.insert(0, 'id', self.identifiers) # Add row index
        distances_df = distances_df.reset_index(drop=True).set_index('id')
        distances_df.to_csv(filename, sep='\t')
            
    def upload_msa(self, fna):
        '''
        Uploads multiple sequence alignment.
        '''

        self.msa = []
        self.identifiers = []
        # Parse sequences and identifiers
        fna_obj = list(SeqIO.parse(fna, 'fasta'))
        for i, seq in enumerate(fna_obj):
            self.msa.append(seq.seq)
            self.identifiers.append(seq.id)
            
    def set_msa(self, msa):
        '''
        Sets given MSA.
        '''
        self.msa = msa
    
    def set_identifiers(self, identifiers):
        '''
        Sets given identifiers.
        '''
        self.identifiers = identifiers
            
    def parse_D(self):
        '''
        Calculates the distance matrix from the MSA. 
        '''

        distances = []
        # For each pair, calculate their pairwise dissimilarity
        for pair in itertools.product(self.msa, repeat=2): # repeat=2? ... and we are doing too much work here but its ok for now, can just do half of the work and flip the matrix across the diagonal
            distances.append(self.calculate_dissimilarity_score(pair[0], pair[1]))
        
        # Reshape distances array into matrix based on number of sequences
        self.D = np.reshape(distances, (len(self.msa), len(self.msa)))
    
    def calculate_dissimilarity_score(self, a, b):
        '''
        Calculates pairwise dissimilarity of sequences a, b.
        If letters do not match, add 1. Otherwise add 0. 
        Sequences a and b are the same length.
        Returns dissimilarity score as float. 
        '''

        score = 0.00
        length = len(a)

        for i, char_a in enumerate(a):
            char_b = b[i]
            if char_a != char_b:
                score += 1

        return score/length
    
    def calculate_q_value(self, distance, sums_a, sums_b):
        '''
        Calculates value of a cell for Q matrix.
        n = number of sequences left in D
        distance — dissimilarity between sequences a, b
        sums_a = sum for all dissimilarities between a and the other sequences 
-       sums_b = sum for all dissimilarities between b and the other sequences 
        '''
        return (self.n - 2) * distance - sums_a - sums_b
            
    def calculate_Q(self):
        '''
        Calculates Q matrix, which uses the following variables:
        '''

        # print('Calculating Q matrix...')
        Q = np.zeros_like(self.D, dtype=float)
        
        for i, row in enumerate(self.D):
            sums_a = self.sums[i] # sum of distances from a to all other nodes
            for j, distance in enumerate(row):
                if i == j:
                    Q[j, i] = np.inf # assign diagonal to positive infinity so it's never the minimum
                else:
                    sums_b = self.sums[j] # sum of distances from b to all other nodes
                    Q[j, i] = self.calculate_q_value(distance, sums_a, sums_b)
        
        return Q
    
    def calculate_branch_length(self, distance, sums_f, sums_g):
        '''
        Calculates branch length between f,g.
        '''
        return distance/2 + abs((sums_f - sums_g))/(2*(self.n-2))

    
    def calculate_next_branch_length(self, distance, first_branch_length):
        '''
        Calculates branch length given the first branch length, 
        to be used as the second branch length.
        '''

        return abs(distance - first_branch_length) 
        
    
    def get_branches(self, f, g, h):
        '''
        Sets up variables and calculations to obtain branch lengths
        between f, g, and h. Usually it is only f, g, until we reach the root.
        '''

        # Get node ids for the pair
        f_id = self.node_ids[f]
        g_id = self.node_ids[g]
        
        # Track new, ancestral node id
        u_id = self.internal_node_id_tracker
        self.internal_node_id_tracker += 1
        
        # Calculate lengths of branches joining f and g to u
        distance_fg = self.D[f, g]
        delta_fu = self.calculate_branch_length(distance_fg, self.sums[f], self.sums[g])
        delta_gu = self.calculate_next_branch_length(distance_fg, delta_fu)
        
        # Format and return branches 
        branch_fu = [u_id, f_id, delta_fu]
        branch_gu = [u_id, g_id, delta_gu]
        branch_hu = None
        
        # If h has been specified, we need to get the third branch
        if h:
            h_id = self.node_ids[h]
            distance_fh = self.D[f, h]
            delta_hu = self.calculate_next_branch_length(distance_fh, delta_fu)
            branch_hu = [u_id, h_id, delta_hu]
            
        # Insert new node (if it goes above, it interferes with h assignment)
        self.node_ids = np.insert(self.node_ids, 0, u_id, axis=0)
        
        return branch_fu, branch_gu, branch_hu
    
    def calculate_uk(self, fk, gk, fg):
        '''
        Calculates and returns distance of node u to node k. 

        dist_fk — distance of node f to node k
        dist_gk — distance of node g to node k
        dist_fg — distance of node f to node g

        where f and g are members of the pair just joined.
        '''
        return (fk + gk - fg)/2
    
    def update_D(self, f, g):
        '''
        Updates distance matrix to include 
        '''

        row_u = np.zeros_like(self.D[0]) 
        for k,_ in enumerate(self.D):
            row_u[k] = self.calculate_uk(self.D[f, k], self.D[g, k], self.D[f, g])

        # Update D matrix with these distances
        column_u = row_u[:, np.newaxis] 
        column_u = np.vstack([0, column_u]) # Add its own diagonal value
        self.D = np.vstack((row_u, self.D))
        self.D = np.hstack((column_u, self.D))
        
    def preorder_traversal(self, root):
        '''
        Performs preorder traversal on the tree, starting at given root node.
        NOTE: Only works if we use the node selected as root in neighbor joining.
        Preorder traversal is a DFS search where we 
            1. Visit the root.
            2. Traverse the left subtree, i.e., call Preorder(left-subtree)
            3. Traverse the right subtree, i.e., call Preorder(right-subtree) 
        '''

        traversed_tree = []
        children = np.where(self.tree[:,0] == root)[0] 

        # If root has children,
        if children.size > 0:
            # Traverse its children
            for child in children: 
                traversed_tree.append(self.tree[child])
                traversed_tree += self.preorder_traversal(self.tree[child][1])
    
        return traversed_tree
    
    def tree_to_newick(self, root):
        '''
        Format tree as NEWICK, which involves postorder traversal.
        '''

        newick_tree = ''
        children = np.where(self.tree[:,0] == root)[0]

        # If we have children, traverse them
        if children.size > 0:
            newick_tree += '('

            for i, node in enumerate(children):
                # Get its id and branch length to root
                node_id = int(self.tree[node][1])
                branch_length = self.tree[node][2]

                newick_tree += self.tree_to_newick(node_id)
                # If node is internal, 
                if node_id > len(self.msa):
                    newick_tree += ':' + str(branch_length) + ','
                else:
                    taxonomic_id = self.identifiers[node_id - 1]
                    newick_tree += str(taxonomic_id) + ':' + str(branch_length) + ','

            newick_tree = newick_tree[:-1] + ')' # removes final comma

        return newick_tree
            
    def run(self):
        '''
        Executes NJT. 
        '''

        print('Joining tree...')
        self.initialize_other_variables()
        
        while len(self.D) >= 3:
            # Calculate the join score of each pair
            self.Q = self.calculate_Q()

            # Find the pair f,g with the minimum join score 
            if len(self.D) == 3:
                # If on the last iteration, we have 3 nodes to join
                f, g, h = 0, 1, 2
            else:
                pairs = np.where(self.Q == np.amin(self.Q))
                f = pairs[1][0]
                g = pairs[0][0]
                h = None
                
            
            # Get branch lengths for f,g to their ancestral node u
            branch_fu, branch_gu, branch_hu = self.get_branches(f, g, h)
            self.tree = np.vstack((self.tree, branch_fu))
            self.tree = np.vstack((self.tree, branch_gu))
            if branch_hu:
                self.tree = np.vstack((self.tree, branch_hu))
                break # since we're on the last iteration, we can quit
            
            # Update D to include distances between u and every other node k
            self.update_D(f, g)

            # Remove f,g from D 
            self.D = np.delete(self.D, [f+1, g+1], axis=1)
            self.D = np.delete(self.D, [f+1, g+1], axis=0)
            
            # update variables
            self.sums = self.D.sum(axis=0)
            self.node_ids = np.delete(self.node_ids, [f+1, g+1])
            self.n = len(self.D)
            
        self.tree = np.delete(self.tree, 0, axis=0) # remove that initialized beginning [0,0,0]
    
    def visualize_tree(self, newick_file):
        '''
        Visualizes the tree using passed file, by calling subprocess to run an external R script.
        (Can be newick or edges, need to update accordingly.)
        '''

        # cmd_edges = 'Rscript hw3-plot-edges.r ' + edges_file + ' tip-labels.txt'
        cmd_newick = 'Rscript hw3-plot-newick.r ' + newick_file + ' tip-labels.txt'
        print(cmd_newick)
        proc = subprocess.Popen(cmd_newick,shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        proc.communicate()

    def print_tree(self, additional_id):
        '''
        Prints tree in edges format [ancestor, descendant, branch_length],
        NEWICK format, and visualizes it using an external R script.
        '''

        print('Printing tree...')

        # Create edges/newick directories if they don't exist
        edges_dir = 'edges/'
        newick_dir = 'newick/'
        if not os.path.exists(edges_dir): 
            os.makedirs(edges_dir)
        if not os.path.exists(newick_dir):
            os.makedirs(newick_dir)
            
        # Make and save edges text file
        traversed_tree = self.preorder_traversal(self.root_id)
        edges_file = edges_dir + 'edges' + str(additional_id) + '.txt'
        np.savetxt(edges_file, traversed_tree, fmt='%i\t%i\t%1.10f')
        
        # Make and save newick tree
        newick_tree = self.tree_to_newick(self.root_id) + ';'
        newick_file = newick_dir + 'tree' + str(additional_id) + '.txt'
        with open(newick_file, 'w') as text_file:
            text_file.write(newick_tree)

        # Only visualize original tree (can change)
        if not additional_id:
            self.visualize_tree(newick_file)

    
