"""
Centroid Identifier

This program takes in a FASTA file containing aligned sequences as its
input and outputs a FASTA file with one or more centroid sequences.

Usage: centroid-identifier.py input output [-p POPULATION] [-g GENERATIONS] [-c CENTROIDS] 

Positional arguments:
    input: input file in FASTA format
    output: output file in FASTA format

Optional arguments:
    -p: Population size (default = 100)
    -g: Number of generations (default = 50)
    -c: Number of centroids to evolve (default = 1)
    -h: Show help message
"""

# Import modules
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from deap import creator, base, tools
from collections import defaultdict
import random

# Command line arguments
parser = argparse.ArgumentParser(description='Centroid identifier')
parser.add_argument('input', type=str, help='Input file (FASTA)')
parser.add_argument('output', type=str, help='Output file (FASTA)')
parser.add_argument('-p', '--population', type=int, default=100, 
                     help='Population size (default = 100)' )

parser.add_argument('-g', '--generations', type=int, default=50,
                    help='Number of generations (default = 50)')

parser.add_argument('-c', '--centroids', type=int, default=1,
                    help='Number of centroids to evolve (default = 1)')

args = parser.parse_args()

# Read in sequences
handle = open(args.input, 'r')
seqrecords = list(SeqIO.parse(handle, 'fasta'))
seq_list = list(set([str(seqrec.seq) for seqrec in seqrecords]))
handle.close()

# Function to generate individual sequences
def generate_seq():
    sequence = ''
    for column in range(len(seq_list[0])):
        row = random.randint(0, len(seq_list) - 1)
        sequence += seq_list[row][column]
    return Seq(sequence).tomutable()

# Fitness and individual objects
creator.create('FitnessMax', base.Fitness, weights=(-1.0,))   # Minimizing the fitness value
creator.create('Individual', list, fitness = creator.FitnessMax)

# Register individual and population objects
toolbox = base.Toolbox()
toolbox.register('attr_bool', generate_seq)
toolbox.register('individual', tools.initRepeat, creator.Individual, toolbox.attr_bool, 1)
toolbox.register('population', tools.initRepeat, list, toolbox.individual)

# Get mutation positions
indexes = []
for seq_1 in seq_list:
    for seq_2 in seq_list:
        if seq_list.index(seq_1) < seq_list.index(seq_2):
            for i in range(len(seq_1)):
                if seq_1[i] != seq_2[i]:
                    indexes.append(i)
indexes = list(set(indexes))

# Dictionary of amino acids and their positions
# Key: position
# Value: amino acids found at that position
amino_acids = defaultdict(list)
for seq in seq_list:
    for idx in indexes:
        amino_acids[idx].append(seq[idx])
        amino_acids[idx] = list(set(amino_acids[idx]))

#Mutation function
def mutation(ind):
    mpoint = random.choice(indexes)
    ind[0][mpoint] = random.choice(amino_acids[mpoint])
    return ind

# Helper function
def count_diff(string1, string2):
    count = 0
    for idx, char in enumerate(string1):
        if char != string2[idx]:
            count += 1
    return count

# Fitness function
# Returns maximum distance between pairs
def eval_distance(ind):
    distances = {}
    for seq in seq_list:
        distances[(str(seq), str(ind[0]))] = count_diff(seq, ind[0])
    return max(distances.values()), # Returned as tuple (DEAP requirement)


# Register the genetic operators and fitness function
toolbox.register('evaluate', eval_distance)
toolbox.register('mutate', mutation)
toolbox.register('select', tools.selTournament, tournsize=3)
MUTPB = 0.7  # the probability for mutating an individual

# List of centroid sequences
centroids = []

# Begin evolution

for i in range(args.centroids):
    pop = toolbox.population(args.population)     # Create population
    fitnesses = list(map(toolbox.evaluate, pop))  # Evaluate individuals
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    fits = [ind.fitness.values[0] for ind in pop] # List of fitnesses

    # Best fitness before evolution
    best_before = min(fits)

    # Evolution
    for generation in range(1, args.generations):
       
        # Select the next generation of individuals
        pop2 = toolbox.select(pop, len(pop))
        # Clone selected individuals
        pop2 = list(map(toolbox.clone, pop2))

        # Apply mutation
        for mutant in pop2:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                # Invalidate fitness of modified offspring
                del mutant.fitness.values
        
        # Reevaluate fitnesses
        invalid_ind = [ind for ind in pop2 if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Replace the old population with the offspring population
        pop[:] = pop2

        # Gather all the fitnesses in one list
        fits = [ind.fitness.values[0] for ind in pop]
    
    # Get the best sequence and append it to the centroid list
    best_seq = [ind for ind in pop if ind.fitness.values[0] == min(fits)]
    if i < 9:
        centroid_seq = SeqRecord(best_seq[0][0].toseq(), id='C' + '0' + str(i + 1),
                                 name='C' + '0' + str(i + 1), description='C' + '0' + str(i + 1))
    else:
        centroid_seq = SeqRecord(best_seq[0][0].toseq(), id= 'C' + str(i + 1),
                                 name='C' + str(i + 1), description = 'C' + str(i + 1))
    centroids.append(centroid_seq)

    # Print out the max distances before and after evolution
    best_after = min(fits)
    print('Centroid ' + str(i + 1) + ' --')
    print('  Max distance before: ', best_before)
    print('  Max distance after: ', best_after)

# Write a new FASTA file
SeqIO.write(centroids, args.output, 'fasta')





