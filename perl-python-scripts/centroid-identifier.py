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
from progress.bar import Bar
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

def readSequences(seq_input):
    '''
        Function that reads an input sequence alignment file (FASTA)
        and returns a list of sequence strings

        Parameters: 
            input (FASTA): Amino acid sequence alignment in FASTA format

        Returns:
            seq_list (list): A list of protein sequence strings
    '''
    seqrecords = SeqIO.parse(seq_input, 'fasta')
    seq_list = list(set([str(seqrec.seq) for seqrec in seqrecords]))
    return seq_list


def generateSeq(reference_sequences):
    '''
        Function that takes in a list of aligned amino acid sequences and
        generates a new sequence based on input sequence alignment. It goes through
        each column in the alignment and randomly selects a row in that column 
        to choose a residue to be added to the new generated sequence.

        Parameters:
            reference_sequences (list): A list of amino acid sequence strings 
                                        obtained from an alignment file in 
                                        FASTA format
        Returns:
            sequence (Seq): A Biopython mutable Seq object (amino acid sequence) 
    '''
    sequence = ''
    # Go through each column in the alignment
    for column in range(len(reference_sequences[0])):
        # Select a random row in the column
        row = random.randint(0, len(reference_sequences) - 1)
        # Append the residue in that position to the new sequence
        sequence += reference_sequences[row][column]
    # Return a mutable (easier to mutate later) Seq object
    return Seq(sequence).tomutable()

def getVariableSites(sequences):
    '''
        Function that takes in a list of aligned amino acid sequence strings
        and returns a dictionary of variable column positions in the sequence
        alignment.

        Parameters:
            sequences (list): A list of aligned amino acid sequence strings

        Returns:
            variable_sites (dict): A dictionary of variable column positions
                                   in the sequence alignment.
                                   Key: variable site in the alignment
                                   Value: amino acids present at that site

    '''
    # Get the reference sequence in the alignment
    reference = sequences[0]
    variable_sites = []
    # Go through each column in the alignment
    for column in range(len(reference)):
        # Go through each row in the alignment
        for row in range(len(sequences[1:])):
            # Check if the amino acid at the position is not the same
            # as the amino acid in reference sequence
            if sequences[row][column] != reference[column]:
                # Add to variable sites list and iterate the next column
                variable_sites.append(column)
                break
    # Initiate a dictionary which will hold lists as values
    amino_acids = defaultdict(list)
    for sequence in sequences:
        for site in variable_sites:
            # Check if amino acid is already present in the list
            if sequence[site] not in amino_acids[site]:
                amino_acids[site].append(sequence[site])
    return amino_acids

def mutate(individual, positions):
    '''
        Function mutate individual sequences in a population.

        Parameters:
            individual (Seq): Amino acid sequence
            positions (dict): A dictionary for variable positions in the
                              amino acid sequence alignment. 
                              Key: position in a sequence
                              Value: Amino acids present in the position
        
        Returns:
            individual (Seq): Amino acid sequence (mutated)
    '''
    # Randomly choose a mutation point from provided values
    mutation_site = random.choice(list(positions.keys()))
    # Randomly select amino acid present at that position and mutate sequence
    individual[mutation_site] = random.choice(positions[mutation_site])
    return individual

def distance(string1, string2):
    '''
        A function that calculates Hamming distance between two given strings.

        Parameters:
            string1, string2 (str): Strings.
        
        Returns:
            distance (int): Calculated Hamming distance between two input
                            strings.
    '''
    distance = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            distance += 1
    return distance

def maxDistance(individual, sequences):
    '''
        Function that calculates the Hamming distance between an individual
        sequence and every sequence in the provided sequence list and returns
        the maximum distance value (value to be minimized)

        Parameters:
            individual (Seq): amino acid sequence
            sequences (list): A list of amino acid sequences from alignment
        
        Returns:
            max_distance (int): Maximum Hamming distance

    '''
    distances = {}
    for seq in sequences:
        distances[(str(seq), str(individual))] = distance(seq, individual)
    max_distance = max(distances.values())
    # Return max distance as tuple (DEAP requirement)
    return max_distance,

def centroidEvolution():
    '''
        A function that uses the Genetic Algorithm to generate
        one or more centroid sequences.
    '''

    # Read in sequences
    sequences = readSequences(args.input)

    # Variable sites
    variable_sites = getVariableSites(sequences)
    
    # Fitness and individual objects
    # # Minimizing the fitness value
    creator.create('FitnessMax', base.Fitness, weights=(-1.0,))
    creator.create('Individual', list, fitness=creator.FitnessMax)

    # Register individual and population objects
    toolbox = base.Toolbox()
    toolbox.register('attr_bool', generateSeq, sequences)
    toolbox.register('individual', tools.initRepeat, creator.Individual, toolbox.attr_bool, 1)
    toolbox.register('population', tools.initRepeat, list, toolbox.individual)
    toolbox.register('evaluate', maxDistance)
    toolbox.register('mutate', mutate)
    toolbox.register('select', tools.selTournament, tournsize=3)

    # List of centroid sequences
    centroids = []

    # Evolution
    for i in range(args.centroids):
        population = toolbox.population(args.population)

        fitnesses = [toolbox.evaluate(ind[0], sequences) for ind in population]
        for ind, fit in zip(population, fitnesses):
            ind.fitness.values = fit
        fits = [ind.fitness.values[0] for ind in population]

        # Best fitness before evolution
        best_before = min(fits)
        print('Evolving Centroid {}'.format(i + 1))
        bar = Bar('Generation', max=args.generations)  # Progress bar
        # Begin evolution
        for generation in range(1, args.generations):
            
            offspring = toolbox.select(population, len(population))
            offspring = list(map(toolbox.clone, offspring))
            
            for mutant in offspring:
                if random.random() < 0.7:
                    toolbox.mutate(mutant[0], variable_sites)
                    del mutant.fitness.values
                    
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = [toolbox.evaluate(ind[0], sequences) for ind in invalid_ind]
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            
            population[:] = offspring
            fits = [ind.fitness.values[0] for ind in population]
            bar.next()
        bar.finish()

        # Get the best sequence and append it to the centroid list
        best_seq = [ind for ind in population if ind.fitness.values[0] == min(fits)]
        if i < 9:
            centroid_seq = SeqRecord(best_seq[0][0].toseq(), id='C' + '0' + str(i + 1),
                    name='C' + '0' + str(i + 1), description='C' + '0' + str(i + 1))
        else:
            centroid_seq = SeqRecord(best_seq[0][0].toseq(), id='C' + str(i + 1),
            name='C' + str(i + 1), description='C' + str(i + 1))
        centroids.append(centroid_seq)

        # Print out the max distances before and after evolution
        best_after = min(fits)
        print('Centroid ' + str(i + 1) + ' --')
        print('  Max distance before: ', best_before)
        print('  Max distance after: ', best_after)
        print()

    SeqIO.write(centroids, args.output, 'fasta')

centroidEvolution()
