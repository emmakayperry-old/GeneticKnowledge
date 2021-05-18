##############################
## Genetic Sequence Checker ##
##############################

import os
from Bio import SeqIO, Entrez

########################################################################
###                              CONFIG                              ###
########################################################################

# set filepath and filename to the path and name of the relevant clingo 
# file (expects GeneticKnowledgebase.lp in the current working directory by default)
FILEPATH, FILENAME = '', 'GeneticKnowledgebase.lp'

# sets the name of the predicate returning the answer for Q3 part 1
ANS1_PRED_NAME = 'answer'

# store length of the predicate name
ANS1_PRED_LENGTH = len(ANS1_PRED_NAME)

# sets the name of the predicate returning possible sequences
PRED_NAME = 'seqResult'

# store length of the predicate name
PRED_LENGTH = len(PRED_NAME)

# sets the accession & email for accessing the database
EMAIL = '' #insert your email here
ACCESSION = 'MT072688.1' #set the accession of the relevant sequence
DATABASE = 'nucleotide' #set the database for the relevant sequence
RETTYPE = 'gb' 
RETMODE = 'text'

########################################################################

# download the record from genbank to current working directory
gbk_file = ACCESSION + '.gbk'
if not os.path.isfile(gbk_file):
    # Downloading...
    print('Saving genbank file to current working directory...')
    Entrez.email = EMAIL
    net_handle = Entrez.efetch(db=DATABASE, id=ACCESSION, 
                               rettype=RETTYPE, retmode=RETMODE)
    out_handle = open(gbk_file, 'w')
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print('File saved!')

# load the sequence associated with the accession
print('Parsing genbank sequence...')
record = SeqIO.read(gbk_file, 'genbank')
print('Sequence parsed!')

# get output from run of clingo module by opening a pipe to the CLI
# and executing the file with clingo, then reads the stream to get
# the output from the clingo run
print('Running clingo ' + FILENAME + '...')
stream = os.popen('clingo ' + FILEPATH + FILENAME)
output = stream.read()
print('Successfully executed clingo ' + FILENAME + '!')

print('Adding clingo output to search list...')
# initialize a list to store possible seqs
possible_seqs = []
found_seqs = {}

# iterate through the output to add possible sequences to the list
for string in  output.split(' '):
    # if it's the correct predicate result
    if string[:PRED_LENGTH] == PRED_NAME:
        # then add the possible sequence without the leading
        # predicate or leading and trailing parentheses or quotation
        # marks to the list of possible sequences
        possible_seqs += [string[PRED_LENGTH + 2:-2]]
    elif string[:ANS1_PRED_LENGTH] == ANS1_PRED_NAME:
        answer = string.split('\n')[0]
print('Output successfully added!')

print('Searching for possible sequences in ' + gbk_file + '...')
for seq in possible_seqs:
    position = record.seq.find(seq) + 1 # add 1 to the position to reflect the genbank position
    if position > 0:                    # conventions rather than pythonic indexing
        found_seqs[position] = seq

print('Search finished!')
print('\n===================')
print('A1 Q3 Part 1 Answer')
print('===================')
print(answer)

print('\n====================')
print('A1 Q3 Part 2 Answer:')
print('====================')
for position in found_seqs.keys():
    print('Sequence found! Position: ' + str(position) + '\nSequence: ' + found_seqs[position])