#! python
# -*- coding: utf-8 -*-
"""
bzipscore-all.py

wrapper to score all the interaction pairs in a given fasta file
"""

def get_ids_from_fasta(fasta_file):
    "Returns all the ids (names) from a fasta file."
    from Bio import SeqIO

    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    ids = []
    for fasta in fasta_sequences:
        ids.append(fasta.id)
    
    return ids   


def get_pairs(ids):
    """
    Returns all possible (symetric) pairs from a list of ids.
    
    For example get_pairs([A B C]) returns:
      AA, BB, CC, AB, AC, CB""" 
    import itertools
    return itertools.combinations_with_replacement(ids, 2)
    
def get_interaction_string(ids):
    """
    Returns a string of all symetric pairs, comma separated. 
    This string is sutable for writing in a interactions file.
    """
    pairs = get_pairs(ids)
    s = ""
    for p in pairs:
      s = s + ",".join(p)+'\n'
    return s        

def get_interaction_lines(ids):
    for p in get_pairs(ids):
        yield ",".join(p)
            
import string
import random

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))

         

def write_interaction_file(ids, file_name=None):
    if file_name is None:
        #get temporary name, so nothing is overwritten
        file_name = id_generator(16) + '.int'

    #int_str = get_interaction_string(ids)
    # 'wb -- don't convert new lines lines'
    with open(file_name, 'w') as f:
        #f.write(bytes(int_str, 'ASCII'))
        for line in get_interaction_lines(ids):
            f.write(line+'\n')
        
    return file_name


def do_bzip_scoring(fasta_file, sfunc='complete', script = "bzipscore.pl"):
    """Writes the interaction file and runs the bzipscore. The result is printed to stdout."""
    import os
    import subprocess
    if not os.path.isfile(fasta_file):
        raise NameError("File "+fasta_file+" does not exist!")
    
    ids = get_ids_from_fasta(fasta_file)
    #print ids
    interactions_file = write_interaction_file(ids)
    
    out_str = subprocess.check_output(['perl',script, '--sfunc',sfunc, interactions_file, fasta_file ], universal_newlines=True) 
    
    print (out_str)
    os.remove(interactions_file)    
    return out_str

import argparse

parser = argparse.ArgumentParser(
                    description="""Wrapper to score all the interaction pairs in a given fasta file.
                    VERSION: 0.1
                    WARNING: All sequences should start at the f position!""")
parser.add_argument('fasta_file', type=str,
                   help='path to fasta file (file with sequences)')
parser.add_argument('--sfunc', choices=['complete','rfe','vinson_ce','fong_svm','bcipa'], default='complete',
                   help='function used for scoring the interactions (default is complete)')

args = parser.parse_args()
#print args.fasta_file
#print args.sfunc

do_bzip_scoring(args.fasta_file, args.sfunc);

