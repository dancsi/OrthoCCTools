"""Test for the solver.
Must be run from the root directory of the project"""

import os
from subprocess import call

def remove_output():
    if os.path.isfile("output.txt"):
        os.remove("output.txt")
        
def ensure_output():
    assert os.path.isfile("output.txt")

def load_output():
    with open ("output.txt", "r") as f:
        data=f.read()   
    return data
    
def load_file(file):
    with open(file, "r") as f:
        data=f.read()
    return data  

def load_set_file(file_name):
    """Loads a score file of comma delimited pairs into a set of sets"""
    str = load_file(file_name)
    return load_set_str(str) 

def load_set_str(str):
    lines = str.split("\n")
    res = list()
    for line in lines:
        line = line.strip()
        if (line=="") or (line[0]=="#"):
            continue
        sl = line.split(",")        
        res.append(list((sl[0],sl[1])))
    return res        
    
def pairlist_to_set(pairlist):
    """Converts a list of pairs to a set of sets"""
    res = set()
    for p in pairlist:
        res.add(tuple(sorted(p)))
    return res    

def pairs_equal(pairlist1, pairlist2):
    ps1=pairlist_to_set(pairlist1)
    ps2=pairlist_to_set(pairlist2)
    return ps1==ps2
    
def test_homo5():
    call("solver data/tests/homo5.out --binding-cutoff=-8.5 --nonbinding-cutoff=-7.5")
    ensure_output()   
    pairs = load_set_file("output.txt")     
    target_pairs = [["P1","P1"],
                    ["P2","P2"],
                    ["P3","P3"],
                    ["P4","P4"],
                    ["P5","P5"]]
    assert pairs_equal(pairs, target_pairs)   

def test_homoP1P5():
    call("solver data/tests/homoP1P5.out --binding-cutoff=-8.5 --nonbinding-cutoff=-7.5")
    ensure_output()   
    pairs = load_set_file("output.txt")     
    assert len(pairs)==4
    
def test_homo2hetero2():
    call("solver data/tests/homo2-hetero2.out --binding-cutoff=-8.5 --nonbinding-cutoff=-7.5")
    ensure_output()   
    pairs = load_set_file("output.txt")     
    target_pairs =[["P1","P2"],
                   ["P3","P4"],
                   ["P5","P5"],
                   ["P6","P6"]]    
    assert pairs_equal(pairs, target_pairs)