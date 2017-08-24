#!/usr/bin/env python

import ROOT
import math
from array import array

import os
import errno

def loader(infile, tree_name):
    chain = ROOT.TChain(tree_name)
    if "root" in infile:
        chain.Add(infile)
        return chain

    # infile is a list of input root files
    ncounter = 0
    with open(infile, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            chain.Add(line[:-1])
            ncounter += 1

    print "total events:",chain.GetEntries(),"in",ncounter,"files"
    return chain

def get_eff(a, b):
    eff = a*1.0/b
    error = eff * math.sqrt(1./a + 1./b)
    print "{:.4f} {:.4f}".format(eff,error)
    return eff,error

def get_eff_error(a, ae, b, be):
    eff = a*1.0/b
    error = eff*math.sqrt((ae/a)**2 + (be/b)**2)
    return eff,error

def column(matrix, i):
    return [row[i] for row in matrix]

def flatten(matrix):
    return [x for y in matrix for x in y]

def significance(s, b):
    if s < 0 or b < 0:
        return 0
    if abs(b) < 1E-6:
        return 0
    return math.sqrt(2*((s+b)*math.log(1+s/b) - s))

def read_np_list(np_list_name):
    np_list = []
    with open(np_list_name) as handle:
        for line in handle:
            if line.strip()[0] == '#':
                continue
            np_list.append(line[:-1].strip())
    
    print "Total nuisance parameters:",len(np_list)
    return np_list

def mkdir_p(path):
    try:
        os.makedirs(path)
        print path,"is created"
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
