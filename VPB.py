#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
VANAPRABHAVA v.1.5
Created on Wed Dec 23 17:33:31 2015
@author: Remi Tournebize
@contact: remi (dot) tournebize (at) gmail (dot) com
@sotd: "Sabre a finances, corne de ma gidouille, madame la financiere, j'ai des oneilles pour parler et vous une bouche pour m'entendre."

@v1.2: more profiling than v1.1 resulting in a significant gain of speed
@v1.3_ms: few debuggings on the SFS (duplicates)
@v1.4: debugged issues related to time decimal precision in Newick strings + modified force_ultrametricity module
@v1.5 31122015: changed binomial to poisson random (large branch lengths throw an error on C long type for binomial)
"""

from ete3 import Tree
import numpy
import time
import sys
import argparse
import os.path

# perspectives: add a protracted speciation mode
# should need some more testing if **multiple** coalescent trees (beta...)

# to speed up runtime we made an extensive use of the set()
# data structure. It implies that THERE MUST NOT BE ANY DUPLICATE
# TIP LABELS IN THE INPUT TREE!

# assumes coalescent, ie ultrametric trees


################
## PARAMETERS ##
################

# set argument possibilities
parser = argparse.ArgumentParser(description='options to VPB algorithm')

parser.add_argument('n', type=str,
                    action="store", help="input Newick-formatted democoalescent tree")
parser.add_argument('-m', '--mutationRate', type=float, required=True,
                    action="store", help="mean mutation rate per lineage per generation")
parser.add_argument('-s', '--scalingFactor', type=float, required=True,
                    action="store", help="a factor to rescale branch lengths (set it to 1 if you do not need to rescale),  should be mandatory for MS input trees which are scaled by No")
parser.add_argument('-ophylo', '--PHYLOoutput', type=str, required=True, 
                    action="store", help="path for the output phylogeny")
parser.add_argument('-osfs', '--SFSoutput', type=str, required=True, 
                    action="store", help="path for the output SFS")
parser.add_argument('-I', '--islands', type=int, nargs="+", default=-9,
                    action="store", help="Number of individuals in each deme separated by a comma, must be in the same order as in your MS command")
parser.add_argument('--forceUltrametric',
                    action="store_true", help="add this option to force the input tree to be ultrametric (useful if you rescale branch lengths)")
parser.add_argument('--drawTrees',
                    action="store_true", help="add this option to draw the analyzed trees in external files (same dir as PHYLOoutput)")
parser.add_argument('--seed', type=int, default=-1,
                    action="store", help="random seed")
parser.add_argument('-q', '--quiet',
                    action="store_true", help="add this option to prompt minimal messages")

# parse actual arguments
args = parser.parse_args()

# export arguments
t = args.n #ok
if os.path.isfile(t) == False:
	sys.exit("The input tree file does not exist apparently!")
mu = args.mutationRate #ok
scale = args.scalingFactor #ok
ophylo = args.PHYLOoutput #ok
osfs = args.SFSoutput #ok
force_ultrametric = args.forceUltrametric #ok
plot_trees = args.drawTrees #ok
seed = args.seed #ok
quiet = args.quiet #ok

# prompting
if not quiet: print "===================================================="
if args.islands != -9:
    ms_islands = args.islands #ok
    ms_input = True #ok
    if not quiet: print "Number of islands     "+str(len(ms_islands))
    if not quiet: print "Island deme sizes:    "+", ".join("Isl_"+str(x)+": "+str(ms_islands[x]) for x in range(len(ms_islands)))
else:
    ms_input = False
    sys.exit("Sorry, exiting... VPB does not know how to deal with input trees other than MS for the moment! should be fixed soon!")
if not quiet:
	print "Scaling factor        "+str(scale)
	print "Mutation rate         "+str(mu)
	print "Force ultrametric     "+str(force_ultrametric)
        print "...................................................."
	print "Reading               "+t
	print "PHYLO output          "+ophylo
	print "SFS output            "+osfs
	print "Draw trees            "+str(plot_trees)
	if seed != -1:
		print "Random seed           "+str(seed)
		numpy.random.seed(seed)
	print "===================================================="



###############
## FUNCTIONS ##
###############

def ubranch_mutation(node, mu):
    lambd = node.dist * mu
    rb = numpy.random.poisson(lambd)
    if rb >= 1:
        return True
    else:
        return False
    
def leaf_species_dictionary(t):
    SPdict = {}
    SPcounts = {}
    for leaf in t:
        sp = leaf.sp
        if sp in SPdict:
            SPdict[sp].append(leaf)
            SPcounts[sp] += 1
        else:
            SPdict[sp] = [leaf]
            SPcounts[sp] = 1
    return SPdict, SPcounts
    
def range_by_length(start, end, length):
    step = ( end - start ) / (1.*length)
    L = list(numpy.arange(start, end, step))
    for l in range(0, len(L)):
        L[l] = int(round(L[l]))
    return L
    
def get_deme(tip_id, ms_islands):
    tip_id = int(tip_id)
    pop = numpy.where(tip_id <= ms_islands)
    if len(pop[0]) == 0:
        sys.exit("1. The MS island structure you provided through -I does not fit the number of tips in the demography. The option might have been mispecified. Check it out!")
        return(None)
    else:
        return(pop[0][0] + 1)



############
## SCRIPT ##
############
  
init = time.time()
sys.stdout.write('[')

#======================================================#
# READ Newick-formatted tree
t = Tree(t, format=5)
sys.stdout.write('R') # Read


#======================================================#
# SET at once: "names" of inner nodes && "species" attributes
if ms_input:
    ms_islands = numpy.array(ms_islands)
    ms_islands = numpy.cumsum(ms_islands)
    nPops = len(ms_islands)

innode = 0
largest_id = 0
nIndsORI = 0
for node in t.traverse():
    node.add_features(sp=1)
    if not node.is_leaf():
        node.name = "iN%d" %innode
        innode += 1
    else:
        nIndsORI += 1
        if ms_input == True:
           inn = int(node.name)
           if inn > largest_id: largest_id = inn
           pop = get_deme(inn, ms_islands)
           node.name = str(pop)+"_"+node.name
    if not node.is_root():
        node.dist *= scale

if (ms_input == True) and (ms_islands[-1] > largest_id):
    sys.exit("2. The MS island structure you provided through -I does not fit the number of tips in the demography. The option might have been mispecified. Check it out!")

# set here the input tree to ultrametric?
# absolute ultrametricity does not matter for topological manipulations.
# It could possibly affect the speciational spreading but, in general, timescales are so big
# that a lack in ultrametricity due to F * precision would eventually have insignificant impacts

if plot_trees: t.render(ophylo+"_1ORIGINAL.png", w=183, units="mm")
sys.stdout.write('T') # Transform


#======================================================#
# SPREAD speciation events
spID = 0
for node in t.traverse("preorder"):
    # if leaf: NO speciation
    if not node.is_leaf():
        umut = ubranch_mutation(node, mu)
        if umut == True:
            spID += 1
            node.sp = spID
            for leaf in node:
                leaf.sp = spID
sys.stdout.write('SS') # SpreadSpeciation

# if requested, print speciational tree
if plot_trees:
    tmut = t.copy()
    for leaf in tmut:
        leaf.name = "["+str(leaf.sp)+"]"+leaf.name
    tmut.render(ophylo+"_2MUT.png", w=183, units="mm")


#======================================================#
# GET all the species IDs
SPdict, SPcounts = leaf_species_dictionary(t)
spIDs = sorted(SPcounts, key=lambda k: SPcounts[k], reverse=True)
sys.stdout.write('G') # Get species


#======================================================#
# CONVERT demography to phylogeny using a homemade traversing method

# also REFORMAT SFS listing in case of persistent paraphylies
# ie when a supposed parent is found as a child of another parent

#> the traversing could be sped up maybe
#> eg by computing the max distance between matching leaves
#> and pruning them out FIRST (ie pruning the largest paraphylies)

iter = 0
SFS = {}
print2iter = set(range_by_length(0, len(spIDs), 5))

for spid in spIDs:
    if spid in SPdict.keys():

	    matchleaves = SPdict[spid]

	    # populate the SFS
	    if len(matchleaves) == 1:
                SFS[matchleaves[0].name] = set()
                del SPdict[spid]
	    
            elif len(matchleaves) >= 2:

                firstleaf = matchleaves[0]
                MRCA = t.get_common_ancestor(matchleaves)
                sisters = MRCA.get_leaves()
                ssp = set()
                for ml in range(0, len(sisters)):
                    ssp.add(sisters[ml].name)

                ssp = ssp - set([firstleaf.name])

                # clear SP dict
                del SPdict[spid]
                for spd in SPdict.keys():
                    if SPdict[spd][0].name in ssp:
                       del SPdict[spd]
                
                # check if any [sisters] is already present as key in the SFS dictionary
                if len(SFS) > 0:
                    prevKeys = set(SFS.keys())
                    intersect = prevKeys & set(ssp)
                    if len(intersect) > 0:
                        #print "KATAM!"
                        if firstleaf.name in prevKeys:
                            SFS[firstleaf.name] = SFS[firstleaf.name].update(SFS[ssp])
                        else:
                            SFS[firstleaf.name] = ssp
                            for iSect in intersect:
                                SFS[firstleaf.name].update(SFS[iSect])
                                del SFS[iSect]
                    else:
                        SFS[firstleaf.name] = ssp
                else:
                    SFS[firstleaf.name] = ssp

	    # detach duplicate lineages at the root of their MRCA
	    if len(matchleaves) >= 2:
		
		if not MRCA.is_root():
		    upMRCA = MRCA.up
		    newdist = t.get_distance(upMRCA, firstleaf)
		    MRCA.detach()
		    upMRCA.add_child(firstleaf, firstleaf.name, newdist)        
		else:
		    newdist = t.get_distance(MRCA, firstleaf)
		    children = MRCA.get_children()
		    for child in children:
		        child.detach()
		    MRCA.add_child(firstleaf, firstleaf.name, newdist)

    iter += 1
    #if iter in print2iter: sys.stdout.write('C') # Convert

sys.stdout.write('CCCC')
#======================================================#

# ultrametricity in the output tree is defined in reference to the furthest leaf
if force_ultrametric:
    tree_dist = t.get_farthest_leaf()[1]
    for l in t:
        dst = t.get_distance(l)
        if dst != tree_dist:
            l.dist += tree_dist - dst
    
#======================================================#
# COMPUTE the SFS

SFS2 = SFS.copy()
nIndsSFS = 0
for i in SFS:
    spSFS = list(SFS[i])
    spSFS.append(i)
    sfs = [0] * nPops
    for j in range(0, len(spSFS)):
        pop = int(spSFS[j].split("_", 1)[0])
        sfs[pop-1] += 1
        nIndsSFS += 1
    SFS[i] = sfs

#======================================================#
# PRINT
## print SFS2
f = open(osfs+"_allinds.txt", 'w+')
for i in SFS2:
    f.write(i + "\t")
    f.write("\t".join(SFS2[i]))
    f.write("\n")
f.close()
## print SFS
f = open(osfs, 'w+')
start = True
for i in SFS:
    if start: f.write("Tip_label\t" + "\t".join("Isl_"+str(x) for x in range(0, len(SFS[i]))) + "\n")
    f.write(i + "\t")
    f.write("\t".join(str(x) for x in SFS[i]))
    f.write("\n")
    start = False
f.close()
## print PHYLO 
t.write(format=5, outfile=ophylo, dist_formatter='%0.20f')
if plot_trees:
    for leaf in t:
        leaf.name = leaf.sp
    t.render(ophylo+"_3PHYLO.png", w=183, units="mm")
sys.stdout.write('P] ') # Print
if not quiet: print "            >> " + str(time.time() - init) + " sec"
if not quiet: print "Total species count      >> " + str(len(SFS))
if (len(SFS) == len(t.get_leaves())) and (nIndsSFS == nIndsORI):
	print "Exit status              >> OK!"
else: sys.exit("ERROR. A bug must be crawling around... Please, please, contact the programmer.")




