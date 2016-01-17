# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 17:33:31 2015
@author: Remi Tournebize
@contact: remi (dot) tournebize (at) gmail (dot) com
@sotd: Sabre a finances, corne de ma gidouille, madame la financiere, j'ai des oneilles pour parler et vous une bouche pour m'entendre.
"""

from ete3 import Tree
import numpy
import time
import sys
import argparse

# more profiling than v1.1
# perspectives: add a protracted speciation mode
# should need some more testing if **multiple** coalescent trees (beta...)

# to speed up runtime we made an extensive use of the set
# data structure. It implies that THERE MUST NOT BE ANY DUPLICATE
# TIP LABELS IN THE INPUT TREE

# assumes coalescent, ie ultrametric trees
# but in case you doubt abt it, set "force_ultrametric" option to "True"


################
## PARAMETERS ##
################

parser = argparse.ArgumentParser(description='Example with simples options')

parser.add_argument('-m', '--mutationRate', type=float, 
                    required=True, action="store", help="Mean mutation rate per lineage per generation")
parser.add_argument('-s', '--scalingFactor', type=float, 
                    required=True, action="store", help="A factor to rescale branch lengths")
parser.add_argument('-ophylo', '--PHYLOoutput', type=str, required=True, 
                    action="store", help="Path of the output phylogeny")
parser.add_argument('-osfs', '--SFSoutput', type=str, required=True, 
                    action="store", help="Path of the output SFS")

parser.add_argument('-I', '--islands', type=int, nargs="+",
                    action="store", help="Number of individuals in each deme separated by a comma, must be in the same order as in your MS command")


parser.add_argument('--forceUltrametric',
                    action="store_true", help="Whether to force the input tree to be ultrametric (useful if you rescale branch lengths)")
parser.add_argument('--drawTrees',
                    action="store_true", help="Whether to plot trees")
parser.add_argument('--seed', type=int, default=-1,
                    action="store", help="Random seed")

args = parser.parse_args()

print "===================================================="
mu = args.m #ok
print "Mutation rate: "+str(mu)
s = args.s #ok
print "Scaling factor: "+str(s)
ophylo = args.ophylo #ok
osfs = args.osfs #ok
if len(args.I) >= 1:
    ms_islands = [args.I] #ok
    ms_input = True
    print "Number of individuals in islands: "+" ".join(ms_islands)
else:
    ms_input = False
force_ultrametric = args.forceUltrametric #ok
if force_ultrametric: print "Will force ultrametricity"
plot_trees = args.drawTrees #ok
if plot_trees: print "Will draw trees"
seed = args.seed #ok
if seed != -1:
    print "Will set random seed at "+str(seed)
    numpy.random.seed(seed)

print "Will output phylo to: "+ophylo
print "Will output SFS to: "+osfs
print "===================================================="

###############
## FUNCTIONS ##
###############

def ubranch_mutation(node, mu):
    rb = numpy.random.binomial(node.dist, mu)
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
        sys.exit("1. The MS island structure you provided through -ms does not fit the number of tips in the demography. The option might have been mispecified. Check it out!")
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
t = Tree("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/test.txt", format=1)
if force_ultrametric:
    for node in t.traverse():
        node.img_style["size"] = 0
    tree_dist = t.get_farthest_leaf()[1]    
    t.convert_to_ultrametric(tree_length=tree_dist, strategy='balanced')
sys.stdout.write('R') # Read


#======================================================#
# SET at once: "names" of inner nodes && "species" attributes
if ms_input:
    ms_islands = numpy.array(ms_islands)
    ms_islands = numpy.cumsum(ms_islands)
    nPops = len(ms_islands)

innode = 0
largest_id = 0
for node in t.traverse():
    node.add_features(sp=1)
    if not node.is_leaf():
        node.name = "iN%d" %innode
        innode += 1
    elif ms_input == True:
        inn = int(node.name)
        if inn > largest_id: largest_id = inn
        pop = get_deme(inn, ms_islands)
        node.name = str(pop)+"_"+node.name

if ms_input == True & ms_islands[-1] > largest_id:
    sys.exit("2. The MS island structure you provided through -ms does not fit the number of tips in the demography. The option might have been mispecified. Check it out!")

if plot_trees: t.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/1ORIGINAL.png", w=183, units="mm")
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
        leaf.name = str(leaf.sp) + "/" + leaf.name
    t.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/2MUT_v11.png", w=183, units="mm")


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
    matchleaves = SPdict[spid]
    # populate the SFS
    if len(matchleaves) >= 1:
        ssp = set()
        for ml in range(1, len(matchleaves)):
            ssp.add(matchleaves[ml].name)
        firstleaf = matchleaves[0]
        # check if any [matchleaves] is already present as key in the SFS dictionary
        if len(SFS) > 0:
            prevKeys = set(SFS.keys())
            intersect = prevKeys & set(matchleaves)
            if len(intersect) > 0:
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
        MRCA = t.get_common_ancestor(matchleaves)
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

    del SPdict[spid]
    iter += 1
    if iter in print2iter: sys.stdout.write('C') # Convert


#======================================================#

if force_ultrametric: 
    tree_dist = t.get_farthest_leaf()[1]    
    t.convert_to_ultrametric(tree_length=tree_dist, strategy='balanced')
    
#======================================================#
# COMPUTE the SFS

SFS2 = SFS.copy()
for i in SFS:
    spSFS = list(SFS[i])
    sfs = [0] * nPops
    for j in range(0, len(spSFS)):
        pop = int(spSFS[j].split("_", 1)[0])
        sfs[pop-1] += 1
    SFS[i] = sfs

#======================================================#
# PRINT
## print SFS2
f = open('C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/SFS2.txt', 'w+')
for i in SFS2:
    f.write(i + "\t")
    f.write("\t".join(SFS2[i]))
    f.write("\n")
f.close()
## print SFS
f = open('C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/SFS2.txt', 'w+')
for i in SFS:
    f.write(i + "\t")
    f.write("\t".join(str(x) for x in SFS[i]))
    f.write("\n")
f.close()
## print PHYLO

for leaf in t: # TMP
    leaf.name = leaf.sp  # TMP
    
t.write(format=5, outfile="C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/V1.2.txt")
if plot_trees:
    for leaf in t:
        leaf.name = leaf.sp
        t.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/3_FINAL_NEW.png", w=183, units="mm")
sys.stdout.write('P] ') # Print

print ">> " + str(time.time() - init) + " sec"
print "Total species count >> " + str(len(SFS))




