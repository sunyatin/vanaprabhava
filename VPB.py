#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
VANAPRABHAVA v.2

Created on Wed Dec 23 17:33:31 2015
@author: Remi Tournebize
@contact: remi (dot) tournebize (at) gmail (dot) com
sotd: "Sabre a finances, corne de ma gidouille, madame la financiere, j'ai des oneilles pour parler et vous une bouche pour m'entendre."

v1.2: more profiling than v1.1 resulting in a significant gain of speed
v1.3_ms: few debuggings on the SFS (duplicates)
v1.4: debugged issues related to time decimal precision in Newick strings + modified force_ultrametricity module
v1.5 31122015: changed binomial to poisson random (large branch lengths throw an error on C long type for binomial)
v1.6 07012015: prints out the number of species 1 line AFTER the newick-formatted phylo tree in the .newick output files
v2.0 08012015: MAJOR CHANGE in the traversing function to convert demo into phylo: huge gain of speed for megademographies
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

"""
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
"""


t="C:/Users/Windows/Desktop/TRAVAIL FM/VANAPRABHAVA/VPB_DELL/temp/treeORI.txt"
mu= 0.000000067256260017746 #0.00000067256260017746
scale=2e+10
ophylo="C:/Users/Windows/Desktop/TRAVAIL FM/VANAPRABHAVA/VPB_DELL/phylo.txt"
osfs="C:/Users/Windows/Desktop/TRAVAIL FM/VANAPRABHAVA/VPB_DELL/sfs.txt"
force_ultrametric=True
plot_trees=False
quiet=False
ms_islands=[20000,20000,20000] #[20000,20000,20000]
ms_input=True
numpy.random.seed(10)

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

""" #v1.6
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
"""

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
# SET at once: "names" of inner nodes && initialize "species" attributes
if ms_input:
    ms_islands = numpy.array(ms_islands)
    ms_islands = numpy.cumsum(ms_islands)
    nPops = len(ms_islands)

innode = 0
largest_id = 0
nIndsORI = 0
for node in t.traverse():
    node.add_features(sp=1)
    node.img_style["size"] = 0
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
           node.add_features(mergedInd=node.name)
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
# CONVERT demography to phylogeny using a traversing method

# __/!\__ to be modified for non dichotomic trees (eg Lambda coalescent)
traversedNodes = set()

for node in t.traverse("preorder"):

	if node not in traversedNodes:

		if not node.is_leaf():

			children = node.get_children()
			if len(children) != 2:
				sys.exit("The algorithm does not know how to deal with non dichotomic trees for the moment!")

			# get child.1's tip species labels
			csp1 = set()
			for i in children[0]:
				csp1.add(i.sp)
			# get child.2's tip species labels
			csp2 = set()
			for i in children[1]:
				csp2.add(i.sp)

			# compare species label between children
			common = csp1.intersection(csp2)

			# if paraphyletic
			if len(common) > 0:

				if not node.is_root():

					# populate "mergedInd" feature for future SFS
					mergedLeaves = ""
					for l in node.iter_leaves():
						mergedLeaves = mergedLeaves+" "+l.name

					# collapse the subtree
					upNode = node.up
					newLeaf = node.get_farthest_leaf()
					newDist = newLeaf[1] + node.dist

					for childnode in node.traverse():
						traversedNodes.add(childnode)
					bad = node.detach()
					upNode.add_child(newLeaf[0], newLeaf[0].name, newDist)

					# the node.prune() function is too lenghty

					# actualize "mergedInd" feature of new leaf
					newLeaf[0].mergedInd = mergedLeaves

				else:

					# populate "mergedInd" feature for future SFS
					mergedLeaves = ""
					for l in node.iter_leaves():
						mergedLeaves = mergedLeaves+" "+l.name

					# collapse the subtree
					newLeaf = t.get_farthest_leaf()
					for child in t.get_children():
						child.detach()
					node.add_child(newLeaf[0], newLeaf[0].name, newLeaf[1])

					# actualize "mergedInd" feature of new leaf
					newLeaf[0].mergedInd = mergedLeaves

sys.stdout.write('CCCC')
nTrueSpecies = len(t)

#======================================================#
# COMPUTE & PRINT the SFS (full names & numerical)

f = open(osfs+"_allinds.txt", 'w+')
fn = open(osfs, 'w+')
fn.write("Tip_label\t" + "\t".join("Isl_"+str(x) for x in range(0, nPops)) + "\n")

for leaf in t.iter_leaves():
	inds = leaf.mergedInd.lstrip(" ")
	# write the SFS with full names
	f.write(leaf.name+" : "+inds+"\n")
	# compute the per-deme numerical SFS
	sfs = [0] * nPops
	inds = inds.split(" ")
	for i in inds:
		pop = int(i.split("_")[0])
		sfs[pop-1] += 1
		nIndsSFS += 1
	fn.write(leaf.name+"\t"+"\t".join(str(x) for x in sfs)+"\n")

f.close()
fn.close()


#======================================================#

# ultrametricity in the output tree is defined in reference to the furthest leaf
if force_ultrametric:
	tree_dist = t.get_farthest_leaf()[1]
	for l in t:
		dst = t.get_distance(l)
		if dst != tree_dist:
			l.dist += tree_dist - dst

#======================================================#

## print PHYLO
t.write(format=5, outfile=ophylo, dist_formatter='%0.20f')
# write the count of resultant species
f = open(ophylo, 'a')
f.write("\n"+str(nTrueSpecies)+"\n")
f.close()

if plot_trees:
	for leaf in t:
		leaf.name = leaf.sp
	t.render(ophylo+"_3PHYLO.png", w=183, units="mm")
"""
sys.stdout.write('P] ') # Print
if not quiet: print "            >> " + str(time.time() - init) + " sec"
if not quiet: print "Total species count      >> " + str(len(SFS))
if (len(SFS) == nTrueSpecies) and (nIndsSFS == nIndsORI):
	print "Exit status              >> OK!"
else: sys.exit("ERROR. A bug must be crawling around... Please, please, contact the programmer.")

"""

