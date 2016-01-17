
from ete3 import Tree
import numpy
import time

# assumes coalescent, ie ultrametric trees

numpy.random.seed(2)

init = time.time()

mu = 1e-4
plot_trees = False

# on peut se passer du search_leaves_by_sp
# si on traverse en preorder
# ou alors en postorder?
       
def ubranch_mutation(node, mu):
    rb = numpy.random.binomial(node.dist, mu)
    if rb >= 1:
        return True
    else:
        return False
        
def get_all_sp(t):
    spIDs = []
    counts = []
    for leaf in t:
        spid = leaf.sp
        if spid not in spIDs:
            spIDs.append(spid)
            counts.append(1)
        else:
            counts[spIDs.index(spid)] += 1
    spIDs = [x for y, x in sorted(zip(counts, spIDs), reverse=True)]
    return spIDs
    
def search_leaves_by_sp(t, value):
    leaves = []
    for leaf in t.iter_leaves():
        if leaf.sp == value:
            leaves.append(leaf)
    return(leaves)


t = Tree("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/test.txt", format=1)

#print t

print "Initializing traversing"
# set: "names" of inner nodes && "species" attribute at once
inode = 0
for node in t.traverse():
    node.add_features(sp=1)
    if not node.is_leaf():
        node.name = "N_%d" %inode
        inode += 1
    
if plot_trees: t.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/1ORIGINAL.png", w=183, units="mm")


print "Spreading speciation events"
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

# print the mutational tree
#t2 = t.copy()
#for leaf in t2:
#    leaf.name = leaf.sp
#if plot_trees: t2.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/2MUT.png", w=183, units="mm")

# print the mutational tree
for leaf in t:
    leaf.name = "_"+str(leaf.sp)+"_" +leaf.name
t.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/2MUT_v11.png", w=183, units="mm")


print "Get (decr) sorted species IDs"
spIDs = get_all_sp(t)
print len(spIDs)

print "Converting to phylogeny"
# the traversing could be sped up maybe
# eg by computing the max distance between matching leaves
# and pruning them out FIRST (ie pruning the largest paraphylies)

# DO IT BY: PREORDER / LEVELORDER / POSTORDER
# PREORDER : serait bcp trop long
# LEVELORDER : idem
# without addition of an external function
# en post order il y a moyen de perdre tres vite de la complexite
# en 5 levels a peu pres isomorphe, on gagne environ
# n/5 diversite

#TreeNode.up.remove_child('x')

iter = 0
SFS = []
outSP = []
for spid in spIDs:
    matchleaves = search_leaves_by_sp(t, spid)
    
    if len(matchleaves) >= 1:
        ssp = []
        for ml in matchleaves:
            ssp.append(ml.name)
        SFS.append(ssp)
        outSP.append(matchleaves[0].name)
    
    if len(matchleaves) >= 2:
        
        # note 4 future profiling: use a list comprehension
        
        MRCA = t.get_common_ancestor(matchleaves)
        if not MRCA.is_root():
            upMRCA = MRCA.up
            keepleaf = matchleaves[0]
            newdist = t.get_distance(upMRCA, keepleaf)
            MRCA.detach()
            upMRCA.add_child(keepleaf, keepleaf.name, newdist)        
        else:
            keepleaf = matchleaves[0]
            newdist = t.get_distance(MRCA, keepleaf)
            MRCA.children[0].detach()
            MRCA.children[0].detach()
            MRCA.add_child(keepleaf, keepleaf.name, newdist)
    iter += 1
    if iter % 100 == 0: print iter,

#print t
#print SFS


def get_duplicates(l):
  seen = set()
  seen_add = seen.add
  overseen = set( x for x in l if x in seen or seen_add(x) )
  return list( overseen )
  
  
# reformat the SFS list of lists
#flatSFS = [item for sublist in SFS for item in sublist]
#dupSP = set([x for x in flatSFS if flatSFS.count(x) > 1])
#if len(dupSP) > 0:
#    for dsp in dupSP:
dupRows = []
for i in range(len(SFS)):
    for j in range(1, len(SFS[i])):
        if SFS[i][j] in outSP:
            print str(i) + "  " + str(j)
            print "REDUNDANCY!!!"
            dupRows.append(i)            
for i in sorted(dupRows, reverse=True):
    del SFS[i]
    
f = open('C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/SFS.txt', 'w+')
for i in range(len(SFS)):
    f.write("\t".join(SFS[i]))
    f.write("\n")
f.close()
    
    
# TMP
for leaf in t:
    leaf.name = leaf.sp    
    
#if plot_trees: 
t.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/3_FINAL_NEW.png", w=183, units="mm")
t.write(format=5, outfile="C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/V1.1.txt")

print ">>>> " + str(time.time() - init)

