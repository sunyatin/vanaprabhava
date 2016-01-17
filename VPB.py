
from ete3 import Tree
from ete3 import PhyloTree
import numpy
import time

# version with a matrix-intermediate
# guess it will be longer than v1.1

# assumes coalescent, ie ultrametric trees

numpy.random.seed(111)

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


t = PhyloTree("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/test.txt", format=1)

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

# cristallize the species info
def set_sp_ID(node):
    return node.sp
t.set_species_naming_function(set_sp_ID)

#t.show()

# print the mutational tree
for leaf in t:
    leaf.name = "_"+str(leaf.sp)+"_" +leaf.name
print "here"
t.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/2MUT.png", w=183, units="mm")

print "CONVERTING"

t.get_descendant_evol_events()

# collapser tous les evenements

#t2 = t.collapse_lineage_specific_expansions(species=None)
#print t2.get_ascii(show_internal=True)

#f = open('C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/SFS.txt', 'w+')

for node in t.traverse("preorder"): # better use inner nodes
    if not node.is_leaf():
        if node.evoltype == 'D':
            leaves = node.get_leaves()
            
            # populate the SFS
            #for leaf in leaves:
                #f.write(leaf.name+"\t")
            #f.write("\n")            
            
            upMRCA = node.up
            leaf = leaves[0]
            newDist = t.get_distance(node.up, leaf)
            leaf.detach()
            upMRCA.add_child(leaf, leaf.name, newDist)
            node.detach()

#f.close()

if plot_trees: t.render("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/3FINAL.png", w=183, units="mm")

print ">>>> " + str(time.time() - init)

for leaf in t:
    leaf.name = leaf.species
t.write(format=5, outfile="C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/V3.txt")

t = PhyloTree("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/V3.txt", format=5)
t2 = PhyloTree("C:/Users/Windows/Desktop/TRAVAIL FM/PythonicWork/Vanaprabhava/V1.1.txt", format=5)
 
comp = t2.compare(t)

print comp






