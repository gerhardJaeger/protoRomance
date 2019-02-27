import numpy as np
from ete3 import Tree


bi = 2500

trees = []
with open('albanoRomance.run1.tre') as f:
    for i, ln in enumerate(f.readlines()):
        if i > bi:
            tr = Tree(ln.strip())
            trees.append(tr)
with open('albanoRomance.run2.tre') as f:
    for i, ln in enumerate(f.readlines()):
        if i > bi:
            t = Tree(ln.strip())
            trees.append(t)
with open('albanoRomance.run3.tre') as f:
    for i, ln in enumerate(f.readlines()):
        if i > bi:
            t = Tree(ln.strip())
            trees.append(t)
with open('albanoRomance.run4.tre') as f:
    for i, ln in enumerate(f.readlines()):
        if i > bi:
            t = Tree(ln.strip())
            trees.append(t)


taxa = np.array(trees[0].get_leaf_names())

romance = np.array([x for x in taxa if 'ALBANIAN' not in x])

for t in trees:
    t.set_outgroup('ALBANIAN')
    t.set_outgroup(t.get_common_ancestor([t & l for l in romance]))

with open('albanoRomance.posterior.tree', 'w') as f:
    f.write('')
with open('albanoRomance.posterior.tree', 'a') as f:
    for t in trees:
        f.write(t.write(format=1)+'\n')


for t in trees:
    t.prune([t & l for l in romance])

with open('romance.posterior.tree', 'w') as f:
    f.write('')
with open('romance.posterior.tree', 'a') as f:
    for t in trees:
        f.write(t.write(format=1)+'\n')
