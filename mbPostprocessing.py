import pandas as pd
from numpy import *
import subprocess,os
from ete2 import Tree



p1 = pd.read_table('albanoRomance.run1.p',
                   sep='\t',skiprows=1)

p2 = pd.read_table('albanoRomance.run2.p',
                   sep='\t',skiprows=1)

pr = pd.concat([p1[p1.Gen>1000000],p2[p2.Gen>1000000]])

trees = []
with open('albanoRomance.run1.tre') as f:
    for i,ln in enumerate(f.readlines()):
        if i >100:
            t = Tree(ln.strip())
            trees.append(t)
with open('albanoRomance.run2.tre') as f:
    for i,ln in enumerate(f.readlines()):
        if i >100:
            t = Tree(ln.strip())
            trees.append(t)


taxa = array(trees[0].get_leaf_names())

romance = array([x for x in taxa if not 'ALBANIAN' in x])

for t in trees:
    t.set_outgroup('ALBANIAN')
    t.set_outgroup(t.get_common_ancestor([t&l for l in romance]))

with open('albanoRomance.posterior.tree','w') as f: f.write('')
with open('albanoRomance.posterior.tree','a') as f:
    for t in trees:
        f.write(t.write(format=1)+'\n')


for t in trees:
    t.prune([t&l for l in romance])

with open('romance.posterior.tree','w') as f: f.write('')
with open('romance.posterior.tree','a') as f:
    for t in trees:
        f.write(t.write(format=1)+'\n')
