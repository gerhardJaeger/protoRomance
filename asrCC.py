import numpy as np
import pandas as pd
from subprocess import Popen
import os


cc = pd.read_csv('albanoRomanceCCbin.csv',
                 index_col=0,
                 dtype='str')

romance = np.array([x for x in cc.index if 'ALBANIAN' not in x])

cc = cc.loc[romance]

with open('romanceCC.tsv', 'w') as f:
    for i in cc.index:
        f.write(i+'\t')
        f.write('\t'.join(cc.ix[i].values)+'\n')

paupCommands = """#Nexus
Begin Paup;
set incr=auto;
gettrees file = romance.posterior.tree;
savetrees file = romance.posterior.nex.tree replace=yes brlen=user;
q;
End;
"""

with open('convertRomancePosterior.paup', 'w') as f:
    f.write(paupCommands)

p = Popen('paup4 convertRomancePosterior.paup>/dev/null', shell=True)
os.waitpid(p.pid, 0)

btCommands = """1
1
seed 12345;
mlt 1;
ga 4;
run"""

with open('asrCC.bt', 'w') as f:
    f.write(btCommands)

cmd = 'BayesTraits romance.posterior.nex.tree '
cmd += 'romanceCC.tsv < asrCC.bt>/dev/null'

p = Popen(cmd, shell=True)
os.waitpid(p.pid, 0)


results = pd.read_table('romanceCC.tsv.log.txt',
                        sep='\t',
                        skiprows=25)

cl = [x for x in results.columns if 'P(1)' in x]

results = results[cl]
results.columns = cc.columns

concepts = np.unique([x.split(':')[0] for x in results.columns])

winners = []
for c in concepts:
    cChars = [x for x in results.columns if x.split(':')[0] == c]
    winners.append(results.mean()[cChars].argmax())

winners = pd.Series(winners, index=concepts)

winners.to_csv('asrCC.csv')
