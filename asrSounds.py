import numpy as np
import pandas as pd
from ete3 import Tree
from alignment import tCoffee
from subprocess import Popen
import os


def nexCharOutput(chMtx, names, outfile,
                  datatype='STANDARD',
                  symbols=[],
                  gap='?',
                  missing='-'):
    f = open(outfile, 'w')
    f.write('#NEXUS\n\n')
    f.write('BEGIN DATA;\n')
    f.write('DIMENSIONS ntax=' + str(len(chMtx))+' NCHAR=' +
            str(len(chMtx.T))+';\n')
    ln = 'FORMAT DATATYPE='+datatype+' GAP='+gap+' MISSING='+missing
    if len(symbols) > 0:
        ln += ' symbols=\"' + ''.join(symbols) + '\"'
    ln += ';\n'
    f.write(ln)
    f.write('MATRIX\n\n')
    txLgth = max(map(len, names))
    for i in range(len(chMtx)):
        f.write(names[i].ljust(txLgth+2))
        for ch in chMtx[i]:
            if ch == -1:
                ch = '-'
            else:
                ch = str(ch)
            f.write(ch)
        f.write('\n')
    f.write('\n;\n\nEND;\n')
    f.close()


guideTree = Tree('romance.mcc.nwk')

asrCC = pd.read_csv('asrCC.csv', index_col=0, squeeze=True)


romance = np.array(guideTree.get_leaf_names())

data = pd.read_csv('albanoRomanceCC.csv', index_col=0)

data = data[(data.language.isin(romance)) & (data.cc.isin(asrCC.values))]

gp1 = -2.49302792222
gp2 = -1.70573165621

pmi = pd.read_csv('pmi-albanoRomance.csv', index_col=0)
sounds = np.array(pmi.index)

pmiDict = {(s1, s2): pmi[s1][s2]
           for s1 in sounds for s2 in sounds}


concepts = asrCC.index

aBlocks = pd.DataFrame()
for c in concepts:
    cData = data[data.cc == asrCC[c]]
    if len(cData) > 1:
        alg = list(tCoffee(cData.language.values, cData.word.values,
                           guideTree, pmiDict, gp1, gp2, sounds))
        cAlg = pd.DataFrame([list(x[1]) for x in alg],
                            index=[x[0] for x in alg])
        cAlg[cAlg == '-'] = '0'
    else:
        cAlg = pd.DataFrame(map(list, cData.word.values),
                            index=cData.language.values)
    cAlg = cAlg.reindex(romance, fill_value='-')
    cAlg.columns = [c+':'+str(i) for i in cAlg.columns]
    aBlocks = pd.concat([aBlocks, cAlg], axis=1)


aBlocks.to_csv('romanceAlignments.csv')

nexCharOutput(aBlocks.values, aBlocks.index,
              'romanceAlignments.nex',
              symbols=list(sounds)+['0'],
              gap='?',
              missing='-')



