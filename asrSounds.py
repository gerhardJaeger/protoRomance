from numpy import *
import pandas as pd
from ete2 import Tree
from alignment import tCoffee
from subprocess import Popen
import os

def nexCharOutput(chMtx,names,outfile,datatype='STANDARD'):
    f = open(outfile,'w')
    f.write('#NEXUS\n\n')
    f.write('BEGIN DATA;\n')
    f.write('DIMENSIONS ntax='+str(len(chMtx))+' NCHAR='+str(len(chMtx.T))+';\n')
    f.write('FORMAT DATATYPE='+datatype+' GAP=? MISSING=- interleave=yes;\n')
    f.write('MATRIX\n\n')
    txLgth = max(map(len,names))
    for i in xrange(len(chMtx)):
        f.write(names[i].ljust(txLgth+2))
        for ch in chMtx[i]:
            if ch==-1: ch='-'
            else:
                ch = str(ch)
            f.write(ch)
        f.write('\n')
    f.write('\n;\n\nEND;\n')
    f.close()


guideTree = Tree('albanoRomance.mcc.nwk').get_children()[0]

asrCC = pd.read_csv('asrCC.csv',index_col=0,header=None)[1]


romance = array(guideTree.get_leaf_names())

data = pd.read_csv('albanoRomanceCC.csv',index_col=0)

data = data[(data.language.isin(romance))&(data.cc.isin(asrCC.values))]

gp1=-2.49302792222
gp2=-1.70573165621

pmi = pd.read_csv('pmi-albanoRomance.csv',index_col=0)
sounds = array(pmi.index)

pmiDict = {(s1,s2):pmi[s1][s2]
           for s1 in sounds for s2 in sounds}


concepts = asrCC.index

aBlocks = pd.DataFrame()
for c in concepts:
    cData = data[data.cc==asrCC[c]]
    if len(cData)>1:
        alg = tCoffee(cData.language.values,cData.word.values,
                      guideTree,pmiDict,gp1,gp2,sounds)
        cAlg = pd.DataFrame([list(x[1]) for x in alg],
                            index = [x[0] for x in alg])
        cAlg[cAlg=='-'] = '0'
    else:
        cAlg = pd.DataFrame(map(list,cData.word.values),index=cData.language.values)
    cAlg = cAlg.reindex(romance,fill_value='-')
    cAlg.columns = [c+':'+str(i) for i in cAlg.columns]
    aBlocks = pd.concat([aBlocks,cAlg],axis=1)


binMtx = pd.DataFrame(index=romance)
for i in aBlocks.columns:
    cl = aBlocks[i].values
    states = unique([x for x in cl if x!='-'])
    clMtx = pd.DataFrame([cl==s for s in states],dtype=int,
                         columns=romance,
                         index=[i+':'+s for s in states]).T
    clMtx[cl=='-'] = '-'
    binMtx = pd.concat([binMtx,clMtx],axis=1)

nexCharOutput(binMtx.values,binMtx.index,'romanceAlignments.nex')

binMtx.to_csv('romanceAlignments.tsv',sep='\t',header=None)

with open('btAlignments.txt','w') as f:
    f.write('1\n1\n')
    f.write('seed 12345;\n')
    f.write('mlt 1;\n')
    f.write('run\n')



p = Popen('BayesTraits romance.posterior.nex.tree romanceAlignments.tsv < btAlignments.txt>/dev/null',
          shell=True)
os.waitpid(p.pid,0)

results = pd.read_csv('romanceAlignments.tsv.log.txt',
                      skiprows=23,sep='\t')


idc = [x for x in results.columns if 'P(1)' in x]

results = results[idc]

results.columns = binMtx.columns

resultsMean = results.mean()

def recon(res):
    asr = []
    for x in aBlocks.columns:
        idc = [y for y in res.index if ':'.join(y.split(':')[:2])==x]
        asr.append(res[idc].sort_values().index[-1].split(':')[-1])
    asr = pd.Series(asr,index=aBlocks.columns)
    reconstruction = []
    for c in concepts:
        idc = [x for x in asr.index if x.split(':')[0]==c]
        reconstruction.append(''.join(asr[idc].values).replace('0',''))
    reconstruction = pd.Series(reconstruction,index=concepts)
    return reconstruction


reconstruction = pd.DataFrame(recon(resultsMean))


asjp = pd.read_table('dataset.tab',index_col=0,
                     sep='\t')

reconstruction['Latin'] = asjp.ix['LATIN'][concepts].values

reconstruction.columns = ['reconstruction','Latin']
reconstruction['concept'] = reconstruction.index

reconstruction[['concept','Latin','reconstruction']].to_csv('reconstruction.csv',
                                                            index=False)
