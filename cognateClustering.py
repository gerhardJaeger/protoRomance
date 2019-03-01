import random as pyrandom
pyrandom.seed(12345)
from numpy import *
random.seed(12345)
import pandas as pd
from Bio import pairwise2
from sklearn.metrics import precision_recall_curve
from sklearn.linear_model import LogisticRegression
import igraph
from scipy import stats
import tempfile
from ete2 import Tree
import subprocess
import os

gp1=-2.49302792222
gp2=-1.70573165621


# Function: nexCharOutput
# Description:
## This function takes a character array, a list of rownames
## and the name of the output nexus file as input
## and writes the character matrix into a nexus file.
## Missing entries are assumed to be coded as "-1"
                   
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


def sscore(a,b,pmiDict,gp1,gp2):
    """a,b: ASJP strings
    pmiDict: logodds dictionary
    gp1,gp2: gap penalties
    return PMI score of a/b
    """
    out = pairwise2.align.globalds(a,b,pmiDict,gp1,gp2)
    if len(out)==0: return nan
    return out[0][2]

def scoreNW(x,y,pmiDict,gp1,gp2):
    """x,y: sequences of ASJP strings, separated by '-'
    pmiDict: logodds dictionary
    gp1,g2: gap penalties
    returns maximal PMI score for the Cartesian product of x and y"""
    if '0' in [x,y]: return nan
    x1=x.split('-')
    y1=y.split('-')
    return max([sscore(xx,yy,pmiDict,gp1,gp2) for xx in x1 for yy in y1])


data = pd.read_csv('albanoRomanceASJP.csv')
data['ID'] = range(len(data))

pmi = pd.read_csv('pmi-albanoRomance.csv',index_col=0)
sounds = array(pmi.index)
pmiDict = {(s1,s2):pmi[s1][s2]
           for s1 in sounds for s2 in sounds}


taxa = data.language.unique()


lpairs = pd.DataFrame([(l1,l2)
                       for i,l1 in enumerate(taxa)
                       for j,l2 in enumerate(taxa)
                       if i<j])


wpairs = pd.DataFrame()
for l1,l2 in lpairs.ix[random.permutation(lpairs.index)].values:
    l1Data = data[data.language==l1]
    l2Data = data[data.language==l2]
    lpPairs = pd.DataFrame([list(l1Data[['concept','language','word','ID']].ix[i])+
                            list(l2Data[['concept','language','word','ID']].ix[j])
                            for i in l1Data.index
                            for j in l2Data.index],
                           columns=['concept1','language1','word1','ID1',
                                    'concept2','language2','word2','ID2'])
    wpairs = pd.concat([wpairs,lpPairs])

wpairs['target'] = array(wpairs.concept1==wpairs.concept2,int)

wpairs['PMI'] = [sscore(a,b,pmiDict,gp1,gp2)
                 for (a,b) in wpairs[['word1','word2']].values]

wpairs.to_csv('albanoRomance.wordpairs.csv',index=False)

lr = LogisticRegression()
lr.fit(c_[wpairs.PMI.values],wpairs.target.values)




synpairs = wpairs[wpairs.target==1][['concept1',
                                     'language1', 'language2',
                                     'word1','word2',
                                     'ID1','ID2','PMI']]

synpairs.columns = ['concept']+list(synpairs.columns[1:])

concepts = data.concept.unique()

synpairs['prediction'] = lr.predict_proba(c_[synpairs.PMI.values])[:,1]




ccData = pd.DataFrame()
th = .5
for c in concepts:
    cData = data[data.concept==c].copy()
    cPairs = synpairs[synpairs.concept==c]
    cIDs = cData.ID.values
    simMtx = zeros((len(cIDs),len(cIDs)))
    simMtx[pd.match(cPairs.ID1.values,cIDs),
           pd.match(cPairs.ID2.values,cIDs)] = cPairs.prediction.values
    simMtx[pd.match(cPairs.ID2.values,cIDs),
           pd.match(cPairs.ID1.values,cIDs)] = cPairs.prediction.values
    simMtx[simMtx<th]=0
    G = igraph.Graph.Weighted_Adjacency(list(simMtx))
    clusters = G.community_label_propagation(weights='weight')
    ccDict = {cIDs[x]:i for i,cl in enumerate(clusters)
              for x in cl}
    cData['cc'] = [c+':'+str(ccDict[i]) for i in cData.ID.values]
    ccData = pd.concat([ccData,cData])


taxa = ccData.language.unique()

ccMtx = pd.DataFrame(index=taxa)
for c in concepts:
    cData = ccData[ccData.concept==c]
    cMtx = pd.crosstab(cData.language,cData.cc)
    cMtx[cMtx>1]=1
    cMtx = cMtx.reindex(taxa,fill_value='-')
    ccMtx = pd.concat([ccMtx,cMtx],axis=1)



ccMtx.to_csv('albanoRomanceCCbin.csv')

nexCharOutput(ccMtx.values,ccMtx.index,'albanoRomanceCC.nex')


ccData.to_csv('albanoRomanceCC.csv',index='False')
