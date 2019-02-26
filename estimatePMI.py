import random as pyrandom
import numpy as np
import pandas as pd
import Levenshtein
import re
from Bio import pairwise2
from multiprocessing import Process, Manager

ncores = 6
np.random.seed(12345)
pyrandom.seed(12345)

# Function: nexCharOutput
# Description:
## This function takes a character array, a list of rownames
## and the name of the output nexus file as input
## and writes the character matrix into a nexus file.
## Missing entries are assumed to be coded as "-1"


def nexCharOutput(chMtx, names, outfile, datatype='STANDARD'):
    f = open(outfile, 'w')
    f.write('#NEXUS\n\n')
    f.write('BEGIN DATA;\n')
    f.write('DIMENSIONS ntax=' +
            str(len(chMtx))+' NCHAR='+str(len(chMtx.T))+';\n')
    f.write('FORMAT DATATYPE='+datatype+' GAP=? MISSING=- interleave=yes;\n')
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


data = pd.read_csv('dataset.tab',
                   index_col=0, na_filter=False, sep='\t')

data = data[data.wls_gen.isin(['ROMANCE', 'ALBANIAN'])]
data = data[data.index != 'LATIN']

concepts100 = np.array(data.columns[9:])

nEntries = data[concepts100].apply(lambda x:
                                   sum(x != '')).sort_values()

concepts = nEntries.index[-40:]

data = data[concepts]

taxa = np.array(data.index)


def cleanASJP(word):
    """takes an ASJP string as argument
    and returns the string with all diacritics removed."""
    word = re.sub(r", ", "-", word)
    word = re.sub(r"\%", "", word)
    word = re.sub(r"\*", "", word)
    word = re.sub(r"\"", "", word)
    word = re.sub(r".~", "", word)
    word = re.sub(r"(.)(.)(.)\$", r"\2", word)
    word = re.sub(r"\$", "", word)
    word = re.sub(r"\s+", "", word)
    return word.replace('~', '')


dataWL = pd.DataFrame([(c, l, cleanASJP(w))
                       for c in concepts
                       for l in taxa
                       for w in data[c][l].split(',')
                       if data[c][l] != ''],
                      columns=['concept', 'language', 'word'])

dataWL = dataWL.drop_duplicates(['concept', 'language'])

training = pd.DataFrame()
for c in concepts:
    cData = dataWL[dataWL.concept == c]
    cTraining = pd.DataFrame([cData.loc[[i, j]].word.values
                              for i in cData.index
                              for j in cData.index
                              if i < j])
    training = training.append(cTraining)


sounds = np.unique(np.concatenate(list(map(list, dataWL.word.values))))


def levalign(w):
    """takes a pair of strings as input
    and returns the Levenshtein alignment
    in column format. Gaps are removed."""
    x, y = w
    algn = np.zeros((0, 2))
    if '0' in [x, y]:
        return algn
    e = Levenshtein.opcodes(x, y)
    for a in e:
        if a[0] in ['replace', 'equal']:
            x_a, x_e = a[1], a[2]
            y_a, y_e = a[3], a[4]
            ag = [list(x[x_a:x_e]),  list(y[y_a:y_e])]
            algn = np.concatenate([algn, np.transpose(ag)])
    return algn


manager = Manager()
return_dict = manager.dict()

packages = np.array_split(training.values, ncores)


def doWork(i, pck):
    return_dict[i] = np.vstack([levalign(p) for p in pck])


jobs = []
for i, pck in enumerate(packages):
    p = Process(target=doWork, args=(i, pck))
    p.start()
    jobs.append(p)

for p in jobs:
    p.join()

alg0 = np.vstack(return_dict.values())


# count alignment frequencies
sFreqs = pd.crosstab(alg0[:, 0], alg0[:, 1])
sFreqs = sFreqs.reindex(sounds, fill_value=0).T
sFreqs = sFreqs.reindex(sounds, fill_value=0).T
# symmetrize
sFreqs = sFreqs.copy()+sFreqs.copy().T

# add-1 smoothing
sFreqs += 1

sProbs = sFreqs/sFreqs.sum().sum()

# extract relative sound frequencies
soundOccurrences = np.concatenate([list(w) for w in dataWL.word.values])

soundProbabilities = pd.value_counts(soundOccurrences, normalize=True)[sounds]


pmi0 = (np.log(sProbs).copy() -
        np.log(soundProbabilities)).T-np.log(soundProbabilities)


def sscore(a, b, pmiDict, gp1, gp2):
    """a, b: ASJP strings
    pmiDict: logodds dictionary
    gp1, gp2: gap penalties
    return PMI score of a/b
    """
    out = pairwise2.align.globalds(a, b, pmiDict, gp1, gp2)
    if len(out) == 0:
        return np.nan
    return out[0][2]


def scoreNW(x, y, pmiDict, gp1, gp2):
    """x, y: sequences of ASJP strings,  separated by '-'
    pmiDict: logodds dictionary
    gp1, g2: gap penalties
    returns maximal PMI score for the Cartesian product of x and y"""
    if '0' in [x, y]:
        return np.nan
    x1 = x.split('-')
    y1 = y.split('-')
    return max([sscore(xx, yy, pmiDict, gp1, gp2) for xx in x1 for yy in y1])


pmi0Dict = {(s1, s2): pmi0[s1][s2]
            for s1 in sounds for s2 in sounds}


def nw(x, y, pmiDict, gp1, gp2):
    """wrapper for Bio.pairwise2.align.globalds"""
    return pairwise2.align.globalds(x, y, pmiDict, gp1, gp2)


def nwalign(w, pmiDict, gp1, gp2, th=-np.Inf):
    """w: pair of ASJP strings
    pmiDict: dictionary of logodds (=PMI scores)
    gp1, gp2: gap penalties (non-positive)
    th: threshold; all pairs with a PMI-score <th will be ignored
    returns: np.array of pairwise alignment,  with gaps removed"""
    x, y = w
    a = nw(x, y, pmiDict, gp1, gp2)
    if len(a) == 0:
        return np.zeros((0, 2))
    algn = []
    if a[0][2] < th:
        return np.zeros((0, 2))
    for aa in a:
        ll = len(aa[0])
        aaa = [[aa[0][i], aa[1][i]] for i in range(ll)]
        algn += [x for x in aaa if '-' not in x]
    return np.array(algn)


def nwalignStar(crp, pmiDict, gp1, gp2, th):
    packages = np.array_split(crp, ncores)
    manager = Manager()
    return_dict = manager.dict()

    def doWork(i, pck):
        return_dict[i] = np.vstack([nwalign(w, pmiDict, gp1, gp2, th)
                                    for w in pck])
    jobs = []
    for i, pck in enumerate(packages):
        p = Process(target=doWork, args=(i, pck))
        p.start()
        jobs.append(p)
    for p in jobs:
        p.join()
    return np.vstack(return_dict.values())


gp1 = -2.49302792222
gp2 = -1.70573165621
th = 4.4451


pmiDict = pmi0Dict.copy()
for i in range(10):
    alg = nwalignStar(training.values, pmiDict, gp1, gp2, th)
    sFreqs = pd.crosstab(alg[:, 0], alg[:, 1])
    sFreqs = sFreqs.reindex(sounds, fill_value=0).T
    sFreqs = sFreqs.reindex(sounds, fill_value=0).T
    sFreqs = sFreqs.copy()+sFreqs.copy().T
    sFreqs += 1
    sProbs = sFreqs/sFreqs.sum().sum()
    pmi = (np.log(sProbs).copy() -
           np.log(soundProbabilities)).T - np.log(soundProbabilities)
    pmiDict = {(s1, s2): pmi[s1][s2]
               for s1 in sounds for s2 in sounds}


pmi.to_csv('pmi-albanoRomance.csv')


dataWL.to_csv('albanoRomanceASJP.csv', index=False)


sc = pd.DataFrame(index=taxa)
for c in concepts:
    cData = dataWL[dataWL.concept == c]
    cTaxa = cData.language.unique()
    cWords = pd.Series([''.join(cData[cData.language == l].word.values)
                        for l in cTaxa],
                       index=cTaxa)
    cMtx = pd.DataFrame([[int(s in cWords[l]) for s in sounds]
                         for l in cTaxa],
                        index=cTaxa)
    cMtx = cMtx.reindex(taxa, fill_value='-')
    sc = pd.concat([sc, cMtx], axis=1)

nexCharOutput(sc.values, sc.index, 'albanoRomanceSC.nex', 'restriction')


#########


def levalignFull(w):
    """takes a pair of strings as input
    and returns the Levenshtein alignment
    in column format."""
    x, y = w
    algn = np.zeros((0, 2))
    e = Levenshtein.opcodes(x, y)
    for a in e:
        if a[0] in ['replace', 'equal']:
            x_a, x_e = a[1], a[2]
            y_a, y_e = a[3], a[4]
            ag = [list(x[x_a:x_e]), list(y[y_a:y_e])]
        elif a[0] == 'delete':
            x_a, x_e = a[1], a[2]
            ag = [list(x[x_a:x_e]), ['-']*(x_e-x_a)]
        else:
            y_a, y_e = a[3], a[4]
            ag = [['-'] * (y_e-y_a), list(y[y_a:y_e])]
        algn = np.concatenate([algn, np.transpose(ag)])
    return algn
