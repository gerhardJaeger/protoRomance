import numpy as np
from ete3 import Tree
import re
from Bio import pairwise2


def match(x, y):
    return [np.argwhere(y == z)[0][0] if z in y else None for z in x]


def nwBio(x, y, lodict, gp1, gp2):
    al = pairwise2.align.globalds(x, y, lodict, gp1, gp2)[0]
    return al[2], np.array(al[:2])


def nw(x, y, lodict, gp1, gp2):
    """
    Needleman-Wunsch algorithm for pairwise string alignment
    with affine gap penalties.
    'lodict' must be a dictionary with all symbol pairs as keys
    and match scores as values.
    gp1 and gp2 are gap penalties for opening/extending a gap.
    Returns the alignment score and one optimal alignment.
    """
    n, m = len(x), len(y)
    dp = np.zeros((n+1, m+1))
    pointers = np.zeros((n+1, m+1), int)
    for i in range(1, n+1):
        dp[i, 0] = dp[i-1, 0] + (gp2 if i > 1 else gp1)
        pointers[i, 0] = 1
    for j in range(1, m+1):
        dp[0, j] = dp[0, j-1] + (gp2 if j > 1 else gp1)
        pointers[0, j] = 2
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = dp[i-1, j-1]+lodict[x[i-1], y[j-1]]
            insert = dp[i-1, j] + (gp2 if pointers[i-1, j] == 1 else gp1)
            delet = dp[i, j-1]+(gp2 if pointers[i, j-1] == 2 else gp1)
            dp[i, j] = max([match, insert, delet])
            pointers[i, j] = np.argmax([match, insert, delet])
    alg = []
    i, j = n, m
    while(i > 0 or j > 0):
        pt = pointers[i, j]
        if pt == 0:
            i -= 1
            j -= 1
            alg = [[x[i], y[j]]]+alg
        if pt == 1:
            i -= 1
            alg = [[x[i], '-']]+alg
        if pt == 2:
            j -= 1
            alg = [['-', y[j]]]+alg
    return dp[-1, -1], np.array([''.join(x) for x in np.array(alg).T])


def algnMtx(al, sounds):
    """
    Takes a pairwise alignment (i.e. a pair of gapped strings
    with identical length)
    as input and returns a matrix representation M as output.
    The matrix M is defined as M[i, j] = 1 if x[i] is matched with y[j]
    in the alignment, 0 else (where x, y are the two ungapped
    strings to be aligned).
    """
    w1 = ''.join(np.array([s for s in al[0] if s != '-']))
    w2 = ''.join(np.array([s for s in al[1] if s != '-']))
    dm = np.zeros((len(w1), len(w2)), int)
    i, j = 0, 0
    for s1, s2 in np.array([list(w) for w in al]).T:
        if s1 in sounds:
            if s2 in sounds:
                dm[i, j] += 1
                i += 1
                j += 1
            else:
                i += 1
        else:
            j += 1
    return dm


def createLibrary(words, lodict, gp1, gp2, sounds):
    """
    Takes a list of sequences and returns a library in the sense of the
    T-Coffee algorithm. A library is a dictionary with sequence pairs
    as keys and pairwise alignments in matrix format as columns.
    """
    library = dict()
    for w1 in words:
        for w2 in words:
            if (w2, w1) in library:
                x = library[w2, w1]
                library[w1, w2] = x[0].T, x[1]
            else:
                a1, a2 = nw(w1, w2, lodict, gp1, gp2)[1]
                library[w1, w2] = algnMtx([a1, a2],
                                          sounds), (1-sHamming(a1, a2))
    return library


def sHamming(x, y):
    """
    Takes two gapped strings and returns the hamming distance between
    them. Positions containing a gap in at least one string are ignored.
    """
    w1, w2 = np.array(list(x)), np.array(list(y))
    r = np.mean(w1[(w1 != '-') * (w2 != '-')] != w2[(w1 != '-') * (w2 != '-')])
    if np.isnan(r):
        r = 1.
    return r


def wHamming(w1, w2, m):
    """
    Takes two words and an alignment matrix and returns the hamming distance
    between them. Positions containing a gap in at least one string
    are ignored.
    """
    id1, id2 = np.argwhere(m == 1).T
    return np.mean(np.array(list(w1))[id1] != np.array(list(w2))[id2])


def createExtendedLibrary(words, lodict, gp1, gp2, sounds):
    """
    Takes a list of sequences and returns an extended library in the
    sense of the T-Coffee algorithm. An extended library is a dictionary with
    sequence pairs as keys and a score matrix as values.
    For a pair of sequences x, y and a corresponding score matrix M,
    M[i, j] is the score for aligning x[i] with y[j].
    """
    library = createLibrary(words, lodict, gp1, gp2, sounds)
    extLibrary = dict()
    for w1 in words:
        for w2 in words:
            dm = np.zeros((len(w1), len(w2)))
            for w3 in words:
                a1, s1 = library[w1, w3]
                a2, s2 = library[w3, w2]
                dm += (s1+s2) * np.dot(a1, a2)
            extLibrary[w1, w2] = dm
    return extLibrary


# pointers[i, j] = argmax([match, insert, delet])

def nwBlock(b1, b2, lib):
    """
    Needleman-Wunsch alignment of two aligned blocks b1 and b2,
    using the scores in the extended library lib.
    """
    def pos(gappedString, i):
        """
        Returns the index of gappedString[i] in the
        ungapped version thereof.
        If gappedString[i] is a gap, returns -1
        """
        if gappedString[i] != '-':
            return (0 if i == 0
                    else sum(np.array(list(gappedString[:i])) != '-'))
        else:
            return -1
    words1 = np.array([re.sub("-", "", w) for w in b1])
    words2 = np.array([re.sub("-", "", w) for w in b2])
    n, m = len(b1[0]), len(b2[0])
    dp = np.zeros((n+1, m+1))
    pointers = np.zeros((n+1, m+1), int)
    pointers[0, 1:] = 2
    pointers[1:, 0] = 1
    for i in range(1, n+1):
        for j in range(1, m+1):
            insert = dp[i-1, j]
            delet = dp[i, j-1]
            match = dp[i-1, j-1] + sum([0 if '-' in [gs1[i-1], gs2[j-1]]
                                        else lib[w1, w2][pos(gs1, i-1),
                                                         pos(gs2, j-1)]
                                        for (w1, gs1) in zip(words1, b1)
                                        for (w2, gs2) in zip(words2, b2)])
            dp[i, j] = max([match, insert, delet])
            pointers[i, j] = np.argmax([match, insert, delet])
    al1 = np.transpose(np.array([list(w) for w in b1]))
    al2 = np.transpose(np.array([list(w) for w in b2]))
    alCombined = []
    while max(i, j) > 0:
        p = pointers[i, j]
        if p == 0:
            alCombined = [list(al1[i-1]) + list(al2[j-1])] + alCombined
            i -= 1
            j -= 1
        elif p == 1:
            alCombined = [list(al1[i-1])+['-']*len(b2)] + alCombined
            i -= 1
        else: 
            alCombined = [['-']*len(b1)+list(al2[j-1])] + alCombined
            j -= 1
    return np.array([''.join(x) for x in np.array(alCombined).T]), dp[-1, -1]


def tCoffee(wTaxa, words, tree, lodict, gp1, gp2, sounds):
    """tree must be a binary branching tree"""
    lib = createExtendedLibrary(words, lodict, gp1, gp2, sounds)
    wTree = Tree(tree.write(format=9))
    wTree.prune([wTree & l for l in wTaxa])
    wDict = dict(zip(wTaxa, words))
    for nd in wTree.traverse('postorder'):
        if nd.is_leaf():
            nd.add_feature('algn', np.array([wDict[nd.name]]))
            nd.add_feature('nTaxa', np.array([nd.name]))
        else:
            dl, dr = nd.get_children()
            b1, b2 = dl.algn, dr.algn
            a = nwBlock(b1, b2, lib)[0]
            nd.add_feature('algn', a)
            nd.add_feature('nTaxa', np.concatenate([dl.nTaxa, dr.nTaxa]))
    return zip(wTree.nTaxa, wTree.algn)

