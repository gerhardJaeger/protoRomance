import numpy as np
import pandas as pd
import Levenshtein
import re


def cleanASJP(word):
    """takes an ASJP string as argument
    and returns the string with all diacritics removed."""
    word = re.sub(r",", "-", word)
    word = re.sub(r"\%", "", word)
    word = re.sub(r"\*", "", word)
    word = re.sub(r"\"", "", word)
    word = re.sub(r".~", "", word)
    word = re.sub(r"(.)(.)(.)\$", r"\2", word)
    word = re.sub(r"\$", "", word)
    word = re.sub(r"\s+", "", word)
    return word.replace('~', '')


asjp = pd.read_table('../dataset.tab', index_col=0,
                     sep='\t', na_filter=False)

romance = np.array([x for x in asjp[asjp.wls_gen == 'ROMANCE'].index
                    if x != 'LATIN'])


reconstruction = pd.read_csv('../reconstruction.csv', index_col=0)

concepts = np.array(reconstruction.index)


def ldn(a, b):
    return min([1.*Levenshtein.distance(x, y)/max(len(x), len(y))
                for x in a.split('-') for y in b.split('-')])


romanceCleaned = pd.DataFrame([[cleanASJP(x).split('-')[0] for x in y]
                               for y in asjp.loc[romance][concepts].values],
                              index=romance,
                              columns=concepts)

latinCleaned = pd.Series([cleanASJP(x) for x in reconstruction.Latin.values],
                         index=reconstruction.index)[concepts]


reconEval = np.mean([ldn(x, y)
                     for x, y in zip(reconstruction.reconstruction.values,
                                     latinCleaned.values)])


romanceEval = pd.Series([np.mean([ldn(romanceCleaned.loc[l][c],
                                      latinCleaned[c])
                                  for c in concepts
                                  if romanceCleaned.loc[l][c] != ''])
                         for l in romance],
                        index=romance)

romanceEval.loc['Proto-Romance'] = reconEval

romanceEval.to_csv('romanceEvaluationCorrected.csv', header=False)
