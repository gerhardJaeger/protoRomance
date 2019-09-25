import pandas as pd
import re
import lingpy as lp
from lingpy.evaluate.acd import bcubes


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


d1 = pd.read_csv('IELex+ASJP.csv', sep='\t')
d1 = d1[~d1.ASJP_transcription.isnull()]
d1 = d1[d1.loan == 0]
d1 = d1[d1['class'].notnull()]
d1['word'] = [cleanASJP(x) for x in d1.ASJP_transcription.values]
d1['language'] = [x.replace('-', '_') for x in d1.ASJP_language.values]

d2 = pd.read_csv('../albanoRomanceCC.csv', index_col=0)

d = pd.merge(d1, d2)[['concept', 'language', 'word', 'class', 'cc']]

d.to_csv('mergedData.tsv', sep='\t')

wl = lp.Wordlist('mergedData.tsv')

bcubes(wl, gold="class", test="cc")

# *************************
# * B-Cubed-Scores        *
# * --------------------- *
# * Precision:     0.9934 *
# * Recall:        0.6194 *
# * F-Scores:      0.7630 *
# *************************'
