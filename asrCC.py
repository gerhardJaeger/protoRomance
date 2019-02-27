import numpy as np
import pandas as pd


cc = pd.read_csv('armenoRomanceCCbin.csv',
                 index_col=0,
                 dtype='str')

romance = np.array([x for x in cc.index if 'ARMENIAN' not in x])

cc = cc.loc[romance]

with open('romanceCC.tsv', 'w') as f:
    for i in cc.index:
        f.write(i+'\t')
        f.write('\t'.join(cc.loc[i].values)+'\n')

p1 = pd.read_csv('armenoRomance.run1.p',
                 sep='\t', skiprows=1)

p2 = pd.read_csv('armenoRomance.run2.p',
                 sep='\t', skiprows=1)

p3 = pd.read_csv('armenoRomance.run3.p',
                 sep='\t', skiprows=1)

p4 = pd.read_csv('armenoRomance.run4.p',
                 sep='\t', skiprows=1)


bi = p1.Gen.max() // 2
pr = pd.concat([p1[p1.Gen > bi],
                p2[p2.Gen > bi],
                p3[p3.Gen > bi],
                p4[p4.Gen > bi]])

asrCharacters = [x for x in pr.columns if 'p(1)' in x]

asr = pr[asrCharacters].mean()

ccConcepts = pd.DataFrame([x.split(':') for x in cc.columns],
                          columns=['concept', 'class'])

ccConcepts['posterior'] = asr.values

concepts = np.sort(ccConcepts.concept.unique())

asr = []
for c in concepts:
    cData = ccConcepts[ccConcepts.concept == c].copy()
    cData.index = range(len(cData))
    i = cData.posterior.argmax()
    asr.append(c+':' + cData['class'][i])

asr = pd.DataFrame(asr, index=concepts)

asr.to_csv('asrCC.csv', header=None)
