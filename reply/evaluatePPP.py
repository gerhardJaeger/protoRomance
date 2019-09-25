import numpy as np
import pandas as pd
import os

np.random.seed(12345)


def nexCharOutput(chMtx, names, outfile, datatype='STANDARD'):
    f = open(outfile, 'w')
    f.write('#NEXUS\n\n')
    f.write('BEGIN DATA;\n')
    f.write('DIMENSIONS ntax=' + str(len(chMtx)) + ' NCHAR=' +
            str(len(chMtx.T))+';\n')
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


trace = pd.read_csv('output/ctmc_posterior_run_1.var',
                    sep="\t")

empData = '../albanoRomanceCC.nex'


id0 = np.array(trace.index[1:])


simIndices = id0[::10]


with open(empData, 'r') as f:
    nex = f.readlines()

df1 = pd.DataFrame({x.split()[0]:
                    list(x.split()[1])
                    for x in nex[7:-4]}).T


simDirs = ['output/ctmc_post_sims/posterior_predictive_sim_'+str(i+1)
           for i in range(len(simIndices))]


def scoreF(df):
    return (df1.T.apply(lambda x:
                        sum(x == '1')) / (df1.T.apply(lambda x:
                                                      sum(x == '0')) +
                                          df1.T.apply(lambda x:
                                                      sum(x == '1')))).std()


scores = [scoreF(df1)]

spl = np.random.permutation(1000)

for i in spl:
    simData = os.getcwd() + '/' + simDirs[i]+'/phyloSeq[1].nex'
    with open(simData, 'r') as f:
        nex = f.readlines()
    df1 = pd.DataFrame({x.split()[0]:
                       list(x.split()[1])
                       for x in nex[6:-3]}).T
    simScore = scoreF(df1)
    scores.append(simScore)

with open('ctmc_std.csv', 'w') as f:
    for x in scores:
        f.write(str(x)+'\n')
