import pandas as pd


def nexCharOutput(chMtx, names, outfile, datatype='STANDARD'):
    f = open(outfile, 'w')
    f.write('#NEXUS\n\n')
    f.write('BEGIN DATA;\n')
    f.write('DIMENSIONS ntax=' + str(len(chMtx)) +
            ' NCHAR=' + str(len(chMtx.T))+';\n')
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


sc = pd.read_csv('albanoRomanceSC.nex',
                 skiprows=7,
                 skipfooter=4, engine='python',
                 index_col=0,
                 sep='\s+', header=None)

cc = pd.read_csv('albanoRomanceCC.nex',
                 skiprows=7,
                 skipfooter=4, engine='python',
                 index_col=0,
                 sep='\s+', header=None)

cc = cc.loc[sc.index]

scCC = pd.DataFrame([list(sc[1][l]+cc[1][l]) for l in sc.index],
                    index=sc.index)

nexCharOutput(scCC.values, scCC.index, 'albanoRomance_sc_cc.nex',
              datatype='restriction')

n = len(sc.values[0][0])
m = len(cc.values[0][0])

romance = [l for l in cc.index if 'ALBANIAN' not in l]
albanian = [l for l in cc.index if 'ALBANIAN' in l]

mbCommands = """#NEXUS
begin MrBayes;
      execute albanoRomance_sc_cc.nex;
      charset sc = 1-"""+str(n)+""";
      charset cc = """+str(n+1)+"""-"""+str(n+m)+""";
      partition dtype = 2:sc, cc;
      set partition = dtype;
      unlink Statefreq=(all) shape=(all);
      lset applyto=(all) rates=gamma;
      lset applyto=(1) coding=all;
      lset applyto=(2) coding=noabsencesites;
      prset brlenspr = clock:uniform;
      prset clockvarpr = igr;
      constraint romance = """+' '.join(romance)+""";
      constraint albanian = """+' '.join(albanian)+""";
      prset topologypr = constraints(romance, albanian);
      report applyto=(2) ancstates=yes;
      mcmcp stoprule=no stopval = 0.01 filename = albanoRomance nruns=4;
      mcmcp mcmcdiagn=yes diagnfreq=10000 samplefreq=10000 burninfrac=.5;
      set seed=12345;
      set swapseed=12345;
      mcmc ngen = 50000000;
      sumt;
      sump;
end;"""

with open('albanoRomance.mb.nex', 'w') as f:
    f.write(mbCommands)
