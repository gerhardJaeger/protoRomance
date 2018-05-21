# protoRomance

This directory contains all data and scripts used and produced for Section 3 of the manuscript *Computational Historical Linguistics*.



## Table of contents

### original data:

dataset.tab (ASJP v. 17, accessed from asjp.clld.org on 8/2/16)

### Script files

#### Python

- alignment.py  
- estimatePMI.py  
- cognateClustering.py  
- combinePartitions.py  
- mbPostprocessing.py
- asrCC.py  
- asrSounds.py  
- evaluation.py  

####R (for conversion between Newick and Nexus format for phylogenies)

- convert_mccTree.r  
- convertTrees.r

#### MrBayes

- albanoRomance.mb.nex

#### BayesTraits

- asrCC.bt
- btAlignments.txt

#### Paup*

- convertRomancePosterior.paup

## Intermediate results

- pmi-albanoRomance.csv (PMI scores)
  albanoRomanceASJP.csv (word lists used)
  albanoRomanceSC.nex (character matrix)
  albanoRomance.wordpairs.csv (word pairs with PMI similarities, used for cognate
       clustering)
  albanoRomanceCC.nex (character matrix)
  albanoRomanceCCbin.csv (character matrix)
  albanoRomanceCC.csv (word lists with inferred class labels)
  albanoRomance_sc_cc.nex (combined character matrix)
- albanoRomance.run1.p (MrBayes output)
  albanoRomance.run1.t (MrBayes output)
  albanoRomance.run2.p (MrBayes output)
  albanoRomance.run2.t (MrBayes output)
  albanoRomance.ckp (MrBayes output)
  albanoRomance.con.tre (MrBayes output)
  albanoRomance.parts (MrBayes output)
  albanoRomance.mcmc (MrBayes output)
  albanoRomance.trprobs (MrBayes output)
  albanoRomance.tstat (MrBayes output)
  albanoRomance.vstat (MrBayes output)

- 
  albanoRomance.run1.tre (MrBayes output, converted to Newick)
  albanoRomance.run2.tre (MrBayes output, converted to Newick)

- albanoRomance.posterior.tree (posterior sample)
  romance.posterior.tree (posterior sample, rooted and pruned to Romance)
  romance.posterior.nex.tree (posterior sample, rooted and pruned to Romance, Nexus format)

- albanoRomance.mcc.tre (maximum clade credibility tree, Nexus format)
  albanoRomance.mcc.nwk (maximum clade credibility tree, Newick format)

- romanceCC.tsv (character matrix for BayesTraits)

- romanceCC.tsv.log.txt (BayesTraits log file)

- asrCC.csv (ASR for cognate classes)

- romanceAlignments.nex (character matrix)
  romanceAlignments.tsv (character matrix for BayesTraits)
  romanceAlignments.tsv.log.txt (BayesTraits log file)

## Final results

- reconstruction.csv (reconstructed word list)
- romanceEvaluation.csv (evaluation results)

## Documentation

- workflow.txt
- README.md