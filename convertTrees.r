library(ape)

t1 <- read.nexus('albanoRomance.run1.t')
write.tree(t1,'albanoRomance.run1.tre')

t2 <- read.nexus('albanoRomance.run2.t')
write.tree(t2,'albanoRomance.run2.tre')
