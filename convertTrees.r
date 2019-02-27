library(ape)

for (i in 1:4) {
    tr <- read.nexus(paste0('albanoRomance.run',i,'.t'))
    write.tree(tr,paste0('albanoRomance.run',i,'.tre'))
}
