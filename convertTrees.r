library(ape)

for (i in 1:4) {
    tr <- read.nexus(paste0('armenoRomance.run',i,'.t'))
    write.tree(tr,paste0('armenoRomance.run',i,'.tre'))
}
