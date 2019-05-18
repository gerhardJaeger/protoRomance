library(ggplot2)

dat <- read.table('romanceEvaluation.csv',
                  sep=',')

colnames(dat) <- c('doculect','LDN')

pr <- dat[dat$doculect=='Proto-Romance',]

extant = dat[dat$doculect!='Proto-Romance',]

p <- ggplot(data=extant,aes(x=LDN)) +
    geom_histogram(colour='black',fill='white',binwidth=.01) +
    geom_vline(aes(xintercept=pr$LDN),colour='red',
               linetype='dashed',size=1) +
    xlab('Average normalized Levenshtein distance to Latin')
p

ggsave('evaluation.svg',p)
