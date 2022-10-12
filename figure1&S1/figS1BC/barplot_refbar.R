
library(stringdist)
library(ggpubr)
b0='CTGCACCGACGCATT'
b1='GAAGGCACAGACTTT'
b2='GTTGGCTGCAATCCA'
b3='AAGACGCCGTCGCAA'
b4='TATTCGGAGGACGAC'
barcodes=c(b0,b1,b2,b3,b4)
tag.dists=data.frame(matrix(ncol = 2,nrow = 10))
colnames(tag.dists)=c('hd','group')
tag.dists$group=c('b0b1','b0b2','b0b3','b0b4',
                  'b1b2','b1b3','b1b4',
                  'b2b3','b2b4',
                  'b3b4')

tag.dists[1,1]=stringdistmatrix(a=barcodes[1], b=barcodes[2],method = 'hamming')
tag.dists[2,1]=stringdistmatrix(a=barcodes[1], b=barcodes[3],method = 'hamming')
tag.dists[3,1]=stringdistmatrix(a=barcodes[1], b=barcodes[4],method = 'hamming')
tag.dists[4,1]=stringdistmatrix(a=barcodes[1], b=barcodes[5],method = 'hamming')
tag.dists[5,1]=stringdistmatrix(a=barcodes[2], b=barcodes[3],method = 'hamming')
tag.dists[6,1]=stringdistmatrix(a=barcodes[2], b=barcodes[4],method = 'hamming')
tag.dists[7,1]=stringdistmatrix(a=barcodes[2], b=barcodes[5],method = 'hamming')
tag.dists[8,1]=stringdistmatrix(a=barcodes[3], b=barcodes[4],method = 'hamming')
tag.dists[9,1]=stringdistmatrix(a=barcodes[3], b=barcodes[5],method = 'hamming')
tag.dists[10,1]=stringdistmatrix(a=barcodes[4], b=barcodes[5],method = 'hamming')

# j=1
# for (i in 1:5) {
#   for (n in 1:5)  {
#     ifelse(i!=n, tag.dists[j] <- stringdistmatrix(a=barcodes[i], b=barcodes[n],method = 'hamming'), print('skip'))
#     j=j+1
#   }
# }

ggbarplot(tag.dists, x='group', y='hd')

temp=as.data.frame(table(tag.dists$hd))
ggbarplot(temp, x='Var1', y='Freq',label = TRUE, label.pos = "out",
          xlab='HD', ylab='Frequency')
ggsave('hd across five ref barcodes.pdf',width = 4,height = 3)

