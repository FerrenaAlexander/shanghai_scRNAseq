library(tidyverse)

samps = as.data.frame(readxl::read_excel('data/samplemetadata.xlsx'))

files = as.data.frame(readxl::read_excel('data/samplemetadata.xlsx', sheet = 2))


#size GB
agg <- aggregate(SizeGB ~ Sample, files, sum)
agg <- agg[match(str_sort(agg$Sample, numeric = T), agg$Sample),]

samps$SizeGB <- agg$SizeGB


#num read pears
agg <- aggregate(NumReadPairs ~ Sample, files, sum)
agg <- agg[match(str_sort(agg$Sample, numeric = T), agg$Sample),]

samps$NumReadPairs <- agg$NumReadPairs


write.csv(samps[,11:12], '~/Desktop/sampsagg.csv', row.names = F, quote = F)
