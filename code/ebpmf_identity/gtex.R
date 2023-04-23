############### gtex #################
sample.df <- read.delim("/project2/mstephens/dongyue/gtex/V8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", as.is=TRUE, header=TRUE, row.names=1)
subject.df<- read.delim('/project2/mstephens/dongyue/gtex/V8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', as.is=TRUE, header=TRUE, row.names=1)

tpm.df <- read.delim("/project2/mstephens/dongyue/gtex/V8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
                     as.is=T, row.names=1, check.names=FALSE, skip=2)
dim(res)

table(sample.df['SMAFRZE'])
rnaseq.sample.df <- sample.df[sample.df['SMAFRZE']=='RNASEQ', ]
as.matrix(sort(table(rnaseq.sample.df['SMTSD']), decreasing=TRUE))


gene_names = read.table('http://stephenslab.github.io/count-clustering/project/utilities/gene_names_all_gtex.tx')
