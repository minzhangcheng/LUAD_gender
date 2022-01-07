library(TCGAbiolinks)
library(SummarizedExperiment)

clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
save(clinical, file="/home/minzhang/TCGA/clinical.RData")
write.table(clinical,file = "/home/minzhang/TCGA/clinical.txt", sep='\t')

query <- GDCquery(project = "TCGA-BRCA",
				  data.category = "Transcriptome Profiling",
				  data.type = "Gene Expression Quantification",
				  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk = 20, directory="/home/minzhang/TCGA")
expr <- GDCprepare(query = query, 
                   save = TRUE, 
                   directory =  "/home/minzhang/TCGA",
                   save.filename = "/home/minzhang/TCGA/BRCA_rnaseq_counts.RData")
write.table(expr,file = "/home/minzhang/TCGA/BRCA_rnaseq_counts.txt", sep='\t')
save.image(file = "/home/minzhang/TCGA/BRCA_rnaseq_counts_all.RData")

