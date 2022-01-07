library(edgeR)
library(TCGAbiolinks)

setwd('/home/minzhang/LUAD_sex')

expr.all <- read.csv('expr.csv')
clinic <- read.csv('clinic.csv')


# 1	11-21	Male vs Female in Non-Smoking
nonsmoking <- clinic[clinic$smoking == 0, ]
ids = row.names(nonsmoking)
expr <- expr.all[, ids]
expr[expr<0] <- 0
groupList <- nonsmoking$gender
design = model.matrix(~ 0 + factor(groupList))
colnames(design) = c('Male', 'Female')
DGElist <- DGEList(counts = expr,group = groupList)
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
write.csv(nrDEG_edgeR, '11-21.csv')


# 2	12-22	Male vs Female in Smoking
smoking <- clinic[clinic$smoking == 1, ]
ids = row.names(smoking)
expr <- expr.all[, ids]
expr[expr<0] <- 0
groupList <- smoking$gender
design = model.matrix(~ 0 + factor(groupList))
colnames(design) = c('Male', 'Female')
DGElist <- DGEList(counts = expr,group = groupList)
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
write.csv(nrDEG_edgeR, '12-22.csv')


# 3	11-12	Non-Smoking vs Smoking in Male
male <- clinic[clinic$gender == 1, ]
ids = row.names(male)
expr <- expr.all[, ids]
expr[expr<0] <- 0
groupList <- male$smoking
design = model.matrix(~ 0 + factor(groupList))
colnames(design) = c('non-smoking', 'smoking')
DGElist <- DGEList(counts = expr,group = groupList)
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
write.csv(nrDEG_edgeR, '11-12.csv')


# 4	21-22	Non-Smoking vs Smoking in Female
female <- clinic[clinic$gender == 0, ]
ids = row.names(female)
expr <- expr.all[, ids]
expr[expr<0] <- 0
groupList <- female$smoking
design = model.matrix(~ 0 + factor(groupList))
colnames(design) = c('non-smoking', 'smoking')
DGElist <- DGEList(counts = expr,group = groupList)
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
write.csv(nrDEG_edgeR, '21-22.csv')


# 5	10-20	Male vs Female
ids = row.names(clinic)
expr <- expr.all[, ids]
expr[expr<0] <- 0
groupList <- clinic$gender
design = model.matrix(~ 0 + factor(groupList))
colnames(design) = c('Male', 'Female')
DGElist <- DGEList(counts = expr,group = groupList)
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
write.csv(nrDEG_edgeR, '10-20.csv')


# 6	01-02	Non-Smoking vs Smoking
ids = row.names(clinic)
expr <- expr.all[, ids]
expr[expr<0] <- 0
groupList <- clinic$smoking
design = model.matrix(~ 0 + factor(groupList))
colnames(design) = c('non-smoking', 'smoking')
DGElist <- DGEList(counts = expr,group = groupList)
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
write.csv(nrDEG_edgeR, '01-02.csv')

