##############################################################################
#
# Copyright (C) 2021 Minzhang Cheng
# Contact: minzhangcheng@gmail.com
#
# GNU Lesser proberal Public License Usage
# This file may be used under the terms of the GNU Lesser proberal Public
# License version 3 as published by the Free Software Foundation and
# appearing in the file LICENSE included in the packaging of this file.
# Please review the following information to ensure the GNU Lesser
# proberal Public License version 3 requirements will be met:
# http://www.gnu.org/licenses/gpl-3.0.html
#
##############################################################################


import random
import os
import pandas as pd


def randomSamples(clinicFile, count, key='', outputFile='', seed=0,
                  stratification='equal'):
    results = list()
    clinic = pd.read_csv(clinicFile, index_col=0)
    if not key:
        key = clinic.columns[-1]
    clinic = pd.DataFrame(clinic, columns=[key])
    clinic.dropna(axis=0, how='any')

    if stratification == 'random':
        random.seed(seed)
        sampleSelected = random.sample(list(clinic.index), count)
    else:
        groups = set(clinic[key])
        groupCount = {
            i: len(clinic[clinic[key] == i])
            for i in groups
        }
        allCount = len(clinic)
        countSelected = dict()
        if stratification == 'ratio':
            countSelected = {
                i: round(groupCount[i] / allCount * count)
                for i in groupCount
            }
        elif stratification == 'equal':
            countLeft = count
            groupLeft = len(groups)
            countSorted = sorted(groupCount.items(), key=lambda x: x[1])
            for i in countSorted:
                avg = countLeft / groupLeft
                k, value = i
                if value < avg:
                    countSelected[k] = value
                    countLeft = countLeft - value
                else:
                    countSelected[k] = round(avg)
                    countLeft = countLeft - round(avg)
                groupLeft = groupLeft - 1
        sampleSelected = list()
        for i in groups:
            random.seed(seed)
            population = list(clinic[clinic[key] == i].index)
            selected = random.sample(population, countSelected[i])
            sampleSelected.extend(selected)

    clinicSelected = pd.DataFrame(clinic, index=sampleSelected)
    if outputFile:
        clinicSelected.to_csv(outputFile)

    return clinicSelected


def exprSelected(exprFile, probeSelected, sampleSelected, outputFile=''):
    expr = pd.read_csv(exprFile, index_col=0)
    exprSelected = pd.DataFrame(expr,
                                index=probeSelected,
                                columns=sampleSelected)
    if outputFile:
        exprSelected.to_csv(outputFile)
    return exprSelected


def randomHeatmap(clinicFile, exprFile, probeSelected, runningDir, count,
                  symbolSelected=[], key='', seed=0, stratification='equal',
                  runningR=False):
    clinicSelectedFile = '{}/clinic_random_{}.csv'.format(runningDir, seed)
    exprSelectedFile = '{}/expr_random_{}.csv'.format(runningDir, seed)

    clinicSelected = randomSamples(clinicFile, count, key,
                                   clinicSelectedFile, seed, stratification)
    sampleSelected = list(clinicSelected.index)
    expr = exprSelected(exprFile, probeSelected,
                        sampleSelected, exprSelectedFile)

    rScript = 'library(gplots)\nlibrary(pheatmap)\nlibrary(RColorBrewer)\n\n'
    rScript += '\nexpr <- read.csv("{}", row.names=1)\n'.format(exprSelectedFile)
    if symbolSelected:
        rScript += 'symbols = c("{}")\n'.format('", "'.join(symbolSelected))
        rScript += 'row.names(expr) <- symbols\n'
    rScript += 'expr <- log(expr + 1, base=2)\n'
    rScript += 'x <- as.matrix(expr)\n'
    rScript += 'x <- t(scale(t(x)))\n'

    rScript += 'png(file="{}/pheatmap_random_{}.png", width=1024, height=1024, bg="transparent")\n' \
        .format(runningDir, seed)
    rScript += 'pheatmap(x, cutree_rows=2, cutree_cols=2, color=greenred(75), border_color=NA)\n'
    rScript += 'dev.off()\n'
    rScript += 'png(file="{}/heatmap2_random_{}.png", width=1024, height=1024, bg="transparent")\n' \
        .format(runningDir, seed)
    rScript += 'heatmap.2(x, col=greenred, scale="row", trace="none")\n'
    rScript += 'dev.off()\n'

    rScript += '\nclinic <- read.csv("{}", head=T, row.names=1)\n'.format(clinicSelectedFile)
    rScript += 'c <- clinic\n'
    rScript += 'annotation_c <- data.frame(c)\n'
    rScript += 'rownames(annotation_c) <- colnames(x)\n'
    for i in ['correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski']:
        rScript += 'png(file="{}/pheatmap2_random_{}_{}.png", width=1024, height=1024, bg="transparent")\n' \
            .format(runningDir, seed, i)
        rScript += 'pheatmap(as.matrix(x), annotation_col=annotation_c, color=bluered(200), border_color=NA, cutree_rows=2, cutree_cols=2, clustering_distance_cols="{}", scale="column")\n' \
            .format(i)
        rScript += 'dev.off()\n'

    rFile = '{}/heatmap_{}.R'.format(runningDir, seed)
    with open(rFile, 'w') as wf:
        print(rScript, file=wf)

    if runningR:
        if runningR:
            cmd = 'Rscript {} > {}.log'.format(rFile, rFile)
            r = os.system(cmd)


def randomHeatmapMultiple(clinicFile, exprFile, probeSelected, runningDir,
                          count, symbolSelected=[], repeat=10, key='', seed=[],
                          stratification='equal', runningR=False):
    if not seed:
        seed = list(range(0, repeat))
    for s in seed:
        randomHeatmap(clinicFile, exprFile, probeSelected, runningDir,
                      count, symbolSelected, key, s, stratification, runningR)


probes_11v21 = ['ENSG00000067646', 'ENSG00000114374', 'ENSG00000183878', 'ENSG00000198692', 'ENSG00000129824', 'ENSG00000067048', 'ENSG00000012817', 'ENSG00000131002', 'ENSG00000176728', 'ENSG00000099725', 'ENSG00000234854', 'ENSG00000165246', 'ENSG00000260197', 'ENSG00000262117', 'ENSG00000273906', 'ENSG00000206159', 'ENSG00000092377', 'ENSG00000241859', 'ENSG00000109132', 'ENSG00000105141', 'ENSG00000215580', 'ENSG00000227234', 'ENSG00000002745', 'ENSG00000169297', 'ENSG00000260458', 'ENSG00000179914', 'ENSG00000229807', 'ENSG00000154478', 'ENSG00000166049', 'ENSG00000139219', 'ENSG00000170370', 'ENSG00000175766', 'ENSG00000183166', 'ENSG00000242366', 'ENSG00000238269', 'ENSG00000157851', 'ENSG00000229847', 'ENSG00000136110', 'ENSG00000177108', 'ENSG00000248300', 'ENSG00000142515', 'ENSG00000125255', 'ENSG00000230805', 'ENSG00000270816', 'ENSG00000161055', 'ENSG00000092054', 'ENSG00000253593', 'ENSG00000152092', 'ENSG00000169594', 'ENSG00000204248']
symbol_11v21 = ['ZFY', 'USP9Y', 'UTY', 'EIF1AY', 'RPS4Y1', 'DDX3Y', 'KDM5D', 'TXLNGY', 'TTTY14', 'PRKY', 'LINC00676', 'NLGN4Y', 'AC010889.1', 'BCAR4', 'AC011297.1', 'GYG2P1', 'TBL1Y', 'ANOS2P', 'PHOX2B', 'CASP14', 'BCORP1', 'SPANXB1', 'WNT16', 'NR0B1', 'KCNJ18', 'ITLN1', 'XIST', 'GPR26', 'PASD1', 'COL2A1', 'EMX2', 'EIF4E1B', 'CALN1', 'UGT1A8', 'PAGE2B', 'DPYSL5', 'EMX2OS', 'CNMD', 'ZDHHC22', 'LINC02360', 'KLK3', 'SLC10A2', 'AL132709.1', 'LINC00221', 'SCGB3A1', 'MYH7', 'AC110741.1', 'ASTN1', 'BNC1', 'COL11A2']
clinicFile = 'C:/Users/minzhang/Desktop/11v21/clinic.csv'
exprFile = 'C:/Users/minzhang/Desktop/11v21/expr.csv'
runningDir = 'C:/Users/minzhang/Desktop/11v21'
randomHeatmapMultiple(clinicFile, exprFile, probes_11v21, runningDir, 50,
                      symbol_11v21, repeat=10, stratification='equal',
                      runningR=True)

probes_12v22 = ['ENSG00000114374', 'ENSG00000131002', 'ENSG00000012817', 'ENSG00000183878', 'ENSG00000273906', 'ENSG00000176728', 'ENSG00000067646', 'ENSG00000198692', 'ENSG00000067048', 'ENSG00000260197', 'ENSG00000129824', 'ENSG00000099725', 'ENSG00000206159', 'ENSG00000164816', 'ENSG00000165246', 'ENSG00000092377', 'ENSG00000231535', 'ENSG00000215580', 'ENSG00000241859', 'ENSG00000131482', 'ENSG00000267793', 'ENSG00000239893', 'ENSG00000081051', 'ENSG00000233070', 'ENSG00000168757', 'ENSG00000204941', 'ENSG00000122585', 'ENSG00000233052', 'ENSG00000254647', 'ENSG00000229807', 'ENSG00000244468', 'ENSG00000159516', 'ENSG00000105198', 'ENSG00000124467', 'ENSG00000242221', 'ENSG00000231924', 'ENSG00000203786', 'ENSG00000170162', 'ENSG00000203857', 'ENSG00000175535', 'ENSG00000243137', 'ENSG00000236731', 'ENSG00000268864', 'ENSG00000145242', 'ENSG00000146352', 'ENSG00000178287', 'ENSG00000236858', 'ENSG00000221826', 'ENSG00000257193', 'ENSG00000264706']
symbol_12v22 = ['USP9Y', 'TXLNGY', 'KDM5D', 'UTY', 'AC011297.1', 'TTTY14', 'ZFY', 'EIF1AY', 'DDX3Y', 'AC010889.1', 'RPS4Y1', 'PRKY', 'GYG2P1', 'DEFA5', 'NLGN4Y', 'TBL1Y', 'LINC00278', 'BCORP1', 'ANOS2P', 'G6PC', 'AC009977.1', 'ZNF736P9Y', 'AFP', 'ZFY-AS1', 'TSPY2', 'PSG5', 'NPY', 'AL513304.1', 'INS', 'XIST', 'AC093001.1', 'SPRR2G', 'LGALS13', 'PSG8', 'PSG2', 'PSG1', 'KPRP', 'VGLL2', 'HSD3B1', 'PNLIP', 'PSG4', 'AL135929.2', 'AC011487.2', 'EPHA5', 'CLVS2', 'SPAG11A', 'AL008638.2', 'PSG3', 'LINC02385', 'RN7SL217P']
clinicFile = 'C:/Users/minzhang/Desktop/12v22/clinic.csv'
exprFile = 'C:/Users/minzhang/Desktop/12v22/expr.csv'
runningDir = 'C:/Users/minzhang/Desktop/12v22'
randomHeatmapMultiple(clinicFile, exprFile, probes_12v22, runningDir, 50,
                      symbol_12v22, repeat=10, stratification='equal',
                      runningR=True)


