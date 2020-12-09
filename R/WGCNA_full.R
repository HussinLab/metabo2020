#!/usr/bin/env Rscript
#adapted from https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

## Data input and cleaning

# Load the WGCNA package
library(WGCNA);
options(stringsAsFactors = FALSE)

# Read expression data
datExpr.ori = read.table("~/Dropbox/Projects/Metabolomics/MI_project/all_scripts_and_inputs_for_paper/WGCNA/data_normalized_all_BL_DC.csv", sep="\t", header=TRUE)
rownames(datExpr.ori) = datExpr.ori$metabolites
datExpr.ori$metabolites = NULL
datExpr = datExpr.ori

# Read metadata file
datDrug = read.table("~/Dropbox/Projects/Metabolomics/MI_project/all_scripts_and_inputs_for_paper/WGCNA/metadata_BL_DC.csv", sep="\t", header=TRUE)
row.names(datDrug) = datDrug$SAMPLE_NAME
datDrug$SAMPLE_NAME = NULL

# Match drug and exp data
datExpr = as.data.frame(t(datExpr))

# Match drug samples to filtered exp data
names(datDrug) = "CONDITION"


## 2.a. Choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# pickSoftThreshold: will use block size 280.
#  pickSoftThreshold: calculating connectivity for given powers...
#    ..working on genes 1 through 280 of 280
#    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1   0.3160 -0.656         0.4200  40.500  3.67e+01  85.90
# 2      2   0.8090 -1.140         0.7540  12.000  7.97e+00  44.90
# 3      3   0.8220 -1.170         0.7910   5.350  2.38e+00  29.30
# 4      4   0.8170 -1.180         0.8030   3.030  9.79e-01  21.50
# 5      5   0.7780 -1.190         0.7550   1.980  4.30e-01  17.00
# 6      6   0.7150 -1.140         0.7330   1.420  2.17e-01  14.10
# 7      7   0.7660 -1.130         0.8280   1.080  1.17e-01  12.10
# 8      8   0.6600 -1.090         0.7700   0.851  6.20e-02  10.60
# 9      9   0.7580 -1.080         0.8400   0.694  3.55e-02   9.37
# 10    10   0.6720 -1.030         0.8590   0.579  2.06e-02   8.41
# 11    12   0.6960 -1.010         0.8850   0.423  7.14e-03   6.92
# 12    13   0.6730 -0.978         0.8740   0.369  4.40e-03   6.33
# 13    14   0.0397 -0.784         0.0943   0.324  2.58e-03   5.81
# 14    15   0.0396 -0.764         0.0783   0.287  1.60e-03   5.36
# 15    16   0.7430 -0.947         0.9380   0.255  1.01e-03   4.95
# 16    17   0.7320 -0.920         0.9350   0.229  6.41e-04   4.58
# 17    18   0.7600 -0.941         0.9680   0.206  4.16e-04   4.25
# 18    19   0.7730 -0.919         0.9790   0.186  2.67e-04   3.94
# 19    20   0.7360 -0.901         0.9700   0.168  1.72e-04   3.67

# Plot the results:
#sizeGrWindow(9, 5)
png("WGCNA/2__soft-thresholding_power.png", width=1500, height=1500, res=200)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"), ylim=c(-0.6,1));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off();


## Construction of the gene network and identification of modules

# Choose a set of soft-thresholding powers
#first step
#power = seq(from = 2, to=12, by=1)
#mergeCutHeight = seq(from = 0.05, to=0.5, by=0.05)

#last step for correlation
power = 4
mergeCutHeight = 0.35


#In WGCNA FAQ:
#Number of samples	Unsigned and signed hybrid networks	Signed networks
#Less than 20	        9	                                18
#20-30	                8	                                16
#30-40	                7	                                14
#more than 40	        6	                                12

for (i in power){
    for (y in mergeCutHeight){ 
        
        # here we define the adjacency matrix using soft thresholding with beta=6
        ADJ1=abs(cor(datExpr,use="p"))^i
        # When you have relatively few genes (<5000) use the following code
        k=as.vector(apply(ADJ1,2,sum, na.rm=T))
        # When you have a lot of genes use the following code
        k=softConnectivity(datE=datExpr,power=i)
        # Plot a histogram of k and a scale free topology plot
        fig = paste("WGCNA/2__network_connectivities.power_", i, ".mergeCutHeight", y, ".png", sep="")
        png(fig, width=1500, height=1500, res=200)
        #pdf(fig)
        par(mfrow=c(1,2))
        hist(k)
        scaleFreePlot(k, main="Check scale free topology\n")
        dev.off();

        ## 2.a.2 One-step network construction and module detection

        # Network construction and module detection
        nethybrid = blockwiseModules(datExpr, power = i, #maxBlockSize = 20000,
                                        TOMType = "signed", 
                                        minModuleSize = 5,
                                        reassignThreshold = 0, mergeCutHeight = y,
                                        numericLabels = TRUE,
                                        saveTOMs = TRUE, networkType = "signed",
                                        #pamStage = TRUE, 
                                        pamRespectsDendro = FALSE,
                                        #pamRespectsDendro = TRUE,
                                        verbose = 5, 
                                        #corType = "bicor", maxPOutliers = 0.1,
                                        pearsonFallback = "individual", deepSplit=3)
                                    
        fig = paste("~/Dropbox/Projects/Metabolomics/MI_project/all_scripts_and_inputs_for_paper/WGCNA/2__dendrograms.geneblocks.auto.nethybrid.power_", i, ".mergeCutHeight_", y, ".png", sep="")
        png(fig, width=1500, height=1500, res=200)
        # Convert labels to colors for plotting
        mergedColors = labels2colors(nethybrid$colors)
        # Plot the dendrogram and the module colors underneath
        plotDendroAndColors(nethybrid$dendrograms[[1]], mergedColors[nethybrid$blockGenes[[1]]],
        #"Module colors",
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05)
        dev.off();

        #note that if the user would like to change some of the tree cut, module
        #membership, and module merging criteria, the package provides the function recutBlockwiseTrees that can apply
        #modified criteria without having to recompute the network and the clustering dendrogram. 

        # save the module assignment and module eigengene information necessary for subsequent analysis
        moduleLabels = nethybrid$colors
        moduleColors = labels2colors(nethybrid$colors)

        # output the number of gene per module created
        genes.mod = as.matrix(table(moduleColors))
        file = paste("WGCNA/2__nbGenesPerModule.power_", i, ".mergeCutHeight_", y, ".csv", sep="")
        colnames(genes.mod) = "gene_nb"
        write.csv(genes.mod, file = file)

        MEs = nethybrid$MEs;
        geneTree = nethybrid$dendrograms[[1]];
        file = paste("WGCNA/2-networkConstruction-auto.power_", i, ".mergeCutHeight_", y, ".RData", sep="")
        save(nethybrid, MEs, moduleLabels, moduleColors, geneTree, file = file)

        
        ### PART3
        
        nGenes = ncol(datExpr)
        nSamples = nrow(datExpr)
        # Recalculate MEs with color labels
        MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
        MEs = orderMEs(MEs0)
        modNames = substring(names(MEs), 3)
        geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))

        ## 3.i Summary output of network analysis results (3.d FemaleLiver-03-relateModsToExt.pdf)

        # output genes membership per module
        membership.df <- data.frame()

        # output TPM per sample per module
        module_exp_all.df <- data.frame()

        # output avg and median TPM per sample per module
        module_exp_avg.df <- data.frame()                            

        for (module in modNames){
            column = match(module, modNames);
            moduleGenes = moduleColors==module;
            # gene name annotation
            #probes2annot = match(colnames(datExpr[,moduleColors==module ]), annot$gene_id)
            membership_module.df <- data.frame(module,
                                            colnames(datExpr[,moduleColors==module ]),
            #                                annot$Gene[probes2annot],
                                            abs(geneModuleMembership[moduleGenes, column])) 
            membership.df <- rbind(membership.df, membership_module.df)

            module_exp.df <- as.data.frame(t(datExpr[,(names(datExpr) %in% colnames(datExpr[,moduleColors==module ]))]))
            module_exp_all.df <- rbind(module_exp_all.df, membership_module.df$module,module_exp.df)

            # variance
                for (sample in names(module_exp.df))
                {
                module_exp_avg_sample.df <- data.frame(sample,
                                                module,
                                                mean(module_exp.df[[sample]]),
                                                median(module_exp.df[[sample]]))
                module_exp_avg.df <- rbind(module_exp_avg.df, module_exp_avg_sample.df)
                }

        }
        colnames(membership.df)=c("module","id","membership")
        colnames(module_exp_all.df)=c(rownames(datExpr[,moduleColors==module ]))
        colnames(module_exp_avg.df)=c("sample","module","mean","median")
        file = paste("WGCNA/3__geneMembershipPerModule.power_", i, ".mergeCutHeight_", y, ".csv", sep="")
        write.csv(membership.df, file = file)
        file = paste("WGCNA/3__TpmPerSamplePerModule.power_", i, ".mergeCutHeight_", y, ".csv", sep="")
        write.csv(module_exp_all.df, file = file)
        file = paste("WGCNA/3__avgTpmPerSamplePerModule.power_", i, ".mergeCutHeight_", y, ".csv", sep="")
        write.csv(module_exp_avg.df, file = file)

        # output eigengene expression per module per sample
        eigengene_exp.df <- data.frame(module=character(),
                                    sample=character(),
                                    eigengene_exp=character(),
                                    stringsAsFactors=FALSE)
        for (module in modNames)
        {
        ME=MEs0[, paste("ME",module, sep="")]
        eigengene_exp_module.df <- data.frame(module,
                                            rownames(datExpr),
                                            ME) 
        eigengene_exp.df <- rbind(eigengene_exp.df, eigengene_exp_module.df)
        }                 
        colnames(eigengene_exp.df)=c("module","sample","eigengene_exp")
        file = paste("WGCNA/3__sampleEigengeneExpPerModule.power_", i, ".mergeCutHeight_", y, ".csv", sep="")
        write.csv(eigengene_exp.df, file = file)

    }
}


# Calculate MEs with color labels
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs)
# method = spearman ! see also line 98, 102 and 115
moduleDrugCor = cor(MEs, datDrug, use = "p", method = "spearman")
moduleDrugPvalue = corPvalueStudent(moduleDrugCor, nrow(datExpr))
moduleDrugCorPvalue = cbind(moduleDrugCor, moduleDrugPvalue)
moduleDrugCorPvalue = as.data.frame(moduleDrugCorPvalue)
names(moduleDrugCorPvalue) = c("corr_BL_DC","corr_AGE", "corr_SEX", "corr_DIABETES","corr_DYSLIPIDEMIA", "corr_TIMEPOINTS","corr_TREATMENT","pval_BL_DC","pval_AGE","pval_SEX","pval_DIABETES","pval_DYSLIPIDEMIA","pval_TIMEPOINTS","pval_TREATMENT") 
write.csv(moduleDrugCorPvalue, file = "WGCNA/moduleDrugCorPvalue.csv")

# Output eigengene expression per module per sample
eigengene_exp.df <- data.frame( module=character(),
                                sample=character(),
                                eigengene_exp=character(),
                                stringsAsFactors=FALSE)

modNames = substring(names(MEs), 3)
for (module in modNames){
    ME=MEs[, paste("ME",module, sep="")]
    eigengene_exp_module.df <- data.frame(module,
                                          rownames(datExpr),
                                          ME) 
    eigengene_exp.df <- rbind(eigengene_exp.df, eigengene_exp_module.df)
}  
colnames(eigengene_exp.df)=c("module","sample","eigengene_exp")
write.csv(eigengene_exp.df, file = "WGCNA/sampleEigengeneExpPerModule.csv")


## Gene relationship to drug

# Get module membership
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p", method = "spearman"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneDrugSignificance = as.data.frame(cor(datExpr, datDrug, use = "p", method = "spearman")) #interesting table
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneDrugSignificance), nrow(datExpr)))
names(geneDrugSignificance) = paste("GS.", names(datDrug), sep="")
names(GSPvalue) = paste("p.GS.", names(datDrug), sep="")

# Gene name annotation
geneInfo0 = data.frame( geneSymbol = names(datExpr),
                        moduleColor = moduleColors,
                        geneDrugSignificance,
                        GSPvalue)

# Order modules by their significance for datDrug
modOrder = order(-abs(cor(MEs, datDrug, use = "p", method = "spearman")))

# Add module membership information
for (mod in 1:ncol(geneModuleMembership)){
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[, modOrder[mod]])
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
    paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneDrugSignificance
GS.drug=paste("GS.", names(datDrug), sep="")
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0[[GS.drug]]))
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "WGCNA/geneInfo.csv")


## FIGURES

# Recalculate MEs with color labels
# plot p-val for each module-trait assoc
textMatrix = paste(signif(moduleDrugCor, 2), " (",signif(moduleDrugPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleDrugCor)
png("WGCNA/3__moduleTraitAssociations.png", width=1800, height=2500, res=200)
par(mar = c(6, 5.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleDrugCor,
xLabels = names(datDrug),
yLabels = names(MEs[,1-3]),
ySymbols = gsub('ME','',colnames(MEs)),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.8,
cex.lab.x = 0.8,
zlim = c(-1,1),
main = "Module-trait relationships in all normalized data")
dev.off();

# Recalculate MEs with color labels
# plot p-val for each module-trait assoc
textMatrix = paste(signif(moduleDrugCor[1:3,], 2), " (",signif(moduleDrugPvalue[1:3,], 1), ")", sep = "")
dim(textMatrix) = dim(moduleDrugCor[1:3,])
png("WGCNA/3__moduleTraitAssociations_noGrey.png", width=1800, height=2500, res=200)
par(mar = c(6, 5.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleDrugCor[1:3,],
xLabels = names(datDrug),
yLabels = names(MEs[,1:3]),
ySymbols = gsub('ME','',colnames(MEs[,1:3])),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.8,
cex.lab.x = 0.8,
zlim = c(-1,1),
main = "Module-trait relationships in all normalized data")
dev.off();

# Convert traits to a color representation: white means low, red means high, grey means missing entry
drugColors = numbers2colors(datDrug, signed = FALSE, colors=blueWhiteRed(100), naColor = "grey")
keep = c("blue","turquoise","brown")
for (module in keep)
{
    fig = paste("WGCNA/3__module_heatmap_eigengene.", module, "_orderedForBL_DC.png", sep="")
    png(fig, width=5000, height=1500, res=200)
    ME=MEs[, paste("ME",module, sep="")]
    par(mfrow=c(2,1), mar=c(0.5, 9.5, 7, 6))
    plotMat(t(scale(datExpr[,moduleColors==module ]) ),
    nrgcols=20,rlabels=F,
    clabels=rownames(datExpr),
    ccols="black", cex.main=1)
    par(mar=c(.7, 5.2, 3, 1.6))
    barplot(ME, col=module, main="", cex.main=1,
    ylab="eigenmetabolite expression",xlab="")
    dev.off();
}

## Hub-metabolite
hubmetabo <- chooseTopHubInEachModule(
   datExpr, 
   moduleColors, 
   omitColors = "grey", 
   power = 4, 
   type = "signed"
   )

#hubmetabo
#>                  blue                    brown                turquoise 
#         "phenylalanine"         "pelargonate 90" "10-nonadecenoate 191n9" 
