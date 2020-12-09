suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(maptools))
suppressMessages(library(plyr))
#38761d dark green
#ff9900 orange

#LOAD R DATA FILE OR TABLES

#fig2A
load(file="fig2_figS2_pca.RData")
pdf(file = "pca_imputed_norm_all_newColors_larger_noLegend.pdf", width = 7, height = 7)
ggplot(as.data.frame(pca_norm_all$x),aes(PC1,PC2))+geom_point(aes(shape=as.factor(meta_norm_all$TREATMENT),colour=as.factor(meta_norm_all$TIME),size=1.4))+ scale_color_manual(values=c("#38761d", "#ff9900")) +xlab("PC1 (18.12%)")+ylab("PC2 (6.23%)")+ theme_bw()+theme(legend.position = "none")
dev.off()
pdf(file = "pca_imputed_norm_all_newColors_larger.pdf", width = 7, height = 7)
ggplot(as.data.frame(pca_norm_all$x),aes(PC1,PC2))+geom_point(aes(shape=as.factor(meta_norm_all$TREATMENT),colour=as.factor(meta_norm_all$TIME),size=1.4))+ scale_color_manual(values=c("#38761d", "#ff9900")) +xlab("PC1 (18.12%)")+ylab("PC2 (6.23%)")+ theme_bw()
dev.off()

#fig2B
table_v <- read.table(file="input/FC-lmer_pval-anova_input_reversedFC.csv", sep = "\t", header = TRUE)
turquoise <- read.table(file="input/turquoise_module.csv", sep="\t", header=T)
blue <- read.table(file="input/blue_module.csv", sep="\t", header=T)
grey <- read.table(file="input/grey_module.csv", sep="\t", header=T)
brown <- read.table(file="input/brown_module.csv", sep="\t", header=T)
brown_plot <- merge(table_v, brown, by="Metabolites")
turquoise_plot <- merge(table_v, turquoise, by="Metabolites")
blue_plot <- merge(table_v, blue, by="Metabolites")
grey_plot <- merge(table_v, grey, by="Metabolites")
pdf("fig2B.pdf", width=7, height=7)
plot(table_v$log2_FC, table_v$log_pval, xlab="log2(FC)", ylab="-log10(p-value)", type="n", xlim=c(-2,2))
points(table_v$log2_FC, table_v$log_pval, col="black", bg="gray48", pch=21, cex=1)
sel1 <- which((table_v$log2_FC<=-0.6) & (table_v$fdr_BL_DC<0.01)) 
sel2 <- which((table_v$log2_FC>=0.6) & (table_v$fdr_BL_DC<0.01))
sel3 <- brown_plot[brown_plot$Metabolites %in% table_v$Metabolites, ] #output the brown metabolites
sel4 <- blue_plot[blue_plot$Metabolites %in% table_v$Metabolites, ] #output the blue metabolites
sel5 <- turquoise_plot[turquoise_plot$Metabolites %in% table_v$Metabolites, ] #output the turquoise metabolites
sel7 <- grey_plot[grey_plot$Metabolites %in% table_v$Metabolites, ] #output the grey metabolites
points(sel4$log2_FC, sel4$log_pval_3way_BL_DC, col="black", bg="blue", pch=21, cex=1)
points(sel5$log2_FC, sel5$log_pval_3way_BL_DC, col="black", bg="turquoise", pch=21, cex=1)
points(sel3$log2_FC, sel3$log_pval_3way_BL_DC, col="black", bg="brown", pch=21, cex=1)
abline(v=-0.6,lty=2)
abline(v=0.6,lty=2)
abline(h=2.099, lty=2)
dev.off()


#fig3A
table <- read.table(file="input/output_logit.csv", sep = "\t", header = TRUE)
pdf("glm_ticaVSclopi_noLOG_withGreenMets.pdf", width=7, height=7)
plot(table$log2_FC, table$log_pval, xlab="log2(FC)", ylab="-log10(p-value)", main="Ticagrelor vs Clopidogrel", type="n", xlim=c(-2,2))
points(table$log2_FC, table$log_pval, col="gray48", bg="gray48", pch=21, cex=1.3)
sel2 <- which((table$log2_FC>=0.6) & (table$fdr<0.05))
#age x state
sel4 <- which(table$Metabolites == "pregnen-diol_disulfate" | table$Metabolites == "21-hydroxypregnenolone_disulfate" | table$Metabolites == "alpha-hydroxyisovalerate" | table$Metabolites == "gamma-glutamylglutamine" | table$Metabolites == "glutamine" | table$Metabolites == "glycerophosphorylcholine GPC" | table$Metabolites == "lysine")
#time x state
sel5 <- which(table$Metabolites == "1-arachidonylglycerol" | table$Metabolites == "1-docosahexaenoylglycerol_1-monodocosahexaenoin" | table$Metabolites == "5-oxoproline" | table$Metabolites == "cortisol" | table$Metabolites == "1-palmitoylglycerophosphate")
#treatment x state put symbol (square) but not color red nor green
sel6 <- which(table$Metabolites == "andro_steroid_monosulfate_2" | table$Metabolites == "2-hydroxybutyrate_AHB")
#green mets in red (omega) and in treatment (square)
sel7 <- which(table$Metabolites == "eicosapentaenoate_EPA_205n3" | table$Metabolites == "linolenate_alpha_or_gamma_183n3_or_6" | table$Metabolites == "arachidonate_204n6" | table$Metabolites == "dihomo-linolenate_203n3_or_n6" | table$Metabolites == "docosahexaenoate_DHA_226n3" | table$Metabolites == "stearidonate_184n3" )
#green mets in red (omega) but circle (because not in treatment)
sel8 <- which(table$Metabolites=="dihomo-linoleate_202n6" | table$Metabolites=="linoleate_182n6" | table$Metabolites=="docosapentaenoate_n6_DPA_225n6" |  table$Metabolites=="docosapentaenoate_n3_DPA_225n3" | table$Metabolites=="adrenate_224n6" )
#green mets in time (diamond)
sel9 <- which(table$Metabolites == "glycerol" | table$Metabolites == "margarate_170" | table$Metabolites == "myristate_140" | table$Metabolites == "myristoleate_141n5" | table$Metabolites == "nonadecanoate_190") 
#green mets in age (triangle)
sel10 <- which(table$Metabolites =="1-arachidonoylglycerophosphocholine" | table$Metabolites == "1-eicosatrienoylglycerophosphocholine" | table$Metabolites == "1-linoleoylglycerophosphocholine") 
#rest of green metabolites (gray circle green contour)
sel11 <- which(table$Metabolites == "1-arachidonoylglycerophosphoinositol" | table$Metabolites == "1-docosahexaenoylglycerophosphocholine" | table$Metabolites == "1-eicosadienoylglycerophosphocholine" | table$Metabolites == "1-heptadecanoylglycerophosphocholine" | table$Metabolites == "1-myristoylglycerophosphocholine" | table$Metabolites == "1-oleoylglycerophosphocholine" | table$Metabolites == "1-palmitoleoylglycerophosphocholine" | table$Metabolites == "1-palmitoylglycerophosphocholine" | table$Metabolites == "1-stearoylglycerophosphocholine" | table$Metabolites == "10-heptadecenoate_171n7" | table$Metabolites == "10-nonadecenoate_191n9" | table$Metabolites == "13-HODE_9-HODE" | table$Metabolites == "16-hydroxypalmitate" | table$Metabolites == "17-methylstearate" | table$Metabolites == "2-palmitoylglycerophosphocholine" | table$Metabolites == "2-stearoylglycerophosphocholine" | table$Metabolites == "5-dodecenoate_121n7" | table$Metabolites == "arachidate_200" | table$Metabolites == "caprate_100" | table$Metabolites == "cis-vaccenate_181n7" | table$Metabolites == "eicosenoate_201n9_or_11" | table$Metabolites == "laurate_120" | table$Metabolites == "oleate_181n9" | table$Metabolites == "palmitate_160" | table$Metabolites == "palmitoleate_161n7" | table$Metabolites == "pentadecanoate_150" | table$Metabolites == "stearate_180")
points(table$log2_FC[sel4], table$log_pval[sel4], col="black",pch=24, cex=1.4) #age x state triangle
points(table$log2_FC[sel5], table$log_pval[sel5], col="black",pch=23, cex=1.4) #time x state diamond
points(table$log2_FC[sel6], table$log_pval[sel6], col = "black", pch=22, cex=1.4) #treatment x state square
points(table$log2_FC[sel7], table$log_pval[sel7], col= "darkgreen", bg = "red", pch=22, cex=1.4) 
points(table$log2_FC[sel8], table$log_pval[sel8], col="darkgreen", bg="red", pch=21, cex=1.4) 
points(table$log2_FC[sel9], table$log_pval[sel9], col="darkgreen", bg="gray48", pch=23, cex=1.4) 
points(table$log2_FC[sel10], table$log_pval[sel10], col="darkgreen", bg="gray48", pch=24, cex=1.4) 
points(table$log2_FC[sel11], table$log_pval[sel11], col="darkgreen", bg="gray48", pch=21, cex=1.4) 
abline(v=-0.6,lty=2)
abline(v=0.6,lty=2)
abline(h=3.267, lty=2, col = "blue") #5%
abline(h=2.094, lty=2, col = "darkgreen") #15%
legend("topleft", inset=c(0.03,0.03), legend=c("Omegas 3 and 6 family", "tSNE green cluster", "Age x State interaction", "Timepoint x State interaction", "Treatment x State interaction"), col=c("red", "darkgreen", "black", "black", "black"), pch=c(21,21,24,23,22), pt.bg=c("red", "white", "white", "white", "white"), cex=1, box.lty=0)
dev.off()

pdf("fig3A_noLegend.pdf", width=7, height=7)
plot(table$log2_FC, table$log_pval, xlab="log2(FC)", ylab="-log10(p-value)", type="n", xlim=c(-2,2))
points(table$log2_FC, table$log_pval, col="gray48", bg="gray48", pch=21, cex=1.3)
sel2 <- which((table$log2_FC>=0.6) & (table$fdr<0.05))
#age x state
sel4 <- which(table$Metabolites == "pregnen-diol_disulfate" | table$Metabolites == "21-hydroxypregnenolone_disulfate" | table$Metabolites == "alpha-hydroxyisovalerate" | table$Metabolites == "gamma-glutamylglutamine" | table$Metabolites == "glutamine" | table$Metabolites == "glycerophosphorylcholine GPC" | table$Metabolites == "lysine")
#time x state
sel5 <- which(table$Metabolites == "1-arachidonylglycerol" | table$Metabolites == "1-docosahexaenoylglycerol_1-monodocosahexaenoin" | table$Metabolites == "5-oxoproline" | table$Metabolites == "cortisol" | table$Metabolites == "1-palmitoylglycerophosphate")
#treatment x state put symbol (square) but not color red nor green
sel6 <- which(table$Metabolites == "andro_steroid_monosulfate_2" | table$Metabolites == "2-hydroxybutyrate_AHB")
#green mets in red (omega) and in treatment (square)
sel7 <- which(table$Metabolites == "eicosapentaenoate_EPA_205n3" | table$Metabolites == "linolenate_alpha_or_gamma_183n3_or_6" | table$Metabolites == "arachidonate_204n6" | table$Metabolites == "dihomo-linolenate_203n3_or_n6" | table$Metabolites == "docosahexaenoate_DHA_226n3" | table$Metabolites == "stearidonate_184n3" )
#green mets in red (omega) but circle (because not in treatment)
sel8 <- which(table$Metabolites=="dihomo-linoleate_202n6" | table$Metabolites=="linoleate_182n6" | table$Metabolites=="docosapentaenoate_n6_DPA_225n6" |  table$Metabolites=="docosapentaenoate_n3_DPA_225n3" | table$Metabolites=="adrenate_224n6" )
#green mets in time (diamond)
sel9 <- which(table$Metabolites == "glycerol" | table$Metabolites == "margarate_170" | table$Metabolites == "myristate_140" | table$Metabolites == "myristoleate_141n5" | table$Metabolites == "nonadecanoate_190") 
#green mets in age (triangle)
sel10 <- which(table$Metabolites =="1-arachidonoylglycerophosphocholine" | table$Metabolites == "1-eicosatrienoylglycerophosphocholine" | table$Metabolites == "1-linoleoylglycerophosphocholine") 
#rest of green metabolites (gray circle green contour)
sel11 <- which(table$Metabolites == "1-arachidonoylglycerophosphoinositol" | table$Metabolites == "1-docosahexaenoylglycerophosphocholine" | table$Metabolites == "1-eicosadienoylglycerophosphocholine" | table$Metabolites == "1-heptadecanoylglycerophosphocholine" | table$Metabolites == "1-myristoylglycerophosphocholine" | table$Metabolites == "1-oleoylglycerophosphocholine" | table$Metabolites == "1-palmitoleoylglycerophosphocholine" | table$Metabolites == "1-palmitoylglycerophosphocholine" | table$Metabolites == "1-stearoylglycerophosphocholine" | table$Metabolites == "10-heptadecenoate_171n7" | table$Metabolites == "10-nonadecenoate_191n9" | table$Metabolites == "13-HODE_9-HODE" | table$Metabolites == "16-hydroxypalmitate" | table$Metabolites == "17-methylstearate" | table$Metabolites == "2-palmitoylglycerophosphocholine" | table$Metabolites == "2-stearoylglycerophosphocholine" | table$Metabolites == "5-dodecenoate_121n7" | table$Metabolites == "arachidate_200" | table$Metabolites == "caprate_100" | table$Metabolites == "cis-vaccenate_181n7" | table$Metabolites == "eicosenoate_201n9_or_11" | table$Metabolites == "laurate_120" | table$Metabolites == "oleate_181n9" | table$Metabolites == "palmitate_160" | table$Metabolites == "palmitoleate_161n7" | table$Metabolites == "pentadecanoate_150" | table$Metabolites == "stearate_180")
points(table$log2_FC[sel4], table$log_pval[sel4], col="black",pch=24, cex=1.4) #age x state triangle
points(table$log2_FC[sel5], table$log_pval[sel5], col="black",pch=23, cex=1.4) #time x state diamond
points(table$log2_FC[sel6], table$log_pval[sel6], col = "black", pch=22, cex=1.4) #treatment x state square
points(table$log2_FC[sel7], table$log_pval[sel7], col= "darkgreen", bg = "red", pch=22, cex=1.4) 
points(table$log2_FC[sel8], table$log_pval[sel8], col="darkgreen", bg="red", pch=21, cex=1.4) 
points(table$log2_FC[sel9], table$log_pval[sel9], col="darkgreen", bg="gray48", pch=23, cex=1.4) 
points(table$log2_FC[sel10], table$log_pval[sel10], col="darkgreen", bg="gray48", pch=24, cex=1.4) 
points(table$log2_FC[sel11], table$log_pval[sel11], col="darkgreen", bg="gray48", pch=21, cex=1.4) 
abline(v=-0.6,lty=2)
abline(v=0.6,lty=2)
abline(h=3.267, lty=2, col = "blue") #5%
abline(h=2.094, lty=2, col = "darkgreen") #15%
dev.off()

#fig3B
median_sous <- read.table(file="input/norm_values_MEDIAN_sous_tig_clo_n3n6_topFA_AA_reorder_abs.csv", sep="\t", header=T)
rownT<-median_sous$subject_BL_DC
median_sous <- as.data.frame(median_sous, row.names=as.character(rownT))
median_sous$subject_BL_DC <- NULL


coul =c("#00BFC4","#F8766D") #blue and red
colors_border <- coul
colors_in <- alpha(coul,0.1)

png("median_ticaVSclo_n3n6_topFA_AA_absoluteValue_reorder_red_blue_v2.png", width=2500, height=2000, res=200)
radarchart( median_sous  , axistype=1 , centerzero = TRUE, seg = 4,
             pcol=colors_border , 
             pfcol=colors_in , plwd=4 , plty=1,
             cglcol="black", cglty=1, axislabcol="black", caxislabels=c("0", "0.5", "1", "1.5", "2"), cglwd=0.8, 
             vlcex=0.8 
 )
legend(x=1, y=1, legend = rownames(median_sous[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=1)
dev.off()

png("median_ticaVSclo_n3n6_topFA_AA_absoluteValue_reorder_red_blue_v2_noLegend.png", width=2500, height=2000, res=200)
radarchart( median_sous  , axistype=1 , centerzero = TRUE, seg = 4,
             pcol=colors_border , 
             pfcol=colors_in , plwd=4 , plty=1,
             cglcol="black", cglty=1, axislabcol="black", caxislabels=c("0", "0.5", "1", "1.5", "2"), cglwd=0.8, 
             vlcex=0.8 
 )
dev.off()



#figS2D
BL <- read.delim(file="input/onlyBL_notImputed_woOutliers.csv",sep="\t", header=T)
DC <- read.delim(file="input/onlyDC_notImputed_woOutliers.csv",sep="\t", header=T)
#https://stackoverflow.com/questions/37801338/r-count-nas-per-row-in-dataframe
BL$na_count <- apply(BL, 1, function(x) sum(is.na(x)))
DC$na_count <- apply(DC, 1, function(x) sum(is.na(x)))
#now we have the number of missing values (NA) in each df, for each metabolite
BL$percent <- (BL$na_count/175)*100
DC$percent <- (DC$na_count/175)*100
png("percent_BL_vs_DC_missingValues.png", width=1500, height=1500, res=200)
plot(BL$percent,DC$percent,type="n")
sel1 <- which((BL$percent<10)&(DC$percent<10))
points(BL$percent,DC$percent,col="gray48")
points(BL$percent[sel1],DC$percent[sel1],col="black")
abline(h=10,lty=2,col="red")
abline(v=10,lty=2,col="red")
dev.off()

cor(BL$percent,DC$percent,method="pearson")
#[1] 0.9428679
cor(BL$percent,DC$percent,method="spearman")
#[1] 0.9285264


#figS4C
library(ggplot2)
library(reshape)
all_M <- read.table(file="input/data_normalized_all_BL_DC_withmetadata.csv", sep="\t", header=T)
rown<-all_M$SAMPLE_NAME
#age*status
metabos <- c("X1.arachidonoylglycerophosphocholine","X1.eicosatrienoylglycerophosphocholine","X1.linoleoylglycerophosphocholine","X21.hydroxypregnenolone.disulfate","alpha.hydroxyisovalerate","gamma.glutamylglutamine","glutamine","glycerophosphorylcholine.GPC","lysine","pregnen.diol.disulfate")

all_M$bins <- cut(all_M$AGE, 5)
t <- as.data.frame(all_M[,metabos], row.names=as.character(rown))
t$TIME <- all_M$TIME
t$TREATMENT <- all_M$TREATMENT
t$bins <- all_M$bins
t_melt <- melt(t)
t_melt$type <- paste(t_melt$TIME,t_melt$bins,sep="_")

#with dots
png("ageXstatus_metabolites_withDots_bins.png", width=3500, height=1700, res=200)
p <- ggplot(t_melt, aes(x=variable, y=value, fill=type)) + geom_boxplot(alpha=0.80) + geom_point(aes(fill=type),size=2, shape=21, position=position_jitterdodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) +ggtitle("Normalized values of age:status significant metabolites")+xlab("Metabolites")
p
dev.off()
#without dots
pdf("figS4C.pdf", width=10, height=7)
p <- ggplot(t_melt, aes(x=variable, y=value, fill=type)) + geom_boxplot(alpha=0.80) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) +xlab("Metabolites")
p
dev.off()
pdf("figS4C_noLegend.pdf", width=10, height=7)
p <- ggplot(t_melt, aes(x=variable, y=value, fill=type)) + geom_boxplot(alpha=0.80) +xlab("Metabolites") + theme_bw() + theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1, size=10))
p
dev.off()

#figS5
load(file="interactions.RData")
plot_diff_ticaclo <- function(tica_22_6, clo_22_6, col1="#63C5CA", col2="#e79289",my_y=c(-4,4),name="metabolite",yshf=0.2)
{
plot(c(0,1),tica_22_6[1,c(1,2)],type="l",xlim=c(-yshf,1+yshf), ylim=my_y,col="white",lwd=1,lty=2,axes=FALSE,ylab=name,xlab="")

for(i in 1:nrow(tica_22_6)){
	lines(c(0,1),tica_22_6[i,c(1,2)],col=col1,lwd=0.1,lty=1)
}
for(i in 1:nrow(clo_22_6)){
	lines(c(0,1),clo_22_6[i,c(1,2)],col=col2,lwd=0.1,lty=1)
}

for(i in 1:nrow(tica_22_6)){
	if(tica_22_6$diff[i] >0){
		lines(c(0,1),tica_22_6[i,c(1,2)],col=col1,lwd=0.5,lty=2)
	}
}
lines(c(0,1),c(mean(tica_22_6[,1]),mean(tica_22_6[,2])),col=col1,lwd=3)
for(i in 1:nrow(clo_22_6)){
	if(clo_22_6$diff[i] >0){
		lines(c(0,1),clo_22_6[i,c(1,2)],col=col2,lwd=0.5,lty=2)
	}
}
lines(c(0,1),c(mean(clo_22_6[,1]),mean(clo_22_6[,2])),col=col2,lwd=3) 

axis(1,at=c(0,1),labels=c("BL","DC"),pos=my_y[1])
axis(2, pos=-yshf)
}
pdf("figS5_n3.pdf", width=20, height=20)
par(mfrow=c(1,3))
plot_diff_ticaclo(tica_22_6, clo_22_6,name="C22:6",my_y=c(-3,4),yshf=0.3)
plot_diff_ticaclo(tica_20_5, clo_20_5,name="C20:5",my_y=c(-3,4),yshf=0.3)
plot_diff_ticaclo(tica_18_4, clo_18_4,name="C18:4",my_y=c(-3,4),yshf=0.3)
dev.off()

all_20_4 <- read.table(file="input/arachidonate_204n6.csv",sep="\t",header=T,row.names=1)
tica_20_4 <- all_20_4[(all_20_4$treatment=="1"),]
tica_20_4$time <- NULL
tica_20_4$pairing <- NULL
tica_20_4$treatment <- NULL
tica_20_4$diabetes <- NULL
tica_20_4$dyslipidemia <- NULL
tica_20_4$sex <- NULL
tica_20_4$age <- NULL
tica_20_4$diff <- tica_20_4$DC_value - tica_20_4$BL_value
clo_20_4 <- all_20_4[(all_20_4$treatment=="0"),]
clo_20_4$time <- NULL
clo_20_4$pairing <- NULL
clo_20_4$treatment <- NULL
clo_20_4$diabetes <- NULL
clo_20_4$dyslipidemia <- NULL
clo_20_4$sex <- NULL
clo_20_4$age <- NULL
clo_20_4$diff <- clo_20_4$DC_value - clo_20_4$BL_value

pdf("figS5_n6_full.pdf", width=20, height=20)
par(mfrow=c(1,3))
plot_diff_ticaclo(tica_22_4, clo_22_4,name="C22:4",my_y=c(-3,4),yshf=0.3)
plot_diff_ticaclo(tica_20_2, clo_20_2,name="C20:2",my_y=c(-3,4),yshf=0.3)
plot_diff_ticaclo(tica_18_2, clo_18_2,name="C18:2",my_y=c(-3,4),yshf=0.3)
plot_diff_ticaclo(tica_20_4, clo_20_4,name="C20:4",my_y=c(-3,4),yshf=0.3)
dev.off()


#figS7A
load(file="fig4_eigenvalues.RData")
#### timpoints T0 vs T1 vs T2 ####
blue_BL_T0 <- blue_info[(blue_info$TIME=="BL")&(blue_info$ONSET=="T0"),]
blue_BL_T1 <- blue_info[(blue_info$TIME=="BL")&(blue_info$ONSET=="T1"),]
blue_BL_T2 <- blue_info[(blue_info$TIME=="BL")&(blue_info$ONSET=="T2"),]
blue_DC_T0 <- blue_info[(blue_info$TIME=="DC")&(blue_info$ONSET=="T0"),]
blue_DC_T1 <- blue_info[(blue_info$TIME=="DC")&(blue_info$ONSET=="T1"),]
blue_DC_T2 <- blue_info[(blue_info$TIME=="DC")&(blue_info$ONSET=="T2"),]
turquoise_BL_T0 <- turquoise_info[(turquoise_info$TIME=="BL")&(turquoise_info$ONSET=="T0"),]
turquoise_BL_T1 <- turquoise_info[(turquoise_info$TIME=="BL")&(turquoise_info$ONSET=="T1"),]
turquoise_BL_T2 <- turquoise_info[(turquoise_info$TIME=="BL")&(turquoise_info$ONSET=="T2"),]
turquoise_DC_T0 <- turquoise_info[(turquoise_info$TIME=="DC")&(turquoise_info$ONSET=="T0"),]
turquoise_DC_T1 <- turquoise_info[(turquoise_info$TIME=="DC")&(turquoise_info$ONSET=="T1"),]
turquoise_DC_T2 <- turquoise_info[(turquoise_info$TIME=="DC")&(turquoise_info$ONSET=="T2"),]
brown_BL_T0 <- brown_info[(brown_info$TIME=="BL")&(brown_info$ONSET=="T0"),]
brown_BL_T1 <- brown_info[(brown_info$TIME=="BL")&(brown_info$ONSET=="T1"),]
brown_BL_T2 <- brown_info[(brown_info$TIME=="BL")&(brown_info$ONSET=="T2"),]
brown_DC_T0 <- brown_info[(brown_info$TIME=="DC")&(brown_info$ONSET=="T0"),]
brown_DC_T1 <- brown_info[(brown_info$TIME=="DC")&(brown_info$ONSET=="T1"),]
brown_DC_T2 <- brown_info[(brown_info$TIME=="DC")&(brown_info$ONSET=="T2"),]

plot_colors <- c("turquoise", "turquoise", "turquoise", "blue", "blue", "blue", "brown", "brown", "brown")
plot_type <- c("l","l","l","l","l","l","l","l","l")
plot_lty <- c(1,2,3,1,2,3,1,2,3)
noms <- c("T0: 0-2h", "T1: 2h-4h", "T2: 4h-12h", "T0: 0-2h", "T1: 2h-4h", "T2: 4h-12h", "T0: 0-2h" , "T1: 2h-4h", "T2: 4h-12h")

#BL
pdf("figS7A.pdf", width=7, height=7)
plot(density(turquoise_BL_T0$eigenvalue), type=plot_type[1], col=plot_colors[1], lty=plot_lty[1], lwd =2, xlab = c("Eigenvalues"), xlim=c(-0.20,0.20), ylim=c(0,18)) 
lines(density(turquoise_BL_T1$eigenvalue), type=plot_type[2], col=plot_colors[2], lty=plot_lty[2], lwd =2) 
lines(density(turquoise_BL_T2$eigenvalue), type=plot_type[3], col=plot_colors[3], lty=plot_lty[3], lwd =2) 
lines(density(blue_BL_T0$eigenvalue), type=plot_type[4], col=plot_colors[4], lty=plot_lty[4], lwd =2) 
lines(density(blue_BL_T1$eigenvalue), type=plot_type[5], col=plot_colors[5], lty=plot_lty[5], lwd =2)
lines(density(blue_BL_T2$eigenvalue), type=plot_type[6], col=plot_colors[6], lty=plot_lty[6], lwd =2)
lines(density(brown_BL_T0$eigenvalue), type=plot_type[7], col=plot_colors[7], lty=plot_lty[7], lwd =2) 
lines(density(brown_BL_T1$eigenvalue), type=plot_type[8], col=plot_colors[8], lty=plot_lty[8], lwd =2) 
lines(density(brown_BL_T2$eigenvalue), type=plot_type[9], col=plot_colors[9], lty=plot_lty[9], lwd =2)   
legend("topleft",noms, border="black",cex=0.8, col=plot_colors, lty =plot_lty, lwd =2)
dev.off()

#DC
png("blue_brown_turquoise_T0vsT1vsT2_DC.png", width=1800, height=1500, res=200)
plot(density(turquoise_DC_T0$eigenvalue), type=plot_type[1], col=plot_colors[1], lty=plot_lty[1], lwd =2, main=c("Blue, Turquoise and Brown modules - T0 + T1 vs T2 at DC"), xlab = c("Eigenvalues"), xlim=c(-0.20,0.20), ylim=c(0,16)) 
lines(density(turquoise_DC_T1$eigenvalue), type=plot_type[2], col=plot_colors[2], lty=plot_lty[2], lwd =2) 
lines(density(turquoise_DC_T2$eigenvalue), type=plot_type[3], col=plot_colors[3], lty=plot_lty[3], lwd =2) 
lines(density(blue_DC_T0$eigenvalue), type=plot_type[4], col=plot_colors[4], lty=plot_lty[4], lwd =2) 
lines(density(blue_DC_T1$eigenvalue), type=plot_type[5], col=plot_colors[5], lty=plot_lty[5], lwd =2)
lines(density(blue_DC_T2$eigenvalue), type=plot_type[6], col=plot_colors[6], lty=plot_lty[6], lwd =2)
lines(density(brown_DC_T0$eigenvalue), type=plot_type[7], col=plot_colors[7], lty=plot_lty[7], lwd =2) 
lines(density(brown_DC_T1$eigenvalue), type=plot_type[8], col=plot_colors[8], lty=plot_lty[8], lwd =2) 
lines(density(brown_DC_T2$eigenvalue), type=plot_type[9], col=plot_colors[9], lty=plot_lty[9], lwd =2)   
legend("topleft",noms, border="black",cex=0.8, col=plot_colors, lty =plot_lty, lwd =2)
dev.off()


#figS7B
library(ggplot2)
library(reshape)
all_M <- read.table(file="input/data_normalized_all_BL_DC_withmetadata.csv", sep="\t", header=T)
rown<-all_M$SAMPLE_NAME
metabos <- c("X1.arachidonylglycerol", "X1.docosahexaenoylglycerol.1.monodocosahexaenoin", "X1.palmitoylglycerophosphate", "X5.oxoproline", "cortisol", "eicosapentaenoate.EPA.205n3", "glycerol", "linolenate.alpha.or.gamma.183n3.or.6", "margarate.170", "myristate.140", "myristoleate.141n5", "nonadecanoate.190", "pregnen.diol.disulfate")
t <- as.data.frame(all_M[,metabos], row.names=as.character(rown))
t$TIME <- all_M$TIME
t$TREATMENT <- all_M$TREATMENT
t$ONSET <- all_M$ONSET
t_melt <- melt(t)
t_melt$type <- paste(t_melt$TIME,t_melt$ONSET,sep="_")

#with dots
png("timeXstatus_metabolites_withDots.png", width=3500, height=1700, res=200)
p <- ggplot(t_melt, aes(x=variable, y=value, fill=type)) + geom_boxplot(alpha=0.80) + geom_point(aes(fill=type),size=2, shape=21, position=position_jitterdodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) +ggtitle("Normalized values of time:status significant metabolites")+xlab("Metabolites")
p
dev.off()
#without dots
pdf("figS7B.pdf", width=10, height=7)
p <- ggplot(t_melt, aes(x=variable, y=value, fill=type)) + geom_boxplot(alpha=0.80) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) +xlab("Metabolites")+ theme_bw()
p
dev.off()

pdf("figS7B_noLegend.pdf", width=10, height=7)
p <- ggplot(t_melt, aes(x=variable, y=value, fill=type)) + geom_boxplot(alpha=0.80) +xlab("Metabolites")+ theme_bw() + theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1, size=10))
p
dev.off()


#figS7C
all_M <- read.table(file="input/data_normalized_all_BL_DC_withmetadata.csv", sep="\t", header=T)
rown<-all_M$SAMPLE_NAME
#age*status
metabos <- c("pregnen.diol.disulfate")
#all_M$bins <- cut(all_M$AGE, 3)
t <- as.data.frame(all_M[,metabos], row.names=as.character(rown))
t$TIME <- all_M$TIME
t$TREATMENT <- all_M$TREATMENT
t$sex <- all_M$SEX
t$onset <- all_M$ONSET
t_melt <- melt(t)
t_melt$type <- paste(t_melt$TIME,t_melt$onset,sep="_")
t_melt$variable2 <- paste(t_melt$variable,t_melt$sex,sep="_")

#with dots
png("sex_timeXstatus_metabolites_withDots.png", width=3500, height=1700, res=200)
p <- ggplot(t_melt, aes(x=variable2, y=value, fill=type)) + geom_boxplot(alpha=0.80) + geom_point(aes(fill=type),size=2, shape=21, position=position_jitterdodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) +xlab("Metabolites") + theme_bw()
p
dev.off()
#without dots
pdf("figS7C_noLegend.pdf", width=7, height=7)
p <- ggplot(t_melt, aes(x=variable2, y=value, fill=type)) + geom_boxplot(alpha=0.80)+theme_bw() + theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1, size=10))+xlab("Metabolites") 
p
dev.off()

pdf("figS7C.pdf", width=7, height=7)
p <- ggplot(t_melt, aes(x=variable2, y=value, fill=type)) + geom_boxplot(alpha=0.80)+theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10))+xlab("Metabolites") 
p
dev.off()
