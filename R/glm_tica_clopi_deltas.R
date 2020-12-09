library("QuantPsyc")
library(tidyr)
library(car)
library(lme4)


#for all the files
metabos <- c("10-heptadecenoate_171n7" , "10-nonadecenoate_191n9" , "10-undecenoate_111n1" , "1-2-propanediol" , "1-3-dipalmitoylglycerol" , "13-HODE_9-HODE" , "1-5-anhydroglucitol_1-5-AG" , "16-hydroxypalmitate" , "17-methylstearate" , "1-arachidonoylglycerophosphocholine" , "1-arachidonoylglycerophosphoethanolamine" , "1-arachidonoylglycerophosphoinositol" , "1-arachidonoyl_LPA" , "1-arachidonylglycerol" , "1-docosahexaenoylglycerol_1-monodocosahexaenoin" , "1-docosahexaenoylglycerophosphocholine" , "1-eicosadienoylglycerophosphocholine" , "1-eicosatrienoylglycerophosphocholine" , "1-heptadecanoylglycerophosphocholine" , "1-linoleoylglycerol_1-monolinolein" , "1-linoleoylglycerophosphocholine" , "1-linoleoylglycerophosphoethanolamine" , "1-myristoylglycerol_1-monomyristin" , "1-myristoylglycerophosphocholine" , "1-oleoylglycerol_1-monoolein" , "1-oleoylglycerophosphocholine" , "1-oleoylglycerophosphoethanolamine" , "1-palmitoleoylglycerophosphocholine" , "1-palmitoylglycerol_1-monopalmitin" , "1-palmitoylglycerophosphate" , "1-palmitoylglycerophosphocholine" , "1-palmitoylglycerophosphoethanolamine" , "1-palmitoylglycerophosphoinositol" , "1-palmitoylplasmenylethanolamine" , "1-pentadecanoylglycerophosphocholine" , "1-stearoylglycerol_1-monostearin" , "1-stearoylglycerophosphocholine" , "1-stearoylglycerophosphoethanolamine" , "1-stearoylglycerophosphoinositol" , "21-hydroxypregnenolone_disulfate" , "2-aminobutyrate" , "2-ethylhexanoate" , "2-hydroxybutyrate_AHB" , "2-hydroxyglutarate" , "2-hydroxyisobutyrate" , "2-hydroxypalmitate" , "2-hydroxystearate" , "2-linoleoylglycerol_2-monolinolein" , "2-linoleoylglycerophosphoethanolamine" , "2-methylbutyrylcarnitine_C5" , "2-oleoylglycerophosphoethanolamine" , "2-palmitoylglycerophosphocholine" , "2-palmitoylglycerophosphoethanolamine" , "2-stearoylglycerophosphocholine" , "3-4-hydroxyphenyllactate" , "3-carboxy-4-methyl-5-propyl-2-furanpropanoate_CMPF" , "3-dehydrocarnitine" , "3-hydroxy-2-ethylpropionate" , "3-hydroxybutyrate_BHBA" , "3-hydroxyisobutyrate" , "3-hydroxypropanoate" , "3-indoxyl_sulfate" , "3-methoxytyrosine" , "3-methyl-2-oxobutyrate" , "3-methyl-2-oxovalerate" , "4-androsten-3beta-17beta-diol_disulfate_1" , "4-androsten-3beta-17beta-diol_disulfate_2" , "4-methyl-2-oxopentanoate" , "4-methylcatechol_sulfate" , "5alpha-androstan-3beta-17beta-diol_disulfate" , "5-dodecenoate_121n7" , "5-methyluridine_ribothymidine" , "5-oxoproline" , "7-alpha-hydroxy-3-oxo-4-cholestenoate7-Hoca" , "7-methylguanine" , "acetylcarnitine" , "acetylphosphate" , "adrenate_224n6" , "ADSGEGDFXAEGGGVR" , "alanine" , "allantoin" , "alpha-hydroxyisovalerate" , "alpha-ketobutyrate" , "alpha-ketoglutarate" , "alpha-tocopherol" , "andro_steroid_monosulfate_1" , "andro_steroid_monosulfate_2" , "androsterone_sulfate" , "arabitol" , "arabonate" , "arachidate_200" , "arachidonate_204n6" , "arginine" , "asparagine" , "aspartate" , "aspartylphenylalanine" , "benzoate" , "beta-alanine" , "beta-hydroxyisovalerate" , "betaine" , "beta-tocopherol" , "bilirubin_E_E" , "bilirubin_Z_Z" , "biliverdin" , "butyrylcarnitine" , "caffeine" , "caprate_100" , "caproate_60" , "caprylate_80" , "carnitine" , "catechol_sulfate" , "C-glycosyltryptophan" , "cholesterol" , "choline" , "cis-4-decenoyl_carnitine" , "cis-vaccenate_181n7" , "citrate" , "citrulline" , "cortisol" , "cortisone" , "creatine" , "creatinine" , "cysteine" , "cysteine-glutathione_disulfide" , "decanoylcarnitine" , "dehydroisoandrosterone_sulfate_DHEA-S" , "deoxycarnitine" , "dihomo-linoleate_202n6" , "dihomo-linolenate_203n3_or_n6" , "dimethylarginine_SDMA_ADMA" , "dimethylglycine" , "docosahexaenoate_DHA_226n3" , "docosapentaenoate_n3_DPA_225n3" , "docosapentaenoate_n6_DPA_225n6" , "DSGEGDFXAEGGGVR" , "eicosapentaenoate_EPA_205n3" , "eicosenoate_201n9_or_11" , "epiandrosterone_sulfate" , "erythritol" , "erythronate" , "erythrulose" , "fructose" , "fucose" , "fumarate" , "gamma-glutamylalanine" , "gamma-glutamylglutamate" , "gamma-glutamylglutamine" , "gamma-glutamylleucine" , "gamma-glutamylmethionine" , "gamma-glutamylphenylalanine") #1:150
metabos2 <- c("gamma-glutamyltyrosine" , "gamma-glutamylvaline" , "gamma-tocopherol" , "gluconate" , "glucose" , "glucuronate" , "glutamate" , "glutamine" , "glycerate" , "glycerol_2-phosphate" , "glycerol_3-phosphate_G3P" , "glycerol" , "glycerophosphorylcholine_GPC" , "glycine" , "glycochenodeoxycholate" , "glycocholate" , "glycocholenate_sulfate" , "glycolate_hydroxyacetate" , "glycolithocholate_sulfate" , "glycoursodeoxycholate" , "glycylvaline" , "heme" , "heptanoate_70" , "hexadecanedioate" , "hexanoylcarnitine" , "hippurate" , "histidine" , "histidyltryptophan" , "homocitrulline" , "homostachydrine" , "HWESASXX" , "hypoxanthine" , "imidazole_lactate" , "indoleacetate" , "indolelactate" , "indolepropionate" , "inositol_1-phosphate_I1P" , "isobutyrylcarnitine" , "isoleucine" , "isovalerylcarnitine" , "kynurenine" , "lactate" , "laurate_120" , "leucine" , "linoleate_182n6" , "linolenate_alpha_or_gamma_183n3_or_6" , "lysine" , "malate" , "mannitol" , "mannose" , "margarate_170" , "methionine" , "methylphosphate" , "methyl_stearate" , "myo-inositol" , "myristate_140" , "myristoleate_141n5" , "N1-Methyl-2-pyridone-5-carboxamide" , "N1-methyladenosine" , "N6-acetyllysine" , "N-acetylalanine" , "N-acetyl-beta-alanine" , "N-acetylglycine" , "N-acetylneuraminate" , "N-acetylornithine" , "N-acetylserine" , "N-acetylthreonine" , "nicotinamide" , "nonadecanoate_190" , "octadecanedioate" , "octanoylcarnitine" , "oleate_181n9" , "oleoylcarnitine" , "ornithine" , "palmitate_160" , "palmitate-_methyl_ester" , "palmitoleate_161n7" , "palmitoylcarnitine" , "palmitoyl_sphingomyelin" , "pantothenate" , "p-cresol_sulfate" , "pelargonate_90" , "pentadecanoate_150" , "phenol_sulfate" , "phenylacetylglutamine" , "phenylalanine" , "phenylalanylphenylalanine" , "phenylalanyltryptophan" , "phosphate" , "pipecolate" , "pregnanediol-3-glucuronide" , "pregnen-diol_disulfate" , "pregnenolone_sulfate" , "pregn_steroid_monosulfate" , "pro-hydroxy-pro" , "proline" , "propionylcarnitine" , "pseudouridine" , "pyroglutamine" , "pyruvate" , "ribitol" , "riboflavin_Vitamin_B2" , "ribose" , "scyllo-inositol" , "serine" , "serotonin_5HT" , "S-methylcysteine" , "stachydrine" , "stearate_180" , "stearidonate_184n3" , "stearoyl_sphingomyelin" , "succinate" , "taurocholenate_sulfate" , "taurolithocholate_3-sulfate" , "theobromine" , "threitol" , "threonate" , "threonine" , "trans-4-hydroxyproline" , "tryptophan" , "tyrosine" , "undecanoate_110" , "urate" , "urea" , "uridine" , "valine" , "xanthine" , "xylitol" , "xylonate" , "xylose") #1:130 or #151:280
#split in two because of some R limitations

data <- list.files(path = "all_metabos_delta/noLog/first_part", full.names = TRUE, pattern="*.csv")
data2 <- list.files(path = "all_metabos_delta/noLog/second_part", full.names = TRUE, pattern="*.csv")


table_time = NULL;
for (i in 1:150){
    f <- read.table(file=data[i], sep="\t", header=T)
    f$treatment <- factor(f$treatment)
    glm_drug <- glm(treatment ~ value, family=binomial(link='logit'),data=f)
    summary_drug <- summary(glm_drug)
    line = c(metabos[i],summary_drug$coefficients[2], summary_drug$coefficients[8])
    #metabo_name, FC, pval
    table_time = rbind(table_time,line)
}

for (i in 1:130){
    f <- read.table(file=data2[i], sep="\t", header=T)
    f$treatment <- factor(f$treatment)
    glm_drug <- glm(treatment ~ value, family=binomial(link='logit'),data=f)
    summary_drug <- summary(glm_drug)
    line2 = c(metabos2[i],summary_drug$coefficients[2], summary_drug$coefficients[8])
    #metabo_name, FC, pval
    table_time = rbind(table_time,line2)
}
#Warning messages:
#1: glm.fit: fitted probabilities numerically 0 or 1 occurred 
#2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
#3: glm.fit: fitted probabilities numerically 0 or 1 occurred 

table_time <- as.data.frame(table_time)
table_time <- plyr::rename(table_time, c("V1"="Metabolites", "V2"="FC", "V3"="pval"))
rownames(table_time) <- table_time$Metabolites
table_time$Metabolites <- NULL

table_time$FC<-as.numeric(levels(table_time$FC))[table_time$FC]
table_time$pval<-as.numeric(levels(table_time$pval))[table_time$pval]

table_time$fdr <-p.adjust(table_time$pval,"fdr")


write.table(table_time, file="all_metabos_delta/noLog/output_logit.csv", sep="\t")

