library(MAGeCKFlute)
library(ggplot2)
library(dplyr)
library(tidyverse)

#setwd("/blue/ferrallm/sneupane/mageck/arielcrispr3/workflow2020/results/")


#data("rra.gene_summary")
#data("rra.sgrna_summary")
#file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),"count/all.countsummary.txt")


countsummary = read.delim("count/all.countsummary.txt", check.names = FALSE)
countsummary$Label <- sub("_R1_001$", "", countsummary$Label)
head(countsummary)
colnames(countsummary)
MapRatesView(countsummary)

IdentBarView(countsummary, x = "Label", y = "GiniIndex", 
             ylab = "Gini index", main = "Evenness of sgRNA reads")



countsummary$Missed = log10(countsummary$Zerocounts)
IdentBarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
             ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")




BarView(countsummary, x = "Label", y = "GiniIndex", ylab = "Gini index", main = "Evenness of sgRNA reads")
countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")



library(ggplot2)
library(reshape2)

countsummary$Unmapped <- countsummary$Reads - countsummary$Mapped
gg <- melt(countsummary[, c("Label", "Mapped", "Unmapped")], id.vars = "Label")
gg$variable <- factor(gg$variable, levels = c("Unmapped", "Mapped"))

p <- ggplot(gg, aes(x = Label, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = NULL, y = "Reads", fill = NULL) +
  ggtitle("Map ratio") +
  scale_fill_manual(values = c("#9BC7E9", "#1C6DAB")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(p)




countsummary$Unmapped = countsummary$Reads - countsummary$Mapped
gg = reshape2::melt(countsummary[, c("Label", "Mapped", "Unmapped")], id.vars = "Label")
gg$variable = factor(gg$variable, levels = c("Unmapped", "Mapped"))
p = BarView(gg, x = "Label", y = "value", fill = "variable", 
            position = "stack", xlab = NULL, ylab = "Reads", main = "Map ratio")
p + scale_fill_manual(values = c("#9BC7E9", "#1C6DAB"))


gdata = ReadRRA("test/mle.gene_summary.txt")
genesummary = read.delim("test/mle.gene_summary.txt")

head(genesummary)
head(gdata)

p1 = VolcanoView(gdata, x = "Score", y = "FDR", Label = "id")
print(p1)

gdatabeta = ReadBeta("test/mle.gene_summary.txt")
head(gdatabeta)
tail(gdatabeta)


gdatabeta_DMSO <- gdatabeta %>% select(Gene,DMSO_day7, DMSO_day14, DMSO_day21)
ctrlnameD = "DMSO_day7"
treatnameD1 = "DMSO_day14"
treatnameD2 = "DMSO_day21"
ViolinView(gdatabeta_DMSO, samples = c(ctrlnameD, treatnameD1, treatnameD2))

gdatabeta_AZA <- gdatabeta %>% select(Gene,AZA_day7, AZA_day14, AZA_day21)
head(gdatabeta_AZA)
ctrlnameA = "AZA_day7"
treatnameA1 = "AZA_day14"
treatnameA2 = "AZA_day21"
ViolinView(gdatabeta_AZA, samples = c(ctrlnameA, treatnameA1, treatnameA2))

gdatabeta_Combo <- gdatabeta %>% select(Gene,Combo_day7, Combo_day14, Combo_day21)
ctrlnameC = "Combo_day7"
treatnameC1 = "Combo_day14"
treatnameC2 = "Combo_day21"
ViolinView(gdatabeta_Combo, samples = c(ctrlnameC, treatnameC1, treatnameC2))


colnames(genesummary)
##########################################################################
gdatabeta_7DC7 <- gdatabeta %>% select(Gene,DMSO_day7,Combo_day7)
head(gdatabeta_7DC7)
ctrlname = "DMSO_day7"
treatname = "Combo_day7"
gdata_7DC7 = NormalizeBeta(gdatabeta_7DC7, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_7DC7)
ViolinView(gdata_7DC7, samples=c(ctrlname, treatname))
DensityView(gdata_7DC7, samples=c(ctrlname, treatname))
DensityDiffView(gdata_7DC7, ctrlname, treatname)
ConsistencyView(gdata_7DC7, ctrlname, treatname)
MAView(gdata_7DC7, ctrlname, treatname)
#CellCycleView(gdata_7DC7, ctrlname, treatname)
gdata_7DC7$Control = rowMeans(gdata_7DC7[,ctrlname, drop = FALSE])
gdata_7DC7$Treatment = rowMeans(gdata_7DC7[,treatname, drop = FALSE])
p4 = ScatterView(gdata_7DC7, ctrlname, treatname)
print(p4)



gdata_7DC7$Diff = gdata_7DC7$Treatment - gdata_7DC7$Control
gdata_7DC7$Rank = rank(gdata_7DC7$Diff)
p5 = ScatterView(gdata_7DC7, x = "Diff", y = "Rank", label = "Gene", 
                 top = 20, model = "rank", max.overlaps = 50)
print(p5)

head(gdata_7DC7)




rankdata = gdata_7DC7$Treatment - gdata_7DC7$Control
names(rankdata) = gdata_7DC7$Gene
RankView(rankdata,
         top = 20, max.overlaps = 50)

rankdata = gdata_7DC7$Treatment - gdata_7DC7$Control
names(rankdata) = gdata_7DC7$Gene
RankView(rankdata)

p6 = ScatterView(gdata_7DC7, x = "DMSO_day7", y = "Combo_day7", label = "Gene", 
                 model = "ninesquare", top = 15, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p6)


p7 = SquareView(gdata_7DC7, label = "Gene", 
                x_cutoff = CutoffCalling(gdata_7DC7$Control, 1), 
                y_cutoff = CutoffCalling(gdata_7DC7$Treatment, 1), model = "ninesquare", top = 15, display_cut = TRUE, force = 1, max.overlaps = 50)
print(p7)

p1 = ScatterView(gdata_7DC7, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)


#install.packages("pathview")
library(pathview)
pathview::pathview(genedata, pathway.id = "hsa04115", species = "hsa", kegg.file = "hsa04115.xml")


# 9-square groups
Square9 = p7$data
idx=Square9$group=="topcenter"
geneList = Square9$Diff
head(geneList)
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene


Square8 = p1$data
idx=Square8$group=="top"
geneList1 = Square8$Diff
head(geneList1)
names(geneList1) = Square8$Gene[idx]
universe = Square8$Gene
kegg1 = EnrichAnalyzer(geneList = geneList1, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)






# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)

genedata = gdata_7DC7[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa03013", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa05166", organism = "hsa", sub = NULL)
hgtRes1 = EnrichAnalyzer(geneList[1:200], method = "HGT", 
                         type = "Pathway", organism = "hsa")
head(hgtRes1@result)

ortRes1 = EnrichAnalyzer(geneList[1:200], method = "ORT", 
                         type = "KEGG", organism = "hsa")
head(ortRes1@result)

gseRes1 = EnrichAnalyzer(geneList, method = "GSEA", type = "Pathway", organism = "hsa")

df = hgtRes1@result
df$logFDR = -log10(df$p.adjust)
p = BarView(df[1:10,], "Description", 'logFDR')
p = p + labs(x = NULL) + coord_flip()
p

# Or use function barplot from enrichplot package
barplot(hgtRes1, showCategory = 10)

EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 1)
EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 2)
dotplot(hgtRes1, showCategory = 10)

library(enrichplot)
hgtRes1@result$geneID = hgtRes1@result$geneName
cnetplot(hgtRes1, 10)
heatplot(hgtRes1, showCategory = 10, foldChange = geneList)
tmp <- pairwise_termsim(hgtRes1)
emapplot(tmp, showCategory = 5, layout = "kk")
# show GSEA results of one pathway
idx = which(gseRes1$NES>0)[1]
gseaplot(gseRes1, geneSetID = idx, title = gseRes1$Description[idx])
# show GSEA results of multiple pathways
gseaplot2(gseRes1, geneSetID = which(gseRes1$NES>0)[1:3])
enrichComb = EnrichAnalyzer(geneList[1:200], type = "KEGG+REACTOME", organism = "hsa")
EnrichedView(enrichComb, top = 20)
## All pathways
enrich = EnrichAnalyzer(geneList = geneList[1:200], type = "REACTOME", organism = "hsa")
EnrichedView(enrich, top = 20)
## All gene ontology
enrichGo = EnrichAnalyzer(geneList[1:200], type = "GOBP", organism = "hsa")
EnrichedView(enrichGo, top = 20)
enrichPro = EnrichAnalyzer(geneList[1:200], type = "Complex", organism = "hsa")
EnrichedView(enrichPro, top = 20)
enrich = EnrichAnalyzer(geneList[1:200], type = "GOBP", limit = c(2, 80), organism = "hsa")
EnrichedView(enrich, top = 20)

enrich1 = EnrichAnalyzer(geneList[1:200], type = "Pathway", organism = "hsa")
enrich2 = EnrichedFilter(enrich1)
EnrichedView(enrich1, top = 15)
EnrichedView(enrich2, top = 15)



# ## Add column of 'diff'
# gdata_7DC7$Control = rowMeans(gdata_7DC7[,ctrlname, drop = FALSE])
# gdata_7DC7$Treatment = rowMeans(gdata_7DC7[,treatname, drop = FALSE])
# 
# rankdata = gdata_7DC7$Treatment - gdata_7DC7$Control
# names(rankdata) = gdata_7DC7$Gene
# p2 = RankView(rankdata)
# print(p2)





################################################################
gdatabeta_DA7 <- gdatabeta %>% select(Gene,DMSO_day7,AZA_day7)
head(gdatabeta_DA7)
ctrlname = "DMSO_day7"
treatname = "AZA_day7"
gdata_DA7 = NormalizeBeta(gdatabeta_DA7, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_DA7)
ViolinView(gdata_DA7, samples=c(ctrlname, treatname))
DensityView(gdata_DA7, samples=c(ctrlname, treatname))
DensityDiffView(gdata_DA7, ctrlname, treatname)
ConsistencyView(gdata_DA7, ctrlname, treatname)
MAView(gdata_DA7, ctrlname, treatname)
#CellCycleView(gdata_7DC7, ctrlname, treatname)
gdata_DA7$Control = rowMeans(gdata_DA7[,ctrlname, drop = FALSE])
gdata_DA7$Treatment = rowMeans(gdata_DA7[,treatname, drop = FALSE])
p4 = ScatterView(gdata_DA7, ctrlname, treatname)
print(p4)
p1 = ScatterView(gdata_DA7, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)

gdata_DA7$Diff = gdata_DA7$Treatment - gdata_DA7$Control
gdata_DA7$Rank = rank(gdata_DA7$Diff)
p5 = ScatterView(gdata_DA7, x = "Diff", y = "Rank", label = "Gene", 
                 top = 10, model = "rank", max.overlaps = 50)
print(p5)

rankdata = gdata_DA7$Treatment - gdata_DA7$Control
names(rankdata) = gdata_DA7$Gene
RankView(rankdata,
         top = 20, max.overlaps = 50)

rankdata = gdata_DA7$Treatment - gdata_DA7$Control
names(rankdata) = gdata_DA7$Gene
RankView(rankdata)

p6 = ScatterView(gdata_DA7, x = "DMSO_day7", y = "AZA_day7", label = "Gene", 
                 model = "ninesquare", top = 10, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p6)
p8 = SquareView(gdata_DA7, label = "Gene", 
                x_cutoff = CutoffCalling(gdata_DA7$Control, 2), 
                y_cutoff = CutoffCalling(gdata_DA7$Treatment, 2), model = "ninesquare", top = 15, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p8)

p1 = ScatterView(gdata_DA7, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)


#install.packages("pathview")
library(pathview)
pathview::pathview(genedata, pathway.id = "hsa04115", species = "hsa", kegg.file = "hsa04115.xml")


Square8 = p1$data
idx=Square8$group=="top"
geneList1 = Square8$Diff
head(geneList1)
names(geneList1) = Square8$Gene[idx]
universe = Square8$Gene
kegg1 = EnrichAnalyzer(geneList = geneList1, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)

# 9-square groups
Square9 = p8$data
idx=Square9$group=="topcenter"
geneList = Square9$Diff
head(geneList)
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene

# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 10, bottom = 0)

genedata = gdata_DA7[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa03013", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa05166", organism = "hsa", sub = NULL)


# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)

genedata = gdata_DA7[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)

hgtRes1 = EnrichAnalyzer(geneList[1:200], method = "HGT", 
                         type = "Pathway", organism = "hsa")
head(hgtRes1@result)

ortRes1 = EnrichAnalyzer(geneList[1:200], method = "ORT", 
                         type = "KEGG", organism = "hsa")
head(ortRes1@result)

gseRes1 = EnrichAnalyzer(geneList, method = "GSEA", type = "Pathway", organism = "hsa")

df = hgtRes1@result
df$logFDR = -log10(df$p.adjust)
p = BarView(df[1:10,], "Description", 'logFDR')
p = p + labs(x = NULL) + coord_flip()
p

# Or use function barplot from enrichplot package
barplot(hgtRes1, showCategory = 10)

EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 1)
EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 2)
dotplot(hgtRes1, showCategory = 10)

library(enrichplot)
hgtRes1@result$geneID = hgtRes1@result$geneName
cnetplot(hgtRes1, 10)
heatplot(hgtRes1, showCategory = 10, foldChange = geneList)
tmp <- pairwise_termsim(hgtRes1)
emapplot(tmp, showCategory = 5, layout = "kk")
# show GSEA results of one pathway
idx = which(gseRes1$NES>0)[1]
gseaplot(gseRes1, geneSetID = idx, title = gseRes1$Description[idx])
# show GSEA results of multiple pathways
gseaplot2(gseRes1, geneSetID = which(gseRes1$NES>0)[1:3])
enrichComb = EnrichAnalyzer(geneList[1:200], type = "KEGG+REACTOME", organism = "hsa")
EnrichedView(enrichComb, top = 20)
## All pathways
enrich = EnrichAnalyzer(geneList = geneList[1:200], type = "REACTOME", organism = "hsa")
EnrichedView(enrich, top = 20)
## All gene ontology
enrichGo = EnrichAnalyzer(geneList[1:200], type = "GOBP", organism = "hsa")
EnrichedView(enrichGo, top = 20)
enrichPro = EnrichAnalyzer(geneList[1:200], type = "Complex", organism = "hsa")
EnrichedView(enrichPro, top = 20)
enrich = EnrichAnalyzer(geneList[1:200], type = "GOBP", limit = c(2, 80), organism = "hsa")
EnrichedView(enrich, top = 20)

enrich1 = EnrichAnalyzer(geneList[1:200], type = "Pathway", organism = "hsa")
enrich2 = EnrichedFilter(enrich1)
EnrichedView(enrich1, top = 15)
EnrichedView(enrich2, top = 15)





########################################################

gdatabeta_DA14 <- gdatabeta %>% select(Gene,DMSO_day14,AZA_day14)
head(gdatabeta_DA14)
ctrlname = "DMSO_day14"
treatname = "AZA_day14"
gdata_DA14 = NormalizeBeta(gdatabeta_DA14, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_DA14)
ViolinView(gdata_DA14, samples=c(ctrlname, treatname))
DensityView(gdata_DA14, samples=c(ctrlname, treatname))
DensityDiffView(gdata_DA14, ctrlname, treatname)
ConsistencyView(gdata_DA14, ctrlname, treatname)
MAView(gdata_DA14, ctrlname, treatname)
#CellCycleView(gdata_7DC7, ctrlname, treatname)
gdata_DA14$Control = rowMeans(gdata_DA14[,ctrlname, drop = FALSE])
gdata_DA14$Treatment = rowMeans(gdata_DA14[,treatname, drop = FALSE])
p4 = ScatterView(gdata_DA14, ctrlname, treatname)
print(p4)
p1 = ScatterView(gdata_DA14, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)

gdata_DA14$Diff = gdata_DA14$Treatment - gdata_DA14$Control
gdata_DA14$Rank = rank(gdata_DA14$Diff)
p5 = ScatterView(gdata_DA14, x = "Diff", y = "Rank", label = "Gene", 
                 top = 10, model = "rank", max.overlaps = 50)
print(p5)

rankdata = gdata_DA14$Treatment - gdata_DA14$Control
names(rankdata) = gdata_DA14$Gene
RankView(rankdata,
         top = 10, max.overlaps = 50)

rankdata = gdata_DA14$Treatment - gdata_DA14$Control
names(rankdata) = gdata_DA14$Gene
RankView(rankdata)

p6 = ScatterView(gdata_DA14, x = "DMSO_day14", y = "AZA_day14", label = "Gene", 
                 model = "ninesquare", top = 10, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p6)

p8 = SquareView(gdata_DA14, label = "Gene", 
                x_cutoff = CutoffCalling(gdata_DA14$Control, 2), 
                y_cutoff = CutoffCalling(gdata_DA14$Treatment, 2), model = "ninesquare", top = 15, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p8)


p1 = ScatterView(gdata_DA14, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)


#install.packages("pathview")
library(pathview)
pathview::pathview(genedata, pathway.id = "hsa04115", species = "hsa", kegg.file = "hsa04115.xml")


Square8 = p1$data
idx=Square8$group=="top"
geneList1 = Square8$Diff
head(geneList1)
names(geneList1) = Square8$Gene[idx]
universe = Square8$Gene
kegg1 = EnrichAnalyzer(geneList = geneList1, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)

# 9-square groups
Square9 = p8$data
idx=Square9$group=="topcenter"
geneList = Square9$Diff
head(geneList)
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene


# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 10, bottom = 0)

genedata = gdata_DA14[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa03013", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa05166", organism = "hsa", sub = NULL)


# 9-square groups
Square9 = p8$data
idx=Square9$group=="midleft"
geneList = Square9$Diff
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene

# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 10, bottom = 0)

genedata = gdata_DA14[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)


# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)

genedata = gdata_DA14[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)

hgtRes1 = EnrichAnalyzer(geneList[1:200], method = "HGT", 
                         type = "Pathway", organism = "hsa")
head(hgtRes1@result)

ortRes1 = EnrichAnalyzer(geneList[1:200], method = "ORT", 
                         type = "KEGG", organism = "hsa")
head(ortRes1@result)

gseRes1 = EnrichAnalyzer(geneList, method = "GSEA", type = "Pathway", organism = "hsa")

df = hgtRes1@result
df$logFDR = -log10(df$p.adjust)
p = BarView(df[1:10,], "Description", 'logFDR')
p = p + labs(x = NULL) + coord_flip()
p

# Or use function barplot from enrichplot package
barplot(hgtRes1, showCategory = 10)

EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 1)
EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 2)
dotplot(hgtRes1, showCategory = 10)

library(enrichplot)
hgtRes1@result$geneID = hgtRes1@result$geneName
cnetplot(hgtRes1, 10)
heatplot(hgtRes1, showCategory = 10, foldChange = geneList)
tmp <- pairwise_termsim(hgtRes1)
emapplot(tmp, showCategory = 5, layout = "kk")
# show GSEA results of one pathway
idx = which(gseRes1$NES>0)[1]
gseaplot(gseRes1, geneSetID = idx, title = gseRes1$Description[idx])
# show GSEA results of multiple pathways
gseaplot2(gseRes1, geneSetID = which(gseRes1$NES>0)[1:3])
enrichComb = EnrichAnalyzer(geneList[1:200], type = "KEGG+REACTOME", organism = "hsa")
EnrichedView(enrichComb, top = 20)
## All pathways
enrich = EnrichAnalyzer(geneList = geneList[1:200], type = "REACTOME", organism = "hsa")
EnrichedView(enrich, top = 20)
## All gene ontology
enrichGo = EnrichAnalyzer(geneList[1:200], type = "GOBP", organism = "hsa")
EnrichedView(enrichGo, top = 20)
enrichPro = EnrichAnalyzer(geneList[1:200], type = "Complex", organism = "hsa")
EnrichedView(enrichPro, top = 20)
enrich = EnrichAnalyzer(geneList[1:200], type = "GOBP", limit = c(2, 80), organism = "hsa")
EnrichedView(enrich, top = 20)

enrich1 = EnrichAnalyzer(geneList[1:200], type = "Pathway", organism = "hsa")
enrich2 = EnrichedFilter(enrich1)
EnrichedView(enrich1, top = 15)
EnrichedView(enrich2, top = 15)
###################################################################################



gdatabeta_DC14 <- gdatabeta %>% select(Gene,DMSO_day14,Combo_day14)
head(gdatabeta_DC14)
ctrlname = "DMSO_day14"
treatname = "Combo_day14"
gdata_DC14 = NormalizeBeta(gdatabeta_DC14, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_DC14)
ViolinView(gdata_DC14, samples=c(ctrlname, treatname))
DensityView(gdata_DC14, samples=c(ctrlname, treatname))
DensityDiffView(gdata_DC14, ctrlname, treatname)
ConsistencyView(gdata_DC14, ctrlname, treatname)
MAView(gdata_DC14, ctrlname, treatname)
#CellCycleView(gdata_7DC7, ctrlname, treatname)
gdata_DC14$Control = rowMeans(gdata_DC14[,ctrlname, drop = FALSE])
gdata_DC14$Treatment = rowMeans(gdata_DC14[,treatname, drop = FALSE])
p4 = ScatterView(gdata_DC14, ctrlname, treatname)
print(p4)
p1 = ScatterView(gdata_DC14, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("P53", "NF2", "EP300"))
print(p1)

gdata_DC14$Diff = gdata_DC14$Treatment - gdata_DC14$Control
gdata_DC14$Rank = rank(gdata_DC14$Diff)
p5 = ScatterView(gdata_DC14, x = "Diff", y = "Rank", label = "Gene", 
                 top = 10, model = "rank", max.overlaps = 50)
print(p5)

rankdata = gdata_DC14$Treatment - gdata_DC14$Control
names(rankdata) = gdata_DC14$Gene
RankView(rankdata,
         top = 10, max.overlaps = 50)

rankdata = gdata_DC14$Treatment - gdata_DC14$Control
names(rankdata) = gdata_DC14$Gene
RankView(rankdata)

p6 = ScatterView(gdata_DC14, x = "DMSO_day14", y = "Combo_day14", label = "Gene", 
                 model = "ninesquare", top = 10, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p6)

p8 = SquareView(gdata_DC14, label = "Gene", 
                x_cutoff = CutoffCalling(gdata_DC14$Control, 2), 
                y_cutoff = CutoffCalling(gdata_DC14$Treatment, 2), model = "ninesquare", top = 15, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p8)

p1 = ScatterView(gdata_DC14, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)


#install.packages("pathview")
library(pathview)
pathview::pathview(genedata, pathway.id = "hsa04115", species = "hsa", kegg.file = "hsa04115.xml")


Square8 = p1$data
idx=Square8$group=="top"
geneList1 = Square8$Diff
head(geneList1)
names(geneList1) = Square8$Gene[idx]
universe = Square8$Gene
kegg1 = EnrichAnalyzer(geneList = geneList1, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)

# 9-square groups
Square9 = p8$data
idx=Square9$group=="topcenter"
geneList = Square9$Diff
head(geneList)
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene


# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 10, bottom = 0)

genedata = gdata_DC14[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa03013", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa05166", organism = "hsa", sub = NULL)

hgtRes1 = EnrichAnalyzer(geneList[1:200], method = "HGT", 
                         type = "Pathway", organism = "hsa")
head(hgtRes1@result)

ortRes1 = EnrichAnalyzer(geneList[1:200], method = "ORT", 
                         type = "KEGG", organism = "hsa")
head(ortRes1@result)

gseRes1 = EnrichAnalyzer(geneList, method = "GSEA", type = "Pathway", organism = "hsa")

df = hgtRes1@result
df$logFDR = -log10(df$p.adjust)
p = BarView(df[1:10,], "Description", 'logFDR')
p = p + labs(x = NULL) + coord_flip()
p

# Or use function barplot from enrichplot package
barplot(hgtRes1, showCategory = 10)

EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 1)
EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 2)
dotplot(hgtRes1, showCategory = 10)

library(enrichplot)
hgtRes1@result$geneID = hgtRes1@result$geneName
cnetplot(hgtRes1, 10)
heatplot(hgtRes1, showCategory = 10, foldChange = geneList)
tmp <- pairwise_termsim(hgtRes1)
emapplot(tmp, showCategory = 5, layout = "kk")
# show GSEA results of one pathway
idx = which(gseRes1$NES>0)[1]
gseaplot(gseRes1, geneSetID = idx, title = gseRes1$Description[idx])
# show GSEA results of multiple pathways
gseaplot2(gseRes1, geneSetID = which(gseRes1$NES>0)[1:3])
enrichComb = EnrichAnalyzer(geneList[1:200], type = "KEGG+REACTOME", organism = "hsa")
EnrichedView(enrichComb, top = 20)
## All pathways
enrich = EnrichAnalyzer(geneList = geneList[1:200], type = "REACTOME", organism = "hsa")
EnrichedView(enrich, top = 20)
## All gene ontology
enrichGo = EnrichAnalyzer(geneList[1:200], type = "GOBP", organism = "hsa")
EnrichedView(enrichGo, top = 20)
enrichPro = EnrichAnalyzer(geneList[1:200], type = "Complex", organism = "hsa")
EnrichedView(enrichPro, top = 20)
enrich = EnrichAnalyzer(geneList[1:200], type = "GOBP", limit = c(2, 80), organism = "hsa")
EnrichedView(enrich, top = 20)

enrich1 = EnrichAnalyzer(geneList[1:200], type = "Pathway", organism = "hsa")
enrich2 = EnrichedFilter(enrich1)
EnrichedView(enrich1, top = 15)
EnrichedView(enrich2, top = 15)

####################################################

gdatabeta_DA21 <- gdatabeta %>% select(Gene,DMSO_day21,AZA_day21)
head(gdatabeta_DA21)
ctrlname = "DMSO_day21"
treatname = "AZA_day21"
gdata_DA21 = NormalizeBeta(gdatabeta_DA21, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_DA21)
ViolinView(gdata_DA21, samples=c(ctrlname, treatname))
DensityView(gdata_DA21, samples=c(ctrlname, treatname))
DensityDiffView(gdata_DA21, ctrlname, treatname)
ConsistencyView(gdata_DA21, ctrlname, treatname)
MAView(gdata_DA21, ctrlname, treatname)
#CellCycleView(gdata_7DC7, ctrlname, treatname)
gdata_DA21$Control = rowMeans(gdata_DA21[,ctrlname, drop = FALSE])
gdata_DA21$Treatment = rowMeans(gdata_DA21[,treatname, drop = FALSE])
p4 = ScatterView(gdata_DA21, ctrlname, treatname)
print(p4)
p1 = ScatterView(gdata_DA21, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)

# Calculate the difference between Control and Treatment
gdata_DA21$Diff <- gdata_DA21$Treatment - gdata_DA21$Control

# Select the top genes based on the Diff column
top_genes <- 6
top_data <- gdata_DA21[order(-gdata_DA21$Diff), ][1:top_genes, ]

# Create the ScatterView plot
p1 <- ScatterView(gdata_DA21, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE)

# Add labels to the top genes
p1 <- p1 + geom_text(data = top_data, aes(label = Gene), color = "red", size = 2)

print(p1)


# Calculate the difference between Control and Treatment
gdata_DA21$Diff <- gdata_DA21$Treatment - gdata_DA21$Control

# Calculate the difference between Control and Treatment
gdata_DA21$Diff <- gdata_DA21$Treatment - gdata_DA21$Control

# Create the ScatterView plot for all genes
library(ggplot2)
library(ggdist)

# Create the ScatterView plot for all genes
p1 <- ScatterView(gdata_DA21, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE)

# Add labels to the XPO1 gene only
xpo1_data <- gdata_DA21[gdata_DA21$Gene == "XPO1", ]
p1 <- p1 + geom_text(data = xpo1_data, aes(label = Gene), color = "red", size = 8)

print(p1)



gdata_DA21$Diff = gdata_DA21$Treatment - gdata_DA21$Control
gdata_DA21$Rank = rank(gdata_DA21$Diff)
p5 = ScatterView(gdata_DA21, x = "Diff", y = "Rank", label = "Gene", 
                 top = 10, model = "rank", max.overlaps = 50)
print(p5)

rankdata = gdata_DA21$Treatment - gdata_DA21$Control
names(rankdata) = gdata_DA21$Gene
RankView(rankdata,
         top = 10, max.overlaps = 50)

rankdata = gdata_DA21$Treatment - gdata_DA21$Control
names(rankdata) = gdata_DA21$Gene
RankView(rankdata)

p6 = ScatterView(gdata_DA21, x = "DMSO_day21", y = "AZA_day21", label = "Gene", 
                 model = "ninesquare", top = 10, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p6)

p8 = SquareView(gdata_DA21, label = "Gene", 
               x_cutoff = CutoffCalling(gdata_DA21$Control, 2), 
               y_cutoff = CutoffCalling(gdata_DA21$Treatment, 2), model = "ninesquare", top = 15, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p8)


head(gdata_DA21)

write.csv(gdata_DA21, "gdata_DA21.csv")
gdata_DC21_100rem <- read.csv("gdata_DC21_top100AZAremoved.csv", header = TRUE, sep = ",")





#install.packages("pathview")
library(pathview)
pathview::pathview(genedata, pathway.id = "hsa04115", species = "hsa", kegg.file = "hsa04115.xml")


Square8 = p1$data
idx=Square8$group=="top"
geneList1 = Square8$Diff
head(geneList1)
names(geneList1) = Square8$Gene[idx]
universe = Square8$Gene
kegg1 = EnrichAnalyzer(geneList = geneList1, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)

# 9-square groups
Square9 = p8$data
idx=Square9$group=="topcenter"
geneList = Square9$Diff
head(geneList)
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene


# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 10, bottom = 0)

genedata = gdata_DA21[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa03013", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa05166", organism = "hsa", sub = NULL)

hgtRes1 = EnrichAnalyzer(geneList[1:200], method = "HGT", 
                         type = "Pathway", organism = "hsa")
head(hgtRes1@result)

ortRes1 = EnrichAnalyzer(geneList[1:200], method = "ORT", 
                         type = "KEGG", organism = "hsa")
head(ortRes1@result)

gseRes1 = EnrichAnalyzer(geneList, method = "GSEA", type = "Pathway", organism = "hsa")

df = hgtRes1@result
df$logFDR = -log10(df$p.adjust)
p = BarView(df[1:10,], "Description", 'logFDR')
p = p + labs(x = NULL) + coord_flip()
p

# Or use function barplot from enrichplot package
barplot(hgtRes1, showCategory = 10)

EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 1)
EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 2)
dotplot(hgtRes1, showCategory = 10)

library(enrichplot)
hgtRes1@result$geneID = hgtRes1@result$geneName
cnetplot(hgtRes1, 10)
heatplot(hgtRes1, showCategory = 10, foldChange = geneList)
tmp <- pairwise_termsim(hgtRes1)
emapplot(tmp, showCategory = 5, layout = "kk")
# show GSEA results of one pathway
idx = which(gseRes1$NES>0)[1]
gseaplot(gseRes1, geneSetID = idx, title = gseRes1$Description[idx])
# show GSEA results of multiple pathways
gseaplot2(gseRes1, geneSetID = which(gseRes1$NES>0)[1:3])
enrichComb = EnrichAnalyzer(geneList[1:200], type = "KEGG+REACTOME", organism = "hsa")
EnrichedView(enrichComb, top = 20)
## All pathways
enrich = EnrichAnalyzer(geneList = geneList[1:200], type = "REACTOME", organism = "hsa")
EnrichedView(enrich, top = 20)
## All gene ontology
enrichGo = EnrichAnalyzer(geneList[1:200], type = "GOBP", organism = "hsa")
EnrichedView(enrichGo, top = 20)
enrichPro = EnrichAnalyzer(geneList[1:200], type = "Complex", organism = "hsa")
EnrichedView(enrichPro, top = 20)
enrich = EnrichAnalyzer(geneList[1:200], type = "GOBP", limit = c(2, 80), organism = "hsa")
EnrichedView(enrich, top = 20)

enrich1 = EnrichAnalyzer(geneList[1:200], type = "Pathway", organism = "hsa")
enrich2 = EnrichedFilter(enrich1)
EnrichedView(enrich1, top = 15)
EnrichedView(enrich2, top = 15)
##################################################################
gdatabeta = ReadBeta("test/mle.gene_summary.txt")
head(gdatabeta)
tail(gdatabeta)

gdatabeta_DC21 <- gdatabeta %>% select(Gene,DMSO_day21,Combo_day21)
gdata_DC21_100rem <- read.csv("gdata_DC21_top100AZAremoved.csv", header = TRUE, sep = ",")
head(gdatabeta_DC21)
ctrlname = "DMSO_day21"
treatname = "Combo_day21"
gdata_DC21 = NormalizeBeta(gdatabeta_DC21, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_DC21)
ViolinView(gdata_DC21, samples=c(ctrlname, treatname))
DensityView(gdata_DC21, samples=c(ctrlname, treatname))
DensityDiffView(gdata_DC21, ctrlname, treatname)
ConsistencyView(gdata_DC21, ctrlname, treatname)
MAView(gdata_DC21, ctrlname, treatname)
#CellCycleView(gdata_7DC7, ctrlname, treatname)
gdata_DC21$Control = rowMeans(gdata_DC21[,ctrlname, drop = FALSE])
gdata_DC21$Treatment = rowMeans(gdata_DC21[,treatname, drop = FALSE])
p4 = ScatterView(gdata_DC21, ctrlname, treatname)
print(p4)
p1 = ScatterView(gdata_DC21, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("NF1", "NF2", "EP300"))
print(p1)

gdata_DC21$Diff = gdata_DC21$Treatment - gdata_DC21$Control
gdata_DC21$Rank = rank(gdata_DC21$Diff)
p5 = ScatterView(gdata_DC21, x = "Diff", y = "Rank", label = "Gene", 
                 top = 15, model = "rank", max.overlaps = 100)
print(p5)

# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", top = 15, model = "rank", max.overlaps = 100)

p5

# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", top = 5, model = "rank", max.overlaps = 100)

# Add label for the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(x = Diff, y = Rank, label = Gene), color = "red", size = 5)

print(p5)








# Calculate the difference between Control and Treatment
gdata_DC21$Diff = gdata_DC21$Treatment - gdata_DC21$Control

# Calculate ranks for the differences in the "Diff" column and store it in the "Rank" column of the gdata_DC21 dataframe using the rank function
gdata_DC21$Rank = rank(gdata_DC21$Diff)

# Create the ScatterView plot for all genes
library(ggplot2)
library(ggdist)

# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50)

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "red", size = 5)

# Add a dashed line at y = 0 to represent the cutoff
p5 <- p5 + geom_vline(xintercept = 0, linetype = "dashed", color = "blue")


print(p5)

# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50)

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "red", size = 5)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "blue")



# Add a square box with a red boundary color after negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 2)

print(p5)




# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50)

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "red", size = 3, nudge_x = 0.1, nudge_y = 0.1)

# Add a red dot for the "XPO1" gene
p5 <- p5 + geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 1)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "blue")

# Add a square box with a red boundary color after the negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 2)

print(p5)



# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", top = 15, model = "rank", max.overlaps = 50)

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "blue", size = 4, nudge_x = 0.5, nudge_y = 5)

# Add a red dot for the "XPO1" gene
p5 <- p5 + geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 1)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "blue")

# Add a square box with a red boundary color after the negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1)

print(p5)

# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", top = 5, model = "rank", max.overlaps = 50)

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "blue", size = 4, nudge_x = 0.4, nudge_y = 0)

# Add a red dot for the "XPO1" gene
p5 <- p5 + geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 1)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", size= 1.5, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box with a red boundary color after the negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5)

# Add "Negative Selection" text inside the square box
p5 <- p5 + annotate("text", x = -3.4, y = 15000, label = "Negative selection", color = "black", size = 4)

print(p5)

##########
# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50)

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold", nudge_x = 0.2, nudge_y = 500)

# Add a red dot for the "XPO1" gene
p5 <- p5 + geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 2)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box with a red boundary color after the negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5)

# Add "Negative Selection" text inside the square box
p5 <- p5 + annotate("text", x = -3.4, y = 15000, label = "Negative selection", color = "black", size = 4)

print(p5)

#######################
library(mageckflute)
library(ggplot2)

# Assuming gdata_DC21_100rem is already loaded with the gene data

# Create the ScatterView plot for all genes with a square aspect ratio
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1)  # Set aspect ratio to 1 for a square plot

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold", nudge_x = 0.2, nudge_y = 500)

# Add a red dot for the "XPO1" gene
p5 <- p5 + geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box with a red boundary color after the negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5)

# Add "Negative Selection" text inside the square box
p5 <- p5 + annotate("text", x = -3.2, y = 15000, label = "Negative selection", color = "black", size = 4.5, fontface = "bold")

# Show the plot with a square aspect ratio
print(p5)

ggsave("scatterview_plot.png", plot = p5, width = 10, height = 10, units = "in")


##############
library(mageckflute)
library(ggplot2)

# Assuming gdata_DC21_100rem is already loaded with the gene data

# Create the ScatterView plot for all genes with a square aspect ratio
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", top=5, model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1,  # Set aspect ratio to 1 for a square plot
        axis.title = element_text(size = 16), # Increase size of axis labels
        axis.text = element_text(size = 16))  # Increase size of axis tick labels

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold", nudge_x = 0.2, nudge_y = 500)

# Add a red dot for the "XPO1" gene
p5 <- p5 + geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box with a red boundary color after the negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5)

# Add "Negative Selection" text inside the square box
p5 <- p5 + annotate("text", x = -3.2, y = 15000, label = "Negative selection", color = "black", size = 4.5, fontface = "bold")

# Show the plot with a square aspect ratio
print(p5)



# Assuming gdata_DC21_100rem is your data frame and it contains the columns 'Gene', 'Diff', and 'Rank'

# Your existing plot code remains the same
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", top=5, model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1, axis.title = element_text(size = 16), axis.text = element_text(size = 16))

# Add label and red dot for the gene "XPO1"
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + 
  geom_text(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold", nudge_x = 0.2, nudge_y = 500) +
  geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Add dashed lines
p5 <- p5 + 
  geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add square box and "Negative Selection" text
p5 <- p5 + 
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5) +
  annotate("text", x = -3.2, y = 15000, label = "Negative selection", color = "black", size = 4.5, fontface = "bold")

# Identify top genes excluding XPO1
top_genes_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene %in% head(sort(gdata_DC21_100rem$Rank, decreasing = TRUE), 5) & gdata_DC21_100rem$Gene != "XPO1", ]

# Modify fonts and colors for other top genes
p5 <- p5 + geom_text(data = top_genes_data, aes(label = Gene), color = "desired_color", size = 5, fontface = "italic")

# Show the plot
print(p5)



# Load necessary library
library(ggrepel)

# Assuming gdata_DC21_100rem is your data frame and it contains the columns 'Gene', 'Diff', and 'Rank'

# Create the initial scatter plot
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1, axis.title = element_text(size = 16), axis.text = element_text(size = 16))

# Customize the appearance for the gene "XPO1"
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + 
  geom_text_repel(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold",  nudge_x = 0.2, nudge_y = 1000) +
  geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Customize the appearance for other specific genes
specific_genes <- c('MSN', 'QTRTD1', 'HSH2D', 'CTBP1', 'LILRB3', 'UCK2', 'MAGT1', 'ACSL4')
specific_genes_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene %in% specific_genes, ]
p5 <- p5 + 
  geom_text_repel(data = specific_genes_data, aes(label = Gene), color = "red", size = 4, fontface = "bold")

# Add dashed lines at specific x-intercepts
p5 <- p5 + 
  geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box and annotate with "Negative Selection"
p5 <- p5 + 
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5) +
  annotate("text", x = -3.2, y = 15000, label = "Negative selection", color = "black", size = 4.5, fontface = "bold")

# Show the final plot
print(p5)



# Load necessary library
library(ggrepel)

# Assuming gdata_DC21_100rem is your data frame and it contains the columns 'Gene', 'Diff', and 'Rank'

# Create the initial scatter plot
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 16, face = "bold"),  # Bold axis titles
        axis.text = element_text(size = 16, face = "bold"))   # Bold axis text

# Customize the appearance for the gene "XPO1"
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + 
  geom_text_repel(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold", nudge_x = 0.2, nudge_y = 1000) +
  geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Customize the appearance for other specific genes
specific_genes <- c('MSN', 'QTRTD1', 'HSH2D', 'CTBP1', 'LILRB3', 'UCK2', 'MAGT1', 'ACSL4')
specific_genes_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene %in% specific_genes, ]
p5 <- p5 + 
  geom_text_repel(data = specific_genes_data, aes(label = Gene), color = "red", size = 4, fontface = "bold")

# Add dashed lines at specific x-intercepts
p5 <- p5 + 
  geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box and annotate with vertical "Negative Selection" text
p5 <- p5 + 
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5) +
  annotate("text", x = -3.2, y = 15000, label = "Negative selection", angle = 90, color = "black", size = 5, fontface = "bold")

# Show the final plot
print(p5)


# Load necessary library
library(ggrepel)

# Assuming gdata_DC21_100rem is your data frame and it contains the columns 'Gene', 'Diff', and 'Rank'

# Create the initial scatter plot
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 16, face = "bold"),  # Bold axis titles
        axis.text = element_text(size = 16, face = "bold"))   # Bold axis text

# Customize the appearance for the gene "XPO1"
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + 
  geom_text_repel(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold", nudge_x = 0.2, nudge_y = 1000) +
  geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Customize the appearance for other specific genes, including adjusting MSN position
specific_genes <- c('MSN', 'QTRTD1', 'HSH2D', 'CTBP1', 'LILRB3', 'UCK2', 'MAGT1', 'ACSL4')
specific_genes_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene %in% specific_genes, ]
p5 <- p5 + 
  geom_text_repel(data = specific_genes_data, aes(label = Gene), color = "red", size = 4, fontface = "bold", 
                  nudge_y = ifelse(gdata_DC21_100rem$Gene == "MSN", -1000, -500)) # Adjust position of MSN

# Add dashed lines at specific x-intercepts
p5 <- p5 + 
  geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box and annotate with a larger, vertically centered "Negative Selection" text
average_rank <- mean(range(gdata_DC21_100rem$Rank))
p5 <- p5 + 
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5) +
  annotate("text", x = -4, y = average_rank, label = "Negative selection", angle = 90, color = "black", size = 6, fontface = "bold")

# Show the final plot
print(p5)


################################################
# Load necessary library
library(ggrepel)

# Assuming gdata_DC21_100rem is your data frame and it contains the columns 'Gene', 'Diff', and 'Rank'

# Create the initial scatter plot
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 18, face = "bold"),  # Bold axis titles
        axis.text = element_text(size = 18, face = "bold"))   # Bold axis text

# Customize the appearance for the gene "XPO1"
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + 
  geom_text_repel(data = xpo1_data, aes(label = Gene), color = "blue", size = 7, fontface = "bold", nudge_x = 0.2, nudge_y = 1000) +
  geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Customize the appearance for other specific genes, including adjusting MSN and QTRTD1 position
specific_genes <- c('MSN', 'QTRTD1', 'HSH2D', 'CTBP1', 'LILRB3', 'UCK2', 'MAGT1', 'ACSL4')
specific_genes_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene %in% specific_genes, ]
p5 <- p5 + 
  geom_text_repel(data = specific_genes_data, aes(label = Gene), color = "red", size = 5, fontface = "bold", 
                  nudge_y = ifelse(gdata_DC21_100rem$Gene == "MSN", -100, 
                                   ifelse(gdata_DC21_100rem$Gene == "QTRTD1", 500, -500))) # Adjust position of MSN and QTRTD1

# Add dashed lines at specific x-intercepts
p5 <- p5 + 
  geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box and annotate with a larger, more right-positioned "Negative Selection" text
average_rank <- mean(range(gdata_DC21_100rem$Rank))
p5 <- p5 + 
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5) +
  annotate("text", x = -3, y = average_rank, label = "Negative selection", angle = 90, color = "black", size = 6, fontface = "bold")

# Show the final plot
print(p5)


library(ggrepel)

# Assuming gdata_DC21_100rem is your data frame and it contains the columns 'Gene', 'Diff', and 'Rank'

# Create the initial scatter plot with Arial font
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1, 
        text = element_text(family = "Arial"), # Set Arial font for all text
        axis.title = element_text(size = 18, face = "bold", family = "Arial"),  # Bold axis titles in Arial
        axis.text = element_text(size = 18, face = "bold", family = "Arial"))   # Bold axis text in Arial

# Customize the appearance for the gene "XPO1" with Arial font
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + 
  geom_text_repel(data = xpo1_data, aes(label = Gene), family = "Arial", color = "blue", size = 10, fontface = "bold", nudge_x = 0.2, nudge_y = 1500) +
  geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Customize the appearance for other specific genes with Arial font, including adjusting MSN and QTRTD1 position
specific_genes <- c('MSN', 'QTRTD1', 'HSH2D', 'CTBP1', 'LILRB3', 'UCK2', 'MAGT1', 'ACSL4')
specific_genes_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene %in% specific_genes, ]
p5 <- p5 + 
  geom_text_repel(data = specific_genes_data, aes(label = Gene), family = "Arial", color = "red", size = 7, fontface = "bold", 
                  nudge_y = ifelse(gdata_DC21_100rem$Gene == "MSN", -1000, 
                                   ifelse(gdata_DC21_100rem$Gene == "QTRTD1", 500, -600))) # Adjust position of MSN and QTRTD1

# Add dashed lines at specific x-intercepts
p5 <- p5 + 
  geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box and annotate with a larger, more right-positioned "Negative Selection" text in Arial font
average_rank <- mean(range(gdata_DC21_100rem$Rank))
p5 <- p5 + 
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5) +
  annotate("text", x = -3, y = average_rank, label = "Negative selection", angle = 90, color = "black", size = 8, fontface = "bold", family = "Arial")

# Show the final plot
print(p5)


# Save the plot with a square aspect ratio and larger axis titles
ggsave("scatterview_plot.png", plot = p5, width = 10, height = 10, units = "in")


############################
library(mageckflute)
library(ggplot2)

# Assuming gdata_DC21_100rem is already loaded with the gene data

# Create the ScatterView plot for all genes with a square aspect ratio
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", top= 5, model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1,  # Set aspect ratio to 1 for a square plot
        axis.title = element_text(size = 16, face = "bold"), # Increase size of axis labels
        axis.text = element_text(size = 16, face = "bold"))  # Increase size of axis tick labels

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold", nudge_x = 0.2, nudge_y = 500)

# Add a red dot for the "XPO1" gene
p5 <- p5 + geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box with a red boundary color after the negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5)

# Add "Negative Selection" text inside the square box
p5 <- p5 + annotate("text", x = -3.2, y = 15000, label = "Negative selection", color = "black", size = 4.5, fontface = "bold")

# Show the plot with a square aspect ratio
print(p5)

# Save the plot with a square aspect ratio and larger axis titles
ggsave("scatterview_plot.png", plot = p5, width = 10, height = 10, units = "in")
ggsave("scatterview_plot.pdf", plot = p5, width = 10, height = 10, units = "in")

###################################
library(mageckflute)
library(ggplot2)

# Assuming gdata_DC21_100rem is already loaded with the gene data

# Create the ScatterView plot for all genes with a square aspect ratio
p5 <- ScatterView(gdata_DC21_100rem, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50) +
  theme(aspect.ratio = 1,  # Set aspect ratio to 1 for a square plot
        axis.title = element_text(size = 16, face = "bold"), # Increase size of axis labels
        axis.text = element_text(size = 16, face = "bold"))  # Increase size of axis tick labels

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21_100rem[gdata_DC21_100rem$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "blue", size = 5, fontface = "bold", nudge_x = 0.2, nudge_y = 500)

# Add a red dot for the "XPO1" gene
p5 <- p5 + geom_point(data = xpo1_data, aes(label = Gene), color = "red", size = 4)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", size= 1.0, color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", size=1.5, color = "blue")

# Add a square box with a red boundary color after the negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 1.5)

# Add "Negative Selection" text inside the square box
p5 <- p5 + annotate("text", x = -3.2, y = 15000, label = "Negative selection", color = "black", size = 4.5, fontface = "bold")

# Show the plot with a square aspect ratio
print(p5)

# Save the plot with a square aspect ratio and larger axis titles
ggsave("scatterview_plot.png", plot = p5, width = 10, height = 10, units = "in")
ggsave("scatterview_plot.pdf", plot = p5, width = 10, height = 10, units = "in")











#gdata_DC21_100rem <- read.csv("gdata_DC21_top100AZAremoved.csv", header = TRUE, sep = ",")









# Create the ScatterView plot for all genes
p5 <- ScatterView(gdata_DC21, x = "Diff", y = "Rank", label = "Gene", model = "rank", max.overlaps = 50)

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21[gdata_DC21$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "red", size = 8)

# Add dashed lines at x = 2 and x = -2
p5 <- p5 + geom_vline(xintercept = 2, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "blue")



# Add a square box with a red boundary color after negative 2 cutoff
p5 <- p5 + annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf, fill = NA, color = "red", size = 2)

print(p5)



# Create the ScatterView plot for all genes with "Diff" on the y-axis and "Rank" on the x-axis
p5 <- ScatterView(gdata_DC21, x = "Rank", y = "Diff", label = "Gene", model = "rank", max.overlaps = 50)

# Add label to the gene "XPO1" only
xpo1_data <- gdata_DC21[gdata_DC21$Gene == "XPO1", ]
p5 <- p5 + geom_text(data = xpo1_data, aes(label = Gene), color = "red", size = 5)

# Add a dashed line at y = 0 to represent the cutoff
p5 <- p5 + geom_hline(yintercept = 0, linetype = "dashed", color = "blue")

print(p5)


rankdata = gdata_DC21$Treatment - gdata_DC21$Control
names(rankdata) = gdata_DC21$Gene
RankView(rankdata,
         top = 100, max.overlaps = 50)




# Calculate the difference between Control and Treatment
gdata_DC21$Diff = gdata_DC21$Treatment - gdata_DC21$Control

# Calculate ranks for the differences in the "Diff" column and store it in the "Rank" column of the gdata_DC21 dataframe using the rank function
gdata_DC21$Rank = rank(gdata_DC21$Diff)

# Create the RankView for all genes
library(ggplot2)
library(ggdist)

# Calculate rankdata using all genes
rankdata <- gdata_DC21$Treatment - gdata_DC21$Control
names(rankdata) <- gdata_DC21$Gene

# Create the RankView plot for all genes
RankView(rankdata, top = 100, max.overlaps = 50, labels = c("XPO1"))





rankdata = gdata_DC21$Treatment - gdata_DC21$Control
names(rankdata) = gdata_DC21$Gene
RankView(rankdata)

p6 = ScatterView(gdata_DC21, x = "DMSO_day21", y = "Combo_day21", label = "Gene", 
                 model = "ninesquare", top = 20, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p6)

p8 = SquareView(gdata_DC21, label = "Gene", 
               x_cutoff = CutoffCalling(gdata_DC21$Control, 2), 
               y_cutoff = CutoffCalling(gdata_DC21$Treatment, 2), model = "ninesquare", top = 15, display_cut = TRUE, force = 2, max.overlaps = 50)
print(p8)


p1 = ScatterView(gdata_DC21, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, toplabels = c("XPO1"))
print(p1)

colnames(gdata_DC21)



head(gdata_DC21)
print(p1)


library(ggplot2)
library(ggdist)
gdata_DC21_100rem <- read.csv("gdata_DC21_top100AZAremoved.csv", header = TRUE, sep = ",")

# Create the ScatterView plot for all genes
p1 <- ScatterView(gdata_DC21, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = FALSE)

# Add labels to the XPO1 gene only
xpo1_data <- gdata_DC21[gdata_DC21$Gene == "XPO1"]
p1 <- p1 + geom_text(data = xpo1_data, aes(label = Gene), color = "red", size = 5)

# Add a rectangular box with red border around data points below 0 in the "Treatment" column
p1 <- p1 + geom_rect(data = gdata_DC21[gdata_DC21$Treatment < 0, ],
                     aes(xmin = Control, xmax = Treatment, ymin = -Inf, ymax = Inf),
                     color = "red", fill = NA)


print(p1)


# Create the ScatterView plot for all genes
library(ggplot2)
library(ggdist)

# Create the ScatterView plot for all genes
p1 <- ScatterView(gdata_DC21, "Control", "Treatment", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE)

# Add labels to the selected genes (XPO1 and MDM4)
selected_genes <- c("XPO1")
selected_data <- gdata_DC21[gdata_DC21$Gene %in% selected_genes, ]
p1 <- p1 + geom_text(data = selected_data, aes(label = Gene), color = "red", size = 2)

print(p1)

head(gdata_DC21)
write.csv(gdata_DC21, "gdata_DC21_0803.csv")

#install.packages("pathview")
library(pathview)
pathview::pathview(genedata, pathway.id = "hsa04115", species = "hsa", kegg.file = "hsa04115.xml")


Square8 = p1$data
idx=Square8$group=="top"
geneList1 = Square8$Diff
head(geneList1)
names(geneList1) = Square8$Gene[idx]
universe = Square8$Gene
kegg1 = EnrichAnalyzer(geneList = geneList1, universe = universe)
EnrichedView(kegg1, top = 20, bottom = 0)

# 9-square groups
Square9 = p8$data
idx=Square9$group=="topcenter"
geneList = Square9$Diff
head(geneList)
names(geneList) = Square9$Gene[idx]
universe = Square9$Gene


# Enrichment analysis
kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
EnrichedView(kegg1, top = 10, bottom = 0)

genedata = gdata_DC21[, c("Control","Treatment")]
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa04115", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa03013", organism = "hsa", sub = NULL)
arrangePathview(genedata, pathways = "hsa05166", organism = "hsa", sub = NULL)

hgtRes1 = EnrichAnalyzer(geneList1[1:200], method = "HGT", 
                         type = "Pathway", organism = "hsa")
head(hgtRes1@result)

ortRes1 = EnrichAnalyzer(geneList1[1:200], method = "ORT", 
                         type = "KEGG", organism = "hsa")
head(ortRes1@result)

gseRes1 = EnrichAnalyzer(geneList1, method = "GSEA", type = "Pathway", organism = "hsa")

df = hgtRes1@result
df$logFDR = -log10(df$p.adjust)
p = BarView(df[1:10,], "Description", 'logFDR')
p = p + labs(x = NULL) + coord_flip()
p

# Or use function barplot from enrichplot package
barplot(hgtRes1, showCategory = 10)

EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 1)
EnrichedView(hgtRes1, top = 10, bottom = 0, mode = 2)
dotplot(hgtRes1, showCategory = 10)

library(enrichplot)
hgtRes1@result$geneID = hgtRes1@result$geneName
cnetplot(hgtRes1, 10)
heatplot(hgtRes1, showCategory = 10, foldChange = geneList)
tmp <- pairwise_termsim(hgtRes1)
emapplot(tmp, showCategory = 5, layout = "kk")
# show GSEA results of one pathway
idx = which(gseRes1$NES>0)[1]
gseaplot(gseRes1, geneSetID = idx, title = gseRes1$Description[idx])
# show GSEA results of multiple pathways
gseaplot2(gseRes1, geneSetID = which(gseRes1$NES>0)[1:3])
enrichComb = EnrichAnalyzer(geneList1[1:200], type = "KEGG+REACTOME", organism = "hsa")
EnrichedView(enrichComb, top = 20)
## All pathways
enrich = EnrichAnalyzer(geneList = geneList1[1:200], type = "REACTOME", organism = "hsa")
EnrichedView(enrich, top = 20)
## All gene ontology
enrichGo = EnrichAnalyzer(geneList1[1:200], type = "GOBP", organism = "hsa")
EnrichedView(enrichGo, top = 20)
enrichPro = EnrichAnalyzer(geneList1[1:200], type = "Complex", organism = "hsa")
EnrichedView(enrichPro, top = 20)
enrich = EnrichAnalyzer(geneList1[1:200], type = "GOBP", limit = c(2, 80), organism = "hsa")
EnrichedView(enrich, top = 20)

enrich1 = EnrichAnalyzer(geneList1[1:200], type = "Pathway", organism = "hsa")
enrich2 = EnrichedFilter(enrich1)
EnrichedView(enrich1, top = 15)
EnrichedView(enrich2, top = 15)
##################################################################################################################################################################################

subset_df <- genesummary[, c("DMSO_day7.beta", "DMSO_day14.beta", "DMSO_day21.beta")]
subset_df_long <- gather(subset_df, key = "Condition", value = "beta")
subset_df_long$Condition <- factor(subset_df_long$Condition, levels = c("DMSO_day7.beta", "DMSO_day14.beta", "DMSO_day21.beta"))
ggplot(subset_df_long, aes(x = Condition, y = beta)) +
  geom_boxplot() +
  labs(x = "Condition", y = "Beta") +
  theme_minimal()
head(genesummary)
subset_df <- genesummary[, c("DMSO_day7.beta", "DMSO_day14.beta", "DMSO_day21.beta")]
head(subset_df)
colnames(subset_df)
subset_df_long <- gather(subset_df, key = "Condition", value = "beta")
subset_df_long$Condition <- factor(subset_df_long$Condition, levels = c("DMSO_day7.beta", "DMSO_day14.beta", "DMSO_day21.beta"))

ggplot(subset_df_long, aes(x = Condition, y = beta)) +
  geom_boxplot(width = 0.2) +
  labs(x = "Condition", y = "Beta") +
  theme_minimal()

subset_df <- genesummary[, c("AZA_day7.beta", "AZA_day14.beta", "AZA_day21.beta")]

subset_df_long <- gather(subset_df, key = "Condition", value = "beta")
subset_df_long$Condition <- factor(subset_df_long$Condition, levels = c("AZA_day7.beta", "AZA_day14.beta", "AZA_day21.beta"))

ggplot(subset_df_long, aes(x = Condition, y = beta)) +
  geom_boxplot(width = 0.2) +
  labs(x = "Condition", y = "Beta") +
  theme_minimal()

subset_df <- genesummary[, c("Combo_day7.beta", "Combo_day14.beta", "Combo_day21.beta")]

subset_df_long <- gather(subset_df, key = "Condition", value = "beta")
subset_df_long$Condition <- factor(subset_df_long$Condition, levels = c("Combo_day7.beta", "Combo_day14.beta", "Combo_day21.beta"))

ggplot(subset_df_long, aes(x = Condition, y = beta)) +
  geom_boxplot(width = 0.2) +
  labs(x = "Condition", y = "Beta") +
  theme_minimal()

