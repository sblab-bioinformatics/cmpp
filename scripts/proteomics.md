Co-authors: Kamal Kishore (CRUK-CI) and Sergio Martinez Cuesta (CRUK-CI and Deparment of Chemistry, Unversity of Cambridge)



### G4IPDB

```python
import csv

g4ipdb_file = csv.reader(open("20201120-G4IPDB_full_list.csv", "rU"), delimiter=',')
g4ipdb_lines = [row for row in g4ipdb_file]
g4ipdb_file.close()

g4ipdb_uniprot_all = []
g4ipdb_uniprot_validated = []

for line in g4ipdb_lines:
  if 'HUMAN' in line[7]:
    g4ipdb_uniprot_all.append(line[6].split()[0])
    if (line[12] == 'N/A' or line[12] == '') and (line[13] == 'N/A' or line[13] == '') and (line[14] == 'N/A' or line[14] == ''):
      continue
    else:
      g4ipdb_uniprot_validated.append(line[6].split()[0])


len(set(g4ipdb_uniprot_all)) # 79
len(set(g4ipdb_uniprot_validated)) # 39

g4ipdb_uniprot_all_output = open("20201120-G4IPDB_full_list_uniprot_all.txt", "w")
g4ipdb_uniprot_all_output.write("\n".join(list(set(g4ipdb_uniprot_all)))+'\n')
g4ipdb_uniprot_all_output.close()

g4ipdb_uniprot_validated_output = open("20201120-G4IPDB_full_list_uniprot_validated.txt", "w")
g4ipdb_uniprot_validated_output.write("\n".join(list(set(g4ipdb_uniprot_validated)))+'\n')
g4ipdb_uniprot_validated_output.close()
```



### Venn, volcano and tables

```r
library(qPLEXanalyzer)
library(gridExtra)
library(pander)
library(readxl)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(ggrepel)

# human annotation
data(human_anno)

# load data
metadata <- read_excel("../data/labelfree.xlsx", sheet = "metadata")
intensities <- read_excel("../data/labelfree.xlsx", sheet = "peptide_intensities")

# Convert intensities and metadata to MSnset
MSnset_data <- convertToMSnset(intensities, metadata = metadata, indExpData = c(7:30), Sequences = 2, Accessions = 6, rmMissing = FALSE)

# Remove values missing in all samples
torm <- apply(exprs(MSnset_data), 1 , function(x) all(is.na(x)))
MSnset_data <- MSnset_data[!torm,]

# Keep peptides identified in one protein and one protein group
tokeep <- which(fData(MSnset_data)[,4] == 1 & fData(MSnset_data)[,5] == 1)
MSnset_data <- MSnset_data[tokeep,]

# Convert SampleGroup to a factor
pData(MSnset_data)$SampleGroup <- as.factor(pData(MSnset_data)$SampleGroup)

# Impute missing value while replacing in control by its minimum value in the group
# In the other groups we keep peptides detected in at least half of the sample (here 2 in 4)
imputeMissingControl <- function(MSnSetObj,control=NULL)
{
  allgrps <- split(MSnSetObj,"SampleGroup")
  aind <- 1:length(allgrps)
  if(!is.null(control))
  {
    cind <- which(names(allgrps)==control)
    exprs(allgrps@x[[cind]])[is.na(exprs(allgrps@x[[cind]]))==TRUE] <- min(exprs(allgrps@x[[cind]]),na.rm=TRUE)
    }
  MSnSetObj_imputed <- unsplit(allgrps, pData(MSnSetObj)$SampleGroup)
  return(MSnSetObj_imputed)
}

MSnset_data <- imputeMissingControl(MSnset_data, control="DMSO_pUV")
MSnset_data <- imputeMissingControl(MSnset_data, control="NC_pUV")

imputeMissingExp <- function(MSnSetObj,noS=NULL)
{
  allgrps <- split(MSnSetObj,"SampleGroup")
  aind <- 1:length(allgrps)
  grpMissing <- lapply(allgrps, function(x) {apply(exprs(x),1,function(x) length(which(is.na(x)==TRUE)))})
  torm <- unique(unlist(lapply(grpMissing,function(x) which(x > noS))))
  MSnSetObj_imputed <- unsplit(allgrps, pData(MSnSetObj)$SampleGroup)
  MSnSetObj_imputed <- MSnSetObj_imputed[-torm,]
  return(MSnSetObj_imputed)
}

MSnset_imputed <- imputeMissingExp(MSnset_data, noS = 2)
MSnset_imputed <- impute(MSnset_imputed, method = "knn")

# The data is normalized via within group median scaling.
MSnset_norm <- groupScaling(MSnset_imputed, median)

# Select NC_pUV, PEG_Photo_PDS_pUV and Photo_PDS_pUV and change SampleGroup
MSnset_norm <- MSnset_norm[, c(5:12,17:20)]
pData(MSnset_norm)$SampleGroup <- factor(c(rep("control", 4), rep("photoPDS_2", 4), rep("photoPDS_1", 4)))


MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)


###########################
# Differential Expression #
###########################

# A **limma** based statistical analysis is carried out to find differentially-expressed proteins.
# The plot below shows the differentially expressed protein (in color).

contrasts <- c(photoPDS_1_vs_control = "photoPDS_1 - control", photoPDS_2_vs_control = "photoPDS_2 - control")

diffstats <- computeDiffStats(MSnset_Pnorm, contrasts=contrasts)

diffexp <- list()
for (i in 1:length(contrasts))
  diffexp[[i]] <- getContrastResults(diffstats = diffstats, contrast = contrasts[i], writeFile = TRUE)


# photoPDS_1_vs_control   diffexp[[1]]
photoPDS_1_vs_control <- data.table(diffexp[[1]])
photoPDS_1_vs_control[Accessions=="Q9UNY4", GeneSymbol := "TTF2"]

# photoPDS_2_vs_control   diffexp[[2]]
photoPDS_2_vs_control <- data.table(diffexp[[2]])
photoPDS_2_vs_control[Accessions=="Q9UNY4", GeneSymbol := "TTF2"]

length(union(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions, photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions)) # 250 proteins
write(union(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions, photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions), "figures/20201121_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union.txt", sep = "\n")


#################
# Venn diagrams #
#################

# photoPDS_1_vs_control and photoPDS_2_vs_control and all G4IPDB
all <- fread("20201120-G4IPDB_full_list_uniprot_all.txt", header = FALSE)

venn.plot <- draw.triple.venn(
  area1 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1]),
  area2 = nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]),
  area3 = nrow(all),
  n12 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions]),
  n23 = nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% all$V1]),
  n13 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% all$V1]),
  n123 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions][Accessions %in% all$V1]),
  category = c("photoPDS-1\n (240)", "photoPDS-2\n (185)", "G4IPDB interactors\n (79)"),
  fill = c("green3", "green4", "darkgoldenrod2"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.12, 0.12, 0.12),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.075)

pdf("figures/20201121_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_venn_known_g4ipdb_all.pdf", width = 14/2.54, height = 14/2.54)
g <- grid.draw(venn.plot)
dev.off()



# photoPDS_1_vs_control and photoPDS_2_vs_control and validated G4IPDB
validated <- fread("20201120-G4IPDB_full_list_uniprot_validated.txt", header = FALSE)

venn.plot <- draw.triple.venn(
  area1 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1]),
  area2 = nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]),
  area3 = nrow(validated),
  n12 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions]),
  n23 = nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% validated$V1]),
  n13 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% validated$V1]),
  n123 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions][Accessions %in% validated$V1]),
  category = c("photoPDS-1\n (240)", "photoPDS-2\n (185)", "G4IPDB validated interactors\n (39)"),
  fill = c("green3", "green4", "darkgoldenrod2"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.12, 0.12, 0.12),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.075)

pdf("figures/20201121_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_venn_known_g4ipdb_validated.pdf", width = 14/2.54, height = 14/2.54)
g <- grid.draw(venn.plot)
dev.off()



###########
# Volcano #
###########

# photoPDS_1_vs_control and all G4IPDB
gg <- ggplot(photoPDS_1_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green3") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("1 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-12, 12)) +
annotate("text", x = 10, y = 8, label = sprintf("n = %s", nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave('figures/20201121_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all.pdf', width = 12, height = 12, units = "cm")


# photoPDS_1_vs_control and all G4IPDB with labels
gg <- ggplot(photoPDS_1_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green3") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("1 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 8, label = sprintf("n = %s", nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_label_repel(data = photoPDS_1_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX1", "DDX24", "HMGB2")], aes(label = GeneSymbol), force = 10)

ggsave('figures/20201203_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all.pdf', width = 12, height = 12, units = "cm")


# photoPDS_1_vs_control and validated G4IPDB
gg <- ggplot(photoPDS_1_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green3") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% validated$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("1 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-12, 12)) +
annotate("text", x = 10, y = 8, label = sprintf("n = %s", nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave('figures/20201121_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_validated.pdf', width = 12, height = 12, units = "cm")


# photoPDS_2_vs_control and all G4IPDB
gg <- ggplot(photoPDS_2_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green4") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("2 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-12, 12)) +
annotate("text", x = 10, y = 8, label = sprintf("n = %s", nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave('figures/20201121_PEG_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all.pdf', width = 12, height = 12, units = "cm")


# photoPDS_2_vs_control and all G4IPDB with labels
gg <- ggplot(photoPDS_2_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green4") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("2 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 8, label = sprintf("n = %s", nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_label_repel(data = photoPDS_2_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX1", "DDX24", "HMGB2")], aes(label = GeneSymbol), force = 10)

ggsave('figures/20201203_PEG_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all.pdf', width = 12, height = 12, units = "cm")


# photoPDS_2_vs_control and validated G4IPDB
gg <- ggplot(photoPDS_2_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green4") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% validated$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("2 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-12, 12)) +
annotate("text", x = 10, y = 8, label = sprintf("n = %s", nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave('figures/20201121_PEG_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_validated.pdf', width = 12, height = 12, units = "cm")



##########
# Tables #
##########
write.table(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][,.(Accessions, GeneSymbol, log2FC, adj.P.Val)], "20201203_Photo_PDS_pUV_vs_NC_pUV.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1][,.(Accessions, GeneSymbol, log2FC, adj.P.Val)], "20201203_PEG_Photo_PDS_pUV_vs_NC_pUV_union.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
```



### Barplots molecular function

Prepare file:

```python
ifile = open("20201121_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot.txt", "r")
ilines = ifile.readlines()
ifile.close()

uniprot_keyword = []
uniprot_go = []

for l in ilines[1:]:
  fields = l.replace("\n", "").split("\t")
  uniprot = fields[0]
  ks = fields[1].split(";")
  gos = fields[2].split("; ")
  for k in ks:
    uniprot_keyword.append((uniprot, k))
  for go in gos:
    go_m = go.split(" [")[0]
    uniprot_go.append((uniprot, go_m))

ofile_k = open("20201121_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot_keyword.txt", "w")
ofile_k.write("\n".join(["\t".join(e) for e in uniprot_keyword]))
ofile_k.close()

ofile_g = open("20201121_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot_go.txt", "w")
ofile_g.write("\n".join(["\t".join(e) for e in uniprot_go]))
ofile_g.close()
```

Barplots uniprot keyword and go:

```r
library(data.table)
library(ggplot2)

# keyword
data_keyword <- fread("20201121_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot_keyword.txt", header = FALSE)
setnames(data_keyword, c("uniprot", "keyword"))

## mf
data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])
#13:               RNA-binding  80
#15:               DNA-binding  58
#22:               ATP-binding  42
#23:             Metal-binding  41
#26:         Ribonucleoprotein  34
#29:                 Hydrolase  29
#33:                 Repressor  25
#42:                 Activator  17
#44:                  Helicase  17
#48:       Chromatin regulator  15

data_keyword_table <- data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])[V1 %in% c("RNA-binding", "DNA-binding", "ATP-binding", "Metal-binding", "Ribonucleoprotein", "Hydrolase", "Repressor", "Activator", "Helicase", "Chromatin regulator")]

gg <- ggplot(data = data_keyword_table, aes(x = N, y = reorder(V1, N))) +
geom_bar(stat="identity") +
xlab("Number of proteins") +
ylab("") +
ggtitle("UniprotKB Keyword Molecular Function") +
theme_bw() +
theme(axis.title = element_text(size=16), axis.text = element_text(size=12, color = "black"), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave("figures/20201121_union_uniprot_keyword_molecular_function.pdf", height = 4, width = 6)


# bp
data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])
#16:           mRNA processing  54
#17:             Transcription  53
#18:  Transcription regulation  52
#20:             mRNA splicing  47
#28:                 Transport  30
#31:    Host-virus interaction  27
#32:                Cell cycle  26
#41:                DNA damage  18
#45:             Cell division  16
#46:                DNA repair  16

data_keyword_table <- data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])[V1 %in% c("mRNA processing", "Transcription", "Transcription regulation", "mRNA splicing", "Transport", "Host-virus interaction", "Cell cycle", "DNA damage", "Cell division", "DNA repair", "Immunity")]
data_keyword_table[V1=="Transcription", N:=53+52]

gg <- ggplot(data = data_keyword_table[V1 != "Transcription regulation"], aes(x = N, y = reorder(V1, N))) +
geom_bar(stat="identity") +
xlab("Number of proteins") +
ylab("") +
ggtitle("UniprotKB Keyword Biological Process") +
theme_bw() +
theme(axis.title = element_text(size=16), axis.text = element_text(size=12, color = "black"), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1), plot.title = element_text(face="bold", size=15, hjust = 0.5))

ggsave("figures/20201121_union_uniprot_keyword_biological_process.pdf", height = 4, width = 6)



# go
data_go <- fread("20201121_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot_go.txt", header = FALSE)
setnames(data_go, c("uniprot", "go"))

data.table(sort(table(data_go$go), decreasing = T)[1:50])
# 1:                                                           RNA binding 144
# 2:                                                           DNA binding  47
# 6:                                                     metal ion binding  26
# 7:                                                     chromatin binding  23
# 8:                                                      cadherin binding  18
# 9:                                                        enzyme binding  16
#10:                                          transcription factor binding  15
#11:                                                protein kinase binding  14
#12:                                                       ATPase activity  12
#13:                                           double-stranded RNA binding  12

data_go_table <- data.table(sort(table(data_go$go), decreasing = T)[1:50])[V1 %in% c("RNA binding", "DNA binding", "metal ion binding", "chromatin binding", "cadherin binding", "enzyme binding", "transcription factor binding", "protein kinase binding", "ATPase activity", "double-stranded RNA binding")]

gg <- ggplot(data = data_go_table, aes(x = N, y = reorder(V1, N))) +
geom_bar(stat="identity") +
xlab("Number of proteins") +
ylab("") +
ggtitle("UniprotKB GO Molecular Function") +
theme_bw() +
theme(axis.title = element_text(size=16), axis.text = element_text(size=12, color = "black"), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1), plot.title = element_text(face="bold", size=15, hjust = 0.5))

ggsave("figures/20201121_union_uniprot_go_molecular_function.pdf", height = 4, width = 6)
```



### Venn, volcano and tables (enhanced)

```r
library(data.table)
library(ggplot2)
library(VennDiagram)
library(ggrepel)


# photoPDS_1_vs_control   
photoPDS_1_vs_control <- fread("Photo_PDS_pUV_vs_NC_pUV.txt")

# photoPDS_2_vs_control
photoPDS_2_vs_control <- fread("PEG_Photo_PDS_pUV_vs_NC_pUV.txt")

length(union(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions, photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions)) # 256 proteins
write(union(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions, photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions), "20210112_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union.txt", sep = "\n")


#################
# Venn diagrams #
#################

# photoPDS_1_vs_control and photoPDS_2_vs_control and all G4IPDB
all <- fread("../20201120-G4IPDB_full_list_uniprot_all.txt", header = FALSE)

venn.plot <- draw.triple.venn(
  area1 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1]),
  area2 = nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]),
  area3 = nrow(all),
  n12 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions]),
  n23 = nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% all$V1]),
  n13 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% all$V1]),
  n123 = nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][Accessions %in% photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1]$Accessions][Accessions %in% all$V1]),
  category = c("photoPDS-1\n (248)", "photoPDS-2\n (209)", "G4IPDB interactors\n (79)"),
  fill = c("green3", "green4", "darkgoldenrod2"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.12, 0.12, 0.12),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.075)

pdf("figures/20210112_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_venn_known_g4ipdb_all.pdf", width = 14/2.54, height = 14/2.54)
g <- grid.draw(venn.plot)
dev.off()


###########
# Volcano #
###########

# photoPDS_1_vs_control and all G4IPDB
gg <- ggplot(photoPDS_1_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green3") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("1 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-12, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave('figures/20210112_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all.pdf', width = 12, height = 12, units = "cm")


# photoPDS_1_vs_control and all G4IPDB with labels
gg <- ggplot(photoPDS_1_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green3") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("1 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_label_repel(data = photoPDS_1_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX1", "DDX24", "HMGB2")], aes(label = GeneSymbol), force = 10)

ggsave('figures/20210112_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all_labels.pdf', width = 12, height = 12, units = "cm")


# photoPDS_1_vs_control and all G4IPDB with geom_text_repel
gg <- ggplot(photoPDS_1_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green3") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("1 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_text_repel(data = photoPDS_1_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX1", "DDX24", "HMGB2")], aes(label = GeneSymbol), min.segment.length = 0, seed = 42, box.padding = 0.5, force = 5, nudge_x = 4)

ggsave('figures/20210112_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all_geom_text_repel.pdf', width = 12, height = 12, units = "cm")


# photoPDS_1_vs_control and all G4IPDB with geom_text_repel and nonudge
gg <- ggplot(photoPDS_1_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green3") +
geom_point(data = photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("1 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_text_repel(data = photoPDS_1_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX1", "DDX24", "HMGB2")], aes(label = GeneSymbol), min.segment.length = 0, seed = 42, box.padding = 0.5, force = 5)

ggsave('figures/20210112_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all_geom_text_repel_nonudge.pdf', width = 12, height = 12, units = "cm", useDingbats = FALSE)


# photoPDS_2_vs_control and all G4IPDB
gg <- ggplot(photoPDS_2_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green4") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("2 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-12, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave('figures/20210112_PEG_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all.pdf', width = 12, height = 12, units = "cm")


# photoPDS_2_vs_control and all G4IPDB with labels
gg <- ggplot(photoPDS_2_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green4") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("2 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_label_repel(data = photoPDS_2_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX1", "DDX24", "HMGB2")], aes(label = GeneSymbol), force = 10)

ggsave('figures/20210112_PEG_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all_labels.pdf', width = 12, height = 12, units = "cm")


# photoPDS_2_vs_control and all G4IPDB with labels and no DDX1
gg <- ggplot(photoPDS_2_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green4") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("2 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_label_repel(data = photoPDS_2_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX24", "HMGB2")], aes(label = GeneSymbol), force = 10)

ggsave('figures/20210112_PEG_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all_labels_noDDX1.pdf', width = 12, height = 12, units = "cm")


# photoPDS_2_vs_control and all G4IPDB with labels and no DDX1 with geom_text_repel
gg <- ggplot(photoPDS_2_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green4") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("2 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_text_repel(data = photoPDS_2_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX24", "HMGB2")], aes(label = GeneSymbol), min.segment.length = 0, seed = 42, box.padding = 0.5, force = 5, nudge_x = 4)

ggsave('figures/20210112_PEG_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all_labels_noDDX1_geom_text_repel.pdf', width = 12, height = 12, units = "cm")


# photoPDS_2_vs_control and all G4IPDB with labels and no DDX1 with geom_text_repel and nonudge
gg <- ggplot(photoPDS_2_vs_control, aes(x = log2FC, y = -log10(adj.P.Val))) +
geom_point(size = 1, alpha = 1, col = "lightgray") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "green4") +
geom_point(data = photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1 & Accessions %in% all$V1], aes(x = log2FC, y = -log10(adj.P.Val)), col = "darkgoldenrod2") +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("2 over 3") +
geom_vline(xintercept = 0, linetype="longdash") +
theme_classic() +
coord_cartesian(xlim = c(-6, 12)) +
annotate("text", x = 10, y = 10, label = sprintf("n = %s", nrow(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1])), size = 6) +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5)) +
geom_text_repel(data = photoPDS_2_vs_control[GeneSymbol %in% c("SMARCA4", "UHRF1", "RBM22", "TTF2", "DDX24", "HMGB2")], aes(label = GeneSymbol), min.segment.length = 0, seed = 42, box.padding = 0.5, force = 5)

ggsave('figures/20210112_PEG_Photo_PDS_pUV_vs_NC_pUV_volcano_g4ipdb_all_labels_noDDX1_geom_text_repel_nonudge.pdf', width = 12, height = 12, units = "cm", useDingbats = FALSE)


##########
# Tables #
##########
write.table(photoPDS_1_vs_control[adj.P.Val < 0.05 & log2FC > 1][,.(Accessions, GeneSymbol, log2FC, adj.P.Val)], "20210112_Photo_PDS_pUV_vs_NC_pUV.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(photoPDS_2_vs_control[adj.P.Val < 0.05 & log2FC > 1][,.(Accessions, GeneSymbol, log2FC, adj.P.Val)], "20210112_PEG_Photo_PDS_pUV_vs_NC_pUV.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
```



### Barplots molecular function

Prepare file:

```python
ifile = open("20210112_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot.txt", "r")
ilines = ifile.readlines()
ifile.close()

uniprot_keyword = []
uniprot_go = []

for l in ilines[1:]:
  fields = l.replace("\n", "").split("\t")
  uniprot = fields[0]
  ks = fields[1].split(";")
  gos = fields[2].split("; ")
  for k in ks:
    uniprot_keyword.append((uniprot, k))
  for go in gos:
    go_m = go.split(" [")[0]
    uniprot_go.append((uniprot, go_m))

ofile_k = open("20210112_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot_keyword.txt", "w")
ofile_k.write("\n".join(["\t".join(e) for e in uniprot_keyword]))
ofile_k.close()

ofile_g = open("20210112_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot_go.txt", "w")
ofile_g.write("\n".join(["\t".join(e) for e in uniprot_go]))
ofile_g.close()
```

Barplots uniprot keyword and go:

```r
library(data.table)
library(ggplot2)

# keyword
data_keyword <- fread("20210112_Photo_PDS_pUV_vs_NC_pUV_vs_PEG_Photo_PDS_pUV_vs_NC_pUV_union_uniprot_keyword.txt", header = FALSE)
setnames(data_keyword, c("uniprot", "keyword"))

data_keyword_table <- data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])[V1 %in% c("RNA-binding", "DNA-binding", "ATP-binding", "Metal-binding", "Ribonucleoprotein", "Hydrolase", "Repressor", "Activator", "Helicase", "Chromatin regulator")]

gg <- ggplot(data = data_keyword_table, aes(x = N, y = reorder(V1, N))) +
geom_bar(stat="identity") +
xlab("Number of proteins") +
ylab("") +
ggtitle("UniprotKB Keyword Molecular Function") +
theme_bw() +
theme(axis.title = element_text(size=16), axis.text = element_text(size=12, color = "black"), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave("figures/20210112_union_uniprot_keyword_molecular_function.pdf", height = 4, width = 6)


# bp
data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])
#17:             Transcription  52
#18:  Transcription regulation  50
#19:           mRNA processing  49
#24:             mRNA splicing  41
#27:                 Transport  31
#28:    Host-virus interaction  29
#30:                Cell cycle  27
#42:                DNA damage  17
#44:             Cell division  16
#49:                DNA repair  15


data_keyword_table <- data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])[V1 %in% c("mRNA processing", "Transcription", "Transcription regulation", "mRNA splicing", "Transport", "Host-virus interaction", "Cell cycle", "DNA damage", "Cell division", "DNA repair")]

gg <- ggplot(data = data_keyword_table, aes(x = N, y = reorder(V1, N))) +
geom_bar(stat="identity") +
xlab("Number of proteins") +
ylab("") +
ggtitle("UniprotKB Keyword Biological Process") +
theme_bw() +
theme(axis.title = element_text(size=16), axis.text = element_text(size=12, color = "black"), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1), plot.title = element_text(face="bold", size=15, hjust = 0.5))

ggsave("figures/20210112_union_uniprot_keyword_biological_process.pdf", height = 4, width = 6)
```
