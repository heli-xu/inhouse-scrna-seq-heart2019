library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)


setwd("C:/Temple/scRNA/sc_2019")

#############Import data (4 samples)###############################
WT_SHAM.data <- Read10X(data.dir = "raw/WT-SHAM/")
WT_SHAM <- CreateSeuratObject(counts=WT_SHAM.data,project = "Heart2019",
                              min.cells = 3,min.features = 200)

WT_MI.data <- Read10X(data.dir="raw/WT-MI/")
WT_MI <- CreateSeuratObject(counts=WT_MI.data,project = "Heart2019",
                            min.cells = 3,min.features = 200)

KO_SHAM.data <- Read10X(data.dir="raw/KO-SHAM/")
KO_SHAM <- CreateSeuratObject(counts=KO_SHAM.data,project = "Heart2019",
                              min.cells = 3,min.features = 200)

KO_MI.data <- Read10X(data.dir = "raw/KO-MI/")
KO_MI <- CreateSeuratObject(counts=KO_MI.data,project = "Heart2019",
                            min.cells = 3,min.features = 200)

#adding column to metadata
WT_SHAM$sample <- sample("WT_SHAM", size = ncol(WT_SHAM), replace = TRUE)
WT_MI$sample <- sample("WT_MI", size = ncol(WT_MI), replace = TRUE)
KO_MI$sample <- sample("KO_MI", size=ncol(KO_MI),replace = TRUE)
KO_SHAM$sample <- sample("KO_SHAM", size = ncol(KO_SHAM),replace=TRUE)

head(x = WT_MI@meta.data, 5)

##merged datasets to a Seurat object
merged_MI <- merge(x=WT_MI, KO_MI,
                   add.cell.id = c("WT", "KO"))

head(merged_MI@meta.data)

#

####Quality Control##############

load("data/mito_gene_name_mouse.rdata")

# Number of UMIs assigned to mitochondrial genes
metadata <- merged_MI@meta.data
metadata$cells <- rownames(metadata)
##keep rownames, add another column SO AS TO incorporate it back to 
#Seurat Object
#in contrast to mutate() or rownames to column() rownames taken out

#rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


counts <- GetAssayData(merged_MI, slot = "counts")
##Number of UMIs assigned to mitochondrial genes
metadata$mtUMI <- 
Matrix::colSums(counts[which(rownames(counts) %in% mt),], na.rm = T)
#if don't understand, used code below instead. 

# #### Ran's code(subset for understanding) ####
# # subset data to small matrix
# matrix_tmp = KO_MI.data[,1:5]
# # get rownames (gene name)
# rownames_tmp  = rownames(matrix_tmp)
# # convert matrix to tibble and add gene name
# df = matrix_tmp %>% 
#   as_tibble() %>% 
#   mutate(gene = rownames_tmp) %>% 
#   dplyr::select(gene, everything())
# 
# # tidy data by pivot longer
# df2 = df %>% pivot_longer(cols = -gene, names_to= "cell.id") %>% 
#   dplyr::select(cell.id,gene,value)
# 
# # categorize gene value into mito or not mito
# df3 = df2 %>% 
#   mutate(mito = ifelse(gene%in%mt,"yes mito","not mito"))
# 
# # aggrgate mito gene counts by cell id
# mito_results = df3 %>% 
#   group_by(cell.id, mito) %>% 
#   summarize(mtUMI=sum(value)) %>% 
#   ungroup()
# 
# 
# #### Cleaning Mito stuff  ####
# rownames_WT  = rownames(WT_MI.data)
# df1 = WT_MI.data%>% 
#   as_tibble() %>% 
#   mutate(gene = rownames_WT) %>% 
#   dplyr::select(gene, everything())
# 
# df2 = df1 %>% pivot_longer(cols = -gene, names_to= "cell.id") %>% 
#   dplyr::select(cell.id,gene,value)
# 
# df3 = df2 %>% 
#   mutate(mito = ifelse(gene%in%mt,"yes","no"))
# 
# WT_mt=df3 %>% 
#   group_by(cell.id, mito) %>% 
#   summarize(mtUMI=sum(value)) %>% 
#   ungroup()%>% 
#   filter(mito=="yes") %>% 
#   mutate(cell.id=paste0("WT_",cell.id))
# 
# rownames_KO = rownames(KO_MI.data)
# df1 = KO_MI.data%>% 
#   as_tibble() %>% 
#   mutate(gene = rownames_KO) %>% 
#   dplyr::select(gene, everything())
# 
# ##takes a while might crash
# df2 = df1 %>% pivot_longer(cols = -gene, names_to= "cell.id") %>% 
#   dplyr::select(cell.id,gene,value)
# 
# df3 = df2 %>% 
#   mutate(mito = ifelse(gene%in%mt,"yes","no"))
# 
# KO_mt=df3 %>% 
#   group_by(cell.id, mito) %>% 
#   summarize(mtUMI=sum(value)) %>% 
#   ungroup() %>% 
#   filter(mito=="yes") %>% 
#   mutate(cell.id=paste0("KO_",cell.id))
# 
# mt_df=bind_rows(WT_mt,KO_mt) %>% 
#   select(-mito)
# 
# metadata = left_join(metadata,mt_df,by= c("cells"="cell.id")) 
# ##this takes the rownames out, so let's add it back.
# metadata <- column_to_rownames(metadata,var = "cells")

# Calculate of mitoRatio per cell
metadata$mitoRatio <- metadata$mtUMI/metadata$nUMI

#add number of genes per UMI for each cell
# this metric with give us an idea of the complexity of our dataset
#(more genes detected per UMI, more complex our data)
metadata$log10GenesPerUMI <- 
  log10(metadata$nGene) / log10(metadata$nUMI)


##save the metadata back to Seurat project
merged_MI@meta.data <- metadata

save(merged_MI,file="data/merged_MI.Rdata")

load("data/merged_MI.rdata")

######some QC metrics
##Cell Count (expected 12k?)
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")#+geom_vline(xintercept =1, color = 'green', linetype =5)
#+geom_hline(yintercept = 12000, col = 'red', lty = 2)

# metadata %>% count(sample)
# 
# hchart(metadata %>% count(sample),"bar",hcaes(x=sample,y = n, group =sample))


##UMI counts (transcripts) per cell (>500)
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

##Genes detected per cell (>300, single peak)
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

#the number of genes versus the numnber of UMIs coloured by 
#the fraction of mitochondrial reads. 
#Mitochondrial read fractions are only high 
#in particularly low count cells with few detected genes.
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

##Distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

#overall complexity of the gene expression by genes detected per UMI
#Sometimes we can detect contamination with low complexity cell types like RBC
#Generally, we expect the novelty score to be above 0.80.
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

##Cell-level, gene-level filtering:
#consider the joint effects of these metrics when setting thresholds
#and set them to be as permissive as possible to 
#avoid filtering out viable cell populations unintentionally.



########Cell cycle scoring######

load("data/cell_cycle_mouse.rdata")

s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")


merged_MI_phase <- NormalizeData(merged_MI) %>% 
  CellCycleScoring(g2m.features = g2m_genes, 
                   s.features = s_genes) %>% 
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000, 
                       verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA() 

save(merged_MI_phase,file = "data/merged_MI_phase.rdata")

load("data/merged_MI_phase.rdata")

DimPlot(merged_MI_phase,
        reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")
##We do not see large differences due to cell cycle phase.--no need to regress out

###Integration or not?####
merged_MI_no_integ <- NormalizeData(merged_MI) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dim=1:20)

save(merged_MI_no_integ,file = "data/merged_MI_no_integ.rdata")

load("data/merged_MI_no_integ.rdata")

DimPlot(merged_MI_no_integ, reduction = "umap", group.by = "sample")

##there's condition-specific clustering (WT vs KO)==>needs integration


##different normalization to use variable feature for cell type id
###SCtransform##############
# Split seurat object by condition to 
#perform cell cycle scoring and SCT on all samples
split_MI <- SplitObject(merged_MI, split.by = "sample")

split_MI <- split_MI[c("WT_MI", "KO_MI")] %>%  
  map(~.x %>% 
        NormalizeData() %>% 
        CellCycleScoring(g2m.features = g2m_genes, 
                         s.features = s_genes) %>% 
        SCTransform(verbose=FALSE))

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_MI, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_MI <- PrepSCTIntegration(object.list = split_MI, 
                                 anchor.features = integ_features)

# CCA: Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_MI, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features,
                                        verbose=FALSE)

# Integrate across conditions
MI_integrated <- IntegrateData(anchorset = integ_anchors, 
                                  normalization.method = "SCT",
                                  verbose=FALSE)


MI_integrated<- RunPCA(object = MI_integrated,verbose = FALSE)


MI_integrated <- RunUMAP(MI_integrated,
                            dims = 1:30,
                            reduction = "pca")

save(MI_integrated, file ="results/integrated_MI.rdata")

load("results/integrated_MI.rdata")

DimPlot(MI_integrated,reduction = "umap", group.by = "sample")
##no clustering by group

ElbowPlot(object = MI_integrated, ndims = 40)
##where the elbow appears is usually the threshold for 
##identifying the majority of the variation


MI_integrated <- FindNeighbors(object = MI_integrated, 
                                  dims = 1:30) %>% 
  FindClusters(resolution=0.8)

DimPlot(MI_integrated,
        reduction = "umap",
        split.by = "Phase")

DimPlot(MI_integrated,
        reduction = "umap",
        group.by = "sample")

save(MI_integrated,file = "results/MI_integrated_cluster.rdata")

load("results/MI_integrated_cluster.rdata")

###Exploration of the PCs driving the different clusters####
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(MI_integrated, 
                     vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(MI_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
#don't run, show figure
library(cowplot)
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

##We can see how the clusters are represented by the different PCs. 

print(MI_integrated[["pca"]], dims = 1:5, nfeatures = 5)


####Exploring known cell type markers#####

# Select the RNA counts slot to be the default assay
DefaultAssay(MI_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
MI_integrated <- NormalizeData(MI_integrated, verbose = FALSE)

##Endothelial cells (marker protein CD31 and CD102)
FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Pecam1","Icam2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#looks like 5 is EC


##Erythrocytes (Hemoglobin)
FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Hbb-bs", "Hba-a1", "Hba-a2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#looks like 1 is RBC

##Leukocytes (CD45, CD11b)
FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Ptprc", "Itgam"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#2,3,4,6 seem like it


##resident mesenchymal cells (fibroblast, VSMC, mural cells)
#no Thy1(CD90), very little Atxn1 (Sca-1)
#VSMC: Cxcl12+,Eln+,Fbln5+
#MS:Postn+,Vcan+,(Sparc, stromal cells)

FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Cxcl12","Eln","Fbln5","Postn","Vcan"), 
            sort.cell = TRUE,
            split.by = "sample",
            min.cutoff = 'q10', 
            label = TRUE)


FeaturePlot(MI_integrated,
            reduction = "umap",
            features = c("Postn","Sparc" ),
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

#########Identification of conserved markers in all conditions###########
##as opposed to FindAllMarkers()--used in a single sample/group
##FindConservedMarkers() use the original counts and not the integrated data
# BiocManager::install('multtest')
# install.packages('metap')

DefaultAssay(MI_integrated) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(MI_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)


##Look for markers with large differences in expression between pct.1 and pct.2 
##and larger fold changes.

load("data/gene_annotations_mouse.rdata")

# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(MI_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    mutate(cluster_id = cluster)
}

conserved_markers <- map_dfr(c(0:6), get_conserved)
##run the clusters that are UNIDENTIFIED, otherwise very long 
##if not enough cells for this function, use FindAllMarkers()


# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (WT_MI_avg_logFC + KO_MI_avg_logFC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc) %>% 
  ungroup()

write.csv(top10,file="results/markers for clusters.csv")

top10_filter <- top10 %>%
  select(1,14,15,4,5,9,10)

#take a closer look at cluster 0-3, should be macro, mono and DC.
#Lyz2: LyzM all myeloid lineage, Macrophage markers: Adgre1: F4/80; Cd14; 
#Fcgr3(CD16),Fcgr1 (CD64) everywhere, Cd68 everywhere, Tfrc(CD71) nowhere;
#Ccr5
##inf.Macro:Itgam,Ccr2;Resident: Adgre1; Fcgr1 everywhere
FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Itgam","Ccr2","Cd14" ), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Adgre1","Fcgr3","Ccr5","Ly6c2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Ccr2","Cd14","Ly6c2", "Ccr5", "Itgam"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

DotPlot(MI_integrated,
        features = c("Ccr2","Cd14","Ly6c2", "Ccr5", "Itgam"), 
        cols = c("lightgrey","blue"),
        scale = TRUE, 
        scale.by = "radius")

VlnPlot(MI_integrated, 
            features = c("Ccr2"),
        split.by = "sample")



# No Itgax (CD11c); No Cd103;sigh....
##Plasmacytoid DC: Gzmb, Serpinf1, itm2c
FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Cd8a", "Gzmb","Il3ra","Itm2c"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#T cell marker; CD8 very little, CD4 nothing;
FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Crem","Cd69","Sell"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#Megakaryocyte
#Ppbp,Itga2b(CD41),Gp5 (CD42d) nothing, Mpl (Tpo R) nothing,.
FeaturePlot(MI_integrated, 
            reduction = "umap", 
            features = c("Ppbp","Lifr","Cxcr4"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# Rename all identities
MI_integrated <- RenameIdents(object = MI_integrated,
                              "0"="Dendritic cells/Monocytes",
                              "1"="Erythrocytes",
                              "2"="Resident macrophages",
                              "3"="Inflammatory macrophages",
                              "4"="Granulocytes",
                              "5"="Endothelial cells/Mesenchymal cells",
                              "6"="M2-polarized macrophages")

MI_integrated <- RenameIdents(object = MI_integrated,
                              "0"="DCs/Mono",
                              "1"="RBC",
                              "2"="Res.M",
                              "3"="Inf.M",
                              "4"="Granu",
                              "5"="EC/MSC",
                              "6"="M2")


DimPlot(object = MI_integrated,
        reduction = "umap",
        label = T,
        label.size = 3,
        repel = T)

save(MI_integrated, file = "results/MI_labelled.rdata")

load("results/MI_labelled.rdata")


n_cells <- FetchData(MI_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample) %>% 
  group_by(sample) %>% 
  mutate(percentage=n/sum(n)*100) %>% 
  pivot_longer(cols = c(n, percentage), names_to = "cell_count") %>% 
  pivot_wider(names_from = ident, values_from = value)

#######heatmap#########
DefaultAssay(MI_integrated) <- "integrated"

#copy top10 code and make it top20

heatmap_genes=top10$gene
a = top10[,"gene"]

##why top20[,"gene"]doesn't work? it worked on conserved_markers
##use class() to check, top10 is tbl_df, conserved_marker is traditional df. 
##somehow tbl_df[]subset returns df, $subset returns characters.

DoHeatmap(MI_integrated, features = heatmap_genes,
          group.by = "ident", slot = "scale.data")


VlnPlot(MI_integrated, 
        features = c("Ccl2"),
        split.by = "sample")

VlnPlot(MI_integrated,
features = c("Ccr2"),
split.by = "sample")

#################
VlnPlot(MI_integrated,
features = c("Ccl4"),
split.by = "sample")

VlnPlot(MI_integrated,
features = c("Ccl3"),
split.by = "sample")


VlnPlot(MI_integrated,
features = c("Ccr5"),
split.by = "sample")
####################

VlnPlot(MI_integrated,
features = c("Ccl7"),
split.by = "sample")
###Ccr2

VlnPlot(MI_integrated,
features = c("Cxcr4"),
split.by = "sample")
##Cxcl12

###############
VlnPlot(MI_integrated,
features = c("Cxcl10"),
split.by = "sample")

VlnPlot(MI_integrated,
features = c("Cxcr3"),
split.by = "sample")
#####################
VlnPlot(MI_integrated,
features = c("Ccl6"),
split.by = "sample")
