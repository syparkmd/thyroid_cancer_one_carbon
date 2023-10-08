
#loae library
if(T){
  library(tidyverse)
  library(Seurat)
  library(ggsci)
  library(scales)
  library(fgsea)
  library(cowplot)
  library(mgsub)
  library(ggsignif)
  library(ggrepel)
}

#output path
out_dir="/path/to/out_dir"

source("/path/to/sypark_scRNAseq_plot_utils.R")


#design
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), panel.border = element_blank(), axis.line = element_line(), plot.title = element_text(size=20, face="bold"))
cancer_type_pal<- pal_d3("category10")(10)[c(1,3,2)]
names(cancer_type_pal) <- c('NT','PT','AT')

#gene_sets
if(T){
  h_list <- GSEABase::getGmt('/path/to/h.all.v7.5.1.symbols.gmt') %>% GSEABase::geneIds()
  c5_list <- GSEABase::getGmt('/path/to/c5.all.v7.4.symbols.gmt') %>% GSEABase::geneIds()
  c2_list <- GSEABase::getGmt('/path/to/c2.all.v7.4.symbols.gmt') %>% GSEABase::geneIds()
  c7_list <- GSEABase::getGmt('/path/to/c7.all.v2022.1.Hs.symbols.gmt') %>% GSEABase::geneIds()
  c8_list <- GSEABase::getGmt('/path/to/c8.all.v2022.1.Hs.symbols.gmt') %>% GSEABase::geneIds()
  
  gobp_names <- names(c5_list)[grepl('GOBP',names(c5_list))]
  gobp_list <- c5_list[gobp_names]
  gocc_names <- names(c5_list)[grepl('GOCC',names(c5_list))]
  gocc_list <- c5_list[gocc_names]
  kegg_names <- names(c2_list)[grepl('KEGG',names(c2_list))]
  kegg_list <- c2_list[kegg_names]
  gse22886_names <- names(c7_list)[grepl("GSE22886", names(c7_list))]
  gse22886_list <- c7_list[gse22886_names]
  
  #oc_gene_list <- c5_list["GOBP_ONE_CARBON_METABOLIC_PROCESS"] %>% unlist()
  oc_gene_list <- c('PHGDH','PSAT1','PSPH','SHMT1','SHMT2','MTHFD1','MTHFD2','MTHFD1L','MTHFD2L','TYMS',
                    'MTHFR','ALDH1L1','ALDH1L2','SFXN1','SFXN2','SFXN3')
  tds_genes = c('DIO1','DIO2','DUOX1','DUOX2','FOXE1','GLIS3','NKX2-1','PAX8','SLC26A4','SLC5A5','SLC5A8','TG','THRA','THRB','TPO','TSHR')
  
  fibro_endo_thyro_genes = c('DCN','COL1A1','ACTA2','COL3A1','CDH5','PECAM1','TG','PAX8','TPO')
  immune_genes=c('PTPRC','CD68','CD14','CD4','CD40LG','FOXP3','CD8A','NKG7','CD79A','CD79B','IGHG1','TPSAB1')  
  
}



#preprocessing
if(T){
  #load data
  PT3_dt <- readRDS('/path/to/PT3_RDS/tc_prefi_ct.rds')
  PT5_dt <- readRDS('/path/to/PT5_RDS/tc_prefi_ct.rds')
  PT7_dt <- readRDS('/path/to/PT7_RDS/tc_prefi_ct.rds')
  PT8_dt <- readRDS('/path/to/PT8_RDS/tc_prefi_ct.rds')
  PT9_dt <- readRDS('/path/to/PT9_RDS/tc_prefi_ct.rds')
  PT10_dt <- readRDS('/path/to/PT10_RDS/tc_prefi_ct.rds')
  PT12_dt <- readRDS('/path/to/PT12_RDS/tc_prefi_ct.rds')
  
  AT9_dt <-readRDS('/path/to/AT9_RDS/tc_prefi_ct.rds')
  AT13_dt <-readRDS('/path/to/AT13_RDS/tc_prefi_ct.rds')
  AT16_dt <-readRDS('/path/to/AT16_RDS/tc_prefi_ct.rds')
  AT17_dt <-readRDS('/path/to/AT17_RDS/tc_prefi_ct.rds')
  AT20_dt <-readRDS('/path/to/AT20_RDS/tc_prefi_ct.rds')
  
  NT3_dt <- readRDS('/path/to/NT3_RDS/tc_prefi_ct.rds')
  NT5_dt <- readRDS('/path/to/NT5_RDS/tc_prefi_ct.rds')
  NT10_dt <- readRDS('/path/to/NT10_RDS/tc_prefi_ct.rds')

  #extract thyroid cells
  # no filtering anymore
  f_PT3 <- PT3_dt
  f_PT5 <- PT5_dt
  f_PT7 <- PT7_dt
  f_PT8 <- PT8_dt
  f_PT9 <- PT9_dt
  f_PT10 <- PT10_dt
  f_PT12 <- PT12_dt
  
  f_AT9 <- AT9_dt
  f_AT13 <- AT13_dt
  f_AT16 <- AT16_dt
  f_AT17 <- AT17_dt
  f_AT20 <- AT20_dt
  
  f_NT3 <- NT3_dt
  f_NT5 <- NT5_dt
  f_NT10 <- NT10_dt
  
  #rename
  f_PT3 <- RenameCells(object = f_PT3, add.cell.id = "PT3")
  f_PT3$orig.ident <- 'PT3'
  f_PT5 <- RenameCells(object = f_PT5, add.cell.id = "PT5")
  f_PT5$orig.ident <- 'PT5'
  f_PT7 <- RenameCells(object = f_PT7, add.cell.id = "PT7")
  f_PT7$orig.ident <- 'PT7'
  f_PT8 <- RenameCells(object = f_PT8, add.cell.id = "PT8")
  f_PT8$orig.ident <- 'PT8'
  f_PT9 <- RenameCells(object = f_PT9, add.cell.id = "PT9")
  f_PT9$orig.ident <- 'PT9'
  f_PT10 <- RenameCells(object = f_PT10, add.cell.id = "PT10")
  f_PT10$orig.ident <- 'PT10'
  f_PT12 <- RenameCells(object = f_PT12, add.cell.id = "PT12")
  f_PT12$orig.ident <- 'PT12'
  
  f_AT9 <- RenameCells(object = f_AT9, add.cell.id = "AT9")
  f_AT9$orig.ident <- 'AT9'
  f_AT13 <- RenameCells(object = f_AT13, add.cell.id = "AT13")
  f_AT13$orig.ident <- 'AT13'
  f_AT16 <- RenameCells(object = f_AT16, add.cell.id = "AT16")
  f_AT16$orig.ident <- 'AT16'
  f_AT17 <- RenameCells(object = f_AT17, add.cell.id = "AT17")
  f_AT17$orig.ident <- 'AT17'
  f_AT20 <- RenameCells(object = f_AT20, add.cell.id = "AT20")
  f_AT20$orig.ident <- 'AT20'
  
  f_NT3 <- RenameCells(object = f_NT3, add.cell.id = "NT3")
  f_NT3$orig.ident <- 'NT3'
  f_NT5 <- RenameCells(object = f_NT5, add.cell.id = "NT5")
  f_NT5$orig.ident <- 'NT5'
  f_NT10 <- RenameCells(object = f_NT10, add.cell.id = "NT10")
  f_NT10$orig.ident <- 'NT10'
 
  f_AT9$nCount_RNA %>% length()
  f_AT13$nCount_RNA %>% length()
  f_AT16$nCount_RNA %>% length()
  f_AT17$nCount_RNA %>% length()
  f_AT20$nCount_RNA %>% length()
  
  f_NT1$nCount_RNA %>% length()
  
  # normalize and identify variable features for each dataset independently
  tc.list = list(f_PT3, f_PT5, f_PT7, f_PT8, f_PT9, f_PT10, f_PT12,
                 f_AT9, f_AT13, f_AT16, f_AT17, f_AT20, 
                 f_NT3,f_NT5,f_NT10)
  tc.list <- lapply(X = tc.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = tc.list)
  
  tc.anchors <- FindIntegrationAnchors(object.list = tc.list, anchor.features = features, k.filter=70)
  
  # this command creates an 'integrated' data assay
  tc.combined <- IntegrateData(anchorset = tc.anchors)
  
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(tc.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  tc.combined <- ScaleData(tc.combined, verbose = FALSE)
  tc.combined <- RunPCA(tc.combined, npcs = 300, verbose = FALSE)
  
  ElbowPlot(tc.combined , ndims=300)
  pc_selection=50
  
  
  #Cluster
  DefaultAssay(tc.combined) <- "integrated"
  tc.combined <- FindNeighbors(tc.combined, dims = 1:pc_selection)
  tc.combined <- FindClusters(tc.combined, resolution = 0.5)
  tc.combined <- FindClusters(tc.combined, resolution = 1)
  tc.combined <- FindClusters(tc.combined, resolution = 1.5)
  tc.combined <- FindClusters(tc.combined, resolution = 2)
  tc.combined <- FindClusters(tc.combined, resolution = 3)
  tc.combined <- FindClusters(tc.combined, resolution = 4)
  
  tc.combined <- RunUMAP(tc.combined, dims = 1:pc_selection)
  draw_umap_plot(tc.combined, 'integrated_snn_res.1.5')
  draw_umap_plot(tc.combined, 'integrated_snn_res.2')
  draw_umap_plot(tc.combined, 'integrated_snn_res.3')
  draw_umap_plot(tc.combined, 'integrated_snn_res.4')
}

if(T){
  # set default assay as RNA  
  DefaultAssay(tc.combined) <- "RNA"
  
  # filter out non follicular cell in NT samples
  draw_umap_plot(tc.combined, 'cell_type')
  tc.combined@meta.data %>% nrow()  #132490 
  tmp <- tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% as_tibble()
  rm_ids <- tmp %>% filter(grepl("NT", orig.ident) & cell_type == "unknown")  %>% pull(cell_id)
  f.tc.combined <- subset(tc.combined, cells = setdiff(tmp$cell_id, rm_ids))
  f.tc.combined@meta.data %>% nrow()  #109988
  draw_umap_plot(f.tc.combined, 'cell_type')
  
  #Draw scatter plot by cancer type
  tmp_meta <- f.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>%
    as_tibble()
  tmp_meta <- tmp_meta %>% rowwise() %>% mutate(cancer_type = ifelse(grepl('PT',orig.ident),'ptc',ifelse(grepl('AT',orig.ident),'atc','normal')))
  f.tc.combined@meta.data <- tmp_meta %>% as.data.frame() %>% column_to_rownames('cell_id')
  DimPlot(f.tc.combined, group.by="cancer_type", reduction='umap')
  
  #QC plots
  FeaturePlot(f.tc.combined, features=c('percent.mt', 'nFeature_RNA','MKI67','GAPDH'))
  FeaturePlot(f.tc.combined, feature=c("nCount_RNA"))
  
  #cell cycle plot
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  f.tc.combined <- CellCycleScoring(f.tc.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  DimPlot(f.tc.combined, group.by="Phase", reduction='umap')
  
  # Draw cluster plot
  draw_umap_plot(f.tc.combined, 'integrated_snn_res.1.5')
  #low-quality clusters: 8,9,12,20
  
  #Filter out low-qual clusters
  Idents(f.tc.combined) <- 'integrated_snn_res.1.5' 
  ff.tc.combined <- subset(f.tc.combined, idents = setdiff(unique(Idents(f.tc.combined)), c(8,9,12,20)))
  ff.tc.combined@meta.data %>% nrow()  #97156
  draw_umap_plot(ff.tc.combined, 'cell_type')
  
  #re clustering
  DefaultAssay(ff.tc.combined) <- "integrated"
  ff.tc.combined <- FindNeighbors(ff.tc.combined, dims = 1:pc_selection)
  ff.tc.combined <- FindClusters(ff.tc.combined, resolution = 1.5)
  ff.tc.combined <- FindClusters(ff.tc.combined, resolution = 2)
  
  ff.tc.combined <- RunUMAP(ff.tc.combined, dims = 1:pc_selection)
  draw_umap_plot(ff.tc.combined, 'integrated_snn_res.1.5')
  draw_umap_plot(ff.tc.combined, 'cell_type')
  draw_umap_plot(ff.tc.combined, 'orig.ident')
  FeaturePlot(ff.tc.combined, features=c('percent.mt', 'nFeature_RNA','MKI67','GAPDH'))
  FeaturePlot(ff.tc.combined, feature=c("nCount_RNA"))
  
  #cell cycle plot
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  ff.tc.combined <- CellCycleScoring(ff.tc.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  DimPlot(ff.tc.combined, group.by="Phase", reduction='umap')
  
  # Cell type re confirm
  # check marker genes
  DefaultAssay(ff.tc.combined) <- "RNA"
  FeaturePlot(ff.tc.combined, features=fibro_endo_thyro_genes, slot='data')
  FeaturePlot(ff.tc.combined, features=immune_genes, slot='data')
  FeaturePlot(tc.combined, tds_genes, slot='data')
  
#chekc DEG, GSEA to identify cell type
  t_obj <- ff.tc.combined
  t_ident <- 'integrated_snn_res.1.5'
  Idents(t_obj) <- t_ident
  #check: 20 vs others, 34 vs c(2,3,4,5,6,17,21,23,24,26,29,36) Tcell,45 vs others, 38 vs others, 38vs16
  # 38vs c(7,8,11,12,15,18,22)M@, 38 vs c(2,4) CD4Tcell, 45 vs M@, 41 vs M@, 39 vs (c25,31), 19vs 37
  # c(7,8,11,12,15,18,22)M@ vs 37, 
  g1=c(2,3,4,5,6,17,21,23,24,26,29,36)
  g2=NULL
  draw_umap_mark_plot(t_obj,t_ident, g1)
  draw_umap_mark_plot(t_obj,t_ident, g2)
 
  mk <- FindMarkers(t_obj, ident.1 = g1, ident.2=g2, min.pct = 0.25, max.cells.per.ident=2000, random.seed=1) 
  DefaultAssay(obj) <- "RNA"
  #umap plot for g1 high
  tmp_feature <- mk  %>%
    as.data.frame() %>% rownames_to_column("gene_name") %>% as_tibble() %>%
    filter(avg_logFC > 0) %>% arrange(p_val) %>% slice(1:9) %>% pull(gene_name)
  FeaturePlot(t_obj, features=tmp_feature, slot='data')
  #umap plot for g2 high
  tmp_feature <- mk  %>%
    as.data.frame() %>% rownames_to_column("gene_name") %>% as_tibble() %>%
    filter(avg_logFC < 0) %>% arrange(p_val) %>% slice(1:9) %>% pull(gene_name)
  FeaturePlot(t_obj, features=tmp_feature, slot='data')
  #volcano and gsea plot together
  draw_volcano_gsea_plot_seurat(mk, gobp_list, name1=g1, name2=g2, repel_overlap=20)
  draw_volcano_gsea_plot_seurat(mk, kegg_list, name1=g1, name2=g2)
  #check RNA amount difference
  compare_meta_two_groups(t_obj, t_ident, "nCount_RNA", g1, g2)
  compare_meta_two_groups(t_obj, t_ident, "nFeature_RNA", g1, g2)
  compare_meta_two_groups(t_obj, t_ident, "percent.mt", g1, g2)
  
  # Cell type assign
  
  "
Fibroblast: 14, 42
Fibro-myofibro: 48
Myofibroblast: 9,32,43
Endothelial : 13,28,30,35,51,54
Follicular: 0,1,19,37,40,47,49,50
Macrophage:7,8,11,12,15,18,22, 33,44
Treg: 3
CD4 Tcell: 2,4,23,29
CD8 Tcell: 5,6,17,21,36
ISG Tcell: 34
NK cell: 24,26
Bcell: 10
Plasmacell : 25,31
mastcell: 46
fibro_cell_cycle:48
immune_cell_cycle:16,27,38
low_qual: 20
AT20_specific: 45,53
"
  
  assign_cell_type <- function(numb){
    if (numb %in% c(9,32,43)){
      clust = 'myofibroblast'
    } else if (numb %in% c(14,42)){
      clust = 'fibroblast'
    }else if (numb %in% c(13,28,30,35,51,54)){
      clust = 'endothelial_cell'
    }else if (numb %in% c(0,1,19,47,49,50)){
      clust = 'thyroid_cell'
    } else if (numb %in% c(7,8,11,12,15,18,22, 33,44)){
      clust = 'macrophage'
    } else if (numb %in% c(55)){
      clust = "CCL17_macrophage"
    } else if (numb %in% c(24,26)){
      clust = 'NK_cell'
    } else if (numb %in% c(5,6,17,21,36)){
      clust = 'CD8_T_cell'
    } else if (numb %in% c( 2,4,23,29)){
      clust = 'CD4_T_cell'
    } else if (numb %in% c(3)){
      clust = 'Treg_cell'
    } else if (numb %in% c()){
        clust = 'naive_T_cell'
    } else if (numb %in% c()){
      clust = 'naive_T_cell'
    } else if (numb %in% c(10)){
      clust = 'B_cell'
    } else if (numb %in% c(25,31)){
      clust = 'plasma_cell'
    } else if (numb %in% c(46)){
      clust = 'mast_cell'
    } else if (numb %in% c(34)){
      clust = "ISG_T_cell"
    } else if (numb %in% c(41)){
      clust = "dendritic_cell"
    } else if (numb %in% c(48)){
      clust = "fibro_cell_cycle"
    } else if (numb %in% c(47)){
      clust = "thyroid_cell_cycle"
    } else if (numb %in% c(40)){
      clust = "thyroid_cell_cycle"
    }else if (numb %in% c(16,27)){
      clust = "immune_cell_cycle"
    } else if (numb %in% c(20, 53)){
      clust = "low_qual"
    } else if (numb %in% c(37,38,39,45)){
      clust="possible_doublet"
    } else if (numb %in% c(52)){
      clust="KIT_high_T_cell"
    }else{
      clust='unknown'
    }
    return(clust)
  }
  
  tmp_meta <- ff.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id')
  tmp_meta <- tmp_meta %>% rowwise() %>% mutate(cell_type = assign_cell_type(integrated_snn_res.1.5))
  ff.tc.combined@meta.data <- tmp_meta %>% as.data.frame() %>% column_to_rownames('cell_id')
  draw_umap_plot(ff.tc.combined, 'cell_type')
}


#additional filture using cluster
# --- filter: possible_double, cell_cycle-related, low_qual

#Filter out low-qual clusters
Idents(ff.tc.combined) <- 'cell_type' 
f3.tc.combined <- subset(ff.tc.combined, 
                          idents = setdiff(unique(Idents(ff.tc.combined)), 
                                           c("immune_cell_cycle","fibro_cell_cycle","thyroid_cell_cycle","possible_doublet","low_qual")))
f3.tc.combined@meta.data %>% nrow()  #97156 -> 88493
draw_umap_plot(f3.tc.combined, 'cell_type')

#re clustering wifh f3
DefaultAssay(f3.tc.combined) <- "integrated"
pc_selection=50
f3.tc.combined <- FindNeighbors(f3.tc.combined, dims = 1:pc_selection)
f3.tc.combined <- FindClusters(f3.tc.combined, resolution = 1.5)
f3.tc.combined <- FindClusters(f3.tc.combined, resolution = 2)

f3.tc.combined <- RunUMAP(f3.tc.combined, dims = 1:pc_selection)
draw_umap_plot(f3.tc.combined, 'integrated_snn_res.1.5')
draw_umap_plot(f3.tc.combined, 'cell_type')
draw_umap_plot(f3.tc.combined, 'orig.ident')
FeaturePlot(f3.tc.combined, features=c('percent.mt', 'nFeature_RNA','MKI67','GAPDH'))
FeaturePlot(f3.tc.combined, feature=c("nCount_RNA"))

#cell cycle plot
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
f3.tc.combined <- CellCycleScoring(f3.tc.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(f3.tc.combined, group.by="Phase", reduction='umap')

# Cell type re confirm
# check marker genes
DefaultAssay(f3.tc.combined) <- "RNA"
FeaturePlot(f3.tc.combined, features=fibro_endo_thyro_genes, slot='data')
FeaturePlot(f3.tc.combined, features=immune_genes, slot='data')
FeaturePlot(f3.tc.combined, tds_genes, slot='data')

#chekc DEG, GSEA to identify cell type
t_obj <- f3.tc.combined
t_ident <- 'integrated_snn_res.1.5'
Idents(t_obj) <- t_ident
#check 39 vs 0, 38 vs 0, 40 vs 0, 41 vs 15, 11 vs 15, 43 vs 26
g1=c(43)
g2=c(26)
draw_umap_mark_plot(t_obj,t_ident, g1)
draw_umap_mark_plot(t_obj,t_ident, g2)

mk <- FindMarkers(t_obj, ident.1 = g1, ident.2=g2, min.pct = 0.25, max.cells.per.ident=2000, random.seed=1) 
DefaultAssay(obj) <- "RNA"
#umap plot for g1 high
tmp_feature <- mk  %>%
  as.data.frame() %>% rownames_to_column("gene_name") %>% as_tibble() %>%
  filter(avg_logFC > 0) %>% arrange(p_val) %>% slice(1:9) %>% pull(gene_name)
FeaturePlot(t_obj, features=tmp_feature, slot='data')
#umap plot for g2 high
tmp_feature <- mk  %>%
  as.data.frame() %>% rownames_to_column("gene_name") %>% as_tibble() %>%
  filter(avg_logFC < 0) %>% arrange(p_val) %>% slice(1:9) %>% pull(gene_name)
FeaturePlot(t_obj, features=tmp_feature, slot='data')
#volcano and gsea plot together
draw_volcano_gsea_plot_seurat(mk, gobp_list, name1=g1, name2=g2, repel_overlap=50)
draw_volcano_gsea_plot_seurat(mk, kegg_list, name1=g1, name2=g2)
#check RNA amount difference
compare_meta_two_groups(t_obj, t_ident, "nCount_RNA", g1, g2)
compare_meta_two_groups(t_obj, t_ident, "nFeature_RNA", g1, g2)
compare_meta_two_groups(t_obj, t_ident, "percent.mt", g1, g2)

# cell type2 assign
#39: TPO high follicular cells in TME
#38: possible doublet between endo + thyroid
#40: possible doublet with myofibro 
#41: possible doublet (endo + myofibro)

tmp_meta <- f3.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id')
tmp_meta <- tmp_meta %>% rowwise() %>% mutate(cell_type2 = case_when(`integrated_snn_res.1.5` == 39 ~ "TPO_high_thyroid_cell",
                                                                     `integrated_snn_res.1.5` %in% c(38,40,41) ~ "possible_doublet",
                                                                     TRUE ~ as.character(`cell_type`)))
f3.tc.combined@meta.data <- tmp_meta %>% as.data.frame() %>% column_to_rownames('cell_id')
draw_umap_plot(f3.tc.combined, 'cell_type2')


#filter and make f4
#Filter out low-qual clusters
Idents(f3.tc.combined) <- 'cell_type2' 
f4.tc.combined <- subset(f3.tc.combined, 
                         idents = setdiff(unique(Idents(f3.tc.combined)), 
                                          c("possible_doublet")))
f4.tc.combined@meta.data %>% nrow()  #88493 -> 88030
draw_umap_plot(f4.tc.combined, 'cell_type2')

#re clustering wifh f4
DefaultAssay(f4.tc.combined) <- "integrated"
pc_selection=50
f4.tc.combined <- FindNeighbors(f4.tc.combined, dims = 1:pc_selection)
f4.tc.combined <- FindClusters(f4.tc.combined, resolution = 1.5)
f4.tc.combined <- FindClusters(f4.tc.combined, resolution = 2)

f4.tc.combined <- RunUMAP(f4.tc.combined, dims = 1:pc_selection)
draw_umap_plot(f4.tc.combined, 'integrated_snn_res.1.5')
draw_umap_plot(f4.tc.combined, 'cell_type2')
draw_umap_plot(f4.tc.combined, 'orig.ident')
DefaultAssay(f4.tc.combined) <- "RNA"
FeaturePlot(f4.tc.combined, features=c('percent.mt', 'nFeature_RNA','MKI67','GAPDH'))
FeaturePlot(f4.tc.combined, feature=c("nCount_RNA"))

#cell cycle plot
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
f4.tc.combined <- CellCycleScoring(f4.tc.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(f4.tc.combined, group.by="Phase", reduction='umap')


# Cell type re confirm
# check marker genes
DefaultAssay(f4.tc.combined) <- "RNA"
FeaturePlot(f4.tc.combined, features=fibro_endo_thyro_genes, slot='data')
FeaturePlot(f4.tc.combined, features=immune_genes, slot='data')
FeaturePlot(f4.tc.combined, features=c("KIT","TPO","CCL17","LAMP3"), slot='data')
FeaturePlot(f4.tc.combined, tds_genes, slot='data')
FeaturePlot(f4.tc.combined, oc_gene_list, slot='data')

#Remove NT cells in non follicular cells again
draw_umap_plot(f4.tc.combined, 'cell_type2')
f4.tc.combined@meta.data %>% nrow()  #88030
unique(f4.tc.combined@meta.data$cell_type2)
tmp <- f4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% as_tibble()
rm_ids <- tmp %>% filter(grepl("NT", orig.ident) & cell_type != "thyroid_cell")  %>% pull(cell_id)
length(rm_ids) #22
f4.tc.combined <- subset(f4.tc.combined, cells = setdiff(tmp$cell_id, rm_ids))
f4.tc.combined@meta.data %>% nrow()  #88008
draw_umap_plot(f4.tc.combined, 'cell_type2')

#cell type3 assign for publication figure
tmp_meta <- f4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id')
tmp_meta <- tmp_meta %>% mutate(cell_type3 = case_when(cell_type2 == "TPO_high_thyroid_cell" ~ "Thyroid cell",
                                                       cell_type2 == "thyroid_cell" ~ "Thyroid cell",
                                                       cell_type2 == "fibroblast" ~ "Fibroblast",
                                                       cell_type2 == "myofibroblast" ~ "Fibroblast",
                                                       cell_type2 == "macrophage" ~ "Macrophage",
                                                       cell_type2 == "CCL17_macrophage" ~ "Macrophage",
                                                       cell_type2 == "endothelial_cell" ~ "Endothelial cell",
                                                       cell_type2 == "dendritic_cell" ~"Dendritic cell",
                                                       cell_type2 == "plasma_cell" ~ "Plasma cell",
                                                       cell_type2 == "B_cell" ~ "B cell",
                                                       cell_type2 == "mast_cell" ~ "Mast cell",
                                                       cell_type2 == "Treg_cell" ~ "Regulatory T cell",
                                                       cell_type2 == "CD4_T_cell" ~ "CD4 T cell",
                                                       cell_type2 == "CD8_T_cell" ~ "CD8 T cell",
                                                       cell_type2 == "ISG_T_cell" ~ "IFN-stimulated T cell",
                                                       cell_type2 == "KIT_high_T_cell" ~ "CD8 T cell",
                                                       cell_type2 == "NK_cell" ~ "NK cell"))
f4.tc.combined@meta.data <- tmp_meta %>% as.data.frame() %>% column_to_rownames('cell_id')


#clean up done


#save
if(T){
  f4.tc.combined %>% saveRDS(paste0(out_dir,'/tc_combined_pc',pc_selection,'_f4.rds'))
}


