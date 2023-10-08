
library(tidyverse)
library(Seurat)
library(DoubletFinder)
library(ggrepel)
library(fgsea)
library(corrplot)
library(cowplot)



#gene_sets

if(T){
  h_list <- GSEABase::getGmt('/path/to/h.all.v7.5.1.symbols.gmt') %>% GSEABase::geneIds()
  c5_list <- GSEABase::getGmt('/path/to/c5.all.v7.4.symbols.gmt') %>% GSEABase::geneIds()
  c2_list <- GSEABase::getGmt('/path/to/c2.all.v7.4.symbols.gmt') %>% GSEABase::geneIds()
  
  gobp_names <- names(c5_list)[grepl('GOBP',names(c5_list))]
  gobp_list <- c5_list[gobp_names]
  kegg_names <- names(c2_list)[grepl('KEGG',names(c2_list))]
  kegg_list <- c2_list[kegg_names]
  
  #oc_gene_list <- c5_list["GOBP_ONE_CARBON_METABOLIC_PROCESS"] %>% unlist()
  oc_gene_list <- c('PHGDH','PSAT1','PSPH','SHMT1','SHMT2','MTHFD1','MTHFD2','MTHFD1L','MTHFD2L','TYMS',
                    'MTHFR','ALDH1L1','ALDH1L2','SFXN1','SFXN2','SFXN3')
  tds_genes = c('DIO1','DIO2','DUOX1','DUOX2','FOXE1','GLIS3','NKX2-1','PAX8','SLC26A4','SLC5A5','SLC5A8','TG','THRA','THRB','TPO','TSHR')
  
  fibro_endo_thyro_genes = c('DCN','COL1A1','ACTA2','COL3A1','CDH5','PECAM1','TG','PAX8','TPO')
  immune_genes=c('PTPRC','CD68','CD14','CD4','CD40LG','FOXP3','CD8A','NKG7','CD79A','CD79B','IGHG1','TPSAB1')  
  
}

#functions
draw_umap_plot <- function(tc, ident, title=NULL){
  Idents(tc) <- ident
  umap_dt <- tc@reductions$umap@cell.embeddings
  umap_dt <- umap_dt %>% as.data.frame %>% rownames_to_column('cell_id')
  dim(umap_dt)
  ident_dt <- Idents(tc) %>% as.data.frame %>% rownames_to_column('cell_id') 
  colnames(ident_dt) <- c('cell_id','cluster')
  umap_dt <- left_join(umap_dt, ident_dt)
  #max_clust <- umap_dt$cluster %>% as.character() %>% as.numeric() %>% max()
  #umap_dt$cluster <- factor(umap_dt$cluster, levels = 0:max_clust)
  umap_dt <- umap_dt %>% separate(cell_id, c('barcode','orig.ident'),sep='-', remove=F) %>% as_tibble()
  umap_pos <- umap_dt %>% group_by(cluster) %>% summarise(med1=median(UMAP_1), med2=median(UMAP_2))
  g <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(aes(color=cluster), alpha=0.5, size=0.5)+
    geom_text_repel(data=umap_pos, aes(x=med1, y=med2, label=cluster), size=5)+
    guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
    theme_bw()+theme(panel.grid = element_blank())+
    ggtitle(title)
  print(g)
}

make_short_name <- function(x, max_n){  #x= string, n=maximum character
    if (nchar(x) < max_n){
      return(x)
    } else{
      return(substr(x, 1,max_n))
    }
  }

#raw data processing
if(T){
  out_dir="/path/to/NT5_RDS"
  system(paste0("mkdir -p ",out_dir))
  tc.data <- Read10X(data.dir = "/path/to/Thy05_GEX/outs/filtered_feature_bc_matrix/")
  tc <- CreateSeuratObject(counts = tc.data, project='thyroid_ca_e',min.cells = 3, min.features=200)
  tc[["percent.mt"]] <- PercentageFeatureSet(tc, pattern = "^MT-")
  
  
  
  if(F){
    #data QC
    VlnPlot(tc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    plot1 <- FeatureScatter(tc, feature1 = "nFeature_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(tc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot1 + plot2
  }
  
  #cell filter
  tc <- subset(tc, subset = nFeature_RNA > 200 & percent.mt < 20)
  
  if(T){ #conventional normalization and scailing
    #normalize
    tc <- NormalizeData(tc)
    #identificaiton of highly variable features
    tc <- FindVariableFeatures(tc, selection.method = "vst", nfeatures = 2000)
    
    if(F){
      # Identify the 10 most highly variable genes
      top10 <- head(VariableFeatures(tc), 10)
      # plot variable features with and without labels
      plot1 <- VariableFeaturePlot(tc)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
      plot2
    }
    
    #scailing the data
    all.genes <- rownames(tc)
    tc <- ScaleData(tc, features = all.genes)
  }

  #PCA
  tc <- RunPCA(tc, features = VariableFeatures(object = tc), npcs=100)
  
  if(T){
    #VizDimLoadings(tc, dims = 1:2, reduction = "pca")
    #DimPlot(tc, reduction='pca')
    
    DimHeatmap(tc, dims = 1:15, cells = 500, balanced = TRUE)
    DimHeatmap(tc, dims = 16:30, cells = 500, balanced = TRUE)
    DimHeatmap(tc, dims = 31:45, cells = 500, balanced = TRUE)
    DimHeatmap(tc, dims = 46:60, cells = 500, balanced = TRUE)
    
    ElbowPlot(tc, ndims=100)
  }
  
  pc_selection=30
  
  
  #Cluster
  tc <- FindNeighbors(tc, dims = 1:pc_selection)
  tc <- FindClusters(tc, resolution = 0.5)
  tc <- FindClusters(tc, resolution = 1)
  tc <- FindClusters(tc, resolution = 1.5)
  
  #UMAP
  tc <- RunUMAP(tc, dims = 1:pc_selection)
  if(F){
    
    draw_umap_plot(tc,'RNA_snn_res.0.5' ,title='res0.5' )
    draw_umap_plot(tc,'RNA_snn_res.1' , title='res1')
    draw_umap_plot(tc,'RNA_snn_res.1.5', title='res1.5')

    
    Idents(tc) <- 'RNA_snn_res.1.5'  #final resolution selection
   
    FeaturePlot(tc, features=c('percent.mt', 'nFeature_RNA','MKI67','GAPDH'))
    
    tmp_meta <- tc@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% as_tibble()
    tmp1 <- tmp_meta %>% group_by(`RNA_snn_res.1.5`) %>% summarise(mean_feature = mean(nFeature_RNA))
    tmp2 <- tmp_meta %>% group_by(`RNA_snn_res.1.5`) %>% summarise(mean_mito = mean(percent.mt))
    
    g1 <- ggplot(tmp2, aes(x='',y=mean_mito))+
      geom_boxplot()+
      geom_point(color='red')+
      geom_text_repel(aes(label=`RNA_snn_res.1.5`), size=7, max.overlaps=30)+
      xlab('')
    
    g2 <- ggplot(tmp1, aes(x='',y=mean_feature))+
      geom_boxplot()+
      geom_point(color='red')+
      geom_text_repel(aes(label=`RNA_snn_res.1.5`), size=7, max.overlaps=30)+
      xlab('')
    
    plot_grid(g1,g2,nrow=1)
    
  }
  
  #Cluster identification
  if(F){
   
    DefaultAssay(tc) <- "RNA"
    FeaturePlot(tc, features=fibro_endo_thyro_genes, slot='data')
    FeaturePlot(tc, features=immune_genes, slot='data')
    FeaturePlot(tc, tds_genes, slot='data')
    FeaturePlot(tc, oc_gene_list, slot='data')
    
  
  }
 
  
  if(T){
    #DoubleFinder
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    print(out_dir)
    if(T){
      sweep.res.list_tc <- paramSweep_v3(tc, PCs = 1:pc_selection, sct = FALSE)
    }
    sweep.stats_tc <- summarizeSweep(sweep.res.list_tc, GT = FALSE)
    bcmvn_tc <- find.pK(sweep.stats_tc)
    dev.off()
    tmp_dt <- bcmvn_tc %>% as_tibble()
    ggplot(tmp_dt, aes(x=pK, y=BCmetric))+
      geom_point()
    pK_selection=0.4
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- tc@meta.data$RNA_snn_res.1
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(tc@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    tc <- doubletFinder_v3(tc, PCs = 1:pc_selection, pN = 0.25, pK = pK_selection, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
    cols <- tc@meta.data %>% colnames()
    DF_name <- cols[grepl('DF',cols)]
    Idents(tc) <- DF_name
    DimPlot(tc, reduction = "umap")
  }
}

#additional cluster define

g1=c(12)
g2=c()

Idents(tc) <- 'RNA_snn_res.1.5' 
mk <- FindMarkers(tc, ident.1 = g1, ident.2=g2, min.pct = 0.25, max.cells.per.ident=2000, random.seed=1)
tmp_feature = rownames(mk)[1:9]
tmp_feature2 = rownames(mk)[10:18]
FeaturePlot(tc, features=tmp_feature, slot='data')
FeaturePlot(tc, features=tmp_feature2, slot='data')
if(T){
  #draw volcano
  tmp_dt <- mk %>% as.data.frame() %>% rownames_to_column('gene_name')
  gg1<- ggplot(tmp_dt, aes(x=avg_logFC,y=-log10(p_val_adj)))+
    geom_point()+
    geom_text_repel(aes(label=gene_name))+
    xlab(paste0(paste(g2,collapse=','),'<----->',paste(g1, collapse=',')))
  #GSEA
  tmp_dt <- tmp_dt %>% mutate(nlpval = -log10(p_val_adj))
  tmp_dt$nlpval[is.finite(tmp_dt$nlpval) == F] <- 300
  tmp_dt <- tmp_dt %>% mutate(dist= sqrt(avg_logFC^2 + nlpval ^2))
  tmp_dt <-tmp_dt %>% mutate(dist = ifelse(avg_logFC < 0, (-1)*dist, dist))
  dds_stats <- tmp_dt$dist
  names(dds_stats) <- tmp_dt$gene_name
  gsea_res <- fgsea(pathways = gobp_list, stats = dds_stats, nper=1000)
  #gsea_res <- fgsea(pathways = h_list, stats = dds_stats, nper=1000)
  #gsea_res <- fgsea(pathways = kegg_list, stats = dds_stats, nper=1000)
  m_res <- gsea_res %>% as_tibble()
  f_res <- m_res %>% filter(pval < 0.05) 
  f_res$rep_genes <-map_chr(f_res$leadingEdge, function(x) paste(x[1:5], collapse=','))
  f_res$short_name <- make.unique(map_chr(f_res$pathway, make_short_name, 50))
  x1 <- f_res %>% arrange(NES) %>% head(n=20) %>% pull(short_name)
  x2 <- f_res %>% arrange(NES) %>% tail(n=20) %>% pull(short_name)
  x_order = unique(c(x1,x2))
  gg2 <- ggplot(f_res, aes(x=short_name, y=NES, fill= pval))+
    geom_bar(stat='identity')+
    geom_text(aes(y=5, label=rep_genes), hjust=0)+
    scale_y_continuous(limits=c(-5,15))+
    scale_x_discrete(limits = x_order)+
    coord_flip()+
    ylab(paste0(paste(g2,collapse=','),'<----->',paste(g1, collapse=',')))
  plot_grid(gg1,gg2,nrow=1, rel_widths=c(1,1.5))
  
}

#assign cell type
draw_umap_plot(tc, "RNA_snn_res.1.5")


assign_cell_type <- function(numb){
  if (numb %in% c()){
    clust = 'myofibroblast'
  } else if (numb %in% c()){
    clust = 'fibroblast'
  }else if (numb %in% c()){
    clust = 'fibro_myofibro'
  }else if (numb %in% c()){
    clust = 'endothelial_cell'
  } else if (numb %in% c()){
    clust= 'atc'
  }else if (numb %in% c(4)){
    clust = 'normal_thyroid'
  } else if (numb %in% c()){
    clust = 'macrophage'
  } else if (numb %in% c()){
    clust = 'NK_cell'
  }else if (numb %in% c()){
    clust = 'unknown_T_cell'
  } else if (numb %in% c()){
    clust = 'CD8_T_cell'
  } else if (numb %in% c()){
    clust = 'CD4_T_cell'
  } else if (numb %in% c()){
    clust = 'Treg_cell'
  } else if (numb %in% c()){
    clust = 'naive_T_cell'
  } else if (numb %in% c()){
    clust = 'B_cell'
  } else if (numb %in% c()){
    clust = 'plasma_cell'
  } else if (numb %in% c()){
    clust = 'mast_cell'
  } else if (numb %in% c()){
    clust = 'low_qual'
  } else if (numb %in% c()){
    clust = 'doublet' 
  } else{
    clust='unknown'
  }
  return(clust)
}

tmp_meta <- tc@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id')
tmp_meta <- tmp_meta %>% rowwise() %>% mutate(cell_type = assign_cell_type(RNA_snn_res.1.5))
tc@meta.data <- tmp_meta %>% as.data.frame() %>% column_to_rownames('cell_id')
draw_umap_plot(tc, 'cell_type')

if(T){
  print(out_dir)
  tc %>% saveRDS(paste0(out_dir,"/tc_prefi_ct.rds"))  
}
