#DEG and GSEA per cell type after combat-seq 

#load libraries
if(T){
  library(tidyverse)
  library(Seurat)
  library(ggsci)
  library(scales)
  library(ggsignif)
  library(corrplot)
  library(cowplot)
  library(ggrepel)
  library(ComplexHeatmap)
  library(DESeq2)
  library(fgsea)
}

#output path-----
out_dir = '/path/to/output'
system(paste0('mkdir -p ',out_dir))

#design
theme_syp = theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=12), 
        panel.border = element_blank(), 
        axis.line = element_line(), 
        plot.title = element_text(size=14, face="bold"))

cancer_type_pal<- pal_d3("category10")(10)[c(1,3,2)]
names(cancer_type_pal) <- c('NT','PTC','ATC')

#script specific functions-----
draw_hm_thy <- function(mx, seurat_obj){
  f_mx <- mx
  f4.tc.combined <- seurat_obj
  #column annotation
  annot_cols <- tibble(cell_id = colnames(f_mx))
  annot_cols <- annot_cols %>% rowwise() %>% mutate(orig_ident = unlist(strsplit(cell_id,'_'))[1])
  tmp <- f4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column("cell_id") %>%
    as_tibble() %>% select(cell_id, nCount_RNA, percent.mt, cancer_type)
  annot_cols <- left_join(annot_cols, tmp)
  annot_cols <- annot_cols %>% column_to_rownames('cell_id')
  all(rownames(annot_cols) == colnames(f_mx))
  #row annotation
  annot_rows <- tibble(gene_name = rownames(f_mx))
  annot_rows <- annot_rows %>% mutate(group = case_when(grepl("^MT-", gene_name) ~ "MT",
                                                        grepl("^RPL", gene_name) ~ "RPL",
                                                        grepl("^HLA-",gene_name) | grepl("^HSPA",gene_name) ~ "AP",
                                                        TRUE ~ "other"))
  annot_rows <- annot_rows %>% column_to_rownames('gene_name')
  all(rownames(annot_rows) == rownames(f_mx))
  # color for annotation
  orig_pal = pal_igv('default')(length(unique(annot_cols$orig_ident)))
  names(orig_pal) <- unique(annot_cols$orig_ident)
  top_anno <- ComplexHeatmap::HeatmapAnnotation(df = annot_cols, col=list(orig_ident = orig_pal))
  row_anno <- ComplexHeatmap::rowAnnotation(df = annot_rows)
  #rowwise scaling
  if(T){
    scaled_dt <- t(apply(f_mx, 1, scale))
    colnames(scaled_dt) <- colnames(f_mx)
  }
  
  #draw plot
  chm <- ComplexHeatmap::Heatmap(scaled_dt,
                                 top_annotation = top_anno,
                                 right_annotation = row_anno,
                                 show_row_names = T, 
                                 show_column_names = F, 
                                 clustering_distance_columns = "pearson", 
                                 clustering_method_columns = 'ward.D2',
                                 clustering_method_rows = 'ward.D2', 
                                 cluster_columns = T,
                                 row_names_gp = grid::gpar(fontsize = 8)
  )
  return(chm)
}

make_short_name <- function(x, max_n){  #x= string, n=maximum character
  if (nchar(x) < max_n){
    return(x)
  } else{
    return(substr(x, 1,max_n))
  }
}


#load Seurat object------

f4.tc.combined <- readRDS("/path/to/Seurat_obj.rds")

# assign cell type2  -----
tmp_meta <- f4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id')
tmp_meta <- tmp_meta %>% mutate(cell_type4 = case_when(cell_type2 == "TPO_high_thyroid_cell" ~ "Thyroid cell",
                                                       cell_type2 == "thyroid_cell" ~ "Thyroid cell",
                                                       cell_type2 == "fibroblast" ~ "Fibroblast",
                                                       cell_type2 == "myofibroblast" ~ "Fibroblast",
                                                       cell_type2 == "macrophage" ~ "Macrophage",
                                                       cell_type2 == "CCL17_macrophage" ~ "Macrophage",
                                                       cell_type2 == "endothelial_cell" ~ "Endothelial cell",
                                                       cell_type2 == "dendritic_cell" ~"Dendritic cell",
                                                       cell_type2 %in% c("B_cell", "plasma_cell") ~ "B/Plasma cell",
                                                       cell_type2 == "mast_cell" ~ "Mast cell",
                                                       cell_type2 == "Treg_cell" ~ "T/NK cell",
                                                       cell_type2 == "CD4_T_cell" ~ "T/NK cell",
                                                       cell_type2 == "CD8_T_cell" ~ "T/NK cell",
                                                       cell_type2 == "ISG_T_cell" ~ "T/NK cell",
                                                       cell_type2 == "KIT_high_T_cell" ~ "T/NK cell",
                                                       cell_type2 == "NK_cell" ~ "T/NK cell"))
f4.tc.combined@meta.data <- tmp_meta %>% as.data.frame() %>% column_to_rownames('cell_id')


#do random sampling,batch correction, and draw heatmap for each cell type------
##for loop for celltype already done and rds saved-----
cell_type_list <- unique(f4.tc.combined@meta.data$cell_type4)
#exclude mast cell, due to lack of sufficient cells
cell_type_list <- setdiff(cell_type_list,c("Mast cell"))
for (t_cell_type in cell_type_list){
  print(t_cell_type)
  t_cell_type2= tolower(gsub("/","_",gsub(" ","_", t_cell_type)))

  #filter in only assigned cell type
  Idents(f4.tc.combined) <- "cell_type4"
  thy.combined <- subset(f4.tc.combined, idents = t_cell_type)
  
  if(file.exists(paste0(out_dir,'/r_merged_mx.',t_cell_type2,'.rds'))==F){
    #generate individual sample objects
    thy.combined@meta.data$orig.ident %>% unique()
    Idents(thy.combined) <- "orig.ident"
    fr_PT3_dt <- subset(thy.combined, idents="PT3")
    fr_PT5_dt <- subset(thy.combined, idents="PT5")
    fr_PT7_dt <- subset(thy.combined, idents="PT7")
    fr_PT8_dt <- subset(thy.combined, idents="PT8")
    fr_PT9_dt <- subset(thy.combined, idents="PT9")
    fr_PT10_dt <- subset(thy.combined, idents="PT10")
    fr_PT12_dt <- subset(thy.combined, idents="PT12")
    
    fr_AT9_dt <- subset(thy.combined, idents="AT9")
    fr_AT13_dt <- subset(thy.combined, idents="AT13")
    fr_AT16_dt <- subset(thy.combined, idents="AT16")
    fr_AT17_dt <- subset(thy.combined, idents="AT17")
    fr_AT20_dt <- subset(thy.combined, idents="AT20")
    
    #generate individual sample matrix
    pt3_mx <- fr_PT3_dt@assays$RNA@counts %>% as.matrix()
    pt5_mx <- fr_PT5_dt@assays$RNA@counts %>% as.matrix()
    pt7_mx <- fr_PT7_dt@assays$RNA@counts %>% as.matrix()
    pt8_mx <- fr_PT8_dt@assays$RNA@counts %>% as.matrix()
    pt9_mx <- fr_PT9_dt@assays$RNA@counts %>% as.matrix()
    pt10_mx <- fr_PT10_dt@assays$RNA@counts %>% as.matrix()
    pt12_mx <- fr_PT12_dt@assays$RNA@counts %>% as.matrix()
    
    at9_mx <- fr_AT9_dt@assays$RNA@counts %>% as.matrix()
    at13_mx <- fr_AT13_dt@assays$RNA@counts %>% as.matrix()
    at16_mx <- fr_AT16_dt@assays$RNA@counts %>% as.matrix()
    at17_mx <- fr_AT17_dt@assays$RNA@counts %>% as.matrix()
    at20_mx <- fr_AT20_dt@assays$RNA@counts %>% as.matrix()
    
    #generate merged matrix  (do not include normal thyroid cells)
    common_genes <- Reduce(intersect, list(rownames(pt3_mx), rownames(pt5_mx), 
                                           rownames(pt7_mx), 
                                           rownames(pt8_mx), rownames(pt9_mx), 
                                           rownames(pt10_mx), rownames(pt12_mx),
                                           rownames(at9_mx),
                                           rownames(at13_mx), rownames(at16_mx),
                                           rownames(at17_mx), rownames(at20_mx)
    ))
    merged_mx <- cbind(pt3_mx[common_genes,,drop=F],pt5_mx[common_genes,,drop=F], 
                       pt7_mx[common_genes,,drop=F], pt8_mx[common_genes,,drop=F],
                       pt9_mx[common_genes,,drop=F], 
                       pt10_mx[common_genes,,drop=F], pt12_mx[common_genes,,drop=F],
                       at9_mx[common_genes,,drop=F], at13_mx[common_genes,,drop=F],
                       at16_mx[common_genes,,drop=F],
                       at17_mx[common_genes,,drop=F],at20_mx[common_genes,,drop=F]
    )
    dim(merged_mx) 

    #remove cells with readcount < 10000
    low_rc_cells <- colnames(merged_mx)[colSums(merged_mx) < 10000]
    f_merged_mx <- merged_mx[,!colnames(merged_mx) %in% low_rc_cells]
    dim(f_merged_mx)
    
    #remove cells with mt.prop > 10%
    hi_mt_cells <- thy.combined@meta.data %>% as.data.frame() %>% 
      rownames_to_column("cell_id") %>%
      as_tibble() %>% filter(percent.mt > 10) %>% pull(cell_id)
    f_merged_mx <- f_merged_mx[,!colnames(f_merged_mx) %in% hi_mt_cells]
    dim(f_merged_mx)
    
    ##random sampling merged_mx to balanced comoparison-----
    #sample proportion adjust before random sampling (1st sampling)
    orig_f <- map_chr(colnames(f_merged_mx), function(x) unlist(strsplit(x,'_'))[1])
    tmp <- tibble(orig_name=names(table(orig_f)), count=as.numeric(table(orig_f)))
    tmp <- tmp %>% mutate(prop = count/sum(count))
    total_cell=ncol(f_merged_mx)
    tmp <- tmp %>% mutate(sample_count = case_when(prop < 1/nrow(tmp) ~ count,
                                                   prop >= 1/nrow(tmp) ~ round(total_cell/nrow(tmp))))
    selected_ids <- c()
    for (t_orig in tmp$orig_name){
      print(t_orig)
      set.seed(1);t_ids <- sample(colnames(f_merged_mx)[which(orig_f==t_orig)], tmp$sample_count[tmp$orig_name==t_orig])
      selected_ids <- c(selected_ids, t_ids)
    }
    ff_merged_mx <- f_merged_mx[,selected_ids]
    dim(ff_merged_mx)
    
    #random sampling up to 500 cells
    sample_n=500
    if(ncol(ff_merged_mx) > sample_n ){
      set.seed(1);selected_columns <- sample(ncol(ff_merged_mx),sample_n)
      r_merged_mx <- ff_merged_mx[,selected_columns]
    }else{
      r_merged_mx <- ff_merged_mx
    }
    dim(r_merged_mx)
    r_merged_mx %>% saveRDS(paste0(out_dir,'/r_merged_mx.',t_cell_type2,'.rds'))
  }else {
    r_merged_mx <- readRDS(paste0(out_dir,'/r_merged_mx.',t_cell_type2,'.rds'))
  }  
 
  if(file.exists(paste0(out_dir,'/fsr_bc_mx.',t_cell_type2,'.rds'))==F){
    #remove orig if it has only single sample since this raise error in combat
    orig <- map_chr(colnames(r_merged_mx), function(x) unlist(strsplit(x,'_'))[1])
    tmp_ids <- dimnames(table(orig))$orig[table(orig) <2]
    sr_merged_mx <- r_merged_mx[,!orig %in% tmp_ids] 
  
    #remove frequently zero count genes
    zero_count <- apply(sr_merged_mx, 1, function(x) sum(x==0))
    tmp <- tibble(gene_name=names(zero_count), zero_count=zero_count, zero_prop=zero_count/ncol(sr_merged_mx))
    zero_genes <- tmp %>% filter(zero_prop > 0.6) %>% pull(gene_name)
    fsr_merged_mx <- sr_merged_mx[!rownames(sr_merged_mx) %in% zero_genes,] 
  
    ##batch correction using ComBat_seqs-----
      batch <- map_chr(colnames(fsr_merged_mx), function(x) unlist(strsplit(x,'_'))[1])
      fsr_bc_mx <- sva::ComBat_seq(fsr_merged_mx, batch=batch, group=NULL)
      fsr_bc_mx %>% saveRDS(paste0(out_dir,'/fsr_bc_mx.',t_cell_type2,'.rds'))
  } else{
    fsr_bc_mx <- readRDS(paste0(out_dir,'/fsr_bc_mx.',t_cell_type2,'.rds'))
  }
  
  #draw heatmaps
  ##heamtap before batch correction using cpm_mx  ------
  if(file.exists(paste0(out_dir,'/heatmap_before_bc.',t_cell_type2,'.png'))==F){ 
    # convert into CPM matrix
    cpm_mx <- apply(fsr_merged_mx, 2, function(x) (x/sum(x))*1000000)
    
    #do not select genes for heatmap

    chm <- draw_hm_thy(cpm_mx, f4.tc.combined)
  
    png(paste0(out_dir,'/heatmap_before_bc.',t_cell_type2,'.png'),
        width=12, height=12, units="in", res=300)
    ComplexHeatmap::draw(chm, heatmap_legend_side = "left", annotation_legend_side = "left")
    dev.off()
  }
  
  ##heamtap after batch correction using cpm_mx------
  if(file.exists(paste0(out_dir,'/heatmap_after_bc.',t_cell_type2,'.png'))==F){ 
    #change bc_mx to CPM mx
    cpm_bc_mx <- apply(fsr_bc_mx, 2, function(x) (x/sum(x))*1000000)
    
    #do not select genes for heatmap

    chm <- draw_hm_thy(cpm_bc_mx, f4.tc.combined)
    chm %>% saveRDS(paste0(out_dir,'/chm_cpm_bc.',t_cell_type2,'.rds'))
    png(paste0(out_dir,'/heatmap_after_bc.',t_cell_type2,'.png'),
        width=12, height=12, units="in", res=300)
    ComplexHeatmap::draw(chm, heatmap_legend_side = "left", annotation_legend_side = "left")
    dev.off()
  }
}


#RPL GENE AND MT GENE FILTER to make balanced quality -----
if(T){
  #Manually define cluster number which is associated with RPL genes
  rpl_cl_tbl <- tribble(~t_cell_type2, ~cl_num,
                        "thyroid_cell", 5,
                        "t_nk_cell",9,
                        "b_plasma_cell",1,
                        "macrophage",5,
                        "fibroblast",8,
                        "endothelial_cell",6,
                        "dendritic_cell",1
  )
  #mast cells, b_cell will be excluded due to lack of high-quality cells
  
  #heatmap split to select quality-related genes
  for (t_cell_type in setdiff(unique(f4.tc.combined$cell_type4),c("Mast cell"))){
    t_cell_type2= tolower(gsub("/","_",gsub(" ","_", t_cell_type)))
    print(t_cell_type2)
    if(file.exists(paste0(out_dir,'/heatmap_after_bc_deg_1en6.',t_cell_type2,'.pdf'))==F){
      chm <- readRDS(paste0(out_dir,'/chm_cpm_bc.',t_cell_type2,'.rds'))
      ro <- row_order(chm)
      rpl_cl_num <- rpl_cl_tbl$cl_num[rpl_cl_tbl$t_cell_type2==t_cell_type2]
      rpl_genes <- rownames(chm@matrix)[ro[[rpl_cl_num]]]
      tibble(rpl_genes = rpl_genes) %>% 
        write_tsv(paste0(out_dir,'/rpl_genes.',t_cell_type2,'.tsv'))
      r_bc_mx <- readRDS(paste0(out_dir,'/fsr_bc_mx.',t_cell_type2,'.rds'))
      #remove rpl related genes and ^MT- genes
      f_bc_mx <- r_bc_mx[!(rownames(r_bc_mx) %in% rpl_genes | grepl("^MT-", rownames(r_bc_mx))),]
      
      ##deg and gsea using cpm_bc_mx -----
      #run DEG
      dim(f_bc_mx)
      ff_bc_mx <- f_bc_mx[rowSums(f_bc_mx)>0,] 
      dim(ff_bc_mx)
      tmp <- f4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column("cell_id") %>%
        as_tibble() %>% select(cell_id, cancer_type)
      pheno_dt<- tibble(cell_id=colnames(ff_bc_mx))
      pheno_dt <- left_join(pheno_dt, tmp) %>% as.data.frame() %>% 
        column_to_rownames("cell_id")
      all(colnames(ff_bc_mx) == rownames(pheno_dt))
      #add pseudocount 1 to all data
      mff_bc_mx <- ff_bc_mx+1
      dds <- DESeqDataSetFromMatrix(countData= mff_bc_mx, 
                                    colData = pheno_dt,
                                    design = ~ cancer_type)
      dds <- DESeq(dds)
      res <- results(dds, contrast = c("cancer_type","atc","ptc"))
      mcols(res, use.names=TRUE)
      
      #volcano plot
      res_dt <- res %>% as.data.frame() %>% rownames_to_column("gene_name") %>% as_tibble()
      res_dt %>% write_tsv(paste0(out_dir,'/deg.',t_cell_type2,'.tsv'))
      label_dt <- res_dt %>% filter(-log10(padj)> 10)
      g <- ggplot(res_dt, aes(x=log2FoldChange,y=-log10(padj)))+
        geom_point()+
        geom_text_repel(data=label_dt, aes(label=gene_name), max.overlaps = 100)+
        theme_syp
      png(paste0(out_dir,'/volcano.',t_cell_type2,'.png'), width=8, height=6,units="in", res=300)
      print(g)
      dev.off()
      
      #heatmap with deg p< 1e-3
      degs <- res_dt %>% filter(padj < 1e-3) %>% pull(gene_name)
      #change bc_mx to CPM mx
      cpm_bc_mx <- apply(r_bc_mx, 2, function(x) (x/sum(x))*1000000)
      # selected only DEGs
      f_mx <- cpm_bc_mx[degs,]
   
      chm <- draw_hm_thy(f_mx, f4.tc.combined)
      chm %>% saveRDS(paste0(out_dir,'/chm_deg_bc_1en3.',t_cell_type2,'.rds'))
      png(paste0(out_dir,'/heatmap_after_bc_deg_1en3.',t_cell_type2,'.png'),
          width=12, height=12, units="in", res=300)
      ComplexHeatmap::draw(chm, heatmap_legend_side = "left", annotation_legend_side = "left")
      dev.off()
    
      pdf(paste0(out_dir,'/heatmap_after_bc_deg_1en3.',t_cell_type2,'.pdf'),
          width=12, height=12)
      ComplexHeatmap::draw(chm, heatmap_legend_side = "left", annotation_legend_side = "left")
      dev.off()
      
      #heatmap with deg padj < 1e-6
      degs <- res_dt %>% filter(padj < 1e-6) %>% pull(gene_name)
      #change bc_mx to CPM mx
      cpm_bc_mx <- apply(r_bc_mx, 2, function(x) (x/sum(x))*1000000)
      # selected only DEGs
      f_mx <- cpm_bc_mx[degs,]
      
      chm <- draw_hm_thy(f_mx, f4.tc.combined)
      chm %>% saveRDS(paste0(out_dir,'/chm_deg_bc_1en6.',t_cell_type2,'.rds'))
      png(paste0(out_dir,'/heatmap_after_bc_deg_1en6.',t_cell_type2,'.png'),
          width=12, height=12, units="in", res=300)
      ComplexHeatmap::draw(chm, heatmap_legend_side = "left", annotation_legend_side = "left")
      dev.off()
      
      pdf(paste0(out_dir,'/heatmap_after_bc_deg_1en6.',t_cell_type2,'.pdf'),
          width=12, height=12)
      ComplexHeatmap::draw(chm, heatmap_legend_side = "left", annotation_legend_side = "left")
      dev.off()
      
      #run GSEA
      h_list <- GSEABase::getGmt('/path/to/h.all.v7.5.1.symbols.gmt') %>% GSEABase::geneIds()
      c5_list <- GSEABase::getGmt('/path/to/c5.all.v7.4.symbols.gmt') %>% GSEABase::geneIds()
      c2_list <- GSEABase::getGmt('/path/to/c2.all.v7.4.symbols.gmt') %>% GSEABase::geneIds()
      
      gobp_names <- names(c5_list)[grepl('GOBP',names(c5_list))]
      gobp_list <- c5_list[gobp_names]
      kegg_names <- names(c2_list)[grepl('KEGG',names(c2_list))]
      kegg_list <- c2_list[kegg_names]
      
      #make distance value
      tmp_dt <-res_dt %>% mutate(nlpval = -log10(pvalue))
      tmp_dt <- tmp_dt %>% mutate(dist= sqrt(log2FoldChange^2 + nlpval ^2))
      tmp_dt <-tmp_dt %>% mutate(dist = ifelse(log2FoldChange < 0, (-1)*dist, dist))
      
      #run gsea with GOBP
      dds_stats <- tmp_dt$dist
      names(dds_stats) <- tmp_dt$gene_name
      gsea_res <- fgsea(pathways = gobp_list, stats = dds_stats, nper=1000)
      m_res <- gsea_res %>% as.tibble()
      m_res$rep_genes <-map_chr(m_res$leadingEdge, function(x) paste(x[1:length(x)], collapse=','))
      m_res %>% select(-leadingEdge) %>%  write_tsv(paste0(out_dir,'/gsea.gobp.',t_cell_type2,'.tsv'))

      f_res <- m_res %>% filter(pval < 0.05) 
      f_res$short_name <- make.unique(map_chr(f_res$pathway, make_short_name, 50))
      x1 <- f_res %>% filter(NES < 0) %>% arrange(NES) %>% head(n=20) %>% pull(short_name)
      x2 <- f_res %>% filter(NES > 0) %>% arrange(NES) %>% tail(n=20) %>% pull(short_name)
      x_order = unique(c(x1,x2))
      g<- ggplot(f_res, aes(x=short_name, y=NES, fill= pval))+
        geom_bar(stat='identity')+
        geom_text(aes(y=5, label=rep_genes), hjust=0)+
        scale_y_continuous(limits=c(-5,15))+
        scale_x_discrete(limits = x_order)+
        coord_flip()+
        theme_syp
      png(paste0(out_dir,'/gsea.gobp.',t_cell_type2,'.png'), width=12, height=6,units="in", res=300)
      print(g)
      dev.off()
      
      #run gsea with KEGG
      gsea_res <- fgsea(pathways = kegg_list, stats = dds_stats, nper=1000)
      m_res <- gsea_res %>% as.tibble()
      m_res$rep_genes <-map_chr(m_res$leadingEdge, function(x) paste(x[1:length(x)], collapse=','))
      m_res %>% select(-leadingEdge) %>%  write_tsv(paste0(out_dir,'/gsea.kegg.',t_cell_type2,'.tsv'))

      f_res <- m_res %>% filter(pval < 0.05) 
      f_res$short_name <- make.unique(map_chr(f_res$pathway, make_short_name, 50))
      x1 <- f_res %>% filter(NES < 0) %>% arrange(NES) %>% head(n=20) %>% pull(short_name)
      x2 <- f_res %>% filter(NES > 0) %>% arrange(NES) %>% tail(n=20) %>% pull(short_name)
      x_order = unique(c(x1,x2))
      g<- ggplot(f_res, aes(x=short_name, y=NES, fill= pval))+
        geom_bar(stat='identity')+
        geom_text(aes(y=5, label=rep_genes), hjust=0)+
        scale_y_continuous(limits=c(-5,15))+
        scale_x_discrete(limits = x_order)+
        coord_flip()+
        theme_syp
      png(paste0(out_dir,'/gsea.kegg.',t_cell_type2,'.png'), width=12, height=6,units="in", res=300)
      print(g)
      dev.off()
     }
  }
}


# make a merged GSEA figure-----
if(T){
  #draw GSEA plot per cell_type and put in the list
  #used ggplot2 default color gradient for 0.05 and 0.001
  g_list=list()
  n=1
  cell_type_list <- unique(f4.tc.combined@meta.data$cell_type4)
  #exclude mast cell, due to lack of cells, this raise error
  cell_type_list <- setdiff(cell_type_list,c("Mast cell","Thyroid cell"))
  for (t_cell_type in cell_type_list){
    print(t_cell_type)
    t_cell_type2= tolower(gsub("/","_",gsub(" ","_", t_cell_type)))
    
    kegg_res <- read_tsv(paste0(out_dir,'/gsea.kegg.',t_cell_type2,'.tsv'))
    f_res <-  kegg_res %>% filter(pval < 0.05) 
    #f_res$short_name <- make.unique(map_chr(f_res$pathway, make_short_name, 50))
    #f_res$short_name <- str_replace_all(str_replace(f_res$short_name,"KEGG_",""),"_"," ")
    x1 <- f_res %>% arrange(NES) %>% head(n=5) %>% pull(pathway)
    x2 <- f_res %>% arrange(NES) %>% tail(n=5) %>% pull(pathway)
    x_order = unique(c(x1,x2))
    g_list[[n]] <- ggplot(f_res, aes(x=pathway, y=NES, fill= log10(pval)))+
      geom_bar(stat='identity')+
      scale_x_discrete(limits = x_order)+
      scale_fill_gradient2(low="#132b42",mid="#56b1f7",high="white",midpoint=log10(0.05),limits=c(log10(0.001),log10(1)),
                           breaks=log10(c(0.001,0.01,0.05, 0.1)),
                           labels=c(0.001,0.01,0.05,0.1))+
      theme_bw() + theme(axis.title =element_blank(), legend.position = "none")+
      ggtitle(t_cell_type)+
      coord_flip()
    n=n+1
  }
  #get legend
  lgd_title <- expression(paste(italic("p"), "-value"))
  g <- ggplot(f_res, aes(x=pathway, y=NES, fill= log10(pval)))+
    geom_bar(stat='identity')+
    scale_x_discrete(limits = x_order)+
    scale_fill_gradient2(low="#132b42",mid="#56b1f7",high="white",midpoint=log10(0.05),limits=c(log10(0.001),log10(1)),
                         breaks=log10(c(0.001,0.01,0.05,1)),
                         labels=c(0.001,0.01,0.05,1))+
    labs(fill=lgd_title)+
    theme_bw() + theme(axis.title =element_blank())+
    ggtitle(t_cell_type_name)+
    coord_flip()
  lgd <- get_legend(g)

  png(paste0(out_dir,"/gsea_tme_cell_type.png"), height=6, width=16, units="in", res=300)
  plot_grid(plotlist = g_list, ncol=3) %>% print()
  dev.off()
    
  pdf(paste0(out_dir,"/gsea_tme_cell_type.pdf"), height=6, width=16)
  plot_grid(plotlist = g_list, ncol=3) %>% print()
  dev.off()

  png(paste0(out_dir,"/gsea_tme_cell_type_lgd.png"), height=4, width=4,units="in", res=300)
  ggdraw(lgd) %>%print()
  dev.off()
  
  pdf(paste0(out_dir,"/gsea_tme_cell_type_lgd.pdf"), height=4, width=4)
  ggdraw(lgd) %>% print()
  dev.off()
}

#Redraw clean plot of thyroid cell volcano and GSEA -----
#volcano
if(T){
  t_cell_type2="thyroid_cell"
  res_dt <- read_tsv(paste0(out_dir,'/deg.',t_cell_type2,'.tsv'))
  kegg_res <- read_tsv(paste0(out_dir,'/gsea.kegg.',t_cell_type2,'.tsv'))
  gobp_res <- read_tsv(paste0(out_dir,'/gsea.gobp.',t_cell_type2,'.tsv'))
  
  mito_paths <- gobp_res %>% arrange(desc(NES)) %>% pull(pathway) %>% head(n=20) %>% .[c(1,2,3,4,6,7,8,9,11,12,14)]
  
  mito_genes <- c()
  for (path in mito_paths){
    mito_genes <- c(mito_genes, unlist(gobp_list[path]))
  }
  res_dt <- res_dt%>% mutate(mito_meta_gene = ifelse(gene_name %in% mito_genes, "yes","no"))
  g <- ggplot(res_dt, aes(x=log2FoldChange,y=-log10(padj)))+
    geom_point(aes(color=mito_meta_gene, size=mito_meta_gene))+
    scale_size_manual(values=c("yes"=2, "no"=0.2))+
    scale_color_manual(values=c("yes"="red", "no"="grey30"))+
    theme_syp
  print(g)
  #not that good for publication
  
  f_res <-  kegg_res %>% filter(pval < 0.05) 
  #f_res$short_name <- make.unique(map_chr(f_res$pathway, make_short_name, 50))
  #f_res$short_name <- str_replace_all(str_replace(f_res$short_name,"KEGG_",""),"_"," ")
  x1 <- f_res %>% arrange(NES) %>% head(n=10) %>% pull(pathway)
  x2 <- f_res %>% arrange(NES) %>% tail(n=10) %>% pull(pathway)
  x_order = unique(c(x1,x2))
  g1 <- ggplot(f_res, aes(x=pathway, y=NES, fill= log10(pval)))+
    geom_bar(stat='identity')+
    scale_x_discrete(limits = x_order)+
    scale_fill_gradient2(low="#132b42",mid="#56b1f7",high="white",midpoint=log10(0.05),limits=c(log10(0.001),log10(1)),
                         breaks=log10(c(0.001,0.01,0.05, 0.1)),
                         labels=c(0.001,0.01,0.05,0.1))+
    theme_bw() + theme(axis.title =element_blank(), legend.position = "none")+
    ggtitle("KEGG")+
    coord_flip()
  
  f_res <-  gobp_res %>% filter(pval < 0.05) 
  #f_res$short_name <- make.unique(map_chr(f_res$pathway, make_short_name, 50))
  #f_res$short_name <- str_replace_all(str_replace(f_res$short_name,"KEGG_",""),"_"," ")
  x1 <- f_res %>% arrange(NES) %>% head(n=10) %>% pull(pathway)
  x2 <- f_res %>% arrange(NES) %>% tail(n=10) %>% pull(pathway)
  x_order = unique(c(x1,x2))
  g2 <- ggplot(f_res, aes(x=pathway, y=NES, fill= log10(pval)))+
    geom_bar(stat='identity')+
    scale_x_discrete(limits = x_order)+
    scale_fill_gradient2(low="#132b42",mid="#56b1f7",high="white",midpoint=log10(0.05),limits=c(log10(0.001),log10(1)),
                         breaks=log10(c(0.001,0.01,0.05, 0.1)),
                         labels=c(0.001,0.01,0.05,0.1))+
    theme_bw() + theme(axis.title =element_blank(), legend.position = "none")+
    ggtitle("Gene Ontology Biological Process")+
    coord_flip()
  
  gg <- plot_grid(g2,g1,nrow=1)
  png(paste0(out_dir,"/gsea_thyroid.png"), height=8, width=18, units="in", res=300)
  print(gg)
  dev.off()
  
  pdf(paste0(out_dir,"/gsea_thyroid.pdf"), height=8, width=18)
  print(gg)
  dev.off()
  
}
