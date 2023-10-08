
library(Seurat)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(fgsea)
theme_syp = theme_bw() + 
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=12), 
        panel.border = element_blank(), 
        axis.line = element_line(), plot.title = element_text(size=12, face="bold"))

make_short_name <- function(x, max_n){  #x= string, n=maximum character
  if (nchar(x) < max_n){
    return(x)
  } else{
    return(substr(x, 1,max_n))
  }
}

draw_volcano_gsea_plot_seurat <- function(mk, geneset_list, name1=NULL, name2=NULL,
                                          n_label=20, fc.weight=1){
  if(F){
    geneset_list=kegg_list
    name1=c('a','b','c')
    name2=c('d','e')
  }
  if(is.null(name1)==F){
    name1 = paste(name1, collapse=",")
  }
  if(is.null(name2)==F){
    name2 = paste(name2, collapse=",")
  }
  
  #draw volcano
  tmp_dt <- mk %>% as.data.frame() %>% rownames_to_column("gene_name")
  tmp_dt <- tmp_dt %>% mutate(nlpval = -log10(p_val_adj))
  tmp_dt$nlpval[is.finite(tmp_dt$nlpval) == F] <- 300
  tmp_dt <- tmp_dt %>% mutate(dist= sqrt((avg_logFC*fc.weight)^2 + nlpval ^2))
  tmp_dt <-tmp_dt %>% mutate(dist = ifelse(avg_logFC < 0, (-1)*dist, dist))
  tmp1 <- tmp_dt %>% filter(avg_logFC <0) %>%  arrange(dist) %>% head(n_label)
  tmp2 <- tmp_dt %>% filter(avg_logFC > 0) %>% arrange(dist) %>% tail(n_label)
  label_dt <- bind_rows(tmp1, tmp2)
  g1<- ggplot(tmp_dt, aes(x=avg_logFC,y=-log10(p_val_adj)))+
    geom_point()+
    geom_text_repel(data=label_dt, aes(label=gene_name), max.overlaps = n_label)+
    theme_syp
  if(is.null(name1)==F | is.null(name2)==F){
    g1 <- g1 + ggtitle(paste0(name2, "<-->", name1))
  }
  
  #GSEA
  dds_stats <- tmp_dt$dist
  names(dds_stats) <- tmp_dt$gene_name
  gsea_res <- fgsea(pathways =geneset_list, stats = dds_stats, nper=1000)
  m_res <- gsea_res %>% as_tibble()
  f_res <- m_res %>% filter(pval < 0.05) 
  f_res$rep_genes <-map_chr(f_res$leadingEdge, function(x) paste(x[1:min(length(x),5)], collapse=','))
  f_res$short_name <- make.unique(map_chr(f_res$pathway, make_short_name, 50))
  x1 <- f_res %>% arrange(NES) %>% head(n=20) %>% pull(short_name)
  x2 <- f_res %>% arrange(NES) %>% tail(n=20) %>% pull(short_name)
  x_order = unique(c(x1,x2))
  g2 <- ggplot(f_res, aes(x=short_name, y=NES, fill= pval))+
    geom_bar(stat='identity')+
    geom_text(aes(y=5, label=rep_genes), hjust=0)+
    scale_y_continuous(limits=c(-5,15))+
    scale_x_discrete(limits = x_order)+
    coord_flip()+
    theme_syp+theme(axis.text.y = element_text(size=8))
  if(is.null(name1)==F | is.null(name2)==F){
    g2 <- g2 + ggtitle(paste0(name2, "<-->", name1))
  }
  return(list(plot=plot_grid(g1,g2,nrow=1, rel_widths=c(1,1.5)), gsea=m_res))

}

draw_volcano_gsea_plot_seurat2 <- function(mk, geneset_list, 
                                           name1=NULL, name2=NULL,
                                          n_label=20, norm.quant=0.9, 
                                          plot="both"){
  #v2 apply quantile normalization to distance
  # plot =c("both","volcano","gsea")
  if(F){
    geneset_list=kegg_list
    name1=c('a','b','c')
    name2=c('d','e')
  }
  if(is.null(name1)==F){
    name1 = paste(name1, collapse=",")
  }
  if(is.null(name2)==F){
    name2 = paste(name2, collapse=",")
  }
  
  #draw volcano
  tmp_dt <- mk %>% as.data.frame() %>% rownames_to_column("gene_name")
  tmp_dt <- tmp_dt %>% mutate(nlpval = -log10(p_val_adj))
  tmp_dt$nlpval[is.finite(tmp_dt$nlpval) == F] <- 300
  pos_dt <- tmp_dt %>% filter(avg_logFC >= 0)
  pos_q_nlpval <- quantile(pos_dt$nlpval, norm.quant)
  print(pos_q_nlpval)
  pos_q_fc <- quantile(pos_dt$avg_logFC, norm.quant)
  print(pos_q_fc)
  pos_dt <- pos_dt %>% mutate(dist= sqrt((avg_logFC/pos_q_fc)^2 + (nlpval/pos_q_nlpval)^2))
  neg_dt <- tmp_dt %>% filter(avg_logFC < 0)
  neg_q_nlpval <- quantile(neg_dt$nlpval, norm.quant)
  print(neg_q_nlpval)
  neg_q_fc <- quantile(abs(neg_dt$avg_logFC), norm.quant)
  print(neg_q_fc)
  neg_dt <- neg_dt %>% mutate(dist= -1*sqrt((avg_logFC/neg_q_fc)^2 + (nlpval/neg_q_nlpval)^2))
  tmp_dt <- bind_rows(pos_dt, neg_dt)
  #tmp_dt <- tmp_dt %>% mutate(dist= sqrt((avg_logFC*fc.weight)^2 + nlpval ^2))
  #tmp_dt <-tmp_dt %>% mutate(dist = ifelse(avg_logFC < 0, (-1)*dist, dist))
  tmp1 <- tmp_dt %>% filter(avg_logFC <0) %>%  arrange(dist) %>% head(n_label)
  tmp2 <- tmp_dt %>% filter(avg_logFC > 0) %>% arrange(dist) %>% tail(n_label)
  label_dt <- bind_rows(tmp1, tmp2)
  g1<- ggplot(tmp_dt, aes(x=avg_logFC,y=-log10(p_val_adj)))+
    geom_point()+
    geom_text_repel(data=label_dt, aes(label=gene_name), max.overlaps = n_label)+
    theme_syp
  if(is.null(name1)==F | is.null(name2)==F){
    g1 <- g1 + ggtitle(paste0(name2, "<-->", name1))
  }
  
  #GSEA
  dds_stats <- tmp_dt$dist
  names(dds_stats) <- tmp_dt$gene_name
  gsea_res <- fgsea(pathways =geneset_list, stats = dds_stats, nper=1000)
  m_res <- gsea_res %>% as_tibble()
  f_res <- m_res %>% filter(pval < 0.05) 
  f_res$rep_genes <-map_chr(f_res$leadingEdge, function(x) paste(x[1:min(length(x),5)], collapse=','))
  f_res$short_name <- make.unique(map_chr(f_res$pathway, make_short_name, 50))
  x1 <- f_res %>% arrange(NES) %>% head(n=20) %>% pull(short_name)
  x2 <- f_res %>% arrange(NES) %>% tail(n=20) %>% pull(short_name)
  x_order = unique(c(x1,x2))
  min_nes <- min(f_res$NES)
  max_nes <- max(f_res$NES)
  g2 <- ggplot(f_res, aes(x=short_name, y=NES, fill= pval))+
    geom_bar(stat='identity')+
    geom_text(aes(y=max_nes+0.5, label=rep_genes),size=3, hjust=0)+
    scale_y_continuous(limits=c(min_nes,max_nes + 5))+
    scale_x_discrete(limits = x_order)+
    coord_flip()+
    theme_syp+theme(axis.text.y = element_text(size=8))
  if(is.null(name1)==F | is.null(name2)==F){
    g2 <- g2 + ggtitle(paste0(name2, "<-->", name1))
  }
  if(plot == "both"){
    return(list(plot=plot_grid(g1,g2,nrow=1, rel_widths=c(1,1.5)), gsea=m_res))
  }else if (plot == "volcano"){
    return(list(plot=g1, gsea=m_res))
  }else if (plot == "gsea"){
    return(list(plot=g2, gsea=m_res))
  }
}

#functions
draw_umap_plot <- function(tc, ident,set_alpha=0.5, remove_na=F,
                           txt_annot=TRUE, title=NULL){
  if(F){
    tc = hcc
    ident="CTNNB1"
    set_alpha=0.5
  }
  library(ggrepel)
  Idents(tc) <- ident
  umap_dt <- tc@reductions$umap@cell.embeddings
  umap_dt <- umap_dt %>% as.data.frame %>% rownames_to_column('cell_id')
  dim(umap_dt)
  ident_dt <- Idents(tc) %>% as.data.frame %>% rownames_to_column('cell_id') 
  colnames(ident_dt) <- c('cell_id','cluster')
  umap_dt <- left_join(umap_dt, ident_dt)
  #max_clust <- umap_dt$cluster %>% as.character() %>% as.numeric() %>% max()
  umap_dt$cluster <- factor(umap_dt$cluster, levels = sort(as.character(unique(umap_dt$cluster))))
  umap_dt <- umap_dt %>% separate(cell_id, c('barcode','orig.ident'),sep='-', remove=F) %>% as_tibble()
  umap_pos <- umap_dt %>% group_by(cluster) %>% summarise(med1=median(UMAP_1), med2=median(UMAP_2))
  if(remove_na == T){
    umap_dt <- umap_dt %>% filter(is.na(cluster)==F)
  }
  g <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(aes(color=cluster), alpha=set_alpha, size=0.5)+
    guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
    labs(color="")+
    theme_bw()+theme(panel.grid = element_blank())
  if(txt_annot == T){
    g <- g +geom_text_repel(data=umap_pos, aes(x=med1, y=med2, label=cluster), size=5)
  }
  if(is.null(title)==F){
    g <- g+ ggtitle(title)
  }
  print(g)
}

draw_orig_plot <- function(tc, orig, set_alpha=0.5){
  library(ggrepel)
  orig_dt <- tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>%
    select(cell_id, orig.ident)
  umap_dt <- tc@reductions$umap@cell.embeddings
  umap_dt <- umap_dt %>% as.data.frame %>% rownames_to_column('cell_id')
  dim(umap_dt)
  umap_dt <- left_join(umap_dt, orig_dt)
  umap_dt <- umap_dt %>% mutate(group = ifelse(orig.ident == orig, 'target', 'others'))
  g <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(aes(color=group), alpha=set_alpha, size=0.5)+
    guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
    scale_color_manual(values = c('target'="red", "others"="grey90"))+
    theme_bw()+theme(panel.grid = element_blank())+
    ggtitle(orig)
  print(g)
}

draw_umap_mark_plot <- function(tc.obj, ident, t_cluster, set_alpha=0.5, remove_other=F){
  if(F){
    tc.obj <- ff.tc.combined
    ident = "integrated_snn_res.1.5"
    t_cluster = c(7,8,11,12,15,18,22)
  }
  library(ggrepel)
  Idents(tc.obj) <- ident
  umap_dt <- tc.obj@reductions$umap@cell.embeddings
  umap_dt <- umap_dt %>% as.data.frame %>% rownames_to_column('cell_id')
  dim(umap_dt)
  ident_dt <- Idents(tc.obj) %>% as.data.frame %>% rownames_to_column('cell_id') 
  colnames(ident_dt) <- c('cell_id','cluster')
  umap_dt <- left_join(umap_dt, ident_dt)
  umap_dt <- umap_dt %>% 
    mutate(group = ifelse(cluster %in% t_cluster, 'target', 'others'))
  if(remove_other==T){
    umap_dt <- umap_dt %>% filter(group == 'target')
  }
  g <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(aes(color=group), alpha=set_alpha, size=0.5)+
    guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
    scale_color_manual(values = c('target'="red", "others"="grey90"))+
    theme_bw()+theme(panel.grid = element_blank())+
    ggtitle(paste0("marking = ",paste(t_cluster, collapse=",")))
  print(g)
}


compare_meta_two_groups<- function(t_obj, t_ident, t_item, g1, g2){
  tmp_meta <- t_obj@meta.data %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
  tmp <- tmp_meta %>% select(t_ident, t_item)
  colnames(tmp)<- c("c_ident","c_item")
  tmp <- tmp %>% mutate(group=ifelse(c_ident %in% g1,"g1",ifelse(c_ident %in% g2, "g2","other")))
  g <- ggplot(tmp, aes(x=group, y=c_item))+
    geom_boxplot()+
    labs(y=t_item)
  print(g)
}

