#compare SHMT2 MTHFD2 TDS after batch correction

if(T){
  library(tidyverse)
  library(Seurat)
  library(ggsci)
  library(scales)
  library(ggsignif)
  library(corrplot)
  library(cowplot)
  library(ggrepel)
}


#output path
out_dir="/path/to/output_dir"
system(paste0('mkdir -p ',out_dir))

#design
theme_syp = theme_bw() + theme(axis.text = element_text(size=12), axis.title = element_text(size=12), panel.border = element_blank(), axis.line = element_line(), plot.title = element_text(size=12, face="bold"))
cancer_type_pal<- pal_d3("category10")(10)[c(1,3,2)]
names(cancer_type_pal) <- c('NT','PTC','ATC')

#gene list
oc_gene_list <- c('PHGDH','PSAT1','PSPH','SHMT1','SHMT2','MTHFD1','MTHFD2','MTHFD1L','MTHFD2L','TYMS',
                  'MTHFR','ALDH1L1','ALDH1L2','SFXN1','SFXN2','SFXN3',"MTFMT","SLC1A4","SLC1A5","FOLH1",
                  "FOLR1","FOLR2","DHFR","SLC2A1","SLC2A3","HK3","GAPDH","PKM")

tds_genes = c('DIO1','DIO2','DUOX1','DUOX2','FOXE1','GLIS3','NKX2-1','PAX8','SLC26A4','SLC5A5','SLC5A8','TG','THRA','THRB','TPO','TSHR')
upr_mt_genes = c("CLPP", "LONP1", "TIMM23", "TOMM40", "ATF5", "HSPD1", "HSPE1", "NRF1", "GDF15", "ATF6", "ATF4", "EIF2A", "IRE1A","XBP1", "SHMT2", "MTHFD2")
fibro_endo_thyro_genes = c('DCN','COL1A1','ACTA2','COL3A1','CDH5','PECAM1','TG','PAX8','TPO')
immune_genes=c('PTPRC','CD68','CD14','CD4','CD40LG','FOXP3','CD8A','NKG7','CD79A','CD79B','IGHG1','TPSAB1')  

h_list <- GSEABase::getGmt('/path/to/h.all.v7.5.1.symbols.gmt') %>% GSEABase::geneIds()
c5_list <- GSEABase::getGmt('/path/to/c5.all.v7.4.symbols.gmt') %>% GSEABase::geneIds()
c2_list <- GSEABase::getGmt('/path/to/c2.all.v7.4.symbols.gmt') %>% GSEABase::geneIds()

gobp_names <- names(c5_list)[grepl('GOBP',names(c5_list))]
gobp_list <- c5_list[gobp_names]
gocc_names <- names(c5_list)[grepl('GOCC',names(c5_list))]
gocc_list <- c5_list[gocc_names]
kegg_names <- names(c2_list)[grepl('KEGG',names(c2_list))]
kegg_list <- c2_list[kegg_names]

#Script specific functions-----
make_signif_mark <- function(pval){
  if (pval < 0.001){
    mk='***'
  }else if (pval < 0.01){
    mk='**'
  }else if (pval < 0.05){
    mk='*'
  }else {
    mk='NS'
  }
  return(mk)
}


draw_violin_plot <- function(tmp_dt, target_gene){
  print(target_gene)
  nt <- tmp_dt %>% filter(cancer_type == 'NT') %>% pull(target_gene)
  pt <- tmp_dt %>% filter(cancer_type == 'PTC') %>% pull(target_gene)
  at <- tmp_dt %>% filter(cancer_type == 'ATC') %>% pull(target_gene)
  ca <- tmp_dt %>% filter(cancer_type %in% c('PTC','ATC')) %>% pull(target_gene)
  nt_ca <- t.test(nt, ca)
  print('normal vs tumor')
  print(paste0('p value=',nt_ca$p.value))
  nt_ca_mk <- make_signif_mark(nt_ca$p.value)
  if(nt_ca_mk == "NS"){
    nt_ca_size=5
  }else{
    nt_ca_size=10
  }
  pt_at <- t.test(pt, at)
  print('PTC vs ATC')
  print(paste0('p value=',pt_at$p.value))
  pt_at_mk <- make_signif_mark(pt_at$p.value)
  if(pt_at_mk == "NS"){
    pt_at_size=5
  }else{
    pt_at_size=10
  }
  max_v=max(c(nt, pt, at))
  
  g <-  ggplot(tmp_dt, aes(x=cancer_type, y=get(target_gene), fill=cancer_type))+
    geom_violin()+
    geom_segment(x=0.8, xend=1.2, y=mean(nt), yend=mean(nt), color="red")+
    geom_segment(x=1.8, xend=2.2, y=mean(pt), yend=mean(pt), color="red")+
    geom_segment(x=2.8, xend=3.2, y=mean(ca), yend=mean(ca), color="red")+
    geom_text(aes(x=1,y=max_v*1.1), label=round(mean(nt),3))+
    geom_text(aes(x=2,y=max_v*1.1), label=round(mean(pt),3))+
    geom_text(aes(x=3,y=max_v*1.1), label=round(mean(ca),3))+
    geom_signif(annotation=nt_ca_mk, xmin=1, xmax=2.5, y_position=max_v*1.4, tip_length=c(0.01,0.01), textsize=nt_ca_size)+
    geom_signif(annotation=pt_at_mk, xmin=2, xmax=3, y_position=max_v*1.2, tip_length=c(0.01,0.01), textsize=pt_at_size)+
    scale_y_continuous(limits=c(0,max_v*1.5))+
    scale_fill_manual(values = cancer_type_pal)+
    ylab("Expresion level (log CPM)")+
    theme_syp+theme(axis.title.x=element_blank(), legend.position='none', panel.grid=element_blank())+
    ggtitle(target_gene)
  return(g)
}

draw_tds_correlation <- function(tmp_dt, target_gene){
  print(target_gene)
  
  g <-  ggplot(data=tmp_dt, aes(x=get(target_gene), y=tds_score))+
    geom_point(aes(color=cancer_type),alpha=0.3)+
    geom_smooth(method="lm", color="red")+
    scale_color_manual(values=cancer_type_pal)+
    xlab("Expression level (log CPM)")+ylab("TDS score")+
    theme_syp+theme(panel.grid=element_blank(), legend.position = 'none')+
    #coord_cartesian(ylim=c(0,3), xlim=c(0,700))+
    #guides(colour = guide_legend(override.aes = list(alpha = 1, size=3), title=element_blank()))+
    ggtitle(target_gene)
  return(g)
}


draw_umap_plot <- function(tc, ident){
  library(ggrepel)
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
    theme_bw()+theme(panel.grid = element_blank())
  print(g)
}

draw_umap_color_plot <- function(umap_dt, gene_name){
  umap_dt <- umap_dt %>% dplyr::rename(target=gene_name)
  umap_dt$target[umap_dt$target > 3] <- 3
  g <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(aes(color=target), alpha=0.5, size=1)+
    scale_color_gradientn(colors = c("grey","orangered", "red1", "red3", "red4", "darkred"), 
                          limits=c(0,3), labels = c(0,1,2,">3"))+
    theme_bw()+theme(panel.grid = element_blank())+
    labs(x="UMAP1",y="UMAP2",color="Expression\nlevel")+
    ggtitle(gene_name)
  print(g)
}

draw_orig_plot <- function(tc, orig){
  library(ggrepel)
  orig_dt <- tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>%
    select(cell_id, orig.ident)
  umap_dt <- tc@reductions$umap@cell.embeddings
  umap_dt <- umap_dt %>% as.data.frame %>% rownames_to_column('cell_id')
  dim(umap_dt)
  umap_dt <- left_join(umap_dt, orig_dt)
  umap_dt <- umap_dt %>% mutate(group = ifelse(orig.ident == orig, 'target', 'others'))
  g <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(aes(color=group), alpha=0.5, size=0.5)+
    guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
    scale_color_manual(values = c('target'="red", "others"="grey90"))+
    theme_bw()+theme(panel.grid = element_blank())+
    ggtitle(orig)
  print(g)
}


make_short_name <- function(x, max_n){  #x= string, n=maximum character
  if (nchar(x) < max_n){
    return(x)
  } else{
    return(substr(x, 1,max_n))
  }
}


#load Seurat object-----
if(T){
  f4.tc.combined <- readRDS("/path/to/Seurat_object.rds")
}

#Include only thyroid cells
if(T){
  f4.tc.combined@meta.data$cell_type2 %>% unique()
  Idents(f4.tc.combined) <- "cell_type2"
  thy.combined <- subset(f4.tc.combined, idents = c('thyroid_cell','TPO_high_thyroid_cell'))
  
}

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

fr_NT3_dt <- subset(thy.combined, idents="NT3")
fr_NT5_dt <- subset(thy.combined, idents="NT5")
fr_NT10_dt <- subset(thy.combined, idents="NT10")


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

nt3_mx  <- fr_NT3_dt@assays$RNA@counts %>% as.matrix()
nt5_mx  <- fr_NT5_dt@assays$RNA@counts %>% as.matrix()
nt10_mx  <- fr_NT10_dt@assays$RNA@counts %>% as.matrix()
  
  
#generate merged matrix  
  if(T){
    common_genes <- Reduce(intersect, list(rownames(pt3_mx), rownames(pt5_mx), 
                                           rownames(pt7_mx), 
                                           rownames(pt8_mx), rownames(pt9_mx), 
                                           rownames(pt10_mx), rownames(pt12_mx),
                                           rownames(at9_mx),
                                           rownames(at13_mx), rownames(at16_mx),
                                           rownames(at17_mx), rownames(at20_mx),
                                           rownames(nt3_mx),
                                           rownames(nt5_mx),
                                           rownames(nt10_mx)))
    merged_mx <- cbind(pt3_mx[common_genes,],pt5_mx[common_genes,], 
                       pt7_mx[common_genes,], pt8_mx[common_genes,],pt9_mx[common_genes,], 
                       pt10_mx[common_genes,], pt12_mx[common_genes,],
                       at9_mx[common_genes,], at13_mx[common_genes,],at16_mx[common_genes,],
                       at17_mx[common_genes,],at20_mx[common_genes,],
                       nt3_mx[common_genes,],
                       nt5_mx[common_genes,],
                       nt10_mx[common_genes,])
    dim(merged_mx) #28745 17412
  }
  

#batch correction using ComBat_seqs-----
if(file.exists(paste0(out_dir,'/ptc7_atc5_nt3_bc_mx.rds'))==F){
  batch <- map_chr(colnames(merged_mx), function(x) unlist(strsplit(x,'_'))[1])
  unique(batch)
  bc_mx <- sva::ComBat_seq(merged_mx, batch=batch, group=NULL)
  bc_mx %>% saveRDS(paste0(out_dir,'/ptc7_atc5_nt3_bc_mx.rds'))
} else {
  bc_mx <- readRDS(paste0(out_dir,'/ptc7_atc5_nt3_bc_mx.rds'))
}

#intersect genelist
oc_gene_list <- intersect(oc_gene_list, rownames(bc_mx))
tds_genes <- intersect(tds_genes, rownames(bc_mx))
upr_mt_genes <- intersect(upr_mt_genes, rownames(bc_mx))


### make data table for expression comparison
# convert into CPM matrix-----
cpm_mx <- apply(bc_mx, 2, function(x) (x/sum(x))*1000000)
# make log CPM matrix
l_cpm_mx <- log10(cpm_mx+1)
# edit logo CPM matrix 
if(T){
  fl_cpm_mx <- l_cpm_mx[c(oc_gene_list, tds_genes),] %>% as.matrix()
  fl_cpm_dt <- fl_cpm_mx %>% t() %>% as.data.frame() %>% rownames_to_column('cell_id') %>% as_tibble()
  #TDS score calculation by scaled mean
  scaled_tds <- t(apply(fl_cpm_mx[tds_genes,], 1, scale))
  colnames(scaled_tds) <- colnames(fl_cpm_mx)
  tds_score <- colMeans(scaled_tds)
  #make Mean value of SHMT2 and MTHFD2
  SM_mean <- colMeans(fl_cpm_mx[c("SHMT2","MTHFD2"),])
  all(names(tds_score) == names(SM_mean))
  fl_cpm_dt$tds_score <- tds_score
  fl_cpm_dt$SM_mean <- SM_mean
  #make orig.ident
  fl_cpm_dt <- fl_cpm_dt %>% rowwise() %>% mutate(orig.ident=unlist(strsplit(cell_id,'_'))[1])
  # make cancer type
  fl_cpm_dt <- fl_cpm_dt %>% mutate(cancer_type=substr(cell_id,1,2))
  fl_cpm_dt$cancer_type[fl_cpm_dt$cancer_type == 'PT'] <- 'PTC'
  fl_cpm_dt$cancer_type[fl_cpm_dt$cancer_type == 'AT'] <- 'ATC'
  fl_cpm_dt$cancer_type <- factor(fl_cpm_dt$cancer_type, levels=c('NT','PTC','ATC'))  
}

#draw pca plot after batch correction-----
#filter protein coding genes
vars <- apply(l_cpm_mx, 1, function(x) var(x))
means <- rowMeans(l_cpm_mx)
tmp <- tibble(gene_name = rownames(l_cpm_mx), var=vars, mean=means)
tmp <- tmp %>% mutate(CV = var/mean)
ggplot(tmp, aes(x=mean, y=log10(CV)))+
  geom_point()+
  geom_smooth(method="loess")
res <- loess(log10(CV)~mean, data=tmp)
tmp2 <- tibble(loess_residual=res$residuals, row_number=as.integer(names(res$residuals)))
tmp$row_number <- 1:nrow(tmp)
tmp <- left_join(tmp, tmp2)
top2000 <- tmp %>% filter(mean > 0.1) %>% arrange(desc(loess_residual)) %>%
  head(n=2000) %>% pull(gene_name)
tmp <- tmp %>% mutate(selected=case_when(gene_name %in% top2000 ~ "Y",
                                  TRUE ~ "N"))
setdiff(tds_genes, top2000)
g <- ggplot(tmp, aes(x=mean, y=log10(CV)))+
  geom_point(aes(color=selected))+
  geom_smooth(method="loess")+
  theme_syp
png(paste0(out_dir,'/selection_top2000_genes.png'), width=6, height=6,
    units="in", res=300)
print(g)
dev.off()

#random select of cells
ct_tbl <- tmp_meta %>% filter(cell_type3 == "Thyroid cell") %>% group_by(orig.ident) %>% count()
n_random=100
f_tmp1 <- tmp_meta %>% filter(cell_type3 == "Thyroid cell" & orig.ident %in% ct_tbl$orig.ident[ct_tbl$n >n_random])
f_tmp2 <- tmp_meta %>% filter(cell_type3 == "Thyroid cell" & orig.ident %in% ct_tbl$orig.ident[ct_tbl$n <=n_random])
set.seed(1); rf_tmp1 <- f_tmp1 %>% 
  group_by(orig.ident) %>% 
  sample_n(n_random) %>% ungroup()
rf_tmp <- bind_rows(rf_tmp1, f_tmp2)
rf_tmp %>% group_by(orig.ident) %>% count()

pc <- prcomp(t(l_cpm_mx[top2000,rf_tmp$cell_id]), center=TRUE)
pc %>% saveRDS(paste0(out_dir,'/thyroid_combat_pc.rds'))
pc <- readRDS(paste0(out_dir,'/thyroid_combat_pc.rds'))

pc_dt <- pc$x %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
tmp_meta <- f4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column("cell_id") %>%
  as_tibble()
pc_dt <- left_join(pc_dt, tmp_meta %>% select(cell_id, cancer_type, orig.ident))
pc_dt <- left_join(pc_dt, fl_cpm_dt %>% select(cell_id, tds_score, SM_mean, 
                                               SHMT2, MTHFD2, TG, TPO, TSHR,
                                               PAX8, SHMT1, PHGDH))
pc_dt <- pc_dt %>% mutate(cancer_type2 = case_when(cancer_type == 'normal' ~'NT',
                                          cancer_type == "ptc" ~ 'PTC',
                                          cancer_type == 'atc' ~ 'ATC'))

point_size=3
g1 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=cancer_type2), alpha=0.5, size=point_size)+
  scale_color_manual(values=cancer_type_pal)+
  theme_syp
g2 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=orig.ident), alpha=0.5, size=point_size)
g3 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=tds_score), alpha=0.5, size=point_size)
g4 <-ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=SM_mean), alpha=0.5, size=point_size)
g5 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=SHMT2), alpha=0.5, size=point_size)
g6 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=MTHFD2), alpha=0.5, size=point_size)
gg <- plot_grid(g1,g2,g3,g4,g5,g6,nrow=2)
png(paste0(out_dir,'/pc1_pc2_r100.png'), width=16, height=8,
    units="in", res=300)
print(gg)
dev.off()

pc_dt$cancer_type2 <- factor(pc_dt$cancer_type2, levels = c("NT","PTC","ATC"))
g <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=cancer_type2), alpha=0.5, size=point_size)+
  scale_color_manual(values=cancer_type_pal)+
  theme_syp+theme(panel.grid = element_blank())+
  guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
  labs(color="")

png(paste0(out_dir,'/pca_caner_type.png'), width=8, height=8, units="in",res=300)
print(g)
dev.off()

pdf(paste0(out_dir,'/pca_caner_type.pdf'), width=8, height=8)
print(g)
dev.off()

g1 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=tds_score), alpha=0.5, size=3)+
  scale_color_gradientn(colors = c("light yellow","yellow","orangered", "red", "dark red"))+
  theme(panel.grid = element_blank())+
  labs(x="PC1",y="PC2",color="Expression\nlevel")+
  ggtitle("TDS score")
g1

g2 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=SM_mean), alpha=0.5, size=3)+
  scale_color_gradientn(colors =c("light yellow","yellow","yellow","orangered","red", "red", "red","dark red"))+
  theme(panel.grid = element_blank())+
  labs(x="PC1",y="PC2",color="Expression\nlevel")+
  ggtitle("Mean of SHMT2 and MTHFD2")
g2

gg <- plot_grid(g1,g2,ncol=1)
png(paste0(out_dir,'/pca_tds_sm2.png'), width=6, height=8, units="in",res=300)
print(gg)
dev.off()

pdf(paste0(out_dir,'/pca_tds_sm2.pdf'), width=6, height=8)
print(gg)
dev.off()


g1 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=TG), alpha=0.5, size=point_size)+
  ggtitle("TG")+
  theme_syp+theme(panel.grid = element_blank())
g2 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=TPO), alpha=0.5, size=point_size)+
  ggtitle("TPO")+
  theme_syp+theme(panel.grid = element_blank())
g3 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=TSHR), alpha=0.5, size=point_size)+
  ggtitle("TSHR")+
  theme_syp+theme(panel.grid = element_blank())
g4 <-ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=PAX8), alpha=0.5, size=point_size)+
  ggtitle("PAX8")+
  theme_syp+theme(panel.grid = element_blank())
g5 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=SHMT2), alpha=0.5, size=point_size)+
  ggtitle("SHMT2")+
  theme_syp+theme(panel.grid = element_blank())
g6 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=MTHFD2), alpha=0.5, size=point_size)+
  ggtitle("MTHFD2")+
  theme_syp+theme(panel.grid = element_blank())
g7 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=SHMT1), alpha=0.5, size=point_size)+
  ggtitle("SHMT1")+
  theme_syp+theme(panel.grid = element_blank())
g8 <- ggplot(pc_dt, aes(x=PC1, y=PC2))+
  geom_point(aes(color=PHGDH), alpha=0.5, size=point_size)+
  ggtitle("PHGDH")+
  theme_syp+theme(panel.grid = element_blank())
gg <- plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,nrow=2)
print(gg)
png(paste0(out_dir,'/pc1_pc2_r100_gene.png'), width=12, height=6,
    units="in", res=300)
print(gg)
dev.off()

pdf(paste0(out_dir,'/pc1_pc2_r100_gene.pdf'), width=12, height=6)
print(gg)
dev.off()

g1 <- ggplot(pc_dt, aes(x=PC2, y=PC3))+
  geom_point(aes(color=cancer_type), alpha=0.5)
g2 <- ggplot(pc_dt, aes(x=PC2, y=PC3))+
  geom_point(aes(color=orig.ident), alpha=0.5)
g3 <- ggplot(pc_dt, aes(x=PC2, y=PC3))+
  geom_point(aes(color=tds_score), alpha=0.5)
g4 <- ggplot(pc_dt, aes(x=PC2, y=PC3))+
  geom_point(aes(color=SM_mean), alpha=0.5)
g5 <- ggplot(pc_dt, aes(x=PC2, y=PC3))+
  geom_point(aes(color=SHMT2), alpha=0.5)
g6 <- ggplot(pc_dt, aes(x=PC2, y=PC3))+
  geom_point(aes(color=MTHFD2), alpha=0.5)
gg <- plot_grid(g1,g2,g3,g4,g5,g6,nrow=2)
print(gg)


#write tds score and SM mean to meta data-----
tmp_meta <- thy.combined@meta.data %>% as.data.frame() %>% 
  rownames_to_column("cell_id") %>% as_tibble()
if ("tds_score" %in% colnames(tmp_meta)){
  tmp_meta <- tmp_meta %>% select(-tds_score)
}
if ("SM_mean" %in% colnames(tmp_meta)){
  tmp_meta <- tmp_meta %>% select(-SM_mean)
} 
tmp_meta <- left_join(tmp_meta,fl_cpm_dt %>% select(cell_id, tds_score, SM_mean))
thy.combined@meta.data <-  tmp_meta %>% as.data.frame() %>% column_to_rownames("cell_id")


#single cell expression comparison among cancer type -----
g1 <- draw_violin_plot(fl_cpm_dt, "SHMT2")
g2 <- draw_violin_plot(fl_cpm_dt, "MTHFD2")
g3 <- draw_violin_plot(fl_cpm_dt, "PHGDH")
g4 <- draw_violin_plot(fl_cpm_dt, "SHMT1")
gg <- plot_grid(g1,g2,g3,g4, nrow=2)

png(paste0(out_dir,'/exp_violin.png'), width=8, height=8, units="in",res=300)
print(gg)
dev.off()

pdf(paste0(out_dir,'/exp_violin.pdf'), width=8, height=8)
print(gg)
dev.off()
#Other genes
g1 <- draw_violin_plot(fl_cpm_dt, "MTFMT")
g2 <- draw_violin_plot(fl_cpm_dt, "MTHFD1L")
g3 <- draw_violin_plot(fl_cpm_dt, "MTHFD2L")
g4 <- draw_violin_plot(fl_cpm_dt, "MTHFD1")
g5 <- draw_violin_plot(fl_cpm_dt, "MTHFR")
g6 <- draw_violin_plot(fl_cpm_dt, "PSPH")
g7 <- draw_violin_plot(fl_cpm_dt, "SFXN1")
g8 <- draw_violin_plot(fl_cpm_dt, "SFXN2")
g9 <- draw_violin_plot(fl_cpm_dt, "SFXN3")
g10 <- draw_violin_plot(fl_cpm_dt, "SLC1A4")
g11 <- draw_violin_plot(fl_cpm_dt, "SLC1A5")
g12 <- draw_violin_plot(fl_cpm_dt, "FOLH1")
gg <- plot_grid(g1, g2, g3,
                g4, g5, g6,
                g7, g8, g9,
                g10,g11,g12,
                nrow=2)
png(paste0(out_dir,'/exp_violin_other.png'), width=20, height=8, units="in",res=300)
print(gg)
dev.off()

pdf(paste0(out_dir,'/exp_violin_other.pdf'), width=16, height=12)
print(gg)
dev.off()

#other genes2
#"FOLR1","FOLR2","DHFR","SLC2A1","SLC2A3","HK3","GAPDH","PKM"
g1 <- draw_violin_plot(fl_cpm_dt, "FOLR1")
g2 <- draw_violin_plot(fl_cpm_dt, "FOLR2")
g3 <- draw_violin_plot(fl_cpm_dt, "DHFR")
g4 <- draw_violin_plot(fl_cpm_dt, "SLC2A1")
g5 <- draw_violin_plot(fl_cpm_dt, "SLC2A3")
g6 <- draw_violin_plot(fl_cpm_dt, "HK3")
g7 <- draw_violin_plot(fl_cpm_dt, "GAPDH")
g8 <- draw_violin_plot(fl_cpm_dt, "PKM")
g9 <- draw_violin_plot(fl_cpm_dt, "PSAT1")
g10 <- draw_violin_plot(fl_cpm_dt, "TYMS")
g11 <- draw_violin_plot(fl_cpm_dt, "ALDH1L2")
gg <- plot_grid(g1, g2, g3,
                g4, g5, g6,
                g7, g8, g9,
                g10,g11,
                nrow=3)
png(paste0(out_dir,'/exp_violin_other2.png'), width=16, height=12, units="in",res=300)
print(gg)
dev.off()

pdf(paste0(out_dir,'/exp_violin_other2.pdf'), width=16, height=12)
print(gg)
dev.off()


# TDS gene expression correlation in single cell level ------
g1 <- draw_tds_correlation(fl_cpm_dt, "SHMT2")
g2 <- draw_tds_correlation(fl_cpm_dt, "MTHFD2")
g3 <- draw_tds_correlation(fl_cpm_dt, "PHGDH")
g4 <- draw_tds_correlation(fl_cpm_dt, "SHMT1")
gg <- plot_grid(g1,g2,g3,g4, nrow=1)

png(paste0(out_dir,'/tds_gene_cor_cell.png'), width=20, height=5, units="in",res=300)
print(gg)
dev.off()

pdf(paste0(out_dir,'/tds_gene_cor_cell.pdf'), width=20, height=5)
print(gg)
dev.off()

#  TDS gene expression correlation in tissue level-----
mean_dt <- tmp_dt  %>% group_by(orig.ident) %>% summarise(m_tds_score = mean(tds_score),
                                                      m_SHMT2 = mean(SHMT2), 
                                                      m_MTHFD2 = mean(MTHFD2), 
                                                      m_PHGDH = mean(PHGDH), 
                                                      m_SHMT1 = mean(SHMT1)
)
mean_dt <- mean_dt %>% mutate(cancer_type = substr(orig.ident,1,2))
mean_dt$cancer_type[mean_dt$cancer_type == "PT"]<- "PTC"
mean_dt$cancer_type[mean_dt$cancer_type == "AT"]<- "ATC"

print("SHMT2")
res = cor.test(mean_dt$m_SHMT2, mean_dt$m_tds_score)
cor_r=res$estimate
cor_p=res$p.value
print(cor_r)
print(cor_p)

g1 <- ggplot(mean_dt, aes(x=m_SHMT2, y=m_tds_score))+
  geom_point(aes(color=cancer_type), size=5, alpha=0.7)+
  geom_smooth(method="lm", color="red", se=F)+
  labs(x="Expression level (log CPM)", y="TDS score", color="")+
  scale_color_manual(values = cancer_type_pal)+
  theme_syp+theme(legend.position='none')+
  ggtitle("SHMT2")

print("MTHFD2")
res = cor.test(mean_dt$m_MTHFD2, mean_dt$m_tds_score)
cor_r=res$estimate
cor_p=res$p.value
print(cor_r)
print(cor_p)

g2 <- ggplot(mean_dt, aes(x=m_MTHFD2, y=m_tds_score))+
  geom_point(aes(color=cancer_type), size=5, alpha=0.7)+
  geom_smooth(method="lm", color="red", se=F)+
  labs(x="Expression level (log CPM)", y="TDS score", color="")+
  scale_color_manual(values = cancer_type_pal)+
  theme_syp+theme(legend.position='none')+
  ggtitle("MTHFD2")

print("PHGDH")
res = cor.test(mean_dt$m_PHGDH, mean_dt$m_tds_score)
cor_r=res$estimate
cor_p=res$p.value
print(cor_r)
print(cor_p)

g3 <- ggplot(mean_dt, aes(x=m_PHGDH, y=m_tds_score))+
  geom_point(aes(color=cancer_type), size=5, alpha=0.7)+
  geom_smooth(method="lm", color="red", se=F)+
  labs(x="Expression level (log CPM)", y="TDS score", color="")+
  scale_color_manual(values = cancer_type_pal)+
  theme_syp+theme(legend.position='none')+
  ggtitle("PHGDH")

print("SHMT1")
res = cor.test(mean_dt$m_SHMT1, mean_dt$m_tds_score)
cor_r=res$estimate
cor_p=res$p.value
print(cor_r)
print(cor_p)

g4 <- ggplot(mean_dt, aes(x=m_SHMT1, y=m_tds_score))+
  geom_point(aes(color=cancer_type), size=5, alpha=0.7)+
  geom_smooth(method="lm", color="red", se=F)+
  labs(x="Expression level (log CPM)", y="TDS score", color="")+
  scale_color_manual(values = cancer_type_pal)+
  theme_syp+theme(legend.position='none')+
  ggtitle("SHMT1")
gg <- plot_grid(g1, g2, g3, g4, nrow=2)       

png(paste0(out_dir,'/tds_gene_cor_tissue.png'), width=8, height=8, units="in",res=300)
print(gg)
dev.off()

pdf(paste0(out_dir,'/tds_gene_cor_tissue.pdf'), width=8, height=8)
print(gg)
dev.off()

#########correlation between genes
target_mx1 <- l_cpm_mx[c(oc_gene_list, tds_genes),]
res <- cor(t(target_mx1))
png(paste0(out_dir,'/oc_tds_gene_cor_plot.png'), width=8, height=8, units="in",res=300)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, hclust.method="ward.D2")
dev.off()


######single cell plotting
#umap cancer type -----
tmp_meta <- thy.combined@meta.data %>% as.data.frame() %>% 
  rownames_to_column("cell_id") %>% as_tibble()

tmp_mx <- l_cpm_mx[c(oc_gene_list, tds_genes,"PTPRC","DCN","EPCAM","KRT18","KRT19"),] %>% as.matrix()
tmp_dt <- tmp_mx %>% t() %>% as.data.frame() %>% rownames_to_column('cell_id') %>% as_tibble()

umap_dt <- thy.combined@reductions$umap@cell.embeddings
umap_dt <- umap_dt %>% as.data.frame %>% rownames_to_column('cell_id')
umap_dt <- left_join(umap_dt, tmp_meta %>% select(cell_id, tds_score, SM_mean))
umap_dt <- left_join(umap_dt, tmp_dt %>% select(cell_id, SHMT2, MTHFD2,PHGDH, SHMT1))
umap_dt <- umap_dt %>% mutate(cancer_type = substr(cell_id, 1,2))
umap_dt$cancer_type[umap_dt$cancer_type == 'PT'] <- 'PTC'
umap_dt$cancer_type[umap_dt$cancer_type == 'AT'] <- 'ATC'
umap_dt$cancer_type <- factor(umap_dt$cancer_type, levels=c("NT","PTC","ATC"))

g <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
  geom_point(aes(color=cancer_type), alpha=0.5, size=1)+
  scale_color_manual(values = cancer_type_pal)+
  theme_bw()+theme(panel.grid = element_blank())+
  guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
  labs(x="UMAP1",y="UMAP2", color="")
png(paste0(out_dir,'/umap_caner_type.png'), width=8, height=8, units="in",res=300)
print(g)
dev.off()

pdf(paste0(out_dir,'/umap_caner_type.pdf'), width=8, height=8)
print(g)
dev.off()


#SHmean, TDS score on UMAP plot-----
umap_dt <- umap_dt %>% filter(!grepl('NT',cell_id))
umap_dt <- umap_dt %>% mutate(SM=ifelse(SHMT2 > 0 | MTHFD2 > 0, 'y','n'))

g1 <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
  geom_point(aes(color=tds_score), alpha=0.5, size=2)+
  scale_color_gradientn(colors = c("light yellow","yellow","orangered", "red", "dark red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x="UMAP1",y="UMAP2",color="Expression\nlevel")+
  ggtitle("TDS score")
g1

g2 <- ggplot(umap_dt, aes(x=UMAP_1, y=UMAP_2))+
  geom_point(aes(color=SM_mean), alpha=0.5, size=2)+
  scale_color_gradientn(colors =c("light yellow","yellow","yellow","orangered","red", "red", "red","dark red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x="UMAP1",y="UMAP2",color="Expression\nlevel")+
  ggtitle("Mean of SHMT2 and MTHFD2")
g2

gg <- plot_grid(g1,g2,ncol=1)
png(paste0(out_dir,'/umap_tds_sm2.png'), width=6, height=8, units="in",res=300)
print(gg)
dev.off()

pdf(paste0(out_dir,'/umap_tds_sm2.pdf'), width=6, height=8)
print(gg)
dev.off()

#UMAP individual TDS OC genes----
thy.combined@meta.data %>% nrow() #17412
ff4.tc.combined <- subset(thy.combined, subset = cancer_type != 'normal')
ff4.tc.combined@meta.data %>% nrow() #14208
ff4.tc.combined@meta.data$cancer_type %>% unique()
DefaultAssay(ff4.tc.combined) <- "RNA"
FeaturePlot(ff4.tc.combined, features=tds_genes, slot='data')
png(paste0(out_dir,'/umap_tds_sm_each.png'), width=12, height=6, units="in",res=300)
FeaturePlot(ff4.tc.combined, features=c("TG","TPO","TSHR","PAX8",
                                        "SHMT2","MTHFD2","PHGDH","SHMT1"), 
            slot='data', ncol=4)
dev.off()

pdf(paste0(out_dir,'/umap_tds_sm_each.pdf'), width=12, height=6)
FeaturePlot(ff4.tc.combined, features=c("TG","TPO","TSHR","PAX8",
                                        "SHMT2","MTHFD2","PHGDH","SHMT1"), 
            slot='data', ncol=4)
dev.off()


#Grouping before DEG and GSEA -----
cfl_cpm_dt <- fl_cpm_dt %>% filter(cancer_type != 'NT')
med_tds_score <- cfl_cpm_dt$tds_score %>% median()
cfl_cpm_dt <- cfl_cpm_dt %>% mutate(tds_sm_group = ifelse(tds_score > med_tds_score, 
                                                          ifelse(SM_mean > 0, 'T_hi_SM_hi', 'T_hi_SM_lo'),
                                                          ifelse(SM_mean > 0, 'T_lo_SM_hi','T_lo_SM_lo')))

cfl_cpm_dt %>% group_by(tds_sm_group) %>% dplyr::count()
cfl_cpm_dt %>% select(cell_id, tds_sm_group) %>%write_tsv(paste0(out_dir,'/tds_sm_group.tsv'))
tmp_meta <- thy.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% as_tibble()
thy.combined@meta.data <- left_join(tmp_meta, cfl_cpm_dt %>% select(cell_id, tds_sm_group)) %>% as.data.frame() %>% column_to_rownames('cell_id')
DimPlot(thy.combined, group.by="tds_sm_group", reduction='umap')
ff4.tc.combined <- subset(thy.combined, subset = cancer_type != 'normal')

umap_dt <- ff4.tc.combined@reductions$umap@cell.embeddings
umap_dt <- umap_dt %>% as.data.frame %>% rownames_to_column('cell_id')
dim(umap_dt)
tmp_meta <- ff4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id') %>% as_tibble()
umap_dt <- left_join(umap_dt, tmp_meta %>% select(cell_id, tds_sm_group))

# DEG between TDS hi SM low vs TDS low SM high
Idents(ff4.tc.combined) <- 'tds_sm_group'
g1=c('T_hi_SM_lo')
g2=c('T_lo_SM_hi')  
DefaultAssay(ff4.tc.combined) <- "RNA"
mk <- FindMarkers(ff4.tc.combined, ident.1 = g1, ident.2=g2, min.pct = 0.25, 
                  max.cells.per.ident=2000, random.seed=1,
                  logfc.threshold = 0.15)
if(T){
  #volcano-----
  mk_dt <- mk %>% as.data.frame() %>% rownames_to_column('gene_name') %>% as_tibble()
  mk_dt <- mk_dt %>% mutate(nlpval = -log10(p_val_adj))
  mk_dt$nlpval[is.finite(mk_dt$nlpval) == F] <- 300
  RP_genes <- mk_dt %>% filter((grepl("^RPS",gene_name) | grepl("^RPL",gene_name))& p_val_adj < 0.05) %>% pull(gene_name)
  oxi_genes <- c("COX5A", "NDUFB9", "COX6B1", "UQCR10", "COX6A1", "COX7A2", 
                 "ATP5F1E", "MRPS15", "SHMT2", "GLRX", "NME2", "LGALS3", "SLC25A5",
                 "MTHFD2")
  
  gobp_names[grepl("_THYROID_HORMONE",gobp_names)]
  thyroid_genes <- unlist(gobp_list[gobp_names[grepl("_THYROID_HORMONE",gobp_names)]])
  
  
  mk_dt <- mk_dt %>% mutate(gene_name_group = ifelse(gene_name %in% RP_genes, 'RP_genes',
                                                     ifelse(gene_name %in% oxi_genes,'oxi_genes',
                                                            ifelse(gene_name %in% thyroid_genes,'thyroid_genes','others'))))
  gene_list1 <- mk_dt %>% filter(p_val_adj < 0.05) %>% arrange(avg_logFC) %>% head(n=10) %>% pull(gene_name)
  gene_list2 <- mk_dt %>% filter(p_val_adj < 0.05) %>% arrange(avg_logFC) %>% tail(n=10) %>% pull(gene_name)
  gene_list3 <- mk_dt %>% filter(p_val_adj < 0.05 & avg_logFC > 0) %>% arrange(p_val_adj) %>% head(n=10) %>% pull(gene_name)
  gene_list4 <- mk_dt %>% filter(p_val_adj < 0.05 & avg_logFC < 0) %>% arrange(p_val_adj) %>% head(n=10) %>% pull(gene_name)
  gene_name_list <- unique(union(union(union(gene_list1, gene_list2), gene_list3), gene_list4))
  g <- ggplot(mk_dt, aes(x=avg_logFC,y=nlpval))+
    geom_point(alpha=0.7)+
    geom_text_repel(data=subset(mk_dt, gene_name %in% gene_name_list), aes(label=gene_name), max.overlaps = 100)+
    geom_hline(yintercept=-log10(0.05), color="red")+
    labs(x="Log fold change", y="-Log10(adjusted P value)")+
    #geom_text_repel(data=subset(mk_dt, gene_name_group != 'others'),aes(label=gene_name, color=gene_name_group), max.overlaps = 100)+
    #xlab(paste0(paste(g2,collapse=','),'<----->',paste(g1, collapse=',')))+
    theme_syp
  png(paste0(out_dir,'/volcano_tds_sm_subgroup.png'), width=12, height=8, units="in",res=300)
  print(g)
  dev.off()
  pdf(paste0(out_dir,'/volcano_tds_sm_subgroup.pdf'), width=12, height=8)
  print(g)
  dev.off()
  
  #GSEA-----
  mk_dt <- mk_dt %>% mutate(nlpval = -log10(p_val_adj))
  mk_dt$nlpval[is.finite(mk_dt$nlpval) == F] <- 300
  mk_dt <- mk_dt %>% mutate(dist= sqrt(avg_logFC^2 + nlpval ^2))
  mk_dt <-mk_dt %>% mutate(dist = ifelse(avg_logFC < 0, (-1)*dist, dist))
  dds_stats <- mk_dt$dist
  names(dds_stats) <- mk_dt$gene_name
  gsea_res <- fgsea::fgsea(pathways = gobp_list, stats = dds_stats, nper=1000)
  gsea_res_collapse <- fgsea::collapsePathways(gsea_res, gobp_list,dds_stats)
  #gsea_res <- fgsea(pathways = gocc_list, stats = dds_stats, nper=1000)
  #gsea_res <- fgsea(pathways = h_list, stats = dds_stats, nper=1000)
  #gsea_res <- fgsea(pathways = kegg_list, stats = dds_stats, nper=1000)
  
  m_res <- gsea_res %>% as_tibble()
  #with collapse
  f_res <- m_res %>% filter(pval < 0.05 & pathway %in% gsea_res_collapse$parentPathways ) 
  f_res$rep_genes <-map_chr(f_res$leadingEdge, function(x) paste(x[1:5], collapse=','))
  x_order <- f_res %>% arrange(pval) %>% head(n=20) %>% arrange(NES) %>% pull(pathway)
  #without collapse
  if(F){
    f_res <- m_res %>% filter(pval < 0.05) 
    x1 <- f_res %>% arrange(NES) %>% head(n=10) %>% pull(pathway)
    x2 <- f_res %>% arrange(NES) %>% tail(n=10) %>% pull(pathway)
    x_order = unique(c(x1,x2))
  }

  g <- ggplot(f_res, aes(x=pathway, y=NES, fill= pval))+
    geom_bar(stat='identity')+
    #geom_text(aes(y=5, label=rep_genes), hjust=0)+
    #scale_y_continuous(limits=c(-5,15))+
    scale_x_discrete(limits = x_order)+
    coord_flip()+
    labs(y="Normalized enrichment score", x="")+
    #ylab(paste0(paste(g2,collapse=','),'<----->',paste(g1, collapse=',')))+
    theme_syp
  png(paste0(out_dir,'/gsea_tds_sm_subgroup.png'), width=16, height=8, units="in",res=300)
  print(g)
  dev.off()
  pdf(paste0(out_dir,'/gsea_tds_sm_subgroup.pdf'), width=16, height=8)
  print(g)
  dev.off()
  
  #GSEA enrichment plot-----
  f_res %>% filter(pathway %in% c("GOBP_THYROID_HORMONE_GENERATION",
                                  "GOBP_OXIDATIVE_PHOSPHORYLATION",
                                  "GOBP_MITOCHONDRIAL_MEMBRANE_ORGANIZATION"))
 
  g1<- fgsea::plotEnrichment(gobp_list[["GOBP_THYROID_HORMONE_GENERATION"]],
                        dds_stats) + labs(title="Thyroid hormone generation")

  g2 <- fgsea::plotEnrichment(gobp_list[["GOBP_OXIDATIVE_PHOSPHORYLATION"]],
                        dds_stats) + labs(title="Oxidative phosphorylation")

  g3 <- fgsea::plotEnrichment(gobp_list[["GOBP_MITOCHONDRIAL_MEMBRANE_ORGANIZATION"]],
                        dds_stats) + labs(title="Mitochondrail membrane organization")
  png(paste0(out_dir,'/gsea_enrichment.png'), width=9, height=12, units="in",res=300)
  plot_grid(g1,g2,g3, ncol=1) %>% print()
  dev.off()
  pdf(paste0(out_dir,'/gsea_enrichment.png'), width=9, height=12)
  plot_grid(g1,g2,g3, ncol=1) %>% print()
  dev.off()
}

### plot other thyroid genes 
if(T){
  tmp_mx <- l_cpm_mx[c(oc_gene_list, tds_genes,"PTPRC","DCN","EPCAM","KRT18","KRT19"),] %>% as.matrix()
  tmp_dt <- tmp_mx %>% t() %>% as.data.frame() %>% rownames_to_column('cell_id') %>% as_tibble()
  umap_dt <- thy.combined@reductions$umap@cell.embeddings
  umap_dt <- umap_dt %>% as.data.frame %>% rownames_to_column('cell_id')
  umap_dt <- left_join(umap_dt, tmp_meta %>% select(cell_id, tds_score, SM_mean))
  umap_dt <- left_join(umap_dt, tmp_dt)
  
  g1 <- draw_umap_color_plot(umap_dt, "TG")
  g2 <- draw_umap_color_plot(umap_dt, "TPO")
  g3 <- draw_umap_color_plot(umap_dt, "PAX8")
  g4 <- draw_umap_color_plot(umap_dt, "TSHR")
  g5 <- draw_umap_color_plot(umap_dt, "DUOX1")
  g6 <- draw_umap_color_plot(umap_dt, "EPCAM")
  g7 <- draw_umap_color_plot(umap_dt, "KRT18")
  g8 <- draw_umap_color_plot(umap_dt, "KRT19")
  
  gg <- plot_grid(g1,g2,g3,g4,g5,g6,g7,g8, nrow=2)
  png(paste0(out_dir,'/umap_other_genes.png'), width=20, height=8, units="in",res=300)
  print(gg)
  dev.off()
  
}


