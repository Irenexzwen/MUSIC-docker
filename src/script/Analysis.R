# Analysis.R

library(Rsamtools)
require(tidyverse)
require(data.table)
require(magrittr)
require(stringi)
require(seqinr)
require(tidyr)
require(chromstaR)
require(plyranges)
require(ggpubr)
options(stringsAsFactors = FALSE)

source("/dataOS/wenxingzhao/project/Rproj/sc_SP_margi/combo_snakemake/script/sciMARGI_utils.R")

MASTER_OUTDIR = "/dataOS/wenxingzhao/project/14____scMARGI/ZF____H1_E14_mix_sciMARGI_lib3_20220702/outputs/"



## ---- mixed species studies ----------
cell_dna_human_reads <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_DNA_human.sort.cell_reads.csv")) %>% set_names(c("CB","human_dna_reads"))
cell_dna_mouse_reads <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_DNA_mouse.sort.cell_reads.csv")) %>% set_names(c("CB","mouse_dna_reads"))

# 95% separation
navy_blue = "#3c5870"
rose_red = "#a5565c"
grey = "grey80"
lemon_yellow = "#f4d15c"

purity_standard = 0.9
low_threshold = 10

cell_purity <- cell_dna_human_reads %>% full_join(cell_dna_mouse_reads, by=c("CB"="CB")) %>% 
  mutate(human_dna_reads = replace_na(human_dna_reads, 0)) %>% 
  mutate(mouse_dna_reads = replace_na(mouse_dna_reads, 0)) %>% 
  dplyr::mutate(human_ratio = human_dna_reads/(human_dna_reads+mouse_dna_reads)) %>% 
  dplyr::mutate(mouse_ratio = mouse_dna_reads / (mouse_dna_reads+human_dna_reads)) %>% 
  rowwise() %>% 
  dplyr::mutate(species = ifelse(max(human_dna_reads, mouse_dna_reads)<low_threshold, "low",
                          ifelse(human_ratio>=purity_standard, "human",
                          ifelse(mouse_ratio>=purity_standard, "mouse", "mix")))) %>% 
  dplyr::mutate(colors_ = case_when(species=="ambient" ~ "grey80",
                                    species=="human" ~ rose_red,
                                    species=="mouse" ~ navy_blue,
                                    species=="mix" ~ lemon_yellow))
  
perc_labels <- cell_purity$species %>% janitor::tabyl() %>% 
  rename(.,Type = `.`) %>% janitor::adorn_pct_formatting(digits = 2, ) %>% 
  mutate(label = paste0(Type,":",n," (", percent,")")) %>% 
  column_to_rownames("Type")


g_mix_dna <- ggplot(cell_purity)+geom_point(aes(x=human_dna_reads, y=mouse_dna_reads, color = species), 
                    size=2, alpha=0.8)+ 
  xlab("Single cell DNA reads\nmapped to human genome") + ylab("Single cell DNA reads\nmapped to mouse genome") + #xlim(c(0,850))+
  ggtitle(paste0("Cell clumping using DNA reads\ncell purity (",scales::percent(purity_standard),")")) +
  scale_color_manual(values = c("human"=rose_red, "mouse"=navy_blue, "mix"=lemon_yellow, "low"="grey80"),
                     labels = c(perc_labels["human","label"], perc_labels["mouse","label"],
                                perc_labels["mix","label"], perc_labels["low","label"]),
                     name = "Type")+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill="white",color="black",size=1),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=14),
    plot.title = element_text(size=16, hjust=0.5),
    legend.position = c(0.7, 0.8),
    legend.title = element_text(size=14), legend.text = element_text(size=12)
  ) 


## RNA ends
cell_rna_human_reads <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_RNA_human.sort.cell_reads.csv")) %>% set_names(c("CB","human_rna_reads"))
cell_rna_mouse_reads <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_RNA_mouse.sort.cell_reads.csv")) %>% set_names(c("CB","mouse_rna_reads"))

g_mix_RNA <- cell_rna_human_reads %>% full_join(cell_rna_mouse_reads, by=c("CB"="CB")) %>% 
  mutate(human_rna_reads = replace_na(human_rna_reads, 0)) %>% 
  mutate(mouse_rna_reads = replace_na(mouse_rna_reads, 0)) %>% 
  inner_join(cell_purity %>% select(CB, species)) %>% 
  
  ggplot()+geom_point(aes(x=human_rna_reads, y=mouse_rna_reads, color = species), 
                                 size=2, alpha=0.8)+ 
  xlab("Single cell RNA reads\nmapped to human genome") + ylab("Single cell RNA reads\nmapped to mouse genome") + #xlim(c(0,850))+
  ggtitle("RNA reads per cell") +
  scale_color_manual(values = c("human"=rose_red, "mouse"=navy_blue, "mix"=lemon_yellow, "low"="grey80"))+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill="white",color="black",size=1),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=14),
    plot.title = element_text(size=16, hjust=0.5),
    legend.position = c(0.7, 0.8),
    legend.title = element_text(size=14), legend.text = element_text(size=12)
  ) 

require(patchwork)
g_mix_DNA | g_mix_RNA
  
## ---- cluster size distribution barplot----------

human_cells <- cell_purity %>% filter(species == "human") %>% .[['CB']]
mouse_cells <- cell_purity %>% filter(species == "mouse") %>% .[['CB']]

human_clusters <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_human.sort.cluster_rna_dna.csv")) 
mouse_clusters <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_mouse.sort.cluster_rna_dna.csv")) 


plot_clu_size_group <- function(species_cluster){
  
  clu_size_df <- 
    species_cluster %>% 
    mutate(cluster_size_group = cut(cluster_size, c(0,1,10,100,1000,Inf))) %>% 
    .[['cluster_size_group']] %>% 
    janitor::tabyl() %>% 
    janitor::adorn_pct_formatting(digits = 2) %>% 
    set_colnames(c("cluster_size_group","group_size","lab"))
  
  ggplot(clu_size_df, aes(x=cluster_size_group, y=group_size, fill=cluster_size_group,
                          label = lab))+geom_bar(stat = 'identity')+
    scale_y_log10() +
    annotate("text", x=1:5, y=clu_size_df$group_size*1.3, label=clu_size_df$lab, size=4)+
    xlab("Cluster size group") + ylab("Number of clusters") +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill="white"),
      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
      axis.text.x = element_text(size=12, angle = 15, hjust=1),
      axis.text.y = element_text(size=12),
      axis.title = element_text(size=14), legend.position = "None"
    ) 
}


g_human_cluster_size_barplot <- plot_clu_size_group(human_clusters)
g_mouse_cluster_size_barplot <- plot_clu_size_group(mouse_clusters)


## ---- cell reads and clusters dot plot ----------

human_cluster_reads <- 
  human_clusters %>% dplyr::filter(CB %in% human_cells) %>% 
  dplyr::group_by(CB) %>%
  dplyr::summarize(total_clusters = length(CBMB), total_reads=sum(cluster_size))

mouse_cluster_reads <- 
  mouse_clusters %>% dplyr::filter(CB %in% mouse_cells) %>% 
  dplyr::group_by(CB) %>%
  dplyr::summarize(total_clusters = length(CBMB), total_reads=sum(cluster_size))


cluster_reads_per_cell_dotplot <- function(df){
  
  df %>% 
    arrange(total_reads) %>% mutate(ID = seq(1, nrow(.))) %>% 
    pivot_longer(., cols = c("total_clusters", "total_reads"), names_to = "num") %>% 
    ggscatter(x="ID", y="value", color='num', alpha=0.5)+ yscale('log10', .format = T) +
    scale_color_manual(name="Type", labels=c("Total DNA clusters",'Total DNA reads'),
                       values = c("#5a9554",'#9e778e'))+
    ggtitle("Human+Mouse: DNA+RNA")
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill="white",color="black",size=1),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      axis.title = element_text(size=14),
      plot.margin = unit(c(0.5,0.8,0.5,0.5), "cm"),
      legend.text = element_text(size=12), legend.position = c(0.3, 0.85),
      plot.title = element_text(size=14, hjust=0.5)
    ) +xlab("Cell Rank") +ylab("Total number per cell")
  
}

g_human_clusters_reads_dotplot <- cluster_reads_per_cell_dotplot(human_cluster_reads)
g_mouse_clusters_reads_dotplot <- cluster_reads_per_cell_dotplot(mouse_cluster_reads)



## ---- DNA-DNA contact map ----------

split_CB_MB <- function(DNA_human_in){
  CB <- gsub(".*\\|(BC3.*-BC2.*-BC1.*)-([ATCG]{16}-lib.*)#([ATCG]{12})$", "\\1",
             DNA_human_in$qname)
  CBMB <- gsub(".*\\|(BC3.*-BC2.*-BC1.*)#[ATCG]{12}$", "\\1",
               DNA_human_in$qname)
  DNA_human_in %<>% plyranges::mutate(CB = CB,CBMB = CBMB) %>% 
    plyranges::select(-mapq, -qname)
  
  return(DNA_human_in)
}

DNA_human <- 
  chromstaR::readBamFileAsGRanges(paste0(MASTER_OUTDIR, "5_merge/merge_DNA_human.sort.bam"),what=c("qname","mapq")) %>% 
  keepSeqlevels(., paste0("chr", c(1:22, 'X')), pruning.mode = "coarse") %>% 
  split_CB_MB %>% dplyr::filter(CB %in% human_cells)

RNA_human <- 
  chromstaR::readBamFileAsGRanges(paste0(MASTER_OUTDIR, "5_merge/merge_RNA_human.sort.bam"),what=c("qname","mapq")) %>% 
  keepSeqlevels(., paste0("chr", c(1:22, 'X')), pruning.mode = "coarse") %>% 
  split_CB_MB %>% dplyr::filter(CB %in% human_cells)

tiles <- gen_windows(window.size = 1e6, species = 'hg38')


### all cells: using different cluster size

plot_four_cluster_condition <- function(h_gr, tile){
  
  plots <- lapply(c("dna_reads>1 & dna_reads<=10",
                    "dna_reads>10 & dna_reads<=100",
                    "dna_reads>100 & dna_reads<=1000",
                    "dna_reads>1000"), function(cond){
                      
                      sub_clu_size_df <- human_clusters %>%  
                        dplyr::filter(eval(rlang::parse_expr(cond))) 
                      
                      test <- gen_DD_2d_contact(dna_gr = h_gr,
                                                clu_size_df = sub_clu_size_df,
                                                DNA_tiles = tile, 
                                                CB_selected=NULL, clu_selected=NULL,
                                                return_sparse=F, diag_rm = F)
                      
                      g <- single_hm_from_full_mtx(contact_matrix_full = test,log_values = T, legend_position = 'None')+
                        ggtitle(cond)+theme(plot.title = element_text(size=10, hjust=0.5))
                    })
  return(plots)
}

cluster_size_plots <- plot_four_cluster_condition(
                            h_gr = DNA_human,
                            tile = tiles %>% plyranges::filter(seqnames == "chr1"))

cowplot::plot_grid(plotlist = cluster_size_plots, nrow=1, byrow = T)



### DNA-DNA contact map on one chr using 10 - 100 - all cells 

plot_DD_by_cell_num <- function(n_rand_cells, chrs){
  
  rdn_cell_CB <- human_clusters %>% 
    filter(dna_reads>1 & dna_reads<1000) %>%
    .[['CB']] %>% unique %>% sample(., n_rand_cells)
  
  test <- gen_DD_2d_contact(dna_gr = DNA_human,
                            clu_size_df = human_clusters %>% filter(dna_reads>1, dna_reads<1000),
                            DNA_tiles = tiles %>% plyranges::filter(seqnames==chrs) , 
                            CB_selected=rdn_cell_CB, clu_selected=NULL,
                            return_sparse=F, diag_rm = F)
  
  g <- single_hm_from_full_mtx(contact_matrix_full = test, log_values = T, legend_position = 'None')+
    ggtitle(paste0(chrs,":", n_rand_cells, " cells")) + theme(plot.title = element_text(hjust=0.5))
  return(g)
}

set.seed(1227)
g_DD_cels <- 
map2(.x = expand.grid(c(10, 100, 500, 900), c("chr1","chr2","chrX"))$Var1,
     .y = expand.grid(c(10, 100, 500, 900), c("chr1","chr2","chrX"))$Var2 %>% as.character,
     function(x, y){plot_DD_by_cell_num(n_rand_cells = x, y)}) %>% 
  cowplot::plot_grid(plotlist = ., nrow=3, byrow = T)

g_DD_cels




### all cells, all chros ----------
allchrs <- gen_DD_2d_contact(dna_gr = DNA_human,
                          clu_size_df = human_clusters %>% filter(dna_reads>1, dna_reads<1000),
                          DNA_tiles = tiles %>% plyranges::filter(seqnames %ni% c("chrM", "chrY")), 
                          CB_selected=NULL, clu_selected=NULL,
                          return_sparse=F, diag_rm = F)

g_all_chrs_map <- single_hm_from_full_mtx(contact_matrix_full = allchrs, log_values = T, legend_position = 'None')+
  ggtitle("Contacts from all chromosomes") + theme(plot.title = element_text(hjust=0.5))



## ---- caRNA abundance ----------

# compare with H1 iMARGI

# transcription correlation with iMARGI at gene level
gtfgr36 <- readRDS("/dataOS/wenxingzhao/database/Human/genome/Robj/gtf36_gene_granges.rds")

ctrl_gi <- readRDS("/dataOS/wenxingzhao/project/13_rna-dna/IGM_iMARGI/MAPQ_30_bedpe/iMARGI_H1_control_merge.mapq30.1k_all_gi.rds")

sciMARGI_bulk_tx_gene <- countOverlaps(gtfgr36, RNA_human, ignore.strand=F)
iMARGI_tx_gene <- countOverlaps(gtfgr36, anchors(ctrl_gi)$first, ignore.strand=F)

# png("/dataOS/wenxingzhao/project/Rproj/sc_SP_margi/ZF____H1_E14_new_design_20220511/Figures/RNA_abundance_cor_with_iMARGI_gene_level.png",
#     width = 4, height = 3, units = 'in', res = 400)
ggscatter(data.frame(iMARGI = iMARGI_tx_gene, sciMARGI_ensemble = sciMARGI_bulk_tx_gene) %>% 
            dplyr::filter(iMARGI*sciMARGI_ensemble>0),
          x = "sciMARGI_ensemble", y="iMARGI", add = 'reg.line', title = "RNA abundance in each gene\n(reads)",
          conf.int = TRUE, alpha=0.5, color='grey10', cor.coef = T, cor.method = 'spearman')+
  xscale('log10', .format = T) + yscale('log10', .format = T) +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=14),
    plot.title = element_text(size=16, hjust=0.5)
  ) 
# dev.off()


## ---- RD interaction 2D map ----------


clu_2to1000 <- human_clusters %>% dplyr::filter(dna_reads>0 & rna_reads>0 & cluster_size<1000) %>% 
  .[['CBMB']]

plot_RD_2D_chrs <- function(chrs){
  rd_chr_mtx_full <- gen_RD_2d_contact(dna_gr = DNA_human, rna_gr = RNA_human,
                                       clu_size_df = human_clusters, clu_selected = clu_2to1000,
                                       RNA_tiles = tiles %>% plyranges::filter(seqnames==chrs),
                                       DNA_tiles = tiles %>% plyranges::filter(seqnames==chrs),
                                       return_sparse = F)
  
  g <- single_hm_from_full_mtx(contact_matrix_full = rd_chr_mtx_full, log_values = T, legend_position = 'None')+
    ggtitle(chrs)+theme(plot.title = element_text(size=12, hjust = 0.5))
  
  return(g)
}

g_RD_2D <- lapply(paste0("chr",1:4), plot_RD_2D_chrs) %>% 
  cowplot::plot_grid(plotlist = ., nrow = 1)


## ---- TAD are 2D confinement of DD RD ----------
