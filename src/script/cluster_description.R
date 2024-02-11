# 
MASTER_OUTDIR = "/dataOS/wenxingzhao/project/14____scMARGI/ZF____brain_tissue_lib3_20220822/outputs/"
source("/dataOS/wenxingzhao/project/Rproj/sc_SP_margi/sciMARGI_utils.R")

DD_col = "#5cc4ef"
RR_col = "#f3af42"
RD_col = "#bebdbe"

gtfgr36 <- readRDS("/dataOS/wenxingzhao/database/Human/genome/Robj/gtf36_gene_granges.rds")

human_clusters <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_human.sort.cluster_rna_dna.csv")) 

human_clusters %<>% mutate(type = ifelse(cluster_size==1, "singleton",
                                         ifelse(dna_reads>=2 & rna_reads==0, "DD", 
                                                ifelse(rna_reads>=2 & dna_reads==0, "RR", "RD"))))

cell_dna_human_reads <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_DNA_human.sort.cell_reads.csv")) %>% set_names(c("CB","human_dna_reads"))
cell_rna_human_reads <- fread(paste0(MASTER_OUTDIR, "6_stats/merge_RNA_human.sort.cell_reads.csv")) %>% set_names(c("CB","human_rna_reads"))

human_cell_reads <- cell_dna_human_reads %>% 
  left_join(cell_rna_human_reads, by=c("CB"="CB")) %>% 
  mutate(human_rna_reads = replace_na(human_rna_reads, 0)) %>% 
  mutate(total_reads = human_dna_reads + human_rna_reads)

human_tol_contacts <- human_clusters %>%
  mutate(DD_con = dna_reads*(dna_reads-1)/2, 
         RR_con = rna_reads*(rna_reads-1)/2, 
         RD_con = as.numeric(dna_reads)*as.numeric(rna_reads)) %>% 
  group_by(CB) %>% summarize(DD_contacts = sum(DD_con),
                             RR_contacts = sum(RR_con),
                             RD_contacts = sum(RD_con))

## ---- Per cell DD, RD, RR contacts | DNA reads, RNA reads | DD,RD,RR clusters ----------
human_clusters_summary <- human_clusters %>% 
  group_by(CB) %>% summarize(DD_num = sum(type=="DD"),
                             RD_num = sum(type=="RD"),
                             RR_num = sum(type=="RR")) %>% 
  mutate(total_clusters = DD_num+RD_num+RR_num) %>% 
  mutate(DD_prop = DD_num/total_clusters,
         RD_prop = RD_num/total_clusters,
         RR_prop = RR_num/total_clusters)


human_comb_df <- human_cell_reads %>% left_join(human_clusters_summary, by=("CB"="CB")) %>% 
  left_join(human_tol_contacts, by=c("CB"="CB")) %>% 
  dplyr::select(CB, human_dna_reads, human_rna_reads, DD_num, RD_num, RR_num,
                DD_contacts, RR_contacts, RD_contacts) 

h_cell_order <- human_cell_reads %>% dplyr::filter(human_dna_reads>10) %>% 
  arrange(desc(human_dna_reads)) %>% .[['CB']]

total_cell <- h_cell_order %>% length

x_labels_at <- quantile(1:total_cell, c(0, 0.33, 0.66, 1)) %>% as.integer()
x_labs = rep("", total_cell)
for(i in x_labels_at){
  x_labs[i] <- paste0("cell #",i)
}

g_all_h <- human_comb_df %>% dplyr::filter(CB %in% h_cell_order) %>% 
  melt() %>% mutate(group = ifelse(grepl("reads",variable), "reads", 
                                   ifelse(grepl("contacts", variable), "contact", "clusters"))) %>% 
  mutate(group=factor(group, levels=c("reads","clusters","contact"))) %>% 
  ggplot()+geom_point(aes(x=factor(CB, levels = h_cell_order),
                          y=value, color=variable, shape=group), alpha=0.5) + 
  facet_wrap(~group, ncol = 1, strip.position = "right", scales = "free_y")+
  scale_y_log10(labels = scales::scientific) + 
  scale_color_manual(labels = c("DNA reads", "RNA reads",
                                "DD clusters", "RD clusters", "RR cluster",
                                "DD contacts", "RD contacts", "RR contacts"),
                     values=c("#354d90", "#e8552b",DD_col, RR_col, RD_col, DD_col, RR_col, RD_col))+
  scale_x_discrete(labels=x_labs)+
  xlab("Ranked cells") + ylab("Number") + ggtitle("Human") +
  theme(axis.ticks.x = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=14), plot.title = element_text(hjust=0.5, size=14)) 

g_all_h


