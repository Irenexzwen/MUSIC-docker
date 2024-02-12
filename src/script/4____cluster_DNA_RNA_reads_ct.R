# 4____cluster_DNA_RNA_reads_ct.R

require(tidyverse)
require(data.table)
require(magrittr)
require(tidyr)
options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)
input_dir <- args[1]


cluster_human_dna <- fread(paste0(input_dir,"merge_DNA.sort.clusters_size.csv")) %>% 
  tidyr::separate(., col = "V1", into = c("CB", "CBMB"), sep = ",") %>% rename(dna_reads = cluster_num_reads)

cluster_human_rna <- fread(paste0(input_dir,"merge_RNA.sort.clusters_size.csv")) %>% 
  tidyr::separate(., col = "V1", into = c("CB", "CBMB"), sep = ",") %>% rename(rna_reads = cluster_num_reads)

human_clusters <- cluster_human_dna %>% dplyr::select(-CB) %>% 
  full_join(cluster_human_rna %>% select(-CB), by=c("CBMB"="CBMB")) %>% 
  mutate(dna_reads = replace_na(dna_reads, 0), rna_reads = replace_na(rna_reads, 0)) %>% 
  mutate(cluster_size = dna_reads + rna_reads) %>% 
  mutate(CB = gsub("^(BC3.*-BC2.*-BC1.*)-([ATCG]{16}-lib.*)$", "\\1",CBMB)) %>% 
  relocate(CB)

fwrite(human_clusters, paste0(input_dir, "merge.sort.cluster_rna_dna.csv"), quote = F, row.names = F)



