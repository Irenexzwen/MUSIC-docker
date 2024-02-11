# notice the processed_fastq has to be in the directory of /dataOS/wenxingzhao/Lab_sequencing/Lab_miniseq_processed/


# design 1 
ls /mnt/extraids/SDSC_NFS/wenxingzhao/Lab_sequencing/IGM/igm-storage2.ucsd.edu/221116_A00953_0655_AHLY3FDSX5/FASTQ/Brain_Lib_10_*|xargs -n2|sed 's/ /,/'|parallel "echo {}|awk -F'[/_]' '{print \$20\",\"\$0}'"|parallel "echo lib{}"|/dataOS/wenxingzhao/software/csvtk/csvtk add-header -n Sample,R1,R2|/dataOS/wenxingzhao/software/csvtk/csvtk csv2json -k Sample > samples.json

