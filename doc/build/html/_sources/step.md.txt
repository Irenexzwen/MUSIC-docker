# Data processing step-by-step

## Demultiplexing

<figure><img src="../.gitbook/assets/Library_configuration.png" alt=""><figcaption><p>MUSIC library configuration</p></figcaption></figure>

Demultiplexing is the step to extract the following important informations from the sequencing reads and output insertions with barcode, UMI information as the read name.&#x20;

| Name                             | Length      | Total combinations | Start location   |
| -------------------------------- | ----------- | ------------------ | ---------------- |
| Cell Barcode # 3                 | 14-17 bp    | 96                 | Read2 bp 1       |
| Cell Barcode # 2                 | 14bp        | 96                 | Read2 bp 21 - 24 |
| Cell Barcode # 1                 | 14bp        | 96                 | Read2 bp 52 - 55 |
| RNA/DNA adaptors                 | 20bp / 24bp | NA                 | Read2 bp 63 - 66 |
| Molecular Barcode (10x barcode)  | 16bp        | 3.5 million        | Read1 bp 1       |
| UMI                              | 12bp        | NA                 | Read1 bp 13      |

The cell barcode is designed such that the [levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein\_distance) between any two of the cell barcode is larger than 2. There is 7bp linker sequence between cell barcodes, and between cell barcodes and adaptor sequence. For cell barcode or adaptor sequences, we allow at most 2 mismatches to the reference sequence. For 10X barcode ([Chemistry: Single Cell 3' v1](https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-)), we only allow perfect match.



### Demultiplexing program

We used a customized python script to finish this step.&#x20;

| Type   | Description                                                                                                                                                                                                  |
| ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Tool   | 1\_\_\_\_combo\_barcode\_resolve.py                                                                                                                                                                          |
| Input  | Pairend fastq files in `gzip` format. R1 is 28bp and R2 is 150bp.                                                                                                                                            |
| Output | Two single end fastq files. One for DNA and the other is RNA. Each sequence is the DNA/RNA insertion sequence, with the read header contains `raw read name`\| `BC3`\_`BC2`\_`BC1` - `10x barcode` # `UMI`.  |

In the Snakemake pipeline, the command line is:

```
// Snakemake
python3 script/1____combo_barcode_resolve.py \
-R1 {input.all_read1} -R2 {input.all_read2} \
-lib {params.lib_name} \
-outRNA {output.outRNA_fq} \
-outDNA {output.outDNA_fq} \
-log {output.decode_log}
```



Sample Input:

```
# Input Read1 sequences (28bp):
@MN00185:343:000H532W5:1:11102:15153:1026 1:N:0:TCGTCTGA
NTCGTGAGTCACTTCCTGAGCTTAATGG
+
#AFFFFFFFFFFFFFFFFFFFFFFFFFF
@MN00185:343:000H532W5:1:11102:15334:1027 1:N:0:TCGTCTGA
NCTATCGGTGTTAAAGGGTTCGGAGACT
+
#AFFFFFFFFFFFFFFFFFFFFFFFFFF


# Input Read2 sequences (150bp)
@MN00185:343:000H532W5:1:11102:15153:1026 2:N:0:TCGTCTGA
CNNTCCTANCTNNANGTNCTCGTACAATCNAANGCNANGTGCTCANAGNGNNTCCAANACGAGGANTACNNACAACNCACANTGNGNAGTCTNTANANNGTGCNNATNACTAANAAANNAAAAAANANANANAAAANNNNANNNANNAAN
+
A##FFFFF#FF##F#FF#FFFFFFFFFFF#FF#FF#F#FFFFFFF#/F#F##FFFFF#FFFFFFF#F/F##FFFFF#F/F=#FA#/#AFFAF#//#=##//FF##FF#///FA#=FF##AFF/FF#/#/#/#//FF####F###/##//#
@MN00185:343:000H532W5:1:11102:15334:1027 2:N:0:TCGTCTGA
CNNGTTAANACNNANATNTTGTTTTCCCCNGGNTANANAATTGCTNCANTNNTTTTANTAATTTGTATCNNTTCATNCCCANATNTNATTCCNCTNANNTTGTNNAGNTATATNAGGNNGTCTATNANTNTNATCTNNNNANNNANNACN
+
A##FFFFF#FF##F#FF#FFFFFFFFFFF#FF#FF#F#FFFFFFF#FF#F##FFFFF#FFFAFFFFFFF##FFFFF#/==F#AF#/#/FFF6#6F#/##/FFF##/F#/FFFF#/F=##/FFFAF#/#/#/#/FFA####/###/##//#
```

Sample output:

```
# DNA end decoded fastq
@MN00185:343:000H532W5:1:11102:15055:1047|BC3_75-BC2_51-BC1_11-CGGAACCAGGCATCGA-lib1#GGTCTGTACCTG
TCATAAAGAACTCTACCAACTTAAGAAGAAGAAAGCAAGCAACCCCATACAAAAATGGCGA
+
=F/FA=/F=/F=FF/FF///6/A=AAFFF/A/=//F///FF/FFFFF/FF///FFF/F///
@MN00185:343:000H532W5:1:11102:12653:1047|BC3_3-BC2_23-BC1_31-GAACTGTTCCAAGCAT-lib1#ACCATGCCATCA
TCAAATTCCACAAAAAGAGTGTTTCAAATCTGCTCTGTGACTAGACACTGTGCGTATCTATAAA
+
FFFFFFFFF6FAFFFFFFFFFFFFF=FFFFFFF=FFFFFFFFF//F/FFFFFFFAAFA=F/F6/


# RNA end decoded fastq
@MN00185:343:000H532W5:1:11102:14470:1060|BC3_68-BC2_89-BC1_46-TGACGCGGTACACGTT-lib1#GTAAGGGTGAAG
AATGGCAGAATGTTTTGGTCAAAATTTTTTAGCCAGCACCCCAATTAAAAAAAAAAAAAAAAAA
+
FAFFFF==AFFFFFAFAFFAAF6=AFFFFFFA=F=/FAF=FFF=A=FFFF=FFFFF=/FF=//6
@MN00185:343:000H532W5:1:11102:23931:1064|BC3_47-BC2_8-BC1_20-CAATTTCCAACACGAG-lib1#ATGGCATAAAGT
AAGCTTTAGGAACAAAGAGCATGGATAGAGTTACGATAGAGTAGGTTAGCATAAAAACTACTGCC
+
F//=F=/AF//FFFFFFFF/FF/=F=A=A/=6F//F/FF//FA=///F/AAFFFF=F=FFFFF=F
```

Sample stats\_log:

```
// Some code
2023-01-02 02:09:58,776 INFO     ---- start ----
2023-01-02 02:10:02,958 INFO     ---- read in 10X BC ----
2023-01-02 02:10:09,158 INFO     ---- finish constructing barcode ----
2023-01-02 02:13:39,927 INFO     total lines in file /input/data/lib5_S5_R1_001.fastq.gz : 4016146
2023-01-02 02:13:39,929 INFO     incorrect 10x barcode (not perfect match) :134796
2023-01-02 02:13:39,929 INFO     reads incomplete CB :1153844
2023-01-02 02:13:39,929 INFO     Final DNA reads :2227173
2023-01-02 02:13:39,929 INFO     Final RNA reads :116670
2023-01-02 02:13:39,929 INFO     ---- finished ----
```

## Reads cleaning

Based on our experimental design, it's possible that we introduced some artificial sequences into the library even if the reads could have perfect matched three cell barcodes and 10x barcodes. We will remove these artifacts using `cutadapt`

### Reads cleaning program

In Snakemake pipeline we used the following command to remove known artifacts from the demultiplexed reads:

```
// For DNA reads
cutadapt -a 'A{{20}}' -a 'G{{20}}' -m 20 -j {threads} -o {output.trimed_DNA} {input.decode_DNA}

// For RNA reads
cutadapt -a CGAGGAGCGCTT -a 'A{{20}}N{{20}}' -q 15 -m 20 -j {threads} -o {output.trimed_RNA} {input.decode_RNA}
```

## Mapping DNA and RNA ends

After raw reads demultiplexing and reads cleaning, we will map DNA ends and RNA ends to reference genome individually. We chose bowtie2 for DNA ends genome mapping and STAR to map RNA ends to genome in a splice aware manner.&#x20;

Before mapping, we need to construct the indexes for both bowtie2 and STAR. The input file for the index construction is the reference fasta file, and gtf file. We chose to build index using snakemake pipeline instead of accepting pre-built index from users. We will sacrifice some time for building the index but will avoid the software compatibility headache.&#x20;

### Building index program

| Type   | Description                                                                                                 |
| ------ | ----------------------------------------------------------------------------------------------------------- |
| Input  | <ol><li>reference genome <code>fasta</code> file.</li><li>reference genome <code>gtf</code> file.</li></ol> |
| Output | <ol><li>bowtie2 index </li><li>STAR index</li></ol>                                                         |

The command to build index in Snakemake pipeline:

{% code lineNumbers="true" %}
```
// bowtie2
bowtie2-build --threads {threads} {input.fasta_file} {params.basename}

// STAR
STAR --runThreadN {threads} --runMode genomeGenerate --genomeFastaFiles {input.fasta_file} --sjdbGTFfile {input.gtf_file} --genomeDir {output.index_dir}
```
{% endcode %}



### Reads mapping program

| Type   | Description                                     |
| ------ | ----------------------------------------------- |
| Input  | Cleaned DNA/RNA reads in fastq format.          |
| Output | Sorted `bam` file for DNA and RNA respectively. |

The command in Snakemake pipeline:

{% code lineNumbers="true" %}
```
# DNA mapping and sorting using bowtie2
bowtie2 -p {threads} -t --phred33 -x {params.index_prefix} -U {input.DNA_fq} 2> {output.human_bowtie2_stats}| samtools view -bq 20 -F 4 -F 256 - | samtools sort -@ 10 -m 64G -O bam -o {output.DNA_human_bam}

# RNA reads mapping and bam file sorting
STAR --runThreadN {threads} --genomeDir {input.index_dir} --readFilesIn {input.non_rRNA_fq} --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0.5 --outFilterMatchNmin 15 --outFileNamePrefix {config[dir_names][mapped_RNA_dir]}{params}/{params}_human_
samtools view -bq 255 -F 4 -F 256 {config[dir_names][mapped_RNA_dir]}{params}/{params}_human_Aligned.sortedByCoord.out.bam | samtools sort -@ 10 -m 64G -O bam -o {output.STAR_human_bam}
samtools index {output.STAR_human_bam}
```
{% endcode %}

The reads mapping, bam file sorting and indexing are using multi threads to accelerate. User only need to specify the total cpu that is available and all details will be taken care of by Snakemake.&#x20;


## Reads deduplication

After reads mapping, the next step is to remove PCR duplicates. PCR duplicates are often recognized by having the same UMI (unique molecular identifier). Deduplication is often done after reads mapping as PCR duplicates are likely to map to the same location on the genome, at the same time, the sequencing errors in the reads are often taken care of by any alignment software.&#x20;

_Deduplicate reads based on the mapping co-ordinate and the UMI attached to the read_ is a hot topic and has been extensively discussed by many, including [umitools](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html). In our case, we wrote our own customized tool to solve this problem (given our UMI sequence is recorded in the read name that is very different from standard input).&#x20;

The idea is that we sort the bam file based on the mapping coordinates and we scan the bam file once. For any read that 1) maps to the location within 8bp from the previous record, 2) sharing the same cell barcode and molecular barcode, 3) and its UMI is less than 2bp levenshitein distance away from the previous read's UMI will be marked as a duplicate. All PCR duplicates are removed in the downstream analysis. &#x20;

### Dedup program

In Snakemake pipeline, the dedup step is as follows:

```
// Dedup step
python3 script/2____Remove_bam_duplicates.py -inbam {input.sorted_bam} -outbam {output.dedupped_bam}
```

We will dedup both DNA and RNA bam files.&#x20;

Sample Input (only showing the first 20 rows):

```
A00953:623:HHY3CDSX3:3:1150:16785:30874|BC3_60-BC2_10-BC1_31-CCGATCTCACATTCGA-lib1#TAGGGGGCTTTG	0	chr1	31681	22	63M	*	0	0	TCACTCCAGGCTGGGCGACAGAGAGAGACCCTGTCTCAGAAAAAAAAAAAAAAGTAATTTGTA	FFFFFF:FF:F,F,FFFFFFFF,FFFFFF,FF:F:FF:FFFF,FFFF:FFFFF:,:,::,F,F	AS:i:-8	XS:i:-31	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0G55C6	YT:Z:UU
A00953:623:HHY3CDSX3:3:2659:30490:37043|BC3_72-BC2_4-BC1_24-TCATGTTAGGACATCG-lib1#ATGATCTAGTCA	16	chr1	55802	24	63M	*	0	0	AGCAGGACTTTGTCGTGTTCACCTCTATCACATCATACATATAGCAAACAGTAAAACTATTTA	F,:FF::FFFFFF:FFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFF:FFFFFFF:FF,F	AS:i:-12	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:37A23G0C0	YT:Z:UU
A00953:623:HHY3CDSX3:3:1111:16459:1626|BC3_42-BC2_65-BC1_32-ATTACTCAGCGCTGAA-lib1#TTTTGGGATTAC	0	chr1	56275	30	64M	*	0	0	ATTCTCCCTCATTCATTTCAATACGCTGTTCGGCCTGCTACCCCAGTTTCCCACTTAGAACAAT	,:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,FFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:1349:7717:8625|BC3_42-BC2_65-BC1_32-ATTACTCAGCGCTGAA-lib1#TTTTGGGATTAC	0	chr1	56275	30	64M	*	0	0	ATTCTCCCTCATTCATTTCAATACGCTGTTCGGCCTGCTACCCCAGTTTCCCACTTAGAACAAT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:1520:12599:10128|BC3_42-BC2_65-BC1_32-ATTACTCAGCGCTGAA-lib1#TTTTGGGATTAC	0	chr1	56275	30	64M	*	0	0	ATTCTCCCTCATTCATTTCAATACGCTGTTCGGCCTGCTACCCCAGTTTCCCACTTAGAACAAT	,,FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:2561:10854:9330|BC3_42-BC2_65-BC1_32-ATTACTCAGCGCTGAA-lib1#TTTTGGGATTAC	0	chr1	56275	30	64M	*	0	0	ATTCTCCCTCATTCATTTCAATACGCTGTTCGGCCTGCTACCCCAGTTTCCCACTTAGAACAAT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:1410:31693:14544|BC3_52-BC2_92-BC1_41-ATTCCATGTCGATTTG-lib1#AAGATCAACCGT	0	chr1	61473	40	63M	*	0	0	TCAATTTGTTTACCATTATTACTCTTGGTATTTTTAAGAAAAGTCTTTCAATTGTTATTATAA	FFFFFFFFFFFF,FFF:FFFFFFFF:FFFFFFFFFFFFF:FFFFFF:FF:FFFFFFFFFFFFF	AS:i:-9	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0G48C13	YT:Z:UU
A00953:623:HHY3CDSX3:3:2101:15212:13902|BC3_52-BC2_92-BC1_41-ATTCCATGTCGATTTG-lib1#AAGATCAACCGT	0	chr1	61473	24	63M	*	0	0	TCAATTTGTTTACCATTATTACTCTTGGTATTTTTAAGAAAAGTCTGTCCATTGTTCTTATAA	FFFF,FFF,F:F:FF:F:FFFFF,FFFFF:F,FF:FFF,::FFFF:,::,:,,:F,:,FF:FF	AS:i:-12	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:0G45T9A6	YT:Z:UU
A00953:623:HHY3CDSX3:3:1507:31159:31344|BC3_37-BC2_69-BC1_1-TCGATTTAGAGGATCC-lib1#GTACAAATTGGC	0	chr1	62417	30	64M	*	0	0	TCACGTGTTGGTTTGGGGTAGATCATTATAGGCACATGTAGGAAACAGCTTTCAGAGATGCCTT	FFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFF::FFF:FFFF:F:	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:2233:19307:4586|BC3_37-BC2_69-BC1_1-TCGATTTAGAGGATCC-lib1#GTACAAATTGGC	0	chr1	62417	30	64M	*	0	0	TCACGTGTTGGTTTGGGGTAGATCATTATAGGCACATGTAGGAAACAGCTTTCAGAGATGCCTT	FFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFF:FF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:2410:24876:20776|BC3_37-BC2_69-BC1_1-TCGATTTAGAGGATCC-lib1#GTACAAATTGGC	0	chr1	62417	30	64M	*	0	0	TCACGTGTTGGTTTGGGGTAGATCATTATAGGCACATGTAGGAAACAGCTTTCAGAGATGCCTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:2526:13476:6637|BC3_37-BC2_69-BC1_1-TCGATTTAGAGGATCC-lib1#GTACAAATTGGC	0	chr1	62417	30	64M	*	0	0	TCACGTGTTGGTTTGGGGTAGATCATTATAGGCACATGTAGGAAACAGCTTTCAGAGATGCCTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFF:FFFFFFFFFFF:FFFFFFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:2560:15483:1564|BC3_37-BC2_69-BC1_1-TCGATTTAGAGGATCC-lib1#GTACAAATTGGC	0	chr1	62417	30	64M	*	0	0	TCACGTGTTGGTTTGGGGTAGATCATTATAGGCACATGTAGGAAACAGCTTTCAGAGATGCCTT	FFFFF,F,FFFF:,FF:F:FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:2408:20202:14340|BC3_4-BC2_37-BC1_28-TGTTTGTCAATAGTGA-lib1#AGACTCTAGCTT	16	chr1	64627	40	65M	*	0	0	TGCACATGTACCCTAGAACTTAAAGTATAATAAAAAAAAATAGACTCTAGGACTCTGTATTATGA	,FFFF,F,FFFF:FF:FFFFF:,FFFF,FF::FF:FFFFFFFFFFFFF,F,FFFF::FFFFF,FF	AS:i:-8	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:50T13C0	YT:Z:UU
A00953:623:HHY3CDSX3:3:2545:21377:2754|BC3_4-BC2_37-BC1_28-TGTTTGTCAATAGTGA-lib1#AGACTCTAGCTT	16	chr1	64627	40	65M	*	0	0	TGCACATGTACCCTCGAACTTAAAGTATAATAAAAAAAAATAGACTCTAGGACTCTGTATTATGA	,FFFF,F:,::FF,,,:FFF:,:F:FF::FFF:FFFFFFF:FFFF:FFF,,FFFF,FFF:FF::F	AS:i:-11	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:14A35T13C0	YT:Z:UU
A00953:623:HHY3CDSX3:3:1342:2555:31281|BC3_56-BC2_78-BC1_13-AGGTAGGGTCACTTCC-lib1#TAAAGCAACTTT	16	chr1	65955	30	63M	*	0	0	ATTTCTTTCACATAAAGGTACAAATCATACTGCTAGAGTTGTGAGGATTTTTACAGCTTTTGA	FFFFFF:FFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF,FFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:63	YT:Z:UU
A00953:623:HHY3CDSX3:3:1365:15872:24972|BC3_56-BC2_78-BC1_13-AGGTAGGGTCACTTCC-lib1#TAAAGCAACTTT	16	chr1	65955	30	63M	*	0	0	ATTTCTTTCACATAAAGGTACAAATCATACTGCTAGAGTTGTGAGGATTTTTACAGCTTTTGA	FFF:FFFF,F:FF,FFF:::FFFFFFFFFFFFFFFFFFFF::FFFF:FFFFFFFFFFFFFFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:63	YT:Z:UU
A00953:623:HHY3CDSX3:3:1375:25536:31469|BC3_56-BC2_78-BC1_13-AGGTAGGGTCACTTCC-lib1#TAAAGCAACTTT	16	chr1	65955	30	63M	*	0	0	ATTTCTTTCACATAAAGGTACAAATCATACTGCTAGAGTTGTGAGGATTTTTACAGCTTTTGA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:63	YT:Z:UU
A00953:623:HHY3CDSX3:3:1411:10059:26115|BC3_56-BC2_78-BC1_13-AGGTAGGGTCACTTCC-lib1#TAAAGCAACTTT	16	chr1	65955	30	63M	*	0	0	ATTTCTTTCACATAAAGGTACAAATCATACTGCTAGAGTTGTGAGGATTTTTACAGCTTTTGA	FFFFF,,FFFFF:,,FFFFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:63	YT:Z:UU
A00953:623:HHY3CDSX3:3:1423:7148:17879|BC3_56-BC2_78-BC1_13-AGGTAGGGTCACTTCC-lib1#TAAAGCAACTTT	16	chr1	65955	30	63M	*	0	0	ATTTCTTTCACATAAAGGTACAAATCATACTGCTAGAGTTGTGAGGATTTTTACAGCTTTTGA	FFF,FFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:63	YT:Z:UU
```

Sample output:

```
A00953:623:HHY3CDSX3:3:1150:16785:30874|BC3_60-BC2_10-BC1_31-CCGATCTCACATTCGA-lib1#TAGGGGGCTTTG	0	chr1	31681	22	63M	*	0	0	TCACTCCAGGCTGGGCGACAGAGAGAGACCCTGTCTCAGAAAAAAAAAAAAAAGTAATTTGTA	FFFFFF:FF:F,F,FFFFFFFF,FFFFFF,FF:F:FF:FFFF,FFFF:FFFFF:,:,::,F,F	AS:i:-8	XS:i:-31	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0G55C6	YT:Z:UU
A00953:623:HHY3CDSX3:3:2659:30490:37043|BC3_72-BC2_4-BC1_24-TCATGTTAGGACATCG-lib1#ATGATCTAGTCA	16	chr1	55802	24	63M	*	0	0	AGCAGGACTTTGTCGTGTTCACCTCTATCACATCATACATATAGCAAACAGTAAAACTATTTA	F,:FF::FFFFFF:FFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFF:FFFFFFF:FF,F	AS:i:-12	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:37A23G0C0	YT:Z:UU
A00953:623:HHY3CDSX3:3:1111:16459:1626|BC3_42-BC2_65-BC1_32-ATTACTCAGCGCTGAA-lib1#TTTTGGGATTAC	0	chr1	56275	30	64M	*	0	0	ATTCTCCCTCATTCATTTCAATACGCTGTTCGGCCTGCTACCCCAGTTTCCCACTTAGAACAAT	,:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,FFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:1410:31693:14544|BC3_52-BC2_92-BC1_41-ATTCCATGTCGATTTG-lib1#AAGATCAACCGT	0	chr1	61473	40	63M	*	0	0	TCAATTTGTTTACCATTATTACTCTTGGTATTTTTAAGAAAAGTCTTTCAATTGTTATTATAA	FFFFFFFFFFFF,FFF:FFFFFFFF:FFFFFFFFFFFFF:FFFFFF:FF:FFFFFFFFFFFFF	AS:i:-9	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0G48C13	YT:Z:UU
A00953:623:HHY3CDSX3:3:1507:31159:31344|BC3_37-BC2_69-BC1_1-TCGATTTAGAGGATCC-lib1#GTACAAATTGGC	0	chr1	62417	30	64M	*	0	0	TCACGTGTTGGTTTGGGGTAGATCATTATAGGCACATGTAGGAAACAGCTTTCAGAGATGCCTT	FFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFF::FFF:FFFF:F:	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU
A00953:623:HHY3CDSX3:3:2408:20202:14340|BC3_4-BC2_37-BC1_28-TGTTTGTCAATAGTGA-lib1#AGACTCTAGCTT	16	chr1	64627	40	65M	*	0	0	TGCACATGTACCCTAGAACTTAAAGTATAATAAAAAAAAATAGACTCTAGGACTCTGTATTATGA	,FFFF,F,FFFF:FF:FFFFF:,FFFF,FF::FF:FFFFFFFFFFFFF,F,FFFF::FFFFF,FF	AS:i:-8	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:50T13C0	YT:Z:UU
A00953:623:HHY3CDSX3:3:1342:2555:31281|BC3_56-BC2_78-BC1_13-AGGTAGGGTCACTTCC-lib1#TAAAGCAACTTT	16	chr1	65955	30	63M	*	0	0	ATTTCTTTCACATAAAGGTACAAATCATACTGCTAGAGTTGTGAGGATTTTTACAGCTTTTGA	FFFFFF:FFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF,FFFF	AS:i:0	XS:i:-5	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:63	YT:Z:UU
```

## Final output and summary report

### Merge bam files from individual library

Remember that the final library is splitted into 8 aliquots (8 libs) before sequencing to further increasing the combinations of molecular barocdes. As the PCR step is separately applied to each aliquots, our dedup pipeline is also applied individually.&#x20;

After these steps, we merged all 8 libs into one as our final output. The final out hence is in bam format that records the cell barcode, molecular barcode, insert mapping location for DNA and RNA individually (two separate files). Notice that the bam file is also sorted, so you can see reads from different libs can be seen in the final output file.&#x20;

```
>samtools view merge_DNA_human.sort.bam|head -n 10
A00953:623:HHY3CDSX3:3:1442:21585:7059|BC3_39-BC2_35-BC1_20-ACAACCATCCGGGACT-lib8#TCAAAAGTGCAC	0	chr1	17455	23	64M	*	0	0	TCACTGTGGGGTCCCAGGCCTCCCGGGCCGAGCCACCCGTCACCCCCTGGCTCCTGGCGGATGT	:FF:FF:,:,F:FFFFFFFFFF,F:,F:F,:FFF,FFFFFFFFF:FF,,FF:F:FF:F:FFF:,	AS:i:-16	XN:i:0	XM:i:4	XO:i:0	XG:i:0	NM:i:4	MD:Z:0G24A32C0T4	YT:Z:UU
A00953:623:HHY3CDSX3:3:2642:19705:24314|BC3_38-BC2_39-BC1_50-AAACGAACATACAGGG-lib2#AAGATAAATACG	16	chr1	31533	40	64M	*	0	0	TATTGTGAGATTCCATCTCTACAAAAATAAAATTAAATAGCCAGTCATGGTGTCACACACCTGA	FFF,F,FFFF,FFFFFFF:FFFF:FFF,,FF:F,F,FFFF,:FF:FFFFFFFFFFFFFFF:,FF	AS:i:-8	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:33A29T0	YT:Z:UU
A00953:623:HHY3CDSX3:3:1150:16785:30874|BC3_60-BC2_10-BC1_31-CCGATCTCACATTCGA-lib1#TAGGGGGCTTTG	0	chr1	31681	22	63M	*	0	0	TCACTCCAGGCTGGGCGACAGAGAGAGACCCTGTCTCAGAAAAAAAAAAAAAAGTAATTTGTA	FFFFFF:FF:F,F,FFFFFFFF,FFFFFF,FF:F:FF:FFFF,FFFF:FFFFF:,:,::,F,F	AS:i:-8	XS:i:-31	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0G55C6	YT:Z:UU
A00953:623:HHY3CDSX3:3:2504:26106:4993|BC3_2-BC2_25-BC1_2-AGGTTGTTCTGGACCG-lib2#AAGTCATACTAG	16	chr1	33708	22	65M	*	0	0	GATGTTCTCTGCCCCACGTCCAAGCGTTCTCATTGTTCAATTCCCACCTGTGAGTGAGAACATGA	FFFFFF:FF:FFF::FFF::F,FF:FFFFFFFFFFFFFFF,FFF:FFFFFFFFFFFF:FFF:FFF	AS:i:-10	XS:i:-37	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:16T47C0	YT:Z:UU
A00953:623:HHY3CDSX3:3:1572:6524:34585|BC3_69-BC2_36-BC1_48-AGAACCTCATTGCTGA-lib5#GATACTTGATTG	16	chr1	33710	22	63M	*	0	0	TGTTCTCTGCCACATGTCCAAGCGTTCTCATTGTTCAATTCCCACCTGTGAGTGAGAACATGA	FFFFFF:FFFF,F:FFFF:FFFFF:FFFFF,FFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFF	AS:i:-8	XS:i:-32	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:11C50C0	YT:Z:UU
A00953:623:HHY3CDSX3:3:1147:4933:22842|BC3_89-BC2_87-BC1_40-AACAGGGAGTTACGTC-lib7#TAAGGTTTTTTG	0	chr1	35357	40	62M	*	0	0	TCCCAAGGTCGGCCAGGTGCAGTGGCTCATGTCTATAATCCCAACACTTTGGGAGGCTAAGG	F:FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF	AS:i:-10	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0G57G3	YT:Z:UU
A00953:623:HHY3CDSX3:3:2639:25201:8030|BC3_85-BC2_25-BC1_59-ACTTCCGCAACACACT-lib4#TTGCACCAGACC	0	chr1	46640	40	62M	*	0	0	TCCTTGAAAGTCTTAACAATTTTTTTAACCAAAGTCCTCACAAAGTCAGTTTACATTAGCCC	,:FFF:F,F,,,F,F,FF:FFFFF,F:FFF,FFFFF:FFF::FF,F:FFFFFFF,FFF:FF:	AS:i:-9	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:0G8T34T17	YT:Z:UU
A00953:623:HHY3CDSX3:3:2515:30273:18756|BC3_69-BC2_90-BC1_25-GAGCCTGAGCACACAG-lib3#AATCTCCACAGC	16	chr1	46642	24	63M	*	0	0	CTTGAAATTATTAAAAATTTTTTTAACCAAAGTCCTCACAAATTAAGTTTACATTAGCCCTGA	FF:,:::F:,,FFF,,F:FFF,FFFF,,FF::F,F:F:,,,,FF,F,FFF:F,FFFF:,,F,F	AS:i:-14	XN:i:0	XM:i:4	XO:i:0	XG:i:0	NM:i:4	MD:Z:9C4C29C17C0	YT:Z:UU
A00953:623:HHY3CDSX3:3:2426:14470:19758|BC3_28-BC2_2-BC1_55-GATCAGTTCACTCACC-lib2#ATCCTTCATGTG	0	chr1	49050	22	64M	*	0	0	TCATGAGAATGTTCATTGCAGCACTCCTCACAATAGCAGAGACATGGAATCAACTTAAATGCCC	FFF,FFFFFFFFFFF:FFFFFFFFF,F,FFF:FF,:FFFFFFFFFFFFFFFFFFF,FF:FFFFF	AS:i:-8	XS:i:-33	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0G24A38	YT:Z:UU
A00953:623:HHY3CDSX3:3:2133:13548:8046|BC3_51-BC2_11-BC1_32-AGGACTTAGTCTTCCC-lib3#TAATTGCCAATA	16	chr1	51853	23	12M1D51M	*	0	0	AGACTCCGCCTCAAAAAAAAAAAAGAAGATTGATCCGAGAGTACCTCCCCTAAGGGTACATGA	F:F:FFFF,:FFFF::FFFF:FFFFF:F,:FF:FF,:F:::,FF:FF,FF,,FFFFFFFFFFF	AS:i:-16	XN:i:0	XM:i:2	XO:i:1	XG:i:1	NM:i:3	MD:Z:12^A23A26C0	YT:Z:UU
```

The bam output is compressed and can be easily imported and processed by widely used packages in `R` or `python`

### `Summary report`

MUSIC is designed as a single-cell single molecule method. It is highly useful for users to get an idea of the per cell reads, per molecular cluster reads quantity and distribution.  To this end, we derived the following summary files:

| File names                                     | Description                                                                                                                                                                                                                                |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| merge\_DNA(RNA)\_human.sort.cell\_clusters.csv | Number of DNA(RNA) clusters each cell have. Two columns record cell barcode and number of cluster.                                                                                                                                         |
| merge\_DNA(RNA)\_human.sort.cell\_reads.csv    | Number of DNA(RNA) reads each cell have. Two columns record cell barcode and number of cluster.                                                                                                                                            |
| merge\_DNA(RNA)\_human.sort.clusters\_size.csv | Number of DNA(RNA) reads in each molecular complex. Three columns: 1) cell barcode, 2) molecular barcode, 3) number of DNA reads in this barcode complex.                                                                                  |
| merge\_human.sort.cluster\_rna\_dna.csv        | Number of DNA and RNA reads in each molecular complex. Five columns: 1) Cell barcode, 2) molecular barcode, 3) number of DNA reads in this molecular complex, 4) number of RNA reads in this molecular complex. 5) Total reads (DNA+RNA).  |


