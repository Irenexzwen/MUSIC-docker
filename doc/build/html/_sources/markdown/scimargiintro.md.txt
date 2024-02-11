# sciMARGI Introduction
## Motivation of the technology

Nucleus RNA is a vital component of the cell nucleus and plays a crucial role in genome organization and gene regulation. However, current technologies for profiling RNA-chromatin interactions can only be applied to large numbers of cells, which limits our understanding of cellular heterogeneity. By detecting RNA-DNA interactions at the single cell level, researchers can gain a deeper understanding of the unique regulatory landscape of nucleus RNA in individual cells and how it changes under different conditions. This is especially important in the study of complex tissues and disease states.


## Experimental design

The sciMARGI method is a technique for capturing multiway nucleotide interactions within the nucleus at single cell resolution. Single cell resolution is achieved through three rounds of split-pool ligation with different DNA barcodes, resulting in 96^3 (884,7367) total combinations. Nucleotide interactions are captured through crosslinking and digestion, and RNA and DNA ends are differentiated and ligated using specialized adaptors. Crosslinked molecular complexes are then physically separated and barcoded using the 10x droplet platform, with the complexes separated into numerous droplets, each containing one of 3.5 million barcodes. The final library is then split into 8 aliquots before sequencing to further improve the detection resolution at the single molecule level.

<figure><img src="https://lh5.googleusercontent.com/Clm82Am_hyAswqVtbiTKiqkKf1VYwnPEiDLXboS741-UVy_tk5ZBgfZgQvPvGAYKb7-pZIxErR9Tl5oY_awTRkA3gbMsVsV6YCs3BFlAWO2qRFShHcilFA7QNAHzYAN12OLSMm7GK6EEIxFpjJUxOXYbokDe38vfzxiKxVxb--7bwnmPxT1wYaPBMd6m1R5mpzTiLA9fgg" alt=""><figcaption><p>Figure 1: sciMARGI experimental workflow.</p></figcaption></figure>

## Library configuration

The experimental protocols and sequencing library designs is described in the previous chapter.  The sequence configuration is illustrated in Figure2. We adopted a pair-end sequencing strategy. Read1 is 28bp which records 16bp 10x barcode and 12bp UMI sequences. Read2 is 150bp which records three cell barcodes, adaptor sequence and insert sequence. Based on the design, RNA and DNA libraries can be pooled to sequence in the same run, and then be distinguished by the adaptor sequence. We use the 19bp unique sequence (the leftmost 19bp before NNNNN) with up to 2 mismatches to select DNA molecules. The RNA adaptor has a 15bp unique sequence interspersed with 5bp random umi, where we use all 20bp to select RNA molecules by allowing at most 2 mismatches. Given the length of either human genome (3 billion) or mouse genome (2.5 billion), the chance of seeing a random 19mer with the same sequence as a DNA adaptor is around 1%.&#x20;

Based on the experimental design, insertions (regardless RNA or DNA) with the same three cell barcodes are from the same single cell and those with same cell barcodes together with the same 10x barcodes, lib indexes are from the same molecular cluster. Each cell barcode is a random sequence from a set of 96. Each 10x barcode is a random sequence from a set of around 3.6 million. In theory, our design would be able to discriminate single DNA or RNA interactions from 884,736 cells and 28 million single molecules per cell. By mapping the DNA and RNA insert sequence to the reference genome would reveal molecular interactions.



<figure><img src="https://lh4.googleusercontent.com/tVMbZXZB_MyjUP405cZYi2HFbk38YHuPG4MH6_w7l0zkre0Pl-apFSjfefl_gJdyLdqJZ30vpnH8IRteJKohQNf-X_hZsbGHPwTneJBD1Tr7nZ-9zn38WnGIa5JruE2UrMdnGZOu2-ZR0FQBac69r7s8m5PWj6GJEUQad5AdC6DZ6S26FjwkebvrrThIcg" alt=""><figcaption><p>Figure 2: Sequencing library configuration for sciMARGI</p></figcaption></figure>

### Terminologies:

* **CB**: cell barcode. It is a combination of BC3, BC2 and BC1.
* **CBMB**: molecular barcode. Since the complete molecular barcode is composed by cell barcode plus 10x barcode, we use CBMB to indicate molecular barcode.
* **Libs:** Lib index. Each lib is one of the eight aliquots.




