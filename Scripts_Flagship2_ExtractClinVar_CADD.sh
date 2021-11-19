#!/bin/bash -l

####################
######CLINVAR#######
####################

####After downloading the ClinVar master summary file (clinvar_variant_summary.txt.gz), extract those variants of interest#######
########Pathogenic#####
zcat clinvar_variant_summary.txt.gz | awk -v FS="\t" '($7 ~ /^Pathogenic$/ && $17 == "GRCh38" && $2 == "single nucleotide variant") || $17 == "Assembly" {print}' > PathogenicOnly
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $19,$20-1,$21,$5,$7,$16}' PathogenicOnly | sed 's/^/chr/g' | sort -k1,1 -k2,2n | uniq | sed 's/, /_/g' | sed 's/ /_/g' > PathogenicOnly_4Intersect

########Likely pathogenic#####
zcat clinvar_variant_summary.txt.gz | awk -v FS="\t" '($7 ~ /^Likely pathogenic$/ && $17 == "GRCh38" && $2 == "single nucleotide variant") || $17 == "Assembly" {print}' > LikelyPathogenicOnly
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $19,$20-1,$21,$5,$7,$16}' LikelyPathogenicOnly | sed 's/^/chr/g' | sort -k1,1 -k2,2n | uniq | sed 's/, /_/g' | sed 's/ /_/g' > LikelyPathogenicOnly_4Intersect

########Benign#####
zcat clinvar_variant_summary.txt.gz | awk -v FS="\t" '($7 ~ /^Benign$/ && $17 == "GRCh38" && $2 == "single nucleotide variant") || $17 == "Assembly" {print}' > BenignOnly
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $19,$20-1,$21,$5,$7,$16}' BenignOnly | sed 's/^/chr/g' | sort -k1,1 -k2,2n | uniq | sed 's/, /_/g' | sed 's/ /_/g' > BenignOnly_4Intersect

########Likely benign#####
zcat clinvar_variant_summary.txt.gz | awk -v FS="\t" '($7 ~ /^Likely benign$/ && $17 == "GRCh38" && $2 == "single nucleotide variant") || $17 == "Assembly" {print}' > LikelyBenignOnly
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $19,$20-1,$21,$5,$7,$16}' LikelyBenignOnly | sed 's/^/chr/g' | sort -k1,1 -k2,2n | uniq | sed 's/, /_/g' | sed 's/ /_/g' > LikelyBenignOnly_4Intersect

#####concatenate all the classes######
cat PathogenicOnly_4Intersect LikelyPathogenicOnly_4Intersect BenignOnly_4Intersect LikelyBenignOnly_4Intersect | sort | uniq > ClinVar_AllClasses_4Intersect

#####Intersect with phyloP scores
bedtools intersect -a $PathToPhyloPScoresFolder/phyloP_scores.bed.gz -b ClinVar_AllClasses_4Intersect -wo | awk -v FS="\t" -v OFS="\t" '{print $6,$7,$8,$5,$9,$10,$11}' | sed '1ichr\tstart\tend\tPhyloP\tgene_name\tClass\tOrigin' > ClinVarHg38_AllClasses_IntersectedWithPhyloP_AllChromosomes.bed


####################
######CADD#######
####################

####After downloading the gnomAD - CADD  summary file (gnomad.genomes.r3.0.snv.tsv.gz), extract those biallelic variants#######
while read i
do
zcat gnomad.genomes.r3.0.snv.tsv.gz | awk -v var="$i" '$1 == var {print}' | sed 's/^/chr/g' | awk -v OFS="\t" '{print $1,$2,$6}' > gnomad_chr${i}_ToGrepWithin
zcat gnomad.genomes.r3.0.snv.tsv.gz | awk -v var="$i" '$1 == var {print}' | sed 's/^/chr/g' | awk -v OFS="\t" '{print $1,$2}' | uniq -u > gnomad_chr${i}_ToGrepFor
grep -wFf gnomad_chr${i}_ToGrepFor gnomad_chr${i}_ToGrepWithin | sed 's/$/\tNA/g' > gnomad_chr${i}_OnlyUniqueVarsBiallelic &&
rm gnomad_chr${i}_ToGrepFor gnomad_chr${i}_ToGrepWithin
done<Chromosome_list   ####Chromosome_list is the list of chromosomes

####Extract those multiallelic variants and get the mean and SD of the CADD score#######
while read i
do
zcat gnomad.genomes.r3.0.snv.tsv.gz | awk -v var="$i" '$1 == var {print}' | sed 's/^/chr/g' | awk -v OFS="\t" '{print $3,$4,$5,$6,$1,$2}' | uniq -f4 -D | awk -v OFS="\t" '{print $5,$6,$4}' > gnomad_chr${i}_DupVarsMultiallelic

R --vanilla <<EOF
library(dplyr)
File <- read.csv("gnomad_chr${i}_DupVarsMultiallelic", sep="\t", header=F)
a <- File %>% group_by(V1,V2) %>% summarise(mean = mean(V3), sd = sd(V3))
write.table(a, "gnomad_chr${i}_MultiallelicPosition_CADD_Mean_SD.txt", sep="\t", quote=F, row.names=F)
EOF

done<Chromosome_list   ####Chromosome_list is the list of chromosomes

##################################################################################################
#######With the above-described CADD scripts, we extract the CADD_Phred scores. To extract the CADD-raw scores (to be used for teh correlation analysis, we only need to specify a different column in the code##
###line 37: zcat gnomad.genomes.r3.0.snv.tsv.gz | awk -v var="$i" '$1 == var {print}' | sed 's/^/chr/g' | awk -v OFS="\t" '{print $1,$2,$5}' > gnomad_chr${i}_ToGrepWithin
###line 46: zcat gnomad.genomes.r3.0.snv.tsv.gz | awk -v var="$i" '$1 == var {print}' | sed 's/^/chr/g' | awk -v OFS="\t" '{print $3,$4,$5,$6,$1,$2}' | uniq -f4 -D | awk -v OFS="\t" '{print $5,$6,$3}' > gnomad_chr${i}_DupVarsMultiallelic
####THE FOLLOWING ALSO APPLY TO THIS SET############
##################################################################################################

#######Create files to be intersected with phyloP scores####
while read i
do
zcat gnomad.genomes.r3.0.snv.tsv.gz | awk -v var="$i" '$1 == var {print}' | sed 's/^/chr/g' | awk -v OFS="\t" '{print $1,$2-1,$2}' > gnomad_chr${i}_4Intersection
gzip gnomad_chr${i}_4Intersection
done<Chromosome_list   ####Chromosome_list is the list of chromosomes

#####Intersect with phyloP scores######
while read i
do
bedtools intersect -a $PathToPhyloPScoresFolder_PerChrom/chr${i}.bed.gz -b gnomad_chr${i}_4Intersection.gz -wb | awk -v OFS="\t" '{print $1,$3,$5}' | sort -k1,1 -k2,2n | uniq > gnomad_chr${i}_IntersectedWithPhyloP.bed
done<Chromosome_list   ####Chromosome_list is the list of chromosomes


######Merging gnomad bi- and multiallelic variants, and then with phyloP scores######
while read i
do
######concatenate bialleleic and multiallelic variants###########
cat gnomad_chr${i}_MultiallelicPosition_CADD_Mean_SD.txt gnomad_chr${i}_OnlyUniqueVarsBiallelic > chr${i}_ALL_GNOMAD_VARIANTS_CADD

####merge with phyloP scores#####
R --vanilla <<EOF

PhyloP_file <- read.csv("gnomad_chr${i}_IntersectedWithPhyloP.bed", sep="\t", head=F)
Cadd_file <- read.csv("chr${i}_ALL_GNOMAD_VARIANTS_CADD", sep="\t", head=T)

TOT_merged <- merge(PhyloP_file, Cadd_file, by=c("V1","V2"))
colnames(TOT_merged) <- c("chr", "pos","PhyloP","CADD_Phred", "CADD_SD")   ###If CADD-raw is extracted, then the 3rd column header is "CADD_Raw"

write.table(TOT_merged, "Merged_PhyloP_CADD_chr${i}.txt", sep="\t", quote=F, row.names=F)

EOF

done<Chromosome_list   ####Chromosome_list is the list of chromosomes


#####Concatenate all chromosomes#######
R --vanilla <<EOF

file_list <- list.files("$Path", pattern=glob2rx("Merged_PhyloP_CADD_chr*.txt"), full.names = T)
file_list
CompleteFile_AllChrom <- do.call(rbind, lapply(file_list, read.csv, header = T, sep="\t"))

write.table(CompleteFile_AllChrom, "MERGED_PhYLOP_CADD_ALL_CHROMOSOMES.txt", sep="\t", quote=F, row.names=F)

EOF


########Spearman correlation CADD-raw vs PhyloP######
#####Use the files where CADD-raw has been extracted####
R --vanilla <<EOF

file_list <- list.files("$Path", pattern=glob2rx("Merged_PhyloP_CADD_chr*.txt"), full.names = T)
file_list
CompleteFile_AllChrom <- do.call(rbind, lapply(file_list, read.csv, header = T, sep="\t"))

####Run correlation test#####

CorSpearman <- cor.test(CompleteFile_AllChrom\$PhyloP, CompleteFile_AllChrom\$CADD_Raw, method="spearman")
sink("Correlation_PhyloP_CADDRaw_Spearman.txt")
print(CorSpearman)
sink()

EOF


