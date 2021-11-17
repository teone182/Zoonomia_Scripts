#!/bin/bash -l

InFolder='Path/To/Start/file'
WorkFolder='Path'
AnnoFolder='Path/To/All/Gencode/Annotations/and/Cadd/values'  ###this dir contains Gencode and other regulatory annotation bed files#####

##########
cd $WorkFolder

######All PhyloP >= 2.27##########
awk -v FS="\t" -v OFS="\t" 'FNR == 1 {next} $3 >= 2.27 {print $1,$2-1,$2,$3,$4}' $InFolder/MERGED_PHYLOP_CADD_ALL_CHROMOSOMES.txt > $WorkFolder/Variants_WithPhyloPFDR05.bed &&

cd $WorkFolder
while read n
do

cd $WorkFolder

######All CADD >= 10-20##########
awk -v var="${n}" -v FS="\t" -v OFS="\t" 'FNR == 1 {next} $4 >= var {print $1,$2-1,$2,$3,$4}' $InFolder/MERGED_PHYLOP_CADD_ALL_CHROMOSOMES.txt > $WorkFolder/Variants_WithCADDGr${n}.bed

done<$AnnoFolder/List_CADD_values.list   #####list of values from 10 to 20#####

while read n
do
while read class
do

cd $WorkFolder

#####################
###PhyloP

Nr_vars_class=$(cat $AnnoFolder/${class}.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}') &&

bedtools intersect -a $AnnoFolder/${class}.bed -b Variants_WithPhyloPFDR05.bed -wo | awk -v OFS="\t" '{print $5,$7,$8,$9,$4}' > ${class}_PhyloP_FDR05 &&

Nr_ThisClass_PhyloPyes=$(wc -l "${class}_PhyloP_FDR05" | cut -f1 -d' ') &&

if [[ $Nr_ThisClass_PhyloPyes -eq 0 ]]
then
Fraction_PhyloPyes_ThisClass=0
else
Fraction_PhyloPyes_ThisClass=$(echo "scale=3; $Nr_ThisClass_PhyloPyes/$Nr_vars_class" | bc)
fi

echo -e "${class}\t$Fraction_PhyloPyes_ThisClass\tPhyloP>=2.27" > Fraction_PhyloPYES_Nr${n}_${class}.txt &&


#####################
###CADD
Nr_vars_class=$(cat $AnnoFolder/${class}.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}') &&

bedtools intersect -a $AnnoFolder/${class}.bed -b Variants_WithCADDGr${n}.bed -wo | awk -v OFS="\t" '{print $5,$7,$8,$9,$4}' > ${class}_CADD_Gr${n} &&

Nr_ThisClass_CADDyes=$(wc -l "${class}_CADD_Gr${n}" | cut -f1 -d' ') &&

if [[ $Nr_ThisClass_CADDyes -eq 0 ]]
then
Fraction_CADDyes_ThisClass=0
else
Fraction_CADDyes_ThisClass=$(echo "scale=3; $Nr_ThisClass_CADDyes/$Nr_vars_class" | bc)
fi

echo -e "${class}\t$Fraction_CADDyes_ThisClass\tCADD>=${n}" > Fraction_CADDyes_CADD${n}_${class}.txt

done<$AnnoFolder/List_ALL_Genomic_Categories.list  #####list of annotation categories (e.g pickOne_cds, pickAll_intron, etc.#####
done<$AnnoFolder/List_CADD_values.list      #####list of values from 10 to 20#####
