#!/bin/bash
#  
#SBATCH --job-name=cohort_common_snplist
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3 
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/common_snplist.out 
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/common_snplist.err 
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu


# Define SNP list files
# Define SNP list files for Asian ancestry
snp_lists_asian=("asian_filtered_all_chr1.snplist" "asian_filtered_all_chr2.snplist" "asian_filtered_all_chr3.snplist" 
                 "asian_filtered_all_chr4.snplist" "asian_filtered_all_chr5.snplist" "asian_filtered_all_chr6.snplist" 
                 "asian_filtered_all_chr7.snplist" "asian_filtered_all_chr8.snplist" "asian_filtered_all_chr9.snplist" 
                 "asian_filtered_all_chr10.snplist" "asian_filtered_all_chr11.snplist" "asian_filtered_all_chr12.snplist" 
                 "asian_filtered_all_chr13.snplist" "asian_filtered_all_chr14.snplist" "asian_filtered_all_chr15.snplist" 
                 "asian_filtered_all_chr16.snplist" "asian_filtered_all_chr17.snplist" "asian_filtered_all_chr18.snplist" 
                 "asian_filtered_all_chr19.snplist" "asian_filtered_all_chr20.snplist" "asian_filtered_all_chr21.snplist" 
                 "asian_filtered_all_chr22.snplist")

# Define SNP list files for Mixed ancestry
snp_lists_mixed=("mixed_filtered_all_chr1.snplist" "mixed_filtered_all_chr2.snplist" "mixed_filtered_all_chr3.snplist" 
                 "mixed_filtered_all_chr4.snplist" "mixed_filtered_all_chr5.snplist" "mixed_filtered_all_chr6.snplist" 
                 "mixed_filtered_all_chr7.snplist" "mixed_filtered_all_chr8.snplist" "mixed_filtered_all_chr9.snplist" 
                 "mixed_filtered_all_chr10.snplist" "mixed_filtered_all_chr11.snplist" "mixed_filtered_all_chr12.snplist" 
                 "mixed_filtered_all_chr13.snplist" "mixed_filtered_all_chr14.snplist" "mixed_filtered_all_chr15.snplist" 
                 "mixed_filtered_all_chr16.snplist" "mixed_filtered_all_chr17.snplist" "mixed_filtered_all_chr18.snplist" 
                 "mixed_filtered_all_chr19.snplist" "mixed_filtered_all_chr20.snplist" "mixed_filtered_all_chr21.snplist" 
                 "mixed_filtered_all_chr22.snplist")

# Define SNP list files for Black ancestry
snp_lists_black=("black_filtered_all_chr1.snplist" "black_filtered_all_chr2.snplist" "black_filtered_all_chr3.snplist" 
                 "black_filtered_all_chr4.snplist" "black_filtered_all_chr5.snplist" "black_filtered_all_chr6.snplist" 
                 "black_filtered_all_chr7.snplist" "black_filtered_all_chr8.snplist" "black_filtered_all_chr9.snplist" 
                 "black_filtered_all_chr10.snplist" "black_filtered_all_chr11.snplist" "black_filtered_all_chr12.snplist" 
                 "black_filtered_all_chr13.snplist" "black_filtered_all_chr14.snplist" "black_filtered_all_chr15.snplist" 
                 "black_filtered_all_chr16.snplist" "black_filtered_all_chr17.snplist" "black_filtered_all_chr18.snplist" 
                 "black_filtered_all_chr19.snplist" "black_filtered_all_chr20.snplist" "black_filtered_all_chr21.snplist" 
                 "black_filtered_all_chr22.snplist")

# Define SNP list files for Chinese ancestry
snp_lists_chinese=("chinese_filtered_all_chr1.snplist" "chinese_filtered_all_chr2.snplist" "chinese_filtered_all_chr3.snplist" 
                   "chinese_filtered_all_chr4.snplist" "chinese_filtered_all_chr5.snplist" "chinese_filtered_all_chr6.snplist" 
                   "chinese_filtered_all_chr7.snplist" "chinese_filtered_all_chr8.snplist" "chinese_filtered_all_chr9.snplist" 
                   "chinese_filtered_all_chr10.snplist" "chinese_filtered_all_chr11.snplist" "chinese_filtered_all_chr12.snplist" 
                   "chinese_filtered_all_chr13.snplist" "chinese_filtered_all_chr14.snplist" "chinese_filtered_all_chr15.snplist" 
                   "chinese_filtered_all_chr16.snplist" "chinese_filtered_all_chr17.snplist" "chinese_filtered_all_chr18.snplist" 
                   "chinese_filtered_all_chr19.snplist" "chinese_filtered_all_chr20.snplist" "chinese_filtered_all_chr21.snplist" 
                   "chinese_filtered_all_chr22.snplist")

# Initialize the common SNPs file with the first SNP list
cp "${snp_lists_asian[0]}" common_snps_asian.txt
cp "${snp_lists_mixed[0]}" common_snps_mixed.txt
cp "${snp_lists_black[0]}" common_snps_black.txt
cp "${snp_lists_chinese[0]}" common_snps_chinese.txt

# Loop through the remaining SNP lists and find common SNPs
for ((i=1; i<${#snp_lists_asian[@]}; i++)); do
  comm -12 <(sort common_snps_asian.txt) <(sort "${snp_lists_asian[i]}") > temp_common_snps_asian.txt
  mv temp_common_snps_asian.txt common_snps_asian.txt
  done
  
for ((i=1; i<${#snp_lists_mixed[@]}; i++)); do
  comm -12 <(sort common_snps_mixed.txt) <(sort "${snp_lists_mixed[i]}") > temp_common_snps_mixed.txt
  mv temp_common_snps_mixed.txt common_snps_mixed.txt
  done

for ((i=1; i<${#snp_lists_black[@]}; i++)); do
  comm -12 <(sort common_snps_black.txt) <(sort "${snp_lists_black[i]}") > temp_common_snps_black.txt
  mv temp_common_snps_black.txt common_snps_black.txt
  done
  
for ((i=1; i<${#snp_lists_chinese[@]}; i++)); do
  comm -12 <(sort common_snps_chinese.txt) <(sort "${snp_lists_chinese[i]}") > temp_common_snps_chinese.txt
  mv temp_common_snps_chinese.txt common_snps_chinese.txt
  done


  
  