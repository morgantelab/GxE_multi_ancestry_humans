#!/bin/bash
#  
#SBATCH --job-name=cohort_common_snplist
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3 
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/common_snplist.out 
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/common_snplist.err 
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

# Set working directory for SNP list files and output
input_dir="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/snp_lists"
output_dir="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output"

# Define SNP list files for Asian ancestry
snp_lists_asian=("$input_dir/asian_filtered_all_chr1.snplist" "$input_dir/asian_filtered_all_chr2.snplist" "$input_dir/asian_filtered_all_chr3.snplist" 
                 "$input_dir/asian_filtered_all_chr4.snplist" "$input_dir/asian_filtered_all_chr5.snplist" "$input_dir/asian_filtered_all_chr6.snplist" 
                 "$input_dir/asian_filtered_all_chr7.snplist" "$input_dir/asian_filtered_all_chr8.snplist" "$input_dir/asian_filtered_all_chr9.snplist" 
                 "$input_dir/asian_filtered_all_chr10.snplist" "$input_dir/asian_filtered_all_chr11.snplist" "$input_dir/asian_filtered_all_chr12.snplist" 
                 "$input_dir/asian_filtered_all_chr13.snplist" "$input_dir/asian_filtered_all_chr14.snplist" "$input_dir/asian_filtered_all_chr15.snplist" 
                 "$input_dir/asian_filtered_all_chr16.snplist" "$input_dir/asian_filtered_all_chr17.snplist" "$input_dir/asian_filtered_all_chr18.snplist" 
                 "$input_dir/asian_filtered_all_chr19.snplist" "$input_dir/asian_filtered_all_chr20.snplist" "$input_dir/asian_filtered_all_chr21.snplist" 
                 "$input_dir/asian_filtered_all_chr22.snplist")

# Define SNP list files for Mixed ancestry
snp_lists_mixed=("$input_dir/mixed_filtered_all_chr1.snplist" "$input_dir/mixed_filtered_all_chr2.snplist" "$input_dir/mixed_filtered_all_chr3.snplist" 
                 "$input_dir/mixed_filtered_all_chr4.snplist" "$input_dir/mixed_filtered_all_chr5.snplist" "$input_dir/mixed_filtered_all_chr6.snplist" 
                 "$input_dir/mixed_filtered_all_chr7.snplist" "$input_dir/mixed_filtered_all_chr8.snplist" "$input_dir/mixed_filtered_all_chr9.snplist" 
                 "$input_dir/mixed_filtered_all_chr10.snplist" "$input_dir/mixed_filtered_all_chr11.snplist" "$input_dir/mixed_filtered_all_chr12.snplist" 
                 "$input_dir/mixed_filtered_all_chr13.snplist" "$input_dir/mixed_filtered_all_chr14.snplist" "$input_dir/mixed_filtered_all_chr15.snplist" 
                 "$input_dir/mixed_filtered_all_chr16.snplist" "$input_dir/mixed_filtered_all_chr17.snplist" "$input_dir/mixed_filtered_all_chr18.snplist" 
                 "$input_dir/mixed_filtered_all_chr19.snplist" "$input_dir/mixed_filtered_all_chr20.snplist" "$input_dir/mixed_filtered_all_chr21.snplist" 
                 "$input_dir/mixed_filtered_all_chr22.snplist")

# Define SNP list files for Black ancestry
snp_lists_black=("$input_dir/black_filtered_all_chr1.snplist" "$input_dir/black_filtered_all_chr2.snplist" "$input_dir/black_filtered_all_chr3.snplist" 
                 "$input_dir/black_filtered_all_chr4.snplist" "$input_dir/black_filtered_all_chr5.snplist" "$input_dir/black_filtered_all_chr6.snplist" 
                 "$input_dir/black_filtered_all_chr7.snplist" "$input_dir/black_filtered_all_chr8.snplist" "$input_dir/black_filtered_all_chr9.snplist" 
                 "$input_dir/black_filtered_all_chr10.snplist" "$input_dir/black_filtered_all_chr11.snplist" "$input_dir/black_filtered_all_chr12.snplist" 
                 "$input_dir/black_filtered_all_chr13.snplist" "$input_dir/black_filtered_all_chr14.snplist" "$input_dir/black_filtered_all_chr15.snplist" 
                 "$input_dir/black_filtered_all_chr16.snplist" "$input_dir/black_filtered_all_chr17.snplist" "$input_dir/black_filtered_all_chr18.snplist" 
                 "$input_dir/black_filtered_all_chr19.snplist" "$input_dir/black_filtered_all_chr20.snplist" "$input_dir/black_filtered_all_chr21.snplist" 
                 "$input_dir/black_filtered_all_chr22.snplist")

# Define SNP list files for Chinese ancestry
snp_lists_chinese=("$input_dir/chinese_filtered_all_chr1.snplist" "$input_dir/chinese_filtered_all_chr2.snplist" "$input_dir/chinese_filtered_all_chr3.snplist" 
                   "$input_dir/chinese_filtered_all_chr4.snplist" "$input_dir/chinese_filtered_all_chr5.snplist" "$input_dir/chinese_filtered_all_chr6.snplist" 
                   "$input_dir/chinese_filtered_all_chr7.snplist" "$input_dir/chinese_filtered_all_chr8.snplist" "$input_dir/chinese_filtered_all_chr9.snplist" 
                   "$input_dir/chinese_filtered_all_chr10.snplist" "$input_dir/chinese_filtered_all_chr11.snplist" "$input_dir/chinese_filtered_all_chr12.snplist" 
                   "$input_dir/chinese_filtered_all_chr13.snplist" "$input_dir/chinese_filtered_all_chr14.snplist" "$input_dir/chinese_filtered_all_chr15.snplist" 
                   "$input_dir/chinese_filtered_all_chr16.snplist" "$input_dir/chinese_filtered_all_chr17.snplist" "$input_dir/chinese_filtered_all_chr18.snplist" 
                   "$input_dir/chinese_filtered_all_chr19.snplist" "$input_dir/chinese_filtered_all_chr20.snplist" "$input_dir/chinese_filtered_all_chr21.snplist" 
                   "$input_dir/chinese_filtered_all_chr22.snplist")

# Initialize the common SNPs file with the first SNP list
cp "${snp_lists_asian[0]}" "$output_dir/common_snps_asian.txt"
cp "${snp_lists_mixed[0]}" "$output_dir/common_snps_mixed.txt"
cp "${snp_lists_black[0]}" "$output_dir/common_snps_black.txt"
cp "${snp_lists_chinese[0]}" "$output_dir/common_snps_chinese.txt"

# Loop through the remaining SNP lists and find common SNPs for Asian ancestry
for ((i=1; i<${#snp_lists_asian[@]}; i++)); do
  comm -12 <(sort "$output_dir/common_snps_asian.txt") <(sort "${snp_lists_asian[i]}") > "$output_dir/temp_common_snps_asian.txt"
  mv "$output_dir/temp_common_snps_asian.txt" "$output_dir/common_snps_asian.txt"
done

# Loop through the remaining SNP lists and find common SNPs for Mixed ancestry
for ((i=1; i<${#snp_lists_mixed[@]}; i++)); do
  comm -12 <(sort "$output_dir/common_snps_mixed.txt") <(sort "${snp_lists_mixed[i]}") > "$output_dir/temp_common_snps_mixed.txt"
  mv "$output_dir/temp_common_snps_mixed.txt" "$output_dir/common_snps_mixed.txt"
done

# Loop through the remaining SNP lists and find common SNPs for Black ancestry
for ((i=1; i<${#snp_lists_black[@]}; i++)); do
  comm -12 <(sort "$output_dir/common_snps_black.txt") <(sort "${snp_lists_black[i]}") > "$output_dir/temp_common_snps_black.txt"
  mv "$output_dir/temp_common_snps_black.txt" "$output_dir/common_snps_black.txt"
done

# Loop through the remaining SNP lists and find common SNPs for Chinese ancestry
for ((i=1; i<${#snp_lists_chinese[@]}; i++)); do
  comm -12 <(sort "$output_dir/common_snps_chinese.txt") <(sort "${snp_lists_chinese[i]}") > "$output_dir/temp_common_snps_chinese.txt"
  mv "$output_dir/temp_common_snps_chinese.txt" "$output_dir/common_snps_chinese.txt"
done

# Combine the common SNPs across all ancestries into one file
comm -12 <(sort "$output_dir/common_snps_asian.txt") <(sort "$output_dir/common_snps_mixed.txt") | \
comm -12 - <(sort "$output_dir/common_snps_black.txt") | \
comm -12 - <(sort "$output_dir/common_snps_chinese.txt") > "$output_dir/common_snps_all.txt"

# Output the final combined common SNPs
echo "Final common SNPs across all ancestries written to $output_dir/common_snps_all.txt"
