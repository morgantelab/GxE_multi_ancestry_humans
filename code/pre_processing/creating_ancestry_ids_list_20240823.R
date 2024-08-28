## creating ancestry id files ##

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data")

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240821.RData")

# Assuming the dataset dtt has a column 'ID' that uniquely identifies individuals

# Create ID lists for each ethnicity and repeat the IDs in two columns
asian_ids <- data.frame(ID1 = dtt$ID[dtt$ethn1_asian == 1], ID2 = dtt$ID[dtt$ethn1_asian == 1])
white_ids <- data.frame(ID1 = dtt$ID[dtt$ethn1_white == 1], ID2 = dtt$ID[dtt$ethn1_white == 1])
mixed_ids <- data.frame(ID1 = dtt$ID[dtt$ethn1_mixed == 1], ID2 = dtt$ID[dtt$ethn1_mixed == 1])
black_ids <- data.frame(ID1 = dtt$ID[dtt$ethn1_black == 1], ID2 = dtt$ID[dtt$ethn1_black == 1])
chinese_ids <- data.frame(ID1 = dtt$ID[dtt$ethn1_chinese == 1], ID2 = dtt$ID[dtt$ethn1_chinese == 1])

# Save the ID lists to files with two columns of the same IDs
write.table(asian_ids, "asian_ids.txt", row.names = FALSE, col.names = FALSE)
write.table(white_ids, "white_ids.txt", row.names = FALSE, col.names = FALSE)
write.table(mixed_ids, "mixed_ids.txt", row.names = FALSE, col.names = FALSE)
write.table(black_ids, "black_ids.txt", row.names = FALSE, col.names = FALSE)
write.table(chinese_ids, "chinese_ids.txt", row.names = FALSE, col.names = FALSE)
