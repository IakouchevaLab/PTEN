library(tidyverse)

# Loading PTEN isoform info 
input_exon <- read_tsv("./Data/PTEN_final_exon_input_ensembl.tsv", col_names = F)
colnames(input_exon) <- c("Seq_Type", "Exon_start", "Exon_end", "Gene_name",
                          "Transcript_ID", "Exon_ID")
input_exon <- input_exon %>%
  mutate(
    Gene_name = gsub(".*id ", "", gsub('\"', "", Gene_name)),
    Transcript_ID = gsub(".*id ", "", gsub('\"', "", Transcript_ID)),
    Exon_ID = gsub(".*id ", "", gsub('\"', "", Exon_ID))
  )

# Loading Variant Info
PTEN_variants <- readxl::read_xlsx("./Data/Variants_Info/PTEN_Variants_Union.xlsx") %>%
  dplyr::select("Gene_symbol", "Chromosome", "Position_start_hg38",
                "Position_end_hg38", "Source")

# Variants Mapping
Final_result_exon <- lapply(seq(nrow(input_exon)), function(r){
  
  temp_MasterTable_res <- PTEN_variants %>%
    filter(Source == "MasterTable") %>%
    mutate(is_in = ifelse(
      Position_start_hg38 >= input_exon$Exon_start[r] &
        Position_end_hg38 <= input_exon$Exon_end[r],
      yes = 1,
      no = 0
    ))
  MasterTable_sum <- sum(temp_MasterTable_res$is_in)
  
  temp_Ruth_res <- PTEN_variants %>%
    filter(Source == "Ruth_lab") %>%
    mutate(is_in = ifelse(
      Position_start_hg38 >= input_exon$Exon_start[r] &
        Position_end_hg38 <= input_exon$Exon_end[r],
      yes = 1,
      no = 0
    ))
  Ruth_sum <- sum(temp_Ruth_res$is_in)
  
  Total_variants <- MasterTable_sum + Ruth_sum
  
  data.frame(
    Gene_name = "PTEN",
    Transcript_ID = input_exon$Transcript_ID[r],
    Exon_ID = input_exon$Exon_ID[r],
    Source = c("MasterTable", "Ruth_Lab", "Total"),
    Variant_Counts = c(MasterTable_sum, Ruth_sum, Total_variants)
  )
  
}) %>%
  bind_rows()

writexl::write_xlsx(Final_result_exon, 
                    "./Results/PTEN_Variant_Mapping_res.xlsx",
                    col_names = T)

# Wider format results with exon start and end points
Final_result_exon %>%
  pivot_wider(id_cols = c("Gene_name", "Transcript_ID", "Exon_ID"),
              names_from = "Source", values_from = "Variant_Counts") %>%
  left_join(., input_exon %>% dplyr::select(Exon_ID, Exon_start, Exon_end),
            by = "Exon_ID") %>%
  dplyr::select(Gene_name, Transcript_ID, Exon_ID, Exon_start, Exon_end, 
                MasterTable, Ruth_Lab, Total) %>%
  writexl::write_xlsx(., "./Results/PTEN_Variant_Mapping_res_widerFormat.xlsx",
                      col_names = T)







