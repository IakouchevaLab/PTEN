library(tidyverse)
library(Gviz)

# Loading Variant Mapping Results
Final_result_exon <- readxl::read_xlsx("./Results/PTEN_Variant_Mapping_res.xlsx")

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
annotation <- input_exon %>%
  mutate(
    chromosome = "chr10"
  ) %>%
  dplyr::rename(
    gene = Gene_name,
    start = Exon_start,
    end = Exon_end,
    exon = Exon_ID,
    transcript = Transcript_ID
  ) %>%
  dplyr::select(chromosome, gene, transcript, exon, start, end)

# Plot Preparation
dtrack_count <- sapply(setNames(nm = unique(Final_result_exon$Transcript_ID)), function(t){
  
  temp_input <- Final_result_exon %>%
    filter(Transcript_ID == t) %>%
    pivot_wider(names_from = "Source", values_from = "Variant_Counts") %>%
    dplyr::select(-Total) %>%
    left_join(., annotation %>%
                dplyr::select(exon, start, end),
              by = c("Exon_ID" = "exon")) %>%
    mutate(end = ifelse(end - start < 1000, 
                        yes = end+1000, no = end))
  
  DataTrack(
    temp_input,
    name = "Variants\nCount",
    start = temp_input$start,
    end = temp_input$end,
    genome = "hg38",
    chromosome = "chr10",
    backgroud.title = "darkblue",
    cex.title = 6,
    cex.axis = 4,
    cex.legend = 6,
    #fontsize.legend = 100,
    col = c("magenta", "darkgreen"),
    col.histogram = "white"
  )
    
})
  
grtrack <- sapply(setNames(nm = unique(Final_result_exon$Transcript_ID)), function(t){
  
  temp_gr <- GeneRegionTrack(
    annotation %>%
      filter(transcript == t),
    genome = "hg38",
    chromosome = "chr10",
    transcriptAnnotation = "transcript",
    fill = "darkgray",
    col.line = "black",
    lex = 20,
    lty = "solid",
    lwd = 20,
    background.panel = "white",
    fontcolor.group = "black",
    background.title = "white",
    col.title = "white",
    cex.group = 6,
    just.group = "above",
    showId=TRUE,
    geneSymbol=FALSE
  )
  
})
gtrack <- GenomeAxisTrack(cex = 4, fontcolor = "black")
itrack <- IdeogramTrack(genome = "hg19", chromosome = "chr10", 
                        cex = 4, cex.bands = 6, fontcolor = "black")

pdf("./Plots/Variant_Mapping_Plot.pdf", height = 96, width = 80)
plotTracks(list(itrack, gtrack,
                grtrack[[1]], dtrack_count[[1]],
                grtrack[[2]], dtrack_count[[2]],
                grtrack[[3]], dtrack_count[[3]],
                grtrack[[4]], dtrack_count[[4]],
                grtrack[[5]], dtrack_count[[5]],
                grtrack[[6]], dtrack_count[[6]],
                grtrack[[7]], dtrack_count[[7]]
),
from = 87850000, to = 88000000,
groups = rep(c("MasterTable", "Ruth_Lab"), each = 1),
type = c("histogram", legend = T)
)
dev.off()

pdf("./Plots/Variant_Mapping_Plot_123.pdf", height = 36, width = 80)
plotTracks(list(itrack, gtrack,
                grtrack[[1]], dtrack_count[[1]],
                grtrack[[2]], dtrack_count[[2]],
                grtrack[[3]], dtrack_count[[3]]
                ),
           from = 87850000, to = 88000000,
           groups = rep(c("MasterTable", "Ruth_Lab"), each = 1),
           type = c("histogram", legend = T)
           )
dev.off()

pdf("./Plots/Variant_Mapping_Plot_4567.pdf", height = 48, width = 80)
plotTracks(list(itrack, gtrack,
                grtrack[[4]], dtrack_count[[4]],
                grtrack[[5]], dtrack_count[[5]],
                grtrack[[6]], dtrack_count[[6]],
                grtrack[[7]], dtrack_count[[7]]
                ),
          from = 87850000, to = 88000000,
          groups = rep(c("MasterTable", "Ruth_Lab"), each = 1),
          type = c("histogram", legend = T)
)
dev.off()





























