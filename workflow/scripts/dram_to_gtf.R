library(argparse)
library(tidyverse)

parser <- ArgumentParser()

parser$add_argument(
  "-i", "--input",
  dest = "input", help = "input DRAM annotation"
)

parser$add_argument(
  "-o", "--output",
  dest = "output", help = "output GTF file"
)

dram_to_gtf <- function(dram_path) {
  dram_path %>%
    read_tsv() %>%
    select(
      gene_id = `...1`, scaffold, start_position, end_position, strandedness
    ) %>%
    mutate(
      `#seqname` = scaffold,
      source = "DRAM",
      feature = "CDS",
      start = start_position,
      end = end_position,
      score = ".",
      strand = case_when(
        strandedness == "1" ~ "+",
        strandedness == "-1" ~ "-",
        .default = ".",
      ),
      frame = ".",
      attribute = str_glue("gene_id {gene_id}")
    ) %>%
    select(
      `#seqname`, source, feature, start, end, score, strand, frame, attribute
    )
}

args <- parser$parse_args()

args$input %>%
  dram_to_gtf() %>%
  write_tsv(file = args$output)
