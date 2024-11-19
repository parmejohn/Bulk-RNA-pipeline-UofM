#!/usr/local/bin/Rscript

library(argparse)
library(biomaRt)
library(data.table)
library(dplyr)
library(tximportData)
library(DESeq2)
library(readr)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(tximport)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(RColorBrewer)
library(ComplexHeatmap)
library(UpSetR)

set.seed(333)

parser <-
  ArgumentParser(description = 'Perform DESeq2 and GSEA on bulk RNA-seq')
parser$add_argument(
  '-input',
  '--i',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Count matrix in a .tsv format, first two columns should be gene_id and gene_name'
)

parser$add_argument(
  '-normalized',
  '--n',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Normalized count matrix in a .tsv format, first two columns should be gene_id and gene_name'
)

parser$add_argument(
  '-species',
  '--s',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Species name (Mus musculus, Homo sapiens); CASE-SENSITIVE'
)

parser$add_argument(
  '-table',
  '--t',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'csv table for sample and condition'
)

parser$add_argument(
  '-filter',
  '--f',
  type = "character",
  required = TRUE,
  nargs = '*',
)

args <- parser$parse_args()

cnts <- read_tsv(args$i)
normalized.cnts <- read_tsv(args$n)

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/dge_analyses.R"))

if (args$s == "musmusculus") {
  species <- "Mus musculus"
} else if (args$s == "homosapiens") {
  species <- "Homo sapiens"
} else {
  print("bad")
}

sample.table <- read.csv(args$t)


PerformDGETests(cnts, normalized.cnts, species, sample.table, args$filter)
