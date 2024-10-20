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

set.seed(333)

parser <-
  ArgumentParser(description = 'Perform DESeq2 and GSEA on bulk RNA-seq')
parser$add_argument(
  '-input',
  '--i',
  type = "character",
  required = TRUE,
  nargs = '*',
  help = 'All kallisto outputs'
)

parser$add_argument(
  '-species',
  '--s',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Species name (Mus musculus, Homo sapiens); CASE-SENSITIVE'
)
args <- parser$parse_args()

indir <- args$i

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

PerformDGETests(indir, species)
