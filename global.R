# Gerekli kütüphaneleri yükleme
if (!require("shiny")) install.packages("shiny")
if (!require("Biostrings")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
if (!require("seqinr")) install.packages("seqinr")
if (!require("shinythemes")) install.packages("shinythemes")

# Kütüphaneleri yükleme
library(shiny)
library(Biostrings)
library(seqinr)
library(shinythemes) 