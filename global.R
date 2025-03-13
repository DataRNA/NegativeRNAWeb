# Gerekli kütüphaneleri yükleme
if (!require("shiny")) install.packages("shiny")
if (!require("Biostrings")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
if (!require("seqinr")) install.packages("seqinr")
if (!require("shinythemes")) install.packages("shinythemes")
if (!require("writexl")) install.packages("writexl")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("magick")) install.packages("magick")
if (!require("base64enc")) install.packages("base64enc")

# Kütüphaneleri yükleme
library(shiny)
library(Biostrings)
library(seqinr)
library(shinythemes)
library(writexl)
library(ggplot2)
library(magick)
library(base64enc)

# Karakter kodlaması ayarı
Sys.setlocale("LC_ALL", "Turkish")

