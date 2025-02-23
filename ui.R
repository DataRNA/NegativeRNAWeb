library(shiny)
library(shinythemes)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  titlePanel("FASTA Negatif Veri Oluşturucu"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_file", "FASTA Dosyası Yükle",
                accept = c(".fasta", ".fa", ".txt")),
      
      hr(),
      
      downloadButton("downloadNegative", "Negatif Verileri İndir")
    ),
    
    mainPanel(
      h4("Dosya Bilgileri"),
      textOutput("sequence_info"),
      
      hr(),
      
      h4("Kullanım Talimatları"),
      tags$ul(
        tags$li("FASTA formatında DNA sekans dosyanızı yükleyin."),
        tags$li("Sistem otomatik olarak negatif örnekler oluşturacaktır."),
        tags$li("'Negatif Verileri İndir' butonunu kullanarak sonuçları indirebilirsiniz.")
      )
    )
  )
) 