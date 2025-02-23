library(shiny)
library(Biostrings)
library(seqinr)

server <- function(input, output, session) {
  
  # FASTA dosyasını okuma fonksiyonu
  fasta_data <- reactive({
    req(input$fasta_file)
    shiny::runApp()
    inFile <- input$fasta_file
    if (is.null(inFile)) return(NULL)
    
    # FASTA dosyasını okuma
    sequences <- readDNAStringSet(inFile$datapath)
    return(sequences)
  })
  
  # Negatif veri oluşturma fonksiyonu
  create_negative_data <- reactive({
    req(fasta_data())
    
    sequences <- fasta_data()
    negative_seqs <- DNAStringSet()
    
    for(i in 1:length(sequences)) {
      # Orijinal sekansı al
      orig_seq <- sequences[[i]]
      
      # Rastgele karıştırma işlemi
      shuffled_seq <- sample(strsplit(as.character(orig_seq), '')[[1]])
      negative_seqs[[i]] <- DNAString(paste0(shuffled_seq, collapse=''))
    }
    
    names(negative_seqs) <- paste0("negative_", names(sequences))
    return(negative_seqs)
  })
  
  # Sonuçları gösterme
  output$sequence_info <- renderText({
    req(fasta_data())
    paste("Yüklenen FASTA dosyasında", length(fasta_data()), "adet sekans bulunmaktadır.")
  })
  
  # Negatif verileri indirme
  output$downloadNegative <- downloadHandler(
    filename = function() {
      paste0("negative_sequences_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      writeXStringSet(create_negative_data(), file)
    }
  )
}
