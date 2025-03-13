library(shiny)
library(Biostrings)
library(seqinr)
library(writexl)  # Excel dosyası oluşturmak için
library(ggplot2)

server <- function(input, output, session) {
  
  # FASTA dosyasını okuma fonksiyonu
  fasta_data <- reactive({
    req(input$fasta_file)
    
    inFile <- input$fasta_file
    if (is.null(inFile)) return(NULL)
    
    sequences <- readDNAStringSet(inFile$datapath)
    return(sequences)
  })
  
  # Maksimum kaydırma miktarını hesapla
  observe({
    req(fasta_data())
    
    # En kısa sekansın uzunluğunu bul
    min_length <- min(width(fasta_data()))
    max_shift <- max(1, min_length - 1)
    
    # Kaydırma slider'ını güncelle
    updateSliderInput(session, "nernaShift", 
                     max = max_shift,
                     value = min(input$nernaShift, max_shift))
  })
  
  # Negatif Data sayfasına git
  observeEvent(input$goToNegative, {
    updateTabsetPanel(session, "mainTabset", selected = "negative")
  })
  
  # RNAfold ile yapı tahmini yapan fonksiyon
  predict_structure <- reactive({
    req(fasta_data())
    
    # Geçici dosya oluştur
    temp_fasta <- tempfile(fileext = ".fasta")
    writeXStringSet(fasta_data(), temp_fasta)
    
    # RNAfold çalıştır
    results <- system2("D:/DUYGU/Desktop/RNAfold/RNAfold.exe",
                      args = c("--noPS", temp_fasta), 
                      stdout = TRUE,
                      stderr = TRUE)
    
    # Sonuçları parse et
    parsed_results <- data.frame(
      Sequence_Name = character(),
      RNA_Type = character(),
      Structure = character(),
      MFE = numeric(),
      stringsAsFactors = FALSE
    )
    
    i <- 1
    for(line in results) {
      if(startsWith(line, ">")) {
        parsed_results[i, "Sequence_Name"] <- substr(line, 2, nchar(line))
        parsed_results[i, "RNA_Type"] <- input$rnaType
      } else if(grepl("[().]+", line)) {
        structure <- gsub("\\s.*$", "", line)
        mfe <- as.numeric(gsub(".*[(]([-.0-9]+)\\s*[)].*", "\\1", line))
        parsed_results[i, "Structure"] <- structure
        parsed_results[i, "MFE"] <- mfe
        i <- i + 1
      }
    }
    
    # Geçici dosyayı sil
    unlink(temp_fasta)
    
    return(parsed_results)
  })
  
  # Metod 1: NeRNA ile negatif veri oluşturma (kaydırma)
  nernaData <- reactiveVal(NULL)
  
  observeEvent(input$generateNeRNA, {
    req(fasta_data())
    
    # Kaydırma miktarını al
    shift_amount <- input$nernaShift
    
    sequences <- fasta_data()
    negative_seqs <- DNAStringSet()
    
    for(i in 1:length(sequences)) {
      orig_seq <- as.character(sequences[[i]])
      
      # Kaydırma işlemi
      if(nchar(orig_seq) > shift_amount) {
        new_seq <- paste0(
          substr(orig_seq, nchar(orig_seq) - shift_amount + 1, nchar(orig_seq)),
          substr(orig_seq, 1, nchar(orig_seq) - shift_amount)
        )
      } else {
        # Sekans çok kısaysa normal komplementer al
        new_seq <- chartr("ACGT", "TGCA", orig_seq)
      }
      
      negative_seqs[[i]] <- DNAString(new_seq)
    }
    
    names(negative_seqs) <- paste0("neRNA_shift", shift_amount, "_", names(sequences))
    nernaData(negative_seqs)
    
    # Sonuç sekmesine geç
    updateTabsetPanel(session, "negativeResultTabs", selected = "NeRNA Results")
    
    showNotification("NeRNA sequences generated successfully!", type = "message")
  })
  
  # Metod 2: Random Shuffling
  shuffleData <- reactiveVal(NULL)
  
  observeEvent(input$generateShuffle, {
    req(fasta_data())
    
    # Random seed ayarla
    set.seed(input$shuffleSeed)
    
    sequences <- fasta_data()
    negative_seqs <- DNAStringSet()
    
    for(i in 1:length(sequences)) {
      orig_seq <- sequences[[i]]
      
      # Rastgele karıştırma işlemi
      shuffled_seq <- sample(strsplit(as.character(orig_seq), '')[[1]])
      negative_seqs[[i]] <- DNAString(paste0(shuffled_seq, collapse=''))
    }
    
    names(negative_seqs) <- paste0("shuffle_seed", input$shuffleSeed, "_", names(sequences))
    shuffleData(negative_seqs)
    
    # Sonuç sekmesine geç
    updateTabsetPanel(session, "negativeResultTabs", selected = "Random Shuffling Results")
    
    showNotification("Shuffled sequences generated successfully!", type = "message")
  })
  
  # Metod 3: Dinucleotide Shuffling
  dinucleotideData <- reactiveVal(NULL)
  
  observeEvent(input$generateDinucleotide, {
    req(fasta_data())
    
    # Random seed ayarla
    set.seed(input$dinucSeed)
    
    sequences <- fasta_data()
    negative_seqs <- DNAStringSet()
    
    for(i in 1:length(sequences)) {
      orig_seq <- as.character(sequences[[i]])
      
      # Dinucleotide shuffling (basitleştirilmiş)
      if(nchar(orig_seq) >= 4) {
        # Dinükleotidleri oluştur
        dinucs <- character(nchar(orig_seq) - 1)
        for(j in 1:(nchar(orig_seq)-1)) {
          dinucs[j] <- substr(orig_seq, j, j+1)
        }
        
        # Dinükleotidleri karıştır
        shuffled_dinucs <- sample(dinucs)
        
        # Yeni sekansı oluştur (basitleştirilmiş)
        new_seq <- shuffled_dinucs[1]
        for(j in 2:length(shuffled_dinucs)) {
          new_seq <- paste0(new_seq, substr(shuffled_dinucs[j], 2, 2))
        }
      } else {
        # Çok kısa sekanslar için normal shuffling
        new_seq <- paste0(sample(strsplit(orig_seq, '')[[1]]), collapse='')
      }
      
      negative_seqs[[i]] <- DNAString(new_seq)
    }
    
    names(negative_seqs) <- paste0("dinuc_seed", input$dinucSeed, "_", names(sequences))
    dinucleotideData(negative_seqs)
    
    # Sonuç sekmesine geç
    updateTabsetPanel(session, "negativeResultTabs", selected = "Dinucleotide Results")
    
    showNotification("Dinucleotide shuffled sequences generated successfully!", type = "message")
  })
  
  # Sonuçları gösterme
  output$sequence_info <- renderText({
    req(fasta_data())
    paste("Number of sequences in uploaded FASTA file:", length(fasta_data()),
          "| RNA Type:", input$rnaType)
  })
  
  # Yapı tahmin sonuçlarını gösterme
  output$structure_results <- renderTable({
    req(predict_structure())
    predict_structure()
  })
  
  # Özet tablo
  output$structure_summary <- renderTable({
    req(predict_structure())
    data.frame(
      "Total Sequences" = nrow(predict_structure()),
      "RNA Type" = input$rnaType,
      "Avg. MFE" = round(mean(predict_structure()$MFE), 2)
    )
  })
  
  # NeRNA sonuçlarını tablo olarak göster
  output$nernaTable <- renderTable({
    req(nernaData())
    
    # Sekansları tablo olarak göster
    data.frame(
      Sequence_Name = names(nernaData()),
      Original_Sequence = as.character(fasta_data()),
      Negative_Sequence = as.character(nernaData()),
      Shift_Amount = input$nernaShift,
      stringsAsFactors = FALSE
    )
  })
  
  # Shuffle sonuçlarını tablo olarak göster
  output$shuffleTable <- renderTable({
    req(shuffleData())
    
    # Sekansları tablo olarak göster
    data.frame(
      Sequence_Name = names(shuffleData()),
      Original_Sequence = as.character(fasta_data()),
      Negative_Sequence = as.character(shuffleData()),
      Random_Seed = input$shuffleSeed,
      stringsAsFactors = FALSE
    )
  })
  
  # Dinucleotide sonuçlarını tablo olarak göster
  output$dinucleotideTable <- renderTable({
    req(dinucleotideData())
    
    # Sekansları tablo olarak göster
    data.frame(
      Sequence_Name = names(dinucleotideData()),
      Original_Sequence = as.character(fasta_data()),
      Negative_Sequence = as.character(dinucleotideData()),
      Random_Seed = input$dinucSeed,
      stringsAsFactors = FALSE
    )
  })
  
  # Karşılaştırma tablosu
  output$comparisonTable <- renderTable({
    req(fasta_data())
    
    # Hangi yöntemler oluşturuldu?
    methods_generated <- c(
      "Original" = TRUE,
      "NeRNA" = !is.null(nernaData()),
      "Random Shuffling" = !is.null(shuffleData()),
      "Dinucleotide" = !is.null(dinucleotideData())
    )
    
    # Sadece oluşturulan yöntemleri göster
    methods <- names(methods_generated)[methods_generated]
    
    # Karşılaştırma tablosu oluştur
    comparison <- data.frame(
      Method = methods,
      Sequences = c(
        length(fasta_data()),
        if("NeRNA" %in% methods) length(nernaData()) else NULL,
        if("Random Shuffling" %in% methods) length(shuffleData()) else NULL,
        if("Dinucleotide" %in% methods) length(dinucleotideData()) else NULL
      ),
      Parameters = c(
        "N/A",
        if("NeRNA" %in% methods) paste("Shift:", input$nernaShift) else NULL,
        if("Random Shuffling" %in% methods) paste("Seed:", input$shuffleSeed) else NULL,
        if("Dinucleotide" %in% methods) paste("Seed:", input$dinucSeed) else NULL
      ),
      stringsAsFactors = FALSE
    )
    
    return(comparison)
  })
  
  # Grafikler
  output$nernaPlot <- renderPlot({
    req(nernaData())
    
    # Basit bir görselleştirme
    orig_lengths <- width(fasta_data())
    neg_lengths <- width(nernaData())
    
    df <- data.frame(
      Sequence = rep(names(fasta_data()), 2),
      Type = c(rep("Original", length(fasta_data())), rep("NeRNA", length(nernaData()))),
      Length = c(orig_lengths, neg_lengths)
    )
    
    ggplot(df, aes(x = Sequence, y = Length, fill = Type)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(title = "Sequence Length Comparison", x = "Sequence", y = "Length") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$shufflePlot <- renderPlot({
    req(shuffleData())
    
    # Basit bir görselleştirme
    orig_lengths <- width(fasta_data())
    neg_lengths <- width(shuffleData())
    
    df <- data.frame(
      Sequence = rep(names(fasta_data()), 2),
      Type = c(rep("Original", length(fasta_data())), rep("Shuffled", length(shuffleData()))),
      Length = c(orig_lengths, neg_lengths)
    )
    
    ggplot(df, aes(x = Sequence, y = Length, fill = Type)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(title = "Sequence Length Comparison", x = "Sequence", y = "Length") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$dinucleotidePlot <- renderPlot({
    req(dinucleotideData())
    
    # Basit bir görselleştirme
    orig_lengths <- width(fasta_data())
    neg_lengths <- width(dinucleotideData())
    
    df <- data.frame(
      Sequence = rep(names(fasta_data()), 2),
      Type = c(rep("Original", length(fasta_data())), rep("Dinucleotide", length(dinucleotideData()))),
      Length = c(orig_lengths, neg_lengths)
    )
    
    ggplot(df, aes(x = Sequence, y = Length, fill = Type)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(title = "Sequence Length Comparison", x = "Sequence", y = "Length") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$comparisonPlot <- renderPlot({
    # En az bir yöntem oluşturulmuş olmalı
    req(any(!is.null(nernaData()), !is.null(shuffleData()), !is.null(dinucleotideData())))
    
    # Veri çerçevesi oluştur
    df <- data.frame(
      Sequence = character(),
      Method = character(),
      Length = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Orijinal sekansları ekle
    df <- rbind(df, data.frame(
      Sequence = names(fasta_data()),
      Method = "Original",
      Length = width(fasta_data()),
      stringsAsFactors = FALSE
    ))
    
    # NeRNA sekanslarını ekle
    if(!is.null(nernaData())) {
      df <- rbind(df, data.frame(
        Sequence = gsub("^neRNA_shift[0-9]+_", "", names(nernaData())),
        Method = "NeRNA",
        Length = width(nernaData()),
        stringsAsFactors = FALSE
      ))
    }
    
    # Shuffle sekanslarını ekle
    if(!is.null(shuffleData())) {
      df <- rbind(df, data.frame(
        Sequence = gsub("^shuffle_seed[0-9]+_", "", names(shuffleData())),
        Method = "Random Shuffling",
        Length = width(shuffleData()),
        stringsAsFactors = FALSE
      ))
    }
    
    # Dinucleotide sekanslarını ekle
    if(!is.null(dinucleotideData())) {
      df <- rbind(df, data.frame(
        Sequence = gsub("^dinuc_seed[0-9]+_", "", names(dinucleotideData())),
        Method = "Dinucleotide",
        Length = width(dinucleotideData()),
        stringsAsFactors = FALSE
      ))
    }
    
    # Grafik oluştur
    ggplot(df, aes(x = Method, y = Length, fill = Method)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = "Sequence Length Comparison Across Methods", 
           x = "Method", y = "Sequence Length") +
      theme(legend.position = "none")
  })
  
  # Excel olarak indirme
  output$downloadExcel <- downloadHandler(
    filename = function() {
      paste0("rna_analysis_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      # Tüm sonuçları bir liste olarak topla
      all_results <- list(
        "Structure_Predictions" = predict_structure(),
        "Comparison" = as.data.frame(output$comparisonTable)
      )
      
      # NeRNA sonuçları
      if(!is.null(nernaData())) {
        all_results[["NeRNA_Results"]] <- data.frame(
          Sequence_Name = names(nernaData()),
          Original_Sequence = as.character(fasta_data()),
          Negative_Sequence = as.character(nernaData()),
          Shift_Amount = input$nernaShift,
          stringsAsFactors = FALSE
        )
      }
      
      # Shuffle sonuçları
      if(!is.null(shuffleData())) {
        all_results[["Shuffle_Results"]] <- data.frame(
          Sequence_Name = names(shuffleData()),
          Original_Sequence = as.character(fasta_data()),
          Negative_Sequence = as.character(shuffleData()),
          Random_Seed = input$shuffleSeed,
          stringsAsFactors = FALSE
        )
      }
      
      # Dinucleotide sonuçları
      if(!is.null(dinucleotideData())) {
        all_results[["Dinucleotide_Results"]] <- data.frame(
          Sequence_Name = names(dinucleotideData()),
          Original_Sequence = as.character(fasta_data()),
          Negative_Sequence = as.character(dinucleotideData()),
          Random_Seed = input$dinucSeed,
          stringsAsFactors = FALSE
        )
      }
      
      # Excel dosyasına yaz
      writexl::write_xlsx(all_results, path = file)
    }
  )
  
  # CSV olarak indirme
  output$downloadCSV <- downloadHandler(
    filename = function() {
      paste0("structure_predictions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(predict_structure())
      write.csv(predict_structure(), file, row.names = FALSE)
    }
  )
  
  # NeRNA indirme
  output$downloadNeRNA <- downloadHandler(
    filename = function() {
      paste0("neRNA_shift", input$nernaShift, "_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      req(nernaData())
      writeXStringSet(nernaData(), file)
    }
  )
  
  # Shuffle indirme
  output$downloadShuffle <- downloadHandler(
    filename = function() {
      paste0("shuffled_seed", input$shuffleSeed, "_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      req(shuffleData())
      writeXStringSet(shuffleData(), file)
    }
  )
  
  # Dinucleotide indirme
  output$downloadDinucleotide <- downloadHandler(
    filename = function() {
      paste0("dinucleotide_seed", input$dinucSeed, "_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      req(dinucleotideData())
      writeXStringSet(dinucleotideData(), file)
    }
  )
}
