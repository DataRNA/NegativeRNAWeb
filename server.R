library(shiny)
library(Biostrings)
library(seqinr)
library(writexl)  # Excel dosyası oluşturmak için
library(ggplot2)
library(magick)

# Başlangıçta .temp dizinini oluştur ve yolunu sakla
temp_dir <- file.path(getwd(), ".temp")
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}

server <- function(input, output, session) {
  
  # FASTA dosyasını okuma fonksiyonu
  fasta_data <- reactive({
    req(input$fasta_file)
    
    inFile <- input$fasta_file
    if (is.null(inFile)) return(NULL)
    
    # FASTA dosyasını oku
    sequences <- readRNAStringSet(inFile$datapath)
    
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
  
  # RNAfold ile yapı tahmini ve görselleştirme
  predict_structure_with_images <- function(sequences, output_dir) {
    # Geçici dosya oluştur
    temp_fasta <- file.path(temp_dir, paste0("temp_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta"))
    writeXStringSet(sequences, temp_fasta)
    
    # Çıktı dizini oluştur
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # RNAfold'un tam yolunu oluştur
    rnafold_path <- file.path(getwd(), "app/RNAfold.exe")
    message("RNAfold path: ", rnafold_path)
    
    # RNAfold'un varlığını kontrol et
    if (!file.exists(rnafold_path)) {
      stop("RNAfold.exe not found at: ", rnafold_path)
    }
    
    # Mevcut çalışma dizinini kaydet
    original_wd <- getwd()
    
    # Çalışma dizinini temp_dir'e değiştir
    setwd(temp_dir)
    
    # Her sekans için ayrı çıktı dosyası oluştur
    parsed_results <- data.frame(
      Sequence_Name = names(sequences),
      RNA_Type = rep(input$rnaType, length(sequences)),
      Structure = character(length(sequences)),
      MFE = numeric(length(sequences)),
      Image_Path = character(length(sequences)),
      stringsAsFactors = FALSE
    )
    
    for(i in 1:length(sequences)) {
      # Her sekans için ayrı bir FASTA dosyası oluştur
      seq_fasta <- paste0("seq_", i, ".fasta")
      writeXStringSet(sequences[i], seq_fasta)
      
      # Çıktı dosya yolunu belirle
      output_file <- file.path(output_dir, paste0(parsed_results$Sequence_Name[i], "_ss.ps"))
      parsed_results$Image_Path[i] <- output_file
      
      # RNAfold çalıştır
      results <- system2(rnafold_path,
                        args = c(seq_fasta), 
                        stdout = TRUE,
                        stderr = TRUE)
      
      # PS dosyasını doğru konuma taşı
      ps_file <- paste0(parsed_results$Sequence_Name[i], "_ss.ps")
      if (file.exists(ps_file)) {
        file.rename(ps_file, output_file)
      }
      
      # Sonuçları parse et
      for(line in results) {
        if(grepl("[().]+", line)) {
          structure <- gsub("\\s.*$", "", line)
          mfe <- as.numeric(gsub(".*[(]([-.0-9]+)\\s*[)].*", "\\1", line))
          
          parsed_results$Structure[i] <- structure
          parsed_results$MFE[i] <- mfe
          break
        }
      }
      
      # Geçici FASTA dosyasını sil
      unlink(seq_fasta)
    }
    
    # Çalışma dizinini geri al
    setwd(original_wd)
    
    # Ana geçici dosyayı sil
    unlink(temp_fasta)
    
    return(parsed_results)
  }
  
  # Orijinal sekanslar için yapı tahmini
  predict_structure <- reactive({
    req(fasta_data())
    
    # Görsellerle birlikte yapı tahmini yap
    output_dir <- file.path(temp_dir, "original_structures")
    predict_structure_with_images(fasta_data(), output_dir)
  })
  
  # Metod 1: NeRNA ile negatif veri oluşturma (kaydırma)
  nernaData <- reactiveVal(NULL)
  
  observeEvent(input$generateNeRNA, {
    req(fasta_data())
    
    # Kaydırma miktarını al
    shift_amount <- input$nernaShift
    
    sequences <- fasta_data()
    negative_seqs <- RNAStringSet()
    
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
      
      negative_seqs[[i]] <- RNAString(new_seq)
    }
    
    names(negative_seqs) <- paste0("neRNA_shift", shift_amount, "_", names(sequences))
    nernaData(negative_seqs)
    
    # Sonuç sekmesine geç
    updateTabsetPanel(session, "negativeResultTabs", selected = "NeRNA Results")
    
    showNotification("NeRNA sequences generated successfully!", type = "message")
  })
  
  # NeRNA sekansları için yapı tahmini
  nerna_structures <- reactive({
    req(nernaData())
    
    # Görsellerle birlikte yapı tahmini yap
    output_dir <- file.path(temp_dir, "nerna_structures")
    predict_structure_with_images(nernaData(), output_dir)
  })
  
  # Metod 2: Random Shuffling
  shuffleData <- reactiveVal(NULL)
  
  observeEvent(input$generateShuffle, {
    req(fasta_data())
    
    # Random seed ayarla
    set.seed(input$shuffleSeed)
    
    sequences <- fasta_data()
    negative_seqs <- RNAStringSet()
    
    for(i in 1:length(sequences)) {
      orig_seq <- sequences[[i]]
      
      # Rastgele karıştırma işlemi
      shuffled_seq <- sample(strsplit(as.character(orig_seq), '')[[1]])
      negative_seqs[[i]] <- RNAString(paste0(shuffled_seq, collapse=''))
    }
    
    names(negative_seqs) <- paste0("shuffle_seed", input$shuffleSeed, "_", names(sequences))
    shuffleData(negative_seqs)
    
    # Sonuç sekmesine geç
    updateTabsetPanel(session, "negativeResultTabs", selected = "Random Shuffling Results")
    
    showNotification("Shuffled sequences generated successfully!", type = "message")
  })
  
  # Shuffle sekansları için yapı tahmini
  shuffle_structures <- reactive({
    req(shuffleData())
    
    # Görsellerle birlikte yapı tahmini yap
    output_dir <- file.path(temp_dir, "shuffle_structures")
    predict_structure_with_images(shuffleData(), output_dir)
  })
  
  # Metod 3: Dinucleotide Shuffling
  dinucleotideData <- reactiveVal(NULL)
  
  observeEvent(input$generateDinucleotide, {
    req(fasta_data())
    
    # Random seed ayarla
    set.seed(input$dinucSeed)
    
    sequences <- fasta_data()
    negative_seqs <- RNAStringSet()
    
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
      
      negative_seqs[[i]] <- RNAString(new_seq)
    }
    
    names(negative_seqs) <- paste0("dinuc_seed", input$dinucSeed, "_", names(sequences))
    dinucleotideData(negative_seqs)
    
    # Sonuç sekmesine geç
    updateTabsetPanel(session, "negativeResultTabs", selected = "Dinucleotide Results")
    
    showNotification("Dinucleotide shuffled sequences generated successfully!", type = "message")
  })
  
  # Dinucleotide sekansları için yapı tahmini
  dinucleotide_structures <- reactive({
    req(dinucleotideData())
    
    # Görsellerle birlikte yapı tahmini yap
    output_dir <- file.path(temp_dir, "dinucleotide_structures")
    predict_structure_with_images(dinucleotideData(), output_dir)
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
  
  # Karşılaştırma verilerini saklayacak reaktif değer
  comparisonData <- reactiveVal(NULL)
  
  # Karşılaştırma tablosu
  output$comparisonTable <- renderTable({
    # En az bir yöntem oluşturulmuş olmalı
    req(any(!is.null(nernaData()), !is.null(shuffleData()), !is.null(dinucleotideData())))
    
    # Orijinal sekansların MFE değerlerini al
    orig_mfe <- mean(predict_structure()$MFE)
    
    # Karşılaştırma tablosu oluştur
    comparison <- data.frame(
      Method = "Original",
      Avg_Length = mean(width(fasta_data())),
      Avg_MFE = orig_mfe,
      stringsAsFactors = FALSE
    )
    
    # NeRNA sonuçlarını ekle
    if(!is.null(nernaData())) {
      # NeRNA sekansları için MFE hesapla
      nerna_mfe <- calculate_mfe(nernaData())
      
      comparison <- rbind(comparison, data.frame(
        Method = "NeRNA",
        Avg_Length = mean(width(nernaData())),
        Avg_MFE = mean(nerna_mfe$MFE),
        stringsAsFactors = FALSE
      ))
    }
    
    # Shuffle sonuçlarını ekle
    if(!is.null(shuffleData())) {
      # Shuffle sekansları için MFE hesapla
      shuffle_mfe <- calculate_mfe(shuffleData())
      
      comparison <- rbind(comparison, data.frame(
        Method = "Random Shuffling",
        Avg_Length = mean(width(shuffleData())),
        Avg_MFE = mean(shuffle_mfe$MFE),
        stringsAsFactors = FALSE
      ))
    }
    
    # Dinucleotide sonuçlarını ekle
    if(!is.null(dinucleotideData())) {
      # Dinucleotide sekansları için MFE hesapla
      dinuc_mfe <- calculate_mfe(dinucleotideData())
      
      comparison <- rbind(comparison, data.frame(
        Method = "Dinucleotide",
        Avg_Length = mean(width(dinucleotideData())),
        Avg_MFE = mean(dinuc_mfe$MFE),
        stringsAsFactors = FALSE
      ))
    }
    
    # Ortalama değerleri yuvarla
    comparison$Avg_Length <- round(comparison$Avg_Length, 2)
    comparison$Avg_MFE <- round(comparison$Avg_MFE, 2)
    
    # Karşılaştırma verilerini sakla
    comparisonData(comparison)
    
    return(comparison)
  })
  
  # Karşılaştırma grafiği
  output$comparisonPlot <- renderPlot({
    # En az bir yöntem oluşturulmuş olmalı
    req(any(!is.null(nernaData()), !is.null(shuffleData()), !is.null(dinucleotideData())))
    req(comparisonData())
    
    # Karşılaştırma verilerini al
    comparison_data <- comparisonData()
    
    # Veriyi uzun formata dönüştür
    library(reshape2)
    melted_data <- melt(comparison_data, id.vars = "Method", 
                        measure.vars = c("Avg_Length", "Avg_MFE"),
                        variable.name = "Metric", value.name = "Value")
    
    # Metrik isimlerini düzelt
    melted_data$Metric <- factor(melted_data$Metric, 
                                levels = c("Avg_Length", "Avg_MFE"),
                                labels = c("Average Length (nt)", "Average MFE (kcal/mol)"))
    
    # Grafik oluştur
    ggplot(melted_data, aes(x = Method, y = Value, fill = Method)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ Metric, scales = "free_y") +
      theme_minimal() +
      labs(title = "Comparison of Methods", x = "", y = "") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  })
  
  # MFE hesaplama fonksiyonu ekleyelim
  calculate_mfe <- function(sequences) {
    # Geçici dosya oluştur
    temp_fasta <- file.path(temp_dir, paste0("temp_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta"))
    writeXStringSet(sequences, temp_fasta)
    original_wd <- getwd()
    rnafold_path <- file.path(original_wd, "app/RNAfold.exe")

    # RNAfold çalıştır
    results <- system2(rnafold_path,
                       args = c(temp_fasta), 
                       stdout = TRUE,
                       stderr = TRUE)
    
    # Sonuçları parse et
    parsed_results <- data.frame(
      Sequence_Name = character(),
      Structure = character(),
      MFE = numeric(),
      stringsAsFactors = FALSE
    )
    
    i <- 1
    for(line in results) {
      if(startsWith(line, ">")) {
        parsed_results[i, "Sequence_Name"] <- substr(line, 2, nchar(line))
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
  }
  
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
  
  # PS dosyalarını PNG'ye dönüştürme fonksiyonu
  convert_ps_to_png <- function(ps_file) {
    if (!file.exists(ps_file)) {
      return(NULL)
    }
    
    # PNG dosya yolunu oluştur
    png_file <- sub("\\.ps$", ".png", ps_file)
    
    # ImageMagick ile dönüştür
    tryCatch({
      # Windows için
      convert_cmd <- "magick"
      
      system2(convert_cmd,
              args = c("convert", ps_file, png_file),
              stdout = TRUE,
              stderr = TRUE)
      
      if (file.exists(png_file)) {
        return(png_file)
      } else {
        return(NULL)
      }
    }, error = function(e) {
      message("Error converting PS to PNG: ", e$message)
      return(NULL)
    })
  }
  
  # Yapı görselleştirme fonksiyonu (sadece metin)
  output$originalStructureNeRNA <- renderUI({
    req(predict_structure())
    req(nernaData())
    
    # İlk sekansın yapısını göster
    structure_data <- predict_structure()[1, ]
    
    div(
      p(strong("Sequence:"), structure_data$Sequence_Name),
      p(strong("Structure:"), structure_data$Structure),
      p(strong("MFE:"), paste0(structure_data$MFE, " kcal/mol")),
      tags$pre(structure_data$Structure)  # Yapıyı monospace font ile göster
    )
  })
  
  output$nernaStructure <- renderUI({
    req(nerna_structures())
    
    # İlk sekansın yapısını göster
    structure_data <- nerna_structures()[1, ]
    
    div(
      p(strong("Sequence:"), structure_data$Sequence_Name),
      p(strong("Structure:"), structure_data$Structure),
      p(strong("MFE:"), paste0(structure_data$MFE, " kcal/mol"))
    )
  })
  
  # Benzer şekilde diğer yöntemler için de yapıları göster
  output$originalStructureShuffle <- renderUI({
    req(predict_structure())
    req(shuffleData())
    
    structure_data <- predict_structure()[1, ]
    
    div(
      p(strong("Sequence:"), structure_data$Sequence_Name),
      p(strong("Structure:"), structure_data$Structure),
      p(strong("MFE:"), paste0(structure_data$MFE, " kcal/mol"))
    )
  })
  
  output$shuffleStructure <- renderUI({
    req(shuffle_structures())
    
    structure_data <- shuffle_structures()[1, ]
    
    div(
      p(strong("Sequence:"), structure_data$Sequence_Name),
      p(strong("Structure:"), structure_data$Structure),
      p(strong("MFE:"), paste0(structure_data$MFE, " kcal/mol"))
    )
  })
  
  output$originalStructureDinuc <- renderUI({
    req(predict_structure())
    req(dinucleotideData())
    
    structure_data <- predict_structure()[1, ]
    
    div(
      p(strong("Sequence:"), structure_data$Sequence_Name),
      p(strong("Structure:"), structure_data$Structure),
      p(strong("MFE:"), paste0(structure_data$MFE, " kcal/mol"))
    )
  })
  
  output$dinucleotideStructure <- renderUI({
    req(dinucleotide_structures())
    
    structure_data <- dinucleotide_structures()[1, ]
    
    div(
      p(strong("Sequence:"), structure_data$Sequence_Name),
      p(strong("Structure:"), structure_data$Structure),
      p(strong("MFE:"), paste0(structure_data$MFE, " kcal/mol"))
    )
  })
  
  # Excel olarak indirme
  output$downloadExcel <- downloadHandler(
    filename = function() {
      paste0("rna_analysis_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      # Tüm sonuçları bir liste olarak topla
      all_results <- list(
        "Structure_Predictions" = predict_structure()
      )
      
      # Karşılaştırma sonuçları
      if(!is.null(comparisonData())) {
        all_results[["Comparison"]] <- comparisonData()
      }
      
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
  
  # RNAfold test fonksiyonu
  test_rnafold <- function() {
    # Test FASTA dosyası oluştur
    test_seq <- RNAStringSet("ACGACGUACGU")
    names(test_seq) <- "test_sequence"
    
    # Mevcut çalışma dizinini kaydet
    original_wd <- getwd()
    
    # Çalışma dizinini temp_dir'e değiştir
    setwd(temp_dir)
    
    test_file <- "test.fasta"
    writeXStringSet(test_seq, test_file)
    
    # RNAfold'un tam yolunu oluştur
    rnafold_path <- file.path(original_wd, "app/RNAfold.exe")
    output_dir <- file.path(temp_dir, "test_output")
    
    # Çıktı dizinini oluştur
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    message("Running RNAfold test...")
    message("Working directory: ", getwd())
    message("RNAfold path: ", rnafold_path)
    message("RNAfold exists: ", file.exists(rnafold_path))
    
    if (!file.exists(rnafold_path)) {
      stop("RNAfold.exe not found at: ", rnafold_path)
    }
    
    results <- system2(rnafold_path,
                      args = c(test_file), 
                      stdout = TRUE,
                      stderr = TRUE)
    
    # PS dosyasını doğru konuma taşı
    ps_file <- "test_sequence_ss.ps"
    if (file.exists(ps_file)) {
      file.rename(ps_file, file.path(output_dir, ps_file))
    }
    
    message("RNAfold output: ", paste(results, collapse = "\n"))
    
    # Çıktı dizinindeki dosyaları listele
    files <- list.files(output_dir, full.names = TRUE)
    message("Files in output directory: ", paste(files, collapse = ", "))
    
    # Çalışma dizinini geri al
    setwd(original_wd)
    
    return(results)
  }
  
  # Uygulama başladığında RNAfold'u test et
  observe({
    test_rnafold()
  })
  
  # Oturum sonlandığında temp dizinini temizle
  onSessionEnded(function() {
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
      dir.create(temp_dir)  # Yeni bir temiz dizin oluştur
    }
  })
}
