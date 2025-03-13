# Gerekli kütüphaneleri yükle
source("global.R")

# UI ve Server dosyalarını yükle
source("ui.R")
source("server.R")

# Uygulamayı çalıştır
shinyApp(ui = ui, server = server)
