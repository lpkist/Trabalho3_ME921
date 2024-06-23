library(rvest)
library(xml2)
library(tidyverse)

#https://www.dac.unicamp.br/sistemas/catalogos/grad/catalogo2024/disciplinas/
siglas <- c("AC", "AD", "AG", "AI", "AM", "AP", "AR", "AU", "AX", "BA",
            "BB", "BC", "BD", "BE", "BF", "BG", "BH", "BI", "BL", "BM",
            "BP", "BS", "BT", "BV", "BX", "BZ", "CE", "CG", "CP", "CS", 
            "CV", "CX", "DB", "DC", "DE", "DM", "DO", "DS", "EA", "EB", 
            "EE", "EF", "EG", "EL", "EM", "EN", "EP", "EQ", "ER", "ES",
            "ET", "EU", "EX", "F", "FA", "FI", "FL", "FM", "FN", "FR",
            "FT", "FX", "GE", "GF", "GL", "GM", "GN", "GT", "HG", "HH",
            "HL", "HZ", "IX", "LA", "LE", "LG", "LT", "MA", "MC", "MD",
            "ME", "MG", "ML", "MP", "MS", "MU", "MX", "NC", "NT", "PF",
            "PG", "QA", "QE", "QF", "QG", "QI", "QL", "QO", "SI", "SL",
            "ST", "TA", "TL", "TT")

siglas <- str_to_lower(siglas)
pre_reqs <- map(siglas, function(sigla){
  if(sigla == "f"){
    url <- "https://www.dac.unicamp.br/sistemas/catalogos/grad/catalogo2024/disciplinas/f_.html"
  } else{
    url <- paste0("https://www.dac.unicamp.br/sistemas/catalogos/grad/catalogo2024/disciplinas/", sigla, ".html")
  }
  download.file(url, destfile = "pagina_sigla.html", quiet=TRUE)
  a <- read_html("pagina_sigla.html")
  aux <- html_elements(a, 'div[class="small-12 columns pad-content"]')
  pre_req <- map(aux, function(materia){
    html_elements(materia, 'p')[7] %>% html_text() %>% str_extract_all("[A-Z]{2}[0-9]{3}") %>% data.frame
  })
  siglas_pre_req <- map(pre_req, function(x)data.frame(materia = x[[1]])) %>% 
    list_rbind() %>% 
    reframe(siglas_pre_req = str_sub(materia, start = 1, end = 2)) %>%
    unique() %>% 
    drop_na()
  fisica_pre_req <- map(aux, function(materia){
    html_elements(materia, 'p')[7] %>% html_text() %>% str_extract("F [0-9]{3}") %>% data.frame
  })
  siglas_pre_req <- ifelse((fisica_pre_req %>%
                              list_rbind() %>%
                              drop_na() %>%
                              nrow()) > 0,
                           c(siglas_pre_req, "F"), siglas_pre_req)
  siglas_pre_req
})
names(pre_reqs) <- str_to_upper(siglas)
matriz <- matrix(0, length(siglas), length(siglas))
colnames(matriz) <- rownames(matriz) <- str_to_upper(siglas)

for(i in str_to_upper(siglas)){
  pre_reqs_i <- pre_reqs[[i]][[1]]
  for(j in str_to_upper(siglas)){
    matriz[j, i] <- as.numeric(j %in% pre_reqs_i)
  }
}

for(i in str_to_upper(siglas)){
  if((sum(matriz[i,], na.rm = T) + sum(matriz[,i], na.rm = T) - 2*matriz[i,i])==0){
    matriz[i,] <- matriz[,i] <- NA
  }
}


matriz <- matriz[rowSums(is.na(matriz)) != ncol(matriz), rowSums(is.na(matriz)) != ncol(matriz)]
diag(matriz) <- 0
save(matriz, file = "matriz_adj.RData")
