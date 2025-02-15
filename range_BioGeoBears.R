# função para construção de matrizes NEXUS usando extrapolação de áreas
Range_BioGeoBears <- function(pres_abs, meto = NULL){
  if(class(pres_abs)[1] == 'matrix'){
    mat <- pres_abs
  } else if(class(pres_abs)[1] == 'character'){
    mat <- read.table(file = pres_abs, sep = '', dec = '.', row.names = 1) # presence-absence matrix
  }
  
  #extraindo os número das células do raster
  taxa <- row.names(mat)
  
  #produzindo os arquivos *.txt
  n.cara <- ncol(mat)
  n.taxa <- nrow(mat)
  
  # substituindo pontos por sublinhas
  colnames(mat) <- gsub(pattern = '\\.', replacement = '_', x = colnames(mat))
  
  nomes <- NULL
  
  conta <- 0
  
  # capture R output:
  zz <- textConnection('texto', "w")
  sink(zz)
  cat(paste(n.cara, n.taxa, paste0('(', paste(rownames(pres_abs), collapse = ' '), ')'), sep = ' '), sep = '\n')
  
  species <- colnames(mat)
  species <- gsub(pattern = " ", replacement = "_", x = species)
  
  #colocando a matriz de caracteres pres_abs...
  for(j in 1:n.cara){
    # nomes[j] <- paste(c(taxa[j]),collapse = ', ')
    # cat(species[j], mat[, j], fill = T)
    cat(paste(mat[, j], collapse = ''), fill = T,
        labels = paste0(species[j]))
  }
  sink()
  close(zz)
  
  #abrindo, escrevendo e fechando o arquivo...
  dir.create('out_mat')
  if(meto == 'MPC'){
    exte <- 'out_mat/pres_abs_MPC_geog.data'
    } else if(meto == 'BUFF'){
      exte <- 'out_mat/pres_abs_BUFF_geog.data'
    } else if(meto == 'MST'){
      exte <- 'out_mat/pres_abs_MST_geog.data'
    }
  tempor = file(exte, "wt")
  writeLines(texto, con = tempor)
  close(tempor)
}