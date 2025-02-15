# função para construção de matrizes NEXUS usando extrapolação de áreas
Range_nexus <- function(pres_abs, meto = NULL){
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
  cat('#NEXUS',sep = '\n')
  cat('\n')
  cat('begin data;',sep = '\n')
  cat(paste0('dimensions ntax=', n.taxa, ' nchar=', n.cara, ';'),sep = '\n')
  cat('format datatype=standard symbols="01" gap=-;',sep = '\n')
  cat('CHARSTATELABELS',sep = '\n')
  #cat('/* Matrix for a PAE-PCE analysis */', sep = '\n')
  
  # colocando a lista de nomes:
  species <- colnames(mat)
  species <- gsub(pattern = " ", replacement = "_", x = species)
  
  for(k in 1:n.cara - 1){
    cat(sprintf("%s %s,", k, species[k]), sep = '\n')
  }
  cat(sprintf("%s %s;", k, species[n.cara]), sep = '\n')
  cat('\n')
  
  # ROOT:
  # cat(paste0(c('ROOT', ' ', rep(0, n.cara)), collapse = ''), sep = '\n')
  
  #colocando a matriz de caracteres pres_abs...
  cat('matrix', sep = '\n')
  
  for(j in 1:n.cara){
    # nomes[j] <- paste(c(taxa[j]),collapse = ', ')
    # cat(species[j], mat[, j], fill = T)
    cat(paste(mat[, j], collapse = ''), fill = T,
        labels = paste0(species[j]))
  }
  
  cat(';', sep = '\n')
  
  cat('end;')
  
  sink()
  close(zz)
  
  #abrindo, escrevendo e fechando o arquivo...
  dir.create('out_mat')
  if(meto == 'MPC'){
    exte <- 'out_mat/pres_abs_MPC.nex'
    } else if(meto == 'BUFF'){
      exte <- 'out_mat/pres_abs_BUFF.nex'
    } else if(meto == 'MST'){
      exte <- 'out_mat/pres_abs_MST.nex'
    }
  tempor = file(exte, "wt")
  writeLines(texto, con = tempor)
  close(tempor)
}
