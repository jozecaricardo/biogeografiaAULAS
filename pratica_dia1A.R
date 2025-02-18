# PRATICA 1A: databio
rm(list = lm())

install.packages('rgbif') # instalação do rgbif

?rm # help da função

library(rgbif) # para carregar o pacote

usr <- "nono" # nome do usuário

pwd <- "nonon" # senha do usuário

eml <- "nonon@gmail.com" # email do usuário

tax <- occ_search(scientificName = "Altiphylax") # você pode usar
# qualquer outro nome
tax

# taxonkey
tpred <- pred("taxonKey", name_backbone("Altiphylax")$usageKey) # predicate and
# execute
tpred

# baixando os dados no GBIF
dload <- occ_download(tpred, format = "SIMPLE_CSV", pred("hasCoordinate", TRUE),
    user = usr, pwd = pwd, email = eml) # Para fazer o download... Pode demorar
#   um pouco dependendo do número de registros!
dload

# usando os IDs:
ID <- taxize::get_gbifid_("Lethocerus patruelis", method = "backbone") ## get_gbifid_
# retornará uma lista
ID

# carregando o pacote tydiverse
library(tidyverse)
# uma pequena filtragem
ID <- ID %>%
  bind_rows() %>% # convertida em data frame
  filter(matchtype == "EXACT" & status == "ACCEPTED") # filtra os dados para os nomes aceitos
ID

# IDs de grupos taxonômicos supra-específicos
ID_Belos <- taxize::get_gbifid_("Belostomatidae", method="backbone") %>%
  bind_rows() %>%
  filter(matchtype == "EXACT" & status == "ACCEPTED")
ID_Belos

# para saber onde o R está "enxergando":
getwd()

# lista de espécies:
splist <- read_csv("bombacoideae_specieslist.csv") # espécies em .csv de Bombacoidea 
splist
gbif_taxon_keys <-
  splist %>%
  pull("accepted_name") %>% # use poucos nomes se for um teste
  taxize::get_gbifid_(method="backbone") %>%
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # coloque os nomes originais no data frame
  bind_rows() %>% # combine todos os data frames em apenas um
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>%
  filter(kingdom == "Plantae") # evite homonímias

gbif_taxon_keys

# Usando polígonos:
# para usar o taxise, carregar primeiro o devtools e depois instalar o taxsize
library(devtools)
install_github('ropensci/taxize') # atenção: só faça isso uma vez! 
library(sf); library(terra); library(tidyverse); library(taxize)

amz <- st_read("Amazonia.kml") # puxando o arquivo para dentro do R
plot(amz$geometry) # para verificar se está tudo bem
amz <- st_as_text(st_as_sfc(amz)) # para transformar num objeto WKT

ID_Belos <- taxize::get_gbifid_("Belostomatidae", method="backbone") %>% 
  bind_rows() %>%
  filter(matchtype == "EXACT" & status == "ACCEPTED") # para pegar o ID da família Belostomatidae
ID_Belos

dload <- occ_download(pred_in("taxonKey", ID_Belos$usagekey), pred("hasCoordinate",
  TRUE), pred_within(amz), format = "SIMPLE_CSV",
  user = usr, pwd = pwd, email = eml) # o argumento pred_within fará o trabalho

occ_download_wait(dload)
##########################################################
##########################################################
# PBDB:
library(paleobioDB); library(readr); library(dplyr) # carregando pacotes 
# Para pegar as ocorrências do Paleobiology DB, estreite sua busca como abaixo:
dat_tax <- paleobioDB::pbdb_occurrences(base_name = c("Belostomatidae"),
    vocab = "pbdb", limit = 500, show = c("coords", "phylo", "attr", "loc",
     "time", "rem", "ident"))
dat_tax
rownames(dat) <- NULL

# É possível também pegar registro de uma lista de espécies ou gêneros
gen_list <- read_csv("bombacoideae_genuslist.csv") # nome do arquivo .csv
gen_list
dat_gen <- paleobioDB::pbdb_occurrences(base_name = gen_list$gen,
    vocab = "pbdb", limit = 500, show = c("coords", "phylo", "attr", "loc", "time", "rem", "ident"))
dat_gen
# combine os dados em um local só:
dat <- bind_rows(dat_tax, dat_gen)
dat

# Mas tome cuidado com as homonímias!
unique(dat$phylum)
dat <- dat[dat$phylum == "Spermatophyta", ]

# Visualize Belostomatidae
pbdb_map(dat_tax)
pbdb_richness(dat_tax, rank = "species",
      temporal_extent = c(100, 0)) # mapa de riqueza                

# grave no arquivo:
write_csv(dat_tax, "BelosFossil_records.csv")


