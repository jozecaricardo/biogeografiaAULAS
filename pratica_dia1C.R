# prática 1C: paleodatabase

library(paleobioDB); library(readr); library(dplyr) # carregando pacotes 

# Para pegar as ocorrências do Paleobiology DB, estreite sua busca como abaixo:
dat_tax <- paleobioDB::pbdb_occurrences(base_name = c("Belostomatidae"),
          vocab = "pbdb", limit = 500, show = c("coords", "phylo", "attr", "loc",
          "time", "rem", "ident"))
rownames(dat) <- NULL

# É possível também pegar registro de uma lista de espécies ou gêneros
gen_list <- read_csv("bombacoideae_genuslist.csv") # nome do arquivo .csv
dat_gen <- paleobioDB::pbdb_occurrences(base_name = gen_list$gen, vocab = "pbdb", limit = 500, show = c("coords", "phylo", "attr", "loc", "time", "rem", "ident"))

# combine os dados em um local só:
dat <- bind_rows(dat_tax, dat_gen)

# Mas tome cuidado com as homonímias!
unique(dat$phylum)
dat <- dat[dat$phylum == "Spermatophyta", ]

# Visualize Belostomatidae
pbdb_map(dat_tax)
pbdb_richness(dat_tax, rank = "species", temporal_extent = c(100, 0)) # mapa de riqueza                

# grave no arquivo:
write_csv(dat_tax, "BelosFossil_records.csv")