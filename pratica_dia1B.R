# prática 1B:

# Dados de Panthera onca:
# 1) Leitura dos dados a partir de um arquivo .txt ou .csv ou .tsv
dat_teste <- read_tsv("onca_original.csv", guess_max = 25000, quote = "") # lendo, no máximo, 25 mil coordenadas
dat_teste # 6236 registros

library(CoordinateCleaner)

# Pronto! Agora, é só rodar...
flagsA <- clean_coordinates(dat_teste)
flagsA
##########################################################
##########################################################
# de forma alternativa, você pode fazer isso:
dat_teste <- data.frame(dat_teste)
flagsB <- clean_coordinates(x = dat_teste, lon = "decimalLongitude", lat = "decimalLatitude", countries = "countryCode",
                            species = "species", tests = c("equal", "gbif", "zeros",
                                                           "seas"), seas_ref = buffland)
##########################################################

# podemos ver os pontos que ficaram e os que foram retirados:
plot(flagsA, lon = "decimalLongitude", lat = "decimalLatitude")

# podemos ver quantos registros foram excluídos:
sum(flagsA$.summary) # 5.768 pontos

# Agora, é só remover dos dados originais:
dat_teste_novo <- dat_teste[flagsA$.summary, ]

# Salve o dado novo e limpo:
write.csv(dat_teste_novo, 'panthera_occ_limpo.csv')