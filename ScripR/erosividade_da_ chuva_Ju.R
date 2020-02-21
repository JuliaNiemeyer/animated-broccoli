

# aqui aplico a funcao da erosividade para cada regiao do cerrado


library(raster)
library(maptools)
library(sp)
library(dplyr)
library(rgdal)


# carrega raster altitude TerraClass


alt_cerrado <- raster("./Mapas_base/dem_buff.tif")


# PRECISA passar para WGS -- precisa estar em latlong
alt_cerrado <- projectRaster(alt_cerrado, crs = "+proj=longlat +datum=WGS84 +towgs84=0,0,0")
#alt_cerrado

##plot(alt_cerrado)
# carrega shape estados - Cerrado
# tenho que separar os dados de altitude por regiao do Brasil, para aplicar a formula de Mello et al 2013

states <- readOGR('./Mapas_base', "UF_BSF")

# passa para WGS
states <- spTransform(states, crs(alt_cerrado))

##plot(states, add=T)


# separa por regiao

# sudeste
SE <- states[states$NM_REGIAO == "SUDESTE",]
##plot(SE)
#Project

# norte
N <- states[states$NM_REGIAO == "NORTE",]
#plot(N, add = T)

# nordeste
NE <- states[states$NM_REGIAO == "NORDESTE",]
##plot(NE, add = T)

# centro-oeste
CO <- states[states$NM_REGIAO == "CENTRO-OESTE",]
##plot(CO, add = T)

##N?o tem SUL

# corta o raster de altitude para cada regiao

#alt_cerrado # raster de altitude do Cerrado

# sudeste
alt_SE <- crop(alt_cerrado, SE) # coloca no extend do SE
alt_SE <- mask(alt_SE, SE) # coloca NA fora dos estados do SE
###plot(alt_SE)

# norte
alt_N <- crop(alt_cerrado, N) # coloca no extend do SE
alt_N <- mask(alt_N, N) # coloca NA fora dos estados do SE
###plot(alt_N)

# nordeste
alt_NE <- crop(alt_cerrado, NE) # coloca no extend do SE
alt_NE <- mask(alt_NE, NE) # coloca NA fora dos estados do SE
##plot(alt_NE)

# centro-oeste
alt_CO <- crop(alt_cerrado, CO) # coloca no extend do SE
alt_CO <- mask(alt_CO, CO) # coloca NA fora dos estados do SE
##plot(alt_CO)



# Prepara o objeto com xy e altitude de cada regiao
# precisa de lat, lot e altitude calcular a erosividade da chuva

# sudeste
# extrai as coordenadas
xy_SE <- coordinates(alt_SE)
xy_SE <- as.matrix(xy_SE)
#xy_SE
#nrow(xy_SE) # 404790
##head(xy_SE)

# extrai os valores altitude
val_SE <- getValues(alt_SE)
val_SE <- as.matrix(val_SE)
#table(val_SE)
#length(val_SE) # 404790
##head(val_SE)


# junta coord e altitude
xyz_SE <- cbind(xy_SE, val_SE)
#xyz_SE  # objeto com x, y e altitude

##head(xyz_SE)


# norte
# extrai as coordenadas
xy_N <- coordinates(alt_N)
#xy_N
#nrow(xy_N)  ## 32076

# extrai altitude
val_N <- getValues(alt_N)
val_N <- as.matrix(val_N)
#length(val_N) # 32076

# junta coord e altitude
xyz_N <- cbind(xy_N, val_N)
#xyz_N  # objeto com x, y e altitude

#tail(xyz_N)


# nordeste
# extrai as coordenadas
xy_NE <- coordinates(alt_NE)
#xy_NE
#nrow(xy_NE) #  1030083

# extrai altitude
val_NE <- getValues(alt_NE)
val_NE <- as.matrix(val_NE)
#length(val_NE) #  1030083

# junta coord e altitude
xyz_NE <- cbind(xy_NE, val_NE)
#xyz_NE  # objeto com x, y e altitude

#head(xyz_NE)
#nrow(xyz_NE) #1030083


# centro-oeste
# extrai as coordenadas
xy_CO <- coordinates(alt_CO)
#nrow(xy_CO) #72150

# extrai altitude
val_CO <- getValues(alt_CO)
val_CO <- as.matrix(val_CO)
#length(val_CO) # 72150

# junta coord e altitude
xyz_CO <- cbind(xy_CO, val_CO)
#xyz_CO  # objeto com x, y e altitude

#head(xyz_CO)
#nrow(xyz_CO) #72150

# carrega tabela com formulas
table <- read.table('./tabela/formulas_erosividade.txt', header = T)


##head(table)



####
# Faz as funcoes para cada regiao

### sudeste

# define a regiao que vai usar
##head(table)
val <- table$Southeast
# transforma em matriz
val <- as.matrix(val)
##head(val)



# constroi funcao
# mudar o nome da funcao
fun_SE <- function(A, LA, LO){

  # tem que tirar os coeficientes que forem NA

  # pega os valores de cada coeficiente
  intercept <- val[1]
  c_A <- val[2]
  c_LA <- val[3]
  #c_LO <- val[4] #NA
  c_A_qua <- val[5]
  c_LA_qua <- val[6]
  c_LO_qua <- val[7]
  # pula 5 NA
  c_LOxA <- val[13]
  c_LAxLO <- val[14]
  c_LO_quaxA <- val[15]
  c_LO_quadxA_qua <- val[16]
  c_LO_quaxLA_qua <- val[17]
  c_LA_quaxLO_tres <- val[18]

  # pega a conta de cada um
  v_A <- A
  v_LA <- LA
  #v_LO <- LO NA
  v_A_qua <- A^2
  v_LA_qua <- LA^2
  v_LO_qua <- LO^2
  # pula 5 NA
  v_LOxA <- LO*A
  v_LAxLO <- LA*LO
  v_LO_quaxA <- LO^2*A
  v_LO_quadxA_qua <- LO^2*A^2
  v_LO_quaxLA_qua <- LO^2*LA^2
  v_LA_quaxLO_tres <- LA^2*LO^3


  # multiplica o coeficiente pela conta de cada um
  r_A <- c_A * v_A
  r_LA <- c_LA * v_LA
  #r_LO <- c_LO * v_LO
  r_A_qua <- c_A_qua * v_A_qua
  r_LA_qua <- c_LA_qua * v_LA_qua
  r_LO_qua <- c_LO_qua * v_LO_qua
  r_LOxA <- c_LOxA * v_LOxA
  r_LAxLO <- c_LAxLO * v_LAxLO
  r_LO_quaxA <- c_LO_quaxA * v_LO_quaxA
  r_LO_quadxA_qua <- c_LO_quadxA_qua * v_LO_quadxA_qua
  r_LO_quaxLA_qua <- c_LO_quaxLA_qua * v_LO_quaxLA_qua
  r_LA_quaxLO_tres <- c_LA_quaxLO_tres * v_LA_quaxLO_tres

    # faz objeto com valor para salvar
    # retirei  "r_LO" NA
    valor <- intercept + r_A + r_LA + r_A_qua + r_LA_qua +
    r_LO_qua + r_LOxA + r_LAxLO + r_LO_quaxA + r_LO_quadxA_qua +
    r_LO_quaxLA_qua + r_LA_quaxLO_tres

  return(valor)

} # fecha a funcao



# cria argumentos para a funcao

# aqui define os argumentos
# mudar a regiao
##head(xyz_SE)
A_SE <- xyz_SE[,3] # altitude
LA_SE <- xyz_SE[,2] # latitude
LO_SE <- xyz_SE[,1] # longitude

# roda a funcao
erod <- fun_SE(A = A_SE, LA_SE, LO_SE)

#erod
#length(erod) #1337560
##plot(erod_SE)
#unique(erod)



# faz o raster da erodibilidade

# pega xy
# mudar aqui a regiao
xy <- xyz_SE[,-3]
#xy
##head(xy)

# coloca o valor da erodibildiade
xy_erod <- cbind(xy, erod)
#xy_erod
##head(xy_erod)

# cria raster
# mudar o nome do raster
r_erod_SE <- rasterFromXYZ(xy_erod)
#r_erod_SE # 1337560

##plot(r_erod_SE)





### norte

# define a regiao que vai usar
##head(table)
val <- table$North.Midwest
#val

# transforma em matriz
val <- as.matrix(val)
##head(val)



# constroi funcao
# mudar o nome da funcao
fun_N <- function(A, LA, LO){

  # tem que tirar os coeficientes que forem NA

  # pega os valores de cada coeficiente
  intercept <- val[1]
  c_A <- val[2]
  c_LA <- val[3]
  c_LO <- val[4]
  c_A_qua <- val[5]
  c_LA_qua <- val[6]
  c_LO_qua <- val[7]
  # pula 5 NA
  c_LOxA <- val[13]
  c_LAxLO <- val[14]
  c_LO_quaxA <- val[15]
  c_LO_quadxA_qua <- val[16]
  c_LO_quaxLA_qua <- val[17]
  c_LA_quaxLO_tres <- val[18]

  # pega a conta de cada um
  v_A <- A
  v_LA <- LA
  v_LO <- LO
  v_A_qua <- A^2
  v_LA_qua <- LA^2
  v_LO_qua <- LO^2
  # pula 5 NA
  v_LOxA <- LO*A
  v_LAxLO <- LA*LO
  v_LO_quaxA <- LO^2*A
  v_LO_quadxA_qua <- LO^2*A^2
  v_LO_quaxLA_qua <- LO^2*LA^2
  v_LA_quaxLO_tres <- LA^2*LO^3


  # multiplica o coeficiente pela conta de cada um
  r_A <- c_A * v_A
  r_LA <- c_LA * v_LA
  r_LO <- c_LO * v_LO
  r_A_qua <- c_A_qua * v_A_qua
  r_LA_qua <- c_LA_qua * v_LA_qua
  r_LO_qua <- c_LO_qua * v_LO_qua
  r_LOxA <- c_LOxA * v_LOxA
  r_LAxLO <- c_LAxLO * v_LAxLO
  r_LO_quaxA <- c_LO_quaxA * v_LO_quaxA
  r_LO_quadxA_qua <- c_LO_quadxA_qua * v_LO_quadxA_qua
  r_LO_quaxLA_qua <- c_LO_quaxLA_qua * v_LO_quaxLA_qua
  r_LA_quaxLO_tres <- c_LA_quaxLO_tres * v_LA_quaxLO_tres

  # faz objeto com valor para salvar
  # retirei NA
  valor <-
    intercept +

    r_LA +
    r_LO +
    r_A_qua +
    r_LA_qua +
    r_LO_qua +

    r_LAxLO +
    r_LO_quaxA +
    r_LO_quadxA_qua +
    r_LO_quaxLA_qua +
    r_LA_quaxLO_tres

  return(valor)

} # fecha a funcao



# cria argumentos para a funcao

# aqui define a regiao
#head(xyz_N)
A_N <- xyz_N[,3] # altitude
LA_N <- xyz_N[,2] # latitude
LO_N <- xyz_N[,1] # longitude

# roda a funcao
erod <- fun_N(A = A_N, LA_N, LO_N)

#erod
#length(erod) #1628433
#unique(erod)



# faz o raster da erosividade

# pega xy
# mudar aqui a regiao
xy <- xyz_N[,-3]
#xy
#head(xy)

# coloca o valor da erodibildiade
xy_erod <- cbind(xy, erod)

#head(xy_erod)

# cria raster
# mudar o nome do raster
r_erod_N <- rasterFromXYZ(xy_erod)
#r_erod_N # 1337560


##plot(r_erod_N)





### nordeste

# define a regiao que vai usar
##head(table)
val <- table$Northeast
#val

# transforma em matriz
val <- as.matrix(val)
#val



# constroi funcao
# mudar o nome da funcao
fun_NE <- function(A, LA, LO){

  # tem que tirar os coeficientes que forem NA

  # pega os valores de cada coeficiente
  intercept <- val[1]
  c_A <- val[2]
  c_LA <- val[3]
  c_LO <- val[4]
  c_A_qua <- val[5]
  c_LA_qua <- val[6]
  c_LO_qua <- val[7]
  # pula 5 NA
  c_LOxA <- val[13]
  c_LAxLO <- val[14]
  c_LO_quaxA <- val[15]
  c_LO_quadxA_qua <- val[16]
  c_LO_quaxLA_qua <- val[17]
  c_LA_quaxLO_tres <- val[18]

  # pega a conta de cada um
  v_A <- A
  v_LA <- LA
  v_LO <- LO
  v_A_qua <- A^2
  v_LA_qua <- LA^2
  v_LO_qua <- LO^2
  # pula 5 NA
  v_LOxA <- LO*A
  v_LAxLO <- LA*LO
  v_LO_quaxA <- LO^2*A
  v_LO_quadxA_qua <- LO^2*A^2
  v_LO_quaxLA_qua <- LO^2*LA^2
  v_LA_quaxLO_tres <- LA^2*LO^3


  # multiplica o coeficiente pela conta de cada um
  r_A <- c_A * v_A
  r_LA <- c_LA * v_LA
  r_LO <- c_LO * v_LO
  r_A_qua <- c_A_qua * v_A_qua
  r_LA_qua <- c_LA_qua * v_LA_qua
  r_LO_qua <- c_LO_qua * v_LO_qua
  r_LOxA <- c_LOxA * v_LOxA
  r_LAxLO <- c_LAxLO * v_LAxLO
  r_LO_quaxA <- c_LO_quaxA * v_LO_quaxA
  r_LO_quadxA_qua <- c_LO_quadxA_qua * v_LO_quadxA_qua
  r_LO_quaxLA_qua <- c_LO_quaxLA_qua * v_LO_quaxLA_qua
  r_LA_quaxLO_tres <- c_LA_quaxLO_tres * v_LA_quaxLO_tres

  # faz objeto com valor para salvar
  # retirei NA
  valor <-
  intercept + r_A + r_LO + r_A_qua + r_LA_qua +
    r_LOxA + r_LO_quaxA + r_LO_quaxLA_qua + r_LA_quaxLO_tres

  return(valor)

} # fecha a funcao



# cria argumentos para a funcao

# aqui define a regiao
##head(xyz_NE)
A_NE <- xyz_NE[,3] # altitude
LA_NE <- xyz_NE[,2] # latitude
LO_NE <- xyz_NE[,1] # longitude

# roda a funcao
erod <- fun_NE(A = A_NE, LA_NE, LO_NE)

#length(erod) #1175675



# faz o raster da erodibilidade

# pega xy
# mudar aqui a regiao
xy <- xyz_NE[,-3]
##head(xy)

# coloca o valor da erodibildiade
xy_erod <- cbind(xy, erod)
##head(xy_erod)

# cria raster
# mudar o nome do raster
r_erod_NE <- rasterFromXYZ(xy_erod)
#r_erod_NE # 1175675

##plot(r_erod_NE)




### centro-oeste

# define a regiao que vai usar
##head(table)
val <- table$North.Midwest

# transforma em matriz
val <- as.matrix(val)
#val



# constroi funcao
# mudar o nome da funcao
fun_CO <- function(A, LA, LO){

  # tem que tirar os coeficientes que forem NA

  # pega os valores de cada coeficiente
  intercept <- val[1]
  c_A <- val[2]
  c_LA <- val[3]
  c_LO <- val[4]
  c_A_qua <- val[5]
  c_LA_qua <- val[6]
  c_LO_qua <- val[7]
  # pula 5 NA
  c_LOxA <- val[13]
  c_LAxLO <- val[14]
  c_LO_quaxA <- val[15]
  c_LO_quadxA_qua <- val[16]
  c_LO_quaxLA_qua <- val[17]
  c_LA_quaxLO_tres <- val[18]

  # pega a conta de cada um
  v_A <- A
  v_LA <- LA
  v_LO <- LO
  v_A_qua <- A^2
  v_LA_qua <- LA^2
  v_LO_qua <- LO^2
  # pula 5 NA
  v_LOxA <- LO*A
  v_LAxLO <- LA*LO
  v_LO_quaxA <- LO^2*A
  v_LO_quadxA_qua <- LO^2*A^2
  v_LO_quaxLA_qua <- LO^2*LA^2
  v_LA_quaxLO_tres <- LA^2*LO^3


  # multiplica o coeficiente pela conta de cada um
  r_A <- c_A * v_A
  r_LA <- c_LA * v_LA
  r_LO <- c_LO * v_LO
  r_A_qua <- c_A_qua * v_A_qua
  r_LA_qua <- c_LA_qua * v_LA_qua
  r_LO_qua <- c_LO_qua * v_LO_qua
  r_LOxA <- c_LOxA * v_LOxA
  r_LAxLO <- c_LAxLO * v_LAxLO
  r_LO_quaxA <- c_LO_quaxA * v_LO_quaxA
  r_LO_quadxA_qua <- c_LO_quadxA_qua * v_LO_quadxA_qua
  r_LO_quaxLA_qua <- c_LO_quaxLA_qua * v_LO_quaxLA_qua
  r_LA_quaxLO_tres <- c_LA_quaxLO_tres * v_LA_quaxLO_tres

  # faz objeto com valor para salvar
  # retirei NA
  valor <-
    intercept +

    r_LA +
    r_LO +
    r_A_qua +
    r_LA_qua +
    r_LO_qua +

    r_LAxLO +
    r_LO_quaxA +
    r_LO_quadxA_qua +
    r_LO_quaxLA_qua +
    r_LA_quaxLO_tres

  return(valor)

} # fecha a funcao



# cria argumentos para a funcao

# aqui define a regiao
##head(xyz_CO)
A_CO <- xyz_CO[,3] # altitude
LA_CO <- xyz_CO[,2] # latitude
LO_CO <- xyz_CO[,1] # longitude

# roda a funcao
erod <- fun_CO(A = A_CO, LA_CO, LO_CO)

# faz o raster da erodibilidade

# pega xy
# mudar aqui a regiao
xy <- xyz_CO[,-3]
##head(xy)

# coloca o valor da erodibildiade
xy_erod <- cbind(xy, erod)
##head(xy_erod)

# cria raster
# mudar o nome do raster
r_erod_CO <- rasterFromXYZ(xy_erod)
#r_erod_CO # 2514900

##plot(r_erod_CO)


# junta os raster de erodibilidade

r_erod_cerr <- merge(r_erod_SE, r_erod_N, r_erod_NE, r_erod_CO)

#r_erod_cerr #5718582
##plot(r_erod_cerr)

#alt_cerrado #6461380 # Parece que o numero de cells e diferente por causa do extend  '
##plot(alt_cerrado)

##ver o que acontece at? aqui sem transformar de novo (tentar transformar no arcgis)

#Coloca no extent de LULC
tc <- raster("./Mapas_base/LULC_bufferProj.tif")
##plot(tc)
#extent(r_erod_cerr)<-extent(tc) #tc ? o raster de LULC

##Projetar pra WGS
proj4string(r_erod_cerr) <- CRS("+proj=longlat +datum=WGS84 +towgs84=0,0,0")
#crs(r_erod_cerr) <- "+proj=longlat +datum=WGS84 +towgs84=0,0,0"
#erosiv <- alignExtent(r_erod_cerr, tc)
erosiv2 <- projectRaster(r_erod_cerr, tc, method = "bilinear")


# exporta raster

writeRaster(x = erosiv2, "./Output/erosividade_WGS2", format = "GTiff", overwrite = T)

