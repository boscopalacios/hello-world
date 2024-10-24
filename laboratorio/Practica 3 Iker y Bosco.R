install.packages("ape")
install.packages("phangorn")
install.packages("phytools")
library(ape)
library(phangorn)
library(phytools)
fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta", 
                        format = "FASTA", type = "AA")
#Mediante la función dist.aa calculamos una matriz de distancias por pares de las secuencias de aminoácidos partir de una objecto de clase AAbin.
matrizdist <- as.AAbin(fraxatin)
matrizdist <- dist.aa(matrizdist)
matrizdist

#vamos a crear un árbol con el método de grupo de pares no ponderados con media aritmética (UPGMA) mediante la matriz de distancia ya calculada
arbolUPGMA <- upgma(matrizdist)
plot(arbolUPGMA)

# creamos otro arbol pero con el metodo de union de unión de vecinos (NJ) con la misma matriz de distancias
arbolNJ <- nj(matrizdist)
plot(arbolNJ)

# Para modificar los arboles podemos añadir parametros como cex (tamaño letra) y edge (color)
plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="red", font=3)
plot(arbolUPGMA, type= "c", cex=0.8, edge.width=2, edge.color="yellow", font=3)
plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=1, node.pos=2, cex=0.8, edge.width=2, edge.color="purple", font=3)

# Podemos hacer graficos de arboles con plotTree
plotTree(arbolNJ)
#Podemos modificarlo otra vez cambiando los colores como en anteriore casos
plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="orange", lwd=2)
# Tambien podemos cambiar el orden de visualizacion
plotTree(ladderize(arbolNJ))
# Para guardar el árbol podemos usar el siguiente comando y con otro comando para leerlo
write.tree(arbolNJ, file = "file_name.nex")
read.tree(file = "file_name.nex")
#Enraizamos los árboles con la función root
arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)
#Se puede hacer lo mismo con un árbol creado por UPGMA
arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
plot(arbolUPGMAraiz)
#Podemos visualizar los dos árboles a la vez con estos comandos
layout(matrix(c(1,2)), height=c(10,10))
par(mar=c(1,1,1,1))
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)
#Se hace el número de pasos necesarios del árbol
parsimony(arbolUPGMAraiz, fraxatin)
#Hacemos lo mismo pero sin la raíz
parsimony(arbolUPGMA, fraxatin)
#Obtenemos el árbol con mejor parsimonia con el siguiente comando
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)
#Hacemos lo mismo con el árbol NJ
mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin)
#Usamos el algoritmo pratchet que es más complejo para buscar el árbol con mejor parsimonia
fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)
fraxatin_parsimonia
#Comparamos y enraizamos
fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex = 0.6)
# Hacemos un 100% de concenso entre las ramas
estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
plot(estrictode100, cex = .6)
estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
plot(estrictode30, cex = .6)

# Podemos hacer una serie de pseudoréplicas con bootstrapping de una matriz
arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)

# Vamos a usar la funcion pratched y generamos replicas
plot(arbolesbootstrap, cex = .6)

# Creamos un consenso del 60%
estricto60 <- consensus(arbolesbootstrap, p = 0.6)
plot(estricto60, cex = .6)

# Con nuestro objeto fraxatin creamos un árbol al azar de 11 ramas con rtree como punto de partida.
arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
plot(arbolazar, cex = .5)

# Enraizamos con ornitorrinco y escalerizamos
arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex = .5); add.scale.bar()

# Calculamos la verosimilitud del árbol dadas las secuencias y usamos pml para verlo
ajustado <- pml(arbolazarR, fraxatin)
ajustado

# Encontramos un árbol que optimice la verosimilitud
ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")
# Vemos el arbol oculto y lo enraizamos
ajustadoconDay$tree

ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()

#Podemos usar diferentes modelos
ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")
ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")

#Calculamos usando el Criterio de información de Akaike AIC
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)
# JTT modelo de Jones-Taylor-Thornton para evaluar distancias
mejorarbol <- optim.pml(
  object = ajustadoconDay, 
  model = "JTT", 
  rearrangement = "ratchet")
mejorarbol

mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()