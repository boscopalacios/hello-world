library(tidyverse)
data(starwars)
# Elegimos a todas las columnas menos el nombre 
starwars %>% select(-name)
# Cogemos unicamente las columnas que tienen subraya (_)
starwars %>% select(contains("_"))
# Cogemos solo las columnas que empiezan con la "s"
starwars %>% select(starts_with("s"))
# Creamos una data frame con los nombres y con su planeta de origen (homeworld)
homeworld <- starwars %>% select(name, homeworld)
#Filtramos por especies: unicamente en humanos
human <- starwars %>% filter(species == "Human")
#Filtramos por especies: unicamente humanos del planeta Tatooine
starwars %>% filter(species == "Human", homeworld == "Tatooine")
#Creamos un nuevo datframe con todas las especies exceptuando a los Droides
starwars_nodroids <- starwars %>% filter(species != "Droid")
#Utilizamos group_by y tally
starwars %>% group_by(species) %>% tally()
#Añadimos otra variable
starwars %>% group_by(species, gender) %>% tally()
#Si queremos guardarlo hay que ponerle un nombre
table_gender <- starwars %>% group_by(species, gender) %>% tally()
#Con el comando na.rm borramos los valores no asignados
starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T),mean_mass = mean(mass,na.rm = T))
#Hacemos un gráfico de altura y masa de los personajes
ggplot(starwars, aes(height, mass)) + geom_point()
#Cambiamos el color 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red")
#Cambiamos el color y el punto
ggplot(starwars, aes(height, mass)) + geom_point(colour = "purple", pch = 3)
#Cambiamos el fondo y el color
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red") + theme_light()
# Creamos un database donde vamos a filtrar por IMC y excluiremos nuestra persona que sobresale
filtradomasa <- starwars %>% filter(name != "Jabba Desilijic Tiure")
#Hacemos un gráfico de altura y masa de los personajes sin el que difiere
ggplot(filtradomasa, aes(height, mass)) + geom_point()
# Con el comando read podemos leer el archivo que hemos cargado
toy <- read_csv("toy.csv")
# Usamos summarise para encontrar todos los parametros estadisticos deseados
toy %>% group_by(Sex) %>% summarise(mean_height = mean(Height_cm, na.rm = T),mean_mass = mean(Weight_Kg, na.rm = T), mean_imc = mean(IMC, na.rm = T), mean_ias = mean(IAS, na.rm = T), mean_ccintura = mean(Ccintura, na.rm = T))
#Filtramos por sexo: unicamente mujeres y además que no sean de IMC normal
toy %>% filter(Sex == "Women", IMC_clas != "Normal")
#Filtramos por sexo: unicamente mujeres y además que no sean de IMC normal ni tengan obesidad
toy %>% filter(Sex == "Women", IMC_clas != "Normal", IMC_clas != "Obesity")
#Hacemos un gráfico de IMC y peso de las personas
ggplot(toy, aes(IMC, Weight_Kg)) + geom_point()
#Hacemos otro database sin las personas de IMC normal
toyno <- toy %>% filter(IMC_clas != "Normal")
#Hacemos un gráfico de IMC no normal y Peso
ggplot(toyno, aes(IMC, Weight_Kg)) + geom_point()

#Con estos comandos descargamos los paquetes necesarios para la siguiente sesión
install.packages("ape")
install.packages("phangorn")
install.packages("phytools")