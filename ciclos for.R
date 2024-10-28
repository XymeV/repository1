#ciclos

for (i in 1:5) {
  print(i)
}

for(i in 1:1001){
  print("Hola Xymena")
}
for(i in 1:1001){
  print("Hola Alondra y Xymena :)")
}

#2
vector <- c(3, 5, 7, 9)
suma <- 0
for (i in vector) {
  suma <- suma + i
}
print(suma)

#con 3
vector <- c(1:100)
suma<-0
for (i in 1:length(vector)) {
 suma<-suma + i
}

print(suma)
print(vector)

#ejercicio data frame

fila <-data.frame(
  nombre =c("Melissa", "Iván", "Xymena", "Alondra"),
  edad = c(20, 20, 20, 19),
  color = c("morado", "azul", "lila", "naranja")
)
 print(fila) 
 for(i in 1:4){
   print(paste("El color favorito de ", fila[i, 1], "es el", fila[ i, 3], 
               "y tiene", fila[ i, 2], "años"))
 }
 
 