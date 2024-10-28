####BIOSTRINGS####
IUPAC_CODE_MAP #funcion del alfabeto extendido.
dna1 <- DNAString("ATACGGT-N")
dna1
dna2 <- DNAStringSet(c("ACGT", "GTCA", "GCTA"))
dna2

dna2[[2]] #para seleccionar especificamente dentro de una secuencia.
#para seleccionar dentro de la seceuncia especifica: dos letras especificas
dna2[[3]][1:3]
dna2[1][1:2]#esta ya no funciona como las anteriores.
names(dna2) <- paste0("secuencia_", 1:3)
dna2
#ya con nombres puedo seleccionar la secuencias especifica que quiero al ubicarla.

#DNAString[Set]: Secuencias de DNA.
#RNAString[Set]: Secuencias de RNA.
#AAString[Set]: Secuencias de Amino Acids (proteínas).
#BString[Set]: Big sequences, cualquier letra.

###FUNCIONES BASICAS###
width(dna2) #te da el tamaño de las secuencias dentro del vector(cuantos nucelotidos tiene)
sort(dna2) #ordena alfabeticamente las secuencias.
rev(dna2) #reverso de la secuencias.
rev(dna1)

dna1
#reverse: reverso de la secuencia.
#complement, reverseComplement: (reverso) complemento de la secuencia.
dna1 #seq: N-ACGT
reverseComplement(dna1) #seq: ACGT-N
dna3 <- DNAString("CGAACG")
dna3
reverseComplement(dna3) #como los primers de Fausto
#translate: traduce DNA o RNA sequencia a aminoacidos.
#alphabetFrequency, letterFrequency: Calcula la frecuencia de todos los caracteres (alphabetFrequency) 
#   solo de secuencias especificas (letterFrequency).
alphabetFrequency(dna1)
letterFrequency(dna2, "GC") #puedo especificar que letras quiero que cuente 

####EJERCICIOS PRACTICA####
dna4 <- DNAString("TAGGAAATCGGTCAGTAGCTATTCGAAG")
dna4
dna5 <-DNAString("TAATAGCTCCTACATCACTTG")

#CODIGO GENETICO.
GENETIC_CODE

#LA FUNCION TRANSLATE IGNORA LOS CODONDES DE INICIO.

translate(dna4)#funcion para traducir las secuencias de ADN -RNA
translate(dna5)#para checar usa la funcion del codigo genetico.

reverse(dna4) #funcion que ordena la secuencia: 5´-3´y la pone 3´-5´
complement(dna4) #saca la secuencia complementaria de la cadena en ORDEN.
reverseComplement(dna4) #complementaria de 3´-5´
alphabetFrequency(dna4) #el total de las letras que estan, incluye el alfabeto extendido
letterFrequency(dna4, "GT") #hace la suma de las dos letras
#Coincidencia exacta:
matchPattern("CAG", dna4) #ENCUENTRA NUCLEOTIDOS ESPECIFICOS EN LA SECUENCIA
countPattern("CAG", dna4) #EL NUMEOR DE MATCHES QUE SI HAY

#Coincidencia aproximada (con mismatches):
matchPattern("TCG", dna4, max.mismatch=1)
countPattern("TCG", dna4, max.mismatch = 2)

#ejercicio con RNA
rna1 <- RNAString("AUGAGCUGUCGAUGUAG")
reverse(rna1)
translate(rna1) #tambien pasa las secuencias de ARN

###ALINEAMIENTO DE SECUENCIAS.
BiocManager::install("pwalign")
#alineamiento global
align1 <- pairwiseAlignment(dna2, dna3, type="global") #te dicen de toda la secuencia cuantos coinciden 
#toda la secuencia de alinemaiento 
align1
#alineamiento local.
align2 <- pairwiseAlignment(dna4, dna5, type="local") #solo toma pedazos, la region en la que localmente coinciden 
#regiones que pueden ser alineadas dentro de la secuencia
align2

#Encontrar palíndromos:
findPalindromes(dna4)

dinucleotideFrequency(dna4) #la probabilidad de encontrar los nucleotidos
trinucleotideFrequency(dna)

#Lectura de archivos FASTA:
# EJEMPLO: fasta_sequences <- readDNAStringSet("~/Downloads/GCF_000005845.2_ASM584v2_genomic.fna")
#se usa la funcion readDNAStringSet: recuerda usar comillas y el tabulador busca la secuencia
#para buscar la secuencia, descargarla y 
secuencia_fasta <- readDNAStringSet("GCF_000005845.2_ASM584v2_genomic.fna")
rev (dna1)
reverse(dna1)


width(secuencia_fasta) #te da el tamaño de las secuencias dentro del vector(cuantos nucelotidos tiene)
rev(secuencia_fasta) #reverso de la secuencias.
reverse(secuencia_fasta)
dinucleotideFrequency(secuencia_fasta) 
trinucleotideFrequency(secuencia_fasta)
alphabetFrequency(secuencia_fasta) #el total de las letras que estan, incluye el alfabeto extendido
letterFrequency(secuencia_fasta, "AT")
letterFrequency(secuencia_fasta, "CG")
#NUMERO DE MATCHES 
vmatchPattern("AAC",secuencia_fasta) #ENCUENTRA NUCLEOTIDOS ESPECIFICOS EN LA SECUENCIA


#Coincidencia aproximada (con mismatches):
vmatchPattern("TCG", secuencia_fasta, max.mismatch=1)
vcountPattern("TCG", secuencia_fasta, max.mismatch = 2)

#Ejercicio 1: Creación y manipulación de secuencias
#Encuentra el complemento inverso de la secuencia.
reverseComplement(secuencia_fasta)
#Cuenta las ocurrencias del nucleótido “A”.
letterFrequency(secuencia_fasta, "A")
#Extrae la subsecuencia de la posición 3 a la 7.
# [[]] dentro de esta secuencia que selecciones []
secuencia_fasta[[1]][3:7] #para que seleccione la secuencia exacta, no funciona con los corchetes.
subseq(secuencia_fasta, start = 3, end = 7)

#Ejercicio 2: Coincidencia de patrones
#Encuentra todas las coincidencias exactas del patrón “AGC”.
vmatchPattern("AGC", secuencia_fasta)
vcountPattern("AGC",secuencia_fasta)
#Realiza coincidencias aproximadas del patrón “AGC” permitiendo 1 desajuste.
vmatchPattern("AGC", secuencia_fasta, max.mismatch = 1)

#Ejercicio 3: Alineación de secuencias
#Realiza una alineación global entre las secuencias “ACGT” y “AGCT”.
align3 <- pairwiseAlignment(DNAString("ACGT"), DNAString("AGCT"), type="global") 
align3
#Realiza una alineación local entre las secuencias “ACGT” y “CG”.
align4 <- pairwiseAlignment(DNAString("ACGT"), DNAString("CG"), type="local")
align4
#Escribe el resultado de la alineación en un archivo.
writeXStringSet(fasta_sequences, "output.fasta")
Lectura y escritura de alineaciones:
  
  #alignments <- readDNAMultipleAlignment("alignment.aln")
  #write.phylip(alignments, "output.phy")