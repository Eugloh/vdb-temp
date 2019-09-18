install.packages("ape")
library(ape)

#Exo1 
?rtree
for (i in 1:9) plot(rtree(20))
plot(rtree(1000))
plot(rtree(4))
plot(tree,show.tip.label=F,type="u")
#install.packages("ggtree")
#library(ggtree)

#Exo2
install.packages("leaflet")
require("leaflet")
m<-leaflet("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png")
m<-addTiles(m)
m<-addMarkers(m, lng=4.85, lat=45.75, label="Ici Lyon")
mP
mP<-addMarkers(m, lng= 2.35, lat=48.85, label="Ici Paris")

?leaflet

recharge <- function() {
  m <- leaflet()
  m<-addTiles(m,"https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png")
  return(m)
}
recharge()



#Exo3


?data.frame
nom <- c("Lyon","Paris")
lg <- c(4.85,2.35)
lat <- c(45.75,48.85)
Map <- data.frame(Ville = nom, Longitude = lg, Latitude = lat )


Map
m <- leaflet(data=Map)
m <- addTiles(m)
m <- addMarkers(m,lng=~Longitude,lat=~Latitude)
m



#http://lifemap-ncbi.univ-lyon1.fr:8080/solr/#/


m<-leaflet()
m <- addTiles(m," http://lifemap-ncbi.univ-lyon1.fr/osm_tiles/{z}/{x}/{y}.png")
m

recharge <- function(dataframe) {
  m <- leaflet(dataframe)
  m<-addTiles(m,"http://lifemap-ncbi.univ-lyon1.fr/osm_tiles/{z}/{x}/{y}.png")
  return(m)
}
install.packages("jsonlite")
library(jsonlite)
fromJSON("http://lifemap-ncbi.univ-lyon1.fr:8080/solr/taxo/select?q=taxid:(2%209443%202087)&wt=json&rows=1000")

Getcoord <- function(taxids) {
  taxids <- as.character(taxids)
  taxids2<-paste(taxids,collapse="%20")
  url <- paste("http://lifemap-ncbi.univ-lyon1.fr:8080/solr/taxo/select?q=taxid:(",taxids2,")&wt=json&rows=1000")
  info <- fromJSON(url)
  }
Getcoord(9443)

require(jsonlite)
GetCooFromTaxID<-function(taxids) {
  ##taxids is an array that contains taxids.
  ## url cannot be too long, so that we need to cut the taxids (100 max in one chunk)
  ## and make as many requests as necessary.
  taxids<-as.character(taxids) #change to characters.
  DATA<-NULL
  i<-1
  while(i<=length(taxids)) {
    cat(".")
    taxids_sub<-taxids[i:(i+99)]
    taxids_sub<-taxids_sub[!is.na(taxids_sub)]
    taxids_sub<-paste(taxids_sub, collapse="%20") #accepted space separator in url
    url<-paste("http://lifemap-ncbi.univ-lyon1.fr:8080/solr/taxo/select?q=taxid:(",taxids_sub,")&wt=json&rows=1000",sep="", collapse="")
    #do the request :
    data_sub<-fromJSON(url)
    DATA<-rbind(DATA,data_sub$response$docs[,c("taxid","lon","lat", "sci_name","zoom","nbdesc")])
    i<-i+100
  } 
  for (j in 1:ncol(DATA)) DATA[,j]<-unlist(DATA[,j])
  class(DATA$taxid)<-"character"
  return(DATA)
}

##test de la fonction
data<-GetCooFromTaxID(c(2,9443,2087))
data


##RÉCUPÉRER LES DONNÉES
EukGenomeInfo<-read.table("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt", sep="\t", header=T, quote="\"", comment.char="")
## liste unique des taxid
taxids<-unique(EukGenomeInfo$TaxID)

## RÉCUPÉRER LES COORDONNÉES
DF<-GetCooFromTaxID(taxids)

## CALCULER LE NOMBRE DE GÉNOMES SÉQUENCÉS POUR CHAQUE TAXID
nbGenomeSequenced<-table(EukGenomeInfo$TaxID)
## l'ajouter à DF
DF$nbGenomeSequenced<-as.numeric(nbGenomeSequenced[DF$taxid])

##CALCULER LE NB DE GENOMES ENTIEREMENT ASSEMBLES POUR CHAQUE TAXID
##le calcul pour un seul taxid nommé 'tid' serait :
sum(EukGenomeInfo[which(EukGenomeInfo$TaxID=='tid'),]$Status=="Chromosome")
##on peut utiliser la fonction sapply pour le faire pour chaque taxid 
nbGenomeAssembled<-sapply(DF$taxid, function(x,tab) sum(tab[which(tab$TaxID==x),]$Status=="Chromosome"), tab=EukGenomeInfo)
DF$nbGenomeAssembled<-nbGenomeAssembled

##CALCULER LE TAUX de GC MOYEN  
tauxgcmoyen<-sapply(DF$taxid, function(x,tab) mean(as.numeric(as.character(tab[which(tab$TaxID==x),]$GC.)), na.rm=TRUE), tab=EukGenomeInfo)
DF$tauxgcmoyen<-tauxgcmoyen

##CALCULER LA TAILLE MOYENNE DES GÉNOMES EN Mb
SizeGenomeMb<-sapply(DF$taxid, function(x,tab) mean(tab[which(tab$TaxID==x),]$Size..Mb., na.rm=TRUE), tab=EukGenomeInfo)
DF$SizeGenomeMb<-SizeGenomeMb 

head(DF)


#eXO 6
addCircleMarkers(DF$nbGenomeSequenced)

recharge <- function(data=NULL) {
  m <- leaflet(data)
  m<-addTiles(m, url ="http://lifemap-ncbi.univ-lyon1.fr/osm_tiles/{z}/{x}/{y}.png",options = tileOptions(maxZomm=42))
  return(m)
}
m <- recharge(DF)
pal <- colorNumeric("Spectral",DF$tauxgcmoyen)
m <- addCircleMarkers(m , lat=~lat, lng=~lon, label =~ sci_name,
                        radius=~sqrt(nbGenomeAssembled),
                      stroke = FALSE,
                      color= "red",
                      fillOpacity = 0.5,
                      fillColor= ~pal(tauxgcmoyen))
m
m <- addLegend(m,position = "bottomright", pal=pal, 
               values =~tauxgcmoyen, title = "%GC moyen des génomes séquencés")

#options = tileOptions(maxZomm=42))

install.packages("RColorBrewer")
library(RColorBrewer)

