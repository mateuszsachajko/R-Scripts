# zad.1
#losowanie danych


indeksy = 1:38
odpowiedzi = data.frame(indeksy)
dane = data.frame(indeksy)
dane[,2]=sample(1:280,size = 38)


dane
#zadanie1 odpowiedz

zad1 = function(x){
  rogate = 300-x
  q = round(sqrt(rogate/300),4)
  p = 1-q
  p2 = round(p**2,3)
  q2 = round(q**2,3)
  pq = round(2*p*q,3)
  odp = c(p,q,p2,pq,q2,pq*300)
}
zad1_x=dane[,2]
odp1 = lapply(as.numeric(zad1_x),zad1)
for(i in 1:38){
  print(odp1[[i]][1]+odp1[[i]][2])
  print(odp1[[i]][3]+odp1[[i]][4]+odp1[[i]][5])
}

for(i in 38){
  
  
}

# zadanie 2 dane
for(i in 1:38){
  
  pop = sample(c("Aa","AA","aa"),
               prob = c(0.5,0.25,0.25),
               size = 2500, replace = T)
  dane[i,3]=table(pop)[1]
  dane[i,4]=table(pop)[2]
  dane[i,5]=table(pop)[3]
}

chiSq <- function(obs, exp) 
{
  sum((obs-exp)**2/exp)
}
#[,3][,4][,5]
zad2 = function(hoR,het,hoD){
    p = round(sqrt(hoD/2500),4);
    q = 1-p
    p2 = p**2
    pq = 2*p*q
    q2 = q**2
    chi = chiSq(c(p2*2500,pq*2500,q2*2500),c(hoD,het,hoR))
    if(chi<=3.8415) wniosek = "jest w rownowadze"
    else wniosek = "nie jest w rownowadze"
    odp = c(chi,wniosek)
}
a = zad2(dane[1,3],dane[1,4],dane[1,5])
a
###zad 3 krew

for(i in 1:38){
  krew = sample(100:1000,size=4,replace = T)
  dane[i,6]=krew[1]
  dane[i,7]=krew[2]
  dane[i,8]=krew[3]
  dane[i,9]=krew[4]
}


zad3 = function(a,b,ab,o){
  pop = a+b+ab+o
  r = round(sqrt(o/pop),4)
  p = round(sqrt((o/pop)+(a/pop)),4) - r
  q = 1 - r - p
  bb = q*q
  bo = 2*q*r
  odp = c(p,q,r,round(bb,4),round(bo,4))
}

###zad 4 niezalezne

for(i in 1:38){
  p = round(runif(n = 1),2)
  q = 1-p
  r = round(runif(1),2)
  s = 1-r
  dane[i,10] = p
  dane[i,11] = q
  dane[i,12] = r
  dane[i,13] = s
}

zad4 = function(p,q,r,s){
  dzikie = (p*p+2*p*q)*(r*r+2*r*s)*20000
  cynamonowe = (p*p+2*p*q)*s*s*20000
  czarneJedno = q*q*(r*r+2*r*s)*20000
  czeko = q*q*s*s*20000
  odp = c(dzikie, cynamonowe, czarneJedno, czeko)
}

### zad 5 inbred
for(i in 1:38){
  inbred = round(runif(min = 0.05,max = 0.4,n = 1),4)
  dane[i,14] = inbred
}

zad5 = function(x){
  x/2
}

zad6 = function(p,q){
  kuraSre = p
  kuraZlo = q
  kogSre = p*p+2*p*q
  kogZlo = q*q
  odp = c(kuraSre, kuraZlo, kogSre,kogZlo)
}



### zad 6 sprzezone z plcia 

for(i in 1:38){
  p = round(runif(min = 0.05,max = 0.4,n = 1),3)
  q = 1 - p
  dane[i,15] = p
  dane[i,16] = q
}

names(dane)=c("indeksy","Zad1 X","Zad2 X","Zad2 Y","Zad2 Z", "Zad3 A","Zad3 B",
              "Zad3 AB","Zad3 0", "Zad4 p", "Zad4 q", "Zad4 r", "zad4 s","Zad5 X",
              "Zad 6 X", "Zad 6 Y")


write.csv(dane, file="dane_kolokwium2.csv")

#dane = read.csv("dane_kolokwium2.csv")

odpowiedzi = data.frame(indeksy)
dane
dane = dane[-1,]
colnames(dane) = dane[1,]
colnames(dane)
a = zad1(x = as.numeric(dane[1,2]))
lapply(as.numeric(dane[,2]),zad1)

dane1=dane[order(dane$indeksy),]
dane = dane1
dane = dane1[,-1]
dane = dane[,-1]
as.numeric(dane[,2])
for(i in 1:38){
  a = zad1(as.numeric(dane[i,2]))
  odpowiedzi[i,2] = a[1]
  odpowiedzi[i,3] = a[2]
  odpowiedzi[i,4] = a[3]
  odpowiedzi[i,5] = a[4]
  odpowiedzi[i,6] = a[5]
  odpowiedzi[i,7] = a[6]
  
}

for(i in 1:38){
 a = zad2(hoR = as.numeric(dane[i,3]),het = as.numeric(dane[i,4]),hoD = as.numeric(dane[i,5]))
 odpowiedzi[i,8] = a[1]
 odpowiedzi[i,9] = a[2]
}

for(i in 1:38){
  a = zad3(a = as.numeric(dane[i,6]),b = as.numeric(dane[i,7]),
           ab = as.numeric(dane[i,8]),o = as.numeric(dane[i,9]))
  odpowiedzi[i,10] = a[1]
  odpowiedzi[i,11] = a[2]
  odpowiedzi[i,12] = a[3]
  odpowiedzi[i,13] = a[4]
  odpowiedzi[i,14] = a[5]
  
}

for(i in 1:38){
  a = zad4(p = as.numeric(dane[i,10]),q = as.numeric(dane[i,11]),
           r = as.numeric(dane[i,12]),s = as.numeric(dane[i,13]))
  odpowiedzi[i,15] = a[1]
  odpowiedzi[i,16] = a[2]
  odpowiedzi[i,17] = a[3]
  odpowiedzi[i,18] = a[4]
}

for(i in 1:38){
  a = zad5(as.numeric(dane[i,14]))
  odpowiedzi[i,19] = a
}

for(i in 1:38){
  a = zad6(p = as.numeric(dane[i,15]),q = as.numeric(dane[i,16]))
  odpowiedzi[i,20] = a[1]
  odpowiedzi[i,21] = a[2]
  odpowiedzi[i,22] = a[3]
  odpowiedzi[i,23] = a[4]
}

write.csv(odpowiedzi,'odpowiedzi2.csv')
getwd()
