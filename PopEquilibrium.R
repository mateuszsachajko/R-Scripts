chiSq <- function(obs, exp) 
{
  sum((obs-exp)**2/exp)
}


n = 1000
wielkPop = 1000
popChi = c()

for(i in 1:n){
  
  pop = sample(c("Aa","AA","aa"),
               prob = c(0.5,0.25,0.25),
               size = wielkPop, replace = T)
  table(pop)
  popChi = c(popChi, chiSq(c(table(pop)[1],table(pop)[2],table(pop)[3]),c(250,500,250)))   
}

plot(popChi, dchisq(popChi,2))
hist(popChi)
quantile(popChi,0.95)

popChi


plot(seq(0,15,by=0.1),dchisq(seq(0,15,by=0.1),2))
quantile(popChi,0.95)

popChi= c()

for(i in 1:n){
  wielkPop = sample(100:10000,1)
  Ho = wielkPop *0.25
  Ht = wielkPop *0.5
  pop = sample(c("Aa","AA","aa"),
               prob = c(0.5,0.25,0.25),
               size = wielkPop, replace = T)
  popChi = c(popChi, chiSq(c(table(pop)[1],table(pop)[2],table(pop)[3]),c(Ho,Ht,Ho)))   
}
plot(popChi, dchisq(popChi,2))
hist(popChi)
max(popChi)



######## migracje

wielkPop = 10000
populacja = sample(c("Aa","AA","aa"),
                   prob = c(0.5,0.25,0.25),
                   size = wielkPop, replace = T)
table(populacja)
chiSq(c(2501,4981,2518),c(2500,5000,2500))
populacja = c(populacja, rep("aa",3000))
table(populacja)
p = (2518*2 + 4981)/(2*13000)
q = (5501*2 + 4981)/(2*13000)

p + q
p
q

p2 = p*p
pq = 2*p*q
q2 = q*q
p2 + pq + q2

AA = p2*13000
Aa = pq * 13000
aa = q2 * 13000

chiSq(c(5501,4981,2518),c(aa,Aa,AA))

populacja = sample(c("Aa","AA","aa"),
                   prob = c(0.5,0.25,0.25),
                   size = wielkPop, replace = T)

populacja = populacja[-sample(1:length(populacja),5000, replace = F)]
table(populacja)
p = (2*sum(populacja=='AA') + sum(populacja=='Aa'))/(2*length(populacja))
q = (2*sum(populacja=='aa') + sum(populacja=='Aa'))/(2*length(populacja))
p+q

p2 = p*p
pq = 2*p*q
q2 = q*q
p2 + pq + q2

AA = p2 * length(populacja)
Aa = pq * length(populacja)
aa = q2 * length(populacja)
table(populacja)
chiSq(c(1287,2470,1243),c(aa,Aa,AA))


#####Selekcja
wielkPop = 1000
populacja = sample(c("Aa","AA","aa"),
                   prob = c(0.5,0.25,0.25),
                   size = wielkPop, replace = T)
populacja = sort(populacja)
populacja

#zrobic parzyste genotypy!!!!!
table(populacja)

chiSq(c(242,523,234),c(250,500,250))

populacjaF1 = c(rep("AA",table(populacja)[3]),rep("aa",table(populacja)[1]))
populacjaF1 = c(populacjaF1, sample(c("AA","Aa","aa"),size = 491,prob = c(0.25,0.5,0.25),replace = T))
table(populacjaF1)

populacjaF2 = c(rep("AA",table(populacjaF1)[3]),rep("aa",table(populacjaF1)[1]))
populacjaF2 = c(populacjaF2, sample(c("AA","Aa","aa"),size = table(populacjaF1)[2],prob = c(0.25,0.5,0.25),replace = T))
table(populacjaF2)
