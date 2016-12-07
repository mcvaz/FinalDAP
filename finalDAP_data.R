#### FINAL DAP - MARCEL VAZ

### DATA

## Tree architectural traits
tt = read.csv("MAO16DATA_all.csv",header=T,as.is=T)[,c("plantID","genus","epithet",
  "code","cii","lblade","phylotax","phylotax2","btype","d10","dbh",
  "sh","cmean","cd1","cd2","tv","tdm","il1","il2","il3","il4","il5",
  "tl1","tl2","tl3","tl4","tl5","dbh0")]
head(tt)

tt$sp = paste(tt$genus,tt$epithet)

tt$sVol = pi*(tt$d10/2000)^2*tt$sh/3 # stem volume (cm^3)
tt$sDens = (tt$tdm/100)/(tt$tv/100) # stem density (g/cm^3)
tt$sMass = tt$sVol*tt$sDens # stem mass (g)

tt$cArea = pi*tt$cd1*tt$cd2/40000 # crown projected area (m^2)

tt$cPos = tt$cmean/100 # crown position (m)

tt$medInL = apply(tt[,c("il1","il2","il3","il4","il5")],1,median,na.rm=T)
tt$medTrL = apply(tt[,c("tl1","tl2","tl3","tl4","tl5")],1,median,na.rm=T)
tt$cAbs = (tt$medInL-tt$medTrL)/tt$medInL
  
tt$hge = tt$cPos/(tt$sMass/1000) # height gain efficiency (m/kg)
tt$lie = tt$cArea*tt$cAbs/(tt$sMass/1000)

tt = tt[!is.na(tt$hge),]
#tt2 = tt2[tt2$cAbs>0 & !is.na(tt2$lie),]
#tt2[tt2$cAbs==0,]

#plot(log(hge)~log(lie),tt2)
#abline(lm(log(hge)~log(lie),tt2),col=2)
#summary(lm(log(hge)~log(lie),tt2))
#plot(lm(log(hge)~log(lie),tt2))

hist(tt$hge)
hist(log(tt$hge),xlab='log(HGE)',main="")
any(tt$hge<=0)


## Adult size for all tree species
aSize = read.csv("maxsizes.csv",header=T,as.is=T)
aSize$sp = paste(aSize$Genus,aSize$epithet)
aSize$aSize = apply(aSize[,c("q975","q975F")],1,max,na.rm=T)
head(aSize)


unique(sort(tt$sp[!tt$sp%in%aSize$sp]))
tt = tt[tt$sp%in%aSize$sp,c("sp","hge")]
tt$aSize = NA
for(i in 1:nrow(tt)){
  tt[i,"aSize"] = aSize[aSize$sp==tt[i,"sp"],"aSize"]
  }

tt = tt[order(tt$sp),]
head(tt)

write.csv(tt,"HGEvsAS.csv",row.names=F)

length(unique(tt$sp)) # number of species





