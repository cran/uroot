
# CONSOLE FUNCTIONS

ret <- function(vari, k)
{
 N <- length(vari)
 vari.r <- matrix(c(1:k*N),N,k)
 i <- 1
 while(i < k){
  vari.r[1:i,(i+1)] <- NA
  vari.r[,1] <- vari
  vari.r[(i+1):N, (i+1)] <- vari[1:(N-i)]
  i <- i+1
             }
 vret <- vari.r
}

#

BoxCox <- function(vari, lambda)
{
  if(lambda==0)
    bcvari <- log(vari, base=exp(1))
  if(lambda!=0)
    bcvari <- (vari^lambda - 1)/lambda
  bcvari
}

#

vfic <- function(y0fic, s0fic, yNfic, sNfic)
{
   aux     <- ysooys(1, .wvari$t0, .wvari$N, .wvari$s)[[2]][,1]
   years   <- seq(aux[1], aux[length(aux)])
   Mfic    <- rep(0, .wvari$N)
   obsfic0 <- (which(years==y0fic)-1) * .wvari$s + s0fic
   obsficN <- (which(years==yNfic)-1) * .wvari$s + sNfic
   Mfic[obsfic0:obsficN] <- 1
   matrix(Mfic, ncol=1)
}

#

contts <- function(lm, a)
{
  X <- model.matrix(lm)
  XX <- solve(t(X)%*%X)
  uhat <- as.matrix(lm$residuals)
  dff <- df.residual(lm)
  var.u <- t(uhat)%*%uhat/dff
  var.coef <- var.u*XX[a,a]
  se.coef <- sqrt(var.coef)
  et <- lm$coef[a]/se.coef
  list(se.coef=se.coef, t.stat=et)
}

#

regr.frec <- function(vari, s, t0){

 N  <- length(vari)
 ML <- ret(vari, s+1)

 if(s == 12)
 {
   Fil.vari <- filtrar(vari, s, t0, c(0,1,1,1,1,1,1), plot=FALSE)[[1]]
   y1 <- ret(Fil.vari, 2)[,2]

   Fil.vari <- filtrar(vari, s, t0, c(1,0,1,1,1,1,1), plot=FALSE)[[1]]
   y2 <- ret(-Fil.vari, 2)[,2]

   Fil.vari <- filtrar(vari, s, t0, c(1,1,0,1,1,1,1), plot=FALSE)[[1]]
   y3 <- ret(-Fil.vari, 3)[,3]
   y4 <- ret(-Fil.vari, 2)[,2]

   Fil.vari <- filtrar(vari, s, t0, c(1,1,1,0,1,1,1), plot=FALSE)
   y6 <- ret((sqrt(3)/2)*Fil.vari[[1]], 2)[,2]

   fcoefrfil <- crpp(Fil.vari[[2]],Fil.vari[[3]],c(1,2),c(0,1))
   fcoef <- fcoefrfil[[1]]
   rfil  <- fcoefrfil[[2]]
   y5.aux  <- matrix(nrow=N, ncol=length(rfil))
   y5.aux2 <- matrix(nrow=N, ncol=1)
   for(i in 1:length(rfil))
     y5.aux[,i] <- fcoef[i]*ML[,(rfil[i]+1)]
   for(i in 1:N)
     y5.aux2[i,1] <- sum(y5.aux[i,])
   y5 <- -ret((1/2)*(y5.aux2), 2)[,2]

   Fil.vari <- filtrar(vari, s, t0, c(1,1,1,1,0,1,1), plot=FALSE)
   y8 <- ret(-(sqrt(3)/2)*Fil.vari[[1]], 2)[,2]

   fcoefrfil <- crpp(Fil.vari[[2]],Fil.vari[[3]],c(1,-2),c(0,1))
   fcoef <- fcoefrfil[[1]]
   rfil  <- fcoefrfil[[2]]
   y7.aux  <- matrix(nrow=N, ncol=length(rfil))
   y7.aux2 <- matrix(nrow=N, ncol=1)
   for(i in 1:length(rfil))
     y7.aux[,i] <- fcoef[i]*ML[,(rfil[i]+1)]
   for(i in 1:N)
     y7.aux2[i,1] <- sum(y7.aux[i,])
   y7 <- ret((1/2)*y7.aux2, 2)[,2]

   Fil.vari <- filtrar(vari, s, t0, c(1,1,1,1,1,0,1), plot=FALSE)
   fcoefrfil <- crpp(Fil.vari[[2]],Fil.vari[[3]],c(sqrt(3),2),c(0,1))
   fcoef <- fcoefrfil[[1]]
   rfil  <- fcoefrfil[[2]]
   y9.aux  <- matrix(nrow=N, ncol=length(rfil))
   y9.aux2 <- matrix(nrow=N, ncol=1)
   for(i in 1:length(rfil))
     y9.aux[,i] <- fcoef[i]*ML[,(rfil[i]+1)]
   for(i in 1:N)
     y9.aux2[i,1] <- sum(y9.aux[i,])
   y9 <- -ret((1/2)*(y9.aux2), 2)[,2]

   Fil.vari <- filtrar(vari, s, t0, c(1,1,1,1,1,0,1), plot=FALSE)[[1]]
   y10 <- ret((1/2)*Fil.vari, 2)[,2]

   Fil.vari <- filtrar(vari, s, t0, c(1,1,1,1,1,1,0), plot=FALSE)
   fcoefrfil <- crpp(Fil.vari[[2]],Fil.vari[[3]],c(sqrt(3),-2),c(0,1))
   fcoef <- fcoefrfil[[1]]
   rfil  <- fcoefrfil[[2]]
   y11.aux  <- matrix(nrow=N, ncol=length(rfil))
   y11.aux2 <- matrix(nrow=N, ncol=1)
   for(i in 1:length(rfil))
     y11.aux[,i] <- fcoef[i]*ML[,(rfil[i]+1)]
   for(i in 1:N)
     y11.aux2[i,1] <- sum(y11.aux[i,])
   y11 <- ret((1/2)*(y11.aux2), 2)[,2]

   Fil.vari <- filtrar(vari, s, t0, c(1,1,1,1,1,1,0), plot=FALSE)[[1]]
   y12 <- ret(-(1/2)*Fil.vari, 2)[,2]

   Myr <- as.matrix(data.frame(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12))
 }
#
 if(s == 4)
 {
   Myr    <- matrix(nrow=N, ncol=s)
   filtra <- c(1,1,1)
   caux   <- c(1,-1,-1,-1)
   retaux <- c(1,1,2,1)

   for(j in 1:s){
     ifelse(j < 4, filtra[j] <- 0, filtra[3] <- 0)
     Fil.vari <- filtrar(vari, s, t0, filtra, plot=FALSE)[[1]]
     Myr[,j] <- ret(caux[j]*Fil.vari, retaux[j]+1)[,retaux[j]+1]
     filtra <- c(1,1,1)
                }
 }
 Myr
}

#

Fsnd <- function(lm2,nr, ncoef, ccc)
{
  # coeficientes a contrastar conjuntamente, ej: c(3, 8)
 R <- matrix(0, nrow=nr ,ncol=ncoef, ccc)
 for(i in 1:length(ccc))
  R[i, ccc[i]] <- 1

 r <- matrix(nrow=nr, ncol=1)
 r <- rep(0, nr)

 lm1 <- lm2
 best <- matrix(nrow=ncoef, ncol=1)
 for(i in 1:ncoef)
   best[i,1] <- lm1$coef[i]

 Rbestr <- R%*%best - r
 X <- model.matrix(lm1)
 XX <- solve(t(X)%*%X)
 uhat <- lm1$residuals
 dff <- lm1$df
 RXXR <- solve(R%*%XX%*%t(R))

 Fs <- (t(Rbestr)%*%RXXR%*%Rbestr/nr)/(sum(uhat^2)/dff)
 Fs
}

crpp <- function(vcoef1, vrfil1, vcoef2, vrfil2)
{
  n1 <- length(vcoef1)
  n2 <- length(vcoef2)
  fcoef <- rep(NA, n1*n2)
  rfil  <- rep(NA, n1*n2)
  fcoef0 <- rep(0, n1*n2)
  rfil0  <- rep(0, n1*n2)
  k <- 1
  for(i in 1:n1){
    for(j in 1:n2){
       fcoef0[k] <- vcoef1[i]*vcoef2[j]
       rfil0[k]  <- vrfil1[i]+vrfil2[j]
       k <- k+1
                  }
                }
  # Simplificar
  for(i in 1:(n1*n2)){
     simpl <- which(rfil0 == rfil0[i])
     if(length(simpl)>0){
        fcoef[i] <- sum(fcoef0[simpl]); fcoef0[simpl] <- NA
        rfil[i]  <- rfil0[i]  ; rfil0[simpl]  <- NA
                        }
                   }
  fcoef[which(fcoef==0)] <- rfil[which(fcoef==0)] <- NA
  fcoef <- na.omit(fcoef)[1:length(na.omit(fcoef))]
  rfil  <- na.omit(rfil)[1:length(na.omit(rfil))]

  # Ordenar
  fcoef0 <- rfil0 <- c(1:length(rfil))
  for(j in 1:length(rfil)){
    wm        <- which.min(rfil)
    rfil0[j]  <- rfil[wm]
    fcoef0[j] <- fcoef[wm]
    rfil[wm]  <- fcoef[wm] <- Inf
                          }
  rfil <- rfil0; fcoef <- fcoef0
  list(fcoef, rfil)
}

#

# filtra <- c(1,1,0,0,0,0,0)  # input
# Funciona con series trimestrales siempre que se tenga
  # cuidado de poner el vector filtra sin ceros a partir de
  # la tercera posición.

filtrar <- function(vari, s, t0, filtra, plot)
{
 N     <- length(vari)
 ML    <- ret(vari, s+1) # para mensuales basta con ret(vari,s)
 VCOEF <- cbind(c(1,-1,0),c(1,1,0),c(1,1,0),c(1,1,1),c(1,-1,1), c(1,sqrt(3),1),c(1,-sqrt(3),1))
 VRFIL <- cbind(c(0,1,0),c(0,1,0),c(0,2,0),c(0,1,2), c(0,1,2),c(0,1,2),c(0,1,2))

 if(length(which(filtra == 1)) == 1)
 {
   Rfil  <- c(1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2)
   Fcoef <- c(-1,1,1,1,1,-1,1,sqrt(3),1,-sqrt(3),1)

   frfil <- which(filtra == 1)
   Fil.vari <-  ML[,1] + Fcoef[frfil]*ML[,Rfil[frfil]+1]

   auxt0 <- length(vari)-length(na.omit(Fil.vari))
   Fil.vari <- ts(Fil.vari, frequency = s, start=t0)
   fcoef <- Fcoef[frfil]
   rfil <- Rfil[frfil]

   if(plot == TRUE)
     plot.ts(Fil.vari, col = "orange", las=1, xlab="", ylab="", main="Filtered series")
 }

 if(length(which(filtra == 1)) > 1)
 {
    frfil <- which(filtra == 1)
    fcoefrfil <- crpp(VCOEF[,frfil[1]], VRFIL[,frfil[1]], VCOEF[,frfil[2]], VRFIL[,frfil[2]])
    fcoef <- fcoefrfil[[1]]
    rfil  <- fcoefrfil[[2]]
    if(length(which(filtra == 1)) > 2){
      for(i in 3:length(frfil)){
        fcoefrfil <- crpp(fcoef, rfil, VCOEF[,frfil[i]], VRFIL[,frfil[i]])
        fcoef <- fcoefrfil[[1]]
        rfil  <- fcoefrfil[[2]]
                               }
                                      }
    Fil.vari.aux  <- matrix(nrow=N, ncol=length(rfil))
    Fil.vari      <- matrix(nrow=N, ncol=1)

    for(i in 1:length(rfil))
       Fil.vari.aux[,i] <- fcoef[i]*ML[,(rfil[i]+1)]
    for(i in 1:N)
       Fil.vari[i,1] <- sum(Fil.vari.aux[i,])
    Fil.vari <- ts(Fil.vari, frequency = s, start=t0)

   if(plot == TRUE)
     plot.ts(Fil.vari, col = "orange", las=1, xlab="", ylab="", main="Filtered series")
 }
 list(Fil.vari, fcoef, rfil)
}

MVFE <- function(vari, s, t0, tipo)
{
  N    <- length(vari)

  if(tipo == "alg"){        # Empleada en Barsky & Miron (1989)

   auxD <- matrix(0,nrow=N, ncol=s)
   sq   <- seq(1,N,s)
   k    <- 0

   for(j in 1:s){
     ifelse(sq[length(sq)] + k > N, n <- length(sq)-1, n <- N)
     for(i in 1:n)
       auxD[sq[i]+k,j] <- 1
     k <- k+1
   }
   VFE <- auxD
   if(t0[2] != 1){
      VFE <- matrix(nrow=N, ncol=s)
      VFE[,1:(t0[2]-1)] <- auxD[,(s-t0[2]+2):s]; VFE[,t0[2]:s] <- auxD[,1:(s-t0[2]+1)]
                  }
   if(t0[2] == 1){ VFE <- auxD }
  }

  if(tipo == "trg"){        # Empleada en Granger & Newbold (1986)

   # La notación, R1, viene de Canova & Hansen (1995)
     # las filas son f1, f2,... en ese artículo

   qq  <- s/2
   R1 <- matrix(nrow=N, ncol=(s-1))

   sq1 <- seq(1,qq*2,2)
   sq2 <- seq(2,qq*2,2)
   j   <- c(1:(qq-1))

   for(i in 1:N){
     for(k in 1:(s-qq-1)){
         R1[i,sq1[k]] <- cos((j[k]*pi/qq)*i)
         R1[i,sq2[k]] <- sin((j[k]*pi/qq)*i)
     }
     R1[i,(s-1)] <- (-1)^i
   }

   VFE <- R1
  }
  VFE
}

#
# usar selecP == "v1"

ADF.test <- function(label, compdet, selecP, Mvfic, VFEp, showcat)
{
  #vari.env <- new.env(FALSE, NULL)
  #Fvarinfo(label$label, env.out=vari.env)
  #vari <- get("vari", env=vari.env)
  #s    <- get("s", env=vari.env)
  #t0   <- get("t0", env=vari.env)
  #N    <- length(vari)
  vari <- label$vari
  s    <- label$s
  t0   <- label$t0
  N    <- length(vari)

  DGP  <- matrix(vari, nrow=N, ncol=1)
  EtFs <- matrix(nrow=1, ncol=(2+3*s)/2)
  Ct   <- matrix(rep(1, N), ncol=1)
  TD   <- matrix(c(1:N), ncol=1)
  VFE  <- MVFE(DGP[,1], s, t0, "alg")[,1:(s-1)]
  ifelse(Mvfic==0, nvf1<-0, nvf1<-ncol(Mvfic))
  ifelse(VFEp==0, nvf2<-0, nvf2<-ncol(VFEp))
  ldet <- compdet[1]*1 + compdet[2]*1 + compdet[3]*(s-1) + nvf1 + nvf2
  MLy  <- matrix(ret(DGP, 2)[,2], ncol=1)

# Dependent variable

  Deltay <- matrix(ret(DGP, 2)[,1]-ret(DGP, 2)[,2], ncol=1)

# Deterministic regressors

  Mregdet <- as.matrix(data.frame(Ct, TD, VFE, Mvfic, VFEp, MLy))
  aux <- rep(1, ncol(Mregdet))
  for(i in 1:ncol(Mregdet)){
    if(length(which(Mregdet[,i]==0))==N)
    aux[i] <- 0
                           }
  Mregdet <- Mregdet[,which(aux==1)]
  cdf1 <- which(compdet[1:2] == 1)
  if(compdet[3] == 1){ cdf2 <- c(3:(s+1)) }
  if(compdet[3] == 0){ cdf2 <- NULL }
  if(length(which(Mvfic==0))>1 || length(which(VFEp==0))>1)
    cdf3 <- c((s+2):(s+2+nvf1+nvf2-1))
  if(length(which(Mvfic==0))==1 && length(which(VFEp==0))==1)
    cdf3 <- NULL
  cdf4 <- c(s+2+nvf1+nvf2)
  cdf <- c(cdf1, cdf2, cdf3, cdf4)
  lmdet <- lm(Deltay[,1] ~ 0+Mregdet[,cdf])

# Adding lags of Deltay, using selecP function

 Pmax <- rsug <- -1
 if(mode(selecP) == "character"){
 if(selecP == "aic")
    Pmax <- selecPabic(Deltay, s, lmdet)[1]
 if(selecP == "bic")
    Pmax <- selecPabic(Deltay, s, lmdet)[2]
 if(selecP == "aiclb")
   Pmax <- selecPv5(Deltay, s, lmdet)
 if(selecP == "biclb")
   Pmax <- selecPv6(Deltay, s, lmdet)
 if(selecP == "aiclut")
    rsug <- selecPv2(Deltay, s, lmdet)
 if(selecP == "biclut")
    rsug <- selecPv3(Deltay, s, lmdet)
 if(selecP == "signf")
    rsug <- selecPv4(Deltay, s, lmdet)
                              }
 if(mode(selecP) == "numeric"){
  if(length(selecP)==1)
     Pmax <- selecP    # si se indica un orden elegido por uno mismo
  if(length(selecP)>1) # cuidado si se quiere sólo un retardo > 1
     rsug <- selecP
                              }
 if(Pmax > rsug[1] && Pmax > 0)
    Deltayr <- as.matrix(ret(Deltay, Pmax+1)[,2:(Pmax+1)])
 if(Pmax < rsug[1] && length(rsug) > 0)
    Deltayr <- as.matrix(ret(Deltay, round(10*log10(N))+2)[,(rsug+1)])
 ifelse(Pmax > 0 || length(which(rsug!=-1)) > 0,
        nP <- ncol(Deltayr), nP <- 0)

# Selection of the auxiliar regression components and testing

 if(Pmax > 0)
   lmadf <- lm(Deltay[,1] ~ 0+Mregdet[,cdf] + Deltayr[,1:nP])
 if(Pmax == 0)
   lmadf <- lm(Deltay[,1] ~ 0+Mregdet[,cdf])
 if(rsug[1] != -1 && length(rsug) > 0)
   lmadf <- lm(Deltay[,1] ~ 0+Mregdet[,cdf] + Deltayr[,1:nP])
 if(length(rsug) == 0)
   lmadf <- lm(Deltay[,1] ~ 0+Mregdet[,cdf])

# Statistics

 t.adf <- round(contts(lmadf, (ldet+1))[[2]], 2)
 #
 rdodet <- NULL
 if(ldet > 0)
 {
   rdodet <- matrix(nrow=ldet, ncol=2)
   for(i in 1:ldet)
   {
     rdodet[i,1] <- lmadf$coef[i]
     rdodet[i,2] <- contts(lmadf, i)[[2]]
   }
   Coeff   <- rdodet[,1]
   t.stat  <- rdodet[,2]

   if(compdet[3] == 1){
     if(s==4) {Fvfeaux <- c(2,3,4) + compdet[2] + nvf1+nvf2}
     if(s==12){Fvfeaux <- c(2,3,4,5,6,7,8,9,10,11,12) + compdet[2] + nvf1+nvf2}
     F.VFE  <- Fsnd(lmadf, s-1, 1+ldet+nP, Fvfeaux)
     rdodet <- round(data.frame(Coeff, t.stat, F.VFE), 2)
                      }
   if(compdet[3] == 0)
     rdodet <- round(data.frame(Coeff, t.stat), 2)
 }
  ifelse(rsug[1]==-1, retardos <- Pmax, retardos <- rsug)


# Summary

  if(showcat==TRUE){
    ifelse(compdet[1]==1, cd1 <- "Intercept", cd1 <- NA)
    ifelse(compdet[2]==1, cd2 <- "Trend", cd2 <- NA)
    ifelse(compdet[3]==1, cd3 <- c("Seasonal dummys", rep("", s-2)), cd3 <- NA)
    ifelse(nvf1!=0, cd4 <- c("Generic dummys", rep("", nvf1-1)), cd4 <- NA)
    ifelse(nvf2!=0, cd5 <- c("Partial seasonal dummys", rep("", nvf2-1)), cd5 <- NA)
    cdlabel <- na.omit(c(cd1, cd2, cd3, cd4, cd5))
    cdrdo   <- data.frame(Components=cdlabel, rdodet[,1:2])
      cat(c("\n ------ ADF test ------ \n"))
      cat(c("Statistic for the null hipothesis of \n unit root: ", t.adf, "\n\n "))
      cat("Deterministic components estimates: \n")
      print(cdrdo)
      if(compdet[3]==1)
      cat(c("Seasonal dummys F.test: ", rdodet[1,3], "\n"))
      cat(c("\n Number of lags selected:\n", retardos, "\n\n"))
      cat(c("Number of effective observations: ", (N-1-max(retardos)), "\n\n"))
#   if(exists("done")) tclvalue(done) <<- 1  # For MakeADF.test (wait before removing Mvfic, VFEp)
  }
  list(stad.t=t.adf, compdeter=rdodet, Lags=retardos, nobs=(N-1-max(retardos)))
}

#

# G R Á F I C O S  B U Y S  &  B A L L O T

# BUYS & BALLOT QUARTERLY PLOT

quarterg <- function(vari, s, t0, plot)
{
 N    <- length(vari)
 colour=c("SlateBlue","SeaGreen","red","magenta")
 leyenda <- c("Qrtr1","Qrtr2","Qrtr3","Qrtr4")
 #leyenda <- c("Trm1","Trm2","Trm3","Trm4")

 #naux <- length(vari)/4
 #ifelse(length(vari)/4 == as.integer(naux), n <- naux, n <- as.integer(naux)+1)
 n <- ysooys(N, t0, N, s)[[1]][1] - ysooys(1, t0, N, s)[[1]][1] + 1
 anyos <- c(1:n)
 for(i in 0:n)
 anyos[i+1] <- t0[1]+i

 MQrt <- matrix(nrow=n, ncol=4)
 r <- c(1:(n*4-t0[2]+1))
 raux1 <- rep(1,5-t0[2])
 raux2 <- gl(n, 4)
 r[1:length(raux1)] <- raux1
 r[(1+length(raux1)):length(r)] <- raux2[5:length(raux2)]

 i <- 1
 while(i < (length(vari)+1)){
  if(i < 5-t0[2])
   for(j in t0[2]:4){
    MQrt[r[i],j] <- vari[i]
    i <- i+1
   }
  else
   for(j in 1:4){
    MQrt[r[i],j] <- vari[i]
    i <- i+1
   }
 }
 Qrt1 <- ts(MQrt[,1], frequency=1 , start = t0[1])
 Qrt2 <- ts(MQrt[,2], frequency=1 , start = t0[1])
 Qrt3 <- ts(MQrt[,3], frequency=1 , start = t0[1])
 Qrt4 <- ts(MQrt[,4], frequency=1 , start = t0[1])
 MQ <- ts(data.frame(Qrt1,Qrt2,Qrt3,Qrt4), frequency=1, start=t0[1])

# opar <- par(mar=c(8,4.7,5,1.5), ps=18, font=1,tcl=-0.5, cex.axis=0.7, las=1)
 if(plot==TRUE){
 plot.ts(MQ, lty=c(1,2,3,4), plot.type="single",
#        main="Evolucion de cada trimestre", cex.main=0.9,
         col = c("SlateBlue","SeaGreen","red","magenta"),
                 #, "navy", "sienna"),
         xlim=c(t0[1],anyos[length(anyos)]+1),
         xlab="", ylab="")
 for(i in 1:4){
  ord <- length(na.omit(MQ[,i]))
  text(anyos[ord]+1.5, MQ[ord,i],
       leyenda[i], cex=0.7, col = colour[i]) #, font=2)
 }
# par(opar)
               }
 MQ
}

# BUYS-BALLOT MONTH PLOT
bbmp <- function(vari, s, t0, mp, vers, plot)
{
 N <- length(vari)

 #naux <- length(vari)/12
 #ifelse(length(vari)/12 == as.integer(naux), n <- naux, n <- as.integer(naux)+1)
 n <- ysooys(N, t0, N, s)[[1]][1] - ysooys(1, t0, N, s)[[1]][1] + 1
 anyos <- c(1:n)
 for(i in 0:n)
 anyos[i+1] <- t0[1]+i

 r <- c(1:(n*12-t0[2]+1))
 raux1 <- rep(1,13-t0[2])
 raux2 <- gl(n, 12)
 r[1:length(raux1)] <- raux1
 r[(1+length(raux1)):length(r)] <- raux2[13:length(raux2)]

 MMth <- matrix(nrow=n, ncol=12)
 MM <- matrix(nrow=n, ncol=12)
 MMp <- matrix(nrow = n, ncol = length(mp))
 leyenda <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

 i <- 1
 while(i < (length(vari)+1)){
  if(i < 13-t0[2])
   for(j in t0[2]:12){
    MMth[r[i],j] <- vari[i]
    i <- i+1
   }
  else
   for(j in 1:12){
    MMth[r[i],j] <- vari[i]
    i <- i+1
   }
 }
 j <- 1
 while(j < (length(mp)+1)){
  for(i in mp[j])
   MM[,i] <- ts(MMth[,i], frequency = 1, start = t0[1])
  j <- j+1
 }
 data.frame(MM[,1:12])

 i <- 1
 j <- 1
 # n = 6, aquí puede dar problemas si falta la n-ésima
   # observación del mes
 while(i < 13){
  ifelse(MM[6,i] != "NA", MMp[,j] <- MM[,i], MM[,i] <- MM[,i])
  ifelse(MM[6,i] != "NA", j <- j+1, j <- j)
  i <- i+1
 }

 MMp <- ts(MMp, frequency=1, start=t0[1])

 if(plot ==TRUE){
 if(vers == "R"){
  colour = rep(c("SlateBlue", "SeaGreen", "red", "magenta", "navy", "sienna"),2)
#  opar <- par(mar=c(8,4.7,5,1.5), ps=18, font=1, tcl=-0.5,
#              cex.axis=0.7, las=1)
  plot.ts(MMp, col = c("SlateBlue","SeaGreen","red","magenta", "navy", "sienna"),
          xlim = c(anyos[1],anyos[n]+1.3),
#         main = "Evolucion de cada mes", cex.main=0.9,
          lty = c(1,2,3,4), plot.type="single", las= 1,
          xlab="", ylab="", cex.axis=0.7)
  for(i in 1:length(mp)){
   ord <- length(na.omit(MMp[,i]))
   text(anyos[ord]+1.1, MMp[ord,i], leyenda[mp[i]], cex=0.7, col = colour[i])
                        }
#  par(opar)
                }
 if(vers == "Prot"){         # Tres meses por gráfico (panel)
  mp <- c(1,2,3,4,5,6,7,8,9,10,11,12)
  colour = rep(c("SlateBlue", "SeaGreen", "red"), 4) #  "magenta"
  ylimmin <- min(na.omit(MMp))
  ylimmax <- max(na.omit(MMp))
  rem <- (ylimmax - ylimmin)/5
#  opar <- par(mfrow=c(2,2), mar=c(3,3.9,1.5,0.5), ps=18, font=1,
#              tcl=-0.5, cex.axis=0.7, las=1)
  plot.ts(MMp[,1:3], col = c("SlateBlue","SeaGreen","red"),
          xlab="", ylab="", ylim=c(ylimmin, ylimmax+rem),
          xlim = c(anyos[1],anyos[n]+1.3),
#         main = "Evolucion de cada mes",
          lty = c(1,2,3,4), plot.type="single")
  for(i in 1:3){
   ord <- length(na.omit(MMp[,i]))
   text(anyos[ord]+1.1, MMp[ord,i],
        leyenda[mp[i]], col = colour[i], cex=0.7)
               }
  plot.ts(MMp[,4:6], col = c("SlateBlue","SeaGreen","red"),
          xlab="", ylab="", ylim=c(ylimmin, ylimmax+rem),
          xlim = c(anyos[1],anyos[n]+1.3),
#         main = "Evolucion de cada mes",
          lty = c(1,2,3,4), plot.type="single")
  for(i in 4:6){
   ord <- length(na.omit(MMp[,i]))
   text(anyos[ord]+1.1, MMp[ord,i],
        leyenda[mp[i]], col = colour[i], cex=0.7)
               }
  plot.ts(MMp[,7:9], col = c("SlateBlue","SeaGreen","red"),
          xlab="", ylab="", ylim=c(ylimmin, ylimmax+rem),
          xlim = c(anyos[1],anyos[n]+1.3),
#         main = "Evolucion de cada mes",
         lty = c(1,2,3,4), plot.type="single")
  for(i in 7:9){
   ord <- length(na.omit(MMp[,i]))
   text(anyos[ord]+1.1, MMp[ord,i],
        leyenda[mp[i]], col = colour[i], cex=0.7)
               }
  plot.ts(MMp[,10:12], col = c("SlateBlue","SeaGreen","red"),
          xlab="", ylab="", ylim=c(ylimmin, ylimmax+rem),
          xlim = c(anyos[1],anyos[n]+1.3),
#         main = "Evolucion de cada mes",
         lty = c(1,2,3,4), plot.type="single")
  for(i in 10:12){
   ord <- length(na.omit(MMp[,i]))
   text(anyos[ord]+1.1, MMp[ord,i],
        leyenda[mp[i]], col = colour[i], cex=0.7)
  }
#  par(opar)
 }
                }
 MMp
}

bbap <- function(vari, s, t0, yearsp)
{
# anyosp son los anyos a dibujar.
 anyosp <- yearsp
 N    <- length(vari)
 #naux <- length(vari)/s
 #ifelse(length(vari)/s == as.integer(naux), n <- naux, n <- as.integer(naux)+1)
 n <- ysooys(N, t0, N, s)[[1]][1] - ysooys(1, t0, N, s)[[1]][1] + 1
 r <- c(1:(n*s-t0[2]+1))
 raux1 <- rep(1,(s+1)-t0[2])
 raux2 <- gl(n, s)
 r[1:length(raux1)] <- raux1
 r[(1+length(raux1)):length(r)] <- raux2[(s+1):length(raux2)]

 anyos <- c(1:n)
 ap <- c(1:length(anyosp))
 for(i in 0:n)
   anyos[i+1] <- t0[1]+i
 for(i in 1:length(anyosp))
   ap[i] <- which(anyos == anyosp[i])

 MMth <- matrix(nrow=n, ncol=s)
 MA   <- matrix(nrow=n, ncol=s)
 MAp  <- matrix(nrow = length(anyosp), ncol = s) #n ) ???

 colour = rep(c("SlateBlue", "SeaGreen", "red", "magenta", "navy", "sienna"),2)
# colour = as.character(gl(6, 2, label = c("SlateBlue", "SeaGreen",
#                     "red", "magenta", "navy", "sienna")))

 i <- 1
 while(i < (length(vari)+1)){
  if(i < (s+1)-t0[2])
    for(j in t0[2]:s){
      MMth[r[i],j] <- vari[i]
      i <- i+1
    }
  else
    for(j in 1:s){
      MMth[r[i],j] <- vari[i]
      i <- i+1
    }
 }
 anyos <- c(1:n)
 ap <- c(1:length(anyosp))
 for(i in 0:n)
   anyos[i+1] <- t0[1]+i
 for(i in 1:length(anyosp))
   ap[i] <- which(anyos == anyosp[i])

 j <- 1
 while(j < (length(anyosp)+1)){
   for(i in ap[j])
     MA[i,] <- ts(MMth[i,], frequency = 1, start = t0[1])
   j <- j+1
 }
 data.frame(MA[,1:s])

 Na <- c(1:s)
 j <- 1
 for(i in 1:nrow(MA))
 {
   for(k in 1:s)
      Na[k] <- MA[i,k] =="NA"
   if(length(which(Na == 0))>1){
      MAp[j,] <- MA[i,]
      j <- j+1
   }
   else
      j <- j
 }
 if(s==4)
   textx <- c("Qrtr1", "Qrtr2", "Qrtr3", "Qrtr4")
 if(s==12)
   textx <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
# textx <- c("Ene","Feb","Mar","Abr","May","Jun","Jul", "Ago","Sep","Oct","Nov","Dic")
# textx <- c("Trm1", "Trm2", "Trm3", "Trm4")

# opar <- par(mar=c(8,4.7,5,1.5), ps=18, font=1,tcl=-0.5, cex.axis=0.7, las=1)
 plot.ts(t(MAp), lty = c(1,2,3,4), plot.type="single", las=1,
#        main = "Anual path", cex.main=0.9
         col = c("SlateBlue","SeaGreen","red","magenta", "navy", "sienna"),
         xlim = c(1,s+0.5), xlab="", ylab="", xaxt="n", cex.axis=0.7)
 for(i in 1:length(ap)){
  for(k in 1:s)
    Na[k] <- MAp[i,k] == "NA"
  ifelse(i==1, ord <- s, ord <- length(na.omit(MAp[i,])))
  text(ord+0.5, MAp[i,which(Na==0)[length(which(Na==0))]],
      anyosp[i], cex=0.9, col = colour[i])   # al lado del primer mes/trimestre
  # text(ord+0.5, MAp[i,which(Na==0)[1]],   # al lado del último mes/trimestre
  #    anyosp[i], cex=0.7, col = colour[i])  #, font=2)
                       }
 axis(side=1, textx, at=c(1:s), tcl=-0.5)
# par(opar)
}

#
CH.test <- function(label, frec, f0, DetTr, showcat)
{
  #vari.env <- new.env(FALSE, NULL)
  #Fvarinfo(label$label, env.out=vari.env)
  #vari <- get("vari", env=vari.env)
  #s    <- get("s", env=vari.env)
  #t0   <- get("t0", env=vari.env)
  #N    <- length(vari)
  vari <- label$vari
  s    <- label$s
  t0   <- label$t0
  N    <- length(vari)

# Regresión auxiliar

 VFEtrg <- MVFE(vari, s, t0, "trg")
 VFEalg <- MVFE(vari, s, t0, "alg")
 tiempo <- c(1:N)

 if(DetTr ==FALSE){
   if(f0 == 1){
     lm2  <- lm(vari[2:N] ~ vari[1:(N-1)] + VFEtrg[2:N,1:(s-1)])
     ehat <- matrix(lm2$residuals, ncol=1)
              }
   if(f0 == 0){
     lm2  <- lm(vari ~ VFEtrg[,1:(s-1)])
     ehat <- matrix(lm2$residuals, ncol=1)
              }
                  }
 if(DetTr ==TRUE){
   if(f0 == 1){
     lm2  <- lm(vari[2:N] ~ vari[1:(N-1)] + tiempo[1:(N-1)] +
                            VFEtrg[2:N,1:(s-1)])
     ehat <- matrix(lm2$residuals, ncol=1)
              }
   if(f0 == 0){
     lm2  <- lm(vari ~ tiempo + VFEtrg[,1:(s-1)])
     ehat <- matrix(lm2$residuals, ncol=1)
              }
                 }
# R1

 R1 <- VFEtrg
 Fhat    <- matrix(nrow=length(ehat), ncol=(s-1))
 Fhataux <- matrix(nrow=length(ehat), ncol=(s-1))

 # ifelse(f0 == 1, op <- 1, op <- 0)# también cambia varnw[3:13...]
 for(i in 1:length(ehat))
   Fhataux[i,] <- R1[i+f0,]*ehat[i]   # No producto matricial %*%

 for(i in 1:nrow(Fhat)){
   for(n in 1:(s-1))
     Fhat[i,n] <- sum(Fhataux[1:i,n]) # F estimada
}

# Omega supra-f estimada  (Omfhat)

ltrunc <- round(4*(N/100)^0.25)
Omfhat <- Omegaf(N, s, ehat, Fhataux)

# Matriz A
 # frec <- c(0,1) # input (no analizar la frecuencia de largo plazo)
 # frec <- c(1,0,0,0,0,0)  # input mensual

 sq     <- seq(1,11,2)
 frecob <- rep(0,(s-1))

  for(i in 1:length(frec)){
   if(frec[i] == 1 && i == s/2)
      frecob[sq[i]]     <- 1
   if(frec[i] == 1 && i < s/2)
      frecob[sq[i]] <- frecob[sq[i]+1] <- 1
  }

 a <- length(which(frecob == 1))
 A <- matrix(0, nrow=s-1, ncol=a)

 j <- 1
 for( i in 1:(s-1) )
  if(frecob[i] == 1){
      A[i,j] <- 1
      ifelse(frecob[i] == 1, j <- j+1, j <- j)
  }


# Estadístico L (ecuación (13) Canova & Hansen (1995))

 estdL <- (1/N^2)*sum(diag( solve(t(A) %*% Omfhat %*% A)
                           %*%
                           t(A) %*% t(Fhat) %*% Fhat %*% A ))

 # Resumen de resultados
 if(showcat==TRUE){
   cat("\n ------ CH test ------ \n\n")
   cat(c("L-statistic:", round(estdL, 2), " \n\n"))
   cat(c("Lag truncation parameter:", ltrunc, " \n\n"))
                  }
  list(L.Stat=estdL, lag.trunc=ltrunc)
}
#
CHseas.test <- function(label, lmax, seas, showcat)
{
  vari <- label$vari
  s    <- label$s
  t0   <- label$t0
  N    <- length(vari)

  if(s==4)
    Mper <- quarterg(vari, s, t0, plot=FALSE)
  if(s==12)
    Mper <- bbmp(vari, s, t0, c(1:12), "R", plot=FALSE)

  if(length(seas)==1)
  {
    rdoCH <- matrix(nrow=lmax+1, ncol=2)
    for(l in 0:lmax){
      rdoCH[(l+1),1] <- as.numeric(KPSS.test(Mper[,(seas+t0[2]-1)], l, showcat=FALSE)[1])
      rdoCH[(l+1),2] <- as.numeric(KPSS.test(Mper[,(seas+t0[2]-1)], l, showcat=FALSE)[2])
                    }
   # Resumen de resultados
   if(showcat==TRUE){
      lret <- matrix(c(0:lmax), ncol=1)
      cat("\n ------ CH test ------ \n\n")
      cat(" Statistics for the null hipothesis of \n level stationarity: \n")
      print(data.frame(Lags=lret, Statistic=rdoCH[,1]))
      cat("\n Statistics for the null hipothesis of \n trend stationarity:\n ")
      print(data.frame(Lags=lret, Statistic=rdoCH[,2]))
                    }
      list(data.frame(Lags=lret, Season=rdoCH[,2]))

  }
  if(length(seas)==s)
  {
    rdoCH <- matrix(nrow=lmax+1, ncol=2*s)
    for(i in 1:s){
      for(l in 0:lmax){
        rdoCH[(l+1),i]     <- as.numeric(KPSS.test(Mper[,i], l, showcat=FALSE)[1])
        rdoCH[(l+1),(i+s)] <- as.numeric(KPSS.test(Mper[,i], l, showcat=FALSE)[2])
                      }
                 }
    rdoCHaux <- matrix(nrow=lmax+1, ncol=2*s)
    rdoCHaux[,t0[2]:s] <- rdoCH[,1:(s+1-t0[2])]
    if(t0[2] != 1)
      rdoCHaux[,1:(t0[2]-1)] <- rdoCH[,(s+2-t0[2]):s]
    rdoCHaux[,(t0[2]+s):(2*s)] <- rdoCH[,(s+1):(2*s+1-t0[2])]
    if(t0[2] != 1)
      rdoCHaux[,(s+1):(t0[2]-1+s)] <- rdoCH[,(2*s+2-t0[2]):(2*s)]
    rdoCH <- rdoCHaux

   # Resumen de resultados
   if(showcat==TRUE){
     lret <- matrix(c(0:lmax), ncol=1)
     cit1 <- "Statistics for the null hipothesis of \n level stationarity:\n"
     if(s==4)
       cit2 <- data.frame(Lags=lret, Qrtr1=rdoCH[,1],
               Qrtr2=rdoCH[,2], Qrtr3=rdoCH[,3], Qrtr4==rdoCH[,4])
     if(s==12)
       cit2 <- data.frame(Lags=lret,
         January=rdoCH[,1], February=rdoCH[,2], March=rdoCH[,3], April=rdoCH[,4],
         May=rdoCH[,5], June=rdoCH[,6], July=rdoCH[,7], August=rdoCH[,8],
         September=rdoCH[,9], October=rdoCH[,10], November=rdoCH[,11], December=rdoCH[,12])

     cit3 <- "\n Statistics for the null hipothesis of \n trend stationarity:\n"
     if(s==4)
       cit4 <- data.frame(Lags=lret, Qrt1=rdoCH[,5],
               Qrt2=rdoCH[,6], Qrt3=rdoCH[,7], Qrt4==rdoCH[,8])
     if(s==12)
       cit4 <- data.frame(Lags=lret,
         January=rdoCH[,13], February=rdoCH[,14], March=rdoCH[,15], April=rdoCH[,16],
         May=rdoCH[,17], June=rdoCH[,18], July=rdoCH[,19], August=rdoCH[,20],
         September=rdoCH[,21], October=rdoCH[,22], November=rdoCH[,23], December=rdoCH[,24])
     cat("\n ------ CH test ------ \n\n")
     cat(cit1)
     print(cit2)
     cat(cit3)
     print(cit4)
                   }
  }
  rdoCH
}

#

# DGP: Phi(L)*y_t = rho(L)eps_t,  eps_t~ N(0,1)
DGPsim <- function(s, N, phi, retphi, rho, retrho)
{
  if(nrow(phi) == 1){
    fcoef <- phi
    rfil  <- retphi
                    }
  if(nrow(phi) == 2){
    crpp.out <- crpp(phi[1,], retphi[1,], phi[2,], retphi[2,])
    fcoef <- crpp.out[[1]]
    rfil  <- crpp.out[[2]]
                    }
  if(nrow(phi) > 2){
    crpp.out <- crpp(phi[1,], retphi[1,], phi[2,], retphi[2,])
    for(i in 3:nrow(phi))
      crpp.out <- crpp(crpp.out[[1]], crpp.out[[2]], phi[i,], retphi[i,])
    fcoef <- crpp.out[[1]]
    rfil  <- crpp.out[[2]]
                   }

  N     <- N + max(rfil)
  DGP   <- matrix(nrow = N, ncol = 1)
  t0    <- c(0, 1)
  eps   <- matrix(rnorm(N, mean = 0, sd = 1), ncol = 1)
  res   <- rep(0, length(retrho))
  MLeps <- ret(eps, max(retrho) + 1)

  DGPaux <- rep(0, length(rfil))
  DGP[1:max(rfil),] <- 0
  for (i in (max(rfil)+1):N)
  {
     for (j in 2:length(rfil))
       DGPaux[j] <- -fcoef[j]*DGP[(i-rfil[j]),]
     if(rho[1] != 0){
        for (j in 1:length(retrho))
           res[j] <- rho[j]*MLeps[i,(retrho[j] + 1)]
        res <- na.omit(res)
                    }
     DGP[i] <- sum(DGPaux) + eps[i] + sum(res)
  }
  DGP <- ts(DGP[(max(rfil)+1):N], frequency = s, start = t0)
  DGP
}

#
freqg <- function(label)
{
 vari.env <- new.env(FALSE, NULL)
 #Fvarinfo(label$label, env.out=vari.env)
 #vari <- get("vari", env=vari.env)
 #s    <- get("s", env=vari.env)
 #t0   <- get("t0", env=vari.env)
 #N    <- length(vari)
 vari <- label$vari
 s    <- label$s
 t0   <- label$t0
 N    <- length(vari)

 if(s == 12)
 {
  ML <- ret(vari, 12)

  Myt    <- matrix(nrow=N, ncol=7)
  for(f in 1:7){
    filtra    <- rep(1, 7)
    filtra[f] <- 0
    Fil.vari  <- filtrar(vari, s, t0, filtra, plot=FALSE)[[1]]
    Myt[,f]   <- Fil.vari
               }
  Myt <- ts(Myt, frequency = s, start = t0)

# par(bg="cornsilk")
# opar <- par(mfrow=c(4,2), mar=c(3,3,3.5,2),
#             tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)

   plot(diff(vari, lag=s), xlab="", ylab="", main="Seasonal difference")

   plot(Myt[,1], cex.axis=0.7, xlab="", ylab="")
   title("Long run", font.main=2, adj=0.5)
   mtext("Frequency zero", side = 3, line = 0.35, cex=0.7)

   plot(Myt[,2], cex.axis=0.7, xlab="", ylab="")
   title("Bi-monthly", font.main=2, adj=0.5)
   mtext("pi frequency", side = 3, line = 0.35, cex=0.7)

   plot(Myt[,3], cex.axis=0.7, xlab="", ylab="")
   title("Four-monthly", font.main=2, adj=0.5)
   mtext("Frequencies (pi/2; 3pi/2)", side = 3, line = 0.35, cex=0.7)

   plot(Myt[,4], cex.axis=0.7, xlab="", ylab="")
   title("Quarterly", font.main=2, adj=0.5)
   mtext("Frequencies (2pi/3; 4pi/3)", side = 3, line = 0.35, cex=0.7)

   plot(Myt[,5], cex.axis=0.7, xlab="", ylab="")
   title("Semianual", font.main=2, adj=0.5)
   mtext("Frequencies (pi/3; 5pi/3)", side = 3, line = 0.35, cex=0.7)

   plot(Myt[,6], cex.axis=0.7, xlab="", ylab="")
   title(" ", font.main=2, adj=0.5)
   mtext("Frequencies (5pi/6; 7pi/6)", side = 3, line = 0.35, cex=0.7)

   plot(Myt[,7], cex.axis=0.7, xlab="", ylab="")
   title("Anual", font.main=2, adj=0.5, cex.main=1)
   mtext("Frequencies (pi/6; 11pi/6)", side = 3,
         line = 0.35, cex=0.7)
#   par(opar)
 }

# S E R I E S  T R I M E S T R A L E S

 if(s == 4)
 {
 ML <- ret(vari, 5)

 Myt    <- matrix(nrow=N, ncol=7)
 for(f in 1:3){
   filtra    <- rep(1, 3)
   filtra[f] <- 0
   Fil.vari  <- filtrar(vari, s, t0, filtra, plot=FALSE)[[1]]
   Myt[,f]   <- Fil.vari
              }
 Myt <- ts(Myt, frequency = s, start = t0)

# opar <- par(mfrow=c(2,2), mar=c(4,3,4,1.5), font=1,
#             tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
   plot(diff(vari, lag=s), xlab="", ylab="", main="Seasonal difference")

   plot(Myt[,1], cex.axis=0.7, xlab="", ylab="")
   title("Long run", font.main=2, adj=0.5)
   mtext("Frequency zero", side = 3,
         line = 0.35, cex=0.7)

   plot(Myt[,2], cex.axis=0.7, xlab="", ylab="")
   title("Semianual", font.main=2, adj=0.5)
   mtext("Frequency pi", side = 3,
         line = 0.35, cex=0.7)

   plot(Myt[,3], cex.axis=0.7, xlab="", ylab="")
   title("Anual", font.main=2, adj=0.5)
   mtext("Frequencies (pi/2; 3pi/2)", side = 3,
         line = 0.35, cex=0.7)
#  par(opar)
 }
}

bb3D <- function(MR, s, t0, color, x, y)
{
  if(color==TRUE)
     colores <- "lightgoldenrod"
  if(color==FALSE)
     colores <- grey((0:6)/6)
  if(s==12){ xlabel <- "months"; ntic <- s/2}
  if(s==4) { xlabel <- "quarters"; ntic <- s}
  xx <- c(1:ncol(MR))
  yy <- c(t0[1]:(t0[1]+nrow(MR)-1))

  persp(xx, yy, t(MR), # main="Buys-Ballot 3D"
        theta = x, phi = y, expand = 0.5,
        xlab=xlabel, ylab="", zlab="", shade=0.4,
        col = colores, ticktype="detailed", nticks=ntic)
}

bbcn <- function(MR, s, t0, color)  # MMp ó Mq de bbmp() y quarterg()
{
  if(color==TRUE)
     colores <- terrain.colors(200)
  if(color==FALSE)
     colores <- grey((0:32)/32)
  if(s==12){
     xlabel <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
           }
  if(s==4){ xlabel <- c("Qrtr1", "Qrtr2", "Qrtr3", "Qrtr4")}
  x <- c(1:ncol(MR))
  y <- c(t0[1]:(t0[1]+nrow(MR)-1))

  image(x, y, t(MR), # main="Buys-Ballot curvas de nivel"
        las=1, xlab="", ylab="", xaxt="n", col = colores)
  mtext(xlabel[1:s], side=1, line=1, at=c(1:s)) #, cex=0.7)
  contour(x, y, t(MR), add = TRUE, drawlabels = TRUE, col="blue")
}

seasboxplot <- function(label, color)
{
  #vari.env <- new.env(FALSE, NULL)
  #Fvarinfo(label$label, env.out=vari.env)
  #vari <- get("vari", env=vari.env)
  #s    <- get("s", env=vari.env)
  #t0   <- get("t0", env=vari.env)
  #N    <- length(vari)
  vari <- label$vari
  s    <- label$s
  t0   <- label$t0
  N    <- length(vari)

  if(s==4){
    xnames <- c("Qrtr1", "Qrtr2", "Qrtr3", "Qrtr4")
    MS     <- quarterg(.wvari$vai, .wvari$s, .wvari$t0, plot=FALSE)
          }
  if(s==12){
    xnames <- c("January", "February", "March", "April", "May", "June", "July",
               "August", "September", "October", "November", "December")
    MS     <- bbmp(.wvari$vari, .wvari$s, .wvari$t0, c(1:12), "Prot", FALSE)
           }

  if(color==TRUE){color1 <- "lightblue"; color2 <- "SeaGreen2"}
  if(color==FALSE){color1 <- "lightgray"; color2 <- "gray60"}
  opar <- par(las=2) # cex.axis=0.7, cex.main=0.8,
  summ.box <- boxplot(split(vari, cycle(vari)), names=xnames,col=color1)
  #                   main="Gráfico de cajas estacional"
  boxplot(split(vari, cycle(vari)), names=xnames, notch=TRUE, add=TRUE, col=color2)
  par(opar)

  if(s==4){
    stats <- round(data.frame("Qrtr1"=summ.box[[1]][,1], "Qrtr2"=summ.box[[1]][,2],
           "Qrtr3"=summ.box[[1]][,3], "Qrtr4"=summ.box[[1]][,4]), 2)
           }
  if(s==12){
    stats <- round(data.frame("January"=summ.box[[1]][,1], "February"=summ.box[[1]][,2],
           "March"=summ.box[[1]][,3], "April"=summ.box[[1]][,4],
           "May"=summ.box[[1]][,5], "June"=summ.box[[1]][,6],
           "July"=summ.box[[1]][,7], "August"=summ.box[[1]][,8],
           "September"=summ.box[[1]][,9], "October"=summ.box[[1]][,10],
           "November"=summ.box[[1]][,11], "December"=summ.box[[1]][,12]), 2)
           }
  cat("\n ------ Seasonal box-plot ------ \n\n")
  cat("     Statictics (by rows): \n")
  cat("     minimum, quartile1, median, quartile3, maximun. \n \n")
  print(stats)
  cat("\n\n Confidence interval of 90% for each median: \n\n")
  print(round(summ.box$conf),2)
  if(length(summ.box$out)){
  cat("\n\n The oints wich lie beyond the extremes of the whiskers are: de\n\n")
  print(round(summ.box$conf,2))
  cat("\n\n and are refered to the following seasons:")
  print(summ.box$group)
                          }
}

#
HEGY.test <- function(label, compdet, selecP, Mvfic, VFEp, showcat)
{
 vari <- label$vari
 s    <- label$s
 t0   <- label$t0
 N    <- length(vari)

 DGP  <- matrix(vari, nrow=N, ncol=1)
 EtFs <- matrix(nrow=1, ncol=(2+3*s)/2)
 Ct   <- matrix(rep(1, N), ncol=1)
 TD   <- matrix(c(1:N), ncol=1)
 VFE  <- MVFE(DGP[,1], s, t0, "alg")[,1:(s-1)]
 ifelse(Mvfic==0, nvf1<-0, nvf1<-ncol(Mvfic))
 ifelse(VFEp==0, nvf2<-0, nvf2<-ncol(VFEp))
 ldet <- compdet[1]*1 + compdet[2]*1 + compdet[3]*(s-1) + nvf1 + nvf2

# Dependent variable

 Deltay <- matrix(NA, nrow=N, ncol=1)
   Deltay[(s+1):N,1] <- DGP[(s+1):N,1] - DGP[1:(N-s),1]

# Frequency regressors

 Myr <- regr.frec(DGP[,1], s, t0)

# Deterministic components regressors

 Mregdet <- as.matrix(data.frame(Ct, TD, VFE, Mvfic, VFEp, Myr))
 aux <- rep(1, ncol(Mregdet))
 for(i in 1:ncol(Mregdet)){
   if(length(which(Mregdet[,i]==0))==N)
   aux[i] <- 0
                          }
 Mregdet <- Mregdet[,which(aux==1)]
 cdf1 <- which(compdet[1:2] == 1)
 if(compdet[3] == 1){ cdf2 <- c(3:(s+1)) }
 if(compdet[3] == 0){ cdf2 <- NULL }
 if(length(which(Mvfic==0))>1 || length(which(VFEp==0))>1)
   cdf3 <- c((s+2):(s+2+nvf1+nvf2-1))
 if(length(which(Mvfic==0))==1 && length(which(VFEp==0))==1)
   cdf3 <- NULL
 cdf4 <- c((s+2+nvf1+nvf2):(2*s+1+nvf1+nvf2))
 cdf  <- c(cdf1, cdf2, cdf3, cdf4)
 lmdet <- lm(Deltay[,1] ~ 0+Mregdet[,cdf])

# Adding lags of Deltay, using slecP function

 Pmax <- rsug <- -1
 if(mode(selecP)=="character"){
 if(selecP == "aic")
    Pmax <- selecPabic(Deltay, s, lmdet)[1]
 if(selecP == "bic")
    Pmax <- selecPabic(Deltay, s, lmdet)[2]
 if(selecP == "aiclb")
   Pmax <- selecPv5(Deltay, s, lmdet)
 if(selecP == "biclb")
   Pmax <- selecPv6(Deltay, s, lmdet)
 if(selecP == "aiclut")
    rsug <- selecPv2(Deltay, s, lmdet)
 if(selecP == "biclut")
    rsug <- selecPv3(Deltay, s, lmdet)
 if(selecP == "signf")
    rsug <- selecPv4(Deltay, s, lmdet)
                         }
 if(mode(selecP) == "numeric"){
    if(length(selecP)==1)
       Pmax <- selecP     # si se indica un orden elegido por uno mismo
    if(length(selecP)>1)  # cuidado si se quiere sólo un retardo > 1
       rsug <- selecP
                              }
 if(Pmax > rsug[1] && Pmax > 0)
    Deltayr <- as.matrix(ret(Deltay, Pmax+1)[,2:(Pmax+1)])
 if(Pmax < rsug[1] && length(rsug) > 0)
    Deltayr <- as.matrix(ret(Deltay, round(10*log10(N))+2)[,(rsug+1)])
 ifelse(Pmax > 0 || length(which(rsug!=-1)) > 0,
        nP <- ncol(Deltayr), nP <- 0)

# Auxiliar regression estimates

 if(Pmax > 0)
   lmhegy <- lm(Deltay[,1] ~ 0+Mregdet[,cdf] + Deltayr[,1:nP])
 if(Pmax == 0)
   lmhegy <- lm(Deltay[,1] ~ 0+Mregdet[,cdf])

 if(rsug[1] != -1 && length(rsug) > 0)
   lmhegy <- lm(Deltay[,1] ~ 0+Mregdet[,cdf] + Deltayr[,1:nP])
 if(length(rsug) == 0)
   lmhegy <- lm(Deltay[,1] ~ 0+Mregdet[,1:cdf])

# Statistics

 if(s == 4) { Fs1.s <- c(1,2,3,4) + ldet
              Fs2.s <- c(2,3,4) + ldet }
 if(s == 12){ Fs1.s <- c(1,2,3,4,5,6,7,8,9,10,11,12) + ldet
              Fs2.s <- c(2,3,4,5,6,7,8,9,10,11,12) + ldet }

 for(j in (ldet+1):(ldet+s)){
    EtFs[1,(j-ldet)] <- contts(lmhegy, j)[[2]]
                            }
 iaux <- seq((3+ldet), (s+ldet), 2)
 for(j in 1:length(iaux))
   EtFs[1,(j+s)] <- Fsnd(lmhegy, 2, s+ldet+nP, c(iaux[j],iaux[j]+1))
 EtFs[1,ncol(EtFs)-1] <- Fsnd(lmhegy, (s-1), s+ldet+nP, Fs2.s)
 EtFs[1,ncol(EtFs)]   <- Fsnd(lmhegy, s, s+ldet+nP, Fs1.s)

 rdofrec <- matrix(nrow=1, ncol=ncol(EtFs))
 for(n in 1:s)
   rdofrec[,n] <- round(EtFs[,n], 2)
 for(n in (s+1):ncol(EtFs))
   rdofrec[,n] <- round(EtFs[,n], 2)
 rdofrec <- list(tpi= rdofrec[,1:s], Fpi= rdofrec[,(s+1):ncol(EtFs)])
 #
 if(compdet[1] == 0){ rdodet=NULL }
 if(compdet[1] == 1){
   rdodet <- matrix(nrow=ldet, ncol=2)
   for(j in 1:(ldet)){
     rdodet[j,1] <- lmhegy$coef[j]
     rdodet[j,2] <- contts(lmhegy, j)[[2]]
                     }
   Coeff    <- round(rdodet[,1], 2)
   t.stat <- round(rdodet[,2], 2)

 if(compdet[3] == 1){
   if(s==4) {Fvfeaux <- c(2,3,4) + compdet[2]}
   if(s==12){Fvfeaux <- c(2,3,4,5,6,7,8,9,10,11,12) + compdet[2]}
   F.VFE  <- Fsnd(lmhegy, s-1, s+ldet+nP, Fvfeaux)
   rdodet <- round(data.frame(Coeff, t.stat, F.VFE), 2)
                    }
 if(compdet[3] == 0)
   rdodet <- round(data.frame(Coeff, t.stat), 2)
                    }
 ifelse(rsug==-1, retardos <- Pmax, retardos <- rsug)

 # Summary

 if(showcat==TRUE){
  if(s==4){
     frlabel <- c("zero", "pi", "pi/2", "")
     tlabel  <- c("t1", "t2", "t3", "t4")
     Flabel  <- c("F.3,4", "F.2,4", "F.1,4", "")
     Faux    <- c(rdofrec[[2]], "")
           }
  if(s==12){
     frlabel <- c("zero", "pi", "pi/2", "", "2pi/3", "", "pi/3", "", "5pi/6", "", "pi/6", "")
     tlabel  <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11", "t12")
     Flabel  <- c("F.3,4", "F.5,6", "F.7,8", "F.9,10", "F.11,12", "F.2,12", "F.1,12",
                  "", "", "", "", "")
     Faux    <- c(rdofrec[[2]], rep("", 5))
           }
  freqrdo <- data.frame(Frequencies=frlabel, t.=tlabel, stat=rdofrec[[1]],
                        F.=Flabel, test=Faux)

  ifelse(compdet[1]==1, cd1 <- "Intercept", cd1 <- NA)
  ifelse(compdet[2]==1, cd2 <- "Trend", cd2 <- NA)
  ifelse(compdet[3]==1, cd3 <- c("Seasonal dummys", rep("", s-2)), cd3 <- NA)
  ifelse(nvf1!=0, cd4 <- c("Generic dummys", rep("", nvf1-1)), cd4 <- NA)
  ifelse(nvf2!=0, cd5 <- c("Partial seasonal dummys", rep("", nvf2-1)), cd5 <- NA)
  cdlabel <- na.omit(c(cd1, cd2, cd3, cd4, cd5))
  cdrdo   <- data.frame(Components=cdlabel, rdodet[,1:2])

    cat(c("\n ------ HEGY test ------ \n\n"))
    cat(c(" Statistics about the frequencies: \n "))
    cat(c("   t-statistics:  \n "))
    print(freqrdo[,1:3])
    cat(c("\n   F-statistics: \n "))
    print(freqrdo[1:(s/2+1),4:5])
    cat("\n Deterministic components estimates: \n\n")
    print(cdrdo)
    if(compdet[3]==1)
      cat(c("Seasonal dummys F.test: ", rdodet[1,3], "\n"))
    cat(c("\n Number of lags selected:\n", retardos, "\n\n"))
    cat(c(" Number of effective observations: ", (N-s-max(retardos)), "\n\n"))
                  }
  list(rdofrec, compdeter=rdodet, Lags=retardos, nobs=(N-s-max(retardos)))
}

#

KPSS.test <- function(vari, l, showcat)
{
  N    <- length(vari)
  tiempo <- c(1:N)
  ML <- ret(vari, 2)
  #l      <- as.integer(3*sqrt(N)/13)  # ver PP.test, artículo KPSS
   #l      <- as.integer(10*sqrt(N)/14)
   #l      <- as.integer(4*(N/100)^(1/4))
   #l      <- as.integer(12*(N/100)^(1/4))

 lmvari <- lm(ML[,1] ~ tiempo)
 ehat   <- residuals(lmvari)
 Sa     <- cumsum(ehat)
 N <- length(ehat)

 # estimador consistente de la varianza de los residuos:
 if(l==0)
   s.2a <- 1/N*sum(ehat^2)
 if(l>0){
   auxa <- c(1:l)
   for(i in 1:l)
     auxa[i] <- (1-i/(l+1))*sum(ehat[(i+1):N]*ehat[1:(N-i)])
   s.2a <- (1/N)*sum(ehat^2)+(2/N)*sum(auxa)
        }
 Trend <- round(N^(-2)*sum(Sa^2/s.2a), 2)

 # para la hipótesis nula de estacionariedad en tendencia los
 # residuos son ee=var-mean(var),
 # se analiza la serie extrayendo la media y la tendencia.

 vari2 <- vari - mean(na.omit(vari))
 Sb    <- cumsum(vari2)

 # estimador consistente de la varianza de los residuos:
 if(l==0)
   s.2b <- 1/N*sum(vari2[1:N]^2)
 if(l>0){
   auxb <- c(1:l)
   for(i in 1:l)
     auxb[i] <- (1-i/(l+1))*sum(vari2[(i+1):N]*vari2[1:(N-i)])
   s.2b <- (1/N)*sum(vari2[1:N]^2)+(2/N)*sum(auxb)
        }
 Level <- round(N^(-2)*sum(Sb[1:N]^2/s.2b), 2)

 #
 if(showcat==TRUE){
     cat("\n------ KPSS test ------ \n\n")
     cat(c(" Statistic for the null hipothesis of \n level stationarity:", Level," \n\n "))
     cat(c("Statistic for the null hipothesis of \n trend stationarity:", Trend," \n\n "))
                  }
 rdok <- data.frame(Level=Level, Trend=Trend)
 rdok
 # list(rdok)
 #list(rdok, rdo)
}

#

Omegaf <- function(N, s, ehat, Fhataux)
{
 m   <- round(4*(N/100)^0.25)
 wnw <- c(1:(m+1))
 j   <- 1
 for(k in -m:m)
 {
    wnw[j] <- (m + 1 - abs(k))/(m+1)
    j <- j+1
 }
 Ne   <- length(ehat)
 Omnw <- 0
 j    <- 1
 for(k in -m:m)
 {
    Om <-
      (t(Fhataux)[,(m+1+k):(Ne-(m-k))] %*% Fhataux[(m+1):(Ne-m),])
    Omnw <- Omnw + wnw[j]*Om
    j <- j+1
 }
 OmNWfaux <- Omnw/N

 OmNWf <- matrix(nrow=nrow(OmNWfaux), ncol=ncol(OmNWfaux))
 diag(OmNWf) <-  diag(OmNWfaux)
 j <- 1
 for(i in 2:(s-1))
 {
   OmNWf[j,i:(s-1)] <- OmNWf[i:(s-1),j] <- OmNWfaux[i:(s-1),j]
   j <- j+1
 }
 OmNWf
}

#
# Elige el número de retardos en función del AIC ó BIC
# Considero k como el número de parámetros:
  # número de retardos más el término constante

AICBIC <- function(lmabic)       
{
   N <- nrow(model.matrix(lmabic))
   k <- ncol(model.matrix(lmabic))
   sigML <- sum(lmabic$residuals^2)/N
   aic   <- N*log(sigML) + 2*k
   bic   <- N*log(sigML) + k*log(N)
   ABIC  <- data.frame(AIC=aic, BIC=bic)
   ABIC
}

#
AICBICaux <- function(vari, s, lmdet)
{
  N    <- length(vari)
  bp   <- round(10*log10(N))
  abic <- matrix(ncol=2, nrow=bp+1)
  ML   <- ret(vari, bp+1)
  k1   <- ncol(model.matrix(lmdet))
  ndet <- ncol(model.matrix(lmdet))
  miss <- s-length(which(vari[1:s]!="NaN"))+1

  lmabic    <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet])
  sigML     <- sum(lmabic$residuals^2)/N
  abic[1,1] <- N*log(sigML) + 2*k1         # aic
  abic[1,2] <- N*log(sigML) + k1*log(N)    # bic

  for(i in 2:(bp+1)){
    lmabic <- lm(ML[miss:N,1] ~
              0+model.matrix(lmdet)[,1:ndet] + ML[miss:N,2:i])
    sigML     <- sum(lmabic$residuals^2)/N
    abic[i,1] <- N*log(sigML) + 2*(i-1+k1)
    abic[i,2] <- N*log(sigML) + (k1+i-1)*log(N)
                    }
  ABIC <-  data.frame(AIC=abic[,1], BIC=abic[,2])
}
#
selecPabic <- function(vari, s, lmdet)  # selecPv1 y selecPv2
{
  Pmax    <- rep(NA, 2)
  ABIC    <- AICBICaux(vari, s, lmdet)
  Pmax[1] <- which.min(ABIC[,1])-1    # AIC
  Pmax[2] <- which.min(ABIC[,2])-1    # BIC
    # -1 porque el primer valor es con cero retardos
  Pmax
}

# SelecPv2 y SlecPv3 puede dar problemas si el valor de ML[N,j] es cero
  # y el de la serie original en ese retardo también es cero

# SelePv2:
  # Selecciona no sólo el número de retardos sino
    # también el orden de cada uno
    # a partir del AIC y como explica Lütkephoul(1990)
  # Regresión auxiliar sin componentes deterministas

selecPv2 <- function(vari, s, lmdet)
{
   # N  <- length(na.omit(vari))
   N    <- length(vari)
   bp   <- round(10*log10(N))
   ML   <- ret(vari, bp + 1)
   ndet <- ncol(model.matrix(lmdet))
   miss <- s - length(which(vari[1:s] != "NaN")) + 1

   aic    <- rep(NA, 2)
   lmarP  <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet] + ML[miss:N,2:(bp+1)])
   aic[1] <- AICBIC(lmarP)[1]
   lmarP  <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet] + ML[miss:N,2:(bp)])
   aic[2] <- AICBIC(lmarP)[1]

   k <- bp-1
   while(k > 0){
     if(aic[2] < aic[1]){
        aic[1] <- aic[2]
        ML[,(k+2)] <- rep(0, length(vari))
                        }
     else
        aic[1] <- aic[1]
     lmarP <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet] +
                          ML[miss:N,2:k] + ML[miss:N,(k+2):(bp+1)])
     aic[2] <- AICBIC(lmarP)[1]
     k <- k-1
               }
   rsugaux <- matrix(nrow=1, ncol=bp)
   j <- 2
   while(j < bp+2){
     ifelse(ML[N,j] != 0, rsugaux[j-1] <- TRUE,
                          rsugaux[j-1] <- FALSE)
     j <- j+1
                  }
   rsug <- which(rsugaux == TRUE)
   rsug
}

#

# SelePv3:
  # Selecciona no sólo el número de retardos sino
    # también el orden de cada uno
    # a partir del BIC y como explica Lütkephoul(1990)
  # Regresión auxiliar sin componentes deterministas

 selecPv3 <- function(vari, s, lmdet)
 {
   N <- length(vari)
   # N  <- length(na.omit(vari))
   bp   <- round(10*log10(N))
   ML   <- ret(vari, bp + 1)
   ndet <- ncol(model.matrix(lmdet))
   miss <- s-length(which(vari[1:s]!="NaN"))+1

   bic    <- c(1:2)
   lmarP  <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet]+
                                ML[miss:N,2:(bp+1)])
   bic[1] <- AICBIC(lmarP)[2]
   lmarP  <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet] +
                                ML[miss:N,2:(bp)])
   bic[2] <- AICBIC(lmarP)[2]

   k <- bp-1
   while(k > 0){
     if(bic[2] < bic[1]){
        bic[1] <- bic[2]
        ML[,(k+2)] <- rep(0, length(vari))
                        }
     else
        bic[1] <- bic[1]
     lmarP  <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet] +
                          ML[miss:N,2:k] + ML[miss:N,(k+2):(bp+1)])
     bic[2] <- AICBIC(lmarP)[2]
     k <- k-1
               }

   rsugaux <- matrix(nrow=1, ncol=bp)
   j <- 2
   while(j < bp+2){
     ifelse(ML[N,j] != 0, rsugaux[j-1] <- TRUE,
                          rsugaux[j-1] <- FALSE)
     j <- j+1
                  }
  rsug <- which(rsugaux == TRUE)
  rsug
}

#

# Selecciona los retardos quitando los que no son significativos al 10%
  # considerando un máximo de retardos posibles igual a 10*log10(N)
# Se supone N >100 y se compara con la Normal.

selecPv4 <- function(vari, s, lmdet)
{
   N <- length(vari)
   # N    <- length(na.omit(vari))
   bp   <- round(10*log10(N))
   ML   <- ret(vari, bp + 1)
   rsug <- rep(NA, bp)
   ndet <- ncol(model.matrix(lmdet))
   miss <- s-length(which(vari[1:s]!="NaN"))+1

   rsug <- c(2:(bp+1))
   aux <- 1
   while(aux==1)
   {
   aux <- 0
   lmarP <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet] +
                               ML[miss:N,na.omit(rsug)])
   j <- 1
   aux2 <- which(rsug > 0)
   for(i in (ndet+1):(ndet+length(na.omit(rsug))))
   {
     et <- contts(lmarP,i)[[2]]
     if(abs(et) <= 1.645){
             # Se supone N>100, se compara con la Normal
        rsug[aux2[j]] <- NA
        aux <- 1
                         }
     j <- j+1
   }
   }
   rsug <- na.omit(rsug)[1:length(na.omit(rsug))]-1
   rsug
}
# Los retardos sugeridos son los de orden rsug,
  # corregido del efecto ML

#

# Elige el numero de retardos Pmax en función del
  # criterio AIC y Ljung-Box
# Si la elección es cero retardos daría error,
  # no se puede calcular lm(ML[,1] ~ ML[,2:0])
# (arrg  <- ar(vari, aic=TRUE, method="ols"))

selecPv5 <- function(vari, s, lmdet)
{
  ABIC <- AICBICaux(vari, s, lmdet)
  arP  <- which.min(ABIC[,1])-1
  ML   <- ret(vari, nrow(ABIC))
  ndet <- ncol(model.matrix(lmdet))
  miss <- s-length(which(vari[1:s]!="NaN"))+1
  N <- length(vari)

  Pmax <- -1
  while(Pmax == -1){
    if(arP > 0)
    lmar <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet] + ML[miss:N,2:(arP+1)])
    if(arP == 0)
    lmar <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet])

    pvLB <- Box.test(lmar$residuals, lag = s, type="Ljung")[3]
    if(pvLB > 0.05)
      Pmax <- arP
    if(pvLB <= 0.05){
      ABIC[arP+1,1] <- Inf   # +1 porque el primero es cero retardos
      arP <- which.min(ABIC[,1])-1
                    }
                   }
  Pmax
}

#

# Elige el numero de retardos Pmax en función del
  # criterio BIC y Ljung-Box

selecPv6 <- function(vari, s, lmdet)
{
  ABIC <- AICBICaux(vari, s, lmdet)
  arP  <- which.min(ABIC[,2]) - 1
  ML   <- ret(vari, nrow(ABIC))
  ndet <- ncol(model.matrix(lmdet))
  miss <- s - length(which(vari[1:s]!="NaN")) + 1
  N <- length(vari)

  Pmax <- -1
  while(Pmax == -1){
    if(arP > 0)
    lmar <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet] + ML[miss:N,2:(arP+1)])
    if(arP == 0)
    lmar <- lm(ML[miss:N,1] ~ 0+model.matrix(lmdet)[,1:ndet])

    pvLB <- Box.test(lmar$residuals, lag = s, type="Ljung")[3]
    if(pvLB > 0.05)
      Pmax <- arP
    if(pvLB <= 0.05){
      ABIC[arP+1,2] <- Inf   # +1 porque el primero es cero retardos
      arP <- which.min(ABIC[,2])-1
                    }
                   }
  Pmax
}

#

rmg <- function(vari, krmg)
{
 # k es el número de periodos entre los que
   # se calcula el rango y la media

 N <- length(vari)
 media <- c(1:(N-krmg+1))
 rango <- c(1:(N-krmg+1))

 for(i in 1:length(media)){
   media[i]<- mean(vari[i:(i+(krmg-1))])
   rango[i]<- (max(vari[i:(i+(krmg-1))])-min(vari[i:(i+(krmg-1))]))
                          }
 cor.mr <- cor(media, rango, use = "complete")

 plot(media, rango, pch=20, main="Range-mean plot", xlab="Mean", ylab="Range")
 mtext(as.expression(substitute(cor(R,M)==cor.mr,
      list(cor.mr=round(cor.mr, 4)))), side = 3, line = 0.35, cex=0.7)

 cit1 <- " \n Correlation range-mean: "
 cit2 <- "\n"
 cat(cit1)
 cat(list=(round(cor.mr, 4)))
 cat(cit2)
}
#
corrgrm <- function(label, transf)
{
  #vari.env <- new.env(FALSE, NULL)
  #Fvarinfo(label$label, env.out=vari.env)
  #vari <- get("vari", env=vari.env)
  #s    <- get("s", env=vari.env)
  #t0   <- get("t0", env=vari.env)
  #N    <- length(vari)
  vari <- label$vari
  s    <- label$s
  t0   <- label$t0
  N    <- length(vari)

  if(transf == "original"){
#     opar <- par(mfrow=c(2,1), mar=c(4,4.7, 1 ,1.5), ps=18,
#                 font=1, tcl=-0.5, cex.lab=0.8, cex.axis=0.7)
     acf(vari, na.action = na.pass, lag=48, main="ACF", xlab="", ylab="")
     mtext("Original series", side=3, line=0.35, cex=0.7)
     pacf(vari, na.action = na.pass, lag=48, main="FACP", xlab="Lag", ylab="")
     #, xlab="Lag", ylab="PACF")
     mtext("Original series", side=3, line=0.35, cex=0.7)
#     par(opar)
                          }
  if(transf == "delta"){
#     opar <- par(mfrow=c(2,1), mar=c(4,4.7, 1 ,1.5), ps=18,
#                 font=1, tcl=-0.5, cex.lab=0.8, cex.axis=0.7)
     acf(diff(vari, lag=1), na.action = na.pass, lag=48, xlab="", ylab="", main="")
     mtext(expression(Delta(y)), side=3, line=0.35)
     pacf(diff(vari, lag=1), na.action = na.pass, lag=48, xlab="Lag", ylab="", main="")
     mtext(expression(Delta(y)), side=3, line=0.35)
#     par(opar)
                       }
  if(transf == "deltas"){
#     opar <- par(mfrow=c(2,1), mar=c(4,4.7, 1 ,1.5), ps=18,
#                 font=1, tcl=-0.5, cex.lab=0.8, cex.axis=0.7)
     acf(diff(vari, lag=s), na.action = na.pass, lag=48, xlab="", ylab="", main="")
     mtext(expression(Delta^s*(y)), side=3, line=0.35)
     pacf(diff(vari, lag=s), na.action = na.pass, lag=48, xlab="Lag", ylab="", main="")
     mtext(expression(Delta^s*(y)), side=3, line=0.35)
#     par(opar)
                        }
  if(transf == "deltadeltas"){
#     opar <- par(mfrow=c(2,1), mar=c(4,4.7,1,1.5), ps=18,
#                 font=1, tcl=-0.5, cex.lab=0.8, cex.axis=0.7)
     ML      <- ret(vari, s+2)
     ddsvari <- ML[,1]-ML[,2]-ML[,s+1]+ML[,s+2]
     acf(ddsvari, na.action = na.pass, lag=48, xlab="", ylab="", main="")
     mtext(expression(Delta*Delta^s* (y)), side=3, line=0.35)
     pacf(ddsvari, na.action = na.pass, lag=48, xlab="Lag", ylab="", main="")
     mtext(expression(Delta*Delta^s*(y)), side=3, line=0.35)
#     par(opar)
                             }
}
#

ExportGraph <- function()
{
    exgrfile <- tclvalue(tkgetSaveFile(filetypes='{"Text Files" {".ps" ".eps"}} {"All Files" {"*"}}'))
      if(!nchar(exgrfile))
         tkmessageBox(message="No file was chosen.", icon="error")
      else{
        dev.copy2eps(file=exgrfile, device=X11)
          }
      # tkdestroy(ttexgr)
      # tkdestroy(ttplot)   # ver con esto o sin esto
}

#

Transfdet <- function(vari, s)
{
  #t0   <- get("t0", env=vari.env)
  N    <- length(vari)

  VFEm <- function(vari, s)
  {
    N    <- length(vari)
    VFE <- matrix(0, nrow = N, ncol = s-1)
    sq1 <- seq(1, N, s)
    sq2 <- seq(0, N, s)
    k <- 0
    for (j in 1:(s-1)){
       ifelse(sq1[length(sq1)] + k > N, n1 <- length(sq1) -1, n1 <- N)
       ifelse(sq2[length(sq2)] + k > N, n2 <- length(sq2) -1, n2 <- N)
       for (i in 1:n1)
          VFE[sq1[i]+k, j] <- 1
       for (i in 1:n2)
          VFE[sq2[i]+k, j] <- -1
        k <- k + 1
                      }
    VFE
  }
  tiempo  <- matrix(c(1:length(vari)), ncol=1)
  VFE     <- VFEm(vari, s)
  ML      <- ret(vari, 2)
  lmdet   <- lm(ML[,1] ~ tiempo + VFE[,1:(s-1)])  # incluye constante

  VFE2 <- matrix(nrow=length(vari), ncol=s-1)
  for(j in 1:length(vari)){
     for(i in 3:(s+1))
       VFE2[j,(i-2)] <- lmdet$coef[i]*VFE[j,(i-2)]
                          }
   for(i in 1:length(vari))
      VFE2[i,1] <- sum(VFE2[i,])

   variodet <- vari-(lmdet$coef[1] + lmdet$coef[2]*tiempo + VFE2[,1])
   print(summary(lmdet))
   variodet
}

#MakeTransfdet <- function()
#{
#   string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
#   .wvari <<- ExeString(c(string, ""))
#   s <- .wvari$s
#   t0 <- .wvari$t0
#   t02 <- .wvari$t02
#   t0t <- .wvari$t0t <- t0t
#   vari2 <- .wvari$vari2
#   varit <- .wvari$varitvarit
#   N <- length(vari)
#
#   ifelse(.wvari$logvari==TRUE, label <<- "Series without deterministic component",
#                                label <<- "Series in logarithms without deterministic component")
#   varit <- ts(Transfdet(varit, s), frequency=s, start=t0t)
#   ref1 <- ysooys(t0, t02, length(vari2), s)[[1]]
#   ref2 <- ysooys(t0t, t02, length(vari2), s)[[1]]
#   if(ref1 < ref2){
#      t0 <- t0t
#      vari <- ts(varit[1:(N-(length(vari2)-length(varit)))], frequency=s, start=t0)}
#   if(ref1 >= ref2){
#      aux1 <- ysooys(t0, t0t, length(varit), s)[[1]]
#      aux2 <- ysooys(t0, t0, length(vari), s)[[2]][N,1:2]
#      aux3 <- ysooys(aux2, t0t, length(varit), s)[[1]]
#      vari <- ts(varit[aux1:aux3], frequency=s, start=t0)}
#   N <- length(vari)
#
   #change.wvari.label(c("label", "t0", "varit", "vari", "N"))
   #tkmessageBox(message="Original series has changed.", icon="info")
#}

#
TablaFrec <- function(outfile)
{

# ADF
  # Cambiar aquí selecP si se quiere
ADF1 <- ADF.test(.wvari, c(0,0,0), "biclb", 0, 0, showcat=FALSE)
  DF1    <- round(as.numeric(ADF1[1]), 2)
  P1adf  <- ADF1[3]
  adobs1 <- ADF1[4]
ADF2 <- ADF.test(.wvari, c(1,0,0), "biclb", 0, 0, showcat=FALSE)
  DF2    <- round(as.numeric(ADF2[1]), 2)
  P2adf  <- ADF2[3]
  adobs2 <- ADF2[4]
ADF3 <- ADF.test(.wvari, c(1,1,0), "biclb", 0, 0, showcat=FALSE)
  DF3    <- round(as.numeric(ADF3[1]), 2)
  P3adf  <- ADF3[3]
  adobs3 <- ADF3[4]
ADF4 <- ADF.test(.wvari, c(1,0,1), "biclb", 0, 0, showcat=FALSE)
  DF4    <- round(as.numeric(ADF4[1]), 2)
  P4adf  <- ADF4[3]
  adobs4 <- ADF4[4]
ADF5 <- ADF.test(.wvari, c(1,1,1), "biclb", 0, 0, showcat=FALSE)
  DF5    <- round(as.numeric(ADF5[1]), 2)
  P5adf  <- ADF5[3]
  adobs5 <- ADF5[4]

# KPSS

kpss0 <- as.numeric(KPSS.test(.wvari$vari, 0, showcat=FALSE))
kpss1 <- as.numeric(KPSS.test(.wvari$vari, 1, showcat=FALSE))
kpss2 <- as.numeric(KPSS.test(.wvari$vari, 2, showcat=FALSE))
kpss3 <- as.numeric(KPSS.test(.wvari$vari, 3, showcat=FALSE))
kpss4 <- as.numeric(KPSS.test(.wvari$vari, 4, showcat=FALSE))

# CH

rdoCH <- c(1:(.wvari$s/2+1))
frecf <- rep(0, .wvari$s/2)
for(i in 1:(.wvari$s/2))
{
  frecf[i] <- 1
  rdoCH[i] <- CH.test(.wvari, frecf, f0=1, DetTr=FALSE, showcat=FALSE)[1]
  frecf <- rep(0, .wvari$s/2)
}
rdoCH[(.wvari$s/2+1)] <- CH.test(.wvari, rep(1, .wvari$s/2), f0=1, DetTr=FALSE, showcat=FALSE)[1]
lCH <- CH.test(.wvari, rep(1, .wvari$s/2), f0=1, DetTr=FALSE, showcat=FALSE)[2]
CH  <- round(as.numeric(rdoCH), 2)
#
rdoCH <- c(1:(.wvari$s/2+1))
frecf <- rep(0, .wvari$s/2)
for(i in 1:(.wvari$s/2))
{
  frecf[i] <- 1
  rdoCH[i] <- CH.test(.wvari, frecf, f0=1, DetTr=TRUE, showcat=FALSE)[1]
  frecf <- rep(0, .wvari$s/2)
}
rdoCH[(.wvari$s/2+1)] <- CH.test(.wvari, rep(1, .wvari$s/2), f0=1, DetTr=TRUE, showcat=FALSE)[1]
CH2  <- round(as.numeric(rdoCH), 2)

# HEGYBM

HEGYBM1 <- HEGY.test(.wvari, c(0,0,0), "biclb", 0, 0, showcat=FALSE)
  H1    <- HEGYBM1[1]
  PH1   <- HEGYBM1[3]
  hobs1 <- HEGYBM1[4]
HEGYBM2 <- HEGY.test(.wvari, c(1,0,0), "biclb", 0, 0, showcat=FALSE)
  H2    <- HEGYBM2[1]
  PH2   <- HEGYBM2[3]
  hobs2 <- HEGYBM2[4]
HEGYBM3 <- HEGY.test(.wvari, c(1,1,0), "biclb", 0, 0, showcat=FALSE)
  H3    <- HEGYBM3[1]
  PH3   <- HEGYBM3[3]
  hobs3 <- HEGYBM3[4]
HEGYBM4 <- HEGY.test(.wvari, c(1,0,1), "biclb", 0, 0, showcat=FALSE)
  H4    <- HEGYBM4[1]
  PH4   <- HEGYBM4[3]
  hobs4 <- HEGYBM4[4]
HEGYBM5 <- HEGY.test(.wvari, c(1,1,1), "biclb", 0, 0, showcat=FALSE)
  H5    <- HEGYBM5[1]
  PH5   <- HEGYBM5[3]
  hobs5 <- HEGYBM5[4]

#
l49 <- "\\documentclass[12pt]{article}\n"
l50 <- "\\usepackage[ansinew]{inputenc}\n"
l51 <- "\\usepackage{a4}\n"
l52 <- "\\usepackage[spanish]{babel}\n"
l53 <- "\\selectlanguage{spanish}\n"
l54 <- "\\begin{document}\n"
#
l1 <- "\\begin{table}[h] \n"
l2 <- "\\centering \n"
l3 <- "\\caption{Integration versus stationary tests} \\label{Tfrec} \n"
l4 <- "\\begin{tabular}{rrrrrr} \n \\\\ \n"
l5 <- "& \\multicolumn{5}{c}{Deterministic components} \\\\ \n"
l6 <- "\\cline{2-6} \n"
l7 <-
  " & \\multicolumn{1}{c}{None} & \\multicolumn{1}{c}{I} & \n"
l8 <-
  "\\multicolumn{1}{c}{I+Tr} & \\multicolumn{1}{c}{I+SD} & \n"
l9 <- "\\multicolumn{1}{c}{I+Tr+SD} \\\\  \\hline \n"

l10 <-
c("\\hspace{1.1cm}($l=0$)&-& $\\eta_\\mu =",kpss0[1],"$& $\\eta_\\tau =",kpss0[2],"$ &-&- \\\\ \n")
l11 <-
c("\\hspace{1.1cm} ($l=1$)&-&$\\eta_\\mu =",kpss1[1],"$& $\\eta_\\tau =",kpss1[2],"$ &-&- \\\\ \n")
l12 <-
c("   KPSS       ($l=2$)&-&$\\eta_\\mu =",kpss2[1],"$ & $\\eta_\\tau =",kpss2[2],"$&-&- \\\\ \n")
l13 <-
c("\\hspace{1.1cm} ($l=3$) &-& $\\eta_\\mu =",kpss3[1],"$& $\\eta_\\tau =",kpss3[2],"$ &-&- \\\\ \n")
l14 <-
c("\\hspace{1.1cm} ($l=4$)&-& $\\eta_\\mu =",kpss4[1],"$& $\\eta_\\tau =",kpss4[2],"$ &-&- \\\\ \\hline \n")

l15 <- c(" & t=$",DF1[1],"$ & t=$",DF2[1],"$ & t=$",DF3[1],"$ & t=$",DF4[1],"$ & t=$",DF5[1],"$ \\\\ \n")
l16 <- "\\multicolumn{1}{c}{ \n"
l17 <- c("ADF }& $p=",P1adf,"$ & $p=",P2adf,"$ & $p=",P3adf,"$ & $p=",P4adf,"$ & $p=",P5adf,"$ \\\\ \n")
l17b <- c("& n.obs=$",adobs1,"$ & n.obs=$",adobs2,"$ & n.obs=$",adobs3,"$ & n.obs=$",adobs4,"$ & n.obs=$",adobs5,"$ \\\\ \\hline \n")

if(.wvari$s==12){
l18 <- c("& & & & \\multicolumn{1}{l}{$L_{\\pi/6}\\; = ",CH[1],"$} &$", CH2[1],"$ \\\\ \n")
l19 <- c("& & & & \\multicolumn{1}{l}{$L_{\\pi/3} \\; = ",CH[2],"$}&$", CH2[2],"$ \\\\ \n")
         }
l20 <- "\\multicolumn{1}{c}{ \n"
l21 <- c(" CH ($l=",lCH,"$)} &-&-&-& \\multicolumn{1}{l}{$L_{\\pi/2} \\;= ",CH[.wvari$s/4],"$}&$", CH2[.wvari$s/4],"$ \\\\ \n")
if(.wvari$s==12){
l22 <- c("& & & & \\multicolumn{1}{l}{$L_{2\\pi/3} = ",CH[4],"$}&$", CH2[4],"$ \\\\ \n")
l23 <- c("& & & & \\multicolumn{1}{l}{$L_{5\\pi/6} = ",CH[5],"$} &$", CH2[5],"$ \\\\ \n")
         }
l24 <- c("& & & & \\multicolumn{1}{l}{$L_{\\pi} \\;\\;\\;\\, = ",CH[.wvari$s/2],"$} &$", CH2[.wvari$s/2],"$ \\\\ \n")
l25 <- c("& & & & \\multicolumn{1}{l}{$L_{f} \\;\\;\\;\\, =   ",CH[.wvari$s/2+1],"$} &$", CH2[.wvari$s/2+1],"$ \\\\ \\hline \n")

l26a <- c("n.obs & $",hobs1,"$ & $",hobs2,"$ & $",hobs3,"$ & $",hobs4,"$ & $",hobs5,"$ \\\\ \n")
l35 <- "\\multicolumn{1}{c}{ HEGY \n"
l26 <- c("$ p \\; = $} & $",PH1,"$ & $",PH2,"$ & $",PH3,"$ & $",PH4,"$ & $",PH5,"$ \\\\ \n")
l27 <- c("$ t_1 = $ & $",H1[[1]][[1]][1],"$ & $",H2[[1]][[1]][1],"$ &  $",H3[[1]][[1]][1]," $ & $",H4[[1]][[1]][1]," $ & $",H5[[1]][[1]][1]," $ \\\\ \n")
l28 <- c("$ t_2 = $ & $",H1[[1]][[1]][2],"$ & $",H2[[1]][[1]][2]," $ & $",H3[[1]][[1]][2]," $ & $",H4[[1]][[1]][2]," $ & $",H5[[1]][[1]][2]," $ \\\\ \n")
l29 <- c("$ t_3 = $ & $",H1[[1]][[1]][3],"$ &  $",H2[[1]][[1]][3],"$ &  $",H3[[1]][[1]][3]," $ & $",H4[[1]][[1]][3],"$  & $",H5[[1]][[1]][3]," $ \\\\ \n")
l30 <- c("$ t_4 = $ & $",H1[[1]][[1]][4],"$ & $",H2[[1]][[1]][4]," $ & $",H3[[1]][[1]][4]," $ & $",H4[[1]][[1]][4]," $ & $",H5[[1]][[1]][4]," $ \\\\ \n")

if(.wvari$s == 12){
l31 <- c("$ t_5 = $ & $",H1[[1]][[1]][5],"$ & $",H2[[1]][[1]][5]," $ & $",H3[[1]][[1]][5]," $ & $",H4[[1]][[1]][5]," $ & $",H5[[1]][[1]][5]," $ \\\\ \n")
l32 <- c("$ t_6 = $ & $",H1[[1]][[1]][6],"$ & $",H2[[1]][[1]][6]," $ & $",H3[[1]][[1]][6]," $ & $",H4[[1]][[1]][6]," $ & $",H5[[1]][[1]][6]," $ \\\\ \n")
l33 <- c("$ t_7 = $ & $",H1[[1]][[1]][7],"$ & $",H2[[1]][[1]][7]," $ & $",H3[[1]][[1]][7]," $ & $",H4[[1]][[1]][7]," $ & $",H5[[1]][[1]][7]," $ \\\\ \n")
l34 <- c("$ t_8 = $ & $",H1[[1]][[1]][8],"$ & $",H2[[1]][[1]][8]," $ & $",H3[[1]][[1]][8]," $ & $",H4[[1]][[1]][8]," $ & $",H5[[1]][[1]][8]," $ \\\\ \n")
           }
if(.wvari$s == 12){
l36 <- c("$\\;\\;\\: t_9 = $ &$",H1[[1]][[1]][9],"$ & $",H2[[1]][[1]][9],"$ & $",H3[[1]][[1]][9]," $ & $",H4[[1]][[1]][9]," $ & $",H5[[1]][[1]][9]," $ \\\\ \n")
l37 <- c("$ t_{10} = $ & $",H1[[1]][[1]][10],"$ & $",H2[[1]][[1]][10]," $ & $",H3[[1]][[1]][10]," $ & $",H4[[1]][[1]][10]," $ & $",H5[[1]][[1]][10]," $ \\\\ \n")
l38 <- c("$ t_{11} = $ & $",H1[[1]][[1]][11],"$ & $",H2[[1]][[1]][11]," $ & $",H3[[1]][[1]][11]," $ & $",H4[[1]][[1]][11]," $ & $",H5[[1]][[1]][11]," $ \\\\ \n")
l39 <- c("$ t_{12} = $ & $",H1[[1]][[1]][12],"$ & $",H2[[1]][[1]][12]," $ & $",H3[[1]][[1]][12]," $ & $",H4[[1]][[1]][12]," $ & $",H5[[1]][[1]][12]," $ \\\\ \n")
           }
l40 <- c("$ F_{34}\\;\\;= $ & $",H1[[1]][[2]][1],"$ & $",H2[[1]][[2]][1]," $ & $",H3[[1]][[2]][1]," $ & $",H4[[1]][[2]][1]," $ & $",H5[[1]][[2]][1]," $ \\\\ \n")

if(.wvari$s == 12){
l41 <- c("$ F_{56}\\;\\,\\, = $ & $",H1[[1]][[2]][2],"$ & $",H2[[1]][[2]][2]," $ & $",H3[[1]][[2]][2]," $ & $",H4[[1]][[2]][2]," $ & $",H5[[1]][[2]][2]," $ \\\\ \n")
l42 <- c("$ F_{78}\\;\\;\\, = $ & $",H1[[1]][[2]][3],"$ & $",H2[[1]][[2]][3]," $ & $",H3[[1]][[2]][3]," $ & $",H4[[1]][[2]][3]," $ & $",H5[[1]][[2]][3]," $ \\\\ \n")
l43 <- c("$ F_{9\\,10}\\;  = $ & $",H1[[1]][[2]][4],"$ & $",H2[[1]][[2]][4]," $ & $",H3[[1]][[2]][4]," $ & $",H4[[1]][[2]][4]," $ & $",H5[[1]][[2]][4]," $ \\\\ \n")
l44 <- c("$ F_{11\\,12} = $  & $",H1[[1]][[2]][5],"$ & $",H2[[1]][[2]][5]," $ & $",H3[[1]][[2]][5]," $ & $",H4[[1]][[2]][5]," $ & $",H5[[1]][[2]][5]," $ \\\\ \n")
           }
l45 <- c("$ F_{2...",.wvari$s,"} = $  & $",H1[[1]][[2]][.wvari$s/2],"$ & $",H2[[1]][[2]][.wvari$s/2]," $ & $",H3[[1]][[2]][.wvari$s/2]," $ & $",H4[[1]][[2]][.wvari$s/2]," $ & $",H5[[1]][[2]][.wvari$s/2]," $ \\\\ \n")
l46 <- c("$ F_{1...",.wvari$s,"} = $  & $",H1[[1]][[2]][.wvari$s/2+1],"$ & $",H2[[1]][[2]][.wvari$s/2+1]," $ & $",H3[[1]][[2]][.wvari$s/2+1]," $ & $",H4[[1]][[2]][.wvari$s/2+1]," $ & $",H5[[1]][[2]][.wvari$s/2+1]," $ \\\\ \n")
l47 <- "\\hline \n \\multicolumn{6}{l}{\\scriptsize{I: Intercept; Tr: Trend; SD; Seasonal dummys.}} \n"
l48 <- "\\end{tabular} \n \\end{table} \n"
#
l56 <- "\\end{document} \n"
#
if(.wvari$s==12){
com1 <- c(l49,l50,l51,l52,l53,l54,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,
l13,l14,l15,l16,l17,l17b,l18,l19,l20,l21,l22,l23,l24,l25,l26a,l35,l26,
l27,l28,l29,l30,l31,l32,l33,l34,l36,l37,l38,l39,l40,l41,l42,l43,
l44,l45,l46,l47,l48,l56)
         }
if(.wvari$s==4){
com1 <- c(l49,l50,l51,l52,l53,l54,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,
l14,l15,l16,l17,l17b,l20,l21,l24,l25,l26a,l35,l26,l27,
l28,l29,l30,l40,l45,l46,l47,l48,l56)
         }
cat(as.character(com1), file=outfile, append=FALSE)
}

#

TablaDet <- function(outfile)
{
   AR  <- matrix(nrow=.wvari$s+1, ncol=3)
   tiempo <- matrix(c(1:.wvari$N), ncol=1)
   VFE    <- MVFE(.wvari$vari, .wvari$s, .wvari$t0, "alg")
   ML     <- ret(.wvari$vari, round(10*log10(.wvari$N))+1)
   lmdet  <- lm(ML[,1] ~ tiempo + VFE[,1:(.wvari$s-1)])      # incluye constante
   Pmax   <- selecPv6(.wvari$vari, .wvari$s, lmdet)          # método bic-lb
   lmardet <- lm(ML[,1] ~ tiempo + VFE[,1:(.wvari$s-1)] + ML[,2:(Pmax+1)])
      # para métodos con rsug (bic-lut) poner ML[,] como en lmadf
   for(i in 1:(.wvari$s+1)){
     et <- contts(lmardet, i)[[2]]
     AR[i,1] <- round(lmardet$coef[i], 2)
     AR[i,2] <- round(et, 2)
                    }
   if(.wvari$s==4) {Fvfeaux <- c(3,4,5)}
   if(.wvari$s==12){Fvfeaux <- c(3,4,5,6,7,8,9,10,11,12,13)}
   AR[1,3] <- round(Fsnd(lmardet, .wvari$s-1, .wvari$s+1+Pmax, Fvfeaux), 2)

   arobs <- .wvari$N - max(Pmax)
   arP   <- Pmax

 ADF5 <- ADF.test(.wvari, c(1,1,1), "biclb", 0, 0, showcat=FALSE)
   DF5   <- ADF5[2]
   P5adf <- as.numeric(ADF5[3])
   adobs <- ADF5[[4]]

 HEGYBM5 <- HEGY.test(.wvari, c(1,1,1), "biclb", 0, 0, showcat=FALSE)
   H5   <- HEGYBM5[2]
   PH5  <- as.numeric(HEGYBM5[3])
   hobs <- HEGYBM5[[4]]
#

l26 <- "\\documentclass[12pt]{article} \n"
l27 <- "\\usepackage[ansinew]{inputenc} \n"
l28 <- "\\usepackage{a4} \n"
l29 <- "\\usepackage[spanish]{babel} \n"
l30 <- "\\selectlanguage{spanish} \n"
l31 <- "\\begin{document} \n"

l1  <- "\\begin{table}[h] \n"
l2  <- "\\centering  \\caption{Deterministic components}  \\label{Tdet} \n"
l3  <- "\\begin{tabular}{l|rr|rr|rr|} \n"
l4  <- "& \\multicolumn{2}{c|}{AR($y_t$)} & \\multicolumn{2}{c|}{ADF} & \\multicolumn{2}{|c|}{HEGYBM} \\\\ \n"
l5  <- "\\cline{2-7} \n"
l6  <- "\\multicolumn{1}{l|}{} &  Coeff & t-stat & Coeff & t-stat & Coeff & t-stat \\\\ \n"
l7  <- "\\hline \n"
l8  <- c("\\multicolumn{1}{|l|}{Intercept} & $",AR[1,1],"$&$",AR[1,2],"$&$",DF5[[1]][1,1],"$ & $",DF5[[1]][1,2],"$ & $",H5[[1]][1,1],"$ & $",H5[[1]][1,2],"$ \\\\ \n")
l9  <- c("\\multicolumn{1}{|l|}{Trend} & $",AR[2,1],"$&$",AR[2,2],"$&$",DF5[[1]][2,1],"$ & $",DF5[[1]][2,2],"$ & $",H5[[1]][2,1],"$ & $",H5[[1]][2,2],"$ \\\\ \n")
l10 <- c("\\multicolumn{1}{|l|}{SD1}      &$",AR[3,1],"$&$",AR[3,2],"$& $",DF5[[1]][3,1],"$ & $",DF5[[1]][3,2],"$ & $",H5[[1]][3,1],"$ & $",H5[[1]][3,2],"$ \\\\ \n")
l11 <- c("\\multicolumn{1}{|l|}{SD2}   & $",AR[4,1],"$&$",AR[4,2],"$&$",DF5[[1]][4,1],"$ & $",DF5[[1]][4,2],"$ & $",H5[[1]][4,1],"$ & $",H5[[1]][4,2],"$ \\\\ \n")
l12 <- c("\\multicolumn{1}{|l|}{SD3}   &$",AR[5,1],"$&$",AR[5,2],"$& $",DF5[[1]][5,1],"$ & $",DF5[[1]][5,2],"$ & $",H5[[1]][5,1],"$ & $",H5[[1]][5,2],"$ \\\\ \n")

if(.wvari$s == 12){
l13 <- c("\\multicolumn{1}{|l|}{SD4}  & $",AR[6,1],"$&$",AR[6,2],"$&$",DF5[[1]][6,1],"$ & $",DF5[[1]][6,2],"$ & $",H5[[1]][6,1],"$ & $",H5[[1]][6,2],"$ \\\\ \n")
l14 <- c("\\multicolumn{1}{|l|}{SD5}  &$",AR[7,1],"$&$",AR[7,2],"$& $",DF5[[1]][7,1],"$ & $",DF5[[1]][7,2],"$ & $",H5[[1]][7,1],"$ & $",H5[[1]][7,2],"$ \\\\ \n")
l14b <- c("\\multicolumn{1}{|l|}{SD6}  & $",AR[8,1],"$&$",AR[8,2],"$&$",DF5[[1]][8,1],"$ & $",DF5[[1]][8,2],"$ & $",H5[[1]][8,1],"$ & $",H5[[1]][8,2],"$ \\\\ \n")
l15 <- c("\\multicolumn{1}{|l|}{SD7}  & $",AR[9,1],"$&$",AR[9,2],"$&$",DF5[[1]][8,1],"$ & $",DF5[[1]][9,2],"$ & $",H5[[1]][9,1],"$ & $",H5[[1]][9,2],"$ \\\\ \n")
l16 <- c("\\multicolumn{1}{|l|}{SD8}  & $",AR[10,1],"$&$",AR[10,2],"$&$",DF5[[1]][10,1],"$ & $",DF5[[1]][10,2],"$ & $",H5[[1]][10,1],"$ & $",H5[[1]][10,2],"$ \\\\ \n")
l17 <- c("\\multicolumn{1}{|l|}{SD9}  & $",AR[11,1],"$&$",AR[11,2],"$&$",DF5[[1]][11,1],"$ & $",DF5[[1]][11,2],"$ & $",H5[[1]][11,1],"$ & $",H5[[1]][11,2],"$ \\\\ \n")
l18 <- c("\\multicolumn{1}{|l|}{SD10} & $",AR[12,1],"$&$",AR[12,2],"$&$",DF5[[1]][12,1],"$ & $",DF5[[1]][12,2],"$ & $",H5[[1]][12,1],"$ & $",H5[[1]][12,2],"$ \\\\ \n")
l19 <- c("\\multicolumn{1}{|l|}{SD11} & $",AR[13,1],"$&$",AR[13,2],"$&$",DF5[[1]][13,1],"$ & $",DF5[[1]][13,2],"$ & $",H5[[1]][13,1],"$ & $",H5[[1]][13,2],"$ \\\\ \n")
          }

l33 <- c("\\multicolumn{1}{|l|}{$F_{SD}$} &\\multicolumn{2}{|c}{$",AR[1,3],"$}  &\\multicolumn{2}{|c}{$",DF5[[1]][1,3],"$}  &\\multicolumn{2}{|c|}{$",H5[[1]][1,3],"$} \\\\ \n")

l20 <- "\\hline \n"
l21 <- c("\\multicolumn{1}{|l}{Lags} & \\multicolumn{2}{|c}{$",arP,"$} & \\multicolumn{2}{|c}{$",P5adf,"$} & \\multicolumn{2}{|c|}{$",PH5,"$} \\\\ \n")
l22 <- c("\\multicolumn{1}{|l}{n.obs} & \\multicolumn{2}{|c|}{$",arobs,"$}& \\multicolumn{2}{|c|}{$",adobs,"$} & \\multicolumn{2}{|c|}{$",hobs,"$} \\\\ \n")
l23 <- "\\hline \n"
l24 <- "\\end{tabular} \n"
l25 <- "\\end{table} \n"

l32 <- "\\end{document}\n"

if(.wvari$s==4)
com1 <- c(l26,l27,l28,l29,l30,l31,l1,l2,l3,l5,l4,l5,l6,l7,l8,l9,l10,l11,
l12,l33,l20,l21,l22,l23,l24,l25,l32)
if(.wvari$s==12)
com1 <- c(l26,l27,l28,l29,l30,l31,l1,l2,l3,l5,l4,l5,l6,l7,l8,l9,l10,l11,
l12,l13,l14,l14b,l15,l16,l17,l18,l19,l33,l20,l21,l22,l23,l24,l25,l32)

cat(com1, file=outfile, append=FALSE)
}

#

PanelqmCorrg <- function()
{
#   opar <- par(mfrow=c(4,2), mar=c(3,4,3,3.5), # mar=c(3,3,3.5,2),
#               tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
   corrgrm(.wvari, "original")
   corrgrm(.wvari, "delta")
   corrgrm(.wvari, "deltas")
   corrgrm(.wvari, "deltadeltas")
#   par(opar)
}
#

PanelmFreq <- function()
{
#   opar <- par(mfrow=c(4,2), mar=c(3,3,3.5,2),
#               tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
   freqg(.wvari)
#   par(opar)
}

#

PanelmRMBB <- function()
{
#   opar <- par(mfrow=c(3,2), mar=c(4,4,3.5,3),
#               tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
   rmg(.wvari$vari, .wvari$s)
   bbap(.wvari$vari, .wvari$s, .wvari$t0,
         c(.wvari$t0[1], .wvari$t0[1]+2, .wvari$t0[1]+4, .wvari$t0[1]+6))
     # Serie de al menos 8 anyos
   bbmp(.wvari$vari, .wvari$s, .wvari$t0, c(1:s), "Prot")
     # mplot <- rep(1, s)
     # which(mplot==1)
#   par(opar)
}

PanelqRMBBFreq <- function()
{
#  opar <- par(mfrow=c(4,2), mar=c(3,4,3,3),
#              tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
  rmg(.wvari$vari, .wvari$s)
  freqg(.wvari)
   bbap(.wvari$vari, .wvari$s, .wvari$t0, c(.wvari$t0[1], .wvari$t0[1]+2, .wvari$t0[1]+4, .wvari$t0[1]+6))
     # Serie de al menos 8 anyos
   quarterg(.wvari$vari, .wvari$s, .wvari$t0, plot=TRUE)
#  par(opar)
}

perdiff <- function(label)
{
  #vari.env <- new.env(FALSE, NULL)
  #Fvarinfo(label$label, env.out=vari.env)
  #vari  <- get("vari", env=vari.env)
  #vari2 <- get("vari2", env=vari.env)
  #varit <- get("varit", env=vari.env)
  #s     <- get("s", env=vari.env)
  #t0    <- get("t0", env=vari.env)
  #t02   <- get("t02", env=vari.env)
  #t0t   <- get("t0t", env=vari.env)
  #N     <- length(vari)
  vari <- label$vari
  s    <- label$s
  t0   <- label$t0
  N    <- length(vari)

  difpvari <- rep(NA, N)
  ML       <- ret(vari, 2)
  VFE      <- MVFE(vari, s, t0, "alg")
  Y1       <- ML[,2]*VFE

  if(s==4){
    nlscoef <- coef(nls(ML[,1] ~ 0+v1*VFE[,1]+v2*VFE[,2]+v3*VFE[,3]+v4*VFE[,4] +
               a1*Y1[,1] + a2*Y1[,2] + a3*Y1[,3] + (1/(a1*a2*a3))*Y1[,4],
               start=list(a1=0.5, a2=0.5, a3=0.5,
                          v1=1,v2=1,v3=1,v4=1), alg="plinear", trace=FALSE))
    alphas <- c(nlscoef[1:(s-1)], 1/prod(nlscoef[1:(s-1)]))
          }
  if(s==12){
    nlscoef <- coef(nls(ML[,1] ~ 0+v1*VFE[,1]+v2*VFE[,2]+v3*VFE[,3]+
               v4*VFE[,4] + v5*VFE[,5] + v6*VFE[,6] + v7*VFE[,7] +
               v8*VFE[,8] + v9*VFE[,9] + v10*VFE[,10] + v11*VFE[,11] + v12*VFE[,12] +
               a1*Y1[,1]+ a2*Y1[,2]+ a3*Y1[,3] + a4*Y1[,4]+ a5*Y1[,5]+ a6*Y1[,6] +
               a7*Y1[,7]+ a8*Y1[,8]+ a9*Y1[,9] + a10*Y1[,10]+ a11*Y1[,11] +
               (1/(a1*a2*a3*a4*a5*a6*a7*a8*a9*a10*a11))*Y1[,12],
               start=list(a1=0.7,a2=0.7,a3=0.7,a4=1,a5=1,a6=1,a7=0.3,
                          a8=0.5,a9=0.5,a10=0.8,a11=0.8,
                          v1=1,v2=1,v3=1,v4=1,v5=1,v6=1,v7=1,v8=1,v9=1,v10=1,v11=1,v12=1),
                          alg="plinear", trace=FALSE))
    alphas <- c(nlscoef[1:(s-1)], 1/prod(nlscoef[1:(s-1)]))
          }

  ifelse(t0[2] == s, seas <- 1, seas <- t0[2]+1)
  for(i in 2:N){
    difpvari[i] <- vari[i] - alphas[seas]*vari[i-1]
    ifelse(seas == s, seas <- 1, seas <- seas+1)
               }

  difpvari <- ts(difpvari, frequency=s, start=t0)        # usado para gráficos
  # difpvari2 <- ts(difpvari[2:N], frequency=s, start=t0)  # usar con la opción urootmenu=TRUE

  #  if(s==4)
  #    quarterg(difpvari, .wvari$s, .wvari$t0, plot=TRUE)
  #  if(s==12)
  #    bbmp(difpvari, .wvari$s, .wvari$t0, c(1:12), "Prot", plot=TRUE)

  # cambiar variables internas: t0t, varit,...
    # sólo usado por uroot cuando para transformar la serie con la que se está trabajando.

  # print(summary(nlscoef))
  difpvari
}

#
#Fvarinfo <- function(label, env.out)
#{
#  tempf <- tempfile()
#  cat(label, file=tempf)
#  assign("varinfo", eval(source(tempf))[[1]])
#  tkcmd("file", "delete", tempf) #; rm(tempf)
#
#  info <- c("vari", "s", "t0", "N", "logvari", "vari2", "varit", "t02", "t0cp", "t0t", "label")
#  for(i in 1:length(varinfo))
#    assign(info[i], varinfo[[i]], env=env.out)
#}

# Ejecutar como commando una variable de tipo carácter
ExeString <- function(Char)
{
  tempf <- tempfile()
  string <- paste(Char[1], Char[2], sep="")
  if(length(Char) > 2){
    for(i in 3:length(Char))
      string <- paste(string, Char[i], sep="")
  }
  string <- cat(string, file=tempf)
  # cat(paste(char, arg, sep=""), file=tempf)
  char.out <- eval(source(tempf))[[1]]

  char.out
}


# input: Objects <- list(a=a, b=b)
clean.auxobjects <- function(Objects, type)
{
  tt <- tktoplevel()
  tkwm.title(tt, "Info request")

  message <- "Please, press OK when results appear on the screen."
  tkgrid(tklabel(tt, text=message))

  OnOK <- function()
  {
    tkdestroy(tt)
    if(type == "ExString")
      ExeString(c(Objects, ""))
    if(type != "ExString")
    {
      Obaux <- names(Objects)
      aux <- Obaux[1]
      if(length(Objects) > 1)
      for(i in 2:length(Obaux))
        aux <- paste(aux, Obaux[i], sep=",")

      string <- paste("remove(", aux[1], ", envir=.GlobalEnv)", sep="")
      ExeString(c(string, ""))
    }
  }
  OK.but <- tkbutton(tt, text="OK", command=OnOK)
  tkgrid(OK.but)
  tkfocus(tt)
}

variinfo_treeW <- function()
{
  tt  <- tktoplevel()
  tkwm.title(tt, "Series info")
  txt <- tktext(tt, bg="white", height=30, width=70)
  tkgrid(txt)

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  # tclvalue(.label) <- as.character(tclvalue(string))
  #.wvari <<- ExeString(c(tclvalue(.label), ""))
  .wvari <<- ExeString(c(string, ""))

  if(.wvari$s == 4){ per <- "4, quarterly" }
  if(.wvari$s == 12){ per <- "12, monthly" }
  end <- ysooys(.wvari$N, .wvari$t0, .wvari$N, .wvari$s)[[1]]

  tkinsert(txt, "end", paste("\n", .wvari$label, "\n"))
  start <- paste(.wvari$t0[1], .wvari$t0[2], sep=".")
  end  <- paste(end[1], end[2], sep=".")
  tkinsert(txt, "end", paste("\n Period:", start, " - ", end, "\n"))
  tkinsert(txt, "end", paste("\n Periodicity: ", per, "\n"))

  tkconfigure(txt, state="disabled")
  tkfocus(txt)
}

#

# FUCIONES GUI

# MENU ARCHIVO

OpenSourceFile <- function()
{
   sourcefile <- tclvalue(
    tkgetOpenFile(filetypes='{"R Code and text files" {".txt" ".TXT" ".R"}} {"All Files" {"*"}}'))
    if(!nchar(sourcefile))
       tkmessageBox(message="No file was chosen.", icon="error")
    else
       source(sourcefile)
}

SaveAsWorkSpace <- function()
{
  savename <- tclvalue(tkgetSaveFile(filetypes='{"R files" {".RData"}} {"All Files" {"*"}}'))
  if(!nchar(savename))
    tkmessageBox(message="No file was chosen.", icon="error")
  else
    save(file=savename, list=ls(all=TRUE))
}

LoadWorkSpace <- function()
{
   wsfile <- tclvalue(tkgetOpenFile(filetypes='{"R files" {".RData"}} {"All Files" {"*"}}'))
   if(!nchar(wsfile))
      tkmessageBox(message="No file was chosen.", icon="error")
   else
      load(wsfile)
}

closerusea <- function()
{
    response <- tclvalue(tkmessageBox(message="Quit uroot?",
        icon="question", type="okcancel", default="cancel"))
    if (response == "cancel") return(invisible(response))
    else{
      tkdestroy(.tt)
      rm(.tt, .treeWidget, envir=.GlobalEnv)
        }
}

# MENU AYUDA

# MENU CONTRATES

MakeADF.test <- function()
{
  ttadf  <- tktoplevel()
  tkwm.title(ttadf, "ADF test")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  N <- .wvari$N

  tkgrid(tklabel(ttadf, text="  Select deterministic components:", fg="blue"), sticky="w")
  Cte  <- tkcheckbutton(ttadf)
  Tdl  <- tkcheckbutton(ttadf)
  Vfet <- tkcheckbutton(ttadf)

  CteValue  <- tclVar("0")
  TdlValue  <- tclVar("0")
  VfetValue <- tclVar("0")

  tkconfigure(Cte, variable=CteValue)
    tkgrid(tklabel(ttadf, text="Intercept"), Cte)
  tkconfigure(Tdl, variable=TdlValue)
    tkgrid(tklabel(ttadf, text="Trend"), Tdl)
  tkconfigure(Vfet, variable=VfetValue)
    tkgrid(tklabel(ttadf, text="Seasonal dummys"), Vfet)

  rb1 <- tkradiobutton(ttadf)
  rb2 <- tkradiobutton(ttadf)
  rb3 <- tkradiobutton(ttadf)
  rb4 <- tkradiobutton(ttadf)
  rb5 <- tkradiobutton(ttadf)
  rb6 <- tkradiobutton(ttadf)
  rbValue <- tclVar("BIC-LB")
  yself <- tclVar()
  entry.yself <- tkentry(ttadf, width="3", textvariable=yself)

  tkconfigure(rb1, variable=rbValue, value="AIC-LB")
  tkconfigure(rb2, variable=rbValue, value="BIC-LB")
  tkconfigure(rb3, variable=rbValue, value="AIC-Lüt")
  tkconfigure(rb4, variable=rbValue, value="BIC-Lüt")
  tkconfigure(rb5, variable=rbValue, value="Signf")
  tkconfigure(rb6, variable=rbValue, value="Tu mismo")

  tkgrid(tklabel(ttadf, text="  Select the method for choosing lags:", fg="blue"), sticky="w")
  tkgrid(tklabel(ttadf, text="AIC-LB"), rb1)
  tkgrid(tklabel(ttadf, text="BIC-LB"), rb2)
  tkgrid(tklabel(ttadf, text="AIC-top-down"), rb3)
  tkgrid(tklabel(ttadf, text="BIC-top-down"), rb4)
  tkgrid(tklabel(ttadf, text="Significant lags"), rb5)
  tkgrid(tklabel(ttadf, text="By yourself"), rb6, entry.yself)

  Definir1 <- tkbutton(ttadf, text="Define", command=Makevfic)
  Definir2 <- tkbutton(ttadf, text="Define", command=MakeVFEp)
  Mvfic <<- matrix(NA, nrow=N, ncol=6)
  tkgrid(tklabel(ttadf, text="  Include dummys:", fg="blue"), sticky="w")
  tkgrid(tklabel(ttadf, text="  Generic dummy           "), Definir1)
  #VFEp <<- 0; done <<- tclVar(0)
  tkgrid(tklabel(ttadf, text="Partial seasonal dummy      "), Definir2)

  tkconfigure(.treeWidget, cursor="watch")
  OnOK <- function()
  {
      tkdestroy(ttadf)
      ifelse(exists("VFEp"), VFEp <- VFEp, VFEp <<- 0)

      mVal <- as.character(tclvalue(CteValue))
      if(mVal =="1")
         Ct <- 1
      if(mVal =="0")
         Ct <- 0
      mVal <- as.character(tclvalue(TdlValue))
      if(mVal =="1")
         TD <- 1
      if(mVal =="0")
         TD <- 0
      mVal <- as.character(tclvalue(VfetValue))
      if(mVal =="1")
         Vfe <- 1
      if(mVal =="0")
         Vfe <- 0

      # SelecP
      rbVal <- as.character(tclvalue(rbValue))
      if(tclvalue(rbValue) == "AIC-LB")
        selecP <- "aiclb"
      if(tclvalue(rbValue) == "BIC-LB")
        selecP <- "biclb"
      if(tclvalue(rbValue) == "AIC-Lüt")
        selecP <- "aiclut"
       if(tclvalue(rbValue) == "BIC-Lüt")
        selecP <- "biclut"
      if(tclvalue(rbValue) == "Signf")
        selecP <- "signf"
      if(tclvalue(rbValue) == "Tu mismo")
        selecP <- as.numeric(tclvalue(yself))

      # VFIC
      aux <- length(which(Mvfic[1,] >= 0))
      ifelse(aux == 0, Mvfic <<- 0, Mvfic <<- as.matrix(Mvfic[,1:aux]))
      rdoADF <- ADF.test(.wvari, c(Ct,TD,Vfe), selecP, Mvfic, VFEp, showcat=TRUE)

      ###
      # Poner ventana con botón OK, CUANDO SE PULSA DONE <- 1
      #tkwait.variable(done)
      #if(doneVal == 1) rm(Mvfic, VFEp, done)
      #clean.auxobjects()
      #tkwait.variable(done)
      #if(tclvalue(done) == 1) rm(Mvfic, VFEp, done)
      cleanlist <- list(Mvfic=Mvfic, VFEp=VFEp)
      clean.auxobjects(cleanlist, type="")
      tkconfigure(.treeWidget, cursor="xterm")
  }
  OK.but <- tkbutton(ttadf,text="OK",command=OnOK)
  tkgrid(OK.but)
}

#

MakeKPSS.test <- function()
{
  ttkpss <- tktoplevel()
  tkwm.title(ttkpss, "KPSS test")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  ltrunc <- tclVar(as.integer(3*sqrt(.wvari$N)/13))
  tkgrid(tklabel(ttkpss, text="  Introduce the lag truncation parameter: \n (By default 3*sqrt(N)/13)",
       fg="blue"), rowspan=2)
  entry.ltrunc <- tkentry(ttkpss, width="5", textvariable=ltrunc)
  tkgrid(entry.ltrunc)

  tkconfigure(.treeWidget, cursor="watch")
  OnOK <- function()
  {
     tkdestroy(ttkpss)
     rdoKPSS <- KPSS.test(.wvari$vari, as.numeric(tclvalue(ltrunc)), showcat=TRUE)
     rdoKPSS <- rdoKPSS[[1]]
     tkconfigure(.treeWidget, cursor="xterm")
  }
  OK.but <- tkbutton(ttkpss, text="OK", command=OnOK)
  tkgrid(OK.but)
}

#

MakeCH.test <- function()
{
  ttch <- tktoplevel()
  tkwm.title(ttch, "CH test")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  trend  <- tkcheckbutton(ttch);  trendValue  <- tclVar("0")
  lag1   <- tkcheckbutton(ttch);  lag1Value  <- tclVar("1")

  tkgrid(tklabel(ttch, text="  Select the elements to include in the auxiliar regression:",
        fg="blue"), sticky="w")
  tkconfigure(trend, variable=trendValue)
    tkgrid(tklabel(ttch, text="Trend"), trend)
  tkconfigure(lag1, variable=lag1Value)
    tkgrid(tklabel(ttch, text="First order lag"), lag1)

  rb1 <- tkradiobutton(ttch)
  rb2 <- tkradiobutton(ttch)
  rbValue <- tclVar("nsumm")
  tkconfigure(rb1, variable = rbValue, value = "nsumm")
  tkconfigure(rb2, variable = rbValue, value = "summ")

  tkgrid(tklabel(ttch, text = "   ---   ---   --- "))
  tkgrid(tklabel(ttch, text = "  Show a summary.", fg="blue"), rb2, sticky="w")
  tkgrid(tklabel(ttch, text = "  Analyse selected frequencies:", fg="blue"), rb1, sticky="w")

  if(.wvari$s==12){
    pi6  <- tkcheckbutton(ttch)
    pi56 <- tkcheckbutton(ttch)
    pi3  <- tkcheckbutton(ttch)
    pi23 <- tkcheckbutton(ttch)
           }
  pi2  <- tkcheckbutton(ttch)
  pi   <- tkcheckbutton(ttch)

  if(.wvari$s==12){
    pi6Value  <- tclVar("0")
    pi56Value <- tclVar("0")
    pi3Value  <- tclVar("0")
    pi23Value <- tclVar("0")
           }
  pi2Value  <- tclVar("0")
  piValue   <- tclVar("0")

  if(.wvari$s==12){
    tkconfigure(pi6, variable=pi6Value)
    tkgrid(tklabel(ttch, text="pi/6"), pi6)
    tkconfigure(pi56, variable=pi56Value)
    tkgrid(tklabel(ttch, text="5pi/6"), pi56)
    tkconfigure(pi3, variable=pi3Value)
    tkgrid(tklabel(ttch, text="pi/3"), pi3)
    tkconfigure(pi23, variable=pi23Value)
    tkgrid(tklabel(ttch, text="2pi/3"), pi23)
           }
  tkconfigure(pi2, variable=pi2Value)
  tkgrid(tklabel(ttch, text="pi/2"), pi2)
  tkconfigure(pi, variable=piValue)
  tkgrid(tklabel(ttch, text="pi"), pi)

  #
  tkconfigure(.treeWidget, cursor="watch")
  OnOK <- function()
  {
    tkdestroy(ttch)
    trendVal <- as.character(tclvalue(trendValue))
      if(trendVal =="1")
         DetTr <- TRUE
      if(trendVal =="0")
         DetTr <- FALSE
   lag1Val <- as.character(tclvalue(lag1Value))
      if(lag1Val =="1")
         f0 <- 1
      if(lag1Val =="0")
         f0 <- 0

    if (tclvalue(rbValue) == "nsumm")
    {
      frec <- rep(0, .wvari$s/2)
      if(.wvari$s == 12)
      {
        mVal <- as.character(tclvalue(pi6Value))
        if(mVal =="1")
           frec[1] <- 1
        if(mVal =="0")
           frec[1] <- 0
        mVal <- as.character(tclvalue(pi56Value))
        if(mVal =="1")
           frec[5] <- 1
        if(mVal =="0")
           frec[5] <- 0
        mVal <- as.character(tclvalue(pi3Value))
        if(mVal =="1")
          frec[2] <- 1
        if(mVal =="0")
           frec[2] <- 0
        mVal <- as.character(tclvalue(pi23Value))
        if(mVal =="1")
           frec[4] <- 1
        if(mVal =="0")
            frec[4] <- 0
      }
      if(.wvari$s == 12){ aux <- c(3,6) }
      if(.wvari$s == 4){ aux <- c(1,2) }
      mVal <- as.character(tclvalue(pi2Value))
        if(mVal =="1")
           frec[aux[1]] <- 1
        if(mVal =="0")
           frec[aux[1]] <- 0
      mVal <- as.character(tclvalue(piValue))
        if(mVal =="1")
           frec[aux[2]] <- 1
        if(mVal =="0")
           frec[aux[2]] <- 0
      rdoCH <- CH.test(.wvari, frec, f0, DetTr=DetTr, showcat=TRUE) # l)
    }
    if (tclvalue(rbValue) == "summ")
    {
      rdoCH <- c(1:(.wvari$s/2+1))
      frecf <- rep(0, .wvari$s/2)
      for(i in 1:(.wvari$s/2))
      {
         frecf[i] <- 1
         rdoCH[i] <- CH.test(.wvari, frecf, f0, DetTr=DetTr, showcat=FALSE)[1]
         frecf <- rep(0, .wvari$s/2)
     }
     rdoCH[(.wvari$s/2+1)] <-
            CH.test(.wvari, rep(1, .wvari$s/2), f0, DetTr=DetTr, showcat=FALSE)[1]
     lCH <- CH.test(.wvari, rep(1, .wvari$s/2), f0, DetTr=FALSE, showcat=FALSE)[2]
     CH  <- round(as.numeric(rdoCH), 2)
     if(.wvari$s==12)
         rdoch <- data.frame("f.pi.6"=CH[1], "f.pi.3"=CH[2], "f.pi.2"=CH[3], "f.2pi.3"=CH[4],
                             "f.5pi.6"=CH[5], "f.pi"=CH[6], "Joint-test"=CH[7])
      if(.wvari$s==4)
         rdoch <- data.frame("f.pi.2"=CH[5], "f.pi"=CH[6], "Joint-test"=CH[7])
      cat("\n ------ CH test ------ \n\n")
      print(rdoch)
    }
    tkconfigure(.treeWidget, cursor="xterm")
  }
  OK.but <- tkbutton(ttch,text="OK", command=OnOK)
  tkgrid(OK.but)
}

MakeCHseas.test <- function()
{
  ttch <- tktoplevel()
  tkwm.title(ttch, "CH test")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  rb1 <- tkradiobutton(ttch)
  rb2 <- tkradiobutton(ttch)
  rbValue <- tclVar("nsumm")
  tkconfigure(rb1, variable = rbValue, value = "nsumm")
  tkconfigure(rb2, variable = rbValue, value = "summ")

  tkgrid(tklabel(ttch, text = "  Show a summary.", fg="blue"), rb2, sticky="w")
  tkgrid(tklabel(ttch, text = " "))
  tkgrid(tklabel(ttch, text = "  Analyse an individual season.", fg="blue"), rb1, sticky="w")

  scr <- tkscrollbar(ttch, repeatinterval=5, command=function(...)tkyview(tl,...))
  tl<-tklistbox(ttch,height=4,selectmode="single", yscrollcommand= function(...)tkset(scr,...),
                background="white")
  tkgrid(tklabel(ttch, text="    Select the season to analyse:"), sticky="w")
  tkgrid(tl, scr)
  tkgrid.configure(scr, rowspan=4, sticky="nsw")
  if(.wvari$s==4)
    seas <- c("Quarter 1", "Quarter 2 ", "Quarter 3", "Quarter 4")
  if(.wvari$s==12)
    seas <- c("January", "February", "March", "April", "May", "June", "July", "August",
              "September", "October", "November", "December")
  for (i in (1:12))
      tkinsert(tl, "end", seas[i])
  tkselection.set(tl, 1)

  Confirm <- function(){
     seasChoice <<- which(seas==seas[as.numeric(tkcurselection(tl))+1])
     tkmessageBox(title="CH input", message="Season has been selected.", icon="info")
     cat(c("\n Season selected:", seas[seasChoice], " \n\n"))
                       }
  Confirm.but <-tkbutton(ttch, text="Select", command=Confirm)
  tkgrid(Confirm.but)

  OnOK <- function()
  {
    if (tclvalue(rbValue) == "summ")
      rdoCH <- CHseas.test(.wvari, 4, c(1:.wvari$s), showcat=TRUE)
    if (tclvalue(rbValue) == "nsumm"){
      rdoCH <- CHseas.test(.wvari, 4, seasChoice, showcat=TRUE)
      rm(seasChoice)
                                    }
    tkdestroy(ttch)
    tkconfigure(.treeWidget, cursor="xterm")
  }
  OK.but <- tkbutton(ttch,text="OK",command=OnOK)
  tkgrid(tklabel(ttch, text = " "))
  tkgrid(OK.but)
}

#

MakeHEGY.test <- function()
{                                # anyadir selecP <- "v1"
  tthegybm  <- tktoplevel()
  tkwm.title(tthegybm, "HEGY test")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  N <- .wvari$N

  tkgrid(tklabel(tthegybm, text="  Select the deterministic components:", fg="blue"), sticky="w")

  Cte  <- tkcheckbutton(tthegybm)
  Tdl  <- tkcheckbutton(tthegybm)
  Vfet <- tkcheckbutton(tthegybm)

  CteValue  <- tclVar("0")
  TdlValue  <- tclVar("0")
  VfetValue <- tclVar("0")

  tkconfigure(Cte, variable=CteValue)
    tkgrid(tklabel(tthegybm, text="Intercept"), Cte)
  tkconfigure(Tdl, variable=TdlValue)
    tkgrid(tklabel(tthegybm, text="Trend"), Tdl)
  tkconfigure(Vfet, variable=VfetValue)
    tkgrid(tklabel(tthegybm, text="Seasonal dummys"), Vfet)

  rb1 <- tkradiobutton(tthegybm)
  rb2 <- tkradiobutton(tthegybm)
  rb3 <- tkradiobutton(tthegybm)
  rb4 <- tkradiobutton(tthegybm)
  rb5 <- tkradiobutton(tthegybm)
  rb6 <- tkradiobutton(tthegybm)
  rbValue <- tclVar("BIC-LB")
  yself <- tclVar()
  entry.yself <- tkentry(tthegybm, width="3", textvariable=yself)
  tkconfigure(rb1, variable=rbValue, value="AIC-LB")
  tkconfigure(rb2, variable=rbValue, value="BIC-LB")
  tkconfigure(rb3, variable=rbValue, value="AIC-Lüt")
  tkconfigure(rb4, variable=rbValue, value="BIC-Lüt")
  tkconfigure(rb5, variable=rbValue, value="Signf")
  tkconfigure(rb6, variable=rbValue, value="Tu mismo")
  tkgrid(tklabel(tthegybm, text="  Select the method for choosing lags:", fg="blue"), sticky="w")
  tkgrid(tklabel(tthegybm, text="AIC-LB"), rb1)
  tkgrid(tklabel(tthegybm, text="BIC-LB"), rb2)
  tkgrid(tklabel(tthegybm, text="AIC-top-down"), rb3)
  tkgrid(tklabel(tthegybm, text="BIC-top-down"), rb4)
  tkgrid(tklabel(tthegybm, text="Significant lags"), rb5)
  tkgrid(tklabel(tthegybm, text="By yourself"), rb6, entry.yself)

  Definir1 <- tkbutton(tthegybm, text="Define", command=Makevfic)
  Definir2 <- tkbutton(tthegybm, text="Define", command=MakeVFEp)
  Mvfic <<- matrix(NA, nrow=N, ncol=6)
  tkgrid(tklabel(tthegybm, text="  Include dummys:", fg="blue"), sticky="w")
  tkgrid(tklabel(tthegybm, text="  Generic dummy           "), Definir1)
  tkgrid(tklabel(tthegybm, text="Partial seasonal dummy     "), Definir2)

  tkconfigure(.treeWidget, cursor="watch")
  OnOK <- function()
  {
      tkdestroy(tthegybm)
      ifelse(exists("VFEp"), VFEp <- VFEp, VFEp <<- 0)

      mVal <- as.character(tclvalue(CteValue))
      if(mVal =="1")
         Ct <- 1
      if(mVal =="0")
         Ct <- 0
      mVal <- as.character(tclvalue(TdlValue))
      if(mVal =="1")
         TD <- 1
      if(mVal =="0")
         TD <- 0
      mVal <- as.character(tclvalue(VfetValue))
      if(mVal =="1")
         Vfe <- 1
      if(mVal =="0")
         Vfe <- 0

      # SelecP
      rbVal <- as.character(tclvalue(rbValue))
      if(tclvalue(rbValue) == "AIC-LB")
        selecP <- "aiclb"
      if(tclvalue(rbValue) == "BIC-LB")
        selecP <- "biclb"
      if(tclvalue(rbValue) == "AIC-Lüt")
        selecP <- "aiclut"
       if(tclvalue(rbValue) == "BIC-Lüt")
        selecP <- "biclut"
      if(tclvalue(rbValue) == "Signf")
        selecP <- "signf"
      if(tclvalue(rbValue) == "Tu mismo")
        selecP <- as.numeric(tclvalue(yself))

      # VFIC
      aux <- length(which(Mvfic[1,] >= 0))
      ifelse(aux == 0, Mvfic <<- 0, Mvfic <<- as.matrix(Mvfic[,1:aux]))

      rdoHEGYBM <- HEGY.test(.wvari, c(Ct,TD,Vfe), selecP, Mvfic, VFEp, showcat=TRUE)
      cleanlist <- list(Mvfic=Mvfic, VFEp=VFEp)
      clean.auxobjects(cleanlist, type="")
      tkconfigure(.treeWidget, cursor="xterm")
  }
  OK.but <- tkbutton(tthegybm,text="OK",command=OnOK)
  tkgrid(OK.but)
}


# MENU DATOS
# modo menu, "DataInfo" ó ""SDGP"

GetDataInfo <- function()   # Datos-Descripción
{
  ttdinfo  <- tktoplevel()
  tkwm.title(ttdinfo, "Description")

  nombre <- tclVar()
  tkgrid(tklabel(ttdinfo, text="  Series label:", fg="blue"), sticky="w")
  entry.nombre <- tkentry(ttdinfo, width="10", textvariable=nombre)
  tkgrid(entry.nombre)

  rb1     <- tkradiobutton(ttdinfo)
  rb2     <- tkradiobutton(ttdinfo)
  rb3     <- tkradiobutton(ttdinfo)
  rbValue <- tclVar("Trimestral")
  tkconfigure(rb1, variable=rbValue, value="Trimestral")
  tkconfigure(rb2, variable=rbValue, value="Mensual")
  tkconfigure(rb3, variable=rbValue, value="Anual")
  tkgrid(tklabel(ttdinfo, text="  Register the periodicity of the series:",
      fg="blue"), sticky="w")
  tkgrid(tklabel(ttdinfo, text="Quarterly"), rb1)
  tkgrid(tklabel(ttdinfo, text="Monthly"), rb2)
  tkgrid(tklabel(ttdinfo, text="Anual"), rb3)

  s0 <- tclVar()
  a0 <- tclVar()

  entry.a0 <- tkentry(ttdinfo, width="4", textvariable=a0)
  entry.s0 <- tkentry(ttdinfo, width="2", textvariable=s0)
  tkgrid(tklabel(ttdinfo, text="  Introduce the year and season of the"), sticky="w")
  tkgrid(tklabel(ttdinfo, text="     first observation:", fg="blue"),
       entry.a0, entry.s0, sticky="w")

  OnOK <- function()
  {
     label <- as.character(tclvalue(nombre))

     #periodicidad
     rbVal <- as.character(tclvalue(rbValue))
     tkdestroy(ttdinfo)
     if(tclvalue(rbValue) == "Trimestral")
        s <- 4
     if(tclvalue(rbValue) == "Mensual")
        s <- 12
     if(tclvalue(rbValue) == "Anual")
        s <- 1

     t0    <- rep(0, 2)
     t0[1] <- as.numeric(tclvalue(a0))
     t0[2] <- as.numeric(tclvalue(s0))
     vari <- ts(vari, frequency=s, start=c(t0[1],t0[2]))
     N    <- length(vari)
     logvari <- FALSE # para gráficos (títulos),
                      # cuidado si se usan datos ya en logaritmos

     assign(label, list(vari=vari, s=s, t0=t0, N=N, logvari=logvari, label=label), env=.GlobalEnv)
     assign(".wvari", list(vari=vari, s=s, t0=t0, N=N, logvari=logvari, label=label),
            env=.GlobalEnv)
     tkinsert(.treeWidget,"end","root",label,text=label)

     tempf <- tempfile()
     cat("rm(datos, vari)", file=tempf)
     eval(source(tempf))
     tkcmd("file", "delete", tempf)

     msg <- paste("Information about the series has been stored in the object ",
                  label, "\n", sep="")
     tkmessageBox(title="Series info", message=msg, icon="info")
     tkdestroy(ttdinfo)
  }
  OK.but <- tkbutton(ttdinfo, text="OK", command=OnOK)
  tkgrid(OK.but)
}

#
MakeDGPsim <- function()
{
  ttdgp  <- tktoplevel()
  tkwm.title(ttdgp, "Simulate DGP")

  rb1     <- tkradiobutton(ttdgp)
  rb2     <- tkradiobutton(ttdgp)
  rb3     <- tkradiobutton(ttdgp)
  rbValue <- tclVar("Trimestral")
  tkconfigure(rb1, variable=rbValue, value="Trimestral")
  tkconfigure(rb2, variable=rbValue, value="Mensual")
  tkconfigure(rb3, variable=rbValue, value="Anual")
  tkgrid(tklabel(ttdgp, text="  Introduce the periodicity of the series:",
      fg="blue"), sticky="w")
  tkgrid(tklabel(ttdgp, text="Quarterly"), rb1)
  tkgrid(tklabel(ttdgp, text="Monthly"), rb2)
  tkgrid(tklabel(ttdgp, text="Anual"), rb3)

  N <- tclVar()
  entry.N <- tkentry(ttdgp, width="5", textvariable=N)
  tkgrid(tklabel(ttdgp, text="  Introduce the number of observations:",
      fg="blue"), sticky="w")
  tkgrid(entry.N)

  phi1Var    <- tclVar(" ")
  phi2Var    <- tclVar(" ")
  phi3Var    <- tclVar(" ")
  phi4Var    <- tclVar(" ")
  retphi1Var <- tclVar(" ")
  retphi2Var <- tclVar(" ")
  retphi3Var <- tclVar(" ")
  retphi4Var <- tclVar(" ")
  rho1Var    <- tclVar(" ")
  rho2Var    <- tclVar(" ")
  rho3Var    <- tclVar(" ")
  rho4Var    <- tclVar(" ")
  retrho1Var <- tclVar(" ")
  retrho2Var <- tclVar(" ")
  retrho3Var <- tclVar(" ")
  retrho4Var <- tclVar(" ")

  phi1Entry    <- tkentry(ttdgp, width="3", textvariable=phi1Var)
  phi2Entry    <- tkentry(ttdgp, width="3", textvariable=phi2Var)
  phi3Entry    <- tkentry(ttdgp, width="3", textvariable=phi3Var)
  phi4Entry    <- tkentry(ttdgp, width="3", textvariable=phi4Var)
  retphi1Entry <- tkentry(ttdgp, width="3", textvariable=retphi1Var)
  retphi2Entry <- tkentry(ttdgp, width="3", textvariable=retphi2Var)
  retphi3Entry <- tkentry(ttdgp, width="3", textvariable=retphi3Var)
  retphi4Entry <- tkentry(ttdgp, width="3", textvariable=retphi4Var)
  rho1Entry    <- tkentry(ttdgp, width="3", textvariable=rho1Var)
  rho2Entry    <- tkentry(ttdgp, width="3", textvariable=rho2Var)
  rho3Entry    <- tkentry(ttdgp, width="3", textvariable=rho3Var)
  rho4Entry    <- tkentry(ttdgp, width="3", textvariable=rho4Var)
  retrho1Entry <- tkentry(ttdgp, width="3", textvariable=retrho1Var)
  retrho2Entry <- tkentry(ttdgp, width="3", textvariable=retrho2Var)
  retrho3Entry <- tkentry(ttdgp, width="3", textvariable=retrho3Var)
  retrho4Entry <- tkentry(ttdgp, width="3", textvariable=retrho4Var)

  tkgrid(tklabel(ttdgp,
  text="  Fill in the following information.", fg="blue"), sticky="w")

  tkgrid(tklabel(ttdgp, text="Autorregresive lags:"))
  tkgrid(tklabel(ttdgp, text="Coefficient:"), phi1Entry, phi2Entry, phi3Entry, phi4Entry)
  tkgrid(tklabel(ttdgp, text="Lag:"), retphi1Entry, retphi2Entry, retphi3Entry, retphi4Entry)
  tkgrid(tklabel(ttdgp, text="Moving average lags:"))
  tkgrid(tklabel(ttdgp, text="Coefficient:"), rho1Entry, rho2Entry, rho3Entry, rho4Entry)
  tkgrid(tklabel(ttdgp, text="Lag:"), retrho1Entry, retrho2Entry, retrho3Entry, retrho4Entry)

  OnOK <- function()
  {
     rbVal <- as.character(tclvalue(rbValue))
     tkdestroy(ttdgp)
     if(tclvalue(rbValue) == "Trimestral")
        s <- 4
     if(tclvalue(rbValue) == "Mensual")
        s <- 12
     if(tclvalue(rbValue) == "Anual")
        s <- 1
     N <- as.numeric(tclvalue(N))
     phi    <- rbind(na.omit(c( as.numeric(tclvalue(phi1Var)), as.numeric(tclvalue(phi2Var)),
                   as.numeric(tclvalue(phi3Var)), as.numeric(tclvalue(phi4Var)) )))
     retphi <- rbind(na.omit(c( as.numeric(tclvalue(retphi1Var)),
                   as.numeric(tclvalue(retphi2Var)), as.numeric(tclvalue(retphi3Var)),
                   as.numeric(tclvalue(retphi4Var)) )))
     rho    <- na.omit(c( as.numeric(tclvalue(rho1Var)), as.numeric(tclvalue(rho2Var)),
                  as.numeric(tclvalue(rho3Var)), as.numeric(tclvalue(rho4Var)) ))
     retrho <- na.omit(c( as.numeric(tclvalue(retrho1Var)), as.numeric(tclvalue(retrho2Var)),
                  as.numeric(tclvalue(retrho3Var)), as.numeric(tclvalue(retrho4Var)) ))
     vari <- DGPsim(s, N, phi, retphi, rho, retrho); t0 <- c(0,1) ;label <- "DGP"

     assign(label, list(vari=vari, s=s, t0=t0, N=N, logvari=FALSE, label=label), env=.GlobalEnv)
     assign(".wvari", list(vari=vari, s=s, t0=t0, N=N, logvari=FALSE, label=label),
            env=.GlobalEnv)
     tkinsert(.treeWidget,"end","root",label,text=label)

     msg <- paste("Data generating process information has been stored in the object ",
                  label, sep="")
     tkmessageBox(title="Series info", message=msg, icon="info")
     tkdestroy(ttdgp)
  }
  OK.but <- tkbutton(ttdgp, text="OK", command=OnOK)
  tkgrid(OK.but)
}

#

ReadDataTXT <- function()
{
  tttxt <- tktoplevel()
  tkwm.title(tttxt, "Read data from text file")

  espaciosep <- tkradiobutton(tttxt)
  comasep    <- tkradiobutton(tttxt)
  ptcomasep  <- tkradiobutton(tttxt)
  sepValue   <- tclVar(";")

  tkconfigure(espaciosep, variable=sepValue, value=" ")
  tkconfigure(comasep, variable=sepValue, value=",")
  tkconfigure(ptcomasep, variable=sepValue, value=";")
  tkgrid(tklabel(tttxt, text="  Separator character:", fg="blue"), sticky="w")
  tkgrid(tklabel(tttxt, text="White space"), espaciosep)
  tkgrid(tklabel(tttxt, text="Comma"), comasep)
  tkgrid(tklabel(tttxt, text="Semicolon"), ptcomasep)

  comadec  <- tkradiobutton(tttxt)
  puntodec <- tkradiobutton(tttxt)
  decValue <- tclVar(",")
  tkconfigure(comadec, variable=decValue, value=",")
  tkconfigure(puntodec, variable=decValue, value=".")
  tkgrid(tklabel(tttxt, text="  Character for decimal points:", fg="blue"), sticky="w")
  tkgrid(tklabel(tttxt, text="Comma"), comadec)
  tkgrid(tklabel(tttxt, text="Dot"), puntodec)

  tkgrid(tklabel(tttxt, text="  The data file contain a header with the series names:",
         fg="blue"), sticky="w")
  header      <- tkcheckbutton(tttxt)
  headerValue <- tclVar("1")

  tkconfigure(header, variable=headerValue)
  tkgrid(tklabel(tttxt, text="With header"), header)

  datacolumn <- tclVar()
  tkgrid(tklabel(tttxt, text="  Introduce the column that contains the series to analyse:",
         fg="blue"), sticky="w")
  entry.datacolumn <- tkentry(tttxt, width="5", textvariable=datacolumn)
  tkgrid(entry.datacolumn)

  OnOK <- function(){
    tkdestroy(tttxt)
    datafile <- tclvalue(tkgetOpenFile(filetypes='{"Text Files" {".txt" ".TXT" ".dat" ".DAT"}} {"All Files" {"*"}}'))
    if(!nchar(datafile))
       tkmessageBox(message="No file was chosen.", icon="error")
    else
       tkmessageBox(message=paste("The file chosen is", datafile, "and has been stored in the object datos."), icon="info")
     mVal <- as.character(tclvalue(headerValue))
     if(mVal =="1")
        headerarg <- TRUE
     if(mVal =="0")
        headerarg <- FALSE

     datos <<- read.csv(datafile, header=headerarg,
                   sep=tclvalue(sepValue), dec=tclvalue(decValue),
                   na.strings="NA")

     variaux <- datos[,as.numeric(tclvalue(datacolumn))]; # rm(datos)
     vari    <<- as.numeric(as.matrix(variaux))

     msg <- paste("Don't forget describe the series (Menu-Data->Description)", sep="")
     tkmessageBox(title="Series info", message=msg, icon="info")
   }
   OK.but <- tkbutton(tttxt, text="OK", command=OnOK)
   tkgrid(OK.but)
}

#

# Cuando hay  8 NA seguidos considera que se acaba la serie
ReadDataCSV <- function()
{
  tttxt <- tktoplevel()
  tkwm.title(tttxt, "Read data from CSV file")

  comadec  <- tkradiobutton(tttxt)
  puntodec <- tkradiobutton(tttxt)
  decValue <- tclVar(",")
  tkconfigure(comadec, variable=decValue, value=",")
  tkconfigure(puntodec, variable=decValue, value=".")
  tkgrid(tklabel(tttxt, text="  Character for decimal points:", fg="blue"), sticky="w")
  tkgrid(tklabel(tttxt, text="Comma"), comadec)
  tkgrid(tklabel(tttxt, text="Dot"), puntodec)

  tkgrid(tklabel(tttxt, text="  The data file contain a header with the series names:",
       fg="blue"), sticky="w")
  header      <- tkcheckbutton(tttxt)
  headerValue <- tclVar("1")

  tkconfigure(header, variable=headerValue)
  tkgrid(tklabel(tttxt, text="With header"), header)

  datacolumn <- tclVar()
  tkgrid(tklabel(tttxt, text="  Introduce the column that contains the series to analyse:",
      fg="blue"), sticky="w")
  entry.datacolumn <- tkentry(tttxt, width="5", textvariable=datacolumn)
  tkgrid(entry.datacolumn)

 OnOK <- function(){
    tkdestroy(tttxt)
    datafile <- tclvalue(tkgetOpenFile(filetypes='{"CSV Files" {".csv" ".CSV"}} {"All Files" {"*"}}'))
    if(!nchar(datafile))
       tkmessageBox(message="No file was chosen.", icon="error")
    else
       tkmessageBox(message=paste("The file chosen was", datafile, "and has been stored in the object named datos."), icon="info")
     mVal <- as.character(tclvalue(headerValue))
     if(mVal =="1")
        headerarg <- TRUE
     if(mVal =="0")
        headerarg <- FALSE

     datos <<- read.csv(datafile, header=headerarg,
                   sep=" ", dec=tclvalue(decValue), na.strings="NA")

     variaux <- datos[,as.numeric(tclvalue(datacolumn))]; # rm(datos)
     vari    <- as.numeric(as.matrix(variaux))

     aux1 <- sort(c(which(vari>0), which(vari<0)))
     aux2 <- c(1:(length(aux1)-1))
     for(i in 2:length(aux1))
       aux2[i] <- aux1[i]-aux1[i-1]
     if(length(which(aux2 > 1)) == 0)
       N <- aux1[length(aux1)]
     if(length(which(aux2 > 1)) > 0){
       if(length(which(aux2 == 8)) == 0)
         N <- length(vari)
       if(length(which(aux2 == 8)) > 0)
         N <- which(aux2 == 8)[1] - 1
                                    }
     vari <<- vari[1:N]

     msg <- paste("Don't forget describe the series (Menu-Data->Description)", sep="")
     tkmessageBox(title="Series info", message=msg, icon="info")
   }
   OK.but <- tkbutton(tttxt, text="OK", command=OnOK)
   tkgrid(OK.but)
}

# Abrir SPSS

ReadSPSS <- function()
{
  name <- tclvalue(tkgetOpenFile(
                   filetypes="{{SPSS Files} {.sav}} {{All files} *}"))
  if (name=="") return;
  zz <- read.spss(name,use.value.label=T,to.data.frame=T)
  assign("myData", zz, envir=.GlobalEnv)
}

#

ysooys <- function(yso, t0, N, s)
{
  index <- matrix(-9999, ncol=3, nrow=N)
  index[,3] <- c(1:N)
  index[,2] <- c(c(t0[2]:s), rep(1:s,N/s))[1:N]
  index[1:length(t0[2]:s),1] <- rep(t0[1], length(t0[2]:s))
  i <- 1
  while(index[N,1] == -9999)
  {
     iaux <- which(index[,1] == -9999)
     reps <- ifelse(length(iaux) >= s, reps <- s, reps <- length(iaux))
     index[iaux[1]:(iaux[1]+(reps-1)),1] <- rep(t0[1]+i, reps)
     i <- i+1
  }

# year and season to observation
  if(length(yso)==2)
  {
    quest1 <- which(c(index[,1] == yso[1]) == TRUE)
    quest2 <- which(c(index[,2] == yso[2]) == TRUE)
     i <- 1; out <- quest1[1]
    while(length(which(quest2 == quest1[i])) != 1){
      i <- i+1
      out <- quest1[i]
    }
  }

# observation to year and season
  if(length(yso)==1)
    out <- index[which(index[,3]==yso),1:2]

  list(out, index)
}

#

CambiarPeriodo <- function()
{
  ttpmtr <- tktoplevel()
  tkwm.title(ttpmtr, "Change sample period")

  a02 <- tclVar()
  s02 <- tclVar()
  aN2 <- tclVar()
  sN2 <- tclVar()

  tkgrid(tklabel(ttpmtr, text="  New sample period:", fg="blue"), sticky="w")
  entry.a02 <- tkentry(ttpmtr, width="4", textvariable=a02)
  entry.s02 <- tkentry(ttpmtr, width="2", textvariable=s02)
  tkgrid(tklabel(ttpmtr, text="     Year and season of the first observation:"),
         entry.a02, entry.s02)
  entry.aN2 <- tkentry(ttpmtr, width="4", textvariable=aN2)
  entry.sN2 <- tkentry(ttpmtr, width="2", textvariable=sN2)
  tkgrid(tklabel(ttpmtr, text="     Year and season of the last observation:"),
         entry.aN2, entry.sN2)

  OnOK <- function(){
    tkdestroy(ttpmtr)
    string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
    .wvari <<- ExeString(c(string, ""))

    obs1 <- ysooys(c(as.numeric(tclvalue(a02)), as.numeric(tclvalue(s02))), .wvari$t0,
                   .wvari$N, .wvari$s)[[1]]
    obs2 <- ysooys(c(as.numeric(tclvalue(aN2)), as.numeric(tclvalue(sN2))), .wvari$t0,
                   .wvari$N, .wvari$s)[[1]]
    t0   <- c(as.numeric(tclvalue(a02)), as.numeric(tclvalue(s02)))
    tN <- c(as.numeric(tclvalue(aN2)), as.numeric(tclvalue(sN2)))
    vari <- ts(.wvari$vari[obs1:obs2], frequency=.wvari$s, start=t0)

    i <- 1; newstring <- "NULL"
    while(newstring == "NULL")
    {
      if(exists(paste(string, "_subs", i, sep="")))
        i <- i+1
      logic <- exists(paste(string, "_subs", i, sep="")) == TRUE
      if(logic == FALSE)
        newstring <- paste(string, "_subs", i, sep="")
    }
    assign(newstring, list(vari=vari, s=.wvari$s, t0=t0, N=length(vari), logvari=.wvari$logvari,
            label=newstring), env=.GlobalEnv)
    showlabel <- paste(newstring, "     ", t0[1], ".", t0[2], "-", tN[1], ".", tN[2], sep="")
    tkinsert(.treeWidget,"end", string, newstring, text=showlabel)

    #change.wvari.label(c("t0", "t0cp", "vari", "N"))
    tkmessageBox(message="The sample period has changed.", icon="info")
                    }
  OK.but <- tkbutton(ttpmtr, text="OK", command=OnOK)
  tkgrid(OK.but)
}

# MENU GRÁFICOS
Makeplotvari <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  plot(.wvari$vari, xlab="", ylab="", las=1, main=.wvari$label)
  print(.wvari$vari)
}
Makeplotlog    <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  plot(log(.wvari$vari, base=exp(1)), xlab="", ylab="", las=1, main="Logarithm of the series")
  print(log(.wvari$vari, base=exp(1)))
}
Makeplotboxcox <- function()
{
  ttplot <- tktoplevel()

  lambda  <- tclVar(0)
  tkgrid(tklabel(ttplot, text=" Introduce Box-Cox algorithm parameter:"), sticky="w")
  entry.lambda <- tkentry(ttplot, width="4", textvariable=lambda)
  tkgrid(tklabel(ttplot, text="lambda"), entry.lambda)

  OnOK <- function()
  {
    tkdestroy(ttplot)
    string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
    .wvari <<- ExeString(c(string, ""))
    varip  <- ts(BoxCox(.wvari$vari, as.numeric(tclvalue(lambda))),
                 frequency=.wvari$s, start=.wvari$t0)
    plot(varip, xlab="", ylab="", las=1, main="Box-Cox transformation")
    print(varip)
  }
  tkgrid(tkbutton(ttplot,text="OK",command=OnOK))
}
Makeboxcox <- function()
{
  ttbc <- tktoplevel()
  lambda  <- tclVar(0)
  tkgrid(tklabel(ttbc, text=" Introduce Box-Cox algorithm parameter:", fg="blue"), sticky="w")
  entry.lambda <- tkentry(ttbc, width="4", textvariable=lambda)
  tkgrid(tklabel(ttbc, text="lambda"), entry.lambda)
  OnOK <- function()
  {
     tkdestroy(ttbc)

     string    <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
     newstring <- paste("boxcox_", string, sep="")
     varit     <- ts(BoxCox(.wvari$vari, as.numeric(tclvalue(lambda))),
                     frequency=.wvari$s, start=.wvari$t0)
     assign(newstring, list(vari=varit, s=.wvari$s, t0=.wvari$t0, N=length(varit),
            logvari=FALSE, label=newstring), env=.GlobalEnv)
     tkinsert(.treeWidget,"end", string, newstring, text=newstring)

     #change.wvari.label(c("label", "logvari", "t0", "varit", "vari", "N"))
     tkmessageBox(message="The scale has changed.", icon="info")
  }
  tkgrid(tkbutton(ttbc,text="OK",command=OnOK))
}
#
Makeplotdelta <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  s <- .wvari$s
  N <- .wvari$N
  t0 <- .wvari$t0

  eje.m  <- rep(seq(0,(1/s)*(s-1),1/s), as.integer(N/s)+1)
  eje.a  <- c(t0[1]:(t0[1]+as.integer(N/s)))
  eje.am <- c(1:length(eje.m))
  aux    <- as.numeric(gl(length(eje.a), s))
  k <- 1
  for(i in 1:length(eje.m)){
    eje.am[k] <- t0[1]-1+aux[i]+eje.m[i]
    k<- k+1
                           }
  eje.am2 <- eje.am[seq(1,length(eje.am),s/2)]
  plot(diff(.wvari$vari, lag=1), xlab="", ylab="", las=1, main="First differences of the series")
  abline(v=eje.am2, lty=2, col="blue")
  print(diff(.wvari$vari, lag=1))
}

Makeplotdeltas  <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  s <- .wvari$s
  N <- .wvari$N
  t0 <- .wvari$t0

  eje.m  <- rep(seq(0,(1/s)*(s-1),1/s), as.integer(N/s)+1)
  eje.a  <- c(t0[1]:(t0[1]+as.integer(N/s)))
  eje.am <- c(1:length(eje.m))
  aux    <- as.numeric(gl(length(eje.a), s))
  k <- 1
  for(i in 1:length(eje.m)){
    eje.am[k] <- t0[1]-1+aux[i]+eje.m[i]
    k<- k+1
                           }
  eje.am2 <- eje.am[seq(1,length(eje.am),s/2)]
  plot(diff(.wvari$vari, lag=s), xlab="", ylab="",las=1, main="Seasonal differences of the series")
  abline(v=eje.am2, lty=2, col="blue")
  print(diff(.wvari$vari, lag=s))
}
Makeplotddeltas <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  plot(diff(diff(.wvari$vari, lag=.wvari$s), lag=1), xlab="", ylab="", las=1, main="First and sesonal differences of the series")
  #abline(v=eje.am2, lty=2, col="blue")
  print(diff(diff(.wvari$vari, lag=.wvari$s), lag=1))
}
Makeplotperdiff  <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  difpvari <- perdiff(.wvari)
  plot(difpvari, xlab="", ylab="", las=1, main="Periodic differences of the series")
  print(difpvari)
}
Makevariodet  <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  Transfdet(.wvari$vari, .wvari$s)
  plot(variodet, xlab="", ylab="",las=1, main="Series without deterministic components")
  print(variodet)
  rm(variodet)
}
Makespec <- function(dif1)
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  vari <- .wvari$vari
  if(dif1 == FALSE){
     opar <- par(tcl=-0.5, cex.axis=0.8, cex.main=1, las=1)
     spectrum(vari, spans=c(3,5), log="no", ylab="", main="Estimated spectral density")
     par(opar)
                   }
  if(dif1 == TRUE){
     opar <- par(tcl=-0.5, cex.axis=0.8, cex.main=1, las=1)
     spectrum(diff(vari, lag=1), spans=c(3,5), log="no", ylab="",
           main="Estimated spectral density upon the first differences")
     par(opar)
                  }
}
#
Makequarterg <- function(transf)
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  if(.wvari$s != 4)
  tkmessageBox(title="Consult the help file", message="This plot is refered only to quarterly series.", icon="error")
  stopifnot(.wvari$s == 4) # o poner un if(s ==4){ hacer todo lo demás }

  if(transf == "orig")
    varibb <- .wvari$vari
  if(transf == "fdiff"){
    ML <- ret(.wvari$vari, 2); varibb <- ts(ML[,1]-ML[,2], frequency= .wvari$s, start= .wvari$t0) }
  if(transf == "pdiff")
    varibb <- perdiff(.wvari)

  opar <- par(mar=c(8,4.7,5,1.5), ps=18, font=1,tcl=-0.5, cex.axis=0.7, las=1)
  quarterg(varibb, .wvari$s, .wvari$t0, plot=TRUE)
  par(opar)
}

#

Makebbmp <- function(transf)
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Buys-Ballot plot")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  rb1     <- tkradiobutton(ttplot)
  rb2     <- tkradiobutton(ttplot)
  rbValue <- tclVar("Panel de graficos con todos los meses")
  tkconfigure(rb1,variable=rbValue, value="Seleccionar meses en un grafico")
  tkconfigure(rb2,variable=rbValue, value="Panel de graficos con todos los meses")
  tkgrid(tklabel(ttplot, text="  Choose the representation type:", fg="blue"), sticky="w")
  tkgrid(tklabel(ttplot, text="    Graphics panel with all months"), rb2, sticky="w")
  tkgrid(tklabel(ttplot, text="    Select months for a one-graphic plot"), rb1, sticky="w")

  tkgrid(tklabel(ttplot, text="Select the months you wish to plot:", fg="blue"))
  tkgrid(tklabel(ttplot, text="(Only if you have chosen the second option.)"))

  En <- tkcheckbutton(ttplot)
  Fe <- tkcheckbutton(ttplot)
  Mar <- tkcheckbutton(ttplot)
  Ab <- tkcheckbutton(ttplot)
  Ma <- tkcheckbutton(ttplot)
  Ju <- tkcheckbutton(ttplot)
  Jul<- tkcheckbutton(ttplot)
  Ag <- tkcheckbutton(ttplot)
  Se <- tkcheckbutton(ttplot)
  Oc <- tkcheckbutton(ttplot)
  No <- tkcheckbutton(ttplot)
  Di <- tkcheckbutton(ttplot)

  EnValue <- tclVar("0")
  FeValue <- tclVar("0")
  MarValue <- tclVar("0")
  AbValue <- tclVar("0")
  MaValue <- tclVar("0")
  JuValue <- tclVar("0")
  JulValue <- tclVar("0")
  AgValue <- tclVar("0")
  SeValue <- tclVar("0")
  OcValue <- tclVar("0")
  NoValue <- tclVar("0")
  DiValue <- tclVar("0")

  tkconfigure(En, variable=EnValue)
  tkgrid(tklabel(ttplot, text="January"), En)
  tkconfigure(Fe, variable=FeValue)
  tkgrid(tklabel(ttplot, text="February"), Fe)
  tkconfigure(Mar, variable=MarValue)
  tkgrid(tklabel(ttplot, text="March"), Mar)
  tkconfigure(Ab, variable=AbValue)
  tkgrid(tklabel(ttplot, text="April"), Ab)
  tkconfigure(Ma, variable=MaValue)
  tkgrid(tklabel(ttplot, text="May"), Ma)
  tkconfigure(Ju, variable=JuValue)
  tkgrid(tklabel(ttplot, text="June"), Ju)
  tkconfigure(Jul, variable=JulValue)
  tkgrid(tklabel(ttplot, text="July"), Jul)
  tkconfigure(Ag, variable=AgValue)
  tkgrid(tklabel(ttplot, text="August"), Ag)
  tkconfigure(Se, variable=SeValue)
  tkgrid(tklabel(ttplot, text="September"), Se)
  tkconfigure(Oc, variable=OcValue)
  tkgrid(tklabel(ttplot, text="October"), Oc)
  tkconfigure(No, variable=NoValue)
  tkgrid(tklabel(ttplot, text="November"), No)
  tkconfigure(Di, variable=DiValue)
  tkgrid(tklabel(ttplot, text="December"), Di)

  OnOK <- function()
  {
    tkdestroy(ttplot)
    s <- .wvari$s
    if(s != 12)
    tkmessageBox(title="Consult the help file", message="This plot is refered only to monthly series.", icon="error")
    stopifnot(s == 12)

    if(transf == "orig")
      varibb <- .wvari$vari
    if(transf == "fdiff"){
      ML <- ret(.wvari$vari, 2)
      varibb <- ts(ML[,1]-ML[,2], frequency= .wvari$s, start= .wvari$t0) }
    if(transf == "pdiff")
      varibb <- perdiff(.wvari)

    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="Seleccionar meses en un grafico")
       vers <- "R"; mplot <- rep(0,12)
    if (rbVal=="Panel de graficos con todos los meses"){
       vers  <- "Prot"
       mplot <- rep(1, 12)                             }
    if(vers == "R")
    {
       tkdestroy(ttplot)
       mVal <- as.character(tclvalue(EnValue))
       if(mVal =="1")
          mplot[1] <- 1
       if(mVal =="0")
          mplot[1] <- 0
       mVal <- as.character(tclvalue(FeValue))
       if(mVal =="1")
          mplot[2] <- 1
       if(mVal =="0")
          mplot[2] <- 0
       mVal <- as.character(tclvalue(MarValue))
        if(mVal =="1")
          mplot[3] <- 1
       if(mVal =="0")
          mplot[3] <- 0
       mVal <- as.character(tclvalue(AbValue))
      if(mVal =="1")
           mplot[4] <- 1
       if(mVal =="0")
          mplot[4] <- 0
       mVal <- as.character(tclvalue(MaValue))
       if(mVal =="1")
          mplot[5] <- 1
       if(mVal =="0")
          mplot[5] <- 0
       mVal <- as.character(tclvalue(JuValue))
       if(mVal =="1")
          mplot[6] <- 1
       if(mVal =="0")
          mplot[6] <- 0
       mVal <- as.character(tclvalue(JulValue))
       if(mVal =="1")
          mplot[7] <- 1
       if(mVal =="0")
          mplot[7] <- 0
       mVal <- as.character(tclvalue(AgValue))
       if(mVal =="1")
          mplot[8] <- 1
       if(mVal =="0")
          mplot[8] <- 0
       mVal <- as.character(tclvalue(SeValue))
       if(mVal =="1")
          mplot[9] <- 1
       if(mVal =="0")
          mplot[9] <- 0
       mVal <- as.character(tclvalue(OcValue))
       if(mVal =="1")
          mplot[10] <- 1
       if(mVal =="0")
          mplot[10] <- 0
       mVal <- as.character(tclvalue(NoValue))
       if(mVal =="1")
          mplot[11] <- 1
       if(mVal =="0")
          mplot[11] <- 0
       mVal <- as.character(tclvalue(DiValue))
       if(mVal =="1")
          mplot[12] <- 1
       if(mVal =="0")
          mplot[12] <- 0
       mplot <- mplot
    }
    if (rbVal=="Panel de graficos con todos los meses"){
      opar <- par(mfrow=c(2,2), mar=c(3,4.5,1.5,1), ps=18, font=1, tcl=-0.5, cex.axis=0.7, las=1)
      bbmp(varibb, .wvari$s, .wvari$t0, which(mplot==1), vers, plot=TRUE)
      par(opar)
                                                       }
    if (rbVal=="Seleccionar meses en un grafico"){
      opar <- par(mar=c(8,4.7,5,1.5), ps=18, font=1, tcl=-0.5, cex.axis=0.7, las=1)
      bbmp(varibb, .wvari$s, .wvari$t0, which(mplot==1), vers, plot=TRUE)
      par(opar)                                  }
    rm(vers, mplot)
  }
  OK.but <- tkbutton(ttplot,text="OK",command=OnOK)
  tkgrid(OK.but)
}

#

Makebbap <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  ttplot  <- tktoplevel()
  tkwm.title(ttplot, "Buys-Ballot plot")

  ap1Var <- tclVar(" ")
  ap2Var <- tclVar(" ")
  ap3Var <- tclVar(" ")
  ap4Var <- tclVar(" ")
  ap5Var <- tclVar(" ")
  ap6Var <- tclVar(" ")

  ap1Entry <- tkentry(ttplot, width="7", textvariable=ap1Var)
  ap2Entry <- tkentry(ttplot, width="7", textvariable=ap2Var)
  ap3Entry <- tkentry(ttplot, width="7", textvariable=ap3Var)
  ap4Entry <- tkentry(ttplot, width="7", textvariable=ap4Var)
  ap5Entry <- tkentry(ttplot, width="7", textvariable=ap5Var)
  ap6Entry <- tkentry(ttplot, width="7", textvariable=ap6Var)

  tkgrid(tklabel(ttplot,
       text="Enter the years to plot:", fg="blue"))
  tkgrid(tklabel(ttplot, text="Year(s):"), ap1Entry, sticky="e")
  tkgrid(tklabel(ttplot, text="       "), ap2Entry, sticky="e")
  tkgrid(tklabel(ttplot, text="       "), ap3Entry, sticky="e")
  tkgrid(tklabel(ttplot, text="       "), ap4Entry, sticky="e")
  tkgrid(tklabel(ttplot, text="       "), ap5Entry, sticky="e")
  tkgrid(tklabel(ttplot, text="       "), ap6Entry, sticky="e")

  onOK <- function()
  {
      # tkdestroy(ttplot)
      yearsp <- c( ap1 <- as.numeric(tclvalue(ap1Var)),
         ap2 <- as.numeric(tclvalue(ap2Var)),
         ap3 <- as.numeric(tclvalue(ap3Var)),
         ap4 <- as.numeric(tclvalue(ap4Var)),
         ap5 <- as.numeric(tclvalue(ap5Var)),
         ap6 <- as.numeric(tclvalue(ap6Var)) )
      anyosp <- yearsp[which(yearsp != " ")]    # <<-
      opar <- par(mar=c(8,4.7,5,1.5), ps=18, font=1, tcl=-0.5)
      bbap(.wvari$vari, .wvari$s, .wvari$t0, anyosp)
      par(opar)
  }
  OK.but <- tkbutton(ttplot, text="OK", command=onOK)
  tkgrid(OK.but)
}
#
Makebbcn <- function()
{
  ttplot  <- tktoplevel()
  tkwm.title(ttplot, "Buys-Ballot contour")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  rb1     <- tkradiobutton(ttplot)
  rb2     <- tkradiobutton(ttplot)

  rbValue <- tclVar("Color")
  tkconfigure(rb1, variable=rbValue, value="Blanco y Negro")
  tkconfigure(rb2, variable=rbValue, value="Color")
  tkgrid(tklabel(ttplot, text="Select the representation type:", fg="blue"))
  tkgrid(tklabel(ttplot, text="Black and white"), rb1)
  tkgrid(tklabel(ttplot, text="Colour"), rb2)

  OnOK <- function()
  {
     tkdestroy(ttplot)
     rbVal <- as.character(tclvalue(rbValue))
     if(rbVal == "Blanco y Negro"){ selecolor <- FALSE }
     if(rbVal == "Color"){ selecolor <- TRUE }

     if(.wvari$s==4)
        MR <- quarterg(.wvari$vari, .wvari$s, .wvari$t0, plot=FALSE)
     if(.wvari$s==12)
        MR <- bbmp(.wvari$vari, .wvari$s, .wvari$t0, c(1:.wvari$s), "Prot", plot=FALSE)
     opar <- par(mar=c(6,3.5,4,2))
     bbcn(MR, .wvari$s, .wvari$t0, color=selecolor)
     par(opar)
  }
  OK.but <- tkbutton(ttplot, text="OK", command=OnOK)
  tkgrid(OK.but)
}
#
Makebb3D <- function()
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Buys-Ballot 3D")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  ttplot2 <- tktoplevel()
  tkconfigure(ttplot2, cursor="fleur")
  tkwm.title(ttplot2, "Rotate")
  RightClick2 <- function(x,y)
  {
    rootx <- as.integer(tkwinfo("rootx",ttplot2))
    rooty <- as.integer(tkwinfo("rooty",ttplot2))
    xTxt <- as.integer(x)+rootx
    yTxt <- as.integer(y)+rooty
    if(.wvari$s==4)
       MR <- quarterg(.wvari$vari, .wvari$s, .wvari$t0, plot=FALSE)
    if(.wvari$s==12)
       MR <- bbmp(.wvari$vari, .wvari$s, .wvari$t0, c(1:.wvari$s), "Prot", plot=FALSE)
    bb3D(MR, .wvari$s, .wvari$t0, color=selecolor, as.integer(x), as.integer(y))
  }
  tkbind(ttplot2, "<Button-3>", RightClick2)

  rb1     <- tkradiobutton(ttplot)
  rb2     <- tkradiobutton(ttplot)

  rbValue <- tclVar("Color")
  tkconfigure(rb1, variable=rbValue, value="Blanco y Negro")
  tkconfigure(rb2, variable=rbValue, value="Color")
  tkgrid(tklabel(ttplot, text="Select the representation type:", fg="blue"))
  tkgrid(tklabel(ttplot, text="Black and white"), rb1)
  tkgrid(tklabel(ttplot, text="Colour"), rb2)

  OnOK <- function()
  {
    tkdestroy(ttplot)

    rbVal <- as.character(tclvalue(rbValue))
    if(rbVal == "Blanco y Negro"){ selecolor <<- FALSE }
    if(rbVal == "Color"){ selecolor <<- TRUE }

    if(.wvari$s==4)
       MR <- quarterg(.wvari$vari, .wvari$s, .wvari$t0, plot=FALSE)
    if(.wvari$s==12)
       MR <- bbmp(.wvari$vari, .wvari$s, .wvari$t0, c(1:.wvari$s), "Prot", plot=FALSE)
    opar <- par(mar=c(6,3.5,4,2))
    bb3D(MR, .wvari$s, .wvari$t0, color=selecolor, 30, 30)
    par(opar)
  }
  OK.but <- tkbutton(ttplot, text="OK", command=OnOK)
  tkgrid(OK.but)
}

MakeSeasboxplot <- function()
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Seasonal box plot")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  rb1 <- tkradiobutton(ttplot)
  rb2 <- tkradiobutton(ttplot)

  rbValue <- tclVar("Color")
  tkconfigure(rb1, variable=rbValue, value="Blanco y Negro")
  tkconfigure(rb2, variable=rbValue, value="Color")
  tkgrid(tklabel(ttplot, text="Select the representation type:", fg="blue"))
  tkgrid(tklabel(ttplot, text="Black and white"), rb1)
  tkgrid(tklabel(ttplot, text="Colour"), rb2)

  OnOK <- function()
  {
     tkdestroy(ttplot)
     rbVal <- as.character(tclvalue(rbValue))
     if(rbVal == "Blanco y Negro"){ color <- FALSE }
     if(rbVal == "Color"){ color <- TRUE }
     seasboxplot(.wvari, color)
  }
  OK.but <- tkbutton(ttplot, text="OK", command=OnOK)
  tkgrid(OK.but)
}

#

Makefreqg <- function()
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  if(.wvari$s==12){
       opar <- par(mfrow=c(4,2), mar=c(2,3,3.5,2),
               tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
       freqg(.wvari)
       par(opar)
             }
    if(.wvari$s==4){
      opar <- par(mfrow=c(2,2), mar=c(4,3,4,1.5), font=1,
             tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
      freqg(.wvari)
      par(opar)
            }
}
#
MakeFiltrarfrec <- function()
{
  ttplot <- tktoplevel()
  tkgrid(tklabel(ttplot, text="Select the frequencies you wish to filter:", fg="blue"))

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  cero <- tkcheckbutton(ttplot)
  pi   <- tkcheckbutton(ttplot)
  pi2  <- tkcheckbutton(ttplot)
  if(.wvari$s==12){
    pi23 <- tkcheckbutton(ttplot)
    pi3  <- tkcheckbutton(ttplot)
    pi56 <- tkcheckbutton(ttplot)
    pi6  <- tkcheckbutton(ttplot)
           }
  ceroValue <- tclVar("0")
  piValue   <- tclVar("0")
  pi2Value  <- tclVar("0")
  if(.wvari$s==12){
    pi23Value <- tclVar("0")
    pi3Value  <- tclVar("0")
    pi56Value <- tclVar("0")
    pi6Value  <- tclVar("0")
           }
  tkconfigure(cero, variable=ceroValue)
  tkgrid(tklabel(ttplot, text="cero"), cero)
  tkconfigure(pi, variable=piValue)
  tkgrid(tklabel(ttplot, text="pi"), pi)
  tkconfigure(pi2, variable=pi2Value)
  tkgrid(tklabel(ttplot, text="pi/2"), pi2)
  if(.wvari$s==12){
    tkconfigure(pi23, variable=pi23Value)
    tkgrid(tklabel(ttplot, text="2pi/3"), pi23)
    tkconfigure(pi3, variable=pi3Value)
    tkgrid(tklabel(ttplot, text="pi/3"), pi3)
    tkconfigure(pi56, variable=pi56Value)
    tkgrid(tklabel(ttplot, text="5pi/6"), pi56)
    tkconfigure(pi6, variable=pi6Value)
    tkgrid(tklabel(ttplot, text="pi/6"), pi6)
           }
  OnOK <- function()
  {
     tkdestroy(ttplot)
     ff1 <- as.numeric(tclvalue(ceroValue))
     ff2 <- as.numeric(tclvalue(piValue))
     ff3 <- as.numeric(tclvalue(pi2Value))
     if(.wvari$s==12){
       ff4 <- as.numeric(tclvalue(pi23Value))
       ff5 <- as.numeric(tclvalue(pi3Value))
       ff6 <- as.numeric(tclvalue(pi56Value))
       ff7 <- as.numeric(tclvalue(pi6Value))
       filtra <- c(ff1,ff2,ff3,ff4,ff5,ff6,ff7)
              }
     if(.wvari$s==4){ filtra <- c(ff1,ff2,ff3) }
     Fil.vari <- filtrar(.wvari$vari, .wvari$s, .wvari$t0, filtra, plot=TRUE)[[1]]
     #tkmessageBox(message="Object Fil.vari has been created. It contains filtered series data.", icon="info")
  }
  OK.but <- tkbutton(ttplot, text="OK", command=OnOK)
  tkgrid(OK.but)
}

#

Makecorrgrm <- function()
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Correlograms")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  rb1     <- tkradiobutton(ttplot)
  rb2     <- tkradiobutton(ttplot)
  rb3     <- tkradiobutton(ttplot)
  rb4     <- tkradiobutton(ttplot)
  rbValue <- tclVar("Serie original")

  tkconfigure(rb1, variable=rbValue, value="original")
  tkconfigure(rb2, variable=rbValue, value="delta")
  tkconfigure(rb3, variable=rbValue, value="deltas")
  tkconfigure(rb4, variable=rbValue, value="deltadeltas")
  tkgrid(tklabel(ttplot, text="Select a transformation:", fg="blue"), sticky="w")
  tkgrid(tklabel(ttplot, text="  Original series"), rb1, sticky="w")
  tkgrid(tklabel(ttplot, text="  First differences"), rb2, sticky="w")
  tkgrid(tklabel(ttplot, text="  Seasonal differences"), rb3, sticky="w")
  tkgrid(tklabel(ttplot, text="  First and sesonal differences"), rb4, sticky="w")

  OnOK <- function()
  {
     # tkdestroy(ttplot)
     rbVal <- as.character(tclvalue(rbValue))
     opar <- par(mfrow=c(2,1)) #, mar=c(4,4.7, 2 ,2), ps=18, font=1, tcl=-0.5,
     #           cex.lab=0.8, cex.axis=0.7)
     corrgrm(.wvari, tclvalue(rbValue))
     par(opar)
  }
  OK.but <- tkbutton(ttplot, text="OK", command=OnOK)
  tkgrid(OK.but)
}

#

Makermp <- function()
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Range-mean plot")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  krmg <- tclVar(round(sqrt(.wvari$N)))

  entry.krmg <- tkentry(ttplot, width="3", textvariable=krmg)
  tkgrid(tklabel(ttplot, text="  Introduce the interval width for calculating the points:  \n (By default sqrt(N))", fg="blue"), entry.krmg, rowspan=2)
  tkgrid(tklabel(ttplot))

  rb1     <- tkradiobutton(ttplot)
  rb2     <- tkradiobutton(ttplot)
  lambda  <- tclVar(0)
  rbValue <- tclVar("Serie original")
  tkconfigure(rb1,variable=rbValue, value="Serie original")
  tkconfigure(rb2,variable=rbValue, value="Transformacion Box-Cox")
  tkgrid(tklabel(ttplot, text="  Introduce the value of lambda if you select the transformation:",
       fg="blue"), sticky="w")
  tkgrid(tklabel(ttplot, text="Original series"), rb1)
  tkgrid(tklabel(ttplot, text="Box-Cox transformation"), rb2)

  entry.lambda <- tkentry(ttplot, width="4", textvariable=lambda)
  tkgrid(tklabel(ttplot, text="lambda"), entry.lambda)

  OnOK <- function()
  {
     tkdestroy(ttplot)
     if(tclvalue(rbValue)=="Serie original")
       rmvari <- .wvari$vari
     if(tclvalue(rbValue)=="Transformacion Box-Cox")
     {
       if(tclvalue(lambda) == 0)
         rmvari <- log(.wvari$vari, base=exp(1))
       if(tclvalue(lambda) != 0){
         rmvari <- BoxCox(.wvari$vari, as.numeric(tclvalue(lambda)))
         cat(c("\n Box-Cox parameter: ", tclvalue(lambda), "\n\n"))
#         bcvari <- rmvari <- BoxCox(.wvari$vari, as.numeric(tclvalue(lambda)))
#         tkmessageBox(message="Object bcvari has been created. It contains transformed series #data.", icon="info")
                                }
     }
     # opar <- par(mar=c(8,4.7,5,2), ps=18, font=1,tcl=-0.5)
     #opar <- par(mar=c(4,4.7, 4 ,2), ps=18, font=1, tcl=-0.5, cex.lab=0.8, cex.axis=0.7)
     rmg(rmvari, as.numeric(tclvalue(krmg)))
     #par(opar)
  }
  OK.but <- tkbutton(ttplot, text="OK", command=OnOK)
  tkgrid(OK.but)
}

#

Makevfic <- function()
{
  tt    <- tktoplevel()
  tkwm.title(tt, "Generic dummy")

  y0fic <- tclVar()
  s0fic <- tclVar()
  yNfic <- tclVar()
  sNfic <- tclVar()

  entry.y0fic <- tkentry(tt, width="5", textvariable=y0fic)
  entry.s0fic <- tkentry(tt, width="3", textvariable=s0fic)
  entry.yNfic <- tkentry(tt, width="5", textvariable=yNfic)
  entry.sNfic <- tkentry(tt, width="3", textvariable=sNfic)

  tkgrid(tklabel(tt, text="  Specify a dummy:", fg="blue"), sticky="w")
  tkgrid(tklabel(tt, text="    Start (Year-season):"), entry.y0fic, entry.s0fic)
  tkgrid(tklabel(tt, text="    End   (Year-season):"), entry.yNfic, entry.sNfic)

  OnOK <- function()
  {
     tkdestroy(tt)
     n <- which.min(c(which(Mvfic[1,] == 0), which(Mvfic[1,] == 1)))
     ifelse(length(n)==0, n<-0, n<-n)
     Mvfic[,(n+1)] <<- vfic(as.numeric(tclvalue(y0fic)),
                    as.numeric(tclvalue(s0fic)),
                    as.numeric(tclvalue(yNfic)),
                    as.numeric(tclvalue(sNfic)))
     tkmessageBox(message="This dummy has been incorporated to the Mvfic object.", icon="info")
  }
  OK.but <- tkbutton(tt, text="OK", command=OnOK)
  tkgrid(OK.but)
}

#

MakeVFEp <- function()
{
  tt <- tktoplevel()
  tkwm.title(tt, "Partial seasonal dummy")

  vfe1 <- tclVar()
  vfe2 <- tclVar()
  vfe3 <- tclVar()
  vfe4 <- tclVar()
  if(.wvari$s==12){
  vfe5 <- tclVar()
  vfe6 <- tclVar()
                  }
  entry.vfe1 <- tkentry(tt, width="2", textvariable=vfe1)
  entry.vfe2 <- tkentry(tt, width="2", textvariable=vfe2)
  entry.vfe3 <- tkentry(tt, width="2", textvariable=vfe3)
  entry.vfe4 <- tkentry(tt, width="2", textvariable=vfe4)
  if(.wvari$s==12){
    entry.vfe5 <- tkentry(tt, width="2", textvariable=vfe5)
    entry.vfe6 <- tkentry(tt, width="2", textvariable=vfe6)
                  }
  tkgrid(tklabel(tt, text="Introduce the seasons you wish to consider:", fg="blue"), sticky="w")
  if(.wvari$s==4)
     tkgrid(tklabel(tt, text="Seasons:"), entry.vfe1, entry.vfe2, entry.vfe3, entry.vfe4)
  if(.wvari$s==12)
     tkgrid(tklabel(tt, text="Seasons:"),
            entry.vfe1, entry.vfe2, entry.vfe3, entry.vfe4, entry.vfe5, entry.vfe6)

  OnOK <- function()
  {
    tkdestroy(tt)
    vfe1 <- as.numeric(tclvalue(vfe1))
    vfe2 <- as.numeric(tclvalue(vfe2))
    vfe3 <- as.numeric(tclvalue(vfe3))
    vfe4 <- as.numeric(tclvalue(vfe4))
    if(.wvari$s==12){
    vfe5 <- as.numeric(tclvalue(vfe5))
    vfe6 <- as.numeric(tclvalue(vfe6))
             }
    ifelse(.wvari$s==12, aux <- c(vfe1, vfe2, vfe3, vfe4, vfe5, vfe6),
                  aux <- c(vfe1, vfe2, vfe3, vfe4))
    VFEp <<- MVFE(.wvari$vari, .wvari$s, .wvari$t0, "alg")[,na.omit(aux)]
    tkmessageBox(message="The dummy has been saved as the object VFEp.", icon="info")
  }
  OK.but <- tkbutton(tt, text="OK", command=OnOK)
  tkgrid(OK.but)
}

SalidaLaTeXfrec <- function()  # poner .tex al elegir archivo
{
  outfile <- tclvalue(tkgetSaveFile(
    filetypes='{"Text files" {".tex"}} {"All Files" {"*"}}'))

  if(!nchar(outfile))
     tkmessageBox(message="No file was chosen.", icon="error")
  else{
     # outfile <- chartr("/", "\\", outfile)
         # Así no cambia de directorio, lo guarda en el actual
     tkconfigure(.treeWidget, cursor="watch")
     TablaFrec(outfile)
     tkconfigure(.treeWidget, cursor="xterm")
     tkmessageBox(message=paste("Table has been saved in", outfile), icon="info")
      }
}
#
SalidaLaTeXdet <- function()  # poner .tex al elegir archivo
{
  outfile <- tclvalue(tkgetSaveFile(
    filetypes='{"Text files" {".tex"}} {"All Files" {"*"}}'))

  if(!nchar(outfile))
     tkmessageBox(message="No file was chosen.", icon="error")
  else{
     tkconfigure(.treeWidget, cursor="watch")
     TablaDet(outfile)
     tkconfigure(.treeWidget, cursor="xterm")
     tkmessageBox(message=paste("Table has been saved in", outfile), icon="info")
      }
}

#

MakePanelqmCorrg <- function()
{
   ttplot <- tktoplevel()
   tkwm.title(ttplot, "Panel")
   tklabel(ttplot, text="Panel")

   string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
   .wvari <<- ExeString(c(string, ""))

   tkdestroy(ttplot)
   opar <- par(mfrow=c(4,2), mar=c(2,3,3.5,2), # mar=c(3,3,3.5,2),
               tcl=-0.5, cex.axis=0.7, las=1)
   corrgrm(.wvari, "original")
   corrgrm(.wvari, "delta")
   corrgrm(.wvari, "deltas")
   corrgrm(.wvari, "deltadeltas")
   par(opar)
}
#
MakePanelmFreqg <- function()
{
   ttplot <- tktoplevel()
   tkwm.title(ttplot, "Panel")
   tklabel(ttplot, text="Panel")

   string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
   .wvari <<- ExeString(c(string, ""))

   tkdestroy(ttplot)
   opar <- par(mfrow=c(4,2), mar=c(2,3,3.5,2),
              tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
   freqg(.wvari)
   par(opar)
}
#
MakePanelmBB1 <- function()
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Panel")
  tklabel(ttplot, text="Panel")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  tkdestroy(ttplot)
  opar <- par(mfrow=c(2,2), mar=c(2,3,3.5,2), tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
  bbmp(.wvari$vari, .wvari$s, .wvari$t0, c(1:.wvari$s), "Prot", plot=TRUE)
  # bbmp((ret(vari,2)[,1]-ret(vari,2)[,2]), .wvari$s, .wvari$t0,
  #     , c(1:.wvari$s), "Prot", plot=TRUE)
  par(opar)
}
#
MakePanelmBB2 <- function()
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Panel")
  tklabel(ttplot, text="Panel")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  tkdestroy(ttplot)
  opar <- par(mfrow=c(2,2), mar=c(2,3,3.5,2), tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
  bbmp(ret(.wvari$vari,2)[,1]-ret(.wvari$vari,2)[,2],
       .wvari$s, .wvari$t0, c(1:.wvari$s), "Prot", plot=TRUE)
  par(opar)
}
#
MakePanelmSerie <- function()
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Panel")
  tklabel(ttplot, text="Panel")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  rb1     <- tkradiobutton(ttplot)
  rb2     <- tkradiobutton(ttplot)
  rbValue <- tclVar("log")
  tkconfigure(rb1, variable=rbValue, value="log")
  tkconfigure(rb2, variable=rbValue, value="nlog")

  tkgrid(tklabel(ttplot, text="  Indicate the scale of the series:", fg="blue"), sticky="w")
  tkgrid(tklabel(ttplot, text="Logarithm"), rb1)
  tkgrid(tklabel(ttplot, text="Without logarithm"), rb2)

  OnOK <- function()
  {
    tkdestroy(ttplot)
    if(tclvalue(rbValue) == "log"){
      variaux <- log(.wvari$vari, base=exp(1))
      main1 <- "First differences of the logarithms"
      main2 <- "Seasonal and regular differences of the logarithms"
      main3 <- "Estimated espectral density upon the logarithms"
      main4 <- "Spectrum of the first differences in logaritmhs"
                                  }
    if(tclvalue(rbValue) == "nlog"){
      variaux <- .wvari$vari
      main1 <- "First differences"
      main2 <- "Seasonal and regular differences"
      main3 <- "Estimated spectral density"
      main4 <- "Spectrum of the first differences"
                                   }
    opar <- par(mfrow=c(4,2), mar=c(2,3,3.5,2), tcl=-0.5, cex.axis=0.7, cex.main=1, las=1)
    plot(variaux, main=.wvari$label)
    rmg(variaux, round(sqrt(.wvari$N)))
    plot(log(variaux, base=exp(1)), main="Logarithms of the series")
    rmg(log(variaux, base=exp(1)), round(sqrt(.wvari$N)))
    plot(diff(variaux, lag=1), main=main1)
    plot(diff(diff(variaux, lag=.wvari$s), lag=1), main=main2)
    spectrum(variaux, spans=c(3,5), log="no", ylab="", xlab="frequency", main=main3)
    spectrum(diff(variaux, lag=1), spans=c(3,5), log="no", ylab="", xlab="frequency", main=main4)
    par(opar)
  }
  OK.but <- tkbutton(ttplot, text="OK", command=OnOK)
  tkgrid(OK.but)
}
#
MakeQPanel1 <- function()  # cambiar nombre en GUI menu
{
  ttplot <- tktoplevel()
  tkwm.title(ttplot, "Panel")
  tklabel(ttplot, text="Panel")

  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))

  rb1     <- tkradiobutton(ttplot)
  rb2     <- tkradiobutton(ttplot)
  rbValue <- tclVar("log")
  tkconfigure(rb1, variable=rbValue, value="log")
  tkconfigure(rb2, variable=rbValue, value="nlog")

  tkgrid(tklabel(ttplot, text="  Indicate the scale of the series:", fg="blue"), sticky="w")
  tkgrid(tklabel(ttplot, text="Logarithm"), rb1)
  tkgrid(tklabel(ttplot, text="Without logarithm"), rb2)

  OnOK <- function()
  {
    tkdestroy(ttplot)
    if(tclvalue(rbValue) == "log"){
      variaux <- log(.wvari$vari, base=exp(1))
      main1 <- "First differences of the logarithms"
      main2 <- "Seasonal and regular differences of the logarithms"
      main3 <- "Estimated espectral density upon the logarithms"
      main4 <- "Spectrum of the first differences in logaritmhs"
                                  }
    if(tclvalue(rbValue) == "nlog"){
      variaux <- .wvari$vari
      main1 <- "First differences"
      main2 <- "Seasonal and regular differences"
      main3 <- "Estimated spectral density"
      main4 <- "Spectrum of the first differences"
                                   }
    opar <- par(mfrow=c(4,2), mar=c(2,3,3.5,2), tcl=-0.5, cex.axis=0.8, cex.main=1, las=1)
    plot(variaux, main=.wvari$label)
    rmg(variaux, round(sqrt(.wvari$N)))
    plot(log(variaux, base=exp(1)), main="Logarithms of the series")
    rmg(log(variaux, base=exp(1)), round(sqrt(.wvari$N)))
    plot(diff(variaux, lag=1), main=main1)
    plot(diff(diff(variaux, lag=.wvari$s), lag=1), main=main2)
    spectrum(variaux, spans=c(3,5), log="no", ylab="", xlab="frequency", main=main3)
    spectrum(diff(variaux, lag=1), spans=c(3,5), log="no", ylab="", xlab="frequency", main=main4)
    par(opar)
  }
  OK.but <- tkbutton(ttplot, text="OK", command=OnOK)
  tkgrid(OK.but)
}
#
MakeQPanelBB <- function()
{
   string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
   .wvari <<- ExeString(c(string, ""))

   opar <- par(mfrow=c(1,2), mar=c(2,2.5,3.5,1.5), tcl=-0.5, cex.axis=0.8, cex.main=1, las=1)
   MQ <- quarterg(.wvari$vari, .wvari$s, .wvari$t0, plot=TRUE)
   quarterg((ret(.wvari$vari,2)[,1]-ret(.wvari$vari,2)[,2]), .wvari$s, .wvari$t0, plot=TRUE)
   par(opar)
}


# DATA BANK FUNCTIONS

dataIKERBIDE <- function()
{
  ttbd <- tktoplevel()
  tkwm.title(ttbd, "CAPV data bank")

  airp    <- tkradiobutton(ttbd)
  epaact  <- tkradiobutton(ttbd)
  iai     <- tkradiobutton(ttbd)
  ipc     <- tkradiobutton(ttbd)
  ipi     <- tkradiobutton(ttbd)
  ipri    <- tkradiobutton(ttbd)
  licof   <- tkradiobutton(ttbd)
  mtur    <- tkradiobutton(ttbd)
  mvic    <- tkradiobutton(ttbd)
  paroreg <- tkradiobutton(ttbd)
  pcemen  <- tkradiobutton(ttbd)
  pernhot <- tkradiobutton(ttbd)
  valorpr <- tkradiobutton(ttbd)
  vpoinic <- tkradiobutton(ttbd)
  vpoterm <- tkradiobutton(ttbd)

  serieValue <- tclVar(" ")
  tkconfigure(airp, variable=serieValue, value="2")
  tkconfigure(epaact, variable=serieValue, value="3")
  tkconfigure(iai, variable=serieValue, value="4")
  tkconfigure(ipc, variable=serieValue, value="5")
  tkconfigure(ipi, variable=serieValue, value="6")
  tkconfigure(ipri, variable=serieValue, value="7")
  tkconfigure(licof, variable=serieValue, value="8")
  tkconfigure(mtur, variable=serieValue, value="9")
  tkconfigure(mvic, variable=serieValue, value="10")
  tkconfigure(paroreg, variable=serieValue, value="11")
  tkconfigure(pcemen, variable=serieValue, value="12")
  tkconfigure(pernhot, variable=serieValue, value="13")
  tkconfigure(valorpr, variable=serieValue, value="14")
  tkconfigure(vpoinic, variable=serieValue, value="15")
  tkconfigure(vpoterm, variable=serieValue, value="16")

  tkgrid(tklabel(ttbd, text="Select a series:", fg="blue"))
  tkgrid(tklabel(ttbd, text="  Airpassengers"), airp, sticky="w")
  tkgrid(tklabel(ttbd, text="  Actives"), epaact, sticky="w")
  tkgrid(tklabel(ttbd, text="  Industrial activity index"), iai, sticky="w")
  tkgrid(tklabel(ttbd, text="  Consumer price index"), ipc, sticky="w")
  tkgrid(tklabel(ttbd, text="  Industrial production index"), ipi, sticky="w")
  tkgrid(tklabel(ttbd, text="  Industrial price index"), ipri, sticky="w")
  tkgrid(tklabel(ttbd, text="  Official bidding"), licof, sticky="w")
  tkgrid(tklabel(ttbd, text="  Private car registration"), mtur, sticky="w")
  tkgrid(tklabel(ttbd, text="  Industrial cargo vehicles registration"), mvic, sticky="w")
  tkgrid(tklabel(ttbd, text="  Registered unemployment"), paroreg, sticky="w")
  tkgrid(tklabel(ttbd, text="  Cement production"), pcemen, sticky="w")
  tkgrid(tklabel(ttbd, text="  Hotel occupation"), pernhot, sticky="w")
  tkgrid(tklabel(ttbd, text="  Production value"), valorpr, sticky="w")
  tkgrid(tklabel(ttbd, text="  Initiated Council houses"), vpoinic, sticky="w")
  tkgrid(tklabel(ttbd, text="  Finished Council houses"), vpoterm, sticky="w")

  Descrbut <- function()
  {
     mytkpager(file.path(R.home(), "library/uroot/data/source_capv"),
       title="CAPV Data Bank", header="", delete.file=FALSE, wwidth=70, wheight=20, export=FALSE)
    # browseURL(file.path(R.home(), "library/urootcook/data/source_capv.html"),
               # browser=getOption("browser"))
  }
  Descr.but <- tkbutton(ttbd, text="Source", command=Descrbut)
  OnOK <- function()
  {
     if(tclvalue(serieValue)==2){
        s<-12; t0 <- c(1982,1); logvari <- FALSE}
     if(tclvalue(serieValue)==3){
        s<-4; t0<- c(1976,3); logvari <- FALSE}
     if(tclvalue(serieValue)==4){
        s<-12; t0<- c(1984,7); logvari <- FALSE}
     if(tclvalue(serieValue)==5){
        s<-12; t0<- c(1993,1); logvari <- FALSE}
     if(tclvalue(serieValue)==6){
        s<-12; t0<- c(1990,1); logvari <- FALSE}
     if(tclvalue(serieValue)==7){
        s<-12; t0<- c(1993,1); logvari <- FALSE}
     if(tclvalue(serieValue)==8){
        s<-4; t0<- c(1985,1); logvari <- FALSE}
     if(tclvalue(serieValue)==9){
        s<-12; t0<- c(1980,1); logvari <- FALSE}
     if(tclvalue(serieValue)==10){
        s<-12; t0<- c(1980,1); logvari <- FALSE}
     if(tclvalue(serieValue)==11){
        s<-12; t0<- c(1980,1); logvari <- FALSE}
     if(tclvalue(serieValue)==12){
        s<-12; t0<- c(1980,1); logvari <- FALSE}
     if(tclvalue(serieValue)==13){
        s<-12; t0<- c(1981,1); logvari <- FALSE}
     if(tclvalue(serieValue)==14){
        s<-4; t0<- c(1989,1); logvari <- FALSE}
     if(tclvalue(serieValue)==15){
        s<-12; t0<- c(1982,1); logvari <- FALSE}
     if(tclvalue(serieValue)==16){
        s<-12; t0<- c(1982,1); logvari <- FALSE}

     datos <- read.csv(file.path(R.home(), "library/uroot/data/dataIKERBIDE.csv"),
                       header=TRUE, sep=" ", dec=",", na.strings="NA")
     vari <- datos[,as.numeric(tclvalue(serieValue))]
     aux1 <- sort(c(which(vari>=0), which(vari<0)))
     aux2 <- c(1:(length(aux1)-1))
     for(i in 2:length(aux1))
       aux2[i] <- aux1[i]-aux1[i-1]
     if(length(which(aux2 > 1)) == 0)
       N <- aux1[length(aux1)]
     if(length(which(aux2 > 1)) > 0){
       if(length(which(aux2 == 8)) == 0)
         N <- length(vari)
       if(length(which(aux2 == 8)) > 0)
         N <- which(aux2 == 8)[1] - 1
                                    }

     vari <- ts(vari[1:N], frequency=s, start=t0)
     labels <- c("airp", "epaact", "iai", "ipc", "ipi", "ipri", "licof", "mtur", "mvic",
                 "paroreg", "pcemen", "pernhot", "valorpr", "vpoinic", "vpoterm")
     label <- as.character(labels[as.numeric(tclvalue(serieValue))-1])
     assign(label, list(vari=vari, s=s, t0=t0, N=N, logvari=logvari, label=label), env=.GlobalEnv)
     assign(".wvari", list(vari=vari, s=s, t0=t0, N=N, logvari=logvari, label=label),
            env=.GlobalEnv)
     # tclvalue(.label) <<- .wvari$label

     tkinsert(.treeWidget,"end","root",label,text=label)

     msg <- paste("Information about the series has been stored in the object ", label, sep="")
     tkmessageBox(title="Series info", message=msg, icon="info")
     tkdestroy(ttbd)
  }
  OK.but <- tkbutton(ttbd, text="OK", command=OnOK)
  tkgrid(Descr.but)
  tkgrid(OK.but)
}

#

dataINE <- function()
{
  ttbd <- tktoplevel()
  tkwm.title(ttbd, "INE data bank")

  costesal <- tkradiobutton(ttbd)
  indcsal  <- tkradiobutton(ttbd)
  ipc      <- tkradiobutton(ttbd)
  ipi      <- tkradiobutton(ttbd)
  ipri     <- tkradiobutton(ttbd)
  mtur     <- tkradiobutton(ttbd)
  mvic     <- tkradiobutton(ttbd)
  paroreg  <- tkradiobutton(ttbd)
  pernhot  <- tkradiobutton(ttbd)
  valojhot <- tkradiobutton(ttbd)
  vpoinic  <- tkradiobutton(ttbd)
  vpoterm  <- tkradiobutton(ttbd)
  epaact   <- tkradiobutton(ttbd)

  serieValue <- tclVar(" ")
  tkconfigure(costesal, variable=serieValue, value="2")
  tkconfigure(indcsal, variable=serieValue, value="3")
  tkconfigure(ipc, variable=serieValue, value="4")
  tkconfigure(ipi, variable=serieValue, value="5")
  tkconfigure(ipri, variable=serieValue, value="6")
  tkconfigure(mtur, variable=serieValue, value="7")
  tkconfigure(mvic, variable=serieValue, value="8")
  tkconfigure(paroreg, variable=serieValue, value="9")
  tkconfigure(pernhot, variable=serieValue, value="10")
  tkconfigure(valojhot, variable=serieValue, value="11")
  tkconfigure(vpoinic, variable=serieValue, value="12")
  tkconfigure(vpoterm, variable=serieValue, value="13")
  tkconfigure(epaact, variable=serieValue, value="14")

  tkgrid(tklabel(ttbd, text="Select a series:", fg="blue"))
  tkgrid(tklabel(ttbd, text="  Wage cost"), costesal, sticky="w")
  tkgrid(tklabel(ttbd, text="  Wage cost index"), indcsal, sticky="w")
  tkgrid(tklabel(ttbd, text="  Consumer price index"), ipc, sticky="w")
  tkgrid(tklabel(ttbd, text="  Industrial production index"), ipi, sticky="w")
  tkgrid(tklabel(ttbd, text="  Industrial price index"), ipri, sticky="w")
  tkgrid(tklabel(ttbd, text="  Private car registration"), mtur, sticky="w")
  tkgrid(tklabel(ttbd, text="  Industrial cargo vehicles registration"), mvic, sticky="w")
  tkgrid(tklabel(ttbd, text="  Registered unemployment"), paroreg, sticky="w")
  tkgrid(tklabel(ttbd, text="  Hotel occupation"), pernhot, sticky="w")
  tkgrid(tklabel(ttbd, text="  Travellers lodged in hotels"), valojhot, sticky="w")
  tkgrid(tklabel(ttbd, text="  Initiated council houses"), vpoinic, sticky="w")
  tkgrid(tklabel(ttbd, text="  Finished council houses"), vpoterm, sticky="w")
  tkgrid(tklabel(ttbd, text="  Actives"), epaact, sticky="w")

  Descrbut <- function()
  {
    mytkpager(file.path(R.home(), "library/uroot/data/source_es"),
       title="INE Data Bank", header="", delete.file=FALSE, wwidth=60, wheight=10, export=FALSE)
    #browseURL(file.path(R.home(), "library/urootcook/data/source_es.html"),
      #browser=getOption("browser"))
  }
  Descr.but <- tkbutton(ttbd, text="Source", command=Descrbut)
  OnOK <- function()
  {
     if(tclvalue(serieValue)==2){
        s<-4; t0<- c(1982,1); logvari <- FALSE}
     if(tclvalue(serieValue)==3){
        s<-4; t0<- c(1982,1); logvari <- FALSE}
     if(tclvalue(serieValue)==4){
        s<-12; t0<- c(1975,1); logvari <- FALSE}
     if(tclvalue(serieValue)==5){
        s<-12; t0<- c(1975,1); logvari <- FALSE}
     if(tclvalue(serieValue)==6){
        s<-12; t0<- c(1975,1); logvari <- FALSE}
     if(tclvalue(serieValue)==7){
        s<-12; t0<- c(1983,1); logvari <- FALSE}
     if(tclvalue(serieValue)==8){
        s<-12; t0<- c(1983,1); logvari <- FALSE}
     if(tclvalue(serieValue)==9){
        s<-12; t0<- c(1977,1); logvari <- FALSE}
     if(tclvalue(serieValue)==10){
        s<-12; t0<- c(1987,1); logvari <- FALSE}
     if(tclvalue(serieValue)==11){
        s<-12; t0<- c(1982,1); logvari <- FALSE}
     if(tclvalue(serieValue)==12){
        s<-12; t0<- c(1972,1); logvari <- FALSE}
     if(tclvalue(serieValue)==13){
        s<-12; t0<- c(1970,1); logvari <- FALSE}
     if(tclvalue(serieValue)==14){
        s<-4; t0<- c(1987,2); logvari <- FALSE}

     datos <- read.csv(file.path(R.home(), "library/uroot/data/dataINE.csv"),
                       header=TRUE, sep=" ", dec=",", na.strings="NA")
     vari <- datos[,as.numeric(tclvalue(serieValue))]

     aux1 <- sort(c(which(vari>=0), which(vari<0)))
     aux2 <- c(1:(length(aux1)-1))
     for(i in 2:length(aux1))
       aux2[i] <- aux1[i]-aux1[i-1]
     if(length(which(aux2 > 1)) == 0)
       N <- aux1[length(aux1)]
     if(length(which(aux2 > 1)) > 0){
       if(length(which(aux2 == 8)) == 0)
         N <- length(vari)
       if(length(which(aux2 == 8)) > 0)
         N <- which(aux2 == 8)[1] - 1
                                    }
     vari <- ts(vari[1:N], frequency=s, start=t0)
     labels <- c("costesal", "indcsal", "ipc", "ipi", "ipri", "mtur", "mvic", "paroreg",
                  "pernhot", "valojhot", "vpoinic", "vpoterm", "epaact")
     label <- as.character(labels[as.numeric(tclvalue(serieValue))-1])
     assign(label, list(vari=vari, s=s, t0=t0, N=N, logvari=logvari, label=label), env=.GlobalEnv)
     assign(".wvari", list(vari=vari, s=s, t0=t0, N=N, logvari=logvari, label=label),
            env=.GlobalEnv)

     tkinsert(.treeWidget,"end","root",label,text=label)

     msg <- paste("Information about the series has been stored in the object ", label, sep="")
     tkmessageBox(title="Series info", message=msg, icon="info")
     tkdestroy(ttbd)
  }
  OK.but <- tkbutton(ttbd, text="OK", command=OnOK)
  tkgrid(Descr.but)
  tkgrid(OK.but)
}

#

dataFranses <- function()
{
  ttbd <- tktoplevel()
  tkwm.title(ttbd, "Data bank")

  usaipi    <- tkradiobutton(ttbd)
  canun     <- tkradiobutton(ttbd)
  gergnp    <- tkradiobutton(ttbd)
  ukinvest  <- tkradiobutton(ttbd)
  usaipisa  <- tkradiobutton(ttbd)
  canunsa   <- tkradiobutton(ttbd)
  gergnpsa  <- tkradiobutton(ttbd)
  ukgdp     <- tkradiobutton(ttbd)
  ukcons    <- tkradiobutton(ttbd)
  ukndcons  <- tkradiobutton(ttbd)
  ukexp     <- tkradiobutton(ttbd)
  ukimp     <- tkradiobutton(ttbd)
  ukpinvest <- tkradiobutton(ttbd)
  ukwf      <- tkradiobutton(ttbd)
  swndcpc   <- tkradiobutton(ttbd)
  swdipc    <- tkradiobutton(ttbd)

  serieValue <- tclVar(" ")
  tkconfigure(usaipi, variable=serieValue, value="2")
  tkconfigure(canun, variable=serieValue, value="3")
  tkconfigure(gergnp, variable=serieValue, value="4")
  tkconfigure(ukinvest, variable=serieValue, value="5")
  tkconfigure(usaipisa, variable=serieValue, value="6")
  tkconfigure(canunsa, variable=serieValue, value="7")
  tkconfigure(gergnpsa, variable=serieValue, value="8")
  tkconfigure(ukgdp, variable=serieValue, value="9")
  tkconfigure(ukcons, variable=serieValue, value="10")
  tkconfigure(ukndcons, variable=serieValue, value="11")
  tkconfigure(ukexp, variable=serieValue, value="12")
  tkconfigure(ukimp, variable=serieValue, value="13")
  tkconfigure(ukpinvest, variable=serieValue, value="14")
  tkconfigure(ukwf, variable=serieValue, value="15")
  tkconfigure(swndcpc, variable=serieValue, value="16")
  tkconfigure(swdipc, variable=serieValue, value="17")

  tkgrid(tklabel(ttbd, text="Select a series:", fg="blue"))
  tkgrid(tklabel(ttbd, text="  Total Industrial Production Index for the United States"), usaipi,
         sticky="w")
  tkgrid(tklabel(ttbd, text="  Unemployment in Canada"), canun, sticky="w")
  tkgrid(tklabel(ttbd, text="  Real GNP in Germany"), gergnp, sticky="w")
  tkgrid(tklabel(ttbd, text="  Real Total Investment in the United Kindom"), ukinvest, sticky="w")
  tkgrid(tklabel(ttbd, text="  Total Industrial Production Index for the United States (seasonally adjusted)"), usaipisa, sticky="w")
  tkgrid(tklabel(ttbd, text="  Unemployment in Canada (seasonally adjusted)"), canunsa, sticky="w")
  tkgrid(tklabel(ttbd, text="  Real GNP in Germany (seasonally adjusted)"), gergnpsa, sticky="w")
  tkgrid(tklabel(ttbd, text="  United Kingdom gross domestic product"), ukgdp, sticky="w")
  tkgrid(tklabel(ttbd, text="  United Kingdom total consumption"), ukcons, sticky="w")
  tkgrid(tklabel(ttbd, text="  United Kindom nondurables consumption"), ukndcons, sticky="w")
  tkgrid(tklabel(ttbd, text="  United Kindom exports of goods and services"), ukexp, sticky="w")
  tkgrid(tklabel(ttbd, text="  United Kindom imports of goods and services"), ukimp, sticky="w")
  tkgrid(tklabel(ttbd, text="  United Kindom public investment"), ukpinvest, sticky="w")
  tkgrid(tklabel(ttbd, text="  United Kindom workforce"), ukwf, sticky="w")
  tkgrid(tklabel(ttbd, text="  Real per capita non-durables consumption in Sweden (measured in logs)"), swndcpc, sticky="w")
  tkgrid(tklabel(ttbd, text="  Real per capita disposable income in Sweden (measured in logs)"),
         swdipc, sticky="w")

  Descrbut <- function()
  {
    mytkpager(file.path(R.home(), "library/uroot/data/source_franses"),
       title="Franses Data Bank", header="", delete.file=FALSE, wwidth=95, wheight=40, export=FALSE)
    #browseURL(file.path(R.home(), "library/urootcook/data/source_franses.html"),
               #browser=getOption("browser"))
  }
  Descr.but <- tkbutton(ttbd, text="Source", command=Descrbut)
  OnOK <- function()
  {
     if(tclvalue(serieValue)==2){
        s<-4; t0<- c(1960,1)}
     if(tclvalue(serieValue)==3){
        s<-4; t0<- c(1960,1)}
     if(tclvalue(serieValue)==4){
        s<-4; t0<- c(1960,1)}
     if(tclvalue(serieValue)==5){
        s<-4; t0<- c(1955,1)}
     if(tclvalue(serieValue)==6){
        s<-4; t0<- c(1960,1)}
     if(tclvalue(serieValue)==7){
        s<-4; t0<- c(1960,1)}
     if(tclvalue(serieValue)==8){
        s<-4; t0<- c(1960,1)}
     if(tclvalue(serieValue)==9){
        s<-4; t0<- c(1955,1)}
     if(tclvalue(serieValue)==10){
        s<-4; t0<- c(1955,1)}
     if(tclvalue(serieValue)==11){
        s<-4; t0<- c(1955,1)}
     if(tclvalue(serieValue)==12){
        s<-4; t0<- c(1955,1)}
     if(tclvalue(serieValue)==13){
        s<-4; t0<- c(1955,1)}
     if(tclvalue(serieValue)==14){
        s<-4; t0<- c(1962,1)}
     if(tclvalue(serieValue)==15){
        s<-4; t0<- c(1955,1)}
     if(tclvalue(serieValue)==16){
        s<-4; t0<- c(1963,1)}
     if(tclvalue(serieValue)==17){
        s<-4; t0<- c(1963,1)}

     datos <- read.csv(file.path(R.home(), "library/uroot/data/dataFranses.csv"),
                       header=TRUE, sep=" ", dec=",", na.strings="NA")
     vari <- datos[,as.numeric(tclvalue(serieValue))]

     aux1 <- sort(c(which(vari>=0), which(vari<0)))
     aux2 <- c(1:(length(aux1)-1))
     for(i in 2:length(aux1))
       aux2[i] <- aux1[i]-aux1[i-1]
     if(length(which(aux2 > 1)) == 0)
       N <- aux1[length(aux1)]
     if(length(which(aux2 > 1)) > 0){
       if(length(which(aux2 == 8)) == 0)
         N <- length(vari)
       if(length(which(aux2 == 8)) > 0)
         N <- which(aux2 == 8)[1] - 1
                                    }

     vari <- ts(vari[1:N], frequency=s, start=t0)
     labels <-  c("usaipi", "canun", "gergnp", "ukinvest", "usaipisa", "canunsa", "gergnpsa",
                   "ukgdp", "ukcons", "ukndcons", "ukexp", "ukimp", "ukpinvest", "ukwf",
                   "swndcpc", "swdipc")
     label <- as.character(labels[as.numeric(tclvalue(serieValue))-1])
     logvaris <- c("FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE",
                   "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE")
     logvari <- as.character(logvaris[as.numeric(tclvalue(serieValue))-1])

     assign(label, list(vari=vari, s=s, t0=t0, N=N, logvari=logvari, label=label), env=.GlobalEnv)
     assign(".wvari", list(vari=vari, s=s, t0=t0, N=N, logvari=logvari, label=label),
            env=.GlobalEnv)

     tkinsert(.treeWidget,"end","root",label,text=label)

     msg <- paste("Information about the series has been stored in the object ", label, sep="")
     tkmessageBox(title="Series info", message=msg, icon="info")
     tkdestroy(ttbd)
  }
  OK.but <- tkbutton(ttbd, text="OK", command=OnOK)
  tkgrid(Descr.but)
  tkgrid(OK.but)
}


# Funciones para hacer los contrastes de forma recursiva (diferentes periodos muestrales)

# Funciones auxiliares (principalmente para recursivetest)
# leyenda code: ejemplo: BM000t0: tablas de BM, compdet=c(0,0,0), frecuencia=t0

lookuptablas <- function(code)
{
  if(code=="BM000t0")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-2.51, -2.52, -2.57),
                         CV2 = c(-2.18, -2.21, -2.24),
                         CV3 = c(-1.89, -1.91, -1.95),
                         CV4 = c(-1.58, -1.59, -1.62))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="BM100t0")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-3.35, -3.40, -3.41),
                         CV2 = c(-3.06, -3.11, -3.12),
                         CV3 = c(-2.80, -2.85, -2.86),
                         CV4 = c(-2.51, -2.55, -2.57))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="BM101t0")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-3.32, -3.37, -3.41),
                         CV2 = c(-3.02, -3.06, -3.12),
                         CV3 = c(-2.76, -2.81, -2.86),
                         CV4 = c(-2.47, -2.53, -2.57))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="BM110t0")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-3.87, -3.92, -3.97),
                         CV2 = c(-3.58, -3.63, -3.67),
                         CV3 = c(-3.32, -3.37, -3.40),
                         CV4 = c(-3.06, -3.09, -3.12))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="BM111t0")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-3.83, -3.85, -3.97),
                         CV2 = c(-3.54, -3.57, -3.67),
                         CV3 = c(-3.28, -3.32, -3.40),
                         CV4 = c(-2.99, -3.04, -3.12))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
#
  if(code=="BM000tpi")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-2.53, -2.52, -2.57),
                         CV2 = c(-2.16, -2.20, -2.24),
                         CV3 = c(-1.87, -1.91, -1.95),
                         CV4 = c(-1.57, -1.59, -1.62))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="BM100tpi")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-2.48, -2.54, -2.57),
                         CV2 = c(-2.15, -2.20, -2.24),
                         CV3 = c(-1.89, -1.91, -1.95),
                         CV4 = c(-1.57, -1.59, -1.62))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="BM101tpi")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-3.28, -3.37, -3.41),
                         CV2 = c(-3.01, -3.07, -3.12),
                         CV3 = c(-2.76, -2.81, -2.86),
                         CV4 = c(-2.48, -2.52, -2.57))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="BM110tpi")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-2.52, -2.55, -2.57),
                         CV2 = c(-2.18, -2.20, -2.24),
                         CV3 = c(-1.88, -1.93, -1.95),
                         CV4 = c(-1.55, -1.60, -1.62))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="BM111tpi")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(-3.31, -3.40, -3.41),
                         CV2 = c(-3.02, -3.08, -3.12),
                         CV3 = c(-2.75, -2.84, -2.86),
                         CV4 = c(-2.47, -2.54, -2.57))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
#
  if(code=="BM000Foddeven")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(2.34, 2.38, 2.40),
                         CV2 = c(3.03, 3.08, 3.10),
                         CV3 = c(3.71, 3.78, 3.79),
                         CV4 = c(4.60, 4.70, 4.68))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
  if(code=="BM100Foddeven")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(2.32, 2.36, 2.40),
                         CV2 = c(3.01, 3.06, 3.10),
                         CV3 = c(3.68, 3.76, 3.79),
                         CV4 = c(4.60, 4.66, 4.68))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
  if(code=="BM101Foddeven")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(5.27, 5.42, 5.64),
                         CV2 = c(6.26, 6.42, 6.67),
                         CV3 = c(7.19, 7.38, 7.63),
                         CV4 = c(8.35, 8.60, 8.79))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
  if(code=="BM110Foddeven")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(2.30, 2.36, 2.40),
                         CV2 = c(2.97, 3.05, 3.10),
                         CV3 = c(3.64, 3.72, 3.79),
                         CV4 = c(4.53, 4.62, 4.68))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
  if(code=="BM111Foddeven")
  {
    tablaVC <- data.frame(N = c(240, 480, 10000),
                         CV1= c(5.25, 5.44, 5.64),
                         CV2 = c(6.23, 6.43, 6.67),
                         CV3 = c(7.14, 7.35, 7.63),
                         CV4 = c(8.33, 8.52, 8.79))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
#
  if(code=="HEGY000t0")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-2.72, -2.60, -2.62, -2.62),
                         CV2 = c(-2.29, -2.26, -2.25, -2.23),
                         CV3 = c(-1.95, -1.97, -1.93, -1.94),
                         CV4 = c(-1.59, -1.61, -1.59, -1.62))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="HEGY100t0")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-3.66, -3.47, -3.51, -3.48),
                         CV2 = c(-3.25, -3.14, -3.17, -3.13),
                         CV3 = c(-2.96, -2.88, -2.89, -2.87),
                         CV4 = c(-2.62, -2.58, -2.58, -2.57))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="HEGY101t0")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-3.77, -3.55, -3.56, -3.51),
                         CV2 = c(-3.39, -3.22, -3.23, -3.18),
                         CV3 = c(-3.08, -2.95, -2.94, -2.91),
                         CV4 = c(-2.72, -2.63, -2.62, -2.59))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="HEGY110t0")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-4.23, -4.07, -4.09, -4.05),
                         CV2 = c(-3.85, -3.73, -3.75, -3.70),
                         CV3 = c(-3.56, -3.47, -3.46, -3.44),
                         CV4 = c(-3.21, -3.16, -3.16, -3.15))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="HEGY111t0")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-4.46, -4.09, -4.15, -4.05),
                         CV2 = c(-4.04, -3.80, -3.80, -3.74),
                         CV3 = c(-3.71, -3.53, -3.52, -3.49),
                         CV4 = c(-3.37, -3.22, -3.21, -3.18))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
#
  if(code=="HEGY000tpi")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-2.67, -2.61, -2.60, -2.60),
                         CV2 = c(-2.27, -222, -2.23, -2.24),
                         CV3 = c(-1.95, -1.92, -1.94, -1.95),
                         CV4 = c(-1.60, -1.57, -1.61, -1.61))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="HEGY100tpi")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-2.68, -2.61, -2.60, -2.58),
                         CV2 = c(-2.27, -2.24, -2.21, -2.22),
                         CV3 = c(-1.95, -1.95, -1.91, -1.92),
                         CV4 = c(-1.60, -1.60, -1.58, -1.59))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="HEGY101tpi")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-3.75, -3.60, -3.49, -3.50),
                         CV2 = c(-3.37, -3.22, -3.15, -3.16),
                         CV3 = c(-3.04, -2.94, -2.90, -2.89),
                         CV4 = c(-2.69, -2.63, -2.59, -2.60))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="HEGY110tpi")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-2.65, -2.58, -2.65, -2.59),
                         CV2 = c(-2.24, -2.24, -2.25, -2.25),
                         CV3 = c(-1.91, -1.94, -1.96, -1.95),
                         CV4 = c(-1.57, -1.60, -1.63, -1.62))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
  if(code=="HEGY111tpi")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(-3.80, -3.60, -3.57, -3.52),
                         CV2 = c(-3.41, -3.22, -3.18, -3.18),
                         CV3 = c(-3.08, -2.94, -2.93, -2.91),
                         CV4 = c(-2.73, -2.63, -2.61, -2.60))
    tablaVCaux <- t(data.frame(signf = c("0.01", "0.025", "0.05", "0.10"),
                               label = c("***", "**", "*", ".")))
    tail <- "left"
  }
#
  if(code=="HEGY000Foddeven")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(2.45, 2.39, 2.41, 2.42),
                         CV2 = c(3.26, 3.12, 3.14, 3.16),
                         CV3 = c(4.04, 3.89, 3.86, 3.92),
                         CV4 = c(5.02, 4.89, 4.81, 4.81))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
  if(code=="HEGY100Foddeven")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(2.32, 2.35, 2.36, 2.37),
                         CV2 = c(3.04, 3.08, 3.00, 3.12),
                         CV3 = c(3.78, 3.81, 3.70, 3.86),
                         CV4 = c(4.78, 4.77, 4.73, 4.76))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
  if(code=="HEGY101Foddeven")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(5.50, 5.56, 5.56, 5.56),
                         CV2 = c(6.60, 6.57, 6.63, 6.61),
                         CV3 = c(7.68, 7.72, 7.66, 7.53),
                         CV4 = c(9.22, 8.74, 8.92, 8.93))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="HEGY110Foddeven")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(2.23, 2.31, 2.33, 2.34),
                         CV2 = c(2.95, 2.98, 3.04, 3.07),
                         CV3 = c(3.70, 3.71, 3.69, 3.76),
                         CV4 = c(4.64, 4.70, 4.57, 4.66))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
  if(code=="HEGY111Foddeven")
  {
    tablaVC <- data.frame(N = c(48, 100, 136, 200),
                         CV1= c(5.37, 5.52, 5.55, 5.56),
                         CV2 = c(6.55, 6.60, 6.62, 6.57),
                         CV3 = c(7.70, 7.52, 7.59, 7.56),
                         CV4 = c(9.27, 8.79, 8.77, 8.96))
    tablaVCaux <- t(data.frame(signf = c("0.99", "0.975", "0.95", "0.90"),
                               label = c("***", "**", "*", ".")))
    tail <- "right"
  }
#
  if(code=="CHp1")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(0.353), CV2 = c(0.470), CV3 = c(0.593), CV4 = c(0.748))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp2")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(0.610), CV2 = c(0.749), CV3 = c(0.898), CV4 = c(1.070))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp3")
  {
    tablaVC <- data.frame(N = c(200),
                          CV1= c(0.846), CV2 = c(1.010), CV3 = c(1.160), CV4 = c(1.350))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp4")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(1.07), CV2 = c(1.24), CV3 = c(1.39), CV4 = c(1.60))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp5")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(1.28), CV2 = c(1.47), CV3 = c(1.63), CV4 = c(1.88))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp6")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(1.49), CV2 = c(1.68), CV3 = c(1.89), CV4 = c(2.12))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp7")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(1.69), CV2 = c(1.90), CV3 = c(2.10), CV4 = c(2.35))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp8")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(1.89), CV2 = c(2.11), CV3 = c(2.33), CV4 = c(2.59))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp9")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(2.10), CV2 = c(2.32), CV3 = c(2.55), CV4 = c(2.82))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp10")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(2.29), CV2 = c(2.54), CV3 = c(2.76), CV4 = c(3.05))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp11")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(2.49), CV2 = c(2.75), CV3 = c(2.99), CV4 = c(3.27))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }
  if(code=="CHp12")
  {
    tablaVC <- data.frame(N = c(200),
                         CV1= c(2.69), CV2 = c(2.96), CV3 = c(3.18), CV4 = c(2.51))
    tablaVCaux <- t(data.frame(signf = c("0.90", "0.95", "0.975", "0.99"),
                               label = c(".", "*", "**", "***")))
    tail <- "right"
  }

#
  list(tablaVC, tablaVCaux, tail)
}

# Mostrar el contenido de un archivo en un ventana redimensionable
mytkpager <- function (file, header, title, delete.file, wwidth, wheight, export)
{
  title <- paste(title, header)
  for (i in seq(along = file))
  {
     zfile <- file[[i]]
     tt <- tktoplevel()
     tkwm.title(tt, if(length(title))
                       title[(i - 1)%%length(title) + 1]
                    else "")
     txt <- tktext(tt, bg = "white", font = "courier", fg="blue")

     scr <- tkscrollbar(tt, repeatinterval = 5, command = function(...) tkyview(txt,...))
     tkconfigure(txt, yscrollcommand = function(...) tkset(scr,...), width=wwidth, height=wheight)
     tkpack(txt, side = "left", fill = "both", expand = TRUE)
     tkpack(scr, side = "right", fill = "y")

     chn <- tkcmd("open", zfile)
     tkinsert(txt, "end", gsub("_\b", "", tclvalue(tkcmd("read", chn))))
     tkcmd("close", chn)
     tkconfigure(txt, state = "disabled")
     tkmark.set(txt, "insert", "0.0")
     tkfocus(txt)

     if(delete.file)
        tkcmd("file", "delete", zfile)
  }
  topMenu <- tkmenu(tt)
  tkconfigure(tt, menu=topMenu)

  if(export==TRUE)
  {
    ToolsMenu <- tkmenu(topMenu, tearoff=FALSE)
    tkadd(ToolsMenu, "command", label="Export to a LaTeX file", command=function()RecTestLaTeX())
    tkadd(topMenu, "cascade", label="Tools", menu=ToolsMenu)
  }
}

#

RecTestLaTeX <- function()  # poner .tex al elegir archivo
{
  N  <- .wvari$N
  s  <- .wvari$s
  t0 <- .wvari$t0

  Trectest <- function(outfile)
  {
    tN <- as.matrix(ysooys(N, t0, N, s))[[1]]
    testname <- rectest.out[[1]]
    freqname <- rectest.out[[2]]
    tabla <- rectest.out[[3]]
    tabla[,2] <- Tround(tabla, 2, 2)[,2]
    tabla[,5] <- Tround(tabla, 5, 2)[,5]

    cat("\\documentclass[11pt]{article}\n", file=outfile, append=TRUE)
    cat("\\begin{document} \n\n", file=outfile, append=TRUE)

    cat("\\begin{table}[h]\n", file=outfile, append=TRUE)
    cat("\\centering\n", file=outfile, append=TRUE)
    cat(paste("\\caption{Recursive testing. ", testname, " test. ", "Frequency ", freqname, "}\n",
              sep=""), file=outfile, append=TRUE)
    cat(paste("\\label{Trectest.", freqname, "}\n", sep=""), file=outfile, append=TRUE)
    cat("\\begin{tabular}{crlcrl} \\\\\n", file=outfile, append=TRUE)

    cat("\\hline \\hline\n", file=outfile, append=TRUE)
    cat(paste("Sample & Stat. & & Sample & Stat. & \\\\\n", sep=""), file=outfile, append=TRUE)

    cat("\\hline\n", file=outfile, append=TRUE)
    for(i in 1:nrow(tabla))
       cat(paste(paste(t0[1], t0[2], sep="."), "-", tabla[i,1], " & ",
                 "$", tabla[i,2], "$", " & ", tabla[i,3], " & ",
                 tabla[i,4], "-", paste(tN[1],tN[2], sep="."), " & ",
                 "$", tabla[i,5], "$", " & ", tabla[i,6], " \\\\\n", sep=""),
                 file=outfile, append=TRUE)
    cat("\\hline \\hline\n", file=outfile, append=TRUE)

    cat("\\end{tabular}", file=outfile, append=TRUE)
    cat("\\end{table}", file=outfile, append=TRUE)

    cat("\n\\end{document}", file=outfile, append=TRUE)
    rm(rectest.out)
  }
  #
  outfile <- tclvalue(tkgetSaveFile(filetypes='{"Text files" {".tex"}} {"All Files" {"*"}}'))

  if(!nchar(outfile))
     tkmessageBox(message="No file was chosen.", icon="error")
  else{
     Trectest(outfile)
     tkmessageBox(message=paste("Table has been saved in", outfile), icon="info")
  }
}

#

# Función para redondear las tabla que se exportan (usando cat...) permitiendo que el último
# decimal sea cero.
  # tabla es una coulumna de matriz de datos
  # colum es una columna de la tabla
  # digits es el número de decimales (se permite que el último sea cero)

Tround <- function(tabla, colum, digits)
{
  catround <- function(x, digits)
  {
    rx <- as.character(round(x, digits=digits))
    logic <- FALSE; i <- 1
    if(substr(rx, 2, 2) != "")
    {
      while(logic == FALSE){
        logic <- substr(rx, i, i) == "."
        i <- i+1
      }
      ifelse(nchar(substr(rx, i, i+digits)) < digits, rx <- paste(rx, "0", sep=""), rx <- rx)
    }
    rx
  }

  for(i in 1:nrow(tabla))
    tabla[i,colum] <- catround(as.numeric(tabla[i,colum]), digits)
  tabla
}

# lookup: Función auxiliar para lookupinfo
# alpha tipo carácter ("0.05")

lookup <- function(code, x, alpha, estad)
{
  tablas     <- lookuptablas(code)
  tablaVC    <- as.matrix(tablas[[1]])
  tablaVCaux <- tablas[[2]]
  tail       <- tablas[[3]]

  sigcol <- as.numeric(which(tablaVCaux[1,] == alpha)) + 1

  if(x <= as.numeric(tablaVC[1,1]))
    f <- tablaVC[1,sigcol]
  if(x > as.numeric(tablaVC[(nrow(tablaVC)),1]))
    f <- as.numeric(tablaVC[(nrow(tablaVC)),sigcol])

  if(x > as.numeric(tablaVC[1,1]) && x <= tablaVC[nrow(tablaVC),1])
  {
    for(i in 1:nrow(tablaVC))
    {
      if(as.numeric(tablaVC[i,1]) < x)
      {
        if(i==1)
          start <- 1
        else{
          if((i+2) > nrow(tablaVC))
            start <- i-1
          else{
            if(x-tablaVC[(i-1),1] > tablaVC[(i+2),1]-x)
               start <- i
            else
               start <- i-1
              }
            }
      }
    }
#
    x0 <- tablaVC[start,1]
    x1 <- tablaVC[start+1,1]
    x2 <- tablaVC[start+2,1]

    f0 <- tablaVC[start,sigcol]
    f1 <- tablaVC[start+1,sigcol]
    f2 <- tablaVC[start+2,sigcol]
#
    b <- (f1-f0)/(x1-x0)
    c <- ((f2-f0)/(x2-x0)-b)/(x2-x1)
    f <- f0 + (b+c*(x-x1))*(x-x0)
  }
  f
}

# lookupinfo: Función para anyadir label (*,**,***,.) al output de lookup

lookupinfo <- function(code, x, alpha, estad)
{
  tablas     <- lookuptablas(code)
  tablaVC    <- as.matrix(tablas[[1]])
  tablaVCaux <- tablas[[2]]
  tail       <- tablas[[3]]

  vcinterp <- lookup(code, x, alpha, estad) # x es el tamanyo muestral, número de datos
                                            # vcinter: valor crítico interpolado

  if(as.numeric(alpha) < 0.5)
    tableinterp <- c(lookup(code, x, "0.01", estad),
                     lookup(code, x, "0.025", estad),
                     lookup(code, x, "0.05", estad),
                     lookup(code, x, "0.10", estad))
  if(as.numeric(alpha) > 0.5)
    tableinterp <- c(lookup(code, x, "0.99", estad),
                     lookup(code, x, "0.975", estad),
                     lookup(code, x, "0.95", estad),
                     lookup(code, x, "0.90", estad))

  # valores críticos interpolados dado un tamanyo muestral

  if(tail=="right")
  {
     #tableaux <- t(data.frame(signf = c("0.05", "0.025", "0.010"), label = c("**", "*", ".")))
     #tableinterp <- sort(tableinterp)
     #if(estad < tableinterp[1])
     #  out <- "*****"
     #if(estad >= tableinterp[length(tableinterp)])
     #  out <- tableaux[2,length(tableinterp)]
     #if(estad >= tableinterp[1] && estad < tableinterp[length(tableinterp)])
     #  out <- tableaux[2,min(which(tableinterp > estad))-1]

     #if(estad < tableinterp[length(tableinterp)])
     #  out <- " "
     #if(estad >= tableinterp[1])
     #  out <- tablaVCaux[2,1]
     #if(estad >= tableinterp[length(tableinterp)] && estad < tableinterp[1])
     #  out <- tablaVCaux[2,max(which(tableinterp > estad))+1]
     if(estad < tableinterp[1])
       out <- " "
     if(estad >= tableinterp[length(tableinterp)])
       out <- tablaVCaux[2,1]
     if(estad < tableinterp[length(tableinterp)] && estad > tableinterp[1]){
       tablaVCaux <- sort(tablas[[2]][2,])
       out <- tablaVCaux[min(which(tableinterp > estad))-1]
     }
  }
  if(tail=="left")  # table[,sigcol], sigcol es la columna es la de los valores críticos
  {
     if(estad > tableinterp[length(tableinterp)])
       out <- " "
     if(estad <= tableinterp[1])
       out <- tablaVCaux[2,1]
     if(estad <= tableinterp[length(tableinterp)] && estad > tableinterp[1])
       out <- tablaVCaux[2,max(which(tableinterp < estad))+1]
  }
  if(tail=="rightleft")
  {

  }
  out
}

# MakeHEGY.rectest: Cuadro de diálogo describir los elementos de la regresión de HEGY

MakeHEGY.rectest <- function()
{
  # tclRequire("BWidget")
  tthegybm <- tktoplevel()
  tkwm.title(tthegybm, "HEGY test")

  # .wvari definido al ejecutar Recursivetesting()
  #string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  #.wvari <<- ExeString(c(string, ""))

  done <<- tclVar(0) # variable para detectar cuándo se ha completado el cuadro de diálogo

  #tkgrid(tklabel(tthegybm,text=" Significance level", fg="blue"), sticky="w")
  #alphas <- c("0.10", "0.05", "0.025", "0.01")
  #comboBox <- tkwidget(tthegybm, "ComboBox", editable=FALSE, width="5", values=alphas)
  #tkgrid(comboBox)

  tkgrid(tklabel(tthegybm, text = "  Select the deterministic components:",
         fg = "blue"), sticky = "w")
  Cte <- tkcheckbutton(tthegybm)
  Tdl <- tkcheckbutton(tthegybm)
  Vfet <- tkcheckbutton(tthegybm)
  CteValue <- tclVar("0")
  TdlValue <- tclVar("0")
  VfetValue <- tclVar("0")
  tkconfigure(Cte, variable = CteValue)
  tkgrid(tklabel(tthegybm, text = "Intercept"), Cte)
  tkconfigure(Tdl, variable = TdlValue)
  tkgrid(tklabel(tthegybm, text = "Trend"), Tdl)
  tkconfigure(Vfet, variable = VfetValue)
  tkgrid(tklabel(tthegybm, text = "Seasonal dummys"), Vfet)

  rb1 <- tkradiobutton(tthegybm)
  rb2 <- tkradiobutton(tthegybm)
  rb3 <- tkradiobutton(tthegybm)
  rb4 <- tkradiobutton(tthegybm)
  rb5 <- tkradiobutton(tthegybm)
  rb6 <- tkradiobutton(tthegybm)
  rbValue <- tclVar("BIC-LB")
  yself <- tclVar()
  entry.yself <- tkentry(tthegybm, width = "3", textvariable = yself)
  tkconfigure(rb1, variable = rbValue, value = "AIC-LB")
  tkconfigure(rb2, variable = rbValue, value = "BIC-LB")
  tkconfigure(rb3, variable = rbValue, value = "AIC-Lüt")
  tkconfigure(rb4, variable = rbValue, value = "BIC-Lüt")
  tkconfigure(rb5, variable = rbValue, value = "Signf")
  tkconfigure(rb6, variable = rbValue, value = "Tu mismo")
  tkgrid(tklabel(tthegybm, text = "  Select the method for choosing lags:",
         fg = "blue"), sticky = "w")
  tkgrid(tklabel(tthegybm, text = "AIC-LB"), rb1)
  tkgrid(tklabel(tthegybm, text = "BIC-LB"), rb2)
  tkgrid(tklabel(tthegybm, text = "AIC-top-down"), rb3)
  tkgrid(tklabel(tthegybm, text = "BIC-top-down"), rb4)
  tkgrid(tklabel(tthegybm, text = "Significant lags"), rb5)
  tkgrid(tklabel(tthegybm, text = "By yourself"), rb6, entry.yself)

  tkgrid(tklabel(tthegybm, text = "  Select the frequency to analyze.", fg="blue"), sticky = "w")
  tl <- tklistbox(tthegybm, height = eval(.wvari$s/2+1), selectmode = "single", bg = "white")
  tkgrid(tl)
  if(.wvari$s == 4)
      freq <- c(" pi/2", " pi", " 0")
  if(.wvari$s == 12)
      freq <- c(" pi/6", " pi/3", " pi/2", " 2pi/3", " 5pi/6", " pi", " 0")
  for(i in 1:(.wvari$s/2+1)) tkinsert(tl, "end", freq[i])
  tkselection.set(tl, 1)

  Confirm <- function(){
    freqChoice <<- rep(NA, 2)
    freqChoice[1] <<- freq[as.numeric(tkcurselection(tl))+1]   # VER +1
    freqChoice[2] <<- which(freq == freq[as.numeric(tkcurselection(tl))+1])

    tkmessageBox(title = "Info", message = "Frequency has been selected.", icon = "info")
    cat(c("\n Frequency selected:", freqChoice[1], " \n\n"))
  }
  Confirm.but <- tkbutton(tthegybm, text = "Select", command = Confirm)
  tkgrid(Confirm.but)

  OnOK <- function(){
     #vari.env <- new.env(FALSE, NULL)
     #Fvarinfo(label$label, env.out=vari.env)
     s <- .wvari$s   # get("s", env=vari.env)
     t0 <- .wvari$t0 # get("t0", env=vari.env)

     # alpha <- alphas[as.numeric(tclvalue(tkcmd(comboBox,"getvalue")))+1]    ###
     tkdestroy(tthegybm)  # después de borrar la ventana comboBox

     mVal <- as.character(tclvalue(CteValue))
     if (mVal == "1")
         Ct <- 1
     if (mVal == "0")
         Ct <- 0
     mVal <- as.character(tclvalue(TdlValue))
     if (mVal == "1")
         TD <- 1
     if (mVal == "0")
         TD <- 0
     mVal <- as.character(tclvalue(VfetValue))
     if (mVal == "1")
         Vfe <- 1
     if (mVal == "0")
         Vfe <- 0
     rbVal <- as.character(tclvalue(rbValue))
     if (tclvalue(rbValue) == "AIC-LB")
         selecP <- "aiclb"
     if (tclvalue(rbValue) == "BIC-LB")
         selecP <- "biclb"
     if (tclvalue(rbValue) == "AIC-Lüt")
         selecP <- "aiclut"
     if (tclvalue(rbValue) == "BIC-Lüt")
         selecP <- "biclut"
     if (tclvalue(rbValue) == "Signf")
         selecP <- "signf"
     if (tclvalue(rbValue) == "Tu mismo")
         selecP <- as.numeric(tclvalue(yself))

     if(s==4){
       if(freqChoice[2]==1) freq <- "[[1]][[2]][1]"
       if(freqChoice[2]==2) freq <- "[[1]][[1]][2]"
       if(freqChoice[2]==3) freq <- "[[1]][[1]][1]"
     }
     if(s==12){
       if(freqChoice[2]==1) freq <- "[[1]][[2]][5]"
       if(freqChoice[2]==2) freq <- "[[1]][[2]][3]"
       if(freqChoice[2]==3) freq <- "[[1]][[2]][1]"
       if(freqChoice[2]==4) freq <- "[[1]][[2]][2]"
       if(freqChoice[2]==5) freq <- "[[1]][[2]][4]"
       if(freqChoice[2]==6) freq <- "[[1]][[1]][2]"
       if(freqChoice[2]==7) freq <- "[[1]][[1]][1]"
     }
     ifelse(freqChoice == " pi" || freqChoice == " 0", alpha <- "0.05", alpha <- "0.95")
     argHEGY <<- list(compdet=c(Ct, TD, Vfe), selecP=selecP, Mvfic=0, VFEp=0,
                      showcat=FALSE, freq=freq, freqname=freqChoice[1], alpha=alpha)
     tclvalue(done) <<- 1
  }
  OK.but <- tkbutton(tthegybm, text = "OK", command = OnOK)
  tkgrid(OK.but)
}

# MakeCH.rectest: Cuadro de diálogo para describir los elementos de la regresión de CH

MakeCH.rectest <-  function()
{
  #tclRequire("BWidget")
  ttch <- tktoplevel()
  tkwm.title(ttch, "CH test")

  # .wvari definido al ejecutar RecursiveTesting()
  #string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  #.wvari <<- ExeString(c(string, ""))

  done <<- tclVar(0)

  #tkgrid(tklabel(ttch,text="Significance level", fg="blue"), sticky="w")
  alphas <- c("0.90", "0.95", "0.975", "0.99")
  #comboBox <- tkwidget(ttch, "ComboBox", editable=FALSE, width="5", values=alphas)
  #tkgrid(comboBox)

  trend <- tkcheckbutton(ttch)
  trendValue <- tclVar("0")
  lag1 <- tkcheckbutton(ttch)
  lag1Value <- tclVar("1")
  tkgrid(tklabel(ttch, text = "  Select the elements to include in the auxiliar regression:",
      fg = "blue"), sticky = "w")
  tkconfigure(trend, variable = trendValue)
  tkgrid(tklabel(ttch, text = "Trend"), trend)
  tkconfigure(lag1, variable = lag1Value)
  tkgrid(tklabel(ttch, text = "First order lag"), lag1)

  #rb1 <- tkradiobutton(ttch)
  #rb2 <- tkradiobutton(ttch)
  #rbValue <- tclVar("nsumm")
  #tkconfigure(rb1, variable = rbValue, value = "nsumm")
  #tkconfigure(rb2, variable = rbValue, value = "summ")
  #tkgrid(tklabel(ttch, text = "   ---   ---   --- "))
  #tkgrid(tklabel(ttch, text = "  Show a summary.", fg = "blue"), rb2, sticky = "w")
  tkgrid(tklabel(ttch, text = "  Select frequencies to analyse:", fg = "blue"), sticky = "w")
  if (.wvari$s == 12) {
      pi6 <- tkcheckbutton(ttch)
      pi56 <- tkcheckbutton(ttch)
      pi3 <- tkcheckbutton(ttch)
      pi23 <- tkcheckbutton(ttch)
  }
  pi2 <- tkcheckbutton(ttch)
  pi <- tkcheckbutton(ttch)
  if (.wvari$s == 12) {
      pi6Value <- tclVar("0")
      pi56Value <- tclVar("0")
      pi3Value <- tclVar("0")
      pi23Value <- tclVar("0")
  }
  pi2Value <- tclVar("0")
  piValue <- tclVar("0")
  if (.wvari$s == 12) {
      tkconfigure(pi6, variable = pi6Value)
      tkgrid(tklabel(ttch, text = "pi/6"), pi6)
      tkconfigure(pi56, variable = pi56Value)
      tkgrid(tklabel(ttch, text = "5pi/6"), pi56)
      tkconfigure(pi3, variable = pi3Value)
      tkgrid(tklabel(ttch, text = "pi/3"), pi3)
      tkconfigure(pi23, variable = pi23Value)
      tkgrid(tklabel(ttch, text = "2pi/3"), pi23)
  }
  tkconfigure(pi2, variable = pi2Value)
  tkgrid(tklabel(ttch, text = "pi/2"), pi2)
  tkconfigure(pi, variable = piValue)
  tkgrid(tklabel(ttch, text = "pi"), pi)

  OnOK <- function() {
    #vari.env <- new.env(FALSE, NULL)
    #Fvarinfo(label$label, env.out=vari.env)
    s <- .wvari$s  #get("s", env=vari.env)

    # alpha <- alphas[as.numeric(tclvalue(tkcmd(comboBox,"getvalue")))+1]
    tkdestroy(ttch)  # después de borrar la ventana comoBox

    trendVal <- as.character(tclvalue(trendValue))
    if (trendVal == "1")
        DetTr <- TRUE
    if (trendVal == "0")
        DetTr <- FALSE
    lag1Val <- as.character(tclvalue(lag1Value))
    if (lag1Val == "1")
        f0 <- 1
    if (lag1Val == "0")
        f0 <- 0

    frec <- rep(0, s/2)
    if(s == 12)
    {
        mVal <- as.character(tclvalue(pi6Value))
        if (mVal == "1")
          frec[1] <- 1
        if (mVal == "0")
          frec[1] <- 0
        mVal <- as.character(tclvalue(pi56Value))
        if (mVal == "1")
          frec[5] <- 1
        if (mVal == "0")
          frec[5] <- 0
        mVal <- as.character(tclvalue(pi3Value))
        if (mVal == "1")
          frec[2] <- 1
        if (mVal == "0")
          frec[2] <- 0
        mVal <- as.character(tclvalue(pi23Value))
        if (mVal == "1")
          frec[4] <- 1
        if (mVal == "0")
          frec[4] <- 0
    }
    if(s == 12)
       aux <- c(3, 6)
    if (s == 4)
       aux <- c(1, 2)

    mVal <- as.character(tclvalue(pi2Value))
    if (mVal == "1")
        frec[aux[1]] <- 1
    if (mVal == "0")
        frec[aux[1]] <- 0
    mVal <- as.character(tclvalue(piValue))
    if (mVal == "1")
        frec[aux[2]] <- 1
    if (mVal == "0")
        frec[aux[2]] <- 0

    argCH <<- list(s=s, frec=frec, f0=f0, DetTr=DetTr, showcat=FALSE, alpha= "0.95") #alpha)
    tclvalue(done) <<- 1
  }
  OK.but <- tkbutton(ttch, text = "OK", command = OnOK)
  tkgrid(OK.but)
}

# RecursiveTesting: Función para aplicar los contrastes de HEGY y CH recursivamente
   # empleando varios periodos muestrales (del principio y final de la muestra)

RecursiveTesting <- function(testname)
{
  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
  .wvari <<- ExeString(c(string, ""))
  vari <- .wvari$vari
  s    <- .wvari$s
  t0   <- .wvari$t0
  N    <- length(vari)

  Nsb <- (N-7*s)/s
  tablaPMs <- matrix(nrow=as.integer(Nsb)+2, ncol=6)
  ifelse(t0[2] == 1, tNs2aux <- s, tNs2aux <- t0[2]-1)

# Cuadro de diálogo para describir los argumentos de las funciones de HEGY y CH

   if(testname == "HEGY"){
      MakeHEGY.rectest()
      tempf <- tempfile()}

   if(testname == "CH")
      MakeCH.rectest()

   tkwait.variable(done)
   # freqname <- argHEGY[[7]]
   #if(testname == "HEGY"){
   #   ifelse(substr(argHEGY[[6]], 12, 12) == "2",
   #          alpha <- as.character(1-as.numeric(argHEGY[[8]])),
   #          alpha <- argHEGY[[8]])
   #   freqname <- argHEGY[[7]]
   #}
   #if(testname == "CH")
   #   alpha <- as.character(1-as.numeric(alpha))

# Construir el argumento "code" para las funciones lookup y lookupinfo.

  if(testname == "HEGY")
  {
    freqname <- argHEGY[[7]]
    if(s==4) { codeaux1 <- "HEGY"}
    if(s==12){ codeaux1 <- "BM"  }

    codeaux2 <- paste(argHEGY[[1]][1], argHEGY[[1]][2], argHEGY[[1]][3], sep="")

    if(substr(argHEGY[[6]], 8, 8) == "1"){
      if(substr(argHEGY[[6]], 12, 12) == "1")
        codeaux3 <- "t0"
      if(substr(argHEGY[[6]], 12, 12) == "2")
        codeaux3 <- "tpi"
    }
    if(substr(argHEGY[[6]], 8, 8) == "2")
      codeaux3 <- "Foddeven"

    code <- paste(codeaux1, codeaux2, codeaux3, sep="")
  }
  if(testname == "CH")
  {
    if(s==4){  fnames <- c("pi2", "pi") }
    #if(s==12){ fnames <- c("pi6", "5pi6", "pi3", "2pi3", "pi2", "pi") }
    if(s==12){ fnames <- c("pi6", "pi3", "pi2", "2pi3", "5pi6", "pi") }
    freqnames <- fnames[which(argCH[[2]] == 1)]
    codeaux1  <- "CH"

    sq     <- seq(1,11,2)
    frecob <- rep(0,(s-1))
    for(i in 1:length(argCH[[2]])){
     if(argCH[[2]][i] == 1 && i == s/2)
       frecob[sq[i]]     <- 1
     if(argCH[[2]][i] == 1 && i < s/2)
       frecob[sq[i]] <- frecob[sq[i]+1] <- 1
    }
    codeaux2 <- paste("p", length(which(frecob == 1)), sep="")

    code <- paste(codeaux1, codeaux2, sep="")
  }

# Con la observación inicial fija, el primer dato:

  for(i in 0:as.integer(Nsb))
  {
    tNsubm <- ysooys(c((t0[1]+7+i), t0[2]), t0, N, s)[[1]] -1
    tablaPMs[(i+1),1] <- paste(ysooys(tNsubm, t0, N, s)[[1]][1],
                               ysooys(tNsubm, t0, N, s)[[1]][2], sep=".")

    # Nsvd <- N
    variaux <- ts(vari[1:tNsubm], frequency=.wvari$s, start=.wvari$t0)
    .wvariaux <<- list(vari=variaux, s=.wvari$s, t0=.wvari$t0, N=length(variaux),
                       logvari=.wvari$logvari, label=".wvariaux")

    if(substr(code, 1,2)=="HE" || substr(code, 1,2)=="BM")
    {
      # t0 <<-c(0,1)  ### # #  #  VER
      rdotest <<- HEGY.test(.wvariaux, argHEGY[[1]], argHEGY[[2]], argHEGY[[3]],
                            argHEGY[[4]], argHEGY[[5]])
      cat(paste("rdotest", argHEGY$freq, sep=""), file=tempf, append=FALSE)
      tablaPMs[(i+1),2] <- eval(source(tempf))[[1]]

      # Última fila (todo el periodo muestral)
      if(i==0){
        # N <<- Nsvd;
        rdotest <<- HEGY.test(.wvari, argHEGY[[1]], argHEGY[[2]], argHEGY[[3]],
                              argHEGY[[4]], argHEGY[[5]])
        cat(paste("rdotest", argHEGY$freq, sep=""), file=tempf, append=FALSE)
        tablaPMs[nrow(tablaPMs),2] <- eval(source(tempf))[[1]]
        # N <<- length(variaux)
      }
      tablaPMs[(i+1),3] <- lookupinfo(code, .wvariaux$N, argHEGY[[8]],
                                      as.numeric(tablaPMs[(i+1),2]))
    }
    if(substr(code, 1,2)=="CH")
    {
      tablaPMs[(i+1),2] <- CH.test(.wvariaux, argCH[[2]], argCH[[3]],
                                   argCH[[4]], argCH[[5]])[[1]]
      if(i==0){  # Última fila (todo el periodo muestral)
      #N <<- Nsvd
        tablaPMs[nrow(tablaPMs),2] <- CH.test(.wvari, argCH[[2]], argCH[[3]],
                                              argCH[[4]], argCH[[5]])[[1]]
      #N <<- length(variaux)
      }
      #tablaPMs[(i+1),3] <- lookupinfo(code, .wvariaux$N, "0.95", as.numeric(tablaPMs[(i+1),2]))
      Tvc <- lookuptablas(code)
      vc.aux <- which(Tvc[[1]][2:5] < as.numeric(tablaPMs[(i+1),2]))
      ifelse(length(vc.aux) == 0, tablaPMs[(i+1),3] <- " ",
                                  tablaPMs[(i+1),3] <- Tvc[[2]][2,vc.aux[length(vc.aux)]])
    }
    #if(testname == "HEGY")
    #  tablaPMs[(i+1),3] <- lookupinfo(code, N, argHEGY[[8]], as.numeric(tablaPMs[(i+1),2]))
    #if(testname == "CH")
    #  tablaPMs[(i+1),3] <- lookupinfo(code, N, "0.95", as.numeric(tablaPMs[(i+1),2]))
    #N <<- Nsvd
  }
  # Última fila (todo el periodo muestral)
  tN <- ysooys(N, t0, N, s)[[1]]
  tablaPMs[nrow(tablaPMs),1] <- paste(tN[1], tN[2], sep=".")
  if(testname == "HEGY"){
    aux <- argHEGY[[8]]
    tablaPMs[nrow(tablaPMs),3] <- lookupinfo(code, .wvari$N, aux,
                                             as.numeric(tablaPMs[nrow(tablaPMs),2]))
                        }
  if(testname == "CH"){
    Tvc <- lookuptablas(code)
    vc.aux <- which(Tvc[[1]][2:5] < as.numeric(tablaPMs[nrow(tablaPMs),2]))
    ifelse(length(vc.aux) == 0, tablaPMs[nrow(tablaPMs),3] <- " ",
                                tablaPMs[nrow(tablaPMs),3] <- Tvc[[2]][2,vc.aux[length(vc.aux)]])
                      }
# Fijando el último dato con la última observación

  for(i in 0:as.integer(Nsb))
  {
    t0subm <- ysooys(c((tN[1]-7-i), tN[2]), t0, N, s)[[1]] - 1
    t0aux <- ysooys(t0subm, t0, N, s)[[1]]
    tablaPMs[(i+1),4] <- paste(ysooys(t0subm, t0, N, s)[[1]][1],
                               ysooys(t0subm, t0, N, s)[[1]][2], sep=".")
    # Nsvd <- N
    # 1 cambiado por t0subm, tNsubm cambiado por tN
    #tN <- ysooys(N, t0, N, s)[[1]]  # ver que N es length(vari), no de una submuestra  ####
    variaux <- ts(vari[t0subm:.wvari$N], frequency=.wvari$s, start=t0aux)
    .wvariaux <<- list(vari=variaux, s=.wvari$s, t0=.wvari$t0,
                       N=length(variaux), logvari=.wvari$logvari, label=".wvariaux")
    # N <<- length(variaux)

    if(substr(code, 1,2)=="HE" || substr(code, 1,2)=="BM")
    {
      rdotest <<- HEGY.test(.wvariaux, argHEGY[[1]], argHEGY[[2]], argHEGY[[3]],
                           argHEGY[[4]], argHEGY[[5]])
      cat(paste("rdotest", argHEGY$freq, sep=""), file=tempf, append=FALSE)
      tablaPMs[(i+1),5] <- eval(source(tempf))[[1]]

      # Última fila (todo el periodo muestral)
      if(i==0){
        # N <<- Nsvd;
        rdotest <<- HEGY.test(.wvari, argHEGY[[1]], argHEGY[[2]], argHEGY[[3]],
                              argHEGY[[4]], argHEGY[[5]])
        cat(paste("rdotest", argHEGY$freq, sep=""), file=tempf, append=FALSE)
        tablaPMs[nrow(tablaPMs),5] <- eval(source(tempf))[[1]]
        # N <<- length(variaux)
      }
      tablaPMs[(i+1),6] <- lookupinfo(code, .wvariaux$N, argHEGY[[8]],
                                      as.numeric(tablaPMs[(i+1),5]))
    }
    if(substr(code, 1,2)=="CH")
    {
      tablaPMs[(i+1),5] <- CH.test(.wvariaux, argCH[[2]], argCH[[3]],
                                   argCH[[4]], argCH[[5]])[[1]]
      if(i==0){ # Última fila (todo el periodo muestral)
      #N <<- Nsvd;
        tablaPMs[nrow(tablaPMs),5] <- CH.test(.wvari, argCH[[2]], argCH[[3]],
                                              argCH[[4]], argCH[[5]])[[1]]
      #N <<- length(variaux)
      }
      #tablaPMs[(i+1),6] <- lookupinfo(code, .wvariaux$N, "0.95", as.numeric(tablaPMs[(i+1),5]))
      Tvc <- lookuptablas(code)
      vc.aux <- which(Tvc[[1]][2:5] < as.numeric(tablaPMs[(i+1),5]))
      ifelse(length(vc.aux) == 0, tablaPMs[(i+1),6] <- " ",
                                  tablaPMs[(i+1),6] <- Tvc[[2]][2,vc.aux[length(vc.aux)]])
    }
    #
    # N <<- Nsvd
  }
  # Última fila (todo el periodo muestral)
  tablaPMs[nrow(tablaPMs),4] <- paste(t0[1], t0[2], sep=".")
  if(testname == "HEGY"){
    aux <- argHEGY[[8]]
    tablaPMs[nrow(tablaPMs),6] <- lookupinfo(code, .wvari$N, aux,
                                           as.numeric(tablaPMs[nrow(tablaPMs),5]))
                         }
  if(testname == "CH"){
    Tvc <- lookuptablas(code)
    vc.aux <- which(Tvc[[1]][2:5] < as.numeric(tablaPMs[nrow(tablaPMs),5]))
    ifelse(length(vc.aux) == 0, tablaPMs[nrow(tablaPMs),6] <- " ",
                                tablaPMs[nrow(tablaPMs),6] <- Tvc[[2]][2,vc.aux[length(vc.aux)]])
                      }
  # rm(.wvariaux)
  # clean.auxobjects()

# showcat
  tablaPMs[,2] <- Tround(tablaPMs, 2, 2)[,2]
  tablaPMs[,5] <- Tround(tablaPMs, 5, 2)[,5]

  tempf <- tempfile()
  cat(c("\n ", testname, " recursive testing for unit roots. \n"), file=tempf, append=TRUE)

  if(testname == "HEGY"){ freqlabel <- argHEGY[[7]] }
  if(testname == "CH"){ freqlabel <- freqnames }
  if(testname == "HEGY") cat(c("\n Frequency", argHEGY[[7]],"\n"), file=tempf, append=TRUE)
  if(testname == "CH")   cat(c("\n Frequency", freqnames,"\n"), file=tempf, append=TRUE)

  cat(" ----- ----- ----- -----\n\n", file=tempf, append=TRUE)
  cat("  Forward recursive samples: \n\n", file=tempf, append=TRUE)
  for(i in 0:(nrow(tablaPMs)-1))
  {
    cat(c(" from", paste(t0[1], t0[2], sep="."), "to", tablaPMs[(i+1),1], ": ", tablaPMs[(i+1),2],
          tablaPMs[(i+1),3], "\n"), file=tempf, append=TRUE)
  }
  cat("\n\n  Backward recursive samples: \n\n", file=tempf, append=TRUE)
  tN <- ysooys(N, t0, N, s)[[1]]        # ya está definida arriba, pero no la detecta aquí
  for(i in (nrow(tablaPMs)-1):0)
  {
    cat(c(" from", tablaPMs[(i+1),4], "to", paste(tN[1], tN[2], sep="."),
        ": ", tablaPMs[(i+1),5], tablaPMs[(i+1),6],"\n"), file=tempf, append=TRUE)
  }
  mytkpager(file=tempf, title=testname, header="recursive testing",
            wheight=25, wwidth=60, delete.file=TRUE, export=TRUE)

  # tablaPMs <<- tablaPMs
  rectest.out <<- list(testname, freqlabel, tablaPMs)

  #tkwait.variable(tablaPMs[nrow(tablaPMs),6])
  file.remove(tempf)

  #cleanlist <- list(); cleanlist[[1]] <- list(done=done)
  if(exists("argHEGY")){
    string <- "rm(argHEGY, rdotest, freqChoice, done)"
    clean.auxobjects(string, type="ExString")
  }
  if(exists("argCH")){
    string <- "rm(argCH, done)"
    clean.auxobjects(string, type="ExString")
  }
  #clean.auxobjects(cleanlist)
}

#change.wvari.label <- function(arg)
#{
#   tempf <- tempfile()
#   for(i in 1:length(arg))
#   {
#     cat(paste(".wvari$", arg[i], " <- ", arg[i], "; rm(", arg[i], ")", sep=""), file=tempf)
#     eval(source(tempf))[[1]]
#   }
#
#   tkcmd("file", "delete", tempf)
#   assign(tclvalue(.label), .wvari, env=.GlobalEnv)
#}

# MAKE urca PACKAGE

#Make.ca.jo <- function(type)
#{
#  require(urca)
#  string <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
#
#  names <- rep(NA, 7)      # name of the selected series
#  sep   <- rep(NA, 7)      # at most 7 series can be selected
#  j <- 1; for(i in 1:nchar(string))
#     if(substr(string, i, i) == " "){ sep[j] <- i; j <- j+1 }
#  sep <- na.omit(sep)
#
#  names[1] <- substr(string, 1, sep[1]-1)
#  for(i in 1:(length(sep)-1))
#    names[i+1] <- substr(string, sep[i]+1, sep[i+1]-1)
#  names[length(sep)+1] <- substr(string, sep[length(sep)]+1, nchar(string))
#  names <- na.omit(names)
#
#  Nvec <- length(names)             # number of series
#  Vvari <- rep(list(), Nvec)        # vector of variables (series)
#  for(i in 1:Nvec)
#    Vvari[[i]] <- ExeString(c(names[i], "$vari"))
#  aux <- data.frame(Vvari[[1]], Vvari[[2]])
#  if(Nvec > 2){
#    for(i in 3:Nvec)
#      aux <- data.frame(aux, Vvari[[i]])
#  }
#  Vvari <- aux
#
#  if(type=="eigen") ca.jo(Vvari, constant=TRUE, type="eigen", K=2, spec="longrun", season=4)
#  if(type=="trace") ca.jo(Vvari, constant=TRUE, type="trace", K=2, spec="longrun", season=4)
#}

# GUI SCHEME

urootgui <- function()
{
  require(stats) || stop("The package stats is not available")
  require(tcltk, keep.source=FALSE) || stop("The package tcltk is not available")
  tclRequire("BWidget")
#
  .tt <<- tktoplevel()
  tkwm.title(.tt, "uroot R-GUI 1.0")

  xScr        <- tkscrollbar(.tt, command=function(...)tkxview(treeWidget,...),
                            orient="horizontal")
  yScr        <- tkscrollbar(.tt, command=function(...)tkyview(treeWidget,...))
  .treeWidget <<- tkwidget(.tt, "Tree", xscrollcommand=function(...)tkset(xScr,...),
           yscrollcommand=function(...)tkset(yScr,...), width=45, height=17, bg="aquamarine3")
  tkgrid(.treeWidget, yScr)
  tkgrid.configure(yScr, stick="nsw")
  tkgrid(xScr)
  tkgrid.configure(xScr, stick="new")
#

# Mouse right button
  editPopupMenu <- tkmenu(.tt, tearoff=FALSE)
  # tkadd(editPopupMenu, "command", label="Export graphic", command=ExportGraph)
  tkadd(editPopupMenu, "command", label="Info", command=function() variinfo_treeW())
  tkadd( editPopupMenu, "command", label="Delete", command=function()
         tkdelete(.treeWidget, tclvalue(tkcmd(.treeWidget,"selection","get"))) )
  RightClick <- function(x,y) # x e y son las coordenadas del ratón
  {
    rootx <- as.integer(tkwinfo("rootx", .tt))
    rooty <- as.integer(tkwinfo("rooty", .tt))
    xTxt <- as.integer(x)+rootx
    yTxt <- as.integer(y)+rooty
    .Tcl(paste("tk_popup",.Tcl.args(editPopupMenu,xTxt,yTxt)))
  }
  tkbind(.tt, "<Button-3>", RightClick)

#
topMenu <- tkmenu(.tt)
tkconfigure(.tt, menu=topMenu)

# ARCHIVO

FileMenu <- tkmenu(topMenu, tearoff=FALSE)
OpenMenu <- tkmenu(topMenu, tearoff=FALSE)

#tkadd(OpenMenu, "command", label="workspace", command=function() LoadWorkSpace())
tkadd(OpenMenu, "command", label="source code", command=function() OpenSourceFile())
tkadd(FileMenu, "cascade", label="Open", menu=OpenMenu)

#tkadd(FileMenu, "command", label="Save", command=function() save.image())
#tkadd(FileMenu, "command", label="Save as", command=function() SaveAsWorkSpace())

tkadd(FileMenu, "command", label="Quit", command=function() closerusea())

tkadd(topMenu, "cascade", label="File", menu=FileMenu)

# DATOS

DataMenu      <- tkmenu(topMenu, tearoff=FALSE)
#DataInfoMenu <- tkmenu(topMenu, tearoff=FALSE)
ImportMenu    <- tkmenu(topMenu, tearoff=FALSE)
BDatosMenu    <- tkmenu(topMenu, tearoff=FALSE)
TransfMenu    <- tkmenu(topMenu, tearoff=FALSE)
PMtrMenu      <- tkmenu(topMenu, tearoff=FALSE)

tkadd(ImportMenu, "command", label="from text file",
      command=function() ReadDataTXT())
tkadd(ImportMenu, "command", label="from CSV file",
      command=function() ReadDataCSV())
# tkadd(ImportMenu, "command", label="from SPSS",
#      command=function() ReadSPSS())
tkadd(DataMenu, "cascade", label="Import data", menu=ImportMenu)
tkadd(DataMenu, "command", label="Description", command=function() GetDataInfo())

tkadd(BDatosMenu, "command", label="CAPV data bank", command=function() dataIKERBIDE())
tkadd(BDatosMenu, "command", label="INE data bank", command=function() dataINE())
tkadd(BDatosMenu, "command", label="Franses (1996)", command=function() dataFranses())
tkadd(DataMenu, "cascade", label="Data bank", menu=BDatosMenu)

#tkadd(DataMenu, "command", label="Simulate DGP", command=function() MakeDGPsim())

tkadd(TransfMenu, "command", label="Logarithmic", command=function(){
      string    <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
      .wvari   <<- ExeString(c(string, ""))
      newstring <- paste("log_", string, sep="")

      vari <- log(.wvari$vari, base=exp(1))
      assign(newstring, list(vari=vari, s=.wvari$s, t0=.wvari$t0,
             N=length(vari), logvari=TRUE, label=newstring), env=.GlobalEnv)
      tkinsert(.treeWidget,"end", string, newstring, text=newstring)
})

tkadd(TransfMenu, "command", label="Box-Cox transformation", command=function() Makeboxcox())

tkadd(TransfMenu, "command", label="--- --- ---", command=function(){})

tkadd(TransfMenu, "command", label="First differences", command=function(){
        # ifelse(.wvari$logvari==TRUE, label <<- "First differences of the logarithms",
        #                              label <<- "First differences")
        string    <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
        .wvari   <<- ExeString(c(string, ""))
        newstring <- paste("fdiff_", string, sep="")

        vari  <- diff(.wvari$vari, lag=1)
        t0aux <- ysooys(.wvari$t0, .wvari$t0, .wvari$N, .wvari$s)
        t0    <- t0aux[[2]][(t0aux[[1]]+1),1:2]
        assign(newstring, list(vari=vari, s=.wvari$s, t0=t0,
                N=length(vari), logvari=.wvari$logvari, label=newstring), env=.GlobalEnv)
        tkinsert(.treeWidget,"end", string, newstring, text=newstring)
        # change.wvari.label(c("t0t", "t0", "varit", "vari", "N"))
})

tkadd(TransfMenu, "command", label="Seasonal differences", command=function(){
        # ifelse(.wvari$logvari==TRUE, label <<- "Seasonal differences of the logarithms",
        #                               label <<- "Seasonal differences")
        string    <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
        .wvari   <<- ExeString(c(string, ""))
        newstring <- paste("sdiff_", string, sep="")

        vari  <- diff(.wvari$vari, lag=.wvari$s)
        t0aux <- ysooys(.wvari$t0, .wvari$t0, .wvari$N, .wvari$s)
        t0    <- t0aux[[2]][(t0aux[[1]]+.wvari$s),1:2]
        assign(newstring, list(vari=vari, s=.wvari$s, t0=t0,
                N=length(vari), logvari=.wvari$logvari, label=newstring), env=.GlobalEnv)
        tkinsert(.treeWidget,"end", string, newstring, text=newstring)
        #change.wvari.label(c("label", "t0t", "varit", "vari", "N"))
})

tkadd(TransfMenu, "command", label="First and seasonal differences", command=function(){
        # ifelse(.wvari$logvari==TRUE,
        # label <<- "First and seasonal differences of the logarithms",
        # label <<- "First and seasonal differences")
        string    <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
        .wvari   <<- ExeString(c(string, ""))
        newstring <- paste("fsdiff_", string, sep="")

        vari <- diff(diff(.wvari$vari, lag=.wvari$s), lag=1)
        t0aux <- ysooys(.wvari$t0, .wvari$t0, .wvari$N, .wvari$s)
        t0    <- t0aux[[2]][(t0aux[[1]]+.wvari$s+1),1:2]
        assign(newstring, list(vari=vari, s=.wvari$s, t0=t0, N=length(vari),
                logvari=.wvari$logvari, label=newstring), env=.GlobalEnv)
        tkinsert(.treeWidget,"end", string, newstring, text=newstring)
        # change.wvari.label(c("label", "t0t", "t0", "varit", "vari", "N"))
})

tkadd(TransfMenu, "command", label="Periodic differences", command=function(){
       string    <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
       .wvari   <<- ExeString(c(string, ""))
       newstring <- paste("pdiff_", string, sep="")

       vari  <- perdiff(.wvari)
       t0aux <- ysooys(.wvari$t0, .wvari$t0, .wvari$N, .wvari$s)
       t0    <- t0aux[[2]][(t0aux[[1]]+1),1:2]
       vari  <- ts(vari[2:length(vari)], frequency=.wvari$s, start=t0)
       assign(newstring, list(vari=vari, s=.wvari$s, t0=t0,
               N=length(vari), logvari=.wvari$logvari, label=newstring), env=.GlobalEnv)
       tkinsert(.treeWidget,"end", string, newstring, text=newstring) })

tkadd(TransfMenu, "command", label="--- --- ---", command=function(){})

tkadd(TransfMenu, "command", label="Remove deterministic component", command=function(){
      #MakeTransfdet()
      string    <- tclvalue(tkcmd(.treeWidget, "selection", "get"))
      .wvari   <<- ExeString(c(string, ""))
      newstring <- paste("undet_", string, sep="")

      vari  <- Transfdet(.wvari$vari, .wvari$s)
      vari  <- ts(vari, frequency=.wvari$s, start=.wvari$t0)
      assign(newstring, list(vari=vari, s=.wvari$s, t0=.wvari$t0, N=length(vari),
              logvari=.wvari$logvari, label=newstring), env=.GlobalEnv)
      tkinsert(.treeWidget,"end", string, newstring, text=newstring)
})

tkadd(DataMenu, "cascade", label="Transformation", menu=TransfMenu)

tkadd(PMtrMenu, "command", label="Change",
      command=function() CambiarPeriodo())
#tkadd(PMtrMenu, "command", label="Restore", command=function(){
#        t0 <<- .wvari$t0t
#        vari <<- ts(.wvari$varit, frequency=.wvari$s, start=t0); N <<- length(vari)
#
        #change.wvari.label(c("t0", "vari", "N"))
#        tkmessageBox(message="Original sample period has been restored", icon="info") })

tkadd(DataMenu, "cascade", label="Sample period", menu=PMtrMenu)
tkadd(topMenu, "cascade", label="Data", menu=DataMenu)

# GRÁFICOS

GrafMenu      <- tkmenu(topMenu, tearoff=FALSE)
TransGrafMenu <- tkmenu(topMenu, tearoff=FALSE)
SpecMenu      <- tkmenu(topMenu, tearoff=FALSE)
BBGrafMenu    <- tkmenu(topMenu, tearoff=FALSE)
QuartergMenu  <- tkmenu(topMenu, tearoff=FALSE)
BBmpMenu      <- tkmenu(topMenu, tearoff=FALSE)

tkadd(TransGrafMenu, "command", label="original", command=function() Makeplotvari())
tkadd(TransGrafMenu, "command", label="logarithms", command=function() Makeplotlog())
tkadd(TransGrafMenu, "command", label="Box-Cox transformation", command=function()
      Makeplotboxcox())
tkadd(TransGrafMenu, "command", label="--- --- ---", command=function(){})

tkadd(TransGrafMenu, "command", label="first differences",
      command=function() Makeplotdelta())
tkadd(TransGrafMenu, "command", label="seasonal differences",
      command=function() Makeplotdeltas())
tkadd(TransGrafMenu, "command", label="first and seasonal differences", command=function()
      Makeplotddeltas())
tkadd(TransGrafMenu, "command", label="periodic differences", command=function() Makeplotperdiff())
tkadd(TransGrafMenu, "command", label="--- --- ---", command=function(){})

tkadd(TransGrafMenu, "command", label="without deterministic component",
      command=function(){ variodet <- Transfdet(.wvari$vari, .wvari$s); plot(variodet)})
tkadd(GrafMenu, "cascade", label="Series", menu=TransGrafMenu)

tkadd(GrafMenu, "command", label="Range-mean", command=function() Makermp())
tkadd(GrafMenu, "command", label="Correlograms", command=function() Makecorrgrm())

tkadd(SpecMenu, "command", label="of the original series",
      command=function() Makespec(dif1=FALSE))
tkadd(SpecMenu, "command", label="of the first differences",
      command=function() Makespec(dif1=TRUE))
tkadd(GrafMenu, "cascade", label="Spectrum", menu=SpecMenu)

tkadd(BBGrafMenu, "command", label="Anual path", command=function() Makebbap())

tkadd(QuartergMenu, "command", label="of the original series",
      command=function() Makequarterg("orig"))
tkadd(QuartergMenu, "command", label="of the first differences",
      command=function() Makequarterg("fdiff"))
tkadd(QuartergMenu, "command", label="of the periodic differences",
      command=function() Makequarterg("pdiff"))
tkadd(BBGrafMenu, "cascade", label="Quarterly path", menu=QuartergMenu)

tkadd(BBmpMenu, "command", label="of the original series", command=function() Makebbmp("orig"))
tkadd(BBmpMenu, "command", label="of the first differences", command=function() Makebbmp("fdiff"))
tkadd(BBmpMenu, "command", label="of the periodic differences",
      command=function() Makebbmp("pdiff"))
tkadd(BBGrafMenu, "cascade", label="Monthly path", menu=BBmpMenu)

tkadd(BBGrafMenu, "command", label="Buys-Ballot 3D",
      command=function() Makebb3D())
tkadd(BBGrafMenu, "command", label="Contour",
      command=function() Makebbcn())
tkadd(GrafMenu, "cascade", label="Buys-Ballot",
      menu=BBGrafMenu)

tkadd(GrafMenu, "command", label="Seasonal box plot",
      command=function() MakeSeasboxplot())

tkadd(GrafMenu, "command", label="Filter frequencies", command=function() MakeFiltrarfrec())

tkadd(topMenu, "cascade", label="Graphics", menu=GrafMenu)

# CONTRASTES

TestMenu  <- tkmenu(topMenu, tearoff=FALSE)
CHMenu    <- tkmenu(topMenu, tearoff=FALSE)
CHrecMenu <- tkmenu(topMenu, tearoff=FALSE)
HEGYMenu  <- tkmenu(topMenu, tearoff=FALSE)

tkadd(TestMenu, "command", label="KPSS", command=function() MakeKPSS.test())
tkadd(TestMenu, "command", label="ADF", command=function() MakeADF.test())

#tkadd(CHMenu, "command", label="about the frequencies", command=function() MakeCH.test())
tkadd(CHMenu, "command", label="about the seasons", command=function() MakeCHseas.test())
tkadd(CHrecMenu, "command", label="original sample", command=function() MakeCH.test())
tkadd(CHrecMenu, "command", label="recursive test",
                 command=function() RecursiveTesting("CH"))
tkadd(CHMenu, "cascade", label="about the frequencies", menu=CHrecMenu)
tkadd(TestMenu, "cascade", label="CH", menu=CHMenu)

tkadd(HEGYMenu, "command", label="original sample", command=function() MakeHEGY.test())
tkadd(HEGYMenu, "command", label="recursive test",
                 command=function() RecursiveTesting("HEGY"))

tkadd(TestMenu, "cascade", label="HEGY", menu=HEGYMenu)
tkadd(topMenu, "cascade", label="Tests", menu=TestMenu)

# UTILIDADES

UtilMenu   <- tkmenu(topMenu, tearoff=FALSE)
LatexMenu  <- tkmenu(topMenu, tearoff=FALSE)
PanelMenu  <- tkmenu(topMenu, tearoff=FALSE)
PanelqMenu <- tkmenu(topMenu, tearoff=FALSE)
PanelmMenu <- tkmenu(topMenu, tearoff=FALSE)

#tkadd(UtilMenu, "command", label="Remove all objects", command=function() Makermall())
tkadd(UtilMenu, "command", label="Save current graph", command=function() ExportGraph())

tkadd(LatexMenu, "command", label="Stochastic components", command=function() SalidaLaTeXfrec())
tkadd(LatexMenu, "command", label="Deterministic components", command=function() SalidaLaTeXdet())
tkadd(UtilMenu, "cascade", label="LaTeX output", menu=LatexMenu)

tkadd(PanelMenu, "command", label="Correlograms", command=function() MakePanelqmCorrg())
tkadd(PanelMenu, "command", label="Seasonal frequencies", command=function() Makefreqg())

tkadd(PanelqMenu, "command", label="Selected graphics", command=function() MakeQPanel1())
tkadd(PanelqMenu, "command", label="Buys-Ballot panel", command=function() MakeQPanelBB())

tkadd(PanelmMenu, "command", label="Selected graphics", command=function() MakePanelmSerie())
tkadd(PanelmMenu, "command", label="Buys-Ballot 1", command=function() MakePanelmBB1())
tkadd(PanelmMenu, "command", label="Buys-Ballot 2", command=function() MakePanelmBB2())
tkadd(PanelMenu, "cascade", label="Quarterly Series ", menu=PanelqMenu)
tkadd(PanelMenu, "cascade", label="Monthly series ", menu=PanelmMenu)
tkadd(UtilMenu, "cascade", label="Panel", menu=PanelMenu)

tkadd(topMenu, "cascade", label="Utilities", menu=UtilMenu)

# urca PACKAGE
#urcaMenu  <- tkmenu(topMenu, tearoff=FALSE)
#ca.joMenu <- tkmenu(topMenu, tearoff=FALSE)
#tkadd(ca.joMenu, "command", label="eigen statistic", command=function() Make.ca.jo(type="eigen"))
#tkadd(ca.joMenu, "command", label="trace statistic", command=function() Make.ca.jo(type="trace"))
#tkadd(urcaMenu, "cascade", label="Tests", menu=ca.joMenu)
#tkadd(topMenu, "cascade", label="urca", menu=urcaMenu)

# AYUDA

AyudaMenu <- tkmenu(topMenu, tearoff=FALSE)

tkadd(AyudaMenu, "command", label="About uroot", command=function()
     mytkpager(file.path(R.home(), "library/uroot/DESCRIPTION"), title="About uroot",
               header="", delete.file=FALSE, wwidth=90, wheight=30, export=FALSE))
#browseURL(file.path(R.home(), "library/urootcook/html/about.html"), browser=getOption("browser")))

tkadd(AyudaMenu, "command", label="Html help", command=function()
      browseURL(file.path(R.home(), "library/uroot/html/00Index.html"),
                browser=getOption("browser")))

tkadd(AyudaMenu, "command", label="Maintainer homepage", command=function()
      browseURL("http://www.bl.ehu.es/~jedlobej", browser=getOption("browser")))

tkadd(topMenu, "cascade", label="  Help", menu=AyudaMenu)
}
urootgui()
