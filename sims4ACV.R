mydetectors<-list(grid10x10large = make.grid(10,10, spacex=600, spacey=600),
              grid10x10med=make.grid(10, 10, spacex = 500, spacey=500), 
              grid10x10small=make.grid(10, 10, spacex=250, spacey=250))


scenarios<-make.scenarios (trapsindex = 1:3, noccasions = 10, nrepeats = 1, D=.003, 
                           g0=c(0.1, 0.2, 0.3),
                sigma=c(600, 1000, 1600), detectfn = 0, recapfactor = .3, popindex = 1,
                detindex = 1, fitindex = 1, crosstraps = TRUE)



x<-run.scenarios (500, scenarios, mydetectors, xsigma = 3,
             nx = 60,  fit = FALSE, extractfn =
               NULL, multisession = FALSE, ncores = 1)
scenarios$EN<-NA
scenarios$meancap<-NA
for (i in 1:27){
  scenarios$EN[i]<-0.003*(spacingtraps[scenarios$trapsindex[i]]*10+6*scenarios$sigma[i])
  scenarios$meancap[i]<-mean(x$output[[i]][,1])
}


scenarios2<-make.scenarios (trapsindex = 1:3, noccasions = 10, nrepeats = 1, D=.004, g0=0.1,
                           sigma=600, detectfn = 0, recapfactor = 1, popindex = 1,
                           detindex = 1, fitindex = 1, crosstraps = TRUE)




x2<-run.scenarios (100, scenarios2, mydetectors, , xsigma = 4,
                  nx = 60,  fit = FALSE, extractfn =
                    NULL, multisession = FALSE, ncores = 1, seed = 123)


scenarios3<-make.scenarios (trapsindex = 1:3, noccasions = 10, nrepeats = 1, D=.003, g0=0.1,
                            sigma=600, detectfn = 0, recapfactor = 0.3, popindex = 1,
                            detindex = 1, fitindex = 1, crosstraps = TRUE)


x3<-run.scenarios (100, scenarios3, mydetectors, , xsigma = 4,
                   nx = 60,  fit = FALSE, extractfn =
                     NULL, multisession = FALSE, ncores = 1)

library(scrbook)
simcaps<-function (N, sigma, alpha0, alpha2, K) {
  ntraps=130
  trapx<-rep(seq(.25, 9.75, .25), 39)
  trapy<-rep(seq(.25, 9.75, .25), 39)
  traplocs <- cbind(sort(trapx), trapy)
  traplocstrue<-traplocs[sample(nrow(traplocs), ntraps),]
  Dmat <- e2dist(traplocstrue, traplocstrue)
  Xl <- 0
  Xu <- 10
  Yl <- 0
  Yu <- 10
  sx <- runif(N, Xl, Xu)
  sy <- runif(N, Yl, Yu)
  S <- cbind(sx, sy)
  D <- e2dist(S, traplocstrue)
  alpha1<--1/(2*sigma*sigma)
  Ycat <- matrix(NA, nrow = N, ncol = K)
  Xlag <- matrix(0, nrow = N, ncol = K + 1)
  Xlag[, 1] <- rep(0, N)
  for (i in 1:N) {
    for (k in 1:K) {
      lp <- alpha0 + alpha2 * Xlag[i, k] + alpha1 * D[i, ] * D[i, ]
      cp <- exp(c(lp, 0))
      cp <- cp/sum(cp)
      Ycat[i, k] <- sample(1:(ntraps+1), 1, prob = cp)
      if (Ycat[i, k] <= ntraps) 
        Xlag[i, (k + 1):ncol(Xlag)] <- 1
    }
  }
  captured <- apply(Ycat <= ntraps, 1, sum)
  captured <- captured > 0
  Ycat <- Ycat[captured, ]
  Xlag <- Xlag[captured, ]
  reencounter = Xlag[, 1:K]
  X <- traplocstrue
  K <- ncol(reencounter)
  Ncap <- S[captured, ]
  Nmissed<-S[!captured,]
  return(Ncap = nrow(Ncap))
  }
  
  
sim=rep(NA,500)
for (i in 1:500){
sim[i]<-simcaps(N=60, sigma=1, alpha0=-1.3, alpha2=-3, K = 10, ssbuff = 2)}


sim2=rep(NA,500)
for (i in 1:500){
  sim2[i]<-simcaps(N=60, sigma=1, alpha0=-2, alpha2=-3, K = 10, ssbuff = 3)}



simtotal=rep(NA, 5000)
for (i in 1:5000){
  simtotal[i]<-simcaps(N=60, sigma=runif(1, 0.5, 1.5), alpha0 = runif(1, -3, -1), alpha2=-3, K=10)}

simtotal3=rep(NA, 50000)
for (i in 1:50000){
  simtotal3[i]<-simcaps(N=60, sigma=runif(1, 0.5, 1.5), alpha0 = runif(1, -3.47, -1.09), alpha2=-3, K=10)}

simcaps<-function (N, sigma1, alpha0, alpha2, K) {
  ntraps=130
  trapx<-rep(seq(.25, 9.75, .25), 39)
  trapy<-rep(seq(.25, 9.75, .25), 39)
  traplocs <- cbind(sort(trapx), trapy)
  traplocstrue<-traplocs[sample(nrow(traplocs), ntraps),]
  Dmat <- e2dist(traplocstrue, traplocstrue)
  Xl <- 0
  Xu <- 10
  Yl <- 0
  Yu <- 10
  sx <- runif(N, Xl, Xu)
  sy <- runif(N, Yl, Yu)
  S <- cbind(sx, sy)
  D <- e2dist(S, traplocstrue)
  alpha1<--1/(2*sigma*sigma)
  Ycat <- matrix(NA, nrow = N, ncol = K)
  Xlag <- matrix(0, nrow = N, ncol = K + 1)
  Xlag[, 1] <- rep(0, N)
  for (i in 1:N) {
    for (k in 1:K) {
      lp <- alpha0 + alpha2 * Xlag[i, k] + alpha1 * D[i, ] * D[i, ]
      cp <- exp(c(lp, 0))
      cp <- cp/sum(cp)
      Ycat[i, k] <- sample(1:(ntraps+1), 1, prob = cp)
      if (Ycat[i, k] <= ntraps) 
        Xlag[i, (k + 1):ncol(Xlag)] <- 1
    }
  }
  captured <- apply(Ycat <= ntraps, 1, sum)
  captured <- captured > 0
  Ycat <- Ycat[captured, ]
  Xlag <- Xlag[captured, ]
  reencounter = Xlag[, 1:K]
  X <- traplocstrue
  K <- ncol(reencounter)
  Ncap <- S[captured, ]
  Nmissed<-S[!captured,]
  return(Ncap = nrow(Ncap))
}