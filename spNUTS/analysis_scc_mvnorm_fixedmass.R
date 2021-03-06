dataDir <- 'data/2020Jan/'
mvnorm.spNUTS1.result <- function(eps, N, astar, repNo, M, burnout, adaptive.massInv) {
    filename <- paste(dataDir,'mvnorm_','spNUTS1','_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_N_',N,'_astar_',astar,'_repNo_',repNo,'.RData',sep='')
    load(filename)
    return(c(as.numeric(runtime)[1], moveprob, mean(nfail), mean(unlist(log2len)), mean(unlist(sym)), minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, minESS.mcmcse, minESS.mcmcse.perSec, avgESS.mcmcse, avgESS.mcmcse.perSec, eps.trace[length(eps.trace)]))
}
mvnorm.spNUTS1.results <- function(eps, N, astar, repNoRange, M, burnout, adaptive.massInv) {
    output <- t(sapply(repNoRange, function(x) mvnorm.spNUTS1.result(eps, N, astar, x, M, burnout, adaptive.massInv)))
    colnames(output) <- c('runtime', 'moveprob', 'nfail', 'log2len', 'sym', 'minESS.coda', 'minESSps.coda', 'avgESS.coda', 'avgESSps.coda', 'minESS.mcmcse', 'minESSps.mcmcse', 'avgESS.mcmcse', 'avgESSps.mcmcse', 'final_eps')
    output
}
mvnorm.spNUTS2.result <- function(eps, astar, repNo, M, burnout, adaptive.massInv) {
    filename <- paste(dataDir,'mvnorm_','spNUTS2','_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_astar_',astar,'_repNo_',repNo,'.RData',sep='')
    load(filename)
    return(c(as.numeric(runtime)[1], moveprob, mean(unlist(log2len)), mean(unlist(sym),na.rm=TRUE), minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, minESS.mcmcse, minESS.mcmcse.perSec, avgESS.mcmcse, avgESS.mcmcse.perSec, eps.trace[length(eps.trace)]))
}
mvnorm.spNUTS2.results <- function(eps, astar, repNoRange, M, burnout, adaptive.massInv) {
    output <- t(sapply(repNoRange, function(x) mvnorm.spNUTS2.result(eps, astar, x, M, burnout, adaptive.massInv)))
    colnames(output) <- c('runtime', 'moveprob', 'log2len', 'sym', 'minESS.coda', 'minESSps.coda', 'avgESS.coda', 'avgESSps.coda', 'minESS.mcmcse', 'minESSps.mcmcse', 'avgESS.mcmcse', 'avgESSps.mcmcse', 'final_eps')
    output
}
mvnorm.NUTS.result <- function(eps, astar, repNo, M, burnout, adaptive.massInv) {
    filename <- paste(dataDir,'mvnorm_','NUTS','_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_astar_',astar,'_repNo_',repNo,'.RData',sep='')
    load(filename)
    return(c(as.numeric(runtime)[1], moveprob, mean(unlist(log2len)), minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, minESS.mcmcse, minESS.mcmcse.perSec, avgESS.mcmcse, avgESS.mcmcse.perSec, eps.trace[length(eps.trace)]))
}
mvnorm.NUTS.results <- function(eps, astar, repNoRange, M, burnout, adaptive.massInv) {
    output <- t(sapply(repNoRange, function(x) mvnorm.NUTS.result(eps, astar, x, M, burnout, adaptive.massInv)))
    colnames(output) <- c('runtime', 'moveprob', 'log2len', 'minESS.coda', 'minESSps.coda', 'avgESS.coda', 'avgESSps.coda', 'minESS.mcmcse', 'minESSps.mcmcse', 'avgESS.mcmcse', 'avgESSps.mcmcse', 'final_eps')
    output
}
mvnorm.HMC.result <- function(eps, astar, repNo, M, burnout, adaptive.massInv) {
    filename <- paste(dataDir,'mvnorm_','HMC','_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_astar_',astar,'_M_',M,'_repNo_',repNo,'.RData',sep='')
    load(filename)
    return(c(as.numeric(runtime)[1], moveprob, minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, minESS.mcmcse, minESS.mcmcse.perSec, avgESS.mcmcse, avgESS.mcmcse.perSec, eps.trace[length(eps.trace)]))
}
mvnorm.HMC.results <- function(eps, astar, repNoRange, M, burnout, adaptive.massInv) {
    output <- t(sapply(repNoRange, function(x) mvnorm.HMC.result(eps, astar, x, M, burnout, adaptive.massInv)))
    colnames(output) <- c('runtime', 'moveprob', 'minESS.coda', 'minESSps.coda', 'avgESS.coda', 'avgESSps.coda', 'minESS.mcmcse', 'minESSps.mcmcse', 'avgESS.mcmcse', 'avgESSps.mcmcse', 'final_eps')
    output
}
mvnorm.spHMC.result <- function(eps, astar, repNo, M, burnout, adaptive.massInv) {
    filename <- paste(dataDir,'mvnorm_','spHMC','_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_l_',50,'_astar_',astar,'_M_',M,'_repNo_',repNo,'.RData',sep='')
    load(filename)
    return(c(as.numeric(runtime)[1], moveprob, minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, minESS.mcmcse, minESS.mcmcse.perSec, avgESS.mcmcse, avgESS.mcmcse.perSec, eps.trace[length(eps.trace)]))
}
mvnorm.spHMC.results <- function(eps, astar, repNoRange, M, burnout, adaptive.massInv) {
    output <- t(sapply(repNoRange, function(x) mvnorm.spHMC.result(eps, astar, x, M, burnout, adaptive.massInv)))
    colnames(output) <- c('runtime', 'moveprob', 'minESS.coda', 'minESSps.coda', 'avgESS.coda', 'avgESSps.coda', 'minESS.mcmcse', 'minESSps.mcmcse', 'avgESS.mcmcse', 'avgESSps.mcmcse', 'final_eps')
    output
}


eps <- 0.01; repNoRange <- 1:10; M_HMC <- 200200; M <- 20200; burnout <- 200; adaptive.massInv <- FALSE

mvnorm_HMC_.95 <- mvnorm.HMC.results(eps, .95, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_HMC_.85 <- mvnorm.HMC.results(eps, .85, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_HMC_.75 <- mvnorm.HMC.results(eps, .75, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_HMC_.65 <- mvnorm.HMC.results(eps, .65, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_HMC_.55 <- mvnorm.HMC.results(eps, .55, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_HMC_.45 <- mvnorm.HMC.results(eps, .45, repNoRange, M_HMC, burnout, adaptive.massInv)
##mvnorm_HMC_.35 <- mvnorm.HMC.results(eps, .35, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_spHMC_.95 <- mvnorm.spHMC.results(eps, .95, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_spHMC_.85 <- mvnorm.spHMC.results(eps, .85, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_spHMC_.75 <- mvnorm.spHMC.results(eps, .75, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_spHMC_.65 <- mvnorm.spHMC.results(eps, .65, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_spHMC_.55 <- mvnorm.spHMC.results(eps, .55, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_spHMC_.45 <- mvnorm.spHMC.results(eps, .45, repNoRange, M_HMC, burnout, adaptive.massInv)
##mvnorm_spHMC_.35 <- mvnorm.spHMC.results(eps, .35, repNoRange, M_HMC, burnout, adaptive.massInv)
mvnorm_spNUTS1_.95_1 <- mvnorm.spNUTS1.results(eps, 1, .95, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.95_5 <- mvnorm.spNUTS1.results(eps, 5, .95, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.85_1 <- mvnorm.spNUTS1.results(eps, 1, .85, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.85_5 <- mvnorm.spNUTS1.results(eps, 5, .85, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.75_1 <- mvnorm.spNUTS1.results(eps, 1, .75, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.75_5 <- mvnorm.spNUTS1.results(eps, 5, .75, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.65_1 <- mvnorm.spNUTS1.results(eps, 1, .65, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.65_5 <- mvnorm.spNUTS1.results(eps, 5, .65, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.55_1 <- mvnorm.spNUTS1.results(eps, 1, .55, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.55_5 <- mvnorm.spNUTS1.results(eps, 5, .55, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.45_1 <- mvnorm.spNUTS1.results(eps, 1, .45, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS1_.45_5 <- mvnorm.spNUTS1.results(eps, 5, .45, repNoRange, M, burnout, adaptive.massInv)
##mvnorm_spNUTS1_.35_1 <- mvnorm.spNUTS1.results(eps, 1, .35, repNoRange, M, burnout, adaptive.massInv)
##mvnorm_spNUTS1_.35_5 <- mvnorm.spNUTS1.results(eps, 5, .35, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS2_.95 <- mvnorm.spNUTS2.results(eps, .95, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS2_.85 <- mvnorm.spNUTS2.results(eps, .85, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS2_.75 <- mvnorm.spNUTS2.results(eps, .75, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS2_.65 <- mvnorm.spNUTS2.results(eps, .65, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS2_.55 <- mvnorm.spNUTS2.results(eps, .55, repNoRange, M, burnout, adaptive.massInv)
mvnorm_spNUTS2_.45 <- mvnorm.spNUTS2.results(eps, .45, repNoRange, M, burnout, adaptive.massInv)
##mvnorm_spNUTS2_.35 <- mvnorm.spNUTS2.results(eps, .35, repNoRange, M, burnout, adaptive.massInv)
mvnorm_NUTS_.95 <- mvnorm.NUTS.results(eps, .95, repNoRange, M, burnout, adaptive.massInv)
mvnorm_NUTS_.85 <- mvnorm.NUTS.results(eps, .85, repNoRange, M, burnout, adaptive.massInv)
mvnorm_NUTS_.75 <- mvnorm.NUTS.results(eps, .75, repNoRange, M, burnout, adaptive.massInv)
mvnorm_NUTS_.65 <- mvnorm.NUTS.results(eps, .65, repNoRange, M, burnout, adaptive.massInv)
mvnorm_NUTS_.55 <- mvnorm.NUTS.results(eps, .55, repNoRange, M, burnout, adaptive.massInv)
mvnorm_NUTS_.45 <- mvnorm.NUTS.results(eps, .45, repNoRange, M, burnout, adaptive.massInv)
##mvnorm_NUTS_.35 <- mvnorm.NUTS.results(eps, .35, repNoRange, M, burnout, adaptive.massInv)

minESSps.coda <- c(
mvnorm_HMC_.95[,"minESSps.coda"],
mvnorm_HMC_.85[,"minESSps.coda"],
mvnorm_HMC_.75[,"minESSps.coda"],
mvnorm_HMC_.65[,"minESSps.coda"],
mvnorm_HMC_.55[,"minESSps.coda"],
mvnorm_HMC_.45[,"minESSps.coda"],
##mvnorm_HMC_.35[,"minESSps.coda"],
mvnorm_spHMC_.95[,"minESSps.coda"],
mvnorm_spHMC_.85[,"minESSps.coda"],
mvnorm_spHMC_.75[,"minESSps.coda"],
mvnorm_spHMC_.65[,"minESSps.coda"],
mvnorm_spHMC_.55[,"minESSps.coda"],
mvnorm_spHMC_.45[,"minESSps.coda"],
##mvnorm_spHMC_.35[,"minESSps.coda"],
mvnorm_NUTS_.95[,"minESSps.coda"],
mvnorm_NUTS_.85[,"minESSps.coda"],
mvnorm_NUTS_.75[,"minESSps.coda"],
mvnorm_NUTS_.65[,"minESSps.coda"],
mvnorm_NUTS_.55[,"minESSps.coda"],
mvnorm_NUTS_.45[,"minESSps.coda"],
##mvnorm_NUTS_.35[,"minESSps.coda"],
mvnorm_spNUTS1_.95_1[,"minESSps.coda"],
mvnorm_spNUTS1_.85_1[,"minESSps.coda"],
mvnorm_spNUTS1_.75_1[,"minESSps.coda"],
mvnorm_spNUTS1_.65_1[,"minESSps.coda"],
mvnorm_spNUTS1_.55_1[,"minESSps.coda"],
mvnorm_spNUTS1_.45_1[,"minESSps.coda"],
##mvnorm_spNUTS1_.35_1[,"minESSps.coda"],
mvnorm_spNUTS1_.95_5[,"minESSps.coda"],
mvnorm_spNUTS1_.85_5[,"minESSps.coda"],
mvnorm_spNUTS1_.75_5[,"minESSps.coda"],
mvnorm_spNUTS1_.65_5[,"minESSps.coda"],
mvnorm_spNUTS1_.55_5[,"minESSps.coda"],
mvnorm_spNUTS1_.45_5[,"minESSps.coda"],
##mvnorm_spNUTS1_.35_5[,"minESSps.coda"],
mvnorm_spNUTS2_.95[,"minESSps.coda"],
mvnorm_spNUTS2_.85[,"minESSps.coda"],
mvnorm_spNUTS2_.75[,"minESSps.coda"],
mvnorm_spNUTS2_.65[,"minESSps.coda"],
mvnorm_spNUTS2_.55[,"minESSps.coda"],
mvnorm_spNUTS2_.45[,"minESSps.coda"]
##mvnorm_spNUTS2_.35[,"minESSps.coda"]
)

avgESSps.coda <- c(
mvnorm_HMC_.95[,"avgESSps.coda"],
mvnorm_HMC_.85[,"avgESSps.coda"],
mvnorm_HMC_.75[,"avgESSps.coda"],
mvnorm_HMC_.65[,"avgESSps.coda"],
mvnorm_HMC_.55[,"avgESSps.coda"],
mvnorm_HMC_.45[,"avgESSps.coda"],
##mvnorm_HMC_.35[,"avgESSps.coda"],
mvnorm_spHMC_.95[,"avgESSps.coda"],
mvnorm_spHMC_.85[,"avgESSps.coda"],
mvnorm_spHMC_.75[,"avgESSps.coda"],
mvnorm_spHMC_.65[,"avgESSps.coda"],
mvnorm_spHMC_.55[,"avgESSps.coda"],
mvnorm_spHMC_.45[,"avgESSps.coda"],
##mvnorm_spHMC_.35[,"avgESSps.coda"],
mvnorm_NUTS_.95[,"avgESSps.coda"],
mvnorm_NUTS_.85[,"avgESSps.coda"],
mvnorm_NUTS_.75[,"avgESSps.coda"],
mvnorm_NUTS_.65[,"avgESSps.coda"],
mvnorm_NUTS_.55[,"avgESSps.coda"],
mvnorm_NUTS_.45[,"avgESSps.coda"],
##mvnorm_NUTS_.35[,"avgESSps.coda"],
mvnorm_spNUTS1_.95_1[,"avgESSps.coda"],
mvnorm_spNUTS1_.85_1[,"avgESSps.coda"],
mvnorm_spNUTS1_.75_1[,"avgESSps.coda"],
mvnorm_spNUTS1_.65_1[,"avgESSps.coda"],
mvnorm_spNUTS1_.55_1[,"avgESSps.coda"],
mvnorm_spNUTS1_.45_1[,"avgESSps.coda"],
##mvnorm_spNUTS1_.35_1[,"avgESSps.coda"],
mvnorm_spNUTS1_.95_5[,"avgESSps.coda"],
mvnorm_spNUTS1_.85_5[,"avgESSps.coda"],
mvnorm_spNUTS1_.75_5[,"avgESSps.coda"],
mvnorm_spNUTS1_.65_5[,"avgESSps.coda"],
mvnorm_spNUTS1_.55_5[,"avgESSps.coda"],
mvnorm_spNUTS1_.45_5[,"avgESSps.coda"],
##mvnorm_spNUTS1_.35_5[,"avgESSps.coda"],
mvnorm_spNUTS2_.95[,"avgESSps.coda"],
mvnorm_spNUTS2_.85[,"avgESSps.coda"],
mvnorm_spNUTS2_.75[,"avgESSps.coda"],
mvnorm_spNUTS2_.65[,"avgESSps.coda"],
mvnorm_spNUTS2_.55[,"avgESSps.coda"],
mvnorm_spNUTS2_.45[,"avgESSps.coda"]
##mvnorm_spNUTS2_.35[,"avgESSps.coda"]
)

minESSps.mcmcse <- c(
mvnorm_HMC_.95[,"minESSps.mcmcse"],
mvnorm_HMC_.85[,"minESSps.mcmcse"],
mvnorm_HMC_.75[,"minESSps.mcmcse"],
mvnorm_HMC_.65[,"minESSps.mcmcse"],
mvnorm_HMC_.55[,"minESSps.mcmcse"],
mvnorm_HMC_.45[,"minESSps.mcmcse"],
##mvnorm_HMC_.35[,"minESSps.mcmcse"],
mvnorm_spHMC_.95[,"minESSps.mcmcse"],
mvnorm_spHMC_.85[,"minESSps.mcmcse"],
mvnorm_spHMC_.75[,"minESSps.mcmcse"],
mvnorm_spHMC_.65[,"minESSps.mcmcse"],
mvnorm_spHMC_.55[,"minESSps.mcmcse"],
mvnorm_spHMC_.45[,"minESSps.mcmcse"],
##mvnorm_spHMC_.35[,"minESSps.mcmcse"],
mvnorm_NUTS_.95[,"minESSps.mcmcse"],
mvnorm_NUTS_.85[,"minESSps.mcmcse"],
mvnorm_NUTS_.75[,"minESSps.mcmcse"],
mvnorm_NUTS_.65[,"minESSps.mcmcse"],
mvnorm_NUTS_.55[,"minESSps.mcmcse"],
mvnorm_NUTS_.45[,"minESSps.mcmcse"],
##mvnorm_NUTS_.35[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.95_1[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.85_1[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.75_1[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.65_1[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.55_1[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.45_1[,"minESSps.mcmcse"],
##mvnorm_spNUTS1_.35_1[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.95_5[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.85_5[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.75_5[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.65_5[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.55_5[,"minESSps.mcmcse"],
mvnorm_spNUTS1_.45_5[,"minESSps.mcmcse"],
##mvnorm_spNUTS1_.35_5[,"minESSps.mcmcse"],
mvnorm_spNUTS2_.95[,"minESSps.mcmcse"],
mvnorm_spNUTS2_.85[,"minESSps.mcmcse"],
mvnorm_spNUTS2_.75[,"minESSps.mcmcse"],
mvnorm_spNUTS2_.65[,"minESSps.mcmcse"],
mvnorm_spNUTS2_.55[,"minESSps.mcmcse"],
mvnorm_spNUTS2_.45[,"minESSps.mcmcse"]
##mvnorm_spNUTS2_.35[,"minESSps.mcmcse"]
)

avgESSps.mcmcse <- c(
mvnorm_HMC_.95[,"avgESSps.mcmcse"],
mvnorm_HMC_.85[,"avgESSps.mcmcse"],
mvnorm_HMC_.75[,"avgESSps.mcmcse"],
mvnorm_HMC_.65[,"avgESSps.mcmcse"],
mvnorm_HMC_.55[,"avgESSps.mcmcse"],
mvnorm_HMC_.45[,"avgESSps.mcmcse"],
##mvnorm_HMC_.35[,"avgESSps.mcmcse"],
mvnorm_spHMC_.95[,"avgESSps.mcmcse"],
mvnorm_spHMC_.85[,"avgESSps.mcmcse"],
mvnorm_spHMC_.75[,"avgESSps.mcmcse"],
mvnorm_spHMC_.65[,"avgESSps.mcmcse"],
mvnorm_spHMC_.55[,"avgESSps.mcmcse"],
mvnorm_spHMC_.45[,"avgESSps.mcmcse"],
##mvnorm_spHMC_.35[,"avgESSps.mcmcse"],
mvnorm_NUTS_.95[,"avgESSps.mcmcse"],
mvnorm_NUTS_.85[,"avgESSps.mcmcse"],
mvnorm_NUTS_.75[,"avgESSps.mcmcse"],
mvnorm_NUTS_.65[,"avgESSps.mcmcse"],
mvnorm_NUTS_.55[,"avgESSps.mcmcse"],
mvnorm_NUTS_.45[,"avgESSps.mcmcse"],
##mvnorm_NUTS_.35[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.95_1[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.85_1[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.75_1[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.65_1[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.55_1[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.45_1[,"avgESSps.mcmcse"],
##mvnorm_spNUTS1_.35_1[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.95_5[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.85_5[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.75_5[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.65_5[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.55_5[,"avgESSps.mcmcse"],
mvnorm_spNUTS1_.45_5[,"avgESSps.mcmcse"],
##mvnorm_spNUTS1_.35_5[,"avgESSps.mcmcse"],
mvnorm_spNUTS2_.95[,"avgESSps.mcmcse"],
mvnorm_spNUTS2_.85[,"avgESSps.mcmcse"],
mvnorm_spNUTS2_.75[,"avgESSps.mcmcse"],
mvnorm_spNUTS2_.65[,"avgESSps.mcmcse"],
mvnorm_spNUTS2_.55[,"avgESSps.mcmcse"],
mvnorm_spNUTS2_.45[,"avgESSps.mcmcse"]
##mvnorm_spNUTS2_.35[,"avgESSps.mcmcse"]
)

runtime.perIter <- c(
mvnorm_HMC_.95[,"runtime"]/M_HMC,
mvnorm_HMC_.85[,"runtime"]/M_HMC,
mvnorm_HMC_.75[,"runtime"]/M_HMC,
mvnorm_HMC_.65[,"runtime"]/M_HMC,
mvnorm_HMC_.55[,"runtime"]/M_HMC,
mvnorm_HMC_.45[,"runtime"]/M_HMC,
##mvnorm_HMC_.35[,"runtime"]/M_HMC,
mvnorm_spHMC_.95[,"runtime"]/M_HMC,
mvnorm_spHMC_.85[,"runtime"]/M_HMC,
mvnorm_spHMC_.75[,"runtime"]/M_HMC,
mvnorm_spHMC_.65[,"runtime"]/M_HMC,
mvnorm_spHMC_.55[,"runtime"]/M_HMC,
mvnorm_spHMC_.45[,"runtime"]/M_HMC,
##mvnorm_spHMC_.35[,"runtime"]/M_HMC,
mvnorm_NUTS_.95[,"runtime"]/M,
mvnorm_NUTS_.85[,"runtime"]/M,
mvnorm_NUTS_.75[,"runtime"]/M,
mvnorm_NUTS_.65[,"runtime"]/M,
mvnorm_NUTS_.55[,"runtime"]/M,
mvnorm_NUTS_.45[,"runtime"]/M,
##mvnorm_NUTS_.35[,"runtime"/M],
mvnorm_spNUTS1_.95_1[,"runtime"]/M,
mvnorm_spNUTS1_.85_1[,"runtime"]/M,
mvnorm_spNUTS1_.75_1[,"runtime"]/M,
mvnorm_spNUTS1_.65_1[,"runtime"]/M,
mvnorm_spNUTS1_.55_1[,"runtime"]/M,
mvnorm_spNUTS1_.45_1[,"runtime"]/M,
##mvnorm_spNUTS1_.35_1[,"runtime"]/M,
mvnorm_spNUTS1_.95_5[,"runtime"]/M,
mvnorm_spNUTS1_.85_5[,"runtime"]/M,
mvnorm_spNUTS1_.75_5[,"runtime"]/M,
mvnorm_spNUTS1_.65_5[,"runtime"]/M,
mvnorm_spNUTS1_.55_5[,"runtime"]/M,
mvnorm_spNUTS1_.45_5[,"runtime"]/M,
##mvnorm_spNUTS1_.35_5[,"runtime"]/M,
mvnorm_spNUTS2_.95[,"runtime"]/M,
mvnorm_spNUTS2_.85[,"runtime"]/M,
mvnorm_spNUTS2_.75[,"runtime"]/M,
mvnorm_spNUTS2_.65[,"runtime"]/M,
mvnorm_spNUTS2_.55[,"runtime"]/M,
mvnorm_spNUTS2_.45[,"runtime"]/M
##mvnorm_spNUTS2_.35[,"runtime"]/M
)

log2len <- c( # final value of adaptively tuned leapfrog jump size
mvnorm_NUTS_.95[,"log2len"],
mvnorm_NUTS_.85[,"log2len"],
mvnorm_NUTS_.75[,"log2len"],
mvnorm_NUTS_.65[,"log2len"],
mvnorm_NUTS_.55[,"log2len"],
mvnorm_NUTS_.45[,"log2len"],
##mvnorm_NUTS_.35[,"log2len"],
mvnorm_spNUTS1_.95_1[,"log2len"],
mvnorm_spNUTS1_.85_1[,"log2len"],
mvnorm_spNUTS1_.75_1[,"log2len"],
mvnorm_spNUTS1_.65_1[,"log2len"],
mvnorm_spNUTS1_.55_1[,"log2len"],
mvnorm_spNUTS1_.45_1[,"log2len"],
##mvnorm_spNUTS1_.35_1[,"log2len"],
mvnorm_spNUTS1_.95_5[,"log2len"],
mvnorm_spNUTS1_.85_5[,"log2len"],
mvnorm_spNUTS1_.75_5[,"log2len"],
mvnorm_spNUTS1_.65_5[,"log2len"],
mvnorm_spNUTS1_.55_5[,"log2len"],
mvnorm_spNUTS1_.45_5[,"log2len"],
##mvnorm_spNUTS1_.35_5[,"log2len"],
mvnorm_spNUTS2_.95[,"log2len"],
mvnorm_spNUTS2_.85[,"log2len"],
mvnorm_spNUTS2_.75[,"log2len"],
mvnorm_spNUTS2_.65[,"log2len"],
mvnorm_spNUTS2_.55[,"log2len"],
mvnorm_spNUTS2_.45[,"log2len"]
##mvnorm_spNUTS2_.35[,"log2len"]
)

astarValues <- c(0.95,0.85,0.75,0.65,0.55,0.45)
algorithms <- c('HMC', 'spHMC', 'NUTS','spNUTS1 (N=1)','spNUTS1 (N=5)','spNUTS2')

## min ESS (coda) per second vs. astar and algorithms
gen_pdf <- TRUE; cex <- 1.3
if(gen_pdf) { pdf(file='figures/mvnormModel/minESSperSec_coda_astar_algo_fixedmass_mvnorm.pdf', width=10, height=3) }
par(mar=c(4,4.5,.3,.3))
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps.coda, pch=3, xaxt='n',yaxt='n', xlab='',ylab='', cex=1, ylim=c(0,max(minESSps.coda)))
axis(side=2, cex.axis=cex)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cex)
axis(side=1, line=2, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cex)
mtext(side=2, line=3, text='minimum ESS per second', cex=cex)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
if(gen_pdf) { dev.off() }

## avg ESS (coda) per second vs. astar and algorithms
gen_pdf <- TRUE; cex <- 1.3
if(gen_pdf) { pdf(file='figures/mvnormModel/avgESSperSec_coda_astar_algo_fixedmass_mvnorm.pdf', width=10, height=3) }
par(mar=c(4,4.5,.3,.3))
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),avgESSps.coda, pch=3, xaxt='n',yaxt='n', xlab='',ylab='', cex=1, ylim=c(0,max(avgESSps.coda)))
axis(side=2, cex.axis=cex)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cex)
axis(side=1, line=2, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cex)
mtext(side=2, line=3, text='average ESS per second', cex=cex)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
if(gen_pdf) { dev.off() }

## min and avg ESS (coda) per second vs. astar and algorithms
gen_pdf <- TRUE; cexlab <- 1.2; cexaxis <- 1.6; cexalgo <- 1.8; cexpt <- 6
if(gen_pdf) { pdf(file='figures/mvnormModel/ESSperSec_astar_coda_algo_fixedmass_mvnorm.pdf', width=10, height=6.5) }
par(mar=c(4,4.5,0.3,.3)); par(mfrow=c(3,1), oma=c(2.8,0,2,0))
## minimum ESS per second
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps.coda, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', xaxs="i", cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), ylim=c(min(minESSps.coda[-c((2*length(astarValues)-1)*length(repNoRange)+repNoRange,(3*length(astarValues)-1)*length(repNoRange)+repNoRange)]),max(minESSps.coda)), log='y') # exclude the NUTS minESSps.coda for astar=0.45 (very small) to get a better range of y-axis
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='minimum ESS /sec', cex=cexlab)
axis(side=3, line=-1, outer=TRUE, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cexalgo)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
## average ESS per second
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),avgESSps.coda, pch='.', xaxt='n',yaxt='n', xaxs="i", xlab='',ylab='', cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), log='y')
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='average ESS /sec', cex=cexlab)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
## runtime in seconds
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)), runtime, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', xaxs="i", cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), log='y')
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='runtime', cex=cexlab)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
mtext(side=1, line=5, text='Target acceptance probability', cex=cexlab)
if(gen_pdf) { dev.off() }


## min and avg ESS (coda and mcmcse) per second vs. astar and algorithms
gen_pdf <- TRUE; cexlab <- 1.2; cexaxis <- 1.6; cexalgo <- 1.7; cexpt <- 6
if(gen_pdf) { pdf(file='figures/mvnormModel/ESSperSec_coda_mcmcse_astar_algo_fixedmass_mvnorm_newM.pdf', width=10, height=10) }
par(mar=c(4,4.5,0.3,.3)); par(mfrow=c(5,1), oma=c(2.8,2,2,0))
## minimum ESS per second (coda)
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps.coda, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', xaxs="i", cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), ylim=c(min(minESSps.coda[-c((2*length(astarValues)-1)*length(repNoRange)+repNoRange,(3*length(astarValues)-1)*length(repNoRange)+repNoRange)]),max(minESSps.coda)), log='y') # exclude the NUTS minESSps.coda for astar=0.45 (very small) to get a better range of y-axis
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='min ESS/sec\n(coda)', cex=cexlab)
axis(side=3, line=-1, outer=TRUE, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cexalgo)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
## minimum ESS per second (mcmcse)
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps.mcmcse, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', xaxs="i", cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), ylim=c(min(minESSps.mcmcse[-c((2*length(astarValues)-1)*length(repNoRange)+repNoRange,(3*length(astarValues)-1)*length(repNoRange)+repNoRange)]),max(minESSps.mcmcse)), log='y') # exclude the NUTS minESSps.mcmcse for astar=0.45 (very small) to get a better range of y-axis
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='min ESS/sec\n(batch)', cex=cexlab)
axis(side=3, line=-1, outer=TRUE, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cexalgo)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
## average ESS per second (coda)
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),avgESSps.coda, pch='.', xaxt='n',yaxt='n', xaxs="i", xlab='',ylab='', cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), log='y')
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='avg ESS/sec\n(coda)', cex=cexlab)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
## average ESS per second (mcmcmse)
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),avgESSps.mcmcse, pch='.', xaxt='n',yaxt='n', xaxs="i", xlab='',ylab='', cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), log='y')
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='avg ESS/sec\n(batch)', cex=cexlab)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
## runtime in seconds
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)), runtime.perIter, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', xaxs="i", cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), log='y')
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='runtime/M', cex=cexlab)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
mtext(side=1, line=5, text='Target acceptance probability', cex=cexlab)
if(gen_pdf) { dev.off() }








final_eps <- c( # final value of adaptively tuned leapfrog jump size
mvnorm_HMC_.95[,"final_eps"],
mvnorm_HMC_.85[,"final_eps"],
mvnorm_HMC_.75[,"final_eps"],
mvnorm_HMC_.65[,"final_eps"],
mvnorm_HMC_.55[,"final_eps"],
mvnorm_HMC_.45[,"final_eps"],
##mvnorm_HMC_.35[,"final_eps"],
mvnorm_spHMC_.95[,"final_eps"],
mvnorm_spHMC_.85[,"final_eps"],
mvnorm_spHMC_.75[,"final_eps"],
mvnorm_spHMC_.65[,"final_eps"],
mvnorm_spHMC_.55[,"final_eps"],
mvnorm_spHMC_.45[,"final_eps"],
##mvnorm_spHMC_.35[,"final_eps"],
mvnorm_NUTS_.95[,"final_eps"],
mvnorm_NUTS_.85[,"final_eps"],
mvnorm_NUTS_.75[,"final_eps"],
mvnorm_NUTS_.65[,"final_eps"],
mvnorm_NUTS_.55[,"final_eps"],
mvnorm_NUTS_.45[,"final_eps"],
##mvnorm_NUTS_.35[,"final_eps"],
mvnorm_spNUTS1_.95_1[,"final_eps"],
mvnorm_spNUTS1_.85_1[,"final_eps"],
mvnorm_spNUTS1_.75_1[,"final_eps"],
mvnorm_spNUTS1_.65_1[,"final_eps"],
mvnorm_spNUTS1_.55_1[,"final_eps"],
mvnorm_spNUTS1_.45_1[,"final_eps"],
##mvnorm_spNUTS1_.35_1[,"final_eps"],
mvnorm_spNUTS1_.95_5[,"final_eps"],
mvnorm_spNUTS1_.85_5[,"final_eps"],
mvnorm_spNUTS1_.75_5[,"final_eps"],
mvnorm_spNUTS1_.65_5[,"final_eps"],
mvnorm_spNUTS1_.55_5[,"final_eps"],
mvnorm_spNUTS1_.45_5[,"final_eps"],
##mvnorm_spNUTS1_.35_5[,"final_eps"],
mvnorm_spNUTS2_.95[,"final_eps"],
mvnorm_spNUTS2_.85[,"final_eps"],
mvnorm_spNUTS2_.75[,"final_eps"],
mvnorm_spNUTS2_.65[,"final_eps"],
mvnorm_spNUTS2_.55[,"final_eps"],
mvnorm_spNUTS2_.45[,"final_eps"]
##mvnorm_spNUTS2_.35[,"final_eps"]
)

## average trajectory length
atl <- matrix(log2len,nrow=length(repNoRange)) * matrix(final_eps,nrow=length(repNoRange))[,-(1:(2*length(astarValues)))]
















## producing plots for NESS 2019
minESSps <- c(
mvnorm_HMC_.95[,"minESSps"],
mvnorm_HMC_.85[,"minESSps"],
mvnorm_HMC_.75[,"minESSps"],
mvnorm_HMC_.65[,"minESSps"],
mvnorm_HMC_.55[,"minESSps"],
mvnorm_HMC_.45[,"minESSps"],
mvnorm_spHMC_.95[,"minESSps"],
mvnorm_spHMC_.85[,"minESSps"],
mvnorm_spHMC_.75[,"minESSps"],
mvnorm_spHMC_.65[,"minESSps"],
mvnorm_spHMC_.55[,"minESSps"],
mvnorm_spHMC_.45[,"minESSps"]
)
astarValues <- c(0.95,0.85,0.75,0.65,0.55,0.45)
algorithms <- c('HMC', 'spHMC')
## min ESS per second vs. astar and algorithms
gen_pdf <- TRUE; cexlab <- 1; cexaxis <- 1; cexalgo <- 1
if(gen_pdf) { pdf(file='figures/mvnormModel/minESSperSec_astar_HMC_spHMC_fixedmass.pdf', width=6, height=3) }
par(mar=c(4,4.5,1.5,.3))
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', cex=3, ylim=c(0,max(minESSps)))
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='minimum ESS /sec', cex=cexlab)
axis(side=3, line=-1, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cexalgo)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
mtext(side=1, line=3, text='Target acceptance probability', cex=cexlab)
if(gen_pdf) { dev.off() }

minESSps <- c(
mvnorm_NUTS_.95[,"minESSps"],
mvnorm_NUTS_.85[,"minESSps"],
mvnorm_NUTS_.75[,"minESSps"],
mvnorm_NUTS_.65[,"minESSps"],
mvnorm_NUTS_.55[,"minESSps"],
mvnorm_NUTS_.45[,"minESSps"],
##mvnorm_NUTS_.35[,"minESSps"],
mvnorm_spNUTS1_.95_5[,"minESSps"],
mvnorm_spNUTS1_.85_5[,"minESSps"],
mvnorm_spNUTS1_.75_5[,"minESSps"],
mvnorm_spNUTS1_.65_5[,"minESSps"],
mvnorm_spNUTS1_.55_5[,"minESSps"],
mvnorm_spNUTS1_.45_5[,"minESSps"]
##mvnorm_spNUTS1_.35_5[,"minESSps"],
)
astarValues <- c(0.95,0.85,0.75,0.65,0.55,0.45)
algorithms <- c('NUTS', 'spNUTS')
## min ESS per second vs. astar and algorithms
gen_pdf <- TRUE; cexlab <- 1; cexaxis <- 1; cexalgo <- 1
if(gen_pdf) { pdf(file='figures/mvnormModel/minESSperSec_astar_NUTS_spNUTS1_fixedmass.pdf', width=6, height=3) }
par(mar=c(4,4.5,1.5,.3))
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', cex=3, ylim=c(0,max(minESSps)))
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='minimum ESS /sec', cex=cexlab)
axis(side=3, line=-1, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cexalgo)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
mtext(side=1, line=3, text='Target acceptance probability', cex=cexlab)
if(gen_pdf) { dev.off() }






minESS.coda <- matrix(
    c(
mvnorm_HMC_.95[,"minESS.coda"],
mvnorm_HMC_.85[,"minESS.coda"],
mvnorm_HMC_.75[,"minESS.coda"],
mvnorm_HMC_.65[,"minESS.coda"],
mvnorm_HMC_.55[,"minESS.coda"],
mvnorm_HMC_.45[,"minESS.coda"],
##mvnorm_HMC_.35[,"minESS.coda"],
mvnorm_spHMC_.95[,"minESS.coda"],
mvnorm_spHMC_.85[,"minESS.coda"],
mvnorm_spHMC_.75[,"minESS.coda"],
mvnorm_spHMC_.65[,"minESS.coda"],
mvnorm_spHMC_.55[,"minESS.coda"],
mvnorm_spHMC_.45[,"minESS.coda"],
##mvnorm_spHMC_.35[,"minESS.coda"],
mvnorm_NUTS_.95[,"minESS.coda"],
mvnorm_NUTS_.85[,"minESS.coda"],
mvnorm_NUTS_.75[,"minESS.coda"],
mvnorm_NUTS_.65[,"minESS.coda"],
mvnorm_NUTS_.55[,"minESS.coda"],
mvnorm_NUTS_.45[,"minESS.coda"],
##mvnorm_NUTS_.35[,"minESS.coda"],
mvnorm_spNUTS1_.95_1[,"minESS.coda"],
mvnorm_spNUTS1_.85_1[,"minESS.coda"],
mvnorm_spNUTS1_.75_1[,"minESS.coda"],
mvnorm_spNUTS1_.65_1[,"minESS.coda"],
mvnorm_spNUTS1_.55_1[,"minESS.coda"],
mvnorm_spNUTS1_.45_1[,"minESS.coda"],
##mvnorm_spNUTS1_.35_1[,"minESS.coda"],
mvnorm_spNUTS1_.95_5[,"minESS.coda"],
mvnorm_spNUTS1_.85_5[,"minESS.coda"],
mvnorm_spNUTS1_.75_5[,"minESS.coda"],
mvnorm_spNUTS1_.65_5[,"minESS.coda"],
mvnorm_spNUTS1_.55_5[,"minESS.coda"],
mvnorm_spNUTS1_.45_5[,"minESS.coda"],
##mvnorm_spNUTS1_.35_5[,"minESS.coda"],
mvnorm_spNUTS2_.95[,"minESS.coda"],
mvnorm_spNUTS2_.85[,"minESS.coda"],
mvnorm_spNUTS2_.75[,"minESS.coda"],
mvnorm_spNUTS2_.65[,"minESS.coda"],
mvnorm_spNUTS2_.55[,"minESS.coda"],
mvnorm_spNUTS2_.45[,"minESS.coda"]
##mvnorm_spNUTS2_.35[,"minESS.coda"]
)
, 60)

runtime <- matrix(runtime, 60)

minESSps.coda <- matrix(minESSps.coda, 60)
