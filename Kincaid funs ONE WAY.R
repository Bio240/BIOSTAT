# Kincaid funs ONE WAY.R     28 April 2018    by   Dwight Kincaid, PhD

# library(DescTools); library(BayesFactor); library(VCA) 
# library(DAAG);      library(QuantPsyc);   library(treeperm) 
# library(moments);   library(nortest);     library(pwr)      

# ---------------------------------------------------------------------------
#  At the end of the code please see a DUMMY, SIMULATED one-way DATA SET to 
#  DEMO FUNCTION CALLS. It's 4 groups, total N=78, from normal distribution.
#
#  Also at the end is the 'Age at First Walking' data set, a=4, total N=23.
#
#  DEMO CODE after each function definition should be studied and run. 
# ---------------------------------------------------------------------------
#
#     NEW FUNCTIONS                                                          
#
#  1. R.sq.oneway(y, x)             # R-square% in one-way Anoval; explained & error
#  2. vc.oneway(y, x, verbose)      # traditional Variance Components (VC) only
#                                     appropriate for random effects one-way Anova                                                    
#  3. CohenD.oneway(y, x, verbose)  # Cohen's D, standardized effect size in one-way Anova                                   
#
#  4. boot.dataframe(dframe, n.rows)              # data frame bootstrapped by row index                                
#  5. boot.vc.oneway(y, x, data.name, NS)         # bootstrap CI of VC
#  6. bootCohenD.oneway(y, x, NS, data.name)      # bootstrap CI of Cohen's D
#  7. bootRsq.oneway(y, x, NS, data.name, export.distribution) # bootstrap CI of R-sq
#
#  8. VCci.oneway(y, x)     # Variance Components in Anova by REML; library(VCA)                                         
#
#  9. new.stats(y, digits, ylab, norm.test)  # univariate descriptive stats
# 10. CImeans.normaltheory(y, x, prob)       # per group in a one-way layout
#                                              by normal theory (t-distribution)
# 11. EDA.graphs.one.sample(y, var.name)     # univariate EDA by annotated graphs
# 12. EDA.plots.oneway(y, x)                 # per group in a one-way layout                                                                               
# 13. NormalQQ.oneway(y, x, data.name)       # per group in a one-way layout; normality
#                                              assessed by sim data after John Maindonald
# 14. Normal.sim.kernel.band(y, NS) # univariate assessment of normality by sim kernel band
# 15. BayesFactorAnova.oneway(y, x, NS, fixed.effects) # anovaBF() in library(BayesFactor)                                   
#
# 16. treeperm.oneway.anova(y, x, NS, type, graph.F)  # permutation test for Anova F                                       
# 17. Kincaid.perm.oneway(y, x, NS, describe.null.F, returnF) # permutation test for Anova F                   
#
# 18. powercurve.normaltheory.oneway(y, x, low.n, high.n, by.n, matrix.out)                                                                
# 19. powercurve.Cohen.oneway(y, x, data.name, low.n, high.n, by.n, matrix.out)         
# 20. boot.power.of.test.oneway(y, x, data.name, NS, N.boot, n.breaks, return.distribution)
# 21. boot.powercurve.oneway(y, x, NS, data.name, low.n, high.n, by.n, matrix.out))
#
# 22. sample.size.barplots(y, x, data.name, response.name, factor.name, ... )
# 23. one.way.histograms.lattice(y, x, data.name, response.name, hist.type, ... ) 
# 24. means.barplot(y, x, err.bar, prob, verbose, bar.lwd, bar.length, bar.col, ...)
# 25. Kincaid.stripchart.error.bars(y, x, err.bar, prob, verbose, 
#        bar.lwd, bar.col, bar.length, stat.cex, ...) 
#
# 26. sim.data.oneway.Normal(a, total.n, balanced, grand.mean, SD, seed, 
#        effect.size, precision, verbose, graph)
# 27. sim.toy.oneway.data(a, n, my.seed, effect, graphs)   # just as it says
# 28. anova.by.stats(n, m, s, digits)  # one-way Anova from N, Mean, 
# 29. Kincaid.error.bars.by.stats()    # many args; means with error bars: SD, SEM, CI

# --  other functions --
#
# my.pause()    idle()    thick.line()    thin.line()    mytick()    skw()    krt() 
#
# --  ideas on functions to add --
#
# x. one.way.sensitivity()       x. assumptions.oneway.anova()  x. my.time.stamp()                                                                     
# x. TukeyOutliersOneway()       x. bootBayesFactor.oneway()    x. nonpar.tests.oneway()    
# x. unpaired.t.test.by.stats()  x. perm.nonpar.tests()
      
 

R.sq.oneway <- function(y, x){  

  # ------------------------------------------------------------------------
  #  R-square(%) in one-way Anova 
  #  
  #  ARGS    y    numeric response vector, NAs OK
  #          x    corresponding factor (grouping) vector, NAs OK
  #
  #  Although trivial, it's a gentle, classroom intro to writing functions
  #  and is easily modified and improved by additional arguments; for
  #  instance, one could add an arg, CI=TRUE to enable bootstrap 95%
  #  confidence intervals for explained R-square in one-way Anova.
  # ------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # arg & data checking
    length(y)>1, length(x)>1, length(y) == length(x)) 

  out         <- anova(lm(y ~ x))            # one-way anova object
  total.SS    <- out[1,2] + out[2,2]         # SS groups + SS within
  R.sq.among  <- 100 * out[1,2] / total.SS   # explained R-sq percent
  R.sq.within <- 100 * out[2,2] / total.SS   # error R-sq percent

  A <- rbind(R.sq.among, R.sq.within) # a [2,1] matrix
  rownames(A) <- c("R-square groups", "R-square unexplained")
  colnames(A) <- "percent"
  return(A)  # a 2 row by 1 column matrix that prints nicely

}  # end of function definition


# demo the call


# sim data: a=4, total n=78
#
# R.sq.oneway(my.frame$response, my.frame$drug)
#
#                      percent
# R-square groups      12.91943
# R-square unexplained 87.08057
#

#
# age at first walking data: a=4, total n=23
#
#
# R.sq.oneway(months, treatment)
#                      percent
# R-square groups      25.2753
# R-square unexplained 74.7247


# library(DescTools)  # Note: R-sq = eta-squared in one-way Anova model
# EtaSq( aov(months ~ treatment), anova=TRUE )   # eta-squared = 25.3 %
#             eta.sq eta.sq.part       SS df       MS        F         p
# treatment 0.252753    0.252753 14.77781  3 4.925936 2.142222 0.1285456
# Residuals 0.747247          NA 43.68958 19 2.299452       NA        NA
# 





vc.oneway <- function(y, x, verbose=FALSE) { 

  # -------------------------------------------------------------------
  #  Variance Components in 'random effects' one-way Anova
  #
  #  ARGS    y         numeric response vector, NAs OK
  #          x         corresponding factor vector, NAs OK
  #          verbose   FALSE returns a vector, TRUE also returns text
  #
  #  function RETURNS textbook traditional variance components as 
  #  %, calculated from one-way ANOVA table after Sokal & Rohlf (1995)
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ------------------------------------------------------------------- 

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.logical(verbose)) 

  dfr <- na.omit(data.frame(y, x))  # remove cases (rows) with NA in y or x
  y   <- dfr$y; x <- dfr$x          # reconstitute y & x, for convenience 

  out  <- anova(lm(y ~ x))    # one-way anova object
  MSw  <- out[2,3]            # mean square WITHIN groups
  MSa  <- out[1,3]            # mean square AMONG groups
  a    <- nlevels(x)          # number of groups
  n    <- length(y)           # total N

  n.sq <- sum(tapply(y, x, length)^2) # sum of squared, group sample sizes    
  No   <- (1/(a-1))*(n-(n.sq/n))      # 'average' sample size for VC purposes
  vc   <- (MSa-MSw)/No                   # 'added variance component'
  vc.among    <- 100*(vc/(MSw+vc))       # % variance component AMONG groups
  if(vc.among < 0) vc.among <- 0         # if < 0, set to 0; CUSTOMARY 
  vc.within   <- 100*(MSw/(MSw+vc))      # % variance component WITHIN groups
  if(vc.within > 100) vc.within <- 100   # if > 100, set to 100; CUSTOMARY

  if(verbose==TRUE) cat("\n\nVariance components(%) in one-way Anova\nGROUP & ERROR\t")
  c(vc.among, vc.within)  # RETURNS vector of % explained, % unexplained VC

} # end. of function definition

# demo the call

# sim data set: a=4, total n=78
# vc.oneway(my.frame$response, my.frame$drug)     
# vc.oneway(my.frame$response, my.frame$drug, verbose=TRUE)  
  
 
# age at first walking data: a=4, total n=23      
# vc.oneway(months, treatment)
# vc.oneway(months, treatment, verbose=TRUE)




CohenD.oneway <- function(y, x, verbose=TRUE){  

  # ------------------------------------------------------------------------
  #  Cohen's D in one-way Anova. Standardized effect size as the difference
  #  between means / pooled SD. Cohen, J. 1988. Statistical power analysis 
  #  for the behavioral sciences, 2nd ed. Lawrence Erlbaum, Hillsdale, NJ.
  #
  #   ARGS    y         numeric response vector, NAs OK
  #           x         corresponding factor vector, NAs OK
  #           verbose   FALSE returns Cohen's D, TRUE also returns text
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),            # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.logical(verbose)) 

  dfr <- na.omit(data.frame(y, x)) # remove all cases (rows) with NA in y or x
  y  <- dfr$y; x <- dfr$x          # reconstitute y & x, for convenience 

  out    <- anova(lm(y ~ x))          # one-way anova object
  RMSE   <- sqrt( out[2,3] )          # root mean square error, RMSE, pooled SD
  diff   <- sqrt( out[1,2]/length(y)) # 'raw' effect size
  CohenD <- diff / RMSE               # raw effect size standardized

  if(verbose == TRUE){
    cat("\nCohen's D =", CohenD, "\t(standardized effect size)\n")
    cat("  In one-way Anova and after Cohen (1988): \n")
    cat("  D = 0.1 (small),  D = 0.25 (medium),  D = 0.4 (large effect size)\n\n")
    cat("Cohen, J. 1988. Statistical power analysis for the behavioral sciences,\n")
    cat("2nd ed. Lawrence Erlbaum, Hillsdale, NJ.\n\n")
  }

  if(verbose == FALSE) return(CohenD)  # returns vector of N=1 as needed

} # end of function definition


# demo the call

# sim data: a=4, total n=78

# CohenD.oneway(my.frame$response, my.frame$drug)  # returns nice text to R Console
#
# Cohen's D = 0.3751714   (standardized effect size)

#

#
# age at first walking data: a=4, total n=23
#

# CohenD.oneway(months, treatment, verbose=FALSE)  # returns vector of N=1
#
# [1] 0.5286022
#
# CohenD.oneway(months, treatment, verbose=TRUE) 
#
#Cohen's D = 0.5286022   (standardized effect size)
#  In one-way Anova and after Cohen (1988): 
#  D = 0.1 (small),  D = 0.25 (medium),  D = 0.4 (large effect size)
#
#Cohen, J. 1988. Statistical power analysis for the behavioral sciences,
#2nd ed. Lawrence Erlbaum, Hillsdale, NJ.


#
# Note: across R there are many Cohen's D funs but mostly only for pairwise comparisons
#



boot.dataframe <- function(dframe, n.rows=0 ){

  # ----------------------------------------------------------------------
  #  BOOTSTRAP an entire data frame by resampling its rows
  #
  #  ARGS    dframe   data frame object of any size -- rows or columns
  #          n.rows   desired N for each bootstrap sample; n.rows=0 
  #                   causes N to be the same as N of observed rows
  #
  #  function RETURNS a bootstrapped data frame; NAs not removed. 
  #  Warning: Dealing with data frame objects is slow, if speed is
  #  an issue when function is called iteratively. 
  # 
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------------

  stopifnot(is.data.frame(dframe), is.numeric(n.rows))   # data & arg checking
  n.rows <- as.integer(n.rows)                           # just in case

  indx <- nrow(dframe)          # N of rows, an integer
  ifelse(n.rows == 0, 
    bf <- dframe[ sample(indx, replace=TRUE), ],               # boot the rows
    bf <- dframe[ sample(indx, size=n.rows, replace=TRUE), ])  # boot the rows 
  
  return(bf)  # the bootstrap resampled data frame

} # end of function definition


# demo the call
#
# age at walking data set: a=4, total n=23
#

#summary(DF)  # original data
#     months        treatment
# Min.   : 9.00   Active :6  
# 1st Qu.:10.00   Passive:6  
# Median :11.50   None   :6  
# Mean   :11.35   Control:5  
# 3rd Qu.:12.50              
# Max.   :15.00 

# NOW get a BOOTSTRAPPED data frame from 'Age at walking' data frame
# boot.frame <- boot.dataframe(DF)
# summary(boot.frame)
#     months        treatment
# Min.   : 9.00   Active :8  
# 1st Qu.: 9.50   Passive:3  
# Median :11.50   None   :8  
# Mean   :11.45   Control:4  
# 3rd Qu.:13.12              
# Max.   :15.00   

# boot.frame <- boot.dataframe(DF, n.rows=50) # n.rows specified
# summary(boot.frame)





boot.vc.oneway <- function(y, x, data.name=NULL, NS=1e3){

  # ---------------------------------------------------------------------
  #  Percentile BOOTSTRAP 96% CI of Variance Components in 
  #  random effects, one-way Anova
  #
  #  ARGS   y,x  numeric response vector, paired with x, factor vector
  #         NS   desired N of bootstrap resamples
  #  data.name   descriptive string for data set, e.g., "Age at walking"
  #
  #  function calls boot.dataframe() and then vc.oneway() 
  #  Note: bootstrap not warranted if N per group is too small.
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ---------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.numeric(NS) ) 

  dfr <- na.omit(data.frame(y, x))  # remove cases (rows) with NA in y or x   
  thin.line <- function(N=25) cat(rep("-", N), sep="", "\n")  # prints "-"

  cat("\n"); thin.line(65)
  cat("BOOTSTRAP Variance Components in one-way Anova\n\n")
  cat("Data set: ", data.name, "\n\n")
  cat("NS =", NS, "bootstrap resamples from observed data frame\n\n")
  flush.console()  # R Console buffering - dump text

  boot.time <- proc.time()  # begin timer; boot data frame NS times; save VC each time
    out <- replicate(NS, {bf <- boot.dataframe(dfr); vc.oneway(bf$y, bf$x)})
  boot.time <- ((proc.time() - boot.time)/60)[3]  # CPU minutes
   
  group.boot <- out[1, ]
  error.boot <- out[2, ]

  group.CI <- quantile(group.boot, c(.025, .975), na.rm=TRUE)  # percentile bootstrap
  error.CI <- quantile(error.boot, c(.025, .975), na.rm=TRUE)  # 95% CI of VC
 
  group.vc <- vc.oneway(y, x)[1]  # observed vc
  error.vc <- vc.oneway(y, x)[2]

  # stats for explained, boot distribution
  grp.range <- range(group.boot)  
  grp.m     <- mean(group.boot); grp.sd <- sd(group.boot)
  grp.skew  <- skw(group.boot)  # Kincaid fun
 
  cat("Using the percentile method\nthe 95% confidence intervals are --\n\n")
  cat("Explained VC%\t\tobserved:",group.vc, "\n"); print(group.CI)

  cat("\n STATS\tfor bootstrap distribution of explained VC%\n")
  cat("\tmin, max:", grp.range[1], ",", grp.range[2], "\n")
  cat("\tmean:", grp.m, " SD:", grp.sd, " skew:", grp.skew, "\n")

  cat("\nError VC%\t\tobserved:",error.vc, "\n");      print(error.CI)
  cat("\nCPU minutes:", round(boot.time,3), "\tfor bootstrap resampling\n")

  cat("\nOutput from:  Kincaid.boot.vc.oneway()  by  Dwight Kincaid, PhD\n")
  thin.line(65); cat("\n")

  par(mfrow=c(2,1))  # 2 rows by 1 column of graphs
    hist(group.boot, main="Bootstrap distribution of EXPLAINED VC%", breaks="scott", 
         xlab="VC%", cex.main=1)
    abline(v=group.vc, col="red", lwd=2)
    mtext("red line at observed VC", cex=.8, side=3, col="red")
    abline(v=quantile(group.boot, .025), col="blue", lwd=2)
    abline(v=quantile(group.boot, .975), col="blue", lwd=2)
    mtext("blue lines at .025 and .975 quantiles", cex=.8,side=1,line=4,adj=0,col="blue")
    mtext(paste("NS =", NS), cex=.8, side=3, adj=0)

    hist(error.boot, main="Bootstrap distribution of RESIDUAL VC%", breaks="scott",
         xlab="VC%", cex.main=1)
    abline(v=error.vc, col="red", lwd=2)
    mtext("red line at observed VC", cex=.8, side=3, col="red")
    abline(v=quantile(error.boot, .025), col="blue", lwd=2)
    abline(v=quantile(error.boot, .975), col="blue", lwd=2)
    mtext("blue lines at .025 and .975 quantiles", cex=.8,side=1,line=4,adj=0,col="blue")
    mtext(paste("NS =", NS), cex=.8, side=3, adj=0)
    mtext(paste("Minutes:", round(boot.time, 3)), side=1, adj=1,cex=.8,line=4,col="blue")
    mtext(paste("Data set:", data.name), side=4, cex=.7, col="blue")
  par(mfrow=c(1,1))  # reset to default of 1 graph in window
  
} # end of function definition


# demo the function

# sim data: a=4, total n=78
# boot.vc.oneway(response, drug, data.name="sim data", NS=1e2)
# boot.vc.oneway(response, drug, data.name="sim data", NS=2e3)


#
# age at walking data: a=4, total n=23
#

# boot.vc.oneway(months, treatment, data.name="Age at walking", NS=1e3)





bootCohenD.oneway <- function(y, x, NS, data.name=NULL){

  # ---------------------------------------------------------------
  #  BOOTSTRAP  CI  for  Cohen's d  in one-way Anova
  #
  #   ARGS    y   numeric response vector, NAs OK
  #           x   corresponding factor vector, NAs OK
  #          NS   desired N of bootstrap samples
  #   data.name   descriptive string, e.g., "Age at Walking"
  #
  #  Nonparametric percentile bootstrap confidence intervals.
  #  Bootstrap not warranted if N per group is too small.
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ---------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.numeric(NS) ) 

  dfr <- na.omit( data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x
  n  <- length(y)                    # total sample size

  Dfun <- function(yy, xx){               # Cohen's d from Anova table
    out  <- anova(lm(yy ~ xx))            # one-way Anova
    RMSE <- sqrt(out[2,3])                # root mean square error
    diff <- sqrt(out[1,2] / length(yy))   # 'raw' effect size
    return(diff/RMSE)                     # standardized effect size, Cohen's d
  }

  obs.D <- Dfun(y, x)  # observed Cohen's d

  cat("\n\n"); cat(rep("-", 30), "\n")
  cat("OUTPUT from bootCohenD.oneway( )  by  Dwight Kincaid, PhD\n")
  cat("\nObserved, Cohen's D:", obs.D, "in one-way Anova", "\n\n")
  flush.console()  # R Console buffering - dump text

  boot.time <- proc.time()  # begin timer
    boot.D  <- replicate( NS, { indx <- sample(n, replace=TRUE);
                                 Dfun(y[indx], x[indx]) } )
  boot.time <- ((proc.time() - boot.time)/60)[3]  # minutes
  
  CI90 <- quantile(boot.D, c(.05,  .95))   # 90% lower and upper bounds
  CI95 <- quantile(boot.D, c(.025, .975))  # 95%
  CI99 <- quantile(boot.D, c(.005, .995))  # 99%

  max.freq <- max(hist(boot.D, plot=FALSE)$counts)

  hist(boot.D, main="Bootstrap distribution of Cohen's D in one-way Anova",
    xlab="Cohen's D", cex.main=1.2, col="wheat", breaks="scott",
    font.lab=2, cex.lab=1.25, font.main=4)

  abline(v=obs.D, col="red", lwd=3)
  mtext(paste("red line: observed Cohen's D =", round(obs.D, 5)), 
        cex=1, side=3, col="red")
  mtext(paste("NS =", NS, "resamples"), cex=.8, side=1, line=3, adj=0)
  mtext(paste("Minutes required:", round(boot.time,2) ), cex=.8, side=1, line=3,adj=1)
  lines(x=c(CI95[1], CI95[1]), y=c(0, max.freq/1.6), lty=3, lwd=3)
  lines(x=c(CI95[2], CI95[2]), y=c(0, max.freq/1.6), lty=3, lwd=3)
  mtext(paste("Dotted lines at percentile bootstrap 95% CI:", 
        round(CI95[1],4), ",", round(CI95[2],4)), 
        side=1, adj=0, line=4, cex=.8)
  mtext(data.name, side=4)

  # stats for boot.D
  m <- mean(boot.D); s <- sd(boot.D); rangeD <- range(boot.D); skewD <- skw(boot.D)
  cat("STATS\tfor bootstrap distribution of Cohen's D\n")
  cat("\tn:", NS, "\t min, max:", rangeD[1], ",", rangeD[2], "\n")
  cat("\tmean:", m, " SD:", s, " skew:", skewD, "\n\n")

  cat("NS:", NS, "bootstrap samples\tminutes:", round(boot.time,3), "\n\n")
  cat("Nonparametric, Percentile Bootstrap Confidence Intervals\n")
  cat("for Cohen's D in one-way Anova\n\n")

  A <- cbind(CI90, CI95, CI99)
  colnames(A) <- c("90% CI", "95% CI", "99% CI")
  rownames(A) <- c("lower-bound", "upper-bound")
  print(A)

  cat("\nTime stamp: ", date(), "\n")
  cat(rep("-", 30), "\n\n")

}  # end function definition


# demo the fun
#
# age at first walking data: a=4, total n=23
#
# bootCohenD.oneway(months, treatment, NS=1e3, data.name="Age at Walking")

# bootCohenD.oneway(response, drug, NS=1e2, data.name="sim data")


bootRsq.oneway <- function(y, x, NS=1e3, 
  data.name="response", export.distribution=FALSE){

  # ---------------------------------------------------------------
  #  BOOTSTRAP  CI  for  explained R-square % in one-way Anova
  #
  #   ARGS    y            numeric response vector, NAs OK
  #           x            corresponding factor vector, NAs OK
  #          NS            desired N of bootstrap samples
  #   data.name            descriptive string, e.g., "Age at Walking"
  #   export.distribution  logical; in case you want to process it
  #
  #  Nonparametric percentile bootstrap confidence intervals.
  #  Bootstrap not warranted if N per group is too small.
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ---------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),      # data & arg checking
    length(y) > 1, length(x) > 1, length(y) == length(x), 
    is.numeric(NS), is.character(data.name), is.logical(export.distribution)) 

  dfr <- na.omit( data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x vectors
  n   <- length(y)             # total n       a CONSTANT

  # Note:  .fastRq() is a private function that can be used independently
  # if you know what you are doing; it does no checking, for instance

  .fastRsq <- function(y, x){    # returns explained R-sq% in oneway Anova
    n.g  <- tapply(y, x, length) # n per group
    gd.m <- sum(y)/length(y)     # grand mean
    SS.t <- sum((y - gd.m)^2)    # total SS
    SS.g <- sum(n.g * ((tapply(y, x, sum)/n.g) - gd.m)^2) # group SS
    100 * SS.g/SS.t }

  obs.Rsq <- .fastRsq(y, x) # observed, explained R-square in one-way Anova

  # the 2 variables are bootstrapped by index
  timer <- proc.time()
    out <- replicate(NS, { indx <- sample(n, replace=TRUE);
                    .fastRsq(y[indx], x[indx]) } )
  timer <- proc.time() - timer
  timer <- round(timer[[3]]/60, 3) # minutes elapsed by bootstrap

  CI <- quantile(out, c(.025, .975), na.rm=TRUE)
  CI <- round(CI, 3)

  density <- hist(out, breaks="scott", plot=FALSE)$density
  y.max <- 1.2 * max(density) # only for plotting

  hist(out, freq=FALSE,
    ylim=c(0, y.max), main=NULL,
    breaks="scott", col="lightblue", border="white",
    xlab="R-square in oneway Anova", font.lab=2
  )

  abline(v=obs.Rsq, col=2, lwd=3.5)
 
  temp <- y.max/3
  lines(c(CI[1], CI[1]), c(0, temp), col=1, lwd=3.5)
  lines(c(CI[2], CI[2]), c(0, temp), col=1, lwd=3.5)

  mtext(paste("minutes:", timer), side=1, line=4, adj=0, cex=.8, col=4)
  mtext(paste("NS =", NS, "resamples"), side=1, line=4, adj=1, cex=.8, col=4)
  mtext(paste("red line: observed R-square =", round(obs.Rsq,3), "%"), 
    side=1, line=2, adj=0, col=2, cex=.8)
  mtext("black lines: 95% CI", side=1, line=2, adj=1, cex=.8)
  mtext("percentile bootstrap ", side=1, line=4, cex=.8, col=4)
  mtext(paste("95% CI: ", CI[1], ",", CI[2]), line=0, adj=0, col=2, cex=.9)
  mtext(paste("data.name =", data.name), side=4, cex=.9)

  mtext("Bootstrap distribution of explained R-square %", line=2,adj=0,font=2,cex=1.2) 
  mtext("STATISTICS for the Bootstrap Distribution of R-sq", line=0, adj=1, cex=.65)

  out <- out[!is.na(out)]  # remove missing values
  m  <- round(mean(out), 2); s <- round(sd(out), 2)
  mn <- round(min(out), 2); mx <- round(max(out), 2)

  mtext(paste("mean:", m, " SD:", s, " min:", mn, " max:", mx), line=-1,adj=1,cex=.65)

  if(export.distribution == TRUE) return(out)

} # end function definition


# demo the function

# bootRsq.oneway(months, treatment, NS=1e4, data.name="Age at Walking")

#a <- bootRsq.oneway(montha, treatment, NS=1e4, export.distribution=TRUE)
#str(a); summary(a); quantile(a, .995)




VCci.oneway <- function(y, x){

  # -----------------------------------------------------------------                                                           
  #  Variance Components by REML for random effects, one-way Anova  
  #  
  #  ARGS  y,x  numeric response vector, paired with x: factor vector
  #             NAs OK for y and x
  #                                              
  #  library(VCA) was new in Spring 2015 and uses REML as in SAS.
  #  Note: VCA can handle more complicated Anova models.
  #                       
  #  This convenience function simply calls anovaVCA() and 
  #  VCAinference(), reduces the output and adds explanatory text.
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu                                                                         
  # -----------------------------------------------------------------                                                           

  require(VCA)

  stopifnot(is.numeric(y), is.factor(x),          # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x)) 

  dfr <- na.omit( data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x

  thin.line <- function(N=25) cat(rep("-", N), sep="", "\n")  # prints "-"
  cat("\n"); thin.line(72)
  cat("95% CI for VARIANCE COMPONENTS in random effects one-way Anova.\n")
  cat("Appropriate only for RANDOM EFFECTS model Anova. From library(VCA) which\n")
  cat("was new in 2015; anovaVCA() and VCAinference() are called.\n")
  cat("Calculations not based on resampling but on theory and REML.\n")
  cat("For documentation:  >?anovaVCA  >?VCAinference\n"); thin.line(72)
  flush.console()  # R Console buffering - dump text

  # get observed variance components and their percentages
  VCA.obj <- anovaVCA(y ~ x, Data=dfr ) 
  print(VCA.obj)

  total.VC <- sum(VCA.obj$VCoriginal)

  thin.line(34)
  cat("\n95% confidence intervals for the\n")
  cat("two variance components: groups('x') and residual('error')\n")
  cat("Note: See Kincaid's translation\n")
  cat(paste("total VC =", round(total.VC,6)), "\n\n")

  # calculate 95% CI for variance components, two-sided
  VCA.ci <- VCAinference(VCA.obj, alpha=.05, VarVC=TRUE)
  out <- VCA.ci$ConfInt$VC  # extract from the huge object
  print(out$OneSided)  # ------------------------------------- 

  LB.groups <- 100 * out$OneSided$LCL[2] / total.VC  # 95% CI as % of total VC
  UB.groups <- 100 * out$OneSided$UCL[2] / total.VC

  LB.error <- 100 * out$OneSided$LCL[3] / total.VC  # 95% CI as % of total VC
  UB.error <- 100 * out$OneSided$UCL[3] / total.VC

  if(LB.groups < 0) LB.groups <- 0; if(UB.groups > 100) UB.groups <- 100
  if(LB.error < 0)  LB.error <- 0;  if(UB.error > 100)  UB.error <- 100

  cat("\nKINCAID translation of the above one-sided table.\n\n")

  A <- rbind(LB.groups, UB.groups, LB.error, UB.error)
  colnames(A) <- "95% CI as % of total VC"
  print(A); cat("\n"); thin.line(34)
  print(out$TwoSided)   # ------------------------------------- 

  LB.groups <- 100 * out$TwoSided$LCL[2] / total.VC  # 95% CI as a % of total VC
  UB.groups <- 100 * out$TwoSided$UCL[2] / total.VC

  LB.error <- 100 * out$TwoSided$LCL[3] / total.VC  # 95% CI as a % of total VC
  UB.error <- 100 * out$TwoSided$UCL[3] / total.VC

  if(LB.groups < 0) LB.groups <- 0; if(UB.groups > 100) UB.groups <- 100
  if(LB.error < 0)  LB.error <- 0;  if(UB.error > 100)  UB.error <- 100

  cat("\nKINCAID translation of the above two-sided table.\n\n")

  A <- rbind(LB.groups, UB.groups, LB.error, UB.error)
  colnames(A) <- "95% CI as % of total VC"
  print(A); cat("\n"); thin.line(34)

  # nothing is returned but text to R Console, but that could change...

} # end of function definition


# demo the code
#
# age at walking data: a=4, total n=23
#

# VCci.oneway(months, treatment) 
#


# development code, e.g. how to call the 2 functions in VCA, on their own

# sim data: a=4, total n=78

# VCA.obj <- anovaVCA( response ~ drug, Data=my.frame)
# VCA.obj

# VCA.ci <- VCAinference(VCA.obj, alpha=.05, VarVC=TRUE)
# VCA.ci
# str(VCA.ci) # huge object

# out <- VCA.ci$ConfInt$VC$TwoSided
# out[ , 2:3]



new.stats <- function(y, digits=3, ylab=NULL, norm.test=FALSE ){

  # -----------------------------------------------------------------
  #  Stats and normality tests for single, numeric samples
  #
  #  ARGS     y: numeric vector     digits: roundoff  
  #        ylab: data name       norm.test: logical, normality tests   
  #
  #  requires packages:  moments, nortest -- both for gof tests
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # -----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y)>1,            # arg & data checking
    is.numeric(digits), is.logical(norm.test) ) 

  y <- y[!is.na(y)]                     # remove missing values from numeric vector
  if(length(y) < 8) norm.test=FALSE     # don't test for normality on small sample
  yy <- y + rnorm(length(y), 0, .00001) # white noise to break ties in ks.test()

  skw <- function(y){                  # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)          # assumes no NA
    (sum(y^3) / length(y)) / sd(y)^3}

  krt <- function(y){                  # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)          # assumes no NA
    (sum(y^4) / length(y)) / sd(y)^4 - 3}
  
  require(moments, quietly=TRUE); require(nortest, quietly=TRUE) # GOF tests

  N <- length(y); Mean <- mean(y); SEM <- sd(y)/sqrt(N)
  LB.95 <- Mean - qt(p=.975, df=N-1) * SEM  # 95% CI of mean by normal theory
  UB.95 <- Mean + qt(p=.975, df=N-1) * SEM

  A <- rbind(N, Mean, sd(y), 100*sd(y)/Mean, min(y), max(y), skw(y), krt(y), 
             SEM, LB.95, UB.95, median(y), quantile(y,.25), quantile(y,.75), IQR(y) )
  A <- round(A, digits)
  colnames(A) <- ylab
  rownames(A) <-c("N", "Mean", "SD", "CV%", "Min", "Max", "skew", "kurtosis", 
             "SEM", "LB.95%CI.mean", "UB.95%CI.mean", "Median", "25th", "75th", "IQR" )

   # Kincaid's revision of cvm.test() in library(nortest) to kill warnings on p, etc.
   # Cramer-von Mises test of normality, library(nortest), >?cvm.test
   my.cvm.test <- function (x){ 
     x <- sort(x); n <- length(x)
     if (n < 8) stop("sample size must be greater than 7")
     p <- pnorm((x-mean(x))/sd(x)); W <- (1/(12*n) + sum((p-(2 * seq(1:n)-1)/(2*n))^2))
     WW <- (1+0.5/n)*W
     if (WW < 0.0275) { pval <- 1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)}
       else if (WW < 0.051) { pval <- 1 - exp(-5.903 + 179.546 * WW - 1515.29 * WW^2)}
       else if (WW < 0.092) { pval <- exp(0.886 - 31.62 * WW + 10.897 * WW^2)}
       else if (WW < 1.1)   { pval <- exp(1.111 - 34.242 * WW + 12.832 * WW^2)}
       else { pval <- 7.37e-10 } #warning("p-value < 7.37e-10 cannot be computed")

     return(pval)
   } # end of my.cvm.test() function definition

  if(norm.test == TRUE){

    if(N <= 5000 & N > 7){
       B <- rbind(shapiro.test(y)$p.value, jarque.test(y)$p.value, ad.test(y)$p.value,
            my.cvm.test(y), sf.test(y)$p.value, pearson.test(y)$p.value,
            agostino.test(y)$p.value, anscombe.test(y)$p.value,
            bonett.test(y)$p.value, ks.test(yy, "pnorm", mean(yy), sd(yy))$p.value,
            lillie.test(y)$p.value)

       colnames(B) <- "p-value"
       rownames(B) <- c("Shapiro-Wilk test","Robust Jarque-Berra test",
                     "Anderson-Darling test",
                     "Cramer-von Mises test", "Shapiro-Francia test", "Pearson test",
                     "D'Agostino test", "Anscombe-Glynn test","Bonett-Seier test", 
                     "Kolmogorov-Smirnov test", "Lilifors Kolmogorov-Smirnov test") 
       return(list(statistics=A, goodness.of.fit.tests.for.normality=B)) }


    if(N > 5000 & N <= 46340){
       B <- rbind(jarque.test(y)$p.value, ad.test(y)$p.value, cvm.test(y)$p.value,
            pearson.test(y)$p.value, agostino.test(y)$p.value, anscombe.test(y)$p.value,
            bonett.test(y)$p.value, ks.test(yy, "pnorm", mean(yy), sd(yy))$p.value,
            lillie.test(y)$p.value) 

       colnames(B) <- "p-value"
       rownames(B) <- c("Robust Jarque-Berra test", "Anderson-Darling test",
                     "Cramer-von Mises test", "Pearson test",
                     "D'Agostino test", "Anscombe-Glynn test","Bonett-Seier test", 
                     "Kolmogorov-Smirnov test", "Lilifors Kolmogorov-Smirnov test") 
       return(list(statistics=A, goodness.of.fit.tests.for.normality=B)) }

    if(N > 46340){
       B <- rbind(jarque.test(y)$p.value, ad.test(y)$p.value, cvm.test(y)$p.value,
            pearson.test(y)$p.value, anscombe.test(y)$p.value,
            bonett.test(y)$p.value, ks.test(yy, "pnorm", mean(yy), sd(yy))$p.value,
            lillie.test(y)$p.value)

       colnames(B) <- "p-value"
       rownames(B) <- c("Robust Jarque-Berra test", "Anderson-Darling test", 
                     "Cramer-von Mises test", 
                     "Pearson test", "Anscombe-Glynn test","Bonett-Seier test", 
                     "Kolmogorov-Smirnov test", "Lilifors Kolmogorov-Smirnov test") 
       return(list(statistics=A, goodness.of.fit.tests.for.normality=B)) }
   }
 
   else

   A   # return just the stats

}  # end of function definition



# demo the function

# library(ISwR); data(IgM)

# new.stats(IgM)                           # minimal call
# new.stats(IgM, ylab="IgM in children")
# new.stats(IgM, digits=4, ylab="IgM in children")
# new.stats(IgM, digits=2, norm.test=TRUE)

# new.stats( rnorm(1e5) )
# new.stats( rnorm(1e5), digits=5, norm.test=TRUE )

#
# age at walking data: a=4, total n=23
#

# new.stats(active)
# new.stats(passive)
# new.stats(none)
# new.stats(control)
#
# tapply(months, treatment, new.stats)    # <---------- WORKS ----------

# new.stats(active, norm.test=TRUE)  # N too small for norm tests; must be N > 7






CImeans.normaltheory <- function(y, x, prob=.95){

  # -----------------------------------------------------------------
  #  CIs of means for each group in a one-way layout, calculated by
  #  normal theory using SEM and t-distribution at n-1 df
  #
  #  ARGS  y,x  numeric response vector paired with x: factor vector
  #             NAs OK for x and y
  #       prob  use prob=.95 for 95% CI, prob=.99 for 99% CI, etc.
  #  
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # -----------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),            # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.numeric(prob)) 

  dfr <- na.omit( data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x

  n     <- tapply(y, x, length) # vector of sample sizes per group
  Mean  <- tapply(y, x, mean)   # vector of means        per group
  SD    <- tapply(y, x, sd)     # vector of SD           per group
  SEM   <- SD / sqrt(n)         # vector of SEM          per group

  t.value <- prob + (1-prob)/2               # convert 'prob=.95' to .975, etc.
  LB      <- Mean - qt( t.value, n-1) * SEM  # lower bound of confidence interval
  UB      <- Mean + qt( t.value, n-1) * SEM  # upper bound of confidence interval 

  total.n <- length(y); total.SEM <- sd(y)/sqrt(total.n)
  LB.total <- mean(y) - qt(t.value, total.n - 1) * total.SEM
  UB.total <- mean(y) + qt(t.value, total.n - 1) * total.SEM

  LB <- c(LB, TOTAL=LB.total); UB <- c(UB, UB.total)

  cat(paste("\n\nMeans with", 100*prob, 
            "% confidence intervals using t-distribution\n\n") )
  A <- cbind(LB, UB, Mean=c(Mean, mean(y)), N=c(n, length(y)))
  return(A) # a nicely labeled matrix

} # end of function definition


# demo the call


# sim data: a=4, total n=78
# CImeans.normaltheory(response, drug) # returns nice text to R Console
#



#
# age at walking data: a=4, total n=23
#

# CImeans.normaltheory(months, treatment)
#
#Means with 95 % confidence intervals using t-distribution
#
#               LB       UB     Mean  N
#Active   8.606488 11.64351 10.12500  6
#Passive  9.385565 13.36443 11.37500  6
#None    10.113185 13.30348 11.70833  6
#Control 11.155806 13.54419 12.35000  5
#TOTAL   10.642867 12.05278 11.34783 23


# CImeans.normaltheory(months, treatment, prob=.99)
# CImeans.normaltheory(months, treatment, prob=.90)


# --------------------------------------------------------
# for SINGLE SAMPLES it's easy to use MeanCI() in
# library(DescTools) -------------------------------------
#
#
# age at first walking data: a=4, total n=23
#
# MeanCI(active, method="classic", conf.level=.95) # using normal theory
#      mean    lwr.ci    upr.ci 
# 10.125000  8.606488 11.643512 
# 
# bootstrap CI of mean, library(DescTools)  # very fast bootstrap and bca
# MeanCI(active, method="boot", conf.level=.95, R=1e4, type="perc") # percentile bootstrap
# MeanCI(active, method="boot", conf.level=.95, R=1e4, type="bca")  # BCA bootstrap 





EDA.graphs.one.sample <- function (y=rnorm(sample(20:100, size=1)), data.name=NULL){

  # -----------------------------------  ------------------------------
  #  EDA for a single, numeric sample
  #
  #  ARGS  y,x        numeric vector, a single sample of data, NAs OK
  #        data.name  decriptive string, e.g., "Sim data"
  #
  #  RETURNS window of 4 intensely annotated graphs; it can be 
  #  called with no args. Inspired by eda.uni() in package, QuantPsyc
  #
  #  EDA: Exploratory Data Analysis. Methods developed by 
  #  John Tukey et al. for data exploration & visualization. 
  #  Today, the analyst does EDA to probe and to understand the data.
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # -------------------------------------------  ----------------------

  if( missing(y) ) data.name <- "Simulated random data, N(0,1)"
  stopifnot(is.numeric(y))
  y <- y[!is.na(y)]  # remove NA, missing values
  require(DescTools, quietly=TRUE)  # skew, kurtosis & robust JarqueBeraTest()

  par( mfrow = c(2, 2) )  # 2 rows by 2 columns of graphs

    hist( y, 
        main = "Histogram", cex.main=1.2, border="white",
        col="lightblue", xlab = data.name, cex.axis=.8,
        sub=paste("skew=",     round(DescTools::Skew(y,method=3),3), 
                "  kurtosis=", round(DescTools::Kurt(y,method=3),3)), 
                cex.sub=1, col.sub="blue" )
 
    if(length(y) <= 1e3) rug(y, col="gray") 

    mtext(paste("N=", length(y), " Mean=", round(mean(y),4), 
                " SD=", round(sd(y),4)), cex=.6, col="blue")

    plot( density(y, n=1024), main = "Kernel Density Estimate", 
        cex.main=1.2, cex.axis=.8, col="red", lwd=3 )    
    
    boxplot(y, horizontal = TRUE,
        col="wheat", xlab=data.name, main= "Boxplot", ylab="Response", boxwex=.9,
        sub=paste("IQR=", round(IQR(y),3),
                  " min=", round(min(y),5), " max=", round(max(y),5)), 
        cex.sub=.8, col.sub="blue", cex.axis=.8 )

    mtext(paste("hinges: 25th=", round(quantile(y,.25),3),  
                " 75th=", round(quantile(y,.75),3)), cex=.65, col="blue")

    a <- boxplot(y, plot=FALSE); outliers <- length(a$out)    
    mtext( paste(" N of outliers=", outliers, "  >1.5*IQR from hinge"), 
        side=3, line=-1, adj=0, cex=.6, col="blue")

    mtext(paste("median=", round(median(y),4) ), 
        side=1, line=-1, cex=.7, col="blue")

    a <- JarqueBeraTest(y, robust=T)$p.value  # p-value of normality test
    if(a == 0) a <- "<2.2e-16"      # if it's zero then make it a string
   
    # ------------------------------ Normal QQ plot to visualize normality
    if(is.character(a)){
        qqnorm(y, col="blue", sub=paste("Robust Jarque-Bera test p", a), 
        ylab="Observed Quantiles", 
        cex.sub=.8, col.sub="blue", cex.axis=.8 )}
    else

   {qqnorm(y, col="blue", sub=paste("Robust Jarque-Bera test p=", 
        signif(a,5)), ylab="Observed Quantiles", 
        cex.sub=.8, col.sub="blue", cex.axis=.8 )}
    # ------------------------------

    if(length(y) <  5000) mtext(paste("Shapiro-Wilk test p=", 
        signif(shapiro.test(y)$p.value,5)), cex=.7, col="blue")

    if(length(y) >= 5000) mtext(paste("Lillifors Kolmogorov-Smirnov test p=", 
        signif(LillieTest(y)$p.value,5)), cex=.7, col="blue") 

    qqline(y, col="blue")  # add straight line to Normal QQ plot 
    
    mtext(" Visualization of observed distribution\n relative to a normal distribution",
        side=3, line=-1.5, adj=0, cex=.6, col="lightblue4")

  par(mfrow = c(1, 1))  # reset to default of 1 graph per page

} # end of function definition


# demo the function

# EDA.graphs.one.sample()  # call function with no arguments, to see what it does

# EDA.graphs.one.sample(rnorm(1e2, 50, 10), data.name="Sim data from N(50,10)")

# library(ISwR); data(IgM)
# EDA.graphs.one.sample(IgM, data.name="IgM, g/l")  




EDA.plots.oneway <- function(y, x){

  # ------------------------------------------------------------------
  #  EDA graphs per group, in a one-way layout
  #
  #  ARGS  y,x  numeric response vector, paired with x: factor vector
  #             NAs OK
  #
  #  Needs:   my.pause()  thin.line()   EDA.graphs.one.sample()  
  #
  #  function simply calls EDA.graphs.one.sample() in a 'for' loop
  # 
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),             # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x)) 

  dfr <- na.omit(data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x          # reconstitute y & x

  cat("\n"); thin.line(60)
  cat("To inform the ANALYST,\nview EDA graphs for response at each level of factor x\n")
  cat("and for the entire response vector and\nfor the residuals from lm(y~x)\n")
  thin.line(60)

  a <- levels(x)  # the group names

  for(i in 1:length(a)){my.pause(); EDA.graphs.one.sample(y[x == a[i]], data.name=a[i])}

  my.pause()
  EDA.graphs.one.sample(y, data.name="-- the merged data --") # the ENTIRE response vector
   
  my.pause()
  EDA.graphs.one.sample(residuals(lm(y ~ x)), 
     "Distribution of lm(y~x) Residuals") # distribution of data set residuals
   
  my.pause()
  boxplot(y ~ x, col=gray.colors(length(a)), ylab="Response", 
    boxwex=.5, cex.lab=1.2, cex.axis=1.2)
  cat("end. EDA.plots.oneway()\n"); thin.line(60)

} # end of function definition


# demo the function 

#
# age at first walking data: a=4, total n=23
#
# EDA.plots.oneway(months, treatment)




NormalQQ.oneway <- function(y, x, data.name="my data"){

  # --------------------------------------------------------------------------------
  #  Normal QQ reference plots per level of x using qreference() from library(DAAG)
  #
  #   ARGS    y   numeric response vector, NAs OK
  #           x   corresponding factor vector, NAs OK
  #   data.name   descriptive string, e.g., "Age at Walking"
  #
  #  Superb, visual assessment of normality thanks to John Maindonald
  # --------------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),            # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.character(data.name)) 

  dfr <- na.omit( data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x

  require(DAAG) # by John Maindonald
  cat("\n\n"); thin.line(63)
  cat("Visual assessment of normality by reference to data simulated\n")
  cat("from normal distribution,\n")
  cat("using normal QQ plots in calls to qreference() in library(DAAG)\n")
  cat("by John Maindonald.\n"); thin.line(63)

  # reference, normal QQ plots per group -- per level of factor x
  a <- levels(x)  # levels of factor x

  for(i in 1:length(a)) {    # loop through the levels of factor x
    plot.new()  # create graph window

    qreference( y[x == a[i]], cex.strip=1,
      xlab=paste("Group: ", a[i], "\nfrom: ", data.name ), ylab="" )

      title(main="Assessment of normality by qreference() in DAAG\nby John Maindonald",
            cex.main=1.1, col.main="red", 
        sub="blue: observed data,  y-axis: response,  x-axis: theoretical quantiles", 
             col.sub="blue", cex.sub=.8, font.sub=2 )
    my.pause()
  }

  qreference(y, cex.strip=1,
     xlab="-- all the data --", ylab="")
     title(main="Assessment of normality by qreference() in DAAG\nby John Maindonald",
     cex.main=1.1, col.main="red", 
     sub="blue: residuals,  x-axis: theoretical quantiles", 
     col.sub="blue", cex.sub=.8, font.sub=2 )

  my.pause()

  qreference( residuals(lm(y ~ x)), cex.strip=1,
     xlab="RESIDUALS from lm(y ~ x)", ylab="")
     title(main="Assessment of normality by qreference() in DAAG\nby John Maindonald",
     cex.main=1.1, col.main="red", 
     sub="blue: residuals,  x-axis: theoretical quantiles", 
     col.sub="blue", cex.sub=.8, font.sub=2 )

  cat("\nend. NormalQQ.oneway()\n"); thin.line(60)

} # end fun definition
 

# demo the function 

#
# age at first walking data: a=4, total n=23
#

# NormalQQ.oneway(months, treatment, data.name="Age at first walking")



Normal.sim.kernel.band <- function(y, NS=500){

  # -----------------------------------------------------------------
  #  Univariate visual assessment of normality by kernel density
  #  band from many data sets simulated from normal distribution 
  #  of the same n, mean and sd as observed data, y. 
  #  Returns annotated graph.
  #
  #  ARGS     y   numeric vector, missing values OK  
  #          NS   number of sim data sets from normal distribution
  #               each set generating a kernel density curve
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # -----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1,            # arg & data checking
    is.numeric(NS) ) 

  y <- y[!is.na(y)]     # remove missing values from numeric vector
  n <- length(y); m <- mean(y); s <- sd(y) # stats

  max.y <- 0  # estimate top of kernel band for sim data -- 'max.y'
  for(i in 1:NS/2){
    a <- density(rnorm(n, m, s))
    if(max(a$y) > max.y) max.y <- max(a$y)
  }
  max.y <- 1.2 * max.y

  hist(y, freq=FALSE, ylim=range(c(0, density(y)$y, max.y)),   
    col="slategray4", xlab="Response", main="", border="white", 
    cex.axis=1.2, cex.lab=1.2, font.lab=2
  )

  mtext(paste("Kernel density estimate for each of NS =", NS,
    "random samples from normal distribution"), line=2, cex=.8) 
  mtext("of same N, Mean and SD as observed sample -- visual assessment of normality", 
    line=1, cex=.8)

  for(i in 1:NS){
    a <- density(rnorm(n, m, s), from=min(y), to=max(y), n=512)
    lines(a, col="gray", lwd=1.5)
  }

  lines( density(y, from=min(y), to=max(y)), lwd=2, col=2)  # add kernel dens curve

  mtext("GRAY: kernel density for\nsim data from N( )", col="gray",
    side=1, line=3, adj=0, cex=.65, font=2)
  mtext("RED: kernel density for observed data", col="red",
    side=1, line=4, adj=0, cex=.65, font=2)
  mtext("If observed kernel density is well within", side=1, line=2, adj=1, cex=.65)
  mtext("the band of sim kernel density", side=1, adj=1, line=3, cex=.65)
  mtext("this is evidence of normality", side=1, adj=1, line=4, cex=.65)

} # end function definition


# demo the function

# Normal.sim.kernel.band(y=rnorm(100))

# how to extract single groups, levels, from a data frame
# y <- df[ df$x == "Control", ]; Normal.sim.kernel.band(y) 



BayesFactorAnova.oneway <- function(y, x, NS=1e5, fixed.effects=TRUE){

  # -----------------------------------------------------------------
  #  BAYES FACTOR ANOVA   for fixed treatment models
  #                       and for random effects, as specified
  #                       using anovaBF() in library(BayesFactor)
  #
  #   ARGS    y       numeric response vector, NAs OK
  #           x       corresponding factor vector, NAs OK
  #          NS       number of Monte Carlo simulations
  #    fixed.effects  logical; fixed effects or random effects model
  #
  #  See  BFManual()  to open a detailed manual.
  #  lmBF() does similar analysis.
  #
  #  This is merely a call to anovaBF() with a text header.
  # -----------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), 
    is.numeric(NS), is.logical(fixed.effects) ) 

  dfr <- na.omit( data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x

  require(BayesFactor)
  cat("\n"); thin.line(65)
  cat("BAYES FACTOR ANOVA\nusing library(BayesFactor)\n")
  cat("\tanovaBF() is for fixed treatment or random effects Anova\n\n")
  cat("A Bayes factor is a multiplicative factor of support for H1 vs. Ho\n")
  cat("Documentation: > help(package=BayesFactor)    >?anovaBF \n")
  thin.line(65); cat("\n")

  cat(paste("You specified an Anova model that's 'fixed effects' =", fixed.effects), "\n\n")

  if(fixed.effects == TRUE) {anovaBF(y ~ x, data=dfr, progress=F, whichRandom=NULL, iterations=NS)}
     else                     
  {anovaBF(y ~ x, data=dfr, whichRandom="x", progress=F, iterations=NS)}

} # end of function definition


# demo the function

#
# age at first walking data: a=4, total n=23
#
# BayesFactorAnova.oneway(months, treatment)  # minimal call, 'fixed effects' is default
# BayesFactorAnova.oneway(months, treatment, fixed.effects=TRUE)
# BayesFactorAnova.oneway(months, treatment, fixed.effects=FALSE)  
# BayesFactorAnova.oneway(DF$months, DF$treatment)

#
# call anovaBF() directly -------------------------------
#
# not for Windows OS but works on multicore Mac OS, >?anovaBF

# anovaBF(months ~ treatment, data=DF, whichRandom="treatment", iterations=1e4, multicore=FALSE)
# treatment:  0.552 +- 0%   # random effects model

# anovaBF(months ~ treatment, data=DF, whichRandom=NULL, whichModels="all", iterations=1e4, multicore=FALSE)
# treatment:  0.939 +- 0%   # fixed effects model


#
# call lmBF() directly. It handles Anova and regression.  >?lmBF
#

# lmBF(months ~ treatment, data=DF, whichRandom="treatment", iterations=1e5)
# treatment:  0.552 +- 0%

# double the data as a very crude form of power analysis
# lmBF(months ~ treatment, data=rbind(DF, DF), whichRandom="treatment", iterations=1e5)
# treatment:  5.499 +- 0%




treeperm.oneway.anova <- function(y, x, NS, type="approximate", graph.F=FALSE){

  # ----------------------------------------------------------------------------
  #  A DIFFERENT APPROACH to a K means,                       USE WITH CAUTION
  #    ONE-WAY PERMUTATION TEST for location                  USE WITH CAUTION
  #                                                           USE WITH CAUTION
  #  library(treeperm); a new library on 4/22/15 by           USE WITH CAUTION
  #  Professor Qiao Kang, St. Andrews College, UK             USE WITH CAUTION
  #                                                           USE WITH CAUTION
  #  Runs compiled C code; uses tree algorithms. Supports       unless total n
  #  EXACT and APPROXIMATE (asymptotic) tests of location.        is low
  #
  #   ARGS    y       numeric response vector, NAs OK         USE WITH CAUTION
  #           x       corresponding factor vector, NAs OK        unless total n
  #          NS       number of permutations desired               is low
  #        type       "approximate"  or "exact"
  #     graph.F       logical
  #          
  #  treeperm() in library(treeperm)  Documentation: > help(package=treeperm)
  #  My function is merely a wrapper with a timer that enforces 'approximate' 
  #  permutation; override it by 'exact' for small N but at your own peril; 
  #  it doesn't have a 'time out' feature. Extreme caution in using 
  #  type="exact" unless total N is low.
  #
  #  Null hypothesis is the same as in classical Anova. Indeed, here the
  #  permutation test statistic is the observed F from ordinary Anova.
  #                                                                   
  #  Reference:  Ernst, M.D. 2004. Permutation methods: A basis for   
  #  exact inference. Statistical Science 19(4):676-685.              
  # ----------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),            # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.numeric(NS), 
    is.character(type), is.logical(graph.F) ) 

  dfr <- na.omit( data.frame(y, x))      # delete rows with NA
  y  <- as.numeric(dfr$y) ; x <- dfr$x   # reconstitute y & x

  require(treeperm)

  cat("\n"); thin.line(64)
  cat("PERMUTATION one-way Anova using library(treeperm)\n")
  cat("treeperm() calls compiled C code and is extremely fast.\n")
  cat("Test statistic: observed F from one-way Anova on observed data. \n")
  cat("'approximate' permutation is wisely enforced here but you may\n")
  cat("change fun definition to allow sending 'exact' if N is small,\n")
  cat("but at your own peril; see Kincaid's demo code.\n\n")

  cat("Description: An implementation of permutation tests in R, supporting both\n") 
  cat("exact and asymptotic K sample test of data locations. The p value of exact\n")
  cat("tests is found using tree algorithms. Tree algorithms treat permutations of\n")
  cat("input data as tree nodes and perform constraint depth-first searches for\n") 
  cat("permutations that fall into the critical region of a test systematically.\n") 
  cat("Pruning of tree search and optimisations at C level enable exact tests for\n")
  cat("certain large data sets. (quoted from treeperm package by Qiao Kang) \n")
  thin.line(64); cat("\n"); flush.console()

  if(type == "exact") type="approximate"   # as a safe guard; comment-out as desired

  # treeperm() may interpret numeric 'y' as integer, if it's integers
  y <- as.numeric(y)  # in case response is 'integer' which treeperm() rejects
  flush.console()  # R Console buffering -- dump text

  treeperm.time <- proc.time()
     perm <- treeperm(y ~ x, type=type, size=NS)  # the heavy lifting
  treeperm.time <- proc.time() - treeperm.time
  treeperm.minutes <- treeperm.time/60  # convert seconds to minutes

  cat("Minutes elapsed: \n")
  print(treeperm.minutes[[3]] )  # elapsed minutes
  cat("\n"); thin.line(5); cat("\n")
  cat("Text output from the call to treeperm()\n\n")
  print(perm)

  if(graph.F == TRUE){  # a quirk: plot.treeperm() sends obs F as text to Console
     if(NS >  1e5){ NS.2 <- 1e5; plot(perm, size=NS.2) }
     if(NS <= 1e5) plot(perm, size=NS)
  } 

} # end of function definition


# demo the function

# sim data: a=4, total n=78

# treeperm.oneway.anova(response, drug, NS=1e4, type="approximate")
# treeperm.oneway.anova(response, drug, NS=1e4, graph.F=TRUE)
# treeperm.oneway.anova(response, drug, NS=1e5)
# treeperm.oneway.anova(response, drug, NS=1e2)

#
# age at first walking data: a=4, total n=23
#
# treeperm.oneway.anova(months, treatment, NS=1e5, graph.F=FALSE)
# treeperm.oneway.anova(months, treatment, NS=1e5, graph.F=TRUE)

# ddf <- rbind(DF, DF)  # double the data
# treeperm.oneway.anova(ddf$months, ddf$treatment, NS=1e4, graph.F=FALSE)


# call treeperm() directly ---------------------------------
#
# NS <- 1e5
# permtime <- proc.time()
#   out <- treeperm(months ~ treatment, frame=DF, type="approximate", size=NS)
# permtime <- proc.time() - permtime
# permtime <- permtime/60  # minutes elapsed
# out
# permtime
# plot(out, size=NS)


# ------------------------------------------------------------------------------
# try an 'exact' permutation on a small data set, Age at First Walking a=4, 
# total n=23, by calling treeperm() directly
#
# CAUTION  CAUTION  CAUTION   if total N is > 25 it may FREEZE YOUR MACHINE
#
# perm Anova using "exact"  -- EXACT PERMUTATION TEST for location
#
#
# timer <- proc.time()
#   treeperm(DF$months ~ DF$treatment, DF, type="exact")
# timer <- proc.time() - timer
# timer <- timer/60  
# timer  # minutes

# Total permutations
# [1] 577185873264       # 577 billion permutations required 
# F statistics
# [1] 2.142222           # F stat from traditional Anova
# P value
# [1] 0.1301595          # two-sided, exact 
#
#     user   system  elapsed   # 4 hrs on an old ASUS laptop
# 243.8842   0.0115 244.3655   # Intel Celeron CPU, 4 gb RAM



# --------------------------------------------------------
# NOTE: library(treeperm) has other permutation test funs BUT...
# Permutation( ) does unpaired 2 group permutation test for difference
# in means, two-sided BUT USE IT ONLY for SMALL N, else it locks up
# computer...

# df <- read.csv("Framingham period 1.csv")
# grab only the first 500 subjects, rows
# bmi <- df$BMI[1:500]; sex <- as.factor(df$SEX[1:500])

# dfr <- na.omit(data.frame(bmi, sex))  #n=498  NAs an issue here
# str(dfr)

# t.test(bmi ~ sex, dfr, var.equal=TRUE) # p=0.09292 (2)
# treeperm(dfr$bmi ~ dfr$sex, dfr, type="approximate", size=1e4)
# p=0.09455 (2)

# Kincaid.perm.oneway(dfr$bmi, dfr$sex, NS=1e4) # p=.09499, 1.4 min

# timer <- proc.time()
#   treeperm(dfr$bmi ~ dfr$sex, dfr, type="approximate", size=1e4)
# timer <- proc.time() - timer
# timer <- timer/60 # min
# timer  # 0.19 minutes, p=0.0957 after 1e4 permutatons





Kincaid.perm.oneway <- function(y, x, NS=1e2, describe.null.F=FALSE, returnF=FALSE){

  # ----------------------------------------------------------------------------
  #  Permutation test for K means in a one-way layout. 
  #  Null hypothesis same as in one-way Anova; F is test statistic. 
  #
  #  ARGS     y          numeric response vector, NAs OK
  #           x          corresponding factor vector, NAs OK
  #          NS          number of desired, global shuffles of y
  #                      start low, note time required, go high
  #    describe.null.F   logical
  #        returnF       logical
  #
  #  Results agree well with treeperm() in package treeperm, which in some cases 
  #  is faster but with fewer diagnostics, but if library(treeperm) is abandoned
  #  or deprecated as has happened for some other permutation libraries in the  
  #  past, then it's best to have backup, home brew code using only base R. Also,
  #  treeperm( ) cannot 'handle' large N but my function will, albeit slower.
  # 
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),            # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.numeric(NS), 
    is.logical(describe.null.F), is.logical(returnF) ) 

  dfr <- na.omit( data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x

  cat("\n"); thin.line(78)
  cat("Kincaid.perm.oneway(y, x, NS, describe.null.F, returnF); for a one-way layout\n")
  cat("Test statistic: observed F from one-way Anova on observed data. \n")
  cat("Null hypothesis is the same as in classical one-way Anova.\n")
  cat("Start low on NS (e.g., 1e3), note the time required then go high.\n")
  thin.line(78); cat("\n"); flush.console()

  perm.time <- proc.time()     # begin timer
    #perm.F  <- replicate(NS, anova(lm(sample(y) ~ x, data=dfr))[1,4])  # iterator for F
    perm.F <- replicate(NS, oneway.test(sample(y) ~ x, data=dfr, var.equal=TRUE)$statistic)
  perm.time <- proc.time() - perm.time
  perm.sec  <- perm.time[[3]]  # pull-out elapsed time in seconds
  perm.min  <- round( perm.sec/60, 4) # convert to minutes and round-off
   
  obs.F  <- anova(lm(y ~ x))[1,4]
  NGE    <- length(perm.F[perm.F >= obs.F]) # how many null F-stats are >= obs.F
  perm.p <- (NGE+1)/(NS+1)  # standard formula for p in any permutation test

  cat("TEST STATISTIC:\n\nobserved F =", obs.F, "in classical one-way Anova.\n\n\n")
  cat("RESULTS:\n\n")
  cat("Minutes elapsed =", perm.min, "\n")
  cat("NGE =", NGE, "\n")
  cat("p =", perm.p, "achieved after NS =", NS, "permutations.\n\n")

  if(describe.null.F == TRUE){
    thin.line(20); cat("\n")

    cat("The NULL DISTRIBUTION of F after NS shuffles of the response\n")
    print(new.stats(perm.F, norm.test=F))
    cat("\nQUANTILES of the empirical, null F distribution under Ho\n")
    print(quantile(perm.F, c(.9, .95, .975, .99, .999, .9999)))
    cat("\n"); thin.line(20); cat("\n")

    if(NS <= 1000) bins <- "sturges" else bins <- "scott"

    hist(perm.F, freq=TRUE,
       breaks=bins,    # comment-out or change as needed
       col="wheat", xlab="F in permutation one-way Anova", 
       cex.lab=1.5, font.lab=2, font.main=4, cex.main=1.5, border="white",
       main=paste("Distribution of null F after", NS, "shuffles")
     )

    abline(v=obs.F, col="black", lwd=2)  # vertical line at observed F-stat
      mtext(paste("permutation p =", signif(perm.p, 5),
         ", observed F =", signif(obs.F, 5), "in classical one-way Anova"),
         cex=.8, col="black" )
      mtext("black line at observed F-statistic", side=1, line=4, cex=.8, adj=1) 
      mtext( paste("elapsed minutes:", perm.min), side=4, cex=.8 )
      mtext(paste("NGE =", NGE), side=3, line=-1, cex=.9)
  }
  if(returnF == TRUE) return(perm.F)  # the entire distribution of null F

} # end function definition


# demo the function

#
# age at first walking data: a=4, total n=23
#
# Kincaid.perm.oneway(months, treatment, NS=1e3, describe.null.F=TRUE)
# Kincaid.perm.oneway(months, treatment, NS=1e4, describe.null.F=TRUE)

# if you need or want the NULL distribution of F to be returned but because
# it's so huge, send it to an object and not to R Console screen!
#
# a <- Kincaid.perm.oneway(months, treatment, NS=1e4, returnF=TRUE)
# hist(a, breaks="scott"); summary(a)





powercurve.normaltheory.oneway <- function(y, x, low.n, high.n, by.n, matrix.out=TRUE){

  # --------------------------------------------------------------------------------
  #  Traditional power analysis in one-way Anova
  #
  #  ARGS      y                  numeric response vector, NAs OK
  #            x                  corresponding factor vector, NAs OK
  #         low.n, high.n, by.n   desired sample sizes (df2 in Anova table)         
  #      matrix.out               logical; table of power by alpha and by N
  #
  #  Works for BALANCED and UNBALANCED one-way Anova designs; PV is explained
  #  percent (enter proportion) of variance, usually taken as variance
  #  component but here, R-square is apparently what it expects, from Anova table.
  #
  #  library(QuantPsyc)   function: powerF(PV, df2, df1, alpha)
  #  powerF() DEALS with TOTAL N (i.e., reflected by 'df2' from Anova table
  #  >?powerF   for details & reference: Murphy, K.R. & Myors, B. 2004.
  #  "Statistical power analysis: A simple and general model for traditional 
  #  and modern hypothesis tests" 2nd ed. Mahwah, NJ: Lawrence Erlbaum Associates
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # --------------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),            # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.logical(matrix.out), 
    is.numeric(low.n), is.numeric(high.n), is.numeric(by.n) ) 

  dfr <- na.omit(data.frame(y, x))   # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x

  require(QuantPsyc, quietly=TRUE)

  vc <- R.sq.oneway(y, x)[1]/100  # apparently it wants explained R-sq as proportion

  a   <- length(levels(x))     # number of groups
  v1  <- anova(lm(y ~ x))[1,1] # v1 from anova table
  v2  <- anova(lm(y ~ x))[2,1] # v2 from anova table
  N   <- v2 + a                # total N

  n <- seq(low.n, high.n, by.n) # values of df2 for plotting on x-axis

  obs.pwr <- powerF(PV=vc, df2=v2, df1=v1, alpha=.05)  # power of the test at alpha=.05
  out05   <- powerF(PV=vc, df2=n,  df1=v1, alpha=.05)  # vectorized
  out01   <- powerF(PV=vc, df2=n,  df1=v1, alpha=.01)

  plot(n, out05, pch=16, type="b", 
    main="Power curves\nusing normal theory after Murphy & Myors (2004)",
    cex.main=1, ylim=c(0,1), font.lab=2, cex.lab=1.2,
    xlab="df2 as measure of total N",  # xlab="Total N",
    ylab="Statistical power"
  )

  abline(h=c(seq(0,.7,.1), .9, 1), lty=2, col="gray")
  abline(h=.8, lty=2, lwd=2, col="black")
  abline(v=N, lwd=2, col="red")
  points(n, out05, pch=16, type="b", col="black") # overplot points, as needed

  points(n, out01, pch=1, type="b", cex=.5, col="blue" )
  mtext("red: observed total N", side=1, line=4, adj=0, col="red")
  mtext("blue: power at alpha=.01", side=1, line=4, adj=1, col="blue")
  mtext( paste("power of the test =", round(obs.pwr, 4), "at 0.05"), side=4)

  if(matrix.out == TRUE){
    power.matrix <- cbind(n+a, n, round(out05, 3), round(out01, 3), round(vc, 3))
    colnames(power.matrix) <- c("total N", "df2", "power.05","power.01", "Explained R-sq")
    return(power.matrix)
  } 

} # end function definition


# demo the function


#
# age at walking data: a=4, total n=23
#
# powercurve.normaltheory.oneway(months, treatment, low.n=5, high.n=70, by.n=2)

# powercurve.normaltheory.oneway(months, treatment, low.n=5, high.n=70, by.n=1, matrix.out=TRUE)

# USING powerF() DIRECTLY, from library(QuantPsyc)
# R.sq = .2528, power of the test = .46     ACTUAL DATA NOT NEEDED!
# powerF(PV=.2528, df2=19, df1=3, alpha=.05)  # returns power ; library(QuantPsyc)
# powerF(PV=.2528, df2=46, df1=3, alpha=.05)  # returns power ; library(QuantPsyc)


# USING pwr.anova.test() DIRECTLY, from library(pwr)
# balanced design only, supposedly;  ACTUAL DATA NOT NEEDED!   library(pwr)
# Cohen's D = .5286, power of the test = .46   ASSUMES BALANCED DESIGN!
# obs.pwr <- pwr.anova.test(k=4, n=5.75, f=.5286, sig.level=.05, power=NULL)  # uses Cohen's D
# obs.pwr

# future.N <- pwr.anova.test(k=4, n=NULL, f=.5286, sig.level=.05, power=.8) # Cohen's D
# future.N  # 11 per group

# USING pwr.f2.test()  DIRECTLY, from library(pwr); f2 is Rsq
# library(pwr);  ACTUAL DATA NOT NEEDED!  obs power = .42, if R-square used
# obs.pwr <- pwr.f2.test(u=3, v=19, f2=.2528, sig.level=.05, power=NULL) 
# obs.pwr

# future.N <- pwr.f2.test(u=3, v=NULL, f2=.2528, sig.level=.05, power=.8)
# future.N  # v2 = 43.2234,  so total N = 44

# USING power.anova.test()  DIRECTLY, from library(stats) in base R
# this works! power.anova.test() in base R, and agrees with pwr, QuantPsyc
# m <- tapply(months, treatment, mean)
# varr <- var(m)
# power.anova.test(groups=4, n=23/4, between.var=varr,
#   within.var=2.2995, sig.level=.05, power=NULL) 
# power.anova.test(groups=4, n=NULL, between.var=varr,
#   within.var=2.2995, sig.level=.05, power=.8) 




powercurve.Cohen.oneway <- function(y, x, data.name=NULL, low.n, high.n, by.n, 
   matrix.out=TRUE){

  # -------------------------------------------------------------------------
  #  Traditional power analysis in one-way Anova
  #
  #  ARGS      y                  numeric response vector, NAs OK
  #            x                  corresponding factor vector, NAs OK
  #        data.name              character
  #         low.n, high.n, by.n   desired sample sizes (N per group)         
  #      matrix.out               logical; table of power by alpha and N
  #
  #  using pwr.anova.test( )      after Cohen (1988)
  #  in library(pwr)
  #
  #  pwr.anova.test() is intended for BALANCED one-way Anova designs; but 
  #  unbalance in N per group is probably OK; power analysis is estimation.
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # -------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),            # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.logical(matrix.out), 
    is.numeric(low.n), is.numeric(high.n), is.numeric(by.n) ) 

  dfr <- na.omit( data.frame(y, x))  # delete rows with NA
  y  <- dfr$y ; x <- dfr$x           # reconstitute y & x

  require(pwr)

  D <- CohenD.oneway(y, x, verbose=FALSE)  # Cohen D as effect size

  a   <- length(levels(x))     # number of groups
  v1  <- anova(lm(y ~ x))[1,1] # v1 from anova table
  v2  <- anova(lm(y ~ x))[2,1] # v2 from anova table
  N   <- v2 + a                # total N

  pwr.of.test <- pwr.anova.test(k=a, n=N/a, f=D, sig.level=.05)$power  # power of test
  N.at.8      <- pwr.anova.test(k=a, f=D, sig.level=.05, power=.8)$n   # N per group
  N.at.9      <- pwr.anova.test(k=a, f=D, sig.level=.05, power=.9)$n   # N per group
   
  n.per.group <- seq(low.n, high.n, by.n) # vector of prospective N per group

  out05  <- pwr.anova.test(k=a, n=n.per.group, f=D, sig.level=.05 ) # vectorized,
  out01  <- pwr.anova.test(k=a, n=n.per.group, f=D, sig.level=.01 ) # returns 
  out001 <- pwr.anova.test(k=a, n=n.per.group, f=D, sig.level=.001) # power

  plot(n.per.group, out05$power, pch=16, type="b", 
     main="Power curve at 0.05 alpha\nusing normal theory after Cohen (1988)",
     cex.main=1, ylim=c(0,1), font.lab=2, cex.lab=1.2,
     xlab=paste("Sample size per group (", data.name, ")", sep=""), 
     ylab="Statistical power",
     sub="red: observed N per group;  black: N for power = .8, .9"
  )

  abline(h=c(seq(0,.7,.1), .9, 1), lty=2, col="gray")
  abline(h=.8, lty=2, lwd=2, col="black")
  arrows(N.at.8, .8, N.at.8, 0, code=2, lty=1, lwd=3, col="black")  # pwr=.8 to x-axis
  arrows(N.at.9, .9, N.at.9, 0, code=2, lty=1, lwd=3, col="black")  # pwr=.9 to x-axis
  arrows(N/a, pwr.of.test, N/a, 0, code=1, lwd=3, lty=1, col="red") 
  points(n.per.group, out05$power, pch=16, type="b", col="black") # to overplot points

  my.pause()
  points(n.per.group, out01$power, pch=1, type="b", cex=.7,col="blue") # plot .01 curve
  mtext("blue: 0.01", side=4, col="blue")

  my.pause()
  points(n.per.group, out001$power, pch=1, type="b", cex=.7,col="gray") # plot .001 curve
  mtext("gray: 0.001", side=4, adj=0, col="gray")
  my.pause()
  abline(v=n.per.group, lty=2, col="gray")

  if(matrix.out == TRUE){
     total.N <- a * n.per.group
     power.matrix <- cbind(total.N, n.per.group, 
                           round(out05$power,4), round(out01$power,4), 
                           round(out001$power,4), round(D,3))
     colnames(power.matrix) <- c("total N", "N per group", 
                                 "power.05", "power.01", "power.001", "Cohen's D")
     return(power.matrix)
  }

} # end of function definition


# demo the function


#
# age at first walking data: a=4, total n=23
#
# y <- months;  x <- treatment   # for convenience
# powercurve.Cohen.oneway(y, x, data.name="Age at Walking", low.n=5, high.n=30, by.n=1)




boot.power.of.test.oneway <- function(y, x, data.name=NULL,
   NS=1e3, N.boot=0, n.breaks=100, return.distribution=FALSE){

  # ------------------------------------------------------------------------
  #  Bootstrap power of the test in one-way Anova       
  #
  #  returns power of the test (one-way Anova) at various alpha, plus graph
  #
  #  Note: Usually, the POWER OF THE TEST is much less interesting 
  #        than PROSPECTIVE power analysis.
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ------------------------------------------------------------------------

  # y: numeric response vector, paired with x: factor vector
  # NS: number of bootstrap samples to get; the entire data frame is bootstrapped
  x   <- as.factor(x)      # just to make sure
  dfr <- data.frame(y, x)

  boot.power.one.way <- function(dframe, N.boot=length(dframe[,1]) ){ 
    N  <- length(dframe[,1])
    bf <- dframe[sample(1:N, N.boot, replace=TRUE), ] # bootstrap sample of data frame
    out <- anova(lm(bf[,1] ~ bf[,2]))
    boot.p <- out[1,5]  # anova p-value is returned for bootstrapped data frame
  }

   N <- length(dfr$y)           # number of observed rows
   if(N.boot < 1) N.boot <- N   # N of desired bootstrap rows, usually all rows

   boot.time <- proc.time()
      a <- replicate(NS, boot.power.one.way(dfr))  # iterator
   boot.time <- proc.time() - boot.time
   boot.time <- boot.time/60  # time in minutes

   pwr.10 <- pwr.05 <- pwr.02 <- pwr.01 <- 0  # initialize

   pwr.10 <- length( a[ a <= .10  ] ) / NS  # delivers power as 
   pwr.05 <- length( a[ a <= .05  ] ) / NS  # proportion of boot-
   pwr.02 <- length( a[ a <= .02  ] ) / NS  # strap samples that
   pwr.01 <- length( a[ a <= .01  ] ) / NS  # had p-value <= alpha

   ylim.max <- max(hist(a, n.breaks, plot=FALSE)$counts) * 1.05 # for hist()

   hist(a, 
     breaks=n.breaks, col="wheat", cex.lab=1.2, 
     main="Bootstrap power of the test (one-way anova)",
     cex.main=1.2, font.main=4, col.main="black",
     sub=paste( round(boot.time[[3]],2), "minutes for resampling"), cex.sub=.8,
     xlab=paste("p-value in one-way Anova (", data.name, ")", sep=""), 
     ylab="Frequency",
     cex.lab=1.2, font.lab=2, cex.axis=1,
     ylim=c(0, ylim.max)
   )

   text(max(a)/2, ylim.max/1.5, paste("power at alpha .10 =", pwr.10),   pos=4, cex=.7)
   text(max(a)/2.7, ylim.max/2, 
      paste("power of the test at alpha .05 =", pwr.05), pos=4, cex=.8, col="red",font=2)
   text(max(a)/2, ylim.max/2.5, paste("power at alpha .02 =", pwr.02),   pos=4, cex=.7)
   text(max(a)/2, ylim.max/3,   paste("power at alpha .01 =", pwr.01),   pos=4, cex=.7)
   mtext(paste("achieved by NS =", NS, "bootstrap samples of N =", N.boot,
      " (observed, total N)"), cex=.8, font=2, col="red" )

   cat("\n"); thin.line(65); cat("\n")
   cat("Output from Kincaid's function: boot.power.of.test.oneway()\n\n") 
   cat("Power of the test at alpha of\n.10, .05, .02, .01 =", 
       c(pwr.10, pwr.05, pwr.02, pwr.01), "\n\n")
   cat(boot.time[[3]], "minutes to achieve NS =", NS, "bootstrapped data frames.\n\n")
   cat("See histogram of the bootstrap distribution of Anova p-values\n\n")
   thin.line(65)

   if(return.distribution == TRUE) return(a)

}  # end of function definition


# demo the function

# sim data

# boot.power.of.test.oneway(response, drug, data.name="sim data", NS=1e3, n.breaks=50)

# out <- boot.power.of.test.oneway(response, drug, data.name="sim data", NS=1e3, n.breaks=50, return.distribution=TRUE)
# psych::describe(out)
# pwr.05  <- length( out[ out <= .05  ] ) / 1e3
# pwr.05


#
# age at first walking data: a=4, total n=23
#

# boot.power.of.test.oneway(months, treatment, data.name="Age at walking", NS=50, n.breaks=50)
# boot.power.of.test.oneway(months, treatment, data.name="Age at walking", NS=1e3, n.breaks=50)




boot.powercurve.oneway <- function(y, x, NS, data.name=NULL, low.n, high.n, by.n){

  # ----------------------------------------------------------------------------
  #  The bootstrap used to vary sample size
  #  in prospective, statistical power analysis in one-way Anova
  #
  #  ARGS      y                  numeric response vector, NAs OK
  #            x                  corresponding factor vector, NAs OK
  #           NS                  desired N of bootstrap samples per sample size
  #        data.name              character
  #         low.n, high.n, by.n   desired sample sizes (N per group)         
  #
  #  Note: function uses observed effect size but effect size will
  #        naturally vary some, among bootstrap samples
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.numeric(NS),
    is.numeric(low.n), is.numeric(high.n), is.numeric(by.n) ) 

  dfr <- na.omit( data.frame(y, x))    # delete rows with NA
  y   <- dfr$y ; x <- dfr$x            # reconstitute y & x

  indx      <- nrow(dfr)                       # N of rows
  n         <- seq(low.n, high.n, by.n)        # desired, prospective sample sizes
  bootpwr05 <- bootpwr01 <- numeric(length(n)) # to write-over in 'for' loop

  boot.df   <- function(df, N.boot=NULL){  
    bf      <- df[sample(indx, N.boot, replace=TRUE), ] # resample data frame by its rows
    out     <- anova(lm(bf[,1] ~ bf[,2]))
    boot.p  <- out[1,5]   # one-way Anova p-value returned for bootstrapped data frame
  }

  boot.time <- proc.time()    # begin timer; loop through prospective sample sizes
    for(i in seq_along(n)){
      a <- replicate(NS, boot.df(dfr, N.boot=n[i]))   # p vector
      bootpwr05[i] <-  length(a[a <= .05]) / NS       # statistical power at .05 alpha
      bootpwr01[i] <-  length(a[a <= .01]) / NS       #  "           "    at .01 alpha
    }
  boot.time <- round(((proc.time() - boot.time)/60)[[3]],2) # bootstrap minutes

  cat("\n\nOutput from Kincaid function: boot.powercurve.oneway()\n\n")
  cat("DATA:\t\t", data.name, "\n\nobserved total N:\t", nrow(dfr), "\n") 
  cat("Bootstrap resamples\n  per sample size:\t", NS, "\nMinutes needed:\t",boot.time,"\n\n")

  plot(n, bootpwr05, pch=16, type="b", cex.main=.8, ylim=c(0,1),
    main=paste("Bootstrap power curve after",NS,"bootstrap resamples per sample size"),
    xlab="Total sample size", ylab="Statistical power (one-way Anova)", 
    cex.lab=1.2, font.lab=2, cex.sub=.8,
    sub=paste("bootstrap time:", boot.time, "minutes ", "  ", data.name) )
   
  abline(h=seq(0,1,.1), lty=2, col="gray"); abline(v=n, lty=2, col="gray")
  abline(h=.8, lty=2, lwd=2, col="black")   # horizontal line at power=0.8
  mtext("solid points: alpha .05,  open points: alpha .01", cex=.8, col="black")
  points(n, bootpwr01, type="b", pch=1,  col="black")  # add power curve at .01

  pwr.matrix           <- cbind(n, round(bootpwr05,3), round(bootpwr01,3))
  colnames(pwr.matrix) <- c("total N", "boot.power.05", "boot.power.01")
  return(pwr.matrix)

} # end fun definition


# demo the function

# sim data

# boot.powercurve.oneway(response, drug, NS=1e2, data.name="sim data", low.n=10, high.n=120, by.n=10)

#
# age at walking data  -- despite the fact that sample size is too low to wisely bootstrap
#
# boot.powercurve.oneway(months,treatment,NS=50,data.name="Age at walking",low.n=10,high.n=60,by.n=2)
# boot.powercurve.oneway(months,treatment,NS=1e2,data.name="Age at walking",low.n=10,high.n=60,by.n=4)
# boot.powercurve.oneway(months,treatment,NS=1e3,data.name="Age at walking",low.n=22,high.n=60,by.n=2)




sample.size.barplots <- function(y, x, data.name="My Data", response.name="Response", 
                                 factor.name="Treatment", ... ){
  # ----------------------------------------------------------------------
  #  Bar graphs of sample size, by barchart() in package, lattice
  #                             by barplot() in package, graphics, base R
  #
  #  ARGS      y          numeric response vector, NAs OK
  #            x          corresponding factor vector, NAs OK
  #        data.name      character
  #       response.name   character
  #       factor.name     character
  #           ...         args to pass to funs in base R         
  #
  #  Merely a wrapper around barchart() in package, lattice 
  #  reflecting my poor command of lattice graphics. NEEDS DEVELOPMENT
  #  
  #  For bar graphs of sample sizes, barplot() in base R is preferable.
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.character(data.name), 
    is.character(response.name), is.character(factor.name) ) 

  dfr <- na.omit(data.frame(y, x))  # delete rows with NA
  y   <- dfr$y ; x <- dfr$x            # reconstitute y & x

  require(lattice)
  x           <- factor(x)               # just to make sure x is a factor
  groups      <- length(levels(x))       # number of groups (levels of factor x)
  sample.size <- tapply(y, x, length)    # sample size per group (level)

  plot.new()   # create graph window

  temp <- barchart(sample.size,
    horizontal=TRUE,            # FALSE makes bars vertical
    xlab=NULL,
    ylab=NULL, 
    main=data.name, 
    box.ratio=2,                # 2 is default
    col=gray.colors(groups),
    scales=list(cex=c(1.25, 1.25, 2), font=2, cex.lab=2)  
  )
  print(temp)
  title(main="", sub="Sample Size", cex.sub=1.5, font.sub=2)
  #mtext(data.name, side=3, line=2.5)  # cex=1.5

} # end of function definition


# demo the function

#
# age at first walking data: a=4, total n=23
#

# y1 <- months; x1 <- treatment   # for convenience

# sample.size.barplots(y1, x1)  # minimal call
# sample.size.barplots(y1, x1, data.name="")  
# sample.size.barplots(y1, x1, data.name="Age at First Walking")  
# sample.size.barplots(y1, x1, data.name="Age at First Walking", factor.name="Infants")  
# sample.size.barplots(y1, x1, data.name="Age at First Walking", cex.sub=1.5)  
# sample.size.barplots(y1, x1, data.name="Age at First Walking", cex.sub=1.5, cex=1.5)  
# sample.size.barplots(y1, x1, data.name="Age at First Walking", cex.sub=1.5, cex=2, font=4)  




one.way.histograms.lattice <- function(y, x, 
  data.name="My Data", response.name="Response", hist.type="count", ... ){

  # ---------------------------------------------------------------
  #  Histograms and kernel density estimates for each level of x
  # 
  #  ARGS      y          numeric response vector, NAs OK
  #            x          corresponding factor vector, NAs OK
  #        data.name      character
  #       response.name   character
  #       hist.type       character; "count" or "density"
  #           ...         args to pass to histogram()              ----------      
  #                                                                NEEDS WORK
  #                                                                ----------
  #  requires library(lattice)                                     NEEDS WORK
  #                                                                ----------
  #  Merely a wrapper around histogram() in lattice, reflecting
  #  my poor command of lattice graphics. 
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ---------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.character(data.name), 
    is.character(response.name), is.character(hist.type) ) 

  dfr <- na.omit(data.frame(y, x))  # delete rows with NA
  y   <- dfr$y ; x <- dfr$x         # reconstitute y & x
  x  <- factor(x)                   # just to make sure x is a factor

  require(lattice)  
  plot.new()   # create graph window

  out <- histogram( ~ y | x, 
    type=hist.type,
    #border="white",
    main="",
    ylab="",
    xlab="",
    #col="gray",
    scales=list(cex=c(.9, .9), font=1), ...
  )

  print(out)
  mtext(response.name, side=1, cex=1.3, line=3.6)
  if(hist.type == "count")   mtext("Frequency", side=2, adj=0, cex=1.3, line=3)
  if(hist.type == "density") mtext("Density", side=2, adj=0, cex=1.3, line=3)
  mtext(data.name, side=3, line=2.6, cex=1.5, font=2)

  idle()
  a <- nlevels(x)
  if(a <= 6)
    kernel.out <- densityplot(~ y | x, layout=c(1, a), xlab=response.name, col="black")
  else
    kernel.out <- densityplot(~ y | x, xlab=response.name, col="black")
  print(kernel.out)
   
} # end of function definition


# demo the function


#
# age at first walking data: a=4, total n=23
#
#
# y1 <- months; x1 <- treatment  # for convenience
#
# one.way.histograms.lattice(y1, x1, response.name="Age at first walking, months")  # minimal call
# one.way.histograms.lattice(y1, x1, col="gray")
# one.way.histograms.lattice(y1, x1, data.name="Age at First Walking Experiment", response.name="Months", col="gray", hist.type="density")




means.barplot <- function(y, x, err.bar="CI", prob=.95, verbose=TRUE,
                 bar.lwd=2, bar.length=.075, bar.col="black", ... ){

  # ---------------------------------------------------------------
  #  Graph means with error bars for a one-way layout
  #
  #  ARGS      y         numeric response vector, NAs OK
  #            x         corresponding factor vector, NAs OK
  #         err.bar      character; "CI" , "SEM" or "SD"
  #          prob        numeric; .95 for 95% CI, etc.
  #         verbose      logical; FALSE turns off graph title
  #        bar.lwd       numeric; thickness of error bars
  #        bar.length    numeric; inches
  #        bar.col       character; a color
  #           ...        args to pass to barplot() 
  #
  #  Merely a wrapper around barplot() in base R, but
  #  allows 3 types of ERROR BARS to be plotted: CI, SEM, SD
  #  It accepts most of args recognized by barplot()
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ---------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.character(err.bar), 
    is.numeric(prob), is.logical(verbose), is.numeric(bar.lwd), 
    is.numeric(bar.length), is.character(bar.col) ) 

  dfr <- na.omit(data.frame(y, x))  # delete rows with NA
  y   <- dfr$y ; x <- dfr$x         # reconstitute y & x
  x  <- factor(x)                   # just to make sure x is a factor

  groups <- length(levels(x))    # number of groups (levels of factor)
  N      <- tapply(y, x, length) # vector of sample sizes
  Mean   <- tapply(y, x, mean)   # vector of means per group
  SD     <- tapply(y, x, sd)     # vector of SD per group
  SEM    <- SD / sqrt(N)         # vector of SEM per group

  if(err.bar == "CI"){
    if(verbose == TRUE) 
      err.name=paste("Means with ", round(100*prob), "% CI", sep="") else err.name=""
    t.value <- prob + (1-prob)/2  # e.g., convert prob=.95 to .975 for t
    LB <- Mean - qt( t.value, N-1)*SEM
    UB <- Mean + qt( t.value, N-1)*SEM }

  if(err.bar == "SEM"){
    if(verbose == TRUE) err.name="Means with SEM" else err.name=""
    LB <- Mean - SEM  # vectors
    UB <- Mean + SEM}

  if(err.bar == "SD"){
    if(verbose == TRUE) err.name="Means with SD" else err.name=""
    LB <- Mean - SD  # vectors
    UB <- Mean + SD}

  y.max <- max( c(Mean, UB))  # for plotting
  y.max <- 1.1 * y.max        # fudge factor
   
  barplot(Mean, 
    ylim=c(0, y.max),
    main=err.name, 
    #cex.main=1.2, font.main=2, 
    offset=0, ...
  )

  out <- barplot(Mean, plot=FALSE, ...) # to get bar midpoints along x-axis
  x.coord <- out[ , 1]  # a numeric vector of x-coordinates per bar, L to R
  arrows(x.coord, UB,   # draw error bars
         x.coord, LB, 
         angle=90, code=3, 
         length=bar.length, lwd=bar.lwd, col=bar.col)

  out <- data.frame(N, Mean, SD, SEM, err.bar, LB, UB)
  if(verbose == TRUE) print(out, digits=5, quote=FALSE)

} # end function definition


# demo the function


#
# age at first walking data: a=4, total n=23
#

# y <- months; x <- treatment  # for convenience

# means.barplot(y, x)  # minimal call, default is means with their 95% CI
# means.barplot(y, x, ylim=c(0, 14)) 

# means.barplot(y, x, err.bar="SEM") 
# means.barplot(y, x, err.bar="SEM", ylab="Months", verbose=FALSE)
# means.barplot(y, x, err.bar="SD") 

# means.barplot(y, x, col="lightblue", names=c("A", "B", "C", "D"), verbose=FALSE)  
# means.barplot(y, x, xlim=c(0,7), col.main="white" )   # squish the plot

# means.barplot(y, x, names.arg=c("ACTIVE training", "Passive", "None", "CONTROL" ), col=terrain.colors(4) )   # rainbow()  heat.colors()  topo.colors()  terrain.colors()  gray.colors()
# means.barplot(y, x, col=c("gray", "wheat", "lightblue", "coral" ))  # specify colors bar-by-bar

# means.barplot(y, x, space=.8 )   # bars narrow 
# means.barplot(y, x, space=.2)

# means.barplot(y, x, border=NA )  # no black line around bars
# means.barplot(y, x, width=c(.5, 1, 1.5, 2))

# means.barplot(y, x, err.bar="CI", prob=.95)
# means.barplot(y, x, err.bar="SEM", cex.names=1.5, cex.axis=1.2, cex.lab=1.5, verbose=FALSE, ylab="Age at walking, months")

# a call using many of the args accepted by barplot() in base R
# means.barplot(y, x, err.bar="CI", verbose=FALSE, 
#   sub="Age at First Walking Experiment", cex.sub=1.5, font.sub=4, 
#   ylab="Months", cex.lab=1.3, font.lab=1, cex.axis=1.2, cex.names=1.2,
#   cex.main=1.2, font.main=4, names.arg=c("ACTIVE", "Passive", "None", "CONTROL"),
#   col=heat.colors(4), err.bar.lwd=3, err.bar.length=0 )





Kincaid.stripchart.error.bars <- function(y, x, err.bar="CI", prob=.95, verbose=TRUE,
                            bar.lwd=1, bar.col="red", bar.length=.05, stat.cex=1, ...){

  # ----------------------------------------------------------------------------
  #  Stripcharts of means with error bars are a commonly used  
  #  alternative to boxplots, when n is small        
  #
  #  ARGS      y          numeric response vector, NAs OK
  #            x          corresponding factor vector, NAs OK
  #         err.bar       character; "none","CI" , "SEM" or "SD"
  #          prob         numeric; .95 for 95% CI if err.bar="CI", etc.
  #         verbose       logical; TRUE prints info in graph window
  #  bar.lwd, bar.length  numeric; thickness & length of error bars
  #         bar.col       character; bar color
  #        stat.cex       numeric; size of central point, e.g., the mean
  #           ...         args to pass to stripchart() in base R  
  #                                               
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------------------

  stopifnot(is.numeric(y), is.factor(x),              # data & arg checking
    length(y)>1, length(x)>1, length(y) == length(x), is.character(err.bar), 
    is.numeric(prob), is.logical(verbose), is.numeric(bar.lwd), 
    is.numeric(bar.length), is.character(bar.col), is.numeric(stat.cex) ) 

  dfr <- na.omit(data.frame(y, x))  # delete rows with NA
  y   <- dfr$y ; x <- dfr$x         # reconstitute y & x
  x  <- factor(x)                   # just to make sure x is a factor

  groups <- length(levels(x))    # number of groups (levels of factor)
  N      <- tapply(y, x, length) # vector of sample sizes
  Mean   <- tapply(y, x, mean)   # vector of means per group
  Median <- tapply(y, x, median) # vector of medians per group
  SD     <- tapply(y, x, sd)     # vector of SD per group
  SEM    <- SD / sqrt(N)         # vector of SEM per group

  if(err.bar == "none"){stripchart(y ~ x, method="jitter", jitter=.1, vertical=TRUE, ... )}

  if(err.bar == "Mean"){
    stripchart(y ~ x, method="jitter", jitter=.1, vertical=TRUE, ... )
    for(i in 1:groups) lines(x=c(i-.15, i+.15), y=c(Mean[[i]], Mean[[i]]), lwd=3, 
        col="black")
    if(verbose == TRUE){
        mtext("jittered points with MEAN", side=1, adj=0, line=4, cex=.8, col="black")}}      

  if(err.bar == "Median"){
     stripchart(y ~ x, method="jitter", jitter=.1, vertical=TRUE, ... )
     for(i in 1:groups) lines(x=c(i-.15, i+.15), y=c(Median[[i]], Median[[i]]), 
        lwd=3, col="black")
     if(verbose == TRUE){
        mtext("jittered points with MEDIAN", side=1, adj=0, line=4, cex=.8, col="black")}}      

  if(err.bar == "SEM"){
     stripchart(y ~ x, method="jitter", jitter=.1, vertical=TRUE, ... )
     LB <- Mean - SEM
     UB <- Mean + SEM
     arrows(1:groups,UB, 1:groups, LB, angle=90, code=3, length=bar.length, 
        lwd=bar.lwd, col=bar.col)
     lines(1:groups, Mean, pch=16, type="p", cex=stat.cex, col="black")
     if(verbose == TRUE){
        mtext("jittered points with means and SEM", side=1, adj=0, line=4, 
           cex=.8, col="black")}}      

  if(err.bar == "SD"){
     stripchart(y ~ x, method="jitter", jitter=.1, vertical=TRUE, ... )
     LB <- Mean - SD
     UB <- Mean + SD
     arrows(1:groups,UB, 1:groups, LB, angle=90, code=3, length=bar.length, 
        lwd=bar.lwd, col=bar.col)
     lines(1:groups, Mean, pch=16, type="p", cex=stat.cex, col="black")
     if(verbose == TRUE){
        mtext("jittered points with means and SD", side=1, adj=0, line=4, 
           cex=.8, col="black")}}      

  if(err.bar == "CI"){
     stripchart(y ~ x, method="jitter", jitter=.1, vertical=TRUE, ... )
     t.value <- qt( p=prob+(1 - prob)/2, df=N - 1)  # vector
     LB <- Mean - (SEM * t.value)  # vectors
     UB <- Mean + (SEM * t.value)
     arrows(1:groups,UB, 1:groups, LB, angle=90, code=3, length=bar.length, 
        lwd=bar.lwd, col=bar.col)
     lines(1:groups, Mean, pch=16, type="p", cex=stat.cex, col="black")
     if(verbose == TRUE){
        mtext(paste("jittered points with means and ", round(100*prob), "% CI", sep=""), 
           side=1, adj=0, line=4, cex=.8, col="black")}}      

}  # end function definition

#  demo the function


#
# age at first walking data: a=4, total n=23
#

# Kincaid.stripchart.error.bars(months, treatment)  # minimal call, default 95% CI of mean

# Kincaid.stripchart.error.bars(months, treatment, col=NULL)   # to ERASE points
# Kincaid.stripchart.error.bars(months, treatment, ylim=c(8, 15)) # adjust y-axis plotting range
# Kincaid.stripchart.error.bars(months, treatment, ylim=c(8, 15), xlim=c(0, 5), verbose=F) # adjust x,y range

# Kincaid.stripchart.error.bars(months, treatment, err.bar="none",  
#   pch=21, bg="lightblue", col="black", ylab="Months", xlab="Treatment",
#   main="Age at First Walking Experiment", 
#   cex=1.2, cex.axis=1.2, cex.lab=1.2 )

# Kincaid.stripchart.error.bars(months, treatment, err.bar="Mean", verbose=TRUE,
#   pch=21, bg="lightblue", col="black", ylab="Months", xlab="Treatment", 
#   cex=1.2, cex.axis=1.2, cex.lab=1.2 )

# Kincaid.stripchart.error.bars(months, treatment, err.bar="Median", verbose=TRUE,
#   pch=21, bg="lightblue", col="black", ylab="Months", xlab="Treatment", 
#   cex=1.2, cex.axis=1.2, cex.lab=1.2 )

# Kincaid.stripchart.error.bars(months, treatment, err.bar="SEM", verbose=TRUE,
#   stat.cex=1.2, pch=21, bg="lightblue", col="black", ylab="Months", xlab="Treatment", 
#   cex=1.2, cex.axis=1.2, cex.lab=1.2 )

# Kincaid.stripchart.error.bars(months, treatment, err.bar="SD", verbose=TRUE,
#   pch=21, bg="lightblue", col="black", ylab="Months", xlab="Treatment", 
#   cex=1.2, cex.axis=1.2, cex.lab=1.2, ylim=c(8, 15) )

# Kincaid.stripchart.error.bars(months, treatment, err.bar="CI", prob=.95, verbose=TRUE, 
#   bar.col="black", bar.lwd=2, stat.cex=1.3,
#   pch=21, bg="lightblue", col="black", ylab="Months", xlab="Treatment", 
#   cex=1.2, cex.axis=1.2, cex.lab=1.2, ylim=c(8, 15) )

# Kincaid.stripchart.error.bars(months, treatment, err.bar="CI", prob=.95, verbose=FALSE, 
#   bar.col="black", bar.lwd=3, bar.length=0, stat.cex=1.7, 
#   pch=21, bg="lightblue", col="black", ylab="Months", xlab="Treatment", 
#   cex=1.2, cex.axis=1.2, cex.lab=1.2, ylim=c(8, 15) )
 




sim.data.oneway.Normal <- function(a=4, total.n=92, balanced=TRUE, grand.mean=50, 
  SD=5, seed=sample(1:1e5, size=1), effect.size=.2, precision=.01, 
  verbose=TRUE, graph=TRUE){

  # --------------------------------------------------------------------
  #  Data simulation from N(grand.mean, SD) for one-way data
  #
  #  Call the function with no arguments as text output to R Console
  #  explains it. It simulates from N(grand.mean, SD) one-way data 
  #  of any number of groups ('a') and total N ('total.n'), given
  #  mean ('grand.mean') and sigma ('SD') i.e., RMSE in Anova; all 
  #  informed by observed data, perhaps. Sim data can be balanced 
  #  (equal N per group) or unbalanced. If desired, a random number 
  #  seed argument ('seed') aids reproducibility.
  # 
  #  I decided to use explained R-sq ('effect.size') as acceptance 
  #  criterion to 'break' out of 'repeat' loop, thus grabbing and 
  #  returning a sim data set as data frame object. Other criteria
  #  could be Cohen's d or the number of Tukey outliers.
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # --------------------------------------------------------------------


  timer <- proc.time() # begin timer
  if(effect.size > .3) effect.size <- .3  # put a cap on the target R-sq
  if(balanced == TRUE)if(total.n %% a == 0) n.level <- total.n / a  # modulo %%
  if(balanced == TRUE)if(total.n %% a != 0) n.level <- ceiling(total.n / a)

  if(balanced == FALSE){  # algorithm for reasonable, unequal N per group
    repeat{
      sizes <- floor(total.n / ( a / sample( seq(.5, 1.5, .01), size=a)))
      if(sum(sizes) == total.n) break 
    }
  }  

  # make factor vector of capital letters -----------------------

  if(balanced == TRUE) x <- gl(n=a, k=n.level, length=total.n, labels=LETTERS[1:a])

  if(balanced == FALSE){   
    x    <- character(total.n)  # create character vector to populate in loop
    indx <- 0  # subscript incremented in 'for' loop for 'x'

    for(i in 1:a){ 
      for(j in 1:sizes[i]){
        indx <- indx + 1
        x[indx] <- LETTERS[i]
      }
    }
  } 
  
  # 'repeat' loop ---------------------------------------------------

  set.seed(seed)  # the target effect size: SSgroups/SStotal, the explained R-sq
  counter <- 0    # to count the 'repeat' loops, for the curious

  repeat{                    # the iterator
    counter <- counter + 1
    y <- round(rnorm(total.n, grand.mean, SD), 2) # sample from N(grand.mean, SD)
    SStotal <- (total.n - 1) * sd(y)^2

    m <- tapply(y, x, mean)   # vector of means per group
    n <- tapply(y, x, length) # "      "  n     "   "

    SSg  <- sum( n * (m - mean(y))^2 )  # group SS in Anova
    R.sq <- SSg/SStotal                 # explained R-square; here, our effect size
 
    if( abs( R.sq - effect.size ) <= precision) break  # leave 'repeat' loop if sim OK
  }

  dframe <- data.frame(x, y)  # capture simulated data as data frame
  timer <- round((proc.time() - timer)[[3]]/60, 4)  # minutes elapsed

  if(verbose == TRUE){
    cat("\n\n\nRUNNING 'sim.data.oneway.Normal()' by Dwight Kincaid, PhD\n\n")
    cat("ARGUMENTS in the call:\n")
    cat("a =", a, " total.n =", total.n, " balanced =", balanced, 
      " grand.mean =", grand.mean, " SD =", SD, "\nseed =", seed,
      " effect.size =", effect.size, 
      " precision=", precision, "\nverbose =", verbose, " graph =", graph, "\n\n")
    cat("NOTE: 'effect.size' is the explained R-square in one-way Anova.\n")
    cat("R-square is the proportion of the total variation (total SS) in the\n")
    cat("response variable that is explained by group differences (group SS).\n")
    cat("R-square is a widely used measure of heterogeneity among the means.\n\n")
    cat("The simulated data has the R-square of 'effect.size', within the\n")
    cat("value of 'precision' but change 'precision' in the call as desired.\n") 
    cat("Simulation time depends mostly on 'effect.size' (the smaller the faster)\n")
    cat("and on 'precision' (the larger the faster) but not on total N, 'total.n'\n")
    cat("or on 'a' the number of groups.\n\n") 
    cat("In simulating one-way data, there is probably no good reason for\n") 
    cat("'effect.size' to exceed R-sq=0.25; the function will not let this arg\n")
    cat("exceed .30, for practical computational purposes. Use the 'seed'\n")
    cat("argument (must be an integer) for reproducibility if desired. Remember,\n")
    cat("data simulation has many uses relative to probing and understanding\n")
    cat("observed data; it also is the best way to learn statistical methods.\n\n")
 
    if(balanced == TRUE) if(total.n %% a != 0){
      cat("WARNING: ") 
      cat("You specified a BALANCED design (= n per group). Yet\n",
      "your total.n/a does not allow it, therefore the sim data will be\n", 
      "nearly balanced. Otherwise, change your arguments in the call.\n\n")
    }
   
    cat("Some DIAGNOSTICS of the sim data frame.\n\n")
    cat(str(dframe)); cat("\n")
    print(summary(dframe)); cat("\n")
  
    cat("The Means\n")
    print(tapply(dframe$y, dframe$x, mean)); cat("\n")
    cat("The Standard Deviations\n")
    print(tapply(dframe$y, dframe$x, sd)); cat("\n\n")

    print(anova(lm(y~x, dframe)), signif.stars=FALSE)
    cat("\nExplained R-square = SS.groups/SS.total = ")
    cat(summary((lm(y~x, dframe)))$r.squared, ", for sim data set.\n\n")
    cat("Minutes elapsed:", timer, ",", " interations:", counter, "required in 'repeat' loop\n")
    cat("< End of run >  ", date(), "\n\n") 
  }

  if(graph == TRUE) boxplot(y~x, dframe, col=terrain.colors(nlevels(x)))

  return(dframe)  # return simulated data frame

} # end function definition


# demo the function

# my.frame <- sim.data.oneway.Normal(effect.size=.2)
# str(my.frame)
# summary(my.frame)
# head(my.frame)

# my.frame <- sim.data.oneway.Normal(effect.size=.2, balanced=FALSE)

# use the random number seed arg for reproducibility
# my.frame <- sim.data.oneway.Normal(seed=1234) 

#
# Age at first walking data: a=4, total n=23
#

# my.frame <- sim.data.oneway.Normal(a=4, total.n=23, grand.mean=11.34783, SD=1.51641, effect.size=.2528, precision=.001)
# my.frame <- sim.data.oneway.Normal(a=4, total.n=23, grand.mean=11.34783, SD=1.51641, effect.size=.2528, precision=.001, balanced=FALSE)
# my.frame <- sim.data.oneway.Normal(a=4, total.n=32, grand.mean=11.34783, SD=1.51641, effect.size=.2528, precision=.001, balanced=TRUE)

# view boxplot for many sim data sets for 'Age at Walking' 
# NS <- 20
# for(i in 1:NS){
#   my.frame <- sim.data.oneway.Normal(a=4, total.n=32, grand.mean=11.34783, 
#   SD=1.51641, effect.size=.2528, precision=.001, balanced=FALSE, verbose=FALSE, graph=TRUE)
#   Sys.sleep(1)  # pause 1 second
# }





sim.toy.oneway.data <- function(a=5L, n=5L, my.seed=123L, effect=c(-3,3), graphs=TRUE){

  # -------------------------------------------------------------------
  # Returns a 2 column data frame of simulated and balanced, oneway
  # toy data useful for testing code and for statistical learning. 
  # This sim data tends toward monotonic positive change in means.  
  # Graphs richly annotated. Call without args as a demo.
  #
  # ARGS   a      number of groups, i.e., levels; default 5
  #        n      sample size equal per group; default 5
  #     my.seed   random number seed; default 123
  #      effect   default c(-3, 3) as range of noise added to response
  #               variable; increase range to decrease heterogeneity
  #               and monotonicity; decrease it to increase 
  #               heterogeneity and monotonicity among group means 
  #      graphs   logical; boxplot followed by means with 95% CI 
  # 
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # -------------------------------------------------------------------

  total.n <- a*n        # total sample size
  y       <- 1:total.n  # populate 'response' vector with integers then add noise

  set.seed(my.seed)
  # add noise to response y  
  y <- y + round(runif(total.n, min=effect[1], max=effect[2]), 1) 

  x <- gl(n=a, k=n, labels=c("Control", LETTERS[1:(a-1)])) # grouping vector

  if(graphs == TRUE){
    idle <- function(){
      if(interactive())readline("\nHit<Enter>to continue...");invisible()}

    obj <- oneway.test(y ~ x, var.equal=TRUE)  # oneway Anova
    F   <- signif(obj[[1]], 5)                 #  stats  
    df1 <- obj$parameter[[1]]                  #   for 
    df2 <- obj$parameter[[2]]                  #    annotating
    p   <- signif(obj[[3]], 5)                 #     graphs

    main.txt <- "sim.toy.oneway.data( )"; y.txt <- "Response"
    aov.txt  <- paste("oneway Anova:  F", "(", df1, ",", df2, ") =", F, ",  p =", p )
    call.txt <- paste("args:  a =", a, " n =", n, " random number seed =", my.seed, 
                      " effect =", effect[1], ",", effect[2], sep="  ") 

    if(a < 5) cx.ax <- 1 else cx.ax <- .8; if(a > 10) cx.ax <- .6; if(a > 15) cx.ax <- .45
     
    boxplot(y ~ x, boxwex=.5, col="gray", ylab=y.txt, cex.axis=cx.ax)

    title(main=main.txt, font.main=2, col.main="red",
          sub=aov.txt, font.sub=2, col.sub="red" )
    mtext(call.txt, cex=.8)

    idle(); library(gplots)
    plotmeans(y ~ x, connect=FALSE, ylab=y.txt, xlab="", n.label=FALSE,
      main=main.txt, cex.axis=cx.ax)
  
    mtext(" means with 95% CI", side=3, line=-1, cex=.8, col="blue", adj=0)
    mtext(call.txt, cex=.8)
    title(main=main.txt, font.main=2, col.main="red",
          sub=aov.txt, font.sub=2, col.sub="red" )
  }

  return(data.frame(y, x))

} # end fun definition

# demo the fun

# df <- sim.toy.oneway.data( )
# df <- sim.toy.oneway.data(graphs=FALSE )
# df <- sim.toy.oneway.data(my.seed=45 )
# df <- sim.toy.oneway.data(a=4)
# sim.toy.oneway.data(a=4, n=3 ) 
# df <- sim.toy.oneway.data(a=4, effect=c(-7, 7), n=5 )

# sim.toy.oneway.data(my.seed=8736, a=9, n=5, effect=c(-30, 30) )
# sim.toy.oneway.data(my.seed=8736, a=15, n=5, effect=c(-30, 30) )


# df <- sim.toy.oneway.data(a=20, n=4, my.seed=87356, effect=c(-200, 200)) # p=.0186
# TukeyHSD( aov(df$y ~ df$x), ordered=TRUE )  #all 190 comparisons not significant
# summary(lm(df$y ~ df$x)) # planned comparisons between 'Control' and each group





anova.by.stats <- function( n, m, s, digits=3 ){

  # --------------------------------------------------------------------
  #  One-way Anova given 3 vectors: n, mean and SD for one-way data
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # --------------------------------------------------------------------

  #  n  is vector of sample sizes
  #  m  "  "      "  means
  #  s  "  "      "  standard deviations

  SS      <- s^2 * (n-1)     # SS for each group
  N       <- sum(n)          # total N
  gmean   <- sum(n * m) / N  # grand mean

  a       <- length(n)       # number of groups
  totaldf <- N-1             # total degree of freedom
  df1     <- a - 1           # group   "    "   "
  df2     <- totaldf - df1   # error   "    "   "

  SSgroup <- sum(n * (m - gmean)^2)  # explained sum of squares
  SSerror <- sum(SS)                 # unexplained sum of squares
  SStotal <- SSgroup + SSerror       # total sum of squares

  MSgroup <- SSgroup / df1  # variance due to groups
  MSerror <- SSerror / df2  # variance within groups

  F    <- MSgroup / MSerror                    # the F test statistic
  p    <- pf(F, df1, df2, lower.tail=FALSE)    # p-value of F
  R.sq <- round( 100 * SSgroup / SStotal  , 3) # explained R-sq %

  SSgroup <- round(SSgroup, digits) # round-off before printing
  SSerror <- round(SSerror, digits)
  SStotal <- round(SStotal, digits)
  MSgroup <- round(MSgroup, digits)
  MSerror <- round(MSerror, digits)
  F       <- round(F, digits)

  t.stat   <- qt(.975, n-1)       # 95% CI of means using normal theory
  SEM      <- s/sqrt(n)
  LB.95.CI <- m - t.stat * SEM
  UB.95.CI <- m + t.stat * SEM

  cat("\n\nDescriptive Statistics for each group (N, Mean, SD)\n")
  cat("plus 95% CI of the means using t-distibution\n\n")
  print(cbind(n, m, s, LB.95.CI, UB.95.CI))

  cat("\n\nOne-way Anova table\n\n")
  cat("\t df \t SS \t MS \tF \t P \n")
  cat("-------------------------------------------------\n")
  cat("groups\t", df1, "\t", SSgroup,  MSgroup, F, p,    "\n\n")
  cat("error\t", df2,  "\t", SSerror, MSerror, "\n")
  cat("------------------------\n")
  cat("total\t", totaldf, "\t", SStotal,   "\n\n\n")
  cat(paste( "R-square = ", R.sq, "% explained", sep=""), "\n\n")

} # end of fun definition

# demo the fun

#
# age at first walking data: a=4, total n=23
#
# n    <- c(6,      6,      6,      5    )  # store descriptive stats                    
# Mean <- c(10.125, 11.375, 11.708, 12.35)  # in three
# SD   <- c(1.447,  1.896,  1.52,   .9617)  # vectors

# anova.by.stats(n, Mean, SD, digits=4)
# anova.by.stats(n=tapply(months, treatment, length),
#                M=tapply(months, treatment, mean),
#                s=tapply(months, treatment, sd ))





Kincaid.error.bars.by.stats <- function(n, m, s, err.bar="SEM", prob=.95, x.names=NULL,
  my.pch=16, verbose=TRUE, bar.lwd=1, bar.col="red", bar.length=.05, point.col="black", 
  stat.cex=NULL, connect=FALSE, add.n=FALSE, add.n.cex=.8, bar.buffer=FALSE, bar.gap=.02, ...){

  # -------------------------------------------------------------------
  #  Given 3 vectors: n, mean and SD, then plot means with error bars 
  #
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # -------------------------------------------------------------------
  #
  # ARGUMENTS
  # n           vector of sample sizes         # DESCRIPTIVE
  # m           vector of means                # STATISTICS as
  # s           vector of standard deviations  # numeric VECTORS

  # err.bar     type of error bar: "SEM", "SD", "CI" ; "SEM" is default
  # prob        type of CI using normal theory, .95 for 95%, .99 for 99%, etc.
  # x.names     character vector of group names, or say LETTERS[1:4]; integer default
  # my.pch      plotting symbol for the means, pch=16 default 
  # verbose     a logical, TRUE prints graph title

  # bar.lwd     thickness of error bar
  # bar.col     color of error bar
  # bar.length  length of bar end; use '0' for no bar end

  # point.col   color of plotted symbol for mean
  # stat.cex    size of plotted symbol for mean
  # connect     a logical, TRUE connects the means with a line
  # add.n       a logical, TRUE adds 'n=_' per group, below error bars
  # add.n.cex   font size for 'add.n"   add.n.cex=.8  is default
  # bar.buffer  a logical, TRUE provides a space between error bars and the means
  # bar.gap     size of space between error bars and the means, default is .02 of ylim range
  # ...         Args specific to stripchart() and points() can be used; see demo code

  a      <- length(n)    # number of groups
  N      <- sum(n)       # total n
  sem    <- s / sqrt(n)  # SEM vector
  t.stat <- qt( p=prob+(1 - prob)/2, df=n-1)  # vector of critical t-values

  if(err.bar == "SEM"){ LB <- m - sem;          UB <- m + sem          }
  if(err.bar == "SD") { LB <- m - s;            UB <- m + s            }
  if(err.bar == "CI") { LB <- m - (sem*t.stat); UB <- m + (sem*t.stat) }

  stripchart(m ~ c(1:a), cex=stat.cex, col=point.col,
    vertical=TRUE, 
    pch=my.pch, group.names=x.names, 
    ylim=range(c(m, LB-(.01*LB), UB)), xlim=c(.5, a+.5), ...)

  if(bar.buffer == TRUE){
    gap <- bar.gap * (max(UB) - min(LB))
    arrows(1:a,UB,  1:a,m+gap, angle=90, code=1, length=bar.length, lwd=bar.lwd, col=bar.col)
    arrows(1:a,LB,  1:a,m-gap, angle=90, code=1, length=bar.length, lwd=bar.lwd, col=bar.col)
  } 
   
  if(bar.buffer == FALSE){arrows(1:a,UB, 1:a, LB, angle=90, code=3, length=bar.length,
    lwd=bar.lwd, col=bar.col)}

  if(verbose == TRUE & err.bar != "CI"){   # plot either SEM or SD
    mtext(paste("Means with", err.bar), side=3, cex=1.2, line=1)}

  if(verbose == TRUE & err.bar == "CI"){ 
      mtext(paste("Means with", round(100*prob), "% CI"), side=3, cex=1.2, line=1)} 

  if(connect == TRUE) lines(1:a, m, type="l", col="gray")

  points(x=c(1:a), y=m, col=point.col, pch=my.pch, cex=stat.cex, ...) # overplot means

  if(add.n == TRUE){
    for(i in 1:a) mtext(paste("n=", n[i], sep=""), at=i, side=1, line=-1, cex=add.n.cex)}

}  # end function definition.  Kincaid.error.bars.by.stats()


# demo the function

#
# Age at first walking data: a=4, total n=23
#

# N    <- c(6,      6,      6,      5    )  # store descriptive stats   
# Mean <- c(10.125, 11.375, 11.708, 12.35)  # in three
# SD   <- c(1.447,  1.896,  1.52,   .9617)  # vectors
# my.names <- c("Active", "Passive", "None", "Control")
# y.string <- "Age at first walking, months"

# Kincaid.error.bars.by.stats(N, Mean, SD) # minimal call
# Kincaid.error.bars.by.stats(N, Mean, SD, x.names=my.names, ylab=y.string)
# Kincaid.error.bars.by.stats(N, Mean, SD, err.bar="SD", x.names=my.names)
# Kincaid.error.bars.by.stats(N, Mean, SD, err.bar="CI", prob=.95, x.names=my.names)

#Kincaid.error.bars.by.stats(N, Mean, SD, err.bar="CI", x.names=my.names,
#  verbose=FALSE, xlab="Treatment", ylab="Age at first walking, months", 
#  cex.lab=1.2, bar.length=0, bar.lwd=1.5, bar.col="gray", stat.cex=1.5)

# Kincaid.error.bars.by.stats(N, Mean, SD, x.names=my.names, add.n=TRUE, add.n.cex=1 )

# Kincaid.error.bars.by.stats(N, Mean, SD, err.bar="CI", x.names=my.names, add.n=TRUE,
#   bar.buffer=TRUE, bar.gap=.03, connect=TRUE, ylab=y.string, cex.lab=1.2, cex.axis=1.2)


# although crude, multiply-out the sample size just to see what happens
# Kincaid.error.bars.by.stats(2*N, Mean, SD, err.bar="CI", x.names=my.names, add.n=TRUE)
# Kincaid.error.bars.by.stats(4*N, Mean, SD, err.bar="CI", x.names=my.names, add.n=TRUE)




# ----------------------------------------------------------------------------
#  my.pause()  is called with no arguments
#  idle()      like my.pause() but more compact and with preceeding linefeed
#
#  both pause text output to R Console and pause between graphs
#
#  modified from DAAG and sm
# ----------------------------------------------------------------------------

my.pause <- function(){
    if(interactive()) readline("Hit <Enter> to continue...")
    invisible()
    #cat("\n")
}   # end of function definition

# demo the function

# my.pause()   # no args in the call


# pause between graphs but with leading linefeed

idle <- function(){if(interactive())readline("\nHit<Enter>to continue...");invisible()}

# idle()


# ---------------------------------
#  simple, line drawing functions 
# ---------------------------------

thick.line <- function(N=25) cat(rep("=", N), sep="", "\n")  # prints "="

thin.line <- function(N=25) cat(rep("-", N), sep="", "\n")  # prints "-"

# demo the funs

# thick.line()
# thin.line(80)

# ------------------------------------------------------------
# skew and kurtosis using Kincaid code optimized for speed 

  skw <- function(y){                  # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)          # assumes no missing values
    (sum(y^3) / length(y)) / sd(y)^3}

  krt <- function(y){                  # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)          # assumes no missing values
    (sum(y^4) / length(y)) / sd(y)^4 - 3}




mytick <- function (nx = 2, ny = 2, tick.ratio = 0.5){

  # --------------------------------------------------------------------
  # ARGS    'nx' & 'ny' are N of desired intervals between major ticks
  #                 in x and in y, respectively; default is 2
  #         'tick.ratio' is minor tick length relative to major tick
  # 
  # a mere capture of fun definition for minor.tick() from 
  # library(Hmisc) because students have issues downloading Hmisc
  #
  # Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # --------------------------------------------------------------------

    ax <- function(w, n, tick.ratio) {
        range <- par("usr")[if (w == "x") 
            1:2
        else 3:4]
        tick.pos <- if (w == "x") 
            par("xaxp")
        else par("yaxp")
        distance.between.minor <- (tick.pos[2] - tick.pos[1])/tick.pos[3]/n
        possible.minors <- tick.pos[1] - (0:100) * distance.between.minor
        low.candidates <- possible.minors >= range[1]
        low.minor <- if (any(low.candidates)) 
            min(possible.minors[low.candidates])
        else tick.pos[1]
        possible.minors <- tick.pos[2] + (0:100) * distance.between.minor
        hi.candidates <- possible.minors <= range[2]
        hi.minor <- if (any(hi.candidates)) 
            max(possible.minors[hi.candidates])
        else tick.pos[2]
        axis(if (w == "x") 
            1
        else 2, seq(low.minor, hi.minor, by = distance.between.minor), 
            labels = FALSE, tcl = par("tcl") * tick.ratio)
    }
    if (nx > 1) ax("x", nx, tick.ratio = tick.ratio)
    if (ny > 1) ax("y", ny, tick.ratio = tick.ratio)
    invisible()
} # end fun definition


# demo the fun

# hist(rnorm(1e4))
# mytick(2, 2)

# mytick(4, 4)
# mytick(2, 2, tick.ratio=1)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# ----------------------------------------------------------------
#  SIMULATED dummy data to exercise calls to the functions above
#  4 groups in an unbalanced, one-way layout with total N=78
# ----------------------------------------------------------------

# set.seed(123)
# drug.A  <- rnorm(20, mean=50, sd=5) # random samples from normal distribution
# drug.B  <- rnorm(15, mean=50, sd=5)
# drug.C  <- rnorm(25, mean=50, sd=5)
# control <- rnorm(18, mean=47, sd=4)  # try, mean=40

# response <- c(drug.A, drug.B, drug.C, control)
# drug     <- c( rep("A", 20), rep("B", 15), rep("C", 25), rep("Control", 18))
# drug     <- factor(drug, levels=c("A", "B", "C", "Control"))

# my.frame <- data.frame(response, drug)
# str(my.frame)

# anova(lm(response ~ drug, data=my.frame)) 

#Analysis of Variance Table
#
#Response: response
#          Df  Sum Sq Mean Sq F value  Pr(>F)  
#treatment  3  221.73  73.909  3.6596 0.01614 *
#Residuals 74 1494.50  20.196                  
#

# require(psych)
# print(describeBy(my.frame$response, my.frame$drug, skew=FALSE, mat=TRUE), digits=3)
#
#   item  group1 vars  n mean   sd median trimmed  mad  min  max range    se
#11    1       A    1 20 50.7 4.86   50.6    50.8 4.35 40.2 58.9  18.8 1.087
#12    2       B    1 15 49.5 4.62   48.9    49.6 6.30 41.6 56.3  14.7 1.193
#13    3       C    1 25 50.5 4.38   49.8    50.4 3.14 42.3 60.8  18.6 0.875
#14    4 Control    1 18 46.5 4.10   45.8    46.5 4.15 37.8 55.2  17.4 0.967

 
# boxplot(response ~ drug, data=my.frame, col=gray.colors(nlevels(treatment)))

# R.sq.oneway(response, drug)   
# vc.oneway(response, drug)
# Cohen.D.oneway(response, drug)
# CI.means.normal.theory(response, drug, prob=.95)


# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
#
# The data are AGE measured in MONTHS at FIRST WALKING
# data(zelazo) from package ISwR by Peter Dalgaard but here we enter data manually
# Data source: Zelazo et al. 1972. Walking in the newborn. Science 176:314-315  
# 
# 'active'  N=6 babies; active training on walking & reflexes; 
#           4, 3 min. sessions/day for the 2-8th weeks of life
# 'passive' N=6 babies; same social & gross motor stimulation but not active training
# 'none'    N=6 babies; no training but were tested with babies who got 'active' or 'passive'
# 'control' N=5 babies; no training & were only tested at age 8 wks; 8 week control group

# data.descriptor <- "Age at Walking"         # handy for output
#
# active  <- c(9, 9.5, 9.75, 10, 13, 9.5)     # N=6 babies, age in months at first walking
# passive <- c(11, 10, 10, 11.75, 10.5, 15)   # N=6 
# none    <- c(11.5, 12, 9, 11.5, 13.25, 13)  # N=6 
# control <- c(13.25, 11.5, 12, 13.5, 11.5)   # N=5 

# months    <- c(active, passive, none, control)  # the response column
# treatment <- c( rep("Active", 6), rep("Passive", 6), rep("None", 6), rep("Control", 5))

# specify desired ORDER, i.e., 'levels', for group names in output, else it's A-Z
# treatment <- factor(treatment, levels=c("Active", "Passive", "None", "Control"))

# DF <- data.frame(months, treatment)  # construct data frame of 2 columns
# str(DF)

# obs.anova.table <- anova(lm(months ~ treatment, data=DF))
# obs.anova.table
# Analysis of Variance Table
#
# Response: months
#           Df Sum Sq Mean Sq F value Pr(>F)
# treatment  3 14.778  4.9259  2.1422 0.1285
# Residuals 19 43.690  2.2995   
#
# R.sq.oneway(months, treatment)  # 25.3 and 74.7 %

# end.  Kincaid funs ONE WAY.R      Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
