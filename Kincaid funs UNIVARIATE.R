# Kincaid funs UNIVARIATE.R    24 March 2018    Author: Dwight Kincaid
#                                                  dwight.kincaid@lehman.cuny.edu
#
#  DEMO CODE after each function definition should be studied and run.
#
#  DEMO DATA includes the N=298, data(IgM) in library(ISwR) as well as sampling
#  from the normal distribution, rnorm(n, mean, sd). Sampling from the Weibull 
#  distribution generates realistic demo data: rweibull(n=50, shape=2.5, scale=1). 
#
#  Generally I maximized use of base R functions to reduce reliance on contributed 
#  CRAN packages for several reasons, including code transparency for students and 
#  frankly, for code 'survival' albeit at slightly slower speed for computationally 
#  intensive tasks although this is usually trivial unless N is large and/or 
#  the re/sampling is huge.  
# ----------------------------------------------------------------------------------
#
##
##  FUNCTION  DEFINITIONS  --  for use on a single, numeric sample: a vector
##
#
#      NEW FUNCTIONS                                                          
# 
#   1. freq.table(y, bins, roundoff)       # frequency distribution table                                  
#
#   2. EDA.graphs.one.sample(y, data.name) # 4 EDA graphs, heavily annotated                                                       
#
#   3. CV.percent(y)                       # coefficient of variation as %
#
#   4. classical.CI.mean(y)                # 90, 95, 99% CI of mean using normal theory
#   5. boot.CI.mean(y, NS)                 # bootstrap CI of mean
#
#   6. my.stats(y)                                # some descriptive stats 
#   7. new.stats(y, roundoff, data.name, norm.test) # descriptive stats and more
#   8. normal.QQ.plot(y)                          # graphical assessment of normality
#
#   9. %mc.skew%      # binary operator; Monte Carlo simulation test for SKEW
#  10. %mc.kurt%      # binary operator; Monte Carlo simulation test for KURTOSIS
#
#  11. %boot.skew%    # binary operator; bootstrap confidence interval for SKEW
#  12. %boot.kurt%    # binary operator; bootstrap confidence interval for KURTOSIS
#                       
#  13. Tukey.outliers(y) # identify outliers by Tukey's 1.5 IQR rule in boxplot()
#
#  14. %bootMedian95CI%  # binary operator; percentile bootstrap 95% CI of median
#
#  15. skw(y)            # skewness (type 3); fast Kincaid code in base R
#  16. krt(y)            # kurtosis (type 3); fast Kincaid code in base R
#
#  --  purely convenience functions  --
#
#  17. thick.line(N)      # draws line of N copies of "="                                            
#  18. thin.line(N)       # draws line of N copies of "-"
#  19. my.pause()         # call with no args; pauses between graphs
#  20. mytick(nx=2, ny=2, tick.ratio=0.5)  # add minor ticks to base R graphs
#
#  --  functions to be added  --
#
#  xx. Agostino.test()                 xx. my.banner()
#  xx. leave.one.out.univariate()      xx. my.footer.header()
#  xx. simUnivarNorm()


freq.table <- function(y=rnorm(1e2, 50, 10), bins="Sturges", roundoff=4 ){

  # ----------------------------------------------------------------------
  #  ARGUMENTS   produces frequency distribution table for numeric vector
  #
  #  y         numeric vector to process for frequency distribution table
  #
  #  bins      value for 'breaks' in arg passed to hist()
  #            default: "Sturges" but can be integer or vector
  #  roundoff  default 4 decimal places
  #
  #  Function RETURNS a data frame object to print 
  #  
  #  Note: for the bins, it's   > ... <=   left open, right closed
  #        function may be called devoid of arguments, as a DEMO
  #   
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove NA, missing values

  # query 'bins' as to whether character or single numeric or vector
  if(is.character(bins)) bins="Sturges" # default binning algorithm in hist()
  if(is.numeric(bins))   bins=bins

  a <- hist(y, breaks=bins, plot=FALSE) # get the hist() object, not the graph
  class.width <- a$mids[2] - a$mids[1]  # width of histogram bins
  L.level     <- round( a$mids - (class.width/2), roundoff) # vector of lower bounds
  U.level     <- round( a$mids + (class.width/2), roundoff) # vector of upper bounds
  my.level    <- paste( L.level, ",", U.level, sep="") # character vector for table

  N <- length(y) # N of response variable

  # frequency distribution DATA FRAME by extraction from the hist() object
  dfr <- data.frame(interval  = 1:length(a$mids),
                    class     = my.level,  
                    midpoint  = a$mids, 
                    FREQ      = a$counts,
                    rel.freq  = round(a$counts / N, roundoff),
                    cum.f     = cumsum(a$counts), 
                    cum.rel.f = round(cumsum(a$counts / N), roundoff),
                    density   = round(a$density, roundoff))

  return(dfr)   # data frame object is returned 

} # end function definition


# demo the function  ----------------------------------------

# freq.table( )   # as a DEMO, call the function with no args

# my.freq <- freq.table(bins="Sturges")
# my.freq

# my.freq <- freq.table(bins=5, roundoff=6)
# my.freq

# library(ISwR); data(IgM)   # N=298
# my.freq <- freq.table( IgM, bins=18 )
# print(my.freq, digits=2)


# freq.table(IgM)  
# freq.table(IgM, bins=c(0,1,2,3,4,5) )
#  interval class midpoint FREQ rel.freq cum.f cum.rel.f density    # text output
#1        1   0,1      0.5  237   0.7953   237    0.7953  0.7953
#2        2   1,2      1.5   55   0.1846   292    0.9799  0.1846
#3        3   2,3      2.5    5   0.0168   297    0.9966  0.0168
#4        4   3,4      3.5    0   0.0000   297    0.9966  0.0000
#5        5   4,5      4.5    1   0.0034   298    1.0000  0.0034

# freq.table(y=rnorm(1e7), roundoff=5)



EDA.graphs.one.sample <- function (y=rnorm(sample(20:100, size=1)), data.name=""){

  # ----------------------------------------------------------------------------  
  #  ARGUMENTS
  #
  #  y          numeric vector, a single sample of data, NAs OK
  #  data.name  descriptive string, e.g., "IgM, g/L"
  #
  #  Function RETURNS a graph window of 4 heavily annotated graphs and can be
  #  called with no args, as a DEMO. Requires library(DescTools) for robust, 
  #  Jarque Bera normality test. Function inspired by QuantPsych::eda.uni()
  #
  #  EDA: Exploratory Data Analysis refers to methods developed by John Tukey
  #  et al. for data exploration & visualization. Today the analyst does EDA
  #  to probe and to understand the data and to enhance hypothesis generation.
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------------------

  if( missing(y) ) data.name <- "Simulated random data, N(0,1)"
  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove NA, missing values
  if(nchar(data.name) > 27) data.name <- paste(substr(data.name,1,24),"...",sep="")
  require(DescTools) 

  skw <- function(y){                # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^3) / length(y)) / sd(y)^3}

  krt <- function(y){                # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^4) / length(y)) / sd(y)^4 - 3}


  par( mfrow = c(2, 2) )  # 2 rows by 2 columns of graphs

    hist( y, main = "Histogram", cex.main=1.2, border="white",
      col="lightblue", xlab = data.name, cex.axis=.8,
      sub=paste("skew=", round(skw(y),3),  "  kurtosis=", 
          round(krt(y),3)), cex.sub=1, col.sub="blue" )
 
    if(length(y) <= 1e3) rug(y, col="gray")  # rug under histogram if N is low

    mtext(paste("N=", length(y), " Mean=", round(mean(y),4), 
                " SD=", round(sd(y),4)), cex=.6, col="blue")

    plot(density(y, n=1024), main = "Kernel Density Estimate", 
         cex.main=1.2, cex.axis=.8, col="red", lwd=3 )    
    
    boxplot(y, horizontal = TRUE,
        col="wheat", xlab=data.name, main= "Boxplot", ylab="Response", boxwex=.9,
        sub=paste("IQR=", round(IQR(y),3),
                  " min=", round(min(y),5), " max=", round(max(y),5)), 
        cex.sub=.8, col.sub="blue", cex.axis=.8 )

    mtext(paste("hinges: 25th=", round(quantile(y,.25),3),  
                " 75th=", round(quantile(y,.75),3)), cex=.65, col="blue")

    outliers <- length( boxplot(y, plot=FALSE)$out)  # extract from boxplot object

    mtext(paste(" N of outliers=", outliers, "  >1.5*IQR from hinge"), 
          side=3, line=-1, adj=0, cex=.6, col="blue")

    mtext(paste("median=", round(median(y),4) ), 
          side=1, line=-1, cex=.7, col="blue")

    a <- DescTools::JarqueBeraTest(y, robust=T, method="chisq")$p.value 
    if(a == 0) a <- "<2.2e-16"  # if it's zero then make it a string to print
   
    # ------------------------------ Normal QQ plot to visualize normality
    if(is.character(a)){
        qqnorm(y, col="blue", sub=paste("Robust Jarque-Bera test p", a), 
        ylab="Observed Quantiles", cex.sub=.8, col.sub="blue", cex.axis=.8 )}
    else

   {qqnorm(y, col="blue", sub=paste("Robust Jarque-Bera test p=", 
        signif(DescTools::JarqueBeraTest(y, robust=T, method="chisq")$p.value,5)), 
        ylab="Observed Quantiles", cex.sub=.8, col.sub="blue", cex.axis=.8 )}

    # ------------------------------

    if(length(y) <=  5000) mtext(paste("Shapiro-Wilk test p=", 
        signif(shapiro.test(y)$p.value,5)), cex=.7, col="blue")

    if(length(y) > 5000) mtext(paste("Lillifors Kolmogorov-Smirnov test p=", 
        signif(LillieTest(y)$p.value,5)), cex=.7, col="blue") 

    qqline(y, col="blue")  # add straight line to Normal QQ plot 
    
    mtext(" Visualization of observed distribution\n relative to normal distribution",
        side=3, line=-1.5, adj=0, cex=.6, col="lightblue4")

  par(mfrow = c(1, 1))  # reset to default of 1 graph per page

} # end function definition


# demo the function  ----------------------------------------

# EDA.graphs.one.sample()  # call function with no args, to see what happens

# EDA.graphs.one.sample(rnorm(1e2, 50, 10), data.name="Sim data from N(50,10)")

# library(ISwR); data(IgM)  # N=298
# EDA.graphs.one.sample(IgM)
# EDA.graphs.one.sample(IgM, data.name="IgM, g/L")  
# EDA.graphs.one.sample(IgM, data.name="The quick brown fox jumped over the lazy dog")


CV.percent <- function(y){
 
  # ----------------------------------------------------------------- 
  #  ARGUMENTS   Useful as an easy 'first' function for students
  #
  #  y         numeric vector, missing values OK
  #
  #  Function RETURNS Coefficient of Variation as percent 
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # -----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[ !is.na(y) ]                      # eliminate any missing values, NA
  abs( 100 * sd(y)/mean(y))                # to avoid a negative number

} # end fun definition


# demo the function  ----------------------------------------

# library(ISwR); data(IgM)   # N=298
# CV.percent(IgM)

# yy <- c(1, 2, 4, NA, 88)  # can it handle missing values? YES
# CV.percent(yy)

# yy <- rnorm(1000, 50, 5)
# CV.percent(yy) 





classical.CI.mean <- function(y){

  # -----------------------------------------------------------------   
  #  ARGUMENTS
  #
  #  y         numeric vector, missing values OK
  #
  #  Function RETURNS 90, 95, 99% CI of mean using t-distribution  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # -----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]  # remove NA, missing values

  N  <- length(y); Mean <- mean(y);    Min <- min(y)
  SD <- sd(y);     SEM  <- SD/sqrt(N); Max <- max(y)

  LB.90 <- Mean - qt(p=.95, df=N-1) * SEM    # 90% CI of mean
  UB.90 <- Mean + qt(p=.95, df=N-1) * SEM    # using t-distribution
  CI.90 <- cbind(LB.90, UB.90)

  LB.95 <- Mean - qt(p=.975, df=N-1) * SEM   # 95% CI of mean
  UB.95 <- Mean + qt(p=.975, df=N-1) * SEM
  CI.95 <- cbind(LB.95, UB.95)

  LB.99 <- Mean - qt(p=.995, df=N-1) * SEM   # 99% CI of mean
  UB.99 <- Mean + qt(p=.995, df=N-1) * SEM
  CI.99 <- cbind(LB.99, UB.99)

  out <- rbind(CI.90, CI.95, CI.99)
  colnames(out) <- c("LB", "UB")
  rownames(out) <- c("90% CI", "95% CI", "99% CI")

  # returns a list  
  list(Confidence.Interval.of.Mean.using.t.distribution=out, 
       Descriptive.Statistics=cbind(N, Mean, SD, SEM, Min, Max))

} # end fun definition


# demo the function  ----------------------------------------

# library(ISwR); data(IgM)   # N=298
# classical.CI.mean(IgM)

# y <- rnorm(10, 50, 10)
# classical.CI.mean(y)



boot.CI.mean <- function(y, NS=1e3, data.name=""){   

  # ----------------------------------------------------------------
  #  ARGUMENTS
  #
  #  y          numeric vector, missing values OK
  #  NS         number of bootstrap samples
  #  data.name  decriptive string, e.g., "IgM, g/L"
  #
  #  A simple, percentile bootstrap procedure using in base R,
  #  replicate() as iterator. Returns a graph of the bootstrap 
  #  distribution of the mean, annotated with all results.
  #   
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove missing values
  if(nchar(data.name) > 27) data.name <- paste(substr(data.name,1,24),"...",sep="")

  Mean <- mean(y); N <- length(y); SD <- sd(y); SEM <- SD/sqrt(N) 

  time <- proc.time()
    out <- replicate(NS, mean( sample(y, replace=TRUE) ))
  time <- ((proc.time() - time)/60)[[3]]  # total time in minutes

  # --- percentile confidence intervals         # 95% CI of mean
  LB.95 <- quantile(out, .025)
  UB.95 <- quantile(out, .975)
  CI.95 <- c(LB.95, UB.95)

  # --- traditional CI based on normal theory   # 95% CI of mean
  LB.95t <- Mean - qt(p=.975, df=N-1) * SEM
  UB.95t <- Mean + qt(p=.975, df=N-1) * SEM

  hist(out, col="lightblue", border="white",
    main="Bootstrap distribution of the means",
    xlab="Mean", sub="red line: observed mean",
    col.sub="red", cex.sub=.8, font.lab=2)
  abline(v=Mean, col="red", lwd=2)
  abline(v=CI.95)

  mtext(data.name, side=3, adj=0, cex=.9)

  mtext("percentile bootstrap 95% CI of mean", side=1, line=2, adj=0, cex=.8)   
  mtext(paste("LB: ", signif(LB.95,6)), side=1, line=3, adj=0, cex=.8)
  mtext(paste("UB: ", signif(UB.95,6)), side=1, line=4, adj=0, cex=.8)
  
  mtext(paste("NS: ", NS, "resamples"), side=4, line=0, adj=0, cex=.8)
  mtext(paste("minutes: ", signif(time,6)), side=4, line=0, adj=1, cex=.8)
  
  mtext("95% CI by normal theory", side=1, line=2, adj=1, cex=.8)
  mtext(paste("LB: ", signif(LB.95t,6)), side=1, line=3, adj=1, cex=.8)
  mtext(paste("UB: ", signif(UB.95t,6)), side=1, line=4, adj=1, cex=.8)

  return(CI.95)

} # end function definition


# demo the function  ----------------------------------------

# library(ISwR); data(IgM)   # N=298
# boot.CI.mean(IgM, NS=1e4, data.name="IgM, g/L")
# boot.CI.mean(IgM, NS=1e5)

# boot.CI.mean(rnorm(1e3))
# boot.CI.mean(rnorm(1e3), NS=1e4)



my.stats <- function(y){

  # ----------------------------------------------------------------  
  #  ARGUMENTS
  #
  #  y         numeric vector, missing values OK
  #
  #  Function RETURNS a routine batch of stats in a column
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]  # remove NA, missing values

  skw <- function(y){                # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^3) / length(y)) / sd(y)^3}

  krt <- function(y){                # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^4) / length(y)) / sd(y)^4 - 3}

  # create a 10 x 1 matrix
  A <- rbind(length(y), mean(y), sd(y), 
       min(y), max(y), skw(y), krt(y), median(y),
       sd(y)/sqrt(length(y)), IQR(y) )

  colnames(A) <- "descriptive stats"

  rownames(A) <-c("N", "Mean", "SD", "Min", "Max", 
                   "skew", "kurtosis", "Median", "SEM", "IQR")
  return(A)  

} # end fun definition


# demo the function  ----------------------------------------

# my.stats(rnorm(10))         # 10 from N(0,1)
# my.stats(rnorm(1e3))        # 1000 from N(0,1)
# my.stats(rnorm(1e6, 50, 5)) # a million from N(50,5)

# library(ISwR); data(IgM)   # N=298
# my.stats(IgM)




new.stats <- function(y, roundoff=3, data.name="", norm.test=FALSE ){

  # ----------------------------------------------------------------
  #  ARGUMENTS
  #
  #  y          numeric vector, missing values OK  
  #  roundoff   decimal places, default 3
  #  data.name  short, descriptive string   
  #  norm.test  logical, default FALSE
  #
  #  requires packages:  DescTools, moments, nortest
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                     # remove missing values
  if(length(y) < 8) norm.test=FALSE     # don't test for normality on small sample

  skw <- function(y){                # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^3) / length(y)) / sd(y)^3}

  krt <- function(y){                # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^4) / length(y)) / sd(y)^4 - 3}

  yy <- y + rnorm(length(y), 0, .00001) # white noise to break ties in ks.test()

  require(DescTools); require(moments); require(nortest) # GOF tests
  N <- length(y); Mean <- mean(y); SEM <- sd(y)/sqrt(N)
  LB.95 <- Mean - qt(p=.975, df=N-1) * SEM  # 95% CI of mean
  UB.95 <- Mean + qt(p=.975, df=N-1) * SEM

  A <- rbind(N, Mean, sd(y), 100*sd(y)/Mean, min(y), max(y), skw(y), krt(y), 
             SEM, LB.95, UB.95, median(y), quantile(y, .25), quantile(y, .75), IQR(y) )
  A <- round(A, roundoff)
  if(data.name != "") colnames(A) <- data.name else colnames(A) <- ""
  rownames(A) <-c("N", "Mean", "SD", "CV%", "Min", "Max", "skew", "kurtosis", 
             "SEM", "LB.95%CI.mean", "UB.95%CI.mean", "Median", "25th", "75th", "IQR" )

  # Kincaid's revision of cvm.test() in library(nortest) to kill warnings about p, etc.
  my.cvm.test <- function (x){    # Cramer-von Mises test of normality  >?cvm.test
    x <- sort(x[complete.cases(x)])
    n <- length(x); if (n < 8) stop("sample size must be greater than 7")
    p <- pnorm((x-mean(x))/sd(x)); W <- (1/(12*n) + sum((p-(2 * seq(1:n)-1)/(2*n))^2))
    WW <- (1+0.5/n)*W
    if (WW < 0.0275) { pval <- 1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)}
       else if (WW < 0.051) { pval <- 1 - exp(-5.903 + 179.546 * WW - 1515.29 * WW^2)}
       else if (WW < 0.092) { pval <- exp(0.886 - 31.62 * WW + 10.897 * WW^2)}
       else if (WW < 1.1)   { pval <- exp(1.111 - 34.242 * WW + 12.832 * WW^2)}
       else { pval <- 7.37e-10 } #warning("p-value < 7.37e-10 cannot be computed")
    return(pval)
  }  # end of my.cvm.test() function definition

  if(norm.test == TRUE){

    if(N <= 5000 & N > 7){
      B <- rbind(shapiro.test(y)$p.value, jarque.test(y)$p.value, ad.test(y)$p.value,
             my.cvm.test(y), sf.test(y)$p.value, pearson.test(y)$p.value,
             agostino.test(y)$p.value, anscombe.test(y)$p.value,
             bonett.test(y)$p.value, ks.test(yy, "pnorm", mean(yy), sd(yy))$p.value,
             lillie.test(y)$p.value)

      colnames(B) <- "p-value"
      rownames(B) <- c("Shapiro-Wilk test","Robust Jarque-Berra test","Anderson-Darling test",
                    "Cramer-von Mises test", "Shapiro-Francia test", "Pearson test",
                    "D'Agostino test", "Anscombe-Glynn test","Bonett-Seier test", 
                    "Kolmogorov-Smirnov test", "Lillifors Kolmogorov-Smirnov test") 
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
                      "Kolmogorov-Smirnov test", "Lillifors Kolmogorov-Smirnov test") 
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
                      "Kolmogorov-Smirnov test", "Lillifors Kolmogorov-Smirnov test") 
     return(list(statistics=A, goodness.of.fit.tests.for.normality=B)) }
   }
 
   else

   return(A)   # return just the stats

}  # end function definition



# demo the function  ----------------------------------------

# library(ISwR); data(IgM)
# new.stats(IgM)        # minimal call
# new.stats(IgM, data.name="IgM in children")
# new.stats(IgM, roundoff=4, data.name="IgM in children")
# new.stats(IgM, roundoff=2, norm.test=TRUE)
# new.stats(IgM, roundoff=2, data.name="IgM in children", norm.test=TRUE)

# new.stats( rnorm(1e5) )
# new.stats( rnorm(1e5), roundoff=5, norm.test=TRUE )
# new.stats( rnorm(1e6), roundoff=2, norm.test=TRUE )



normal.QQ.plot <- function(y=rnorm(15)){

  # --------------------------------------------------------------------
  #  ARGUMENTS  graphical assessment of normality of a sample, prevails
  #
  #  y          numeric vector, missing values OK  
  #
  #  requires package: DAAG  for  qreference()
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # --------------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)
  y <- y[!is.na(y)]  # remove NA, missing values
  fun.name <- "normal.QQ.plot"

  y <- sort(y); n <- length(y); m <- mean(y); s <- sd(y)
  quantiles <- ((1:n) - .5) / n      # quantile of each obs value
  obs.z.score <- (y - m) / s         # for graphical convenience
  normal.z.score <- qnorm(quantiles) # normal distribution scores

  skw <- function(y){                # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^3) / length(y)) / sd(y)^3}

  krt <- function(y){                # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^4) / length(y)) / sd(y)^4 - 3}

  my.pause <- function(){
    if(interactive()) readline("Hit <Enter> to continue...");invisible()}

  thin.line <- function(N=25) cat(rep("-", N), sep="", "\n\n")  # prints "-"

  plot(normal.z.score, obs.z.score, 
    ylab="Observed z-scores", xlab="Normal z-scores", 
    cex.lab=1.2, font.lab=2, 
    main="Normal QQ plot with linear fit", cex.main=1.8
  )

  grid(); abline(lm(obs.z.score ~ normal.z.score), lwd=2)

  if(n <= 5000){ 
    p <- signif(shapiro.test(y)$p.value, 4)
    mtext(paste("Shapiro-Wilk test of normality, p =", p), cex=.9)
  }

  mtext(paste("skew =", round(skw(y), 3)), cex=.8, line=-1 )
  mtext(paste("kurtosis =", round(krt(y), 3)), cex=.8, line=-2 )
  mtext(paste("", " n =", n), cex=.8, line=-1, adj=0)
  mtext(paste("", " mean =", round(m, 4)), cex=.8, line=-2, adj=0)
  mtext(paste("", " SD =", round(s, 4)), cex=.8, line=-3, adj=0)

  mtext("z = (y - mean)/SD", side=1, line=4, adj=0, cex=.8)
  mtext("z-scores are the data", side=1, line=2, adj=1, cex=.8)
  mtext("converted into SD", side=1, line=3, adj=1, cex=.8)
  mtext("units above & below the mean", side=1, adj=1, line=4, cex=.8)

  cat("\n\n"); thin.line(70)
  cat("OUTPUT from the R function: ", fun.name, " by Dwight Kincaid, PhD\n\n")
  cat("The Call: "); print(match.call()); cat("\n")

  cat("The sample's stats:\n")
  cat("\tn =", n, ", mean =", m, ", SD =", s, "\n")
  cat("\tskew =", round(skw(y), 4), ", kurtosis =", round(krt(y), 4), "\n\n")
  if(n <= 5000) cat(paste("Shapiro-Wilk test of normality, p =", p, "\n\n"))
  cat("See the graphs.\n")

  my.pause() # --------------------
  hist(y, col="gray", main="Histogram", xlab="Response", cex.lab=1.2, font.lab=2)

  my.pause()
  dens <- density(y, n=2^10 ) 
  plot(dens, main="Kernel density estimate", 
    lwd=.5, font.lab=2, cex.lab=1.2)
  
  polygon(dens, col="lightblue") # color-in the kernel density

  my.pause() # --------------------
  library(DAAG)
  qreference(y, cex.strip=1, xlab="z-scores", ylab="" )
  title(
    main="Assessment of normality by qreference()\nin package DAAG by John Maindonald",
    cex.main=1.1, col.main="red", 
  sub="blue: observed data,  y-axis: response,  x-axis: theoretical quantiles", 
    col.sub="blue", cex.sub=.8, font.sub=2 )

  my.pause() # --------------------
  yy <- seq(min(y), max(y), length.out=500)
  norm.fit <- dnorm(yy, mean(y), sd(y))
  ylim <- range(c(0, norm.fit, hist(y, plot=F)$density ))

  hist(y, freq=F, ylim=ylim, xlab="Response", col="gray", cex.lab=1.2, font.lab=2,
    main="Density histogram with normal fit")
  points( yy, norm.fit, type="l", lwd=2.5)

  cat("Time stamp:", date(), "\n"); thin.line(70)

} # end function definition


# demo the function -----------------------------------


# normal.QQ.plot()  # minimal call

# normal.QQ.plot(IgM)  # library(ISwR); data)IgM)

# y <- rnorm(1e3)
# normal.QQ.plot(y)

# normal.QQ.plot(rweibull(50, 2.5, 1)) # sample from the Weibull distribution






"%mc.skew%" <- function(y, NS){

  # ----------------------------------------------------------------
  #  ARGS    'y'  numeric vector    'NS'  the number of samples
  #                                       simulated by Monte Carlo
  #
  #  binary operator for one-sided, Monte Carlo simulation 
  #  significance test for SKEW, achieved by random sampling from
  #  the normal distribution of the same N as y
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove missing values

  skw <- function(y){                # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^3) / length(y)) / sd(y)^3}

  N <- length(y); M <- mean(y); SD <- sd(y); obs.skew <- skw(y)

  time <- proc.time()
     out <- replicate(NS, skw(rnorm(N)))
  time <- ((proc.time() - time)/60)[[3]]  # total time in minutes

  ifelse(obs.skew > 0, NGE <- length(out[out >= obs.skew]),
     NGE <- length(out[out <= obs.skew]))

  p <- (NGE + 1) / (NS + 1)  # the standard formula, one-sided

  # graph the null distribution of skew under MC simulation test
  hist(out, col="gray", xlab="Skewness", main=NULL,
    font.lab=2, cex.lab=1.2, xlim=range(out, obs.skew) )
  
  mtext("Distribution of SKEW in samples from normal distribution", line=2, cex=1.2) 
  mtext("of same sample size as observed data", cex=1.2, line=1) 
  abline(v=obs.skew, col="red", lwd=2)
  mtext("red line: observed skew", side=3, line=0, col="red", cex=.9)
  mtext(paste("NS =", NS, "Monte Carlo simulatios"), side=1, adj=1, 
    line=4, col="red", cex=.8)
  mtext(paste("observed skew =", round(obs.skew,4)), 
    side=1, adj=1, line=3, col="red", cex=.8)
  mtext("Monte Carlo simulation test of skewness", side=4, line=-0, cex=1)
  mtext(paste("p =", signif(p, 5), "(1)"), side=1, line=3, adj=0, cex=.8)    
  mtext(paste("NGE =", NGE), side=1, line=4, adj=0, cex=.8)    

  # function returns a list to R Console

  list("observed.skew"=obs.skew, "NS"=NS, "NGE"=NGE, 
     "p.one.sided"=p, "minutes"=time)

} # end binary operator definition


# demo binary operator  --------------------------------

# library(ISwR); data(IgM)   # N=298
# IgM%mc.skew%1e4

# rnorm(50)%mc.skew%1e4




"%mc.kurt%" <- function(y, NS){

  # ----------------------------------------------------------------
  #  ARGS    'y'  numeric vector    'NS'  the number of samples
  #                                       simulated by Monte Carlo
  #  
  #  binary operator for one-sided, Monte Carlo simulation 
  #  significance test for KURTOSIS, achieved by random sampling from
  #  the normal distribution of the same N as y
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove missing values

  krt <- function(y){                # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^4) / length(y)) / sd(y)^4 - 3}

  N <- length(y); M <- mean(y); SD <- sd(y); obs.kurt <- krt(y)

  time <- proc.time()
     out <- replicate(NS, krt(rnorm(N)))
  time <- ((proc.time() - time)/60)[[3]]  # total time in minutes

  ifelse(obs.kurt > 0, NGE <- length(out[out >= obs.kurt]),
     NGE <- length(out[out <= obs.kurt]))

  p <- (NGE + 1) / (NS + 1)  # the standard formula, one-sided

  # graph the null distribution of kurtosis under MC simulation test
  hist(out, col="gray", xlab="Kurtosis", main=NULL,
    font.lab=2, cex.lab=1.2, xlim=range(out, obs.kurt) )
  
  mtext("Distribution of KURTOSIS in samples from normal distribution", 
    line=2, cex=1.2) 
  mtext("of same sample size as observed data", cex=1.2, line=1) 
  abline(v=obs.kurt, col="red", lwd=2)
  mtext("red line: observed kurtosis", side=3, line=0, col="red", cex=.9)
  mtext(paste("NS =", NS, "Monte Carlo simulatios"), side=1, adj=1, 
    line=4, col="red", cex=.8)
  mtext(paste("observed kurtosis =", round(obs.kurt,4)), 
    side=1, adj=1, line=3, col="red", cex=.8)
  mtext("Monte Carlo simulation test of kurtosis", side=4, line=-0, cex=1)
  mtext(paste("p =", signif(p, 5), "(1)"), side=1, line=3, adj=0, cex=.8)    
  mtext(paste("NGE =", NGE), side=1, line=4, adj=0, cex=.8)    

  list("observed.kurtosis"=obs.kurt, "NS"=NS, "NGE"=NGE, 
     "p.one.sided"=p, "minutes"=time)

} # end binary operator definition


# demo binary operator  ----------------------

# library(ISwR); data(IgM)   # N=298
# IgM%mc.kurt%1e4
# 

# rnorm(20)%mc.kurt%1e4





"%boot.skew%" <- function(y, NS){

  # ---------------------------------------------------------------------
  #  ARGS    'y'  numeric vector    'NS'  number of bootstrap resamples
  #
  #  binary operator for a percentile bootstrap 95% CI of SKEW
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # ---------------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove missing values

  skw <- function(y){                # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^3) / length(y)) / sd(y)^3}

  N <- length(y); M <- mean(y); SD <- sd(y); obs.skew <- skw(y)

  time <- proc.time()
     out <- replicate(NS, skw( sample(y, replace=TRUE) ))
  time <- ((proc.time() - time)/60)[[3]]  # total time in minutes

  # --- bootstrap percentile confidence interval, 95%
  LB.95 <- quantile(out, .025); UB.95 <- quantile(out, .975)

  # graph the bootstrap distribution of skew 
  hist(out, col="gray", xlab="Skewness", main=NULL,
    font.lab=2, cex.lab=1.2, xlim=range(out, obs.skew) )
  
  mtext("Bootstrap Distribution of SKEW", line=2, font=2, cex=1.2) 
  mtext("Percentile Bootstrap Confidence Intervals", line=1, cex=.9) 
  abline(v=obs.skew, col="red", lwd=2)
  mtext("red line: observed skew", side=3, line=0, col="red", cex=.9)
  mtext(paste("NS =", NS, "Bootstrap resamples"), side=1, adj=1, 
    line=4, col="red", cex=.8)
  mtext(paste("observed skew =", round(obs.skew,4)), 
    side=1, adj=1, line=3, col="red", cex=.8)
  mtext("Percentile Bootstrap, 95% CI", 
    side=1, adj=0, line=2, cex=.8)
  mtext(paste("LB =", signif(LB.95, 5)), side=1, line=3, adj=0, cex=.8)    
  mtext(paste("UB =", signif(UB.95, 5)), side=1, line=4, adj=0, cex=.8)   

  list("observed.skew"=obs.skew, "NS"=NS, "LB.95"=LB.95, "UB.95"=UB.95, 
     "skew.of.boot.distribution.of.skew"=skw(out),
     "minutes"=time)

} # end binary operator definition


# demo binary operator  ---------------------------

# library(ISwR); data(IgM)   #N=298
# IgM%boot.skew%1e4
# IgM%boot.skew%1e5


# rnorm(20)%boot.skew%1e4




"%boot.kurt%" <- function(y, NS){

  # -----------------------------------------------------------------------
  #  ARGS    'y'  numeric vector    'NS'  number of bootstrap resamples
  #
  #  binary operator for a percentile bootstrap 95% CI of KURTOSIS
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # -----------------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove missing values 

  krt <- function(y){                # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^4) / length(y)) / sd(y)^4 - 3}

  N <- length(y); M <- mean(y); SD <- sd(y); obs.kurt <- krt(y)

  time <- proc.time()
     out <- replicate(NS, krt( sample(y, replace=TRUE) ))
  time <- ((proc.time() - time)/60)[[3]]  # total time in minutes

  # --- bootstrap percentile confidence interval, 95%
  LB.95 <- quantile(out, .025); UB.95 <- quantile(out, .975)

  # graph the bootstrap distribution of kurtosis 
  hist(out, col="gray", xlab="Kurtosis", main=NULL,
    font.lab=2, cex.lab=1.2, xlim=range(out, obs.kurt) )
  
  mtext("Bootstrap Distribution of KURTOSIS", line=2, font=2, cex=1.2) 
  mtext("Percentile Bootstrap Confidence Intervals", line=1, cex=.9) 
  abline(v=obs.kurt, col="red", lwd=2)
  mtext("red line: observed kurtosis", side=3, line=0, col="red", cex=.9)
  mtext(paste("NS =", NS, "Bootstrap resamples"), side=1, adj=1, 
    line=4, col="red", cex=.8)
  mtext(paste("observed kurtosis =", round(obs.kurt,4)), 
    side=1, adj=1, line=3, col="red", cex=.8)
  mtext("Percentile Bootstrap, 95% CI", 
    side=1, adj=0, line=2, cex=.8)
  mtext(paste("LB =", signif(LB.95, 5)), side=1, line=3, adj=0, cex=.8)    
  mtext(paste("UB =", signif(UB.95, 5)), side=1, line=4, adj=0, cex=.8)   

  list("observed.kurt"=obs.kurt, "NS"=NS, "LB.95"=LB.95, "UB.95"=UB.95, 
     "skew.of.boot.distribution.of.kurtosis"=krt(out),
     "minutes"=time)

} # end binary operator definition


# demo binary operator  ---------------------------

# library(ISwR); data(IgM)   # N=298
# IgM%boot.kurt%1e4
# IgM%boot.kurt%1e5

# rnorm(20)%boot.kurt%1e4
# rnorm(1e3)%boot.kurt%1e3





Tukey.outliers <- function(y){

  # ----------------------------------------------------------------
  #  ARGS     Tukey's rule:  points beyond 1.5 IQR beyond a hinge 
  #                                        are potential OUTLIERS
  #  y    numeric vector, missing values OK
  #
  #  Returns a list of the trimmed vector and
  #  matrix of row numbers of outliers and the values
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # ----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove missing values 

  outlier.values <- boxplot(y, plot=FALSE)$out   # outlier values
  outlier.rows   <- which(y %in% outlier.values) # rows of outliers
  trimmed <- y[-outlier.rows]                    # delete outliers
  list( trimmed.vector = trimmed,    
        outliers = cbind(outlier.rows, outlier.values))

}  # end fun definition


# demo the function  ----------------------------------------

# library(ISwR); data(IgM)   # N=298
# Tukey.outliers(IgM)
#
# out <- Tukey.outliers(IgM)
# str(out)
# new.IgM <- out$trimmed.vector  # create new vector devoid of Tukey outliers
# hist(new.IgM)

# new.IgM%mc.skew%1e4    # compare to original data,    IgM%mc.skew%1e4
# new.IgM%mc.kurt%1e4    # compare to original data,    IgM%mc.kurt%1e4
# shapiro.test(new.IgM); DAAG::qreference(new.IgM)
# normal.QQ.plot(new.IgM)




"%bootMedian95CI%" <- function(y, NS){   

  # -----------------------------------------------------------------------
  #  ARGS    'y'  numeric vector    'NS'  number of bootstrap resamples
  #
  #  binary operator for a percentile bootstrap 95% CI of MEDIAN
  #
  #  Needs histogram of bootstrap distribution of the median.
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # -----------------------------------------------------------------------

  stopifnot(is.numeric(y), length(y) > 1)  # argument checking
  y <- y[!is.na(y)]                        # remove missing values
  quantile(replicate(NS, median(sample(y, replace=TRUE))), c(.025, .5, .975) )

} # end fun definition


# demo the function  ----------------------------------------

# library(ISwR); data(IgM)
# IgM%bootMedian95CI%1e4
# y <- rnorm(1e4); y%bootMedian95CI%1e3
# median(y)




skw <- function(y){                # SKEWNESS, type 3, fast Kincaid code
  y <- y - sum(y)/length(y)
  (sum(y^3) / length(y)) / sd(y)^3}

krt <- function(y){                # KURTOSIS, type 3, fast Kincaid code
  y <- y - sum(y)/length(y)
  (sum(y^4) / length(y)) / sd(y)^4 - 3}


# skw(IgM)
# krt(IgM)

# skw( rnorm(1e6) )
# krt( rnorm(1e6) )



################################################
#  definitions for some convenience functions  #
#  the names of which are self-explanatory     #
################################################

thick.line <- function(N=50) cat(rep("=", N), sep="", "\n")  # prints "="

thin.line <- function(N=50) cat(rep("-", N), sep="", "\n")  # prints "-"

my.pause <- function(){if(interactive()) readline("Hit <Enter> to continue...");invisible()}


mytick <- function (nx = 2, ny = 2, tick.ratio = 0.5){
  # a shameless, mere capture of the fun definition for minor.tick() from 
  # library(Hmisc) because some students have issues with Hmisc on old machines

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


# demo the function  ----------------------------------------

# hist(rnorm(1e4))
# mytick(2, 2)

# mytick(4, 4)
# mytick(2, 2, tick.ratio=1)

# hist(IgM); mytick(2, 2)

# ------------------------------
# END of FUNCTION DEFINITIONS    for single numeric samples
# ------------------------------

# end. Kincaid funs UNIVARIATE.R