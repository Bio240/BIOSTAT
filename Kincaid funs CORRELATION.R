# Kincaid funs CORRELATION.R    May 19, 2018  by Dwight Kincaid, PhD  

#library(MASS); library(hexbin);  library(boot); library(pwr); library(DAAG)
#library(moments); library(psych); library(nortest); library(MVN); library(ellipse) 

##
## DEMO x,y data suitable for demos of Kincaid's correlation functions
##
# sim data
# my.xy <- bivar.rnorm(rho=.6, N=20, txt.out=FALSE, graph=FALSE, empirical=TRUE)
# x <- my.xy[, "x"] 
# y <- my.xy[, "y"] 
#
##
## DEMO x,y data explainedat end of this file. Testosterone vs. Age, N=15,  
#    testosterone level, ng/dL
## variable names:   'testosterone'    'age'
##
# testosterone <- c(1305,1000,1175,1495,1060,800,1005,710,1150,605,690,700,625,610,450)
#
# age <- c(11, 12, 13, 14, 15, 16, 16, 17, 18, 20, 21, 23, 24, 27, 30)
#
#
# ---------------------------------------------------------------------------------
#  
#  DEMO FUNCTION CALLS listed after function definitions should be studied and run
#
# ---------------------------------------------------------------------------------
#
#     NEW FUNCTIONS                                                          
#
#  1. bivar.rnorm(rho, N, Mx, My,    # random x,y samples from the bivariate normal
#     SDx, SDy, rndf, txt.out, graph, ...)
#    
#  2. bivar.graphs(x, y, mvn.tests)  # IN PROGRESS; 2D/3D graphs; bivar normality tests
#
#  3. bootCIr(x, y, NS, CPUcore, OS, ...)  # bootstrap CI for Pearson r; percentile & BCA
#
#  4. boot.pwr.of.test(x, y, NS, traditional)   # bootstrap power of the test
#                                                 two-sided, for Pearson r
#
#  5. bootpowercurve.r(x, y, NS, low.n, high.n, by.n)   # bootstrap power curve
#  
#  6. new.stats(y, digits, ylab, norm.test)  # univariate stats & normality tests
#
#  7. perm.r(x, y, NS)                       # permutation test of Pearson r
#
#  8. confidence.ellipse(x, y, prob, pts,  # confidence ellipses by bivar. normal dist.
#         ell.col, ell.lwd, time.stamp, time.name, margin.txt, graph, return.mat,...)
#
#  9. leave.one.out.cor(x, y, graph=FALSE) # each x,y point deleted in turn &   
#                                          # influence on correlation measured
#
# --  other functions --
#
#  my.pause()      idle()    thick.line(N=25)    thin.line(N=25) 
#
# -- ideas on functions to add --
#
# xx. bivar.norm.tests( )        xx. perm.Spearman.r( )
# xx. %PearsonRperm%             xx. %SpearmanRperm%
# xx. bootCIr.ztrans( )          xx. 2D.kernel( )
# xx. bayesian.correlation( )    xx. bootpowercurve.Spearman.r( )
# xx. traditional.powercurve.r() xx. my.tick( )             
#
#


bivar.rnorm <- function(rho=runif(1, -.999, .999), 
  N=sample(c(10:40, seq(1e3, 1e4, 1e3)), size=1), Mx=50, My=100, 
  SDx=5, SDy=10, rndf=2, txt.out=TRUE, graph=TRUE, ...){

  #  --------------------------------------------------------------------
  #  RANDOM SAMPLING from the bivariate normal probability distribution
  #
  #  A convenience wrapper for mvrnorm() in library(MASS)
  #  Requires packages: MASS, MVN, hexbin
  #
  #  ARGUMENTS
  #
  #  rho       parametric correlation; default, rho=.8
  #  N         desired sample size for random x,y from 
  #            bivariate normal distribution
  #  Mx        parametric mean of x variable; default Mx=50
  #  My        parametric mean of y variable; default My=100
  #  SDx       parametric SD   of x variable; default SDx=5
  #  SDy       parametric SD   of y variable; default SDy=10
  #  rndf      decimal places for x,y round-off; default rndf=2
  #  txt.out   logical; text to R Console; default TRUE
  #  graph     logical; annotated scatter plot; default TRUE
  #  ...       arguments passed to mvrnorm() e.g., 'empirical=TRUE'
  #
  #  Function RETURNS x,y in data frame; can be called as arg free demo
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  #  --------------------------------------------------------------------

  my.pause <- function(){
    if(interactive()) readline("Hit <Enter> to continue...")
    invisible()}

  require(MASS); require(hexbin); require(MVN)

  g <- mvrnorm(n=N, mu=c(0, 0), Sigma=matrix(c(1,rho,rho,1),2,2), ... )
  g[,1] <- (g[,1] * SDx) + Mx   # mvrnorm() in library(MASS)
  g[,2] <- (g[,2] * SDy) + My

  x <- round(g[,1], rndf)
  y <- round(g[,2], rndf)
  obs.cor <- cor(x, y)

  if(txt.out == TRUE){
    cat("PARAMETERS specified in function call\n")
    print(rbind("rho"=rho, "N"=N, "Mean x"=Mx, "SD x"=SDx,
      "Mean y"=My, "SD y"=SDy))

    cat("\nDescriptive stats for random x,y sample\n")
    print(rbind("obs.r"= obs.cor, "N"=N, "Mean x"=mean(x), "SD x"=sd(x), 
      "Mean y"=mean(y), "SD y"=sd(y)))
  }

  if(N < 1e3) {cex.N <- 1.2; col.N <- "blue"} 
    else 
  {cex.N <- .8; col.N <- "lightblue"}

  if(graph == TRUE)
    plot(x, y, pch=16, main="Random sample from bivariate normal distribution",
      cex.lab=1.2, font.lab=2, col=col.N, cex=cex.N,
      sub=paste("N =", N, "points,  observed Pearson r =", round(cor(x,y),3),
                ", ", "p =", signif(cor.test(x,y)$p.value, 5), "(2)")
  )

  data.frame(x, y)  # returns a data frame: x,y random sample from the bivar. normal

}  # end function definition


# demo the function

# 
# my.sim <- bivar.rnorm()  # minimal call, taking all default args
# 
# plot(my.sim[, "x"], my.sim[, "y"])
# my.sim <- bivar.rnorm(empirical=TRUE)  # caution in use of 'empirical=TRUE'
# my.sim <- bivar.rnorm(empirical=FALSE) # 'empirical=FALSE' is default in mvrnorm()
# my.sim <- bivar.rnorm(N=3e2)
# my.sim <- bivar.rnorm(rho=-.2, N=2e3)


# my.sim <- bivar.rnorm(rho=-.9, N=300, Mx=50, SDx=2, My=75, SDy=3)
# str(my.sim)
# y1 <- my.sim[, "x"]
# y2 <- my.sim[, "y"]
# plot(y1, y2)


# rho=0    # zero correlation
# my.sim <- bivar.rnorm(rho=0, N=1e3, Mx=0, My=0, SDx=1, SDy=1 )


# DEMO
#
# demo x,y data suitable for correlation demos of Kincaid and other funs
# my.xy <- bivar.rnorm(rho=.6, N=200, txt.out=FALSE, graph=FALSE, empirical=TRUE)
# x <- my.xy[, "x"] 
# y <- my.xy[, "y"]  

# require(MVN)
# mardiaTest(cbind(x,y))
# mvnPlot( mardiaTest(cbind(x,y)), type="contour")
# mtext("contour plot with mvnPlot() of a mardiaTest() object, library(MVN", cex=.9)
#
# mvnPlot( mardiaTest(cbind(x,y)), type="persp")
# mtext("mesh plot with mvnPlot() of a mardiaTest() object, library(MVN", cex=.9)
#
# mvnPlot( mardiaTest(cbind(x,y)), type="persp", default=FALSE, col="gray", shade=.05, 
#   axes=F, expand=.8, theta=-20, r=20, phi=40, box=F)
# mtext("mesh plot with mvnPlot() of a mardiaTest() object, library(MVN", cex=.9)
#
# mvnPlot( mardiaTest(cbind(x,y)), type="persp", default=F, shade=.5, border=F, 
#   box=F, expand=1, col="wheat")



bivar.graphs <- function(x, y, sim=FALSE){

##  NON-FUNCTIONAL development code -- do not call!

    plot(x, y, pch=1, main="Random sample from bivariate normal distribution",
      sub=paste("N =", N, "points, obs.r =", round(cor(x,y),3),
         ",", "p =", signif(cor.test(x,y)$p.value, 6), "(2)"))
    my.pause()
    out <- kde2d(x, y, n=80)                     # library(MASS)
    persp(out, xlab="x", ylab="y", zlab="Count") # base R
    my.pause()
    persp(out, shade=.1, expand=.8, theta=10, d=.1, r=200, zlim=c(0, max(out$z)), 
      xlab="x", ylab="y", zlab="Count")
    my.pause()
    persp(out, shade=.05, axes=FALSE, box=FALSE, expand=.8, theta=10, d=.1, r=200 )
    my.pause()
    image(out)
    my.pause()
    contour(out)                                 # base R
    my.pause()
    contour(out, axes=FALSE, col=adjustcolor(1:3, .5))
    my.pause()
    contour(out, nlevels=15, axes=TRUE, col="black", lwd=1.5,
      drawlabels=TRUE, xlab="x", ylab="y", vfont=c("sans serif", "plain"),
      labcex=1)    # Hershey font used!
    my.pause()

    if(N > 200){
      plot( hexbin(x, y), main="hexagonal binning" )   # library(hexbin)
      my.pause()
      plot( hexbin(x, y), colramp=BTY, main="hexagonal binning") }
    my.pause()

    mvnPlot( mardiaTest(cbind(x, y)), type="contour")    # library(MVN)
    mtext("contour plot with mvnPlot() of a mardiaTest() object, library MVN", cex=.9)
    my.pause()

    mvnPlot( mardiaTest(cbind(x, y)), type="persp")
    mtext("mesh plot with mvnPlot() of a mardiaTest() object, library MVN", cex=.9)
    my.pause()

    mvnPlot( mardiaTest(cbind(x, y)), type="persp", default=FALSE, col="gray",
      shade=.05, axes=FALSE, expand=.8, theta=-20, r=20, phi=40, box=FALSE)
    mtext("mesh plot with mvnPlot() of a mardiaTest() object, library MVN", cex=.9)
    my.pause()

    mvnPlot( mardiaTest(cbind(x, y)), type="persp", default=FALSE,
      shade=.5, border=FALSE, box=FALSE, expand=1, col="wheat")

} # end function definition

# demo the function

# bivar.graphs(age, testosterone, sim=FALSE)


bootCIr <- function(x, y, NS=1e3, CPUcore=FALSE, OS="Win", ...){

  #  ---------------------------------------------------------------------------
  #  Nonparametric bootstrap two-sided CI for Pearson r
  #
  #  ARGUMENTS    Requires library(boot) for percentile & BCA bootstrap methods
  #
  #  x        numeric, random variables; paired. Missing values OK.
  #  y        numeric, random variable, paired with 'x'
  #  NS       desired N of bootstrap samples
  #  CPUcore  logical, default FALSE; TRUE for multicore CPU parallel processing
  #  OS       character, either "Win" or "Mac"
  #  ...      args to pass to boot.ci() & abc.ci() 
  #             e.g., conf=c(.9, .95, .99), or conf=.99 ; .95 is default
  #
  #  RETURNS text to R Console & graph of bootstrap distribution of r.
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  #  ---------------------------------------------------------------------------

  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y), 
    is.numeric(NS), is.logical(CPUcore), is.character(OS))  # data & arg checking

  run.time <- proc.time()          # begin timing the entire run
  D <- na.omit(cbind(x, y))        # matrix, delete rows with missing values
  x  <- D[, "x"] ; y  <- D[, "y"]  # reconstitute x,y vectors
  cores <- parallel::detectCores() # returns N of CPU cores, Win or Mac

  require(boot); require(snow); require(parallel)
  thin.line <- function(N=30) cat(rep("-", N), sep="", "\n\n")  # prints "-"
  my.skew <- function(x) sum((x-sum(x)/length(x))^3)/(length(x)*sd(x)^3) # type 3

  N <- length(y); r <- cor(x, y); p <- cor.test(x, y)$p.value
  LB.trad.95.CI <- signif(cor.test(x, y)$conf.int[1], 3)
  UB.trad.95.CI <- signif(cor.test(x, y)$conf.int[2], 3)

  cat("\n\n"); thin.line(66)
  cat("OUTPUT from Kincaid's:  bootCIr()", "\t", date(), "\n\n")
  cat("OS:\t\t", sessionInfo()$running, "\n")
  cat("R:\t\t", sessionInfo()$R.version$version.string, " www.r-project.org\n")
  cat("Platform:\t", sessionInfo()$platform, "\n\n")

  cat("CALL:\n"); print(sys.call()); cat("\n")
  cat("CPU cores available:", cores, "\n")
  if(CPUcore == TRUE) cat("Note: PARALLEL, multicore processing in progress...\n\n") 
    else cat("Note: SEQUENTIAL, nonparallel processing in progress...\n\n")
  
  cat("STATISTIC:  Pearson correlation via cor(x, y)\n")
  cat(paste("NS =", NS, "bootstrap samples taken from observed x,y data."), "\n")
  cat("From library(boot) the functions boot(), boot.ci() and\n")
  cat("abc.ci() were called.", "\n\n")
  cat(paste("observed N:", N, "  observed Pearson r:", round(r,5)), "\n") 
  cat(paste("p-value: ", signif(p,4), "(2) using normal theory (Fisher's z, t-dist)",sep=""), "\n\n")
  flush.console()  # Console buffering

  # establish Windows CPU multicore CLUSTER if needed
  if(OS == "Win" & CPUcore == TRUE) cl=makeCluster(cores, type="SOCK") 

  # configure string to send to 'parallel' arg in boot()
  if(OS == "Win") para.type = "snow"
  if(OS == "Mac") para.type = "multicore"

  my.fun    <- function(data, indices){cor(data[indices,1], data[indices,2])}

  boot.time <- proc.time()

    if(CPUcore == FALSE) 
      out <- boot(data=D, statistic=my.fun, R=NS, parallel="no" )
    if(CPUcore == TRUE & OS == "Win")  
      out <- boot(data=D, statistic=my.fun, R=NS, parallel=para.type, ncpus=cores, cl=cl) 
    if(CPUcore == TRUE & OS == "Mac")  
      out <- boot(data=D, statistic=my.fun, R=NS, parallel=para.type) 

  boot.time <- proc.time() - boot.time  # seconds
  boot.time <- boot.time/60             # minutes
  boot.time <- round(boot.time[[3]], 4) # extract minutes then round

  if(OS == "Win" & CPUcore == TRUE) stopCluster(cl) # stop Windows parallel cluster

  boot.bias <- signif(mean(out$t) - cor(x,y), 3) # in the bootstrap distribution of r
  boot.SE   <- signif(sd(out$t), 3)              #     we seek small bias and
  boot.skew <- signif(my.skew(out$t), 3)         #         low skewness

  cat(paste("MINUTES for bootstrap resampling: ", boot.time), "\n\n"); flush.console()

  thin.line(5); print(out); cat("\n\n"); thin.line(5); flush.console()

  if(N < 1e3) boot.type=c("perc", "bca") else boot.type="perc" # else boot.ci() chokes

  ci.time  <- proc.time()
    ci.out <- boot.ci(out, type=boot.type, ...)  # library(boot)
  ci.time  <- proc.time() - ci.time
  ci.time  <- ci.time[[3]]  # seconds

  print(ci.out); cat("\n")
  cat("SECONDS to calculate CIs from bootstrap distribution:", ci.time, "\n\n")
  thin.line(5); flush.console()

  if(N < 5000){  # bca chokes at high n 
    cat("BCA bootstrap CIs by numerical differentiation after R. Tibshirani.\n\n")
    cat("Two-sided nonparametric approximate, BCA CIs for Pearson r\n")
    cat("using boot::abc.ci() ; for references see documentation >?abc.ci \n\n")

    abc.time  <- proc.time()
      abc.out <- abc.ci(D, boot::corr, ...)        # library(boot)
    abc.time  <- proc.time() - abc.time; abc.time <- abc.time[[3]]  # seconds

    cat("\t\tLB\t UB \n")
    print(abc.out)
    cat("\n"); cat("SECONDS to approximate BCA CIs of r by ")
    cat("numerical differentiation:", abc.time, "\n\n"); thin.line(5)
  }

  run.time <- (proc.time() - run.time)/60
  cat("MINUTES for entire run:", round(run.time[[3]], 3), " ending:",date(),"\n\n")
  cat("AUTHOR:\tDwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu\n")
  cat("\tContact me for suggestions on this wrapper function.\n")
  thin.line(66); flush.console()  # Console buffering

  hist(out$t, main=NULL, col="lightblue", xlab="Correlation", 
    font.lab=2, cex.lab=1.2, cex.main=1, font.main=4, col.main="blue")
  abline(v=r, col="red", lwd=2) # vertical red line at observed r

  mtext(paste("Bootstrap distribution of Pearson r achieved by", NS, "resamples"), 
    side=3, line=2.5, col="blue", cex=1, font=2)
  mtext(paste("Boot minutes: ", boot.time), side=1, line=4, adj=1, cex=.9)
  mtext("Red line:  observed r", side=1, line=4, adj=0, col="red", cex=.9)

  mtext(paste("observed r =", signif(r, 3), ", p =", signif(p, 4), "(2)"), 
    side=3, line=1.6, adj=0, cex=.85)
  mtext(paste("traditional CI 95%: ", LB.trad.95.CI, ",", UB.trad.95.CI),
    side=3, line= .8, adj=0, cex=.85)
  mtext(paste("N =", N, "x,y points"), side=3, line=0, adj=0, cex=.85)
  
  cx <- .6  # font size for mtext()
  mtext("STATS for Bootstrap Distribution", side=3, line=1.65, adj=1, cex=.75, font=2)
  mtext(paste("bootstrap mean, bias: ", signif(mean(out$t), 3), ",", boot.bias), 
    side=3, line=1, adj=1, cex=cx)
  mtext(paste( "range =", signif(range(out$t)[1], 3), " to ", 
    signif(range(out$t)[2], 3) ), side=3, line=.5, adj=1, cex=cx)
  mtext(paste("bootstrap SE =", boot.SE), side=3, line=0, adj=1, cex=cx)
  mtext(paste("skewness =", boot.skew), side=3, line=-.5, adj=1, cex=cx)

  mtext(paste("Percentile 95% CI:", signif(ci.out$percent[4], 4), ",", 
    signif(ci.out$percent[5], 4)), side=1, adj=0, line=2.2, cex=.7)

  if(N < 5000){ 
    mtext(paste("BCA 95% CI:", signif(abc.out[2], 4), ",", signif(abc.out[3], 4)),
      side=1, adj=1, line=2.2, cex=.7)} # by differentiation after Tibshirani

  if(N >= 1e3) warning(paste("\n\nCalculating BCA bootstrap when observed N >= 1e3 may",
    "\nthrow vector memory allocation error in boot::boot.ci() or simply",
    "\ntake too long. Here, view BCA results from boot::abc.ci()\n" ))

}  # end function definition


#  demo the function

##
## testosterone and age, N=15
##
#
# bootCIr(age, testosterone, NS=1e4)           # default is 95%CI of r
#   compare to normal theory
#     cor.test(age, testosterone) # 95% CI = -.936, -.517 (2)
#
# bootCIr(age, testosterone, NS=1e4, CPUcore=TRUE, OS="Win")
# bootCIr(age, testosterone, NS=1e5, CPUcore=TRUE, OS="Win", conf=c(.9, .95, .99))
#
#
# bootCIr(age, testosterone, NS=1e4, conf=.99) # 99%CI of r
# bootCIr(age, testosterone, NS=1e4, conf=.9)  # 90%CI of r
# bootCIr(age, testosterone, NS=1e4, conf=c(.9, .95, .99))
 

# sim data 

# bootCIr(rnorm(1e3), rnorm(1e3), NS=1e5)
# bootCIr(rnorm(1e3), rnorm(1e3), NS=1e5, CPUcore=TRUE, OS="Win")
#
# bootCIr(rnorm(2e3), rnorm(2e3), NS=1e3, conf=c(.9, .95, .99))



#
# change the BOOT object...useful... THIS  CODE  WORKS   << DEVELOPMENT CODE >>
#
#a <- out  # 'out' is object from a call to boot()
#new.n <- 1e3 # select the first 1e3 boot samples
#a[3] <- new.n ; str(a)
#a$t <- matrix(data=a$t[1:new.n], nrow=new.n, ncol=1)
#str(a)

#hist(a$t)
#a.ci.out <- boot.ci(a, type=c("perc", "bca"), conf=c(.9, .95, .99) )
#Level     Percentile            BCa          
#90%   (-0.9160, -0.6976 )   (-0.9027, -0.6656 )   
#95%   (-0.9303, -0.6639 )   (-0.9168, -0.6302 )   
#99%   (-0.9537, -0.5817 )   (-0.9433, -0.5565 )  
 
#aa <- out
#new.n <- 1e3 # select the next 1e3 boot samples
#aa[3] <- new.n ; str(aa)
#aa$t <- matrix(data=aa$t[10001:20000], nrow=new.n, ncol=1)
#str(aa)

#hist(aa$t)
#aa.ci.out <- boot.ci(aa, type=c("perc", "bca"), conf=c(.9, .95, .99) )
#Level     Percentile            BCa          
#90%   (-0.9164, -0.6988 )   (-0.9004, -0.6651 )   
#95%   (-0.9301, -0.6664 )   (-0.9153, -0.6299 )   
#99%   (-0.9535, -0.5986 )   (-0.9424, -0.5381 )  
# -------------------------------------------------------------

#jack.after.boot(out, index=1, alpha=c(.01, .05, 0, .95, .99) )
#empinf(out, type="jack" )  # where 'out' is a boot() object




boot.pwr.of.test <- function(x, y, NS=1e3, traditional=FALSE){

  #  ----------------------------------------------------------------------------
  #  Power of the test for Pearson r by bootstrapping observed x,y data
  #  Requires packages:  boot, pwr
  #
  #  ARGUMENTS
  #
  #  x            numeric, random variable
  #  y            numeric, random variable, paired with 'x'
  #  NS           number of desired bootstrap samples
  #  traditional  also delivers traditional power using normal theory (Cohen)
  #
  #  Function RETURNS text to R Console. 
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  #  ----------------------------------------------------------------------------

  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))

  D  <- na.omit(cbind(x, y))         # delete rows with missing values
  x  <- D[, "x"]  # reconstitute x
  y  <- D[, "y"]  # reconstitute y

  require(boot); require(pwr)
  thin.line <- function(N=25) cat(rep("-", N), sep="", "\n\n")  # prints "-"

  obs.r <- cor(y, x)
  obs.p <- cor.test(y, x)$p.value
  N     <- length(y)

  my.fun <- function(data, indices){cor.test(data[indices,1], data[indices,2])$p.value}
  #my.fun <- function(data, indices){cor.test(data$x, data$y)$p.value}

  boot.time   <- proc.time()  # a 'matrix' boots faster than a 'data frame' object
    out       <- boot(data=D, statistic=my.fun, R=NS, parallel="no") 
  boot.time   <- (proc.time() - boot.time)/60  # minutes

  p <- out$t  # extract p-value vector
  power.05  <- length(p[p <= .05]) /NS  # power
  power.01  <- length(p[p <= .01]) /NS
  power.02  <- length(p[p <= .02]) /NS
  power.005 <- length(p[p <= .005])/NS
  power.001 <- length(p[p <= .001])/NS

  cat("\n\n"); thin.line(62)
  cat(paste("Output from Kincaid's:  boot.pwr.of.test()\n\n\t", date()), "\n\n")
  cat(paste("NS =", NS, "bootstrap samples were taken."), "\n")
  cat("Statistic:\tp-value, two-sided, of Pearson correlation\n")
  cat("\t\tcor.test(y,x)$p.value\n\n")
  cat("From library(boot) the function boot() was called.", "\n\n")
  cat(paste("obs. N:", N, "\tobs. r:", round(obs.r,5), "\tobs. p:", signif(obs.p,4)), "\n\n")
  cat(paste("Minutes required for bootstrap: ", round(boot.time[[3]],4)), "\n\n")
  temp <- rbind(power.05, power.01, power.02, power.005, power.001)
  colnames(temp) <- "BOOTSTRAP power of the test"
  print(temp); cat("\n\n")

  if(traditional == TRUE){
    pwr.05  <- pwr.r.test(n=N, r=obs.r, sig.level=.05,  alt="two.sided", power=NULL)[[4]]
    pwr.01  <- pwr.r.test(n=N, r=obs.r, sig.level=.01,  alt="two.sided", power=NULL)[[4]]
    pwr.02  <- pwr.r.test(n=N, r=obs.r, sig.level=.02,  alt="two.sided", power=NULL)[[4]]
    pwr.005 <- pwr.r.test(n=N, r=obs.r, sig.level=.005, alt="two.sided", power=NULL)[[4]]
    pwr.001 <- pwr.r.test(n=N, r=obs.r, sig.level=.001, alt="two.sided", power=NULL)[[4]]
    temp <- rbind(pwr.05, pwr.01, pwr.02, pwr.005, pwr.001)
    cat("\n'traditional' power of the test, using", "\n")
    colnames(temp) <- "pwr.r.test() in library(pwr)"
    print(temp); cat("\n\n")}

  cat("Note: All power calculated from two-sided p-values.", "\n\n")
  cat("AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu\n")
  thin.line(62)

}  # end function definition


# demo the function

##
## testosterone and age, N=15  which is probably not enough data to warrant bootstrap
##
#
# boot.pwr.of.test(age, testosterone, NS=1e4)
# boot.pwr.of.test(age, testosterone, NS=1e4, traditional=TRUE)


##
## simulated data
##
#
# sim.df <- bivar.rnorm(rho=.8, N=30, graph=FALSE, txt.out=FALSE, empirical=TRUE)
# x <- sim.df[,1]
# y <- sim.df[,2]
#
# boot.pwr.of.test(x, y, NS=1e4)
# boot.pwr.of.test(x, y, NS=1e4, traditional=TRUE)




bootpowercurve.r <- function(x, y, NS=1e3, low.n, high.n, by.n){

  #  ----------------------------------------------------------------------
  #  bootstrap power curves for Pearson r
  #
  #  ARGUMENTS
  #
  #  x         numeric, random variable
  #  y         numeric, random variable, paired with 'x'
  #  NS        desired number of bootstrap samples, per sample size step
  #  low.n, high.n, by.n    self explanatory, prospective sample sizes
  #
  #  Function RETURNS text to R Console and graph. 
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  #  ----------------------------------------------------------------------

  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))

  d  <- na.omit(data.frame(x, y))  # delete rows with missing values
  x  <- d[, "x"]  # reconstitute x
  y  <- d[, "y"]  # reconstitute y

  N       <- length(y)
  new.n   <- seq(low.n, high.n, by.n)  # vector of prospective sample sizes
  obs.r   <- cor(x, y)
  obs.p   <- cor.test(x, y)$p.value
  power05 <- power01 <- power001 <- numeric(length(new.n)) # initialize vectors

  timer <- proc.time()     # should try boot() instead of replicate() for speed

  for(i in seq_along(new.n)){
    rows <- new.n[i]
    p <- replicate(NS, { 
                        my.rows <- sample(1:N, size=rows, replace=TRUE);
                        cor.test(x[my.rows], y[my.rows])$p.value
                         }
         ) # only need p-value!

    power05[i]  <- length( p[p <= .05]) /NS  
    power01[i]  <- length( p[p <= .01]) /NS
    power001[i] <- length( p[p <= .001])/NS
  }

  timer <- ((proc.time() - timer)/60)[3]  # minutes
 
  cat("\nTotal minutes:", round(timer,4), "\n\n")
  print(cbind("bootstrap.NS"=NS, new.n, power05, power01, power001)); cat("\n")

  # create empty graph to add to
  plot(c(min(new.n), max(new.n)), c(0,1), type="n", 
    xlab="Sample Size", ylab="Statistical Power", font.lab=2, cex.main=.9,
    main="Bootstrap power curves for Pearson correlation", cex.lab=1.2  )
  abline(h=.8, lty=2)  # dotted line at power=.8

  points(new.n, power05,  type="b", pch=16, col="red",   lwd=1)
  points(new.n, power01,  type="b", pch=1,  col="black", lwd=.5, cex=.8)
  points(new.n, power001, type="b", pch=1,  col="blue",  lwd=.5, cex=.8)

  mtext("red     .05 alpha",  side=1, line=2, adj=0, col="red", font=2) 
  mtext("black  .01",  side=1, line=3, adj=0, col="black", font=2) 
  mtext("blue  .001", side=1, line=4, adj=0, col="blue", font=2) 
  mtext(paste(NS,"bootstrap samples per sample size increment"),side=3,adj=0,cex=.8)
  mtext(paste("minutes =", round(timer,3)), side=1, line=4, adj=1, cex=.9)
  mtext("vertical line at obs. N", side=3, adj=1, cex=.8)
  mtext(paste("obs. p =", signif(obs.p, 4)), side=1, line=3, adj=1, cex=.9)
  mtext(paste("obs. r =", round(obs.r, 3)), side=1, line=2, adj=1, cex=.9)
  abline(v=N); grid()

}  # end function definition


# demo the function


##
## testosterone and age, N=15
##
#
# bootpowercurve.r(age, testosterone, NS=1e3, low.n=9, high.n=16, by.n=1)




new.stats <- function(y, digits=3, ylab=NULL, norm.test=TRUE ){

  # -----------------------------------------------------------------
  #  Stats and normality tests for single, numeric samples
  #
  #  ARGS     y: numeric vector     digits: roundoff  
  #        ylab: data name       norm.test: logical, normality tests   
  #
  #  requires packages:  moments, nortest -- both for GOF tests
  #
  #  AUTHOR: Dwight Kincaid, PhD    dwight.kincaid@lehman.cuny.edu
  # -----------------------------------------------------------------

  stopifnot(is.numeric(y), length(y)>1,            # arg & data checking
    is.numeric(digits), is.logical(norm.test) ) 

  y <- y[!is.na(y)]                     # remove missing values from numeric vector
  if(length(y) < 8) norm.test=FALSE     # don't test for normality on small sample
  yy <- y + rnorm(length(y), 0, .00001) # white noise to break ties in ks.test()

  skw <- function(y){                   # SKEWNESS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
    (sum(y^3) / length(y)) / sd(y)^3}

  krt <- function(y){                   # KURTOSIS, type 3, fast Kincaid code
    y <- y - sum(y)/length(y)
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


##
## testosterone and age, N=15
##
#
# The so-called 'marginal distributions' of x and y should be normally distributed for p
# and CI to be reliable in Pearson correlation. And while testing marginal distributions
# for normality is not a test for bivariate normality, it is an indicator of same.
#
# new.stats(testosterone, ylab="Testosterone level", norm.test=TRUE)
# new.stats(age, ylab="Age", norm.test=TRUE)




perm.r <- function(x, y, NS=1e3){

  #  -------------------------------------------------------------------------
  #  permutation test of Pearson r; one-sided and two-sided p-values
  #
  #  ARGUMENTS
  #
  #  x         numeric, random variable
  #  y         numeric, random variable, paired with 'x'
  #  NS        desired N of shuffles of one variable, holding the other fixed
  #
  #  Function RETURNS a fully annotated graph of the null distribution of r. 
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  #  --------------------------------------------------------------------------

  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y), is.numeric(NS))

  D  <- na.omit(cbind(x, y))  # delete rows with missing values
  x  <- D[, "x"]  # reconstitute x
  y  <- D[, "y"]  # reconstitute y

  obs.r <- cor(x, y)  # Pearson r

  perm.time <- proc.time()
    null.r  <- replicate(NS, cor(sample(y), x))  # embarassingly simple in R
  perm.time <- (proc.time() - perm.time)/60      # minutes

  # Tally null r, for N Greater than or Equal to observed r; one-sided & two-sided
  NGE1 <- NGE2 <- 0  

  if(obs.r >  0) NGE1 <- length(null.r[null.r >= obs.r])    # one-sided NGE
  if(obs.r <= 0) NGE1 <- length(null.r[null.r <= obs.r])    # one-sided NGE
                  
  abs.null <- abs(null.r)
  NGE2     <- length(abs.null[abs.null >= abs(obs.r)])  # two-sided NGE for any obs.r

  p.value1 <- (NGE1+1) / (NS+1)    # one-sided p-value
  p.value2 <- (NGE2+1) / (NS+1)    # two-sided "  "

  # histogram of distribution of null r, under permutation

  hist(null.r, main="Permutation, null distribution of r",
    xlab="Correlation", font.lab=2, col="lightgray",
    xlim=range(null.r, obs.r))

  abline(v=obs.r, col="red", lwd=2) # vertical line at observed r

  mtext("red line at observed Pearson r", adj=0, col="red", cex=.9)
  mtext(paste("observed r =", round(obs.r,5), " (Pearson)"), 
    side=1, adj=0, line=4, col="red", cex=.9)
  mtext(paste("permutation p=", signif(p.value1, 5), " (1)", sep=""), 
    side=1, adj=1, line=3, col="black", cex=.85)
  mtext(paste("permutation p=", signif(p.value2, 5), " (2)", sep=""), 
    side=1, adj=1, line=4, col="black", cex=.85)
  mtext(paste("NS =", NS, "shuffles"), adj=1, col="black", cex=.9)
  mtext(paste("CPU minutes =", round(perm.time[[3]], 4)), 
    side=1, adj=0, line=2, cex=.8)
  mtext(paste("NGE(1)= ", NGE1, " , ", "p= (NGE+1)/(NS+1)", sep=""), 
    side=1, adj=1, line=2, cex=.75)

} # end function definition


# demo the function


##
## testosterone and age, N=15
##
#
# perm.r(age, testosterone, NS=1e4)
# cor.test(age, testosterone, alternative="less") # compare to result by normal theory





confidence.ellipse <- function(x, y, prob=.95, pts=200, 
  ell.col="blue", ell.lwd=2, time.stamp=FALSE, time.name=NULL, 
  margin.txt=FALSE, graph=TRUE, return.mat=FALSE,...){

  #  ---------------------------------------------------------------------------------
  #  confidence ellipse around points in a scatter plot, based on bivariate normality 
  #
  #  ARGUMENTS   requires library 'ellipse'
  #
  #  x           numeric, random variable
  #  y           numeric, random variable, paired with 'x'
  #  prob        default is .95 for 95% probability ellipse around points
  #  pts         number of points around ellipse
  #  ell.col     color of ellipse
  #  ell.lwd     thickness of ellipse line
  #  time.stamp  logical
  #  time.name   string; name to paste to time stamp
  #  margin.txt  logical
  #  graph       logical
  #  return.mat  logical
  #  ...         args passed to plot()
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  #  ---------------------------------------------------------------------------------

  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))

  d  <- na.omit(data.frame(x, y))  # delete rows with missing values
  x  <- d[, "x"]                   # reconstitute x
  y  <- d[, "y"]                   # reconstitute y

  require(ellipse)

  elps <- ellipse( cor(x, y), npoints=pts, # ellipse() returns the x,y
    scale=c(sd(x), sd(y)), level=prob,     # matrix of 200 rows for
    centre=c(mean(x),mean(y)))             # column 'x' and column 'y'

  if(graph == TRUE){
    range.x <- range(x, elps[,1])   # useful for plotting to
    range.y <- range(y, elps[,2])   # expand axis scales for ellipse 

    N     <- length(x)
    obs.r <- round(cor(x, y), 4)
    p     <- signif(cor.test(x, y)$p.value, 4)  # p-value, 4 digits 

    plot(x, y, xlim=range.x, ylim=range.y, ...) # user may use plot() args in call

    lines(elps, type="l", col=ell.col, lwd=ell.lwd)

    if(margin.txt == TRUE){ 
      mtext(paste("Pearson r =", obs.r, ",", " N =", N,
        ",",  " p =", p, "(2)"), side=1, line=4, adj=0, cex=.7) 
      mtext(paste(prob*100, "% confidence ellipse of points", sep=""),
        side=1, line=4, adj=1, cex=.7)}   

    if(time.stamp == TRUE){ 
      mtext( paste(time.name, Sys.Date()), side=1, adj=1, line=3, 
        cex=.75, col="gray")}
  }

  if(return.mat == TRUE) return(elps) # returns matrix of x,y points around ellipse

}  # end of function definition


# demo the function

#
# testosterone and age data set
#
# confidence.ellipse(age, testosterone)
# confidence.ellipse(age, testosterone, time.stamp=TRUE, time.name="Cain")
# confidence.ellipse(x=age, y=testosterone)

# confidence.ellipse(age, testosterone, margin.txt=TRUE, prob=.9)

# confidence.ellipse(age, testosterone, pch=20, 
#    xlab="Age at First Conviction", ylab="Testosterone, ng/dL")

# confidence.ellipse(age, testosterone, pts=50, margin.txt=TRUE, pch=1, 
#   ell.col="red", ell.lwd=1, font.lab=2, cex.lab=1.2, cex.main=1.5,
#   xlab="Age at First Conviction", ylab="Testosterone, ng/dL",
#   main="95% confidence ellipse of points")


#
# send ellipse x,y coordinates matrix to an object for use in plotting
# your way, with absolute control!

# out <- confidence.ellipse(age, testosterone, pts=500, graph=FALSE, return.mat=TRUE)
# str(out)

# plot(age, testosterone)
# lines(out)



leave.one.out.cor <- function(x, y, graph=FALSE){

  #  -------------------------------------------------------------------
  #  each x,y point is left out in turn, and its influence on 
  #  correlation analysis reported as text
  #
  #  ARGUMENTS
  #
  #  x            numeric, random variable
  #  y            numeric, random variable, paired with 'x'
  #  graph        logical; NOT USED but to be added
  #  
  #  AUTHOR: Dwight Kincaid, PhD   dwight.kincaid@lehman.cuny.edu
  #  --------------------------------------------------------------------

  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  d   <- na.omit(data.frame(x, y))  # delete rows with missing values
  x   <- d[, "x"]  # reconstitute x
  y   <- d[, "y"]  # reconstitute y

  my.pause <- function(){
    if(interactive()) readline("Hit <Enter> to continue...")
    invisible()}
  thin.line <- function(N=25) cat(rep("-", N), sep="", "\n\n")  # prints "-"

  N        <- length(y)         # stats for observed data set
  cor.obj  <- cor.test(x, y)    #   for
  obs.r    <- cor.obj$estimate  #     observed data
  obs.Rsq  <- obs.r^2           #       set
  obs.p    <- cor.obj$p.value   # two-sided p-value

  r <- p <- LB95 <- UB95 <- numeric(N) # initialize vectors to populate in 'for' loop

  # 'for' loop omits each x,y data point in turn, by negative integer indexing

  for(i in 1:N){
     out     <- cor.test(x[-i], y[-i], method="pearson", 
                alternative="two.sided", conf.level=.95)

     r[i]    <- out$estimate     # Pearson r
     p[i]    <- out$p.value      # two-sided p-value
     LB95[i] <- out$conf.int[1]  # Lower Bound of 95% CI of r
     UB95[i] <- out$conf.int[2]  # Upper Bound "  "   "  "  "
  }

  Rsq      <- r^2           # vectorized R-square per omitted point
  Rsq.diff <- Rsq - obs.Rsq
  dfr      <- data.frame(point.rm=1:N, x=x, y=y, Pearson.r=r, 
                p=p, LB95CI= LB95, UB95CI=UB95, Rsq.diff, new.R.sq=Rsq) 
  dfr      <- dfr[order(dfr$new.R.sq), ]   # sort by increasing R-square

  list(obs.stats=cor.obj, dfr=signif(dfr, 4)) # return list of cor.test obj & sorted data frame

}  # end function definition



# demo the function


##
## sim data
##
#
# a <- leave.one.out.cor(x, y)
# str(a)
# print(a$obs.stats)
# print(a$dfr, digits=3)

# N <- length(x)
# indx.dwn <- a$dfr[1, 1]  # index of point at top of output table
# indx.up  <- a$dfr[N, 1] # index of point at bottom of output table
# obs.r <- cor(x, y)
#p <- cor.test(x,y)$p.value

# plot(x, y)
#confidence.ellipse(x,y)
#abline(lm(y~x), col="gray")
# points(x[indx.dwn], y[indx.dwn], pch=16, lwd=2, col="red")
# points(x[indx.up ], y[indx.up ], pch=16, lwd=2, col="blue")
#
# add ellipse for RED point removed, for BLUE point removed 
#
# mtext(paste("RED  point", "(#", indx.dwn, ") gone, correlation FIT drops the most"), 
#   col="red", side=3, adj=0, cex=1)
# mtext(paste("BLUE point", "(#", indx.up, ") gone, correlation FIT increases the most"), 
#   col="blue", side=3, adj=0, line=1, cex=1)
#
# mtext("OLS line only drawn to help visualize linear association",
#   side=1, adj=0, line=4, cex=.6)

#mtext(paste("Pearson r=", round(obs.r,4), ", entire data set", sep=""),
#   side=1, adj=1, line=4, cex=.8)

# plot(a$dfr$new.R.sq, a$df$p, 
#  xlab="R-sq, after point removed",
#  ylab="p-value, after point removed",
#  main="Leave-one-out analysis in Pearson correlation")
#
# abline(v=obs.r^2, lwd=2, col="red"); abline(h=p, lwd=2)
# mtext("RED line at observed R-sq (all x,y points included)", 
#   side=3, line=.7, col="red", cex=.85)
# mtext("BLACK line at observed p-value(2) (all x,y points included)", 
#   side=3, col="black", cex=.85)
#
# locator(2)


##
## testosterone data
##
#
# a <- leave.one.out.cor(age, testosterone)
# str(a)
# print(a$obs.stats)
#          Pearson's product-moment correlation
#
#  data:  x and y
#  t = -5.0468, df = 13, p-value = 0.0002237
#  alternative hypothesis: true correlation is not equal to 0
#  95 percent confidence interval:
#   -0.9358598 -0.5168499
#  sample estimates:
#         cor 
#  -0.8136803 

# print(a$dfr, digits=3)
#   point.rm  x    y Pearson.r        p LB95CI UB95CI Rsq.diff new.R.sq
#15       15 30  450    -0.776 0.001102 -0.926 -0.417 -0.05982    0.602
#1         1 11 1305    -0.784 0.000913 -0.928 -0.433 -0.04805    0.614
#3         3 13 1175    -0.800 0.000596 -0.934 -0.468 -0.02248    0.640
#13       13 24  625    -0.801 0.000579 -0.934 -0.470 -0.02087    0.641
#14       14 27  610    -0.805 0.000520 -0.936 -0.478 -0.01469    0.647
#12       12 23  700    -0.807 0.000484 -0.937 -0.484 -0.01070    0.651
#5         5 15 1060    -0.809 0.000461 -0.937 -0.487 -0.00795    0.654
#11       11 21  690    -0.811 0.000435 -0.938 -0.492 -0.00470    0.657
#7         7 16 1005    -0.811 0.000426 -0.938 -0.493 -0.00363    0.658
#10       10 20  605    -0.825 0.000281 -0.943 -0.524  0.01858    0.681
#2         2 12 1000    -0.829 0.000244 -0.944 -0.534  0.02581    0.688
#9         9 18 1150    -0.832 0.000224 -0.945 -0.540  0.03006    0.692
#6         6 16  800    -0.833 0.000216 -0.946 -0.542  0.03190    0.694
#8         8 17  710    -0.840 0.000172 -0.948 -0.557  0.04287    0.705
#4         4 14 1495    -0.851 0.000116 -0.952 -0.583  0.06138    0.723

# OPTIONAL code to visualize text output of correlation influence
# perhaps thess 2 graphs should be added to the function... 
#  
# N <- length(age)
# a <- leave.one.out.cor(age, testosterone)
# indx.fit.down <- a$dfr[1, 1]
# indx.fit.up   <- a$dfr[N, 1]

# plot(a$dfr$new.R.sq, a$dfr$p, xlab="new.R.sq", ylab="p-value", col="gray60")
# points(a$dfr$new.R.sq[1], a$dfr$p[1], cex=1.5, pch=16, lwd=2, col="red")
# points(a$dfr$new.R.sq[N], a$dfr$p[N], cex=1.5, pch=16, lwd=2, col="blue")
# mtext("RED point gone, fit drops most", col="red", side=3, adj=0, cex=.8)
# mtext("BLUE point gone, fit increases most", col="blue", side=3, adj=1, cex=.8)

# A new graph --
# plot(age, testosterone, col="gray60")
# points(age[indx.fit.down], testosterone[indx.fit.down], pch=16, lwd=2, col="red")
# points(age[indx.fit.up], testosterone[indx.fit.up], pch=16, lwd=2, col="blue")
# mtext("RED point gone, fit drops most", col="red", side=3, adj=0, cex=.8)
# mtext("BLUE point gone, fit increases most", col="blue", side=3, adj=1, cex=.8)


# ------------------------------------------------------------------------------
#  my.pause()  is called with no arguments to pause between text & graph output
# ------------------------------------------------------------------------------

my.pause <- function(){if(interactive())readline("Hit<Enter> to continue...");invisible()}  

# demo the function

# my.pause()   # no args in the call

# ---------------------------------
#  simple, line drawing functions 
# ---------------------------------

thick.line <- function(N=25) cat(rep("=", N), sep="", "\n")  # prints "="
thin.line  <- function(N=25) cat(rep("-", N), sep="", "\n")  # prints "-"

# demo the funs

# thick.line()
# thin.line(); thin.line()
# thin.line(80)


#
#
# ---------------  D E M O   D A T A  ---------------------
#
# Daniel (1991) p. 413. "...plasma testosterone levels and age at
# first conviction for violent and aggressive crimes collected on
# a sample of young male prisoners."
#
# TESTOSTERONE    AGE at FIRST    TESTOSTERONE    AGE at FIRST
# LEVEL ng/dL       CONVICTION      LEVEL           CONVICTION
# ----------------------------------------------------------------
#  1305               11             710              17     
#  1000               12            1150              18
#  1175               13             605              20
#  1495               14             690              21
#  1060               15             700              23
#   800               16             625              24
#  1005               16             610              27
#                                    450              30
# ----------------------------------------------------------------
# testosterone level, ng/dL
# testosterone <- c(1305, 1000, 1175, 1495, 1060, 800, 1005, 710, 1150,  605,  690,  700,  625, 610,  450)  
# age <- c(11, 12, 13, 14, 15, 16, 16, 17, 18, 20, 21, 23, 24, 27, 30)
#
# plot(age, testosterone)  # quick look at minimal scatterplot
#
# cor.test(age, testosterone) # default two-sided, Pearson r with its 95% CI by normal theory
# r = - 0.81368, p = 0.0002237 (2), 95% CI  LB, UB = -.936, -.517 (2)

# end.  Kincaid funs CORRELATION.R    AUTHOR: Dwight Kincaid, PhD 
