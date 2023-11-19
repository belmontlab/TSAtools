#' @title PRESS
#' @author Thomas Hopper
#' @description Returns the PRESS statistic (predictive residual sum of squares).
#'              Useful for evaluating predictive power of regression models.
#' @param linear.model A linear regression model (class 'lm'). Required.
#'
PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)
  return(PRESS)
}


### This calculate the Predictive r-squared

#' @title Predictive R-squared
#' @author Thomas Hopper
#' @description returns the predictive r-squared. Requires the function PRESS(), which returns
#'              the PRESS statistic.
#' @param linear.model A linear regression model (class 'lm'). Required.
#'
pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)
  return(pred.r.squared)
}



### This calculate the R-squared, the adj-R-squared and the Predictive R-squared (from the functions above)

#' @title Model Fit Statistics
#' @description Returns lm model fit statistics R-squared, adjucted R-squared,
#'      predicted R-squared and PRESS.
#'      Thanks to John Mount for his 6-June-2014 blog post, R style tip: prefer functions that return data frames" for
#'      the idea \link{http://www.win-vector.com/blog/2014/06/r-style-tip-prefer-functions-that-return-data-frames}
#' @return Returns a data frame with one row and a column for each statistic
#' @param linear.model A \code{lm()} model.
model_fit_stats <- function(model) {
  r.sqr <- summary(model)$r.squared
  adj.r.sqr <- summary(model)$adj.r.squared
  ratio.adjr2.to.r2 <- (adj.r.sqr/r.sqr)
  pre.r.sqr <- pred_r_squared(model)
  press <- PRESS(model)
  return.df <- data.frame("R-sq" = r.sqr, "Adj-R-sq" = adj.r.sqr,
                          "Adj.R2-to-R2" = ratio.adjr2.to.r2, "Pred-R-sq" = pre.r.sqr, PRESS = press)
  return(round(return.df,3))
}
#'TSA smthr
#'
#'This approximates TSA signal (or any other continous genomic signal) by a cubic polynomial.
#'The way it works is that it move a window along the genome and fit cubic polynomial in that window and reads the value at the centeral bin.
#' @param BW input BigWig file - GRanges object
#' @param bin bin window. Defaults to 20000 bp
#' @param smth_wind number of bins to include in polynomial smoothing
#' @param p_deg polynomial degree
#' @return smooth GRanges
#' @seealso
#' @export
#' @example x <- BW_smthr_zoo(bigwig, bin=2e4, smth_wind=20)
BW_smthr_zoo<-function(BW,bin,smth_wind,p_deg=3){
  require(rtracklayer)
  require(GenomicRanges)
  require(zoo)
  library(pracma)
  BW<-keepStandardChromosomes(BW, species = "Homo_sapiens",pruning.mode = "coarse")

  ## Defining central bin index
  if (is.numeric(smth_wind) && smth_wind%%1!=0){
    stop("Smoothing window should be an integer. you have provided: ",smth_wind)
  } else if (smth_wind%%2==1){
    center_ind<-(smth_wind%/%2)+1
    message("Using window size of ",smth_wind,". Center index is ", center_ind)
  } else if(smth_wind%%2==0){
    smth_wind<-smth_wind+1
    center_ind<-(smth_wind%/%2)+1
    message("Window size is even. Using window size of ",smth_wind," instead. Center index is ", center_ind)
  }

  BW_bnr<-function(BW,bin){
    bins<-tileGenome(seqinfo(BW),tilewidth = bin,cut.last.tile.in.chrom = TRUE)
    BW_Rle<-coverage(BW,weight="score")
    BW_bin<-binnedAverage(bins,BW_Rle,"score")
    return(BW_bin)
  }

  if (!is.na(bin)){
    BW_bn<-BW_bnr(BW,bin)
  } else{
    BW_bn<-BW
  }

  #
  ## NEED to generelize this function. To make it faster, maybe using as.function.polynomial or simplify predict.
  smth<-function(x){
    #fits a degree 3 polynomial
    if (any(is.na(x))|any(x==0)){
      p_vals<-as.data.frame(t(rep(0,p_deg)))
      colnames(p_vals)<-paste0("pval-c",seq(p_deg))
      model.stats<- data.frame("R-sq" = 0, "Adj-R-sq" = 0,
                    "Adj.R2-to-R2" = 0, "Pred-R-sq" = 0, PRESS = 0)
      return (data.frame(smth=0,
                         slp=0,
                         p_vals,
                         model.stats

      )
      )
    } else{

      poly_fit<-lm(x~poly(1:length(x),p_deg,raw = TRUE))
      fit_sum<-summary(poly_fit)
      coef_mat<-fit_sum$coefficients #Finding the coefs
      smth_coef<-rev(coef_mat[,"Estimate"])

      smth_coef_pval<-t(coef_mat[,'Std. Error'][-1]) # removing the intercept p-val
      colnames(smth_coef_pval)<-paste0("pval-c",seq(dim(smth_coef_pval)[2])) # changing the column names
      slp_coef<-polyder(smth_coef)

      model.stats<-model_fit_stats(poly_fit)

      return(data.frame(smth=polyval(smth_coef,center_ind), #evaluating poly at centeral index
                        slp=abs(polyval(slp_coef,center_ind)), #evaluating derivative at centeral index
                        # r_sq=fit_sum$r.squared, #getting r squared value of the fit
                        smth_coef_pval,
                        model.stats

      )
      )
    }

  }
  fits<-data.frame()
  for (chr in unique(seqnames(BW_bn))){
    message("Smoothing: ",chr)
    chr_sub<-keepSeqlevels(BW_bn,value = chr,pruning.mode = "coarse")
    fits<-rbind(fits,rollapply(data = chr_sub$score,width=smth_wind,FUN=smth, fill=0, align="center"))

  }
  mcols(BW_bn)<-fits
  valids<-which(apply(mcols(BW_bn),MARGIN = 1,FUN = sum)!=0) # filter the ones with all zero
  return(BW_bn[valids])
}
