#'#' This is some descriptio of this function.
#' @title realGARCH_select
#'
#' @description select distribution of innovation along with the realizedVol in realGARCH based on specified methods
#'
#' @details
#'
#' @param gp is high frequency data
#'
#' @return text



realGARCH_select <- function(gp,insample=F) {
  rm(list=ls())

  library(rugarch)
  library(stargazer)
  # ---------read data----------------#
  library(readxl)
  # gp<- read_excel("2019_2022_5min_000016SH.xlsx", col_types = c("date", "numeric", "numeric", "numeric", "numeric", "numeric"))

  # gp$date=paste(substr(gp$time,1,4),substr(gp$time,6,7), substr(gp$time,9,10), sep = "")
  gp$date=substr(gp$time,1,10)
  dates=sort(unique(gp$date))
  gp$rtn=log(gp$pct_chg+1)
  gp=gp[!is.na(gp$time),]
  gp=gp[,c("time","date","open","close","high","low")]
  head(gp)
  library(stargazer)
  stargazer(head(gp),summary = F,type = 'text')

  daily=data.frame(dates)
  colnames(daily)=c('date')

  daily$lgrtn=c()
  daily$RV=c()
  daily$RRV=c()
  daily$RBV=c()
  daily$RRBV=c()
  daily$close=c()
  for(i in 2:length(dates)){
    tmp=gp[gp$date==dates[i],]
    tn=length(tmp)
    daily$RV[i]=sum(diff(log(tmp$close))^2)
    daily$RRV[i]=sum((log(tmp$high)-log(tmp$low))^2/4/log(2))
    daily$RBV[i]=pi/2*tn/(tn-1)*sum(abs(diff(log(tmp$close))[1:(tn-2)])*abs(diff(log(tmp$close))[2:(tn-1)]))
    daily$RRBV[i]=pi/8*sum(abs(log(tmp$high[1:(tn-1)])-log(tmp$low[1:(tn-1)]))*abs(log(tmp$high[2:tn])-log(tmp$low[2:tn])))
    daily$close[i]=tmp$close[tn]
    #     if(i==1){
    #       daily$rtn[i]=0
    #   }else{
    #       daily$rtn[i]=(tmp$close[tn]-tmp$open[1])/tmp$open[1]
    #   }
  }

  # ---------------------------generate log-rtn---------------------------------------#
  for(i in 2:length(dates)){
    cur=gp[gp$date==dates[i],]
    pre=gp[gp$date==dates[i-1],]
    tn=length(cur)
    daily$lgrtn[i]=log(cur$close[tn])-log(pre$close[tn])
  }

  # ---------------------------generate xts object---------------------------------------#
  library(xts)
  daily=daily[2:dim(daily)[1],]
  nd=xts(daily[,c('lgrtn','RV','RRV','RBV','RRBV','close')],order.by = as.Date(daily$date) )
  nd_df=as.data.frame(nd)
  head(nd)


  # ---------------------------return modeling---------------------------------------#
  library(rugarch)
  rtnfit=autoarfima(nd[,1], ar.max = 3, ma.max = 3, criterion = "AIC", include.mean = F,arfima = T,method = 'full')

  rtnfit

  kk=row.names(data.frame(rtnfit$fit@fit$coef))

  # ---------------------------in-sample way of section--------------------------------------#
  if(insample==T){
    sec_dtb=c('snorm','std','sstd','ged','sged')
    tab_in=matrix(0,nrow = 4,ncol = 5)
    for(i in 1:4){
      for(j in 1:5){
        spec_tmp <- ugarchspec(mean.model = list(armaOrder=c(3,2),arfima=T, include.mean = F),
                               variance.model = list(model = 'realGARCH',garchOrder = c(1, 1)),
                               fixed.pars=list(arfima = 0.046141),
                               distribution.model = sec_dtb[j])
        fit_tmp <- ugarchfit(spec=spec_tmp, nd[,1], solver = 'hybrid', realizedVol=nd[,i+1])

        tab_in[i,j]=infocriteria(fit_tmp)[1]
      }
    }
  }else{
    library(GAS)
    tab_out=matrix(0,nrow = 4,ncol = 5)
    for(i in 1:4){
      for(j in 1:5){
        i=4
        j=3
        spec_tmp <- ugarchspec(mean.model = list(armaOrder=c(3,2),arfima=T, include.mean = F),
                               variance.model = list(model = 'realGARCH',garchOrder = c(1, 1)),
                               fixed.pars=list(arfima = 0.046141),
                               distribution.model = sec_dtb[j])
        pf_show=ugarchroll(spec = spec_tmp,data = nd[,1], realizedVol=nd[,i+1],
                           n.ahead = 1,n.start = 400,  refit.every = 50, refit.window = "recursive",fit.control = list(),solver = "hybrid",
                           calculate.VaR = TRUE, VaR.alpha = 0.05,
                           keep.coef = F)
        df=data.frame(pf_show@forecast$VaR)
        if(sum(dim(df))==0) {
          tab_out[i,j]=0
          next
        }
        colnames(df)=c('a005','real')
        t=BacktestVaR(df$real, df$a005, 0.05)

        tab_out[i,j]=(t$LRcc[2]+t$DQ$pvalue)/2

      }
    }
  }

  # ---------------------------out-sample way of section--------------------------------------#



}




