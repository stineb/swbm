plotting <- function (list) {

  melting_default <- 3
  whc_default <- 220
  exp_runoff_default <- 6.4
  exp_et_default <- 0.06
  beta_default <- 0.66

  ## have default parameters been modified?
  if ( list$exp_runoff != exp_runoff_default || list$exp_et != exp_et_default || list$beta != beta_default || list$whc != whc_default || list$melting != melting_default ){
    modified <- TRUE
  } else { 
    modified <- FALSE
  }

  ######### plot functions of normalized ET and normalized runoff versus soil moisture ########################
  pdf( "./fig/vars_vs_sm.pdf", width=7, height=4.5 )

    ## compute and plot default functions for site Payerne
    help1 <- ((seq(1,1.05*whc_default,1))/whc_default)^(exp_runoff_default)
    help2 <- beta_default*((seq(1,1.05*list$whc,1))/list$whc)^(exp_et_default)
    help1[1] <- 0
    help2[1] <- 0
    help1[which(help1>1)] <- 1
    help2[which(help2>1)] <- 1

    par( mar=c(3,3,0.5,0.5) ,mgp=c(1.5,0.3,0), tck=-0.02 )

    plot( help1, type='l', col=2, ylim=c(-0.04,1), lwd=2, ylab='Q/P and ET/Rnet', xlab='Soil Moisture (mm)' )
    lines( help2, lwd=2 )
    abline( h=0 )

    if (modified){
      legend( 'topleft', c("Q/P default","Q/P current",'ET/Rnet default','ET/Rnet current'), col=c(2,2,1,1), text.col='black', lwd=c(1.5,1.5,1.5,1.5), lty=c(1,2,1,2), cex=0.85 )
      legend( 80, 1.04, c("Soil Moisture Histogram Default","Soil Moisture Histogram Current"), col=c('grey85',1), text.col='black', lwd=c(10,1.5), lty=c(1,2), cex=0.85 )
    } else {
      legend('topleft', c("Q/P",'ET/Rnet'), col=c(2,1), text.col='black', lwd=c(1.5,1.5), cex=0.85 )
      legend( 80, 1.04, c("Soil Moisture Histogram Default"), col=c('grey85'), text.col='black', lwd=c(10), lty=c(1), cex=0.85 )
    }

    help3 <- help1
    help4 <- help2

    ## if changed, plot actual functions for site Payerne
    if (modified){
      help3 <- ((seq(1,1.05*list$whc,1))/list$whc)^(list$exp_runoff)
      help4 <- list$beta*((seq(1,1.05*list$whc,1))/list$whc)^(list$exp_et)
      help3[1] <- 0
      help4[1] <- 0
      help3[which(help3>1)] <- 1
      help4[which(help4>1)] <- 1

      lines(help3,col=2,lwd=2,lty=2)
      lines(help4,lwd=2,lty=2)
    }
    #############################################################################################################

    ########### add soil moisture histogram and display SM-ET coupling strength ################################
    summer <- array(NaN,c(list$length))
    years <- list$length/365

    for(l in 0:(years-1)){
      for(m in 121:273){
        summer[m+l*365] <- 1
      }
    }

    ## add histogram into functions plot
    soilm_default <- read.table('./Data/SoilMoisture_SimpleModel_Payerne_1998-2012')[,4]
    
    histo <- hist( soilm_default*summer, plot=F, breaks=11 )
    barplot( c(0,histo$counts/(1.5*max(histo$counts))), c(histo$breaks[1],rep(20,13)), space=0, add=T, col='grey85', border=F )

    sm_et_coupling <- array(NaN,c(4))

    for (i in 1:11){
      counts=length(which(list$soilmoisture[which(list$soilmoisture>((i-1)*20))]*summer[which(list$soilmoisture>((i-1)*20))]<((i)*20)))/(1.5*max(histo$counts))
      counts_next=length(which(list$soilmoisture[which(list$soilmoisture>((i)*20))]*summer[which(list$soilmoisture>((i)*20))]<((i+1)*20)))/(1.5*max(histo$counts))

      if (i == 1){ lines(c(0,0),c(0,counts),lty=2) }

      lines( c(0+(i-1)*20, 20+(i-1)*20), rep(counts,2), lty=2 )
      lines( c(20+(i-1)*20,20+(i-1)*20), c(counts,counts_next), lty=2 )
    }

    lines( help1, col=2, lwd=2 )
    lines( help2, lwd=2 )
    lines( help3, col=2, lwd=2, lty=2 )
    lines( help4, lwd=2, lty=2 )

    ## add information on SM-ET coupling strength 
    sm_et_coupling <- array(NaN,c(4))

    for (i in 1:4){
      help <- which(list$soilmoisture[which(list$soilmoisture>((i-1)*50))]*summer[which(list$soilmoisture>((i-1)*50))]<((i)*50))
      sm_et_coupling[i] <- cor(list$soilmoisture[which(list$soilmoisture>((i-1)*50))[help]],list$evapotranspiration[which(list$soilmoisture>((i-1)*50))[help]])
    }

    par( mgp=c(1.5,-1.2,0), options(warn=-1) )
    axis( 1, at=c(seq(25,175,50),210), labels=c(round(sm_et_coupling,digits=2),'cor (SM,ET)'), lwd=2, tick=F, col.axis=4 )
    par( mgp=c(1.5,0.3,0) )

  dev.off()
  #############################################################################################################


  ################# plot modeled versus observed soil moisture, including climatology ##########################
  soilm_clim <- array(NaN,c(list$length))
  for (i in 1:365){
    soilm_clim[i] <- mean(list$soilmoisture[i+(0:(years-1))*365])
  }
  soilm_clim <- rep(soilm_clim[1:365],years)


  pdf( "./fig/sm_tseries.pdf", width=9, height=4 ) 

    par( mar=c(1.6,2.8,0.7,0.5), mgp=c(1.5,0.3,0), tck=-0.02 )
    
    plot( soilm_default[3651:list$length], type='l', xaxt='n', xlab='', ylab='Soil Moisture (mm)', ylim=c(-20,240) )

    axis( 1, at=c(0:5*365), labels=F)
    par(tck=0)
    axis( 3, at=c(0:5*365), labels=F)
    axis( 1, at=c(seq(1,9,2)*182.5), labels=c('2008','2009','2010','2011','2012') )
    par(tck=-0.02)

    lines( soilm_clim[3651:list$length], lwd=2.5, col='grey70' )
    lines( soilm_default[3651:list$length], lwd=2.5 )

    if ( modified ){ lines( list$soilmoisture[3651:list$length], lwd=2, col=3 ) }

    ## read observed soil moisture
    soilm_obs=array(NaN,c(list$length))

    s <- read.table('./Data/PAY_SM_2008-12.txt',skip=16)
    co <- 3870
    for (i in seq(1,1606*24,24)){
      soilm_obs[co] <- (7.5*mean(s[i:(i+23),7])+12.5*mean(s[i:(i+23),8]) + 20*mean(s[i:(i+23),9]) + 25*mean(s[i:(i+23),10])+ 30*mean(s[i:(i+23),11]))/95
      co <- co+1
    }

    ## scale observed soil moisture
    soilm_obs <- soilm_obs*(sd(soilm_default[3651:list$length],na.rm=T)/sd(soilm_obs[3651:list$length],na.rm=T))
    soilm_obs <- soilm_obs-(mean(soilm_obs[3651:list$length],na.rm=T)-mean(soilm_default[3651:list$length],na.rm=T))

    ## plot observed soil moisture
    lines(soilm_obs[3651:list$length],col=2)

    if (modified){
      legend('bottomleft', c("Climatology Default",'Modelled Default','Modelled current',
      'Observed'),col = c('grey70',1,3,2),lty=c(1,1,1,1),
      text.col = 'black', lwd=c(1.5,1.5,1.5,1.5),cex=0.85)
    } else {
      legend('bottomleft', c("Climatology",'Modelled','Observed'),col = c('grey70',1,2),
      text.col = 'black', lwd=c(1.5,1.5,1.5),cex=0.85)
    }

  dev.off()
  #############################################################################################################

}
