######################################################################################################
################################## SIMPLE WATER BALANCE MODEL ########################################
######################################################################################################
# Code by Rene Orth, passed on to Beni Stocker, July 2016


simplemodel <- function ( exp_runoff, exp_et, beta, whc, melting ) {
  
  library('fields')


  ##################### read forcing data for site Payerne ####################################################
  length <- 5475
  temp   <- array(NaN,c(length))
  precip <- array(NaN,c(length))
  rad    <- array(NaN,c(length))

  temp[]   <- t(read.table('./Data/pay_meteo_1998-2012.txt',skip=10,header=T)[1:5475,7])
  precip[] <- t(read.table('./Data/pay_meteo_1998-2012.txt',skip=10,header=T)[1:5475,9])
  rad[]    <- t(read.table('./Data/pay_meteo_1998-2012.txt',skip=10,header=T)[1:5475,8])

  # Scaling solar radiation to approximate net radiation (xxx ref? xxx)
  rad <- rad * 0.65 - 35.0

  # convert radiation into mm
  rad <- rad * 86400 / 2260000

  #############################################################################################################

  ## define arrays
  soilm <- array(NaN,c(length))
  runoff <- array(NaN,c(length))
  et <- array(NaN,c(length))
  et_corr <- array(NaN,c(length))
  snow <- array(NaN,c(length))
  infiltration <- array(NaN,c(length))
  infiltration_corr <- array(NaN,c(length)) 
  #runoffsum <- array(0,c(length,100))

  ############# ACCOUNTING FOR SNOW & DEW, converting precipitation to rainfall+snow melt+dew melt #############
  ## threshold temperature in deg C
  temp_grenze <- 1.0

  for (iterate in 1:2){

    ## add dew to precipitation
    precip[ which(rad<0) ] <- precip[ which(rad<0) ] + (-1) * beta * rad[ which(rad<0) ]
    rad[ which(rad<0) ] <- 0

    ## first start from zero snow, in the second loop start from average 31 December snow
    if (iterate == 1){
      snow_current <- 0
    } else {
      snow_current <- mean( snow[ seq(365,length,365) ] )
    }

    ## i is the day
    for (i in 1:length){

      ## form snow and decrease precip accordingly if precip falls at less than 1 deg C
      if ( precip[i]>0 && temp[i]<(temp_grenze+1.0) ){

        ## depending on temperature, all precipitation or only some fraction is converted to snow
        if (temp[i]<((temp_grenze-1.0))){
          snow_current <- snow_current + precip[i]
          precip[i]    <- 0.0
        } else {
          fsnow        <- ( temp[i] - ( temp_grenze - 1.0 ) ) / 2.0
          snow_current <- snow_current + precip[i] * ( 1.0 - fsnow )
          precip[i]    <- precip[i] * fsnow
        }

      }

      ## melt snow and increase precip accrodingly if snow is present and temperature is above 1 deg C
      if (snow_current > 0 && temp[i]>temp_grenze){
        melt <- min( snow_current, melting * ( temp[i] - temp_grenze ) )
        snow_current <- snow_current - melt
        precip[i] <- precip[i] + melt
      }

      snow[i] <- snow_current

    }

  }
  #############################################################################################################



  ##################### compute delayed runoff ################################################################
  # runoffsum[] <- 0
  # 
  # # to compute delayed runoff, precip from previous is days is considered if 
  # # (i) it occurred max 100 days ago, and
  # # (ii) more than 0.5% of that daily precip would still contribute to the runoff on the current day 
  # # --> this way, 99.5 % of all precip will eventually be captured
  # maxdays=min(100,which(exp((-1)*runoffdelaytimescale*1:1000)<0.005)[1])
  # 
  # for (i in 2:le){
  # for (j in c(which(precip[i-1:(min((i-1),maxdays))]>0))){
  # runoffsum[i,j]=precip[i-j]*(exp((-1)*runoffdelaytimescale*(j-1))-exp((-1)*runoffdelaytimescale*j))
  # }
  # }
  #############################################################################################################



  ############################# SPIN UP, 5 years ##############################################################

  ## initializing soil moisture to 90% of water holding capacity
  soilm[1] <- 0.9 * whc

  for (i in 2:1825){

    ## calculate ET from net radiation and Eq. 2 in Orth et al., 2013, limited to <=1
    et[i-1] <- rad[i-1] * beta * min( 1.0, ( soilm[i-1] / whc ) ^ exp_et )

    ## calculate derivative of ET w.r.t. soil moisture
    et_corr[i-1] <- rad[i-1] * beta * min( max( 0.0, whc - soilm[i-1]), ( exp_et / whc ) ) * ( soilm[i-1] / whc ) ^ (exp_et - 1.0)

    ## calculate infiltration (P-Q) from Eq. 3 in Orth et al., 2013
    infiltration[i-1] <- ( 1.0 - min( 1.0,( ( soilm[i-1] / whc) ^ exp_runoff ) ) ) * precip[i-1]

    ## calculate derivative of infiltration w.r.t. soil moisture
    infiltration_corr[i-1] <- (-1) * min( max( 0, whc - soilm[i-1] ), ( exp_runoff / whc ) ) * ( ( soilm[i-1] / whc ) ^ ( exp_runoff - 1.0 ) ) * precip[i-1]

    ## XXX is 5.0 a permanent wilting point parameter?
    et[i-1] <- min( et[i-1], soilm[i-1] - 5.0 )

    ## implicit solution, see Eq. 7 in Orth et al., 2013
    soilm[i] <- soilm[i-1] + ( (infiltration[i-1] - et[i-1]) / ( 1.0 + et_corr[i-1] - infiltration_corr[i-1]) )
  
  }
  #############################################################################################################



  ############################## ACTUAL MODEL RUN #############################################################

  ## initializing soil moisture to mean of 31 December values from spin-up
  soilm[1] <- mean(c(soilm[730],soilm[1095],soilm[1460],soilm[1825]))

  for (i in 2:length){

    ## calculate ET from net radiation and Eq. 2 in Orth et al., 2013, limited to <=1
    et[i-1] <- rad[i-1] * beta * min( 1.0, ( soilm[i-1] / whc ) ^ exp_et )

    ## calculate derivative of ET w.r.t. soil moisture
    et_corr[i-1] <- rad[i-1] * beta * min( max( 0.0, whc - soilm[i-1]), ( exp_et / whc ) ) * ( soilm[i-1] / whc ) ^ (exp_et - 1.0)

    ## calculate infiltration (P-Q) from Eq. 3 in Orth et al., 2013
    infiltration[i-1] <- ( 1.0 - min( 1.0,( ( soilm[i-1] / whc) ^ exp_runoff ) ) ) * precip[i-1]

    ## calculate derivative of infiltration w.r.t. soil moisture
    infiltration_corr[i-1] <- (-1) * min( max( 0, whc - soilm[i-1] ), ( exp_runoff / whc ) ) * ( ( soilm[i-1] / whc ) ^ ( exp_runoff - 1.0 ) ) * precip[i-1]

    ## XXX is 5.0 a permanent wilting point parameter?
    et[i-1] <- min( et[i-1], soilm[i-1] - 5.0 )

    ## implicit solution, see Eq. 7 in Orth et al., 2013
    soilm[i] <- soilm[i-1] + ( (infiltration[i-1] - et[i-1]) / ( 1.0 + et_corr[i-1] - infiltration_corr[i-1]) )

    runoff[i-1] <- ( min( 1.0, ( ( soilm[i-1] / whc ) ^ exp_runoff ) ) ) * precip[i-1]

    # in the case of delayed runoff, use instead:
    # for (j in c(which(precip[i-1:(min((i-1),maxdays))]>0)) ){
    # runoff[i-1] <- runoff[i-1]+min(1,((soilm[i-j+1]/whc)^exp_runoff))*runoffsum[i,j]
    # }

    et[i-1] <- et[i-1] + ( soilm[i] - soilm[i-1] ) * et_corr[i-1]

  }
  
  #############################################################################################################

  list <- list( soilm, runoff, et, snow, exp_runoff, exp_et, beta, whc, melting, length )
  names(list) <- c( 'soilmoisture', 'runoff', 'evapotranspiration', 'snow', 'exp_runoff', 'exp_et', 'beta', 'whc', 'melting', 'length' )

  return(list)

}
