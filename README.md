# Documentation of Simple Water Balance Model (SWBM)
### Rene Orth, 14 July 2016


- R code of the model is provided in simplemodel.R
  - SWBM set up to compute hydrology at site Payerne
  - sections concerning delayed runoff are commented out
  - the delyed runoff requires an additional parameter and takes more time to compute
  - it is required if daily runoff should be analyzed
  - it is not required if the focus is on soil moisture, ET, or monthly runoff

- R code of to plot simple model output is provided in figures.R

- The model is tested with R Versions >= 2.10

- Meteorological input data file required to run simplemodel.R:
```
./Data/pay_meteo_1998-2012.txt
```

- Input data files required to run figures.R:
```
./SoilMoisture_SimpleModel_Payerne_1998-2012
./Data/PAY_SM_2008-12.txt
```


## Example to run SWBM at site Payerne

Load functions of SWBM
```R
source('simplemodel.R')
source('plotting.R')
```

Execute model with default parameters and create standard plots (make sure to create directory `./fig` beforehand)
```R
output = simplemodel( exp_runoff=6.4, exp_et=0.06, beta=0.66, whc=220, melting=3 )
str(output)
plotting(output)
```

Execute model with  changed parameters and create standard plots
```R
output=simplemodel( exp_runoff=2.4, exp_et=0.06, beta=0.86, whc=220, melting=3 )
plotting(output)
```


## General

- Model is introduced in:
  * Orth, Koster, Seneviratne 2013, J.Hydrometeorol. (based on Koster and Mahanama, 2012, J.Hydrometeorol.)

- Modifications to the model described in 
  * Orth et al. 2015, J.Hydrol [update to snow module]
  * Orth and Seneviratne, 2015, ERL [accounting for dew formation]

### Input for SWBM
- Precipitation
- Net Radiation
- Temperature

### Output of SWBM
- Soil Moisture
- Runoff
- Evapotranspiration
- Snow Water Equivalent

### Methods to calibrate SWBM:

- local-scale: using streamflow observations from near-natural (usually rather small) catchments
  * applied in Orth et al. 2015, J.Hydrol & Orth and Seneviratne, 2014, Clim.Dyn & Orth et al. 2013, J.Hydrometeorol.
- continental-scale: using multiple hydrological reference datasets
  * applied in Orth and Seneviratne, 2015, ERL & Orth et al., 2016, Scientific Reports 

- continental-scale calibration (applicable all over Europe): 
  * exp_runoff = 3.89, exp_et = 1.14, beta = 0.674 (maximum of normalized ET function), whc = 970.5 [mm], melting = 10.39 [mm snow water equivalent per degreee C above 1 deg C]

