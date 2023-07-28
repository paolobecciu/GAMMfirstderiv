# GAMMfirstderiv
Example data and script to calculate first derivative from GAMM as in Becciu et al. 2023 "Soaring migrants flexibly respond to sea-breeze in a migratory bottleneck: using first derivatives to identify behavioural adjustments over time" Movement Ecology (https://doi.org/10.1186/s40462-023-00402-4)


## The data
Find description of the data in the Method section of the manuscript. 
For practical purposes here a list of description of variables (skipping the obvious ones or those not used in the paper):

- n.days: julian date each year starting from 1 on the first of January
- h.from.sunset: hours between sunset and the recording time (in negative decimal hours; i.e. -9:30 = -9.5)
- h.sset: round of the previous variable (i.e. -9.89 = -10)
- n.tracks: number of tracks recorded in that radar session (radar sessions are every 15/30 minutes)
- gspeed: average groundspeed of the radar tracks recorded in that specific radar session
- Wrad: circular mean wind direction in radians
- Brad: circular mean track direction in radians
- Ws: wind speed
- dist.coast_km: average distance to the coast
- CW: crosswind component of the wind relative to the mean migration direction (187.7 degrees)
- TW: tailwind component of the wind relative to the mean migration direction (187.7 degrees)
- aspeed: airspeed calculated relative to the CW and TW components of the wind
- s.aspeed187: lateral component of airspeed (LcA) 
