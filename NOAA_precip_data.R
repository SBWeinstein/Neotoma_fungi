#get annual precipitation data for each site
#adapted from code written by Marc T. Mayes

setwd("<directory_location>")

library(rnoaa)

#establish noaa token
noaakey <- "<NOAA_Key>"  

#vector of station IDs to retrieve data
stn <- c("GHCND:USC00029542",
		"GHCND:USC00042319",
		"GHCND:USC00047836",
		"GHCND:USC00042327",
		"GHCND:USC00044211",
		"GHCND:USC00046635",
		"GHCND:USC00044405",
		"GHCND:USC00047306",
		"GHCND:USC00046386",
		"GHCND:USC00262243",
		"GHCND:USC00265400",
		"GHCND:USC00421241",
		"GHCND:USC00422150",
		"GHCND:USC00423348",
		"GHCND:USC00422558",
		"GHCND:USC00425252")

out <- ncdc(datasetid='NORMAL_ANN', datatypeid='ANN-PRCP-NORMAL', 
	startdate = '2010-01-01', enddate = '2010-01-01', 
	token=noaakey, stationid=stn, limit=25, add_units=TRUE)

#see https://www1.ncdc.noaa.gov/pub/data/normals/1981-2010/readme.txt
#Note for units: "1" is 0.01 inches

out1<-data.frame(out$data)

write.csv(out1, file="ANN-PRCP-NORMAL_sites.csv")
