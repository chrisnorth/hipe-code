###plot aperture photometry on sources

obsid=1342213459

###read in maps from files
file_PSW='/data/Herschel/Calibration/RelGains/Obs/mapPSW_'+str(hex(obsid))[2:10]+'.fits'
print file_PSW
file_PMW='/data/Herschel/Calibration/RelGains/Obs/mapPMW_'+str(hex(obsid))[2:10]+'.fits'
file_PLW='/data/Herschel/Calibration/RelGains/Obs/mapPLW_'+str(hex(obsid))[2:10]+'.fits'

mapPSW = fitsReader(file = file_PSW)
mapPMW = fitsReader(file = file_PMW)
mapPLW = fitsReader(file = file_PLW)

###set radius range
radius=[]
rad=1.0
rstep=1.0
nrad=20

###set-up lists
apphotPSW_tar=[]
apphotPSW_bg=[]
apphotPSW_sub=[]

#apphotPMW_tar=[]
#apphotPMW_bg=[]
#apphotPMW_sub=[]
#apphotPLW_tar=[]
#apphotPLW_bg=[]
#apphotPLW_sub=[]

for i in range(nrad):
	radius.append(rad)
	result_PSW = annularSkyAperturePhotometry(image=mapPSW, centerX=141.0, centerY=144.0, radiusPixels=rad, innerPixels=rad, outerPixels=100.0, fractional=1, algorithm=0)
	apphotPSW_tar.append(result_PSW["Results table"][0].data[0])
	apphotPSW_bg.append(result_PSW["Results table"][0].data[1])
	apphotPSW_sub.append(result_PSW["Results table"][0].data[2])
	
	#result_PMW = annularSkyAperturePhotometry(image=mapPMW, centerX=87., centerY=88.5, radiusPixels=rad, innerPixels=rad, outerPixels=100.0, fractional=1, algorithm=0)
	#apphotPMW_tar.append(result_PMW["Results table"][0].data[0])
	#apphotPMW_bg.append(result_PMW["Results table"][0].data[1])
	#apphotPMW_sub.append(result_PMW["Results table"][0].data[2])
	
	#result_PLW = annularSkyAperturePhotometry(image=mapPLW, centerX=61.5, centerY=62.5, radiusPixels=rad, innerPixels=rad, outerPixels=100.0, fractional=1, algorithm=0)
	#apphotPLW_tar.append(result_PLW["Results table"][0].data[0])
	#apphotPLW_bg.append(result_PLW["Results table"][0].data[1])
	#apphotPLW_sub.append(result_PLW["Results table"][0].data[2])
	
	print 'radius: %d arcsec'%(int(rad))
	rad=rad + rstep

###make plots

radP=Double1d(radius)
tar_PSW=Double1d(apphotPSW_tar)
bg_PSW=Double1d(apphotPSW_bg)
sub_PSW=Double1d(apphotPSW_sub)

plot=PlotXY()
plot.autoBoxAxes=1
ltar_PSW=LayerXY(radP,tar_PSW)
lbg_PSW=LayerXY(radP,bg_PSW)
lsub_PSW=LayerXY(radP,sub_PSW)

ltar_PSW.color=java.awt.Color.red
lbg_PSW.color=java.awt.Color.blue
lsub_PSW.color=java.awt.Color.red
lsub_PSW.line=Style.DASHED
plot.addLayer(ltar_PSW)
plot.addLayer(lbg_PSW)
plot.addLayer(lsub_PSW)
###end of program
