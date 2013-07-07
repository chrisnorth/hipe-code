## Construct tables (FITS files) of beam areas for aperture photometry

def makedat(dirPath,profFile,filePSW,filePMW,filePLW):

	########################################################
	##Make table files
	########################################################
	from herschel.ia.toolbox.fit import FitterFunction
	from herschel.ia.numeric.toolbox import RealFunction
	from herschel.ia.numeric.toolbox.integr import TrapezoidalIntegrator
	from herschel.ia.numeric import String1d,Double1d,Double2d
	import math
	from herschel.ia.numeric.toolbox.interp import LinearInterpolator
	from herschel.ia.io.ascii import AsciiTableTool
	from herschel.ia.dataset import TableDataset,Column,StringParameter
	from herschel.ia.toolbox.image import Wcs,SimpleImage
	from herschel.ia.toolbox.util import simpleFitsWriter
	########################################################
	##Read in azimuthal beam profiles (made from beam maps)
	########################################################

	print dirPath+profFile
	f=open(dirPath+profFile,'r')
	lines=String1d(f.readlines())
	nl=lines.size
	
	rad=Double1d(nl)
	profPSW=Double1d(nl)
	profPMW=Double1d(nl)
	profPLW=Double1d(nl)
	
	minProf=1.e-8
	
	for l in range (nl):
		line=lines[l].split(',')
		rad[l]=float(line[0])
		profPSW[l]=max(float(line[1]),minProf)
		profPMW[l]=max(float(line[2]),minProf)
		profPLW[l]=max(float(line[3]),minProf)
	
	#####################################
	## Calculate beam areas
	#####################################
	
	## Calculate Integrand 2*pi*rad*B(rad)
	intargPSW=profPSW * 2. * math.pi * rad
	intargPMW=profPMW * 2. * math.pi * rad
	intargPLW=profPLW * 2. * math.pi * rad
	
	## Make interpolation of integrands
	intpPSW=LinearInterpolator(rad,intargPSW)
	intpPMW=LinearInterpolator(rad,intargPMW)
	intpPLW=LinearInterpolator(rad,intargPLW)
	
	areaPSW=Double1d(nl)
	areaPMW=Double1d(nl)
	areaPLW=Double1d(nl)
	#areaPSW2=Double1d(nl)
	#areaPMW2=Double1d(nl)
	#areaPLW2=Double1d(nl)
	print 'Calculating beam areas...'
	
	## Integrate using interpolation (takes ~10-20s)
	for r in range(nl):
		if r == 0:
			areaPSW[r]=0.
			areaPMW[r]=0.
			areaPLW[r]=0.
		else:
			areaPSW[r]=areaPSW[r-1] + profPSW[r]*2.*math.pi*rad[r]
			areaPMW[r]=areaPMW[r-1] + profPMW[r]*2.*math.pi*rad[r]
			areaPLW[r]=areaPLW[r-1] + profPLW[r]*2.*math.pi*rad[r]
	
		##BUILT IN INTEGRATOR MISBEHAVES NEAR zero-values!!!
		#intLim=TrapezoidalIntegrator(0.,rad[r])
		#areaPSW2[r]=intLim.integrate(intpPSW)
		#areaPMW2[r]=intLim.integrate(intpPMW)
		#areaPLW2[r]=intLim.integrate(intpPLW)
	#	print rad_r,areaPSW[r],areaPMW[r],areaPLW[r]
	
	#comp_PSW=ABS(areaPSW2-areaPSW)/areaPSW2
	#comp_PMW=ABS(areaPMW2-areaPMW)/areaPMW2
	#comp_PLW=ABS(areaPLW2-areaPLW)/areaPLW2
	
	maxAreaPSW=max(areaPSW)
	maxAreaPMW=max(areaPMW)
	maxAreaPLW=max(areaPLW)
	
	normPSW=areaPSW/maxAreaPSW
	normPMW=areaPMW/maxAreaPMW
	normPLW=areaPLW/maxAreaPLW
	
	ascii=AsciiTableTool()
	
	areaTable=TableDataset(description="Integrated beam areas")
	
	#ascii.formatter.header=1
	ascii.formatter.commented=1
	
	from herschel.share.unit.Angle import SECONDS_ARC
	from herschel.share.unit.Scalar import ONE
	
	areaTable['Radius']=Column(rad,unit=SECONDS_ARC,description='Integrated radius')
	areaTable['PSW']=Column(data=normPSW,unit=ONE,description='PSW beam integral')
	areaTable['PMW']=Column(normPMW,unit=ONE,description='PMW beam integral')
	areaTable['PLW']=Column(normPLW,unit=ONE,description='PLW beam integral')
	
	ascii.save('area.dat',areaTable)
	
	bgArrPSW=Double2d(nl,nl)
	bgArrPMW=Double2d(nl,nl)
	bgArrPLW=Double2d(nl,nl)
	#bgArrPSW[:]=None
	#bgArrPMW[:]=None
	#bgArrPLW[:]=None
	for r1 in range(nl):
		bgArrPSW[r1,r1:]=normPSW[r1:]-normPSW[r1]
		bgArrPMW[r1,r1:]=normPMW[r1:]-normPMW[r1]
		bgArrPLW[r1,r1:]=normPLW[r1:]-normPLW[r1]
		if r1 > 0:
			bgArrPSW[r1,0:r1-1]=normPSW[r1]-normPSW[0:r1-1]
			bgArrPMW[r1,0:r1-1]=normPMW[r1]-normPMW[0:r1-1]
			bgArrPLW[r1,0:r1-1]=normPLW[r1]-normPLW[0:r1-1]
	
	wcs=Wcs(crval1=0.,crval2=0.,crpix1=0.,crpix2=0.,cdelt1=1.,cdelt2=1.,\
		cunit1='arcsec',cunit2='arcsec')
	bgImgPSW=SimpleImage(description="PSW background subtraction",\
		image=bgArrPSW,wcs=wcs,instrument='SPIRE')
	bgImgPSW.meta["USAGE1"]=StringParameter("This is the PSW integrated beam between radius1 and radius2")
	bgImgPSW.meta["USAGE2"]=StringParameter("The beam between r1 and r2 is DATA[r1,r2] (order not important)")
	bgImgPSW.meta["wavelength"]=StringParameter("PSW")
	simpleFitsWriter(bgImgPSW,file=dirPath+filePSW)
	
	bgImgPMW=SimpleImage(description="PSW background subtraction",\
		image=bgArrPMW,wcs=wcs,instrument='SPIRE')
	bgImgPMW.meta["USAGE1"]=StringParameter("This is the PSW integrated beam between radius1 and radius2")
	bgImgPMW.meta["USAGE2"]=StringParameter("The beam between r1 and r2 is DATA[r1,r2] (order not important)")
	bgImgPMW.meta["wavelength"]=StringParameter("PMW")
	simpleFitsWriter(bgImgPMW,file=dirPath+filePMW)
	
	bgImgPLW=SimpleImage(description="PSW background subtraction",\
		image=bgArrPLW,wcs=wcs,instrument='SPIRE')
	bgImgPLW.meta["USAGE1"]=StringParameter("This is the PSW integrated beam between radius1 and radius2")
	bgImgPLW.meta["USAGE2"]=StringParameter("The beam between r1 and r2 is DATA[r1,r2] (order not important)")
	bgImgPLW.meta["wavelength"]=StringParameter("PLW")
	simpleFitsWriter(bgImgPLW,file=dirPath+filePLW)


def plotbeams(rad,profPSW,profPMW,profPLW,areaPSW,areaPWM,areaPLW):

	#####################################
	## Plot beam functions
	#####################################
	
	plot=PlotXY()
	plot.autoBoxAxes=1
	plot.title.setText('Beam Profiles')
	
	## add PSW in blue
	lyrPSW=LayerXY(rad,profPSW,color=java.awt.Color(0,0,150))
	plot.addLayer(lyrPSW)
	
	## set axes ranges
	plot.yaxis.range = [1.e-8,1.]
	plot.yaxis.type = Axis.LOG
	plot.xaxis.range = [0,1000]
	
	## PMW in green
	lyrPMW=LayerXY(rad,profPMW,color=java.awt.Color(0,150,0))
	plot.addLayer(lyrPMW)
	
	## PLW in red
	lyrPLW=LayerXY(rad,profPLW,color=java.awt.Color(150,0,0))
	plot.addLayer(lyrPLW)
	
	#####################################
	## Plot beam areas
	#####################################
	
	plota=PlotXY()
	plota.autoBoxAxes=1
	plota.title.setText('Beam Areas')
	
	## PSW in blue
	lyrAreaPSW=LayerXY(rad,areaPSW,color=java.awt.Color(0,0,150))
	plota.addLayer(lyrAreaPSW)
	
	## Set axes ranges
	#plota.yaxis.range = [1.e-8,1.]
	#plota.yaxis.type = Axis.LOG
	plota.xaxis.range = [0,500.]
	
	## PMW in green
	lyrAreaPMW=LayerXY(rad,areaPMW,color=java.awt.Color(0,150,0))
	plota.addLayer(lyrAreaPMW)
	
	##PLW in red
	lyrAreaPLW=LayerXY(rad,areaPLW,color=java.awt.Color(150,0,0))
	plota.addLayer(lyrAreaPLW)
	
#	#####################################
#	## Plot beam area comparisons
#	#####################################
#	
#	plotc=PlotXY()
#	plotc.autoBoxAxes=1
#	plotc.title.setText('Beam Area Comparison')
#	
#	## PSW in blue
#	lyrCompPSW=LayerXY(rad,comp_PSW,color=java.awt.Color(0,0,150))
#	plotc.addLayer(lyrCompPSW)
#	
#	## Set axes ranges
#	#plotc.yaxis.range = [1.e-8,1.]
#	plotc.yaxis.type = Axis.LOG
#	plotc.xaxis.range = [0,500.]
#	
#	## PMW in green
#	lyrCompPMW=LayerXY(rad,comp_PMW,color=java.awt.Color(0,150,0))
#	plotc.addLayer(lyrCompPMW)
#	
#	##PLW in red
#	lyrCompPLW=LayerXY(rad,comp_PLW,color=java.awt.Color(150,0,0))
#	plotc.addLayer(lyrCompPLW)
