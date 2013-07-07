import sys
from herschel.ia.gui.plot.renderer.PComponentEngine import HAlign
from herschel.ia.gui.plot.renderer.PComponentEngine import VAlign

from herschel.ia.numeric.toolbox.basic import Histogram
from herschel.ia.numeric.toolbox.basic import BinCentres

## Plot aperture photometry stuff

pc=False
pf=False


dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/Obs/'
dirPlot = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/Plots/'
##Get list of ObsIDs

fileObsName='AlphaBoo_obs_All.dat'
tarName='AlphaBoo'

#fileObsName='GammaDra_obs_All.dat'
#tarName='GammaDra'

fileObsList=dirPath+fileObsName

bandStr=["PSW","PMW","PLW"]

fObs=open(fileObsList,'r')
obsidsIn=fObs.readlines()
fObs.close()
obsids=[]
raSrcs=[]
decSrcs=[]
nobIn=len(obsidsIn)
for ob in range(nobIn):
	line=obsidsIn[ob]
	if line.find('#') < 0:
		obsids.append(int(line.split(',')[0],16))
		raSrcs.append(float(line.split(',')[1]))
		decSrcs.append(float(line.split(',')[2]))
obsids=Int1d(obsids)
raSrcs=Double1d(raSrcs)
decSrcs=Double1d(decSrcs)
nob=len(obsids)

#radiusArcsec=[[22.,30.,45.,50.,60.,100.],\
#	[30.,45.,55.,65.,80.,100.],\
#	[42.,50.,60.,90.,100.,120.]]
##Max source radii that are valid
#radMax=[[5,6,6,6],[3,6,6,6],[3,4,6,6]]

radiusArcsec=[[22.],[30.],[42.]]
radMax=[[1],[1],[1]]

#innerArcsec=[60.,100.,200.,300.]
#outerArcsec=[90.,150.,250.,350.]
##Min BR radii that are valid
#bgMin=[[0,0,0,0,0,1],[0,0,0,1,1,1],[0,0,0,1,1,2]]

innerArcsec=[60.]
outerArcsec=[90.]
bgMin=[[0],[0],[0]]

nRad=len(radiusArcsec[0])
nBg=len(innerArcsec)

bgStr=[]
for bg in range(nBg):
	bgStr.append('Background: %i"-%i"'%(int(innerArcsec[bg]),int(outerArcsec[bg])))

#print bgStr

#File columns
#0:srcRad
#1:bgInner
#2:bgOuter
#3:TimelineFlux
#4:TimelineErr
#5:MapFlux
#6:MapErr
#7:MapRgFlux
#8:MapRgErr
#9:SrcCorr
#10:ApCorr
#11:SrcFlux
#12:SrcErr
#13:ApFlux
#14:ApErr
#15:ApRgFlux
#16:ApRgErr

###Plot correction parameters
corrFact=Double4d(nRad,nBg,3,2)
##[srcRad , bgRad , band , [src|ap] ]
fluxArr=Double4d(nRad,nBg,3,6)
##[srcRad, bgRad , band , [timeline|map|mapRg|src|ap|apRg]
ob=0
myObsid=obsids[ob]
raSrc=raSrcs[ob]
decSrc=decSrcs[ob]
for b in range(3):
	band=bandStr[b]
	#print '%s Band'%band
	fileDat='%s0x%x_SrcFlux1run_%s.dat'%(dirPath,myObsid,band)
	fDat=open(fileDat,'r')
	lines=fDat.readlines()
	nLine=len(lines)
	for l in range(nLine):
		if lines[l].find('#') < 0:
			line=lines[l].split(',')
			iRad=radiusArcsec[b].index(float(line[0]))
			iBg=innerArcsec.index(float(line[1]))
			corrFact[iRad,iBg,b,0] = float(line[9]) #srcCorr
			corrFact[iRad,iBg,b,1] = float(line[10]) #apCorr
			fluxArr[iRad,iBg,b,0] = float(line[3]) #timeline
			fluxArr[iRad,iBg,b,1] = float(line[5]) #map
			fluxArr[iRad,iBg,b,2] = float(line[7]) #mapRg
			fluxArr[iRad,iBg,b,3] = float(line[11]) #srcFlux
			fluxArr[iRad,iBg,b,4] = float(line[13]) #apFlux
			fluxArr[iRad,iBg,b,5] = float(line[15]) #apRgFlux

##############################
##Plot correction factors
##############################

if pc:
    for b in range(3):
    	pCorr=PlotXY()
    	pCorr.batch = 1
    	ptStyles=[Style.SQUARE,Style.DIAMOND,Style.DCROSS,Style.CIRCLE]
    	lCols=[java.awt.Color.black,java.awt.Color.red,java.awt.Color.blue,java.awt.Color.green]
    	for bg in range(nBg):
    		#lSrcLine=LayerXY(Double1d(radiusArcsec[b][0:radMax[b][bg]]),Double1d(corrFact[0:radMax[b][bg],bg,b,0]))
    		#lSrcLine.style.stroke=1.0
    		#lSrcLine.line=Style.MARK_DASHED
    		#lSrcLine.symbol=Style.VCROSS
    		#lSrcLine.color=lCols[bg]
    		#lSrcLine.setInLegend(False)
    
    		lApLine=LayerXY(Double1d(radiusArcsec[b][0:radMax[b][bg]]),Double1d(corrFact[0:radMax[b][bg],bg,b,1]))
    		lApPt=LayerXY(Double1d(radiusArcsec[b][0:radMax[b][bg]]),Double1d(corrFact[0:radMax[b][bg],bg,b,1]))
    		lApLine.style.stroke=2.0
    		lApLine.symbol=Style.NONE
    		lApLine.color=lCols[bg]
    		lApLine.name=bgStr[bg]
    		
    		lApPt.line=Style.NONE
    		lApPt.style.stroke=2.0
    		lApPt.symbol=Style.DCROSS
    		#lApPt.symbol=ptStyles[bg]
    		lApPt.color=lCols[bg]
    		lApPt.setInLegend(False)
    		
    		#pCorr.addLayer(lSrcLine)	
    		pCorr.addLayer(lApLine)	
    		pCorr.addLayer(lApPt)
    
    
    	pCorr.batch = 0
    	pCorr.legend.visible=1
    	pCorr.legend.setColumns(2)
    	pCorr.legend.setBorderVisible(False)
    	#pCorr.legend.setTitleText('Background annulus radius [arcsec]')
    	#pCorr.legend.setPosition(PlotLegend.RIGHT)
    	#pCorr.legend.setLocation(0.9,0.5)
    	pCorr.xaxis.setTitleText('Source Aperture [arcsec]')
    	pCorr.yaxis.setTitleText('Correction Factor')
    	
    	pCorr.yaxis.setRange(1.,1.3)
    	pCorr.xaxis.setRange(20.,130.)
    	yr=pCorr.yaxis.getRange()
    	xr=pCorr.xaxis.getRange()
    
    
    	aBand=Annotation(0.98*xr[1],yr[0]+0.9*(yr[1]-yr[0]),bandStr[b])
    	aBand.setHalign(HAlign.LEFT)
    	aBand.setFontSize(12)
    	pCorr.addAnnotation(aBand)
    
    	plotFile='%s/Plots/ApCorr_%s.png'%(dirPath,bandStr[b])
    	pCorr.saveAsPNG(plotFile)

###########################
##Plot flux (1obsid)
###########################

if pf:
    for b in range(3):
    	pFlux=PlotXY()
    	pFlux.batch=1
    	ptStyles=[Style.SQUARE,Style.DIAMOND,Style.DCROSS,Style.CIRCLE]
    	lCols=[java.awt.Color.black,\
    		java.awt.Color.red,\
    		java.awt.Color.blue,\
    		java.awt.Color.green,\
    		java.awt.Color.orange,\
    		java.awt.Color.magenta]
    	for bg in range(nBg):
    		timeLine=LayerXY(Double1d(radiusArcsec[b][0:radMax[b][bg]]),Double1d(fluxArr[0:radMax[b][bg],bg,b,0]))
    		timeLine.line=Style.MARKED
    		timeLine.style.stroke=2.0
    		timeLine.name=bgStr[bg]
    		timeLine.color=lCols[0]
    		timeLine.symbol=ptStyles[bg]
    
    		pFlux.addLayer(timeLine)
    
    	for bg in range (nBg):
    		mapLine=LayerXY(Double1d(radiusArcsec[b][0:radMax[b][bg]]),Double1d(fluxArr[0:radMax[b][bg],bg,b,2]))
    		mapLine.line=Style.MARK_DASHED
    		if bg == 0:
    			mapLine.name='MapRg Flux'
    		else:
    			mapLine.setInLegend(False)
    		mapLine.color=lCols[1]
    		mapLine.symbol=ptStyles[bg]
    	
    		pFlux.addLayer(mapLine)
    
    		srcLine=LayerXY(Double1d(radiusArcsec[b][0:radMax[b][bg]]),Double1d(fluxArr[0:radMax[b][bg],bg,b,3]))
    		srcLine.line=Style.MARK_DASHED
    		if bg == 0:
    			srcLine.name='SrcCorr Flux'
    		else:
    			srcLine.setInLegend(False)
    		srcLine.color=lCols[2]
    		srcLine.symbol=ptStyles[bg]
    	
    		pFlux.addLayer(srcLine)
    
    		apLine=LayerXY(Double1d(radiusArcsec[b][0:radMax[b][bg]]),Double1d(fluxArr[0:radMax[b][bg],bg,b,5]))
    		apLine.line=Style.MARKED
    		if bg == 0:
    			apLine.name='ApCorr Flux'
    		else:
    			apLine.setInLegend(False)
    		apLine.color=lCols[3]
    		apLine.symbol=ptStyles[bg]
    	
    		pFlux.addLayer(apLine)
    	
    	pFlux.batch=0
    
    	pFlux.legend.visible=1
    	pFlux.legend.setColumns(2)
    	pFlux.legend.setBorderVisible(False)
    
    	#pCorr.yaxis.setRange(1.,1.3)
    	pFlux.xaxis.setRange(20.,130.)
    	yr=pFlux.yaxis.getRange()
    	xr=pFlux.xaxis.getRange()
    	aBand=Annotation(0.98*xr[1],yr[0]+0.9*(yr[1]-yr[0]),bandStr[b])
    	aBand.setHalign(HAlign.LEFT)
    	aBand.setFontSize(12)
    	pFlux.addAnnotation(aBand)
    	obsstr='0x%x'%myObsid
    	aObs=Annotation(0.98*xr[1],yr[0]+0.85*(yr[1]-yr[0]),obsstr)
    	aObs.setHalign(HAlign.LEFT)
    	aObs.setFontSize(8)
    	pFlux.addAnnotation(aObs)
    	pFlux.xaxis.setTitleText('Source Aperture [arcsec]')
    	pFlux.yaxis.setTitleText('Timeline Flux')

#######################################
## Read in all obsids
#######################################

fluxArrOb=Double5d(nob,nRad,nBg,3,6)
##[obs, srcRad, bgRad , band , [timeline|map|mapRg|src|ap|apRg]
CorrFact=Double4d(nRad,nBg,3,2)
##[srcRad, bgRad , band , [src|app]


maxFlux=0.
minFlux=999.

for ob in range(nob):
	myObsid=obsids[ob]
	raSrc=raSrcs[ob]
	decSrc=decSrcs[ob]
	for b in range(3):
		band=bandStr[b]
		#print '%s Band'%band
		fileDat='%s0x%x_SrcFlux1run_%s.dat'%(dirPath,myObsid,band)
		fDat=open(fileDat,'r')
		lines=fDat.readlines()
		nLine=len(lines)
		for l in range(nLine):
			if lines[l].find('#') < 0:
				line=lines[l].split(',')
				iRad=radiusArcsec[b].index(float(line[0]))
				iBg=innerArcsec.index(float(line[1]))
				corrFact[iRad,iBg,b,0]=float(line[9])
				corrFact[iRad,iBg,b,1]=float(line[10])
				fluxArrOb[ob,iRad,iBg,b,0] = float(line[3]) #timeline
				fluxArrOb[ob,iRad,iBg,b,1] = float(line[5]) #map
				fluxArrOb[ob,iRad,iBg,b,2] = float(line[7]) #mapRg
				fluxArrOb[ob,iRad,iBg,b,3] = float(line[11]) #srcFlux
				fluxArrOb[ob,iRad,iBg,b,4] = float(line[13]) #apFlux
				fluxArrOb[ob,iRad,iBg,b,5] = float(line[15]) #apRgFlux
				if min(fluxArrOb[ob,iRad,iBg,b,:]) < minFlux:
					minFlux=min(fluxArrOb[ob,iRad,iBg,b,:])
				if max(fluxArrOb[ob,iRad,iBg,b,:]) > maxFlux:
					maxFlux=max(fluxArrOb[ob,iRad,iBg,b,:])

		fDat.close()
#		for r in range(nRad):
#			for bg in range(nBg):
#				for i in range(6):
#					if fluxArrOb[ob,r,bg,b,i]==0:
#						fluxArrOb[ob,r,bg,b,i]=None

nbin=30
#minFlux=0.
#minFlux=1.
minFlux=Double1d(3)
maxFlux=Double1d(3)
binsize=Double1d(3)

bins=Double3d(nbin,3,3,0) #bin , band , [min|max|centre]
hist=Double3d(nbin,3,6,0) #bin , band , [time|map|mapRg|src|ap|apRg]

lCols=[java.awt.Color.black,\
    		java.awt.Color.red,\
    		java.awt.Color.blue,\
    		java.awt.Color.green,\
    		java.awt.Color.orange,\
    		java.awt.Color.magenta]

plotX=Double2d(nbin*2,3) #bin , band
plotY=Double3d(nbin*2,3,6) #bin , band , [time|map|mapRg|src|ap|apRg]

medians=Double2d(3,6) #band , [time|map|mapRg|src|ap|apRg]
means=Double2d(3,6) #band , [time|map|mapRg|src|ap|apRg]
stdevs=Double2d(3,6) #band , [time|map|mapRg|src|ap|apRg]

fileSum='%s%s_r%i_bg%i-%i_Summary.dat'%\
	(dirPath,fileObsName[:-4],int(radiusArcsec[b][r]),int(innerArcsec[r]),int(outerArcsec[r]))
fSum=open(fileSum,'w')
#binsize=1.
for b in range(3):
	minFlux[b]=0.95*min(fluxArrOb[:,:,:,b,:])
	maxFlux[b]=1.05*max(fluxArrOb[:,:,:,b,:])
	binsize[b]=(maxFlux[b]-minFlux[b])/nbin
	#for r in range(nRad):
	#for bg in range(nBg):
	r=0
	bg=0

	for f in range(6):
		medians[b,f] = MEDIAN(fluxArrOb[:,r,bg,b,f])
		means[b,f] = MEAN(fluxArrOb[:,r,bg,b,f])
		stdevs[b,f] = STDDEV(fluxArrOb[:,r,bg,b,f])

	line='%s:\n'%bandStr[b]
	print line[:-2]
	fSum.write(line)

	line='  Src Corr Fact:     %.4f\n'%corrFact[0,0,b,0]
	print line[:-2]
	fSum.write(line)
	line='  Ap Corr Fact:      %.4f\n'%corrFact[0,0,b,1]
	print line[:-2]
	fSum.write(line)
	line='  Timeline Flux:     %.4f +/- %.4f Jy\n'%(means[b,0],stdevs[b,0])
	print line[:-2]
	fSum.write(line)
	line='  Map Flux:          %.4f +/- %.4f Jy\n'%(means[b,1],stdevs[b,1])
	print line[:-2]
	fSum.write(line)
	line='  Map Flux (RG):     %.4f +/- %.4f Jy\n'%(means[b,2],stdevs[b,2])
	print line[:-2]
	fSum.write(line)
	line='  Src corr Flux:     %.4f +/- %.4f Jy\n'%(means[b,3],stdevs[b,3])
	print line[:-2]
	fSum.write(line)
	line='  Ap corr Flux:      %.4f +/- %.4f Jy\n'%(means[b,4],stdevs[b,4])
	print line[:-2]
	fSum.write(line)
	line='  Ap corr Flux (RG): %.4f +/- %.4f Jy\n'%(means[b,5],stdevs[b,5])
	print line[:-2]
	fSum.write(line)

	for bin in range(nbin):
		bins[bin,b,0]=minFlux[b] + bin*binsize[b]
		bins[bin,b,1]=minFlux[b] + (bin+1.)*binsize[b]
		bins[bin,b,2]=minFlux[b] + (bin+0.5)*binsize[b]
		for ob in range(nob):
			for f in range(6):
				if abs(fluxArrOb[ob,r,bg,b,f] - bins[bin,b,2]) < binsize[b]/2.:
					hist[bin,b,f]=hist[bin,b,f]+1
		plotX[2*bin,b]=bins[bin,b,0]
		plotX[2*bin+1,b]=bins[bin,b,1]
		for f in range(6):
			plotY[2*bin,b,f]=hist[bin,b,f]
			plotY[2*bin+1,b,f]=hist[bin,b,f]
	pH=PlotXY()

	lTime=LayerXY(plotX[:,b],plotY[:,b,0])
	lTime.name='Timeline'
	lTime.style.stroke=3.0
	lTime.color=lCols[0]
	pH.addLayer(lTime)

	lMap=LayerXY(plotX[:,b],plotY[:,b,1])
	lMap.name='Map'
	lMap.style.stroke=1.0
	lMap.line=Style.DASHED
	lMap.color=lCols[1]
	pH.addLayer(lMap)

	lMapRg=LayerXY(plotX[:,b],plotY[:,b,2])
	lMapRg.name='Map (RG)'
	lMapRg.style.stroke=2.0
	lMapRg.color=lCols[1]
	pH.addLayer(lMapRg)

	lSrc=LayerXY(plotX[:,b],plotY[:,b,3])
	lSrc.name='Src Corr only'
	lSrc.color=lCols[2]
	pH.addLayer(lSrc)

	lAp=LayerXY(plotX[:,b],plotY[:,b,4])
	lAp.name='Ap Corr'
	lAp.style.stroke=1.0
	lAp.line=Style.DASHED
	lAp.color=lCols[3]
	pH.addLayer(lAp)

	lApRg=LayerXY(plotX[:,b],plotY[:,b,5])
	lApRg.name='Ap Corr (RG)'
	lApRg.color=lCols[3]
	lApRg.style.stroke=2.0
	pH.addLayer(lApRg)

	yr=pH.yaxis.getRange()
	pH.yaxis.setRange(0.,yr[1])
	yr=pH.yaxis.getRange()
	xr=pH.xaxis.getRange()
	aBand=Annotation(xr[0]+0.98*(xr[1]-xr[0]),yr[0]+0.9*(yr[1]-yr[0]),bandStr[b])
	aBand.setHalign(HAlign.LEFT)
	aBand.setFontSize(12)
	pH.addAnnotation(aBand)

	aTar=Annotation(xr[0]+0.02*(xr[1]-xr[0]),yr[0]+0.9*(yr[1]-yr[0]),tarName)
	aTar.setHalign(HAlign.RIGHT)
	aTar.setFontSize(12)
	pH.addAnnotation(aTar)

	aNob=Annotation(xr[0]+0.02*(xr[1]-xr[0]),yr[0]+0.85*(yr[1]-yr[0]),'%i obs'%nob)
	aNob.setHalign(HAlign.RIGHT)
	aNob.setFontSize(8)
	pH.addAnnotation(aNob)

	medTimeStr='Timeline: %.4f +/- %.4f Jy'%(means[b,0],stdevs[b,0])
	aMedTime=Annotation(xr[0]+0.98*(xr[1]-xr[0]),yr[0]+0.6*(yr[1]-yr[0]),medTimeStr)
	aMedTime.setHalign(HAlign.LEFT)
	aMedTime.setFontSize(8)
	#pH.addAnnotation(aMedTime)

	medApStr='ApCorr: %.4f +/- %.4f Jy'%(means[b,0],stdevs[b,0])
	aApTime=Annotation(xr[0]+0.98*(xr[1]-xr[0]),yr[0]+0.55*(yr[1]-yr[0]),medTimeStr)
	aApTime.setHalign(HAlign.LEFT)
	aApTime.setFontSize(8)
	aApTime.setColor(java.awt.Color.green)
	#pH.addAnnotation(aApTime)

	pH.legend.visible=1
	pH.xaxis.setTitleText('Source Flux [Jy]')
	pH.yaxis.setTitleText('# Obs')

	plotFile='%s%s_r%i_bg%i-%i_%s.png'%\
		(dirPlot,fileObsName[:-4],int(radiusArcsec[b][r]),int(innerArcsec[r]),int(outerArcsec[r]),bandStr[b])

	pH.saveAsPNG(plotFile)
	
fSum.close()
	#pH.style.setChartType(Style.HISTOGRAM)

