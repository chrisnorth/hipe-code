import sys
from herschel.ia.gui.plot.renderer.PComponentEngine import HAlign
from herschel.ia.gui.plot.renderer.PComponentEngine import VAlign

from herschel.ia.numeric.toolbox.basic import Histogram
from herschel.ia.numeric.toolbox.basic import BinCentres

## Plot aperture photometry stuff

pc=True
pf=True


dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/'

##Get list of ObsIDs
fileObsList=dirPath+'Obs/AlphaBoo_obs_All.dat'
tarName='AlphaBoo'
#fileObsList=dirPath+'Obs/GammaDra_obs_All.dat'
#tarName='GammaDra'

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

radiusArcsec=[[22.,30.,45.,50.,60.,100.],\
	[30.,45.,55.,65.,80.,100.],\
	[42.,50.,60.,90.,100.,120.]]

#Max source radii that are valid
radMax=[[5,6,6,6],[3,6,6,6],[3,4,6,6]]

innerArcsec=[60.,100.,200.,300.]
outerArcsec=[90.,150.,250.,350.]

#Min BR radii that are valid
bgMin=[[0,0,0,0,0,1],[0,0,0,1,1,1],[0,0,0,1,1,2]]

nRad=len(radiusArcsec[0])
nBg=len(innerArcsec)

bgStr=[]
for bg in range(nBg):
	bgStr.append('Background: %i"-%i"'%(int(innerArcsec[bg]),int(outerArcsec[bg])))

print bgStr

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
	print '%s Band'%band
	fileDat='%sObs/0x%x_SrcFlux_%s.dat'%(dirPath,myObsid,band)
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
    
    	plotFile='%s/Plots/ApCorr_%s_%s.png'%(dirPath,tarName,bandStr[b])
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
    		mapLine.style.stroke=2.0
    		mapLine.symbol=ptStyles[bg]
    	
    		pFlux.addLayer(mapLine)
    
    		srcLine=LayerXY(Double1d(radiusArcsec[b][0:radMax[b][bg]]),Double1d(fluxArr[0:radMax[b][bg],bg,b,3]))
    		srcLine.line=Style.MARK_DASHED
    		srcLine.style.stroke=2.0
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
    		apLine.style.stroke=2.0
    		apLine.symbol=ptStyles[bg]
    	
    		pFlux.addLayer(apLine)
    	
    	pFlux.batch=0
    
    	pFlux.legend.visible=1
    	pFlux.legend.setColumns(2)
    	pFlux.legend.setBorderVisible(False)
    
    	#pCorr.yaxis.setRange(1.,1.3)
    	pFlux.xaxis.setRange(20.,130.)
    	yr=pFlux.yaxis.getRange()
    	pFlux.yaxis.setRange(yr[0],yr[0]+1.05*(yr[1]-yr[0]))
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

    	plotFile='%s/Plots/ApCorrFlux_%s_%s.png'%(dirPath,tarName,bandStr[b])
    	pFlux.saveAsPNG(plotFile)

#print hdr
