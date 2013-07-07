import sys
from herschel.ia.gui.plot.renderer.PComponentEngine import HAlign
from herschel.ia.gui.plot.renderer.PComponentEngine import VAlign

from herschel.ia.numeric.toolbox.basic import Histogram
from herschel.ia.numeric.toolbox.basic import BinCentres

## Plot aperture photometry stuff

pc=False
pf=False

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

dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/Obs/'
dirPlot = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/Plots/'
##Get list of ObsIDs
files=[\
	'10_Hygeia_obs_All.dat',\
	'173_Ino_obs_All.dat',\
	'19_Fortuna_obs_All.dat',\
	'1_Ceres_obs_All.dat',\
	'20_Massalia_obs_All.dat',\
	'21_Lutetia_obs_All.dat',\
	'253_Mathilde_obs_All.dat',\
	'29_Amphitrite_obs_All.dat',\
	'2_Pallas_obs_All.dat',\
	'372_Palma_obs_All.dat',\
	'37_Fides_obs_All.dat',\
	'3_Juno_obs_All.dat',\
	'40_Harmonia_obs_All.dat',\
	'47_Aglaja_obs_All.dat',\
	'4_Vesta_obs_All.dat',\
	'511_Davida_obs_All.dat',\
	'52_Europa_obs_All.dat',\
	'54_Alexandra_obs_All.dat',\
	'65_Cybele_obs_All.dat',\
	'6_Hebe_obs_All.dat',\
	'7_Iris_obs_All.dat',\
	'88_Thisbe_obs_All.dat',\
	'8_Flora_obs_All.dat',\
	'93_Minerva_obs_All.dat',\
	]

nFiles=len(files)
tarNames=[]
nObsFiles=Int1d(nFiles)
obsIds=[]

for f in range(nFiles):

	fileObsList=dirPath+files[f]

	bandStr=["PSW","PMW","PLW"]

	fObs=open(fileObsList,'r')
	obsidsIn=fObs.readlines()
	fObs.close()
	obsids=[]
	raSrcs=[]
	decSrcs=[]
	nobIn=len(obsidsIn)
	nameFound=False
	for ob in range(nobIn):
		line=obsidsIn[ob]
		if line.find('#') < 0:
			obsids.append(int(line.split(',')[0],16))
			raSrcs.append(float(line.split(',')[1]))
			decSrcs.append(float(line.split(',')[2]))
		else:
			if line.find('Name:') > 0:
				tarNames.append(line.split(':')[1][1:-1])
				nameFound=True
	if nameFound==False:
		tarNames.append['UNKNOWN']
	
	obsIds.append(Int1d(obsids))
	#raSrcs=Double1d(raSrcs)
	#decSrcs=Double1d(decSrcs)
	nob=len(obsids)
	nObsFiles[f]=nob

	tarNames[f]='%s (%d)'%(tarNames[f],nObsFiles[f])
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
	

	#######################################
	## Read in all obsids
	#######################################

fluxMeanFiles=Double5d(nFiles,nRad,nBg,3,6)
fluxRelFiles=Double5d(nFiles,nRad,nBg,3,6)
fluxErrFiles=Double5d(nFiles,nRad,nBg,3,6)
fluxMeanAll=Double4d(nRad,nBg,3,6)
fluxErrAll=Double4d(nRad,nBg,3,6)
##[object, srcRad, bgRad , band , [timeline|map|mapRg|src|ap|apRg]
#CorrFact=Double4d(nRad,nBg,3,2)
###[srcRad, bgRad , band , [src|app]

maxFlux=0.
minFlux=999.

for f in range(nFiles):
	nObs=nObsFiles[f]
	obsList=obsIds[f]
	fluxArrOb=Double5d(nObs,nRad,nBg,3,6)
	fluxErrOb=Double5d(nObs,nRad,nBg,3,6)
	##[object, srcRad, bgRad , band , [timeline|map|mapRg|src|ap|apRg]
	for b in range(3):
		band=bandStr[b]
		for ob in range(nObs):
			myObsid=obsList[ob]
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

		print iRad
		print iBg
		for r in [iRad]:
			for bg in [iBg]:
				for p in range(6):
					fluxMeanFiles[f,r,bg,b,p]=MEAN(fluxArrOb[:,r,bg,b,p])
					print 'File',f,'param',p,'MEAN:',MEAN(fluxArrOb[:,r,bg,b,p])
					fluxErrFiles[f,r,bg,b,p]=STDDEV(fluxArrOb[:,r,bg,b,p])
					print 'File',f,'param',p,'STDDEV:',STDDEV(fluxArrOb[:,r,bg,b,p])
					fluxRelFiles[f,r,bg,b,p]=fluxMeanFiles[f,r,bg,b,p]/fluxMeanFiles[f,r,bg,b,0]

for b in range(3):
	for r in [iRad]:
		for bg in [iBg]:
			for p in range(6):
				fluxMeanAll[r,bg,b,p]=MEAN(fluxMeanFiles[:,r,bg,b,p])
				fluxErrAll[r,bg,b,p]=STDDEV(fluxMeanFiles[:,r,bg,b,p])

rFiles=Float1d(range(nFiles)) + 0.5

#pFluxPSW=PlotXY()
#lTime=LayerXY(rFiles,fluxMeanFiles[:,0,0,0,0])
#lTime.line=Style.NONE
#lTime.name='Timeline'
#lTime.symbol=Style.DCROSS
#pFluxPSW.addLayer(lTime)
#
#lMap=LayerXY(rFiles,fluxMeanFiles[:,0,0,0,1])
#lMap.line=Style.NONE
#lMap.name='Map'
#lMap.symbol=Style.VCROSS
#pFluxPSW.addLayer(lMap)
#
#lAp=LayerXY(rFiles,fluxMeanFiles[:,0,0,0,5])
#lAp.line=Style.NONE
#lAp.name='ApCorr'
#lAp.symbol=Style.DIAMOND
#pFluxPSW.addLayer(lAp)
#
#pFluxPSW.legend.visible=1

pErrPSW=PlotXY()

lTimeRel=LayerXY(rFiles,fluxRelFiles[:,0,0,0,0])
lTimeRel.line=Style.NONE
lTimeRel.name='Timeline'
lTimeRel.symbol=Style.DCROSS
lTimeRel.setErrorY(Double1d(fluxErrFiles[:,0,0,0,0]),Double1d(fluxErrFiles[:,0,0,0,0]))
lTimeRel.setColor(java.awt.Color.black)
lTimeRel.style.stroke=2.
pErrPSW.addLayer(lTimeRel)

#lTimeRelUp=LayerXY(rFiles-0.2,fluxRelFiles[:,0,0,0,0]+fluxErrFiles[:,0,0,0,0])
#lTimeRelUp.setColor(java.awt.Color.black)
#pErrPSW.addLayer(lTimeRelUp)
#lTimeRelDown=LayerXY(rFiles-0.2,fluxRelFiles[:,0,0,0,0]-fluxErrFiles[:,0,0,0,0])
#lTimeRelDown.setColor(java.awt.Color.black)
#pErrPSW.addLayer(lTimeRelDown)

lMapRel=LayerXY(rFiles-0.2,fluxRelFiles[:,0,0,0,1])
lMapRel.line=Style.NONE
lMapRel.name='Map'
lMapRel.symbol=Style.VCROSS
#lMapRel.setErrorY(Double1d(fluxErrFiles[:,0,0,0,1]),Double1d(fluxErrFiles[:,0,0,0,1]))
lMapRel.setColor(java.awt.Color.red)
lMapRel.style.stroke=2.
pErrPSW.addLayer(lMapRel)

lApRel=LayerXY(rFiles+0.2,fluxRelFiles[:,0,0,0,5])
lApRel.line=Style.NONE
lApRel.name='Corrected'
lApRel.symbol=Style.DIAMOND
lApRel.setErrorY(Double1d(fluxErrFiles[:,0,0,0,5]),Double1d(fluxErrFiles[:,0,0,0,5]))
lApRel.setColor(java.awt.Color.blue)
lApRel.style.stroke=2.
pErrPSW.addLayer(lApRel)

#lApRelUp=LayerXY(rFiles-0.2,fluxRelFiles[:,0,0,0,5]+fluxErrFiles[:,0,0,0,5])
#lApRelUp.setColor(java.awt.Color.blue)
#pErrPSW.addLayer(lApRelUp)
#lApRelDown=LayerXY(rFiles-0.2,fluxRelFiles[:,0,0,0,5]-fluxErrFiles[:,0,0,0,5])
#lApRelDown.setColor(java.awt.Color.blue)
#pErrPSW.addLayer(lApRelDown)

pErrPSW.yaxis.setRange(0,2)
pErrPSW.yaxis.setTitleText('Flux relative to timeline flux density')
pErrPSW.legend.visible=1

pErrPSW.xaxis.setRange(0,nFiles)
pErrPSW.xaxis.tick.setFixedValues(rFiles)
pErrPSW.xaxis.tick.label.setFixedStrings(tarNames)
pErrPSW.xaxis.tick.label.setOrientation(1)
pErrPSW.xaxis.setTitleText('Object (# Obs)')