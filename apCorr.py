## Run timeline fitter and aperture correction
from os.path import isfile
import sys
sys.path.append("/home/astrog82/spxcen/Herschel/Calibration/RelGains")
import apcorrfn

dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/'

##Get list of ObsIDs
fileObsList=dirPath+'Obs/GammaDra_obs_All.dat'
#fileDb=dirPath+'Obs/AlphaBoo_ApCorr_LargeScan.csv'

#raSrc=213.91175352734098
#decSrc=19.176395668750974

bandStr=["PSW","PMW","PLW"]
beamAreasArcsec=[433.1,776.9,1628.3]

##Set ApPhot radii
radiusArcsec=[[22.,30.,45.,50.,60.,100.],\
	[30.,45.,55.,65.,80.,100.],\
	[42.,50.,60.,90.,100.,120.]]

innerArcsec=[60.,100.,200.,300.]
outerArcsec=[90.,150.,250.,350.]

myDataPool="default"
fOb=open(fileObsList,'r')
obsids_in=String1d(fOb.readlines())

#nob=len(obsids_in)
#obsids=Int1d(nob)
#raSrcs=Double1d(nob)
#decSrcs=Double1d(nob)
#for ob in range(nob):
#	obsids[ob]=int(obsids_in[ob],16)
#	raSrcs[ob]=raSrc
#	decSrcs[ob]=decSrc

obsids=[]
raSrcs=[]
decSrcs=[]
nobIn=len(obsids_in)
for ob in range(nobIn):
	line=obsids_in[ob]
	if line.find('#') < 0:
		obsids.append(int(line.split(',')[0],16))
		raSrcs.append(float(line.split(',')[1]))
		decSrcs.append(float(line.split(',')[2]))
obsids=Int1d(obsids)
raSrcs=Double1d(raSrcs)
decSrcs=Double1d(decSrcs)
nob=len(obsids)

## Read in existing database
#if isfile(fileDb):
#	#don't do anything - just overwrite the file
#	print 'Output file exists. Overwriting.'
#else:
#	print 'Output file does not exist. Making a new one.'

##read in L1 data

for ob in range(nob):
	myObsid=obsids[ob]
	raSrc=raSrcs[ob]
	decSrc=decSrcs[ob]

	obs=getObservation(myObsid,poolName=myDataPool)
	print "Processing observation %i (0x%X) from data pool %s."%(myObsid, myObsid, myDataPool)

	scans = obs.level1
	
	# Remove baselines
	#scans=baselineRemovalMedian(scans)
	
	scansInput = []
	for i in range(scans.count):
		#print 'scan %d/%d'%(i+1,scans.count)
		tempPSP=scans.getProduct(i)
		tempPSP.setType('PST')
		scansInput.append(tempPSP)

	scansContext = ScanContext(scansInput)

	for b in range(3):
		band=bandStr[b]
		print '%s Band'%band
		fileOut='%sObs/0x%x_SrcFlux_%s.dat'%(dirPath,myObsid,band)
		fOut=open(fileOut,'w')
		line='#srcRad,bgInner,bgOuter,TimelineFlux,TimelineErr,MapFlux,MapErr,MapRgFlux,MapRgErr,SrcCorr,ApCorr,SrcFlux,SrcErr,ApFlux,ApErr,ApRgFlux,ApRgErr\n'
		mapId=band
		mapIdRg=band+'Rg'

		##Read in maps
		mapIn = obs.level2.refs[mapId].product
		mapInRg = obs.level2.refs[mapIdRg].product

		beamAreaArcsec=beamAreasArcsec[b]
		##Pixel sizes (arcsec)
		pixelArcsec=(mapIn.wcs.cdelt1*3600.)**2
		##Beam areas (pixels)
		beamAreaPixel=beamAreaArcsec/pixelArcsec

		##Convert map units
		mapInArcsec = convertImageUnit(image=mapIn,newUnit='Jy/arcsec^2',beamArea=beamAreaPixel)
		mapInRgArcsec = convertImageUnit(image=mapInRg,newUnit='Jy/arcsec^2',beamArea=beamAreaPixel)
		
		nRad=len(radiusArcsec[b])
		nBg=len(innerArcsec)

		for r in range(nRad):
			srcRad=radiusArcsec[b][r]
			for bg in range(nBg):
				bgInner=innerArcsec[bg]
				bgOuter=outerArcsec[bg]
				if srcRad <= bgInner:
					print 'Source: %.1f ; Background: [%.1f,%.1f]'%(srcRad,bgInner,bgOuter)
					##Aperture photometry
					result = annularSkyAperturePhotometry(image=mapInArcsec, fractional=1, \
						centerRA=str(raSrc), centerDec=str(decSrc), \
						radiusArcsec=srcRad, \
						innerArcsec=bgInner, outerArcsec=bgOuter)

					resultRg = annularSkyAperturePhotometry(image=mapInRgArcsec, fractional=1, \
						centerRA=str(raSrc), centerDec=str(decSrc), \
						radiusArcsec=srcRad, \
						innerArcsec=bgInner, outerArcsec=bgOuter)

					mapFlux=result["Results table"][0].data[2]
					mapErr=result["Results table"][2].data[2]
					mapRgFlux=resultRg["Results table"][0].data[2]
					mapRgErr=resultRg["Results table"][2].data[2]

					(srcCorr,bgCorr,apCorr)=apcorrfn.getApCorr(srcRad,bgInner,bgOuter,band=band)
					print 'Map Flux: %.5f +/- %.5f ; %.5f +/- %.5f (RG)'%(mapFlux,mapErr,mapRgFlux,mapRgErr)
					print 'Aperture correction: %.5f, %.5f, %.5f'%(srcCorr,bgCorr,apCorr)
					srcFlux=mapFlux*srcCorr
					srcErr=mapErr*srcCorr
					srcRgFlux=mapRgFlux*srcCorr
					srcRgErr=mapRgErr*srcCorr
					apFlux=mapFlux*apCorr
					apErr=mapErr*apCorr
					apRgFlux=mapRgFlux*apCorr
					apRgErr=mapRgErr*apCorr
					print 'SrcCorr Flux: %.5f +/- %.5f ; %.5f +/- %.5f (RG)'%(srcFlux,srcErr,srcRgFlux,srcRgErr)
					print 'ApCorr Flux:  %.5f +/- %.5f ; %.5f +/- %.5f (RG)'%(apFlux,apErr,apRgFlux,apRgErr)

					#srcPos=Double1d[raSrc,decSrc]
					fitter=TimelineSourceFitterTask()
					outputTimeline=fitter(input=scansContext,array=band,\
						sourcePositionEstimate=Double1d([raSrc,decSrc]),\
						rPeak=srcRad,rBackground=Double1d([bgInner,bgOuter]),\
						fitEllipticalGaussian=True,\
						useBackInFit=True,allowVaryBackground=True)
					#sys.exit()

					timelineFlux=outputTimeline.meta["fitFlux"].getValue()
					timelineErr=outputTimeline.meta["fitFluxErr"].getValue()
					print 'Timeline flux: %.5f +/- %.5f'%(timelineFlux,timelineErr)
					print'----'
					line='%.1f, %.1f, %.1f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n'%\
						(srcRad,bgInner,bgOuter,timelineFlux,timelineErr,\
						mapFlux,mapErr,mapRgFlux,mapRgErr,srcCorr,apCorr,\
						srcFlux,srcErr,apFlux,apErr,apRgFlux,apRgErr)
					fOut.write(line)
		fOut.close()
