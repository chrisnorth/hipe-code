## Run timeline fitter and aperture correction
from os.path import isfile
import sys
sys.path.append("/home/astrog82/spxcen/Herschel/Calibration/RelGains")
import apcorrfn

dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/'
dirPathObs = dirPath + 'Obs/'
files=[\
#	'AlphaBoo_obs_LargeScan.dat', 'AlphaBoo_obs_SmallScan.dat', 'AlphaBoo_obs_ParallelMode.dat',\
#	'GammaDra_obs_LargeScan.dat', 'GammaDra_obs_SmallScan.dat', 'GammaDra_obs_ParallelMode.dat',\
#	'AlphaTau_obs_LargeScan.dat', 'AlphaTau_obs_SmallScan.dat',\
#	'BetaPeg_obs_LargeScan.dat', 'BetaPeg_obs_SmallScan.dat',\
#	'BetaUmi_obs_LargeScan.dat', 'BetaUmi_obs_SmallScan.dat', 'BetaUmi_obs_ParallelMode.dat',\
#	'Neptune_obs_LargeScan.dat', 'Neptune_obs_SmallScan.dat',\
#	'Uranus_obs_LargeScan.dat', 'Uranus_obs_SmallScan.dat',\
#	'1_Ceres_obs_LargeScan.dat', '1_Ceres_obs_SmallScan.dat',\
#	'2_Pallas_obs_LargeScan.dat', '2_Pallas_obs_SmallScan.dat',\
#	'3_Juno_obs_LargeScan.dat', '3_Juno_obs_SmallScan.dat',\
#	'4_Vesta_obs_LargeScan.dat', '4_Vesta_obs_SmallScan.dat',\
#	'10_Hygeia_obs_LargeScan.dat', '10_Hygeia_obs_SmallScan.dat'\
#
#
#
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

##Get list of ObsIDs
#fileObsList=dirPath+'Obs/GammaDra_obs_All.dat'
#fileDb=dirPath+'Obs/AlphaBoo_ApCorr_LargeScan.csv'

#raSrc=213.91175352734098
#decSrc=19.176395668750974

bandStr=["PSW","PMW","PLW"]
beamAreasArcsec=[433.1,776.9,1628.3]

##Set ApPhot radii
#radiusArcsec=[[22.,30.,45.,50.,60.,100.],\
#	[30.,45.,55.,65.,80.,100.],\
#	[42.,50.,60.,90.,100.,120.]]
#innerArcsec=[60.,100.,200.,300.]
#outerArcsec=[90.,150.,250.,350.]

radiusArcsec=[[22.],[30.],[42.]]
innerArcsec=[60.]
outerArcsec=[90.]

myDataPool="default"
#fOb=open(fileObsList,'r')
#obsids_in=String1d(fOb.readlines())

#nob=len(obsids_in)
#obsids=Int1d(nob)
#raSrcs=Double1d(nob)
#decSrcs=Double1d(nob)
#for ob in range(nob):
#	obsids[ob]=int(obsids_in[ob],16)
#	raSrcs[ob]=raSrc
#	decSrcs[ob]=decSrc

nFile=len(files)
obsids=[]
obsidsX=[]
names=[]
modes=[]
raSrcs=[]
decSrcs=[]
for file in files:
	fileObsList=dirPathObs+file
	print file
	fOb=open(fileObsList,'r')
	obsids_in=String1d(fOb.readlines())
	nob=len(obsids_in)

	nobIn=len(obsids_in)
	tarName='UNKNOWN'
	tarMode='UNKNOWN'
	tarRA=-1
	tarDec=-1
	for ob in range(nobIn):
		if obsids_in[ob].find('#') >= 0:
			if obsids_in[ob].find('Name:') > 0:
				tarName=obsids_in[ob].split(':')[1][1:-1]
			if obsids_in[ob].find('Mode:') > 0:
				tarMode=obsids_in[ob].split(':')[1][1:-1]
			if obsids_in[ob].find('RA:') > 0:
				tarRA=obsids_in[ob].split(':')[1][1:-1]
			if obsids_in[ob].find('Dec:') > 0:
				tarDec=obsids_in[ob].split(':')[1][1:-1]

	for ob in range(nobIn):
		obRA=-1
		obDec=-1
		if obsids_in[ob].find('#') < 0:
			obsids.append(int(obsids_in[ob].split(',')[0],16))
			if len(obsids_in[ob].split(',')) > 1:
				obRA=float(obsids_in[ob].split(',')[1])
				if obRA == "Variable":
					obRA=tarRA
				obDec=float(obsids_in[ob].split(',')[2])
				if obDec == "Variable":
					obDec=tarDec
			obsidsX.append(obsids_in[ob].split(',')[0])
			names.append(tarName)
			modes.append(tarMode)
			raSrcs.append(obRA)
			decSrcs.append(obDec)

obsids=Int1d(obsids)
names=String1d(names)
modes=String1d(modes)
obsidsX=String1d(obsidsX)
raSrcs=Double1d(raSrcs)
decSrcs=Double1d(decSrcs)
nob=len(obsids)

for ob in range(nob):
	print '0x%x %s (%s): [%.3f, %.3f]'%(obsids[ob],names[ob],modes[ob],float(raSrcs[ob]),float(decSrcs[ob]))

#obsids=[]
#raSrcs=[]
#decSrcs=[]
#nobIn=len(obsids_in)
#for ob in range(nobIn):
#	line=obsids_in[ob]
#	if line.find('#') < 0:
#		obsids.append(int(line.split(',')[0],16))
#		raSrcs.append(float(line.split(',')[1]))
#		decSrcs.append(float(line.split(',')[2]))
#obsids=Int1d(obsids)
#raSrcs=Double1d(raSrcs)
#decSrcs=Double1d(decSrcs)
#nob=len(obsids)

## Read in existing database
#if isfile(fileDb):
#	#don't do anything - just overwrite the file
#	print 'Output file exists. Overwriting.'
#else:
#	print 'Output file does not exist. Making a new one.'

##read in L1 data

for ob in range(nob):
	myObsid=obsids[ob]
	#raSrc=raSrcs[ob]
	#decSrc=decSrcs[ob]

	obs=getObservation(myObsid,poolName=myDataPool)
	print "Processing observation %i (0x%X) [%s] from data pool %s."%(myObsid, myObsid, names[ob],myDataPool)

	try:
		decSrc=obs.refs["level2"].meta["decNominal"].value
		raSrc=obs.refs["level2"].meta["raNominal"].value
	except:
		raSrc=raSrcs[ob]
		decSrc=decSrcs[ob]

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
		fileOut='%sObs/0x%x_SrcFlux1run_%s.dat'%(dirPath,myObsid,band)
		fOut=open(fileOut,'w')
		line='#srcRad,bgInner,bgOuter,TimelineFlux,TimelineErr,MapFlux,MapErr,MapRgFlux,MapRgErr,SrcCorr,ApCorr,SrcFlux,SrcErr,ApFlux,ApErr,ApRgFlux,ApRgErr\n'
		fOut.write(line)
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

					## Timeline fitter
					try:
						fitter=TimelineSourceFitterTask()
						outputTimeline=fitter(input=scansContext,array=band,\
							sourcePositionEstimate=Double1d([raSrc,decSrc]),\
							rPeak=srcRad,rBackground=Double1d([bgInner,bgOuter]),\
							fitEllipticalGaussian=True,\
							useBackInFit=True,allowVaryBackground=True)
						timelineFlux=outputTimeline.meta["fitFlux"].getValue()
						timelineErr=outputTimeline.meta["fitFluxErr"].getValue()
						print 'Timeline flux: %.5f +/- %.5f'%(timelineFlux,timelineErr)
						#sys.exit()

						## Get real coordinates
						raSrc=outputTimeline.meta["fitRa"].getValue()
						decSrc=outputTimeline.meta["fitDec"].getValue()

					except:
						timelineFlux=-999.
						timelineErr=-999.
						print 'Timeline flux: ERROR'

					try:
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
					except:
						mapFlux=-999.
						mapErr=-999.
						mapRgFlux=-999.
						mapRgErr=-999.
						print 'Map Flux: ERROR'
						srcCorr=-999.
						apCorr=-999.
						bgCorr=-999.
						print 'Aperture correction: ERROR'
						srcFlux=-999
						srcErr=-999.
						srcRgFlux=-999.
						srcRgErr=-999.
						apFlux=-999
						apErr=-999.
						apRgFlux=-999.
						apRgErr=-999.
						print 'SrcCorr Flux: ERROR'
						print 'Aperture correction: ERROR'
					print'----'

					#srcPos=Double1d[raSrc,decSrc]
					line='%.1f, %.1f, %.1f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n'%\
						(srcRad,bgInner,bgOuter,timelineFlux,timelineErr,\
						mapFlux,mapErr,mapRgFlux,mapRgErr,srcCorr,apCorr,\
						srcFlux,srcErr,apFlux,apErr,apRgFlux,apRgErr)
					fOut.write(line)
		fOut.close()
	#sys.exit()