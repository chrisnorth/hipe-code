
from herschel.ia.pal.pool.http import HttpClientPool
import datetime
import sys

dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/Obs/'

pool = HttpClientPool ("http://wakefield.bnsc.rl.ac.uk/hcss/pal", "ops")

#pool.setAuthentication(Configuration.getProperty("hcss.access.authentication.authstring"))
store = ProductStorage (pool)	# Just the Data
store.authenticate()

#### INPUT FILE(S)
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

#### SET TO OVERWRITE OBS
overwrite=False

nFile=len(files)
obsids=[]
obsidsX=[]
names=[]
modes=[]
RAs=[]
Decs=[]
for file in files:
	fileObsList=dirPath+file
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
			RAs.append(obRA)
			Decs.append(obDec)

obsids=Int1d(obsids)
names=String1d(names)
modes=String1d(modes)
obsidsX=String1d(obsidsX)
RAs=Double1d(RAs)
Decs=Double1d(Decs)
nob=len(obsids)

for ob in range(nob):
	print '0x%x %s (%s): [%.3f, %.3f]'%(obsids[ob],names[ob],modes[ob],float(RAs[ob]),float(Decs[ob]))

#sys.exit()
fileLog=dirPath+'ImportLog.log'
fLog=open(fileLog,'a')
fLog.write('\n')
date=str(datetime.datetime.now())+'\n'
fLog.write(date)
fLog.close()

cal = spireCal(jarFile='/home/user/spxcen/dp-spire/cal/spire_cal_8_phot_27Sep2011.jar')

for ob in range(nob):
	myObsid=obsids[ob]
	thisName=names[ob]
	thisMode=modes[ob]
	store_out=ProductStorage("default")
	
	print 'Searching for obsid 0x%x...'%myObsid
	
	#infile=
	
	query=MetaQuery(ObservationContext,"p","p.meta['obsid'].value==%iL"%myObsid)

	#search for local observations
	refs_local=store_out.select(query)
	nObLocal=len(refs_local)
	doObs=True
	if nObLocal > 0 and not overwrite:
		#obsid exists in local storage
		logLine='0x%x [%s: %s] Already processed (skipped)\n'%(myObsid,thisName,thisMode)
		fLog=open(fileLog,'a')
		fLog.write(logLine)
		fLog.close()
		print 'ObsID 0x%x [%s: %s] already processed. Skipping ObsID'%(myObsid,thisName,thisMode)
		doObs=False
	else:
		#get remote observations
		refs=store.select(query)
		nObRemote=len(refs)
		if nObRemote < 1:
			#no remote observation found
			logLine='0x%x [%s: %s] NOT FOUND (skipped)\n'%(myObsid,thisName,thisMode)
			fLog=open(fileLog,'a')
			fLog.write(logLine)
			fLog.close()
			print 'ObsID 0x%x [%s: %s] not found on remote server. Skipping ObsID'%(myObsid,thisName,thisMode)
			doObs=False
	if doObs:
		print 'Processing ObsID 0x%x [%s: %s]...'%(myObsid,thisName,thisMode)
		#sys.exit()
		obs=refs[0].product
		
		level1Ral = obs.level1
		
		print 'Getting scans...'
		level0_5Ral = obs.level0_5
		bbids=level0_5Ral.getBbids(0xa103)
		nlines=len(bbids)
		print "Total number of scan lines: ",nlines
		
		bsmPos		   = cal.phot.bsmPos
		lpfPar		   = cal.phot.lpfPar
		detAngOff	   = cal.phot.detAngOff
		chanTimeConst	   = cal.phot.chanTimeConst
		chanNum		   = cal.phot.chanNum
		fluxConvList	   = cal.phot.fluxConvList
		tempDriftCorrList  = cal.phot.tempDriftCorrList
		relGains		   = cal.phot.chanRelGain
		
		try:
			hpp	 = obs.auxiliary.pointing
		except:
			print 'ObsID 0x%x [%s: %s] missing Pointing Product. Skipping ObsID'%(myObsid,thisName,thisMode)
			logLine='0x%x [%s: %s] MISSING HPP (skipped)\n'%(myObsid,thisName,thisMode)
			fLog=open(fileLog,'a')
			fLog.write(logLine)
			fLog.close()
			continue
		siam = obs.auxiliary.siam
		timeCorr = obs.auxiliary.timeCorrelation
	
		level1=Level1Context(obs.meta['obsid'].value)
		level1Rg=Level1Context(obs.meta['obsid'].value)
			
		count=1
		for bbid in bbids:
			print "Starting BBID=0x%x: scan %i / %i"%(bbid,count,nlines)
			pdt  = level0_5Ral.get(bbid).pdt
			nhkt = level0_5Ral.get(bbid).nhkt
			# record the calibration tree version used by the pipeline
			pdt.calVersion = cal.version
			#
			# -----------------------------------------------------------
			# (1) join all scan legs and turnarounds together
			bbCount=bbid & 0xFFFF
			pdtLead=None
			nhktLead=None
			pdtTrail=None
			nhktTrail=None
			if bbCount >1:
				blockLead=obs.level0_5.get(0xaf000000L+bbCount-1)
				pdtLead=blockLead.pdt
				nhktLead=blockLead.nhkt
				if pdtLead != None and pdtLead.sampleTime[-1] < pdt.sampleTime[0]-3.0:
					pdtLead=None
					nhktLead=None
			if bbid < MAX(Long1d(bbids)):
				blockTrail=obs.level0_5.get(0xaf000000L+bbCount)
				pdtTrail=blockTrail.pdt
				nhktTrail=blockTrail.nhkt
				if pdtTrail != None and pdtTrail.sampleTime[0] > pdt.sampleTime[-1]+3.0:
					pdtTrail=None
					nhktTrail=None
			pdt=joinPhotDetTimelines(pdt,pdtLead,pdtTrail)
			nhkt=joinNhkTimelines(nhkt,nhktLead,nhktTrail)
			#
			# -----------------------------------------------------------
			# (2) Convert BSM timeline to angles on sky (constant for scan map)
			bat=calcBsmAngles(nhkt,bsmPos=bsmPos)
			# -----------------------------------------------------------
			# (3) Create the Spire Pointing Product for this observation
			spp=createSpirePointing(detAngOff=detAngOff,bat=bat,hpp=hpp,siam=siam)
			# -----------------------------------------------------------
			# (4) Detect jumps in the Thermistor timelines that occasionally occur,
			# leading to map artefacts introduced in the temperature drift correction
			# Also requires the Temperature Drift Correct Calibration File.
			tempDriftCorr=tempDriftCorrList.getProduct(pdt.meta["biasMode"].value,pdt.startDate)
			if pdt.meta["biasMode"].value == "nominal":
				pdt=signalJumpDetector(pdt,tempDriftCorr=tempDriftCorr, kappa=2.0,gamma=6.0,\
					gapWidth=1.0,windowWidth=40.0, filterType="DISCRETEDERIVATIVE",glitchinfo="NULL")
			# -----------------------------------------------------------
			# (5) Run the concurrent deglitcher on the timeline data
			pdt=concurrentGlitchDeglitcher(pdt,chanNum=chanNum,kappa=2.0, size = 15, correctGlitches = True)
			# -----------------------------------------------------------
			# (6) Run the wavelet deglitcher on the timeline data
			obsMode=obs.meta["obsMode"].value
			if obsMode.find("Parallel") >= 0:
				#Deglitch with Parallel Mode parameters
				pdt=waveletDeglitcher(pdt, scaleMin=1.0, scaleMax=8.0, scaleInterval=5, holderMin=-1.9,\
					holderMax=-0.3, correlationThreshold=0.69, optionReconstruction='linearAdaptive20',\
					reconstructionPointsBefore=1, reconstructionPointsAfter=3)
			else:
				#Deglitch with SPIRE-only parameters
				pdt=waveletDeglitcher(pdt, scaleMin=1.0, scaleMax=8.0, scaleInterval=7, holderMin=-3.0,\
						holderMax=-0.3, correlationThreshold=0.3, optionReconstruction='linear',\
						reconstructionPointsBefore=1, reconstructionPointsAfter=6)
			# -----------------------------------------------------------
			# (7) Apply the Low Pass Filter response correction
			pdt=lpfResponseCorrection(pdt,lpfPar=lpfPar)
			# -----------------------------------------------------------
			# (8) Apply the flux conversion 
			fluxConv=fluxConvList.getProduct(pdt.meta["biasMode"].value,pdt.startDate)
			pdt=photFluxConversion(pdt,fluxConv=fluxConv)
			# -----------------------------------------------------------
			# (9) Make the temperature drift correction
			pdt=temperatureDriftCorrection(pdt,tempDriftCorr=tempDriftCorr)
			# -----------------------------------------------------------
			# (10) Apply the bolometer time response correction
			pdt=bolometerResponseCorrection(pdt,chanTimeConst=chanTimeConst)
			# -----------------------------------------------------------
			# (11) Add pointing timelines to the data
			psp=associateSkyPosition(pdt,spp=spp)
			# -----------------------------------------------------------
			# (12) Cut the timeline back into individual scan lines .
			# NOTE: If you want include turnaround data in map making, call the following
			# task with the option "extend=True"
			psp=cutPhotDetTimelines(psp,extend=False)
			# -----------------------------------------------------------
			# (13) Apply the time correlation 
			tempPsp=timeCorrelation(psp,timeCorr)
						#tempPsp = level1Ral.getProduct(i)
			#tempPsp.setType('PSP')
			level1.addProduct(tempPsp)
			
			tempPspRg=applyRelativeGains(tempPsp, relGains)
			level1Rg.addProduct(tempPspRg)
				
			count = count+1
			
		scans=baselineRemovalMedian(level1)
			
		obs.level1=scans
		scansRg=baselineRemovalMedian(level1Rg)
			
		print 'beginning Map Maker...'
		mapPlw=naiveScanMapper(scans, array="PLW")
		mapPmw=naiveScanMapper(scans, array="PMW")
		mapPsw=naiveScanMapper(scans, array="PSW")
			
		mapPlwRg=naiveScanMapper(scansRg, array="PLW")
		mapPmwRg=naiveScanMapper(scansRg, array="PMW")
		mapPswRg=naiveScanMapper(scansRg, array="PSW")
			
		#level2=Level2Context(obs.meta.['obsid'].value)
		#level2Rg=Level2Context(obs.meta.['obsid'].value)
		#level2.setProduct("PSW", mapPsw)
		#level2.setProduct("PMW", mapPmw)
		#level2.setProduct("PLW", mapPlw)
		#level2Rg.setProduct("PSW", mapPswRg)
		#level2Rg.setProduct("PMW", mapPmwRg)
		#level2Rg.setProduct("PLW", mapPlwRg)
			
		obs.level2.setProduct("PSW", mapPsw)
		obs.level2.setProduct("PMW", mapPmw)
		obs.level2.setProduct("PLW", mapPlw)
		obs.level2.setProduct("PSWRg", mapPswRg)
		obs.level2.setProduct("PMWRg", mapPmwRg)
		obs.level2.setProduct("PLWRg", mapPlwRg)
			
		print 'Writing maps to FITS files...'
		simpleFitsWriter(mapPsw, "%smapPSW_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
		simpleFitsWriter(mapPmw, "%smapPMW_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
		simpleFitsWriter(mapPlw, "%smapPLW_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
			
		simpleFitsWriter(mapPswRg, "%smapPSW_RG_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
		simpleFitsWriter(mapPmwRg, "%smapPMW_RG_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
		simpleFitsWriter(mapPlwRg, "%smapPLW_RG_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
	
		obs.calibration.spec.refs.clear()
		obs.calibration.phot.refs.clear()
	
		store_out.save(obs)
	
		if nObLocal > 0:
			logLine='0x%x [%s: %s] Reprocessed\n'%(myObsid,thisName,thisMode)
		else:
			logLine='0x%x [%s: %s] Processed\n'%(myObsid,thisName,thisMode)
		fLog=open(fileLog,'a')
		fLog.write(logLine)
		fLog.close()
