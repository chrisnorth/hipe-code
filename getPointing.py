##Get asteroid locations
from herschel.ia.pal.pool.http import HttpClientPool
from herschel.share.fltdyn.ephem.horizons import BasicHorizons
from herschel.share.fltdyn.math import Direction
from herschel.ia.obs.auxiliary.fltdyn import Ephemerides
from herschel.ia.obs.auxiliary.fltdyn import Horizons
import sys

##Set up pool
pool = HttpClientPool ("http://wakefield.bnsc.rl.ac.uk/hcss/pal", "ops")
#pool.setAuthentication(Configuration.getProperty("hcss.access.authentication.authstring"))
store = ProductStorage (pool)	# Just the Data
store.authenticate()


#files of objects without pointing info
dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/Obs/'
filesIn=[\
	'Neptune_obs_LargeScan_noPt.dat', 'Neptune_obs_SmallScan_noPt.dat',\
	'Uranus_obs_LargeScan_noPt.dat', 'Uranus_obs_SmallScan_noPt.dat',\
	'1_Ceres_obs_LargeScan_noPt.dat', '1_Ceres_obs_SmallScan_noPt.dat',\
	'2_Pallas_obs_LargeScan_noPt.dat', '2_Pallas_obs_SmallScan_noPt.dat',\
	'3_Juno_obs_LargeScan_noPt.dat', '3_Juno_obs_SmallScan_noPt.dat',\
	'4_Vesta_obs_LargeScan_noPt.dat', '4_Vesta_obs_SmallScan_noPt.dat',\
	'10_Hygeia_obs_LargeScan_noPt.dat', '10_Hygeia_obs_SmallScan_noPt.dat'\
	]

for fileIn in filesIn:
	fileOut=fileIn[:-9]+'_withPt.dat'

	#Opne files
	fIn=open(dirPath+fileIn,'r')
	fOut=open(dirPath+fileOut,'w')
	fOut.close()
	linesIn=String1d(fIn.readlines())
	nLine=len(linesIn)
	#tarName='UNKNOWN'

	for l in range(nLine):
		if linesIn[l].find('#') >=0:
			#HEADER line, copy it directly to outfile
			fOut=open(dirPath+fileOut,'a')
			fOut.write(linesIn[l])
			fOut.close()
			if linesIn[l].find('Name') > 0:
				tarName=linesIn[l].split(':')[1][1:-1]
		else:
			myObsid=int(linesIn[l].split(',')[0],16)
			print linesIn[l]
			print 'Searching for obsid 0x%x...'%myObsid

			query=MetaQuery(ObservationContext,"p","p.meta['obsid'].value==%iL"%myObsid)

			#get remote observations
			refs=store.select(query)
			nObRemote=len(refs)

			if nObRemote < 1:
				#no remote observation found
				print 'ObsID 0x%x [%s: %s] not found on remote server. Skipping ObsID'%(myObsid,thisName,thisMode)
				continue
			
			#read in obsid
			print 'Reading obsid 0x%x...'%myObsid
			obs=refs[0].product

			level0=obs.level0
			nScan=level0.getCount()
			ppt0=level0.refs[0].product
			ppt1=level0.refs[nScan-1].product
			try: 
				ephemProduct = obs.calibration.orbit
				horizonsProduct = obs.calibration.horizons
				#horizonsProduct  = FitsArchive().load('/home/alp/Spire_work/SPIRE_3615/hauxauxhorizons.fits')
				#ephemProduct = FitsArchive().load('/home/alp/Spire_work/SPIRE_3615/hauxauxorbitp.fits')
				naifid = ppt0.meta["naifId"].value

				ephemObs    = Ephemerides(ephemProduct)
				horizonsObs = Horizons(horizonsProduct, ephemObs)
	
				###GEOMETRIC state (No correction)
				#usePointingCorrection = BasicHorizons.Correction.NONE 
	
				###APPARENT state (Light Time and Stellar aberration correction)
				### This is the Ephemeride position  that Tanya used and reported in SPIRE-3615
				#usePointingCorrection = BasicHorizons.Correction.LTS 
	
				## ASTROMETRIC state (Light Time correction only)
				## It has been agreed that all 3 Herschel instruments will apply this correction to SSO poinintg. 
				## Spire Map positions of SSOs should be compared to Ephemeride positions obtained in the
				## way shown below
				usePointingCorrection = BasicHorizons.Correction.LT
	
				dir1 = Direction(horizonsObs.spacecraftVectorTo(int(naifid), ppt0.startDate, usePointingCorrection))
				dir2 = Direction(horizonsObs.spacecraftVectorTo(int(naifid), ppt1.endDate, usePointingCorrection))
	
				posEphemStart = [dir1.raDegrees,dir1.decDegrees]
				posEphemEnd   = [dir2.raDegrees,dir2.decDegrees]
	
				print '0x%x (%s) Start (%s): [%f , %f]'%(obsid,tarName,ppt0.startDate,posEphemStart[0],posEphemStart[1])
				print '0x%x (%s) End (%s):   [%f , %f]'%(obsid,tarName,ppt1.endDate,posEphemEnd[0],posEphemEnd[1])
				fOut=open(dirPath+fileOut,'a')
				line='0x%x , %f , %f\n'%(myObsid,posEphemStart[0],posEphemStart[1])
				fOut.write(line)
				fOut.close()
			except:
				print 'No horizons/ephem file available for 0x%x (%s)'%(myObsid,tarName)
				line='0x%x , %i , %i, #no pointing available\n'%(myObsid,-1, -1)
				fOut=open(dirPath+fileOut,'a')
				fOut.write(line)
				fOut.close()
	#fOut.close()
	fIn.close()