
from herschel.ia.pal.pool.http import HttpClientPool

dirPath = '/data/Herschel/Calibration/Inputs/'

pool = HttpClientPool ("http://wakefield.bnsc.rl.ac.uk/hcss/pal", "ops")

#pool.setAuthentication(Configuration.getProperty("hcss.access.authentication.authstring"))
store = ProductStorage (pool)    # Just the Data
store.authenticate()


store_out=ProductStorage("Mars")
cal = spireCal(jarFile='/home/chris/dp-spire/cal/spire_cal_8_phot_27Sep2011.jar')

print 'Accessing data...'

myObsid	=0x5000532c
query=MetaQuery(ObservationContext,"p","p.meta['obsid'].value==%iL"%myObsid)

refs=store.select(query)
obs=refs[0].product
relGains= cal.phot.chanRelGain

level1Ral = obs.level1

print 'Getting scans...'

level1=Level1Context(obs.meta['obsid'].value)
#level1Rg=Level1Context(obs.meta['obsid'].value)
for i in range (level1Ral.getCount()):
	tempPsp = level1Ral.getProduct(i)
	tempPsp.setType('PSP')
#	tempPspRg = level1Ral.getProduct(i)
#	tempPspRg.setType('PSP')
	level1.addProduct(tempPsp)
#	applyRelativeGains(tempPspRg, relGains)
#	level1Rg.addProduct(tempPspRg)

scans=baselineRemovalMedian(level1)
#scansRg=baselineRemovalMedian(level1Rg)

print 'beginning Map Maker...'
mapPlw=naiveScanMapper(scans, array="PLW",resolution=1.)
mapPmw=naiveScanMapper(scans, array="PMW",resolution=1.)
mapPsw=naiveScanMapper(scans, array="PSW",resolution=1.)

#mapPlwRg=naiveScanMapper(scansRg, array="PLW")
#mapPmwRg=naiveScanMapper(scansRg, array="PMW")
#mapPswRg=naiveScanMapper(scansRg, array="PSW")

#level2=Level2Context(obs.meta.['obsid'].value)
#level2Rg=Level2Context(obs.meta.['obsid'].value)
#level2.setProduct("PSW", mapPsw)
#level2.setProduct("PMW", mapPmw)
#level2.setProduct("PLW", mapPlw)
#level2Rg.setProduct("PSW", mapPswRg)
#level2Rg.setProduct("PMW", mapPmwRg)
#level2Rg.setProduct("PLW", mapPlwRg)

print 'Writing maps to FITS files...'
simpleFitsWriter(mapPsw, "%sMars_mapPSW_%s_1arcsec.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
simpleFitsWriter(mapPmw, "%sMars_mapPMW_%s_1arcsec.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
simpleFitsWriter(mapPlw, "%sMars_mapPLW_%s_1arcsec.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))

#simpleFitsWriter(mapPswRg, "%smapPSW_RG_%s_1arcsec.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
#simpleFitsWriter(mapPmwRg, "%smapPMW_RG_%s_1arcsec.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
#simpleFitsWriter(mapPlwRg, "%smapPLW_RG_%s_1arcsec.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
