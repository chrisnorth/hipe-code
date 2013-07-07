import sys
sys.path.append("/data/Herschel/Calibration/RelGains")
#from herschel.ia.toolbox.image import bendoSourceFit_v0.9
from herschel.ia.pal.pool.http import HttpClientPool

dirPath = '/data/Herschel/Calibration/RelGains/Obs/'

pool = HttpClientPool ("http://wakefield.bnsc.rl.ac.uk/hcss/pal", "ops")

#pool.setAuthentication(Configuration.getProperty("hcss.access.authentication.authstring"))
store = ProductStorage (pool)    # Just the Data
store.authenticate()

myObsid=0x50005984
#store_out=ProductStorage("Neptune")
cal = spireCal(jarFile='/home/user/spxcen/dp-spire/cal/spire_cal_8_phot_27Sep2011.jar')

print 'Accessing data...'

#infile=

query=MetaQuery(ObservationContext,"p","p.meta['obsid'].value==%iL"%myObsid)
refs=store.select(query)
obs=refs[0].product
relGains= cal.phot.chanRelGain

level1Ral = obs.level1

scans=baselineRemovalMedian(level1Ral)

scansInput = []
for i in range(scans.count):
	print 'scan %d/%d'%(i+1,scans.count)
	tempPSP=scans.getProduct(i)
	scansInput.append(tempPSP)
sys.exit()
scansContext = ScanContext(scansInput)
raSrc=269.1508446821
decSrc=51.4885644774
radiusArcsec=22.
innerArcsec=250.
outerArcsec=300.

fitter=TimelineSourceFitterTask()
output=fitter(input=scansContext,array='PSW',\
	sourcePositionEstimate=Double1d([raSrc,decSrc]),\
	rPeak=radiusArcsec,rBackground=Double1d([innerArcsec,outerArcsec]),\
	fitEllipticalGaussian=True,\
	useBackInFit=True,allowVaryBackground=True)



fitter2=bendoSourceFit(scans)
fitter2.setModelGauss2DRot()
fitter2.setUseBackInFit(Boolean.TRUE)
fitter2.setVaryBackInFit(Boolean.TRUE)
fitter2.setPlotFitCenter()
param=fitter2.fit("PSW",raSrc,decSrc,radiusArcsec,innerArcsec,outerArcsec)
paramErr=fitter2.getStandardDeviation()

print param[0]-output.meta["fitFlux"].getValue()
print param[1]-output.meta["fitRa"].getValue()
print param[2]-output.meta["fitDec"].getValue()
print param[3]-output.meta["fitSigma1"].getValue()
print param[6]-output.meta["fitBackgroundParm1"].getValue()
