###Runs George's PSF fitter

obsid=1342213459

myObsid	=  0x50008d53
myDataPool = "default"
outDir	 = "/data/Herschel/Calibration/RelGains"

obs=getObservation(myObsid,poolName=myDataPool)
bbids=obs.level0_5.getBbids(0xa103)
nscan=len(bbids)

dirPath='/data/Herschel/Calibration/RelGains/Obs/'
count=0
for i in range(nscan):
	count+=1
	filename = 'psp_'+hex(obs.meta['obsid'].value)[2:10]+"_scan"+str(counter)+'_BLsub.fits'
	FitsArchive.load(dirPath+filename,scan)
	asdfsdaf
