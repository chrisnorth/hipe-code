## Construct tables (FITS files) of beam areas for aperture photometry

import sys

###plot aperture photometry on sources
#sys.exit()

def makeApCorr(dirPath,profFile,fileOut,band):

	########################################################
	##Make table files
	########################################################
	from herschel.ia.toolbox.fit import FitterFunction
	from herschel.ia.numeric.toolbox import RealFunction
	from herschel.ia.numeric.toolbox.integr import TrapezoidalIntegrator
	from herschel.ia.numeric import String1d,Double1d,Double2d
	import math
	import sys.exit
	from herschel.ia.numeric.toolbox.interp import LinearInterpolator
	from herschel.ia.io.ascii import AsciiTableTool
	from herschel.ia.dataset import TableDataset,Column,StringParameter
	from herschel.ia.toolbox.image import Wcs,SimpleImage
	from herschel.ia.toolbox.util import simpleFitsWriter
	########################################################
	##Read in azimuthal beam profiles (made from beam maps)
	########################################################

	print dirPath+profFile
	f=open(dirPath+profFile,'r')
	lines=String1d(f.readlines())
	nl=lines.size
	
	rad=Double1d(nl)
	profPSW=Double1d(nl)
	profPMW=Double1d(nl)
	profPLW=Double1d(nl)
	
	minProf=1.e-8
	
	for l in range (nl):
		line=lines[l].split(',')
		rad[l]=float(line[0])
		profPSW[l]=max(float(line[1]),minProf)
		profPMW[l]=max(float(line[2]),minProf)
		profPLW[l]=max(float(line[3]),minProf)

	if band=='PSW':
		profBeam=profPSW
		fileDat='beam_frac_PSW.dat'
	elif band=='PMW':
		profBeam=profPMW
		fileDat='beam_frac_PMW.dat'
	elif band=='PLW':
		profBeam=profPLW
		fileDat='beam_frac_PLW.dat'
	else:
		print '"%s" is unknown band. Must be [PSW|PMW|PLW]'%(band)
		sys.exit()
	
	#####################################
	## Calculate beam areas
	#####################################
	
	## Calculate Integrand 2*pi*rad*B(rad)
	intArgBeam=profBeam * 2. * math.pi * rad
	
	## Make interpolation of integrands
	intpBeam=LinearInterpolator(rad,intArgBeam)
	
	areaBeam=Double1d(nl)

	#areaBeam2=Double1d(nl)
	print 'Calculating beam areas...'
	
	## Integration using interpolation takes ~10-20s
	## Simple integrator much faster
	for r in range(nl):
		if r == 0:
			areaBeam[r]=0.
		else:
			areaBeam[r]=areaBeam[r-1] + profBeam[r]*2.*math.pi*rad[r]
	
		##BUILT IN INTEGRATOR MISBEHAVES NEAR zero-values!!!
		#intLim=TrapezoidalIntegrator(0.,rad[r])
		#areaBeam2[r]=intLim.integrate(intpBeam)
	#	print rad_r,areaBeam[r]
		
	maxAreaBeam=max(areaBeam)
	normBeam=areaBeam/maxAreaBeam

	## Make array
	bgArrBeam=Double2d(nl,nl)
	
	for r1 in range(nl):
		bgArrBeam[r1,r1:]=normBeam[r1:]-normBeam[r1]
		if r1 > 0:
			bgArrBeam[r1,0:r1-1]=normBeam[r1]-normBeam[0:r1-1]
	
	wcs=Wcs(crval1=0.,crval2=0.,crpix1=0.,crpix2=0.,cdelt1=1.,cdelt2=1.,\
		cunit1='arcsec',cunit2='arcsec')
	desc='%s background subtraction'%(band)
	bgImgBeam=SimpleImage(description=desc,\
		image=bgArrBeam,wcs=wcs,instrument='SPIRE')
	bgImgBeam.meta["USAGE1"]=StringParameter("This is the %s integrated beam between radius1 and radius2"%(band))
	bgImgBeam.meta["USAGE2"]=StringParameter("The beam between r1 and r2 is DATA[r1,r2] (order not important)")
	bgImgBeam.meta["wavelength"]=StringParameter(band)
	simpleFitsWriter(bgImgBeam,file=dirPath+fileOut)

	##Write to data file
	fDat=open(dirPath+fileDat,'w')
	line='#Fraction of beam (normalised area) at a given radius for %s band.\n'%(band)
	fDat.write(line)
	line='#Beam is assumed to be >=0, so area monotonically increases\n'
	fDat.write(line)
	line='#Beam area at rad=%.1f arcsec: %.6f arcsec^2\n'%(rad[nl-1],areaBeam[nl-1])
	fDat.write(line)
	line='#\n'
	fDat.write(line)
	line='#radius [arcsec], Fraction (%s)\n'%(band)
	fDat.write(line)
	
	for l in range(nl):
		line='%.1f, %.6f\n'%(rad[l],normBeam[l])
		fDat.write(line)
	
	print 'Written to %s'%(dirPath+fileDat)
	fDat.close()
	
def getApCorr(radSrc,radIn,radOut,band="PSW"):
	from os.path import isfile
	from herschel.ia.toolbox.util import simpleFitsReader
	import sys.exit
	import signal
	
	###############################################
	##Set path and filenames
	###############################################
	
	dirPath='/home/astrog82/spxcen/Herschel/Calibration/RelGains/'
	profFile='beamprof_rad.dat'
	filePSW='bg_im_PSW.fits'
	filePMW='bg_im_PMW.fits'
	filePLW='bg_im_PLW.fits'

	if band=='PSW':
		fileIn=filePSW
	elif band=='PMW':
		fileIn=filePSW
	elif band=='PLW':
		fileIn=filePLW
	else:
		print '"%s" is unknown band. Must be [PSW|PMW|PLW]'%(band)
		sys.exit()
	##Check if file exists
	
	if not isfile(dirPath + fileIn):
		print 'File does not exist...'
		## Make new FITS files
		
		makeApCorr(dirPath,profFile,fileIn,band)
		print 'Now it does...'
	#else:
		#print 'File exists...'
	
	apCorrIm=simpleFitsReader(dirPath+fileIn)
	
	srcCorr=apCorrIm["image"].data[0,int(radSrc)]
	print 'srcCorr=',srcCorr
	bkgCorr=apCorrIm["image"].data[int(radIn),int(radOut)]
	#print 'bkgCorr=',bkgCorr
	bkgCorr2=apCorrIm["image"].data[0,int(radOut)]-apCorrIm["image"].data[0,int(radIn)]
	print 'bkgCor',bkgCorr,bkgCorr2
	apCorr=srcCorr - bkgCorr
	print 'apCorr=',apCorr
	
	return(1./srcCorr,1./bkgCorr,1./apCorr)


#if __name__=='__main__':
#	maincode()
