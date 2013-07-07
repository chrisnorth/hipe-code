##regrid beam maps to (6,10,12)arsec
import sys

def reGrid():

    dirPath='/home/astrog82/spxcen/Herschel/Calibration/spire_beams_measured/'
    fileInPSW='psw_beam_1arcsec.fits'
    fileInPMW='psw_beam_1arcsec.fits'
    fileInPLW='psw_beam_1arcsec.fits'
    
    ##Read in maps
    mapInPSW=simpleFitsReader(dirPath+fileInPSW)
    mapInPMW=simpleFitsReader(dirPath+fileInPMW)
    mapInPLW=simpleFitsReader(dirPath+fileInPLW)
    
    ##Extract WCS
    wPSW=mapInPSW.wcs
    wPMW=mapInPMW.wcs
    wPLW=mapInPLW.wcs
    
    ##Regrid maps
    wPSW2=Wcs(crpix1=wPSW.crpix1,crpix2=wPSW.crpix2,crval1=wPSW.crval1,crval2=wPSW.crval2,\
    	cdelt1=wPSW.cdelt1*6.,cdelt2=wPSW.cdelt2*6.,ctype1=wPSW.ctype1,ctype2=wPSW.ctype2,\
    	equinox=wPSW.equinox,naxis2=wPSW.naxis2,naxis1=wPSW.naxis1,crota2=wPSW.crota2)
    mapOutPSW = regrid(source=mapInPSW, wcs=wPSW2)
    
    
    wPMW2=Wcs(crpix1=wPMW.crpix1,crpix2=wPMW.crpix2,crval1=wPMW.crval1,crval2=wPMW.crval2,\
    	cdelt1=wPMW.cdelt1*10.,cdelt2=wPMW.cdelt2*10.,ctype1=wPMW.ctype1,ctype2=wPMW.ctype2,\
    	equinox=wPMW.equinox,naxis2=wPMW.naxis2,naxis1=wPMW.naxis1,crota2=wPMW.crota2)
    mapOutPMW = regrid(source=mapInPMW, wcs=wPMW2)
    
    wPLW2=Wcs(crpix1=wPLW.crpix1,crpix2=wPLW.crpix2,crval1=wPLW.crval1,crval2=wPLW.crval2,\
    	cdelt1=wPLW.cdelt1*12.,cdelt2=wPLW.cdelt2*12.,ctype1=wPLW.ctype1,ctype2=wPLW.ctype2,\
    	equinox=wPLW.equinox,naxis2=wPLW.naxis2,naxis1=wPLW.naxis1,crota2=wPLW.crota2)
    mapOutPLW = regrid(source=mapInPLW, wcs=wPLW2)
    
    ##write out to files
    fileOutPSW='psw_beam_6arcsec.fits'
    fileOutPMW='pmw_beam_10arcsec.fits'
    fileOutPLW='plw_beam_12arcsec.fits'
    
    simpleFitsWriter(mapOutPSW,file=dirPath+fileOutPSW)
    simpleFitsWriter(mapOutPMW,file=dirPath+fileOutPMW)
    simpleFitsWriter(mapOutPLW,file=dirPath+fileOutPLW)
    
    #sys.exit()


dirPath='/home/astrog82/spxcen/Herschel/Calibration/spire_beams_measured/'
fileIn1=['psw_beam_1arcsec.fits',\
	'pmw_beam_1arcsec.fits',\
	'plw_beam_1arcsec.fits']

fileIn2=['psw_beam_6arcsec.fits',\
	'pmw_beam_10arcsec.fits',\
	'plw_beam_12arcsec.fits']

radiusArcsec=[[22.,30.,45.,50.,60.,100.,500.],\
	[30.,45.,55.,65.,80.,100.,500.],\
	[42.,50.,60.,90.,100.,120.,500.]]

#Max source radii that are valid
radMax=[[5,6,6,6,7],[3,6,6,6,7],[3,4,6,6,7]]

innerArcsec=[60.,100.,200.,300.,500]
outerArcsec=[90.,150.,250.,350.,550]

#Min BR radii that are valid
bgMin=[[0,0,0,0,0,1,4],[0,0,0,1,1,1,4],[0,0,0,1,1,2,4]]

nRad=len(radiusArcsec[0])
nBg=len(innerArcsec)

fileDat=['beam_regridding_PSW.dat',\
	'beam_regridding_PMW.dat',\
	'beam_regridding_PLW.dat']

for b in range(3):
	fDat=open(fileDat[b],'w')
	hdr='#bgInner, bgOuter, BeamArea (small), Error, BeamArea (large), Error'
	fDat.write(hdr)
	mapIn1=simpleFitsReader(dirPath+fileIn1[b])
	mapIn2=simpleFitsReader(dirPath+fileIn2[b])
	for r in range(nRad):
		radSrc=radiusArcsec[b][r]
		bgInner=innerArcsec[bgMin[b][r]]
		bgOuter=outerArcsec[bgMin[b][r]]

		result1 = annularSkyAperturePhotometry(image=mapIn1,\
			fractional=1, centerRA='0', centerDec='0', \
			radiusArcsec=radSrc, innerArcsec=bgInner, outerArcsec=bgOuter)
		result2 = annularSkyAperturePhotometry(image=mapIn2,\
			fractional=1, centerRA='0', centerDec='0', \
			radiusArcsec=radSrc, innerArcsec=bgInner, outerArcsec=bgOuter)

		#sys.exit()
		beamArea1=result1["Results table"][0].data[0]
		beamErr1=result1["Results table"][3].data[0]
		beamArea2=result2["Results table"][0].data[0]
		beamErr2=result2["Results table"][3].data[0]

		line='%.1f , %.1f , %.3f , %.5f , %.5f , %.5f'%\
			(0.,radSrc,beamArea1,beamErr1,beamArea2,beamErr2)
		fDat.write(line)
		print line

	for bg in range(nBg):
		radSrc=radiusArcsec[b][0]
		bgInner=innerArcsec[bg]
		bgOuter=outerArcsec[bg]

		result1 = annularSkyAperturePhotometry(image=mapIn1,\
			fractional=1, centerRA='0', centerDec='0', \
			radiusArcsec=radSrc, innerArcsec=bgInner, outerArcsec=bgOuter)
		result2 = annularSkyAperturePhotometry(image=mapIn2,\
			fractional=1, centerRA='0', centerDec='0', \
			radiusArcsec=radSrc, innerArcsec=bgInner, outerArcsec=bgOuter)
		beamArea1=result1["Results table"][0].data[1]
		beamErr1=result1["Results table"][3].data[1]
		beamArea2=result2["Results table"][0].data[1]
		beamErr2=result2["Results table"][3].data[1]

		line='%.1f , %.1f , %.5f , %.5f , %.3f , %.5f'%\
			(bgInner,bgOuter,beamArea1,beamErr1,beamArea2,beamErr2)
		fDat.write(line)
		print line

	fDat.close()