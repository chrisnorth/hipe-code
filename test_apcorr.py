import sys
sys.path.append("/data/Herschel/Calibration/RelGains")
from herschel.ia.toolbox.image import sourceFitting
import apcorrfn

#def maincode():

obsid=0x5000AE55

dirPath='/data/Herschel/Calibration/RelGains/'
tlFile='timelinefit.csv'

f=open(dirPath+tlFile,'r')
lines=String1d(f.readlines())
nl=lines.size

for l in range(nl):
	if not lines[l][0] == '#':
		line=lines[l].split(',')
		obsIDx=int(line[0],16)
		if obsIDx == obsid:
			obsName=line[1]
			centerRA=float(line[2])
			centerDec=float(line[3])
			tlFluxPSW=float(line[4])
			tlFluxPMW=float(line[5])
			tlFluxPLW=float(line[6])
			tlErrPSW=float(line[7])
			tlErrPMW=float(line[8])
			tlErrPLW=float(line[9])

###Use timeline fitter

print 'Using obsID %X [%s] (RA,Dec)=(%.3f,%.3f)'%(obsIDx,obsName,centerRA,centerDec)

#centerRA='269.15084468207004'
#centerDec='51.48856447735525'

###read in maps from files
filePSWrg='Obs/mapPSW_rg_'+str(hex(obsid))[2:10]+'.fits'
filePMWrg='Obs/mapPMW_rg_'+str(hex(obsid))[2:10]+'.fits'
filePLWrg='Obs/mapPLW_rg_'+str(hex(obsid))[2:10]+'.fits'

filePSW='Obs/mapPSW_'+str(hex(obsid))[2:10]+'.fits'
filePMW='Obs/mapPMW_'+str(hex(obsid))[2:10]+'.fits'
filePLW='Obs/mapPLW_'+str(hex(obsid))[2:10]+'.fits'

##read in files (Jy/beam)
mapPSWin = fitsReader(file = dirPath+filePSW)
mapPMWin = fitsReader(file = dirPath+filePMW)
mapPLWin = fitsReader(file = dirPath+filePLW)
mapPSWinrg = fitsReader(file = dirPath+filePSWrg)
mapPMWinrg = fitsReader(file = dirPath+filePMWrg)
mapPLWinrg = fitsReader(file = dirPath+filePLWrg)

print mapPSWin.getUnit()

###Caluclate beam areas
#Beam areas (arcsec)
beamAreasArcsec=[433.1,776.9,1628.3]
##Pixel sizes (arcsec)
pixelArcsec=[(mapPSWin.wcs.cdelt1*3600.)**2,(mapPMWin.wcs.cdelt1*3600.)**2,(mapPLWin.wcs.cdelt1*3600.)**2]
##Beam areas (pixels)
beamAreasPixel=[beamAreasArcsec[0]/pixelArcsec[0],beamAreasArcsec[1]/pixelArcsec[1],beamAreasArcsec[2]/pixelArcsec[2]]

##Convert to Jy/arcsec^2
mapPSWArcsec = convertImageUnit(image=mapPSWin,newUnit='Jy/arcsec^2',beamArea=beamAreasPixel[0])
mapPMWArcsec = convertImageUnit(image=mapPMWin,newUnit='Jy/arcsec^2',beamArea=beamAreasPixel[1])
mapPLWArcsec = convertImageUnit(image=mapPLWin,newUnit='Jy/arcsec^2',beamArea=beamAreasPixel[2])
mapPSWArcsecrg = convertImageUnit(image=mapPSWinrg,newUnit='Jy/arcsec^2',beamArea=beamAreasPixel[0])
mapPMWArcsecrg = convertImageUnit(image=mapPMWinrg,newUnit='Jy/arcsec^2',beamArea=beamAreasPixel[1])
mapPLWArcsecrg = convertImageUnit(image=mapPLWinrg,newUnit='Jy/arcsec^2',beamArea=beamAreasPixel[2])

print mapPSWArcsec.getUnit()
print beamAreasArcsec
print pixelArcsec
print beamAreasPixel

##Set ApPhot radii
radiusArcsecPSW=22.0
innerArcsecPSW=60.
outerArcsecPSW=90.

radiusArcsecPMW=30.0
innerArcsecPMW=60.
outerArcsecPMW=90.

radiusArcsecPLW=42.0
innerArcsecPLW=60.
outerArcsecPLW=90.

##do Aperture photometry
resultPSW = annularSkyAperturePhotometry(image=mapPSWArcsec, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPSW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPSW)

resultPMW = annularSkyAperturePhotometry(image=mapPMWArcsec, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPMW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPMW)

resultPLW = annularSkyAperturePhotometry(image=mapPLWArcsec, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPLW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPLW)

resultPSWrg = annularSkyAperturePhotometry(image=mapPSWArcsecrg, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPSW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPSW)

resultPMWrg = annularSkyAperturePhotometry(image=mapPMWArcsecrg, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPMW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPMW)

resultPLWrg = annularSkyAperturePhotometry(image=mapPLWArcsecrg, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPLW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPLW)

##Get pixel locations of source locations
(yPSW,xPSW)=mapPSWin.wcs.getPixelCoordinates(centerRA,centerDec)
(yPMW,xPMW)=mapPMWin.wcs.getPixelCoordinates(centerRA,centerDec)
(yPLW,xPLW)=mapPLWin.wcs.getPixelCoordinates(centerRA,centerDec)

##Extract source flux (in Jy)
mapFluxPSW=resultPSW["Results table"][0].data[2]#*pixelArcsec[0]
mapFluxPMW=resultPMW["Results table"][0].data[2]#*pixelArcsec[0]
mapFluxPLW=resultPLW["Results table"][0].data[2]#*pixelArcsec[0]
mapFluxPSWrg=resultPSWrg["Results table"][0].data[2]#*pixelArcsec[0]
mapFluxPMWrg=resultPMWrg["Results table"][0].data[2]#*pixelArcsec[0]
mapFluxPLWrg=resultPLWrg["Results table"][0].data[2]#*pixelArcsec[0]

##Extract source flux error
mapErrPSW=resultPSW["Results table"][2].data[2]
mapErrPMW=resultPMW["Results table"][2].data[2]
mapErrPLW=resultPLW["Results table"][2].data[2]
mapErrPSWrg=resultPSWrg["Results table"][2].data[2]
mapErrPMWrg=resultPMWrg["Results table"][2].data[2]
mapErrPLWrg=resultPLWrg["Results table"][2].data[2]

#print results
print 'Map value (PSW): %.4f +/- %.4f'%(mapFluxPSW,mapErrPSW)
print 'Map value (PMW): %.4f +/- %.4f'%(mapFluxPMW,mapErrPMW)
print 'Map value (PLW): %.4f +/- %.4f'%(mapFluxPLW,mapErrPLW)
print '--'

##Get aperturecorrection parameters
(srcCorrPSW,bgCorrPSW,apCorrPSW)=apcorrfn.getApCorr(radiusArcsecPSW,innerArcsecPSW,outerArcsecPSW,band='PSW')
(srcCorrPMW,bgCorrPMW,apCorrPMW)=apcorrfn.getApCorr(radiusArcsecPMW,innerArcsecPMW,outerArcsecPMW,band='PMW')
(srcCorrPLW,bgCorrPLW,apCorrPLW)=apcorrfn.getApCorr(radiusArcsecPLW,innerArcsecPLW,outerArcsecPLW,band='PLW')

##Print numbers
print '--'
print 'Source Correction (PSW):  %.3f'%(1./srcCorrPSW)
print 'Source Correction (PMW):  %.3f'%(1./srcCorrPMW)
print 'Source Correction (PLW):  %.3f'%(1./srcCorrPLW)
print '--'
print 'Background Correction (PSW):  %.3f'%(1./bgCorrPSW)
print 'Background Correction (PMW):  %.3f'%(1./bgCorrPMW)
print 'Background Correction (PLW):  %.3f'%(1./bgCorrPLW)
print '--'
print 'Correction (PSW):  %.3f'%(apCorrPSW)
print 'Correction (PMW):  %.3f'%(apCorrPMW)
print 'Correction (PLW):  %.3f'%(apCorrPLW)

mapCorrFluxPSW=mapFluxPSW * apCorrPSW
mapCorrFluxPMW=mapFluxPMW * apCorrPMW
mapCorrFluxPLW=mapFluxPLW * apCorrPLW
mapCorrFluxPSWrg=mapFluxPSWrg * apCorrPSW
mapCorrFluxPMWrg=mapFluxPMWrg * apCorrPMW
mapCorrFluxPLWrg=mapFluxPLWrg * apCorrPLW

mapCorrErrPSW=mapErrPSW * apCorrPSW
mapCorrErrPMW=mapErrPMW * apCorrPMW
mapCorrErrPLW=mapErrPLW * apCorrPLW
mapCorrErrPSWrg=mapErrPSWrg * apCorrPSW
mapCorrErrPMWrg=mapErrPMWrg * apCorrPMW
mapCorrErrPLWrg=mapErrPLWrg * apCorrPLW

badmapFluxPSW=mapFluxPSW * srcCorrPSW
badmapFluxPMW=mapFluxPMW * srcCorrPMW
badmapFluxPLW=mapFluxPLW * srcCorrPLW
badmapFluxPSWrg=mapFluxPSWrg * srcCorrPSW
badmapFluxPMWrg=mapFluxPMWrg * srcCorrPMW
badmapFluxPLWrg=mapFluxPLWrg * srcCorrPLW

badmapErrPSW=mapErrPSW * srcCorrPSW
badmapErrPMW=mapErrPMW * srcCorrPMW
badmapErrPLW=mapErrPLW * srcCorrPLW
badmapErrPSWrg=mapErrPSWrg * srcCorrPSW
badmapErrPMWrg=mapErrPMWrg * srcCorrPMW
badmapErrPLWrg=mapErrPLWrg * srcCorrPLW

print '--'
print 'Map Bad Corr. (PSW):   %.4f +/- %.4f'%(badmapFluxPSW,badmapErrPSW)
print 'Map Bad Corr. (PMW):   %.4f +/- %.4f'%(badmapFluxPMW,badmapErrPMW)
print 'Map Bad Corr. (PLW):   %.4f +/- %.4f'%(badmapFluxPLW,badmapErrPLW)
print '--'
print 'Map Corrected (PSW):   %.4f +/- %.4f'%(mapCorrFluxPSW,mapCorrErrPSW)
print 'Map Corrected (PMW):   %.4f +/- %.4f'%(mapCorrFluxPMW,mapCorrErrPMW)
print 'Map Corrected (PLW):   %.4f +/- %.4f'%(mapCorrFluxPLW,mapCorrErrPLW)
print '--'
print 'Timeline (PSW):    %.4f +/- %.4f'%(tlFluxPSW,tlErrPSW)
print 'Timeline (PMW):    %.4f +/- %.4f'%(tlFluxPMW,tlErrPMW)
print 'Timeline (PLW):    %.4f +/- %.4f'%(tlFluxPLW,tlErrPLW)

#fileOut=obsName+'_'+hex(obsid)+'.dat'
fileOut='%s_0x%X_rg_%.1f_%.1f_%.1f_%.1f-%.1f.dat'%\
	(obsName,obsid,radiusArcsecPSW,radiusArcsecPMW,radiusArcsecPLW,\
	innerArcsecPSW,outerArcsecPSW)
fileOut='%s_0x%X_%.1f_%.1f_%.1f_%.1f-%.1f.dat'%\
	(obsName,obsid,radiusArcsecPSW,radiusArcsecPMW,radiusArcsecPLW,\
	innerArcsecPSW,outerArcsecPSW)
print 'Writing to '+fileOut+'...'
fOut=open(dirPath+fileOut,'w')
line='#ObsId: 0x%X\n'%obsid
fOut.write(line)
line='#Name: %s\n'%obsName
fOut.write(line)
line='#RA/Dec [deg]: %.7f,%.7f\n'%(centerRA,centerDec)
fOut.write(line)
line='#Radius [arcsec]: %.1f, %.1f, %.1f\n'%(radiusArcsecPSW,radiusArcsecPMW,radiusArcsecPLW)
fOut.write(line)
line='#Bg inner [arcsec]: %.1f, %.1f, %.1f\n'%(innerArcsecPSW,innerArcsecPMW,innerArcsecPLW)
fOut.write(line)
line='#Bg outer [arcsec]: %.1f, %.1f, %.1f\n'%(outerArcsecPSW,outerArcsecPMW,outerArcsecPLW)
fOut.write(line)
line='\n'
fOut.write(line)
line='#Parameter: PSW, ErrPSW, PMW, ErrPMW, PLW, ErrPLW\n'
fOut.write(line)

line='Timeline Flux: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(tlFluxPSW, tlErrPSW, tlFluxPMW, tlErrPMW, tlFluxPLW, tlErrPLW)
fOut.write(line)

line='Map Flux: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(mapFluxPSW, mapErrPSW, mapFluxPMW, mapErrPMW, mapFluxPLW, mapErrPLW)
fOut.write(line)

line='Map Flux RG: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(mapFluxPSWrg, mapErrPSWrg, mapFluxPMWrg, mapErrPMWrg, mapFluxPLWrg, mapErrPLWrg)
fOut.write(line)

line='SrcCorr Fact: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(srcCorrPSW, 0., srcCorrPMW, 0., srcCorrPLW, 0.)
fOut.write(line)

line='SrcCorr Flux: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(badmapFluxPSW, badmapErrPSW, badmapFluxPMW, badmapErrPMW, badmapFluxPLW, badmapErrPLW)
fOut.write(line)

line='SrcCorr Flux RG: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(badmapFluxPSWrg, badmapErrPSWrg, badmapFluxPMWrg, badmapErrPMWrg, badmapFluxPLWrg, badmapErrPLWrg)
fOut.write(line)

line='ApCorr Fact: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(apCorrPSW, 0., apCorrPMW, 0., apCorrPLW, 0.)
fOut.write(line)

line='ApCorr Flux: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(mapCorrFluxPSW, mapCorrErrPSW, mapCorrFluxPMW, mapCorrErrPMW, mapCorrFluxPLW, mapCorrErrPLW)
fOut.write(line)

line='ApCorr Flux RG: %.7f, %.7f, %.7f, %.7f, %.7f,%.7f\n'%\
	(mapCorrFluxPSWrg, mapCorrErrPSWrg, mapCorrFluxPMWrg, mapCorrErrPMWrg, mapCorrFluxPLWrg, mapCorrErrPLWrg)
fOut.write(line)

fOut.close()

#if __name__=='__main__':
#	maincode()
