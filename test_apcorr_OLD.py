import sys
sys.path.append("/data/Herschel/Calibration/RelGains")
from herschel.ia.toolbox.image import sourceFitting
import apcorrfn

#def maincode():

obsid=0x50005984

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
			tlSbPSW=float(line[4])
			tlSbPMW=float(line[5])
			tlSbPLW=float(line[6])
			tlErrPSW=float(line[7])
			tlErrPMW=float(line[8])
			tlErrPLW=float(line[9])

print 'Using obsID %X [%s] (RA,Dec)=(%.3f,%.3f)'%(obsIDx,obsName,centerRA,centerDec)

#centerRA='269.15084468207004'
#centerDec='51.48856447735525'

###read in maps from files
filePSW='Obs/mapPSW_'+str(hex(obsid))[2:10]+'.fits'
filePMW='Obs/mapPMW_'+str(hex(obsid))[2:10]+'.fits'
filePLW='Obs/mapPLW_'+str(hex(obsid))[2:10]+'.fits'

mapPSWin = fitsReader(file = dirPath+filePSW)
mapPMWin = fitsReader(file = dirPath+filePMW)
mapPLWin = fitsReader(file = dirPath+filePLW)

beamAreas=[433.1,776.9,1628.3]

##convert to Jy/arcsec^2
mapPSW = convertImageUnit(image=mapPSWin,newUnit='Jy/pixel',beamArea=beamAreas[0])
mapPLW = convertImageUnit(image=mapPMWin,newUnit='Jy/pixel',beamArea=beamAreas[1])
mapPLW = convertImageUnit(image=mapPLWin,newUnit='Jy/pixel',beamArea=beamAreas[2])

radiusArcsecPSW=22.0
innerArcsecPSW=100.0
outerArcsecPSW=200.0

radiusArcsecPMW=30.0
innerArcsecPMW=100.0
outerArcsecPMW=200.0

radiusArcsecPLW=42.0
innerArcsecPLW=100.0
outerArcsecPLW=200.0

resultPSW = annularSkyAperturePhotometry(image=mapPSW, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPSW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPSW)

resultPMW = annularSkyAperturePhotometry(image=mapPMW, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPMW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPMW)

resultPLW = annularSkyAperturePhotometry(image=mapPLW, fractional=1, \
	centerRA=str(centerRA), centerDec=str(centerDec), \
	radiusArcsec=radiusArcsecPLW, \
	innerArcsec=innerArcsecPSW, outerArcsec=outerArcsecPLW)

(yPSW,xPSW)=mapPSW.wcs.getPixelCoordinates(centerRA,centerDec)
(yPMW,xPMW)=mapPMW.wcs.getPixelCoordinates(centerRA,centerDec)
(yPLW,xPLW)=mapPLW.wcs.getPixelCoordinates(centerRA,centerDec)

mapSbPSW=mapPSW.image[int(yPSW),int(xPSW)]
mapSbPMW=mapPMW.image[int(yPMW),int(xPMW)]
mapSbPLW=mapPLW.image[int(yPLW),int(xPLW)]

mapErrPSW=mapPSW.error[int(yPSW),int(xPSW)]
mapErrPMW=mapPMW.error[int(yPMW),int(xPMW)]
mapErrPLW=mapPLW.error[int(yPLW),int(xPLW)]

srcPSW=sourceFitting(mapPSW,minX=xPSW-7.5,minY=yPSW-7.5,width=15.,height=15.,elongated=False)
srcPMW=sourceFitting(mapPMW,minX=xPMW-7.5,minY=yPMW-7.5,width=15.,height=15.,elongated=False)
srcPLW=sourceFitting(mapPLW,minX=xPLW-5,minY=yPLW-5,width=10.,height=10.,elongated=False)

fitSbPSW=srcPSW.getPeak()
fitSbPMW=srcPMW.getPeak()
fitSbPLW=srcPLW.getPeak()

#fitWidPSW=srcPSW.getsigmaArcsec()
#fitWidPMW=srcPSW.getsigmaArcsec()
#fitWidPLW=srcPSW.getsigmaArcsec()

#use same error as map
fitErrPSW=mapPSW.error[int(yPSW),int(xPSW)]
fitErrPMW=mapPMW.error[int(yPMW),int(xPMW)]
fitErrPLW=mapPLW.error[int(yPLW),int(xPLW)]

#print resultPSW["Results table"][0].data[2])
print 'Map value (PSW): %.4f +/- %.4f'%(mapSbPSW,mapErrPSW)
print 'Map value (PMW): %.4f +/- %.4f'%(mapSbPMW,mapErrPMW)
print 'Map value (PLW): %.4f +/- %.4f'%(mapSbPLW,mapErrPLW)
print '--'

(srcCorrPSW,bgCorrPSW,apCorrPSW)=apcorrfn.getApCorr(radiusArcsecPSW,innerArcsecPSW,outerArcsecPSW,band='PSW')
(srcCorrPMW,bgCorrPMW,apCorrPMW)=apcorrfn.getApCorr(radiusArcsecPMW,innerArcsecPMW,outerArcsecPMW,band='PMW')
(srcCorrPLW,bgCorrPLW,apCorrPLW)=apcorrfn.getApCorr(radiusArcsecPLW,innerArcsecPLW,outerArcsecPLW,band='PLW')

print '--'
print 'Correction (PSW):  %.3f (%.3f)'%(apCorrPSW,srcCorrPSW)
print 'Correction (PMW):  %.3f (%.3f)'%(apCorrPMW,srcCorrPMW)
print 'Correction (PLW):  %.3f (%.3f)'%(apCorrPLW,srcCorrPLW)

mapCorrSbPSW=mapSbPSW * apCorrPSW
mapCorrSbPMW=mapSbPMW * apCorrPMW
mapCorrSbPLW=mapSbPLW * apCorrPLW
fitCorrSbPSW=fitSbPSW * apCorrPSW
fitCorrSbPMW=fitSbPMW * apCorrPMW
fitCorrSbPLW=fitSbPLW * apCorrPLW

mapCorrErrPSW=mapErrPSW * apCorrPSW
mapCorrErrPMW=mapErrPMW * apCorrPMW
mapCorrErrPLW=mapErrPLW * apCorrPLW

fitCorrErrPSW=fitErrPSW * apCorrPSW
fitCorrErrPMW=fitErrPMW * apCorrPMW
fitCorrErrPLW=fitErrPLW * apCorrPLW

badmapSbPSW=mapSbPSW * srcCorrPSW
badmapSbPMW=mapSbPMW * srcCorrPMW
badmapSbPLW=mapSbPLW * srcCorrPLW
badmapErrPSW=mapErrPSW * srcCorrPSW
badmapErrPMW=mapErrPMW * srcCorrPMW
badmapErrPLW=mapErrPLW * srcCorrPLW

badfitSbPSW=fitSbPSW * srcCorrPSW
badfitSbPMW=fitSbPMW * srcCorrPMW
badfitSbPLW=fitSbPLW * srcCorrPLW
badfitErrPSW=fitErrPSW * srcCorrPSW
badfitErrPMW=fitErrPMW * srcCorrPMW
badfitErrPLW=fitErrPLW * srcCorrPLW

print '--'
print 'Map Bad Corr. (PSW):   %.4f +/- %.4f'%(badmapSbPSW,badmapErrPSW)
print 'Map Bad Corr. (PMW):   %.4f +/- %.4f'%(badmapSbPMW,badmapErrPMW)
print 'Map Bad Corr. (PLW):   %.4f +/- %.4f'%(badmapSbPLW,badmapErrPLW)
print '--'
print 'Fit Bad Corr. (PSW):   %.4f +/- %.4f'%(badfitSbPSW,badfitErrPSW)
print 'Fit Bad Corr. (PMW):   %.4f +/- %.4f'%(badfitSbPMW,badfitErrPMW)
print 'Fit Bad Corr. (PLW):   %.4f +/- %.4f'%(badfitSbPLW,badfitErrPLW)
print '--'
print 'Map Corrected (PSW):   %.4f +/- %.4f'%(mapCorrSbPSW,mapCorrErrPSW)
print 'Map Corrected (PMW):   %.4f +/- %.4f'%(mapCorrSbPMW,mapCorrErrPMW)
print 'Map Corrected (PLW):   %.4f +/- %.4f'%(mapCorrSbPLW,mapCorrErrPLW)
print '--'
print 'Fit Corrected (PSW):   %.4f +/- %.4f'%(fitCorrSbPSW,fitCorrErrPSW)
print 'Fit Corrected (PMW):   %.4f +/- %.4f'%(fitCorrSbPMW,fitCorrErrPMW)
print 'Fit Corrected (PLW):   %.4f +/- %.4f'%(fitCorrSbPLW,fitCorrErrPLW)
print '--'
print 'Timeline (PSW):    %.4f +/- %.4f'%(tlSbPSW,tlErrPSW)
print 'Timeline (PMW):    %.4f +/- %.4f'%(tlSbPMW,tlErrPMW)
print 'Timeline (PLW):    %.4f +/- %.4f'%(tlSbPLW,tlErrPLW)

fileOut=obsName+'_'+hex(obsid)+'.dat'
fileOut='%s_0x%X_%.1f_%.1f_%.1f_%.1f-%.1f.dat'%\
	(obsName,obsid,radiusArcsecPSW,radiusArcsecPMW,radiusArcsecPLW,\
	innerArcsecPSW,outerArcsecPSW)
print 'Writing to '+fileOut+'...'
fOut=open(dirPath+fileOut,'w')
line='#ObsId: 0x%X\n'%obsid
fOut.write(line)
line='#Name: %s\n'%obsName
fOut.write(line)
line='#RA/Dec [deg]: %.4f,%.4f\n'%(centerRA,centerDec)
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
line='Map SB: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(mapSbPSW, mapErrPSW, mapSbPMW, mapErrPMW, mapSbPLW, mapErrPLW)
fOut.write(line)

line='Fit SB: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(fitSbPSW, fitErrPSW, fitSbPMW, fitErrPMW, fitSbPLW, fitErrPLW)
fOut.write(line)

line='Timeline SB: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(tlSbPSW, tlErrPSW, tlSbPMW, tlErrPMW, tlSbPLW, tlErrPLW)
fOut.write(line)

line='SrcCorr Fact: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(srcCorrPSW, 0., srcCorrPMW, 0., srcCorrPLW, 0.)
fOut.write(line)

line='SrcCorr Map Flux: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(badmapSbPSW, badmapErrPSW, badmapSbPMW, badmapSbPMW, badmapSbPLW, badmapSbPLW)
fOut.write(line)

line='SrcCorr Fit Flux: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(badfitSbPSW, badfitErrPSW, badfitSbPMW, badfitErrPMW, badfitSbPLW, badfitErrPLW)
fOut.write(line)

line='ApCorr Fact: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(apCorrPSW, 0., apCorrPMW, 0., apCorrPLW, 0.)
fOut.write(line)

line='ApCorr Map Flux: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(mapCorrSbPSW, mapCorrErrPSW, mapCorrSbPMW, mapCorrErrPMW, mapCorrSbPLW, mapCorrErrPLW)
fOut.write(line)

line='ApCorr Fit Flux: %.4f, %.4f, %.4f, %.4f, %.4f,%.4f\n'%\
	(fitCorrSbPSW, fitCorrErrPSW, fitCorrSbPMW, fitCorrErrPMW, fitCorrSbPLW, fitCorrErrPLW)
fOut.write(line)

fOut.close()

#if __name__=='__main__':
#	maincode()
