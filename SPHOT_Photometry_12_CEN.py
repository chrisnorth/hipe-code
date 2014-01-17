###############################################################
#
# Title: SPHOT_Photometry_12
#
# Scope: this script presents different ways of doing
#        photometry with SPIRE - in this case of a 
#        point source with assumed spectral index alpha = +2.0
#
# Luca Conversi - 01/11/2013 - ver. 2
#
###############################################################


######################### Import data #########################

# Loading an observation of Gamma Dra from the HSA and its level1
try:
	obs
except:
	obs    = getObservation(1342195871, useHsa=1)
mapPsw = obs.level2.refs["psrcPSW"].product
mapPmw = obs.level2.refs["psrcPMW"].product
mapPlw = obs.level2.refs["psrcPLW"].product


#################### Correction parameters ####################

# Beam areas, FWHM rPeak and aperture correcton for PSW, PMW, and PLW, respectively.
# Values are obtained from the SPIRE calibration tree assuming a point source with spectral index alpha = 2
# N.B.: spire_cal_12_0 or later version is needed

alpha  = 2.0
arrays = ['PSW','PMW','PLW']
#import cal from working directory
import os
try:
	cal
except:
	cal = spireCal(jarFile=os.path.join(Configuration.getProperty('var.hcss.workdir'),'spire_cal_12_0.jar'))
#cal    = spireCal(pool='spire_cal_12_0')

beamCorrTable  = cal.phot.refs["ColorCorrBeam"].product
aperCorrTable  = cal.phot.colorCorrApertureList.refs[0].product
kCorrPsrcTable = cal.phot.colorCorrKList.refs[0].product
kCorrExtdTable = cal.phot.colorCorrKList.refs[1].product

#fwhm      = [18.2, 24.9, 36.9]
fwhm      = [17.0, 23.9, 35.2] #geometric mean
peak      = [22, 30, 42]
kPtoE    = [90.681, 51.432, 23.908]
beamArea  = [beamCorrTable.meta['beamPswArc'].double,\
	beamCorrTable.meta['beamPmwArc'].double, \
	beamCorrTable.meta['beamPlwArc'].double]
beamCorr  = []
kCorrPsrc = [] 
kCorrExtd = [] 
aperCorr  = []

for array in arrays:
	beamCorr.append(beamCorrTable.getAlphaCorrection(alpha, array))
	kCorrPsrc.append(kCorrPsrcTable.getAlphaCorrection(alpha, array))
	kCorrExtd.append(kCorrPsrcTable.getAlphaCorrection(alpha, array))
	aperCorr.append(aperCorrTable.getApertColorCorrection(alpha, array))

beamAreaCorr = Double1d(3)
for i in range(3):
	beamAreaCorr[i] = beamArea[i] / beamCorr[i]


############### Running Sussextractor on PSW map ##############

srcSussex = sourceExtractorSussextractor(image=mapPsw, detThreshold=5.0, fwhm=fwhm[0])

# Colour correcting results for a source with alpha = +2.0 (default being -1)
srcSussex['sources']['flux'].data = srcSussex['sources']['flux'].data * kCorrPsrc[0]


mapPsrc2 = obs.level2.refs["extdPSW"].product
mapPsrc2["image"].data = mapPsrc2["image"].data / kPtoE[0]
mapPsrc2.setUnit('Jy/beam')

#################### Running DAOphot on PSW ####################

# By default since HIPE 9, results are already corrected for aperture correction
srcDao = sourceExtractorDaophot(image=mapPsrc2, detThreshold=5.0, fwhm=fwhm[0], beamArea=beamArea[0])
srcDaoCorr = sourceExtractorDaophot(image=mapPsrc2, detThreshold=5.0, fwhm=fwhm[0], beamArea=beamAreaCorr[0])
srcDaoNoAp = sourceExtractorDaophot(image=mapPsrc2, detThreshold=5.0, fwhm=fwhm[0], beamArea=beamArea[0],doApertureCorrection=False)
srcDaoCorrNoAp = sourceExtractorDaophot(image=mapPsrc2, detThreshold=5.0, fwhm=fwhm[0], beamArea=beamAreaCorr[0],doApertureCorrection=False)

# Colour correcting results for a source with alpha = +2.0 (default being -1)
srcDao['sources']['flux'].data = srcDao['sources']['flux'].data * kCorrPsrc[0] * beamCorr[0]
srcDaoCorr['sources']['flux'].data = srcDaoCorr['sources']['flux'].data * kCorrPsrc[0]
srcDaoNoAp['sources']['flux'].data = srcDaoNoAp['sources']['flux'].data * kCorrPsrc[0] * aperCorr[0] * beamCorr[0]
srcDaoCorrNoAp['sources']['flux'].data = srcDaoCorrNoAp['sources']['flux'].data * kCorrPsrc[0] * aperCorr[0]

#################    Running TimelineFitter    #################
################# using Sussextractor as prior #################

srcTimeline = sourceExtractorTimeline(input=obs.level1, array='PSW', rPeak=peak[0], inputSourceList=srcSussex)

# Colour correcting results for a source with alpha = +2.0 (default being -1)
srcTimeline['sources']['flux'].data = srcTimeline['sources']['flux'].data * kCorrPsrc[0]


############# Comparison with aperture photometry ##############
########## starting from point source calibrated maps ##########

# First, convert the map to Jy/pixel dividing by the pipeline beam area 
mapPsrc = convertImageUnit(image=mapPsw, newUnit='Jy/pixel', beamArea=beamArea[0])

mapPsrc2Ext = convertImageUnit(image=mapPsrc2, newUnit='Jy/pixel', beamArea=beamArea[0])
# Select/highlight mapPsrc in the "Variables" view and then double click on
# "annularSkyAperturePhotometry" under "Applicable" in the "Tasks" view.
# Click on the central source, select "Centroiding", select some values for the radii
# (e.g. 22, 60, 90" for PSW) and then run the task. This is equivalent to running the following line
# (centerX and centerY will change depending on where you clicked with the mouse!):
resultPsrc = annularSkyAperturePhotometry(image=mapPsrc, centerX=91, centerY=104,\
	centroid=True, fractional=1, radiusArcsec=peak[0], innerArcsec=60.0, outerArcsec=90.0)
resultPsrc2 = annularSkyAperturePhotometry(image=mapPsrc2Ext, centerX=91, centerY=104,\
	centroid=True, fractional=1, radiusArcsec=peak[0], innerArcsec=60.0, outerArcsec=90.0)
#
## Alternatively, you can also specify a [RA, Dec] position,
## e.g. the one obtained above using "Get coordinates"
#raCentre  = str(displayWorldCoords[0])
#decCentre = str(displayWorldCoords[1])
#
#resultPsrc = annularSkyAperturePhotometry(image=mapPsrc, centerRA=raCentre, centerDec=decCentre,\
#	centroid=True, fractional=1, radiusArcsec=peak[0], innerArcsec=60.0, outerArcsec=90.0)

# Correcting result for aperture and colour corrections for a source with alpha = +2.0 (default being -1)
print 'Source flux (in mJy) for PSW array is: %f +/- %f'%(\
	resultPsrc['Results table'][0].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0] * beamCorr[0],\
	resultPsrc['Results table'][3].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0] * beamCorr[0])

print 'Source flux 2 (in mJy) for PSW array is: %f +/- %f'%(\
	resultPsrc2['Results table'][0].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0] * beamCorr[0],\
	resultPsrc2['Results table'][3].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0] * beamCorr[0])
############# Comparison with aperture photometry ##############
####### starting from extended emission calibrated maps ########

# First, convert the map from MJy/sr to Jy/pixel (no need of the beam area so far...)
mapExtd = convertImageUnit(image=obs.level2.refs["extdPSW"].product, newUnit='Jy/pixel')

# Select/highlight mapConverted in the "Variables" view and then double click on
# "annularSkyAperturePhotometry" under "Applicable" in the "Tasks" view.
# Click on the central source, select "Centroiding", select some values for the radii
# (e.g. 30, 40, 50" for PSW) and then run the task. This is equivalent to running the following line
# (centerX and centerY will change depending on where you clicked with the mouse!):
resultExtd = annularSkyAperturePhotometry(image=mapExtd, centerX=91, centerY=104,\
	centroid=True, fractional=1, radiusArcsec=peak[0], innerArcsec=60.0, outerArcsec=90.0)

## Alternatively, you can also specify a [RA, Dec] position,
## e.g. the one obtained above using "Get coordinates"
#raCentre  = str(displayWorldCoords[0])
#decCentre = str(displayWorldCoords[1])
#
#resultExtd = annularSkyAperturePhotometry(image=mapExtd, centerRA=raCentre, centerDec=decCentre,\
#	centroid=True, fractional=1, radiusArcsec=peak[0], innerArcsec=60.0, outerArcsec=90.0)

# Correcting result for aperture and colour corrections for a source with alpha = +2.0 (default being -1)
# N.B.: results must also be multiplied by the beam correction for alpha = +2.0,
# as the original maps in MJy/sr were calculated under the assumption of alpha = -1.0
print 'Source flux (in mJy) for PSW array is: %f +/- %f'%(\
	resultExtd['Results table'][0].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0] * beamCorr[0],\
	resultExtd['Results table'][3].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0] * beamCorr[0])

########################## End of script #########################
##################################################################


##Test numbers
print '------Results------'
print 'Timeline flux (in mJy): %.4f +/- %.4f'%(\
	srcTimeline['sources']['flux'].data[0],srcTimeline['sources']['fluxPlusErr'].data[0])
print 'Sussextractor flux (in mJy): %.4f +/- %.4f'%(\
	srcSussex['sources']['flux'].data[0],srcSussex['sources']['fluxPlusErr'].data[0])
print 'DaoPhot flux (in mJy): %.4f +/- %.4f'%(\
	srcDao['sources']['flux'].data[0],srcDao['sources']['fluxPlusErr'].data[0])
print 'DaoPhot flux [corr beam] (in mJy): %.4f +/- %.4f'%(\
	srcDaoCorr['sources']['flux'].data[0],srcDaoCorr['sources']['fluxPlusErr'].data[0])
print 'DaoPhot flux [no auto apcorr] (in mJy): %.4f +/- %.4f'%(\
	srcDaoNoAp['sources']['flux'].data[0],srcDaoNoAp['sources']['fluxPlusErr'].data[0])
print 'DaoPhot flux [corr beam, no auto apcorr] (in mJy): %.4f +/- %.4f'%(\
	srcDaoCorrNoAp['sources']['flux'].data[0],srcDaoCorrNoAp['sources']['fluxPlusErr'].data[0])
print 'ApCorr flux [psrc] (in mJy): %.4f +/- %.4f'%(\
	resultPsrc2['Results table'][0].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0] * beamCorr[0],\
	resultPsrc2['Results table'][3].data[2] * 1.e3 * kCorrPsrc[0])
print 'ApCorr flux [extd] (in mJy): %.4f +/- %.4f'%(\
	resultExtd['Results table'][0].data[2] * 1.e3 * aperCorr[0] * kCorrExtd[0] * beamCorr[0],\
	resultExtd['Results table'][3].data[2] * 1.e3 * kCorrExtd[0])
