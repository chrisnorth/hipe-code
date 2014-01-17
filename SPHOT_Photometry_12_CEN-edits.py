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
cal    = spireCal(pool='spire_cal_12_0')

beamCorrTable  = cal.phot.refs["ColorCorrBeam"].product
aperCorrTable  = cal.phot.colorCorrApertureList.refs[0].product
kCorrPsrcTable = cal.phot.colorCorrKList.refs[0].product
kCorrExtdTable = cal.phot.colorCorrKList.refs[1].product

#fwhm      = [18.2, 24.9, 36.9]
fwhm      = [17.0, 23.9, 35.2] #geometric mean ###CEN
peak      = [22, 30, 42]
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


#################### Running DAOphot on PSW ####################

# By default since HIPE 9, results are already corrected for aperture correction
srcDao = sourceExtractorDaophot(image=mapPsw, detThreshold=5.0, fwhm=fwhm[0], beamArea=beamAreaCorr[0])

# Colour correcting results for a source with alpha = +2.0 (default being -1)
srcDao['sources']['flux'].data = srcDao['sources']['flux'].data * kCorrPsrc[0]

#################    Running TimelineFitter    #################
################# using Sussextractor as prior #################

srcTimeline = sourceExtractorTimeline(input=obs.level1, array='PSW', rPeak=peak[0], inputSourceList=srcSussex)

# Colour correcting results for a source with alpha = +2.0 (default being -1)
srcTimeline['sources']['flux'].data = srcTimeline['sources']['flux'].data * kCorrPsrc[0]


############# Comparison with aperture photometry ##############
########## starting from point source calibrated maps ##########

# First, convert the map to Jy/pixel dividing by the beam area
mapPsrc = convertImageUnit(image=mapPsw, newUnit='Jy/pixel', beamArea=beamAreaCorr[0])

# Select/highlight mapConverted in the "Variables" view and then double click on
# "annularSkyAperturePhotometry" under "Applicable" in the "Tasks" view.
# Click on the central source, select "Centroiding", select some values for the radii
# (e.g. 30, 40, 50" for PSW) and then run the task. This is equivalent to running the following line
# (centerX and centerY will change depending on where you clicked with the mouse!):
resultPsrc = annularSkyAperturePhotometry(image=mapPsrc, centerX=91, centerY=104,\
	centroid=True, fractional=1, radiusArcsec=peak[0], innerArcsec=60.0, outerArcsec=90.0)

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
	resultPsrc['Results table'][3].data[2] * 1.e3 * kCorrPsrc[0] * beamCorr[0]) ##CEN - added beamCorr


############# Comparison with aperture photometry ##############
####### starting from extended emission calibrated maps ########

# First get the extended map from the obs (###CEN - new method)
mapExtd = obs.level2.refs["extdPSW"].product
# Then get the KPtoE values (from OM/DRG)
# These convert from MJy/sr (extended calibration) to Jy/beam (point src calibration)
kPtoE     = [91.6814, 51.4315, 23.9076]
# Divide map by that KPtoE and set units to Jy/beam
mapExtd["image"].data = mapExtd["image"].data / kPtoE[0]
mapExtd.setUnit('Jy/beam')
# Convert to Jy/pixel
mapExtd = convertImageUnit(mapExtd, newUnit='Jy/pixel', beamArea=beamAreaCorr[0])

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
	resultExtd['Results table'][3].data[2] * 1.e3 * kCorrPsrc[0] * beamCorr[0]) ##CEN - added beamCorr
########################## End of script #########################
##################################################################

print '------Results------'
print 'Timeline flux (in mJy): %.4f +/- %.4f'%(\
	srcTimeline['sources']['flux'].data[0],srcTimeline['sources']['fluxPlusErr'].data[0])
print 'Sussextractor flux (in mJy): %.4f +/- %.4f'%(\
	srcSussex['sources']['flux'].data[0],srcSussex['sources']['fluxPlusErr'].data[0])
print 'DaoPhot flux (in mJy): %.4f +/- %.4f'%(\
	srcDao['sources']['flux'].data[0],srcDao['sources']['fluxPlusErr'].data[0])
print 'ApCorr flux [psrc] (in mJy): %.4f +/- %.4f'%(\
	resultPsrc['Results table'][0].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0]* beamCorr[0],\
	resultPsrc['Results table'][3].data[2] * 1.e3 * kCorrPsrc[0]* beamCorr[0])
print 'ApCorr flux [extd] (in mJy): %.4f +/- %.4f'%(\
	resultExtd['Results table'][0].data[2] * 1.e3 * aperCorr[0] * kCorrPsrc[0]* beamCorr[0],\
	resultExtd['Results table'][3].data[2] * 1.e3 * kCorrPsrc[0]* beamCorr[0])
