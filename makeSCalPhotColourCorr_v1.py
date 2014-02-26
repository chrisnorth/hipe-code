#===============================================================================
# 
#  This routine computes the following calibration datasets for Herschel-SPIRE
#  * Radial beam profiles (RadialCorrBeam):
#    Radial profiles of the beam for each band, separated into a freqeuncy-dependent
#    "core" section and a freqeuncy-independent "constant" section. The cumulative
#    solid angle as a function of radius, normalised to the final value is also
#    output. Metadata values are produced for the solid angle as measured on
#    Neptune, the spectral index of Neptune, the beam FWHM frequency-dependence,
#    the effective frequency at which these beam profiles can be assumed to apply,
#    and the beam solid angle used in the pipeline.
#
#  * Colour correction tables (ColorCorrK):
#    Colour correction tables to convert from the standard pipeline flux densities,
#    (which are quoted for a nu*F_nu=const. spectrum) into flux densities for a
#    range of source spectra. These produce monochromatic flux densities at the
#    SPIRE reference wavelengths of 250, 350 and 500 micron. The source spectra
#    include power law spectra, and modified black body spectra with a range of
#    emissivities. Metadata values are produced for the flux conversion parameters
#    applied to the pipeline (the "K4P" and "K4E" parameters).
#
#  * Beam colour correction tables (ColorCorrBeam):
#    Colour correction tables to convert from the pipeline beam area (which assumes
#    a source spectrum of nu*F_nu=const.) into the effective beam solid angle for
#    a range of source spectra. The source spectra include power law spectra,
#    and modified black body spectra with a range of emissivities. The effective
#    beam area is calculated as [pipeline beam area] / [beam correction factor].
#
#  * Aperture correction tables (ColorCorrAperture):
#    Aperture correction tables to be applied after aperture photometry, assuming
#    the standard aperture and annulus radii. These aperture corrections are provided
#    for a range of source spectra. The sourcce spectra include power law spectra
#    and modified black body spectra with a range of emissivities.
#
#  * HFI-SPIRE Colour Correction factors (ColorCorrHfi):
#    A table of colour correction factors for a range of
#    realistic celestial background spectra. The factors transfer standard 
#    monochromatic fluxdensities from Planck-HFI all-sky maps at 857 and 545 GHz
#    which are quoted for a nu*F_nu = const. spectrum into Herschel-SPIRE 
#    monochromatic flux densities with the same definition but reference 
#    wavelengths at 250, 350, and 500 micron. It is assumed that the actual 
#    celestial background is described by the Planck-Law multiplied by frequency 
#    to the power of beta which is set to 1.8. In addition the flux ratios as 
#    measured by HFI between the 545 and 857 Gz filters is given such that the 
#    colour corrections can be determined directly as a function of the flux 
#    ratios in those maps.
# 
#  Inputs:
#    General:
#      An existing calibration tree.
#      RSRF profiles from SPIRE calibration tree
#      New version number
#      Input and Output directories
#      Which products to update
#    HFI RIMO files (containing HFI RSRF profiles), downloadable at:
#        http://pla.esac.esa.int/pla/aio/product-action?DOCUMENT.DOCUMENT_ID=HFI_RIMO_R1.10.fits
#    RSRF profiles from SPIRE calibration tree
#    Exponent for modified Blackbody to model dust emission background colour
#    Name and version of output file
# 
#  Output:
#    FITS binary table file containing 5 columns:
#      'Temp'         dust temperature of model in Kelvin
#      'ratio545_857' ratio between HFI 545 GHz and 857 GHz maps
#      'k545toPLW'    colour correction factor 545 GHz map to PLW map
#      'k857toPMW'    colour correction factor 857 GHz map to PMW map
#      'k857toPSW'    colour correction factor 857 GHz map to PSW map
#    K4 factors for point sources as meta data parameters
# 
#  Calculations:
#   K-Correction factors to go from HFI to SPIRE maps. The colour correction 
#   computed here:
#   1. starts from HFI maps computed under the assumption of having extended 
#      source with alpha=-1 spectrum
#   2. "converts" them to sky flux density assuming a grey body spectrum
#   3. "converts" them again to SPIRE wavebands using the same grey body 
#      assumption
#   4. finally changes them to the SPIRE extended emission calibration 
#      assuming an alpha=-1 spectrum
#  
#   In other terms: 
#   kHFItoSPIRE = 
#      1/K4E_HFI * K4E_MOD_HFI / K4E_MOD_SPIRE * K4E_SPIRE * fSky_SPIRE / fSky_HFI
#  
#   We also calculate the ratio of the Planck maps as measured by the pipeline,
#   i.e. under the assumption alpha=-1, but assuming a grey body of temperature T 
#   and spectral index beta. Hence:
#   R_HFI = fPipe_545 / fPipe_857 =
#         = (fSky_GB_545 / k545 * k4E_545) / (fSky_GB_857 / k857 * k4E_857) 
#  
# 
# 
#  Edition History
#  Luca Conversi   - 13/Jul/2012 - 0.1: First test version
#  Bernhard Schulz - 18/Aug/2012 - 0.2: Separated out and added plots
#  Bernhard Schulz - 22/Aug/2012 - 0.3: Added comments, reformatted, file output added,
#                                       renamed kCorr to hpXcalKcorr
#  Bernhard Schulz - 24/Aug/2012 - 0.4: Addded point source RSRF to plot and calculate 
#                                       K4 factors for point sources
#  Bernhard Schulz - 27/Aug/2012 - 0.5: Cleaned up code and added more comments, added 
#                                       point source and extended source colour correction 
#                                       factors for standard spectrum as metadata
#  Bernhard Schulz - 28/Aug/2012 - 1.0: Renamed ratioPMWoverPLW to kPMWtoPLW and included 
#                                       into output file as column PMWtoPLW
#  Luca Conversi   - 16/Oct/2012 - 1.1: Change column names and order;
#                                       computing k857toPSW instead of kPMWtoPSW
#  Luca Conversi   - 21/Mar/2013 - 2.0: Adding the possibility of constructing a
#                                       colour correction table with T fixed and variable beta
#                                       Other minor changes (e.g. toggles plots on/off)
#  Bernhard Schulz - 05/Apr/2013 - 2.1: Change to tabular integration, 
#                                       added gamma parameter to header.
#  Bernhard Schulz - 12/Apr/2013 - 2.2: Common frequency interpolation and aperture 
#                                       efficiencies included
#  Luca Conversi   - 30/May/2013 - 2.3: Making file indipendent from external librabry
#                                       Changing script to read HFI RSRF from "RIMO" file
#                                       (Table of the performance and instrumental characteristics
#                                       as downloaded from the Planck legacy archive)
#                                       Cleaning of filenames, paths, etc.
#  Luca Conversi   - 18/Jun/2013 - 2.4: Update Ks computation following C. North inputs:
#                                       now gamma is kept separatated from SPIRE RSRFs definition
#                                       and included in the denominator of the Ks calculations
#                                       Correct typo in table metadata when fixBeta = False
#  Chris North     - 15/Jan/2014 - 2.5: Added full beam model treatment
#                                       Changed integration method to TrapezoidalIntegrator
#                                       Changed RSRF/Aperture Efficiency interpolation to Linear
#
#===============================================================================

import os
scriptVersionString = "makeSCalPhotColorCorr.py $Revision: 1.1 $"

# define the start and end dates for the product
# starting at beginning of PFM1
df  = java.text.SimpleDateFormat("yyyy.MM.dd/HH:mm:ss/z")
startDate = df.parse("2005.02.22/00:00:00/GMT")
endDate   = df.parse("2020.01.01/00:00:00/GMT")
dateNow = java.util.Date().toString()

# Input directory
dataDir = Configuration.getProperty('var.hcss.workdir')
#dataDir = "//disks//winchester2//calibration_data//"
# Output directory
directory = Configuration.getProperty('var.hcss.workdir')
#directory = "..//..//..//..//..//..//data//spire//cal//SCal"

dirRadialCorrBeam=os.path.join(directory,'Phot/SCalPhotRadialCorrBeam')
dirColorCorrK=os.path.join(directory,'Phot/SCalPhotColorCorrK')
dirColorCorrBeam=os.path.join(directory,'Phot/SCalPhotColorCorrBeam')
dirColorCorrHfi=os.path.join(directory,'Phot/SCalPhotColorCorrHfi')
dirColorCorrAperture=os.path.join(directory,'Phot/SCalPhotColorCorrAperture')
dirFluxConv=os.path.join(directory,'Phot/SCalPhotFluxConv')

#-------------------------------------------------------------------------------
# Set version numbers

# Old calibration version
oldCalVersion = "spire_cal_12_0"

# New table version
version = "2.8.1"
formatVersion="1.0"

#set whether to print more info
verbose=True
# Choose if you want to plot the computed RSRFs and colour correction parameters
plot = True

#-------------------------------------------------------------------------------
# Set which colour correction tables to re-compute
calcRadialCorrBeam =    True	# recalculate new beam profiles
calcSpireEffFreq =      True	# recalculate SPIRE effective frequencies
calcColorCorrK =        True	# recalculate K-params and beam colour correction parameters
calcColorCorrAperture = False	# recalculate Aperture Correction
calcColorCorrHfi =      True	# recalculate HFI X-cal parameters

#-------------------------------------------------------------------------------
# Check dependencies
#overrideDepCheck = True #set to override dependency check
#depCheck=True
#if calcRadialCorrBeam:
#	if not calcSpireEffFreq:
#		depCheck=False
#		print "***WARNING: new RadialCorrBeam will change effective frequencies***"
#	if not calcColorCorrK:
#		depCheck=False
#		print "***WARNING: new RadialCorrBeam will change ColorCorrK_extended and ColorCorrBeam***"
#	if not calcColorCorrHfi:
#		depCheck=False
#		print "***WARNING: new RadialCorrBeam will change ColorCorrHfi***"
#	if not calcColorCorrAperture:
#		depCheck=False
#		print "***WARNING: new RadialCorrBeam will change ColorCorrAperture***"
#
#if calcSpireEffFreq:
#	if not calcRadialCorrBeam:
#		depCheck=False
#		print "***WARNING: new effective frequencies will change RadialCorrBeam metadata***"
#	if not calcColorCorrK:
#		depCheck=False
#		print "***WARNING: new effective frequencies will change ColorCorrK_extended and ColorCorrBeam***"
#	if not calcColorCorrHfi:
#		depCheck=False
#		print "***WARNING: new effective frequencies will change ColorCorrHfi***"
#	if not calcColorCorrAperture:
#		depCheck=False
#		print "***WARNING: new effective frequencies will change ColorCorrAperture***"
#
#if calcColorCorrK:
#	if not calcColorCorrHfi:
#		depCheck=False
#		print "***WARNING: new ColorCorrBeam & ColorCorrK will change ColorCorrHfi***"
#	if not calcRadialCorrBeam:
#		depCheck=False
#		print "***WARNING: new ColorCorrBeam will change RadialCorrBeam metadata***"
#	if not calcColorCorrAperture:
#		depCheck=False
#		print "***WARNING: new ColorCorrBeam will change ColorCorrAperture***"
#
##stop if not overriden
#if not overrideDepCheck:
#	assert depCheck==True,"***ERROR: Update dependencies not met. Stopping."
#else:
#	if depCheck==False:
#		print "***WARNING: Update dependencies not met. Overriding."

#
#-------------------------------------------------------------------------------
#===============================================================================
#=====                         SET INPUT PARAMETERS                        =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Input parameters for Beam

#beam profile version
beamNewVersion = "1.1"
beamNewFileConstant = 'beamProfs_constant_v%s.csv'%beamNewVersion
beamNewFileCore = 'beamProfs_core_v%s.csv'%beamNewVersion

# Beam areas at effective frequency in Sr
# Due to beam model, this is the same as the value measured on Neptune.
# These values are online at https://nhscsci.ipac.caltech.edu/sc/index.php/Spire/PhotBeamProfileDataAndAnalysis
arcsec2Sr = (Math.PI/(60.*60.*180))**2
spireAreaEffFreq = {"PSW":450.*arcsec2Sr, "PMW":795.*arcsec2Sr, "PLW":1665.*arcsec2Sr}
# Neptune spectral index (from ESA4 model)
alphaNep={"PSW":1.29, "PMW":1.42, "PLW":1.47}
nepVersion="ESA4"

# Exponent of powerlaw describing FWHM dependence on frequency
# FWHM ~ freq**gamma
gamma = -0.85 

# Use full beam treatment
simpleBeam=False

#-------------------------------------------------------------------------------
# Input parameters for colour correction
# range of alphas to compute colour corrections for
alphaK=[-4.,-3.5,-3.,-2.5,-2.,-1.5,-1.,-0.5,0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.]
#alphaK=[-4,-1,0,2]
# range of beta and temp to calculate colour corrections for
betaK=[0.,0.5,1.,1.25,1.5,1.75,2.,2.5,3.]
tempK=range(3,300)

#-------------------------------------------------------------------------------
# Input parameters for aperture correction
apPhotRad={"PSW":22.,"PMW":30.,"PLW":45.}
apPhotBGRad=[60.,90.]

#-------------------------------------------------------------------------------
# Input parameters for ColorCorrHfi

# Choose if you want to keep constant either:
# - grey body spectral index beta: fixBeta = True [default]
# - or grey body temperature: fixBeta = False
fixBeta = True

# Grey body spectral index used for simulated background dust emission spectrum
# Only used if fixBeta = True
beta0 = 1.8  

# Grey body temeprature (in K) used for simulated background dust emission spectrum
# Only used if fixBeta = False
temp0 = 20.0

# HFI RIMO file
hfiRIMOFileName = 'HFI_RIMO_R1.10.fits'
hfiRIMOFile = os.path.join(dataDir,hfiRIMOFileName)

# Name of output table
#outFilename = os.path.join(dataDir,'SpireHfiColourCorrTab_v%3.1f.fits'%version)

#-------------------------------------------------------------------------------
# Loading physical constants
from herschel.share.unit import *
h = Constant.H_PLANCK.value
k = Constant.K_BOLTZMANN.value
c = Constant.SPEED_OF_LIGHT.value

# Frequency raster of common frequency grid
deltaNu = 0.1e9		# 0.1 GHz
nuMin   = 300.e9
nuMax   = 1800.e9
nNu     = FIX((nuMax-nuMin)/deltaNu)
freq    = Double1d(range(nNu)) * deltaNu + nuMin

# Define the temperatures range and number of rows in output table
if fixBeta == True:
	Tmin = 5.		# Min. modified blackbody temperature
	Tmax = 300.		# Max. modified blackbody temperature
	ndiv = 500		# Number of elements in temperature vector

# Define the spectral index range and number of rows in output table
if fixBeta == False:
	Bmin = 1.		# Min. modified blackbody temperature
	Bmax = 3.		# Max. modified blackbody temperature
	ndiv = 200		# Number of elements in temperature vector

# SPIRE band names
spireBands=["PSW","PMW","PLW"]

# Three SPIRE filter reference frequencies for PSW, PMW, PLW respectively
spireRefWl = {"PSW":250.*1e-6, "PMW":350.*1.e-6, "PLW":500.*1.e-6}
spireRefFreq = {}
for band in spireBands:
	spireRefFreq[band] = c/spireRefWl[band]

# Two HFI filter reference frequencies for the 857 and 545 GHz filters respectively
hfiBands= ["545","857"]
hfiRefFreq = Double1d([857.,545.])*1e9       

#-------------------------------------------------------------------------------
#===============================================================================
#=====                LOAD CALIBRATION AND FILTER FUNCTIONS                =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# load calibration
calibration = spireCal(jarFile=os.path.join(dataDir,'%s.jar'%oldCalVersion))
#calibration = SpireCal.getInstance()
calVersion=calibration.getVersion()
print 'Using SPIRE Calibration Tree version %s'%(calVersion)
#-------------------------------------------------------------------------------
# Load SPIRE filter functions

rsrf        = calibration.phot.getProduct("Rsrf")
rsrfVersion = rsrf.getVersion()
print 'Using SPIRE RSRF version %s from %s'%(rsrfVersion,calVersion)
spireFreq   = rsrf.getFrequency()*1e9	# Frequency in Hz
#indexes of freq in rsrf
ixR = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))

# Read aperture efficiency table from calibration tree
spireApEffCal = calibration.phot.getProduct("ApertureEfficiency")
apEffVersion=spireApEffCal.getVersion()
print 'Using SPIRE Aperture Efficiency version %s from %s'%(apEffVersion,calVersion)
spireApEffFreq = spireApEffCal.getApertEffTable()["frequency"].data * 1e9 #comes in [GHz]
#indexes of freq in apEff
ixA = freq.where((freq>=MIN(spireApEffFreq)) & (freq<=MAX(spireApEffFreq)))

# spire RSRF only
spireFiltOnly={}
# spire RSRF * ApEff
spireFilt={}

#interpolate to freq array
for band in spireBands:
	#create Rsrf and ApEff interpolation objects
	interpRsrf = LinearInterpolator(spireFreq, rsrf.getRsrf(band))
	interpAp = LinearInterpolator(spireApEffFreq, spireApEffCal.getApertEffTable()[band].data)
	#make arrays for final objects
	spireFiltOnly[band] = Double1d(nNu)
	spireFilt[band] = Double1d(nNu)
	#interpolate Rsrf to freq array
	spireFiltOnly[band][ixR] = interpRsrf(freq[ixR])
	#copy into Rsrf*ApEff array
	spireFilt[band] = spireFiltOnly[band].copy()
	spireFilt[band][ixA] = spireFilt[band][ixA] * interpAp(freq[ixA])

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           DEFINE FUNCTIONS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

# Set up functions to calculate beam profile, effective frequency and effective area
# Function list:
#   * spireMonoBeam: Calculate monochromatic beam profile and area at a given frequency
#   * spireMonoAreas: Calculate monochromatic beam areas over a range of frequencies
#   * spireMonoAreasSimple: Calculate monochromatic beam areas over a range of frequencies
#        using simple beam model [beam area ~ (freq/effFreq)**2gamma]
#   * spireEffArea: Calculate the effective beam area for a given spectrum
#   * spireFindEffFreq: Calculate the effective frequency for SPIRE
#   * spireEffBeam: Calculate the effective beam profile, area and beam map for SPIRE
#   * hpXcalKcorr: Calculate K-correction parameters for given spectrum & source type

#-------------------------------------------------------------------------------
# Calculate monochromatic beam profile and area at a given frequency
def spireMonoBeam(freqx,beamRad,beamProfs,beamConst,effFreq,gamma,array):
	"""
	========================================================================
	Implements the full beam model to generate the monochromatic beam profile
	and corresponding monochromatic beam solid angle at a given frequency.

	Inputs:
	  freqx:     (float) frequency [Hz] for which monochromatic beam
	               profile and area should be calculated
	  beamRad:   (array float) radius vector from the input beam profiles
	  beamProfs: (dataset) PhotRadialCorrBeam object
	               [used to retrieve core beam profile]
	  beamConst: (array float) constant part of beam profile for "array"
                       [passed to prevent repeated calls]
	  effFreq:   (float) effective frequency [Hz] of array
	  gamma:     (float) Exponent of powerlaw describing FWHM dependence
	               on frequency
	  array:     (string) spire array ('Psw'|'Pmw'|'Plw')

	Outputs:     (list of objects)
	  [0]:       (float) Beam area [arcsec^2] at frequency freqx
	  [1]:       (array float) Monochromatic beam profile at frequency freqx

	Calculation:
          Scales the core beam profile width as (freqx/effFreq)^gamma.
	  Queries the calibration file to generate new core beam profile.
	  Uses constant beam profile where it is larger than core profile.
	  Integrates over radius to caluclate beam area.

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
	  
	2013/12/19  C. North  initial version

	"""

	#calculate the "scaled" radius, as nu^-gamma
	radNew=beamRad*(freqx/effFreq)**-gamma
	maxRad=max(beamRad)
	nRad=len(beamRad)
	#ensure it doesn't go out of range
	radNew[radNew.where(radNew > maxRad)]=maxRad
	#get the corresponding beam profiles
	beamNew=Double1d(nRad)
	for r in range(nRad):
		beamNew[r]=beamProfs.getCoreCorrection(radNew[r],array)
	#apply the "constant" beam where appropriate
	#beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data
	isConst=beamNew.where(beamNew < beamConst)
	beamNew[isConst]=beamConst[isConst]

	#integrate to get solid angle (in arcsec^2)
	
	beamInterp=LinearInterpolator(beamRad,beamNew * 2. * Math.PI * beamRad)
	integrator=TrapezoidalIntegrator(0,maxRad)
	beamMonoArea=integrator.integrate(beamInterp)

	return(beamMonoArea,beamNew)

#-------------------------------------------------------------------------------
# Calculate monochromatic beam areas over a range of frequencies
def spireMonoAreas(freq,beamProfs,effFreq,gamma,array,freqFact=100):

	"""
	========================================================================
	Generates array of monochromatic beam areas over frequency range by
	calculating over a sparser array and interpolating

	Inputs:
	  freq:      (array float) frequency vector [Hz] for which monochromatic
	               beams areas should be Beam Radial Profiles calculated
	  beamProfs: (dataset) PhotRadialCorrBeam object from calibration tree
	  effFreq:   (float) effective frequency [Hz] of array
	  gamma:     (float) Exponent of powerlaw describing FWHM dependence
	               on frequency
	  array:     (string) spire array ('Psw'|'Pmw'|'Plw')
	  freqFact:  (int) Factor by which to reduce size of freq.
	               OPTIONAL. Default=100.

	Outputs:     
	             (array float) Monochromatic Beam area [sr] at frequencies
	                corresponding to freq

	Calculation:
          Geneates sparse frequency array of full range
	  Uses spireMonoBeam to calculate monochromatic beam area at sparse freqs
	  Interpolates to full frequency grid

	Dependencies:
	  spireMonoBeam
	  herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
	  
	2013/12/19  C. North  initial version

	"""

	#set up a sparser range of frequencies (otherwise it takes too long)
	nNu=len(freq)
	nNuArea=nNu/freqFact + 1
	#array of indices of full frequency array to use
	iNuArea=Int1d(range(nNuArea))*freqFact
	iNuArea[-1]=nNu-1

	#set up arrays
	beamMonoFreqSparse=Double1d(nNuArea)
	beamMonoAreaSparse=Double1d(nNuArea)

	#get beam radius array from calibration table
	beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
	#get constant beam profile from calibration table
	beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data

	# calculate at sparse frequencies
	for fx in range(nNuArea):
		#get corresponding index in full frequency array
		f=iNuArea[fx]
		#populate frequency array
		beamMonoFreqSparse[fx]=freq[f]
		#populate beam area array
		beamMonoAreaSparse[fx]=spireMonoBeam(freq[f],beamRad,beamProfs,beamConst,effFreq,gamma,array)[0]

	# interpolate to full frequency array and convert to Sr
	beamInterp=CubicSplineInterpolator(beamMonoFreqSparse,beamMonoAreaSparse)
	beamMonoArea=beamInterp(freq)*arcsec2Sr #in sr
	
	return(beamMonoArea)

##-------------------------------------------------------------------------------
## Calculate monochromatic beam areas over a range of frequencies using simple
##   beam model (beam area ~ (freq/effFreq)**2gamma]
#def spireMonoAreasSimple(freq,effFreq,areaEffFreq,gamma):
#
#	"""
#	========================================================================
#	Generates array of monochromatic beam areas over frequency range by
#	calculating over a sparser array and interpolating
#
#	Inputs:
#	  freq:        (array float) frequency vector [Hz] for which monochromatic
#	                 beams areas should be calculated
#	  effFreq:     (float) effective frequency [Hz] of array
#	  areaEffFreq: (float) beam area at effective frequency 
#	  gamma:       (float) Exponent of powerlaw describing FWHM dependence
#	                 on frequency
#	Outputs:     
#	             (array float) Monochromatic Beam area at frequencies
#	                corresponding to freq [same units as areaEffFreq]
#
#	Calculation:
#	  Uses simple power law to scale beam area with frequency
#
#	Dependencies:
#	  None
#	
#	2014/01/07  C. North  initial version
#
#	"""
#
#	beamMonoArea = areaEffFreq * (freq/effFreq)**(2.*gamma)
#	return(beamMonoArea)

#-------------------------------------------------------------------------------
# Calculate the effective beam area for a given spectrum
def spireEffArea(freq, transm, monoArea, BB=False, temp=20.0, beta=1.8, alpha=-1.0):
	"""
	========================================================================
	Calculate the effective beam area for a source of a given spectrum

	Inputs:
	  freq:       (array float) frequency vector corresponding to RSRF values [Hz]
	  transm:     (array float) relative spectral response (RSRF) corresponding to freq
	                Note that this should *not* include the aperture efficiency
	  monoArea:   (array float) monochromatic beam solid angle corresponding
	                to frequencies in freq
	  BB:         (boolean) spectral function to use for source spectrum:
	                'True': a modified black body
	                'False' a power-law with exponent alpha=-1
	                OPTIONAL. Default=False
	  temp:       (float) Dust/sky temperature (if BB=True)
	                OPTIONAL. Default=20.0
	  beta:       (float) Dust/sky spectral index (if BB=True)
	                OPTIONAL. Default=1.8
	  alpha:      (float) Exponent of power-law sky background model (if BB=False)
	                OPTIONAL. Default=-1

	Outputs:     
	            (float) Beam area for given spectrum, in same units as monoArea

	Calculation:
	  Calculates the source spectrum (either modifies black body or power law)
          Multiplies the monochromatic beam area by RSRF and source spectrum
	  Integrates over frequency
	  Normalises by integral over frequency of RSRF and source spectrum

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

	2013/12/19  C. North  initial version

	"""	
	#
	# Calculate sky background model
	#
	if BB == 1:
		#print temp,beta,c,h,k
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha

	# Integrate monochromatic area over frequency, weighted by rsrf and fSky
	numInterp=LinearInterpolator(freq,transm * fSky * monoArea)
	denomInterp=LinearInterpolator(freq,transm * fSky)
	minFreq=min(freq)
	maxFreq=max(freq)
	integrator=TrapezoidalIntegrator(minFreq,maxFreq)
	numInteg=integrator.integrate(numInterp)
	denomInteg=integrator.integrate(denomInterp)
	effArea = numInteg / denomInteg

	return(effArea)


#-------------------------------------------------------------------------------
# Calculate the effective frequency for SPIRE
def spireFindEffFreq(freq, rsrf, beamProfs, effFreqInit, gamma,
  areaNep, alphaNep, array, simpleBeam=False, freqFact=500,
  initRange=0.01, reqPrec=1.e-6, maxIter=5, verbose=False):
	"""
	========================================================================
	Derive effective frequency for spire bands. This is the frequency at
	which the monochromatic beam area is equal to the area as measured on
	Neptune.

	Inputs:
	  freq:        (array float) frequency vector corresponding to RSRF values [Hz]
	  rsrf:        (array float) relative spectral response (RSRF) corresponding to freq
	                 Note that this should *not* include the aperture efficiency
	  beamProfs:   (dataset) PhotRadialCorrBeam object from calibration tree
	  effFreqInit: (float) Initial estimate of effective frequency [Hz]
	  gamma:       (float) Exponent of powerlaw describing FWHM dependence
	                 on frequency
	  areaNep:     (float) Solid angle of beam measured on Neptune [sr]
	  alphaNep:    (float) Spectral index of Neptune frequency spectrum
	  array:       (string) spire array ('Psw'|'Pmw'|'Plw')
	  simpleBeam:  (boolean) set to use simple beam area model, which
                         scales as freq^(2.gamma)
	  freqFact:    (int) Factor by which to reduce size of freq for calculations.
	                 Only applicable if SimpleBeam=False
	                 OPTIONAL. Default=500.
	  initRange:   (float) Fractional intial range to use in calculations
	                 OPTIONAL. Default=0.01
	  reqPrec:     (float) Required relative precision for convergence.
	                 OPTIONAL. Default=1.e-6
          maxIter:     (int) Maximum interations to try
	                 OPTIONAL. Default=5
	  verbose:     (boolean) set to print more detailed info
                         OPTIONAL. Default=False

	Outputs:
	               (float) Effective frequency [Hz]

	Calculation:
	  Uses effFreqInit +/- initRange to calculate monochromatic Areas
	  Calculates effective beam area for Neptune spectrum (alphaNep)
 	  Compares with measured Neptune beam area (areaNep)
	  Adjusts effFreqInit and iterates until maximum iterations (maxIter) or
	  required precicion (reqPrec) is reached

	Dependencies:
	  spireMonoAreasSimple
	  spireMonoAreas
	  spireEffArea

	2014/01/07  C. North  initial version
	"""
		
	#calculate initial estimates of effective Frequency
	#parameterised be offset from original estimate
	relEff=Double1d([1.-initRange,1.+initRange])
	effFreqs=effFreqInit*relEff
	if verbose:
		print 'Calculating effective frequency for %s'%array
		print '  Initial %s Effective Frequencies: [%.2f : %.2f] GHz'%(array,effFreqs[0]/1.e9,effFreqs[1]/1.e9)
	#calculate effective beam area for initial estimates of effFreq
	if simpleBeam:
		beamMonoArea0=spireMonoAreasSimple(freq,effFreqs[0],areaNep,gamma)
		beamMonoArea1=spireMonoAreasSimple(freq,effFreqs[1],areaNep,gamma)
	else:
		beamMonoArea0=spireMonoAreas(freq,beamProfs,effFreqs[0],gamma,array,freqFact=freqFact)
		beamMonoArea1=spireMonoAreas(freq,beamProfs,effFreqs[1],gamma,array,freqFact=freqFact)
	beamMonoDiff=beamMonoArea1-beamMonoArea0
	effAreas=Double1d(2)
	effAreas[0]=spireEffArea(freq,rsrf, beamMonoArea0, BB=False, alpha=alphaNep)
	effAreas[1]=spireEffArea(freq,rsrf, beamMonoArea1, BB=False, alpha=alphaNep)

	iter=0
	done=False
	while ((done==False) and (iter <= maxIter)):
		iter=iter+1
		#difference from measured beam area
		diffAreas=effAreas-areaNep
		relAreas=diffAreas/areaNep
		#calculate new esitmate of rel
		grad=(diffAreas[1]-diffAreas[0])/(relEff[1]-relEff[0])
		relEffNew=relEff[1] - diffAreas[1]/grad

		#move values in arrays
		relEff[0]=relEff[1]
		effAreas[0]=effAreas[1]

		#calculate new effective beam area
		relEff[1]=relEffNew
		effFreqs=effFreqInit*relEff

		if simpleBeam:
			beamMonoNew=spireMonoAreasSimple(freq,effFreqs[1],areaNep,gamma)
		else:
			beamMonoNew=spireMonoAreas(freq,beamProfs,effFreqs[1],gamma,array,freqFact=freqFact)
		effAreas[1]=spireEffArea(freq,rsrf, beamMonoNew, BB=False, alpha=alphaNep)

		diffAreas=effAreas-areaNep
		relAreas=diffAreas/areaNep
		if verbose:
			print '    iter %d: %.4f [effFreq %.2f GHz], Area=%.2f, RelDiff=%.4g'%(iter,relEff[1],effFreqs[1]/1.e9,effAreas[1]/arcsec2Sr,relAreas[1])
		if (Math.abs(relAreas[1]) < reqPrec):
			done=True
	if ((iter > maxIter) and (done==False)):
		print "  Warning: maximum iterations [%d] exceeded without conversion [%g]"%(maxIter,reqPrec)

	if verbose:
		print '  Final %s effFreq: %.4f'%(array,effFreqs[1]/1.e9)
		print '  Resulting %s Neptune area: %.2f [rel Diff: %.3g]'%(array,effAreas[1]/arcsec2Sr,relAreas[1])
	return(effFreqs[1])

#-------------------------------------------------------------------------------
# Calculate the effective beam profile, area and beam map for SPIRE
def spireEffBeam(freq, transm, beamProfs, effFreq, gamma, array, beamRadMap,
  BB=False,temp=20.0,beta=1.8,alpha=-1.0,verbose=False):

	"""
	========================================================================
	Computes an effective beam profile for a given source spectrum
	***
	N.B. Computing the beam area this way integrates over frequency *then* radius,
	  while spireEffArea integrates over radius then frequency.
	  The method used here produces areas which area lower by ~0.1%
	***

	Inputs:
	  freq:       (array float) frequency vector [Hz] for which monochromatic
	                beams areas should be calculated
	  transm:     (array float) relative spectral response (RSRF) corresponding to freq
	                Note that this should *not* include the aperture efficiency
	  beamProfs:  (dataset) PhotRadialCorrBeam object from calibration tree
	  effFreq:    (float) effective frequency [Hz] of array
	  gamma:      (float) Exponent of powerlaw describing FWHM dependence
	                on frequency
	  array:      (string) spire array ('Psw'|'Pmw'|'Plw')
	  beamRadMap: (Simple Image) image containing radius at each point (in arcsec)
	  BB:         (boolean) spectral function to use for source spectrum:
	                'True': a modified black body
	                'False' a power-law with exponent alpha=-1
	                OPTIONAL. Default=False
	  temp:       (float) Dust/sky temperature (if BB=True)
	                OPTIONAL. Default=20.0
	  beta:       (float) Dust/sky spectral index (if BB=True)
	                OPTIONAL. Default=1.8
	  alpha:      (float) Exponent of power-law sky background model (if BB=False)
	                OPTIONAL. Default=-1

	Outputs:
	              (array float) radialised beam profile

	Calculation:
	  Calculates the source spectrum (either modifies black body or power law)
	  Loops over beam radius, and calculates scaled radii for full frequency range
	  For that radius, gets the values of profile at those scaled radii
	  Integrates profile values over frequency, weighted by RSRF and source spectrum

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

	2014/01/16  C. North  initial version

	"""

	#
	# Calculate sky background model
	#
	if BB == 1:
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha

	#integrate transm*fSky over frequency for nomalisation
	integrator=TrapezoidalIntegrator(min(freq),max(freq))
	denomInterp=CubicSplineInterpolator(freq,transm*fSky)
	denomInteg=integrator.integrate(denomInterp)

	#get beam radius list from calibration table
	beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
	#get core beam profile from calibration table
	beamCore=beamProfs.getCoreCorrectionTable().getColumn(array).data
	#create interpolation object
	beamCoreInt=CubicSplineInterpolator(beamRad,beamCore)

	#make array for new beam
	nRad=len(beamRad)
	maxRad=max(beamRad)
	effBeam=Float1d(beamRad)
	#loop over radius
	for r in range(nRad):
		#calculate the "scaled" radius for range of frequencies
		radFreq=beamRad[r]*(freq/effFreq)**-gamma
		#ensure it doesn't fo beyong maximum radius
		radFreq[radFreq.where(radFreq > maxRad)]=maxRad
		#compute value beam profile at each scaled radius
		beamCoreFreq=beamCoreInt(radFreq)
		#apply constant beam profile value where appropriate
		beamConstRad=beamProfs.getConstantCorrection(beamRad[r],array)
		isConst=beamCoreFreq.where(beamCoreFreq < beamConstRad)
		beamCoreFreq[isConst]=beamConstRad

		#integrate beamCoreFreq*transm*fSky over frequency
		numInterp=CubicSplineInterpolator(freq,beamCoreFreq*transm*fSky)
		numInteg = integrator.integrate(numInterp)

		#write value into table
		effBeam[r]=numInteg/denomInteg		

	#integrate over radius to get solid angle (in arcsec^2)
	
	effBeamAreaInterp=LinearInterpolator(beamRad,effBeam * 2. * Math.PI * beamRad)
	integrator=TrapezoidalIntegrator(0,maxRad)
	effBeamArea=integrator.integrate(effBeamAreaInterp)

	#create beam image by copying beamRad
	effBeamMap=beamRadMap.copy()
	#remove data in image
	effBeamMap.setUnit('Jy/beam')
	effBeamMap['image'].data[:,:]=0

	effBeamInterp=LinearInterpolator(beamRad,effBeam)
	nxMap=int(effBeamMap['image'].data[:,0].size)
	nyMap=int(effBeamMap['image'].data[0,:].size)
	if verbose:
		print 'making beam map'
	for x in range(nxMap):
		for y in range(nyMap):
			if beamRadMap['image'].data[x,y] <= maxRad-1:
				effBeamMap['image'].data[x,y]= \
				    effBeamInterp(beamRadMap['image'].data[x,y])

	effBeamDict={'area':effBeamArea,'profile':effBeam,'map':effBeamMap}
	return(effBeamDict)

#-------------------------------------------------------------------------------
# OLD METHOD: Calculate K-correction parameters for given spectrum & source type
def OLD_hpXcalKcorr(freq0, freq, transm, BB=True, temp=20.0, beta=1.8, alpha=-1.0, gamma=0.0):
	"""
	================================================================================
	Calculation of the K-correction factor from isophotal flux to a monochromatic 
	flux-density at a given reference frequency (data to be multiplied!)
	This routine is needed by hpXcalColorCorr.py
	
	Inputs:
	  freq0:    (float) waveband reference frequency [Hz] for which monochromatic
	            flux-density is given
	  freq:     (array float) frequency vector corresponding to RSRF values [Hz]
	  transm:   (array float) relative spectral response (RSRF) corresponding to freq
	  BB:       (boolean) spectral function to use for source spectrum:
	                'True': a modified black body
	                'False' a power-law with exponent alpha=-1
	  temp:     (float) Dust/sky temperature [optional; default is 20K; 
	                only for modified black body]
	  beta:     (float) Dust/sky spectral index [optional; default is 1.8; 
	                only for modified black body]
	  alpha:    (float) Exponent of power-law sky background model
	  gamma:    (float) Exponent of powerlaw describing FWHM dependence on frequency
	
	Outputs:
	 (list)     1st item: K-correction factor
	            2nd item: Sky emission at reference fequency (fSky0)
	
	Calculation:
	  Depending on the state of the input parameter BB, either the spectrum of a
	  Planck function multiplied by frequency to the power of beta, or a power-law
	  spectrum with spectral index alpha is calculated for all values in the vector 
	  frequ. In addition the same value is calculated at the discrete frequency 
	  freq0. Then the product of this value and the integral of the RSRF over all
	  frequencies, divided by the integral over all products of frequency and RSRF	
	  is calculated. Note that the integrals are coded as simple sums as the 
	  HIPE numeric integral doesn't allow too may discrete points and the RSRF
	  is sampled to quite some detail.
	
	
	2012/04/18  L. Conversi  initial version in file hp_xcal_Total_v1.6.py.txt
	2012/08/22  B. Schulz    added comments, reformatted, renamed kCorr to hpXcalKcorr
	2012/08/24  B. Schulz    fixed defaults for temp and beta, changed inputs to frequency,
	                         updated header and comments, renamed inputFreq to freq
	                         and inputFilt to transm added powerlaw sky background
	2012/08/27  B. Schulz    brush -up on header and comments
	2013/03/28  B. Schulz    removed implicit limitation to fixed frequency interval
	2013/04/05  B. Schulz    implemented proper tabulated integration function
	2013/06/18  L. Conversi  implemented gamma parameter in the denominator
	                         (previsuly applie directly to RSRFs, i.e. in the numerator too)
	
	================================================================================
	"""
	#
	# Calculate sky background model
	#
	# 1) As a modified Blackbody
	if BB == 1:
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
		fSky0 = 2*h * freq0**3 / c**2 / (EXP(h*freq0/k/temp) - 1.) * freq0**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha
		fSky0 = freq0**alpha
	#
	# Using the K-correction as defined in th OM (eq. 5.5, 5.15 & 5.16)
	kWave = fSky0 * IntTabulated(freq)(transm) / IntTabulated(freq)(transm * fSky * (freq/freq0)**(2*gamma))
	#
	# Return the result as a 2-element list of K-correction and flux at freq0
	return (kWave, fSky0)

#-------------------------------------------------------------------------------
# Calculate K-correction parameters for given spectrum & source type
def hpXcalKcorr(freq0, freq, transm, BB=True, temp=20.0, beta=1.8, alpha=-1.0,
  ext=False, monoArea=None):
	"""
	================================================================================
	Calculation of the K-correction factor from isophotal flux to a monochromatic 
	flux-density at a given reference frequency (data to be multiplied!)
	This routine is needed by hpXcalColorCorr.py
	
	Inputs:
	  freq0:     (float) waveband reference frequency [Hz] for which monochromatic
	               flux-density is given
	  freq:      (array float) frequency vector corresponding to RSRF values [Hz]
	  transm:    (array float) relative spectral response (RSRF) corresponding to freq
	  BB:        (boolean) spectral function to use for source spectrum:
	                'True': a modified black body
	                'False' a power-law with exponent alpha=-1
	                OPTIONAL. Default=False
	  temp:      (float) Dust/sky temperature [K] (only for modified black body)
			OPTIONAL. Deafult=20K; 
	                only for modified black body]
	  beta:      (float) Dust/sky spectral index (only for modified black body]
			OPTIONAL. Default=1.8
	  alpha:     (float) Exponent of power-law sky background model (only for
			power-law spectrum)
			OPTIONAL. Default=-1
	  ext:       (boolean) calculating for extended source
	                OPTIONAL. Default=False
	  monoArea:  (array float) Monochromatic Beam solid angle [Sr] corresponding
	                to freq.
	                OPTIONAL. Only required if ext=True

	Outputs:
	 (list)     [0]: K-correction factor
	            [1]: Sky emission at reference fequency (fSky0)
	
	Calculation:
	  Depending on the state of the input parameter BB, either the spectrum of a
	  Planck function multiplied by frequency to the power of beta, or a power-law
	  spectrum with spectral index alpha is calculated for all values in the vector 
	  frequ. In addition the same value is calculated at the discrete frequency 
	  freq0. Then the product of this value and the integral of the RSRF over all
	  frequencies, divided by the integral over all products of frequency and RSRF	
	  is calculated. Note that the integrals are coded as simple sums as the 
	  HIPE numeric integral doesn't allow too may discrete points and the RSRF
	  is sampled to quite some detail.

	  N.B. If ext=True, and monoBeam is in units sr, then this procedure outputs
	    K-correction factor in [Jy/sr per Jy/beam] and Sky emission in Jy/sr.
	    Units will change if a different input unit is used.
	
	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

	2012/04/18  L. Conversi  initial version in file hp_xcal_Total_v1.6.py.txt
	2012/08/22  B. Schulz    added comments, reformatted, renamed kCorr to hpXcalKcorr
	2012/08/24  B. Schulz    fixed defaults for temp and beta, changed inputs to frequency,
	                         updated header and comments, renamed inputFreq to freq
	                         and inputFilt to transm added powerlaw sky background
	2012/08/27  B. Schulz    brush -up on header and comments
	2013/03/28  B. Schulz    removed implicit limitation to fixed frequency interval
	2013/04/05  B. Schulz    implemented proper tabulated integration function
	2013/06/18  L. Conversi  implemented gamma parameter in the denominator
	                         (previsuly applie directly to RSRFs, i.e. in the numerator too)
	2013/12/03  C. North     corrected procedure include area where appropriate
			 	 NB: if ext=True , note output units
	
	================================================================================
	"""
	#
	# Calculate sky background model
	#
	# 1) As a modified Blackbody
	if BB == 1:
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
		fSky0 = 2*h * freq0**3 / c**2 / (EXP(h*freq0/k/temp) - 1.) * freq0**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha
		fSky0 = freq0**alpha
	#
	# Using the K-correction as defined in th OM (eq. 5.5, 5.15 & 5.16)

	if ext == True:
		area = monoArea
	else:
		#don't use area for point sources
		area = Float1d(len(freq))
		area[:] = 1.0

        # monoArea is the monochromatic beam solid angle at effFreq

	# integrate over frequency
	numInterp=LinearInterpolator(freq,transm)
	denomInterp=LinearInterpolator(freq,transm * fSky * area)
	minFreq=min(freq)
	maxFreq=max(freq)

	integrator=TrapezoidalIntegrator(minFreq,maxFreq)
	numInteg = integrator.integrate(numInterp)
	denomInteg = integrator.integrate(denomInterp)

	kWave = fSky0 * numInteg / denomInteg

	#
	# Return the result as a 2-element array of K-correction and flux at freq0
	return (Double1d([kWave, fSky0]))

def logMeta(name, meta, type='double'):
	"""
	Produce a metadata string for the log file
	
	Inputs:
	  name: (string) name of table in obj
	  meta: (Mtadata) metadata object containing key [name]
	  type: (string) type of metadate [double|string|int]
	        OPTIONAL. Default is double.

	Outputs:
	  	(string) String containing name, value, unit and description

	2014/01/20  C. North  initial version

	"""

	if type=='double':
		#check if there's a unit
		if herschel.share.util.StringUtil.asString(meta[name].unit)=='null':
			unit=''
		else:
			unit=meta[name].unit.getName()
		metaStr='%s: %.9g [%s] (Type=Double, Desc="%s")'%(name, meta[name].double, unit, meta[name].description)

	elif type=='int':
		#check if there's a unit
		if herschel.share.util.StringUtil.asString(meta[name].unit)=='null':
			unit=''
		else:
			unit=meta[name].unit.getName()
		metaStr='%s: %d [%s] (Type=Integer, Desc="%s")'%(name, meta[name].int, unit, meta[name].description)

	elif type=='string':
		metaStr='%s: "%s" (Type=String, Desc="%s")'%(name, meta[name].string, meta[name].description)

	return(metaStr)

def logTable(name, file, obj):
	"""
	Produce a table string for the log file
	
	Inputs:
	  name: (string) name of table in obj
	  file: (string) filename of ascii file containing table
	  obj:  (Product) calibration product containing table [name]

	Outputs:
	  	(string) String containing name, file and description

	2014/01/20  C. North  initial version

	"""
	tableStr='%s: %s (Desc="%s")'%(name, file, obj[name].getDescription())
	return(tableStr)

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           END OF FUNCTIONS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#===============================================================================
#=====                        CALCULATE BEAM PROFILES                      =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
if calcRadialCorrBeam:
	# copy old calibration product
	beamProfsOld = calibration.phot.getProduct("RadialCorrBeam")
	beamVersionOld = beamProfsOld.getVersion()
	# Set log file
	beamLogFile=os.path.join(dirRadialCorrBeam,'RadialCorrBeam_log_v%s.dat'%version)
	#create product and set basic metadata
	beamProfs = herschel.spire.ia.dataset.PhotRadialCorrBeam()
	beamProfs.setDescription('Spire radial beam profile product')
	beamProfs.meta["creator"].value = scriptVersionString
	beamProfs.meta["creationDate"].value = FineTime(java.util.Date())
	beamProfs.meta["startDate"].value = FineTime(startDate)
	beamProfs.meta["endDate"].value = FineTime(endDate)
	beamProfs.setVersion(version)
	beamProfs.setFormatVersion(formatVersion)
	beamProfsFileList='%s,%s'%(beamNewFileConstant,beamNewFileCore)
	beamProfs.meta['fileOrigin'] = StringParameter(beamProfsFileList,description="Origin of the data")

	beamProfs["constant"] = asciiTableReader(os.path.join(dataDir,beamNewFileConstant))
	beamProfs["core"] = asciiTableReader(os.path.join(dataDir,beamNewFileCore))
	#update labels
	beamProfs["core"].setDescription("Frequency-dependent core part of the radial beam profile")
	beamProfs["constant"].setDescription("Frequency-independent constant part of the radial beam profile")
	#generate file list (for metadata

	print "\nGenerating new RadialCorrBeam (version %s) [replaces version %s from %s]"% \
	(version,beamVersionOld,calVersion)

	# Re-compute normalised beam areas for new beams
	# get beam radius
	beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
	nRad=len(beamRad)

	# Create normArea table
	beamNormArea=TableDataset(description="Area as a function of radius, normalised by final value")
	# add radius column
	beamNormArea.addColumn("radius",Column(Float1d(beamRad)))

	print 'Calculating normalised beam area...'
	for band in spireBands:
		# add column to table
		beamNormArea.addColumn(band,Column(Float1d(nRad)))
		# get core beam
		beamComb=beamProfs.getCoreCorrectionTable().getColumn(band).data.copy()
		# get const beam
		beamConst=beamProfs.getConstantCorrectionTable().getColumn(band).data
		# work out where constant beam applies
		isConst = beamComb.where(beamComb < beamConst)
		# apply constant beam where applicable
		beamComb[isConst] = beamConst[isConst]
		# interpolate and integrate
		beamInterp=CubicSplineInterpolator(beamRad,beamComb*2.*Math.PI*beamRad)
		integrator=TrapezoidalIntegrator(0,max(beamRad))
		beamTotArea=integrator.integrate(beamInterp)
		for r in range(nRad):
			integrator=TrapezoidalIntegrator(0,beamRad[r])
			beamNormArea[band].data[r]=integrator.integrate(beamInterp)/beamTotArea
	# write table to beamProfs object
	beamProfs["normArea"]=beamNormArea

	# add gamma to metadata
	beamProfs.meta['gamma']=DoubleParameter(gamma,description='Exponent describing FWHM dependence on frequency')
	# add Neptune beam area to metaData (square arseconds)
	beamProfs.meta['beamNeptunePswArc']= DoubleParameter(spireAreaEffFreq['PSW']/arcsec2Sr, \
		unit=SolidAngle.SQUARE_SECONDS_ARC,description='PSW beam area as measured on Neptune')
	beamProfs.meta['beamNeptunePmwArc']= DoubleParameter(spireAreaEffFreq['PMW']/arcsec2Sr, \
		unit=SolidAngle.SQUARE_SECONDS_ARC,description='PMW beam area as measured on Neptune')
	beamProfs.meta['beamNeptunePlwArc']= DoubleParameter(spireAreaEffFreq['PLW']/arcsec2Sr, \
		unit=SolidAngle.SQUARE_SECONDS_ARC,description='PLW beam area as measured on Neptune')
	# add Neptune beam area to metaData (steradians)
	beamProfs.meta['beamNeptunePswSr']= DoubleParameter(spireAreaEffFreq['PSW'], \
		unit=SolidAngle.STERADIANS,description='PSW beam area as measured on Neptune')
	beamProfs.meta['beamNeptunePmwSr']= DoubleParameter(spireAreaEffFreq['PMW'], \
		unit=SolidAngle.STERADIANS,description='PMW beam area as measured on Neptune')
	beamProfs.meta['beamNeptunePlwSr']= DoubleParameter(spireAreaEffFreq['PLW'], \
		unit=SolidAngle.STERADIANS,description='PLW beam area as measured on Neptune')
	#add Neptune spectral index to metadata (steradians)
	beamProfs.meta['alphaNeptunePsw']= DoubleParameter(alphaNep['PSW'],\
		description='Neptune spectral used for PSW')
	beamProfs.meta['alphaNeptunePmw']= DoubleParameter(alphaNep['PMW'],\
		description='Neptune spectral used for PMW')
	beamProfs.meta['alphaNeptunePlw']= DoubleParameter(alphaNep['PLW'],\
		description='Neptune spectral used for PLW')

	#write tables to ascii files 
	asciiFileCore='RadialCorrBeam_core_v%s.csv'%version
	asciiFileConstant='RadialCorrBeam_constant_v%s.csv'%version
	asciiFileNormArea='RadialCorrBeam_normArea_v%s.csv'%version
	asciiTableWriter(table=beamProfs['core'],file=os.path.join(dirRadialCorrBeam,asciiFileCore))
	asciiTableWriter(table=beamProfs['constant'],file=os.path.join(dirRadialCorrBeam,asciiFileConstant))
	asciiTableWriter(table=beamProfs['normArea'],file=os.path.join(dirRadialCorrBeam,asciiFileNormArea))
	if verbose:
		print 'RadialCorrBeam written to ascii files'
	# set FITS filename
	beamProfsFits = java.io.File(r"%s//SCalPhotRadialCorrBeam_v%s.fits"%(dirRadialCorrBeam, version))
	beamProfs.meta['fileName'] = herschel.ia.dataset.StringParameter(value=beamProfsFits.name,\
	  description="Name of file when exported")

	# write log file
	beamLog=open(beamLogFile,'w')
	beamLog.write('Beam Radial Profiles Version %s\n'%version)
	beamLog.write('Creation Date: %s\n'%java.util.Date())
	beamLog.write('\nTables:\n')
	beamLog.write('  %s\n'%(logTable('core',asciiFileCore,beamProfs)))
	beamLog.write('  %s\n'%(logTable('constant',asciiFileConstant,beamProfs)))
	beamLog.write('  %s\n'%(logTable('normArea',asciiFileNormArea,beamProfs)))
	beamLog.write('\nMetadata:\n')
	beamLog.write('  %s\n'%(logMeta('gamma',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('beamNeptunePswArc',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('beamNeptunePmwArc',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('beamNeptunePlwArc',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('beamNeptunePswSr',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('beamNeptunePmwSr',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('beamNeptunePlwSr',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('alphaNeptunePsw',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('alphaNeptunePmw',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('alphaNeptunePlw',beamProfs.meta,'double')))
	beamLog.write('  %s\n'%(logMeta('fileOrigin',beamProfs.meta,'string')))
	beamLog.close()

else:
	beamProfs = calibration.phot.getProduct('RadialCorrBeam')
	beamVersion=beamProfs.getVersion()
	print "Using RadialCorrBeam version %s from calibration %s"%(beamVersion,calVersion)
	beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
	nRad=len(beamRad)
#-------------------------------------------------------------------------------
#===============================================================================
#=====                    CALCULATE EFFECTIVE FREQUENCIES                  =====
#===============================================================================
#-------------------------------------------------------------------------------
#load existing effective frequencies from calibration tree
spireEffFreq = {"PSW":beamProfsOld.meta["freqEffPsw"].double*1.e9,\
	"PMW":beamProfsOld.meta["freqEffPmw"].double*1.e9,\
	"PLW":beamProfsOld.meta["freqEffPlw"].double*1.e9}
#-------------------------------------------------------------------------------
if calcSpireEffFreq:
	# Find spire Effective frequencies
	print '\nGenerating new effective frequencies...'
	#uses existing numbers as initial guess
	for band in spireBands:
		spireEffFreq[band] = spireFindEffFreq(freq, spireFiltOnly[band],
		  beamProfs,spireEffFreq[band], gamma, spireAreaEffFreq[band],
		  alphaNep[band], band , simpleBeam=simpleBeam, verbose=verbose)

	if calcRadialCorrBeam:
		# update RadialBeamCorr metadata (in GHz)
		beamProfs.meta['freqEffPsw']=\
		  DoubleParameter(spireEffFreq['PSW']/1.e9,\
		  unit=Frequency.GIGAHERTZ,\
		  description='Effective frequency at which the measured PSW beam profile applies')
		beamProfs.meta['freqEffPmw']=\
		  DoubleParameter(spireEffFreq['PMW']/1.e9,\
		  unit=Frequency.GIGAHERTZ,\
		  description='Effective frequency at which the measured PMW beam profile applies')
		beamProfs.meta['freqEffPlw']=\
		  DoubleParameter(spireEffFreq['PLW']/1.e9,\
		  unit=Frequency.GIGAHERTZ,\
		  description='Effective frequency at which the measured PLW beam profile applies')
	
		#if not calcRadialCorrBeam:
		#	beamVersion = version
		#	print ' *** WARNING: RadialCorrBeam v%s uses v%s profiles from calibration %s'% \
		#		(beamVersion,beamVersionOld,calibration)
	
		# update log file
		beamLog=open(beamLogFile,'a')
		beamLog.write('  %s\n'%(logMeta('freqEffPsw',beamProfs.meta,'double')))
		beamLog.write('  %s\n'%(logMeta('freqEffPmw',beamProfs.meta,'double')))
		beamLog.write('  %s\n'%(logMeta('freqEffPlw',beamProfs.meta,'double')))
		beamLog.close()
else:
	beamProfs = calibration.phot.getProduct('RadialCorrBeam')
	beamVersion = beamProfs.getVersion()
	# Get effective frequencies from Calibration tree and convert to Hz
	print 'Using effective frequencies from RadialCorrBeam version %s from calibration %s'% \
	(beamVersion,calVersion)
	if calcRadialCorrBeam:
		print ' *** WARNING: Using new RadialCorrBeam (v%s) and old effective frequencies (from v%s)'%\
		(beamVersion,beamVersionOld)


#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CALCULATE MONOCHROMATIC BEAM AREAS                 =====
#===============================================================================
#-------------------------------------------------------------------------------

#calculate monochromatic beam areas using full or simple beam treatment
beamMonoArea={}
beamRefArea={}
print '\nCalculating monochromatic beam areas...'
for band in spireBands:
	beamMonoArea[band] = spireMonoAreas(freq, beamProfs, 
	  spireEffFreq[band], gamma, band)
	beamConst=beamProfs.getConstantCorrectionTable().getColumn(band).data
	beamRefArea[band] = spireMonoBeam(spireRefFreq[band],beamRad,\
	  beamProfs,beamConst,spireEffFreq[band],gamma,band)[0] * arcsec2Sr

#-------------------------------------------------------------------------------
#===============================================================================
#=====                      CALCULATE COLOR CORRECTIONS                    =====
#===============================================================================
#-------------------------------------------------------------------------------
# Description of parameters:
# K4P: * Scalar value per band
#      * Converts from broadband RSRF-weighted flux density (Jy/beam) to
#        monochromatic flux density (Jy/beam) at ref frequency a for point source
#        with alpha=-1 power law spectrum (nu*Fnu = const)
#      * This is applied in pipeline for the point source (PxWpsrc) L2 products
#      * These are stored in FluxConv metadata
#
# beamAreaPipSr: * Scalar value per band.
#                * Effective beam area (in sr) for a source with alpha=-1 spectrum
#                * These are stored in RadialCorrBeam and ColorCorrBeam metadata
#
# beamAreaPipArc: * Scalar value per band.
#                 * Effective beam area (in arcsec^2) for a source with alpha=-1 spectrum
#                 * These are stored in RadialCorrBeam and ColorCorrBeam metadata
#
# K4E_Tot: * Scalar value per band.
#          * Converts from broadband RSRF-weighted flux density (Jy/beam) to
#            monochromatic surface brightness (Jy/sr) at reference frequency for an
#            extended source with alpha=-1 spectrum (i.e. nu*Fnu = const)
#          * This is essentially what is applied in pipeline for the extended
#            source (PxWextd) L2 products. In actualitly, point source flux 
#            (with K4P applied) is multiplied by K4E * beamAreaPipSr / K4P.
#          * These values are not stored in the calibration tree
#
# K4E: * Scalar value per band.
#      * Defined as K4E_Tot * beamAreaPipSr.
#      * Converts from broadband RSRF-weighted flux density (Jy/beam) to
#        monochromatic flux density (Jy/beam) at reference frequency for an
#        extended source with alpha=-1 spectrum (i.e. nu*Fnu = const)
#      * These are stored in FluxConv metadata
#
# KPtoE: * Scalar value per band.
#        * Defined as K4E_Tot / K4P.
#        * Converts from monochromatic flux density (Jy/beam) of a point source
#          with alpha=-1 spectrum to the monochromatic surface brightness (Jy/sr)
#          of an extended source with alpha=-1 spectrum.
#        * These values are not stored in the calibration tree
#
# KcP: * Set of tables per band for various source spectra
#      * Colour correction to convert from monochromatic flux density (Jy/beam)
#        of a point source with alpha=-1 spectrum to the monochromatic flux density
#        (Jy/beam) of a point source with a given spectrum.
#      * These values are stored in the ColorCorrK_point tables
#
# KcE: * Set of tables per band for various source spectra
#      * Colour correction to convert from monochromatic flux density (Jy/beam)
#        of an extended source with alpha=-1 spectrum to the monochromatic flux 
#        density (Jy/beam) of an extended source with a given spectrum.
#      * These values are stored in the ColorCorrK_extended tables
#
# effBeamSr: * Set of tables per band for various source spectra
#            * Effective beam solid angle (in sr) for a source with a given spectrum
#            * These values are not stored in the calibration tree
#
# kBeam: * Set of tables per band for various source spectra
#        * Defined as beamAreaPipSr / effBeamSr
#        * Colour correction for beam solid angle
#        * These values are stores in the ColorCorrBeam tables
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Compute colour corrections for range of alpha and temp/beta
#set up dictionary for pipeline values
k4P = {}
k4E_Tot = {}
k4E = {}
k4P_Prev = {}
k4E_Prev = {}
k4P_Old = {}
k4E_Old = {}
kPtoE = {}
beamAreaPipSr  = {}
beamAreaPipArc  = {}
beamAreaPipSr_Old  = {}
beamAreaPipArc_Old  = {}
#create dictionaries for filenames
fileKCorrBeam = {}
fileKCorrPsrc = {}
fileKCorrExtd = {}

if calcColorCorrK:
	#-----------------------------------------------------------------------
	#=======================================================================
	#=====                 CALCULATE PIPELINE PARAMETERS               =====
	#=======================================================================
	#-----------------------------------------------------------------------
	#-----------------------------------------------------------------------
	# Calculate pipeline colour correction parameters
	print '\nGenerating new pipeline beam areas for ColorCorrBeam version %s'%version
	fluxConvOld=calibration.phot.getProduct('FluxConvList')[0]
	fluxConvOldVersion = fluxConvOld.getVersion()
	#set up log files
	fluxConvLogFile=os.path.join(dirFluxConv,'FluxConv_log_v%s.dat'%version)	

	#make dummy for fluxConv metadata
	fluxConv = fluxConvOld.copy()
	fluxConv.setDescription(fluxConvOld.getDescription())
	fluxConv.meta["creator"].value = scriptVersionString
	fluxConv.meta["creationDate"].value = FineTime(java.util.Date())
	fluxConv.meta["startDate"].value = FineTime(startDate)
	fluxConv.meta["endDate"].value = FineTime(endDate)
	fluxConv.setVersion(version)

	for band in spireBands:
		k4P_Old[band]=fluxConvOld[band].meta['k4P_%s'%band].double
		k4E_Old[band]=fluxConvOld[band].meta['k4E_%s'%band].double
		#pipeline beam areas
		beamAreaPipSr[band]=spireEffArea(freq, spireFiltOnly[band], \
		  beamMonoArea[band], BB=False, alpha=-1)
		beamAreaPipArc[band]=beamAreaPipSr[band]/arcsec2Sr
		#pipeline point source K4P conversion
		k4P[band] = hpXcalKcorr(spireRefFreq[band], freq, spireFilt[band], \
		  BB=False, ext=False)[0]
		#pipeline extended source K4E conversion
		k4E_Tot[band] = hpXcalKcorr(spireRefFreq[band], freq, spireFilt[band],
		  BB=False, ext=True, monoArea=beamMonoArea[band])[0]
		k4E[band] = k4E_Tot[band] * beamAreaPipSr[band]
		#calculate "K_PtoE" conversion
		kPtoE[band] = k4E_Tot[band] / k4P[band]
		#update metadata
		fluxConv.meta['K4P_%s'%band] = DoubleParameter(k4P[band],\
		  description='%s colour correction parameter for standard nu*F_nu=const. reference SED of a point source'%band)
		fluxConv.meta['K4E_%s'%band] = DoubleParameter(k4E[band],\
		  description='%s colour correction parameter for standard nu*F_nu=const. reference SED of an extended source'%band)

		# use previous method
		k4P_Prev[band] = OLD_hpXcalKcorr(spireRefFreq[band], freq, spireFilt[band], \
		  BB=False, alpha=-1.0)[0]
		k4E_Prev[band] = OLD_hpXcalKcorr(spireRefFreq[band], freq, spireFilt[band], \
		  BB=False, alpha=-1.0, gamma=gamma)[0]


	# set FITS filename
	fluxConvFits = java.io.File(r"%s//SCalPhotFluxConv_0_v%s.fits"%(dirFluxConv, version))
	fluxConv.meta['fileName'] = herschel.ia.dataset.StringParameter(value=fluxConvFits.name,\
	  description="Name of file when exported")

	# update FluxConv log file
	fluxConvVersion=fluxConv.getVersion()
	fluxConvLog=open(fluxConvLogFile,'w')
	fluxConvLog.write('Flux Conv Correction Version %s\n'%version)
	fluxConvLog.write('Description: %s\n'%fluxConv.getDescription())
	fluxConvLog.write('Creation Date: %s\n'%java.util.Date())
	fluxConvLog.write('\nTables: [UNCHANGED from v%s (cal %s)]\n'%(fluxConvOldVersion,calVersion))
	fluxConvLog.write('  ...\n')
	fluxConvLog.write('\nMetadata:\n')
	for band in spireBands:
		fluxConvLog.write('  %s\n'%(logMeta('K4P_%s'%band,fluxConv.meta,'double')))
		fluxConvLog.write('  %s\n'%(logMeta('K4E_%s'%band,fluxConv.meta,'double')))
	fluxConvLog.write('  %s\n'%(logMeta('fileName',fluxConv.meta,'string')))
	fluxConvLog.close()

	#-----------------------------------------------------------------------
	# write to FITS file
	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(fluxConvFits.toString(), fluxConv)

	#-----------------------------------------------------------------------
	#update RadialBeamCorr metadata
	if calcRadialCorrBeam:
		print 'Updating pipeline beam areas in RadialCorrBeam version %s'%version
		beamProfs.meta['beamPipelinePswSr']= DoubleParameter(beamAreaPipSr['PSW'], \
			unit=SolidAngle.STERADIANS,description='PSW beam area for spectral index alpha=-1 (as assumed in pipeline)')
		beamProfs.meta['beamPipelinePmwSr']= DoubleParameter(beamAreaPipSr['PMW'], \
			unit=SolidAngle.STERADIANS,description='PMW beam area for spectral index alpha=-1 (as assumed in pipeline)')
		beamProfs.meta['beamPipelinePlwSr']= DoubleParameter(beamAreaPipSr['PLW'], \
			unit=SolidAngle.STERADIANS,description='PLW beam area for spectral index alpha=-1 (as assumed in pipeline)')
		beamProfs.meta['beamPipelinePswArc']= DoubleParameter(beamAreaPipArc['PSW'], \
			unit=SolidAngle.SQUARE_SECONDS_ARC,description='PSW beam area for spectral index alpha=-1 (as assumed in pipeline)')
		beamProfs.meta['beamPipelinePmwArc']= DoubleParameter(beamAreaPipArc['PMW'], \
			unit=SolidAngle.SQUARE_SECONDS_ARC,description='PMW beam area for spectral index alpha=-1 (as assumed in pipeline)')
		beamProfs.meta['beamPipelinePlwArc']= DoubleParameter(beamAreaPipArc['PLW'], \
			unit=SolidAngle.SQUARE_SECONDS_ARC,description='PLW beam area for spectral index alpha=-1 (as assumed in pipeline)')
	
		#write to beamProfs log file
		beamLog=open(beamLogFile,'a')
		beamLog.write('  %s\n'%(logMeta('beamPipelinePswSr',beamProfs.meta,'double')))
		beamLog.write('  %s\n'%(logMeta('beamPipelinePmwSr',beamProfs.meta,'double')))
		beamLog.write('  %s\n'%(logMeta('beamPipelinePlwSr',beamProfs.meta,'double')))
		beamLog.write('  %s\n'%(logMeta('beamPipelinePswArc',beamProfs.meta,'double')))
		beamLog.write('  %s\n'%(logMeta('beamPipelinePmwArc',beamProfs.meta,'double')))
		beamLog.write('  %s\n'%(logMeta('beamPipelinePlwArc',beamProfs.meta,'double')))
		beamLog.close()
	
	#-----------------------------------------------------------------------
	#=======================================================================
	#=====               CALCULATE BEAM COLOR CORRECTIONS              =====
	#=======================================================================
	#-----------------------------------------------------------------------
	#-----------------------------------------------------------------------

	#get old files from calibration
	kCorrBeamOld=calibration.phot.getProduct('ColorCorrBeam')
	kCorrBeamOldVersion = kCorrBeamOld.getVersion()
	#get old piipeline beam areas
	beamAreaPipSr_Old["PSW"]=kCorrBeamOld.meta['beamPswSr'].double
	beamAreaPipSr_Old["PMW"]=kCorrBeamOld.meta['beamPmwSr'].double
	beamAreaPipSr_Old["PLW"]=kCorrBeamOld.meta['beamPlwSr'].double
	beamAreaPipArc_Old["PSW"]=kCorrBeamOld.meta['beamPswArc'].double
	beamAreaPipArc_Old["PMW"]=kCorrBeamOld.meta['beamPmwArc'].double
	beamAreaPipArc_Old["PLW"]=kCorrBeamOld.meta['beamPlwArc'].double
	effBeamSrOld=kCorrBeamOld.copy()
	#compute old beam area
	for name in kCorrBeamOld.getSets():
		for band in spireBands:
			effBeamSrOld[name][band].data = \
			  beamAreaPipSr_Old[band] / kCorrBeamOld[name][band].data 

	#set up log files
	kCorrBeamLogFile=os.path.join(dirColorCorrBeam,'ColorCorrBeam_log_v%s.dat'%version)
	#create product and set basic metadata
	kCorrBeam=herschel.spire.ia.dataset.PhotColorCorrBeam()
	kCorrBeam.setDescription('Spire beam corrections with spectral index and temperature product')
	kCorrBeam.meta["creator"].value = scriptVersionString
	kCorrBeam.meta["creationDate"].value = FineTime(java.util.Date())
	kCorrBeam.meta["startDate"].value = FineTime(startDate)
	kCorrBeam.meta["endDate"].value = FineTime(endDate)
	kCorrBeam.setVersion(version)
	kCorrBeam.setFormatVersion(formatVersion)

	# make copy of kBeam to contain raw effective beam areas
	effBeamSr=herschel.spire.ia.dataset.PhotColorCorrBeam()
	effBeamSr.setDescription('Spire beam solid angle with spectral index and temperature product')

	#-----------------------------------------------------------------------
	# Compute beam corrections for range of alpha and temp/beta
	# Create tables for alpha arrays for KBeam
	effBeamSr['alpha']=TableDataset(description='Beam Solid Angle (Spectral Index)')
	kCorrBeam['alpha']=TableDataset(description='Beam Colour Correction (Spectral Index)')

	#add alpha column
	effBeamSr['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	kCorrBeam['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	
	print 'Calculating beam correction parameters over alpha...'
	for band in spireBands:
		# add columns to tables
		effBeamSr['alpha'].addColumn(band,Column(Float1d(len(alphaK)),unit=SolidAngle.STERADIANS,description=''))
		kCorrBeam['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
	
		for a in range(len(alphaK)):
			#effective area
			effBeamSr['alpha'][band].data[a]=spireEffArea(freq, spireFiltOnly[band],
			  beamMonoArea[band], BB=False, alpha=alphaK[a])
			#beam correction factor
			kCorrBeam['alpha'][band].data[a] = \
			  beamAreaPipSr[band] / effBeamSr['alpha'][band].data[a]

	#-----------------------------------------------------------------------
	print 'Calculating beam colour correction parameters over beta & temp...'
	for b in range(len(betaK)):
		#create format version of beta text (beta=x.yz -> beta_x_yz)
		betaTxt='beta_%d_%d%d'%(int(betaK[b]),int(10*(betaK[b]%1)),int(10*((10*betaK[b])%1)))
		if verbose:
			print '  %s: beta=%f'%(betaTxt,betaK[b])
		#create tables
		effBeamSr[betaTxt]=TableDataset(description='Beam Colour Correction (Modified Black Body, beta=%.2f)'%betaK[b])
		kCorrBeam[betaTxt]=TableDataset(description='Beam Colour Correction (Modified Black Body, beta=%.2f)'%betaK[b])
		# add temp column
		effBeamSr[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		kCorrBeam[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		for band in spireBands:
			effBeamSr[betaTxt].addColumn(band,Column(Float1d(len(tempK)),unit=SolidAngle.STERADIANS,description=''))
			kCorrBeam[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			for t in range(len(tempK)):
				#effective area
				effBeamSr[betaTxt][band].data[t]=spireEffArea(freq, spireFiltOnly[band],
				  beamMonoArea[band], BB=True, beta=betaK[b], temp=tempK[t])
				#beam correction factor
				kCorrBeam[betaTxt][band].data[t] = \
				  beamAreaPipSr[band] / effBeamSr[betaTxt][band].data[t]

	#write to ascii table files
	kCorrBeamFileList=''
	for name in kCorrBeam.getSets():
		#set filename
		fileKCorrBeam[name]='ColorCorrBeam_%s_v%s.csv'%(name,version)
		#write to ascii table
		asciiTableWriter(table=kCorrBeam[name],file=os.path.join(dirColorCorrBeam,fileKCorrBeam[name]))
		#update filelist (for metadata)
		if kCorrBeamFileList=='':
			kCorrBeamFileList='%s'%fileKCorrBeam[name]
		else:
			kCorrBeamFileList='%s,%s'%(kCorrBeamFileList,fileKCorrBeam[name])
	kCorrBeam.meta["fileOrigin"] = herschel.ia.dataset.StringParameter(value="%s"%kCorrBeamFileList,\
	  description="Origin of the data")

	# set FITS filename
	kCorrBeamFits = java.io.File(r"%s//SCalPhotColorCorrBeam_v%s.fits"%(dirColorCorrBeam, version))
	kCorrBeam.meta['fileName'] = herschel.ia.dataset.StringParameter(value=kCorrBeamFits.name,\
	  description="Name of file when exported")

	#update ColorCorrBeam metadata
	kCorrBeam.meta['beamPswSr']=DoubleParameter(beamAreaPipSr['PSW'],\
	 unit=SolidAngle.STERADIANS,description='PSW beam area for spectral index alpha=-1 (as assumed in the pipeline')
	kCorrBeam.meta['beamPmwSr']=DoubleParameter(beamAreaPipSr['PMW'],\
	 unit=SolidAngle.STERADIANS,description='PMW beam area for spectral index alpha=-1 (as assumed in the pipeline')
	kCorrBeam.meta['beamPlwSr']=DoubleParameter(beamAreaPipSr['PLW'],\
	 unit=SolidAngle.STERADIANS,description='PLW beam area for spectral index alpha=-1 (as assumed in the pipeline')
	kCorrBeam.meta['beamPswArc']=DoubleParameter(beamAreaPipArc['PSW'],\
	 unit=SolidAngle.SQUARE_SECONDS_ARC,description='PSW beam area for spectral index alpha=-1 (as assumed in the pipeline')
	kCorrBeam.meta['beamPmwArc']=DoubleParameter(beamAreaPipArc['PMW'],\
	 unit=SolidAngle.SQUARE_SECONDS_ARC,description='PMW beam area for spectral index alpha=-1 (as assumed in the pipeline')
	kCorrBeam.meta['beamPlwArc']=DoubleParameter(beamAreaPipArc['PLW'],\
	 unit=SolidAngle.SQUARE_SECONDS_ARC,description='PLW beam area for spectral index alpha=-1 (as assumed in the pipeline')

	# write log file
	kCorrBeamLog=open(kCorrBeamLogFile,'w')
	kCorrBeamLog.write('Beam Colour Correction Version %s\n'%version)
	kCorrBeamLog.write('Description: %s\n'%kCorrBeam.getDescription())
	kCorrBeamLog.write('Creation Date: %s\n'%java.util.Date())
	kCorrBeamLog.write('\nTables:\n')
	for name in kCorrBeam.getSets():
		kCorrBeamLog.write('  %s\n'%(logTable(name,fileKCorrBeam[name],kCorrBeam)))
	kCorrBeamLog.write('\nMetadata:\n')
	kCorrBeamLog.write('  %s\n'%(logMeta('beamPswSr',kCorrBeam.meta,'double')))
	kCorrBeamLog.write('  %s\n'%(logMeta('beamPmwSr',kCorrBeam.meta,'double')))
	kCorrBeamLog.write('  %s\n'%(logMeta('beamPlwSr',kCorrBeam.meta,'double')))
	kCorrBeamLog.write('  %s\n'%(logMeta('beamPswArc',kCorrBeam.meta,'double')))
	kCorrBeamLog.write('  %s\n'%(logMeta('beamPmwArc',kCorrBeam.meta,'double')))
	kCorrBeamLog.write('  %s\n'%(logMeta('beamPlwArc',kCorrBeam.meta,'double')))
	kCorrBeamLog.write('  %s\n'%(logMeta('fileOrigin',kCorrBeam.meta,'string')))
	kCorrBeamLog.write('  %s\n'%(logMeta('fileName',kCorrBeam.meta,'string')))
	kCorrBeamLog.close()

	#-----------------------------------------------------------------------
	# write to FITS file
	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(kCorrBeamFits.toString(), kCorrBeam)

	#-----------------------------------------------------------------------
	# plot comparisons
	if plot:
		cols={'PSW':java.awt.Color.BLUE,\
		  'PMW':java.awt.Color.GREEN,\
		  'PLW':java.awt.Color.RED}
		#plot beam correction against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kCorrBeam['alpha']['alpha'].data,\
			  kCorrBeam['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New'))
			p.addLayer(LayerXY( \
			  kCorrBeamOld['alpha']['alpha'].data,\
			  kCorrBeamOld['alpha'][band].data,\
			  color=cols[band],name=band+' Old'))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Beam correction"
		p.setTitleText("ColorCorrBeam value")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrBeamOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrBeam,'SCalPhotColorCorrBeam_v%s_v%s.png'%(kCorrBeamOldVersion,version)))

		#plot beam correction change against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kCorrBeam['alpha']['alpha'].data,\
			  kCorrBeam['alpha'][band].data/kCorrBeamOld['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New/Old'))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Beam correction w.r.t. v%s"%(kCorrBeamOldVersion)
		p.setTitleText("ColorCorrBeam relative change")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrBeamOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrBeam,'SCalPhotColorCorrBeam_RelDiff_v%s_v%s.png'%(kCorrBeamOldVersion,version)))

		#plot beam area against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  effBeamSr['alpha']['alpha'].data,\
			  effBeamSr['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New'))
			p.addLayer(LayerXY( \
			  effBeamSrOld['alpha']['alpha'].data,\
			  effBeamSrOld['alpha'][band].data,\
			  color=cols[band],name=band+' Old'))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Beam area (sr)"
		p.setTitleText("Effective Beam Area value")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrBeamOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrBeam,'effArea_v%s_v%s.png'%(kCorrBeamOldVersion,version)))

		#plot beam area change against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  effBeamSr['alpha']['alpha'].data,\
			  effBeamSr['alpha'][band].data/effBeamSrOld['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New/Old'))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Effective Beam area w.r.t. v%s"%(kCorrBeamOldVersion)
		p.setTitleText("Effective Beam Area relative change")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrBeamOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrBeam,'effArea_RelDiff_v%s_v%s.png'%(kCorrBeamOldVersion,version)))

		#plot beam area error against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  effBeamSr['alpha']['alpha'].data,\
			  effBeamSr['alpha'][band].data/beamRefArea[band],\
			  color=cols[band],stroke=2.0,name=band+' Eff/Ref'))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Effective Beam area w.r.t Ref. Freq"
		p.setTitleText("Beam Area error")
		p.setSubtitleText("New method vs. Old method")
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrBeam,'effArea_Error_v%s.fits'%version))
	#-----------------------------------------------------------------------
	#=======================================================================
	#=====           CALCULATE POINT SOURCE COLOR CORRECTIONS          =====
	#=======================================================================
	#-----------------------------------------------------------------------
	#-----------------------------------------------------------------------
	# Compute point source colour correction values
	#get old files from calibration
	kCorrPsrcOld=calibration.phot.getProduct('ColorCorrKList')[1]
	kConvPsrcOld=kCorrPsrcOld.copy()
	for name in kCorrPsrcOld.getSets():
		for band in spireBands:
			kConvPsrcOld[name][band].data = \
			 kCorrPsrcOld[name][band].data * k4P_Old[band]
	kCorrPsrcOldVersion = kCorrPsrcOld.getVersion()
	#set up log files
	kCorrPsrcLogFile=os.path.join(dirColorCorrK,'ColorCorrK_point_log_v%s.dat'%version)

	# create product and update basic metadata
	kCorrPsrc = herschel.spire.ia.dataset.PhotColorCorrK()
	kCorrPsrc.setDescription('Spire color corrections with spectral index and temperature product')
	kCorrPsrc.meta["creator"].value = scriptVersionString
	kCorrPsrc.meta["creationDate"].value = FineTime(java.util.Date())
	kCorrPsrc.meta["startDate"].value = FineTime(startDate)
	kCorrPsrc.meta["endDate"].value = FineTime(endDate)
	kCorrPsrc.setVersion(version)
	kCorrPsrc.setFormatVersion(formatVersion)
	
	#make dummy product for point source conversion
	kConvPsrc=herschel.spire.ia.dataset.PhotColorCorrK()
	kConvPsrc.setDescription('Spire color conversions with spectral index and temperature product')

	#-----------------------------------------------------------------------
	# Compute point source colour corrections for range of alpha and temp/beta
	beamAreaPipSr_Old["PSW"]=kCorrBeamOld.meta['beamPswSr'].double
	beamAreaPipSr_Old["PMW"]=kCorrBeamOld.meta['beamPmwSr'].double
	beamAreaPipSr_Old["PLW"]=kCorrBeamOld.meta['beamPlwSr'].double
	beamAreaPipArc_Old["PSW"]=kCorrBeamOld.meta['beamPswArc'].double
	beamAreaPipArc_Old["PMW"]=kCorrBeamOld.meta['beamPmwArc'].double
	beamAreaPipArc_Old["PLW"]=kCorrBeamOld.meta['beamPlwArc'].double

	print '\nCalculating point source colour correction parameters over alpha...'
	# Create tables for alpha arrays for K (point source)
	kCorrPsrc['alpha']=TableDataset(description='Point Source Colour Correction (Spectral Index)')
	kConvPsrc['alpha']=TableDataset(description='Point Source Conversion (Spectral Index)')
	# add alpha column
	kCorrPsrc['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	kConvPsrc['alpha'].addColumn('alpha',Column(Float1d(alphaK)))

	for band in spireBands:
		# add band columns to tables
		kCorrPsrc['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		kConvPsrc['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		for a in range(len(alphaK)):
			#point source conversion for current alpha
			kConvPsrc['alpha'][band].data[a]=hpXcalKcorr(spireRefFreq[band],\
			   freq, spireFilt[band], BB=False, alpha=alphaK[a])[0]
			#point source colour correction for current alpha
			kCorrPsrc['alpha'][band].data[a] = \
			  kConvPsrc['alpha'][band].data[a]/k4P[band]

	#-----------------------------------------------------------------------
	#use previous method
	kCorrPsrcPrev=herschel.spire.ia.dataset.PhotColorCorrK()
	kCorrPsrcPrev['alpha']=TableDataset()
	kConvPsrcPrev=herschel.spire.ia.dataset.PhotColorCorrK()
	kConvPsrcPrev['alpha']=TableDataset()
	kCorrPsrcPrev['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	kConvPsrcPrev['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	for band in spireBands:
		kCorrPsrcPrev['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		kConvPsrcPrev['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		for a in range(len(alphaK)):
			#point source conversion for current alpha (previous method)
			kConvPsrcPrev['alpha'][band].data[a]=\
			  OLD_hpXcalKcorr(spireRefFreq[band],\
			  freq, spireFilt[band], BB=False, alpha=alphaK[a])[0]
			#point source colour correction for current alpha (previous method)
			kCorrPsrcPrev['alpha'][band].data[a] = \
			  kConvPsrcPrev['alpha'][band].data[a] / k4P_Prev[band]
	
	#-----------------------------------------------------------------------
	print 'Calculating point source colour correction parameters over beta & temp...'
	for b in range(len(betaK)):
		#create format version of beta text (beta=x.yz -> beta_x_yz)
		betaTxt='beta_%d_%d%d'%(int(betaK[b]),int(10*(betaK[b]%1)),int(10*((10*betaK[b])%1)))
		if verbose:
			print '  %s: beta=%f'%(betaTxt,betaK[b])
		# Create tables for beta arrays for K (point source)
		kCorrPsrc[betaTxt]=TableDataset(description='Point Source Colour Correction (Modified Black Body, beta=%.2f)'%betaK[b])
		kConvPsrc[betaTxt]=TableDataset(description='Point Source Conversion (Modified Black Body, beta=%.2f)'%betaK[b])
		# add temp column
		kCorrPsrc[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		kConvPsrc[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		for band in spireBands:
			kCorrPsrc[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			kConvPsrc[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			for t in range(len(tempK)):
				#point source conversion for current beta,temp
				kConvPsrc[betaTxt][band].data[t]=hpXcalKcorr(spireRefFreq[band],
				  freq, spireFilt[band], BB=True, beta=betaK[b], temp=tempK[t])[0]
				#point source colour correction for current beta,temp
				kCorrPsrc[betaTxt][band].data[t] = \
				  kConvPsrc[betaTxt][band].data[t]/k4P[band]

		#---------------------------------------------------------------
		#use previous method
		kCorrPsrcPrev[betaTxt]=TableDataset()
		kConvPsrcPrev[betaTxt]=TableDataset()
		kCorrPsrcPrev[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		kConvPsrcPrev[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		for band in spireBands:
			kCorrPsrcPrev[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			kConvPsrcPrev[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			for t in range(len(tempK)):
				#point source conversion for current beta,temp (previous method)
				kConvPsrcPrev[betaTxt][band].data[t]=\
				  hpXcalKcorr(spireRefFreq[band],
				  freq, spireFilt[band], BB=True, \
				  beta=betaK[b], temp=tempK[t])[0]
				#point source colour correction for current beta,temp (previous method)
				kCorrPsrcPrev[betaTxt][band].data[t] = \
				  kConvPsrcPrev[betaTxt][band].data[t] / k4P_Prev[band]

	#write to ascii table files
	kCorrPsrcFileList=''
	for name in kCorrPsrc.getSets():
		#set filename
		fileKCorrPsrc[name]='ColorCorrK_point_%s_v%s.csv'%(name,version)
		#write to ascii table
		asciiTableWriter(table=kCorrPsrc[name],file=os.path.join(dirColorCorrK,fileKCorrPsrc[name]))
		#update filelist (for metadata)
		if kCorrPsrcFileList=='':
			kCorrPsrcFileList='%s'%fileKCorrPsrc[name]
		else:
			kCorrPsrcFileList='%s,%s'%(kCorrPsrcFileList,fileKCorrPsrc[name])
	kCorrPsrc.meta["fileOrigin"] = herschel.ia.dataset.StringParameter(value="%s"%kCorrPsrcFileList,\
	  description="Origin of the data")

	# set FITS filename
	kConvPsrcFits = java.io.File(r"%s//SCalPhotConvK_point_v%s.fits"%(dirColorCorrK, version))
	kCorrPsrcFits = java.io.File(r"%s//SCalPhotColorCorrK_point_v%s.fits"%(dirColorCorrK, version))
	kCorrPsrc.meta['fileName'] = herschel.ia.dataset.StringParameter(value=kCorrPsrcFits.name,\
	  description="Name of file when exported")

	# write log file
	kCorrPsrcLog=open(kCorrPsrcLogFile,'w')
	kCorrPsrcLog.write('Point Source Colour Correction Version %s\n'%version)
	kCorrPsrcLog.write('Description: %s\n'%kCorrPsrc.getDescription())
	kCorrPsrcLog.write('Creation Date: %s\n'%java.util.Date())
	kCorrPsrcLog.write('\nTables:\n')
	for name in kCorrPsrc.getSets():
		kCorrPsrcLog.write('  %s\n'%(logTable(name,fileKCorrPsrc[name],kCorrPsrc)))
	kCorrPsrcLog.write('\nMetadata:\n')
	kCorrPsrcLog.write('  %s\n'%(logMeta('fileOrigin',kCorrPsrc.meta,'string')))
	kCorrPsrcLog.write('  %s\n'%(logMeta('fileName',kCorrPsrc.meta,'string')))
	kCorrPsrcLog.close()

	#-----------------------------------------------------------------------
	# write to FITS file
	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(kCorrPsrcFits.toString(), kCorrPsrc)

	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(kConvPsrcFits.toString(), kConvPsrc)

	#-----------------------------------------------------------------------
	# plot comparisons
	if plot:
		cols={'PSW':java.awt.Color.BLUE,\
		  'PMW':java.awt.Color.GREEN,\
		  'PLW':java.awt.Color.RED}
		#plot colour correction against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kCorrPsrc['alpha']['alpha'].data,\
			  kCorrPsrc['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New'))
			p.addLayer(LayerXY( \
			  kCorrPsrcOld['alpha']['alpha'].data,\
			  kCorrPsrcOld['alpha'][band].data,\
			  color=cols[band],name=band+' Old'))
			p.addLayer(LayerXY( \
			  kCorrPsrcPrev['alpha']['alpha'].data,\
			  kCorrPsrcPrev['alpha'][band].data,\
			  color=cols[band],name=band+' Prev',line=Style.DASHED))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Colour correction"
		p.setTitleText("ColorCorrK_point value")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrPsrcOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrK,'SCalPhotColorCorrK_point_v%s_v%s.png'%(kCorrPsrcOldVersion,version)))

		#plot colour correction change against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kCorrPsrc['alpha']['alpha'].data,\
			  kCorrPsrc['alpha'][band].data/kCorrPsrcOld['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New/Old'))
			p.addLayer(LayerXY( \
			  kCorrPsrcPrev['alpha']['alpha'].data,\
			  kCorrPsrcPrev['alpha'][band].data/kCorrPsrcOld['alpha'][band].data,\
			  color=cols[band],name=band+' Prev/Old',line=Style.DASHED))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Colour correction w.r.t. v%s"%(kCorrPsrcOldVersion)
		p.setTitleText("ColorCorrK_point relative change")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrPsrcOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrK,'SCalPhotColorCorrK_point_RelDiff_v%s_v%s.png'%(kCorrPsrcOldVersion,version)))
		#plot conversion against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kConvPsrc['alpha']['alpha'].data,\
			  kConvPsrc['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New'))
			p.addLayer(LayerXY( \
			  kConvPsrcOld['alpha']['alpha'].data,\
			  kConvPsrcOld['alpha'][band].data,\
			  color=cols[band],name=band+' Old'))
			p.addLayer(LayerXY( \
			  kConvPsrcPrev['alpha']['alpha'].data,\
			  kConvPsrcPrev['alpha'][band].data,\
			  color=cols[band],name=band+' Prev',line=Style.DASHED))
		p.xaxis.titleText = "Spectral Index (alpha)"
		p.yaxis.titleText = "Point Source Conversion"
		p.setTitleText("K4P x ColorCorrK_point value")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrPsrcOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrK,'KConv_point_v%s_v%s.png'%(kCorrPsrcOldVersion,version)))
		#plot conversion change against alpha
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kConvPsrc['alpha']['alpha'].data,\
			  kConvPsrc['alpha'][band].data/kConvPsrcOld['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+ 'New/Old'))
			p.addLayer(LayerXY( \
			  kConvPsrcPrev['alpha']['alpha'].data,\
			  kConvPsrcPrev['alpha'][band].data/kConvPsrcOld['alpha'][band].data,\
			  color=cols[band],name=band+' Prev/Old',line=Style.DASHED))
		p.xaxis.titleText = "Spectral Index (alpha)"
		p.yaxis.titleText = "Point Source Conversion w.r.t. v%s"%kCorrPsrcOldVersion
		p.setTitleText("K4P x ColorCorrK_point relative change")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrPsrcOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrK,'KConv_point_RelDiff_v%s_v%s.png'%(kCorrPsrcOldVersion,version)))

	#-----------------------------------------------------------------------
	#=======================================================================
	#=====         CALCULATE EXTENDED SOURCE COLOR CORRECTIONS         =====
	#=======================================================================
	#-----------------------------------------------------------------------
	#-----------------------------------------------------------------------
	# Compute extended source colour correction values
	#-----------------------------------------------------------------------
	#get old files from calibration
	kCorrExtdOld=calibration.phot.getProduct('ColorCorrKList')[0]
	kCorrExtdOldVersion = kCorrExtdOld.getVersion()
	kConvExtdOld=kCorrExtdOld.copy()
	for name in kCorrExtdOld.getSets():
		for band in spireBands:
			kConvExtdOld[name][band].data = \
			  kCorrExtdOld[name][band].data * k4E_Old[band]
	#set up log files
	kCorrExtdLogFile=os.path.join(dirColorCorrK,'ColorCorrK_extended_log_v%s.dat'%version)

	# update basic metadata
	kCorrExtd = herschel.spire.ia.dataset.PhotColorCorrK() 
	kCorrExtd.setDescription('Spire color corrections with spectral index and temperature product')
	kCorrExtd.meta["creator"].value = scriptVersionString
	kCorrExtd.meta["creationDate"].value = FineTime(java.util.Date())
	kCorrExtd.meta["startDate"].value = FineTime(startDate)
	kCorrExtd.meta["endDate"].value = FineTime(endDate)
	kCorrExtd.setVersion(version)
	kCorrExtd.setFormatVersion(formatVersion)

	#create dummy for extended conversion parameter
	kConvExtd = herschel.spire.ia.dataset.PhotColorCorrK()
	kConvExtd.setDescription('Spire conversions with spectral index and temperature product')
	#-----------------------------------------------------------------------
	# Compute extended source colour corrections for range of alpha and temp/beta
	# Create tables for alpha arrays for K (extended source)
	kCorrExtd['alpha']=TableDataset(description='Extended Source Colour Correction (Spectral Index)')
	kConvExtd['alpha']=TableDataset(description='Extended Source Conversion (Spectral Index)')
	# add alpha column
	kCorrExtd['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	kConvExtd['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	#-----------------------------------------------------------------------
	print '\nCalculating extended source colour correction parameters over alpha...'
	for band in spireBands:

		# add band columns to tables
		kCorrExtd['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		kConvExtd['alpha'].addColumn(band,Column(Float1d(len(alphaK))))

		for a in range(len(alphaK)):
			#total extended source conversion for current alpha
			k4EaTot_x=hpXcalKcorr(spireRefFreq[band], freq,
			 spireFilt[band], BB=False, alpha=alphaK[a],
			 ext=True, monoArea=beamMonoArea[band])[0]
			kConvExtd['alpha'][band].data[a] = \
			  k4EaTot_x * effBeamSr['alpha'][band].data[a]
			#extended source colour correction for current alpha
			kCorrExtd['alpha'][band].data[a]=k4EaTot_x / k4E_Tot[band]

	#-----------------------------------------------------------------------
	#use previous method
	kCorrExtdPrev=herschel.spire.ia.dataset.PhotColorCorrK()
	kCorrExtdPrev['alpha']=TableDataset()
	kConvExtdPrev=herschel.spire.ia.dataset.PhotColorCorrK()
	kConvExtdPrev['alpha']=TableDataset()
	kCorrExtdPrev['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	kConvExtdPrev['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	for band in spireBands:
		kCorrExtdPrev['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		kConvExtdPrev['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		for a in range(len(alphaK)):
			#point source conversion for current alpha (previous method)
			kConvExtdPrev['alpha'][band].data[a]=\
			  OLD_hpXcalKcorr(spireRefFreq[band],\
			  freq, spireFilt[band], BB=False,\
			  alpha=alphaK[a],gamma=gamma)[0]
			#point source colour correction for current alpha (previous method)
			kCorrExtdPrev['alpha'][band].data[a] = \
			  kConvExtdPrev['alpha'][band].data[a] / k4E_Prev[band]

	#-----------------------------------------------------------------------
	print 'Calculating extended source colour correction parameters over beta & temp...'
	for b in range(len(betaK)):
		#create format version of beta text (beta=x.yz -> beta_x_yz)
		betaTxt='beta_%d_%d%d'%(int(betaK[b]),int(10*(betaK[b]%1)),int(10*((10*betaK[b])%1)))
		if verbose:
			print '  %s: beta=%f'%(betaTxt,betaK[b])
		# Create tables for beta for K (extended source)
		kCorrExtd[betaTxt]=TableDataset(description='Extended Source Colour Correction (Modified Black Body, beta=%.2f)'%betaK[b])
		kConvExtd[betaTxt]=TableDataset(description='Extended Source Conversion (Modified Black Body, beta=%.2f)'%betaK[b])
		# add temp column
		kCorrExtd[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		kConvExtd[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		for band in spireBands:
			kCorrExtd[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			kConvExtd[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			for t in range(len(tempK)):
				#extended source conversion for current beta,temp
				k4EbTot_x=hpXcalKcorr(spireRefFreq[band],
				  freq, spireFilt[band], BB=True, beta=betaK[b], temp=tempK[t],
				  ext=True, monoArea=beamMonoArea[band])[0]
				kConvExtd[betaTxt][band].data[t] = \
				  k4EbTot_x * effBeamSr[betaTxt][band].data[t]
				#extended source solour correction for current beta,temp
				kCorrExtd[betaTxt][band].data[t]=k4EbTot_x / k4E_Tot[band]

		#-----------------------------------------------------------------------
		#use previous method
		kCorrExtdPrev[betaTxt]=TableDataset()
		kConvExtdPrev[betaTxt]=TableDataset()
		kCorrExtdPrev[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		kConvExtdPrev[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
		for band in spireBands:
			kCorrExtdPrev[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			kConvExtdPrev[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
			for t in range(len(tempK)):
				#point source conversion for current beta,temp (previous method)
				kConvExtdPrev[betaTxt][band].data[t]=\
				  OLD_hpXcalKcorr(spireRefFreq[band],
				  freq, spireFilt[band], BB=True, \
				  beta=betaK[b], temp=tempK[t], gamma=gamma)[0]
				#point source colour correction for current beta,temp (previous method)
				kCorrExtdPrev[betaTxt][band].data[t] = \
				  kConvExtdPrev[betaTxt][band].data[t] / k4E_Prev[band]

	#-----------------------------------------------------------------------
	#write to ascii table files
	kCorrExtdFileList=''
	for name in kCorrExtd.getSets():
		#set filename
		fileKCorrExtd[name]='ColorCorrK_extended_%s_v%s.csv'%(name,version)
		#write to ascii table
		asciiTableWriter(table=kCorrExtd[name],file=os.path.join(dirColorCorrK,fileKCorrExtd[name]))
		#add to list of filenames (for metadata
		if kCorrExtdFileList=='':
			kCorrExtdFileList='%s'%fileKCorrExtd[name]
		else:
			kCorrExtdFileList='%s,%s'%(kCorrExtdFileList,fileKCorrExtd[name])
	kCorrExtd.meta["fileOrigin"] = herschel.ia.dataset.StringParameter(value="%s"%kCorrExtdFileList,\
	  description="Origin of the data")

	# set FITS filename
	kCorrExtdFits = java.io.File(r"%s//SCalPhotColorCorrK_extended_v%s.fits"%(dirColorCorrK, version))
	kCorrExtd.meta['fileName'] = herschel.ia.dataset.StringParameter(value=kCorrExtdFits.name,\
	  description="Name of file when exported")
	kConvExtdFits = java.io.File(r"%s//SCalPhotConvK_extended_v%s.fits"%(dirColorCorrK, version))
	
	#-----------------------------------------------------------------------
	# write log file
	kCorrExtdLog=open(kCorrExtdLogFile,'w')
	kCorrExtdLog.write('Extended Source Colour Correction Version %s\n'%version)
	kCorrExtdLog.write('Description: %s\n'%kCorrExtd.getDescription())
	kCorrExtdLog.write('Creation Date: %s\n'%java.util.Date())
	kCorrExtdLog.write('\nTables:\n')
	for name in kCorrExtd.getSets():
		kCorrExtdLog.write('  %s\n'%(logTable(name,fileKCorrExtd[name],kCorrExtd)))
	kCorrExtdLog.write('\nMetadata:\n')
	kCorrExtdLog.write('  %s\n'%(logMeta('fileOrigin',kCorrExtd.meta,'string')))
	kCorrExtdLog.write('  %s\n'%(logMeta('fileName',kCorrExtd.meta,'string')))
	kCorrExtdLog.close()

	#-----------------------------------------------------------------------
	# write to FITS file
	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(kCorrExtdFits.toString(), kCorrExtd)

	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(kConvExtdFits.toString(), kConvExtd)

	#-----------------------------------------------------------------------
	# plot comparisons
	if plot:
		#plot colour correction
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kCorrExtd['alpha']['alpha'].data,\
			  kCorrExtd['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New'))
			p.addLayer(LayerXY( \
			  kCorrExtdOld['alpha']['alpha'].data,\
			  kCorrExtdOld['alpha'][band].data,\
			  color=cols[band],name=band+' Old'))
			p.addLayer(LayerXY( \
			  kCorrExtdPrev['alpha']['alpha'].data,\
			  kCorrExtdPrev['alpha'][band].data,\
			  color=cols[band],name=band+' Prev',line=Style.DASHED))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Colour correction"
		p.setTitleText("ColorCorrK_extended value")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrExtdOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrK,'SCalPhotColorCorrK_extended_v%s_v%s.png'%(kCorrExtdOldVersion,version)))

		#plot colour correction relative change
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kCorrExtd['alpha']['alpha'].data,\
			  kCorrExtd['alpha'][band].data/kCorrExtdOld['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New/Old'))
			p.addLayer(LayerXY( \
			  kCorrExtdPrev['alpha']['alpha'].data,\
			  kCorrExtdPrev['alpha'][band].data/kCorrExtdOld['alpha'][band].data,\
			  color=cols[band],name=band+' Prev/Old'))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Colour correction w.r.t. v%s"%kCorrExtdOldVersion
		p.setTitleText("ColorCorrK_extended relative change")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrExtdOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrK,'SCalPhotColorCorrK_extended_RelDiff_v%s_v%s.png'%(kCorrExtdOldVersion,version)))

		#plot conversion value
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kConvExtd['alpha']['alpha'].data,\
			  kConvExtd['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New'))
			p.addLayer(LayerXY( \
			  kConvExtdOld['alpha']['alpha'].data,\
			  kConvExtdOld['alpha'][band].data,\
			  color=cols[band],name=band+' Old'))
			p.addLayer(LayerXY( \
			  kConvExtdPrev['alpha']['alpha'].data,\
			  kConvExtdPrev['alpha'][band].data,\
			  color=cols[band],name=band+' Prev',line=Style.DASHED))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Colour conversion"
		p.setTitleText("K4E x ColorCorrK_extended value")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrExtdOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrK,'KConv_extended_v%s_v%s.png'%(kCorrExtdOldVersion,version)))

		#plot conversion relative change
		p=PlotXY()
		for band in spireBands:
			p.addLayer(LayerXY( \
			  kConvExtd['alpha']['alpha'].data,\
			  kConvExtd['alpha'][band].data/kConvExtdOld['alpha'][band].data,\
			  color=cols[band],stroke=2.0,name=band+' New/Old'))
			p.addLayer(LayerXY( \
			  kConvExtdPrev['alpha']['alpha'].data,\
			  kConvExtdPrev['alpha'][band].data/kConvExtdOld['alpha'][band].data,\
			  color=cols[band],name=band+' Prev/Old'))
			p.addLayer(LayerXY( \
			  kConvExtdPrev['alpha']['alpha'].data,\
			  kConvExtd['alpha'][band].data/kConvExtdPrev['alpha'][band].data,\
			  color=cols[band],name=band+' New/Prev',line=Style.DASHED))
		p.xaxis.titleText = "Spectral index (alpha)"
		p.yaxis.titleText = "Colour conversion w.r.t v%s"%kCorrExtdOldVersion
		p.setTitleText("K4E x ColorCorrK_extended relative change ")
		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrExtdOldVersion,version))
		p.legend.visible = 1
		p.saveAsPNG(os.path.join(dirColorCorrK,'KConv_extended_RelDiff_v%s_v%s.png'%(kCorrExtdOldVersion,version)))

	#-----------------------------------------------------------------------
	#=======================================================================
	#=====                 PRINT TO TERMINAL AS CHECK                  =====
	#=======================================================================

	if verbose:
		# Print Spire K4 factors for point source calibration for verification
		print 'v%s SPIRE K4P factors for point source: %6.4f, %6.4f, %6.4f' % \
		  (kCorrPsrcOldVersion,k4P_Old["PSW"],k4P_Old["PMW"],k4P_Old["PLW"])
		print 'v%s SPIRE K4P factors for extended source (old method): %6.4f, %6.4f, %6.4f' % \
		  (version,k4P_Prev["PSW"],k4P_Prev["PMW"],k4P_Prev["PLW"])
		print 'v%s SPIRE K4P factors for point source: %6.4f, %6.4f, %6.4f' % \
		  (version,k4P["PSW"],k4P["PMW"],k4P["PLW"])
	
		# Print Spire K4 factors for extended source calibration for verification
		print 'v%s SPIRE K4E factors for extended source: %6.4f, %6.4f, %6.4f' % \
		  (kCorrExtdOldVersion,k4E_Old["PSW"],k4E_Old["PMW"],k4E_Old["PLW"])
		print 'v%s SPIRE K4E factors for extended source (old method): %6.4f, %6.4f, %6.4f' % \
		  (version,k4E_Prev["PSW"],k4E_Prev["PMW"],k4E_Prev["PLW"])
		print 'v%s SPIRE K4E factors for extended source: %6.4f, %6.4f, %6.4f' % \
		  (version,k4E["PSW"],k4E["PMW"],k4E["PLW"])

		# Print SPIRE beam areas
		print 'v%s pipeline beam areas : %.4f, %.4f, %.4f' %\
		  (kCorrBeamOldVersion,beamAreaPipArc_Old["PSW"],beamAreaPipArc_Old["PMW"],beamAreaPipArc_Old["PLW"])
		print 'v%s pipeline beam areas : %.4f, %.4f, %.4f' %\
		  (version,beamAreaPipArc["PSW"],beamAreaPipArc["PMW"],beamAreaPipArc["PLW"])
		print 'v%s beam areas at reference frequencies: %.4f, %.4f, %.4f' %\
		  (version,beamRefArea["PSW"]/arcsec2Sr,beamRefArea["PMW"]/arcsec2Sr,beamRefArea["PLW"]/arcsec2Sr)



else:
	# use old values
	kCorrBeam = calibration.phot.getProduct('ColorCorrBeam')
	kCorrBeamVersion = kCorrBeam.getVersion()
	print 'Using pipeline beam areas from ColorCorrBeam version %s from calibration %s)'% \
	  (kCorrBeamVersion,calVersion)
	beamAreaPipSr["PSW"]=kCorrBeam.meta['beamPswSr'].double
	beamAreaPipSr["PMW"]=kCorrBeam.meta['beamPmwSr'].double
	beamAreaPipSr["PLW"]=kCorrBeam.meta['beamPlwSr'].double
	beamAreaPipArc["PSW"]=kCorrBeam.meta['beamPswArc'].double
	beamAreaPipArc["PMW"]=kCorrBeam.meta['beamPmwArc'].double
	beamAreaPipArc["PLW"]=kCorrBeam.meta['beamPlwArc'].double
	
	#calculate effBeamSr	
	effBeamSr=kCorrBeam.copy()
	for name in kCorrBeam.getSets():
		for band in spireBands:
			effBeamSr[name][band].data = \
			  beamAreaPipSr[band] / Float1d(kCorrBeam[name][band].data)
	
	kCorrPsrc = calibration.phot.getProduct('ColorCorrKList')[1]
	kCorrExtd = calibration.phot.getProduct('ColorCorrKList')[0]
	kCorrPsrcVersion=kCorrPsrc.getVersion()
	kCorrExtdVersion=kCorrExtd.getVersion()
	print 'Using point source colour corrections from ColorCorrK_point version %s from calibration %s'%\
	  (kCorrPsrcVersion,calVersion)
	print 'Using extended source colour corrections from ColorCorrK_point version %s from calibration %s'%\
	  (kCorrExtdVersion,calVersion)

	fluxConv = calibration.phot.getProduct('FluxConvList')[0]
	fluxConvVersion=fluxConv.getVersion()
	#get old K4P & K4E pipeline values from calibration file
	print 'Using pipeline K4-parameters from FluxConvList version %s from calibration %s'%\
	  (fluxConvVersion,calVersion)
	for band in spireBands:
		k4P[band]=fluxConv[band].meta['k4P_%s'%band].double
		k4E[band]=fluxConv[band].meta['k4P_%s'%band].double
		k4E_Tot[band]=k4E[band]/beamAreaPipSr[band]
		k4P_Prev[band]=k4P[band]
		k4E_Prev[band]=k4E[band]

#-------------------------------------------------------------------------------
#===============================================================================
#=====                    WRITE BEAM PROFILES TO FITS FILE                 =====
#===============================================================================
#-------------------------------------------------------------------------------
# delayed from earlier due to additional metadata
if calcRadialCorrBeam:
	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(beamProfsFits.toString(), beamProfs)

#-------------------------------------------------------------------------------
#===============================================================================
#=====                     CALCULATE APERTURE CORRECTIONS                  =====
#===============================================================================
#-------------------------------------------------------------------------------

if calcColorCorrAperture:
	#set aperture photometry radius and annulus size
	apPhotRad={"PSW":22.,"PMW":30.,"PLW":45.}
	apPhotBGRad={'in':60.,'out':90.}
	#read in map from file
	beamNames = {"PSW":"0x5000241aL_PSW_pmcorr_1arcsec_norm_beam.fits",
		"PMW":"0x5000241aL_PMW_pmcorr_1arcsec_norm_beam.fits",
		"PLW":"0x5000241aL_PLW_pmcorr_1arcsec_norm_beam.fits"}
	try:
		beamIn = fitsReader(file = os.path.join(dataDir,"0x5000241aL_PSW_pmcorr_1arcsec_norm_beam.fits"))
	except:
		#download if not available
		urllib.urlretrieve ("https://nhscsci.ipac.caltech.edu/spire/data/beam_profiles/"+beamNames["PSW"],\
		    os.path.join(dataDir,beamNames[band]))
		beamIn = fitsReader(file = os.path.join(dataDir,beamNames[band]))

	# make map of radius (speeds up processing later)
	print 'Making beam Radius map'
	beamRadMap=SimpleImage()
	beamRadMap['image']=beamIn['image']
	nxMap=beamIn['image'].data[:,0].size
	nyMap=beamIn['image'].data[0,:].size
	bcenter=[int(nxMap/2.),int(nyMap/2.)]
	for x in range(nxMap):
		for y in range(nyMap):
			beamRadMap['image'].data[x,y]= \
				Math.sqrt((x-bcenter[0])**2 + (y-bcenter[1])**2)

	# set up product for no background
	apCorrNoBG = herschel.spire.ia.dataset.PhotColorCorrAperture() 
	apCorrNoBG.setDescription('Spire aperture correction product')
	apCorrNoBG.meta["creator"].value = scriptVersionString
	apCorrNoBG.meta["creationDate"].value = FineTime(java.util.Date())
	apCorrNoBG.meta["startDate"].value = FineTime(startDate)
	apCorrNoBG.meta["endDate"].value = FineTime(endDate)
	apCorrNoBG.setVersion(version)
	apCorrNoBG.setFormatVersion(formatVersion)

	#copy product for including background
	apCorrIncBG = apCorrNoBG.copy()

	#set up log files
	apCorrNoBGLogFile=os.path.join(dirColorCorrAperture,'ColorCorrAperture_noBG_log_v%s.dat'%version)
	apCorrIncBGLogFile=os.path.join(dirColorCorrAperture,'ColorCorrAperture_incBG_log_v%s.dat'%version)

	#-----------------------------------------------------------------------
	# Compute aperture corrections for range of alpha
	# Create tables for alpha arrays for aperture corrections
	apCorrNoBG['alpha']=TableDataset(description='Aperture Correction without background (Spectral Index)')
	apCorrIncBG['alpha']=TableDataset(description='Aperture Correction including background (Spectral Index)')
	apCorrNoBG['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	apCorrIncBG['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
	#-----------------------------------------------------------------------
	print '\nCalculating aperture corrections corrections over alpha...'
	for band in spireBands:
		apCorrNoBG['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		apCorrIncBG['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
		for a in range(len(alphaK)):
			print '%s alpha=%.1f'%(band,alphaK[a])
			#calculate beam areas
			#effBeamPSW contains (beam area, beam profile, beam map)
			effBeam_x=spireEffBeam(freq,spireFiltOnly[band],beamProfs,spireEffFreq[band],\
			  gamma,band,beamRadMap,BB=False,alpha=alphaK[a],verbose=verbose)
			#perform aperture photometry
			apPhot_x = annularSkyAperturePhotometry(image=effBeam_x['map'], \
			  fractional=1, centerX=bcenter[0], centerY=bcenter[1], \
			  radiusArcsec=apPhotRad[band], \
			  innerArcsec=apPhotBGRad['in'], outerArcsec=apPhotBGRad['out'])
			#get result of aperture correction procedure
			apPhotIncBG_x = apPhot_x.getTargetTotal()
			apPhotNoBG_x = apPhot_x.getTargetPlusSkyTotal()
			#compute correction to provide actual beam area
			apCorrIncBG['alpha'][band].data[a]=effBeam_x['area']/apPhotIncBG_x
			apCorrNoBG['alpha'][band].data[a]=effBeam_x['area']/apPhotNoBG_x


	#-----------------------------------------------------------------------
	#write to ascii table files for no background
	apCorrNoBGFileList=''
	fileApCorrNoBG={}
	for name in apCorrNoBG.getSets():
		#set filename
		fileApCorrNoBG[name]='ColorCorrAperture_noBG_%s_v%s.csv'%(name,version)
		#write to ascii table
		asciiTableWriter(table=apCorrNoBG[name],file=os.path.join(dirColorCorrAperture,fileApCorrNoBG[name]))
		#add to list of filenames (for metadata
		if apCorrNoBGFileList=='':
			apCorrNoBGFileList='%s'%fileApCorrNoBG[name]
		else:
			apCorrNoBGFileList='%s,%s'%(apCorrNoBGFileList,fileApCorrNoBG[name])
	apCorrNoBG.meta["fileOrigin"] = herschel.ia.dataset.StringParameter(value="%s"%apCorrNoBGFileList,\
	  description="Origin of the data")

	#-----------------------------------------------------------------------
	#write to ascii table files for including background
	apCorrIncBGFileList=''
	fileApCorrIncBG={}
	for name in apCorrIncBG.getSets():
		#set filename
		fileApCorrIncBG[name]='ColorCorrAperture_incBG_%s_v%s.csv'%(name,version)
		#write to ascii table
		asciiTableWriter(table=apCorrIncBG[name],file=os.path.join(dirColorCorrAperture,fileApCorrIncBG[name]))
		#add to list of filenames (for metadata
		if apCorrIncBGFileList=='':
			apCorrIncBGFileList='%s'%fileApCorrIncBG[name]
		else:
			apCorrIncBGFileList='%s,%s'%(apCorrIncBGFileList,fileApCorrIncBG[name])
	apCorrIncBG.meta["fileOrigin"] = herschel.ia.dataset.StringParameter(value="%s"%apCorrIncBGFileList,\
	  description="Origin of the data")

	# set FITS filenames
	apCorrNoBGFits = java.io.File(r"%s//SCalPhotColorCorrAperture_noBG_v%s.fits"%(dirColorCorrAperture, version))
	apCorrIncBGFits = java.io.File(r"%s//SCalPhotColorCorrAperture_incBG_v%s.fits"%(dirColorCorrAperture, version))
	apCorrNoBG.meta['fileName'] = herschel.ia.dataset.StringParameter(value=apCorrNoBGFits.name,\
	  description="Name of file when exported")
	apCorrIncBG.meta['fileName'] = herschel.ia.dataset.StringParameter(value=apCorrIncBGFits.name,\
	  description="Name of file when exported")
	
	#-----------------------------------------------------------------------
	# write log file for No BG
	apCorrNoBGLog=open(apCorrNoBGLogFile,'w')
	apCorrNoBGLog.write('Aperture Correction Version %s\n'%version)
	apCorrNoBGLog.write('Description: %s\n'%apCorrNoBG.getDescription())
	apCorrNoBGLog.write('Creation Date: %s\n'%java.util.Date())
	apCorrNoBGLog.write('\nTables:\n')
	for name in apCorrNoBG.getSets():
		apCorrNoBGLog.write('  %s\n'%(logTable(name,fileApCorrNoBG[name],apCorrNoBG)))
	apCorrNoBGLog.write('\nMetadata:\n')
	apCorrNoBGLog.write('  %s\n'%(logMeta('fileOrigin',apCorrNoBG.meta,'string')))
	apCorrNoBGLog.write('  %s\n'%(logMeta('fileName',apCorrNoBG.meta,'string')))
	apCorrNoBGLog.close()

	#-----------------------------------------------------------------------
	# write log file for Inc BG
	apCorrIncBGLog=open(apCorrIncBGLogFile,'w')
	apCorrIncBGLog.write('Aperture Correction Version %s\n'%version)
	apCorrIncBGLog.write('Description: %s\n'%apCorrIncBG.getDescription())
	apCorrIncBGLog.write('Creation Date: %s\n'%java.util.Date())
	apCorrIncBGLog.write('\nTables:\n')
	for name in apCorrIncBG.getSets():
		apCorrIncBGLog.write('  %s\n'%(logTable(name,fileApCorrIncBG[name],apCorrIncBG)))
	apCorrIncBGLog.write('\nMetadata:\n')
	apCorrIncBGLog.write('  %s\n'%(logMeta('fileOrigin',apCorrIncBG.meta,'string')))
	apCorrIncBGLog.write('  %s\n'%(logMeta('fileName',apCorrIncBG.meta,'string')))
	apCorrIncBGLog.close()

	#-----------------------------------------------------------------------
	# write to FITS file for No BG
	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(apCorrNoBGFits.toString(), apCorrNoBG)

	#-----------------------------------------------------------------------
	# write to FITS file for Inc BG
	fitsWriter = FitsArchive()
	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
	fitsWriter.save(apCorrIncBGFits.toString(), apCorrIncBG)

##-------------------------------------------------------------------------------
##===============================================================================
##=====                     CALCULATE SPIRE-HFI CONVERSION                  =====
##===============================================================================
##-------------------------------------------------------------------------------
#if calcColorCorrHfi:
#	# copy old calibration product
#	kCorrHfiOld = calibration.phot.getProduct('ColorCorrHfi')
#	kCorrHfiOldVersion = kCorrHfiOld.getVersion()
#
#	print "Calculating HFI-SPIRE cross-calibration conversion factors"
#	#-------------------------------------------------------------------------------
#	# Load HFI filter functions
#	hfiRIMO = fitsReader(file = hfiRIMOFile)
#	
#	# Read and parse HFI filter files into tables 
#	hfiFreq545 = 100. * c * hfiRIMO['BANDPASS_F545']['WAVENUMBER'].data[1:-1]
#	hfiFreq857 = 100. * c * hfiRIMO['BANDPASS_F857']['WAVENUMBER'].data[1:-1]
#	
#	hfiTrans545 = Double1d(hfiRIMO['BANDPASS_F545']['TRANSMISSION'].data[1:-1])
#	hfiTrans857 = Double1d(hfiRIMO['BANDPASS_F857']['TRANSMISSION'].data[1:-1])
#	
#	
#	# Excluding data points that are not monotonically increasing for HFI-545
#	diff = Double1d(len(hfiFreq545)-1)
#	for i in range(len(hfiFreq545)-1):
#		diff[i] = hfiFreq545[i+1]-hfiFreq545[i]
#	
#	ind = diff.where(diff <= 0.).toInt1d()
#	for i in ind:
#		hfiFreq545.delete(i+1,1)
#		hfiTrans545.delete(i+1,1)
#	
#	
#	# Excluding data points that are not monotonically increasing for HFI-545
#	diff = Double1d(len(hfiFreq857)-1)
#	for i in range(len(hfiFreq857)-1):
#		diff[i] = hfiFreq857[i+1] - hfiFreq857[i]
#	
#	ind = diff.where(diff <= 0.).toInt1d()
#	for i in ind:
#		hfiFreq857.delete(i+1,1)
#		hfiTrans857.delete(i+1,1)
#	
#	
#	# Interpolate HFI transmissions to common frequency grid 
#	hfiFilt545 = Double1d(nNu)
#	hfiFilt857 = Double1d(nNu)
#	
#	ix545     = freq.where((freq>=MIN(hfiFreq545)) & (freq<=MAX(hfiFreq545)))
#	interp545 = CubicSplineInterpolator(hfiFreq545, hfiTrans545)
#	hfiFilt545[ix545] = interp545(freq[ix545])
#	
#	ix857     = freq.where((freq>=MIN(hfiFreq857)) & (freq<=MAX(hfiFreq857)))
#	interp857 = CubicSplineInterpolator(hfiFreq857, hfiTrans857)
#	hfiFilt857[ix857] = interp857(freq[ix857])
#	
#	# Normalize transmissions
#	hfiFilt545 = hfiFilt545 / MAX(hfiFilt545)
#	hfiFilt857 = hfiFilt857 / MAX(hfiFilt857)
#	
#	# Calculate K-correction factors for extended source assuming alpha=-1
#	# Note that since HFI beams do not vary with frequency, extended sources
#	#   are treated like point sources in this context
#	k4E_857 = hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857,  BB=False)[0]
#	k4E_545 = hpXcalKcorr(hfiRefFreq[1],   freq, hfiFilt545,  BB=False)[0]
##
##	#-------------------------------------------------------------------------------
##	# Plot relative spectral response functions
##	
##	if plot:
##		# Make linear plot of all filter bands including those for point sources
##		p = PlotXY()
##		p.addLayer(LayerXY(c/freq*1e6, hfiFilt857, name='HFI 857GHz',\
##		  color=java.awt.Color.MAGENTA))
##		p.addLayer(LayerXY(c/freq*1e6, hfiFilt545, name='HFI 545GHz',\
##		  color=java.awt.Color.CYAN))
##		#
##		p.addLayer(LayerXY(c/freq*1e6, spireFiltOnly['PSW'], name='PSW RSRF',\
##		  color=java.awt.Color.BLUE,line=Style.SOLID))
##		p.addLayer(LayerXY(c/freq*1e6, spireFiltOnly['PMW'], name='PMW RSRF',\
##		  color=java.awt.Color.GREEN,line=Style.SOLID))
##		p.addLayer(LayerXY(c/freq*1e6, spireFiltOnly['PLW'], name='PLW RSRF',\
##		  color=java.awt.Color.RED,line=Style.SOLID))
##		#
##		p.addLayer(LayerXY(c/freq*1e6, spireFilt['PSW'], name='PSW RSRFxApEff',\
##		  color=java.awt.Color.BLUE,line=Style.DASHED))
##		p.addLayer(LayerXY(c/freq*1e6, spireFilt['PMW'], name='PMW RSRFxApEff',\
##		  color=java.awt.Color.GREEN,line=Style.DASHED))
##		p.addLayer(LayerXY(c/freq*1e6, spireFilt['PLW'], name='PLW RSRFxApEff',\
##		  color=java.awt.Color.RED,line=Style.DASHED))
##		#
##		p.xaxis.range = [100,800]
##		p.yaxis.range = [-0.1,1.1]
##		p.xaxis.titleText = "Wavelength [micron]"
##		p.yaxis.titleText = "Relative Spectral Response"
##		p.legend.visible = 1
##		#	
#	#-------------------------------------------------------------------------------
#	# Calculate and tabulate colour correction parameters for a range of temperatures
#	k545toPLW = Double1d()  # K-correction from HFI-545 to PLW
#	k857toPMW = Double1d()  # K-correction from HFI-857 to PMW
#	k857toPSW = Double1d()  # K-correction from HFI-857 to PSW
#	ratio545over857 = Double1d()  # Ratio 545GHz to 857GHz filter
#	#create arrays for previous method
#	k545toPLWPrev = Double1d()  # K-correction from HFI-545 to PLW
#	k857toPMWPrev = Double1d()  # K-correction from HFI-857 to PMW
#	k857toPSWPrev = Double1d()  # K-correction from HFI-857 to PSW
#	#loat old versions
#	kCorrHfiOldTemp = Double1d(kCorrHfiOld['colorCorr']['Temperature'].data)
#	k545toPLWOld = Double1d(kCorrHfiOld['colorCorr']['k545toPLW'].data)
#	k857toPMWOld = Double1d(kCorrHfiOld['colorCorr']['k857toPMW'].data)
#	k857toPSWOld = Double1d(kCorrHfiOld['colorCorr']['k857toPSW'].data)
#	ratio545over857Old = Double1d(kCorrHfiOld['colorCorr']['ratio545_857'].data)
#	if fixBeta == True:
#		#
#		# Make temperature vector with log distances between tiers
#		#tvect = Double1d(range(ndiv)) * (Tmax - Tmin) / (ndiv-1) + Tmin
#		tvect = Double1d(tempK)
#		#
#		for temp in tvect:
#			#
#			# K-correction from HFI-545 to PLW
#			# calculate K-correction for HFI 545
#			k545 = Float1d(hpXcalKcorr(hfiRefFreq[1], freq, hfiFilt545,
#			  BB=True, temp=temp, beta=beta0))
#			# calculate SPIRE PLW extended conversion
#			kPLW = Float1d(hpXcalKcorr(spireRefFreq["PLW"], freq, spireFilt["PLW"], 
#			  BB=True, temp=temp, beta=beta0, ext=True,
#			  monoArea=beamMonoArea["PLW"]))
#			# multiply by effective beam area
#			kPLW[0] = kPLW[0] * spireEffArea(freq, spireFiltOnly["PLW"],
#			  beamMonoArea["PLW"], BB=True, temp=temp, beta=beta0)
#			
#			k545toPLW.append(k545[0] / kPLW[0] * k4E["PLW"] / k4E_545 * kPLW[1] / k545[1])
#			
#			# K-correction from HFI-857 to PMW
#			# calculate K-correction for HFI 857
#			k857 = Float1d(hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857, 
#			  BB=True, temp=temp, beta=beta0))
#			# calculate SPIRE PMW extended conversion
#			kPMW = Float1d(hpXcalKcorr(spireRefFreq["PMW"], freq, spireFilt["PMW"], 
#			  BB=True, temp=temp, beta=beta0, ext=True,
#			  monoArea=beamMonoArea["PMW"]))
#			# multiply by effective beam area
#			kPMW[0] = kPMW[0] * spireEffArea(freq, spireFiltOnly["PMW"],
#			  beamMonoArea["PMW"], BB=True, temp=temp, beta=beta0)
#			
#			k857toPMW.append(k857[0] / kPMW[0] * k4E["PMW"] / k4E_857 * kPMW[1] / k857[1])
#			
#			#
#			# K-correction from HFI-857 to PSW
#			# calculate SPIRE PMW extended conversion
#			kPSW = Float1d(hpXcalKcorr(spireRefFreq["PSW"], freq, spireFilt["PSW"], 
#			  BB=True, temp=temp, beta=beta0, ext=True,
#			  monoArea=beamMonoArea["PSW"]))
#			# multiply by effective beam area
#			kPSW[0] = kPSW[0] * spireEffArea(freq, spireFiltOnly["PSW"],
#			  beamMonoArea["PSW"], BB=True, temp=temp, beta=beta0)
#			k857toPSW.append(k857[0] / kPSW[0] * k4E["PSW"] / k4E_857 * kPSW[1] / k857[1])
#			#
#			# Ratio of 545 and 845 GHz filters
#			ratio545over857.append(k545[1] / k857[1] * k4E_545 / k4E_857 * k857[0] / k545[0])
#
#			#-------------------------------------------------------
#			# use previous method
#			k545 = OLD_hpXcalKcorr(hfiRefFreq[1], freq, hfiFilt545,
#			  BB=True, temp=temp, beta=beta0)
#			kPLW = OLD_hpXcalKcorr(spireRefFreq["PLW"], freq, \
#			  spireFilt["PLW"], BB=True, temp=temp, beta=beta0, gamma=gamma)
#			k545toPLWPrev.append(k545[0] / kPLW[0] * k4E_Prev["PLW"] / k4E_545 * kPLW[1] / k545[1])
#
#			k857 = OLD_hpXcalKcorr(hfiRefFreq[0], freq, hfiFilt857, 
#			  BB=True, temp=temp, beta=beta0)
#			kPMW = OLD_hpXcalKcorr(spireRefFreq["PMW"], freq, \
#			  spireFilt["PMW"], BB=True, temp=temp, beta=beta0, gamma=gamma)
#			k857toPMWPrev.append(k857[0] / kPMW[0] * k4E_Prev["PMW"] / k4E_857 * kPMW[1] / k857[1])
#
#			kPSW = OLD_hpXcalKcorr(spireRefFreq["PSW"], freq, \
#			  spireFilt["PSW"], BB=True, temp=temp, beta=beta0, gamma=gamma)
#			k857toPSWPrev.append(k857[0] / kPSW[0] * k4E_Prev["PSW"] / k4E_857 * kPSW[1] / k857[1])
#
#		#
#	else:
#		#
#		# Make temperature vector with log distances between tiers
#		bvect = Double1d(range(ndiv)) * (Bmax - Bmin) / (ndiv-1) + Bmin
#		#
#		for beta in bvect:
#			#
#			# K-correction from HFI-545 to PLW
#			k545 = Double1d(hpXcalKcorr(hfiRefFreq[1],   freq, hfiFilt545, 
#			  BB=True, temp=temp0, beta=beta))
#			# calculate SPIRE PLW extended conversion
#			kPLW = Double1d(hpXcalKcorr(spireRefFreq["PLW"], freq, spireFilt["PLW"], 
#			  BB=True, temp=temp0, beta=beta, ext=True,
#			  monoArea=beamMonoArea["PLW"]))
#			# multiply by effective beam area
#			kPLW[0] = kPLW[0] * spireEffArea(freq, spireFiltOnly["PLW"],
#			  beamMonoArea["PLW"], BB=True, temp=temp0, beta=beta)
#			k545toPLW.append(k545[0] / kPLW[0] * k4E["PLW"] / k4E_545 * kPLW[1] / k545[1])
#			#
#			# K-correction from HFI-857 to PMW
#			k857 = Double1d(hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857,
#			  BB=True, temp=temp0, beta=beta))
#			# calculate SPIRE PMW extended conversion
#			kPMW = Double1d(hpXcalKcorr(spireRefFreq["PMW"], freq, spireFilt["PMW"], 
#			  BB=True, temp=temp0, beta=beta, ext=True,
#			  monoArea=beamMonoArea["PMW"]))
#			# multiply by effective beam area
#			kPMW[0] = kPMW[0] * spireEffArea(freq, spireFiltOnly["PMW"],
#			  beamMonoArea["PMW"], BB=True, temp=temp0, beta=beta)
#			k857toPMW.append(k857[0] / kPMW[0] * k4E["PMW"] / k4E_857 * kPMW[1] / k857[1])
#			#
#			# K-correction from HFI-857 to PSW
#			# calculate SPIRE PMW extended conversion
#			kPSW = Double1d(hpXcalKcorr(spireRefFreq["PSW"], freq, spireFilt["PSW"], 
#			  BB=True, temp=temp0, beta=beta, ext=True,
#			  monoArea=beamMonoArea["PSW"]))
#			# multiply by effective beam area
#			kPSW = kPSW * spireEffArea(freq, spireFiltOnly["PSW"],
#			  beamMonoArea["PSW"], BB=True, temp=temp0, beta=beta)
#			k857toPSW.append(k857[0] / kPSW[0] * k4E["PSW"] / k4E_857 * kPSW[1] / k857[1])
#			#
#			# Ratio of 545 and 845 GHz filters
#			ratio545over857.append(k545[1] / k857[1] * k4E_545 / k4E_857 * k857[0] / k545[0])
#
#			#-------------------------------------------------------
#			# use previous method
#			k545 = OLD_hpXcalKcorr(hfiRefFreq[1], freq, hfiFilt545, 
#			  BB=True, temp=temp0, beta=beta)
#			kPLW = OLD_hpXcalKcorr(spireRefFreq["PLW"], freq, 
#			  spireFilt["PLW"], BB=True, temp=temp0, beta=beta, gamma=gamma)
#			k545toPLWPrev.append(k545[0] / kPLW[0] * k4E_Prev["PLW"] / k4E_545 * kPLW[1] / k545[1])
#
#			k857 = OLD_hpXcalKcorr(hfiRefFreq[0], freq, hfiFilt857,
#			  BB=True, temp=temp0, beta=beta)
#			kPMW = OLD_hpXcalKcorr(spireRefFreq["PMW"], freq,
#			  spireFilt["PMW"], BB=True, temp=temp0, beta=beta, gamma=gamma)
#			k857toPMWPrev.append(k857[0] / kPMW[0] * k4E_Prev["PMW"] / k4E_857 * kPMW[1] / k857[1])
#
#			kPSW = OLD_hpXcalKcorr(spireRefFreq["PSW"], freq, 
#			  spireFilt["PSW"], BB=True, temp=temp0, beta=beta, gamma=gamma)
#			k857toPSWPrev.append(k857[0] / kPSW[0] * k4E_Prev["PSW"] / k4E_857 * kPSW[1] / k857[1])
#
#	#-------------------------------------------------------------------------------
#		
#	#-------------------------------------------------------------------------------
#	# Create calibration product with tabulated colour correction factors
#	#-------------------------------------------------------------------------------
#
#	#create blank calibration product
#	kCorrHfi = herschel.spire.ia.dataset.PhotColorCorrHfi()
#	kCorrHfi.meta["creator"].value   = scriptVersionString
#	kCorrHfi.meta["modelName"].value = "FM"
#	kCorrHfi.meta["creationDate"].value = FineTime(java.util.Date())
#	kCorrHfi.meta["startDate"].value = FineTime(startDate)
#	kCorrHfi.meta["endDate"].value   = FineTime(endDate)
#	kCorrHfi.meta["fileOrigin"]  = herschel.ia.dataset.StringParameter(value="%s"%hfiRIMOFileName, description="Origin of the data")
#	kCorrHfi.setVersion(version)
#	kCorrHfi.setFormatVersion(formatVersion)
#	
#	# Save beta dependent colour correction tables into binary table
#	if fixBeta == True:
#		kCorrHfi.meta["beta"] = DoubleParameter(beta0, "Modified black-body spectral index used")
#	        kCorrHfi.setTempVals(tvect)      # Temperature of modified BB
#	        #
#	# Save temperature dependent colour correction tables into binary table
#	else:
#		kCorrHfi.meta['temperature'] = DoubleParameter(temp0, "Modified black-body temperature used")
#		kCorrHfi.meta['temperature'].unit = Temperature.KELVIN
#		#
#		kCorrHfi.setTempVals(bvect)          # Spectral index of modified BB
#	
#	# Save standard point & extended source correction factor into meta data
#	# *** NOT IN NEW VERSIONS ***
#	#for band in spireBands:
#	#	kCorrHfi.meta["k4P_%s"%band] = DoubleParameter(k4P[band], "%s point source K-factor")
#	#	kCorrHfi.meta["k4E_%s"%band] = DoubleParameter(k4E[band], "%s extended source K-factor")
#	
#	# save  gamma to metadata
#	# ***NOT IN NEW VERSIONS ***
#	#kCorrHfi.meta["gamma"] = DoubleParameter(gamma, "Exponent describing FWHM dependence on frequency used")
#	
#	# Save colour correction tables into binary table
#	kCorrHfi.setRatio545_857CorrVals(ratio545over857)  # Ratio 545/857 GHz filter
#	kCorrHfi.setK545toPLWCorrVals(k545toPLW)           # 545GHz to PLW K-factor
#	kCorrHfi.setK857toPMWCorrVals(k857toPMW)           # 857GHz to PMW K-factor
#	kCorrHfi.setK857toPSWCorrVals(k857toPSW)           # 857GHz to PSW K-factor
#
#	#make dummy table for previous method
#	kCorrHfiPrev = herschel.spire.ia.dataset.PhotColorCorrHfi()
#	kCorrHfiPrev.setRatio545_857CorrVals(ratio545over857)  # Ratio 545/857 GHz filter
#	kCorrHfiPrev.setK545toPLWCorrVals(k545toPLWPrev)           # 545GHz to PLW K-factor
#	kCorrHfiPrev.setK857toPMWCorrVals(k857toPMWPrev)           # 857GHz to PMW K-factor
#	kCorrHfiPrev.setK857toPSWCorrVals(k857toPSWPrev)           # 857GHz to PSW K-factor
#	if fixBeta:
#		kCorrHfiPrev.setTempVals(tvect)
#	else:
#		kCorrHfiPrev.setTempVals(bvect)
#
#	#get table from product
#	kCorrHfiTable=kCorrHfi['colorCorr']
#	#write ascii files
#	asciiFileHfi='ColorCorrHfi_v%s.csv'%version
#	asciiTableWriter(table=kCorrHfiTable,file=os.path.join(dirColorCorrHfi,asciiFileHfi))
#	
#	#write log file
#	kCorrHfiLogFile=os.path.join(dirColorCorrK,'ColorCorrHfi_log_v%s.dat'%version)
#	# write log file
#	kCorrHfiLog=open(kCorrHfiLogFile,'w')
#	kCorrHfiLog.write('HFI-SPIRE Colour Correction Version %s\n'%version)
#	kCorrHfiLog.write('Description: %s\n'%kCorrExtd.getDescription())
#	kCorrHfiLog.write('Creation Date: %s\n'%java.util.Date())
#	kCorrHfiLog.write('\nTables:\n')
#	kCorrHfiLog.write('  %s\n'%(logTable('colorCorr',asciiFileHfi,kCorrHfi)))
#	kCorrHfiLog.write('\nMetadata:\n')
#	if fixBeta:
#		kCorrHfiLog.write('  %s\n'%(logMeta('beta',kCorrHfi.meta,'double')))
#	else:
#		kCorrHfiLog.write('  %s\n'%(logMeta('temp',kCorrHfi.meta,'double')))
#	kCorrHfiLog.close()
#
#	##########################################
#	# ******* write to FITS *******
#	filename = java.io.File(r"%s//Phot//SCalPhotColorCorrHfi//SCalPhotColorCorrHfi_v%s.fits"%(directory, version))
#	#filename = java.io.File(r"%s/SCalPhotColorCorrHfi_v%s.fits"%(directory, version))
#	kCorrHfi.meta['fileName'] = herschel.ia.dataset.StringParameter(value=filename.name)
#	fitsWriter = FitsArchive()
#	fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
#	fitsWriter.save(filename.toString(), kCorrHfi)
#
#	if plot:
#		# Plot SPIRE-HFI conversion values
#		p=PlotXY()
#		colsHfi={'k545toPLW':java.awt.Color.RED,\
#		  'k857toPMW':java.awt.Color.GREEN,\
#		  'k857toPSW':java.awt.Color.BLUE}
#		for col in ['k545toPLW','k857toPMW','k857toPSW']:
#			p.addLayer(LayerXY(kCorrHfi['colorCorr']['ratio545_857'].data,\
#			  kCorrHfi['colorCorr'][col].data,name=col+' New',\
#			  color=colsHfi[col],line=Style.SOLID,stroke=2.0))
#			p.addLayer(LayerXY(kCorrHfiOld['colorCorr']['ratio545_857'].data,\
#			  kCorrHfiOld['colorCorr'][col].data,name=col+' Old',\
#			  color=colsHfi[col],line=Style.SOLID))
#			p.addLayer(LayerXY(kCorrHfiPrev['colorCorr']['ratio545_857'].data,\
#			  kCorrHfiPrev['colorCorr'][col].data,name=col+' Prev',\
#			  color=colsHfi[col],line=Style.DASHED))
#
#		p.xaxis.titleText = "ratio 545/857"
#		p.yaxis.titleText = "SPIRE-HFI conversion"
#		p.setTitleText("ColorCorrHfi value ")
#		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrHfiOldVersion,version))
#		p.legend.visible = 1
#		p.saveAsPNG(os.path.join(dirColorCorrHfi,'SCalPhotColorCorrHfi_v%s_v%s.png'%(kCorrHfiOldVersion,version)))
#		#
#		# Plot SPIRE-HFI conversion relative change
##		p=PlotXY()
##
##		for col in ['k545toPLW','k857toPMW','k857toPSW']:
##			p.addLayer(LayerXY(kCorrHfi['colorCorr']['ratio545_8beamNewFileConstant57'].data,\
##			  kCorrHfi['colorCorr'][col].data/kCorrHfiOld['colorCorr'][col].data,\
##			  name=col+' New/Old',color=colsHfi[col],stroke=2.0))
##			p.addLayer(LayerXY(kCorrHfiPrev['colorCorr']['ratio545_857'].data,\
##			  kCorrHfiPrev['colorCorr'][col].data/kCorrHfiOld['colorCorr'][col].data,\
##			  name=col+' Prev/Old',color=colsHfi[col],line=Style.DASHED))
##
##		p.xaxis.titleText = "ratio 545/857"
##		p.yaxis.titleText = "SPIRE-HFI conversion w.r.t v%s"%kCorrHfiOldVersion
##		p.setTitleText("ColorCorrHfi value ")
##		p.setSubtitleText("Old:v%s vs New:v%s"%(kCorrHfiOldVersion,version))
##		p.legend.visible = 1
##		p.saveAsPNG(os.path.join(dirColorCorrHfi,'SCalPhotColorCorrHfi_RelDiff_v%s_v%s.png'%(kCorrHfiOldVersion,version)))
##		#
#
##-------------------------------------------------------------------------------
## End of File
