#===============================================================================
# 
#  Planck-HFI / Herschel-SPIRE Colour Correction factors 
# 
#  This routine calculates a table of colour correction factors for a range of
#  realistic celestial background spectra. The factors transfer standard 
#  monochromatic fluxdensities from Planck-HFI all-sky maps at 857 and 545 GHz
#  which are quoted for a nu*F_nu = const. spectrum into Herschel-SPIRE 
#  monochromatic flux densities with the same definition but reference 
#  wavelengths at 250, 350, and 500 micron. It is assumed that the actual 
#  celestial background is described by the Planck-Law multiplied by frequency 
#  to the power of beta which is set to 1.8. In addition the flux ratios as 
#  measured by HFI between the 545 and 857 Gz filters is given such that the 
#  colour corrections can be determined directly as a function of the flux 
#  ratios in those maps.
# 
#  Input:
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


#-------------------------------------------------------------------------------
# Input parameters

# Colour correction table version
ver = 2.5

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
import os
workDir = Configuration.getProperty('var.hcss.workdir')
hfiRIMOFile = os.path.join(workDir,'HFI_RIMO_R1.10.fits')

simpleBeam=False
calcEffFreq=False

# Name of output table
outFilename = os.path.join(workDir,'SpireHfiColourCorrTab_v%3.1f.fits'%ver)

# Choose if you want to plot the computed RSRFs and colour correction parameters
plot = False


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

# Exponent of powerlaw describing FWHM dependence on frequency
# FWHM ~ frequ**gamma
gamma = -0.85 

# Three SPIRE filter reference frequencies for PSW, PMW, PLW respectively
spireRefFreq = c/Double1d([250.,350.,500.])*1e6 

# Two HFI filter reference frequencies for the 857 and 545 GHz filters respectively
hfiRefFreq = Double1d([857.,545.])*1e9       


#-------------------------------------------------------------------------------
# Load SPIRE filter functions

# Get SPIRE calibration tree
#calibration = SpireCal.getInstance()
try:
	calibration
except:
	calibration = spireCal(jarFile=os.path.join(workDir,'spire_cal_12_0.jar'))
	#calibration = SpireCal.getInstance()
rsrf        = calibration.phot.rsrf						# Load relative spectral response function
spireFreq   = rsrf['rsrf']['frequency'].data*1e9	# Frequency in Hz

ix = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))

#-------------------------------------------------------------------------------
# Get pipeline beam areas from Calibration tree
beamCorrTable  = calibration.phot.refs["ColorCorrBeam"].product
beamAreaPipArc  = [beamCorrTable.meta["beamPswArc"].double,
	beamCorrTable.meta["beamPmwArc"].double,
	beamCorrTable.meta["beamPlwArc"].double]
beamAreaPipSr  = [beamCorrTable.meta["beamPswSr"].double,
	beamCorrTable.meta["beamPmwSr"].double,
	beamCorrTable.meta["beamPlwSr"].double]

#-------------------------------------------------------------------------------
# Get beam profiles from Calibration tree
beamProfs = calibration.phot.refs["RadialCorrBeam"].product

#############################    TEMPORARY FIX    ##############################
# Add in core and constant tables for maxrad 700 (***TEMPORARY FIX***)
beamProfs["constant"] = asciiTableReader(file=os.path.join(workDir,'beamProfs_constant_r700.csv'))
beamProfs["core"] = asciiTableReader(file=os.path.join(workDir,'beamProfs_core_r700.csv'))

# Re-compute normalised beam areas for new beams
beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
nRad=len(beamRad)
beamNormArea=TableDataset(description="Area as a function of radius, normalised by final value")
beamNormArea.addColumn("radius",Column(Float1d(beamRad)))

for array in ["PSW","PMW","PLW"]:
	#add column to table
	beamNormArea.addColumn(array,Column(Float1d(nRad)))
	#get core beam
	beamComb=beamProfs.getCoreCorrectionTable().getColumn(array).data
	#get const beam
	beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data
	#work out where constant beam applies
	isConst = beamComb.where(beamComb < beamConst)
	#apply constant beam where applicable
	beamComb[isConst] = beamConst[isConst]
	beamInterp=CubicSplineInterpolator(beamRad,beamComb*2.*Math.PI*beamRad)
	integrator=TrapezoidalIntegrator(0,max(beamRad))
	beamTotArea=integrator.integrate(beamInterp)
	for r in range(nRad):
		integrator=TrapezoidalIntegrator(0,beamRad[r])
		beamNormArea[array].data[r]=integrator.integrate(beamInterp)/beamTotArea

beamProfs["normArea"]=beamNormArea
################################################################################

#-------------------------------------------------------------------------------
# Get effective frequencies from Calibration tree and convert to Hz
spireEffFreq = [beamProfs.meta["freqEffPsw"].double*1.e9,
	beamProfs.meta["freqEffPmw"].double*1.e9,
	beamProfs.meta["freqEffPlw"].double*1.e9]

#####    TEMPORARY FIX    #####
# Correct effective frequencies (***TEMPORARY FIX***)
spireEffFreq=[1217.1400e9,865.4342e9,610.1790e9] #for maxRad=700
###############################

#-------------------------------------------------------------------------------
# Beam areas at effective frequency in Sr
# Due to beam model, this is the same as the value measured on Neptune.
# These values are online at https://nhscsci.ipac.caltech.edu/sc/index.php/Spire/PhotBeamProfileDataAndAnalysis
arcsec2Sr = (Math.PI/(60.*60.*180))**2
spireAreaEffFreq = [450.*arcsec2Sr, 795.*arcsec2Sr, 1665.*arcsec2Sr]
# Neptune spectral index (from ESA4 model)
alphaNep=[1.29,1.42,1.47]

# Get RSRF 
# Interpolate to common frequency grid
interpPLW = LinearInterpolator(spireFreq, rsrf['rsrf']['plw'].data)
interpPMW = LinearInterpolator(spireFreq, rsrf['rsrf']['pmw'].data)
interpPSW = LinearInterpolator(spireFreq, rsrf['rsrf']['psw'].data)

spireFiltPLW = Double1d(nNu)
spireFiltPMW = Double1d(nNu)
spireFiltPSW = Double1d(nNu)

spireFiltPLW[ix] = interpPLW(freq[ix])
spireFiltPMW[ix] = interpPMW(freq[ix])
spireFiltPSW[ix] = interpPSW(freq[ix])

# Read aperture efficiency table from calibration tree
SpireApEff = calibration.refs["Phot"].product.refs["ApertureEfficiency"].product["frequency"]
spireApEffFreq = SpireApEff["frequency"].data * 1e9 #comes in [GHz]
spireApEffPsw  = SpireApEff["PSW"].data
spireApEffPmw  = SpireApEff["PMW"].data
spireApEffPlw  = SpireApEff["PLW"].data

# Fold in and interpolate aperture efficiency
ix = freq.where((freq>=MIN(spireApEffFreq)) & (freq<=MAX(spireApEffFreq)))
interpPLW = LinearInterpolator(spireApEffFreq, spireApEffPlw)
interpPMW = LinearInterpolator(spireApEffFreq, spireApEffPmw)
interpPSW = LinearInterpolator(spireApEffFreq, spireApEffPsw)

# Keep copy of RSRF without aperture efficiency
spireFiltOnlyPSW=spireFiltPSW.copy()
spireFiltOnlyPMW=spireFiltPMW.copy()
spireFiltOnlyPLW=spireFiltPLW.copy()

# Also for point source filter profiles
spireFiltPLW[ix] = interpPLW(freq[ix]) * spireFiltPLW[ix]
spireFiltPMW[ix] = interpPMW(freq[ix]) * spireFiltPMW[ix]
spireFiltPSW[ix] = interpPSW(freq[ix]) * spireFiltPSW[ix]

#for f in range(len(freq)):
#	if f % 1000 == 0:
#		print freq[f]/1.e12, spireFiltOnlyPSW[f],interpPSW(freq[f]),spireFiltPSW[f]

#-------------------------------------------------------------------------------
#===============================================================================
#-------------------------------------------------------------------------------

# Set up functions to calculate beam profile, effective frequency and effective area


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
	  
	2013/12/19  C. North  initial version

	"""

	import herschel.ia.toolbox.spectrum.utils.integration.TrapezoidIntegrator as TrapezoidIntegrator

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

#	if verbose:
#		for f in range(len(freq)):
#			if (f<10) or (f % 1000 ==0):
#				intf=TrapezoidalIntegrator(minFreq,freq[f])
#				print f, freq[f]/1.e12,transm[f],fSky[f],monoArea[f],numInterp(freq[f]),denomInterp(freq[f]),intf.integrate(numInterp),intf.integrate(denomInterp)

	return(effArea)

def spireMonoAreas(freq,beamProfs,effFreq,gamma,array,freqFact=100):

	"""
	========================================================================
	Generates array of monochromatic beam areas over frequency range by
	calculating over a sparser array and interpolating

	Inputs:
	  freq:      (array float) frequency vector [Hz] for which monochromatic
	               beams areas should be calculated
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

def spireMonoAreasSimple(freq,effFreq,areaEffFreq,gamma):

	"""
	========================================================================
	Generates array of monochromatic beam areas over frequency range by
	calculating over a sparser array and interpolating

	Inputs:
	  freq:        (array float) frequency vector [Hz] for which monochromatic
	                 beams areas should be calculated
	  effFreq:     (float) effective frequency [Hz] of array
	  areaEffFreq: (float) beam area at effective frequency 
	  gamma:       (float) Exponent of powerlaw describing FWHM dependence
	                 on frequency
	Outputs:     
	             (array float) Monochromatic Beam area at frequencies
	                corresponding to freq [same units as areaEffFreq]

	Calculation:
	  Uses simple power law to scale beam area with frequency
	  
	2014/01/07  C. North  initial version

	"""

	beamMonoArea = areaEffFreq * (freq/effFreq)**(2.*gamma)
	return(beamMonoArea)

def spireFindEffFreq(freq,rsrf,beamProfs,effFreqInit,gamma,areaNep,alphaNep,array,simpleBeam=False,freqFact=500,initRange=0.01,reqPrec=1.e-6,maxIter=5,verbose=False):
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

	2014/01/07  C. North  initial version
	"""
		
	#calculate initial estimates of effective Frequency
	#parameterised be offset from original estimate
	relEff=Double1d([1.-initRange,1.+initRange])
	effFreqs=effFreqInit*relEff
	if verbose:
		print 'Initial Effective Frequencies: [%.2f : %.2f] GHz'%(effFreqs[0]/1.e9,effFreqs[1]/1.e9)
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
			print 'iter %d: %.4f [effFreq %.2f GHz], Area=%.2f, RelDiff=%.4g'%(iter,relEff[1],effFreqs[1]/1.e9,effAreas[1]/arcsec2Sr,relAreas[1])
		if (Math.abs(relAreas[1]) < reqPrec):
			done=True
	if ((iter > maxIter) and (done==False)):
		print "Warning: maximum iterations [%d] exceeded without conversion [%g]"%(maxIter,reqPrec)

	if verbose:
		print 'Final effFreq: %.4f'%(effFreqs[1]/1.e9)
		print 'Resulting Neptune area: %.2f [rel Diff: %.3g]'%(effAreas[1]/arcsec2Sr,relAreas[1])
	return(effFreqs[1])

#-------------------------------------------------------------------------------
#===============================================================================
#-------------------------------------------------------------------------------

##Find spire Effective frequencies
if calcEffFreq:
	print 'Calculating effective frequencies...'
	spireEffFreq[0] = spireFindEffFreq(freq,spireFiltOnlyPSW,beamProfs,spireEffFreq[0],gamma,spireAreaEffFreq[0],alphaNep[0],'PSW',simpleBeam=simpleBeam,verbose=False)
	spireEffFreq[1] = spireFindEffFreq(freq,spireFiltOnlyPMW,beamProfs,spireEffFreq[1],gamma,spireAreaEffFreq[1],alphaNep[1],'PMW',simpleBeam=simpleBeam,verbose=False)
	spireEffFreq[2] = spireFindEffFreq(freq,spireFiltOnlyPLW,beamProfs,spireEffFreq[2],gamma,spireAreaEffFreq[2],alphaNep[2],'PLW',simpleBeam=simpleBeam,verbose=False)

#calculate monochromatic beam areas using full beam treatment
if simpleBeam:
	print 'Calculating monochromatic beam areas for simple beam...'
	beamMonoAreaPsw = spireMonoAreasSimple(freq,spireEffFreq[0],spireAreaEffFreq[0],gamma)
	beamMonoAreaPmw = spireMonoAreasSimple(freq,spireEffFreq[1],spireAreaEffFreq[1],gamma)
	beamMonoAreaPlw = spireMonoAreasSimple(freq,spireEffFreq[2],spireAreaEffFreq[2],gamma)
else:
	print 'Calculating monochromatic beam areas for simple beam...'
	beamMonoAreaPsw = spireMonoAreas(freq,beamProfs,spireEffFreq[0],gamma,'PSW')
	beamMonoAreaPmw = spireMonoAreas(freq,beamProfs,spireEffFreq[1],gamma,'PMW')
	beamMonoAreaPlw = spireMonoAreas(freq,beamProfs,spireEffFreq[2],gamma,'PLW')

#-------------------------------------------------------------------------------
# Optional check of effective beam area
#  calculate effective beam area for Neptune (alphaNep)
#  should be equal to measured area spireAreaEffFreq
#effNepArea=Double1d(3)
#effNepArea[0]=spireEffArea(freq, spireFiltOnlyPSW, beamMonoAreaPsw, BB=False, alpha=alphaNep[0])
#effNepArea[1]=spireEffArea(freq, spireFiltOnlyPMW, beamMonoAreaPmw, BB=False, alpha=alphaNep[1])
#effNepArea[2]=spireEffArea(freq, spireFiltOnlyPLW, beamMonoAreaPlw, BB=False, alpha=alphaNep[2])
#print 'Measured area: [%.2f, %.2f, %.2f] arcsec^2'%
#	(spireAreaEffFreq[0]/arcsec2Sr,spireAreaEffFreq[1]/arcsec2Sr,spireAreaEffFreq[2]/arcsec2Sr)
#print 'Effective area: [%.2f, %.2f, %.2f] arcsec^2'%
#	(effNepArea[0]/arcsec2Sr,effNepArea[1]/arcsec2Sr,effNepArea[2]/arcsec2Sr)
#print 'Relative difference: [%.2f, %.2f, %.2f] %%'%
#	(100.*(1.-effNepArea[0]/spireAreaEffFreq[0]), 100.*(1.-effNepArea[1]/spireAreaEffFreq[1]), 100.*(1.-effNepArea[02/spireAreaEffFreq[2]))

#-------------------------------------------------------------------------------
# Load HFI filter functions
hfiRIMO = fitsReader(file = hfiRIMOFile)

# Read and parse HFI filter files into tables 
hfiFreq545 = 100. * c * hfiRIMO['BANDPASS_F545']['WAVENUMBER'].data[1:-1]
hfiFreq857 = 100. * c * hfiRIMO['BANDPASS_F857']['WAVENUMBER'].data[1:-1]

hfiTrans545 = Double1d(hfiRIMO['BANDPASS_F545']['TRANSMISSION'].data[1:-1])
hfiTrans857 = Double1d(hfiRIMO['BANDPASS_F857']['TRANSMISSION'].data[1:-1])


# Excluding data points that are not monotonically increasing for HFI-545
diff = Double1d(len(hfiFreq545)-1)
for i in range(len(hfiFreq545)-1):
	diff[i] = hfiFreq545[i+1]-hfiFreq545[i]

ind = diff.where(diff <= 0.).toInt1d()
for i in ind:
	hfiFreq545.delete(i+1,1)
	hfiTrans545.delete(i+1,1)


# Excluding data points that are not monotonically increasing for HFI-545
diff = Double1d(len(hfiFreq857)-1)
for i in range(len(hfiFreq857)-1):
	diff[i] = hfiFreq857[i+1] - hfiFreq857[i]

ind = diff.where(diff <= 0.).toInt1d()
for i in ind:
	hfiFreq857.delete(i+1,1)
	hfiTrans857.delete(i+1,1)


# Interpolate HFI transmissions to common frequency grid 
hfiFilt545 = Double1d(nNu)
hfiFilt857 = Double1d(nNu)

ix545     = freq.where((freq>=MIN(hfiFreq545)) & (freq<=MAX(hfiFreq545)))
interp545 = CubicSplineInterpolator(hfiFreq545, hfiTrans545)
hfiFilt545[ix545] = interp545(freq[ix545])

ix857     = freq.where((freq>=MIN(hfiFreq857)) & (freq<=MAX(hfiFreq857)))
interp857 = CubicSplineInterpolator(hfiFreq857, hfiTrans857)
hfiFilt857[ix857] = interp857(freq[ix857])

# Normalize transmissions
hfiFilt545 = hfiFilt545 / MAX(hfiFilt545)
hfiFilt857 = hfiFilt857 / MAX(hfiFilt857)

#-------------------------------------------------------------------------------
# Plot relative spectral response functions

if plot:
	# Make linear plot of all filter bands including those for point sources
	p = PlotXY()
	p.addLayer(LayerXY(c/freq*1e6, hfiFilt857, name='hfiFilt857'))
	p.addLayer(LayerXY(c/freq*1e6, hfiFilt545, name='hfiFilt545'))
	#
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPSW * (freq/spireRefFreq[0])**(2*gamma), name='spireXFiltPSW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPMW * (freq/spireRefFreq[1])**(2*gamma), name='spireXFiltPMW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPLW * (freq/spireRefFreq[2])**(2*gamma), name='spireXFiltPLW'))
	#
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPSW, name='spireFiltPSW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPMW, name='spireFiltPMW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPLW, name='spireFiltPLW'))
	#
	p.xaxis.range = [100,800]
	p.yaxis.range = [-0.1,1.1]
	p.xaxis.titleText = "Wavelength [micron]"
	p.yaxis.titleText = "Relative Spectral Response"
	p.legend.visible = 1
	#
	# Make log plot in y direction of all extended source filter bands
	p = PlotXY()
	p.addLayer(LayerXY(c/freq*1e6, hfiFilt857, name='hfiFilt857'))
	p.addLayer(LayerXY(c/freq*1e6, hfiFilt545, name='hfiFilt545'))
	#
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPSW * (freq/spireRefFreq[0])**(2*gamma), name='spireXFiltPSW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPMW * (freq/spireRefFreq[1])**(2*gamma), name='spireXFiltPMW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPLW * (freq/spireRefFreq[2])**(2*gamma), name='spireXFiltPLW'))
	#
	p.xaxis.range = [0,2000]
	p.yaxis.range = [1e-8,1.1]
	p.xaxis.titleText = "Wavelength [micron]"
	p.yaxis.titleText = "Relative Spectral Response"
	p.yaxis.type = Axis.LOG
	p.legend.visible = 1


# hpXcalKcorr function definition
def hpXcalKcorr(freq0, freq, transm, BB=True, temp=20.0, beta=1.8, alpha=-1.0, ext=False, monoArea=None):
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
	  intMethod: (string) Integration method to use.
	                [Trapezoidal | Spectrum | Tabulated | Simpson | Romberg | Rectangular]
	                OPTIONAL. Default='Trapezoidal'	
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

	  N.B. If ext=True, and monoBeam is in units sr, then this procedure outputs
	    K-correction factor in [Jy/sr per Jy/beam] and Sky emission in Jy/sr.
	    Units will change if a different input unit is used.
	
	
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
	# Return the result as a 2-element list of K-correction and flux at freq0
	return (kWave, fSky0)

#-------------------------------------------------------------------------------
# Calculate K-correction factors for extended source assuming alpha=-1
# Compute the total conversion, then divide by the effective beam area
# SPIRE PSW (alpha=-1)
k4E_PSW_Tot = hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, False, ext=True, monoArea=beamMonoAreaPsw)[0]
beamAreaPipSr[0] = spireEffArea(freq, spireFiltOnlyPSW, beamMonoAreaPsw, BB=False, alpha=-1)
k4E_PSW = k4E_PSW_Tot * beamAreaPipSr[0]

# SPIRE PMW (alpha=-1)
k4E_PMW_Tot = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, False, ext=True, monoArea=beamMonoAreaPmw)[0]
beamAreaPipSr[1] = spireEffArea(freq, spireFiltOnlyPMW, beamMonoAreaPmw, BB=False, alpha=-1)
k4E_PMW = k4E_PMW_Tot * beamAreaPipSr[1]

# SPIRE PLW (alpha=-1)
k4E_PLW_Tot = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, False, ext=True, monoArea=beamMonoAreaPlw)[0]
beamAreaPipSr[2] = spireEffArea(freq, spireFiltOnlyPLW, beamMonoAreaPlw, BB=False, alpha=-1)
k4E_PLW = k4E_PLW_Tot * beamAreaPipSr[2]

# HFI bands (alpha=-1)
k4E_857 = hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857,   False)[0]
k4E_545 = hpXcalKcorr(hfiRefFreq[1],   freq, hfiFilt545,   False)[0]

# Calculate SPIRE K-correction factors for point source assuming alpha=-1
k4P_PSW = hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, False, ext=False)[0]
k4P_PMW = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, False, ext=False)[0]
k4P_PLW = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, False, ext=False)[0]

# Pipeline conversion (point source to extended source)
kPtoE_PSW = k4E_PSW_Tot/k4P_PSW
kPtoE_PMW = k4E_PMW_Tot/k4P_PMW
kPtoE_PLW = k4E_PLW_Tot/k4P_PLW

# Print Spire K4 factors for point source calibration for verification
print 'SPIRE k4 factors for point source: %6.4f, %6.4f, %6.4f' % (k4P_PSW,k4P_PMW,k4P_PLW)

# Print Spire K4 factors for extended source calibration for verification
print 'SPIRE k4 (Total) factors for extended source: %6.4f, %6.4f, %6.4f' % (k4E_PSW_Tot,k4E_PMW_Tot,k4E_PLW_Tot)
print 'SPIRE k4 factors for extended source: %6.4f, %6.4f, %6.4f' % (k4E_PSW,k4E_PMW,k4E_PLW)

print 'SPIRE kPtoE conversion factor: %6.4f, %6.4f, %6.4f' % (kPtoE_PSW,kPtoE_PMW,kPtoE_PLW)

#-------------------------------------------------------------------------------
# Compute colour corrections for range of alpha and temp/beta

ColorCorrKList=herschel.spire.ia.cal.PhotColorCorrKList(description='Set of source-type correction products')
ColorCorrKPsrc=herschel.spire.ia.dataset.PhotColorCorrK(description='Spire Colour Corrections product (Point source)')
ColorCorrKExtd=herschel.spire.ia.dataset.PhotColorCorrK(description='Spire Colour Corrections product (Point source)')
#ColorCorrApertureList=herschel.spire.ia.cal.PhotColorCorrApertureList()
ColorCorrBeam=herschel.spire.ia.dataset.PhotColorCorrBeam(description='Spire Beam Corrections with Spectral Index and Temperature product')

#-------------------------------------------------------------------------------
doKcorr=True
if doKcorr:
	# Compute colour corrections for alpha
	alpha=[-4.,-3.5,-3.,-2.5,-2.,-1.5,-1.,-0.5,0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.]
	
	nA=len(alpha)
	
	k4Pa=TableDataset(description='Point Source Conversion (Spectral Index)')
	k4Pa.addColumn('alpha',Column(Float1d(alpha)))
	k4Pa.addColumn('PSW',Column(Float1d(nA)))
	k4Pa.addColumn('PMW',Column(Float1d(nA)))
	k4Pa.addColumn('PLW',Column(Float1d(nA)))
	
	kcPa=k4Pa.copy()
	kcPa.setDescription('Point Source Colour Correction (Spectral Index)')
	
	k4Ea=k4Pa.copy()
	k4Ea.setDescription('Extended Source Conversion (Spectral Index)')
	k4Ea.meta['effFreqPsw']=DoubleParameter(spireEffFreq[0],unit=Frequency.HERTZ,description='Effective frequency for PSW')
	k4Ea.meta['effFreqPmw']=DoubleParameter(spireEffFreq[1],unit=Frequency.HERTZ,description='Effective frequency for PSW')
	k4Ea.meta['effFreqPlw']=DoubleParameter(spireEffFreq[2],unit=Frequency.HERTZ,description='Effective frequency for PSW')
	
	k4EaTot=k4Ea.copy()
	k4EaTot.setDescription('Total Extended Source Conversion (Spectral Index)')
	
	kcEa=k4Ea.copy()
	kcEa.setDescription('Point Source Colour Correction (Spectral Index)')
	
	effBeamaSr=TableDataset(description='Beam Solid Angle (Spectral Index)')
	effBeamaSr.addColumn('alpha',Column(Float1d(alpha)))
	effBeamaSr.addColumn('PSW',Column(Float1d(nA),unit=SolidAngle.STERADIANS,description=''))
	effBeamaSr.addColumn('PMW',Column(Float1d(nA),unit=SolidAngle.STERADIANS,description=''))
	effBeamaSr.addColumn('PLW',Column(Float1d(nA),unit=SolidAngle.STERADIANS,description=''))
	
	effBeamaSr2=TableDataset(description='Beam Solid Angle (Spectral Index)')
	effBeamaSr2.addColumn('alpha',Column(Float1d(alpha)))
	effBeamaSr2.addColumn('PSW',Column(Float1d(nA),unit=SolidAngle.STERADIANS,description=''))
	effBeamaSr2.addColumn('PMW',Column(Float1d(nA),unit=SolidAngle.STERADIANS,description=''))
	effBeamaSr2.addColumn('PLW',Column(Float1d(nA),unit=SolidAngle.STERADIANS,description=''))
	
	kBeama=TableDataset(description='Beam Solid Angle Colour Correction (Spectral Index)')
	kBeama.addColumn('alpha',Column(Float1d(alpha)))
	kBeama.addColumn('PSW',Column(Float1d(nA),unit=Scalar.ONE,description=''))
	kBeama.addColumn('PMW',Column(Float1d(nA),unit=Scalar.ONE,description=''))
	kBeama.addColumn('PLW',Column(Float1d(nA),unit=Scalar.ONE,description=''))
	
	for a in range(nA):
		#Beam from calibration file
		#point source conversion
		k4Pa['PSW'].data[a]=hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, False, alpha=alpha[a])[0]
		k4Pa['PMW'].data[a]=hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, False, alpha=alpha[a])[0]
		k4Pa['PLW'].data[a]=hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, False, alpha=alpha[a])[0]
	
		#point source colour correction
		kcPa['PSW'].data[a] = k4Pa['PSW'].data[a]/k4P_PSW
		kcPa['PMW'].data[a] = k4Pa['PMW'].data[a]/k4P_PMW
		kcPa['PLW'].data[a] = k4Pa['PLW'].data[a]/k4P_PLW
	
		#total extended source conversion
		k4EaTot['PSW'].data[a]=hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, False, alpha=alpha[a], ext=True, monoArea=beamMonoAreaPsw)[0]
		k4EaTot['PMW'].data[a]=hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, False, alpha=alpha[a], ext=True, monoArea=beamMonoAreaPmw)[0]
		k4EaTot['PLW'].data[a]=hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, False, alpha=alpha[a], ext=True, monoArea=beamMonoAreaPlw)[0]
	
		#effective area
		effBeamaSr['PSW'].data[a]=spireEffArea(freq, spireFiltOnlyPSW, beamMonoAreaPsw, BB=False, alpha=alpha[a])
		effBeamaSr['PMW'].data[a]=spireEffArea(freq, spireFiltOnlyPMW, beamMonoAreaPmw, BB=False, alpha=alpha[a])
		effBeamaSr['PLW'].data[a]=spireEffArea(freq, spireFiltOnlyPLW, beamMonoAreaPlw, BB=False, alpha=alpha[a])
	
		#beam correction factor
		kBeama['PSW'].data[a]=effBeamaSr['PSW'].data[a] / beamAreaPipSr[0]
		kBeama['PMW'].data[a]=effBeamaSr['PMW'].data[a] / beamAreaPipSr[1]
		kBeama['PLW'].data[a]=effBeamaSr['PLW'].data[a] / beamAreaPipSr[2]
	
		#extended source conversion
		k4Ea['PSW'].data[a]=k4EaTot['PSW'].data[a] * effBeamaSr['PSW'].data[a]
		k4Ea['PMW'].data[a]=k4EaTot['PMW'].data[a] * effBeamaSr['PMW'].data[a]
		k4Ea['PLW'].data[a]=k4EaTot['PLW'].data[a] * effBeamaSr['PLW'].data[a]
	
		#extended source conversion
		kcEa['PSW'].data[a]=k4EaTot['PSW'].data[a] / k4E_PSW_Tot
		kcEa['PMW'].data[a]=k4EaTot['PMW'].data[a] / k4E_PMW_Tot
		kcEa['PLW'].data[a]=k4EaTot['PLW'].data[a] / k4E_PLW_Tot
	
	asciiTableWriter(table=kBeama,file=os.path.join(workDir,'kBeam_alpha.csv'))
	asciiTableWriter(table=kcPa,file=os.path.join(workDir,'kcP_alpha.csv'))
	asciiTableWriter(table=kcEa,file=os.path.join(workDir,'kcE_alpha.csv'))
	
	
	chkPrint=True
	if chkPrint:
		axPrint=[0,6,12,18] #indices of alpha to print to terminal
		nAp=len(axPrint)
		print 'KmonP'
		for ax in range(nAp):
			row=k4Pa.getRow(axPrint[ax])
			print "%.1f   %.5g   %.5g   %.5g"%(row[0],row[1],row[2],row[3])
	
		print 'KmonE total'
		for ax in range(nAp):
			row=k4EaTot.getRow(axPrint[ax])
			print "%.1f   %.9g   %.9g   %.9g"%(row[0],row[1],row[2],row[3])
	
		print'KmonE'
		for ax in range(nAp):
			row=k4Ea.getRow(axPrint[ax])
			print "%.1f   %.9f   %.9f   %.9f"%(row[0],row[1],row[2],row[3])
	
		print 'effArea'
		for ax in range(nAp):
			row=effBeamaSr.getRow(axPrint[ax])
			print "%.1f   %.9f   %.9f   %.9f"%(row[0],row[1]/arcsec2Sr,row[2]/arcsec2Sr,row[3]/arcsec2Sr)
	
		print 'MonoArea'
		for n in (0,nNu/4,nNu/2,3*nNu/4,nNu-1):
			print "%.1f   %.2f   %.2f   %.2f"%(freq[n]/1.e9,beamMonoAreaPsw[n]/arcsec2Sr,beamMonoAreaPmw[n]/arcsec2Sr,beamMonoAreaPlw[n]/arcsec2Sr)


	#-----------------------------------------------------------------------
	# compute colour corrections for betas
	beta=[0.,0.5,1.,1.25,1.5,1.75,2.,2.5,3.]
	temp=range(3,41)
	
	nB=len(beta)
	nT=len(temp)
	
	k4Pb=TableDataset(description='Point Source Conversion (Modified Black Body)')
	k4Pb.addColumn('temp',Column(Float1d(temp),unit=Temperature.KELVIN,description=''))
	k4Pb.addColumn('PSW',Column(Float1d(nT)))
	k4Pb.addColumn('PMW',Column(Float1d(nT)))
	k4Pb.addColumn('PLW',Column(Float1d(nT)))
	
	kcPb=k4Pb.copy()
	kcPb.setDescription('Point Source Colour Correction (Modified Black Body)')
	
	k4Eb=k4Pb.copy()
	k4Eb.setDescription('Extended Source Conversion (Modified Black Body)')
	#k4Eb.meta['effFreqPsw']=DoubleParameter(spireEffFreq[0],unit=Frequency.HERTZ,description='Effective frequency for PSW')
	#k4Eb.meta['effFreqPmw']=DoubleParameter(spireEffFreq[1],unit=Frequency.HERTZ,description='Effective frequency for PSW')
	#k4Eb.meta['effFreqPlw']=DoubleParameter(spireEffFreq[2],unit=Frequency.HERTZ,description='Effective frequency for PSW')
	
	k4EbTot=k4Eb.copy()
	k4EbTot.setDescription('Total Extended Source Conversion (Modified Black Body)')
	
	kcEb=k4Eb.copy()
	kcEb.setDescription('Point Source Colour Correction (Modified Black Body)')
	
	effBeambSr=TableDataset(description='Beam Solid Angle (Modified Black Body)')
	effBeambSr.addColumn('temp',Column(Float1d(temp),unit=Temperature.KELVIN,description=''))
	effBeambSr.addColumn('PSW',Column(Float1d(nT),unit=SolidAngle.STERADIANS,description=''))
	effBeambSr.addColumn('PMW',Column(Float1d(nT),unit=SolidAngle.STERADIANS,description=''))
	effBeambSr.addColumn('PLW',Column(Float1d(nT),unit=SolidAngle.STERADIANS,description=''))
	
	kBeamb=TableDataset(description='Beam Solid Angle Colour Correction (Modified Black Body)')
	kBeamb.addColumn('temp',Column(Float1d(temp),unit=Temperature.KELVIN,description=''))
	kBeamb.addColumn('PSW',Column(Float1d(nT),unit=Scalar.ONE,description=''))
	kBeamb.addColumn('PMW',Column(Float1d(nT),unit=Scalar.ONE,description=''))
	kBeamb.addColumn('PLW',Column(Float1d(nT),unit=Scalar.ONE,description=''))
	
	
	for b in range(len(beta)):
		for t in range(nT):
			#Beam from calibration file
			#point source conversion
			k4Pb['PSW'].data[t]=hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, BB=True, beta=beta[b], temp=temp[t])[0]
			k4Pb['PMW'].data[t]=hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, BB=True, beta=beta[b], temp=temp[t])[0]
			k4Pb['PLW'].data[t]=hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, BB=True, beta=beta[b], temp=temp[t])[0]
		
			#point source colour correction
			kcPb['PSW'].data[t] = k4Pb['PSW'].data[t]/k4P_PSW
			kcPb['PMW'].data[t] = k4Pb['PMW'].data[t]/k4P_PMW
			kcPb['PLW'].data[t] = k4Pb['PLW'].data[t]/k4P_PLW
		
			#total extended source conversion
			k4EbTot['PSW'].data[t]=hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, BB=True, beta=beta[b], temp=temp[t], ext=True, monoArea=beamMonoAreaPsw)[0]
			k4EbTot['PMW'].data[t]=hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, BB=True, beta=beta[b], temp=temp[t], ext=True, monoArea=beamMonoAreaPmw)[0]
			k4EbTot['PLW'].data[t]=hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, BB=True, beta=beta[b], temp=temp[t], ext=True, monoArea=beamMonoAreaPlw)[0]
		
			#effective area
			effBeambSr['PSW'].data[t]=spireEffArea(freq, spireFiltOnlyPSW, beamMonoAreaPsw, BB=True, beta=beta[b], temp=temp[t])
			effBeambSr['PMW'].data[t]=spireEffArea(freq, spireFiltOnlyPMW, beamMonoAreaPmw, BB=True, beta=beta[b], temp=temp[t])
			effBeambSr['PLW'].data[t]=spireEffArea(freq, spireFiltOnlyPLW, beamMonoAreaPlw, BB=True, beta=beta[b], temp=temp[t])
		
			kBeamb['PSW'].data[t]=effBeambSr['PSW'].data[t] / beamAreaPipSr[0]
			kBeamb['PMW'].data[t]=effBeambSr['PMW'].data[t] / beamAreaPipSr[1]
			kBeamb['PLW'].data[t]=effBeambSr['PLW'].data[t] / beamAreaPipSr[2]
		
			#extended source conversion
			k4Eb['PSW'].data[t]=k4EbTot['PSW'].data[t] * effBeambSr['PSW'].data[t]
			k4Eb['PMW'].data[t]=k4EbTot['PMW'].data[t] * effBeambSr['PMW'].data[t]
			k4Eb['PLW'].data[t]=k4EbTot['PLW'].data[t] * effBeambSr['PLW'].data[t]
		
			#extended source conversion
			kcEb['PSW'].data[t]=k4EbTot['PSW'].data[t] / k4E_PSW_Tot
			kcEb['PMW'].data[t]=k4EbTot['PMW'].data[t] / k4E_PMW_Tot
			kcEb['PLW'].data[t]=k4EbTot['PLW'].data[t] / k4E_PLW_Tot
		
		#produce string for beta filename suffix
		betaTxt='%d_%d'%(int(beta[b]),int(100*(beta[b])%1))
	
		#write to ascii files	
		asciiTableWriter(table=kBeamb,file=os.path.join(workDir,'kBeam_beta_%s.csv'%betaTxt))
		asciiTableWriter(table=kcPb,file=os.path.join(workDir,'kcP_beta_%s.csv'%betaTxt))
		asciiTableWriter(table=kcEb,file=os.path.join(workDir,'kcE_beta_%s.csv'%betaTxt))

def spireEffBeam(freq,transm,beamProfs,effFreq,gamma,array,beamRadMap,BB=False,temp=20.0,beta=1.8,alpha=-1.0):

	"""
	========================================================================
	Computes an effective beam profile for a given source spectrum
	***
	N.B. Computing the beam area this way integrates over frequency *then* radius,
	  while beamMonoAreas->spireEffArea integrates over radius then frequency.
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
	print 'making beam map'
	for x in range(nxMap):
		for y in range(nyMap):
			if beamRadMap['image'].data[x,y] <= maxRad-1:
				effBeamMap['image'].data[x,y]= \
				    effBeamInterp(beamRadMap['image'].data[x,y])

	return(effBeamArea,effBeam,effBeamMap)

doApCorr=True
if doApCorr:
	#set aperture photometry radius and annulus size
	apPhotRad=[22.,30.,45.]
	apPhotBGRad=[60.,90.]
	#read in map from file
	beamNames = {"PSW":"0x5000241aL_PSW_pmcorr_1arcsec_norm_beam.fits",
		"PMW":"0x5000241aL_PMW_pmcorr_1arcsec_norm_beam.fits",
		"PLW":"0x5000241aL_PLW_pmcorr_1arcsec_norm_beam.fits"}
	try:
		beamIn = fitsReader(file = os.path.join(workDir,"0x5000241aL_PSW_pmcorr_1arcsec_norm_beam.fits"))
	except:
		#download if not available
		urllib.urlretrieve ("https://nhscsci.ipac.caltech.edu/spire/data/beam_profiles/"+beamNameNorm,\
		    os.path.join(workDir,beamNames[band]))
		beamIn = fitsReader(file = os.path.join(workDir,beamNames[band]))

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

	apCorrIncBGa=TableDataset(description='Aperture Correction Factor including background annulus, against spectral index')
	apCorrIncBGa.addColumn('alpha',Column(Float1d(alpha)))
	apCorrIncBGa.addColumn('PSW',Column(Float1d(nA),unit=Scalar.ONE,description=''))
	apCorrIncBGa.addColumn('PMW',Column(Float1d(nA),unit=Scalar.ONE,description=''))
	apCorrIncBGa.addColumn('PLW',Column(Float1d(nA),unit=Scalar.ONE,description=''))

	apCorrNoBGa=apCorrIncBGa.copy()
	apCorrNoBGa.setDescription('Aperture Correction Factor not including background, against spectral index')

	for a in range(nA):
		print 'alpha=%.1f'%alpha[a]
		(effBeamAreaPSW,effBeamProfPSW,effBeamMapPSW)=spireEffBeam(freq,spireFiltOnlyPSW,beamProfs,spireEffFreq[0],gamma,'PSW',beamRadMap,BB=False,alpha=alpha[a])
		apPhotPSW = annularSkyAperturePhotometry(image=effBeamMapPSW, fractional=1, centerX=700.0, centerY=700.0, radiusArcsec=apPhotRad[0], innerArcsec=apPhotBGRad[0], outerArcsec=apPhotBGRad[1])
		apCorrIncBGa['PSW'].data[a]=effBeamAreaPSW/apPhotPSW.getTargetPlusSkyTotal()
		apCorrNoBGa['PSW'].data[a]=effBeamAreaPSW/apPhotPSW.getTargetTotal()
		
		(effBeamAreaPMW,effBeamProfPMW,effBeamMapPMW)=spireEffBeam(freq,spireFiltOnlyPMW,beamProfs,spireEffFreq[1],gamma,'PMW',beamRadMap,BB=False,alpha=alpha[a])
		apPhotPMW = annularSkyAperturePhotometry(image=effBeamMapPMW, fractional=1, centerX=700.0, centerY=700.0, radiusArcsec=apPhotRad[1], innerArcsec=apPhotBGRad[0], outerArcsec=apPhotBGRad[1])
		apCorrIncBGa['PMW'].data[a]=effBeamAreaPMW/apPhotPMW.getTargetPlusSkyTotal()
		apCorrNoBGa['PMW'].data[a]=effBeamAreaPMW/apPhotPMW.getTargetTotal()
		
		(effBeamAreaPLW,effBeamProfPLW,effBeamMapPLW)=spireEffBeam(freq,spireFiltOnlyPLW,beamProfs,spireEffFreq[2],gamma,'PLW',beamRadMap,BB=False,alpha=alpha[a])
		apPhotPSW = annularSkyAperturePhotometry(image=effBeamMapPLW, fractional=1, centerX=700.0, centerY=700.0, radiusArcsec=apPhotRad[2], innerArcsec=apPhotBGRad[0], outerArcsec=apPhotBGRad[1])
		apCorrIncBGa['PLW'].data[a]=effBeamAreaPLW/apPhotPLW.getTargetPlusSkyTotal()
		apCorrNoBGa['PLW'].data[a]=effBeamAreaPLW/apPhotPLW.getTargetTotal()

sastdast

#-------------------------------------------------------------------------------
# Calculate and tabulate colour correction parameters for a range of temperatures
k545toPLW = Double1d()				# K-correction from HFI-545 to PLW
k857toPMW = Double1d()				# K-correction from HFI-857 to PMW
k857toPSW = Double1d()				# K-correction from HFI-857 to PSW
ratio545over857 = Double1d()		# Ratio 545GHz to 857GHz filter

if fixBeta == True:
	#
	# Make temperature vector with log distances between tiers
	tvect = Double1d(range(ndiv)) * (Tmax - Tmin) / (ndiv-1) + Tmin
	#
	for temp in tvect:
		#
		# K-correction from HFI-545 to PLW
		k545 = hpXcalKcorr(hfiRefFreq[1],   freq, hfiFilt545,   True, temp, beta0)
		kPLW = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, True, temp, beta0, gamma=gamma)
		k545toPLW.append(k545[0] / kPLW[0] * k4E_PLW / k4E_545 * kPLW[1] / k545[1])
		#
		# K-correction from HFI-857 to PMW
		k857 = hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857,   True, temp, beta0)
		kPMW = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, True, temp, beta0, gamma=gamma)
		k857toPMW.append(k857[0] / kPMW[0] * k4E_PMW / k4E_857 * kPMW[1] / k857[1])
		#
		# K-correction from HFI-857 to PSW
		kPSW = hpXcalKcorr(spireRefFreq[0], freq,  spireFiltPSW, True, temp, beta0, gamma=gamma)
		k857toPSW.append(k857[0] / kPSW[0] * k4E_PSW / k4E_857 * kPSW[1] / k857[1])
		#
		# Ratio of 545 and 845 GHz filters
		ratio545over857.append(k545[1] / k857[1] * k4E_545 / k4E_857 * k857[0] / k545[0])
	#
else:
	#
	# Make temperature vector with log distances between tiers
	bvect = Double1d(range(ndiv)) * (Bmax - Bmin) / (ndiv-1) + Bmin
	#
	for beta in bvect:
		#
		# K-correction from HFI-545 to PLW
		k545 = hpXcalKcorr(hfiRefFreq[1],   freq, hfiFilt545,    True, temp0, beta, gamma=gamma)
		kPLW = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, True, temp0, beta, gamma=gamma)
		k545toPLW.append(k545[0] / kPLW[0] * k4E_PLW / k4E_545 * kPLW[1] / k545[1])
		#
		# K-correction from HFI-857 to PMW
		k857 = hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857,    True, temp0, beta, gamma=gamma)
		kPMW = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, True, temp0, beta, gamma=gamma)
		k857toPMW.append(k857[0] / kPMW[0] * k4E_PMW / k4E_857 * kPMW[1] / k857[1])
		#
		# K-correction from HFI-857 to PSW
		kPSW = hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, True, temp0, beta, gamma=gamma)
		k857toPSW.append(k857[0] / kPSW[0] * k4E_PSW / k4E_857 * kPSW[1] / k857[1])
		#
		# Ratio of 545 and 845 GHz filters
		ratio545over857.append(k545[1] / k857[1] * k4E_545 / k4E_857 * k857[0] / k545[0])
	#

#-------------------------------------------------------------------------------

if plot:
	# Plot colour correction factors
	if fixBeta == True:
		# First plot of parameters against dust temperature of assumed background
		p = PlotXY()
		p.addLayer(LayerXY(tvect,k545toPLW,name='k545toPLW'))
		p.addLayer(LayerXY(tvect,k857toPMW,name='k857toPMW'))
		p.addLayer(LayerXY(tvect,k857toPSW,name='k857toPSW'))
		p.addLayer(LayerXY(tvect,ratio545over857,name='ratio545over857'))
		p.xaxis.titleText = "T [K]"
		p.yaxis.titleText = "K-Factor"
		p.xaxis.type = Axis.LOG
		p.xaxis.range = [1.,100.]
		p.legend.visible = 1
	#
	else:
		# First plot of parameters against dust spectral index of assumed background
		p = PlotXY()
		p.addLayer(LayerXY(bvect,k545toPLW,name='k545toPLW'))
		p.addLayer(LayerXY(bvect,k857toPMW,name='k857toPMW'))
		p.addLayer(LayerXY(bvect,k857toPSW,name='k857toPSW'))
		p.addLayer(LayerXY(bvect,ratio545over857,name='ratio545over857'))
		p.xaxis.titleText = "Beta"
		p.yaxis.titleText = "K-Factor"
		p.legend.visible = 1
	#
	# Second plot of parameters against the ratio of the fluxes found in the 
	# all-sky maps of HFI
	p = PlotXY()
	p.addLayer(LayerXY(ratio545over857,k545toPLW,name='k545toPLW'))
	p.addLayer(LayerXY(ratio545over857,k857toPMW,name='k857toPMW'))
	p.addLayer(LayerXY(ratio545over857,k857toPSW,name='k857toPSW'))
	p.xaxis.titleText = "ratio 545/857"
	p.yaxis.titleText = "K-Factor"
	p.yaxis.range = [-0.1,4]
	p.legend.visible = 1

#-------------------------------------------------------------------------------
# Create calibration product with tabulated colour correction factors

description = 'Colour correction factors between HFI and SPIRE'
hpKcorrTable = TableDataset(description = description) 
hpKcorrTable.meta["Version"] = DoubleParameter(ver, "Colour correction table version")

# Save standard point source correction factor into meta data
hpKcorrTable.meta["k4P_PSW"] = DoubleParameter(k4P_PSW, "PSW point source K-factor")
hpKcorrTable.meta["k4P_PMW"] = DoubleParameter(k4P_PMW, "PMW point source K-factor")
hpKcorrTable.meta["k4P_PLW"] = DoubleParameter(k4P_PLW, "PLW point source K-factor")

hpKcorrTable.meta["k4E_PSW"] = DoubleParameter(k4E_PSW, "PSW extended source K-factor")
hpKcorrTable.meta["k4E_PMW"] = DoubleParameter(k4E_PMW, "PMW extended source K-factor")
hpKcorrTable.meta["k4E_PLW"] = DoubleParameter(k4E_PLW, "PLW extended source K-factor")

# Save beta dependent colour correction tables into binary table
if fixBeta == True:
	hpKcorrTable.meta["Beta"] = DoubleParameter(beta0, "Modified black-body spectral index used")
	#
	hpKcorrTable['Temperature'] = Column(tvect)      # Temperature of modified BB
	hpKcorrTable['Temperature'].setDescription('Temperature of modified black-body')
	#
# Save temperature dependent colour correction tables into binary table
else:
	hpKcorrTable.meta['Temperature'] = DoubleParameter(temp0, "Modified black-body temperature used")
	hpKcorrTable.meta['Temperature'].unit = Temperature.KELVIN
	#
	hpKcorrTable['Beta'] = Column(bvect)          # Spectral index of modified BB
	hpKcorrTable['Beta'].setDescription('Spectral index of modified black-body')

hpKcorrTable.meta["Gamma"] = DoubleParameter(gamma, "Exponent describing FWHM dependence on frequency used")

# Save colour correction tables into binary table
hpKcorrTable['ratio545_857'] = Column(ratio545over857)  # Ratio 545/857 GHz filter
hpKcorrTable['k545toPLW']    = Column(k545toPLW)      # 545GHz to PLW K-factor
hpKcorrTable['k857toPMW']    = Column(k857toPMW)      # 857GHz to PMW K-factor
hpKcorrTable['k857toPSW']    = Column(k857toPSW)      # 857GHz to PSW K-factor

hpKcorrTable['ratio545_857'].setDescription('Ratio of 545/857 GHz filters')
hpKcorrTable['k545toPLW'].setDescription('545GHz to PLW K-factor')
hpKcorrTable['k857toPMW'].setDescription('857GHz to PMW K-factor')
hpKcorrTable['k857toPSW'].setDescription('857GHz to PSW K-factor')

# Embed table dataset into product
hpKcorr = Product(creator='HIPE hpXcalColourCorr.py', description = description)
hpKcorr.setType('HFI to SPIRE colour correction table')
hpKcorr.setInstrument('Planck HFI')
hpKcorr.setStartDate(FineTime(1620997954000000))
hpKcorr.setEndDate(FineTime(1956528035000000))
hpKcorr.meta["version"] = StringParameter(value=str(ver), description='Colour correction table version')
hpKcorr['HDU_1'] = hpKcorrTable

# Save product into FITS file
simpleFitsWriter(product=hpKcorr, file=outFilename, warn=True)

print "Table written in %s"%(outFilename)

#-------------------------------------------------------------------------------
# End of File
