#
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2013 Herschel Science Ground Segment Consortium
#
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
#
# $Id: makeSCalPhotColorCorrHfi.py,v 1.1 2013/10/28 14:52:42 epoleham Exp $
#
#===============================================================================
# 
#  Edition History
#  Chris North   - 18/Feb/2014 - 1.0: First test version
#===============================================================================

import os
scriptVersionString = "makeSCalPhotRadialCorrBeam.py $Revision: 1.0 $"

#directory = "..//..//..//..//..//..//data//spire//cal//SCal"
#dataDir = "//disks//winchester2//calibration_data//"
#LOCAL VERSION
directory = Configuration.getProperty('var.hcss.workdir')
dataDir = Configuration.getProperty('var.hcss.workdir')

df  = java.text.SimpleDateFormat("yyyy.MM.dd/HH:mm:ss/z")

#-------------------------------------------------------------------------------
# Input parameters

#beam profile version
version = "3.0"

beamNewVersion="1.1"
beamNewFileConstant = 'beamProfs_constant_v%s.csv'%beamNewVersion
beamNewFileCore = 'beamProfs_core_v%s.csv'%beamNewVersion

formatVersion = "1.0"

verbose = True

#-------------------------------------------------------------------------------
arcsec2Sr = (Math.PI/(60.*60.*180))**2
spireAreaEffFreq = {"PSW":450.*arcsec2Sr, "PMW":795.*arcsec2Sr, "PLW":1665.*arcsec2Sr}
# Neptune spectral index (from ESA4 model)
alphaNep={"PSW":1.29, "PMW":1.42, "PLW":1.47}

# Exponent of powerlaw describing FWHM dependence on frequency
# FWHM ~ freq**gamma
gamma = -0.85 

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

#-------------------------------------------------------------------------------
# Load SPIRE filter functions

# SPIRE band names
spireBands=["PSW","PMW","PLW"]

# SPIRE Photometer RSRF calibration product from cal tree
rsrfVersion = "2"
rsrf = fitsReader("%s//Phot//SCalPhotRsrf//SCalPhotRsrf_v%s.fits"%(directory, rsrfVersion))
# Photometer RSRF
spireFreq   = rsrf['rsrf']['frequency'].data*1e9  # Frequency in Hz
#indexes of freq in rsrf
ixR = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))

# SPIRE aperture efficiency product from cal tree
apertureEfficiencyVersion = "1"
apertureEfficiency = fitsReader("%s//Phot//SCalPhotApertureEfficiency//SCalPhotApertureEfficiency_v%s.fits"%(directory, apertureEfficiencyVersion))
spireApEffFreq = apertureEfficiency.getApertEffTable()["frequency"].data * 1e9 #comes in [GHz]
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
	interpAp = LinearInterpolator(spireApEffFreq, apertureEfficiency.getApertEffTable()[band].data)
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

#create product and set basic metadata
# define the start and end dates for the product
# starting at beginning of PFM1
startDate = df.parse("2005.02.22/00:00:00/GMT")
endDate   = df.parse("2020.01.01/00:00:00/GMT")

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

#read in beams from file
beamProfs["constant"] = asciiTableReader(os.path.join(dataDir,beamNewFileConstant))
beamProfs["core"] = asciiTableReader(os.path.join(dataDir,beamNewFileCore))
#update labels
beamProfs["core"].setDescription("Frequency-dependent core part of the radial beam profile")
beamProfs["constant"].setDescription("Frequency-independent constant part of the radial beam profile")
#generate file list (for metadata

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

#-------------------------------------------------------------------------------
# Add metadata
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

#-------------------------------------------------------------------------------
#===============================================================================
#=====                    CALCULATE EFFECTIVE FREQUENCIES                  =====
#===============================================================================
#-------------------------------------------------------------------------------
#uses existing numbers as initial guess
spireEffFreq = {"PSW":1217.27*1.e9,\
	"PMW":867.75*1.e9,\
	"PLW":610.87*1.e9}
print '\nGenerating new effective frequencies...'

for band in spireBands:
	spireEffFreq[band] = spireFindEffFreq(freq, spireFiltOnly[band],
	  beamProfs,spireEffFreq[band], gamma, spireAreaEffFreq[band],
	  alphaNep[band], band, verbose=verbose)

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

#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CALCULATE MONOCHROMATIC BEAM AREAS                 =====
#===============================================================================
#-------------------------------------------------------------------------------

#calculate monochromatic beam areas using full or simple beam treatment
beamMonoArea={}
beamAreaPipSr  = {}
beamAreaPipArc  = {}
print '\nCalculating monochromatic beam areas...'
for band in spireBands:
	#monochromatic beam areas
	beamMonoArea[band] = spireMonoAreas(freq, beamProfs, 
	  spireEffFreq[band], gamma, band)
	#pipeline beam areas
	beamAreaPipSr[band]=spireEffArea(freq, spireFiltOnly[band], \
	  beamMonoArea[band], BB=False, alpha=-1)
	beamAreaPipArc[band]=beamAreaPipSr[band]/arcsec2Sr

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

filename = java.io.File(r"%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_v%s.fits"%(directory, version))
beamProfs.meta['fileName'] = herschel.ia.dataset.StringParameter(value=filename.name)
fitsWriter = FitsArchive()
fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
fitsWriter.save(filename.toString(), beamProfs)
