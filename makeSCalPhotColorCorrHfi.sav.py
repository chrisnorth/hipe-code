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
#    SCalPhotColorCorrHfi product
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
#  Ed Polehampton  - 23/Oct/2013 - 2.5  Adapt to produce the new SCalPhotColorCorrHfi calibration
#                                       product (SPCAL-77)
#                                       Removed K4 parameters from metadata
#===============================================================================

scriptVersionString = "makeSCalPhotColorCorrHfi.py $Revision: 1.1 $"

directory = "..//..//..//..//..//..//data//spire//cal//SCal"
dataDir = "//disks//winchester2//calibration_data//"

df  = java.text.SimpleDateFormat("yyyy.MM.dd/HH:mm:ss/z")

#-------------------------------------------------------------------------------
# Input parameters

# Colour correction table version
version = "1"
formatVersion = "1.0"


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
hfiRIMOFile = dataDir+hfiRIMOFileName

# SPIRE aperture efficiency product from cal tree
apertureEfficiencyVersion = "1"
apertureEfficiency = fitsReader("%s//Phot//SCalPhotApertureEfficiency//SCalPhotApertureEfficiency_v%s.fits"%(directory, apertureEfficiencyVersion))
 
# SPIRE Photometer RSRF calibration product from cal tree
rsrfVersion = "2"
rsrf = fitsReader("%s//Phot//SCalPhotRsrf//SCalPhotRsrf_v%s.fits"%(directory, rsrfVersion))

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
nuMin   = 150.e9
nuMax   = 1500.e9
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

# Photometer RSRF
spireFreq   = rsrf['rsrf']['frequency'].data*1e9  # Frequency in Hz

ix = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))

# Get RSRF for normal point sources
# Interpolate to common frequency grid
interpPLW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['plw'].data)
interpPMW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['pmw'].data)
interpPSW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['psw'].data)

spireFiltPLW = Double1d(nNu)
spireFiltPMW = Double1d(nNu)
spireFiltPSW = Double1d(nNu)

spireFiltPLW[ix] = interpPLW(freq[ix])
spireFiltPMW[ix] = interpPMW(freq[ix])
spireFiltPSW[ix] = interpPSW(freq[ix])

# Aperture efficiency table

spireApEffFreq = apertureEfficiency['frequency']['frequency'].data * 1e9 #comes in [GHz]
spireApEffPsw  = apertureEfficiency['frequency']["PSW"].data
spireApEffPmw  = apertureEfficiency['frequency']["PMW"].data
spireApEffPlw  = apertureEfficiency['frequency']["PLW"].data

# Fold in and interpolate aperture efficiency
ix = freq.where((freq>=MIN(spireApEffFreq)) & (freq<=MAX(spireApEffFreq)))
interpPLW = CubicSplineInterpolator(spireApEffFreq, spireApEffPlw)
interpPMW = CubicSplineInterpolator(spireApEffFreq, spireApEffPmw)
interpPSW = CubicSplineInterpolator(spireApEffFreq, spireApEffPsw)

# Also for point source filter profiles
spireFiltPLW[ix] = interpPLW(freq[ix]) * spireFiltPLW[ix]
spireFiltPMW[ix] = interpPMW(freq[ix]) * spireFiltPMW[ix]
spireFiltPSW[ix] = interpPSW(freq[ix]) * spireFiltPSW[ix]


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
def hpXcalKcorr(freq0, freq, transm, BB=True, temp=20.0, beta=1.8, alpha=-1.0, gamma=0.0):
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
# Calculate K-correction factors for extended source assuming alpha=-1
k4E_PSW = hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, False, gamma=gamma)[0]
k4E_PMW = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, False, gamma=gamma)[0]
k4E_PLW = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, False, gamma=gamma)[0]
k4E_857 = hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857,   False)[0]
k4E_545 = hpXcalKcorr(hfiRefFreq[1],   freq, hfiFilt545,   False)[0]

# Calculate K-correction factors for point source assuming alpha=-1
k4P_PSW = hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, False)[0]
k4P_PMW = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, False)[0]
k4P_PLW = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, False)[0]

# Print Spire K4 factors for point source calibration for verification
print 'SPIRE k4 factors for point source: %s, %s, %s' % (k4P_PSW,k4P_PMW,k4P_PLW)

# Print Spire K4 factors for extended source calibration for verification
print 'SPIRE k4 factors for extended source: %s, %s, %s' % (k4E_PSW,k4E_PMW,k4E_PLW)

#-------------------------------------------------------------------------------
# Calculate and tabulate colour correction parameters for a range of temperatures
k545toPLW = Double1d()       # K-correction from HFI-545 to PLW
k857toPMW = Double1d()       # K-correction from HFI-857 to PMW
k857toPSW = Double1d()       # K-correction from HFI-857 to PSW
ratio545over857 = Double1d() # Ratio 545GHz to 857GHz filter

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

# define the start and end dates for the product
# starting at beginning of PFM1
startDate = df.parse("2005.02.22/00:00:00/GMT")
endDate   = df.parse("2020.01.01/00:00:00/GMT")

photColorCorrHfi = herschel.spire.ia.dataset.PhotColorCorrHfi()

photColorCorrHfi.meta["creator"].value   = scriptVersionString
photColorCorrHfi.meta["modelName"].value = "FM"
photColorCorrHfi.meta["creationDate"].value = FineTime(java.util.Date())
photColorCorrHfi.meta["startDate"].value = FineTime(startDate)
photColorCorrHfi.meta["endDate"].value   = FineTime(endDate)
photColorCorrHfi.meta["fileOrigin"]  = herschel.ia.dataset.StringParameter(value="%s"%hfiRIMOFileName, description="Origin of the data")
photColorCorrHfi.setVersion(version)
photColorCorrHfi.setFormatVersion(formatVersion)

# Save standard point source correction factor into meta data
#photColorCorrHfi.meta["k4P_PSW"] = DoubleParameter(k4P_PSW, "PSW point source K-factor")
#photColorCorrHfi.meta["k4P_PMW"] = DoubleParameter(k4P_PMW, "PMW point source K-factor")
#photColorCorrHfi.meta["k4P_PLW"] = DoubleParameter(k4P_PLW, "PLW point source K-factor")
#photColorCorrHfi.meta["k4E_PSW"] = DoubleParameter(k4E_PSW, "PSW extended source K-factor")
#photColorCorrHfi.meta["k4E_PMW"] = DoubleParameter(k4E_PMW, "PMW extended source K-factor")
#photColorCorrHfi.meta["k4E_PLW"] = DoubleParameter(k4E_PLW, "PLW extended source K-factor")

# Save beta dependent colour correction tables into binary table
if fixBeta == True:
	photColorCorrHfi.meta["beta"] = DoubleParameter(beta0, "Modified black-body spectral index used")
        photColorCorrHfi.setTempVals(tvect)      # Temperature of modified BB
        #
# Save temperature dependent colour correction tables into binary table
else:
	photColorCorrHfi.meta['temperature'] = DoubleParameter(temp0, "Modified black-body temperature used")
	photColorCorrHfi.meta['temperature'].unit = Temperature.KELVIN
	#
	photColorCorrHfi.setTempVals(bvect)          # Spectral index of modified BB

photColorCorrHfi.meta["gamma"] = DoubleParameter(gamma, "Exponent describing FWHM dependence on frequency used")

# Save colour correction tables into binary table
photColorCorrHfi.setRatio545_857CorrVals(ratio545over857)  # Ratio 545/857 GHz filter
photColorCorrHfi.setK545toPLWCorrVals(k545toPLW)           # 545GHz to PLW K-factor
photColorCorrHfi.setK857toPMWCorrVals(k857toPMW)           # 857GHz to PMW K-factor
photColorCorrHfi.setK857toPSWCorrVals(k857toPSW)           # 857GHz to PSW K-factor


##########################################
# ******* write to FITS *******
filename = java.io.File(r"%s//Phot//SCalPhotColorCorrHfi//SCalPhotColorCorrHfi_v%s.fits"%(directory, version))
#filename = java.io.File(r"%s/SCalPhotColorCorrHfi_v%s.fits"%(directory, version))
photColorCorrHfi.meta['fileName'] = herschel.ia.dataset.StringParameter(value=filename.name)
print
fitsWriter = FitsArchive()
fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
fitsWriter.save(filename.toString(), photColorCorrHfi)
print "written: %s"%filename.toString()
print

