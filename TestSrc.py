
from herschel.spire.all import *
from herschel.spire.util import *
from herschel.ia.all import *
from herschel.ia.task.mode import *
from herschel.ia.pg import ProductSink
from java.lang import *
from java.util import *
from os.path import join as pjoin
import herschel.spire.ia.pipeline.common.util
import math
import os
from os.path import join as pj

#obsid = 0x5001301E #Neptune shadow obs 1
#obsid = 0x5001301F] #Neptune shadow obs 2
obsid = 0x500060E6   #HeVICS V3
objName='HeViCS_v3'
pool = HttpClientPool ("http://wakefield.bnsc.rl.ac.uk/hcss/pal", "ops")
pools   = [pool,pool]
outDir     = "/home/astrog82/spxcen/Herschel/Calibration/ArtificialSources"

fits=FitsArchive()

###########################################################################
#### Set source options ####
########################################################################### 
# Select input peak in Jy/Beam
peak = 1.0

alphaFWHM= [[-4.0 , 17.0964582 , 23.2158027 , 33.4193693 ], \
	[-3.5 , 17.1587462 , 23.2954734 , 33.6264675  ], \
	[-3.0 , 17.2210894 , 23.3765151 , 33.8366761  ], \
	[-2.5 , 17.2839022 , 23.458792 , 34.0493927  ], \
	[-2.0 , 17.3473485 , 23.5421564 , 34.2639343  ], \
	[-1.5 , 17.4114569 , 23.6264504 , 34.4795555  ], \
	[-1.0 , 17.4761811 , 23.7115072 , 34.695467  ], \
	[-0.5 , 17.5414331 , 23.7971529 , 34.9108563  ], \
	[0.0 , 17.6071016 , 23.8832086 , 35.1249068  ], \
	[0.5 , 17.6730623 , 23.9694922 , 35.3368182  ], \
	[1.0 , 17.7391841 , 24.0558209 , 35.5458238  ], \
	[1.5 , 17.8053334 , 24.1420136 , 35.7512068  ], \
	[2.0 , 17.8713775 , 24.227893 , 35.9523134  ], \
	[2.5 , 17.937186 , 24.313288 , 36.1485637  ], \
	[3.0 , 18.0026338 , 24.3980365 , 36.3394589  ], \
	[3.5 , 18.0676021 , 24.4819875 , 36.5245852  ], \
	[4.0 , 18.1319799 , 24.5650037 , 36.7036159  ], \
	[4.5 , 18.1956657 , 24.6469637 , 36.8763098  ], \
	[5.0 , 18.2585678 , 24.7277653 , 37.0425076  ]]

print alphaFWHM[0]

alpha=alphaFWHM[12][0]
fwhmPxW= {'PSW':alphaFWHM[12][1], 'PMW': alphaFWHM[12][2], 'PLW': alphaFWHM[12][3]}
extent = 5.0 #number of beam-widths to apply source to

# Do you want to set the pixel size
manPix = True
if manPix:
	# Set pixel size
	pixsize = {'PSW':6, 'PMW':10, 'PLW':14}

# calculate sigma of gaussian
sigmaPxW = {}
sigmaPxW['PSW'] = fwhmPxW['PSW'] / (2.0 * math.sqrt(2.0*math.log(2.0)) * 60.0 * 60.0)
sigmaPxW['PMW'] = fwhmPxW['PMW'] / (2.0 * math.sqrt(2.0*math.log(2.0)) * 60.0 * 60.0)
sigmaPxW['PLW'] = fwhmPxW['PLW'] / (2.0 * math.sqrt(2.0*math.log(2.0)) * 60.0 * 60.0)

# calculate integral of gaussian
lengthPxW = {}
lengthPxW['PSW'] = extent * fwhmPxW['PSW'] / (2.0*60.0*60.0)
lengthPxW['PMW'] = extent * fwhmPxW['PMW'] / (2.0*60.0*60.0)
lengthPxW['PLW'] = extent * fwhmPxW['PLW'] / (2.0*60.0*60.0)
#integral = 2.0 * peak * math.pi * sigma**2.0 * (1.0 - math.exp(-length**2.0/(2.0*sigma**2.0)))
#integral = integral / (math.pi * length**2.0)
stepPxW={}
stepPxW['PSw'] = lengthPxW['PSW']/100.0
stepPxW['PMw'] = lengthPxW['PMW']/100.0
stepPxW['PLw'] = lengthPxW['PLW']/100.0

sum = 0.0
#for j in range(-2000,2001):
#	for i in range(-2000,2001):
#		if (i * step)**2.0 + (j*step)**2.0 > length**2.0:
#			continue
#		sum += peak * math.exp(-((i*step)**2.0+(j*step)**2.0)/(2.0*sigma**2.0))
#PSWsource = sum * (step/(pixsize["PSW"]/3600.0))**2.0
#PMWsource = sum * (step/(pixsize["PMW"]/3600.0))**2.0
#PLWsource = sum * (step/(pixsize["PLW"]/3600.0))**2.0
#print "Total Source flux in PSW is ", PSWsource, " Jy/Beam"
#print "Total Source flux in PMW is ", PMWsource, " Jy/Beam"
#print "Total Source flux in PLW is ", PLWsource, " Jy/Beam"


print "Reading in ObsID =", obsid,"("+hex(obsid)+") from Pool ",pool

###########################################################################
#### Apply source to observation ####
########################################################################### 

###get location of objects from obsid
query=MetaQuery(ObservationContext,"p","p.meta['obsid'].value==%iL"%obsid)
store = ProductStorage (pool)    # Just the Data
store.authenticate()
refs=store.select(query)
obs_ref=refs[0].product
level1=obs_ref.level1

##make initial map
##basic baseline removal
scans=level1
#scans=baselineRemovalMedian(level1)

scansout=[]
first=True
for pdtref in scans.refs:

	pdt=pdtref.product
	bolometers=pdt.getChannelNames()
	
	for bolo in bolometers:
		# if first time then calculate where the centre is best
		if first:
			# create centre
			centre = [0.0,0.0]
			centre[0] = (max(pdt.getRa("PSWE8")) - min(pdt.getRa("PSWE8"))) / 2.0 + min(pdt.getRa("PSWE8"))
			centre[1] = (max(pdt.getDec("PSWE8")) - min(pdt.getDec("PSWE8"))) / 2.0 + min(pdt.getDec("PSWE8"))
			first = False
		
		# skip anything that is not a bolometer
		if bolo[3:4] == "T" or bolo[3:5] == "DP" or bolo[3:4] == "R":
			continue
				
		if bolo.find('PSW')>=0:
			band='PSW'
		elif bolo.find('PMW')>=0:
			band='PMW'
		elif bolo.find('PLW')>=0:
			band='PLW'

		FWHM=fwhmPxW[band]
		sigma=sigmaPxW[band]
		length=lengthPxW[band]
		# locate all samples that are in the range required
		boloRA = pdt.getRa(bolo)
		boloDEC = pdt.getDec(bolo)
		boloSignal = pdt.getSignal(bolo)
		selection = boloRA.where(((boloRA - centre[0]) * math.cos(centre[1] * math.pi / 180.0))**2.0 + (boloDEC - centre[1])**2.0 <= (FWHM/(2.0*60.0*60.0)*extent)**2.0)
		objRA = boloRA[selection]
		objDEC = boloDEC[selection]
		objSignal = boloSignal[selection]
		
		# add the source to the timeline
		for i in range(0,len(objSignal)):
			radSqu = ((objRA[i] -centre[0])*math.cos(objDEC[i] * math.pi/180.0))**2.0 + (objDEC[i]-centre[1])**2.0
			objSignal[i] = objSignal[i] + peak * math.exp(-radSqu / (2.0*sigma**2.0))
		# add modified values back to full array
		boloSignal[selection] = objSignal
		
		# set modified array back to pdt
		pdt.setSignal(bolo, boloSignal)
	
	scansout.append(pdt)

scanConCal=ScanContext(scansout)

if manPix:
	PSWmap = naiveScanMapper(scanConCal, array='PSW', resolution = pixsize['PSW'])
	PMWmap = naiveScanMapper(scanConCal, array='PMW',resolution = pixsize['PMW'])
	PLWmap = naiveScanMapper(scanConCal, array='PLW',resolution = pixsize['PLW'])
else:
	PSWmap = naiveScanMapper(scanConCal, array='PSW')
	PMWmap = naiveScanMapper(scanConCal, array='PMW')
	PLWmap = naiveScanMapper(scanConCal, array='PLW')

# for PSW
meta = PSWmap.getMeta()
meta.set('author',StringParameter('C. North (Cardiff) - artificial source insertion'))
PSWmap.setMeta(meta)
# for PMW
meta = PMWmap.getMeta()
meta.set('author',StringParameter('C. North (Cardiff) - artificial source insertion'))
PMWmap.setMeta(meta)
# for PLW
meta = PLWmap.getMeta()
meta.set('author',StringParameter('C. North (Cardiff) - artificial source insertion'))
PLWmap.setMeta(meta)

# Save maps
fits.save(pj(outDir,objName+"_"+obsid+"_PSWmap-src-"+alpha+".fits"),PSWmap)
fits.save(pj(outDir,objName+"_"+obsid+"_PMWmap-src-"+alpha+".fits"),PMWmap)
fits.save(pj(outDir,objName+"_"+obsid+"_PLWmap-src-"+alpha+".fits"),PLWmap)



## Save Maps to output directory
#simpleFitsWriter(psw_combined, outDir+"MarBeam_psw_1arcsec_combined.fits")
#simpleFitsWriter(pmw_combined, outDir+"MarBeam_pmw_1arcsec_combined.fits")
#simpleFitsWriter(plw_combined, outDir+"MarBeam_plw_1arcsec_combined.fits")
#print "Map saved as FITS files to %s"%(outDir)

# End of Script
###########################################################################



