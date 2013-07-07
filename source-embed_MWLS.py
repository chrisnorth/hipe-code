### Program to embed artificial sources into the timelines ###
# Created October 2010

#### import stuff #####
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

# Import the script tasks.py that contains the task definitions
from herschel.spire.ia.pipeline.scripts.POF5.POF5_tasks import *

# Input definition:
from herschel.spire.ia.pipeline.scripts.POF5.POF5_input import *

### Decide on inputs

# select folders
pdtfolder = '/home/herdata/spxmws/mega-drive/test/non-detec/pdt'
outfolder = '/home/herdata/spxmws/mega-drive/test/non-detec/test'

# Select input peak in Jy/Beam
peak = 0.5
FWHM = 60.0
extent = 5.0

# Do you want to set the pixel size
manPix = True
if manPix:
	# Set pixel size
	pixsize = {'PSW':6, 'PMW':8, 'PLW':12}

version = "P05F1"

################################################################################

# calculate sigma of gaussian
sigma = FWHM / (2.0 * math.sqrt(2.0*math.log(2.0)) * 60.0 * 60.0)

# calculate integral of gaussian
length  = extent * FWHM / (2.0*60.0*60.0)
#integral = 2.0 * peak * math.pi * sigma**2.0 * (1.0 - math.exp(-length**2.0/(2.0*sigma**2.0)))
#integral = integral / (math.pi * length**2.0)
step = length/100.0
sum = 0.0
for j in range(-2000,2001):
	for i in range(-2000,2001):
		if (i * step)**2.0 + (j*step)**2.0 > length**2.0:
			continue
		sum += peak * math.exp(-((i*step)**2.0+(j*step)**2.0)/(2.0*sigma**2.0))
PSWsource = sum * (step/(pixsize["PSW"]/3600.0))**2.0
PMWsource = sum * (step/(pixsize["PMW"]/3600.0))**2.0
PLWsource = sum * (step/(pixsize["PLW"]/3600.0))**2.0
print "Total Source flux in PSW is ", PSWsource, " Jy/Beam"
print "Total Source flux in PMW is ", PMWsource, " Jy/Beam"
print "Total Source flux in PLW is ", PLWsource, " Jy/Beam"

pause()

# Create FITS archive
fits = FitsArchive()

# Search folder for FITS files
allfiles = os.listdir(pdtfolder)
pdtfiles = [] 
for i in allfiles:
	if i[-4:] == 'fits':
		pdtfiles.append(i)

# locate objects... HRS only have 1 file, HeViCS have 8
objects = {}
for file in pdtfiles:
	name = file.split('-')[0]
	if objects.has_key(name):
		objects[name].append(file)
	else:
		objects[name] = [file]

# loop over all the objects
for object in objects.keys():
	# load in each observation and add in new source
	first = True
	
	# create empty array for map-making pdts
	obs = []
	
	for file in objects[object]:
		# read in data
		pdt = fits.load(pj(pdtfolder,file))
		
		# create a list of bolometers
		bolometers = pdt.getChannelNames()
		
		# if first time then calculate where the centre is best
		if first:
			# create centre
			centre = [0.0,0.0]
			centre[0] = (max(pdt.getRa("PSWE8")) - min(pdt.getRa("PSWE8"))) / 2.0 + min(pdt.getRa("PSWE8"))
			centre[1] = (max(pdt.getDec("PSWE8")) - min(pdt.getDec("PSWE8"))) / 2.0 + min(pdt.getDec("PSWE8"))
			first = False
		
		# now apply source to the timelines
		for bolo in bolometers:
			# skip anything that is not a bolometer
			if bolo[3:4] == "T" or bolo[3:5] == "DP" or bolo[3:4] == "R":
				continue
			
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
		
		# add pdt to obs list
		obs.append(pdt)
		
	# Make scan context for map making
	scanConCal=ScanContext(obs)
		
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
	meta.set('author',StringParameter('M.Smith (Cardiff) BriGAdE L1-L2 Pipeline'))
	PSWmap.setMeta(meta)
	# for PMW
	meta = PMWmap.getMeta()
	meta.set('author',StringParameter('M.Smith (Cardiff) BriGAdE L1-L2 Pipeline'))
	PMWmap.setMeta(meta)
	# for PLW
	meta = PLWmap.getMeta()
	meta.set('author',StringParameter('M.Smith (Cardiff) BriGAdE L1-L2 Pipeline'))
	PLWmap.setMeta(meta)
	
	# Save maps
	fits.save(pj(outfolder,object+"-PSWmap-mosaic_MS-"+version+".fits"),PSWmap)
	fits.save(pj(outfolder,object+"-PMWmap-mosaic_MS-"+version+".fits"),PMWmap)
	fits.save(pj(outfolder,object+"-PLWmap-mosaic_MS-"+version+".fits"),PLWmap)
	
	# delete large memory values
	del(pdt)
	del(obs)

print "Program finished successfully"
