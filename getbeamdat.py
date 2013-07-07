
from os.path import isfile
from herschel.ia.toolbox.util import simpleFitsWriter

##Import make_beamdat user module
import sys
sys.path.append('/data/Herschel/Calibration/RelGains/')
import make_beamdat

###############################################
##Set path and filenames
###############################################

dirPath='/data/Herschel/Calibration/RelGains/'
profFile='beamprof_rad.dat'
filePSW='bg_im_PSW.fits'
filePMW='bg_im_PMW.fits'
filePLW='bg_im_PLW.fits'

##Check if file exists

if not isfile(dirPath + filePSW):
	print 'File does not exist...'
	## Make new FITS files
	
	make_beamdat.makedat(dirPath,profFile,filePSW,filePMW,filePLW)
	print 'Now it does...'
else:
	print 'File exists...'

bgImPSW=simpleFitsReader(filePSW)
bgImPMW=simpleFitsReader(filePMW)
bgImPLW=simpleFitsReader(filePLW)
