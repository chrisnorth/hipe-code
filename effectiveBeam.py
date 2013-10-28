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
import sys

def readBeam(band,fromFile=True,maxRad=600):
	
	radArr=Double1d(maxRad)
	if (fromFile):
		##read from CSV file
		fileCore=band+'_Core.csv'
		fileConstant=band+'_Constant.csv'
		print 'need to write reading code'
		beamCore=Double1d(maxRad)
		beamCore=Double1d(maxRad)
	else:
		print 'Cal table input not possible'
		sys.exit()
		##read from Cal table

	return(maxRad,beamCore,beamConst)

def combBeam(beamCore,beamConst):
	nrad=len(beamCore)
	if len(beamConst) != nrad:
		print 'Core and Constant beams of different lengths [%d , %d]'%(len(beamCore),len(beamConst))

	beamComb=Double1d(nrad)
	for i in range(nrad):
		beamComb[i]=max(beamCore[i],beamConst[i])
	
	return(beamComb)

def effectiveBeam(beamRad,beamCore,beamConst,band='PSW',effNu=None,nuArr=None,rsrf=None,alpha=-1.0,gamma=0.85):
	
	if effNu==None:
		effNu={'PSW':1217.27e9,'PMW':867.75e9,'PLW':610.87e9}

	effNuBand=effNu[band]

	