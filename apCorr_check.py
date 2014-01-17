from herschel.spire.all import *
from herschel.spire.util import *
from herschel.ia.all import *
from herschel.ia.task.mode import *
from herschel.ia.pg import ProductSink
from herschel.ia.toolbox.astro import Planck
from java.lang import *
from java.util import *
from os.path import join as pjoin
import herschel.spire.ia.pipeline.common.util
from herschel.share.unit import Frequency
from herschel.share.unit import Angle
from herschel.share.unit import SolidAngle

import math
import os
import urllib
from os.path import join as pj
import sys

#Working directory
workDir = Configuration.getProperty('var.hcss.workdir')

##Read cal data from file
try:
	#check whether cal exists
	cal
except:
	print "Reading calibration data from file"
	cal = spireCal(jarFile='/home/chris/dp_spire/cal/spire_cal_11_0.jar')

#BAApCorr_alpha=Float2d([-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0])
#BSApCorr_noBG=Float2d([[1.277,1.227,1.216],[1.275,1.226,1.214],[1.274,1.226,1.213],[1.273,1.226,1.212],[1.272,1.225,1.21],[1.271,1.225,1.209],[1.27,1.224,1.208],[1.269,1.224,1.207],[1.269,1.223,1.205],[1.268,1.223,1.204],[1.267,1.222,1.203],[1.267,1.222,1.202],[1.266,1.222,1.202],[1.265,1.221,1.2],[1.264,1.221,1.199],[1.264,1.220,1.198],[1.263,1.22,1.197],[1.262,1.219,1.196],[1.262,1.219,1.195]])
#BSApCorr_BG=Float2d([[1.282,1.234,1.236],[1.28,1.234,1.234],[1.279,1.234,1.232],[1.277,1.233,1.231],[1.276,1.233,1.229],[1.276,1.232,1.227],[1.275,1.232,1.226],[1.274,1.231,1.224],[1.273,1.231,1.222],[1.272,1.231,1.221],[1.272,1.23,1.219],[1.271,1.23,1.218],[1.27,1.229,1.216],[1.269,1.229,1.215],[1.269,1.229,1.213],[1.268,1.228,1.212],[1.267,1.228,1.211],[1.266,1.227,1.209],[1.266,1.227,1.208]])

BSApCorr_alpha=Float1d([-4.0,-3.5,-3.0,-2.0,-2.5,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0])
BSApCorr_noBG=Float2d([[1.277,1.227,1.216],[1.275,1.226,1.214],[1.274,1.226,1.213],[1.272,1.225,1.21],[1.273,1.226,1.212],[1.271,1.225,1.209],[1.27,1.224,1.208],[1.269,1.224,1.207],[1.269,1.223,1.205],[1.268,1.223,1.204],[1.267,1.222,1.203],[1.267,1.222,1.202],[1.266,1.222,1.202],[1.265,1.221,1.2],[1.264,1.221,1.199],[1.264,1.220,1.198],[1.263,1.22,1.197],[1.262,1.219,1.196],[1.262,1.219,1.195]])
BSApCorr_BG=Float2d([[1.282,1.234,1.236],[1.28,1.234,1.234],[1.279,1.234,1.232],[1.276,1.233,1.229],[1.277,1.233,1.231],[1.276,1.232,1.227],[1.275,1.232,1.226],[1.274,1.231,1.224],[1.273,1.231,1.222],[1.272,1.231,1.221],[1.272,1.23,1.219],[1.271,1.23,1.218],[1.27,1.229,1.216],[1.269,1.229,1.215],[1.269,1.229,1.213],[1.268,1.228,1.212],[1.267,1.228,1.211],[1.266,1.227,1.209],[1.266,1.227,1.208]])

fileSCalApCorr='SCalPhotColorCorrAperture_noBG_v1.fits'
fileSCalApCorrBg='SCalPhotColorCorrAperture_incBG_v1.fits'
try:
	ApCorrSCalIn=fitsReader(file=os.path.join(workDir,fileSCalApCorr))
except:
	#urllib.urlretrieve ("ftp://www.spire.rl.ac.uk/newCalTree/"+fileSCalApCorr,\
	#	    os.path.join(workDir,fileSCalApCorr))
	#ApCorrSCalIn=fitsReader(file=os.path.join(workDir,fileSCalApCorr))
	print 'no file found'

try:
	ApCorrBgSCalIn=fitsReader(file=os.path.join(workDir,fileSCalApCorrBg))
except:
	#urllib.urlretrieve ("ftp://www.spire.rl.ac.uk/newCalTree/"+fileSCalApCorrBg,\
	#	    os.path.join(workDir,fileSCalApCorrBg))
	#ApCorrBgSCalIn=fitsReader(file=os.path.join(workDir,fileSCalApCorrBg))
	print 'no file found'

ApCorrSCal={'noBg': \
	{'alpha':ApCorrSCalIn['alpha']['alpha'].data, \
	'PSW':ApCorrSCalIn['alpha']['PSW'].data, \
	'PMW':ApCorrSCalIn['alpha']['PMW'].data, \
	'PLW':ApCorrSCalIn['alpha']['PLW'].data}, \
	'withBg': \
	{'alpha':ApCorrBgSCalIn['alpha']['alpha'].data, \
	'PSW':ApCorrBgSCalIn['alpha']['PSW'].data, \
	'PMW':ApCorrBgSCalIn['alpha']['PMW'].data, \
	'PLW':ApCorrBgSCalIn['alpha']['PLW'].data}}

fileApCorr='/data/Herschel/Calibration/CalProducts/HIPE12/apcorr_alpha_Cal_131031.csv'
chkApCorr={'noBg':
	{'alpha':asciiTableReader(file=fileApCorr, tableType='CSV')['c0'].data, \
	'PSW':asciiTableReader(file=fileApCorr, tableType='CSV')['c1'].data, \
	'PMW':asciiTableReader(file=fileApCorr, tableType='CSV')['c2'].data, \
	'PLW':asciiTableReader(file=fileApCorr, tableType='CSV')['c3'].data}, \
	'withBg':
	{'alpha':asciiTableReader(file=fileApCorr, tableType='CSV')['c0'].data, \
	'PSW':asciiTableReader(file=fileApCorr, tableType='CSV')['c4'].data, \
	'PMW':asciiTableReader(file=fileApCorr, tableType='CSV')['c5'].data, \
	'PLW':asciiTableReader(file=fileApCorr, tableType='CSV')['c6'].data}}

pac=PlotXY()
lac1=LayerXY(BSApCorr_alpha,BSApCorr_noBG[:,0],color=java.awt.Color(0,0,150))
lac2=LayerXY(BSApCorr_alpha,BSApCorr_noBG[:,1],color=java.awt.Color(0,150,0))
lac3=LayerXY(BSApCorr_alpha,BSApCorr_noBG[:,2],color=java.awt.Color(150,0,0))
pac.addLayers([lac1,lac2,lac3])
lac4=LayerXY(ApCorrSCal['noBg']['alpha'],ApCorrSCal['noBg']['PSW'],color=java.awt.Color(0,0,255))
lac4.line=Style.DASHED
lac5=LayerXY(ApCorrSCal['noBg']['alpha'],ApCorrSCal['noBg']['PMW'],color=java.awt.Color(0,255,0))
lac5.line=Style.DASHED
lac6=LayerXY(ApCorrSCal['noBg']['alpha'],ApCorrSCal['noBg']['PLW'],color=java.awt.Color(255,0,0))
lac6.line=Style.DASHED
pac.addLayers([lac4,lac5,lac6])
pac.setTitleText('Aperture Corrections Check (no BG)')
pac.setSubtitleText('Python vs SCal')
pac.setXtitle('Spectral Index (alpha)')
pac.setYtitle('Aperture Correction')

pacd=PlotXY()
pacd.addLayer(LayerXY(BSApCorr_alpha,1.-BSApCorr_noBG[:,0]/ApCorrSCal['noBg']['PSW'],color=java.awt.Color(0,0,150)))
pacd.addLayer(LayerXY(BSApCorr_alpha,1.-BSApCorr_noBG[:,1]/ApCorrSCal['noBg']['PMW'],color=java.awt.Color(0,150,0)))
pacd.addLayer(LayerXY(BSApCorr_alpha,1.-BSApCorr_noBG[:,2]/ApCorrSCal['noBg']['PLW'],color=java.awt.Color(150,0,0)))
pacd.setTitleText('Aperture Corrections Difference (no BG)')
pacd.setSubtitleText('SCal - BS')
pacd.setXtitle('Spectral Index (alpha)')
pacd.setYtitle('Aperture Correction Difference')

pacb=PlotXY()
lacb1=LayerXY(BSApCorr_alpha,BSApCorr_BG[:,0],color=java.awt.Color(0,0,150))
lacb2=LayerXY(BSApCorr_alpha,BSApCorr_BG[:,1],color=java.awt.Color(0,150,0))
lacb3=LayerXY(BSApCorr_alpha,BSApCorr_BG[:,2],color=java.awt.Color(150,0,0))
pacb.addLayers([lacb1,lacb2,lacb3])
lacb4=LayerXY(ApCorrSCal['withBg']['alpha'],ApCorrSCal['withBg']['PSW'],color=java.awt.Color(0,0,255))
lacb4.line=Style.DASHED
lacb5=LayerXY(ApCorrSCal['withBg']['alpha'],ApCorrSCal['withBg']['PMW'],color=java.awt.Color(0,255,0))
lacb5.line=Style.DASHED
lacb6=LayerXY(ApCorrSCal['withBg']['alpha'],ApCorrSCal['withBg']['PLW'],color=java.awt.Color(255,0,0))
lacb6.line=Style.DASHED
pacb.addLayers([lacb4,lacb5,lacb6])
pacb.setTitleText('Aperture Corrections Check (with BG)')
pacb.setSubtitleText('SCal va BS')
pacb.setXtitle('Spectral Index (alpha)')
pacb.setYtitle('Aperture Correction')

pacbd=PlotXY()
pacbd.addLayer(LayerXY(BSApCorr_alpha,1.-BSApCorr_BG[:,0]/ApCorrSCal['withBg']['PSW'],color=java.awt.Color(0,0,150)))
pacbd.addLayer(LayerXY(BSApCorr_alpha,1.-BSApCorr_BG[:,1]/ApCorrSCal['withBg']['PMW'],color=java.awt.Color(0,150,0)))
pacbd.addLayer(LayerXY(BSApCorr_alpha,1.-BSApCorr_BG[:,2]/ApCorrSCal['withBg']['PLW'],color=java.awt.Color(150,0,0)))
pacbd.setTitleText('Aperture Corrections Difference (with BG)')
pacbd.setSubtitleText('SCal - BS')
pacbd.setXtitle('Spectral Index (alpha)')
pacbd.setYtitle('Aperture Correction Difference')

