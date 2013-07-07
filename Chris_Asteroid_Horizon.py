from herschel.share.fltdyn.ephem.horizons import BasicHorizons
from herschel.share.fltdyn.math import Direction
from herschel.ia.obs.auxiliary.fltdyn import Ephemerides
from herschel.ia.obs.auxiliary.fltdyn import Horizons

#### Load the horizons and orbit files using your favourite way. 
#### Alternative loading could be done with obs.calibration.horizons and obs.calibration.orbit
horizonsProduct  = FitsArchive().load('/home/alp/Spire_work/SPIRE_3615/hauxauxhorizons.fits')
ephemProduct = FitsArchive().load('/home/alp/Spire_work/SPIRE_3615/hauxauxorbitp.fits')
naifid = ppt.meta["naifId"].value

ephemObs    = Ephemerides(ephemProduct)
horizonsObs = Horizons(horizonsProduct, ephemObs)


### Comment out accordingly ####

##GEOMETRIC state (No correction)
usePointingCorrection = BasicHorizons.Correction.NONE 

##APPARENT state (Light Time and Stellar aberration correction)
## This is the Ephemeride position  that Tanya used and reported in SPIRE-3615
usePointingCorrection = BasicHorizons.Correction.LTS 

## ASTROMETRIC state (Light Time correction only)
## It has been agreed that all 3 Herschel instruments will apply this correction to SSO poinintg. 
## Spire Map positions of SSOs should be compared to Ephemeride positions obtained in the
## way shown below
usePointingCorrection = BasicHorizons.Correction.LT


### pptList[0].startDate: The start time of the earlies ppt in the obsid
### pptList[7].endDate: The start time of the latest ppt in the obsid
dir1 = Direction(horizonsObs.spacecraftVectorTo(int(naifid), pptList[0].startDate, usePointingCorrection))
dir2 = Direction(horizonsObs.spacecraftVectorTo(int(naifid), pptList[7].endDate, usePointingCorrection))

print dir1
print dir2

posEphemStart = [dir1.raDegrees,dir1.decDegrees]
posEphemEnd   = [dir2.raDegrees,dir2.decDegrees]



