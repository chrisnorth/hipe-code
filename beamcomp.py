##compare Calibration and CEN beams

cal = spireCal(jarFile='/home/chris/dp-spire/cal/spire_cal_8_phot_27Sep2011.jar')

beamPswCal = cal.getPhot().getBeamProf("PSW")
beamPmwCal = cal.getPhot().getBeamProf("PMW")
beamPlwCal = cal.getPhot().getBeamProf("PLW")

wcsPswCal = beamPswCal.wcs
wcsPmwCal = beamPmwCal.wcs
wcsPlwCal = beamPlwCal.wcs

beamPath = '/data/Herschel/Calibration/Inputs/spire_beams_measured/'

beamPswCEN = simpleFitsReader(beamPath+'psw_beam_1arcsec_withblanks.fits')
beamPmwCEN = simpleFitsReader(beamPath+'pmw_beam_1arcsec_withblanks.fits')
beamPlwCEN = simpleFitsReader(beamPath+'plw_beam_1arcsec_withblanks.fits')

wcsPswCEN = beamPswCEN.wcs
wcsPmwCEN = beamPmwCEN.wcs
wcsPlwCEN = beamPlwCEN.wcs

beamPswCalRegrid = regrid(source=beamPswCal, wcs=wcsPswCEN)
beamPmwCalRegrid = regrid(source=beamPmwCal, wcs=wcsPmwCEN)
beamPlwCalRegrid = regrid(source=beamPlwCal, wcs=wcsPlwCEN)

beamPswDiff = beamPswCEN
beamPmwDiff = beamPmwCEN
beamPlwDiff = beamPlwCEN

beamPswDiff.image[:,:] = beamPswCEN.image[:,:] - beamPswCalRegrid.image[:,:]
beamPmwDiff.image[:,:] = beamPmwCEN.image[:,:] - beamPmwCalRegrid.image[:,:]
beamPlwDiff.image[:,:] = beamPlwCEN.image[:,:] - beamPlwCalRegrid.image[:,:]

simpleFitsWriter(beamPswDiff,beamPath+'psw_beamDiff_1arcsec.fits')
simpleFitsWriter(beamPmwDiff,beamPath+'pmw_beamDiff_1arcsec.fits')
simpleFitsWriter(beamPlwDiff,beamPath+'plw_beamDiff_1arcsec.fits')
