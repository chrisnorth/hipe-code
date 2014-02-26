#check calibration versions

ver1="2."
ver2="3.1"

directory=Configuration.getProperty('var.hcss.workdir')

spireBands=["PSW","PMW","PLW"]

#-------------------------------------------------------------------------------
# RadialCorrBeam
beamProf1=fitsReader('%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_v%s.fits'%(directory,ver1))
beamProf2=fitsReader('%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_v%s.fits'%(directory,ver2))

print 'RadialCorrBeam (1-v%s/v%s):'%(ver2,ver1)
for band in spireBands:
	print '  core     (%s): min=%g, max=%g'%\
		(band,min(1.-beamProf2['core'][band].data/beamProf1['core'][band].data),\
		max(1.-beamProf2['core'][band].data/beamProf1['core'][band].data))
	print '  constant (%s): min=%g, max=%g'%\
		(band,min(1.-beamProf2['constant'][band].data/beamProf1['constant'][band].data),\
		max(1.-beamProf2['constant'][band].data/beamProf1['constant'][band].data))
	print '  normArea (%s): min=%g, max=%g'%\
		(band,min(1.-beamProf2['normArea'][band].data/beamProf1['normArea'][band].data),\
		max(1.-beamProf2['normArea'][band].data/beamProf1['normArea'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrBeam
Kbeam1=fitsReader('%s//Phot//SCalPhotColorCorrBeam//SCalPhotColorCorrBeam_v%s.fits'%(directory,ver1))
Kbeam2=fitsReader('%s//Phot//SCalPhotColorCorrBeam//SCalPhotColorCorrBeam_v%s.fits'%(directory,ver2))

print 'ColorCorrBeam (1-v%s/v%s):'%(ver2,ver1)
for band in spireBands:
	print '  alpha     (%s): min=%g, max=%g'%\
		(band,min(1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data),\
		max(1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data))
	print '  beta_2_00 (%s): min=%g, max=%g'%\
		(band,min(1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data),\
		max(1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrK_point
KPsrc1=fitsReader('%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_point_v%s.fits'%(directory,ver1))
KPsrc2=fitsReader('%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_point_v%s.fits'%(directory,ver2))

print 'ColorCorrK_point (1-v%s/v%s):'%(ver2,ver1)
for band in spireBands:
	print '  alpha     (%s): min=%g, max=%g'%\
		(band,min(1.-KPsrc2['alpha'][band].data/KPsrc1['alpha'][band].data),\
		max(1.-KPsrc2['alpha'][band].data/KPsrc1['alpha'][band].data))
	print '  beta_2_00 (%s): min=%g, max=%g'%\
		(band,min(1.-KPsrc2['beta_2_00'][band].data/KPsrc1['beta_2_00'][band].data),\
		max(1.-KPsrc2['beta_2_00'][band].data/KPsrc1['beta_2_00'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrK_extended
KExtd1=fitsReader('%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_extended_v%s.fits'%(directory,ver1))
KExtd2=fitsReader('%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_extended_v%s.fits'%(directory,ver2))

print 'ColorCorrK_extended (1-v%s/v%s):'%(ver2,ver1)
for band in spireBands:
	print '  alpha     (%s): min=%g, max=%g'%\
		(band,min(1.-KExtd2['alpha'][band].data/KExtd1['alpha'][band].data),\
		max(1.-KExtd2['alpha'][band].data/KExtd1['alpha'][band].data))
	print '  beta_2_00 (%s): min=%g, max=%g'%\
		(band,min(1.-KExtd2['beta_2_00'][band].data/KExtd1['beta_2_00'][band].data),\
		max(1.-KExtd2['beta_2_00'][band].data/KExtd1['beta_2_00'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrAperture_noBG
apCorrNoBG1=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_noBG_v%s.fits'%(directory,ver2))
apCorrNoBG2=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_Analytical_noBG_v%s.fits'%(directory,ver2))
#apCorrNoBG2=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_noBG_v%s.fits'%(directory,ver2))

print 'ColorCorrAperture_noBG (1-v%s/v%s):'%(ver2,ver1)
for band in spireBands:
	print '  alpha (%s): min=%g, max=%g'%\
		(band,min(1.-apCorrNoBG2['alpha'][band].data/apCorrNoBG1['alpha'][band].data),\
		max(1.-apCorrNoBG2['alpha'][band].data/apCorrNoBG1['alpha'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrAperture_incBG
apCorrIncBG1=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_incBG_v%s.fits'%(directory,ver2))
apCorrIncBG2=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_Analytical_incBG_v%s.fits'%(directory,ver2))
#apCorrIncBG2=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_incBG_v%s.fits'%(directory,ver2))

print 'ColorCorrAperture_incBG (1-v%s/v%s):'%(ver2,ver1)
for band in spireBands:
	print '  alpha (%s): min=%g, max=%g'%\
		(band,min(1.-apCorrIncBG2['alpha'][band].data/apCorrIncBG1['alpha'][band].data),\
		max(1.-apCorrIncBG2['alpha'][band].data/apCorrIncBG1['alpha'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrHfi
#KHfi1=fitsReader('%s//Phot//SCalPhotColorCorrHfi//SCalPhotColorCorrHfi_v%s.fits'%(directory,ver1))
KHfi1=cal12.getPhot().getProduct('ColorCorrHfi')
KHfi2=fitsReader('%s//Phot//SCalPhotColorCorrHfi//SCalPhotColorCorrHfi_v%s.fits'%(directory,ver2))

KHfi1i=KHfi1.copy()
KHfi2i=KHfi2.copy()
if len(KHfi2['colorCorr']['Temperature'].data) != len(KHfi1['colorCorr']['Temperature'].data):
	#interpolate ver 1 to match ver2 Temperature grid
	temp1=KHfi1['colorCorr']['Temperature'].data
	temp2=KHfi2['colorCorr']['Temperature'].data
	ix=temp2.where((temp2 >= MIN(temp1)) & (temp2 <= MAX(temp1)))
	KHfi1r_interp=CubicSplineInterpolator(KHfi1['colorCorr']['Temperature'].data,\
		KHfi1['colorCorr']['ratio545_857'].data)
	KHfi1plw_interp=CubicSplineInterpolator(KHfi1['colorCorr']['Temperature'].data,\
		KHfi1['colorCorr']['k545toPLW'].data)
	KHfi1pmw_interp=CubicSplineInterpolator(KHfi1['colorCorr']['Temperature'].data,\
		KHfi1['colorCorr']['k857toPMW'].data)
	KHfi1psw_interp=CubicSplineInterpolator(KHfi1['colorCorr']['Temperature'].data,\
		KHfi1['colorCorr']['k857toPSW'].data)
	KHfi1i['colorCorr']['Temperature'].data=KHfi2['colorCorr']['Temperature'].data[ix]
	KHfi1i['colorCorr']['ratio545_857'].data=KHfi1r_interp(KHfi1i['colorCorr']['Temperature'].data)
	KHfi1i['colorCorr']['k545toPLW'].data=KHfi1plw_interp(KHfi1i['colorCorr']['Temperature'].data)
	KHfi1i['colorCorr']['k857toPMW'].data=KHfi1pmw_interp(KHfi1i['colorCorr']['Temperature'].data)
	KHfi1i['colorCorr']['k857toPSW'].data=KHfi1psw_interp(KHfi1i['colorCorr']['Temperature'].data)

	KHfi2i['colorCorr']['Temperature'].data=KHfi2['colorCorr']['Temperature'].data[ix]
	KHfi2i['colorCorr']['ratio545_857'].data=KHfi2['colorCorr']['ratio545_857'].data[ix]
	KHfi2i['colorCorr']['k545toPLW'].data=KHfi2['colorCorr']['k545toPLW'].data[ix]
	KHfi2i['colorCorr']['k857toPMW'].data=KHfi2['colorCorr']['k857toPMW'].data[ix]
	KHfi2i['colorCorr']['k857toPSW'].data=KHfi2['colorCorr']['k857toPSW'].data[ix]
	print 'ColorCorrHfi (1-v%s/v%s) [v%s interpolated to match v%s]:'%(ver2,ver1,ver1,ver2)
else:
	print 'ColorCorrHfi (1-v%s/v%s):'%(ver2,ver1)

print '  ratio545_857: min=%g, max=%g'%\
	(min(1.-KHfi2i['colorCorr']['ratio545_857'].data/KHfi1i['colorCorr']['ratio545_857'].data),\
	max(1.-KHfi2i['colorCorr']['ratio545_857'].data/KHfi1i['colorCorr']['ratio545_857'].data))
print '  k545toPLW   : min=%g, max=%g'%\
	(min(1.-KHfi2i['colorCorr']['k545toPLW'].data/KHfi1i['colorCorr']['k545toPLW'].data),\
	max(1.-KHfi2i['colorCorr']['k545toPLW'].data/KHfi1i['colorCorr']['k545toPLW'].data))
print '  k857toPMW   : min=%g, max=%g'%\
	(min(1.-KHfi2i['colorCorr']['k857toPMW'].data/KHfi1i['colorCorr']['k857toPMW'].data),\
	max(1.-KHfi2i['colorCorr']['k857toPMW'].data/KHfi1i['colorCorr']['k857toPMW'].data))
print '  k857toPSW   : min=%g, max=%g'%\
	(min(1.-KHfi2i['colorCorr']['k857toPSW'].data/KHfi1i['colorCorr']['k857toPSW'].data),\
	max(1.-KHfi2i['colorCorr']['k857toPSW'].data/KHfi1i['colorCorr']['k857toPSW'].data))


sdfgf

cols={'PSW':java.awt.Color.BLUE,'PMW':java.awt.Color.GREEN,'PLW':java.awt.Color.RED}
pCore=PlotXY()
for band in spireBands:
	pCore.addLayer(LayerXY(beamProf2['core']['radius'].data,\
		1.-beamProf2['core'][band].data/beamProf1['core'][band].data, \
		color=cols[band],name='%s'%(band)))
	pCore.setTitleText('Beam Profile (core)')
	pCore.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pCore.xaxis.titleText = 'radius (arcsec)'
	pCore.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-beamProf2['core'][band].data/beamProf1['core'][band].data),\
		max(1.-beamProf2['core'][band].data/beamProf1['core'][band].data)))

pConst=PlotXY()
for band in spireBands:
	pConst.addLayer(LayerXY(beamProf2['constant']['radius'].data,\
		1.-beamProf2['constant'][band].data/beamProf1['constant'][band].data, \
		color=cols[band],name='%s'%(band)))
	pConst.setTitleText('Beam Profile (constant)')
	pConst.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pConst.xaxis.titleText = 'radius (arcsec)'
	pConst.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-beamProf2['constant'][band].data/beamProf1['constant'][band].data),\
		max(1.-beamProf2['constant'][band].data/beamProf1['constant'][band].data)))

pNorm=PlotXY()
for band in spireBands:
	pNorm.addLayer(LayerXY(beamProf2['normArea']['radius'].data,\
		1.-beamProf2['normArea'][band].data/beamProf1['normArea'][band].data, \
		color=cols[band],name='%s'%(band)))
	pNorm.setTitleText('Beam Profile (normArea)')
	pNorm.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pNorm.xaxis.titleText = 'radius (arcsec)'
	pNorm.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-beamProf2['normArea'][band].data/beamProf1['normArea'][band].data),\
		max(1.-beamProf2['normArea'][band].data/beamProf1['normArea'][band].data)))

pKBeama=PlotXY()
for band in spireBands:
	pKBeama.addLayer(LayerXY(Kbeam2['alpha']['alpha'].data,\
		1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data, \
		color=cols[band],name='%s'%(band)))
	pKBeama.setTitleText('Beam Correction (alpha)')
	pKBeama.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pKBeama.xaxis.titleText = 'Spectral index'
	pKBeama.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data),\
		max(1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data)))

pKBeamt=PlotXY()
for band in spireBands:
	pKBeamt.addLayer(LayerXY(Kbeam2['beta_2_00']['Temperature'].data,\
		1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data, \
		color=cols[band],name='%s'%(band)))
	pKBeamt.setTitleText('Beam Correction (beta=2.0)')
	pKBeamt.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pKBeamt.xaxis.titleText = 'Temperature (K)'
	pKBeamt.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data),\
		max(1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data)))

pKHFIr=PlotXY()
for band in spireBands:
	pKHFIr.addLayer(LayerXY(KHfi2i['colorCorr']['Temperature'].data,\
		1.-KHfi2i['colorCorr']['ratio545_857'].data/KHfi1i['colorCorr']['ratio545_857'].data, \
		color=cols[band],name='%s'%(band)))
	pKHFIr.setTitleText('SPIRE-HFI correction (545/857 ratio)')
	pKHFIr.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pKHFIr.xaxis.titleText = 'Temperature (K)'
	pKHFIr.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KHfi2['colorCorr']['ratio545_857'].data/KHfi1['colorCorr']['ratio545_857'].data),\
		max(1.-KHfi2['colorCorr']['ratio545_857'].data/KHfi1['colorCorr']['ratio545_857'].data)))

pKHFIpsw=PlotXY()
for band in spireBands:
	pKHFIpsw.addLayer(LayerXY(KHfi2i['colorCorr']['Temperature'].data,\
		1.-KHfi2i['colorCorr']['k857toPSW'].data/KHfi1i['colorCorr']['k857toPSW'].data, \
		color=cols[band],name='%s'%(band)))
	pKHFIpsw.setTitleText('SPIRE-HFI correction (857 to PSW)')
	pKHFIpsw.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pKHFIpsw.xaxis.titleText = 'Temperature (K)'
	pKHFIpsw.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KHfi2['colorCorr']['k857toPSW'].data/KHfi1['colorCorr']['k857toPSW'].data),\
		max(1.-KHfi2['colorCorr']['k857toPSW'].data/KHfi1['colorCorr']['k857toPSW'].data)))

pKHFIpmw=PlotXY()
for band in spireBands:
	pKHFIpmw.addLayer(LayerXY(KHfi2i['colorCorr']['Temperature'].data,\
		1.-KHfi2i['colorCorr']['k857toPMW'].data/KHfi1i['colorCorr']['k857toPMW'].data, \
		color=cols[band],name='%s'%(band)))
	pKHFIpmw.setTitleText('SPIRE-HFI correction (857 to PMW)')
	pKHFIpmw.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pKHFIpmw.xaxis.titleText = 'Temperature (K)'
	pKHFIpmw.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KHfi2['colorCorr']['k857toPMW'].data/KHfi1['colorCorr']['k857toPMW'].data),\
		max(1.-KHfi2['colorCorr']['k857toPMW'].data/KHfi1['colorCorr']['k857toPMW'].data)))

pKHFIplw=PlotXY()
for band in spireBands:
	pKHFIplw.addLayer(LayerXY(KHfi2i['colorCorr']['Temperature'].data,\
		1.-KHfi2i['colorCorr']['k545toPLW'].data/KHfi1i['colorCorr']['k545toPLW'].data, \
		color=cols[band],name='%s'%(band)))
	pKHFIplw.setTitleText('SPIRE-HFI correction (857 to PMW)')
	pKHFIplw.yaxis.titleText = '1-v%s/v%s'%(ver2,ver1)
	pKHFIplw.xaxis.titleText = 'Temperature (K)'
	pKHFIplw.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KHfi2['colorCorr']['k545toPLW'].data/KHfi1['colorCorr']['k545toPLW'].data),\
		max(1.-KHfi2['colorCorr']['k545toPLW'].data/KHfi1['colorCorr']['k545toPLW'].data)))