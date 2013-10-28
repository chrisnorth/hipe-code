import math
import os
from os.path import join as pj
import sys

##check output of effBeam.py
doCheck=True

doBeam=True
doArea=True
doApCorr=False

if doArea:
	if doCheck:
		print "Checking..."
		fileNorm='/data/Herschel/Calibration/Beam_plots/EffBeams/Unnorm_areas.csv'
		chkArea={'radius':asciiTableReader(file=fileNorm, tableType='CSV')['radius'].data[0:maxRad],
			'PSW':asciiTableReader(file=fileNorm, tableType='CSV')['PSW'].data[0:maxRad],
			'PMW':asciiTableReader(file=fileNorm, tableType='CSV')['PMW'].data[0:maxRad],
			'PLW':asciiTableReader(file=fileNorm, tableType='CSV')['PLW'].data[0:maxRad]}
	
		diffArea={'radius':beamsIn['radius'],
			'PSW':beamsIn['PSW']['AreaCumul']-chkArea['PSW'],
			'PMW':beamsIn['PMW']['AreaCumul']-chkArea['PMW'],
			'PLW':beamsIn['PLW']['AreaCumul']-chkArea['PLW']}
		relArea={'radius':beamsIn['radius'],
			'PSW':1.-beamsIn['PSW']['AreaCumul']/chkArea['PSW'],
			'PMW':1.-beamsIn['PMW']['AreaCumul']/chkArea['PMW'],
			'PLW':1.-beamsIn['PLW']['AreaCumul']/chkArea['PLW']}
	
		pchk=PlotXY()
		pchk.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(beamsIn['PSW']['AreaCumul'])))
		pchk.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(beamsIn['PMW']['AreaCumul'])))
		pchk.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(beamsIn['PLW']['AreaCumul'])))
		pchk.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(chkArea['PSW'])))
		pchk.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(chkArea['PMW'])))
		pchk.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(chkArea['PLW'])))
		pchk.setTitleText('Absolute Area Check')
		pchk.setSubtitleText('HIPE vs python')
		pchk.setXtitle('Radius [arcsec]')
		pchk.setYtitle('Cumulative Beam Area [sq. arcsec]')
	
		pchkrel=PlotXY()
		pchkrel.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(relArea['PSW'])))
		pchkrel.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(relArea['PMW'])))
		pchkrel.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(relArea['PLW'])))
		pchkrel.setTitleText('Relative Area Check')
		pchkrel.setSubtitleText('HIPE vs python')
		pchkrel.setXtitle('Radius [arcsec]')
		pchkrel.setYtitle('1 - HIPE/python')
	
		pchkrelz=PlotXY()
		pchkrelz.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(relArea['PSW'])))
		pchkrelz.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(relArea['PMW'])))
		pchkrelz.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(relArea['PLW'])))
		pchkrelz.setYrange([0,0.005])
		pchkrelz.setXrange([0,100])
		pchkrelz.setTitleText('Relative Area Check (zoomed)')
		pchkrelz.setSubtitleText('HIPE vs python')
		pchkrelz.setXtitle('Radius [arcsec]')
		pchkrelz.setYtitle('1 - HIPE/python')

		pchkdiff=PlotXY()
		pchkdiff.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(diffArea['PSW'])))
		pchkdiff.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(diffArea['PMW'])))
		pchkdiff.addLayer(LayerXY(Double1d(beamsIn['radius']),Double1d(diffArea['PLW'])))
		pchkdiff.setTitleText('Absolute Area Check')
		pchkdiff.setSubtitleText('HIPE vs python')
		pchkdiff.setXtitle('Radius [arcsec]')
		pchkdiff.setYtitle('HIPE - python')


if doArea:
	if doCheck:
		print "Checking..."
		fileArea='/data/Herschel/Calibration/CalProducts/HIPE11/beamarea_alpha_Cal.csv'
		chkAreaPy={'alpha':asciiTableReader(file=fileArea, tableType='CSV')['c0'].data,
			'PSW':asciiTableReader(file=fileArea, tableType='CSV')['c1'].data,
			'PMW':asciiTableReader(file=fileArea, tableType='CSV')['c2'].data,
			'PLW':asciiTableReader(file=fileArea, tableType='CSV')['c3'].data}
		chkArea={'alpha':[],'PSW':[],'PMW':[],'PLW':[]}
		for n in effBeams:
			if effBeams[n]['type']=='PowerLaw':
				chkArea['alpha'].append(effBeams[n]['alpha'])
				for band in bands:
					chkArea[band].append(effBeams[n][band]['Area'])
		for x in chkArea:
			chkArea[x]=Float1d(chkArea[x])
		pArea=PlotXY()
		lArea1=LayerXY(chkArea['alpha'],chkArea['PSW'],color=java.awt.Color(0,0,150))
		lArea2=LayerXY(chkArea['alpha'],chkArea['PMW'],color=java.awt.Color(0,150,0))
		lArea3=LayerXY(chkArea['alpha'],chkArea['PLW'],color=java.awt.Color(150,0,0))
		pArea.addLayer(lArea1)
		pArea.addLayer(lArea2)
		pArea.addLayer(lArea3)
		lArea4=LayerXY(chkAreaPy['alpha'],chkAreaPy['PSW'],name='PSW (Py)',color=java.awt.Color(0,0,150))
		lArea4.line=Style.DASHED
		lArea5=LayerXY(chkAreaPy['alpha'],chkAreaPy['PMW'],name='PMW (Py)',color=java.awt.Color(0,150,0))
		lArea4.line=Style.DASHED
		lArea6=LayerXY(chkAreaPy['alpha'],chkAreaPy['PLW'],name='PLW (Py)',color=java.awt.Color(150,0,0))
		lArea4.line=Style.DASHED
		pArea.addLayers([lArea4,lArea5,lArea6])
		pArea.setTitleText('Total areas')
		pArea.setSubtitleText('HIPE vs. Python')
		pArea.setXtitle('Spectral Index')
		pArea.setYtitle('Total Area (sq. arcsec)')

if doApCorr:
	if doCheck:
		print "Checking..."
		fileApCorr='/data/Herschel/Calibration/CalProducts/HIPE11/apcorr_Cal.csv'
		chkApCorr={'alpha':asciiTableReader(file=fileApCorr, tableType='CSV')['c0'].data,
			'PSW_noBG':asciiTableReader(file=fileApCorr, tableType='CSV')['c1'].data,
			'PMW_noBG':asciiTableReader(file=fileApCorr, tableType='CSV')['c2'].data,
			'PLW_noBG':asciiTableReader(file=fileApCorr, tableType='CSV')['c3'].data,
			'PSW_BG':asciiTableReader(file=fileApCorr, tableType='CSV')['c4'].data,
			'PMW_BG':asciiTableReader(file=fileApCorr, tableType='CSV')['c5'].data,
			'PLW_BG':asciiTableReader(file=fileApCorr, tableType='CSV')['c6'].data}

		BSApCorr_noBG=Float2d([[1.277,1.227,1.216],[1.275,1.226,1.214],[1.274,1.226,1.213],[1.273,1.226,1.212],[1.272,1.225,1.21],[1.271,1.225,1.209],[1.27,1.224,1.208],[1.269,1.224,1.207],[1.269,1.223,1.205],[1.268,1.223,1.204],[1.267,1.222,1.203],[1.267,1.222,1.202],[1.266,1.222,1.202],[1.265,1.221,1.2],[1.264,1.221,1.199],[1.264,1.220,1.198],[1.263,1.22,1.197],[1.262,1.219,1.196],[1.262,1.219,1.195]])

		BSApCorr_BG=Float2d([[1.282,1.234,1.236],[1.28,1.234,1.234],[1.279,1.234,1.232],[1.277,1.233,1.231],[1.276,1.233,1.229],[1.276,1.232,1.227],[1.275,1.232,1.226],[1.274,1.231,1.224],[1.273,1.231,1.222],[1.272,1.231,1.221],[1.272,1.23,1.219],[1.271,1.23,1.218],[1.27,1.229,1.216],[1.269,1.229,1.215],[1.269,1.229,1.213],[1.268,1.228,1.212],[1.267,1.228,1.211],[1.266,1.227,1.209],[1.266,1.227,1.208]])
		chkplot={'alpha':[],'PSW':[],'PMW':[],'PLW':[],'PSW_noBG':[],'PMW_noBG':[],'PLW_noBG':[]}
		for n in effBeams:
			if effBeams[n]['type']=='PowerLaw':
				chkplot['alpha'].append(effBeams[n]['alpha'])
				for band in bands:
					chkplot[band].append(effBeams[n][band]['apCorrWithBg'])
					chkplot[band+'_noBG'].append(effBeams[n][band]['apCorrNoBg'])
		pAp1=PlotXY()
		lAp1_1=LayerXY(chkApCorr['alpha'],chkApCorr['PSW_noBG'],name='PSW (py)',color=java.awt.Color(0,0,150))
		lAp1_1.line=Style.DASHED
		lAp1_2=LayerXY(chkApCorr['alpha'],chkApCorr['PMW_noBG'],name='PMW (py)',color=java.awt.Color(0,150,0))
		lAp1_2.line=Style.DASHED
		lAp1_3=LayerXY(chkApCorr['alpha'],chkApCorr['PLW_noBG'],name='PLW (py)',color=java.awt.Color(150,0,0))
		lAp1_3.line=Style.DASHED
		pAp1.addLayers([lAp1_1,lAp1_2,lAp1_3])
		
		lAp1_4=LayerXY(Float1d(chkplot['alpha']),Float1d(chkplot['PSW_noBG']),name='PSW (H)',color=java.awt.Color(0,0,150))
		lAp1_4.line=Style.SOLID
		#lAp1_4.symbol=Style.DCROSS
		lAp1_5=LayerXY(Float1d(chkplot['alpha']),Float1d(chkplot['PMW_noBG']),name='PMW (H)',color=java.awt.Color(0,150,0))
		lAp1_5.line=Style.SOLID
		#lAp1_5.symbol=Style.DCROSS
		lAp1_6=LayerXY(Float1d(chkplot['alpha']),Float1d(chkplot['PLW_noBG']),name='PLW (H)',color=java.awt.Color(150,0,0))
		lAp1_6.line=Style.SOLID
		#lAp1_6.symbol=Style.DCROSS
		pAp1.addLayers([lAp1_4,lAp1_5,lAp1_6])
		pAp1.setXtitle('Spectral Index')
		pAp1.setYtitle('Aperture correction')
		pAp1.setTitleText('Aperture Correction')
		pAp1.setSubtitleText('no Background Annulus')
		pAp1.getLegend().setVisible(1)
		
		
		lAp1_7=LayerXY(chkApCorr['alpha'],BSApCorr_noBG[:,0],name='PSW (BS)',color=java.awt.Color(0,0,150))
		lAp1_7.line=Style.DOTTED
		lAp1_8=LayerXY(chkApCorr['alpha'],BSApCorr_noBG[:,1],name='PMW (BS)',color=java.awt.Color(0,150,0))
		lAp1_8.line=Style.DOTTED
		lAp1_9=LayerXY(chkApCorr['alpha'],BSApCorr_noBG[:,2],name='PLW (BS)',color=java.awt.Color(150,0,0))
		lAp1_9.line=Style.DOTTED
		pAp1.addLayers([lAp1_7,lAp1_8,lAp1_9])

		pAp2=PlotXY()
		pAp2.addLayer(LayerXY(chkApCorr['alpha'],chkApCorr['PSW_BG'],name='PSW',color=java.awt.Color(0,0,150)))
		pAp2.addLayer(LayerXY(chkApCorr['alpha'],chkApCorr['PMW_BG'],name='PMW',color=java.awt.Color(0,150,0)))
		pAp2.addLayer(LayerXY(chkApCorr['alpha'],chkApCorr['PLW_BG'],name='PLW',color=java.awt.Color(150,0,0)))
		pAp2.addLayer(LayerXY(Float1d(chkplot['alpha']),Float1d(chkplot['PSW']),name='PSW (H)',color=java.awt.Color(0,0,96)))
		pAp2.addLayer(LayerXY(Float1d(chkplot['alpha']),Float1d(chkplot['PMW']),name='PMW (H)',color=java.awt.Color(0,96,0)))
		pAp2.addLayer(LayerXY(Float1d(chkplot['alpha']),Float1d(chkplot['PLW']),name='PLW (H)',color=java.awt.Color(96,0,0)))
		chkplot_BG={'alpha':{},'PSW':{},'PMW':{},'PLW':{}}
		pAp2.setXtitle('Spectral Index')
		pAp2.setYtitle('Aperture correction')
		pAp2.setTitleText('Aperture Correction')
		pAp2.setSubtitleText('with Background Annulus')
		pAp2.getLegend().setVisible(1)

		pAp2.addLayer(LayerXY(chkApCorr['alpha'],BSApCorr_BG[:,0],name='PSW (BS)'))
		pAp2.addLayer(LayerXY(chkApCorr['alpha'],BSApCorr_BG[:,1],name='PMW (BS)'))
		pAp2.addLayer(LayerXY(chkApCorr['alpha'],BSApCorr_BG[:,2],name='PLW (BS)'))
