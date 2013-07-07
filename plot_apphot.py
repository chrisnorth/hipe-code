from herschel.ia.gui.plot.renderer.PComponentEngine import HAlign

## Plot aperture corrections

dirPath='/data/Herschel/Calibration/RelGains/'

files=String1d(['AlphaBoo_0x500057E4_rg_22.0_30.0_42.0_60.0-90.0.dat',\
	'AlphaBoo_0x500057E5_rg_22.0_30.0_42.0_60.0-90.0.dat',\
	'GammaDra_0x50005983_rg_22.0_30.0_42.0_60.0-90.0.dat',\
	'GammaDra_0x50005984_rg_22.0_30.0_42.0_60.0-90.0.dat',\
	'Neptune_0x5000AE55_rg_22.0_30.0_42.0_60.0-90.0.dat',\
	'Neptune_0x5000B10A_rg_22.0_30.0_42.0_60.0-90.0.dat'])
nf=files.size

names=String1d(nf)
obsidstr=String1d(nf)
obsidi=Int1d(nf)
obs=Float1d.range(nf)+1
bands=Float1d.range(3)+1

radiusArcsecs=Float2d(nf,3)
innerArcsecs=Float2d(nf,3)
outerArcsecs=Float2d(nf,3)

mapFlux=Float2d(nf,3)
mapErr=Float2d(nf,3)

mapFluxrg=Float2d(nf,3)
mapErrrg=Float2d(nf,3)

timeFlux=Float2d(nf,3)
timeErr=Float2d(nf,3)

mapCorrFlux=Float2d(nf,3)
mapCorrErr=Float2d(nf,3)

mapCorrFluxrg=Float2d(nf,3)
mapCorrErrrg=Float2d(nf,3)

srcCorrFlux=Float2d(nf,3)
srcCorrErr=Float2d(nf,3)

srcCorrFluxrg=Float2d(nf,3)
srcCorrErrrg=Float2d(nf,3)

for f in range(nf):
	file=files[f]
	fin=open(dirPath+file,'r')
	lines=String1d(fin.readlines())
	nl=lines.size
	for l in range(nl):
		line=lines[l]
		s0=line.find(':')
		if line.find('Name') >= 0:
			names[f]=line[s0+1:-1]
		if line.find('ObsId') >= 0:
			obsid=int(line[s0+1:-1],16)
			obsidi[f]=obsid
			obsidstr[f]=str(obsid)

		if line.find('Map Flux:') >=0:
			for b in range(3):
				mapFlux[f,b]=float(line[s0+1:-1].split(',')[2*b])
				mapErr[f,b]=float(line[s0+1:-1].split(',')[2*b+1])
		
		if line.find('Map Flux RG:') >=0:
			for b in range(3):
				mapFluxrg[f,b]=float(line[s0+1:-1].split(',')[2*b])
				mapErrrg[f,b]=float(line[s0+1:-1].split(',')[2*b+1])
		
		if line.find('ApCorr Flux:') >=0:
			for b in range(3):
				mapCorrFlux[f,b]=float(line[s0+1:-1].split(',')[2*b])
				mapCorrErr[f,b]=float(line[s0+1:-1].split(',')[2*b+1])
		
		if line.find('ApCorr Flux RG:') >=0:
			for b in range(3):
				mapCorrFluxrg[f,b]=float(line[s0+1:-1].split(',')[2*b])
				mapCorrErrrg[f,b]=float(line[s0+1:-1].split(',')[2*b+1])
		
		if line.find('SrcCorr Flux:') >=0:
			for b in range(3):
				srcCorrFlux[f,b]=float(line[s0+1:-1].split(',')[2*b])
				srcCorrErr[f,b]=float(line[s0+1:-1].split(',')[2*b+1])
		
		if line.find('SrcCorr Flux RG:') >=0:
			for b in range(3):
				srcCorrFluxrg[f,b]=float(line[s0+1:-1].split(',')[2*b])
				srcCorrErrrg[f,b]=float(line[s0+1:-1].split(',')[2*b+1])
		
		if line.find('Timeline Flux:') >=0:
			for b in range(3):
				timeFlux[f,b]=float(line[s0+1:-1].split(',')[2*b])
				timeErr[f,b]=float(line[s0+1:-1].split(',')[2*b+1])

		if line.find('Radius') >=0:
			for b in range(3):
				radiusArcsecs[f,b]=float(line[s0+1:-1].split(',')[b])
		if line.find('inner') >=0:
			for b in range(3):
				innerArcsecs[f,b]=float(line[s0+1:-1].split(',')[b])
		if line.find('outer') >=0:
			for b in range(3):
				outerArcsecs[f,b]=float(line[s0+1:-1].split(',')[b])
		

		filePSW='Obs/mapPSW_'+str(hex(obsid))[2:10]+'.fits'
		filePMW='Obs/mapPMW_'+str(hex(obsid))[2:10]+'.fits'
		filePLW='Obs/mapPLW_'+str(hex(obsid))[2:10]+'.fits'

		#mapPSW = fitsReader(file = dirPath+filePSW)
		#mapPMW = fitsReader(file = dirPath+filePMW)
		#mapPLW = fitsReader(file = dirPath+filePLW)
		
	plot=PlotXY()
	plot.clearLayers()
	plot.autoBoxAxes=1
	ltime=LayerXY(bands-0.5,timeFlux[f,:])
	ltime.line=Style.NONE
	ltime.symbol=Style.DCROSS
	ltime.style.stroke=3.0
	ltime.setErrorY(timeErr[f,:],timeErr[f,:])
	ltime.name='Timeline'
	plot.addLayer(ltime)

	lmap=LayerXY(bands-0.45,mapFlux[f,:])
	lmap.line=Style.NONE
	lmap.symbol=Style.DIAMOND
	lmap.style.stroke=3.0
	lmap.setErrorY(mapErr[f,:],mapErr[f,:])
	lmap.name='Map (uncorr.)'
	plot.addLayer(lmap)
	
	lmaprg=LayerXY(bands-0.4,mapFluxrg[f,:])
	lmaprg.line=Style.NONE
	lmaprg.symbol=Style.DIAMOND
	lmaprg.style.stroke=3.0
	lmaprg.setErrorY(mapErrrg[f,:],mapErrrg[f,:])
	lmaprg.name='Map (uncorr) RG'
	plot.addLayer(lmaprg)
	
	lmapc=LayerXY(bands-0.55,mapCorrFlux[f,:])
	lmapc.line=Style.NONE
	lmapc.symbol=Style.SQUARE
	lmapc.style.stroke=3.0
	lmapc.setErrorY(mapCorrErr[f,:],mapCorrErr[f,:])
	lmapc.name='Map (corr.)'
	plot.addLayer(lmapc)
	
	lmapcrg=LayerXY(bands-0.6,mapCorrFluxrg[f,:])
	lmapcrg.line=Style.NONE
	lmapcrg.symbol=Style.SQUARE
	lmapcrg.style.stroke=3.0
	lmapcrg.setErrorY(mapCorrErrrg[f,:],mapCorrErrrg[f,:])
	lmapcrg.name='Map (corr.) RG'
	plot.addLayer(lmapcrg)
	
	lsrc=LayerXY(bands-0.45,srcCorrFlux[f,:])
	lsrc.line=Style.NONE
	lsrc.symbol=Style.CIRCLE
	lsrc.style.stroke=3.0
	lsrc.setErrorY(srcCorrErr[f,:],srcCorrErr[f,:])
	plot.addLayer(lsrc)
	lsrc.name='Source only'
	
	lsrcrg=LayerXY(bands-0.4,srcCorrFluxrg[f,:])
	lsrcrg.line=Style.NONE
	lsrcrg.symbol=Style.CIRCLE
	lsrcrg.style.stroke=3.0
	lsrcrg.setErrorY(srcCorrErrrg[f,:],srcCorrErrrg[f,:])
	plot.addLayer(lsrcrg)
	lsrcrg.name='Source only (RG)'
	
	plot.legend.visible=1

	plot.xaxis.setRange(0,3)
	plot.xaxis.tick.setFixedValues(Float1d.range(3)+0.5)
	plot.xaxis.tick.label.setFixedStrings(['PSW','PLW','PLW'])
	yr=plot.yaxis.getRange()
	plot.yaxis.setRange(0,yr[1])
	plot.yaxis.setTitleText('Source Flux Density')
	plot.xaxis.setTitleText('Band')
	labname=Annotation(2.9,0.9*yr[1],names[f])
	labobs=Annotation(2.9,0.85*yr[1],str(hex(obsidi[f])))
	labname.setHalign(HAlign.LEFT)
	labobs.setHalign(HAlign.LEFT)
	plot.addAnnotation(labobs)
	plot.addAnnotation(labname)
	
	plotFile='Plots/'+file[:-4]+'.png'
	plot.saveAsPNG(dirPath+plotFile)
