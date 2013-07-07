
beamDir ='/home/astrog82/spxcen/Herschel/Calibration/spire_beams_varpix/'
beamDir0='/home/astrog82/spxcen/Herschel/Calibration/spire_beams_measured/'

obsId='5000241a'
pixSize1=[1,1,1]
pixSize2=[6,10,14]
bands=['Psw','Pmw','Plw']
bands0=['psw','pmw','plw']

files0=[]
files1=[]
files2=[]

for b in range(3):
    files0.append('%s%s_beam_1arcsec.fits'%(beamDir0,bands0[b]))
    files1.append('/home/astrog82/spxcen/Herschel/Calibration/spire_beams_varpix/%s_Beams_multiPixelSize/%s_%sBeam_pxSize_%.1f.fits'%(obsId,obsId,bands[b],pixSize1[b]))
    files2.append('%s%s_Beams_multiPixelSize/%s_%sBeam_pxSize_%.1f.fits'%(beamDir,obsId,obsId,bands[b],pixSize2[b]))

maps0=[]
maps1=[]
maps2=[]
rad0=[]
rad1=[]
rad2=[]

for b in range(1):
	#print 'Map 0:'
	#print files0[b]
	#map0=simpleFitsReader(files0[b],reader=herschel.ia.toolbox.util.SimpleFitsReaderTask.ReaderType.STANDARD)
	#cpx0Est=map0.meta["crpix1"].value
	#cpy0Est=map0.meta["crpix2"].value
	#srcFit = sourceFitting(image=map0["PrimaryImage"], minX=cpx0Est-50., minY=cpy0Est-50., width=100., height=100.)
	#cpx0=srcFit["Column1"].data[1]
	#cpy0=srcFit["Column1"].data[2]
	#print map0[int(cpx0),int(cpy0)]
	#print max(map0)

	print files1[b]
	map1=simpleFitsReader(files1[b])
	nx1=map1.meta["naxis1"].value
	ny1=map1.meta["naxis2"].value
	cpx1Est=map1.meta["crpix1"].value - 50.
	cpy1Est=map1.meta["crpix2"].value - 50.
	srcFit = sourceFitting(image=map1, minX=cpx1Est, minY=cpy1Est, width=100., height=100)
	cpx1=srcFit["Column1"].data[1]
	cpy1=srcFit["Column1"].data[2]
	#pixSize1=ABS(map1.meta["cdelt1"].value) * 3600. #in arcsec
	print max(map1.image)
	print map1.image[int(cpx1),int(cpy1)]
	rad1=Float2d(nx1,ny1)
	print 'Making rad:'
	for x in range(nx1):
		for y in range(ny1):
			rad1[x,y]=SQRT((float(x)-cpx1)**2 + (float(y)-cpy1)**2)
	
			#print x,cpx1,float(x)-cpx1,y,cpy1,float(y)-cpy1,rad1[x,y]

	#result = annularSkyAperturePhotometry(image=map1, centerX=cpx1, centerY=cpx2, fractional=1, radiusArcsec=1.0)


	print files2[b]
	map2=simpleFitsReader(files2[b])
	cpx2Est=map2.meta["crpix1"].value
	cpy2Est=map2.meta["crpix2"].value
	srcFit = sourceFitting(image=map2, minX=cpx2Est-50., minY=cpy2Est-50., width=100, height=100)
	cpx2=srcFit["Column1"].data[1]
	cpy2=srcFit["Column1"].data[2]
	pixSize2=abs(map2.meta["cdelt1"].value) * 3600. #in arcsec

	#for 
