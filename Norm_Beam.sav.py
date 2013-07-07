####
## Normalise beam maps to peak at 1 and have background at zero
####

bands=["psw","pmw","plw"]
inDir="/data/Herschel/Calibration/"

radArcsec=[22,30,42]

files=[]
for b in range(3):
	files.append("%sNepBeam_%s_1arcsec_combined.fits"%(inDir,bands[b]))

mapsIn=[]
mapsNorm=[]
for b in range(3):
	mapIn=simpleFitsReader(files[b])
	mapsIn.append(mapIn)
	srcParam=sourceFitting(image=mapsIn[b], \
		minX=1000-radArcsec[b], minY=1000-radArcsec[b], \
		width=2*radArcsec[b],height=2*radArcsec[b])
	srcPeak=srcParam["Column1"].data[0]
	print '%s: %f'%(bands[b],srcPeak)
	
	DataNAN=IS_NAN(mapsIn[b].image)
	mapsIn[b].image[mapsIn[b].image.where(DataNAN)]=0.
	
	mapNorm=SimpleImage()
	mapNorm.setImage(mapIn.image/srcPeak)
	mapsNorm.append(mapNorm)


#for pix in mapInPsw.image:
#	if IS_NAN(pix):
#		pix=0.
