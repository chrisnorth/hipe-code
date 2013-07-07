
from herschel.ia.pal.pool.http import HttpClientPool

pool = HttpClientPool ("http://wakefield.bnsc.rl.ac.uk/hcss/pal", "ops")

#pool.setAuthentication(Configuration.getProperty("hcss.access.authentication.authstring"))
store = ProductStorage (pool)    # Just the Data
store.authenticate()

#************* EDIT THESE TWO PARAMETERS******************

# Storage where observations will be saved (within lstore)
#store_out=ProductStorage("Neptune_0x5000E665")

# OBSIDs of observations to be downloaded
#obsids=[0x5000E665]

##Neptune
store_out=ProductStorage("Neptune")
#obsids=[0x50002582,0x500048E1,0x50004BC5,0x50004BC7,0x50004CB1,0x50004CB2,0x50004E65,\
#	0x50004E66,0x50007D19,0x500081E5,0x5000820E,0x5000833C,0x5000A6AC,0x5000A8F4,\
#	0x5000AA4B,0x5000AE55,0x5000B10A]
obsids=[0x5000AE55,0x5000B10A]

##AlphaBoo
#store_out=ProductStorage("AlphaBoo")
#obsids=[0x500057E4,0x500057E5]

##GammaDra
#store_out=ProductStorage("GammaDra")
#obsids=[0x50005983,0x50005984]

#*********************************************************


for obsid in obsids:
	query=MetaQuery(ObservationContext,"p","p.meta['obsid'].value==%iL"%obsid)
	refs=store.select(query)
        obs=refs[0].product
        obs.calibration.spec.refs.clear()
        obs.calibration.phot.refs.clear()
	store_out.save(obs)
pass
print
print 'Download finished!'
print
