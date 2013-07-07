from herschel.ia.pal.pool.http import HttpClientPool

pool = HttpClientPool ("http://wakefield.bnsc.rl.ac.uk/hcss/pal", "ops")

#pool.setAuthentication(Configuration.getProperty("hcss.access.authentication.authstring"))
store = ProductStorage (pool)    # Just the Data
store.authenticate()

#************* EDIT THESE TWO PARAMETERS******************

# Storage where observations will be saved (within lstore)
#store_out=ProductStorage("Neptune_0x5000E665")
store_out=ProductStorage("Neptune_0x50002582")


# OBSIDs of observations to be downloaded
#obsids=[0x5000E665]
obsids=[0x50002582]


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
