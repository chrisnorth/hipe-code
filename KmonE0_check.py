calibration = spireCal()
#rsrf = calibration.phot.getProduct('Rsrf')
#apEff=calibration.phot.getProduct('ApertureEfficiency')

freq=calibration.phot.getProduct('Rsrf').getFrequency() * 1.e9 #in Hz

arrays=['PSW','PMW','PLW']

rsrf={}
apEff={}
for array in arrays:
	rsrf[array] = calibration.phot.getProduct('Rsrf').getRsrfColumn(array).data
	apEff[array] = calibration.phot.getProduct('ApertureEfficiency').getApertEffTable().getColumn(array).data

#make combinations with interpolators
rsrfOnly={}
rsrfOmega={}
rsrfApEff={}
rsrfApEffOmega={}
for array in arrays:
	rsrfOnly[array] = CubicSplineInterpolator(freq,rsrf[array])
	rsrfOmega[array] = CubicSplineInterpolator(freq,rsrf[array]*freq**-1.7)
	rsrfApEff[array] = CubicSplineInterpolator(freq,rsrf[array]*apEff[array])
	rsrfApEffOmega[array] = CubicSplineInterpolator(freq,rsrf[array]*apEff[array]*freq**-1.7)

integ = TrapezoidalIntegrator(min(freq),max(freq))
KmonE0={}
for array in arrays:
	KmonE0[array] = integ.integrate(rsrfOmega[array]) * integ.integrate(rsrfApEff[array]) / \
	(integ.integrate(rsrfApEffOmega[array]) * integ.integrate(rsrfOnly[array]))
	print 'KmonE(0) %s: %f (%f)'%(array,KmonE0[array],1-KmonE0[array])

comment = """First, let's deal with the point source

For a point source with spectral index a, KmonP(a) is the conversion from the broad-band RSRF-weighted flux density to the monochromatic flux density (at nu=nu0):

[Eq.1] KmonP(a) = int{ RSRF(nu) * ApEff(nu) dnu} / int{ RSRF(nu) * ApEff(nu) * (nu/nu0)^(a) dnu}

K4P, used in the pipeline for a=-1, is therefore just KmonP(a=-1)

[Eq.2] K4P = int{ RSRF(nu) * ApEff(nu) dnu} / int{ RSRF(nu) * ApEff(nu) * nu^(-1) dnu}

The colour correction for a source of spectrum alpha, KcolP(a), converts from the monochromatic flux density assuming a=-1, to the same for spectral index a:

[Eq.3] KcolP(a) = KmonP(a) / K4P

From [Eq.1] it is clear that KmonP(0)=1, and so:

[Eq.4] KcolP(0) = 1 / K4P

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For an extended source with spectral index a, KmonP(a) is the conversion from the broad-band RSRF-weighted flux density to the monochromatic flux density (at nu=nu0):

[Eq.5] KmonE(a) = Omega_eff(a) * int{ RSRF(nu) * ApEff(nu) dnu} / int{ RSRF(nu) * ApEff(nu) * (nu/nu0)^a * Omega(nu) dnu}

where Omega(nu) is the monochromatic beam area at nu, while Omega_eff(a) is the effective area for spectrum a:

[Eq.6] Omega_eff(a) = int{ RSRF(nu) * (nu/nu0)^a * Omega(nu) dnu} / int{ RSRF(nu) * (nu/nu0)^a }

Note that Eq.6 does not include the aperture efficiency, and so KmonE(0) != 1 (though it is close). The values of KmonE(0) are:

[Eq.7] KmonE(0) = [ int{RSRF(nu) * Omega(nu) dnu} * int{ RSRF(nu) * ApEff(nu) dnu} ] / [ int{RSRF(nu) * ApEff(nu) * Omega(nu) dnu} * int{RSRF(nu)} ]

Assuming that Omega(nu) ~ nu^-1.7 (which is a good approximation, then by integrating the appropriate combinations of RSRF and ApEff in HIPE, I get:
"""