#-------------------------------------------------------------------------------
## Check dependencies
#overrideDepCheck = True #set to override dependency check
#depCheck=True
#if calcRadialCorrBeam:
#	if not calcSpireEffFreq:
#		depCheck=False
#		print "***WARNING: new RadialCorrBeam will change effective frequencies***"
#	if not calcColorCorrK:
#		depCheck=False
#		print "***WARNING: new RadialCorrBeam will change ColorCorrK_extended and ColorCorrBeam***"
#	if not calcColorCorrHfi:
#		depCheck=False
#		print "***WARNING: new RadialCorrBeam will change ColorCorrHfi***"
#	if not calcColorCorrAperture:
#		depCheck=False
#		print "***WARNING: new RadialCorrBeam will change ColorCorrAperture***"
#
#if calcSpireEffFreq:
#	if not calcRadialCorrBeam:
#		depCheck=False
#		print "***WARNING: new effective frequencies will change RadialCorrBeam metadata***"
#	if not calcColorCorrK:
#		depCheck=False
#		print "***WARNING: new effective frequencies will change ColorCorrK_extended and ColorCorrBeam***"
#	if not calcColorCorrHfi:
#		depCheck=False
#		print "***WARNING: new effective frequencies will change ColorCorrHfi***"
#	if not calcColorCorrAperture:
#		depCheck=False
#		print "***WARNING: new effective frequencies will change ColorCorrAperture***"
#
#if calcColorCorrK:
#	if not calcColorCorrHfi:
#		depCheck=False
#		print "***WARNING: new ColorCorrBeam & ColorCorrK will change ColorCorrHfi***"
#	if not calcRadialCorrBeam:
#		depCheck=False
#		print "***WARNING: new ColorCorrBeam will change RadialCorrBeam metadata***"
#	if not calcColorCorrAperture:
#		depCheck=False
#		print "***WARNING: new ColorCorrBeam will change ColorCorrAperture***"
#
##stop if not overriden
#if not overrideDepCheck:
#	assert depCheck==True,"***ERROR: Update dependencies not met. Stopping."
#else:
#	if depCheck==False:
#		print "***WARNING: Update dependencies not met. Overriding."
