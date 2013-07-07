#bendoSourceFit.py - A HIPE class for performing 2D Gaussian fitting to
#unresolved sources in PointedPhotTimelines in a SpireListContext.
#
#Written by George J. Bendo in 2010 with a little assistance from Andreas
#Papageorgiou and Ivan Valtchanov.
#
#Usage (basic):
#
# > fitter=bendoSourceFit(inputContext)
# > param=fitter.fit(array, raCent, decCent, rPeak)
# > param=fitter.fit(array, raCent, decCent, rPeak, rBackIn, rBackOut)
#
#Usage (using optional methods):
#
# > fitter=bendoSourceFit(inputContext)
# > fitter.setMaxIterations(10000)
# > fitter.setModelGauss2D()
# > fitter.setModelGauss2DFixedWidth(0.002069)
# > fitter.setModelGauss2DRot()
# > fitter.setPlot(0)
# > fitter.setPlotRegionCenter()
# > fitter.setPlotFitCenter()
# > fitter.setTiltBack(Boolean.FALSE)
# > fitter.setTolerance(1e-4)
# > fitter.setUseBackInFit(Boolean.TRUE)
# > fitter.setVaryBackInFit(Boolean.TRUE)
# > param=fitter.fit(array, raCent, decCent, rPeak, rBackIn, rBackOut)
# > chiSquared=fitter.getChiSquared()
# > evidence=fitter.getEvidence()
# > nPointBack=fitter.getNPointBack()
# > nPointFit=fitter.getNPointFit()
# > reducedChiSquared=fitter.getReducedChiSquared()
# > paramErr=fitter.getStandardDeviation()
#
#Input:
#  inputContext: A SpireListContext with PointedPhotTimelines.  These
#    data should be processed through a baseline removal tool before this
#    class is used.
#  array: The array (PSW, PMW, or PLW) used for the fit.
#  raCent: The RA (in degrees) of the center of the region to use in the fit.  
#    (The class will derive the center of the soruce.)
#  decCent: The declination (in degrees) of the center of the region to use in 
#    the fit.  (The class will derive the center of the soruce.)
#  rPeak: The radius of the region that will include the peak of the source, in
#    arcsec.  Suggested values are 22, 30, and 42 for PSW, PMW, and PLW,
#    respectively.  (These radii exclude the Airy rings.)
#
#Optional input:
#  rBackIn: Optional inner radius of the annular region to use as the 
#    background, in arcsec.  When provided, a background is measured within the
#    background and subtracted before the fit is performed.  Both rBackIn and 
#    rBackOut must be provided for the background annulus to be used.  
#    Suggested value is 300. 
#  rBackOut: Optional outer radius of the annular region to use as the 
#    background, in arcsec.  When provided, a background is measured within the
#    background and subtracted before the fit is performed.  Both rBackIn and 
#    rBackOut must be provided for the background annulus to be used.  
#    Suggested value is 350.
#
#Output (from the fit method): A Double1d array with the following numbers:
#  If the model used is Gauss2dModel:
#    [0]: Amplitude
#    [1]: RA of central position (degrees)
#    [2]: Dec of central position (degrees)
#    [3]: Sigma (degrees)
#    [4]: Background constant term (if the background is treated as a free 
#         parameter)
#    [5]: Slope in background in declination (if the tilted background is used;
#         signal units/degrees)
#    [6]: Slope in background in right ascension (if the tilted background is 
#         used; signal units/degrees)
#  If the model used is Gauss2dModel with a fixed width:
#    [0]: Amplitude
#    [1]: RA of central position (degrees)
#    [2]: Dec of central position (degrees)
#    [3]: Background constant term (if the background is treated as a free 
#         parameter)
#    [4]: Slope in background in declination (if the tilted background is used;
#         signal units/degrees)
#    [5]: Slope in background in right ascension (if the tilted background is 
#         used; signal units/degrees)
#  If the model used is Gauss2DRotModel:
#    [0]: Amplitude
#    [1]: RA of central position (degrees)
#    [2]: Dec of central position (degrees)
#    [3]: Sigma of major axis (degrees)
#    [4]: Sigma of minor axis (degrees)
#    [5]: Position angle of major axis (radians)
#    [6]: Background constant term (if the background is treated as a free 
#         parameter)
#    [7]: Slope in background in declination (if the tilted background is used;
#         signal units/degrees)
#    [8]: Slope in background in right ascension (if the tilted background is 
#         used; signal units/degrees)
#
#Optional methods.
#  getChiSquared(): Returns the chi squared value for the fit to the data.
#  getEvidence(): Returns the evidence (log Z) value for the fit.  (See the
#    method below for a more verbose description of what this is.)
#  getNPointBack(): Returns the number of data points in the background
#    annulus.
#  getNPointFit(): Returns the number of data points used in the fit.  If the
#    background annulus is included in the fit, then the number given by this
#    method will include the number of data points in the background.
#  getReducedChiSquared(): Returns the reduced chi squared value for the fit 
#    to the data.
#  getStandardDeviation(): Returns the standard deviation for the best-fitting
#    parameters from the last fit performed (in the same format as the 
#    parameters themselves).
#  setMaxInterations(): Sets the maximum number of iterations performed
#    by the LevenbergMarquardtFitter. The default is 10000.
#  setModelGauss2D(): Sets the function to be fitted to the data as a circular
#    2D Gaussian function.  The default is the use an elliptical 2D Gaussian
#    function.
#  setModelGauss2DFixedWidth(sigma): Sets the function to be fitted to the
#    data as a circular 2D Gaussian function with a fixed width sigma 
#    (in degrees) that is specified by the user.   The default is the use 
#    an elliptical 2D Gaussian function (with a variable width).
#  setModelGauss2DRot(): Sets the function to be fitted to the data as a 
#    elliptical 2D Gaussian function.  This is the default.
#  setPlot(input): Sets whether a radial profile is plotted at the end of the 
#    execution of the fitter method.  The options are 0 (no plot), 1 (a plot 
#    with the center set to the center of the target aperture), or 2 (a plot 
#    with the center set to the position given by the fit coordinates).  The 
#    default is not to produce a plot (0).
#  setPlotRegionCenter(): Sets whether a radial profile is plotted at the end 
#    of the execution of the fitter method  with the center position set to 
#    the center of the target aperture.  The default is to produce no plot.
#  setPlotFitCenter(): Sets whether a radial profile is plotted at the end 
#    of the execution of the fitter method  with the center position set to 
#    the center of the Gaussian function from the fit.  The plot also includes
#    the radial profile of the best fitting function.  The default is to 
#    produce no plot.
#  setTiltBack(): Sets whether a tilted plane is used for the background.
#    Setting this to true will set varyBackInFit to true.  The default is
#    false.  Note that tiltBack will bet set to false if varyBackInFit is
#    set to false.
#  setTolerance(value): Sets the tolerance that the LevenbergMarquardtFitter
#    attempts to achieve when performing the fit.  The default is 1e-4.
#  setUseBackInFit(booleanValue): Sets whether the background annulus is
#    included in the data that is input into the fitter.  Note that this
#    is ignored if the background annulus is not specified.  The default is 
#    true.
#  setVaryBackInFit(booleanValue): Sets whether the background is treated as
#    a free parameter in the fit.  The default is true.  This is set to true
#    if the background is set to tilt (i.e. tiltBack=true). 
#
#Notes:
#  * Although this class includes methods for returning chi squared, reduced 
#    chi squared, and a Bayesian evidence value, tests with simulated sources 
#    added to real timeline data showed that these numbers did not necessarily 
#    reflect the accuracy of the resulting fit.  For example, the evidence
#    values for a Gauss2DModel with a fixed width fit to the simulated source
#    was systematically higher (e.g. better) than fit, although the standard
#    Gauss2DModel fit (with the variable width) and the Gauss2DRotModel fit
#    produced more accurate measurements.  Also keep in mind that so many data
#    points are used in the fit that the reduced chi squared values will vary 
#    very little between models with different numbers of parameters. Use 
#    these metrics with extreme caution.
#  * While this class can produce radial profiles of the best fitting functions,
#    these radial profiles may be difficult to interpret if the data are
#    asymmetric (e.g. the source is elliptical or the background is sloped).
#    If the best-fitting radial profile is plotted, it will be an average
#    best-fitting function.
#
#Development history:
#V0-1 (Aug 2010): Created class.
#V0-2 (Aug 2010): Added information in comments on output based on comments
#    from A. Smith.  Added methods for selecting between elliptical and 
#    circular Gaussian functions and for setting the number of iterations.
#V0-3 (Aug 2010): Made small fixes to the part of the class that reads in
#    timeline data.  The fitter's speed improved significantly.  Altered
#    the way that different models are selected for fitting so as to 
#    avoid chashes on different HIPE builds.  Added lines to handle standard
#    deviations in parameters from LevenbergMarquardtFitter (based on
#    information from I. Valtchanov).
#V0-4 (Sep 2010): Added ability to perform fit without using a background
#    annulus.  Added the ability to select a tilted plane for the fit.
#    Added the ability to fix the width of the Gaussian function (based on
#    information from I. Valtchanov).  Fixed bug related to Boolean flags.  
#    Increased the default value for the tolerance, which will have a 0.05% 
#    effect on measured peak flux densities.  Replaced numerical values of pi 
#    with Math.PI.  Added description about getStandardDeviation in the
#    comments.
#V0-5 (Oct 2010): Added some memory-saving lines.  Added a check to determine
#    if the specified array is valid; if not, then an error message is  
#    displayed.  Added a check to determine whether any data fall within the 
#    specified aperture for the fit; if not, then an error message is 
#    displayed.  Added methods to return chi squared, reduced chi squared, 
#    and evidence values along with a warning note.
#V0-6 (Feb 2011): Revised the checks for mask bits so that the metadata
#    for the input are checked for the correct mask bit values.  Also added
#    checks for uncorrected glitches and noisy bolometers.  (Why can't we use
#    the %*&$! master bit to track all of these things?)
#V0-7 (Apr 2011): Fixed an issue with determining whether points had RA values
#    that placed them within the target and background apertures.  Added
#    options for plotting.  In description at top of file, added examples
#    of using getChiSquared, getEvidence, and getReducedChiSquared methods.
#V0-8 (Apr 2011): Added lines to import classes needed for this class to
#    work properly (with thanks to I. Valtchanov for finding this issue
#    and explaining how to fix it).  Changed the name of one of the plotting 
#    methods.  Made some minor adjustments to comment formatting.
#V0_9 (May 2011): Added fix to coordinates that affected plotting routines.
#    Added nPointBack and nPointFit methods.  Now rejecting data where
#    maskL1GlitchDetected bit is set.  Added fix to deal with crash
#    when the maskL1GlitchDetected metadata is not in the input.



#Import classes needed for this class.
from herschel.ia.numeric import *
from herschel.ia.numeric.toolbox.basic.Basic import *
from herschel.ia.numeric.toolbox.fit import *
from herschel.ia.gui.plot import *
from java.lang import Boolean, Math

#Define the class.
class bendoSourceFit:

	#Initialize the class.
	def __init__(self,inputContext):

		#Define universal variables.
		self.inputContext=inputContext

		#Set defaults for how the class performs the fit.
		self.fixedsigma=Double1d()
		self.maxIterations=10000
		self.model="Gauss2DRot"
		self.tiltBack=Boolean(0)
		self.tolerance=1e-4
		self.useBackInFit=Boolean(1)
		self.varyBackInFit=Boolean(1)

		#Set other public variables related to the results.
		self.chiSquared=0.
		self.evidence=0.
		self.nPointFit=0l
		self.nPointBack=0l
		self.param=Double1d()
		self.paramErr=Double1d()
		self.plotFit="No Plot"
		self.reducedChiSquared=0.

	#Create a method to get the chi squared value from the fit.
	def getChiSquared(self):
		return self.chiSquared

	#Create a method to get the evidence (log Z) value for the fit.
	#This is a Bayesian measure of the quality of the fit.  The values are 
	#dimensionless and can only be used to compare two model fits relative 
	#to each other.  In a case where two models are fit to the same data
	#and one has an evidence value that is 1 higher than the other, the 
	#model with the higher value is 10 times more likely to represent the 
	#data. For more information, see "A Software Package for Parameter 
	#Estimation and Model Comparison" by Do Kester (in "Bayesian Inference 
	#and Maximum Entropy Methods in Science and Engineering", Eds: R. 
	#Fischer et al., Garching, 2004, AIP Conference Proceedings 735 2004, 
	#p. 379.  Also, read the warning in the notes section.
	def getEvidence(self):
		return self.evidence

	#Create a method for returning the number of points in the background
	#region.
	def getNPointBack(self):
		return self.nPointBack

	#Create a method for returning the number of points used in the fit.
	def getNPointFit(self):
		return self.nPointFit

	#Create a method to get the reduced chi squared value from the fit.
	def getReducedChiSquared(self):
		return self.reducedChiSquared

	#Create a method to get the standard deviation for the best fitting
	#parameters.
	def getStandardDeviation(self):
		return self.paramErr

	#Create a method that will set whether a radial plot is produced.  A
	#value of 0 produces no plot.  A value of 1 uses the center of the
	#target aperture as the point from which the radial profile is 
	#measured.  A value of 2 uses the center of the fit as the point from
	#which the radial profile is measured.
	def setPlot(self,plotFitInput):
		if plotFitInput==0:
			self.plotFit="No Plot"
		if plotFitInput==1:
			self.plotFit="Region Center"
		if plotFitInput==2:
			self.plotFit="Fit Center"
		if plotFitInput!=0 and plotFitInput!=1 and plotFitInput!=2:
			print "Error: The input for this method is invalid (valid numbers"
			print "    are 0, 1, or 2).  The plotting option has not been"
			print "    changed."

	#Create a method that will set whether a radial plot is produced with
	#the center set as the target aperture center.
	def setPlotRegionCenter(self):
		self.plotFit="Region Center"

	#Create a method that will set whether a radial plot is produced with
	#the center set as the center determined from the fit.
	def setPlotFitCenter(self):
		self.plotFit="Fit Center"

	#Create a method to set the maximum number of iterations performed
	#by the fitting program.
	def setMaxIterations(self,maxIterations):
		self.maxIterations=maxIterations

	#Create a method for selecting a circular Gaussian function to use in
	#in the fit.
	def setModelGauss2D(self):
		self.model="Gauss2D"
		self.fixedsigma=Double1d()

	#Create a method to set a Gaussian function with a fixed width (which
	#automatically sets model="Gauss2D".
	def setModelGauss2DFixedWidth(self,sigma):
		self.model="Gauss2D"
		self.fixedsigma=Double1d([sigma])

	#Create a method for selecting an ellipsoidal Gaussian function to use
	#in the fit.
	def setModelGauss2DRot(self):
		self.model="Gauss2DRot"
		self.fixedsigma=Double1d()

	#Set whether the background is allowed to be a tilted plane in the fit.
	def setTiltBack(self,booleanVal):
		self.tiltBack=Boolean(booleanVal)
		if Boolean(booleanVal)==Boolean(1):
			self.varyBackInFit=Boolean(Boolean.TRUE)

	#Set the tolerance level to be reached by the fitting program.
	def setTolerance(self,tolerance):
		self.tolerance=tolerance

	#Create a method for changing whether the background annulus is used in
	#the fit.
	def setUseBackInFit(self,booleanVal):
		self.useBackInFit=Boolean(booleanVal)

	#Create a method for changing whether the background is varied in the
	#fit.
	def setVaryBackInFit(self,booleanVal):
		self.varyBackInFit=Boolean(booleanVal)
		if Boolean(booleanVal)==Boolean(0):
			self.tiltBack=Boolean(Boolean.FALSE)

	#Determine the parameters of the best fitting function using a 
	#method.  This version uses a background annulus.
	def fit(self,array,raCent,decCent,rPeak,*rBack):

		#Check to see if array is equal to either "PSW", "PMW", or
		#"PLW".  If not, print an error message and exit.
		if array!="PSW" and array!="PMW" and array!="PLW":
			print "Error: An invalid array was given in the input.  Please"
			print "    use either PSW, PMW, or PLW. No fit will be"
			print "    performed."
			return 

		#Determine if rBackIn and rBackOut are given.  If not, then 
		#set backRegionFlag to 0.  If it is set, set backRegionFlag
		#to 1.
		backRegionFlag=0
		if len(rBack)==2:
			backRegionFlag=1

		#Set up storage arrays for the signal and coordinate data.
		raFit=Double1d()
		decFit=Double1d()
		signalFit=Double1d()
		raBack=Double1d()
		decBack=Double1d()
		signalBack=Double1d()

		#Loop through each Product in the Context and each bolometer.
		for ref in self.inputContext.getAllRefs():
			product=ref.getProduct()

			#Determine the mask bits to check.  In cases where
			#these mask bits are set, the data will not be used.
			#For now, we are not using data where the master bit,
			#truncated, noisy, slow, and uncorrected L1 glitches
			#are not removed.  Since not all data include 
			#"maskGlitchL1NotRemoved" mask bits, use a try statement
			#to identify this bit, and substitute the master bit
			#for the value if it does not work.
			mc1=product.meta["maskMaster"].value
			mc2=product.meta["maskTruncated"].value
			mc3=product.meta["maskNoisy"].value
			mc4=product.meta["maskSlow"].value
			try:
				mc5=product.meta["maskGlitchL1Detected"].value
			except:
				mc5=mc1
				print 'Warning: Problem identifying data with glitches because'
				print '    metadata in the input does not identify the mask bit.'
				print '    Uncorrected glitches may be included in the fit.'

			#Reset the variables for counting the number of
			#points used in the fit and in the background.
			self.nPointFit=0l
			self.nPointBack=0l

			for bolom in product.getChannelNames():
				if bolom[0:3]==array and bolom[3]!='T' and bolom[3]!='R' and bolom[3:5]!='DP':

					#Extract the coordinates and mask from 
					#the level 1 Product.
					signalBolom=Double1d(product["signal"].getColumn(bolom).data)
					raBolom=product["ra"].getColumn(bolom).data
					decBolom=product["dec"].getColumn(bolom).data
					maskBolom=product["mask"].getColumn(bolom).data

					#Select only data that fall within the 
					#region specified by the user and that 
					#are not masked out.  
					r=( ((raBolom-raCent)*COS(decCent/180*Math.PI))**2. + (decBolom-decCent)**2. )**0.5*3600.
					if MIN(r)<rPeak:
						iFit=r.where((r<rPeak)\
							.and((maskBolom/mc1)%2==0)\
							.and((maskBolom/mc2)%2==0)\
							.and((maskBolom/mc3)%2==0)\
							.and((maskBolom/mc4)%2==0)\
							.and((maskBolom/mc5)%2==0)\
							.and(IS_FINITE(raBolom))\
							.and(IS_FINITE(decBolom))\
							.and(IS_FINITE(signalBolom)))
					else:
						iFit=Double1d()
					if (backRegionFlag==1) and (MIN(r)<rBack[1]):
						iBack=r.where((r>rBack[0]).and(r<rBack[1])\
							.and((maskBolom/mc1)%2==0)\
							.and((maskBolom/mc2)%2==0)\
							.and((maskBolom/mc3)%2==0)\
							.and((maskBolom/mc4)%2==0)\
							.and((maskBolom/mc5)%2==0)\
							.and(IS_FINITE(raBolom))\
							.and(IS_FINITE(decBolom))\
							.and(IS_FINITE(signalBolom)))
					else:
						iBack=Double1d()

					#Store the coordinate data and signal 
					#data that meet the above criteria.
					if iFit.length()>0:
						raFit.append(raBolom[iFit])
						decFit.append(decBolom[iFit])
						signalFit.append(signalBolom[iFit])
						self.nPointFit+=iFit.length()
					if iBack.length()>0:
						raBack.append(raBolom[iBack])
						decBack.append(decBolom[iBack])
						signalBack.append(signalBolom[iBack])
						self.nPointBack+=iBack.length()

			#Erase product, which otherwise might hog memory.
			del(product)

		#If no data fell within the fit radius, display an error
		#message and return.
		if signalFit.length()==0:
			print "Error: The aperture defined for the source contains no data."
			print "    Please check your input coordinates and radius. No fit"
			print "    will be performed."
			return 

		#Perform the next steps only if signalBack contains data.
		subBack=0.
		if len(signalBack)>0:

			#If useBackInFit is set, append background data to fit 
			#data.
			if self.useBackInFit==Boolean(1):
				raFit.append(raBack)
				decFit.append(decBack)
				signalFit.append(signalBack)
				self.nPointFit+=self.nPointBack

			#Subtract a preliminary background from the data that 
			#will be fit with the 2D Gaussian function.
			subBack=MEDIAN(signalBack)
			signalFit-=subBack
			del(raBack,decBack,signalBack)

		#Set up a model. We assume that the deprojection will set the 
		#source close to raCent and decCent, so the initial guess for 
		#the central coordinates is assumed to be 0.  The width 
		#parameters are in degrees and are within an order of magnitude 
		#for all of the SPIRE photometer arrays.  The final value is 
		#the position angle of the major axis for Gauss2DRotModel is 
		#IN RADIANS AND NOT DEGREES.  This angle is given a non-zero 
		#value just to avoid unforseen issues with the functionality 
		#of the nonlinear inverse analysis.  If varyBackInFit is set, 
		#add the background model.  Note that BinomialModel(0) is 
		#effectively a two-dimensional model of a constant value.
		if self.model=="Gauss2D":
			model=Gauss2DModel()
			if self.varyBackInFit==Boolean(1):
				if self.tiltBack==Boolean(1):
					model+=PolySurfaceModel(1)
					model.setParameters(Double1d([MAX(signalFit),0.00001,0.00001,0.0021,1e-10,0.001,-0.001]))
				else:
					model+=BinomialModel(0)
					model.setParameters(Double1d([MAX(signalFit),0.00001,0.00001,0.0021,1e-10]))
			else:
				model.setParameters(Double1d([MAX(signalFit),0.00001,0.00001,0.0021]))
			if len(self.fixedsigma)==1:
				param=model.getParameters()
				param[3]=self.fixedsigma[0] 
				model.setParameters(param)
				model.keepFixed(Int1d([3]))
		if self.model=="Gauss2DRot":
			model=Gauss2DRotModel()
			if self.varyBackInFit==Boolean(1):
				if self.tiltBack==Boolean(1):
					model+=PolySurfaceModel(1)
					model.setParameters(Double1d([MAX(signalFit),0.00001,0.00001,0.0021,0.0019,0.01,1e-10,0.001,-0.001]))
				else:
					model+=BinomialModel(0)
					model.setParameters(Double1d([MAX(signalFit),0.00001,0.00001,0.0021,0.0019,0.01,1e-10]))
			else:
				model.setParameters(Double1d([MAX(signalFit),0.00001,0.00001,0.0021,0.0019,0.01]))

		#Deproject the coordinates.
		coordinates=Double2d()
		coordinates.append((raFit-raCent)*COS(decCent/180.*Math.PI),1)
		coordinates.append(decFit-decCent,1)
		#del(raFit,decFit)

		#Perform the fit.
		fitter=LevenbergMarquardtFitter(coordinates,model)
		fitter.setMaxIterations(self.maxIterations)
		fitter.setTolerance(self.tolerance)
		self.param=fitter.fit(signalFit)
		self.paramErr=fitter.getStandardDeviation()
		self.chiSquared=fitter.getChiSquared()
		self.reducedChiSquared=self.chiSquared/(len(signalFit)-len(self.param))
		self.evidence=fitter.getEvidence(Double1d(model.getNumberOfParameters(),1000))

		#Add the offset back to the parameters for the center of the 
		#Gaussian function. 
		self.param[1]=self.param[1]/COS(decCent/180.*Math.PI)+raCent
		self.paramErr[1]/=COS(decCent/180.*Math.PI)
		self.param[2]=self.param[2]+decCent

		#Add subBack to the background parameter if the background was 
		#treated as a free parameter.  The value that is returned will 
		#be equivalent to the real background and not just the residual
		#from the fit.
		if self.varyBackInFit==Boolean(1):
			if self.model=="Gauss2D" and len(self.fixedsigma)==0:
				self.param[4]+=subBack
			if self.model=="Gauss2D" and len(self.fixedsigma)==1:
				self.param[3]+=subBack
			if self.model=="Gauss2DRot":
				self.param[6]+=subBack

		#Plot the radial profile if self.plotFit="Fit Center".  This
		#will include a radial profile of the fit as well as the data.
		if self.plotFit=="Fit Center":
			rPlot=( ((raFit-self.param[1])*COS(self.param[2]/180.*Math.PI))**2. + (decFit-self.param[2])**2. )**0.5*3600.
			p=PlotXY(rPlot,signalFit)
			p.setLine(0)
			p.setSymbolSize(2)
			p.setSymbol(5)
			p.setXtitle("Radius (arcsec)")
			p.setYtitle("Signal")
			rProfile=Double1d.range(FIX(CEIL(MAX(rPlot)))[0])
			if self.model=="Gauss2D" and len(self.fixedsigma)==0:
				profileFit=self.param[0]*EXP(-rProfile**2./2./(self.param[3]*3600.)**2.)
				if self.varyBackInFit==Boolean(1):
					profileFit+=self.param[4]
					if self.tiltBack==Boolean(1):
						profileFit+=self.param[5]*(self.param[2]-decCent)+self.param[6]*(self.param[1]-raCent)*COS(decCent/180.*Math.PI)
			if self.model=="Gauss2D" and len(self.fixedsigma)==1:
				profileFit=self.param[0]*EXP(-rProfile**2./2./(self.fixedsigma[0]*3600.)**2.)
				if self.varyBackInFit==Boolean(1):
					profileFit+=self.param[3]
					if self.tiltBack==Boolean(1):
						profileFit+=param[4]*(self.param[2]-decCent)+self.param[5]*(self.param[1]-raCent)*COS(decCent/180.*Math.PI)
			if self.model=="Gauss2DRot":
				profileFit=self.param[0]*EXP(-rProfile**2./2./(self.param[3]*self.param[4]*3600.**2.))
				if self.varyBackInFit==Boolean(1):
					profileFit+=self.param[6]
					if self.tiltBack==Boolean(1):
						profileFit+=self.param[7]*(self.param[2]-decCent)+self.param[8]*(self.param[1]-raCent)*COS(decCent/180.*Math.PI)
			p.addLayer(LayerXY(rProfile,profileFit))
			p[1].setStroke(2)

		#Plot the radial profile if self.plotFit="Region Center".
		if self.plotFit=="Region Center":
			rPlot=( ((raFit-raCent)*COS(decCent/180.*Math.PI))**2. + (decFit-decCent)**2. )**0.5*3600.
			p=PlotXY(rPlot,signalFit)
			p.setLine(0)
			p.setSymbolSize(2)
			p.setSymbol(5)
			p.setXtitle("Radius (arcsec)")
			p.setYtitle("Signal")

		#Return param.
		return self.param
