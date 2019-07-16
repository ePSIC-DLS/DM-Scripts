// Script for automated acquisition of STEM images of nanoparticles

// T. J. A. Slater, September 2018

// version:20181002, v0.3
// thomas.slater@diamond.ac.uk

taggroup minxfield, minyfield, maxxfield, maxyfield, pathfield
taggroup imtags, focustags
	
taggroup lowmagfield, limsizefield, lexposurefield, highmagfield, himsizefield, hexposurefield

taggroup c_check, f_check, c_freqfield, c_stepfield, c_scansizefield, c_exposurefield, c_itfield, f_freqfield, frangefield, f_stepsizefield, f_scansizefield, f_exposurefield

object FocusSettingsDialogObject
object ImageSettingsDialogObject

Image CreateHanning(image img)
	{
	/* 
		Function to create a Hanning window on an image.
		Input: Initial image
		Output: Image with Hanning window applied
	*/
	
	number size, sizeX, sizeY, ii
	image hannX, hannY, hann, avg, hannout
	
	GetSize(img, sizeX, sizeY);
	
	// Create Hanning window.
	ii = 1;
	hannX := CreateFloatImage("", sizeX, sizeY);
	hannX = 0;
	hannX[0, 0, 1, sizeX] = 1 - cos( 2 * Pi() * icol / sizeY);
	while( ii < sizeY )
	{
		hannX[ii, 0, 2*ii, sizeX] = hannX[0, 0, ii, sizeY];
		ii = ii * 2;
	}

	ii = 1;
	hannY := CreateFloatImage("", sizeX, sizeY);
	hannY = 0;
	hannY[0, 0, sizeY, 1] = 1 - cos( 2 * Pi() * irow / sizeY);
	while( ii < sizeX )
	{
		hannY[0, ii, sizeY, 2*ii] = hannY[0, 0, sizeY, ii];
		ii = ii * 2;
	}

	hann = hannX * hannY;
	//hann.showimage()
			
	// Subtract average from image.
	avg = img - Average(img);

	// Multiply with Hanning window.
	hannout = avg * hann;
	//hannout.showimage()
	return(hannout)
	}

number fft_coeff_find(image img)
	{
	/* 
		Function to calculate a measure of the intensity and number of 
		spots in an FFT of an atomic resolution image.
		Input: Initial image
		Output: Coefficient equal to number of spots * spot intensity
	*/
	
	image hanning_img
	hanning_img = CreateHanning(img)
	image fft_im = realfft(hanning_img)
	//showimage(fft_im)
	
	imagechangedatatype(fft_im,2)
	//showimage(fft_im)
	
	number sx, sy
	getsize( img, sx, sy )

	image mask := BinaryImage( "Binary Mask", sx, sy )      // The real input for analysis
	Image foundParticlesImage := BinaryImage( "found", sx, sy ) // Accepted particles mask (result)

	image median = medianfilter(fft_im,0,2)

	mask = median>2*mean(median) ? 1 : 0           // This is a simple Treshold to test...
	//ShowImage(mask)
	
	mask=tert(iradius<100,0,mask)
	
	//ShowImage(mask)
	
	string mFields = "CenterX,CenterY"      // Specify analysis results
	Number minParticleSize = 100                  // Limit accepted
	Number doLabel = 0                              // Add labeled mask to inputImage True/False

	imageDocument ResultsDoc = FindParticles(img,mask,mFields, minParticleSize, doLabel, foundParticlesImage)
	
	if(imagedocumentisvalid(ResultsDoc)!=0)
		{
		image x_img :=ImageDocumentGetImage(ResultsDoc,0)
		number imageid=x_img.imagegetid()
		
		//ImageDocumentShow( ResultsDoc )

		number num_p, num_vals
		getsize(x_img, num_vals, num_p)
		
		number intensity = sum(mask*fft_im)
		result("\n Peaks = "+num_p)
		result(", Intensity = "+intensity)
		result(", Coefficient = "+(num_p*intensity))
		return(num_p*intensity)
		}
	else
		{
		result("\nNo peaks found.")
		return(0)
		}

	}
	
number autocorr_coeff_find(image img)
	{
	/* 
		Function to calculate the coefficient of autocorrelation of an image.
		Input: Initial image
		Output: Autocorrelation coefficient
	*/
	
	number fh, fv, ch, cv, vc
	image autoIm = AutoCorrelate(img)
	// return the size of the image
	GetSize(autoIm, fh, fv)
	
	//calculate central pixel
	ch = fh/2
	cv = fv/2
	
	//set central pixel = 0 to get rid of spike in autocorrolation
	autoIm[ch ,cv] = 0
	
	//get max value of pixels surrounding and set central pixel value
	vc = max (autoIm[cv+2, ch-2, cv-2, ch+2])
	
	return(vc)
	}
	
number autofocus_correlation(string coarsetype, number start_step, number auto_stop, number scan_size, number exposure, number max_it)
	{
	/* 
		Function to determine the correct focus of a microscope using 
		autocorrelation.
	*/
	
	number signalIndex = 0 // The signal on channel zero - set in the microscope interface
	number rotation = 0 //currentframes * 90   // degree of rotation
	number continuous = 0 // have continuous or single frame acquisitions
	number lineSynch = 0 // do not use line synching
	number synchronous = 0 // 0 = return immediately, 1 = return when finished
	number selected=1 // acquire signal
	number bitdepth = 4
	number signed = 0
	number auto_val_previous = 0.5
	number auto_val_current
	number step_size = start_step //Set initial focus step
	number it_num = 0
	number diff_auto = 1
	number now_focus, start_focus
	
	Result("\n Auto focus: autocorrelation step")
	
	if(coarsetype == "Fit")
	{
		number i
		number c_fit_it = 10
		
		now_focus = emgetfocus()
		
		image auto_vals:= realimage( "Auto_vals" ,4, c_fit_it, 1 )
		emsetfocus(now_focus-c_fit_it*start_step*0.5)
		start_focus = emgetfocus()
		
		for(i = 0;i < c_fit_it; i++)
			{
			now_focus = emgetfocus()
			image img:= IntegerImage( "Img " ,bitdepth, signed, scan_size, scan_size )
			number imageid=img.imagegetid()
			
			number paramID = DSCreateParameters(scan_size, scan_size, rotation, exposure, linesynch )
			DSSetParametersSignal( paramID, signalIndex, bitdepth, selected, imageID )
			
			DSStartAcquisition(paramID, 0, synchronous);
			while(DSIsAcquisitionActive( )) yield() // wait for acquistion to finish
			
			auto_vals[i,0] = autocorr_coeff_find(img)
			//Result("\n Autocorr value is: " + autocorr_coeff_find(img))
			deleteimage(img)
			
			emsetfocus(now_focus+start_step)
			}
			
		auto_vals.showimage()
		
		// setup fit:
		Image pars := NewImage("pars", 2, 3)
		Image parsToFit := NewImage("pars to fit", 2, 3)
		pars = 10;          // starting values
		parsToFit = 1;
		Number chiSqr = 1e6
		Number conv_cond = 0.00001
		Image errors := auto_vals.ImageClone()
		
		String formulaStr = "p0 + p1*x + p2*x**2"
		Number ok = FitFormula(formulaStr, auto_vals,errors, pars, parsToFit, chiSqr, conv_cond)
		number zdf  = -GetPixel(pars, 1,0) / (2* GetPixel(pars,2,0))
		
		Image plot := PlotFormula(formulaStr, auto_vals, pars)
		Image compare := NewImage("Compare Fit", 2,c_fit_it, 3)
		compare[icol, 0] = auto_vals        // original data
		compare[icol, 1] = plot         // fit function
		ImageDocument linePlotDoc = CreateImageDocument("Test Fitting")
		ImageDisplay linePlotDsp = linePlotDoc.ImageDocumentAddImageDisplay(compare, 3)
		linePlotDoc.ImageDocumentShow()
		
		emsetfocus(start_focus + (zdf*step_size))
		Result("\n Zero defocus at :" + zdf + " " + (start_focus+(zdf*step_size)) )
	}
	
	if(coarsetype == "Slope")
	{
	while(diff_auto > auto_stop)
		{
		now_focus = emgetfocus()
		
		image img:= IntegerImage( "Img " ,bitdepth, signed, scan_size, scan_size )
		number imageid=img.imagegetid()
		
		number paramID = DSCreateParameters(scan_size, scan_size, rotation, exposure, linesynch )
		DSSetParametersSignal( paramID, signalIndex, bitdepth, selected, imageID )
		
		DSStartAcquisition(paramID, 0, synchronous);
		while(DSIsAcquisitionActive( )) yield() // wait for acquistion to finish
		
		auto_val_current = autocorr_coeff_find(img)
		diff_auto = auto_val_current - auto_val_previous
		auto_val_previous = auto_val_current
		
		Result("\n Current autocorrelation coefficient: " + auto_val_current)
		Result(", difference in autocorr coeff: " + diff_auto)
		Result("\n Focus step size: " + step_size)
		Result(", current focus: " + now_focus)
		Result("\n")
		
		DeleteImage(img)
		
		if (diff_auto < 0)
			{
			step_size = -step_size / (2)
			diff_auto = -diff_auto 
			}
			
		emsetfocus(now_focus+step_size)
		
		it_num = it_num + 1
		if(it_num > max_it)
			{
			break
			}
		
		}
		emsetfocus(now_focus-step_size)
	}
	}
	
number autofocus_fine(string finetype,number focusrange, number step_size, number scan_size, number exposure, number top, number left, number bottom, number right, number xsize, number ysize, image img)
	{
	/* 
		Function to determine the correct focus of a microscope 
		using a measure of quality of FFTs.
	*/
	
	number signalIndex = 0 // The signal on channel zero - set in the microscope interface
	number rotation = 0 //currentframes * 90   // degree of rotation
	number continuous = 0 // have continuous or single frame acquisitions
	number lineSynch = 0 // do not use line synching
	number synchronous = 0 // 0 = return immediately, 1 = return when finished
	number selected=1 // acquire signal
	number bitdepth = 4
	number signed = 0
	number xmax
	number ymax

	number init_focus = emgetfocus()
	number max_it = focusrange / step_size
	image coefficients := realimage("Focus Coefficients",4,max_it,1)
	
	emsetfocus(init_focus-max_it*step_size/2)
	number start_focus = emgetfocus()
	
	number i
	for(i = 0; i<max_it; i++)
		{
		image sub_img:= IntegerImage( "Fine Focus Image",4, 0, xsize, ysize )    
		
		//number paramID = DSCreateParameters(scan_size, scan_size, rotation, exposure, linesynch )
		//DSSetParametersSignal( paramID, signalIndex, bitdepth, selected, imageID )
		
		//DSStartAcquisition(paramID, 0, synchronous);
		
		DSScanSubRegion(img,sub_img,top,left,bottom,right,exposure )
		
		while(DSIsAcquisitionActive( )) yield() // wait for acquistion to finish
		
		//result("\n FFT " + i + ": ")
		if(finetype=="FFT")
			{coefficients[i,0] = fft_coeff_find(sub_img)
			}
			
		if(finetype=="AC")
			{coefficients[i,0] = autocorr_coeff_find(sub_img)
			}
		
		DeleteImage(sub_img)
		emsetfocus(emgetfocus()+step_size)
		}
		
	coefficients.showimage()	
	min(coefficients,xmax,ymax)
	emsetfocus(start_focus+(xmax*step_size))
	Result("\n Found fine focus value:" + (start_focus+(xmax*step_size)))
	}
	
number autostig(number stigrange, number step_size, number scan_size, number exposure, number top, number left, number bottom, number right, number xsize, number ysize, image img)
	{
	/* 
		WARNING: Untested.
		Function to automatically correct image stigmation using FFT 
		coefficients.
	*/
	
	number signalIndex = 0 // The signal on channel zero - set in the microscope interface
	number rotation = 0 //currentframes * 90   // degree of rotation
	number continuous = 0 // have continuous or single frame acquisitions
	number lineSynch = 0 // do not use line synching
	number synchronous = 0 // 0 = return immediately, 1 = return when finished
	number selected=1 // acquire signal
	number bitdepth = 4
	number signed = 0
	number xmax
	number ymax

	number init_xstig, init_ystig,start_xstig, start_ystig
	emgetcondensorstigmation(init_xstig,init_ystig)
	number max_it = stigrange / step_size
	image coefficients := realimage("Focus Coefficients",4,max_it,1)
	
	emsetcondensorstigmation(init_xstig-max_it*step_size/2,init_ystig)
	emgetcondensorstigmation(start_xstig,init_ystig)
	
	number i
	for(i = 0; i<max_it; i++)
		{
		image sub_img:= IntegerImage( "Stig X Image",4, 0, xsize, ysize )    
		
		//number paramID = DSCreateParameters(scan_size, scan_size, rotation, exposure, linesynch )
		//DSSetParametersSignal( paramID, signalIndex, bitdepth, selected, imageID )
		
		//DSStartAcquisition(paramID, 0, synchronous);
		
		DSScanSubRegion(img,sub_img,top,left,bottom,right,exposure )
		
		while(DSIsAcquisitionActive( )) yield() // wait for acquistion to finish
		
		//result("\n FFT " + i + ": ")
		coefficients[i,0] = fft_coeff_find(sub_img)
		
		DeleteImage(sub_img)
		number xstignow
		emgetcondensorstigmation(xstignow,init_ystig)
		emsetcondensorstigmation(xstignow+step_size,init_ystig)
		}
		
	max(coefficients,xmax,ymax)
	number cal_xstig = start_xstig+(xmax*step_size)
	emsetcondensorstigmation(cal_xstig,init_ystig)
	
	emsetcondensorstigmation(cal_xstig,init_ystig-max_it*step_size/2)
	emgetcondensorstigmation(cal_xstig,start_ystig)
	
	number j
	for(j = 0; j<max_it; j++)
		{
		image sub_img:= IntegerImage( "Stig Y Image",4, 0, xsize, ysize )    
		
		//number paramID = DSCreateParameters(scan_size, scan_size, rotation, exposure, linesynch )
		//DSSetParametersSignal( paramID, signalIndex, bitdepth, selected, imageID )
		
		//DSStartAcquisition(paramID, 0, synchronous);
		
		DSScanSubRegion(img,sub_img,top,left,bottom,right,exposure )
		
		while(DSIsAcquisitionActive( )) yield() // wait for acquistion to finish
		
		//result("\n FFT " + i + ": ")
		coefficients[i,0] = fft_coeff_find(sub_img)
		
		DeleteImage(sub_img)
		number ystignow
		emgetcondensorstigmation(cal_xstig,ystignow)
		emsetcondensorstigmation(cal_xstig,ystignow+step_size)
		}
		
	max(coefficients,xmax,ymax)
	number cal_ystig = start_ystig+(xmax*step_size)
	emsetcondensorstigmation(cal_xstig,cal_ystig)
	
	Result("\n Found stigmator values: x=" + cal_xstig + " y=" + cal_ystig)
	}
	
/*number autofocus(string coarseflag, string fineflag, string coarsetype, string finetype, number coarse_step, number auto_stop, number coarse_scansize, number coarse_exposure, number coarse_it, number focusrange, number fine_stepsize, number fine_scansize, number fine_exposure)
	{
	if(coarseflag=="True")
	autofocus_correlation(coarsetype, coarse_step, auto_stop, coarse_scansize, coarse_exposure, coarse_it)
	
	if(fineflag=="True")
	autofocus_fine(finetype, focusrange, fine_stepsize, fine_scansize, fine_exposure)
	}
*/

number autostem(taggroup stagelims, taggroup imsettings, taggroup focussettings)
	{
	/* 
		Function to perform automated STEM imaging of nanoparticle samples.
		Input: Tag groups containing the necessary parameters.
	*/
	
	number im_count = 0
	
	//Get the stage positions from the tag group
	number min_x, max_x, min_y, max_y
	stagelims.taggroupgettagasnumber("min_x", min_x)
	stagelims.taggroupgettagasnumber("max_x", max_x)
	stagelims.taggroupgettagasnumber("min_y", min_y)
	stagelims.taggroupgettagasnumber("max_y", max_y)
	
	number diststage_x = (max_x-min_x)
	number diststage_y = (max_y-min_y)
	
	string dev_loc = dsgetdevicelocation()
	number cal_fov
	
	//Get the image settings from the tag group
	number lowmag, highmag, scansize_low, scansize_high
	imsettings.taggroupgettagasnumber("lowmag", lowmag)
	imsettings.taggroupgettagasnumber("highmag", highmag)
	imsettings.taggroupgettagasnumber("scansize_low", scansize_low)
	imsettings.taggroupgettagasnumber("scansize_high", scansize_high)
	
	//Get the autofocus settings from the tag group
	number c_check, f_check, c_freq, c_step, c_scansize, c_exposure, c_it, f_freq, frange, f_stepsize, f_scansize, f_exposure
	focussettings.taggroupgettagasnumber("c_check", c_check)
	focussettings.taggroupgettagasnumber("f_check", f_check)
	focussettings.taggroupgettagasnumber("c_freq", c_freq)
	focussettings.taggroupgettagasnumber("c_step", c_step)
	focussettings.taggroupgettagasnumber("c_scansize", c_scansize)
	focussettings.taggroupgettagasnumber("c_exposure", c_exposure)
	focussettings.taggroupgettagasnumber("c_it", c_it)
	focussettings.taggroupgettagasnumber("f_freq", f_freq)
	focussettings.taggroupgettagasnumber("frange", frange)
	focussettings.taggroupgettagasnumber("f_stepsize", f_stepsize)
	focussettings.taggroupgettagasnumber("f_scansize", f_scansize)
	focussettings.taggroupgettagasnumber("f_exposure", f_exposure)
	
	//Get lowmag field of view
	EMSetMagnification(LowMag)
	TagGroup state_info = emgetcalibrationstatetags()
	EMGetCalibratedFieldOfView(dev_loc,state_info,cal_fov,0)
	Result("\n FOV = "+cal_fov)
	//Set folder path
	string folder
	imsettings.taggroupgettagasstring("path", folder)
	//Result("\n "+folder)
	
	number signalIndex = 0 // The signal on channel zero - set in the microscope interface
	number rotation =  0 //currentframes * 90   // degree of rotation
	number continuous=0 // have continuous or single frame acquisitions
	number lineSynch = 0 // do not use line synching
	number synchronous = 0 // 0 = return immediately, 1 = return when finished
	number selected=1 // acquire signal
	
	//Set exposures
	number exposure_low, exposure_high
	imsettings.taggroupgettagasnumber("exposure_low",exposure_low)
	imsettings.taggroupgettagasnumber("exposure_high",exposure_high)
	
	number t, s
	
	for(t=0;t<(diststage_x/cal_fov);t++)
		{
		Result("\n Starting scan...")
		emsetstagex((min_x+t*cal_fov))
		EMWaitUntilReady()
		for(s=0;s<(diststage_y/cal_fov);s++)
			{
			emsetstagey((min_y+s*cal_fov))
			EMWaitUntilReady()
			
			if(s%c_freq==0)
			{
				string c_flag
				//if(c_check==1) c_flag = "True"
				
				//string f_flag = "False"
				
				string c_type = "Fit"
				//string f_type = "AC"

				number autostop = 0.001
				Result("Check")
				if(c_check==1) autofocus_correlation(c_type, c_step, autostop, c_scansize, c_exposure, c_it)
				
				//autofocus(c_flag,f_flag,c_type,f_type,c_step,autostop,c_scansize,c_exposure,c_it,frange,f_stepsize,f_scansize,f_exposure)
				
			}
			
			//Result("\n"+t+" "+s)
			number xsize_low = scansize_low
			number ysize_low = scansize_low
			//Result("\n " + ysize_low)
			image img
			img:= IntegerImage( "LowMag Img", 4, 0, xsize_low, ysize_low )   

			number imageid=img.imagegetid()
			number paramID = DSCreateParameters( xsize_low, ysize_low, rotation, exposure_low, linesynch )
			DSSetParametersSignal( paramID, signalIndex, 4, selected, imageID )
			//number waitforstopsignal=0.05 // dwell in seconds
			//StartSignal.WaitOnSignal(waitforstopsignal, StopSignal)
			DSStartAcquisition(paramID, 0, synchronous);
			
			while(DSIsAcquisitionActive( )) yield() // wait for acquistion to finish

			//Check if any particles in each area
			number sx, sy
			getsize( img, sx, sy )

			image mask := BinaryImage( "Binary Mask", sx, sy )      // The real input for analysis
			Image foundParticlesImage := BinaryImage( "found", sx, sy ) // Accepted particles mask (result)

			mask = img>0.6*max(img) ? 1 : 0           // This is a simple Treshold to test...

			string mFields = "CenterX,CenterY"      // Specify analysis results
			Number minParticleSize = 10                   // Limit accepted
			Number doLabel = 0                              // Add labeled mask to inputImage True/False

			imageDocument ResultsDoc = FindParticles(img,mask,mFields, minParticleSize, doLabel, foundParticlesImage)
			ImageDocumentShow( ResultsDoc )
			image x_img :=ImageDocumentGetImage(ResultsDoc,0)
			number num_p, num_vals
			getsize(x_img, num_vals, num_p)
			//Result(num_p)

			number k
			number l

			//Subscans
			number xsize_high = scansize_high
			number ysize_high = scansize_high
			number overlap = 0.1
			number M_ratio = LowMag/HighMag
			//Result("\n"+(1*ysize_low*M_ratio+M_ratio*overlap))

			number i
			number j
			number flag = 0

			for (i=0; i<(HighMag/LowMag); i++)
				{
				for (j=0; j<(HighMag/LowMag); j++)
					{
					flag = 0
					for (k=0; k<num_p; k++)
						{
						if((min(x_img[0,k])!=256.0) && (min(x_img[1,k]) != 256.0))
						if((min(x_img[0,k])>(i*xsize_low*M_ratio)) && (min(x_img[0,k]) < (i+1)*(xsize_low*M_ratio))&& (min(x_img[1,k]) > (j*ysize_low*M_ratio) && (min(x_img[1,k]) < (j+1)*(xsize_low*M_ratio))))
							{
							flag = 1
							//Result("\n Particles found")
						//if((min(x_img[0,k])>(i*xsize_low*M_ratio)) && (min(x_img[0,k]) < (i+1)*(xsize_low*M_ratio))&& (min(x_img[1,k]) > (j*ysize_low*M_ratio) && (min(x_img[1,k]) < (i+1)*(xsize_low*M_ratio))))
							
							//component imgdisp=imagegetimagedisplay(img, 0)
							//component ovalannot=newovalannotation(min(x_img[1,k])- 5, min(x_img[0,k])-5, min(x_img[1,k])+ 5, min(x_img[0,k])+5)
							//ovalannot.componentsetdrawingmode(1) // sets the background to off ie lines are not outlined
							//ovalannot.componentsetforegroundcolor(1,0,1) // sets the foreground colour to magenta
							//imgdisp.componentaddchildatend(ovalannot) // add the annotation to the image display
							//Result("\n "+(i*xsize_low*M_ratio)+" "+(i+1)*(xsize_low*M_ratio))
							
						}
							
					}
					//Result("\n "+flag)
					if(flag == 1)
						{
						
						//Result(j*ysize_low*M_ratio-M_ratio*overlap*ysize_low)
						number top = (j*ysize_low*M_ratio-M_ratio*overlap*ysize_low)
						number left = (i*xsize_low*M_ratio-M_ratio*overlap*xsize_low)
						number bottom = ((j+1)*ysize_low*M_ratio+M_ratio*overlap*ysize_low)
						number right = ((i+1)*xsize_low*M_ratio+M_ratio*overlap*xsize_low)
						
						if(j%f_freq==0)
							{						
							string f_type = "AC"
							
							if(f_check==1) autofocus_fine(f_type, frange, f_stepsize, f_scansize, f_exposure, top, left, bottom, right, xsize_high, ysize_high, img)

							}

						image sub_img:= IntegerImage( "sub_img ",4, 0, xsize_high, ysize_high )
						number sub_imageid=sub_img.imagegetid()
						string filename = folder+"goldtest"+im_count
						im_count = im_count + 1
						
						Showimage( sub_img )
						DSScanSubRegion(img,sub_img,top,left,bottom,right,exposure_high )
						while(DSIsAcquisitionActive( )) yield()
						//sub_img.saveimage(filename)
						saveasgatan(sub_img,filename)
						
						//Result("\n X range:"+(i*xsize_low*M_ratio)+" "+(i+1)*(xsize_low*M_ratio))
						//Result("\n Y range:"+(j*ysize_low*M_ratio)+" "+(j+1)*(ysize_low*M_ratio))
					}
				}
			}
			//deleteimage(img)
		}
	}
}

class FocusSettingsDialogClass : uiframe
	{
	// Creates the subdialog
	
	TagGroup CreateFocusSettingsDialog(object self) 
		{
			TagGroup SubDialog=DLGCreateDialog("Sub-dialog")
			
			TagGroup cbox_items
			Taggroup cBox=DLGCreateBox("Coarse Focus", cbox_items).dlgexternalpadding(3,3).dlginternalpadding(10,10)
			
			taggroup c_freq= DLGCreateBox("Frequency")
			c_freqfield = DLGCreateRealField(3, 12, 1);
			c_freq.dlgaddelement(c_freqfield)
			cbox_items.dlgaddelement(c_freq)

			taggroup c_step= DLGCreateBox("Focus Step (nm)")
			c_stepfield = DLGCreateRealField(50, 12, 1);
			c_step.dlgaddelement(c_stepfield)
			cbox_items.dlgaddelement(c_step)
			
			taggroup c_scansize= DLGCreateBox("Scan size (px)")
			c_scansizefield = DLGCreateRealField(512, 12, 1);
			c_scansize.dlgaddelement(c_scansizefield)
			cbox_items.dlgaddelement(c_scansize)
			
			taggroup c_exposure= DLGCreateBox("Exposure (s)")
			c_exposurefield = DLGCreateRealField(2, 12, 1);
			c_exposure.dlgaddelement(c_exposurefield)
			cbox_items.dlgaddelement(c_exposure)
			
			taggroup c_it= DLGCreateBox("Iterations")
			c_itfield = DLGCreateRealField(20, 12, 1);
			c_it.dlgaddelement(c_itfield)
			cbox_items.dlgaddelement(c_it)
			
			SubDialog.DLGAddElement(cbox)
			
			TagGroup fbox_items
			TagGroup fBox = DLGCreateBox("Fine Focus", fbox_items).dlgexternalpadding(3,3).dlginternalpadding(10,10)
		
			taggroup f_freq= DLGCreateBox("Frequency")
			f_freqfield = DLGCreateRealField(3, 12, 1);
			f_freq.dlgaddelement(f_freqfield)
			fbox_items.dlgaddelement(f_freq)
			
			taggroup frange= DLGCreateBox("Focus Range (nm)")
			frangefield = DLGCreateRealField(20, 12, 1);
			frange.dlgaddelement(frangefield)
			fbox_items.dlgaddelement(frange)
			
			taggroup f_stepsize= DLGCreateBox("Focus Step (nm)")
			f_stepsizefield = DLGCreateRealField(2, 12, 1);
			f_stepsize.dlgaddelement(f_stepsizefield)
			fbox_items.dlgaddelement(f_stepsize)
			
			taggroup f_scansize= DLGCreateBox("Scan size (px)")
			f_scansizefield = DLGCreateRealField(512, 12, 1);
			f_scansize.dlgaddelement(f_scansizefield)
			fbox_items.dlgaddelement(f_scansize)
			
			taggroup f_exposure= DLGCreateBox("Expsoure (s)")
			f_exposurefield = DLGCreateRealField(10, 12, 1);
			f_exposure.dlgaddelement(f_exposurefield)
			fbox_items.dlgaddelement(f_exposure)
			
			SubDialog.DLGAddElement(fbox)
			
			return SubDialog
		}

	FocusSettingsDialogClass(object self)
		{ 
			self.init(self.CreateFocusSettingsDialog())
		}
		
	~FocusSettingsDialogClass(object self) 
		{
		}
}

class ImageSettingsDialogClass : uiframe
	{
	
	TagGroup CreateImageSettingsDialog(object self) 
		{
			TagGroup SubDialog=DLGCreateDialog("Sub-dialog2")
			
			TagGroup lbox_items
			Taggroup lBox=DLGCreateBox("Low Mag Image", lbox_items).dlgexternalpadding(3,3).dlginternalpadding(10,10)
			
			taggroup lowmag= DLGCreateBox("Magnification")
			lowmagfield = DLGCreateRealField(1e6, 12, 1);
			lowmag.dlgaddelement(lowmagfield)
			lbox_items.dlgaddelement(lowmag)
			
			taggroup limsize= DLGCreateBox("Image size (px)")
			limsizefield = DLGCreateRealField(512, 12, 1);
			limsize.dlgaddelement(limsizefield)
			lbox_items.dlgaddelement(limsize)
			
			taggroup lexposure= DLGCreateBox("Exposure (us)")
			lexposurefield = DLGCreateRealField(2, 12, 1);
			lexposure.dlgaddelement(lexposurefield)
			lbox_items.dlgaddelement(lexposure)
			
			SubDialog.DLGAddElement(lbox)
			
			TagGroup hbox_items
			Taggroup hBox=DLGCreateBox("High Mag Image", hbox_items).dlgexternalpadding(3,3).dlginternalpadding(10,10)
			
			taggroup highmag= DLGCreateBox("Magnification")
			highmagfield = DLGCreateRealField(4e6, 12, 1);
			highmag.dlgaddelement(highmagfield)
			hbox_items.dlgaddelement(highmag)
			
			taggroup himsize= DLGCreateBox("Image size (px)")
			himsizefield = DLGCreateRealField(1024, 12, 1);
			himsize.dlgaddelement(himsizefield)
			hbox_items.dlgaddelement(himsize)
			
			taggroup hexposure= DLGCreateBox("Exposure (us)")
			hexposurefield = DLGCreateRealField(10, 12, 1);
			hexposure.dlgaddelement(hexposurefield)
			hbox_items.dlgaddelement(hexposure)
			
			SubDialog.DLGAddElement(hbox)
			
			return SubDialog
		}
		
	ImageSettingsDialogClass(object self)
		{ 
			self.init(self.CreateImageSettingsDialog())
		}
		
	~ImageSettingsDialogClass(object self) 
		{
		}
	}

class CreateAutoImDialog : uiframe
	{
	//Class to create the Auto Imaging Dialog
	
	TagGroup CreateStageTagGroup(object self)
	{
		number index
		taggroup postags = NewTagGroup()
		
		index = postags.TagGroupCreateNewLabeledTag( "min_x" ) 
		postags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(minxfield) )
		index = postags.TagGroupCreateNewLabeledTag( "max_x" ) 
		postags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(maxxfield) )
		index = postags.TagGroupCreateNewLabeledTag( "min_y" ) 
		postags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(minyfield) )
		index = postags.TagGroupCreateNewLabeledTag( "max_y" ) 
		postags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(maxyfield) )
		
		return postags
	}
	
	TagGroup CreateImTagGroup(object self)
		{
		number index
		imtags = NewTagGroup()
		
		index = imtags.TagGroupCreateNewLabeledTag( "lowmag" ) 
		imtags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(lowmagfield) )
		index = imtags.TagGroupCreateNewLabeledTag( "highmag" ) 
		imtags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(highmagfield) )
		index = imtags.TagGroupCreateNewLabeledTag( "scansize_low" ) 
		imtags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(limsizefield) )
		index = imtags.TagGroupCreateNewLabeledTag( "scansize_high" ) 
		imtags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(himsizefield) )
		index = imtags.TagGroupCreateNewLabeledTag( "exposure_low" ) 
		imtags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(lexposurefield) )
		index = imtags.TagGroupCreateNewLabeledTag( "exposure_high" ) 
		imtags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(hexposurefield) )
		index = imtags.TagGroupCreateNewLabeledTag( "path" ) 
		imtags.TagGroupSetIndexedTagAsString( index, dlggetstringvalue(pathfield) )
		
		return imtags
		}
		
	TagGroup CreateFocusTagGroup(object self)
		{
		number index
		focustags = NewTagGroup()
		
		number c_on, f_on
		index = focustags.TagGroupCreateNewLabeledTag( "c_check" )
		self.dlggetvalue("c_check",c_on)
		focustags.TagGroupSetIndexedTagAsNumber( index, c_on )
		index = focustags.TagGroupCreateNewLabeledTag( "f_check" ) 
		self.dlggetvalue("f_check",f_on)
		focustags.TagGroupSetIndexedTagAsNumber( index, f_on )
		
		index = focustags.TagGroupCreateNewLabeledTag( "c_freq" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(c_freqfield) )
		index = focustags.TagGroupCreateNewLabeledTag( "c_step" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(c_stepfield) )
		index = focustags.TagGroupCreateNewLabeledTag( "c_scansize" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(c_scansizefield) )
		index = focustags.TagGroupCreateNewLabeledTag( "c_exposure" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(c_exposurefield) )
		index = focustags.TagGroupCreateNewLabeledTag( "c_it" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(c_itfield) )
		index = focustags.TagGroupCreateNewLabeledTag( "f_freq" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(f_freqfield) )
		index = focustags.TagGroupCreateNewLabeledTag( "frange" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(frangefield) )
		index = focustags.TagGroupCreateNewLabeledTag( "f_stepsize" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(f_stepsizefield) )
		index = focustags.TagGroupCreateNewLabeledTag( "f_scansize" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(f_scansizefield) )
		index = focustags.TagGroupCreateNewLabeledTag( "f_exposure" ) 
		focustags.TagGroupSetIndexedTagAsNumber( index, dlggetvalue(f_exposurefield) )
		
		return focustags
		}

	void PressMinPosButton(object self)
		{
		number xstagepos = emgetstagex()
		number ystagepos = emgetstagey()
		dlgvalue(minxfield,xstagepos)
		dlgvalue(minyfield,ystagepos)
		//Result("X Stage Pos: "+xstagepos)
	}
	
	void PressMaxPosButton(object self)
		{
		number xstagepos = emgetstagex()
		number ystagepos = emgetstagey()
		dlgvalue(maxxfield,xstagepos)
		dlgvalue(maxyfield,ystagepos)
		//Result("X Stage Pos: "+xstagepos)
	}
	
	void PressPathButton(object self)
		{
		String EMFolder

		GetDirectoryDialog( "Select folder" , "X:/SessionData/data" , EMFolder )
		dlgvalue(pathfield,EMFolder)
		
		}
	
	void PressAcquireButton(object self)
		{
		
		taggroup postags = self.CreateStageTagGroup()
		taggroup imsettings = self.CreateImTagGroup()
		taggroup focus_settings = self.CreateFocusTagGroup()
		
		number min_x
		postags.taggroupgettagasnumber("min_x", min_x)
		
		//string test
		//imsettings.taggroupgettagasstring("path", test)
		
		Result("\nAcquiring automated STEM capture within the following spatial limits:")
		Result("\n Min x pos: " + min_x + " Min y pos: " + dlggetvalue(minyfield))
		Result("\n Max x pos: " + dlggetvalue(maxxfield) + " Max y pos: " + dlggetvalue(maxyfield))
		Result("\n Folder: " + dlggetstringvalue(pathfield))

		//Result("\n Test: " + testnum)
		
		autostem(postags,imsettings,focus_settings)
	}
	
	void AutoFocusSettings(object self)
		{
		FocusSettingsDialogObject.display("Auto Focus Settings").WindowSetFramePosition(500, 300 )
		}
		
	void ImageSettings(object self)
		{
		ImageSettingsDialogObject.display("Image Settings").WindowSetFramePosition(500, 300 )
		}
	
	TagGroup CreateDialog(object self)
		{
		TagGroup dialog = DLGCreateDialog("Automated Acquisition")
		
		TagGroup posfields = DLGCreateGroup()
		DLGLayout(posfields, DLGCreateTableLayout(3,2,0))
		dialog.dlgaddelement(posfields)
		
		TagGroup MinPosButton = DLGCreatePushButton("Get Min Stage Position","PressMinPosButton")
		posfields.dlgaddelement(MinPosButton)
		
		minxfield = DLGCreateRealField(emgetstagex(),10,4)
		posfields.dlgaddelement(minxfield)
		minyfield = DLGCreateRealField(emgetstagey(),10,4)
		posfields.dlgaddelement(minyfield)
		
		TagGroup MaxPosButton = DLGCreatePushButton("Get Max Stage Position","PressMaxPosButton")
		posfields.dlgaddelement(MaxPosButton)
		
		maxxfield = DLGCreateRealField(emgetstagex(),10,4)
		posfields.dlgaddelement(maxxfield)
		maxyfield = DLGCreateRealField(emgetstagey(),10,4)
		posfields.dlgaddelement(maxyfield)
		
		TagGroup settingsfields = DLGCreateGroup()
		DLGLayout(settingsfields, DLGCreateTableLayout(2,1,0))
		dialog.dlgaddelement(settingsfields)
		
		tagGroup DisplayImageSettingsButton=DLGCreatePushButton("Image Settings", "ImageSettings")
		settingsfields.dlgaddelement(displayimagesettingsbutton)
		
		tagGroup DisplayFocusSettingsButton=DLGCreatePushButton("Auto Focus Settings", "AutoFocusSettings")
		settingsfields.dlgaddelement(displayfocussettingsbutton)
		
		taggroup c_check= DLGCreateCheckBox("Apply Coarse Focus",1)
		c_check.dlgidentifier("c_check")
		dialog.dlgaddelement(c_check)
		
		taggroup f_check= DLGCreateCheckBox("Apply Fine Focus",1)
		f_check.dlgidentifier("f_check")
		dialog.dlgaddelement(f_check)
		
		TagGroup pathbox_items
		Taggroup pathBox=DLGCreateBox("Set Path for Images", pathbox_items).dlgexternalpadding(3,3).dlginternalpadding(10,10)
		
		TagGroup pathfields = DLGCreateGroup()
		DLGLayout(pathfields, DLGCreateTableLayout(2,1,0))
		pathbox_items.dlgaddelement(pathfields)
		
		pathfield = DLGCreateStringField("X:/SessionData/data",50)
		pathfields.dlgaddelement(pathfield)
		tagGroup PathButton=DLGCreatePushButton("Path", "PressPathButton")
		pathfields.dlgaddelement(pathbutton)
		
		dialog.DLGAddElement(pathbox)
		
		TagGroup AcquireButton = DLGCreatePushButton("Acquire","PressAcquireButton")
		dialog.dlgaddelement(acquirebutton)
		
		return dialog

	}
	
	CreateAutoImDialog(object self)
		{
		self.init(self.CreateDialog())
		self.display("Automated Acquisition");
	}
	
	~CreateAutoImDialog(object self) 
		{
			//result("\nMain Dialog destructed.")
		}
}

void main()
	{

	Alloc(CreateAutoImDialog)
	
	FocusSettingsDialogObject =Alloc(FocusSettingsDialogClass)
	
	ImageSettingsDialogObject =Alloc(ImageSettingsDialogClass)

}

main()
