; docformat = 'rst'
;
; NAME:
;       MrVar_Resample
;
;*****************************************************************************************
;   Copyright (c) 2018, Matthew Argall                                                   ;
;   All rights reserved.                                                                 ;
;                                                                                        ;
;   Redistribution and use in source and binary forms, with or without modification,     ;
;   are permitted provided that the following conditions are met:                        ;
;                                                                                        ;
;       * Redistributions of source code must retain the above copyright notice,         ;
;         this list of conditions and the following disclaimer.                          ;
;       * Redistributions in binary form must reproduce the above copyright notice,      ;
;         this list of conditions and the following disclaimer in the documentation      ;
;         and/or other materials provided with the distribution.                         ;
;       * Neither the name of the <ORGANIZATION> nor the names of its contributors may   ;
;         be used to endorse or promote products derived from this software without      ;
;         specific prior written permission.                                             ;
;                                                                                        ;
;   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY  ;
;   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES ;
;   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT  ;
;   SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,       ;
;   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED ;
;   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR   ;
;   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     ;
;   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN   ;
;   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH  ;
;   DAMAGE.                                                                              ;
;*****************************************************************************************
;
; PURPOSE:
;+
;   Time average one variable onto the time stamps of another.
;
; :Author:
;   Matthew Argall::
;       University of New Hampshire
;       Morse Hall, Room 348
;       8 College Rd.
;       Durham, NH, 03824
;       matthew.argall@unh.edu
;
; :History:
;   Modification History::
;       2018/11/04  -   Written by Matthew Argall
;       2019/01/15  -   Added the T_EDGE keyword. - MRA
;-
;*****************************************************************************************
;+
;   Test MrVar_Time_Avg using artificial data.
;
; :Returns:
;       WIN:        out, required, type=objref
;                   A MrWindow object reference containing a visualization of the results.
;-
FUNCTION MrVar_Resample_Test
	Compile_Opt idl2
	On_Error, 2
	
	;Reference time
	t_ref  = '2017-08-03T15:18:33'
	length = 10           ; s
	f1     = 1.0          ; Hz
	f2     = 3.0          ; Hz
	
;-------------------------------------------
; Create Two Signals ///////////////////////
;-------------------------------------------
	
	;SIGNAL 1
	fs1     = 5.0           ; S/s
	dt1     = 1.0 / fs1     ; s
	N1      = length * fs1  ; Samples
	t1      = dt1 * FIndGen(N1)
	s1      = Sin( 2 * !pi * f1 * t1 ) + RandomU(5, N1) - 0.5

	oT1 = MrTimeVar( t1, 'SSM', $
	                 T_REF = t_ref )
	oT1['UNITS']       = 's'
	oT1['DELTA_MINUS'] = 0.5 * dt1
	oT1['DELTA_PLUS']  = 0.5 * dt1
	
	oS1 = MrTimeSeries( oT1, s1 )
	
	;SIGNAL 2
	fs2     = 20.0          ; S/s
	dt2     = 1.0 / fs2     ; s
	N2      = length * fs2  ; Samples
	t2      = dt2 * FIndGen(N2)
	s2      = Sin( 2 * !pi * f2 * t2 ) + RandomU(5, N2) - 0.5

	oT2 = MrTimeVar( t2, 'SSM', $
	                 T_REF = t_ref )
	oT2['UNITS']       = 's'
	oT2['DELTA_MINUS'] = 0.5 * dt2
	oT2['DELTA_PLUS']  = 0.5 * dt2
	
	oS2 = MrTimeSeries( oT2, s2 )
	
;-------------------------------------------
; Plot Results /////////////////////////////
;-------------------------------------------
	;Time-Average Signal 2
	oS2_Avg = MrVar_Resample(oS2, oS1)
	
;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	;SIGNAL 1
	oS1['AXIS_RANGE'] = [-1.5, 1.5]
	oS1['PLOT_TITLE'] = 'Targeted Signal'
	oS1['SYMBOL']     = 6
	oS1['TITLE']      = 'S1'
	
	;SIGNAL 2
	oS2['AXIS_RANGE'] = [-1.5, 1.5]
	oS2['PLOT_TITLE'] = 'Source Signal'
	oS2['LABEL']      = 'Orignal'
	oS2['SYMBOL']     = 6
	oS2['TITLE']      = 'S2'
	
	;RESULTS
	oS2_Avg['COLOR']     = 'Red'
	oS2_Avg['LABEL']     = 'Avg'
;	oS2_Avg['LINESTYLE'] = 'None'
	oS2_Avg['SYMBOL']    = 7
	
;-------------------------------------------
; Plot Results /////////////////////////////
;-------------------------------------------
	;Create a window in which to display the graphics
	win = MrWindow( OXMARGIN = [10, 8], $
	                REFRESH=0 )
	
	;Plot the results
	p1 = MrVar_Plot( oS1, /CURRENT )
	p2 = MrVar_Plot( oS2, /CURRENT )
	p3 = MrVar_Plot( oS2_Avg, OVERPLOT=p2 )
	
	;Finish
	win -> SetGlobal, XRANGE=[0,10], XTICKINTERVAL=1
	win -> Refresh
	RETURN, win
END


;+
;   Average data from one dataset onto the time tags of another dataset. It is assumed
;   that both the sample times and the target times are center times of the measurement.
;   A nearest-neighbor scheme is used to match sample to target.
;
; :Params:
;       OV:         in, required, type=objref
;                   The MrTimeSeries instance to be averaged.
;       OT:         in, required, type=objref
;                   The MrTimeVar instance containing the target times.
;       Q:          in, required, type=float
;                   The filter factor
;       WIN:        in, required, type=string
;                   Type of tapering window to use when averaging. Choices
;                       include: {'BOXCAR' | 'RECTANGULAR' | 'TRIANGLE' | 'BARTLETT'}.
;
; :Returns:
;       OV_AVG:     out, required, type=objref
;                   Averaged dataset.
;-
FUNCTION MrVar_Resample_BoxCar, oV, oT, q, win
	Compile_Opt idl2
	On_Error, 2
	
	dt_in     = oV['TIMEVAR'] -> GetSI(RATE=fs_in)
	dt_target = oT -> GetSI(RATE=fs_target)
	n_in      = oV -> GetNPts()
	n_target  = oT.N_Elements
	
;-------------------------------------------
; Box Car Filter ///////////////////////////
;-------------------------------------------
	n_filter = Fix(2*q, TYPE=3) + 1
	
	;Allocate memory
	;   - The number of time samples in the out put will be reduced from the input
	nDims    = Size(oV, /N_DIMENSIONS)
	dims     = Size(oV, /DIMENSIONS)
	dims[0]  = n_target
	v_out    = Make_Array( dims, TYPE=Size(self, /TYPE), VALUE=!Values.F_NaN )
	sdev_out = Make_Array( dims, TYPE=Size(self, /TYPE), VALUE=!Values.F_NaN )
	
	;Locate the targets
	iNear = oV['TIMEVAR'] -> Nearest_Neighbor(oT['DATA','TT2000'], 'TT2000')
	
	t_ssm_target = oT['DATA', 'SSM']
	t_ssm        = oV['TIME', 'SSM']
	
	;Resample and filter
	pv = oV['PTR']
	l  = 0
	n  = 0
	FOR i = 0, n_target - 1 DO BEGIN
		;Save previous interval
		lp = l
		np = n
		
		;BoXCar
		IF win EQ 'BOXCAR' || win EQ 'SQUARE' || win EQ 'RECTANGULAR' THEN BEGIN
			l = iNear[i] - Floor(n_filter/2.0)
			n = l + n_filter - 1
			IF l GT 0 && (t_ssm_target[i] - t_ssm[l-1]) LT (t_ssm[n] - t_ssm_target[i]) THEN BEGIN
				l -= 1
				n -= 1
			ENDIF
			
			;Normalized filter weights
			w = Replicate(1.0, n_filter)
		
		;Bartlett
		ENDIF ELSE IF win EQ 'BARTLETT' THEN BEGIN
			l = iNear[i]
			x = Abs(t_ssm_target[i] - t_ssm[l]) / dt_in
			
			;Number of points before and after center
			lm = Fix(2*q - x, TYPE=3)
			nm = Fix(2*q + x, TYPE=3)
			
			;Barlett Window
			n = l + lm + nm
			w = FltArr(lm + nm + 1)
			
			w[0:lm-1]     = Reverse(1 - (IndGen(lm)+1 + x) / (2*q))
			w[lm]         = 1 - Abs(x) / (2*q)
			w[lm+1:lm+nm] = 1 - (IndGen(nm)+1 - x) / (2*q)
		ENDIF ELSE BEGIN
			Message, 'Window name not recognized: "' + win + '".'
		ENDELSE
		
		;Normalize
		dims_out    = dims
		dims_out[0] = lm+nm+1
		w    /= Total(w)
		!Null = Where(w NE 0, nw)
		w     = Rebin(w, dims_out)
		
		;Skip
		IF l LT 0 || n GE n_in THEN CONTINUE
		
		;Amount of overlap
;		overlap = Float(np - l + 1) / (n - l + 1) * 100.0
;		IF overlap LT 50 THEN MrPrintF, 'LogWarn', 'Overlap should be > 50%.'
		
;		dlo = ( (t_ssm[l] + dt_in/2*Fix(2*q, TYPE=3)) - t_ssm_target[i] ) / dt_in ; (-1, 0]
;		dhi = ( (t_ssm[n] - dt_in/2*Fix(2*q, TYPE=3)) - t_ssm_target[i] ) / dt_in ; [0, 1)
;		MrPrintF, 'LogText', dlo, dhi, overlap, n_filter, l, n, n-l+1, FORMAT='(%"dt_lo=%f dt_hi=%f overlap=%f\% N=%i l=%i n=%i n-l+1=%i")'
		
		;Mean
		avg = Total(w*(*pv)[l:n,*,*,*,*,*,*,*], 1, NAN=NaN)
		
		;Standard Deviation
		IF nDims EQ 1 THEN BEGIN
			sdev = Sqrt( Total( w * ((*pv)[l:n,*,*,*,*,*,*,*] - avg)^2 ) $ 
			              / ( (nw-1)*Total(w) / nw ) )
		ENDIF ELSE BEGIN
			sdev = Sqrt( Total( w * ((*pv)[l:n,*,*,*,*,*,*,*] $
			                    - Rebin(Reform(avg, 1, dims[1:*]), dims_out))^2, 1 ) $
			             / ( (nw-1)*Total(w, 1) / nw ) )
		ENDELSE
		
		;Store the data
		ref    = dims_out
		ref[0] = 1
		v_out[i,0,0,0,0,0,0,0]    = Reform(Temporary(avg),  ref)
		sdev_out[i,0,0,0,0,0,0,0] = Reform(Temporary(sdev), ref)
	ENDFOR
	
	;Create the variable
	oV_avg = Obj_New( Obj_Class(oV), oT, v_out, $
	                  CACHE = cache, $
	                  NAME  = name )
	oV_sdev = Obj_New( Obj_Class(oV), oT, sdev_out, $
	                   NAME = 'Time_Avg_SDev(' + oV.name + ')' )
	
	;Attributes
	oV_avg['DELTA_PLUS_VAR'] = oV_sdev
	oV_avg['DELTA_MINUS_VAR'] = oV_sdev
	
	RETURN, oV_avg
END


;+
;   Time average one variable onto the time stamps of another.
;
; :Params:
;       VAR:            in, required, type=string/integer/objref
;                       Name, number, or MrTimeSeries variable to be averaged.
;       TIME:           in, required, type=string/integer/objref
;                       Name, number, MrTimeSeries or MrTimeVar objref of the variable
;                           containing the target time tags.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output will be added to the variable cache.
;       NAN:            in, optional, type=boolean, default=0
;                       If set, NaNs will be ignored during the averaging process.
;       NAME:           in, optional, type=string, default='Cyclotron_Frequency'
;                       Name to be given to the output variable.
;       NO_CLOBBER:     in, optional, type=boolean, default=0
;                       If set, cached variables with name `NAME` are not over-written.
;                           Instead, the output variable will have "_#" appended, were "#"
;                           represents a unique number.
;       TEST:           in, optional, type=boolean, default=0
;                       If set, run MrVar_Time_Avg with a set of synthetic data and return
;                           a MrGraphics window visualizing the result.
;
; :Returns:
;       OV_AVG:         out, required, type=objref
;                       A MrTimeSeries objref containing the averaged data.
;
FUNCTION MrVar_Resample, var, time, $
CACHE=cache, $
FILTER_FACTOR=filter_factor, $
NAN=nan, $
NAME=name, $
NO_CLOBBER=no_clobber, $
REDUCTION_FACTOR=reduction_factor, $
TEST=test, $
WINDOW=window
	Compile_Opt idl2
	On_Error, 2
	
	IF Keyword_Set(test) THEN BEGIN
		win = MrVar_Time_Avg_Test()
		RETURN, win
	ENDIF
	
	win_name = N_Elements(window) EQ 0 ? 'BARTLETT' : StrUpCase(window)
	IF N_Elements(t_type) EQ 0 THEN t_type = 'TT2000'
	
	;Grab the variables
	oV = MrVar_Get(var)
	oT = MrVar_Get(time)
	IF ~Obj_IsA(oV, 'MrTimeSeries') THEN Message, 'VAR must be a MrTimeSeries variable.'
	
	;Time
	IF ~Obj_IsA(oT, 'MrTimeVar') THEN BEGIN
		IF Obj_IsA(oT, 'MrTimeSeries') $
			THEN oT = oT['TIMEVAR'] $
			ELSE Message, 'TIME must be a MrTimeVar or MrTimeSeries variable.'
	ENDIF
	
	IF N_Elements(name) EQ 0 THEN name = 'Resample(' + oV.name + ')'
	
;-------------------------------------------
; Filter Coefficients //////////////////////
;-------------------------------------------
	;Sampling intervals and rates
	dt_in     = oV['TIMEVAR'] -> GetSI(RATE=fs_in)
	dt_target = oT -> GetSI(RATE=fs_target)
	
	;Nyquist Frequencies
	fN_in     = fs_in / 2.0
	fN_target = fs_target / 2.0
	
	;Sampling factor
	s = N_Elements(sample_factor) EQ 0 ? 1.0 : sample_factor
	
	;Reduction factor
	r = fN_in / fN_target
	
	;Filter factor
	q = N_Elements(filter_factor) EQ 0 ? r / s : filter_factor
	
	;Sampling rates are almost equal
	IF r GE 0.5 || R LE 2 THEN MrPrintF, 'LogText', 'r~1. Careful when resampling!'
	
	;Warn about over- or under-sampled data
	IF s GT 1 THEN BEGIN
		MrPrintF, 'LogText', 'Data is over-sampled.'
	ENDIF ELSE IF s LT 1 THEN BEGIN
		MrPrintF, 'LogText', 'Data is under-sampled.'
	ENDIF
	
	;Skip resampling or filtering
	IF s LT 1 THEN BEGIN
		MrPrintF, 'LogText', 'S < 1. Leaving as is.'
		MrPrintF, 'LogText', 'Resampling does not change effects of aliasing already present.'
	ENDIF ELSE IF q/s GT 100 && r/s GT 100 THEN BEGIN
		MrPrintF, 'LogText', 'q >> S and r >> S, so s = q/r.'
		MrPrintF, 'LogText', 'Resampling without filtering as filter extends to ' + $
		                     'frequencies much higher than fN of the input.'
	ENDIF
	
	;Interpolate if
	;   - ∆T >= ∆t   Reducing sampling interval
	;   - F_N < f_N  Increasing Nyquist frequency
	;   - r <= 1     Increasing the number of samples
	tf_interpolate = r LE 1
	IF tf_interpolate THEN MrPrintF, 'LogText', 'r <= 1. Interpolating.'
	 
	;Resultant sample factor
	s_target = q/r * s
	
	MrPrintF, 'StdOut', s, r, q, s_target, FORMAT='(%"S = %f, r = %f, q = %f, s = %f")'
	
;-------------------------------------------
; Box Car Average //////////////////////////
;-------------------------------------------
	IF tf_interpolate THEN BEGIN
		oV_avg = oV -> Interpol(oT)
	ENDIF ELSE BEGIN
		oV_avg = MrVar_Resample_BoxCar(oV, oT, q, win_name)
	ENDELSE
	
	oV_avg -> SetName, name
	IF Keyword_Set(cache) THEN oV_avg -> Cache
	
	RETURN, oV_avg
END