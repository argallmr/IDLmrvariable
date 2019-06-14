; docformat = 'rst'
;
; NAME:
;       MrVar_PoyntingSpectra
;
;*****************************************************************************************
;   Copyright (c) 2016, Matthew Argall                                                   ;
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
;   The purpose of this program is to calculate the Poynting Vector, in micro-Watts per
;   square meter (W/m^2 * 1e-6), in the spectral domain::
;       \vec{S}_{avg} = \frac{1}{4 \mu} \left(\delta \vec{E}^{*} \times \delta \vec{B} +
;                       \delta \vec{E} \times \delta \vec{B}^{*} \right)
;
;   References::
;       Loto'aniu, T, M., et. al., Propagation of electromagnetic ion cyclotron wave
;           energy in the magnetosphere, J. Geophys. Res., 110, 2005.
;           doi:10.1029/2004JA010816
;
; :Categories:
;   MrVariable
;
; :Params:
;       E:              in, required, type=string/integer/objref
;                       Name, number, or objref of a MrVectorTS variable containing the
;                           vector electric field (mV/m). If E and `B` do not have
;                           identical time stamps, `E` will be inteprolated to `B`.
;       B:              in, required, type=string/integer/objref
;                       Name, number, or objref of a MrVectorTS variable containing the
;                           vector magnetic field (nT).
;       NFFT:           in, optional, type=int, default=256
;                       The number of points to use per FFT
;       NSHIFT:         in, optional, type=int. default=`NFFT`/2
;                       The number of points to shift ahead after each FFT.
;
; :Keywords:
;       B0:             in, type=string/integer/objref
;                       Name, number, or objref of a MrVectorTS variable containing the
;                           background magnetic field (nT). If present, `B` and `E` will
;                           be rotating into a field-aligned coordinate system.
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output will be added to the variable cache.
;       NAME:           in, optional, type=string, default='Cyclotron_Frequency'
;                       Name to be given to the output variable.
;       NO_CLOBBER:     in, optional, type=boolean, default=0
;                       If set, cached variables with name `NAME` are not over-written.
;                           Instead, the output variable will have "_#" appended, were "#"
;                           represents a unique number.
;       S_MIN:          in, optional, type=float
;                       Values of `S_poynting` below this value will be set to NaN
;       SPHERE:         in, optional, type=boolean, default=0
;                       If set, poynting flux is returned in spherical coordinates.
;       WINDOW:         in, type=string/fltarr(`NFFT`)
;                       A tapering window or the name of a tapering window.
;                       
; :Returns:
;       S:              out, type=objarr(3)
;                       The poynting vector in the spectral domain. Resulting units
;                           are in micro-Watts per square meter (uW / m^2). Three
;                           MrTimeSeries objects are returned, one for each component.
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
;       2016/12/08  -   Written by Matthew Argall
;-
FUNCTION MrVar_PoyntingSpectra, E, B, nfft, nshift, $
CACHE=cache, $
DRIFT=drift, $
NAME=name, $
NO_CLOBBER=no_clobber, $
S_MIN=s_min, $
B_MIN=b_min, $
SPHERE=sphere, $
WINDOW=win
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	tf_drift  = Keyword_Set(drift)
	tf_sphere = Keyword_Set(sphere)
	IF N_Elements(name)   EQ 0 THEN name   = tf_drift ? 'Drift_Velocity_Spectra' : 'Poynting_Spectra'
	IF N_Elements(nfft)   EQ 0 THEN nfft   = 256
	IF N_Elements(nshift) EQ 0 THEN nshIFt = nfft / 4
	
	;Get variables
	oB = MrVar_Get(B)
	oE = MrVar_Get(E)
	
	;Interpolate
	IF ~oE -> IsTimeIdentical(oB) $
		THEN oE = oE -> Interpol(oB)

;-----------------------------------------------------
; FFT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Take the FFT
	E_fft = oE -> FFT( nfft, nshift, /ONE_SIDED, WINDOW=win )
	B_fft = oB -> FFT( nfft, nshift, /ONE_SIDED, WINDOW=win )
	
	pE = E_fft['PTR']
	pB = B_fft['PTR']

;-----------------------------------------------------
; Calculate Poynting Vector \\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;The coefficient 1 / (4*mu_0)
	;   Poynting Flux strenGThs reported in Loto'aniu ranged from 1-150 uW/m^2, so
	;   convert to microWeber per square meter.
	;
	;   In MKS units, to convert to micro-W/m^2, muLTiply by 1e-6
	;
	;   IF not in MKS units, mV/m * nT resuLTs an answer 1e+12 too big. Therefore, muLTiply
	;   by 1e-12 to convert to W/m^2 and another 1e+6 to convert to micro-W/m^2. In total,
	;   that means we have to muLTiply by 1e-6.
	
	
	IF tf_drift THEN BEGIN
		dims  = Size(*pB, /DIMENSIONS)
		bsqr  = Total(Real_Part(Conj(*pB)*(*pB)), 3)
		coeff = Rebin(1e3 / bsqr, dims)
	ENDIF ELSE BEGIN
		coeff = 1.0 / (4.0 * MrConstants('mu_0') * 1.0e6)
	ENDELSE

	;Sx = 1/(4u) * [(Ey*Bz - Ez*By) + (EyBz* - EzBy*)]
	Sx = Conj((*pE)[*,*,1]) * (*pB)[*,*,2] - Conj((*pE)[*,*,2]) * (*pB)[*,*,1] + $
	     (*pE)[*,*,1] * Conj((*pB)[*,*,2]) - (*pE)[*,*,2] * Conj((*pB)[*,*,1])
	
	;Sy = 1/(4u) * [(Ez*Bx - Ex*Bz) + (EzBx* - ExBz*)]
	Sy = Conj((*pE)[*,*,2]) * (*pB)[*,*,0] - Conj((*pE)[*,*,0]) * (*pB)[*,*,2] + $
	     (*pE)[*,*,2] * Conj((*pB)[*,*,0]) - (*pE)[*,*,0] * Conj((*pB)[*,*,2])

	;Sz = 1/(4u) * [(Ex*By - Ey*Bx) + (ExBy* - EyBx*)]
	Sz = Conj((*pE)[*,*,0]) * (*pB)[*,*,1] - Conj((*pE)[*,*,1]) * (*pB)[*,*,0] + $
	     (*pE)[*,*,0] * Conj((*pB)[*,*,1]) - (*pE)[*,*,1] * Conj((*pB)[*,*,0])
	
	;MuLTiply by the coefficient
	Sx *= coeff
	Sy *= coeff
	Sz *= coeff

;-----------------------------------------------------
; Threshold \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Set values |S_poynting| < S_MIN to !values.f_nan
	IF N_Elements(s_min) GT 0 THEN BEGIN
		iBad = Where(Abs(Sx) LT s_min, nBad)
		IF nBad GT 0 THEN Sx[iBad] = !values.f_nan
		
		iBad = Where(Abs(Sy) LT s_min, nBad)
		IF nBad GT 0 THEN Sy[iBad] = !values.f_nan
		
		iBad = Where(Abs(Sz) LT s_min, nBad)
		IF nBad GT 0 THEN Sz[iBad] = !values.f_nan
	ENDIF
	
	IF N_Elements(b_min) GT 0 THEN BEGIN
		iBad = Where(Rebin(bsqr, dims) LT b_min, nBad)
		IF nBad GT 0 THEN BEGIN
			Sx[iBad] = !values.f_nan
			Sy[iBad] = !values.f_nan
			Sz[iBad] = !values.f_nan
		ENDIF
	ENDIF

;-----------------------------------------------------
; Spherical Coordinates \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	param = tf_drift ? 'V$\downD$' : 'S'
	
	IF tf_sphere THEN BEGIN
		;Allocate memory
		dims   = size(Sx, /DIMENSIONS)
		Sx_hat = fltarr(dims)
		Sy_hat = fltarr(dims)
		Sz_hat = fltarr(dims)
		
		;Magnitude
		Smag = Sqrt( Sx^2 + Sy^2 + Sz^2 )
		
		;Normalize
		iZero = Where(Smag eq 0, nZero, COMPLEMENT=iFinite, NCOMPLEMENT=nFinite)
		IF nFinite GT 0 THEN BEGIN
			Sx_hat[iFinite] = Sx[iFinite] / Smag[iFinite]
			Sy_hat[iFinite] = Sy[iFinite] / Smag[iFinite]
			Sz_hat[iFinite] = Sz[iFinite] / Smag[iFinite]
		ENDIF
		IF nZero GT 0 THEN BEGIN
			Sx_hat[iZero] = !values.f_nan
			Sy_hat[iZero] = !values.f_nan
			Sz_hat[iZero] = !values.f_nan
		ENDIF

		;Radial, Azimuthal, Polar
		Sx = Temporary(Smag)
		Sy = !RaDeg * ATan( Temporary(Sy_hat), Temporary(Sx_hat) )
		Sz = !RaDeg * ACos( Temporary(Sz_hat) )
		
		;Axis labels
		desc   = 'Poynting Flux'
		labels = ['|'+param+'|', param+'$\phi$', param+'$\theta$']
		title  = labels + '!C' + '($\mu$W/m$\up2$)'
		units  = '\mu W/m^2'
	ENDIF ELSE BEGIN
		desc   = 'Drift Velodity'
		labels = param + ['x', 'y', 'z']
		title  = labels + '!C' + '(km/s)'
		units  = 'km/s'
	ENDELSE

;-----------------------------------------------------
; Variables \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Component labels
	comps = tf_sphere ? ['r', 'phi', 'theta'] : ['x', 'y', 'z']

	;Create variables.
	oSx = MrTimeSeries( B_fft['TIMEVAR'], Real_Part(Temporary(Sx)), $
	                    CACHE      = cache, $
	                    NAME       = name + '_' + comps[0], $
	                    NO_CLOBBER = no_clobber  )
	oSy = MrTimeSeries( B_fft['TIMEVAR'], Real_Part(Temporary(Sy)), $
	                    CACHE      = cache, $
	                    NAME       = name + '_' + comps[1], $
	                    NO_CLOBBER = no_clobber  )
	oSz = MrTimeSeries( B_fft['TIMEVAR'], Real_Part(Temporary(Sz)), $
	                    CACHE      = cache, $
	                    NAME       = name + '_' + comps[2], $
	                    NO_CLOBBER = no_clobber )

;-----------------------------------------------------
; Attributes \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;SX
	oSx['DEPEND_1']      = B_fft['DEPEND_1']
	oSx['CATDESC']       = desc + ' as a function of frequency'
	oSx['MISSING_COLOR'] = 'Grey'
	oSx['MISSING_VALUE'] = !values.f_nan
	oSx['PLOT_TITLE']    = desc + ' Spectrogram'
	oSx['RGB_TABLE']     = 67
	oSx['SCALE']         = 1B
	oSx['UNITS']         = units
	oSx['TITLE']         = title[0]

	;SY
	oSy['DEPEND_1']      = B_fft['DEPEND_1']
	oSy['CATDESC']       = desc + ' as a function of frequency'
	oSy['MISSING_COLOR'] = 'Grey'
	oSy['MISSING_VALUE'] = !values.f_nan
	oSy['PLOT_TITLE']    = desc + ' Spectrogram'
	oSy['RGB_TABLE']     = 67
	oSy['SCALE']         = 1B
	oSy['UNITS']         = units
	oSy['TITLE']         = title[1]

	;SZ
	oSz['DEPEND_1']      = B_fft['DEPEND_1']
	oSz['CATDESC']       = desc + ' as a function of frequency'
	oSz['MISSING_COLOR'] = 'Grey'
	oSz['MISSING_VALUE'] = !values.f_nan
	oSz['PLOT_TITLE']    = desc + ' Spectrogram'
	oSz['RGB_TABLE']     = 67
	oSz['SCALE']         = 1B
	oSz['UNITS']         = units
	oSz['TITLE']         = title[2]

	;SPHERE-specific properties
	IF tf_sphere THEN BEGIN
		oSx['LOG']       = 1
		oSx['RGB_TABLE'] = 13
		oSx['TITLE']     = '|' + param + '|' + '!C' + (tf_drift ? '(km/s)' : '($\mu$W/m$\up2$/Hz)')
		
		oSy['AXIS_RANGE']   = [-180.0, 180.0]
		oSy['TICKINTERVAL'] = 180
		oSy['MINOR']        = 3
		oSy['TITLE']        = param + '$\phi$!C(Deg)'
		oSy['UNITS']        = 'degrees'
		
		oSz['AXIS_RANGE']   = [0.0, 180.0]
		oSz['TICKINTERVAL'] = 90.0
		oSz['RGB_TABLE']    = 13
		oSz['TITLE']        = param + '$\theta$!C(Deg)'
		oSz['UNITS']        = 'degrees'
		
		;Angles do not mean anything without power
		IF N_Elements(s_min) EQ 0 THEN BEGIN
			s_min = 1e-6
			iBad = oSx -> Where(s_min, /LEQ, COUNT=nBad)

			IF nBad GT 0 THEN BEGIN
				oSy[iBad]  = !Values.f_nan
				oSz[iBad]  = !Values.f_nan
				oSy['NAN'] = 1B
				oSz['NAN'] = 1B
			ENDIF
		ENDIF
	
	;CARTESIAN-specific properties
	ENDIF ELSE BEGIN
		;On a linear scale, the flux is not clearly visible unless it is scaled significantly
		;   TODO: Take log of absolute value, then negate negative values
		smax   = Abs( Min( [oSx.min, oSy.min, oSz.min] ) ) > Abs( Max( [oSx.max, oSy.max, oSz.max] ) )
		srange = 1e-5 * [-smax, smax]
		oSx['AXIS_RANGE'] = srange
		oSy['AXIS_RANGE'] = srange
		oSz['AXIS_RANGE'] = srange
	ENDELSE

	RETURN, [oSx, oSy, oSz]
END