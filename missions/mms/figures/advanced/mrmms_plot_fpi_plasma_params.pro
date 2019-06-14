; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_PolPoynt
;
;*****************************************************************************************
;   Copyright (c) 2017, Matthew Argall                                                   ;
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
;   Plot the parallel and perpendicular components of the electric field.
;        1. Bx, By, Bz, |B|
;        2. Density
;        3. Temprature
;        4. Plasma Beta
;        5. Alfven speed
;        6. Inertial length
;        7. Ratio of plasma to cyclotron frequencies
;
; :Categories:
;   MMS, EDI, MrVariable
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;       TSTART:     in, optional, type=string/strarr
;                   Date-time string, formatted as YYYY-MM-DDThh:mm:ss.fff indicating the
;                       time of the first distribution. If an array, `TEND` and `TSTRIDE`
;                       are ignored. If not provided, the first time stamp will be used.
;       TEND:       in, optional, type=string/integer, default=1
;                   Date-time string, formatted as YYYY-MM-DDThh:mm:ss[.fff] indicating
;                       the time of the last distribution. If an integer, the total
;                       number of distributions to be plotted.
;       TSTRIDE:    in, optional, type=integer, defualt=1
;                   The number of distributions to skip between successive plots.
;
; :Keywords:
;       DT:                 in, optional, type=string/objref, default=2.5
;                           Duration, in seconds, of each time bin.
;       DGA:                in, optional, type=float, default=11.25
;                           Width, in degrees, of each gyrophase bin.
;       FAC:                in, optional, type=string, default='EXB' for srvy and 'VXB' for brst
;                           Name of the field-aligned coordinate system used
;                               to define the directions perpendicular to B.
;                               Options include 'VXB' and 'EXB'
;       GARANGE:            out, optional, type=string, default=[-180\, 180]
;                           Range of gyrophase angles over which to bin.
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
;       2017/05/04  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_FPI_Plasma_Params, sc, mode, species, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
NO_LOAD=no_load, $
TRANGE=trange
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF N_Elements(win) GT 0 THEN Obj_Destroy, win
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(species) EQ 0 THEN species  = 'e'
	IF N_Elements(coords)  EQ 0 THEN coords   = 'gse'
	IF N_Elements(level)   EQ 0 THEN level    = 'l2'
	IF N_Elements(trange)  GT 0 THEN MrVar_SetTRange, trange
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	;FGM parameters
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	fgm_coords = (coords EQ 'dsl'  || coords EQ 'dbcs') ? 'dmpa' : coords
	fgm_mode   = (mode   EQ 'fast' || mode   EQ 'slow') ? 'srvy' : mode
	
	;FPI parameters
	fpi_instr  = 'd' + species + 's'
	fpi_coords = (coords EQ 'dsl' || coords EQ 'dmpa') ? 'dbcs' : coords

	;Load the data
	IF tf_load THEN BEGIN
		;FPI
		MrMMS_FPI_Load_Data, sc, instr, $
		                     OPTDESC   = 'd' + species + 's-moms', $
		                     VARFORMAT = ['*numberdensity_'+mode, '*temptensor_'+coords+'*', $
		                                  '*tempp*', '*prestensor_'+coords+'*']
		
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = level, $
		                     VARFORMAT = '*b_'+fgm_coords+'_'+fgm_mode+'*'
	ENDIF

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	b_vname      = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, mode, level], '_' )
	bmag_vname   = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, mode, level], '_' )
	bvec_vname   = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, mode, level], '_' )
	n_vname      = StrJoin( [sc, fpi_instr, 'numberdensity',             mode], '_' )
	t_vname      = StrJoin( [sc, fpi_instr, 'temptensor',    fpi_coords, mode], '_' )
	t_par_vname  = StrJoin( [sc, fpi_instr, 'temppar',                   mode], '_' )
	t_perp_vname = StrJoin( [sc, fpi_instr, 'tempperp',                  mode], '_' )
	p_vname      = StrJoin( [sc, fpi_instr, 'prestensor',    fpi_coords, mode], '_' )
	
	;Output names
	t_scl_vname  = StrJoin( [sc, fpi_instr, 'scalartemp',      mode], '_' )
	p_scl_vname  = StrJoin( [sc, fpi_instr, 'scalarpres',      mode], '_' )
	va_vname     = StrJoin( [sc, fpi_instr, 'alfvenspeed',     mode], '_' )
	lambda_vname = StrJoin( [sc, fpi_instr, 'intertiallength', mode], '_' )
	beta_vname   = StrJoin( [sc, fpi_instr, 'beta',            mode], '_' )
	ww_vname     = StrJoin( [sc, fpi_instr, 'wc_wp',           mode], '_' )

;-------------------------------------------
; Calculate Plasma Parameters //////////////
;-------------------------------------------
	mu0    = MrConstants('mu_0')
	c      = MrConstants('c')
	mass   = species EQ 'e' ? 'm_e' : 'm_p'
	oB_mag = MrVar_Get(bmag_vname)
	oB_xyz = MrVar_Get(bvec_vname)
	oN     = MrVar_Get(n_vname)
	oT     = MrVar_Get(t_vname)
	oP     = MrVar_Get(p_vname)
	
	IF species EQ 'e' THEN BEGIN
		oB = oB_mag -> Interpol(oN)
	ENDIF ELSE BEGIN
		oB = oB_mag
		oN = oN -> Interpol(oB_mag)
		oT = oT -> Interpol(oB_mag)
		oP = oP -> Interpol(oB_mag)
	ENDELSE
	
	;Scalar pressure & temperature
	oPp = (oP[*,0,0] + oP[*,1,1] + oP[*,2,2]) / 3.0
	oTp = (oT[*,0,0] + oT[*,1,1] + oT[*,2,2]) / 3.0
	
	;Plasma and cyclotron frequencies
	f_p = MrVar_Freq_Plasma(oN, mass)
	f_c = MrVar_Freq_Cyclotron(oB_mag, mass)
	
	;Plasma Beta
	;   - Beta = P_plasma / P_magnetic
	oPb   = 1e-9 * oB^2 / (2*mu0)  ;Magnetic pressure (nPa)
	oBeta = oPp / oPb               ;Plasma beta
	
	;Inertial length
	;   - Lambda = c / w_pe
	;   - 1e-3 converts to km/s
	oLambda = 1e-3 * c / (2*!pi*f_p)
	
	;Alfven Speed
	;   - 1e-6 to m/s, 1e-3 to km/s
	oVA = 1e-15 * oB / Sqrt( mu0 * MrConstants(mass) * oN['DATA'] )
	
	;Ratio of plasma to cyclotron frequencies
	ofp_fc = f_p / f_c

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	param_type = species EQ 'i' ? 'Ion' : 'Electron'
	
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase(sc) + ' ' + param_type + ' Parameters'
	
	;DENSITY
	oN['PLOT_TITLE']    = 'Density'
	oN['SI_CONVERSION'] = '1e-9>Pa'
	oN['TITLE']         = 'N$\down' + species + '$!C(cm$\up-3$)'
	oN['UNITS']         = 'cm^-3'
	
	;SCALAR PRESSURE
	oPp -> SetName, p_scl_vname
	oPp -> Cache
	oPp['PLOT_TITLE']    = 'Scalar Pressure'
	oPp['SI_CONVERSION'] = '1e-9>Pa'
	oPp['TITLE']         = 'P$\down' + species + '$!C(nPa)'
	oPp['UNITS']         = 'nPa'
	
	;SCALAR TEMPERATURE
	oTp -> SetName, t_scl_vname
	oTp -> Cache
	oTp['PLOT_TITLE']    = 'Scalar Temperature'
	oTp['SI_CONVERSION'] = '>'
	oTp['TITLE']         = 'T$\down' + species + '$!C(eV)'
	oTp['UNITS']         = 'eV'
	
	;ALFVEN SPEED
	oVA -> SetName, va_vname
	oVA -> Cache
	oVA['PLOT_TITLE']    = 'Alfven Speed'
	oVA['SI_CONVERSION'] = '1e3>m/s'
	oVA['TITLE']         = 'v$\downA' + species + '$!C(km/s)'
	oVA['UNITS']         = 'km/s'
	
	;INERTIAL LENGTH
	oLambda -> SetName, lambda_vname
	oLambda -> Cache
	oLambda['PLOT_TITLE']    = 'Inertial Length'
	oLambda['SI_CONVERSION'] = '1e3>m/s'
	oLambda['TITLE']         = '$\lambda$$\down' + species + '$!C(km)'
	oLambda['UNITS']         = 'km'
	
	;PLASMA BETA
	oBeta -> SetName, beta_vname
	oBeta -> Cache
	oBeta['LOG']           = 1B
	oBeta['PLOT_TITLE']    = 'Plasma Beta'
	oBeta['SI_CONVERSION'] = '>'
	oBeta['TITLE']         = '$\beta$$\down' + species + '$'
	oBeta['UNITS']         = ''
	
	;RATIO OF CYCLOTRON & PLASMA FREQUENCIES
	ofp_fc -> SetName, ww_vname
	ofp_fc -> Cache
	ofp_fc['PLOT_TITLE']    = 'Ratio of Plasma to Cyclotron Freqencies'
	ofp_fc['SI_CONVERSION'] = '>'
	ofp_fc['TITLE']         = 'f$\downc' + species + '$/f$\downp' + species + '$'
	ofp_fc['UNITS']         = ''

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [b_vname, n_vname, t_scl_vname, p_scl_vname, beta_vname, va_vname, lambda_vname, ww_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 650, $
	                    YSIZE = 800 )
	
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[14,5]
	win    -> refresh
	return, win
end