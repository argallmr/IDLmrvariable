; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_DSP
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
;   Generate a plot of FGM quantities:
;       1. Bx, By, Bz, |B|
;       2. Bx PSD + fce, fce/2 lines
;       3. By PSD + fce, fce/2 lines
;       4. Bz PSD + fce, fce/2 lines
;       5. B Omni PSD + fce, fce/2 lines
;       6. Ex PSD + fce, fce/2 lines
;       7. Ey PSD + fce, fce/2 lines
;       8. Ez PSD + fce, fce/2 lines
;       9. E Omni PSD + fce, fce/2 lines
;
; :Categories:
;   MMS
;
; :Params:
;       MODE:       in, required, type=string, default='fast'
;                   Data rate mode. Options are {'slow' | 'fast'}
;       INSTR:      in, required, type=string, default='fgm'
;                   FGM strument to use. Options are: {'afg' | 'dfg' | 'fgm'}
;
; :Keywords:
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;       TRANGE:     in, optional, type=string/strarr(2), default=MrVar_GetTRange()
;                   The start and end times of the data interval to be plotted, formatted
;                       as 'YYYY-MM-DDThh:mm:ss'
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
;       2017/01/17  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_Mirror_Overview, sc, instr, mode, $
COORDS=coords, $
LEVEL=level, $
FGM_INSTR=fgm_instr, $
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
	
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(level)     EQ 0 THEN level      = 'l2'
	IF N_Elements(coords)    EQ 0 THEN coords     = 'gse'
	IF N_Elements(mode)      EQ 0 THEN mode       = 'fast'
	IF N_Elements(instr)     EQ 0 THEN instr      = 'scm'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr  = 'fgm'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;FGM
	fgm_mode = 'srvy'
	CASE fgm_instr OF
		'afg': fgm_level = 'l2pre'
		'dfg': fgm_level = 'l2pre'
		'fgm': fgm_level = 'l2'
		ELSE: Message, 'Invalid FGM instrument: "' + instr + '".'
	ENDCASE
	CASE coords OF
		'dsl':  fgm_coords = 'dmpa'
		'dmpa': fgm_coords = 'dmpa'
		'dbcs': fgm_coords = 'dmpa'
		'gse':  fgm_coords = coords
		'gsm':  fgm_coords = coords
		ELSE: Message, 'FPI does not provide data in "' + coords + '" coordinates.'
	ENDCASE
	
	;DSP
	dsp_mode = mode EQ 'brst' ? mode : ['fast', 'slow']
	
	;SCM
	edp_mode = mode EQ 'brst' ? mode : 'srvy'
	IF N_Elements(scm_optdesc) EQ 0 THEN BEGIN
		CASE mode OF
			'slow': scm_optdesc = 'scs'
			'fast': scm_optdesc = 'scf'
			'srvy': scm_optdesc = 'scsrvy'
			'brst': scm_optdesc = 'scb'
			ELSE: Message, 'Invalid value for MODE: "' + mode + '".'
		ENDCASE
	ENDIF
	CASE coords OF
		'gse':  scm_coords = coords
		ELSE: Message, 'SCM does not provide data in "' + coords + '" coordinates.'
	ENDCASE
		
	
	;EDP
	edp_mode = mode EQ 'brst' ? mode : 'fast'
	
	;FPI
	fpi_mode = mode EQ 'brst' ? mode : 'fast'
	CASE coords OF
		'dsl':  fpi_coords = 'dbcs'
		'dmpa': fpi_coords = 'dbcs'
		'dbcs': fpi_coords = coords
		'gse':  fpi_coords = coords
		ELSE: Message, 'FPI does not provide data in "' + coords + '" coordinates.'
	ENDCASE

	;Source names
	b_vname     = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_' )
	bvec_vname  = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_' )
	bmag_vname  = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_' )
	e_edp_vname = ''
	b_scm_vname = StrJoin( [sc, 'scm', 'acb', coords, 'scb', mode, level], '_' )
	bxpsd_vname = StrJoin( [sc, 'dsp', 'bpsd', 'scm1', mode, level], '_' )
	bypsd_vname = StrJoin( [sc, 'dsp', 'bpsd', 'scm2', mode, level], '_' )
	bzpsd_vname = StrJoin( [sc, 'dsp', 'bpsd', 'scm3', mode, level], '_' )
	bpsd_vname  = StrJoin( [sc, 'dsp', 'bpsd', 'omni', mode, level], '_' )
	expsd_vname = StrJoin( [sc, 'dsp', 'epsd', 'x'], '_' )
	eypsd_vname = StrJoin( [sc, 'dsp', 'epsd', 'y'], '_' )
	ezpsd_vname = StrJoin( [sc, 'dsp', 'epsd', 'z'], '_' )
	epsd_vname  = StrJoin( [sc, 'dsp', 'epsd', 'omni'], '_' )
	ne_vname    = StrJoin( [sc, 'des', 'numberdensity', fpi_mode], '_' )
	ni_vname    = StrJoin( [sc, 'dis', 'numberdensity', fpi_mode], '_' )
	
	;Derived names
	omni_hi_vname = StrJoin( [sc, instr,     'psd',     'omni', 'hi', mode, level], '_' )
	omni_lo_vname = StrJoin( [sc, instr,     'psd',     'omni', 'lo', mode, level], '_' )
	fce_vname     = StrJoin( [sc, fgm_instr, 'fce',                   mode, level], '_' )
	fce_2_vname   = StrJoin( [sc, fgm_instr, 'halffce',               mode, level], '_' )
	fch_vname     = StrJoin( [sc, fgm_instr, 'fcH',                   mode, level], '_' )
	fche_vname    = StrJoin( [sc, fgm_instr, 'fcHe',                  mode, level], '_' )
	fco_vname     = StrJoin( [sc, fgm_instr, 'fcO',                   mode, level], '_' )

;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = fgm_level, $
		                     VARFORMAT = b_vname
		
		;Low-Frequency PSD of E-Field
		IF instr EQ 'eqp' THEN BEGIN
			MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
			                 OPTDESC   = edp_optdesc
		ENDIF
		
		;Old naming convention?
		MrVar_Names, names, b_vname
		IF names[0] EQ '' THEN BEGIN
			b_vname     = StrJoin( [sc, fgm_instr,        fgm_mode, fgm_level, fgm_coords], '_' )
			bvec_vname  = StrJoin( [sc, fgm_instr, 'vec', fgm_mode, fgm_level, fgm_coords], '_' )
			bmag_vname  = StrJoin( [sc, fgm_instr, 'mag', fgm_mode, fgm_level, fgm_coords], '_' )
			
			;FGM
			MrMMS_FGM_Load_Data, sc, fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = fgm_level, $
			                     VARFORMAT = b_vname
			
			IF fgm_level NE 'l2' && ~MrVar_IsCached(b_vname) $
				THEN Message, 'FGM L2PRE variable name incorrect.'
		ENDIF
		
		;MAGNETIC FIELD
		IF instr EQ 'scm' THEN BEGIN
			IF mode EQ 'brst' THEN BEGIN
				MrMMS_Load_Data, sc, 'scm', scm_mode, level, $
				                 OPTDESC   = scm_optdesc
			ENDIF ELSE BEGIN
				MrMMS_Load_Data, sc, 'dsp', dsp_mode, level, $
				                 OPTDESC   = 'bpsd'
			ENDELSE
			
		;ELECTRIC FIELD
		ENDIF ELSE IF instr EQ 'edp' THEN BEGIN
			IF mode EQ 'brst' THEN BEGIN
				MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
				                 OPTDESC   = edp_optdesc
			ENDIF ELSE BEGIN
				MrMMS_Load_Data, sc, 'dsp', mode, level, $
				                 OPTDESC   = 'epsd'
			ENDELSE
		
		;OTHER
		ENDIF ELSE BEGIN
			Message, 'Invalid value for INSTR: "' + instr + '".'
		ENDELSE
		
		;FPI
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     LEVEL     = level, $
		                     OPTDESC   = 'des-moms', $
		                     VARFORMAT = ['*density_'+fpi_mode+'*', '*tempperp*', '*temppara*', $
		                                  '*prestensor_'+fpi_coords+'*']
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     LEVEL     = level, $
		                     OPTDESC   = 'dis-moms', $
		                     VARFORMAT = ['*density_'+fpi_mode+'*', '*tempperp*', '*temppara*', $
		                                  '*prestensor_'+fpi_coords+'*']
	ENDIF ELSE BEGIN
		;Old naming convention?
		MrVar_Names, names, b_vname
		IF names[0] EQ '' THEN BEGIN
			b_vname     = StrJoin( [sc, fgm_instr,        fgm_mode, fgm_level, fgm_coords], '_' )
			bvec_vname  = StrJoin( [sc, fgm_instr, 'vec', fgm_mode, fgm_level, fgm_coords], '_' )
			bmag_vname  = StrJoin( [sc, fgm_instr, 'mag', fgm_mode, fgm_level, fgm_coords], '_' )
			
			IF ~MrVar_IsCached(b_vname) $
				THEN Message, 'FGM L2PRE variable name incorrect.'
		ENDIF
	ENDELSE
;-------------------------------------------
; High Frequency ///////////////////////////
;-------------------------------------------
	;BRST
	IF mode EQ 'brst' THEN BEGIN
		;Pick the variables
		IF instr EQ 'edp' $
			THEN vname = e_edp_vname $
			ELSE vname = b_scm_vname
		
		;Compute the Spectrogram
		oV    = MrVar_Get(vname)
		oSpec = oV -> Spectrogram( 4096, 2048, WINDOW='Hamming' )
		
		;Omni-Directional spectrogram
		oOmni  = oSpec[*,*,0] + oSpec[*,*,1] + oSpec[*,*,2]
		oOmni -> SetName, omni_hi_vname
		oOmni -> Cache
		oSpec -> CopyAttrTo, oOmni
		
	;SRVY
	ENDIF ELSE BEGIN
		omni_hi_vname = instr EQ 'edp' ? epsd_vname : bpsd_vname
	ENDELSE
	
;-------------------------------------------
; Low Frequency ////////////////////////////
;-------------------------------------------
	IF mode EQ 'brst' $
		THEN vname = bvec_vname $
		ELSE vname = instr EQ 'edp' ? edp_vname : bvec_vname
	
	;Compute the spectrogram
	oV    = MrVar_Get(vname)
	oSpec = oV -> Spectrogram( 4096, 2048, $
	                           /DETREND, $
	                           WINDOW = 'Hamming' )
	
	;Omni-directional spectrogram
	oOmni  = oSpec[*,*,0] + oSpec[*,*,1] + oSpec[*,*,2]
	oSpec -> CopyAttrTo, oOmni
	oOmni -> SetName, omni_lo_vname
	oOmni -> Cache

;-------------------------------------------
; Gyrofrequency Lines //////////////////////
;-------------------------------------------
	oBmag  = MrVar_Get(bmag_vname)
	
	;FCE
	ofce   = MrVar_Freq_Cyclotron(oBmag, 'm_e',  /CACHE, NAME=fce_vname)
	
	;FCE/2
	ofce_2 = ofce / 2.0
	ofce_2 -> SetName, fce_2_vname
	ofce_2 -> Cache
	
	;FC H+, He+, O+
	ofch  = MrVar_Freq_Cyclotron(oBmag, 'm_p',  /CACHE, NAME=fch_vname)
	ofche = MrVar_Freq_Cyclotron(oBmag, 'm_He', /CACHE, NAME=fche_vname)
	ofco  = MrVar_Freq_Cyclotron(oBmag, 'm_O',  /CACHE, NAME=fco_vname)
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	;BMAG
	oB = MrVar_Get(b_vname)
	oB['AXIS_RANGE'] = [oB.min > (-100), oB.max < 100]
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level], ' ' ) )
	
	;Omin Lo
	oOmni['LOG'] = 1B
	
	;Density
	oNe = MrVar_Get(ne_vname)
	oNi = MrVar_Get(ni_vname)
	range = [ Min( [ oNe['DATA', Where(oNe['DATA'] GT 0)], $
	                 oNi['DATA', Where(oNi['DATA'] GT 0)] ], nmax ), nmax ]
	oNe['AXIS_RANGE'] = range
	oNe['COLOR']      = 'Red'
	oNe['LABEL']      = 'Ne'
	oNe['TITLE']      = 'N!C(cm$\up-3$)'
	
	;PSD
	IF mode EQ 'brst' THEN BEGIN
	
	ENDIF ELSE BEGIN
		IF instr EQ 'scm' THEN BEGIN
			oOmni = MrVar_Get(omni_hi_vname)
			oOmni['AXIS_RANGE'] = [1e-7,1e-1]
			oOmni['LOG']        = 1B
			oOmni['SCALE']      = 1B
		ENDIF
	ENDELSE
	
	;Ion Density
	oNi['AXIS_RANGE'] = range
	oNi['COLOR']      = 'Blue'
	oNi['LABEL']      = 'Ni'
	
	;FCE
	ofce['COLOR']     = 'White'
;	ofce['LINESTYLE'] = '--'
	ofce['NSUM']      = 4
	
	;FCE_2
	ofce_2['COLOR']      = 'Magenta'
;	ofce_2['LINESTYLE'] = '.'
	ofce_2['NSUM']      = 4
	ofce_2['UNITS']     = 'Hz'
	
	;FC H+
	ofch['COLOR']     = 'White'
;	ofce['LINESTYLE'] = '--'
	ofce['NSUM']      = 4
	
	;FC He+
	ofche['COLOR']     = 'Green'
;	ofche['LINESTYLE'] = '--'
	ofche['NSUM']      = 4
	
	;FC O+
	ofco['COLOR']     = 'Magenta'
;	ofco['LINESTYLE'] = '--'
	ofco['NSUM']      = 4

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------

	;Plot data
	win = MrVar_PlotTS( [b_vname, omni_hi_vname, omni_lo_vname, ne_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 680, $
	                    YSIZE = 700 )
	win = MrVar_OPlotTS( omni_hi_vname, [fce_vname, fce_2_vname] )
	win = MrVar_OPlotTS( omni_lo_vname, [fch_vname, fche_vname, fco_vname] )
	win = MrVar_OPlotTS( ne_vname,      ni_vname )

	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 14]
	win    -> Refresh

	RETURN, win
END