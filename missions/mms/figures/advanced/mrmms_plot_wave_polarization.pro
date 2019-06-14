; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_Wave_Polarization
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
;       1. Bx, By, Bz, |B|
;       2. Detrended field (FGM, SCM, or EDP)
;       3. Intnesity
;       4. Percent Polarization
;       5. Ellipticity
;       6. Wave Normal Angle
;       7. Coherency
;
; :Categories:
;   MrVariable, MMS
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;       NFFT:       in, optional, type=integer, default=512
;                   Number of points per FFT.
;       NSHIFT:     in, optional, type=integer, default=`NFFT`/2
;                   Number of points to shift after each FFT.
;
; :Keywords:
;       INSTR:      in, optional, type=string, default='fgm'
;                   Instrument to use. Options are: { 'afg' | 'dfg' | 'fgm', 'edp', 'scm' }
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
;       2017/05/04  -   Written by Matthew Argall
;       2017/06/06  -   Apply window function to FFT. - MRA
;       2018/04/18  -   Added the LEVEL keyword. - MRA
;-
FUNCTION MrMMS_Plot_Wave_Polarization, sc, instr, mode, nfft, nshift, $
LEVEL=level, $
NDETREND=nDetrend, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
P_MIN=p_min, $
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
	IF N_Elements(nfft)     EQ 0 THEN nfft     = 512
	IF N_Elements(nshift)   EQ 0 THEN nshift   = nfft / 2.0
	IF N_Elements(nDetrend) EQ 0 THEN nDetrend = nfft
	IF N_Elements(coords)   EQ 0 THEN coords   = 'gse'
	IF N_Elements(level)    EQ 0 THEN level    = 'l2'
	IF N_Elements(trange)   GT 0 THEN MrVar_SetTRange, trange
	IF N_Elements(p_min)    EQ 0 THEN p_min = 0.0
	
	;Conflicts
	IF Array_Equal(instr EQ ['afg', 'dfg', 'fgm', 'edp', 'scm'], 0) $
		THEN Message, 'Invalid value for instr: "' + instr + '".'
	
;-------------------------------------------
; Variable Name Parameters /////////////////
;-------------------------------------------
	;FGM parameters
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = (instr NE 'edp' && instr NE 'scm') ? instr : 'fgm'
	fgm_coords = (coords EQ 'dsl'  || coords EQ 'dbcs') ? 'dmpa' : coords
	IF level EQ 'ql' THEN fgm_coords = 'dmpa'
	fgm_mode   = (mode   EQ 'fast' || mode   EQ 'slow') ? 'srvy' : mode
	
	;EDP parameters
	edp_coords = (coords EQ 'dmpa' || coords EQ 'dbcs') ? 'dsl'  : coords
	IF mode EQ 'srvy' THEN BEGIN
		MrPrintF, 'LogWarn', 'EDP does not have SRVY data. Using FAST.'
		edp_mode = 'fast'
	ENDIF ELSE BEGIN
		edp_mode = mode
	ENDELSE
	
	;SCM parameters
	IF N_Elements(scm_optdesc) EQ 0 THEN BEGIN
		CASE mode OF
			'slow': scm_optdesc = 'scf'
			'fast': scm_optdesc = 'scf'
			'srvy': scm_optdesc = 'scsrvy'
			'brst': scm_optdesc = 'scb'  ; or 'sch'
			ELSE: Message, 'Invalid MODE: "' + mode + '".'
		ENDCASE
	ENDIF
	scm_coords = coords
	IF scm_coords NE 'gse' THEN Message, 'SCM data is available only in GSE.'
	
	;FPI Parameters
	fpi_mode    = mode EQ 'brst' ? mode : 'fast'
	fpi_optdesc = level EQ 'ql' ? 'des' : 'des-moms'
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	
	;Load the data
	IF tf_load THEN BEGIN
		;EDP
		IF instr EQ 'edp' THEN BEGIN
			MrMMS_Load_Data, sc, instr, edp_mode, level, $
			                 OPTDESC   = 'dce', $
			                 VARFORMAT = '*_dce_'+edp_coords+'_*'
		
		;SCM
		ENDIF ELSE IF instr EQ 'scm' THEN BEGIN
			MrMMS_Load_Data, sc, instr, 'brst', 'l2', $
			                 OPTDESC = scm_optdesc
		
		;FGM
		ENDIF ELSE IF instr NE fgm_instr THEN BEGIN
			MrMMS_FGM_Load_Data, sc, fgm_mode, $
			                     INSTR     = instr, $
			                     LEVEL     = level, $
			                     VARFORMAT = '*b_'+fgm_coords+'_'+fgm_mode
		ENDIF
		
		;FPI
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     LEVEL     = level, $
		                     OPTDESC   = fpi_optdesc, $
		                     VARFORMAT = '*numberdensity_'+fpi_mode
		
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = level, $
		                     VARFORMAT = ['*'+fgm_mode+'_'+fgm_coords,'*b_'+fgm_coords+'_'+fgm_mode+'*']
	ENDIF
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	IF level EQ 'ql' THEN BEGIN
		b_vname    = StrJoin( [sc, fgm_instr,        mode, fgm_coords], '_' )
		bvec_vname = StrJoin( [sc, fgm_instr, 'vec', mode, fgm_coords], '_' )
		bmag_vname = StrJoin( [sc, fgm_instr, 'mag', mode, fgm_coords], '_' )
	ENDIF ELSE BEGIN
		b_vname    = StrJoin( [sc, fgm_instr, 'b',    fgm_coords,              mode, level], '_' )
		bmag_vname = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords,              mode, level], '_' )
		bvec_vname = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords,              mode, level], '_' )
	ENDELSE
	fgm_vname  = bvec_vname
	scm_vname  = StrJoin( [sc, instr,     'acb',  scm_coords, scm_optdesc, mode, level], '_' )
	edp_vname  = StrJoin( [sc, instr,     'dce',  edp_coords,              edp_mode, level], '_' )
	ne_vname   = StrJoin( [sc, 'des',     'numberdensity',                 fpi_mode       ], '_' )
	
	CASE instr OF
		'edp': vin_vname = edp_vname
		'scm': vin_vname = scm_vname
		ELSE:  vin_vname = fgm_vname
	ENDCASE
	
	;Output names
	CASE instr OF
		'edp': dv_vname = StrJoin( [sc, instr, 'dce',  'dfac',               edp_mode, level], '_' )
		'scm': dv_vname = StrJoin( [sc, instr, 'acb',  'dfac',  scm_optdesc, mode,     level], '_' )
		ELSE:  dv_vname = StrJoin( [sc, instr, 'bvec', 'dfac',               fgm_mode, level], '_' )
	ENDCASE
	coherency_vname    = 'Coherency('       + dv_vname + ')'
	intensity_vname    = 'Intensity('       + dv_vname + ')'
	wavenormal_vname   = 'WaveNormalAngle(' + dv_vname + ')'
	ellipticity_vname  = 'Ellipticity('     + dv_vname + ')'
	polarization_vname = 'Polarization('    + dv_vname + ')'
	fce_vname          = StrJoin( [sc, 'plasma', 'fce',         mode, level], '_' )
	fce2_vname         = StrJoin( [sc, 'plasma', 'fce', 'half', mode, level], '_' )
	fci_vname          = StrJoin( [sc, 'plasma', 'fci',         mode, level], '_' )
	flh_vname          = StrJoin( [sc, 'plasma', 'flh',         mode, level], '_' )
	fpe_vname          = StrJoin( [sc, 'plasma', 'fpe',         mode, level], '_' )
	fpi_vname          = StrJoin( [sc, 'plasma', 'fpi',         mode, level], '_' )

;-------------------------------------------
; Detrend & FAC ////////////////////////////
;-------------------------------------------
	oB   = MrVar_Get(bvec_vname)
	oVec = MrVar_Get(vin_vname)
	
	;Remove DC field
	oDVec = oVec -> Detrend(nDetrend, /EDGE_TRUNCATE)
	
	;Background field
	IF StRegEx(instr, '(afg|dfg|fgm)', /BOOLEAN) THEN BEGIN
		oB0 = oVec - oDVec
	ENDIF ELSE BEGIN
		nSmooth = Round( nDetrend * oVec['TIMEVAR'] -> GetSI() / oB['TIMEVAR'] -> GetSI() )
		oB0  = oB -> Copy()
		oB0 -> SetData, Smooth( oB0['DATA'], [nSmooth, 0], /EDGE_TRUNCATE )
		oB0  = oB0 -> Interpol(oVec)
	ENDELSE
	
	;Field-Aligned Coordinates
	oT       = MrVar_FAC(oB0)
	oVec_FAC = oT ## oDVec
	
	;Name and cache
	oVec_FAC -> SetName, dv_vname
	oVec_FAC -> Cache
	
	;Free memory
	Obj_Destroy, [oB0, oDVec, oT]

;-------------------------------------------
; Polarization Analysis ////////////////////
;-------------------------------------------
	MrVar_Polarization, oVec_FAC, nfft, nshift, $
	                    /CACHE, $
	                    NFAVG        = 3, $
	                    COHERENCY    = oCoh, $
	                    ELLIPTICITY  = oEllip, $
	                    INTENSITY    = oIntensity, $
	                    KDOTB_ANGLE  = oWaveNormal, $
	                    POLARIZATION = oPol, $
	                    WINDOW       = 'hanning'

;-------------------------------------------
; Filter By % Polarization /////////////////
;-------------------------------------------
	IF p_min GT 0.0 THEN BEGIN
		iBad = oPol -> Where(p_min, /LESS, COUNT=nBad)
		IF nBad GT 0 THEN BEGIN
			oCoh[iBad]        = !values.f_nan
			oEllip[iBad]      = !values.f_nan
			oIntensity[iBad]  = !values.f_nan
			oWaveNormal[iBad] = !values.f_nan
			oPol[iBad]        = !values.f_nan
		ENDIF
	ENDIF

;-------------------------------------------
; Plasma Frequencies ///////////////////////
;-------------------------------------------
	;Plasma frequencies
	oFci = MrVar_Freq_Cyclotron(bmag_vname, 'm_p', /CACHE, NAME=fci_vname)
	oFce = MrVar_Freq_Cyclotron(bmag_vname, 'm_e', /CACHE, NAME=fce_vname)
	oFpi = MrVar_Freq_Plasma(ne_vname, 'm_p', /CACHE, NAME=fpi_vname)
	oFpe = MrVar_Freq_Plasma(ne_vname, 'm_e', /CACHE, NAME=fpe_vname)
	oFlh = MrVar_Freq_LowerHybrid(bmag_vname, ne_vname, /CACHE, NAME=flh_vname)
	
	;Half the electron cyclotron frequency
	oFce2 = oFce / 2.0
	oFce2 -> SetName, fce2_vname
	oFce2 -> Cache
	
	;Determine which fit in the graph
	oplots = []
	oFreq  = (MrVar_Get(intensity_vname))['DEPEND_1']
	frange = [oFreq.min, oFreq.max]
	IF ~(oFci.min  GE frange[1] || oFci.max  LE frange[0]) THEN oplots = [ oplots, oFci  ]
	IF ~(oFce.min  GE frange[1] || oFce.max  LE frange[0]) THEN oplots = [ oplots, oFce  ]
	IF ~(oFce2.min GE frange[1] || oFce2.max LE frange[0]) THEN oplots = [ oplots, oFce2 ]
	IF ~(oFpi.min  GE frange[1] || oFpi.max  LE frange[0]) THEN oplots = [ oplots, oFpi  ]
	IF ~(oFpe.min  GE frange[1] || oFpe.max  LE frange[0]) THEN oplots = [ oplots, oFpe  ]
	IF ~(oFlh.min  GE frange[1] || oFlh.max  LE frange[0]) THEN oplots = [ oplots, oFlh  ]

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	oB['PLOT_TITLE'] = StrUpCase(sc) + '!C' + String(nfft, nshift, ndetrend, FORMAT='(%"NFFT=%i NSHIFT=%i NDETREND=%i")')

	IF instr EQ 'scm' THEN BEGIN
		oVec['LABEL'] = ['B$\downX$', 'B$\downY$', 'B$\downZ$']
		oVec['TITLE'] = 'B!C(nT)'
		
		oVec_FAC['LABEL'] = ['dB$\down||$', 'dB$\downPerp1$', 'dB$\downPerp2$']
		oVec_FAC['TITLE'] = 'dB!C(nT)'
	ENDIF ELSE IF instr EQ 'edp' THEN BEGIN
		oVec['LABEL'] = ['E$\downX$', 'E$\downY$', 'E$\downZ$']
		oVec['TITLE'] = 'E!C(mV/m)'
		
		oVec_FAC['LABEL'] = ['dE$\down||$', 'dE$\downPerp1$', 'dE$\downPerp2$']
		oVec_FAC['TITLE'] = 'dE!C(mV/m)'
	ENDIF
	
	oFci['COLOR'] = 'White'
	oFci['LABEL'] = 'fci'
	oFci['THICK'] = 2
	
	oFce['COLOR'] = 'Black'
	oFce['LABEL'] = 'fce'
	oFce['THICK'] = 2
	
	oFpi['COLOR'] = 'Cyan'
	oFpi['LABEL'] = 'fpi'
	oFpi['THICK'] = 2
	
	oFpe['COLOR'] = 'Magenta'
	oFpe['LABEL'] = 'fpe'
	oFpe['THICK'] = 2
	
	oFlh['COLOR'] = 'Yellow'
	oFlh['LABEL'] = 'flh'
	oFlh['THICK'] = 2
	
	oFce2['COLOR']     = 'Black'
	oFce2['LABEL']     = 'fce/2'
	oFce2['LINESTYLE'] = '--'
	oFce2['THICK']     = 2

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [bvec_vname, dv_vname, intensity_vname, polarization_vname, $
	                     ellipticity_vname, wavenormal_vname, coherency_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	
	;Overplot plasma frequencies
	win = MrVar_OPlotTS( intensity_vname, oplots )
	FOR i = 0, N_Elements(oplots) - 1 DO oplots[i] -> RemoveAttr, 'LABEL'
	
	win = MrVar_OPlotTS( polarization_vname, oplots )
	win = MrVar_OPlotTS( ellipticity_vname,  oplots )
	win = MrVar_OPlotTS( wavenormal_vname,   oplots )
	win = MrVar_OPlotTS( coherency_vname,    oplots )

	;Set legend properties
	lgds = win -> Get(/ALL, COUNT=nLgds, ISA='MrLegend')
	FOR i = 0, nLgds - 1 DO BEGIN
		IF StRegEx(lgds[i].name, 'f(ce|ci|pe|pi|lh)', /BOOLEAN) THEN BEGIN
			lgds[i] -> SetProperty, ALIGNMENT  ='NE', $
			                        COLOR      = 'Black', $
			                        FILL_COLOR ='Light Grey', $
			                        LINESTYLE  ='-'
		ENDIF
	ENDFOR
	
	;Prettify
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win -> refresh

;-------------------------------------------
; Save Figure //////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
			MrPrintF, 'LogText', 'Saving file to: "' + output_dir + '".'
		ENDIF
		
		;File name
		sn    = String(nfft,     FORMAT='(%"%infft")')
		ss    = String(nshift,   FORMAT='(%"%inshift")')
		sd    = String(ndetrend, FORMAT='(%"%indetrend")')
		fname = StrJoin( [sc, instr, mode, level, 'pol-'+sn+'-'+ss+'-'+sd], '_' )
		fname = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, win
END