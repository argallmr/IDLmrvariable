; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_dBE
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
;   Generate a plot of EDP quantities:
;        1. B barycenter
;        2. dBx from MMS 1, 2, 3, 4
;        3. dBy from MMS 1, 2, 3, 4
;        4. dBz from MMS 1, 2, 3, 4
;        5. dEx from MMS 1, 2, 3, 4
;        6. dEy from MMS 1, 2, 3, 4
;        7. dEz from MMS 1, 2, 3, 4
;
; :Categories:
;   MMS
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' }
;
; :Keywords:
;       COORDS:     in, optional, type=string, default='gse'
;                   Coordinate system in which data is plotted. Options are: {'dsl' | 'gse' | 'fac'}
;       FGM_INSTR:  in, optional, type=string, default='fgm'
;                   FGM instrument to use. Options are: { 'afg' | 'dfg' | 'fgm' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       OPTDESC:    in, optional, type=string, default='dce'
;                   Optional filename descriptor. Options are: {'dce' | 'hmfe'}
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
;       F_PASS:     in, optional, type=fltarr(2), default=[20.0, 0.0]
;                   Frequency range in which the components of E are filtered. A value
;                       of 0.0 in the second element is equivalent to the Nyquist frequency.
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
;       2018/02/23  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_dBE, mode, ffilt, $
COORDS=coords, $
MAG_INSTR=mag_instr, $
LEVEL=level, $
OPTDESC=optdesc, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
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
	IF N_Elements(optdesc)   EQ 0 THEN optdesc   = 'dce'
	IF N_Elements(coords)    EQ 0 THEN coords    = optdesc EQ 'hmfe' ? 'dsl' : 'gse'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(mode)      EQ 0 THEN mode      = 'brst'
	IF N_Elements(mag_instr) EQ 0 THEN mag_instr = mode EQ 'brst' ? 'fsm' : 'scm'
	IF N_Elements(ffilt)     EQ 0 THEN ffilt     = [20.0, 0.0]
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange 
	
;-------------------------------------------
; Dataset Parameters ///////////////////////
;-------------------------------------------
	
	
	;EDP
	;   - OPDESC of 'scpot' is not allowed
	edp_instr  = 'edp'
	edp_coords = (coords EQ 'dmpa' || coords EQ 'dbcs') ? 'dsl' : coords
	IF coords EQ 'fac' THEN edp_coords = 'gse'
	edp_mode   = mode EQ 'brst' ? mode : 'fast'
	IF ~MrIsMember(['dce', 'hmfe'], optdesc) $
		THEN Message, 'Optional descriptor "' + optdesc + '" not allowed. Must be {"dce" | "hmfe"}' $
		ELSE edp_optdesc = optdesc
		
	;FGM
	fgm_coords  = MrIsMember(['dbcs', 'dsl'], coords) ? 'dmpa' : coords
	IF coords EQ 'fac' THEN fgm_coords = 'gse'
	fgm_mode    = mode EQ 'brst' ? mode : 'srvy'
	fgm_optdesc = mag_instr EQ 'fsm' ? '8khz' : ''
	CASE mag_instr OF
		'fsm': fgm_level = 'l3'
		'fgm': fgm_level = 'l2'
		'dfg': fgm_level = 'l2pre'
		'afg': fgm_level = 'l2pre'
		ELSE: Message, 'MAG_INSTR not recognized: "' + mag_instr + '".'
	ENDCASE
	
	;SCM
	scm_coords = 'gse'
	CASE mode OF
		'slow': scm_optdesc = 'scs'
		'fast': scm_optdesc = 'scf'
		'srvy': scm_optdesc = 'scsrvy'
		'brst': scm_optdesc = 'scb'
		ELSE:   Message, 'Invalid operating mode: "' + mode + '".'
	ENDCASE
	IF coords NE 'gse' && coords NE 'fac' THEN BEGIN
		MrPrintF, 'LogWarn', 'SCM data is in GSE coordinates only.'
		scm_coords = 'gse'
	ENDIF
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	sc = 'mms' + ['1', '2', '3', '4']
	sc_colors = ['Black', 'Red', 'Forest Green', 'Blue']
	
	
	b_vnames = sc + '_' + StrJoin([mag_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_')
	IF mag_instr EQ 'fsm' THEN BEGIN
		bvec_vnames = sc + '_' + StrJoin([mag_instr, 'b', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([mag_instr, 'b', 'mag',      fgm_mode, fgm_level], '_')
	ENDIF ELSE IF mag_instr EQ 'scm' THEN BEGIN
		
	ENDIF ELSE BEGIN
		bvec_vnames = sc + '_' + StrJoin([mag_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([mag_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_')
	ENDELSE
	
	;Source names
	e_vnames    = sc + '_' + StrJoin( [edp_instr, edp_optdesc, edp_coords, mode, level], '_' )
	epar_vnames = sc + '_' + StrJoin( [edp_instr, edp_optdesc, 'par', 'epar', mode, level], '_' )
	
	;Inermediate names
	bvec_interp_vnames = sc + '_' + StrJoin( ['fgm', 'binerp', fgm_coords, fgm_mode, fgm_level], '_' )
	e_interp_vnames    = sc + '_' + StrJoin( [edp_instr, 'interp', edp_optdesc, edp_coords, mode, level], '_' )
	
	;Output names
	bmag_bary_vname = StrJoin( ['mms', mag_instr, 'bmag', fgm_coords, 'bary', fgm_mode, fgm_level], '_' )
	bx_vnames    = bvec_vnames + '_x'
	by_vnames    = bvec_vnames + '_y'
	bz_vnames    = bvec_vnames + '_z'
	db_vnames    = sc + '_' + StrJoin( [mag_instr, 'db',  fgm_coords, fgm_mode, fgm_level], '_' )
	dbx_vnames   = db_vnames + '_x'
	dby_vnames   = db_vnames + '_y'
	dbz_vnames   = db_vnames + '_z'
	b_bary_vname = StrJoin( ['mms', 'fgm', 'b', mag_instr, 'bary', fgm_mode, fgm_level], '_' )
	ex_vnames    = e_vnames + '_x'
	ey_vnames    = e_vnames + '_y'
	ez_vnames    = e_vnames + '_z'
	de_vnames    = sc + '_' + StrJoin( [edp_instr, 'de', edp_coords, edp_mode, level], '_' )
	dex_vnames   = de_vnames + '_x'
	dey_vnames   = de_vnames + '_y'
	dez_vnames   = de_vnames + '_z'
	e_bary_vname = StrJoin( ['mms', edp_instr, 'e', edp_coords, 'bary', edp_mode, level], '_' )

;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;B-Field
		IF mag_instr EQ 'fsm' THEN BEGIN
			!MrMMS.dropbox_root = '/nfs/fsm/temp/'
			MrMMS_Load_Data, '', mag_instr, fgm_mode, fgm_level, $
			                 /TEAM_SITE, $
			                 OPTDESC  = fgm_optdesc, $
			                 VARFORMAT = '*_b_'+fgm_coords+'_'+fgm_mode+'*'
		ENDIF ELSE IF mag_instr EQ 'scm' THEN BEGIN
		
		ENDIF ELSE BEGIN
			MrMMS_FGM_Load_Data, '', fgm_mode, $
			                     INSTR     = mag_instr, $
			                     LEVEL     = fgm_level, $
			                     OPTDESC   = fgm_optdesc, $
			                     VARFORMAT = '*_b_*'
		ENDELSE

		;EDP
		MrMMS_Load_Data, '', edp_instr, edp_mode, level, $
		                 OPTDESC   = edp_optdesc, $
		                 VARFORMAT = '*'+edp_optdesc+'_'+[edp_coords, 'par_epar']+'_'+edp_mode+'*'
	ENDIF
	
;-------------------------------------------
; Common Time Basis ////////////////////////
;-------------------------------------------
	;Interpolate B-Field to MMS1
	oB1 = MrVar_Get(bvec_vnames[0])
	oB2 = MrVar_Get(bvec_vnames[1])
	oB3 = MrVar_Get(bvec_vnames[2])
	oB4 = MrVar_Get(bvec_vnames[3])
	oB1 = oB1 -> Copy(bvec_interp_vnames[0], /CACHE)
	oB2 = oB2 -> Interpol(oB1, /CACHE, NAME=bvec_interp_vnames[1])
	oB3 = oB3 -> Interpol(oB1, /CACHE, NAME=bvec_interp_vnames[2])
	oB4 = oB4 -> Interpol(oB1, /CACHE, NAME=bvec_interp_vnames[3])
	
	;Interpolate E-Field to MMS1
	oE1 = MrVar_Get(e_vnames[0])
	oE2 = MrVar_Get(e_vnames[1])
	oE3 = MrVar_Get(e_vnames[2])
	oE4 = MrVar_Get(e_vnames[3])
	oE1 = oE1 -> Copy(e_interp_vnames[0], /CACHE)
	oE2 = oE2 -> Interpol(oE1, /CACHE, NAME=e_interp_vnames[1])
	oE3 = oE4 -> Interpol(oE1, /CACHE, NAME=e_interp_vnames[2])
	oE4 = oE4 -> Interpol(oE1, /CACHE, NAME=e_interp_vnames[3])
	
;-------------------------------------------
; Barycentric Averages /////////////////////
;-------------------------------------------
	;B
	oB_bary = (oB1 + oB2 + oB3 + oB4) / 4.0
	oB_bary -> SetName, b_bary_vname
	oB_bary -> Cache
	oB_bary['CATDESC'] = 'Barycentric average of the vector magnetic field.'
	oB_bary['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oB_bary['LABEL']   = ['Bx', 'By', 'Bz']
	oB_bary['TITLE']   = 'B$\downBC$!C(nT)'
	
	;|B|
	oBmag_bary = oB_bary -> Magnitude(/CACHE, NAME=bmag_bary_vname)
	oBmag_bary['LABEL'] = '|B|'
	oBmag_bary['TITLE'] = '|B|!C(nT)'
	
	;E
	oE_bary = (oE1 + oE2 + oE3 + oE4) / 4.0
	oE_bary -> SetName, e_bary_vname
	oE_bary -> Cache
	oE_bary['COLOR'] = ['Blue', 'Forest Green', 'Red']
	oE_bary['LABEL'] = 'E$\down' + ['X', 'Y', 'Z'] + '$'
	oE_bary['TITLE'] = 'E$\downBC$!C(mV/m^2)'
	
	;Transformation to Field-Aligned Barycentric Coordinates
	oT = MrVar_FAC(oB_bary, !Null, 'CrossX')
	
;-------------------------------------------
; Process Each Spacecraft //////////////////
;-------------------------------------------
	A = 50
	len = 1.0
	FOR i = 0, 3 DO BEGIN
		oB = MrVar_Get(bvec_interp_vnames[i])
		oE = MrVar_Get(e_interp_vnames[i])
		
		;B Filter
		dt     = oB['TIMEVAR'] -> GetSI(RATE=fs)
		fN     = fs / 2.0
		fc     = ffilt / fN
		nTerms = Round(len / dt)
		odB     = oB -> Digital_Filter(fc[0], fc[1], A, nTerms, NAME=db_vnames[i])
		
		;E Filter
		dt     = oE['TIMEVAR'] -> GetSI(RATE=fs)
		fN     = fs / 2.0
		fc     = ffilt / fN
		nTerms = Round(len / dt)
		odE    = oE -> Digital_Filter(fc[0], fc[1], A, nTerms, NAME=de_vnames[i])
		
		;Rotate to FAC
		IF coords EQ 'fac' THEN BEGIN
		
		
		;Keep same coordinates
		ENDIF ELSE BEGIN
			;Split into components
			odB -> Split, odBx, odBy, odBz, /CACHE
			odE -> Split, odEx, odEy, odEz, /CACHE
			
			odBx['TITLE'] = 'dBx!C(nT)'
			odBy['TITLE'] = 'dBy!C(nT)'
			odBz['TITLE'] = 'dBz!C(nT)'
			odBx['COLOR'] = sc_colors[i]
			odBy['COLOR'] = sc_colors[i]
			odBz['COLOR'] = sc_colors[i]
			odBx['LABEL'] = 'mms' + String(i+1, FORMAT='(i1)')
			
			odEx['TITLE'] = 'dEx!c(mV/m)'
			odEy['TITLE'] = 'dEy!c(mV/m)'
			odEz['TITLE'] = 'dEz!c(mV/m)'
			odEx['COLOR'] = sc_colors[i]
			odEy['COLOR'] = sc_colors[i]
			odEz['COLOR'] = sc_colors[i]
		ENDELSE
		
		;Clean up data
		Obj_Destroy, [oB, oE]
	ENDFOR
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	title = String('MMS', mode, len, ffilt, FORMAT='(%"%s %s N=%0.1fs, fc=[%0.1f,%0.1f]Hz")')

	;B_BARY
	oB_bary['PLOT_TITLE'] = title

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [b_bary_vname, $
	                     dbx_vnames[0], dby_vnames[0], dbz_vnames[0], $
	                     dex_vnames[0], dey_vnames[0], dez_vnames[0]], $
	                    /NO_REFRESH, $
	                    XSIZE = 680, $
	                    YSIZE = 700 )
	
	win = MrVar_OPlotTS( b_bary_vname, bmag_bary_vname )
	win = MrVar_OPlotTS( dbx_vnames[0], dbx_vnames[1:3] )
	win = MrVar_OPlotTS( dby_vnames[0], dby_vnames[1:3] )
	win = MrVar_OPlotTS( dbz_vnames[0], dbz_vnames[1:3] )
	win = MrVar_OPlotTS( dex_vnames[0], dex_vnames[1:3] )
	win = MrVar_OPlotTS( dey_vnames[0], dey_vnames[1:3] )
	win = MrVar_OPlotTS( dez_vnames[0], dez_vnames[1:3] )

	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 14]
	win    -> Refresh

;-------------------------------------------
; Save Figure //////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
		ENDIF ELSE IF ~File_Test(output_dir, /DIRECTORY) THEN BEGIN
			MrPrintF, 'LogText', 'Creating directory: "' + output_dir + '".'
			File_MKDir, output_dir
		ENDIF
		
		outdesc = String(optdesc, coords, nfft, nshift, f_pass, FORMAT='(%"%s-de-%s-n%i-s%i-f%i-%i")')
		
		;File name
		fname = StrJoin( [sc, 'edp', mode, level, outdesc], '_' )
		fname = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Finished! ////////////////////////////////
;-------------------------------------------
	RETURN, win
END