; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_Rho
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
;   Generate a plot to provide an overview of reconnection quantities:
;       1. Bx MMS 1-4
;       2. By MMS 1-4
;       3. Bz MMS 1-4
;       4. Ex MMS 1-4
;       5. Ey MMS 1-4
;       6. Ez MMS 1-4
;       7. J = Curl(B)
;       8. Charge density: e0 * Div(E) / q
;
; :Categories:
;   MMS
;
; :Params:
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;
; :Keywords:
;       EPHDESC:    in, optional, type=string, default='ephts04d'
;                   Optional descriptor of the definitive ephemeris datatype to use.
;                       Options are: { 'epht89d' | 'epht89q' | 'ephts04d' | 'defeph' | 'predeph' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
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
;       2018/02/16  -   Written by Matthew Argall
;       2018/10/15  -   Resample data instead of interpolate. Rotate to LMN. Add DiV(E)
;                           visual. - MRA
;-
FUNCTION MrMMS_Plot_4sc_Rho_Div, aoR, aoE, t_ref, $
COMPS=comps, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext
	Compile_Opt idl2
	On_Error, 2
	
	IF N_Elements(cOmps) EQ 0 THEN comps = ['X', 'Y', 'Z']
	sc_colors = ['Black', 'Red', 'Forest Green', 'Blue']
	
	;Which location to display
	it = (aoR[0])['TIMEVAR'] -> Value_Locate(t_ref)

;-------------------------------------------
; Separation ///////////////////////////////
;-------------------------------------------
	
	;Finally, take the mean position
	r1 = Reform( (aoR[0])['DATA',it,*] )
	r2 = Reform( (aoR[1])['DATA',it,*] )
	r3 = Reform( (aoR[2])['DATA',it,*] )
	r4 = Reform( (aoR[3])['DATA',it,*] )

	;First, find the barycenter
	r_bary = (r1 + r2 + r3 + r4) / 4.0
	
	;Next, find separation relative to barycenter
	r1 = r1 - r_bary
	r2 = r2 - r_bary
	r3 = r3 - r_bary
	r4 = r4 - r_bary

;-------------------------------------------
; E-Field //////////////////////////////////
;-------------------------------------------
	E1 = (aoE[0])['DATA',it,*]
	E2 = (aoE[1])['DATA',it,*]
	E3 = (aoE[2])['DATA',it,*]
	E4 = (aoE[3])['DATA',it,*]
	
;-------------------------------------------
; Plot the Positions ///////////////////////
;-------------------------------------------
	win = MrWindow( ASPECT   = 1.0, $
	                LAYOUT   = [3,1], $
	                OXMARGIN = [10,8], $
	                NAME     = '4sc-DivE-Sep', $
	                REFRESH  = 0, $
	                XGAP     = 6, $
	                XSIZE    = 1000, $
	                YGAP     = 0, $
	                YSIZE    = 300 )
	                
	pos = MrLayout( [3,1], $
	                ASPECT   = 1.0, $
	                CHARSIZE = 2.0, $
	                OXMARGIN = [10,8], $
	                WDIMS    = [1000,300], $
	                XGAP     = 6, $
	                YGAP     = 0 )
	
	xrange = [ Min([r1[0], r2[0], r3[0], r4[0]]), Max([r1[0], r2[0], r3[0], r4[0]]) ]
	yrange = [ Min([r1[1], r2[1], r3[1], r4[1]]), Max([r1[1], r2[1], r3[1], r4[1]]) ]
	zrange = [ Min([r1[2], r2[2], r3[2], r4[2]]), Max([r1[2], r2[2], r3[2], r4[2]]) ]
	
	;XY
	range  = [ Min([xrange[0], yrange[0]]), Max([xrange[1], yrange[1]]) ] * 1.15
	zvalue = 3 * ([r1[2], r2[2], r3[2], r4[2]] + abs(zrange[0])) / (zrange[1] - zrange[0]) + 1
	
;-------------------------------------------
; X-Y Plane ////////////////////////////////
;-------------------------------------------
	
	;MMS1
	p1 = MrPlot( r1[0], r1[1], $
	             /CURRENT, $
	             COLOR    = sc_colors[0], $
;	             LAYOUT   = [1,1], $
	             NAME     = StrJoin(comps[[0,1]]) + ' MMS1', $
	             POSITION = pos[*,1], $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[0], $
	             TITLE    = title, $
	             XRANGE   = Reverse(range), $
	             XTITLE   = '$\Delta$' + comps[0] + ' (km)', $
	             YRANGE   = range, $
	             YTITLE   = '$\Delta$' + comps[1] + ' (km)' )
	
	;MMS2
	p2 = MrPlot( r2[0], r2[1], $
	             COLOR    = sc_colors[1], $
	             NAME     = StrJoin(comps[[0,1]]) + ' MMS2', $
	             OVERPLOT = p1, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[1] )
	
	;MMS3
	p3 = MrPlot( r3[0], r3[1], $
	             COLOR    = sc_colors[2], $
	             NAME     = StrJoin(comps[[0,1]]) + ' MMS3', $
	             OVERPLOT = p1, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[2] )
	
	;MMS4
	p4 = MrPlot( r4[0], r4[1], $
	             COLOR    = sc_colors[3], $
	             NAME     = StrJoin(comps[[0,1]]) + ' MMS4', $
	             OVERPLOT = p1, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[3])
	
	;E-vectors
	Ex = [ E1[0], E2[0], E3[0], E4[0] ]
	rx = [ r1[0], r2[0], r3[0], r4[0] ]
	Ey = [ E1[1], E2[1], E3[1], E4[1] ]
	ry = [ r1[1], r2[1], r3[1], r4[1] ]
	v1 = MrVector( Ex, Ey, rx, ry, $
	               NAME          = 'Vector ' + StrJoin(comps[[0,1]]), $
	               OVERPLOT      = p1 )
	
;-------------------------------------------
; X-Z Plane ////////////////////////////////
;-------------------------------------------
	
	;XZ
	range  = [ Min([xrange[0], zrange[0]]), Max([xrange[1], zrange[1]]) ] * 1.15
	zvalue = 3 * ([r1[1], r2[1], r3[1], r4[1]] + abs(yrange[0])) / (yrange[1] - yrange[0]) + 1
	
	p5 = MrPlot( r1[0], r1[2], $
	             /CURRENT, $
	             COLOR    = sc_colors[0], $
;	             LAYOUT   = [1,1], $
	             NAME     = StrJoin(comps[[0,2]]) + ' MMS1', $
	             POSITION = pos[*,0], $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[0], $
	             XRANGE   = Reverse(range), $
	             XTITLE   = '$\Delta$' + comps[0] + ' (km)', $
	             YRANGE   = range, $
	             YTITLE   = '$\Delta$' + comps[2] + ' (km)' )

	p6 = MrPlot( r2[0], r2[2], $
	             COLOR    = sc_colors[1], $
	             NAME     = StrJoin(comps[[0,2]]) + ' MMS2', $
	             OVERPLOT = p5, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[1] )
	
	p7 = MrPlot( r3[0], r3[2], $
	             COLOR    = sc_colors[2], $
	             NAME     = StrJoin(comps[[0,2]]) + ' MMS3', $
	             OVERPLOT = p5, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[2] )
	
	p8 = MrPlot( r4[0], r4[2], $
	             COLOR    = sc_colors[3], $
	             NAME     = StrJoin(comps[[0,2]]) + ' MMS4', $
	             OVERPLOT = p5, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[3] )
	
	;E-vectors
	Ex = [ E1[0], E2[0], E3[0], E4[0] ]
	rx = [ r1[0], r2[0], r3[0], r4[0] ]
	Ez = [ E1[2], E2[2], E3[2], E4[2] ]
	rz = [ r1[2], r2[2], r3[2], r4[2] ]
	v2 = MrVector( Ex, Ez, rx, rz, $
	               NAME          = 'Vector ' + StrJoin(comps[[0,2]]), $
	               OVERPLOT      = p5 )
	
;-------------------------------------------
; Y-Z Plane ////////////////////////////////
;-------------------------------------------
	
	;YZ
	range  = [ Min([yrange[0], zrange[0]]), Max([yrange[1], zrange[1]]) ] * 1.15
	zvalue = 3 * ([r1[0], r2[0], r3[0], r4[0]] + abs(xrange[0])) / (xrange[1] - xrange[0]) + 1
	
	p9 = MrPlot( r1[1], r1[2], $
	             /CURRENT, $
	             COLOR    = sc_colors[0], $
;	             LAYOUT   = [1,1], $
	             NAME     = StrJoin(comps[[1,2]]) + ' MMS1', $
	             POSITION = pos[*,2], $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[0], $
	             XRANGE   = range, $
	             XTITLE   = '$\Delta$' + comps[1] + ' (km)', $
	             YRANGE   = range, $
	             YTITLE   = '$\Delta$' + comps[2] + ' (km)' )

	p10 = MrPlot( r2[1], r2[2], $
	              COLOR    = sc_colors[1], $
	              NAME     = StrJoin(comps[[1,2]]) + ' MMS2', $
	              OVERPLOT = p9, $
	              PSYM     = 'FilledCircle', $
	              SYMSIZE  = zvalue[1], $
	              XRANGE   = range, $
	              YRANGE   = range )
	
	p11 = MrPlot( r3[1], r3[2], $
	              COLOR    = sc_colors[2], $
	              NAME     = StrJoin(comps[[1,2]]) + ' MMS3', $
	              OVERPLOT = p9, $
	              PSYM     = 'FilledCircle', $
	              SYMSIZE  = zvalue[2], $
	              XRANGE   = range, $
	              YRANGE   = range )
	
	p12 = MrPlot( r4[1], r4[2], $
	              COLOR    = sc_colors[3], $
	              NAME     = StrJoin(comps[[1,2]]) + ' MMS4', $
	              OVERPLOT = p9, $
	              PSYM     = 'FilledCircle', $
	              SYMSIZE  = zvalue[3], $
	              XRANGE   = range, $
	              YRANGE   = range )
	;E-vectors
	Ex = [ E1[1], E2[1], E3[1], E4[1] ]
	rx = [ r1[1], r2[1], r3[1], r4[1] ]
	Ey = [ E1[2], E2[2], E3[2], E4[2] ]
	ry = [ r1[2], r2[2], r3[2], r4[2] ]
	v3 = MrVector( Ex, Ey, rx, ry, $
	               NAME          = 'Vector ' + StrJoin(comps[[1,2]]), $
	               OVERPLOT      = p9 )
	
	;Legend
	l1 = MrLegend( ALIGNMENT    = 'NW', $
	               FILL_COLOR   = '', $
	               LABEL        = 'MMS' + ['1', '2', '3', '4'], $
	               LINESTYLE    = 'None', $
	               POSITION     = [1.0, 1.0], $
	               /RELATIVE, $
	               SAMPLE_WIDTH = 0, $
	               SYMBOL       = 'FilledCircle', $
	               /SYM_CENTER, $
	               SYM_COLOR    = sc_colors, $
	               TARGET       = p9, $
	               TEXT_COLOR   = sc_colors )
	
	;Title
	tt     = StrMid(t_ref, 0, 10) + ' ' + StrMid(t_ref, 11)
	title  = 'MMS Tetrahedron ' + tt
	t1     = MrText( 0.5, 0.91, title, $
	                 ALIGNMENT = 0.5, $
	                 CHARSIZE  = 2.0, $
	                 /CURRENT, $
	                 /NORMAL )

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	win -> Refresh
	win -> Refresh ;Needs a second refresh for the vectors to appear

	RETURN, win
END


;+
;
;-
FUNCTION MrMMS_Plot_4sc_Rho, mode, $
COORDS=coords, $
FC = fc, $
FGM_INSTR = fgm_instr, $
EDP_OPTDESC = edp_optdesc, $
EPHDESC=ephdesc, $
LEVEL=level, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
NO_LOAD=no_load, $
T_DIV=t_div, $
TRANGE=trange
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF N_Elements(win) GT 0 THEN Obj_Destroy, win
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	tf_load   = ~Keyword_Set(no_load)
	IF N_Elements(fgm_instr)   EQ 0 THEN fgm_instr   = 'fgm'
	IF N_Elements(edp_optdesc) EQ 0 THEN edp_optdesc = 'dce' 
	IF N_Elements(coords)      EQ 0 THEN coords      = 'gse'
	IF N_Elements(fc)          EQ 0 THEN fc          = 0.0
	IF N_Elements(level)       EQ 0 THEN level       = 'l2'
	IF N_Elements(trange)      GT 0 THEN MrVar_SetTRange, trange
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	
	e0        = MrConstants('epsilon_0')
	mu0       = MrConstants('mu_0')
	q         = MrConstants('q')
	perp      = '!9' + String(120B) + '!X'
	sc_colors = ['Black', 'Red', 'Forest Green', 'Blue']

;-------------------------------------------
; Variable Parameters //////////////////////
;-------------------------------------------
	;FGM
	fgm_mode    = 'srvy' ;mode EQ 'brst' ? mode : 'srvy'
	fgm_coords  = MrIsMember(['dsl', 'dbcs'], coords) ? 'dmpa' : coords
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	CASE coords OF
		'dsl':  fgm_coords = 'dmpa'
		'dbcs': fgm_coords = 'dmpa'
		'gse':  fgm_coords = coords
		'gsm':  fgm_coords = coords
		ELSE:   fgm_coords = 'gse'
	ENDCASE
	CASE fgm_instr OF
		'fsm': fgm_level = 'l3'
		'fgm': fgm_level = 'l2'
		'dfg': fgm_level = 'l2pre'
		'afg': fgm_level = 'l2pre'
		ELSE: Message, 'Invalid FGM instrument: "' + fgm_instr + '".'
	ENDCASE
	
	;EDP
	edp_instr  = 'edp'
	edp_mode   = mode
	IF mode EQ 'srvy' THEN BEGIN
		MrPrintF, 'LogWarn', 'EDP does not have "srvy" data. Using "fast".'
		edp_mode = 'fast'
	ENDIF
	CASE coords OF
		'dmpa': edp_coords = 'dsl'
		'dbcs': edp_coords = 'dsl'
		'gse':  edp_coords = coords
		'gsm':  Message, 'EDP does not have GSM coordinates available.'
		ELSE:   edp_coords = 'gse'
	ENDCASE
	IF edp_optdesc EQ 'hmfe' THEN edp_coords = 'dsl'
	
	;EPH
	CASE coords OF
		'dmpa': eph_coords = 'dsl'
		'dbcs': eph_coords = 'dsl'
		'dsl':  eph_coords = 'dsl'
		'gse':  eph_coords = coords
		'gsm':  eph_coords = coords
		ELSE:   eph_coords = 'gse'
	ENDCASE
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	sc          = 'mms' + ['1', '2', '3', '4']
	b_vnames    = sc + '_' + StrJoin([fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_')
	IF fgm_instr EQ 'fsm' THEN BEGIN
		bvec_vnames = sc + '_' + StrJoin([fgm_instr, 'b', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([fgm_instr, 'b', 'mag',      fgm_mode, fgm_level], '_')
	ENDIF ELSE BEGIN
		bvec_vnames = sc + '_' + StrJoin([fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_')
	ENDELSE
	e_vnames    = sc + '_' + StrJoin([edp_instr, edp_optdesc, edp_coords, edp_mode, level], '_' )
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', eph_coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	bx_vnames = bvec_vnames + '_x'
	by_vnames = bvec_vnames + '_y'
	bz_vnames = bvec_vnames + '_z'
	ex_vnames = e_vnames    + '_x'
	ey_vnames = e_vnames    + '_y'
	ez_vnames = e_vnames    + '_z'
	e_perp_vnames  = sc + '_' + StrJoin([edp_instr, edp_optdesc, 'perp', edp_mode, level], '_' )
	ex_perp_vnames = e_perp_vnames + '_x'
	ey_perp_vnames = e_perp_vnames + '_y'
	ez_perp_vnames = e_perp_vnames + '_z'
	j_vname   = StrJoin( ['mms', fgm_instr,   'j',   fgm_coords, fgm_mode, fgm_level], '_' )
	rho_vname = StrJoin( ['mms', 'dce',     'rho', edp_mode, level], '_' )
		
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;B-Field
		IF fgm_instr EQ 'fsm' THEN BEGIN
			!MrMMS.dropbox_root = '/nfs/fsm/temp/'
			MrMMS_Load_Data, '', fgm_instr, fgm_mode, fgm_level, $
			                 /TEAM_SITE, $
			                 OPTDESC  = fgm_optdesc, $
			                 VARFORMAT = '*_b_'+fgm_coords+'_'+fgm_mode+'*'
		ENDIF ELSE BEGIN
			MrMMS_FGM_Load_Data, 'mms' + ['1', '2', '3', '4'], fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = fgm_level, $
			                     OPTDESC   = fgm_optdesc, $
			                     VARFORMAT = '*_b_'+fgm_coords+'_'+fgm_mode+'_'+fgm_level
		ENDELSE
		
		;E-Field
		MrMMS_Load_Data, '', edp_instr, edp_mode, level, $
		                 OPTDESC   = edp_optdesc, $
		                 VARFORMAT = '*'+edp_optdesc+'_'+edp_coords+'*'
		
		;MEC
		IF mode EQ 'brst' THEN BEGIN
			IF (mrmms_get_filenames('mms1', 'mec', mode, level, OPTDESC=ephdesc))[0] EQ '' $
				THEN mec_mode = 'srvy' $
				ELSE mec_mode = mode
		ENDIF
		
		;Ephemeris
		MrMMS_Load_Data, '', 'mec', mec_mode, level, $
		                 OPTDESC   = ephdesc, $
		                 VARFORMAT = r_vnames
	ENDIF
	
;-------------------------------------------
; Components & Resample ////////////////////
;-------------------------------------------
	;Reference times for resampling
	;   - Resample to B for current density
	;   - Resample to E for charge density
	oB = MrVar_Get(bvec_vnames[0])
;	oE = MrVar_Get(e_vnames[0])
	oTref_B = oB['TIMEVAR']
;	oTref_E = oE['TIMEVAR']
	
	;Axis ranges
	bx_range = [!Values.F_Infinity, -!Values.F_Infinity]
	by_range = [!Values.F_Infinity, -!Values.F_Infinity]
	bz_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ex_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ey_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ez_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ex_perp_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ey_perp_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ez_perp_range = [!Values.F_Infinity, -!Values.F_Infinity]
	
	aoB   = ObjArr(4)
	aoE   = ObjArr(4)
	aoEperp = ObjArr(4)
	aoPos = ObjArr(4)
	aoR_E = ObjArr(4)
	FOR i = 0, 3 DO BEGIN
		;Grab data
		oB   = MrVar_Get(bvec_vnames[i])
		oE   = MrVar_Get(e_vnames[i])
		oPos = MrVar_Get(r_vnames[i])
	
	;-------------------------------------------
	; Resample /////////////////////////////////
	;-------------------------------------------
		;Use MVA coordinates for B so that J is in the new system.
		oB   = MrVar_Resample(oB,   oTref_B)
		oPos = MrVar_Resample(oPos, oTref_B)
		
		;Use orignal coordinates for E because rho is a scalar
		oE   = MrVar_Resample(oE,   oTref_B)
		oR_E = oPos
	
	;-------------------------------------------
	; Perpendicular Electric Field /////////////
	;-------------------------------------------
;		oTx = MrVar_FAC(oB)
;		oE_perp      = oTx ## oE
;		oE_perp[*,2] = 0.0
;		oE_perp      = (oTx -> Transpose()) ## oE_perp
;		oE_perp -> SetName, e_perp_vnames[i]
;		oE_perp -> Cache
	
	;-------------------------------------------
	; Rotate ///////////////////////////////////
	;-------------------------------------------
		IF coords EQ 'mva' THEN BEGIN
			aoB[i]   = MrVar_xForm(oB)
			aoE[i]   = MrVar_xForm(oE)
;			oE_perp  = MrVar_xForm(oE_perp)
			aoPos[i] = MrVar_xForm(oPos)
			aoR_E[i] = aoPos[i]
			comps    = ['N', 'M', 'L']
		ENDIF ELSE BEGIN
			aoB[i]   = oB   ;oB -> Copy()
			aoE[i]   = oE   ;oE -> Copy()
			aoPos[i] = oPos ;oPos -> Copy()
			aoR_E[i] = oPos
			comps    = ['X', 'Y', 'Z']
		ENDELSE
	
	;-------------------------------------------
	; Split Into Components ////////////////////
	;-------------------------------------------
		IF aoB[i] -> HasAttr('MIN_VALUE') THEN aoB[i] -> RemoveAttr, ['MIN_VALUE', 'MAX_VALUE']
		aoB[i] -> Split, oBx, oBy, oBz, /CACHE
		aoE[i] -> Split, oEx, oEy, oEz, /CACHE
;		oE_perp -> Split, oEx_perp, oEy_perp, oEz_perp, /CACHE
		
		;Bx
		oBx -> SetName, bx_vnames[i]
		oBx['COLOR'] = sc_colors[i]
		oBx['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
		oBx['TITLE'] = 'B$\down' + comps[0] + '$!C(nT)'
		bx_range[0] <= oBx.min
		bx_range[1] >= oBx.max
		
		;By
		oBy -> SetName, by_vnames[i]
		oBy['COLOR'] = sc_colors[i]
;		oBy['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
		oBy['TITLE'] = 'B$\down' + comps[1] + '$!C(nT)'
		by_range[0] <= oBy.min
		by_range[1] >= oBy.max
		
		;Bz
		oBz -> SetName, bz_vnames[i]
		oBz['COLOR'] = sc_colors[i]
;		oBz['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
		oBz['TITLE'] = 'B$\down' + comps[2] + '$!C(nT)'
		bz_range[0] <= oBz.min
		bz_range[1] >= oBz.max
		
		;Ex
		oEx -> SetName, ex_vnames[i]
		oEx['COLOR'] = sc_colors[i]
;		oEx['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
		oEx['TITLE'] = 'E$\down' + comps[0] + '$!C(mV/m)'
		ex_range[0] <= oEx.min
		ex_range[1] >= oEx.max
		
		;Ey
		oEy -> SetName, ey_vnames[i]
		oEy['COLOR'] = sc_colors[i]
;		oEy['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
		oEy['TITLE'] = 'E$\down' + comps[1] + '$!C(mV/m)'
		ey_range[0] <= oEy.min
		ey_range[1] >= oEy.max
		
		;Ez
		oEz -> SetName, ez_vnames[i]
		oEz['COLOR'] = sc_colors[i]
;		oEz['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
		oEz['TITLE'] = 'E$\down' + comps[2] + '$!C(mV/m)'
		ez_range[0] <= oEz.min
		ez_range[1] >= oEz.max
		
		;Ex-perp
;		oEx_perp -> SetName, ex_perp_vnames[i]
;		oEx_perp['COLOR'] = sc_colors[i]
;		oEx_perp['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
;		oEx_perp['TITLE'] = 'E$\down' + perp + comps[0] + '$!C(mV/m)'
;		ex_perp_range[0] <= oEx_perp.min
;		ex_perp_range[1] >= oEx_perp.max
		
		;Ey-perp
;		oEy_perp -> SetName, ey_perp_vnames[i]
;		oEy_perp['COLOR'] = sc_colors[i]
;		oEy_perp['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
;		oEy_perp['TITLE'] = 'E$\down' + perp + comps[1] + '$!C(mV/m)'
;		ey_perp_range[0] <= oEy_perp.min
;		ey_perp_range[1] >= oEy_perp.max
		
		;Ez-perp
;		oEz_perp -> SetName, ez_perp_vnames[i]
;		oEz_perp['COLOR'] = sc_colors[i]
;		oEz_perp['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
;		oEz_perp['TITLE'] = 'E$\down' + perp + comps[2] + '$!C(mV/m)'
;		ez_perp_range[0] <= oEz_perp.min
;		ez_perp_range[1] >= oEz_perp.max
	ENDFOR
	
;-------------------------------------------
; Charge Density ///////////////////////////
;-------------------------------------------
	
	;Reciprocal vectors for taking derivatives
	oRecipVec = MrVar_RecipVec(aoR_E[0], aoR_E[1], aoR_E[2], aoR_E[3])
	
	;Charge Density
	;   - 1e10 converts to uC/cm^2
	oRho = (e0 / q * 1e-12) * oRecipVec -> Divergence(aoE[0], aoE[1], aoE[2], aoE[3])
	oRho -> SetName, rho_vname
	oRho -> Cache
	oRho['CATDESC'] = 'Free charge density, computed from the divergence of the electric field.'
	oRho['TITLE']   = '$\rho$/e!C(cm$\up-3$)'
	oRho['UNITS']   = 'cm^-3'
	
	;Free memory
	Obj_Destroy, oRecipVec
	
;-------------------------------------------
; Curlometer ///////////////////////////////
;-------------------------------------------
	
	;Reciprocal vectors for taking derivatives
	oRecipVec = MrVar_RecipVec(aoPos[0], aoPos[1], aoPos[2], aoPos[3])
	
	;Curlometer
	oJ  = (1e-6/mu0) * oRecipVec -> Curl(aoB[0], aoB[1], aoB[2], aoB[3])
	oJ -> SetName, j_vname
	oJ -> Cache
	oJ['CATDESC'] = 'Current density derived from the curlometer technique using MEC & ' + fgm_instr + '.'
	oJ['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oJ['LABEL']   = 'J$\down' + comps + '$'
	oJ['TITLE']   = 'J!C($\mu$A/m^2)'
	oJ['UNITS']   = 'nA/m^2'
	
	;Free Memory
	Obj_Destroy, oRecipVec
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	;B
	oBx = MrVar_Get(bx_vnames[0])
	oBy = MrVar_Get(by_vnames[0])
	oBz = MrVar_Get(bz_vnames[0])
	oBx['PLOT_TITLE'] = 'MMS ' + StrUpCase(mode)
	oBx['AXIS_RANGE'] = bx_range
	oBy['AXIS_RANGE'] = by_range
	oBz['AXIS_RANGE'] = bz_range
	
	;E
	oEx = MrVar_Get(ex_vnames[0])
	oEy = MrVar_Get(ey_vnames[0])
	oEz = MrVar_Get(ez_vnames[0])
	oEx['AXIS_RANGE'] = ex_range
	oEy['AXIS_RANGE'] = ey_range
	oEz['AXIS_RANGE'] = ez_range
	
	;E-perp
;	oEx_perp = MrVar_Get(ex_perp_vnames[0])
;	oEy_perp = MrVar_Get(ey_perp_vnames[0])
;	oEz_perp = MrVar_Get(ez_perp_vnames[0])
;	oEx_perp['AXIS_RANGE'] = ex_perp_range
;	oEy_perp['AXIS_RANGE'] = ey_perp_range
;	oEz_perp['AXIS_RANGE'] = ez_perp_range

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot MMS1 data
	win = MrVar_PlotTS( [ bx_vnames[0], by_vnames[0], bz_vnames[0], $
	                      ex_vnames[0], ey_vnames[0], ez_vnames[0], $
	                      j_vname, rho_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	
	;Overplot MMS 2-4 data
	win = MrVar_OPlotTS( bx_vnames[0], bx_vnames[1:3] )
	win = MrVar_OPlotTS( by_vnames[0], by_vnames[1:3] )
	win = MrVar_OPlotTS( bz_vnames[0], bz_vnames[1:3] )
	win = MrVar_OPlotTS( ex_vnames[0], ex_vnames[1:3] )
	win = MrVar_OPlotTS( ey_vnames[0], ey_vnames[1:3] )
	win = MrVar_OPlotTS( ez_vnames[0], ez_vnames[1:3] )
	
	;Add horizontal lines
	l1 = MrPlotS(win[bx_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Bx', TARGET=win[bx_vnames[0]])
	l2 = MrPlotS(win[by_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: By', TARGET=win[by_vnames[0]])
	l3 = MrPlotS(win[bz_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Bz', TARGET=win[bz_vnames[0]])
	l4 = MrPlotS(win[ex_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Ex', TARGET=win[ex_vnames[0]])
	l5 = MrPlotS(win[ey_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Ey', TARGET=win[ey_vnames[0]])
	l6 = MrPlotS(win[ez_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Ez', TARGET=win[ez_vnames[0]])
	l7 = MrPlotS(win[j_vname].xrange,   [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: J',   TARGET=win[j_vname])
	l8 = MrPlotS(win[rho_vname].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Rho', TARGET=win[rho_vname])
	
	;Pretty-up the window
	win.name = '4sc-Rho'
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 11]
	
;-------------------------------------------
; Position /////////////////////////////////
;-------------------------------------------
	;Default reference time
	IF N_Elements(t_div) EQ 0 THEN BEGIN
		!Null = Min(oRho['DATA'], iMin)
		t_div = oRho['TIME', iMin]
	ENDIF ELSE BEGIN
		iMin = oRho['TIMEVAR'] -> Value_Locate(t_div)
	ENDELSE
	
	;Plot positions and E-field vectors
	w2 = MrMMS_Plot_4sc_Rho_Div( aoR_E, aoE, t_div, $
	                             COMPS      = comps, $
	                             OUTPUT_DIR = output_dir, $
	                             OUTPUT_EXT = output_ext )
	
	;Draw line on divergence graph
	t_ssm = oRho['TIME', iMin, 'SSM']
	l9    = MrPlotS( Replicate(t_ssm, 2), win[rho_vname].yrange, COLOR='Magenta', LINESTYLE='--', NAME='Line: Div', TARGET=win[rho_vname])
	
	;Free memory
	Obj_Destroy, aoR_E

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
		
		;File name
		fname1   = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-rho'], '_' )
		fname2   = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-divE-sep'], '_' )
		fname1   = FilePath( fname1, ROOT_DIR=output_dir )
		fname2   = FilePath( fname2, ROOT_DIR=output_dir )
		
		;Save the figure
		fout1 = MrVar_PlotTS_Save( win, fname1, output_ext )
		fout2 = MrVar_PlotTS_Save( w2, fname2, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	win    -> Refresh
	RETURN, win
END