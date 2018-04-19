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
;       7. Charge density: e0 * Div(E)
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
;-
FUNCTION MrMMS_Plot_4sc_Rho, mode, $
FC = fc, $
FGM_INSTR = fgm_instr, $
EPHDESC=ephdesc, $
LEVEL=level, $
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
	
	tf_load   = ~Keyword_Set(no_load)
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = mode EQ 'brst' ? 'fsm' : 'fgm'
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fc)        EQ 0 THEN fc        = 0.0
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	
	e0        = MrConstants('epsilon_0')
	mu0       = MrConstants('mu_0')
	q         = MrConstants('q')
	sc_colors = ['Black', 'Red', 'Green', 'Blue']

;-------------------------------------------
; Variable Parameters //////////////////////
;-------------------------------------------
	;FGM
	fgm_mode    = mode EQ 'brst' ? mode : 'srvy'
	fgm_coords  = MrIsMember(['dsl', 'dbcs'], coords) ? 'dmpa' : coords
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	CASE fgm_instr OF
		'fsm': fgm_level = 'l3'
		'fgm': fgm_level = 'l2'
		'dfg': fgm_level = 'l2pre'
		'afg': fgm_level = 'l2pre'
		ELSE: Message, 'Invalid FGM instrument: "' + fgm_instr + '".'
	ENDCASE
	
	;EDP
	edp_instr  = 'edp'
	edp_coords = MrIsMember(['dmpa', 'dbcs'], coords) ? 'dsl' : coords
	edp_mode   = mode
	IF mode EQ 'srvy' THEN BEGIN
		MrPrintF, 'LogWarn', 'EDP does not have "srvy" data. Using "fast".'
		edp_mode = 'fast'
	ENDIF
	
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
	e_vnames    = sc + '_' + StrJoin([edp_instr, 'dce',  edp_coords, edp_mode, level], '_' )
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	bx_vnames = bvec_vnames + '_x'
	by_vnames = bvec_vnames + '_y'
	bz_vnames = bvec_vnames + '_z'
	ex_vnames = e_vnames    + '_x'
	ey_vnames = e_vnames    + '_y'
	ez_vnames = e_vnames    + '_z'
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
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = '*dce_'+edp_coords+'*'
		
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
; Split Vectors into Components ////////////
;-------------------------------------------
	bx_range = [!Values.F_Infinity, -!Values.F_Infinity]
	by_range = [!Values.F_Infinity, -!Values.F_Infinity]
	bz_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ex_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ey_range = [!Values.F_Infinity, -!Values.F_Infinity]
	ez_range = [!Values.F_Infinity, -!Values.F_Infinity]
	
	FOR i = 0, 3 DO BEGIN
		oB = MrVar_Get(bvec_vnames[i])
		IF oB -> HasAttr('MIN_VALUE') THEN oB -> RemoveAttr, ['MIN_VALUE', 'MAX_VALUE']
		oB -> Split, oBx, oBy, oBz, /CACHE
		oBx['COLOR'] = sc_colors[i]
		oBy['COLOR'] = sc_colors[i]
		oBz['COLOR'] = sc_colors[i]
		oBx['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
;		oBy['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
;		oBz['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
		oBx['TITLE'] = 'Bx!C(nT)'
		oBy['TITLE'] = 'By!C(nT)'
		oBz['TITLE'] = 'Bz!C(nT)'
		bx_range[0] <= oBx.min
		by_range[0] <= oBy.min
		bz_range[0] <= oBz.min
		bx_range[1] >= oBx.max
		by_range[1] >= oBy.max
		bz_range[1] >= oBz.max
		
		oE = MrVar_Get(e_vnames[i])
		oE -> Split, oEx, oEy, oEz, /CACHE
		oEx['COLOR'] = sc_colors[i]
		oEy['COLOR'] = sc_colors[i]
		oEz['COLOR'] = sc_colors[i]
;		oEx['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
;		oEy['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
;		oEz['LABEL'] = 'MMS' + String(i+1, FORMAT='(i1)')
		oEx['TITLE'] = 'Ex!C(mV/m)'
		oEy['TITLE'] = 'Ey!C(mV/m)'
		oEz['TITLE'] = 'Ez!C(mV/m)'
		ex_range[0] <= oEx.min
		ey_range[0] <= oEy.min
		ez_range[0] <= oEz.min
		ex_range[1] >= oEx.max
		ey_range[1] >= oEy.max
		ez_range[1] >= oEz.max
	ENDFOR
	
;-------------------------------------------
; Charge Density ///////////////////////////
;-------------------------------------------
	;Get E-field
	oE1 = MrVar_Get(e_vnames[0])
	oE2 = MrVar_Get(e_vnames[1])
	oE3 = MrVar_Get(e_vnames[2])
	oE4 = MrVar_Get(e_vnames[3])
	
	;Get Position
	oR1 = MrVar_Get(r_vnames[0])
	oR2 = MrVar_Get(r_vnames[1])
	oR3 = MrVar_Get(r_vnames[2])
	oR4 = MrVar_Get(r_vnames[3])
	
	;Interpolate to MMS1
	oE2 = oE2 -> Interpol(oE1)
	oE3 = oE3 -> Interpol(oE1)
	oE4 = oE4 -> Interpol(oE1)
	oR1 = oR1 -> Interpol(oE1)
	oR2 = oR2 -> Interpol(oE1)
	oR3 = oR3 -> Interpol(oE1)
	oR4 = oR4 -> Interpol(oE1)
	
	;Reciprocal vectors for taking derivatives
	oRecipVec = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	
	;Charge Density
	;   - 1e10 converts to uC/cm^2
	oRho = (e0 / q * 1e-12) * oRecipVec -> Divergence( oE1, oE2, oE3, oE4 )
	oRho -> SetName, rho_vname
	oRho -> Cache
	oRho['CATDESC'] = 'Free charge density, computed from the divergence of the electric field.'
	oRho['TITLE']   = '$\Delta$#!C(#/cm^3)'
	oRho['UNITS']   = '#/cm^3'
	
;-------------------------------------------
; Curlometer ///////////////////////////////
;-------------------------------------------
	;Get B-Field
	oB1 = MrVar_Get(bvec_vnames[0])
	oB2 = MrVar_Get(bvec_vnames[1])
	oB3 = MrVar_Get(bvec_vnames[2])
	oB4 = MrVar_Get(bvec_vnames[3])
	
	;Inerpolate to E on MMS1
	oB1 = oB1 -> Interpol(oE1)
	oB2 = oB2 -> Interpol(oE1)
	oB3 = oB3 -> Interpol(oE1)
	oB4 = oB4 -> Interpol(oE1)
	
	;Curlometer
	oJ  = (1e-6/mu0) * oRecipVec -> Curl(oB1, oB2, oB3, oB4)
	oJ -> SetName, j_vname
	oJ -> Cache
	oJ['CATDESC'] = 'Current density derived from the curlometer technique using MEC & ' + fgm_instr + '.'
	oJ['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oJ['LABEL']   = 'J$\down' + ['X', 'Y', 'Z'] + '$'
	oJ['TITLE']   = 'J!C($\mu$A/m^2)'
	oJ['UNITS']   = 'nA/m^2'
	
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
	l1 = MrPlotS(win[bx_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Bx',  TARGET=win[bx_vnames[0]])
	l2 = MrPlotS(win[by_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: By',  TARGET=win[by_vnames[0]])
	l3 = MrPlotS(win[bz_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Bz',  TARGET=win[bz_vnames[0]])
	l4 = MrPlotS(win[ex_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Ex',  TARGET=win[ex_vnames[0]])
	l5 = MrPlotS(win[ey_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Ey',  TARGET=win[ey_vnames[0]])
	l6 = MrPlotS(win[ez_vnames[0]].xrange, [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Ez',  TARGET=win[ez_vnames[0]])
	l7 = MrPlotS(win[j_vname].xrange,      [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: J',   TARGET=win[j_vname])
	l8 = MrPlotS(win[rho_vname].xrange,    [0,0], COLOR='Magenta', LINESTYLE='--', NAME='Line: Rho', TARGET=win[rho_vname])
		
	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 11]
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
		
		;File name
		fname   = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-rho'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END