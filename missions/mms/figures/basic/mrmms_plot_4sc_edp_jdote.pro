; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FGM_4sc
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
;   Generate a plot to provide an overview of reconnection quantities:
;       1. Ex MMS1-4
;       2. Ey MMS1-4
;       3. Ez MMS1-4
;       4. E Parallel
;       5. Vsc
;       6. 1/e0 Div(E)
;       7. Curl(E)
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
;       2017/01/05  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_EDP_JdotE, mode, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
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
	
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fsm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	
	mu0       = MrConstants('mu_0')
	q         = MrConstants('q')
	sc_colors = ['Black', 'Red', 'Green', 'Blue']
	
;-------------------------------------------
; Data Parameters //////////////////////////
;-------------------------------------------
	;FGM
	CASE fgm_instr OF
		'afg': fgm_level = 'l2pre'
		'dfg': fgm_level = 'l2pre'
		'fgm': fgm_level = 'l2'
		'fsm': fgm_level = 'l3'
		ELSE: Message, 'Invalid FGM instrument: "' + fgm_instr + '".'
	ENDCASE
	fgm_coords  = MrIsMember(['dsl', 'dbcs'], coords) ? 'dmpa' : coords
	IF coords EQ 'fac' THEN fgm_coords = 'gse'
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	tf_team     = fgm_instr EQ 'fsm'
	
	;EDP
	edp_instr   = 'edp'
	edp_coords  = MrIsMember(['dmpa', 'dbcs'], coords) ? 'dsl' : coords
	IF coords EQ 'fac' THEN edp_coords = 'gse'
	edp_mode    = mode EQ 'brst' ? mode : 'fast'
	edp_optdesc = 'dce'
	
	;MEC
	mec_mode   = 'srvy'
	mec_coords = coords
	IF coords EQ 'fac' THEN mec_coords = 'gse'
	
;	IF N_Elements(coords) EQ 0 THEN BEGIN
;		CASE level OF
;			'ql': coords = 'dsl'
;			ELSE: coords = 'gse'
;		ENDCASE
;	ENDIF
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------

	;Source names
	sc = 'mms' + ['1', '2', '3', '4']
	bvec_vnames = sc + '_' + StrJoin( [fgm_instr, 'b',   fgm_coords, mode, fgm_level], '_' )
	bmag_vnames = sc + '_' + StrJoin( [fgm_instr, 'b',   'mag',      mode, fgm_level], '_' )
	e_vnames    = sc + '_' + StrJoin( [edp_instr, 'dce', edp_coords, mode, level], '_' )
	epe_vnames  = sc + '_' + StrJoin( [edp_instr, 'dce',   'par',  'epar', mode, level], '_' )
	
	;MEC
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames = sc + '_' + StrJoin( ['mec', 'r', mec_coords], '_' ) $
		ELSE r_vnames = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	bout_vnames       = bvec_vnames + '_out'
	eout_vnames       = e_vnames + '_out'
	rout_vnames       = r_vnames + '_out'
	ex_vnames         = e_vnames + '_x'
	ey_vnames         = e_vnames + '_y'
	ez_vnames         = e_vnames + '_z'
	b_bary_vname      = StrJoin( ['mms1234', fgm_instr, 'bbarycenter', coords, mode, level], '_' )
	e_par_vnames      = sc + '_' + StrJoin( [edp_instr, 'dce', 'par', mode, level], '_' )
	e_bary_vname      = StrJoin( ['mms1234', edp_instr, 'ebarycenter', coords, mode, level], '_' )
	e_par_bary_vname  = StrJoin( ['mms1234', edp_instr, 'ebarycenter', 'par', mode, level], '_' )
	j_vname           = StrJoin( ['mms1234', fgm_instr, 'j',     coords, mode, level], '_' )
	jdote_vname       = StrJoin( [fgm_instr, 'jdote', coords, mode, level], '_' )
	jdote_par_vname   = StrJoin( [fgm_instr, 'jdote', 'par',  mode, level], '_' )
	jdote_perp_vname  = StrJoin( [fgm_instr, 'jdote', 'perp', mode, level], '_' )
	jdotepar_bary_vname = StrJoin( [fgm_instr, 'jdotepar', 'bary', mode, level], '_' )
	exb_bary_vname    = StrJoin( ['mms1234', edp_instr, 'exbbary', coords, mode, level], '_' )
	exb_perp1_vname   = StrJoin( ['mms1234', edp_instr, 'exbbary', 'perp1', mode, level], '_' )
	exb_perp2_vname   = StrJoin( ['mms1234', edp_instr, 'exbbary', 'perp2', mode, level], '_' )
	v_vname           = StrJoin( [fgm_instr, 'v', coords, mode, level], '_' )
	v_perp1_vname     = StrJoin( [fgm_instr, 'v', 'perp1', mode, level], '_' )
	v_perp2_vname     = StrJoin( [fgm_instr, 'v', 'perp2', mode, level], '_' )
		
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;B-Field
		MrMMS_Load_Data, '', fgm_instr, mode, fgm_level, $
		                 OPTDESC   = fgm_optdesc, $
		                 TEAM_SITE = tf_team, $
		                 VARFORMAT = ['*b_mag*', '*b_'+fgm_coords+'_'+mode+'*']
		
		;E-Field
		MrMMS_Load_Data, '', edp_instr, edp_mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = ['*dce_'+edp_coords+'*', '*dce_par*']
		
		;Ephemeris
		IF mec_mode EQ 'brst' THEN BEGIN
			IF (mrmms_get_filenames(sc, 'mec', mode, level, OPTDESC=ephdesc))[0] EQ '' $
				THEN mec_mode = 'srvy' $
				ELSE mec_mode = mode
		ENDIF
		MrMMS_Load_Data, '', 'mec', mec_mode, level, $
		                 OPTDESC   = ephdesc, $
		                 VARFORMAT = r_vnames
	ENDIF
	
;-------------------------------------------
; Remove Limits ////////////////////////////
;-------------------------------------------
	oB1 = MrVar_Get(bvec_vnames[0])
	IF oB1 -> HasAttr('MIN_VALUE') && oB1['MIN_VALUE'] EQ 0.0 THEN BEGIN
		oB2 = MrVar_Get(bvec_vnames[1])
		oB3 = MrVar_Get(bvec_vnames[2])
		oB4 = MrVar_Get(bvec_vnames[3])
		oB1 -> RemoveAttr, ['MIN_VALUE', 'VALIDMIN']
		oB2 -> RemoveAttr, ['MIN_VALUE', 'VALIDMIN']
		oB3 -> RemoveAttr, ['MIN_VALUE', 'VALIDMIN']
		oB4 -> RemoveAttr, ['MIN_VALUE', 'VALIDMIN']
	ENDIF
	
;-------------------------------------------
; Interpolate Everything to MMS1 ///////////
;-------------------------------------------
	oE1 = MrVar_Get(e_vnames[0])
	oE2 = MrVar_Get(e_vnames[1])
	oE3 = MrVar_Get(e_vnames[2])
	oE4 = MrVar_Get(e_vnames[3])
	oE1 = oE1 -> Copy(eout_vnames[0], /CACHE)
	oE2 = oE2 -> Interpol(oE1, /CACHE, NAME=eout_vnames[1])
	oE3 = oE3 -> Interpol(oE1, /CACHE, NAME=eout_vnames[2])
	oE4 = oE4 -> Interpol(oE1, /CACHE, NAME=eout_vnames[3])
	
	oEpar1 = MrVar_Get(epe_vnames[0])
	oEpar2 = MrVar_Get(epe_vnames[1])
	oEpar3 = MrVar_Get(epe_vnames[2])
	oEpar4 = MrVar_Get(epe_vnames[3])
	oEpar1 = oEpar1 -> Interpol(oE1)
	oEpar2 = oEpar2 -> Interpol(oE2)
	oEpar3 = oEpar3 -> Interpol(oE3)
	oEpar4 = oEpar4 -> Interpol(oE4)
	
	oB1 = MrVar_Get(bvec_vnames[0])
	oB2 = MrVar_Get(bvec_vnames[1])
	oB3 = MrVar_Get(bvec_vnames[2])
	oB4 = MrVar_Get(bvec_vnames[3])
	oB1 = oB1 -> Interpol(oE1, /CACHE, NAME=bout_vnames[0])
	oB2 = oB2 -> Interpol(oE1, /CACHE, NAME=bout_vnames[1])
	oB3 = oB3 -> Interpol(oE1, /CACHE, NAME=bout_vnames[2])
	oB4 = oB4 -> Interpol(oE1, /CACHE, NAME=bout_vnames[3])
	
	oBmag1 = MrVar_Get(bmag_vnames[0])
	oBmag2 = MrVar_Get(bmag_vnames[1])
	oBmag3 = MrVar_Get(bmag_vnames[2])
	oBmag4 = MrVar_Get(bmag_vnames[3])
	oBmag1 = oBmag1 -> Interpol(oE1)
	oBmag2 = oBmag2 -> Interpol(oE1)
	oBmag3 = oBmag3 -> Interpol(oE1)
	oBmag4 = oBmag4 -> Interpol(oE1)
	
	oR1 = MrVar_Get(r_vnames[0])
	oR2 = MrVar_Get(r_vnames[1])
	oR3 = MrVar_Get(r_vnames[2])
	oR4 = MrVar_Get(r_vnames[3])
	oR1 = oR1 -> Interpol(oE1, /CACHE, NAME=rout_vnames[0])
	oR2 = oR2 -> Interpol(oE1, /CACHE, NAME=rout_vnames[1])
	oR3 = oR3 -> Interpol(oE1, /CACHE, NAME=rout_vnames[2])
	oR4 = oR4 -> Interpol(oE1, /CACHE, NAME=rout_vnames[3])
	
;-------------------------------------------
; Rotate to FAC ////////////////////////////
;-------------------------------------------
	perp = '!D!9' + String(120B) + '!X'
	IF coords EQ 'fac' THEN BEGIN
		labels = [perp+'1!N', perp+'2!N', '!D||!X']
		FOR i = 0, N_Elements(e_vnames) - 1 DO BEGIN
			;Get the data
			oE    = MrVar_Get(eout_vnames[i])
			oB    = MrVar_Get(bout_vnames[i])
		
			;Rotate E to FAC
			oT = MrVar_FAC(oB, '', 'CROSSX')
			oE = oT ## oE
			oE -> Split, oEx, oEy, oEz
			
			;Ex
			oEx -> SetName, ex_vnames[i]
			oEx -> Cache
			oEx['CATDESC'] = 'Electric field perp1-component in field-aligned coordinates.'
			oEx['COLOR']   = sc_colors[i]
			oEx['COORDS']  = 'FAC: PAR=B, PERP2=(PAR x [1,0,0]), PERP1=(PERP2 x PAR)'
			oEx['TITLE']   = 'E!D!9' + String(120B) + '!X1!N!C(mV/m)'
			oEx['UNITS']   = 'mV/m'
		
			;Ey
			oEy -> SetName, ey_vnames[i]
			oEy -> Cache
			oEy['CATDESC'] = 'Electric field perp2-component in field-aligned coordinates.'
			oEy['COLOR']   = sc_colors[i]
			oEy['COORDS']  = 'FAC: PAR=B, PERP2=(PAR x [1,0,0]), PERP1=(PERP2 x PAR)'
			oEy['TITLE']   = 'E!D!9' + String(120B) + '!X2!C(mV/m)'
			oEy['UNITS']   = 'mV/m'
		
			;Ez
			oEz -> SetName, ez_vnames[i]
			oEz -> Cache
			oEz['CATDESC'] = 'Electric field parallel-component in field-aligned coordinates.'
			oEz['COLOR']   = sc_colors[i]
			oEz['TITLE']   = 'E$\down||$!C(mV/m)'
			oEz['UNITS']   = 'mV/m'
		ENDFOR
	
	;Split into components without rotating
	ENDIF ELSE BEGIN
		labels = '$\down' + ['X', 'Y', 'Z'] + '$'
		FOR i = 0, N_Elements(e_vnames) - 1 DO BEGIN
			oE  = MrVar_Get(e_vnames[i])
			oE -> Split, oEx, oEy, oEz, /CACHE

			;Ex
			oEx['TITLE'] = 'Ex!C(mV/m)'
	
			;Ey
			oEy['TITLE'] = 'Ey!C(mV/m)'
	
			;Ez
			oEz['TITLE'] = 'Ez!C(mV/m)'
		ENDFOR
	ENDELSE

;-------------------------------------------
; Barycentric Fields ///////////////////////
;-------------------------------------------
	;B
	oB_bary     = (oB1 + oB2 + oB3 + oB4) / 4.0
	oB_bary_hat = oB_bary -> Normalize()
	oT_bary     = MrVar_FAC(oB_bary, 'CROSSX')
	
	;|B|
	oBmag_bary = (oBmag1 + oBmag2 + oBmag3 + oBmag4) / 4.0
	
	;R
	oR_bary = (oR1 + oR2 + oR3 + oR4) / 4.0
	
	;E-Parallel
	oE_par_bary = (oEpar1[*,1] + oEpar2[*,1] + oEpar3[*,1] + oEpar4[*,1]) / 4.0
	oE_par_bary -> SetName, e_par_bary_vname
	oE_par_bary -> Cache
	
	;E
	oE_bary = (oE1 + oE2 + oE3 + oE4) / 4.0
	oE_bary -> SetName, e_bary_vname
	oE_bary -> Cache
	
	;E-Parallel
	oE_bary_par = oE_bary -> Dot(oB_bary_hat)
	oE_bary_par['CATDESC'] = ['Parallel component of the barycentric average of the Electric field.']
	oE_bary_par['TITLE']   = 'E$\downBC,||$!C(mV/m)'
	oE_barY_par['UNITS']   = 'mV/m'
	
	;E-Perpendicular
	oE_bary_par_vec = oE_bary_par * oB_bary_hat
	oE_bary_perp    = oE_bary - oE_bary_par_vec
	oE_bary_perp['CATDESC'] = ['Perpendicular components of the barycentric average of the Electric field.']
	oE_bary_perp['TITLE']   = 'E$\downBC,$' + perp + '!N!C(mV/m)'
	oE_bary_perp['UNITS']   = 'mV/m'

;-------------------------------------------
; Current Density //////////////////////////
;-------------------------------------------
	oRecip = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	oJ     = (1e-3/mu0) * oRecip -> Curl(oB1, oB2, oB3, oB4)
	oJ_fac = oT_bary ## oJ
	oJ_fac -> SetName, j_vname
	oJ_fac -> Cache
	oJ_fac['CATDESC'] = 'Current density derived from the curlometer technique using MEC & ' + fgm_instr + '.'
	oJ_fac['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oJ_fac['LABEL']   = 'J' + labels
	oJ_fac['TITLE']   = 'J!C(nA/m^2)'
	oJ_fac['UNITS']   = 'nA/m^2'
	
	;J Par
	oJ_par = oJ -> Dot(oB_bary_hat)
	oJ_par['CATDESC'] = 'Parallel current density derived from the curlometer technique using MEC & ' + fgm_instr + '.'
	oJ_par['COLOR']   = 'Red'
	oJ_par['LABEL']   = 'J' + labels[2]
	oJ_par['TITLE']   = 'J' + labels[2] + '!C(nA/m^2)'
	oJ_par['UNITS']   = 'nA/m^2'
	
	;J Perp
	oJ_par_vec = oJ_par * oB_bary_hat
	oJ_perp    = oJ - oJ_par_vec
	oJ_perp['CATDESC'] = 'Parallel current density derived from the curlometer technique using MEC & ' + fgm_instr + '.'
	oJ_perp['COLOR']   = 'Red'
	oJ_perp['LABEL']   = 'J' + perp + '!N'
	oJ_perp['TITLE']   = 'J' + perp + '!N!C(nA/m^2)'
	oJ_perp['UNITS']   = 'nA/m^2'

;-------------------------------------------
; Velocities ///////////////////////////////
;-------------------------------------------
	;Drift velocity
	oExB_bary     = 1e3 * oE_bary -> Cross(oB_bary) / oBmag_bary^2
	oExB_bary_fac = oT_bary ## oExB_bary
	oExB_bary_fac -> SetName, exb_bary_vname
	
	oExB_bary_perp1 = 1e-3 * oExB_bary_fac[*,0]
	oExB_bary_perp1 -> SetName, exb_perp1_vname
	oExB_bary_perp1 -> Cache
	oExB_bary_perp1['LABEL'] = '(ExB)'+labels[0]+'/|B|$\up2$!C'
	oExB_bary_perp1['TITLE'] = 'V'+labels[0]+'!C(x10$\up3$km/s)'
	
	oExB_bary_perp2 = 1e-3 * oExB_bary_fac[*,1]
	oExB_bary_perp2 -> Cache
	oExB_bary_perp2 -> SetName, exb_perp2_vname
	oExB_bary_perp2['LABEL'] = '(ExB)'+labels[1]+'/|B|$\up2$'
	oExB_bary_perp2['TITLE'] = 'V'+labels[1]+'!C(x10$\up3$km/s)'
	
	;Elecron velocity
	oV = -oJ_fac / (1e18 * q * 0.03)
	oV -> SetName, v_vname
	oV -> Cache
	
	oV_perp1 = 1e-3 * oV[*,0]
	oV_perp1 -> SetName, v_perp1_vname
	oV_perp1 -> Cache
	oV_perp1['COLOR'] = 'Blue'
	oV_perp1['LABEL'] = 'J'+labels[0]+'/en'

	oV_perp2 = 1e-3 * oV[*,1]
	oV_perp2 -> SetName, v_perp2_vname
	oV_perp2 -> Cache
	oV_perp2['COLOR'] = 'Blue'
	oV_perp2['LABEL'] = 'J'+labels[1]+'/en'
;-------------------------------------------
; J.E //////////////////////////////////////
;-------------------------------------------
	
	;J.E
	oJdotE = 1e-3 * oJ -> Dot(oE_bary)
	oJdotE -> SetName, jdote_vname
	oJdotE -> Cache
	oJdotE['CATDESC'] = 'Energy dissipation, J.E, where J is determined by the ' + $
	                    'curlometer technique and E is the barycentric average ' + $
	                    'of the electric field.'
	oJdotE['LABEL']   = 'J.E'
	oJdotE['TITLE']   = 'J.E!C(nW/m^3)'
	oJdotE['UNITS']   = 'nW/m^3'
	
	;(J.E)_par
	oJdotE_par = 1e-3 * oJ_par * oE_par_bary
	oJdotE_par -> SetName, jdotepar_bary_vname
	oJdotE_par -> Cache
	oJdotE_par['CATDESC'] = 'Parallel energy dissipation, J_||.E_||, where J is determined by the ' + $
	                        'curlometer technique and E is the barycentric average ' + $
	                        'of the electric field.'
	oJdotE_par['COLOR']   = 'Green'
	oJdotE_par['LABEL']   = 'J$\down||$.E$\down||$'
	oJdotE_par['TITLE']   = 'J$\down||$.E$\down||$!C(nW/m^3)'
	oJdotE_par['UNITS']   = 'nW/m^3'
	
	;(J.E)_par
	oJdotE_par = 1e-3 * oJ_par * oE_bary_par
	oJdotE_par -> SetName, jdote_par_vname
	oJdotE_par -> Cache
	oJdotE_par['CATDESC'] = 'Parallel energy dissipation, J_||.E_||, where J is determined by the ' + $
	                        'curlometer technique and E is the barycentric average ' + $
	                        'of the electric field.'
	oJdotE_par['COLOR']   = 'Red'
	oJdotE_par['LABEL']   = 'J$\down||$.E$\down||$'
	oJdotE_par['TITLE']   = 'J$\down||$.E$\down||$!C(nW/m^3)'
	oJdotE_par['UNITS']   = 'nW/m^3'
	
	;(J.E)_perp
;	oEperp      = oE_bary -> Magnitude() - oE_bary_par
;	oJperp      = oJ -> Magnitude() - oJ_par
;	oEperp      = oE_bary_perp -> Magnitude()
;	oJperp      = oJ_perp      -> Magnitude()
;	oJdotE_perp = 1e-3 * oJperp * oEperp
	oJdotE_perp = oJdotE - oJdotE_par
	oJdotE_perp -> SetName, jdote_perp_vname
	oJdotE_perp -> Cache
	oJdotE_perp['CATDESC'] = 'Perpendicular energy dissipation, Jperp.Eperp, where J is determined by the ' + $
	                         'curlometer technique and E is the barycentric average ' + $
	                         'of the electric field.'
	oJdotE_perp['COLOR']   = 'Blue'
	oJdotE_perp['LABEL']   = 'J'+perp+'!N.E'+perp+'!N'
	oJdotE_perp['TITLE']   = 'J'+perp+'!N.E'+perp+'!N!C(nW/m^3)'
	oJdotE_perp['UNITS']   = 'nW/m^3'
	
;-------------------------------------------
; Extract Components ///////////////////////
;-------------------------------------------

	;Split into components
	xrange   = [!values.f_infinity, -!values.f_infinity]
	yrange   = [!values.f_infinity, -!values.f_infinity]
	zrange   = [!values.f_infinity, -!values.f_infinity]
	parrange = [!values.f_infinity, -!values.f_infinity]
	FOR i = 0, 3 DO BEGIN
		;Epar
		;   - Remove the error term
		oEpar = MrVar_Get(epe_vnames[i])
		oEpar = oEpar[*,1]
		oEpar -> RemoveAttr, 'LABEL'
		oEpar -> SetName, e_par_vnames[i]
		oEpar -> Cache
		
		;E_||
		oEpar['COLOR'] = sc_colors[i]
		parrange[0] <= oEpar.min
		parrange[1] >= oEpar.max		
		
		;Ex
		oEx = MrVar_Get(ex_vnames[i])
		xrange[0] <= oEx.min
		xrange[1] >= oEx.max
		
		;Ey
		oEy = MrVar_Get(ey_vnames[i])
		yrange[0] <= oEy.min
		yrange[1] >= oEy.max
		
		;Ez
		oEz = MrVar_Get(ez_vnames[i])
		zrange[0] <= oEz.min
		zrange[1] >= oEz.max
	ENDFOR
	
	;Clamp at +/-100mV/m
;	parrange[0] >= -100
;	parrange[1] <= 100
;	xrange[0]   >= -100
;	xrange[1]   <= 100
;	yrange[0]   >= -100
;	yrange[1]   <= 100
;	zrange[0]   >= -100
;	zrange[1]   <= 100
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	title = StrUpCase(StrJoin(['mms1234', mode, level], ' '))

	oEx1 = MrVar_Get(ex_vnames[0])
	oEx2 = MrVar_Get(ex_vnames[1])
	oEx3 = MrVar_Get(ex_vnames[2])
	oEx4 = MrVar_Get(ex_vnames[3])
	oEx1['PLOT_TITLE'] = title
	oEx1['TITLE']      = 'E' + labels[0] + '!C(mV/m)'
	oEx1['AXIS_RANGE'] = xrange
	oEx1['LABEL']      = 'mms1'
	oEx2['LABEL']      = 'mms2'
	oEx3['LABEL']      = 'mms3'
	oEx4['LABEL']      = 'mms4'
	
	oEy = MrVar_Get(ey_vnames[0])
	oEy['AXIS_RANGE'] = yrange
	oEy['TITLE']      = 'E' + labels[1] + '!C(mV/m)'
	
	oEz = MrVar_Get(ez_vnames[0])
	oEz['AXIS_RANGE'] = zrange
	oEz['TITLE']      = 'E' + labels[2] + '!C(mV/m)'
	
	oEpar = MrVar_Get(e_par_vnames[0])
	oEpar['AXIS_RANGE'] = parrange
	oEpar['TITLE']      = 'E$\down||$!C(mV/m)'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [ ex_vnames[0], ey_vnames[0], ez_vnames[0], e_par_vnames[0], $
	                      j_vname, jdote_vname, exb_perp1_vname, exb_perp2_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	win = MrVar_OPlotTS( ex_vnames[0],    ex_vnames[1:3] )
	win = MrVar_OPlotTS( ey_vnames[0],    ey_vnames[1:3] )
	win = MrVar_OPlotTS( ez_vnames[0],    ez_vnames[1:3] )
	win = MrVar_OPlotTS( e_par_vnames[0], e_par_vnames[1:3] )
	win = MrVar_OPlotTS( jdote_vname,     [jdote_par_vname, jdote_perp_vname, jdotepar_bary_vname] )
	win = MrVar_OPlotTS( exb_perp1_vname, v_perp1_vname )
	win = MrVar_OPlotTS( exb_perp2_vname, v_perp2_vname )
	
	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 10]
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
		fname   = StrJoin( ['mms1234', edp_instr, mode, level, '4sc-jdote-fac'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END