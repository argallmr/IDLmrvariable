; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_EMaxwell
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
;       1. B barycentric average
;       2. Curl(E)x & -dBx/dt
;       3. Curl(E)y & -dBy/dt
;       4. Curl(E)z & -dBz/dt
;       5. Charge density: e0 * Div(E)
;       6. Ex barycenter & -Grad(V)x
;       7. Ey barycenter & -Grad(V)y
;       8. Ez barycenter & -Grad(V)z
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
;       2018/02/10  -   Written by Matthew Argall
;       2018/05/07  -   Filter individual spacecraft fields before gradients/averaging. - MRA
;-
FUNCTION MrMMS_Plot_4sc_EMaxwell, mode, $
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
	fgm_instr = mode EQ 'brst' ? 'fsm' : 'fgm'
	IF N_Elements(coords) EQ 0 THEN coords = 'gse'
	IF N_Elements(fc)     EQ 0 THEN fc     = 0.0
	IF N_Elements(level)  EQ 0 THEN level  = 'l2'
	IF N_Elements(trange) GT 0 THEN MrVar_SetTRange, trange
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	
	e0        = MrConstants('epsilon_0')
	q         = MrConstants('q')
	sc_colors = ['Black', 'Red', 'Green', 'Blue']
	nabla     = '!9'+String(71B)+'!X'
	partial   = '!9'+String(68B)+'!X'
	
;-------------------------------------------
; Variable Parameters //////////////////////
;-------------------------------------------
	;FGM
	fgm_mode    = mode
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
	v_vnames    = sc + '_' + StrJoin([edp_instr, 'scpot',            edp_mode, level], '_' )
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	bf_vnames     = sc + '_' + StrJoin([fgm_instr, 'bvec', fgm_coords, 'lpfilt', fgm_mode, fgm_level], '_')
	ef_vnames     = sc + '_' + StrJoin([edp_instr, 'dce',  edp_coords, 'lpfilt', edp_mode, level], '_' )
	vf_vnames     = sc + '_' + StrJoin([edp_instr, 'scpot', 'lpfilt',            edp_mode, level], '_' )
	bmag_bary_vname = StrJoin( ['mms', 'fgm', 'bmag', fgm_coords, 'bary', fgm_mode, fgm_level], '_' )
	b_bary_vname  = StrJoin( ['mms', 'fgm', 'b', fgm_coords, 'bary', fgm_mode, fgm_level], '_' )
	bx_bary_vname = b_bary_vname + '_x'
	by_bary_vname = b_bary_vname + '_y'
	bz_bary_vname = b_bary_vname + '_z'
;	bf_bary_vname = StrJoin( ['mms', 'fgm', 'b', fgm_coords, 'bary', 'lpfilt', fgm_mode, fgm_level], '_' )
	dBdt_vname    = StrJoin( ['mms', 'fgm', 'bdot', fgm_coords, 'bary', fgm_mode, fgm_level], '_' )
	dBxdt_vname   = dBdt_vname + '_x'
	dBydt_vname   = dBdt_vname + '_y'
	dBzdt_vname   = dBdt_vname + '_z'
	e_bary_vname  = StrJoin( ['mms', edp_instr, 'e', edp_coords, 'bary', edp_mode, level], '_' )
	ex_bary_vname = e_bary_vname + '_x'
	ey_bary_vname = e_bary_vname + '_y'
	ez_bary_vname = e_bary_vname + '_z'
	gradv_vname   = StrJoin( ['mms', 'dce', 'gradv', edp_mode, level], '_' )
	gradvx_vname  = gradv_vname + '_x'
	gradvy_vname  = gradv_vname + '_y'
	gradvz_vname  = gradv_vname + '_z'
	rho_vname     = StrJoin( ['mms', 'dce', 'rho',   edp_mode, level], '_' )
	curlE_vname   = StrJoin( ['mms', 'dce', 'curle', edp_mode, level], '_' )
	curlEx_vname  = curlE_vname + '_x'
	curlEy_vname  = curlE_vname + '_y'
	curlEz_vname  = curlE_vname + '_z'
;	curlEf_vname  = StrJoin( ['mms', 'dce', 'curle', 'lpfilt', edp_mode, level], '_' )
		
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
			MrMMS_FGM_Load_Data, '', fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = fgm_level, $
			                     OPTDESC   = fgm_optdesc, $
			                     VARFORMAT = '*_b_*'
		ENDELSE
		
		;E-Field
		MrMMS_Load_Data, '', instr, mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = '*dce_'+coords+'*'
		
		;Spacecraft Potential
		MrMMS_Load_Data, '', instr, mode, level, $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
		                 
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
; Interpolate //////////////////////////////
;-------------------------------------------
	oBmag1 = MrVar_Get(bmag_vnames[0])
	oBmag2 = MrVar_Get(bmag_vnames[1])
	oBmag3 = MrVar_Get(bmag_vnames[2])
	oBmag4 = MrVar_Get(bmag_vnames[3])
	oB1 = MrVar_Get(bvec_vnames[0])
	oB2 = MrVar_Get(bvec_vnames[1])
	oB3 = MrVar_Get(bvec_vnames[2])
	oB4 = MrVar_Get(bvec_vnames[3])
	oE1 = MrVar_Get(e_vnames[0])
	oE2 = MrVar_Get(e_vnames[1])
	oE3 = MrVar_Get(e_vnames[2])
	oE4 = MrVar_Get(e_vnames[3])
	oV1 = MrVar_Get(v_vnames[0])
	oV2 = MrVar_Get(v_vnames[1])
	oV3 = MrVar_Get(v_vnames[2])
	oV4 = MrVar_Get(v_vnames[3])
	oR1 = MrVar_Get(r_vnames[0])
	oR2 = MrVar_Get(r_vnames[1])
	oR3 = MrVar_Get(r_vnames[2])
	oR4 = MrVar_Get(r_vnames[3])
	
	oB1 = oB1 -> Interpol(oE1)
	oB2 = oB2 -> Interpol(oE1)
	oB3 = oB3 -> Interpol(oE1)
	oB4 = oB4 -> Interpol(oE1)
	oE2 = oE2 -> Interpol(oE1)
	oE3 = oE3 -> Interpol(oE1)
	oE4 = oE4 -> Interpol(oE1)
	oV1 = oV1 -> Interpol(oE1)
	oV2 = oV2 -> Interpol(oE1)
	oV3 = oV3 -> Interpol(oE1)
	oV4 = oV4 -> Interpol(oE1)
	oR1 = oR1 -> Interpol(oE1)
	oR2 = oR2 -> Interpol(oE1)
	oR3 = oR3 -> Interpol(oE1)
	oR4 = oR4 -> Interpol(oE1)
	
;-------------------------------------------
; Filter Data //////////////////////////////
;-------------------------------------------
	dtB = oB1['TIMEVAR'] -> GetSI(RATE=fsB)
	dtE = oE1['TIMEVAR'] -> GetSI(RATE=fsE)
	dtV = oV1['TIMEVAR'] -> GetSI(RATE=fsV)
	
	IF fc GT 0 THEN BEGIN
		fN = fsB / 2.0
		f0 = 0.0
		f1 = fc / fN
		A  = 75
		N  = Round(fsB)
		oB1 = oB1 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=bf_vnames[0])
		oB2 = oB2 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=bf_vnames[1])
		oB3 = oB3 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=bf_vnames[2])
		oB4 = oB4 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=bf_vnames[3])
		
		;Filter
		fN = fsE / 2.0
		f0 = 0.0
		f1 = fc / fN
		A  = 50
		N  = Round(fsE)
		oE1 = oE1 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=ef_vnames[0])
		oE2 = oE2 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=ef_vnames[1])
		oE3 = oE3 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=ef_vnames[2])
		oE4 = oE4 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=ef_vnames[3])
		
		;Filter
		fN = fsE / 2.0
		f0 = 0.0
		f1 = fc / fN
		A  = 50
		N  = Round(fsE)
		oV1 = oV1 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=vf_vnames[0])
		oV2 = oV2 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=vf_vnames[1])
		oV3 = oV3 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=vf_vnames[2])
		oV4 = oV4 -> Digital_Filter(f0, f1, A, N, /CACHE, NAME=vf_vnames[3])
	ENDIF
	
	
;-------------------------------------------
; Barycentric Averages /////////////////////
;-------------------------------------------
	;B
	oB_bary = (oB1 + oB2 + oB3 + oB4) / 4.0
	oB_bary -> SetName, b_bary_vname
	oB_bary -> Cache
	
	;Split into components
	oB_bary -> Split, oBx_bary, oBy_bary, oBz_bary, /CACHE
	oB_bary['COLOR'] = ['Blue', 'Forest Green', 'Red']
	oB_bary['LABEL'] = 'B$\down' + ['X', 'Y', 'Z'] + '$'
	oB_bary['TITLE'] = 'B$\downBC$!C(nT)'
	
	;|B|
	oBmag_bary = oB_bary -> Magnitude(/CACHE, NAME=bmag_bary_vname)
	oBmag_bary['LABEL'] = '|B|'
	oBmag_bary['TITLE'] = '|B|!C(nT)'
	
	;E
	oE_bary = (oE1 + oE2 + oE3 + oE4) / 4.0
	oE_bary -> SetName, e_bary_vname
	oE_bary -> Cache
	
	;Split into components
	oE_bary -> Split, oEx_bary, oEy_bary, oEz_bary, /CACHE
	oE_bary['COLOR']      = ['Blue', 'Forest Green', 'Red']
	oE_bary['LABEL']      = 'E$\down' + ['X', 'Y', 'Z'] + '$'
	oE_bary['TITLE']      = 'E$\downBC$!C(mV/m^2)'
	
	oEx_bary['LABEL'] = 'E$\downBC,X$'
	oEx_bary['TITLE'] = 'E!C(mV/m)'
	oEy_bary['LABEL'] = 'E$\downBC,Y$'
	oEy_bary['TITLE'] = 'E!C(mV/m)'
	oEz_bary['LABEL'] = 'E$\downBC,Z$'
	oEz_bary['TITLE'] = 'E!C(mV/m)'
	
;-------------------------------------------
; dB/dt ////////////////////////////////////
;-------------------------------------------
	;dB/dt
	;   - 1e-3 converts to uT/s = kg/Cs^2 * 1e6
	dB_dt = -1e-3 * ( Double(oB_bary['DATA', 1:-1, *]) - Double(oB_bary['DATA', 0:-2, *]) ) / dtB
	odB_bary_dt = MrVectorTS( (oB_bary['TIMEVAR'])[1:*], dB_dt )
	
	;Set properties
	odB_bary_dt -> SetName, dBdt_vname
	odB_bary_dt -> Cache
	odB_bary_dt['COLOR'] = ['Blue', 'Forest Green', 'Red']
	odB_bary_dt['LABEL'] = '-dB$\down' + ['X', 'Y', 'Z'] + '$/dt'
	odB_bary_dt['TITLE'] = '-dB/dt!C($\mu$T/s)'
	
	;Split into components
	odB_bary_dt -> Split, odBx_dt, odBy_dt, odBz_dt
	
	;dBx/dt
	odBx_dt -> SetName, dBxdt_vname
	odBx_dt -> Cache
	odBx_dt['COLOR'] = 'Blue'
	odBx_dt['LABEL'] = '-dB$\downBC,X$/dt'
	odBx_dt['TITLE'] = '$\mu$T/s'
	
	;dBy/dt
	odBy_dt -> SetName, dBydt_vname
	odBy_dt -> Cache
	odBy_dt['COLOR'] = 'Blue'
	odBy_dt['LABEL'] = '-dB$\downBC,Y$/dt'
	odBy_dt['TITLE'] = '$\mu$T/s'
	
	;dBz/dt
	odBz_dt -> SetName, dBzdt_vname
	odBz_dt -> Cache
	odBz_dt['COLOR'] = 'Blue'
	odBz_dt['LABEL'] = '-dB$\downBC,Z$/dt'
	odBz_dt['TITLE'] = '$\mu$T/s'
	
;-------------------------------------------
; Div(E) ///////////////////////////////////
;-------------------------------------------
	oRecipVec = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	
	;Charge Density
	;   - 1e4 converts to C/cm^2
	oRho = (e0 / q * 1e-12) * oRecipVec -> Divergence( oE1, oE2, oE3, oE4 )
	oRho -> SetName, rho_vname
	oRho -> Cache
	oRho['CATDESC'] = 'Free charge density, computed from the divergence of the electric field.'
	oRho['TITLE']   = '$\Delta$#!C(#/cm^3)'
	oRho['UNITS']   = '#/cm^3'
	
;-------------------------------------------
; Curl(E) //////////////////////////////////
;-------------------------------------------
	
	;Curl of E
	;   - 1e0 converts to uV/m^2 = kg/Cs^2 * 1e6
	oCurl =  oRecipVec -> Curl( oE1, oE2, oE3, oE4 )
	oCurl -> SetName, curlE_vname
	oCurl -> Cache
	oCurl['COLOR'] = ['Blue', 'Forest Green', 'Red']
	oCurl['LABEL'] = '(' + nabla + 'xE)$\down' + ['X', 'Y', 'Z'] + '$'
	oCurl['TITLE'] = nabla + 'xE!C($\mu$V/m^2)'
	
	;Attributes
	oCurl['CATDESC'] = 'Curl(E)'
	oCurl['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oCurl['LABEL']   = '(' + nabla + 'xE)$\down' + ['X', 'Y', 'Z'] + '$'
	oCurl['TITLE']   = nabla + 'xE!C($\mu$V/m^2)'
	oCurl['UNITS']   = ''
	
	;Split Curl(E) into components
	oCurl -> Split, oCurlEx, oCurlEy, oCurlEz
	
	;Curl(E)x
	oCurlEx -> SetName, curlEx_vname
	oCurlEx -> Cache
	oCurlEx['LABEL'] = '(' + nabla + 'xE)$\downX$'
	oCurlEx['TITLE'] = '$\mu$V/m$\up2$'
	
	;Curl(E)y
	oCurlEy -> SetName, curlEy_vname
	oCurlEy -> Cache
	oCurlEy['LABEL'] = '(' + nabla + 'xE)$\downY$'
	oCurlEy['TITLE'] = '$\mu$V/m$\up2$'
	
	;Curl(E)z
	oCurlEz -> SetName, curlEz_vname
	oCurlEz -> Cache
	oCurlEz['LABEL'] = '(' + nabla + 'xE)$\downZ$'
	oCurlEz['TITLE'] = '$\mu$V/m$\up2$'
	
;-------------------------------------------
; -Grad(V) /////////////////////////////////
;-------------------------------------------
	;E = -Grad(V)
	;   - 1e0 converts V/km to mV/m
	oGradV = -oRecipVec -> Gradient( oV1, oV2, oV3, oV4 )
	oGradV -> SetName, gradv_vname
	oGradV -> Cache
	oGradV['CATDESC'] = '-Grad(V)'
	oGradV['COLOR']   = ['Blue', 'Green', 'Red']
	oGradV['LABEL']   = ['X', 'Y', 'Z']
	oGradV['TITLE']   = '-Grad(V)!C(mV/m)'
	oGradV['UNITS']   = 'mV/m'
	
	;Split into components
	oGradV -> Split, oGradVx, oGradVy, oGradVz, /CACHE
	
	oGradVx['COLOR']     = 'Magenta'
	oGradVx['LABEL']     = '(' + nabla + 'V)$\downX$'
;	oGradVx['LINESTYLE'] = '--'
	oGradVx['TITLE']     = 'E$\downX$!C(mV/m)'
	
	oGradVy['COLOR']     = 'Magenta'
	oGradVy['LABEL']     = '(' + nabla + 'V)$\downY$'
;	oGradVy['LINESTYLE'] = '--'
	oGradVy['TITLE']     = 'E$\downY$!C(mV/m)'
	
	oGradVz['COLOR']     = 'Magenta'
	oGradVz['LABEL']     = '(' + nabla + 'V)$\downZ$'
;	oGradVz['LINESTYLE'] = '--'
	oGradVz['TITLE']     = 'E$\downZ$!C(mV/m)'
	
	Obj_Destroy, oRecipVec
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
;	title = StrUpCase(StrJoin(['mms1234', instr, mode, level], ' '))
;
;	oEx1 = MrVar_Get(ex_vnames[0])
;	oEx2 = MrVar_Get(ex_vnames[1])
;	oEx3 = MrVar_Get(ex_vnames[2])
;	oEx4 = MrVar_Get(ex_vnames[3])
;	odVx = MrVar_Get(gradvx_vname)
;	oEx1['PLOT_TITLE'] = title
;	oEx1['TITLE'] = 'Ex!C(mV/m)'
;	oEx1['LABEL'] = 'mms1'
;	oEx2['LABEL'] = 'mms2'
;	oEx3['LABEL'] = 'mms3'
;	oEx4['LABEL'] = 'mms4'
;	odVx['LABEL'] = '-Grad(V)'
;	
;	oEy = MrVar_Get(ey_vnames[0])
;	oEy['TITLE'] = 'Ey!C(mV/m)'
;	
;	oEz = MrVar_Get(ez_vnames[0])
;	oEz['TITLE'] = 'Ez!C(mV/m)'
;	
;	oEpar = MrVar_Get(epar_vnames[0])
;	oEpar['TITLE'] = 'E$\down||$!C(mV/m)'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [ b_bary_vname, curlEx_vname, curlEy_vname, curlEz_vname, rho_vname, $
	                      ex_bary_vname, ey_bary_vname, ez_bary_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	win = MrVar_OPlotTS( b_bary_vname, bmag_bary_vname )
	win = MrVar_OPlotTS( curlEx_vname, dBxdt_vname )
	win = MrVar_OPlotTS( curlEy_vname, dBydt_vname )
	win = MrVar_OPlotTS( curlEz_vname, dBzdt_vname )
	win = MrVar_OPlotTS( ex_bary_vname, gradvx_vname )
	win = MrVar_OPlotTS( ey_bary_vname, gradvy_vname )
	win = MrVar_OPlotTS( ez_bary_vname, gradvz_vname )
	
	;Pretty-up the window
	win.name = 'EMaxwell'
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 11]
	win    -> Refresh

;-------------------------------------------
; E vs. -Grad(V) ///////////////////////////
;-------------------------------------------
	;LADFIT variables
;	px = LADFit(oGradVx['DATA'], oEx_bary['DATA'])
;	py = LADFit(oGradVy['DATA'], oEy_bary['DATA'])
;	pz = LADFit(oGradVz['DATA'], oEz_bary['DATA'])

	;Scatter plots of Grad(V) and E_bary to view gain and offset parameters
	w2 = MrWindow( ASPECT  = 1.0, $
	               LAYOUT  = [3,1], $
	               NAME    = 'EvGradV', $
	               REFRESH = 0, $
	               XSIZE   = 1000 )
	p1 = MrVar_Plot( gradvx_vname, ex_bary_vname, $
	                 XTITLE        = '(' + nabla + 'V)$\downX$ (mV/m)', $
	                 XTICKINTERVAL = 0.04, $
	                 YTITLE        = 'E$\downX$ (mV/m)', $
	                 /CURRENT )
;	p1f = MrPlot( oGradVx['DATA'], px[0] + px[1]*oGradVx['DATA'], $
;	              COLOR    = 'Blue', $
;	              OVERPLOT = p1 )
	
	p2 = MrVar_Plot( gradvy_vname, ey_bary_vname, $
	                 XTITLE = '(' + nabla + 'V)$\downY$ (mV/m)', $
	                 YTITLE = 'E$\downY$ (mV/m)', $
	                 /CURRENT )
	p3 = MrVar_Plot( gradvz_vname, ez_bary_vname, $
	                 XTITLE = '(' + nabla + 'V)$\downZ$ (mV/m)', $
	                 YTITLE = 'E$\downZ$ (mV/m)', $
	                 /CURRENT )
	p1 -> SetLayout, [1,1]
	trange   = MrVar_GetTRange()
	p2.title = StrJoin(StrSplit(StrMid(trange[0], 0, 19), 'T', /EXTRACT), ' ') + ' -- ' + StrMid(trange[1], 11, 8)
	w2 -> TrimLayout
	w2 -> Remove, w2 -> Get(/ALL, ISA='MrLegend')
	w2 -> Refresh
;	w2 -> Save, '/home/argall/figures/20151206/mms_edp_brst_l2_4sc-e-maxwell-gradv-scatter_20151206_233827_233836.png'

;-------------------------------------------
; Curl(E) vs. -dB/dt ///////////////////////
;-------------------------------------------

	;Scatter plots of -dB/dt and Curl(E) to view gain and offset parameters
	oCurlEx = MrVar_Get(curlEx_vname)
	oCurlEy = MrVar_Get(curlEy_vname)
	oCurlEz = MrVar_Get(curlEz_vname)
	w3 = MrWindow( ASPECT  = 1.0, $
	               LAYOUT  = [3,1], $
	               NAME    = 'CurlEvdBdt', $
	               REFRESH = 0, $
	               XSIZE   = 1000 )
	p4 = MrVar_Plot( dBxdt_vname, oCurlEx[1:-1], $
	                 TITLE         = '', $
	                 XRANGE        = [-0.03,0.03], $
	                 XTITLE        = '-'+partial+'B$\downX$/'+partial+'t ($\mu$T/s)', $
	                 XTICKINTERVAL = 0.02, $
	                 YTITLE        = nabla+'E$\downX$ ($\mu$V/m$\up2$)', $
	                 /CURRENT )
	p5 = MrVar_Plot( dBydt_vname, oCurlEy[1:-1], $
	                 XRANGE = [-0.04, 0.01], $
	                 XTITLE = '-'+partial+'B$\downY$/'+partial+'t ($\mu$T/s)', $
	                 YTITLE = nabla+'E$\downZ$ ($\mu$V/m$\up2$)', $
	                 /CURRENT )
	p6 = MrVar_Plot( dBzdt_vname, oCurlEz[1:-1], $
	                 TITLE  = '', $
	                 XRANGE = [-0.04,0.01], $
	                 XTITLE = '-'+partial+'B$\downZ$/'+partial+'t ($\mu$T/s)', $
	                 YTITLE = nabla+'E$\downZ$ ($\mu$V/m$\up2$)', $
	                 /CURRENT )
	p4 -> SetLayout, [1,1]
	trange   = MrVar_GetTRange()
	p5.title = StrJoin(StrSplit(StrMid(trange[0], 0, 19), 'T', /EXTRACT), ' ') + ' -- ' + StrMid(trange[1], 11, 8)
	w3 -> TrimLayout
	w3 -> Remove, w3 -> Get(/ALL, ISA='MrLegend')
	w3 -> Refresh
;	w3 -> Save, '/home/argall/figures/20151206/mms_edp_brst_l2_4sc-e-maxwell-curlE-scatter_20151206_233827_233836.png'

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
		fname   = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-e-maxwell'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END