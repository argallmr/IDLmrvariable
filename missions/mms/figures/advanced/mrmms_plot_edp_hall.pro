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
;       1. Ex, VixB, VexB, JxB/ne
;       2. Ey, VixB, VexB, JxB/ne
;       3. Ez, VixB, VexB, JxB/ne
;       4. ni, ne
;       5. Bxyz, |B|
;       6. J (moments)
;       7. Rho (Div E)
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
FUNCTION MrMMS_Plot_EDP_Hall, sc, mode, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
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
	
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	IF mode EQ 'brst' THEN BEGIN
		IF (mrmms_get_filenames(sc, 'mec', mode, level, OPTDESC=ephdesc))[0] EQ '' $
			THEN mec_mode = 'srvy' $
			ELSE mec_mode = mode
	ENDIF
	
	;Constants
	q     = MrConstants('q')
	e0    = MrConstants('epsilon_0')
	perp  = '!9' + String(120B) + '!X'
	nabla = '!9' + String(71B) + '!X'
	cdot  = '!9' + String(46B) + '!X'
	
;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	;FGM
	CASE fgm_instr OF
		'afg': fgm_level = 'l2pre'
		'dfg': fgm_level = 'l2pre'
		'fgm': fgm_level = 'l2'
		'fsm': fgm_level = 'l3'
		ELSE: Message, 'Invalid FGM instrument: "' + fgm_instr + '".'
	ENDCASE
	CASE coords OF
		'dsl':  fgm_coords = 'dmpa'
		'dbcs': fgm_coords = 'dmpa'
		'gse':  fgm_coords = coords
		'gsm':  fgm_coords = coords
		ELSE:   fgm_coords = 'gse'
	ENDCASE
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	
	;EDP
	edp_mode    = mode EQ 'brst' ? mode : 'fast'
	edp_optdesc = 'dce'
	CASE coords OF
		'dmpa': edp_coords = 'dsl'
		'dbcs': edp_coords = 'dsl'
		'gse':  edp_coords = coords
		'gsm':  edp_coords = coords
		ELSE:   edp_coords = 'gse'
	ENDCASE
	
	;FPI
	fpi_mode    = mode EQ 'brst' ? mode : 'fast'
	des_optdesc = 'des-moms'
	dis_optdesc = 'dis-moms'
	CASE coords OF
		'dmpa': fpi_coords = 'dsl'
		'dsl':  fpi_coords = 'dbcs'
		'gse':  fpi_coords = coords
		'gsm':  fpi_coords = coords
		ELSE:   fpi_coords = 'gse'
	ENDCASE
	
	;EPH
	CASE coords OF
		'dmpa': Message, 'EPHEMERIS products do not have coordinate system "' + coords + '".'
		'dsl':  Message, 'EPHEMERIS products do not have coordinate system "' + coords + '".'
		'dmpa': Message, 'EPHEMERIS products do not have coordinate system "' + coords + '".'
		'gse':  eph_coords = coords
		'gsm':  eph_coords = coords
		ELSE:   eph_coords = 'gse'
	ENDCASE
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	all_sc = 'mms' + ['1', '2', '3', '4']
	nSC    = N_Elements(all_sc)
	
	;FGM
	b_vnames    = all_sc + '_' + StrJoin([fgm_instr, 'b',    fgm_coords, mode, fgm_level], '_')
	bmag_vnames = all_sc + '_' + StrJoin([fgm_instr, 'bmag', fgm_coords, mode, fgm_level], '_')
	bvec_vnames = all_sc + '_' + StrJoin([fgm_instr, 'bvec', fgm_coords, mode, fgm_level], '_')
	
	;EDP
	e_vnames = all_sc + '_' + StrJoin( ['edp', 'dce', edp_coords, edp_mode, level], '_' )
	
	;MEC
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames = all_sc + '_' + StrJoin( ['mec', 'r', eph_coords], '_' ) $
		ELSE r_vnames = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;FPI
	ni_vnames = all_sc + '_' + StrJoin(['dis', 'numberdensity',             fpi_mode], '_')
	vi_vnames = all_sc + '_' + StrJoin(['dis', 'bulkv',         fpi_coords, fpi_mode], '_')
	ne_vnames = all_sc + '_' + StrJoin(['des', 'numberdensity',             fpi_mode], '_')
	ve_vnames = all_sc + '_' + StrJoin(['des', 'bulkv',         fpi_coords, fpi_mode], '_')
	
	;Output names
	CASE StrUpCase(sc) OF
		'MMS1': e_vname = e_vnames[0]
		'MMS2': e_vname = e_vnames[1]
		'MMS3': e_vname = e_vnames[2]
		'MMS4': e_vname = e_vnames[3]
		ELSE: Message, 'Incorrect value for SC: "' + sc + '".'
	ENDCASE
	ex_vname    = e_vname + '_x'
	ey_vname    = e_vname + '_y'
	ez_vname    = e_vname + '_z'
	eperp_vname = StrJoin( [sc, 'edp', 'dce', 'perp', edp_mode, level], '_' )
	eperpx_vname = eperp_vname + '_x'
	eperpy_vname = eperp_vname + '_y'
	eperpz_vname = eperp_vname + '_z' 
	jxb_vnames   = all_sc + '_' + StrJoin( ['fpi', 'jxb', edp_coords, fpi_mode, level], '_' )
	jxbx_vnames  = jxb_vnames + '_x'
	jxby_vnames  = jxb_vnames + '_y'
	jxbz_vnames  = jxb_vnames + '_z'
	vixb_vname  = StrJoin( [sc, 'dis', 'vxb', edp_coords, edp_mode, level], '_' )
	vixbx_vname = vixb_vname + '_x'
	vixby_vname = vixb_vname + '_y'
	vixbz_vname = vixb_vname + '_z'
	vexb_vname  = StrJoin( [sc, 'des', 'vxb', edp_coords, edp_mode, level], '_' )
	vexbx_vname = vexb_vname + '_x'
	vexby_vname = vexb_vname + '_y'
	vexbz_vname = vexb_vname + '_z'
	j_vnames    = all_sc + '_' + StrJoin(['fpi', 'j', fpi_coords, fpi_mode], '_')
	jxb_vnames  = all_sc + '_' + StrJoin(['fpi', 'jxb', fpi_coords, fpi_mode], '_')
	jxbx_vnames = jxb_vnames + '_x'
	jxby_vnames = jxb_vnames + '_y'
	jxbz_vnames = jxb_vnames + '_z'
	rho_vname   = StrJoin( ['mms', 'edp', 'rho', mode, level], '_' )
	div_jxb_vname = StrJoin( ['mms', 'fpi', 'divjxb', mode, level], '_' )
		
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, '', mode, $
		                     INSTR   = fgm_instr, $
		                     LEVEL   = fgm_level, $
		                     OPTDESC = fgm_optdesc, $
		                     VARFORMAT = '*_b_'+fgm_coords+'_'+mode+'*'
		
		;E-Field
		MrMMS_Load_Data, '', 'edp', edp_mode, level, $
		                 OPTDESC   = edp_optdesc, $
		                 VARFORMAT = ['*dce_'+edp_coords+'*', '*dce_par*']

		;DES
		MrMMS_FPI_Load_Data, '', fpi_mode, $
		                     OPTDESC   = des_optdesc, $
		                     VARFORMAT = ['*density_'+fpi_mode, '*bulkv_'+fpi_coords+'_'+fpi_mode]
		                     
		;DIS
		MrMMS_FPI_Load_Data, '', fpi_mode, $
		                     OPTDESC   = dis_optdesc, $
		                     VARFORMAT = ['*density_'+fpi_mode, '*bulkv_'+fpi_coords+'_'+fpi_mode]
		
		;Ephemeris
		MrMMS_Load_Data, '', 'mec', mec_mode, level, $
		                 OPTDESC   = ephdesc, $
		                 VARFORMAT = r_vnames
	ENDIF
	
	idx = Fix(StrMid(sc, 3, 1)) - 1
	
;-------------------------------------------
; Free Charge Density //////////////////////
;-------------------------------------------
	;Grab E & R
	oNe = MrVar_Get(ne_vnames[idx])
	oTref = oNe['DEPEND_0']
	
	oE = ObjArr(4)
	oPos = ObjArr(4)
	FOR i = 0, 3 DO BEGIN
		oE[i]   = MrVar_Resample(e_vnames[i], oTref)
		oPos[i] = MrVar_Resample(r_vnames[i], oTref)
	ENDFOR
	
	;Create the reciprocal vectors
	oRecipVec = MrVar_RecipVec(oPos[0], oPos[1], oPos[2], oPos[3])
	
	;Charge Density
	;   - 1e4 converts to C/cm^2
	oRho = (e0 / q * 1e-12) * oRecipVec -> Divergence( oE[0], oE[1], oE[2], oE[3] )
	oRho -> SetName, rho_vname
	oRho -> Cache
	oRho['CATDESC'] = 'Free charge density, computed from the divergence of the electric field.'
	oRho['LABEL']   = '$\epsilon$$\down0$/e' + nabla + cdot + 'E'
	oRho['TITLE']   = '$\rho$/e!C(cm$\up-3$)'
	oRho['UNITS']   = 'cm^-3'
	
	Obj_Destroy, [oRecipVec, oE, oPos]
	
;-------------------------------------------
; J, JxB, Div(JxB) /////////////////////////
;-------------------------------------------
	oB   = ObjArr(nSC)
	oE   = ObjArr(nSC)
	oNe  = ObjArr(nSC)
	oVe  = ObjArr(nSC)
	oNi  = ObjArr(nSC)
	oVi  = ObjArr(nSC)
	oJ   = ObjArr(nSC)
	oJxB = ObjArr(nSC)
	oPos = ObjArr(nSC)
	FOR i = 0, 3 DO BEGIN
		;Interpolate to DES
		oB[i]   = MrVar_Resample(bvec_vnames[i], oTref)
		oE[i]   = MrVar_Resample(e_vnames[i],    oTref)
		oNe[i]  = MrVar_Resample(ne_vnames[i],   oTref)
		oVe[i]  = MrVar_Resample(ve_vnames[i],   oTref)
		oNi[i]  = MrVar_Resample(ni_vnames[i],   oTref)
		oVi[i]  = MrVar_Resample(vi_vnames[i],   oTref)
		oPos[i] = MrVar_Resample(r_vnames[i],    oTref)
	
		;Compute current density
		;   - 1e15 converts to uA/m^2
		q  = MrConstants('q')
		oJ[i] = q * 1e15 * oNe[i] * (oVi[i] - oVe[i])
		oJ[i] -> SetName, j_vnames[i]
		oJ[i] -> Cache
		
		;Hall field
		;   - 1e-18 converts to mV/m
		oJxB[i] = 1.0/(1e18 * q * oNe[i]) * oJ[i] -> Cross(oB[i])
		oJxB[i] -> SetName, jxb_vnames[i]
		oJxB[i] -> Cache
	ENDFOR
	
	;Create the reciprocal vectors
	oRecipVec = MrVar_RecipVec(oPos[0], oPos[1], oPos[2], oPos[3])
	
	;Charge Density
	;   - 1e4 converts to C/cm^2
	oDivJxB = (e0 / q * 1e-12) * oRecipVec -> Divergence( oJxB[0], oJxB[1], oJxB[2], oJxB[3] )
	oDivJxB -> SetName, div_jxb_vname
	oDivJxB -> Cache
	oDivJxB['CATDESC'] = 'Free charge density, computed from the divergence of JxB.'
	oDivJxB['COLOR']   = 'Blue'
	oDivJxB['LABEL']   = '$\epsilon$$\down0$/e' + nabla + cdot + '(JxB)'
	oDivJxB['TITLE']   = '$\rho$/e!C(cm$\up-3$)'
	oDivJxB['UNITS']   = 'cm^-3'
	
	Obj_Destroy, oRecipVec
	
;-------------------------------------------
; Electric Fields //////////////////////////
;-------------------------------------------
	;Perpendicular field
	;  - Roy's method does not work...
;	oB_hat  = oB[idx] -> Normalize()
;	oE_perp = oE[idx] * (1.0 - oB_hat)
	oTx = MrVar_FAC(oB[idx])
	oE_perp      = oTx ## oE[idx]
	oE_perp[*,2] = 0.0
	oE_perp      = (oTx -> Transpose()) ## oE_perp
	oE_perp -> SetName, eperp_vname
	oE_perp -> Cache
	
	;Convective field
	;   - 1e-3 converts to mV/m
	oVexB = -1e-3 * oVe[idx] -> Cross(oB[idx])
	oVexB -> SetName, vexb_vname
	oVexB -> Cache
	
	;Convective field
	;   - 1e-3 converts to mV/m
	oVixB = -1e-3 * oVi[idx] -> Cross(oB[idx])
	oVixB -> SetName, vixb_vname
	oVixB -> Cache
	
;-------------------------------------------
; Rotate ///////////////////////////////////
;-------------------------------------------
	IF coords EQ 'mva' THEN BEGIN
		oB[idx]   = MrVar_xForm(oB[idx])
		oE_perp   = MrVar_xForm(oE_perp)
		oJxB[idx] = MrVar_xForm(oJxB[idx], /CACHE, NAME=jxb_vnames[idx])
		oVexB     = MrVar_xForm(oVexB,     /CACHE, NAME=vexb_vname)
		oVixB     = MrVar_xForm(oVixB,     /CACHE, NAME=vixb_vname)
		oJ[idx]   = MrVar_xForm(oJ[idx],   /CACHE, NAME=j_vnames[idx])
		comps   = ['N', 'M', 'L']
	ENDIF ELSE BEGIN
		oB[idx] = oB[idx] -> Copy()
		comps   = ['X', 'Y', 'Z']
	ENDELSE
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	;Magnitude
	oBmag = MrVar_Get(bmag_vnames[idx])
	oBmag['AXIS_RANGE'] = [(oB[idx]).min, oBmag.max]
	oBmag['TITLE']      = 'B!C(nT)'
	
	;Magnetic Field
	(oB[idx])['CATDESC'] = 'Vector magnetic field.'
	(oB[idx])['COLOR']   = ['Blue', 'Forest Green', 'Red']
	(oB[idx])['LABEL']   = comps
	(oB[idx])['TITLE']   = 'B!C(nT)'
	(oB[idx])['UNITS']   = 'nT'
	
	;Current Density
	(oJ[idx])['CATDESC'] = 'Current density calculated from particle moments: J=q*ne*(Vi-Ve).'
	(oJ[idx])['COLOR']   = ['Blue', 'Forest Green', 'Red']
	(oJ[idx])['LABEL']   = comps
	(oJ[idx])['TITLE']   = 'J!C($\mu$A/m^2)'
	(oJ[idx])['UNITS']   = 'uA/m^2'
	
	;Pepr E-field
	oE_perp['CATDESC'] = 'Electric field perpendicular to the magnetic field.'
	oE_perp['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oE_perp['LABEL']   = comps
	oE_perp['TITLE']   = 'E$\down' + perp + '$!C(mV/m)'
	oE_perp['UNITS']   = 'mV/m'
	
	;Hall E-field
	(oJxB[idx])['CATDESC'] = 'E = (JxB)/(ne)'
	(oJxB[idx])['COLOR']   = ['Blue', 'Forest Green', 'Red']
	(oJxB[idx])['LABEL']   = comps
	(oJxB[idx])['TITLE']   = 'JxB/ne!C(mV/m)'
	(oJxB[idx])['UNITS']   = 'mV/m'
	
	;Convective E-field
	oVexB['CATDESC'] = 'E = -VexB'
	oVexB['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oVexB['LABEL']   = comps
	oVexB['TITLE']   = '-VexB!C(mV/m)'
	oVexB['UNITS']   = 'mV/m'
	
	;Convective E-field
	oVixB['CATDESC'] = 'E = -VixB'
	oVixB['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oVixB['LABEL']   = comps
	oVixB['TITLE']   = '-VexB!C(mV/m)'
	oVixB['UNITS']   = 'mV/m'
	
;-------------------------------------------
; Split E-Fields into X, Y, Z Components ///
;-------------------------------------------
	oE[idx]   -> Split, oEx, oEy, oEz, /CACHE, NAME=[ex_vname, ey_vname, ez_vname]
	oE_perp   -> Split, oEx_perp, oEy_perp, oEz_perp, /CACHE, NAME=[eperpx_vname, eperpy_vname, eperpz_vname]
	oJxB[idx] -> Split, oJxBx, oJxBy, oJxBz, /CACHE
	oVexB     -> Split, oVexBx, oVexBy, oVexBz, /CACHE
	oVixB     -> Split, oVixBx, oVixBy, oVixBz, /CACHE
	
	;Axis range
	Exrange = [ Min( [oEx_perp.min, oVexBx.min, oVixBx.min, oJxBx.min] ), $
	            Max( [oEx_perp.max, oVexBx.max, oVixBx.max, oJxBx.max] ) ]
	Eyrange = [ Min( [oEy_perp.min, oVexBy.min, oVixBy.min, oJxBy.min] ), $
	            Max( [oEy_perp.max, oVexBy.max, oVixBy.max, oJxBy.max] ) ]
	Ezrange = [ Min( [oEz_perp.min, oVexBz.min, oVixBz.min, oJxBz.min] ), $
	            Max( [oEz_perp.max, oVexBz.max, oVixBz.max, oJxBz.max] ) ]
	
	;E
	oEx['AXIS_RANGE'] = Exrange
	oEx['LABEL']      = 'E'
	oEx['TITLE']      = 'E$\down' + comps[0] + '$!C(mV/m)'
	
	oEy['AXIS_RANGE'] = Eyrange
;	oEy['LABEL']      = 'E'
	oEy['TITLE']      = 'E$\down' + comps[1] + '$!C(mV/m)'
	
	oEz['AXIS_RANGE'] = Ezrange
;	oEz['LABEL']      = 'E'
	oEz['TITLE']      = 'E$\down' + comps[2] + '$!C(mV/m)'
	
	;Eperp
	oEx_perp['AXIS_RANGE'] = Exrange
	oEx_perp['LABEL']      = 'E$\down' + perp + '$'
	oEx_perp['TITLE']      = 'E$\down' + perp + comps[0] + '$!C(mV/m)'
	
	oEy_perp['AXIS_RANGE'] = Eyrange
	oEy_perp['TITLE']      = 'E$\down' + perp + comps[1] + '$!C(mV/m)'
	
	oEz_perp['AXIS_RANGE'] = Ezrange
	oEz_perp['TITLE']      = 'E$\down' + perp + comps[2] + '$!C(mV/m)'
	
	
	;JXB
	oJxBx['COLOR'] = 'Blue'
	oJxBx['LABEL'] = 'JxB/ne'
	
	oJxBy['COLOR'] = 'Blue'
;	oJxBy['LABEL'] = 'JxB/ne'
	
	oJxBz['COLOR'] = 'Blue'
;	oJxBz['LABEL'] = 'JxB/ne'
	
	;VIxB
	oVixBx['COLOR'] = 'Forest Green'
	oVixBx['LABEL'] = '-VixB'
	
	oVixBy['COLOR'] = 'Forest Green'
;	oVixBy['LABEL'] = '-VixB'
	
	oVixBz['COLOR'] = 'Forest Green'
;	oVixBz['LABEL'] = '-VixB'
	
	;VexB
	oVexBx['COLOR'] = 'Red'
	oVexBx['LABEL'] = '-VexB
	
	oVexBy['COLOR'] = 'Red'
;	oVexBy['LABEL'] = '-VexB
	
	oVexBz['COLOR'] = 'Red'
;	oVexBz['LABEL'] = '-VexB
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	oNi = MrVar_Get(ni_vnames[idx])
	oNe = MrVar_Get(ne_vnames[idx])
	oNi['AXIS_RANGE'] = [ Min([oNi.min, oNe.min]), Max([oNi.max, oNe.max]) ]
	oNi['COLOR']      = 'Blue'
	oNi['LABEL']      = 'Ni'

	oNe['COLOR'] = 'Red'
	oNe['LABEL'] = 'Ne'


	title = StrUpCase(StrJoin([sc, mode, level, 'Hall'], ' '))
	oEx_perp['PLOT_TITLE'] = title

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [oEx_perp, oEy_perp, oEz_perp, oNi, oBmag, oJ[idx], oRho], $ ;[ex_vname, ey_vname, ez_vname, ni_vname, bvec_vname, j_vname, rho_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	win = MrVar_OPlotTS( bmag_vnames[idx], oB[idx] )
	win = MrVar_OPlotTS( eperpx_vname, [ vixbx_vname, vexbx_vname, jxbx_vnames[idx] ] )
	win = MrVar_OPlotTS( eperpy_vname, [ vixby_vname, vexby_vname, jxby_vnames[idx] ] )
	win = MrVar_OPlotTS( eperpz_vname, [ vixbz_vname, vexbz_vname, jxbz_vnames[idx] ] )
	win = MrVar_OPlotTS( ni_vnames[idx], ne_vnames[idx] )
	win = MrVar_OPlotTS( rho_vname, div_jxb_vname )
	
	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 12]
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
		fname   = StrJoin( [sc, 'edp', mode, level, 'hall'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END