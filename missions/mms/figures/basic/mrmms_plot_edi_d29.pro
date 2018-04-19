; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_EDI_Q0
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
;   Generate a plot of EDP quantities:
;       1. FGM Bxyz
;       2. EDP Exyz
;       3. GDU Theta
;       4. GDU Phi
;       5. FAC Phi
;       6. FAC Theta
;       7. PAD
;       8. GPD
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
;       FGM_INSTR:  in, optional, type=string, default='fgm'
;                   FGM instrument to use. Options are: { 'afg' | 'dfg' | 'fgm' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'ql' | 'l2'}
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
;       2017/01/22  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_EDI_D29, sc, mode, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
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
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'dfg'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(mode)      EQ 0 THEN mode      = 'srvy'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	instr      = 'edi'
	edi_level  = 'l1a'
	optdesc    = 'efield'
	coords     = 'dbcs'
	
	;FGM 
	CASE fgm_instr OF
		'afg': fgm_level = 'l2pre'
		'dfg': fgm_level = 'l2pre'
		'fgm': fgm_level = 'l2'
		ELSE: Message, 'Invalid FGM instrument: "' + fgm_instr + '".'
	ENDCASE
	fgm_mode   = mode EQ 'brst' ? mode : 'srvy'
	fgm_coords = coords EQ 'dbcs'  ? 'dmpa' : coords
	
	;EDP
	edp_instr  = 'edp'
	edp_mode   = mode eq 'brst' ? mode : 'fast'
	edp_coords = coords EQ 'dbcs'  ? 'dsl'  : coords

	;Source names
	fgm_b_vname       = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_' )
	fgm_bvec_vname    = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_' )
	fgm_bmag_vname    = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_' )
	e_vname           = StrJoin( [sc, edp_instr, 'dce',  edp_coords, edp_mode, level], '_' )
	d29_gd12_vname    = StrJoin( [sc, 'edi', 'data29', 'gd12'], '_')
	d29_gd21_vname    = StrJoin( [sc, 'edi', 'data29', 'gd21'], '_')
	phi_gd12_vname    = StrJoin( [sc, 'edi', 'phi',    'gd12'], '_')
	phi_gd21_vname    = StrJoin( [sc, 'edi', 'phi',    'gd21'], '_')
	theta_gd12_vname  = StrJoin( [sc, 'edi', 'theta',  'gd12'], '_')
	theta_gd21_vname  = StrJoin( [sc, 'edi', 'theta',  'gd21'], '_')
	q_gd12_vname      = StrJoin( [sc, 'edi', 'sq',     'gd12'], '_')
	q_gd21_vname      = StrJoin( [sc, 'edi', 'sq',     'gd21'], '_')
	gdu1_abscal_vname = StrJoin( [sc, 'edi', 'abs',    'gdu1', 'cal', level], '_')
	gdu2_abscal_vname = StrJoin( [sc, 'edi', 'abs',    'gdu2', 'cal', level], '_')

	;Derived names
	q0_gd12_vname        = StrJoin( [sc, 'edi', 'q0',     'gd12'], '_')
	q1_gd12_vname        = StrJoin( [sc, 'edi', 'q1',     'gd12'], '_')
	q2_gd12_vname        = StrJoin( [sc, 'edi', 'q2',     'gd12'], '_')
	q3_gd12_vname        = StrJoin( [sc, 'edi', 'q3',     'gd12'], '_')
	q0_gd21_vname        = StrJoin( [sc, 'edi', 'q0',     'gd21'], '_')
	q1_gd21_vname        = StrJoin( [sc, 'edi', 'q1',     'gd21'], '_')
	q2_gd21_vname        = StrJoin( [sc, 'edi', 'q2',     'gd21'], '_')
	q3_gd21_vname        = StrJoin( [sc, 'edi', 'q3',     'gd21'], '_')
	phi_gdu1_vname       = StrJoin( [sc, instr, 'phi',   'data29', 'gd21', mode, level], '_' )
	phi_gdu2_vname       = StrJoin( [sc, instr, 'phi',   'data29', 'gd12', mode, level], '_' )
	theta_gdu1_vname     = StrJoin( [sc, instr, 'theta', 'data29', 'gd21', mode, level], '_' )
	theta_gdu2_vname     = StrJoin( [sc, instr, 'theta', 'data29', 'gd12', mode, level], '_' )
	ang_diff_vname       = StrJoin( [sc, instr, 'dlook', 'data29',         mode, level], '_' )

;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = fgm_level, $
		                     VARFORMAT = [ '*'+fgm_instr+'_'+fgm_mode+'*'+fgm_coords, $
		                                   '*'+fgm_instr+'_b_*']
		
		;Check which naming convention is in use.
		IF ~MrVar_IsCached(fgm_b_vname) THEN BEGIN
			;Old names
			fgm_b_vname     = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_' )
			fgm_bvec_vname  = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_' )
			fgm_bmag_vname  = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_' )
			
			;Check if old names exist
			IF ~MrVar_IsCached(fgm_b_vname) THEN Message, 'FGM variable names unknown.'
		ENDIF
			
			

		;EDP
		MrMMS_Load_Data, sc, edp_instr, edp_mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = e_vname

		;EDI
		MrMMS_Load_Data, sc, 'edi', mode, edi_level, $
		                 OPTDESC   = optdesc, $
		                 VARFORMAT = [ '*phi*', '*theta*', '*data29*', '*sq*' ]
		
		;EDI Calibrations
		MrMMS_EDI_Load_Cals, sc, $
		                     VARFORMAT='*abs*cal*'
	ENDIF

;-------------------------------------------
; Interpolate Angles to D29 ////////////////
;-------------------------------------------
	;Look-directions in instrument coordinates
	oPhi_GD12   = MrVar_Get( phi_gd12_vname   )
	oPhi_GD21   = MrVar_Get( phi_gd21_vname   )
	oTheta_GD12 = MrVar_Get( theta_gd12_vname )
	oTheta_GD21 = MrVar_Get( theta_gd21_vname )
	oD29_GD12   = MrVar_Get( d29_gd12_vname   )
	oD29_GD21   = MrVar_Get( d29_gd21_vname   )

	;Interpolate angles to DATA29 cadence
	;   - Data29 GD12 & GD21 have the same time tags
	;   - Angles for GD12 & GD21 do not have the same time tags
	oPhi_GDU1   = oPhi_GD21   -> Interpol(oD29_GD21)
	oPhi_GDU2   = oPhi_GD12   -> Interpol(oD29_GD12)
	oTheta_GDU1 = oTheta_GD21 -> Interpol(oD29_GD21)
	oTheta_GDU2 = oTheta_GD12 -> Interpol(oD29_GD12)
	
	;GDU1: Do not extrapolate!
	t_gdu1 = oPhi_GDU1['TIME', 'TT2000']
	t_gd21 = oPhi_GD21['TIME', 'TT2000']
	iLT = Where(t_gdu1 LT t_gd21[0], nLT)
	iGT = Where(t_gdu1 GT t_gd21[-1], nGT)
	IF nLT GT 0 THEN BEGIN
		oPhi_GDU1[iLT]   = !Values.F_NaN
		oTheta_GDU1[iLT] = !Values.F_NaN
	ENDIF
	IF nGT GT 0 THEN BEGIN
		oPhi_GDU1[iGT]   = !Values.F_NaN
		oTheta_GDU1[iGT] = !Values.F_NaN
	ENDIF
	
	;GDU2: Do not extrapolate!
	
	t_gdu2 = oPhi_GDU1['TIME', 'TT2000']
	t_gd12 = oPhi_GD12['TIME', 'TT2000']
	iLT = Where(t_gdu2 LT t_gd12[0], nLT)
	iGT = Where(t_gdu2 GT t_gd12[-1], nGT)
	IF nLT GT 0 THEN BEGIN
		oPhi_GDU2[iLT]   = !Values.F_NaN
		oTheta_GDU2[iLT] = !Values.F_NaN
	ENDIF
	IF nGT GT 0 THEN BEGIN
		oPhi_GDU2[iGT]   = !Values.F_NaN
		oTheta_GDU2[iGT] = !Values.F_NaN
	ENDIF
	
	;Name the variables
	oPhi_GDU1   -> SetName, phi_gdu1_vname
	oPhi_GDU2   -> SetName, phi_gdu2_vname
	oTheta_GDU1 -> SetName, theta_gdu1_vname
	oTheta_GDU2 -> SetName, theta_gdu2_vname
	
	;Stash them away
	oPhi_GDU1   -> Cache
	oPhi_GDU2   -> Cache
	oTheta_GDU1 -> Cache
	oTheta_GDU2 -> Cache

;-------------------------------------------
; Angular Difference ///////////////////////
;-------------------------------------------
	oPhi_GDU1   *= !dtor
	oPhi_GDU2   *= !dtor
	oTheta_GDU1 *= !dtor
	oTheta_GDU2 *= !dtor

	;Convert to cartesian coords
	x1 = Sin( oTheta_GDU1['DATA'] ) * Cos( oPhi_GDU1['DATA'] )
	y1 = Sin( oTheta_GDU1['DATA'] ) * Sin( oPhi_GDU1['DATA'] )
	z1 = Cos( oTheta_GDU1['DATA'] )
	
	;Angles are in GDU coordinates, so one vector must be negated
	;   - X-Axes are aligned. Others are inverted
	x2 =   Sin( oTheta_GDU2['DATA'] ) * Cos( oPhi_GDU2['DATA'] )
	y2 = - Sin( oTheta_GDU2['DATA'] ) * Sin( oPhi_GDU2['DATA'] )
	z2 = - Cos( oTheta_GDU2['DATA'] )

	;Angular difference between GDU1 & GDU2
	angdiff = ACos( x1*x2 + y1*y2 + z1*z2 ) * !radeg
	oDAng   = MrScalarTS( oPhi_GDU1['TIMEVAR'], angdiff, $
	                      /CACHE, $
	                      NAME   = ang_diff_vname )

;-------------------------------------------
; Split Quality by Flag ////////////////////
;-------------------------------------------
	;Grab the data
	oQ_gd21 = MrVar_Get(q_gd21_vname)
	oQ_gd12 = MrVar_Get(q_gd12_vname)
	oT_gd12 = oQ_gd12['TIMEVAR']
	oT_gd21 = oQ_gd21['TIMEVAR']
	
	;
	; GDU1
	;
	iq0 = oQ_gd12 -> Where(0, /EQUAL, COUNT=nq0)
	iq1 = oQ_gd12 -> Where(1, /EQUAL, COUNT=nq1)
	iq2 = oQ_gd12 -> Where(2, /EQUAL, COUNT=nq2)
	iq3 = oQ_gd12 -> Where(3, /EQUAL, COUNT=nq3)
	
;	IF nq0 GT 0 THEN oQ0_gd12 = MrScalarTS( oT_gd12[iq0], oQ_gd12[iq0], /CACHE, NAME=q0_gd12_vname )
;	IF nq1 GT 0 THEN oQ1_gd12 = MrScalarTS( oT_gd12[iq1], oQ_gd12[iq1], /CACHE, NAME=q1_gd12_vname )
;	IF nq2 GT 0 THEN oQ2_gd12 = MrScalarTS( oT_gd12[iq2], oQ_gd12[iq2], /CACHE, NAME=q2_gd12_vname )
;	IF nq3 GT 0 THEN oQ3_gd12 = MrScalarTS( oT_gd12[iq3], oQ_gd12[iq3], /CACHE, NAME=q3_gd12_vname )

	trange = MrVar_GetTRange()
	oQ0_gd12 = nq0 EQ 0 ? MrScalarTS( trange, Replicate(!Values.F_NaN, 2), /CACHE, NAME=q0_gd12_vname ) $
	                    : MrScalarTS( oT_gd12[iq0], oQ_gd12[iq0],          /CACHE, NAME=q0_gd12_vname )
	oQ1_gd12 = nq1 EQ 0 ? MrScalarTS( trange, Replicate(!Values.F_NaN, 2), /CACHE, NAME=q1_gd12_vname ) $
	                    : MrScalarTS( oT_gd12[iq1], oQ_gd12[iq1],          /CACHE, NAME=q1_gd12_vname )
	oQ2_gd12 = nq2 EQ 0 ? MrScalarTS( trange, Replicate(!Values.F_NaN, 2), /CACHE, NAME=q2_gd12_vname ) $
	                    : MrScalarTS( oT_gd12[iq2], oQ_gd12[iq2],          /CACHE, NAME=q2_gd12_vname )
	oQ3_gd12 = nq3 EQ 0 ? MrScalarTS( trange, Replicate(!Values.F_NaN, 2), /CACHE, NAME=q3_gd12_vname ) $
	                    : MrScalarTS( oT_gd12[iq3], oQ_gd12[iq3],          /CACHE, NAME=q3_gd12_vname )

	;
	; GDU2
	;
	iq0 = oQ_gd21 -> Where(0, /EQUAL, COUNT=nq0)
	iq1 = oQ_gd21 -> Where(1, /EQUAL, COUNT=nq1)
	iq2 = oQ_gd21 -> Where(2, /EQUAL, COUNT=nq2)
	iq3 = oQ_gd21 -> Where(3, /EQUAL, COUNT=nq3)
	
	oQ0_gd21 = nq0 EQ 0 ? MrScalarTS( trange, Replicate(!Values.F_NaN, 2), /CACHE, NAME=q0_gd21_vname ) $
	                    : MrScalarTS( oT_gd21[iq0], oQ_gd21[iq0],          /CACHE, NAME=q0_gd21_vname )
	oQ1_gd21 = nq1 EQ 0 ? MrScalarTS( trange, Replicate(!Values.F_NaN, 2), /CACHE, NAME=q1_gd21_vname ) $
	                    : MrScalarTS( oT_gd21[iq1], oQ_gd21[iq1],          /CACHE, NAME=q1_gd21_vname )
	oQ2_gd21 = nq2 EQ 0 ? MrScalarTS( trange, Replicate(!Values.F_NaN, 2), /CACHE, NAME=q2_gd21_vname ) $
	                    : MrScalarTS( oT_gd21[iq2], oQ_gd21[iq2],          /CACHE, NAME=q2_gd21_vname )
	oQ3_gd21 = nq3 EQ 0 ? MrScalarTS( trange, Replicate(!Values.F_NaN, 2), /CACHE, NAME=q3_gd21_vname ) $
	                    : MrScalarTS( oT_gd21[iq3], oQ_gd21[iq3],          /CACHE, NAME=q3_gd21_vname )

;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	;BMAG
	oB = MrVar_Get(fgm_b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, 'DATA29'], ' ' ) )
	
	;E FIELD
	oE = MrVar_Get(e_vname)
	oE['TITLE'] = 'E!C(mV/m)'
	
	;ANGULAR DIFFERENCE
	oDAng['AXIS_RANGE']   = [0,180]
	oDAng['TICKINTERVAL'] = 90
	oDAng['TITLE']        = '$\Delta$$\Theta$!C(deg)'
	oDAng['UNITS']        = 'degrees'
	
	;
	; EFIELD - THETA & PHI
	;
	
	;GDU2
	oPhi_GD12['AXIS_RANGE']   = [-180,180]
	oPhi_GD12['COLOR']        = 'Red'
	oPhi_GD12['LABEL']        = 'GDU2'
	oPhi_GD12['TICKINTERVAL'] = 90.0
	oPhi_GD12['TITLE']        = 'Phi!C(deg)'
	oPhi_GD12['UNITS']        = 'degrees'
	
	;GDU1
	oPhi_GD21['AXIS_RANGE']   = [-180,180]
	oPhi_GD21['COLOR']        = 'Blue'
	oPhi_GD21['LABEL']        = 'GDU1'
	oPhi_GD21['TICKINTERVAL'] = 90.0
	oPhi_GD21['TITLE']        = 'Phi!C(deg)'
	oPhi_GD21['UNITS']        = 'degrees'
	
	;GDU2
	oTheta_GD12['AXIS_RANGE']   = [0,180]
	oTheta_GD12['COLOR']        = 'Red'
	oTheta_GD12['LABEL']        = 'GDU2'
	oTheta_GD12['TICKINTERVAL'] = 90.0
	oTheta_GD12['TITLE']        = 'Theta!C(deg)'
	oTheta_GD12['UNITS']        = 'degrees'
	
	;GDU1
	oTheta_GD21['AXIS_RANGE']   = [0,180]
	oTheta_GD21['COLOR']        = 'Blue'
	oTheta_GD21['LABEL']        = 'GDU1'
	oTheta_GD21['TICKINTERVAL'] = 90.0
	oTheta_GD21['TITLE']        = 'Theta!C(deg)'
	oTheta_GD21['UNITS']        = 'degrees'
	
	;
	; EFIELD - Quality
	;
	
	;GDU2
	oQ0_GD12['AXIS_RANGE'] = [-1,4]
	oQ0_GD12['COLOR']      = 'Purple'
	oQ0_GD12['LINESTYLE']  = 'None'
	oQ0_GD12['SYMBOL']     = 9
	oQ0_GD12['TITLE']      = 'Quality'
	
	oQ1_GD12['AXIS_RANGE'] = [-1,4]
	oQ1_GD12['COLOR']      = 'Blue'
	oQ1_GD12['LINESTYLE']  = 'None'
	oQ1_GD12['SYMBOL']     = 9
	oQ1_GD12['TITLE']      = 'Quality'

	oQ2_GD12['AXIS_RANGE'] = [-1,4]
	oQ2_GD12['COLOR']      = 'Yellow'
	oQ2_GD12['LINESTYLE']  = 'None'
	oQ2_GD12['SYMBOL']     = 9
	oQ2_GD12['TITLE']      = 'Quality'
	
	oQ3_GD12['AXIS_RANGE'] = [-1,4]
	oQ3_GD12['COLOR']      = 'Red'
	oQ3_GD12['LINESTYLE']  = 'None'
	oQ3_GD12['SYMBOL']     = 9
	oQ3_GD12['TITLE']      = 'Quality'
	
	;GDU1
	oQ0_GD21['AXIS_RANGE'] = [-1,4]
	oQ0_GD21['COLOR']      = 'Purple'
	oQ0_GD21['LINESTYLE']  = 'None'
	oQ0_GD21['SYMBOL']     = 9
	oQ0_GD21['TITLE']      = 'Quality'
	
	oQ1_GD21['AXIS_RANGE'] = [-1,4]
	oQ1_GD21['COLOR']      = 'Blue'
	oQ1_GD21['LINESTYLE']  = 'None'
	oQ1_GD21['SYMBOL']     = 9
	oQ1_GD21['TITLE']      = 'Quality'

	oQ2_GD21['AXIS_RANGE'] = [-1,4]
	oQ2_GD21['COLOR']      = 'Yellow'
	oQ2_GD21['LINESTYLE']  = 'None'
	oQ2_GD21['SYMBOL']     = 9
	oQ2_GD21['TITLE']      = 'Quality'
	
	oQ3_GD21['AXIS_RANGE'] = [-1,4]
	oQ3_GD21['COLOR']      = 'Red'
	oQ3_GD21['LINESTYLE']  = 'None'
	oQ3_GD21['SYMBOL']     = 9
	oQ3_GD21['TITLE']      = 'Quality'
	
	;
	; D29 - Counts
	;
	
	;GDU1
	oD29_GD12 = MrVar_Get(d29_gd12_vname)
	oD29_GD12['TITLE'] = 'Counts!CGDU2'
	
	;GDU1
	oD29_GD21 = MrVar_Get(d29_gd21_vname)
	oD29_GD21['TITLE'] = 'Counts!CGDU1'
	
;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [ fgm_b_vname, e_vname, theta_gd21_vname, phi_gd21_vname, ang_diff_vname, $
	                      d29_gd21_vname, q0_gd21_vname, d29_gd12_vname, q0_gd12_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 680, $
	                    YSIZE = 700 )
	
	win = MrVar_OPlotTS( theta_gd21_vname,     theta_gd12_vname )
	win = MrVar_OPlotTS( phi_gd21_vname,       phi_gd12_vname )
	win = MrVar_OPlotTS( q0_gd21_vname,        [q1_gd21_vname, q2_gd21_vname, q3_gd21_vname] )
	win = MrVar_OPlotTS( q0_gd12_vname,        [q1_gd12_vname, q2_gd12_vname, q3_gd12_vname] )

	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 6]
	win    -> Refresh

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
		fname = StrJoin( [sc, instr, mode, level, 'data29'], '_' )
		fname = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, win
END