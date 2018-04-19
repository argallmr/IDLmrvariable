; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_EDP
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
;       2018/02/06  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_EDP, mode, $
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
	IF N_Elements(level)  EQ 0 THEN level = 'l2'
	IF N_Elements(trange) GT 0 THEN MrVar_SetTRange, trange
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	
	e0        = MrConstants('epsilon_0')
	sc_colors = ['Black', 'Red', 'Green', 'Blue']
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	instr = 'edp'
	IF N_Elements(coords) EQ 0 THEN BEGIN
		CASE level OF
			'ql': coords = 'dsl'
			ELSE: coords = 'gse'
		ENDCASE
	ENDIF

	;Source names
	sc         = 'mms' + ['1', '2', '3', '4']
	e_vnames   = sc + '_' + StrJoin( [instr, 'dce',   coords, mode, level], '_' )
	v_vnames   = sc + '_' + StrJoin( [instr, 'scpot',         mode, level], '_' )
	epe_vnames = sc + '_' + StrJoin( [instr, 'dce',   'par',  'epar', mode, level], '_' )
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;MEC
	IF mode EQ 'brst' THEN BEGIN
		IF (mrmms_get_filenames('mms1', 'mec', mode, level, OPTDESC=ephdesc))[0] EQ '' $
			THEN mec_mode = 'srvy' $
			ELSE mec_mode = mode
	ENDIF
	
	;Output names
	ex_vnames    = sc + '_' + StrJoin( ['dce', coords, mode, level, 'x'], '_' )
	ey_vnames    = sc + '_' + StrJoin( ['dce', coords, mode, level, 'y'], '_' )
	ez_vnames    = sc + '_' + StrJoin( ['dce', coords, mode, level, 'z'], '_' )
	epar_vnames  = sc + '_' + StrJoin( ['dce', 'par',  mode, level], '_' )
	gradv_vname  = StrJoin( ['mms', 'dce', 'gradv', mode, level], '_' )
	gradvx_vname = StrJoin( ['mms', 'dce', 'gradv', mode, level, 'x'], '_' )
	gradvy_vname = StrJoin( ['mms', 'dce', 'gradv', mode, level, 'y'], '_' )
	gradvz_vname = StrJoin( ['mms', 'dce', 'gradv', mode, level, 'z'], '_' )
	rho_vname    = StrJoin( ['mms', 'dce', 'rho',   mode, level], '_' )
	curlE_vname  = StrJoin( ['mms', 'dce', 'curle', mode, level], '_' )
		
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;E-Field
		MrMMS_Load_Data, '', instr, mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = ['*dce_'+coords+'*', '*dce_par*']
		
		;Spacecraft Potential
		MrMMS_Load_Data, '', instr, mode, level, $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
		
		;Ephemeris
		MrMMS_Load_Data, '', 'mec', mec_mode, level, $
		                 OPTDESC   = ephdesc, $
		                 VARFORMAT = r_vnames
	ENDIF

;-------------------------------------------
; Interpolate //////////////////////////////
;-------------------------------------------
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
; Curl and Divergence of E /////////////////
;-------------------------------------------
	oRecipVec = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	
	;Charge Density
	;   - 1e4 converts to C/cm^2
	oRho = e0 * 1e4 * oRecipVec -> Divergence( oE1, oE2, oE3, oE4 )
	oRho -> SetName, rho_vname
	oRho -> Cache
	oRho['CATDESC'] = 'Free charge density, computed from the divergence of the electric field.'
	oRho['TITLE']   = '$\rho$!C(C/cm^2)'
	oRho['UNITS']   = 'C/cm^2'
	
	;Curl of E
	;   - 1e4 converts to C/cm^2
	oCurl = e0 * 1e1 * oRecipVec -> Curl( oE1, oE2, oE3, oE4 )
	oCurl -> SetName, curlE_vname
	oCurl -> Cache
	oCurl['CATDESC'] = 'Curl(E)'
	oCurl['COLOR']   = ['Blue', 'Green', 'Red']
	oCurl['LABEL']   = ['X', 'Y', 'Z']
	oCurl['TITLE']   = 'Curl(E)'
	oCurl['UNITS']   = ''
	
;-------------------------------------------
; -Gradient of V ///////////////////////////
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
	oGradVx['LINESTYLE'] = '--'
	oGradVy['COLOR']     = 'Magenta'
	oGradVy['LINESTYLE'] = '--'
	oGradVz['COLOR']     = 'Magenta'
	oGradVz['LINESTYLE'] = '--'

;-------------------------------------------
; Extract Components ///////////////////////
;-------------------------------------------

	;Split into components
	xrange   = [!values.f_infinity, -!values.f_infinity]
	yrange   = [!values.f_infinity, -!values.f_infinity]
	zrange   = [!values.f_infinity, -!values.f_infinity]
	parrange = [!values.f_infinity, -!values.f_infinity]
	vrange   = [!values.f_infinity, -!values.f_infinity]
	FOR i = 0, 3 DO BEGIN
		oV  = MrVar_Get(v_vnames[i])
		oE  = MrVar_Get(e_vnames[i])
		oE -> Split, oEx, oEy, oEz, /CACHE
		
		;Epar
		;   - Remove the error term
		oEpar = MrVar_Get(epe_vnames[i])
		oEpar = oEpar[*,1]
		oEpar -> RemoveAttr, 'LABEL'
		oEpar -> SetName, epar_vnames[i]
		oEpar -> Cache
		
		;Spacecraft color
		CASE i OF
			0: color = 'Black'
			1: color = 'Red'
			2: color = 'Green'
			3: color = 'Blue'
			ELSE: Message, 'Invalid spacecraft number: ' + String(i, FORMAT='(i0)')
		ENDCASE
		
		;E_||
		oEpar['COLOR'] = color
		parrange[0]   <= oEpar.min
		parrange[1]   >= oEpar.max
		
		;V
		oV['COLOR'] = color
		vrange[0]   <= oV.min
		vrange[1]   >= oV.max
		
		
		;Ex
		oEx -> SetName, ex_vnames[i]
		oEx['COLOR'] = color
		xrange[0]   <= oEx.min
		xrange[1]   >= oEx.max
		
		;Ey
		oEy -> SetName, ey_vnames[i]
		oEy['COLOR'] = color
		yrange[0]   <= oEy.min
		yrange[1]   >= oEy.max
		
		;Ez
		oEz -> SetName, ez_vnames[i]
		oEz['COLOR'] = color
		zrange[0]   <= oEz.min
		zrange[1]   >= oEz.max
	ENDFOR
	
	;Clamp at +/-100mV/m
	parrange[0] >= -100
	parrange[1] <= 100
	vrange[0]   >= -100
	vrange[1]   <= 100
	xrange[0]   >= -100
	xrange[1]   <= 100
	yrange[0]   >= -100
	yrange[1]   <= 100
	zrange[0]   >= -100
	zrange[1]   <= 100
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	title = StrUpCase(StrJoin(['mms1234', instr, mode, level], ' '))

	oEx1 = MrVar_Get(ex_vnames[0])
	oEx2 = MrVar_Get(ex_vnames[1])
	oEx3 = MrVar_Get(ex_vnames[2])
	oEx4 = MrVar_Get(ex_vnames[3])
	odVx = MrVar_Get(gradvx_vname)
	oEx1['PLOT_TITLE'] = title
	oEx1['TITLE'] = 'Ex!C(mV/m)'
	oEx1['LABEL'] = 'mms1'
	oEx2['LABEL'] = 'mms2'
	oEx3['LABEL'] = 'mms3'
	oEx4['LABEL'] = 'mms4'
	odVx['LABEL'] = '-Grad(V)'
	
	oEy = MrVar_Get(ey_vnames[0])
	oEy['TITLE'] = 'Ey!C(mV/m)'
	
	oEz = MrVar_Get(ez_vnames[0])
	oEz['TITLE'] = 'Ez!C(mV/m)'
	
	oEpar = MrVar_Get(epar_vnames[0])
	oEpar['TITLE'] = 'E$\down||$!C(mV/m)'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [ ex_vnames[0], ey_vnames[0], ez_vnames[0], epar_vnames[0], v_vnames[0], rho_vname, curlE_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	win = MrVar_OPlotTS( v_vnames[0],    v_vnames[1:3] )
	win = MrVar_OPlotTS( ex_vnames[0],   [ex_vnames[1:3], gradvx_vname] )
	win = MrVar_OPlotTS( ey_vnames[0],   [ey_vnames[1:3], gradvy_vname] )
	win = MrVar_OPlotTS( ez_vnames[0],   [ez_vnames[1:3], gradvz_vname] )
	win = MrVar_OPlotTS( epar_vnames[0], epar_vnames[1:3] )
	
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
		fname   = StrJoin( ['mms1234', instr, mode, level, '4sc'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END