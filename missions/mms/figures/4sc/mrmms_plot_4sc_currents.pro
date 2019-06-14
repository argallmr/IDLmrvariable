; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_Currents
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
;       1. |B| MMS 1-4
;       2. Bx MMS 1-4
;       3. By MMS 1-4
;       4. Bz MMS 1-4
;       5. Jx MMS 1-4 + Curlometer
;       6. Jy MMS 1-4 + Curlometer
;       7. Jz MMS 1-4 + Curlometer
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
;       2018/09/25  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_Currents, mode, $
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
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
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
	nabla     = '!9'+String(71B)+'!X'
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
	
	;FPI
	fpi_mode  = mode
	IF fpi_mode EQ 'srvy' THEN BEGIN
		MrPrintF, 'LogWarn', 'FPI does not have "srvy" data. Using "fast".'
		fpi_mode = 'fast'
	ENDIF
	fpi_coords = (coords EQ 'dmpa' || coords EQ 'dsl') ? 'dbcs' : coords
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	sc  = 'mms' + ['1', '2', '3', '4']
	nSC = N_Elements(sc)
	
	;FGM
	b_vnames = sc + '_' + StrJoin([fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_')
	IF fgm_instr EQ 'fsm' THEN BEGIN
		bvec_vnames = sc + '_' + StrJoin([fgm_instr, 'b', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([fgm_instr, 'b', 'mag',      fgm_mode, fgm_level], '_')
	ENDIF ELSE BEGIN
		bvec_vnames = sc + '_' + StrJoin([fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_')
	ENDELSE
	
	;FPI
	ni_vnames = sc + '_' + StrJoin(['dis', 'numberdensity', fpi_mode], '_')
	ne_vnames = sc + '_' + StrJoin(['des', 'numberdensity', fpi_mode], '_')
	vi_vnames = sc + '_' + StrJoin(['dis', 'bulkv', fpi_coords, fpi_mode], '_')
	ve_vnames = sc + '_' + StrJoin(['des', 'bulkv', fpi_coords, fpi_mode], '_')
	
	;EPHEM
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Derived names
	bx_vnames    = bvec_vnames + '_x'
	by_vnames    = bvec_vnames + '_y'
	bz_vnames    = bvec_vnames + '_z'
	j_vnames     = sc + '_' + StrJoin(['fpi', 'j', fpi_coords, fpi_mode], '_')
	jx_vnames    = j_vnames + '_x'
	jy_vnames    = j_vnames + '_y'
	jz_vnames    = j_vnames + '_z'
	curlb_vname  = StrJoin(['mms', fgm_instr, 'curlb', fgm_coords, fgm_mode, fgm_level], '_')
	curlbx_vname = curlb_vname + '_x'
	curlby_vname = curlb_vname + '_y'
	curlbz_vname = curlb_vname + '_z'
	divb_vname  = StrJoin(['mms', fgm_instr, 'divb',  fgm_mode, fgm_level], '_')
		
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
		
		;FPI
		MrMMS_FPI_Load_Data, '', fpi_mode, $
		                     OPTDESC   = ['des-moms', 'dis-moms'], $
		                     VARFORMAT = '*' + ['numberdensity', 'bulkv_'+fpi_coords] + '_'+fpi_mode
		
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
; Split Vectors; Current Density ///////////
;-------------------------------------------
	;Interpolate DIS to DES
	oNe    = MrVar_Get(ne_vnames[0])
	oB     = MrVar_Get(b_vnames[0])
	oT_fpi = oNe['DEPEND_0']
	oT_fgm = oB['DEPEND_0']
	
	;Ranges
	b_range  = [!values.f_infinity, -!values.f_infinity]
	bx_range = [!values.f_infinity, -!values.f_infinity]
	by_range = [!values.f_infinity, -!values.f_infinity]
	bz_range = [!values.f_infinity, -!values.f_infinity]
	jx_range = [!values.f_infinity, -!values.f_infinity]
	jy_range = [!values.f_infinity, -!values.f_infinity]
	jz_range = [!values.f_infinity, -!values.f_infinity]
	
	oaB  = ObjArr(nSC)
	oaBvec = MrVar_Get(bvec_vnames)
	oaBmag = MrVar_Get(bmag_vnames)
	oaBx = ObjArr(nSC)
	oaBy = ObjArr(nSC)
	oaBz = ObjArr(nSC)
	oaNe = ObjArr(nSC)
	oaNi = ObjArr(nSC)
	oaVe = ObjArr(nSC)
	oaVi = ObjArr(nSC)
	oaJ  = ObjArr(nSC)
	oaJx = ObjArr(nSC)
	oaJy = ObjArr(nSC)
	oaJz = ObjArr(nSC)
	oaR  = ObjArr(nSC)
	FOR i = 0, nSC - 1 DO BEGIN
		scstr = 'mms' + String(i+1, FORMAT='(i1)')
		
		;Extract variables
		b   = MrVar_Get(b_vnames[i])
		bvec = MrVar_Get(bvec_vnames[i])
		n_e = MrVar_Get(ne_vnames[i])
		n_i = MrVar_Get(ni_vnames[i])
		v_e = MrVar_Get(ve_vnames[i])
		v_i = MrVar_Get(vi_vnames[i])
		r   = MrVar_Get(r_vnames[i])
		
	;-------------------------------------------
	; Interpolate //////////////////////////////
	;-------------------------------------------
		oaB[i]  = b -> Interpol(oT_fgm)
		oaBvec[i] = bvec -> Interpol(oT_fgm)
		oaNi[i] = n_i -> Interpol(oT_fpi)
		oaNe[i] = n_e -> Interpol(oT_fpi)
		oaVi[i] = v_i -> Interpol(oT_fpi)
		oaVe[i] = v_e -> Interpol(oT_fpi)
		oaR[i]  = r   -> Interpol(oT_fgm)
		
	;-------------------------------------------
	; Current Density //////////////////////////
	;-------------------------------------------
		oaJ[i] = 1e15 * q * oaNe[i] * (oaVi[i] - oaVe[i])
		jtemp  = (oaJ[i])['DATA']
		FOR j = 0, 2 DO jtemp[0,j] = Convol(jtemp[*,j], [0.25, 0.5, 0.25])
		oaJ[i] -> SetData, Temporary(jtemp)
		oaJ[i] -> SetName, j_vnames[i]
		oaJ[i] -> Cache
		oaJ[i] -> AddAttr, 'CATDESC', 'Current density computed from particle moments.'
		oaJ[i] -> AddAttr, 'TITLE',   'J!C($\mu$A/m$\up2$'
		oaJ[i] -> AddAttr, 'UNITS',   'uA/m^2'
		
	;-------------------------------------------
	; Split B //////////////////////////////////
	;-------------------------------------------
		;Split vector
		oaBvec[i] -> Split, oBx, oBy, oBz, /CACHE
		
		;Ranges
		bx_range[0] <= oBx.min
		bx_range[1] >= oBx.max
		by_range[0] <= oBy.min
		by_range[1] >= oBy.max
		bz_range[0] <= oBz.min
		bz_range[1] >= oBz.max
		
		;Set Attributes
		oBx -> SetName, bx_vnames[i]
		oBx['CATDESC'] = 'X-component of the magnetic field from ' + fgm_instr + '.'
		oBx['COLOR']   = sc_colors[i]
		oBx['TITLE']   = 'Bx!C(nT)'
		oBx['UNITS']   = 'nT'
		
		oBy -> SetName, by_vnames[i]
		oBy['CATDESC'] = 'Y-component of the magnetic field from ' + fgm_instr + '.'
		oBy['COLOR']   = sc_colors[i]
		oBy['TITLE']   = 'By!C(nT)'
		oBy['UNITS']   = 'nT'
		
		oBz -> SetName, bz_vnames[i]
		oBz['CATDESC'] = 'Z-component of the magnetic field from ' + fgm_instr + '.'
		oBz['COLOR']   = sc_colors[i]
		oBz['TITLE']   = 'Bz!C(nT)'
		oBz['UNITS']   = 'nT'
		
		;|B|
		oBmag = oaBmag[i]
		b_range[0]    <= oBmag.min
		b_range[1]    >= oBmag.max
		oBmag['COLOR'] = sc_colors[i]
		oBmag['LABEL'] = scstr
		
		;Store components
		oaBx[i]   = oBx
		oaBy[i]   = oBy
		oaBz[i]   = oBz
		
	;-------------------------------------------
	; Split J //////////////////////////////////
	;-------------------------------------------
		;Split vector
		oaJ[i] -> Split, oJx, oJy, oJz, /CACHE
		
		;Ranges
		jx_range[0] <= oJx.min
		jx_range[1] >= oJx.max
		jy_range[0] <= oJy.min
		jy_range[1] >= oJy.max
		jz_range[0] <= oJz.min
		jz_range[1] >= oJz.max
		
		;Set Attributes
		oJx['CATDESC'] = 'X-component of the current density.'
		oJx['COLOR']   = sc_colors[i]
		oJx['TITLE']   = 'Jx!C($\mu$A/m$\up2$)'
		oJx['UNITS']   = 'uA/m^2'
		
		oJy['CATDESC'] = 'X-component of the current density.'
		oJy['COLOR']   = sc_colors[i]
		oJy['TITLE']   = 'Jy!C($\mu$A/m$\up2$)'
		oJy['UNITS']   = 'uA/m^2'
		
		oJz['CATDESC'] = 'Z-component of the current density.'
		oJz['COLOR']   = sc_colors[i]
		oJz['TITLE']   = 'Jz!C($\mu$A/m$\up2$)'
		oJz['UNITS']   = 'uA/m^2'
		
		;Store components
		oaJx[i] = oJx
		oaJy[i] = oJy
		oaJz[i] = oJz
	ENDFOR
		
;-------------------------------------------
; Current Density //////////////////////////
;-------------------------------------------
	
	;Reciprocal vectors for taking derivatives
	oRecipVec = MrVar_RecipVec(oaR[0], oaR[1], oaR[2], oaR[3])
	
	;Curlometer Current density
	oCurlB = (1e-6/mu0) * oRecipVec -> Curl( oaBvec[0], oaBvec[1], oaBvec[2], oaBvec[3] )
	oCurlB -> SetName, curlb_vname
	oCurlB -> Cache
	oCurlB['CATDESC'] = 'Current density derived from the curlometer technique using MEC & ' + fgm_instr + '.'
	oCurlB['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oCurlB['LABEL']   = 'J$\down' + ['X', 'Y', 'Z'] + '$'
	oCurlB['TITLE']   = nabla + 'xB/$\mu$$\down0$!C($\mu$A/m^2)'
	oCurlB['UNITS']   = 'uA/m^2'
	
	;Split into components
	oCurlB -> Split, oCurlBx, oCurlBy, oCurlBz, /CACHE
		
	;Set Attributes
	oCurlBx['CATDESC'] = 'X-component of the current density computed from the curlometer technique.'
	oCurlBx['COLOR']   = 'Magenta'
	oCurlBx['LABEL']   = nabla + 'xB/$\mu$$\down0$'
	oCurlBx['TITLE']   = 'Jx!C($\mu$A/m$\up2$)'
	oCurlBx['UNITS']   = 'uA/m^2'
	
	oCurlBy['CATDESC'] = 'Y-component of the current density computed from the curlometer technique.' + fgm_instr + '.'
	oCurlBy['COLOR']   = 'Magenta'
	oCurlBy['TITLE']   = 'Jy!C($\mu$A/m$\up2$)'
	oCurlBy['UNITS']   = 'uA/m^2'
	
	oCurlBz['CATDESC'] = 'Z-component of the current density computed from the curlometer technique.' + fgm_instr + '.'
	oCurlBz['COLOR']   = 'Magenta'
	oCurlBz['TITLE']   = 'Jz!C($\mu$A/m$\up2$)'
	oCurlBz['UNITS']   = 'uA/m^2'
	
	;Divergence of B
	oDivB = (1e-6/mu0) * oRecipVec -> Divergence( oaBvec[0], oaBvec[1], oaBvec[2], oaBvec[3] )
	oDivB -> SetName, divb_vname
	oDivB -> Cache
	oDivB['AXIS_RANGE'] = Replicate( Max( Abs([oDivB.min, oDivB.max]) ), 2) * [-1,1]
	oDivB['CATDESC']    = 'Divergence of B derived from the curlometer technique using MEC & ' + fgm_instr + '.'
	oDivB['TITLE']      = nabla + '.B/$\mu$$\down0$!C($\mu$A/m^2)'
	oDivB['UNITS']      = 'uA/m^2'
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	(oaBmag[0])['PLOT_TITLE'] = 'MMS1-4 Currents'
	(oaBmag[0])['AXIS_RANGE'] = b_range
	(oaBx[0])['AXIS_RANGE'] = bx_range
	(oaBy[0])['AXIS_RANGE'] = by_range
	(oaBz[0])['AXIS_RANGE'] = bz_range
	(oaJx[0])['AXIS_RANGE'] = jx_range
	(oaJy[0])['AXIS_RANGE'] = jy_range
	(oaJz[0])['AXIS_RANGE'] = jz_range
	
;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot MMS1 data
	win = MrVar_PlotTS( [ bmag_vnames[0], bx_vnames[0], by_vnames[0], bz_vnames[0], $
	                      jx_vnames[0], jy_vnames[0], jz_vnames[0], divb_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	
	;Overplot MMS 2-4 data
	win = MrVar_OPlotTS( bmag_vnames[0], bmag_vnames[1:3] )
	win = MrVar_OPlotTS( bx_vnames[0],   bx_vnames[1:3] )
	win = MrVar_OPlotTS( by_vnames[0],   by_vnames[1:3] )
	win = MrVar_OPlotTS( bz_vnames[0],   bz_vnames[1:3] )
	win = MrVar_OPlotTS( jx_vnames[0],   [jx_vnames[1:3], curlbx_vname] )
	win = MrVar_OPlotTS( jy_vnames[0],   [jy_vnames[1:3], curlby_vname] )
	win = MrVar_OPlotTS( jz_vnames[0],   [jz_vnames[1:3], curlbz_vname] )
		
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
		fname   = StrJoin( ['mms', fgm_instr, mode, level, '4sc-currents'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END