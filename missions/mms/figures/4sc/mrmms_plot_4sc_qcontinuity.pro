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
FUNCTION MrMMS_Plot_4sc_QContinuity, mode, $
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
	
	;EDP
	edp_instr  = 'edp'
	edp_coords = MrIsMember(['dmpa', 'dbcs'], coords) ? 'dsl' : coords
	edp_mode   = mode
	IF mode EQ 'srvy' THEN BEGIN
		MrPrintF, 'LogWarn', 'EDP does not have "srvy" data. Using "fast".'
		edp_mode = 'fast'
	ENDIF
	
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
	
	;EDP
	e_vnames = sc + '_' + StrJoin([edp_instr, 'dce',  edp_coords, edp_mode, level], '_' )
	
	;FPI
	ni_vnames = sc + '_' + StrJoin(['dis', 'numberdensity', fpi_mode], '_')
	ne_vnames = sc + '_' + StrJoin(['des', 'numberdensity', fpi_mode], '_')
	vi_vnames = sc + '_' + StrJoin(['dis', 'bulkv', fpi_coords, fpi_mode], '_')
	ve_vnames = sc + '_' + StrJoin(['des', 'bulkv', fpi_coords, fpi_mode], '_')
	j_vnames  = sc + '_' + StrJoin(['fpi', 'j',     fpi_coords, fpi_mode], '_')
	
	;EPHEM
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	b_bc_vname    = StrJoin( ['mms', fgm_instr, 'bbc', fgm_coords, fgm_mode, fgm_level], '_' )
	e_bc_vname    = StrJoin( ['mms', 'dce', 'ebc', edp_coords, edp_mode, level], '_' )
	j_bc_vname    = StrJoin( ['mms', 'fpi', 'jbc', fpi_coords, fpi_mode], '_' )
	curlb_vname   = StrJoin( ['mms', fgm_instr,   'j',   fgm_coords, fgm_mode, fgm_level], '_' )
	rho_vname     = StrJoin( ['mms', 'dce',     'rho', edp_mode, level], '_' )
	divj_vname    = StrJoin( ['mms', 'fpi', 'divj', fpi_coords, fpi_mode], '_' )
	drho_dt_vname = StrJoin( ['mms', 'dce', 'dRho_dt', edp_mode, level], '_' )
		
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
		
		;DES
		MrMMS_FPI_Load_Data, '', fpi_mode, $
		                     OPTDESC   = ['des-moms', 'dis-moms'], $
		                     VARFORMAT = '*' + ['numberdensity', 'bulkv_'+fpi_coords] + '_'+fpi_mode
				
		;EDP
		MrMMS_Load_Data, '', 'edp', edp_mode, edp_level, $
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
; Div Current Density //////////////////////
;-------------------------------------------
	;Interpolate DIS to DES
	oNe   = MrVar_Get(ne_vnames[0])
	oTref = oNe['DEPEND_0']
	
	oaB  = ObjArr(nSC)
	oaE  = ObjArr(nSC)
	oaNe = ObjArr(nSC)
	oaNi = ObjArr(nSC)
	oaVe = ObjArr(nSC)
	oaVi = ObjArr(nSC)
	oaJ  = ObjArr(nSC)
	oaR  = ObjArr(nSC)
	FOR i = 0, nSC - 1 DO BEGIN
		;Extract variables
		b   = MrVar_Get(bvec_vnames[i])
		e   = MrVar_Get(e_vnames[i])
		n_e = MrVar_Get(ne_vnames[i])
		n_i = MrVar_Get(ni_vnames[i])
		v_e = MrVar_Get(ve_vnames[i])
		v_i = MrVar_Get(vi_vnames[i])
		r   = MrVar_Get(r_vnames[i])
		
		;Interpolate
		oaB[i]  = b   -> Interpol(oTref)
		oaE[i]  = e   -> Interpol(oTref)
		oaNi[i] = n_i -> Interpol(oTref)
		oaNe[i] = n_e -> Interpol(oTref)
		oaVi[i] = v_i -> Interpol(oTref)
		oaVe[i] = v_e -> Interpol(oTref)
		oaR[i]  = r   -> Interpol(oTref)
		
		;Current density
		oaJ[i] = 1e15 * q * oaNe[i] * (oaVi[i] - oaVe[i])
		jtemp  = (oaJ[i])['DATA']
		FOR j = 0, 2 DO jtemp[0,j] = Convol(jtemp[*,j], [0.25, 0.5, 0.25])
		oaJ[i] -> SetData, Temporary(jtemp)
		oaJ[i] -> SetName, j_vnames[i]
		oaJ[i] -> Cache
		oaJ[i] -> AddAttr, 'CATDESC', 'Current density computed from particle moments.'
		oaJ[i] -> AddAttr, 'TITLE',   'J!C($\mu$A/m$\up2$'
		oaJ[i] -> AddAttr, 'UNITS',   'uA/m^2'
	ENDFOR
		
;-------------------------------------------
; Current Density //////////////////////////
;-------------------------------------------
	
	;Reciprocal vectors for taking derivatives
	oRecipVec = MrVar_RecipVec(oaR[0], oaR[1], oaR[2], oaR[3])
	
	;Current density
	oCurlB = (1e-6/mu0) * oRecipVec -> Curl( oaB[0], oaB[1], oaB[2], oaB[3] )
	oCurlB -> SetName, curlb_vname
	oCurlB -> Cache
	oCurlB['CATDESC'] = 'Current density derived from the curlometer technique using MEC & ' + fgm_instr + '.'
	oCurlB['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oCurlB['LABEL']   = 'J$\down' + ['X', 'Y', 'Z'] + '$'
	oCurlB['TITLE']   = nabla + 'xB/$\mu$$\down0$!C($\mu$A/m^2)'
	oCurlB['UNITS']   = 'uA/m^2'
	
;-------------------------------------------
; Divergence of Current Density ////////////
;-------------------------------------------
	
	;Charge Density
	;   - 1e-3 converts 1/km to 1/m ==> div(J) = uA/m^3
	oDivJ = (1e-15/q) * oRecipVec -> Divergence( oaJ[0], oaJ[1], oaJ[2], oaJ[3] )
	oDivJ -> SetName, divj_vname
	oDivJ -> Cache
	oDivJ['CATDESC'] = 'Divergence of the current density'
	oDivJ['COLOR']   = 'Black'
	oDivJ['LABEL']   = nabla + '.J/q'
	oDivJ['TITLE']   = nabla + '.J/q!C(cm$\up-3$ s$\up-1$)'
	oDivJ['UNITS']   = 'cm^-3 s^-1'
		
;-------------------------------------------
; Charge Density ///////////////////////////
;-------------------------------------------
	
	;Charge Density
	;   - 1e10 converts to uC/cm^2
	oRho = (e0 / q * 1e-12) * oRecipVec -> Divergence( oaE[0], oaE[1], oaE[2], oaE[3] )
	oRho -> SetData, Convol(oRho['DATA'], [0.25, 0.5, 0.25])
	oRho -> SetName, rho_vname
	oRho -> Cache
	oRho['CATDESC'] = 'Free charge density, computed from the divergence of the electric field.'
	oRho['TITLE']   = '$\Delta$$\downe$!C(1/cm^3)'
	oRho['UNITS']   = '1/cm^3'
	
;-------------------------------------------
; Time Derivative of Charge Density ////////
;-------------------------------------------
	
	dt       = oTref['DATA', 1:-1, 'SSM'] - oTref['DATA', 0:-2, 'SSM']
	dRho_dt  = (oRho['DATA',1:-1] - oRho['DATA',0:-2]) / dt
	odRho_dt = MrScalarTS(oTref[1:-1], dRho_dt)
	odRho_dt -> SetName, drho_dt_vname
	odRho_dt -> Cache
	odRho_dt['CATDESC'] = 'Time derivative of the free charge density.'
	odRho_dt['COLOR']   = 'Blue'
	odRho_dt['LABEL']   = 'd$\rho$/dt'
	odRho_dt['TITLE']   = 'd$\rho$/dt!C(cm$\up-3$ s$\up-1$)'
	odRho_dt['UNITS']   = 'cm^-3 s^-1'
		
;-------------------------------------------
; Barycenter ///////////////////////////////
;-------------------------------------------
	;B
	oBbc = (oaB[0] + oaB[1] + oaB[2] + oaB[3]) / 4.0
	oBbc -> SetName, b_bc_vname
	oBbc -> Cache
	oBbc['CATDESC'] = 'Barycentric average of the vector magnetic field.'
	oBbc['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oBbc['LABEL']   = 'B$\down' + ['X', 'Y', 'Z'] + '$'
	oBbc['TITLE']   = 'B!C(nT)'
	oBbc['UNITS']   = 'nT'
	
	;E
	oEbc = (oaE[0] + oaE[1] + oaE[2] + oaE[3]) / 4.0
	oEbc -> SetName, e_bc_vname
	oEbc -> Cache
	oEbc['CATDESC'] = 'Barycentric average of the vector electric field.'
	oEbc['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oEbc['LABEL']   = 'E$\down' + ['X', 'Y', 'Z'] + '$'
	oEbc['TITLE']   = 'E!C(mV/m)'
	oEbc['UNITS']   = 'mV/m'
	
	;J
	oJbc = (oaJ[0] + oaJ[1] + oaJ[2] + oaJ[3]) / 4.0
	oJbc -> SetName, j_bc_vname
	oJbc -> Cache
	oJbc['CATDESC'] = 'Barycentric average of the current density.'
	oJbc['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oJbc['LABEL']   = 'J$\down' + ['X', 'Y', 'Z'] + '$'
	oJbc['TITLE']   = 'J!C($\mu$A/m$\up2$)'
	oJbc['UNITS']   = 'uA/m^2'
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	oJbc['AXIS_RANGE'] = [ Min([oJbc.min, oCurlB.min]), Max([oJbc.max, oCurlB.max]) ]
	
	oCurlB['AXIS_RANGE'] = oJbc['AXIS_RANGE']
;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot MMS1 data
	win = MrVar_PlotTS( [ b_bc_vname, e_bc_vname, j_bc_vname, curlb_vname, $
	                      rho_vname, divj_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	
	;Overplot MMS 2-4 data
	win = MrVar_OPlotTS( divj_vname, drho_dt_vname)
		
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
		fname   = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-q-continuity'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END