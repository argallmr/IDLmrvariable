; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_DES
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
;   Generate a plot of FPI quantities:
;       1. FGM Bxyz
;       2. Energy Spectr Par
;       2. Energy Spectr Perp
;       2. Energy Spectr Anti
;       3. PA Spectr Low
;       3. PA Spectr Mid
;       3. PA Spectr High
;       4. Density
;       5. Vxyz
;       6. Tpara, Tperp
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
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
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
;       2017/01/13  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_DES, sc, mode, $
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
	IF N_Elements(level)   EQ 0 THEN level   = 'l2'
	IF N_Elements(mode)    EQ 0 THEN mode    = 'fast'
	IF N_Elements(trange)  GT 0 THEN MrVar_SetTRange, trange
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;FPI
	species     = 'e'
	fpi_instr   = 'd' + species + 's'
	fpi_mode    = mode EQ 'brst' ? mode : 'fast'
	fpi_optdesc = fpi_instr + '-moms'
	
	;FGM
	IF N_Elements(fgm_instr) EQ 0 THEN BEGIN
		CASE level OF
			'l2': fgm_instr = 'fgm'
			ELSE: fgm_instr = 'dfg'
		ENDCASE
	ENDIF
	IF N_Elements(coords) EQ 0 THEN BEGIN
		CASE level OF
			'ql': coords = 'dbcs'
			ELSE: coords = 'gse'
		ENDCASE
	ENDIF
	fgm_coords = coords EQ 'dbcs' ? 'dmpa' : coords
	fgm_mode   = mode   EQ 'brst' ? mode   : 'srvy'

	;Source names
	b_vname      = StrJoin( [sc, fgm_instr, 'b',     fgm_coords, fgm_mode, level], '_' )
	bvec_vname   = StrJoin( [sc, fgm_instr, 'bvec',  fgm_coords, fgm_mode, level], '_' )
	bmag_vname   = StrJoin( [sc, fgm_instr, 'bmag',  fgm_coords, fgm_mode, level], '_' )
	eanti_vname  = StrJoin( [sc, fpi_instr,     'energyspectr', 'anti',   fpi_mode], '_' )
	epar_vname   = StrJoin( [sc, fpi_instr,     'energyspectr', 'par',    fpi_mode], '_' )
	eperp_vname  = StrJoin( [sc, fpi_instr,     'energyspectr', 'perp',   fpi_mode], '_' )
	palow_vname  = StrJoin( [sc, fpi_instr,     'pitchangdist', 'lowen',  fpi_mode], '_' )
	pamid_vname  = StrJoin( [sc, fpi_instr,     'pitchangdist', 'miden',  fpi_mode], '_' )
	pahigh_vname = StrJoin( [sc, fpi_instr,     'pitchangdist', 'highen', fpi_mode], '_' )
	n_vname      = StrJoin( [sc, fpi_instr,     'numberdensity', fpi_mode], '_' )
	v_vname      = StrJoin( [sc, fpi_instr,     'bulkv', coords, fpi_mode], '_' )
	tpara_vname  = StrJoin( [sc, fpi_instr,     'temppara', fpi_mode], '_' )
	tperp_vname  = StrJoin( [sc, fpi_instr,     'tempperp', fpi_mode], '_' )
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = level, $
		                     VARFORMAT = b_vname

		;DIS
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = fpi_optdesc, $
		                     VARFORMAT = [eanti_vname, epar_vname, eperp_vname, $
		                                  palow_vname, pamid_vname, pahigh_vname, $
		                                  n_vname, v_vname, tpara_vname, tperp_vname]
	ENDIF
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	;BMAG
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, fpi_mode, level, fpi_optdesc], ' ' ) )
	
	;Epar
	oEpar = MrVar_Get(epar_vname)
	oEpar['TITLE'] = 'DEF!CPar'
	
	;EPERP
	oEperp = MrVar_Get(eperp_vname)
	oEperp['TITLE']    = 'DEF!CPerp!C' + oEperp['UNITS'] + ''
	
	;EANTI
	oEanti = MrVar_Get(eanti_vname)
	oEanti['TITLE']    = 'DEF!CAnti'

	;PAHIGH
	oPAhigh = MrVar_Get(pahigh_vname)
	oPAhigh['TITLE'] = 'DEF!CHigh-E'
	
	;PAMID
	oPAmid = MrVar_Get(paMID_vname)
	oPAmid['TITLE'] = 'DEF!CMid-E!C' + oPAhigh['UNITS'] + ''
	
	;PALOW
	oPAlow = MrVar_Get(palow_vname)
	oPAlow['TITLE'] = 'DEF!CLow-E'
	
	;N
	oN = MrVar_Get(n_vname)
	oN['TITLE'] = 'Ne!C(cm$\up-3$)'
	
	;V
	oV = MrVar_Get(v_vname)
	oV['TITLE'] = 'Ve!C(km/s)'
	oV['LABEL'] = ['Vx', 'Vy', 'Vz']

	;TPARA
	oTpara = MrVar_Get(tpara_vname)
	oTpara['COLOR'] = 'Blue'
	oTpara['LABEL'] = 'Tpara'
	oTpara['TITLE'] = 'Te!C(eV)'
	
	;TPERP
	oTperp = MrVar_Get(tperp_vname)
	oTperp['COLOR'] = 'Red'
	oTperp['LABEL'] = 'Tperp'
	oTperp['TITLE'] = 'Te!C(eV)'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [b_vname, epar_vname, eperp_vname, eanti_vname, $
	                     pahigh_vname, pamid_vname, palow_vname, $
	                     n_vname, v_vname, tpara_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 680, $
	                    YSIZE = 700 )
	
	;Temperatures
	win = MrVar_OPlotTS( tpara_vname, tperp_vname )

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
			MrPrintF, 'LogText', 'Saving file to: "' + output_dir + '".'
		ENDIF
		
		;File name
		fname   = StrJoin( [sc, fpi_instr, fpi_mode, level, fpi_optdesc], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, win
END