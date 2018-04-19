; docformat = 'rst'
;
; NAME:
;       MrVar_Plot_CVA
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
;   Create a plot to both assist in and display the results of determining the boundary
;   velocity and normal direction of a 1D structure via the Constant Velocity Approach (CVA).
;
; :Categories:
;   MMS
;
; :Params:
;       T:          in, required, type=FltArr(4)
;                   The times, in seconds, at which each spacecraft encountered the 1D
;                       boundary. Analysis assumes the order is 1, 2, 3, 4.
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;
; :Keywords:
;       COORDS:             in, optional, type=string, default='gse'
;                           Coordinate system of the original data.
;       INSTR:              in, optional, type=string, default='fgm'
;                           Name of the instrument to be used in determining the times
;                               at which each spacecraft crossed the boundary.
;       LEVEL:              in, optional, type=string, default='l2'
;                           Quality level of data to be loaded.
;       NO_LOAD:            in, optional, type=boolean, default=0
;                           If set, data is not loaded from source files. Instead, it
;                               is taken from the variable cache.
;       OUTPUT_DIR:         in, optional, type=string, default=pwd
;                           A directory in which to save the figure. If neither `OUTPUT_DIR`
;                               nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT:         in, optional, type=string, default=pwd
;                           File extensions for the output figure. Options include: 'eps', 'gif',
;                               'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                               `OUTPUT_EXT` are defined, no file is generated.
;       TEST:               in, optional, type=boolean, default=0
;                           If set, no CVA analysis is performed and `T` is ignored.
;                               Instead, data from `INSTR` is plotted and the tools in the
;                               window GUI can be used to select times at which the space-
;                               craft cross the boundary.
;       TRANGE:             in, optional, type=strarr(2), default=MrVar_GetTRange()
;                           Date-time string, formatted as YYYY-MM-DDThh:mm:ss.fff indicating
;                               time interval for which data is loaded. If given, the
;                               global analysis interval is changed via MrVar_SetTRange.
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
;       2018/02/12  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_CVA, t, mode, $
COORDS=coords, $
INSTR=instr, $
LEVEL=level, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
TEST=test, $
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
	tf_test = Keyword_Set(test)
	IF N_Elements(coords)  EQ 0 THEN coords  = 'gse'
	IF N_Elements(ephdesc) EQ 0 THEN ephdesc = ''
	IF N_Elements(instr)   EQ 0 THEN instr   = 'fgm'
	IF N_Elements(level)   EQ 0 THEN level   = 'l2'
	IF N_Elements(trange)  GT 0 THEN MrVar_SetTRange, trange
	IF N_Elements(t_mva)   EQ 0 THEN t_mva   = MrVar_GetTRange()
	
;-------------------------------------------
; Data Parameters //////////////////////////
;-------------------------------------------
	;FGM
	fgm_instr   = instr
	fgm_mode    = (mode   EQ 'fast' || mode   EQ 'slow') ? 'srvy' : mode
	fgm_coords  = (coords EQ 'dsl'  || coords EQ 'dbcs') ? 'dmpa' : coords
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	CASE fgm_instr OF
		'fsm': fgm_level = 'l3'
		'fgm': fgm_level = 'l2'
		'dfg': fgm_level = 'l2pre'
		'afg': fgm_level = 'l2pre'
		ELSE: ;Do nothing -- another instrument was selected.
	ENDCASE
	
	;EDP
	edp_instr  = 'edp'
	edp_coords = (coords EQ 'dmpa' || coords EQ 'dbcs') ? 'dsl'  : coords
	edp_mode   = mode
	IF instr EQ 'edp' && mode EQ 'srvy' THEN BEGIN
		MrPrintF, 'LogWarn', 'EDP does not have "srvy" data. Using "fast".'
		edp_mode = 'fast'
	ENDIF
	
	;FPI
	fpi_coords = (coords EQ 'dmpa' || coords EQ 'dsl') ? 'dbcs' : coords
	fpi_mode   = mode
	IF (instr EQ 'des' || instr EQ 'dis') && (mode EQ 'srvy' || mode EQ 'slow') THEN BEGIN
		MrPrintF, 'LogWarn', 'FPI does not have "srvy" or "slow" data. Using "fast".'
		fpi_mode = 'fast'
	ENDIF
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	sc = 'mms' + ['1', '2', '3', '4']
	
	;FGM
	b_vnames    = sc + '_' + StrJoin([fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_')
	IF fgm_instr EQ 'fsm' THEN BEGIN
		bvec_vnames = sc + '_' + StrJoin([fgm_instr, 'b', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([fgm_instr, 'b', 'mag',      fgm_mode, fgm_level], '_')
	ENDIF ELSE BEGIN
		bvec_vnames = sc + '_' + StrJoin([fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_')
	ENDELSE
	
	;EDP
	e_vnames    = sc + '_' + StrJoin([edp_instr, 'dce',  edp_coords, edp_mode, level], '_' )
	v_vnames    = sc + '_' + StrJoin([edp_instr, 'scpot',            edp_mode, level], '_' )
	
	;FPI
	ne_vnames = sc + '_' + StrJoin(['des', 'numberdensity',      fpi_mode, level], '_' )
	ve_vnames = sc + '_' + StrJoin(['des', 'bulkv',  fpi_coords, fpi_mode, level], '_' )
	ni_vnames = sc + '_' + StrJoin(['dis', 'numberdensity',      fpi_mode, level], '_' )
	vi_vnames = sc + '_' + StrJoin(['dis', 'bulkv',  fpi_coords, fpi_mode, level], '_' )
	
	;EPHEMERIS
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	bx_vnames  = bvec_vnames  + '_x'
	by_vnames  = bvec_vnames  + '_y'
	bz_vnames  = bvec_vnames  + '_z'
	ex_vnames  = e_vnames  + '_x'
	ey_vnames  = e_vnames  + '_y'
	ez_vnames  = e_vnames  + '_z'
	vex_vnames = ve_vnames + '_x'
	vey_vnames = ve_vnames + '_y'
	vez_vnames = ve_vnames + '_z'
	vix_vnames = vi_vnames + '_x'
	viy_vnames = vi_vnames + '_y'
	viz_vnames = vi_vnames + '_z'
	
;-------------------------------------------
; Output Variables /////////////////////////
;-------------------------------------------
		
	;FGM
	IF instr EQ 'fgm' || instr EQ 'dfg' || instr EQ 'afg' THEN BEGIN
		nParams    = 4
		vnames     = StrArr(4, nParams)
		vec_vnames  = bvec_vnames
		vnames[0,*] = [bmag_vnames[0], bx_vnames[0], by_vnames[0], bz_vnames[0]]
		vnames[1,*] = [bmag_vnames[1], bx_vnames[1], by_vnames[1], bz_vnames[1]]
		vnames[2,*] = [bmag_vnames[2], bx_vnames[2], by_vnames[2], bz_vnames[2]]
		vnames[3,*] = [bmag_vnames[3], bx_vnames[3], by_vnames[3], bz_vnames[3]]
	;EDP
	ENDIF ELSE IF instr EQ 'edp' THEN BEGIN
		nParams    = 3
		sc1_vnames = [ex_vnames[0], ey_vnames[0], ez_vnames[0]]
		sc2_vnames = [ex_vnames[1], ey_vnames[1], ez_vnames[1]]
		sc3_vnames = [ex_vnames[2], ey_vnames[2], ez_vnames[2]]
		sc4_vnames = [ex_vnames[3], ey_vnames[3], ez_vnames[3]]
	;DES
	ENDIF ELSE IF instr EQ 'des' THEN BEGIN
	
	;DIS
	ENDIF ELSE IF instr EQ 'dis' THEN BEGIN
	
	;OTHER
	ENDIF ELSE BEGIN
		Message, 'Instrument not recognized: "' + param + '".'
	ENDELSE
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
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
		
		;FGM
		IF instr EQ 'fgm' || instr EQ 'dfg' || instr EQ 'afg' THEN BEGIN
			MrMMS_FGM_Load_Data, sc, fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = level, $
			                     VARFORMAT = '*b_'+fgm_coords+'_'+fgm_mode+'*'
		;EDP
		ENDIF ELSE IF instr EQ 'edp' THEN BEGIN
			MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
			                 OPTDESC   = 'dce', $
			                 VARFORMAT = '*_dce_'+edp_coords+'_*'
		;DES
		ENDIF ELSE IF instr EQ 'des' THEN BEGIN
		
		;DIS
		ENDIF ELSE IF instr EQ 'dis' THEN BEGIN
		
		;OTHER
		ENDIF ELSE BEGIN
			Message, 'Instrument not recognized: "' + param + '".'
		ENDELSE
	ENDIF
	
;-------------------------------------------
; Split Components /////////////////////////
;-------------------------------------------
	IF N_Elements(vec_vnames) GT 0 THEN BEGIN
		FOR i = 0, 3 DO BEGIN
			oV = MrVar_Get(vec_vnames[i])
			oV -> Split, oVx, oVy, oVz, /CACHE
		ENDFOR
	ENDIF
	
;-------------------------------------------
; CVA //////////////////////////////////////
;-------------------------------------------
	oR1 = MrVar_Get(r_vnames[0])
	oR2 = MrVar_Get(r_vnames[1])
	oR3 = MrVar_Get(r_vnames[2])
	oR4 = MrVar_Get(r_vnames[3])
	
	oR2 = oR2 -> Interpol(oR1)
	oR3 = oR3 -> Interpol(oR1)
	oR4 = oR4 -> Interpol(oR1)
	
	IF ~tf_test THEN V = MrVar_CVA(oR1, oR2, oR3, oR4, t[0], t[1], t[2], t[3], N=N)
	
;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	title     = 'Constant Velocity Approach'
	sc_colors = ['Black', 'Red', 'Forest Green', 'Blue']
	yrange    = FltArr(4,nParams,2)
	
	;Color & Label
	FOR i = 0, 3 DO BEGIN
		FOR j = 0, nParams - 1 DO BEGIN
			oV = MrVar_Get(vnames[i,j])
			oV['COLOR']   = sc_colors[i]
			oV['LABEL']   = 'MMS' + String(i+1, FORMAT='(i0)')
			yrange[i,j,*] = [oV.min, oV.max]
		ENDFOR
	ENDFOR
	
	;YRANGE
	FOR j = 0, nParams-1 DO BEGIN
		oV = MrVar_Get(vnames[0,j])
		oV['AXIS_RANGE'] = [Min(yrange[*,j,0]), Max(yrange[*,j,1])]
	ENDFOR
	
	;FGM
	IF instr EQ 'fgm' || instr EQ 'dfg' || instr EQ 'afg' THEN BEGIN
		oBmag = MrVar_Get(bmag_vnames[0])
		oBmag['PLOT_TITLE'] = title
		
		oBx = MrVar_Get(bx_vnames[0])
		oBx['TITLE'] = 'Bx!C(nT)'
		oBy = MrVar_Get(by_vnames[0])
		oBy['TITLE'] = 'By!C(nT)'
		oBz = MrVar_Get(bz_vnames[0])
		oBz['TITLE'] = 'Bz!C(nT)'
	;EDP
	ENDIF ELSE IF instr EQ 'edp' THEN BEGIN
		nParams    = 3
		sc1_vnames = [ex_vnames[0], ey_vnames[0], ez_vnames[0]]
		sc2_vnames = [ex_vnames[1], ey_vnames[1], ez_vnames[1]]
		sc3_vnames = [ex_vnames[2], ey_vnames[2], ez_vnames[2]]
		sc4_vnames = [ex_vnames[3], ey_vnames[3], ez_vnames[3]]
	;DES
	ENDIF ELSE IF instr EQ 'des' THEN BEGIN
	
	;DIS
	ENDIF ELSE IF instr EQ 'dis' THEN BEGIN
	
	;OTHER
	ENDIF ELSE BEGIN
		Message, 'Instrument not recognized: "' + instr + '".'
	ENDELSE
	
;-------------------------------------------
; Plot Time Series /////////////////////////
;-------------------------------------------
	win = MrWindow( OYMARGIN = [9, 2], $
	                REFRESH  = 0, $
	                XSIZE    = 800, $
	                YGAP     = 0.5, $
	                YSIZE    = 700 )
	
	win = MrVar_PlotTS( vnames[0,*], $
	                    /CURRENT, $
	                    /NO_REFRESH )
	
	;Overplot other spacecraft
	FOR i = 1, 3 $
		DO win  = MrVar_OPlotTS( vnames[0,*], vnames[i,*] )
	
	;Draw zero-lines
	FOR j = 0, nParams - 1 $
		DO line = MrPlotS( win[vnames[0,j]].xrange, [0,0], LINESTYLE='--', TARGET=win[vnames[0,j]] )
	

;-------------------------------------------
; Annotate /////////////////////////////////
;-------------------------------------------
	IF ~tf_test THEN BEGIN
		;Outline the MVA interval
		FOR j = 0, nParams - 1 DO BEGIN
			gfx = win[2*j]
			ps1 = MrPlotS( t[[0,0]], gfx.yrange, $
			               COLOR     = 'Blue', $
			               NAME      = 'Line: CVA1', $
			               TARGET    = gfx, $
			               LINESTYLE = 2 )
			ps2 = MrPlotS( t[[1,1]], gfx.yrange, $
			               COLOR     = 'Blue', $
			               NAME      = 'Line: CVA2', $
			               TARGET    = gfx, $
			               LINESTYLE = 2 )
			ps2 = MrPlotS( t[[2,2]], gfx.yrange, $
			               COLOR     = 'Blue', $
			               NAME      = 'Line: CVA3', $
			               TARGET    = gfx, $
			               LINESTYLE = 2 )
			ps2 = MrPlotS( t[[3,3]], gfx.yrange, $
			               COLOR     = 'Blue', $
			               NAME      = 'Line: CVA4', $
			               TARGET    = gfx, $
			               LINESTYLE = 2 )
		ENDFOR
	
		;Print eigenvectors and eigenvalues
		order     = Sort(t)+1
		st        = 't_0 = [' + StrJoin(StrMid(ssm_to_hms(t), 0, 13), ', ') + ']'
		svel      = String( V,   FORMAT='(%"|V| = %9.4f")' )
		sorder    = String( order, FORMAT='(%"Order = [%1i, %1i, %1i, %1i]")' )
		snormal   = String( N,     FORMAT='(%"n = [%7.4f, %7.4f, %7.4f]")' )
		svelocity = String( N*V,   FORMAT='(%"V = [%9.4f, %9.4f, %9.4f]")' )
		txt1 = MrText( 0.50, 0.07, st,        ALIGNMENT=0.5, NAME='Txt: Times', /NORMAL )
		txt1 = MrText( 0.25, 0.04, sorder,    ALIGNMENT=0.5, NAME='Txt: Order', /NORMAL )
		txt2 = MrText( 0.25, 0.01, svel,      ALIGNMENT=0.5, NAME='Txt: |V|', /NORMAL )
		txt3 = MrText( 0.75, 0.04, snormal,   ALIGNMENT=0.5, NAME='Txt: Normal',   /NORMAL )
		txt4 = MrText( 0.75, 0.01, svelocity, ALIGNMENT=0.5, NAME='Txt: Velocity', /NORMAL )
	ENDIF

;-------------------------------------------
; Prettify /////////////////////////////////
;-------------------------------------------
	
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[12,8]
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
		fname   = StrJoin( ['mms', instr, mode, level, 'cva'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, win
END