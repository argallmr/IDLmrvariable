; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_JVA
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
;   Perform a minimum or maximum variance analysis on the electric or magnetic field.
;   Plot both the original and rotated data as well as hodograms in the variance frame.
;
; :Categories:
;   MMS
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;       T_JVA:      in, optional, type=strarr(2), default=MrVar_GetTRange()
;                   Date-time string, formatted as YYYY-MM-DDThh:mm:ss.fff indicating the
;                       time interval over which to apply JVA.
;
; :Keywords:
;       COORDS:             in, optional, type=string, default='gse'
;                           Coordinate system of the original data.
;       FGM_INSTR:          in, optional, type=string, default='fgm'
;                           Name of the FGM instrument used as the background magnetic
;                               field when calculating the wave normal angle. Options are
;                               {'afg' | 'dfg' | 'fgm'}.
;       LEVEL:              in, optional, type=string, default='l2'
;                           Quality level of data to be loaded.
;       NO_LOAD:            in, optional, type=boolean, default=0
;                           If set, data is not loaded from source files. Instead, it
;                               is taken from the variable cache.
;       OUTPUT_EXT:         in, optional, type=string/strarr, default='png'
;                           Type of image to save. Options include 'eps', 'png', 'gif',
;                               'jpeg', 'ps', 'pdf', 'tiff'. If neither `OUTPUT_DIR` nor
;                               `OUTPUT_EXT` are give, no file is produced.
;       OUTPUT_DIR:         in, optional, type=string, default='~/'
;                           Directroy in which to save the image. If neither `OUTPUT_DIR`
;                               nor `OUTPUT_EXT` are give, no file is produced.
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
;       2017/05/16  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_JVA, sc, mode, t_jva, $
COORDS=coords, $
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
	IF N_Elements(instr)     EQ 0 THEN instr     = 'fgm'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	IF N_Elements(t_jva)     EQ 0 THEN t_jva     = MrVar_GetTRange()
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	;FGM parameters
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	fgm_coords = (coords EQ 'dsl'  || coords EQ 'dbcs') ? 'dmpa' : coords
	fgm_mode   = (mode   EQ 'fast' || mode   EQ 'slow') ? 'srvy' : mode
	
	;EDP parameters
	edp_coords = (coords EQ 'dmpa' || coords EQ 'dbcs') ? 'dsl'  : coords
	IF mode EQ 'srvy' THEN BEGIN
		MrPrintF, 'LogWarn', 'EDP does not have SRVY data. Using FAST.'
		edp_mode = 'fast'
	ENDIF ELSE BEGIN
		edp_mode = mode
	ENDELSE
	
	;Load the data
	IF tf_load THEN BEGIN
		;EDP E-Field
		MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = '*_dce_'+edp_coords+'_*'
		
		
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = level, $
		                     VARFORMAT = '*b_'+fgm_coords+'_'+fgm_mode+'*'
	ENDIF

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	b_vname    = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, fgm_mode, level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, fgm_mode, level], '_' )
	e_vname    = StrJoin( [sc, 'edp',     'dce',  edp_coords, edp_mode, level], '_' )
	
	;output names
	b_jva_vname = StrJoin( [sc, fgm_instr, 'b',    'jva', fgm_mode, level], '_' )
	e_jva_vname = StrJoin( [sc, 'edp',     'dce',  'jva', edp_mode, level], '_' )

;-------------------------------------------
; JVA //////////////////////////////////////
;-------------------------------------------

	;Grab the data
	oB = MrVar_Get(bvec_vname)
	oE = MrVar_Get(e_vname)
	
	;Find time interval
	oTb = oB['TIMEVAR']
	itb = oTb -> Value_Locate(t_jva) > 0
	
	oTe = oE['TIMEVAR']
	ite = oTe -> Value_Locate(t_jva) > 0
	
	;Minimum variance
	oEigVec = MrVar_JVA( oB[itb[0]:itb[1], *], oE[ite[0]:ite[1], *], $
	                     EIGENVALS_MAX = oEigVal_Max, $
	                     EIGENVALS_MIN = oEigVal_Min, $
	                     EIGENVEC_MAX  = oEigVec_Max, $
	                     EIGENVEC_MIN  = oEigVec_Min, $
	                     THETA         = theta )
	
	;Rotate field MVA coordiantes
	oB_jva = oEigVec ## oB
	oE_jva = oEigVec ## oE
	
	;Turn them into vectors
	oB_jva = MrVectorTS( oTb, oB_jva, /CACHE, NAME=b_jva_vname )
	oE_jva = MrVectorTS( oTe, oE_jva, /CACHE, NAME=e_jva_vname )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	;B XYZ
	oB['COLOR']      = ['Blue', 'Forest Green', 'Red']
	oB['LABEL']      = 'B$\down' + ['X', 'Y', 'Z'] + '$'
	oB['PLOT_TITLE'] = 'JVA Analysis'
	oB['TITLE']      = 'B!C(nT)'
	
	;E XYZ
	oE['COLOR']      = ['Blue', 'Forest Green', 'Red']
	oE['LABEL']      = 'E$\down' + ['X', 'Y', 'Z'] + '$'
	oE['PLOT_TITLE'] = 'JVA Analysis'
	oE['TITLE']      = 'E!C(mV/m)'
	
	;B JVA
	oB_jva['COLOR'] = ['Blue', 'Forest Green', 'Red']
	oB_jva['LABEL'] = 'B$\down' + ['N', 'M', 'L'] + '$'
	oB_jva['TITLE'] = 'B!C(nT)'
	oB_jva['UNITS'] = 'nT'
	
	;E JVA
	oE_jva['COLOR'] = ['Blue', 'Forest Green', 'Red']
	oE_jva['LABEL'] = 'E$\down' + ['N', 'M', 'L'] + '$'
	oE_jva['TITLE'] = 'E!C(mV/m)'
	oE_jva['UNITS'] = 'mV/m'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	win = MrWindow( LAYOUT   = [2,6], $
	                OXMARGIN = [10,6], $
	                XGAP     = 13, $
	                XSIZE    = 740, $
	                YGAP     = 0.5, $
	                YSIZE    = 775, $
	                REFRESH  = 0 )
	
	layout = MrLayout( [3,2], $
	                   ASPECT   = 1.0, $
	                   OXMARGIN = [12,8], $
	                   OYMARGIN = [15,29], $
	                   XGAP     = 14, $
	                   YGAP     = 6 )
	
	;B XYZ
	p1 = MrVar_Plot( oB, $
	                 /CURRENT, $
	                 LAYOUT = [1,1] )
	
	;E XYZ
	p2 = MrVar_Plot( oE, $
	                 /CURRENT, $
	                 LAYOUT = [1,2] )
	
	;B JVA
	p3 = MrVar_Plot( oB_jva, $
	                 /CURRENT, $
	                 LAYOUT = [2,1] )
	
	;E JVA
	p4 = MrVar_Plot( oE_jva, $
	                 /CURRENT, $
	                 LAYOUT = [2,2] )

;-------------------------------------------
; B Hodogram ///////////////////////////////
;-------------------------------------------
	cgLoadCT, 13

	oB_jva -> RemoveAttr, ['COLOR', 'LABEL']
	brange = [ Min(oB_jva['DATA', itb[0]:itb[1], *], MAX=bmax), bmax ]
	
	;BL-BN
	p5 = MrVar_Plot( oB_jva[itb[0]:itb[1], 0], oB_jva[itb[0]:itb[1], 2], $
	                 /CURRENT, $
	                 POSITION = layout[*,0] )
	x5a = MrPlotS( oB_jva['DATA', itb[0]:itb[1], 0], oB_jva['DATA', itb[0]:itb[1], 2], $
	               COLOR  = BytScl(oTb['DATA', itb[0]:itb[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p5 )
	
	;BL-BM
	p6 = MrVar_Plot( oB_jva[itb[0]:itb[1], 1], oB_jva[itb[0]:itb[1], 2], $
	                 /CURRENT, $
	                 POSITION = layout[*,1] )
	x6 = MrPlotS( oB_jva['DATA', itb[0]:itb[1], 1], oB_jva['DATA', itb[0]:itb[1], 2], $
	               COLOR  = BytScl(oTb['DATA', itb[0]:itb[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p6 )
	
	;BM-BN
	p7 = MrVar_Plot( oB_jva[itb[0]:itb[1], 0], oB_jva[itb[0]:itb[1], 1], $
	                 /CURRENT, $
	                 POSITION = layout[*,2] )
	x7 = MrPlotS( oB_jva['DATA', itb[0]:itb[1], 0], oB_jva['DATA', itb[0]:itb[1], 1], $
	               COLOR  = BytScl(oTb['DATA', itb[0]:itb[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p7 )

;-------------------------------------------
; E Hodogram ///////////////////////////////
;-------------------------------------------
	oE_jva -> RemoveAttr, ['COLOR', 'LABEL']
	erange = [ Min(oE_jva['DATA', ite[0]:ite[1], *], MAX=emax), emax ]
	
	;EL-EN
	p8 = MrVar_Plot( oE_jva[ite[0]:ite[1], 0], oE_jva[ite[0]:ite[1], 2], $
	                 /CURRENT, $
	                 POSITION = layout[*,3] )
	x8a = MrPlotS( oE_jva['DATA', ite[0]:ite[1], 0], oE_jva['DATA', ite[0]:ite[1], 2], $
	               COLOR  = BytScl(oTe['DATA', ite[0]:ite[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p8 )
	
	;EL-EM
	p9 = MrVar_Plot( oE_jva[ite[0]:ite[1], 1], oE_jva[ite[0]:ite[1], 2], $
	                 /CURRENT, $
	                 POSITION = layout[*,4] )
	x9 = MrPlotS( oE_jva['DATA', ite[0]:ite[1], 1], oE_jva['DATA', ite[0]:ite[1], 2], $
	               COLOR  = BytScl(oTe['DATA', ite[0]:ite[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p9 )
	
	;EM-EN
	p10 = MrVar_Plot( oE_jva[ite[0]:ite[1], 0], oE_jva[ite[0]:ite[1], 1], $
	                  /CURRENT, $
	                  POSITION = layout[*,5] )
	x10 = MrPlotS( oE_jva['DATA', ite[0]:ite[1], 0], oE_jva['DATA', ite[0]:ite[1], 1], $
	               COLOR  = BytScl(oTe['DATA', ite[0]:ite[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p10 )

;-------------------------------------------
; Annotate /////////////////////////////////
;-------------------------------------------
	tr = MrVar_GetTRange('SSM')
	tt = oTb[ 'DATA', [itb[0], itb[1]], 'SSM' ]

	;
	;Outline the JVA interval
	;
	
	;B
	ps1 = MrPlotS( tt[[0,0]], p1.yrange, $
	               COLOR     = 'Blue', $
	               NAME      = 'Line: B-JVA0', $
	               TARGET    = p1, $
	               LINESTYLE = 2 )
	ps2 = MrPlotS( tt[[1,1]], p1.yrange, $
	               COLOR     = 'Blue', $
	               NAME      = 'Line: B-JVA1', $
	               TARGET    = p1, $
	               LINESTYLE = 2 )

	;E
	ps3 = MrPlotS( tt[[0,0]], p2.yrange, $
	               COLOR     = 'Blue', $
	               NAME      = 'Line: E-JVA0', $
	               TARGET    = p2, $
	               LINESTYLE = 2 )
	ps4 = MrPlotS( tt[[1,1]], p2.yrange, $
	               COLOR     = 'Blue', $
	               NAME      = 'Line: E-JVA1', $
	               TARGET    = p2, $
	               LINESTYLE = 2 )
	
	;Print eigenvectors and eigenvalues
	jointvec  = String( oEigVec['DATA'], FORMAT='(f6.3)' )
	sjva      = 'Joint Variance:'                              + '!C' + $
	            'N = [' + StrJoin( jointvec[0:2], ', ' ) + ']' + '!C' + $
	            'M = [' + StrJoin( jointvec[3:5], ', ' ) + ']' + '!C' + $
	            'L = [' + StrJoin( jointvec[6:8], ', ' ) + ']'
	
	maxval = String( oEigVal_Max['DATA'], FORMAT='(f0.3)' )
	maxvec = String( oEigVec_Max['DATA'], FORMAT='(f6.3)' )
	smvae  = 'Maximum Variance:'                           + '!C' + $
	         '\lambda = [' + StrJoin( maxval, ', ' ) + ']' + '!C' + $
	         'N = [' + StrJoin( maxvec[0:2], ', ' )  + ']' + '!C' + $
	         'M = [' + StrJoin( maxvec[3:5], ', ' )  + ']' + '!C' + $
	         'L = [' + StrJoin( maxvec[6:8], ', ' )  + ']'
	
	minval = String( oEigVal_Min['DATA'], FORMAT='(f0.3)' )
	minvec = String( oEigVec_Min['DATA'], FORMAT='(f6.3)' )
	smvab  = 'Minimum Variance:'                           + '!C' + $
	         '\lambda = [' + StrJoin( minval, ', ' ) + ']' + '!C' + $
	         'N = [' + StrJoin( minvec[0:2], ', ' )  + ']' + '!C' + $
	         'M = [' + StrJoin( minvec[3:5], ', ' )  + ']' + '!C' + $
	         'L = [' + StrJoin( minvec[6:8], ', ' )  + ']'
	
	;Analysis time range
	st_jva = 't_{JVA}=' + StrMid(t_jva[0], 11) + '-' + StrMid(t_jva[1], 11)

	txt1 = MrText( 0.01,  0.085, sjva,   ALIGNMENT=0.0, NAME='Txt: JVA',  /NORMAL )
	txt2 = MrText( 0.345, 0.085, smvae,  ALIGNMENT=0.0, NAME='Txt: MVAE', /NORMAL )
	txt3 = MrText( 0.68,  0.085, smvab,  ALIGNMENT=0.0, NAME='Txt: MVAB', /NORMAL )
	txt4 = MrText( 0.5,   0.11,  st_jva, ALIGNMENT=0.5, NAME='Txt: Time', /NORMAL )

;-------------------------------------------
; Prettify /////////////////////////////////
;-------------------------------------------
	trange = MrVar_GetTRange('SSM')

	p1 -> SetLayout, [1,1]
	p2 -> SetLayout, [2,1]
	p3 -> SetLayout, [1,2]
	p4 -> SetLayout, [2,2]

	p1  -> SetProperty, XRANGE=trange, XTICKFORMAT='(a1)', XTICKS=3, XTITLE=''
	p2  -> SetProperty, XRANGE=trange, XTICKFORMAT='(a1)', XTICKS=3, XTITLE=''
	p3  -> SetProperty, XRANGE=trange, XTICKS=3
	p4  -> SetProperty, XRANGE=trange, XTICKS=3
	p5  -> SetProperty, XRANGE=brange, XTITLE='B$\downN$', YTITLE='B$\downL$', YRANGE=brange
	p6  -> SetProperty, XRANGE=brange, XTITLE='B$\downM$', YTITLE='B$\downL$', YRANGE=brange
	p7  -> SetProperty, XRANGE=brange, XTITLE='B$\downN$', YTITLE='B$\downM$', YRANGE=brange
	p8  -> SetProperty, XRANGE=erange, XTITLE='E$\downN$', YTITLE='E$\downL$', YRANGE=erange
	p9  -> SetProperty, XRANGE=erange, XTITLE='E$\downM$', YTITLE='E$\downL$', YRANGE=erange
	p10 -> SetProperty, XRANGE=erange, XTITLE='E$\downN$', YTITLE='E$\downM$', YRANGE=erange

	win -> SetProperty, OXMARGIN=[10,6]

;-------------------------------------------
; Save /////////////////////////////////////
;-------------------------------------------
	IF N_Elements(output_ext) GT 0 || N_Elements(output_dir) GT 0 THEN BEGIN
		IF N_Elements(output_ext) EQ 0 THEN output_ext = 'png'
		IF N_Elements(output_dir) EQ 0 THEN output_dir = File_Search('~', /TEST_DIRECTORY)
		
		;Parse analysis time
		t0 = StrJoin( StrSplit(t_jva[0], '-T:', /EXTRACT) )
		t1 = StrJoin( StrSplit(t_jva[1], '-T:', /EXTRACT) )
		date0 = StrMid(t0, 0, 8)
		date1 = StrMid(t1, 0, 8)
		time0 = StrMid(t0, 8, 6) + ( StrLen(t0) LE 14 ? '' : '_p' + StrMid(t0, 15) )
		time1 = StrMid(t1, 8, 6) + ( StrLen(t1) LE 14 ? '' : '_p' + StrMid(t1, 15) )
		
		;Time stamp for file
		IF date0 NE date1 $
			THEN ftime = StrJoin([date0, time0, date1, time1], '_') $
			ELSE ftime = StrJoin([date0, time0, time1], '_')
			
		;File name
		fname = StrJoin([sc, 'edp-'+fgm_instr, mode, level, 'jva', ftime], '_')
		fname = FilePath( fname + '.' + output_ext, $
		                  ROOT_DIR = output_dir )

		;Save
		FOR i = 0, N_Elements(fname) - 1 DO win -> Save, fname[i]
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	win -> refresh
	return, win
end