; docformat = 'rst'
;
; NAME:
;       MrVar_Hodogram
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
;       T_MVA:      in, optional, type=strarr(2), default=MrVar_GetTRange()
;                   Date-time string, formatted as YYYY-MM-DDThh:mm:ss.fff indicating the
;                       time interval over which to apply MVA.
;
; :Keywords:
;       COORDS:             in, optional, type=string, default='gse'
;                           Coordinate system of the original data.
;       FGM_INSTR:          in, optional, type=string, default='fgm'
;                           Name of the FGM instrument used as the background magnetic
;                               field when calculating the wave normal angle. Options are
;                               {'afg' | 'dfg' | 'fgm'}. If FGM_INSTR is the same as
;                               `INSTR`, it is ignored.
;       INSTR:              in, optional, type=string, default='fgm'
;                           Name of the instrument on which minimum variance analysis is
;                               performed. Options are: {'afg', 'dfg', 'edp', 'fgm', 'scm'}.
;       LEVEL:              in, optional, type=string, default='l2'
;                           Quality level of data to be loaded.
;       MAXIMUM:            in, optional, type=boolean, default=0
;                           If set, a maximum variance analysis is performed.
;       NO_LOAD:            in, optional, type=boolean, default=0
;                           If set, data is not loaded from source files. Instead, it
;                               is taken from the variable cache.
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
;       2017/05/04  -   Written by Matthew Argall
;       2017/05/15  -   AFG, DFG, and FGM may be used. Added MAXIMUM keyword. - MRA
;-
FUNCTION MrMMS_Plot_MVA, sc, mode, t_mva, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
MAXIMUM=maximum, $
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
	
	tf_load    = ~Keyword_Set(no_load)
	tf_maximum = Keyword_Set(maximum)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	IF N_Elements(t_mva)     EQ 0 THEN t_mva     = MrVar_GetTRange()
	
;-------------------------------------------
; Data Parameters //////////////////////////
;-------------------------------------------
	;FGM parameters
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
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	b_vname    = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, fgm_mode, level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, fgm_mode, level], '_' )
	bmag_vname = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, fgm_mode, level], '_' )
	e_vname    = StrJoin( [sc, 'edp',     'dce',  edp_coords, edp_mode, level], '_' )
	
	;output names
	b_mva_vname = StrJoin( [sc, fgm_instr, 'b',    'mva', fgm_mode, level], '_' )
	e_mva_vname = StrJoin( [sc, 'edp',     'dce',  'mva', edp_mode, level], '_' )
	
	IF tf_maximum THEN BEGIN
		param     = 'E'
		varname   = e_vname
		mva_vname = e_mva_vname
	ENDIF ELSE BEGIN
		param     = 'B'
		varname   = bvec_vname
		mva_vname = b_mva_vname
	ENDELSE
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;EDP E-Field
		IF tf_maximum THEN BEGIN
			MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
			                 OPTDESC   = 'dce', $
			                 VARFORMAT = '*_dce_'+edp_coords+'_*'
		
		;FGM
		ENDIF ELSE BEGIN
			MrMMS_FGM_Load_Data, sc, fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = level, $
			                     VARFORMAT = '*b_'+fgm_coords+'_'+fgm_mode+'*'
		ENDELSE
	ENDIF
	
;-------------------------------------------
; MVA //////////////////////////////////////
;-------------------------------------------

	;Grab the data
	oVec = MrVar_Get(varname)
	oVec['PLOT_TITLE'] = 'Minimum Variance Analysis'
	
	;Find time interval
	oT = oVec['TIMEVAR']
	it = oT -> Value_Locate(t_mva) > 0
	
	;Minimum variance
	oEigVec = MrVar_MVA( oVec[it[0]:it[1], *], $
	                     EIGENVALUES = oEigVal, $
	                     MAXIMUM     = tf_maximum )
	
	;Rotate field MVA coordiantes
	;   - The rotation returns a MrVariable. Convert it to a MrVectorTS.
	oVec_mva  = oEigVec ## oVec
	oVec_mva = MrVectorTS( oVec['TIMEVAR'], oVec_mva, $
	                       /CACHE, $
	                       NAME = name )
	oVec_mva['COLOR'] = ['Blue', 'Forest Green', 'Red']
	oVec_mva['LABEL'] = param + '$\down' + ['N', 'M', 'L'] + '$'
	oVec_mva['TITLE'] = param + '!C(' + (param EQ 'B' ? 'nT' : 'mV/m') + ')'
	
;-------------------------------------------
; Plot Time Series /////////////////////////
;-------------------------------------------
	cgLoadCT, 13
	
	win = MrWindow( LAYOUT   = [1,2], $
	                OXMARGIN = [10,6], $
	                OYMARGIN = [20,3], $
	                XGAP     = 0.0, $
	                XSIZE    = 740, $
	                YGAP     = 0.5, $
	                YSIZE    = 700, $
	                REFRESH  = 0 )
	
	layout = MrLayout( [3,1], $
	                   ASPECT   = 1.0, $
	                   OXMARGIN = [10,6], $
	                   OYMARGIN = [4,25], $
	                   XGAP     = 15, $
	                   YGAP     = 0.5 )
	
	;XYZ
	p1 = MrVar_Plot( oVec, $
	                 /CURRENT, $
	                 LAYOUT = [1,1] )
	
	;LMN
	p2 = MrVar_Plot( oVec_mva, $
	                 /CURRENT, $
	                 LAYOUT = [1,2] )
	
;-------------------------------------------
; Plot Hodograms ///////////////////////////
;-------------------------------------------
	oVec_mva -> RemoveAttr, ['COLOR', 'LABEL']
	xyrange   = [oVec_mva.min, oVec_mva.max]
	
	;BL-BN
	p3 = MrVar_Plot( oVec_mva[it[0]:it[1], 0], oVec_mva[it[0]:it[1], 2], $
	                 /CURRENT, $
	                 POSITION = layout[*,0], $
	                 XRANGE   = xyrange, $
	                 YRANGE   = xyrange )
	x3a = MrPlotS( oVec_mva['DATA', it[0]:it[1], 0], oVec_mva['DATA', it[0]:it[1], 2], $
	               COLOR  = BytScl(oT['DATA', it[0]:it[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p3 )
	
	;BL-BM
	p4 = MrVar_Plot( oVec_mva[it[0]:it[1], 1], oVec_mva[it[0]:it[1], 2], $
	                 /CURRENT, $
	                 POSITION = layout[*,1], $
	                 XRANGE   = xyrange, $
	                 YRANGE   = xyrange )
	x4 = MrPlotS( oVec_mva['DATA', it[0]:it[1], 1], oVec_mva['DATA', it[0]:it[1], 2], $
	               COLOR  = BytScl(oT['DATA', it[0]:it[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p4 )
	
	;BM-BN
	p5 = MrVar_Plot( oVec_mva[it[0]:it[1], 0], oVec_mva[it[0]:it[1], 1], $
	                 /CURRENT, $
	                 POSITION = layout[*,2], $
	                 XRANGE   = xyrange, $
	                 YRANGE   = xyrange )
	x5 = MrPlotS( oVec_mva['DATA', it[0]:it[1], 0], oVec_mva['DATA', it[0]:it[1], 1], $
	               COLOR  = BytScl(oT['DATA', it[0]:it[1], 'SSM']), $
	               PSYM   = 3, $
	               TARGET = p5 )

;-------------------------------------------
; Annotate /////////////////////////////////
;-------------------------------------------
	tt = oT[ 'DATA', [it[0], it[1]], 'SSM' ]

	;Outline the MVA interval
	ps1 = MrPlotS( tt[[0,0]], p1.yrange, $
	               COLOR     = 'Blue', $
	               NAME      = 'Line: MVA0', $
	               TARGET    = p1, $
	               LINESTYLE = 2 )
	ps2 = MrPlotS( tt[[1,1]], p1.yrange, $
	               COLOR     = 'Blue', $
	               NAME      = 'Line: MVA1', $
	               TARGET    = p1, $
	               LINESTYLE = 2 )
	
	;Print eigenvectors and eigenvalues
	eigvec  = String( oEigVec['DATA'], FORMAT='(f7.4)' )
	seigvec = 'N = [' + StrJoin( eigvec[0:2], ', ' ) + ']' + '!C' + $
	          'M = [' + StrJoin( eigvec[3:5], ', ' ) + ']' + '!C' + $
	          'L = [' + StrJoin( eigvec[6:8], ', ' ) + ']'
	
	eigval  = String( oEigVal['DATA'], FORMAT='(f0.3)' )
	seigval = '\lambda = [' + StrJoin( eigval, ', ' ) + ']'
	
;	sB     = 'B = [' + StrJoin( String(b_avg, FORMAT='(f0.2)'), ', ' ) + ']'
;	sTheta = '\theta_{k} = ' + String(theta, FORMAT='(f0.2)')
	
	txt1 = MrText( 0.50, 0.14, seigval, ALIGNMENT=0.5, NAME='Txt: Eigenvalues',  /NORMAL )
	txt2 = MrText( 0.50, 0.11, seigvec, ALIGNMENT=0.5, NAME='Txt: Eigenvectors', /NORMAL )
;	txt3 = MrText( 0.50, 0.25, sB,      ALIGNMENT=0.5, NAME='Txt: B avg', /NORMAL )
;	txt3 = MrText( 0.50, 0.05, stheta,  ALIGNMENT=0.5, NAME='Txt: ThetaK', /NORMAL )

;-------------------------------------------
; Prettify /////////////////////////////////
;-------------------------------------------
	p1 -> SetProperty, XTICKFORMAT='(a1)', XTITLE=''
	p3 -> SetProperty, XTITLE=param+'$\downN$', YTITLE=param+'$\downL$'
	p4 -> SetProperty, XTITLE=param+'$\downM$', YTITLE=param+'$\downL$'
	p5 -> SetProperty, XTITLE=param+'$\downN$', YTITLE=param+'$\downM$'
	
	
	win[0] -> SetLayout, [1,1]
;	win    -> TrimLayout
	win -> SetProperty, OXMARGIN=[10,6]

;-------------------------------------------
; Save /////////////////////////////////////
;-------------------------------------------
	IF N_Elements(output_ext) GT 0 || N_Elements(output_dir) GT 0 THEN BEGIN
		IF N_Elements(output_ext) EQ 0 THEN output_ext = 'png'
		IF N_Elements(output_dir) EQ 0 THEN output_dir = File_Search('~', /TEST_DIRECTORY)
		
		;Parse analysis time
		t0 = StrJoin( StrSplit(t_mva[0], '-T:', /EXTRACT) )
		t1 = StrJoin( StrSplit(t_mva[1], '-T:', /EXTRACT) )
		date0 = StrMid(t0, 0, 8)
		date1 = StrMid(t1, 0, 8)
		time0 = StrMid(t0, 8, 6) + ( StrLen(t0) LE 14 ? '' : 'p' + StrMid(t0, 15) )
		time1 = StrMid(t1, 8, 6) + ( StrLen(t1) LE 14 ? '' : 'p' + StrMid(t1, 15) )
		
		;Time stamp for file
		IF date0 NE date1 $
			THEN ftime = StrJoin([date0, time0, date1, time1], '_') $
			ELSE ftime = StrJoin([date0, time0, time1], '_')
			
		;File name
		sinstr = tf_maximum ? 'edp' : fgm_instr
		fname  = StrJoin([sc, sinstr, mode, level, 'mva', ftime], '_')
		fname  = FilePath( fname + '.' + output_ext, $
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