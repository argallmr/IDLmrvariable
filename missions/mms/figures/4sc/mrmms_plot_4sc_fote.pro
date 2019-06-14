; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FGM_FOTE
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
;   Generate a plot summarizing the results of the First Order Taylor Expansion method.
;       1. Rx
;       2. Ry
;       3. Rz
;
; :Categories:
;   MMS
;
; :Params:
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;       INSTR:      in, required, type=string, default='fgm'
;                   FGM strument to use. Options are: {'afg' | 'dfg' | 'fgm'}
;
; :Keywords:
;       EPHDESC:    in, optional, type=string, default='ephts04d'
;                   Optional descriptor of the definitive ephemeris datatype to use.
;                       Options are: { 'epht89d' | 'epht89q' | 'ephts04d' | 'defeph' | 'predeph' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       OPTDESC:    in, optional, type=string, default=''
;                   Optional filename descriptor.
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
FUNCTION MrMMS_Plot_4sc_FOTE, mode, instr, $
EPHDESC=ephdesc, $
LEVEL=level, $
OPTDESC=optdesc, $
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
	
	RE2km = 1e-3 * MrConstants('RE')
	sc_colors = ['Black', 'Red', 'Forest Green', 'Blue']
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(trange)  GT 0 THEN MrVar_SetTRange, trange
	
	;INSTR
	IF N_Elements(instr) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN instr = 'fgm' $
			ELSE instr = level EQ 'ql' ? 'dfg' : 'fgm'
	ENDIF
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	
	;LEVEL
	IF N_Elements(level) EQ 0 THEN BEGIN
		CASE instr OF
			'afg': level = 'l2pre'
			'dfg': level = 'l2pre'
			'fgm': level = 'l2'
			ELSE: Message, 'Invalid FGM instrument: "' + instr + '".'
		ENDCASE
	ENDIF
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	IF N_Elements(coords) EQ 0 THEN BEGIN
		CASE level OF
			'ql': coords = 'dmpa'
			ELSE: coords = 'gse'
		ENDCASE
	ENDIF

	;Source names
	sc          = 'mms' + ['1', '2', '3', '4']
	b_vnames    = sc + '_' + StrJoin( [instr, 'b',    coords, mode, level], '_' )
	bvec_vnames = sc + '_' + StrJoin( [instr, 'bvec', coords, mode, level], '_' )
	bmag_vnames = sc + '_' + StrJoin( [instr, 'bmag', coords, mode, level], '_' )
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	r_interp_vnames = sc + '_' + StrJoin( ['mec', 'r', 'interp', coords], '_' )
	b1_vnames     = [ bmag_vnames[0], bvec_vnames[0] + '_' + ['x', 'y', 'z'] ]
	b2_vnames     = [ bmag_vnames[1], bvec_vnames[1] + '_' + ['x', 'y', 'z'] ]
	b3_vnames     = [ bmag_vnames[2], bvec_vnames[2] + '_' + ['x', 'y', 'z'] ]
	b4_vnames     = [ bmag_vnames[3], bvec_vnames[3] + '_' + ['x', 'y', 'z'] ]
	r_null_vname  = StrJoin( ['mms', instr, 'R', 'null', mode, level], '_' )
	rmag_null_vname = r_null_vname + '_magnitude'
	rmin_null_vname = r_null_vname + '_min'
	err_vname     = r_null_vname + '_uncertainty'
	rx_null_vname = r_null_vname + '_x'
	ry_null_vname = r_null_vname + '_y'
	rz_null_vname = r_null_vname + '_z'
	r_null_bary_vname  = StrJoin( ['mms', 'ephem', 'R', 'null', 'barycenter', mode, level], '_' )
	rx_null_bary_vname = r_null_bary_vname + '_x'
	ry_null_bary_vname = r_null_bary_vname + '_y'
	rz_null_bary_vname = r_null_bary_vname + '_z'
	pi_vname      = StrJoin( ['mms', instr, 'pi', mode, level], '_' )
	r_plus_vname  = StrJoin( ['mms', 'ephem', 'dr', 'plus',  mode, level], '_' )
	r_minus_vname = StrJoin( ['mms', 'ephem', 'dr', 'minus', mode, level], '_' )
	dr_mean_vname = StrJoin( ['mms', 'ephem', 'dr', 'mean',  mode, level], '_' )
	r_bary_vnames  = sc + StrJoin( ['mms', 'ephem', 'dr', 'bary',  mode, level], '_' )
	rx_bary_vnames = r_bary_vnames + '_x'
	ry_bary_vnames = r_bary_vnames + '_y'
	rz_bary_vnames = r_bary_vnames + '_z'
	
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		oFGM = MrMMS_FGM_4sc( mode, ephdesc, $
		                      COORD_SYS = coords, $
		                      INSTR     = instr, $
		                      LEVEL     = level )
	ENDIF
	
	;Old naming convention
	IF ~MrVar_IsCached(b_vnames[0]) THEN BEGIN
		b_vnames    = sc + '_' + StrJoin( [instr,        mode, coords], '_' )
		bvec_vnames = sc + '_' + StrJoin( [instr, 'vec', mode, coords], '_' )
		b1_vnames   = bvec_vnames[0] + '_' + ['x', 'y', 'z']
		b2_vnames   = bvec_vnames[1] + '_' + ['x', 'y', 'z']
		b3_vnames   = bvec_vnames[2] + '_' + ['x', 'y', 'z']
		b4_vnames   = bvec_vnames[3] + '_' + ['x', 'y', 'z']
		IF ~MrVar_IsCached(b_vnames[0]) $
			THEN Message, 'Unexpected variable naming convention for FGM.'
	ENDIF
	
	IF Obj_Valid(oFGM) EQ 0 THEN BEGIN
		oFGM = MrMMS_FGM_4sc( bvec_vnames[0], bvec_vnames[1], bvec_vnames[2], bvec_vnames[3], $
		                      r_vnames[0],    r_vnames[1],    r_vnames[2],    r_vnames[3] )
	ENDIF

;-------------------------------------------
; Null Position ////////////////////////////
;-------------------------------------------
	
	;FOTE
	oRnull = oFGM -> FOTE(/CACHE, ERR=oErr, NAME=r_null_vname, RMAG=oRmag)
	
	;Rmin
	Rmin  = Min(oRmag['DATA'], DIMENSION=2)
	oRmin = MrScalarTS(oRmag['TIMEVAR'], Rmin, /CACHE, NAME=rmin_null_vname)
	oRmin['AXIS_RANGE']   = [1, 1.5e3]
	oRmin['LOG']          = 1B
	oRmin['SYMBOL']       = oRmag['SYMBOL']
	oRmin['TICKINTERVAL'] = 500
	oRmin['TITLE']        = 'Rmin!C(km)'
	
	;Rmag
	oRmag['AXIS_RANGE'] = [1e0, 5e3]
	oRmag['LABEL']      = ['mms1', 'mms2', 'mms3', 'mms4']
	oRmag['LOG']        = 1B
	oRmag -> RemoveAttr, 'SYMBOL'
	
	;Rxyz
	oRnull -> Split, oRx_null, oRy_null, oRz_null, /CACHE
	
	;POINCARE INDEX
;	oPI = oFGM -> Poincare(NAME=pi_vname, /CACHE)

;-------------------------------------------
; Separation ///////////////////////////////
;-------------------------------------------
	oR1 = MrVar_Get(r_vnames[0])
	oR1 = oR1 -> Copy(r_interp_vnames[0], /CACHE)
	oT = oR1['TIMEVAR']
	FOR i = 1, N_Elements(sc) - 1 DO BEGIN
		oPos = MrVar_Get(r_vnames[i])
		oPos = oPos -> Interpol(oT, /CACHE, NAME=r_interp_vnames[i])
	ENDFOR
	
	;Position of barycenter (km)
	oR1 = MrVar_Get(r_interp_vnames[0])
	oR2 = MrVar_Get(r_interp_vnames[1])
	oR3 = MrVar_Get(r_interp_vnames[2])
	oR4 = MrVar_Get(r_interp_vnames[3])
	oRbary = (oR1 + oR2 + oR3 + oR4) / 4.0
	
	;Mean separation
	odR12 = (oR2 - oR1) -> Magnitude()
	odR13 = (oR2 - oR1) -> Magnitude()
	odR14 = (oR2 - oR1) -> Magnitude()
	odR23 = (oR3 - oR2) -> Magnitude()
	odR24 = (oR3 - oR2) -> Magnitude()
	odR34 = (oR4 - oR3) -> Magnitude()
	dr_mean = Mean( [ [odR12['DATA']], [odR12['DATA']], [odR12['DATA']], $
	                  [odR12['DATA']], [odR12['DATA']], [odR12['DATA']] ] )

;-------------------------------------------
; Barycenter Coordinates ///////////////////
;-------------------------------------------
	
	;Spacecraft position with respect to barycenter
	FOR i = 0, N_Elements(sc) - 1 DO BEGIN
		oPos     = MrVar_Get(r_interp_vnames[i])
		odR_bary = oPos - oRbary
		odR_bary -> SetName, r_bary_vnames[i]
		odR_bary -> Cache
		
		;Components
		odR_bary -> Split, odRx_bary, odRy_bary, odRz_bary, /CACHE
		
		;Attributes
		odRx_bary['COLOR'] = sc_colors[i]
		odRy_bary['COLOR'] = sc_colors[i]
		odRz_bary['COLOR'] = sc_colors[i]
	ENDFOR
		
	;Null
	oRbary = oRbary -> Interpol(oRnull)
	odRnull_bary = Re2km*oRnull - oRbary
	odRnull_bary -> SetName, r_null_bary_vname
	odRnull_bary -> Cache
	
	;Components
	rbc_range = 1.5 * [-dr_mean, dr_mean]
	odRnull_bary -> Split, odRxnull_bary, odRynull_bary, odRznull_bary, /CACHE
	odRxnull_bary['AXIS_RANGE']    = rbc_range
	odRxnull_bary['CATDESC']       = 'Distance to null from barycenter.'
	odRxnull_bary['LINESTYLE']     = 'None'
	odRxnull_bary['SYMBOL']        = 1
	odRxnull_bary['TITLE']         = 'Rx$\downBC$!C(km)'
	odRxnull_bary['UNITS']         = 'km'
	odRxnull_bary['SI_CONVERSION'] = '1e3>m'
	
	odRynull_bary['AXIS_RANGE']    = rbc_range
	odRynull_bary['CATDESC']       = 'Distance to null from barycenter.'
	odRynull_bary['LINESTYLE']     = 'None'
	odRynull_bary['SYMBOL']        = 1
	odRynull_bary['TITLE']         = 'Ry$\downBC$!C(km)'
	odRynull_bary['UNITS']         = 'km'
	odRynull_bary['SI_CONVERSION'] = '1e3>m'
	
	odRznull_bary['AXIS_RANGE']    = rbc_range
	odRznull_bary['CATDESC']       = 'Distance to null from barycenter.'
	odRznull_bary['LINESTYLE']     = 'None'
	odRznull_bary['SYMBOL']        = 1
	odRznull_bary['TITLE']         = 'Rz$\downBC$!C(km)'
	odRznull_bary['UNITS']         = 'km'
	odRznull_bary['SI_CONVERSION'] = '1e3>m'
	
	
;-------------------------------------------
; Ranges of the Spacecraft Positions //////
;-------------------------------------------
	
	;Deltas on barycenter
	r_plus = Max( [ [[oR1['DATA']]], [[oR2['DATA']]], [[oR3['DATA']]], [[oR4['DATA']]] ], $
	                DIMENSION = 3, $
	                MIN       = r_minus )
	
	;Create variables
	oR_plus = MrVectorTS( oR1['TIMEVAR'], r_plus / (1e-3 * MrConstants('RE')), $
	                      /CACHE, $
	                      NAME = r_plus_vname, $
	                      /NO_COPY )
	
	oR_minus = MrVectorTS( oR1['TIMEVAR'], r_minus / (1e-3 * MrConstants('RE')), $
	                       /CACHE, $
	                       NAME = r_minus_vname, $
	                       /NO_COPY )
	
	;RX
	rx_range = [oR_minus[*,0].min, oR_plus[*,0].max]
	oRx_null['AXIS_RANGE'] = rx_range + [-1,1] * 0.01 * Abs(rx_range)
	oRx_null['TITLE']      = 'Rx!C(RE)'
	
	;RY
	ry_range = [oR_minus[*,1].min, oR_plus[*,1].max]
	oRy_null['AXIS_RANGE'] = ry_range + [-1,1] * 0.01 * Abs(ry_range)
	oRy_null['TITLE']      = 'Ry!C(RE)'
	
	;RZ
	rz_range = [oR_minus[*,2].min, oR_plus[*,2].max]
	oRz_null['AXIS_RANGE'] = rz_range + [-1,1] * 0.01 * Abs(rz_range)
	oRz_null['TITLE']      = 'Rz!C(RE)'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [ rmag_null_vname, rx_null_vname, ry_null_vname, rz_null_vname, $
	                      rx_null_bary_vname, ry_null_bary_vname, rz_null_bary_vname, $
	                      rmin_null_vname, err_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	
	win = MrVar_OPlotTS(rx_null_bary_vname, rx_bary_vnames)
	win = MrVar_OPlotTS(ry_null_bary_vname, ry_bary_vnames)
	win = MrVar_OPlotTS(rz_null_bary_vname, rz_bary_vnames)
	
;	op = MrVar_OPlotTS( rx_null_vname, [oR_minus[*,0], oR_plus[*,0]] )
;	op = MrVar_OPlotTS( ry_null_vname, [oR_minus[*,1], oR_plus[*,1]] )
;	op = MrVar_OPlotTS( rz_null_vname, [oR_minus[*,2], oR_plus[*,2]] )
	
	s = MrPlotS( win[err_vname].xrange, [0.4, 0.4], $
	             LINESTYLE = '--', $
	             TARGET    = win[err_vname] )
	
	s = MrPlotS( win[rx_null_vname].xrange, rx_range[[0,0]], LINESTYLE='--', TARGET=win[rx_null_vname] )
	s = MrPlotS( win[rx_null_vname].xrange, rx_range[[1,1]], LINESTYLE='--', TARGET=win[rx_null_vname] )
	s = MrPlotS( win[ry_null_vname].xrange, ry_range[[0,0]], LINESTYLE='--', TARGET=win[ry_null_vname] )
	s = MrPlotS( win[ry_null_vname].xrange, ry_range[[1,1]], LINESTYLE='--', TARGET=win[ry_null_vname] )
	s = MrPlotS( win[rz_null_vname].xrange, rz_range[[0,0]], LINESTYLE='--', TARGET=win[rz_null_vname] )
	s = MrPlotS( win[rz_null_vname].xrange, rz_range[[1,1]], LINESTYLE='--', TARGET=win[rz_null_vname] )
	s = MrPlotS( win[rmin_null_vname].xrange, [dr_mean, dr_mean], LINESTYLE='--', TARGET=win[rmin_null_vname] )
	
	l = MrLegend( ALIGNMENT    = 'SE', $
	              LABEL        = ['A', 'B', 'As', 'Bs'], $
	              ORIENTATION  = 1, $
	              POSITION     = [1.0, 0.0], $
	              /RELATIVE, $
	              SAMPLE_WIDTH = 0.0, $
	              /SYM_CENTER, $
	              SYM_COLOR    = 'Black', $
	              SYMBOL       = [5, 12, 17, 19], $
	              TARGET       = win[rmin_null_vname], $
	              TEXT_COLOR   = 'Black' )
	
	;Pretty-up the window
	win    -> Refresh, /DISABLE
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[12, 9]
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
		fname   = StrJoin( ['mms1234', instr, mode, level, '4sc-fote'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END