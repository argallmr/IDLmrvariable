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
FUNCTION MrMMS_Plot_Formation, mode, $
INSTR=instr, $
COORDS=coords, $
NO_LOAD=no_load, $
OPTDESC=optdesc, $
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
	
	tf_load   = ~Keyword_Set(no_load)
	eph_instr = N_Elements(instr) EQ 0 ? 'mec' : instr
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(optdesc)   EQ 0 THEN optdesc   = 'epht89d'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
;-------------------------------------------
; Data Parameters //////////////////////////
;-------------------------------------------
	eph_instr = 'mec'
	eph_mode  = mode
		
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	sc = 'mms' + ['1', '2', '3', '4']
	IF Array_Equal(eph_instr EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(optdesc) + '_' + 'R'
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;MEC
		IF eph_instr EQ 'mec' THEN BEGIN
			;For part of the mission (2015-2016), burst mode is not available
			IF ~!MrMMS.offline && eph_mode EQ 'brst' THEN BEGIN
				IF (mrmms_get_filenames('mms1', mec_instr, eph_mode, level, OPTDESC=optdesc))[0] EQ '' $
					THEN eph_mode = 'srvy'
			ENDIF
			
			MrMMS_Load_Data, '', eph_instr, eph_mode, level, $
			                 OPTDESC   = optdesc, $
			                 VARFORMAT = r_vnames
		
		;DEFEPH & PREDEPH
		ENDIF ELSE BEGIN
			Message, 'Definitive and Predictive ephemeris have not been implemented.'
		ENDELSE
	ENDIF
	
;-------------------------------------------
; Separation ///////////////////////////////
;-------------------------------------------
	oR1 = MrVar_Get(r_vnames[0])
	oR2 = MrVar_Get(r_vnames[1])
	oR3 = MrVar_Get(r_vnames[2])
	oR4 = MrVar_Get(r_vnames[3])
	
	;First, find the barycenter
	oR_bary = (oR1 + oR2 + oR3 + oR4) / 4.0
	
	;Next, find separation relative to barycenter
	oR1b = oR1 - oR_bary
	oR2b = oR2 - oR_bary
	oR3b = oR3 - oR_bary
	oR4b = oR4 - oR_bary
	
	;Finally, take the mean position
	r1 = Reform( Mean(oR1b['DATA'], DIMENSION=1) )
	r2 = Reform( Mean(oR2b['DATA'], DIMENSION=1) )
	r3 = Reform( Mean(oR3b['DATA'], DIMENSION=1) )
	r4 = Reform( Mean(oR4b['DATA'], DIMENSION=1) )
	
;-------------------------------------------
; Plot the Positions ///////////////////////
;-------------------------------------------
	win = MrWindow( ASPECT   = 1.0, $
	                LAYOUT   = [3,1], $
	                OXMARGIN = [10,8], $
	                REFRESH  = 0, $
	                XGAP     = 6, $
	                XSIZE    = 1000, $
	                YGAP     = 0, $
	                YSIZE    = 300 )
	                
	pos = MrLayout( [3,1], $
	                ASPECT   = 1.0, $
	                CHARSIZE = 2.0, $
	                OXMARGIN = [10,8], $
	                WDIMS    = [1000,300], $
	                XGAP     = 6, $
	                YGAP     = 0 )
	
	xrange = [ Min([r1[0], r2[0], r3[0], r4[0]]), Max([r1[0], r2[0], r3[0], r4[0]]) ]
	yrange = [ Min([r1[1], r2[1], r3[1], r4[1]]), Max([r1[1], r2[1], r3[1], r4[1]]) ]
	zrange = [ Min([r1[2], r2[2], r3[2], r4[2]]), Max([r1[2], r2[2], r3[2], r4[2]]) ]
	
	;XY
	range  = [ Min([xrange[0], yrange[0]]), Max([xrange[1], yrange[1]]) ] * 1.15
	zvalue = 3 * ([r1[2], r2[2], r3[2], r4[2]] + abs(zrange[0])) / (zrange[1] - zrange[0]) + 1
	
	p1 = MrPlot( r1[0], r1[1], $
	             /CURRENT, $
	             COLOR    = 'Black', $
;	             LAYOUT   = [1,1], $
	             NAME     = 'XY MMS1', $
	             POSITION = pos[*,1], $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[0], $
	             TITLE    = title, $
	             XRANGE   = Reverse(range), $
	             XTITLE   = '$\Delta$X (km)', $
	             YRANGE   = range, $
	             YTITLE   = '$\Delta$Y (km)' )

	p2 = MrPlot( r2[0], r2[1], $
	             COLOR    = 'Red', $
	             NAME     = 'XY MMS2', $
	             OVERPLOT = p1, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[1] )
	
	p3 = MrPlot( r3[0], r3[1], $
	             COLOR    = 'Green', $
	             NAME     = 'XY MMS3', $
	             OVERPLOT = p1, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[2] )
	
	p4 = MrPlot( r4[0], r4[1], $
	             COLOR    = 'Blue', $
	             NAME     = 'XY MMS4', $
	             OVERPLOT = p1, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[3])
	
	;XZ
	range  = [ Min([xrange[0], zrange[0]]), Max([xrange[1], zrange[1]]) ] * 1.15
	zvalue = 3 * ([r1[1], r2[1], r3[1], r4[1]] + abs(yrange[0])) / (yrange[1] - yrange[0]) + 1
	
	p5 = MrPlot( r1[0], r1[2], $
	             /CURRENT, $
	             COLOR    = 'Black', $
;	             LAYOUT   = [1,1], $
	             NAME     = 'XZ MMS1', $
	             POSITION = pos[*,0], $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[0], $
	             XRANGE   = Reverse(range), $
	             XTITLE   = '$\Delta$X (km)', $
	             YRANGE   = range, $
	             YTITLE   = '$\Delta$Z (km)' )

	p6 = MrPlot( r2[0], r2[2], $
	             COLOR    = 'Red', $
	             NAME     = 'XZ MMS2', $
	             OVERPLOT = p5, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[1] )
	
	p7 = MrPlot( r3[0], r3[2], $
	             COLOR    = 'Green', $
	             NAME     = 'XZ MMS3', $
	             OVERPLOT = p5, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[2] )
	
	p8 = MrPlot( r4[0], r4[2], $
	             COLOR    = 'Blue', $
	             NAME     = 'XZ MMS4', $
	             OVERPLOT = p5, $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[3] )
	
	;YZ
	range  = [ Min([yrange[0], zrange[0]]), Max([yrange[1], zrange[1]]) ] * 1.15
	zvalue = 3 * ([r1[0], r2[0], r3[0], r4[0]] + abs(xrange[0])) / (xrange[1] - xrange[0]) + 1
	
	p9 = MrPlot( r1[1], r1[2], $
	             /CURRENT, $
	             COLOR    = 'Black', $
;	             LAYOUT   = [1,1], $
	             NAME     = 'YZ MMS1', $
	             POSITION = pos[*,2], $
	             PSYM     = 'FilledCircle', $
	             SYMSIZE  = zvalue[0], $
	             XRANGE   = range, $
	             XTITLE   = '$\Delta$Y (km)', $
	             YRANGE   = range, $
	             YTITLE   = '$\Delta$Z (km)' )

	p10 = MrPlot( r2[1], r2[2], $
	              COLOR    = 'Red', $
	              NAME     = 'YZ MMS2', $
	              OVERPLOT = p9, $
	              PSYM     = 'FilledCircle', $
	              SYMSIZE  = zvalue[1], $
	              XRANGE   = range, $
	              YRANGE   = range )
	
	p11 = MrPlot( r3[1], r3[2], $
	              COLOR    = 'Green', $
	              NAME     = 'YZ MMS3', $
	              OVERPLOT = p9, $
	              PSYM     = 'FilledCircle', $
	              SYMSIZE  = zvalue[2], $
	              XRANGE   = range, $
	              YRANGE   = range )
	
	p12 = MrPlot( r4[1], r4[2], $
	              COLOR    = 'Blue', $
	              NAME     = 'YZ MMS4', $
	              OVERPLOT = p9, $
	              PSYM     = 'FilledCircle', $
	              SYMSIZE  = zvalue[3], $
	              XRANGE   = range, $
	              YRANGE   = range )
	
	;Legend
	l1 = MrLegend( ALIGNMENT    = 'NW', $
	               FILL_COLOR   = '', $
	               LABEL        = 'MMS' + ['1', '2', '3', '4'], $
	               LINESTYLE    = 'None', $
	               POSITION     = [1.0, 1.0], $
	               /RELATIVE, $
	               SAMPLE_WIDTH = 0, $
	               SYMBOL       = 'FilledCircle', $
	               /SYM_CENTER, $
	               SYM_COLOR    = ['Black', 'Red', 'Green', 'Blue'], $
	               TARGET       = p9, $
	               TEXT_COLOR   = ['Black', 'Red', 'Green', 'Blue'] )
	
	;Title
	tt     = MrVar_GetTRange()
	tt_ssm = MrVar_GetTRange('SSM')
	tt     = StrMid(tt[0], 0, 10) + ' ' + StrMid(ssm_to_hms(Mean(tt_ssm)), 0, 8)
	title  = 'MMS Tetrahedron ' + tt
	t1     = MrText( 0.5, 0.91, title, $
	                 ALIGNMENT = 0.5, $
	                 CHARSIZE  = 2.0, $
	                 /CURRENT, $
	                 /NORMAL )

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
		fname   = StrJoin( ['mms', eph_instr, eph_mode, level, 'formation'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	win -> Refresh

	RETURN, win
END