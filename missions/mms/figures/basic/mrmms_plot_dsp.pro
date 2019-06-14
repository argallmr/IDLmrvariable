; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_DSP
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
;   Generate a plot of FGM quantities:
;       1. Bx, By, Bz, |B|
;       2. Bx PSD + fce, fce/2 lines
;       3. By PSD + fce, fce/2 lines
;       4. Bz PSD + fce, fce/2 lines
;       5. B Omni PSD + fce, fce/2 lines
;       6. Ex PSD + fce, fce/2 lines
;       7. Ey PSD + fce, fce/2 lines
;       8. Ez PSD + fce, fce/2 lines
;       9. E Omni PSD + fce, fce/2 lines
;
; :Categories:
;   MMS
;
; :Params:
;       MODE:       in, required, type=string, default='fast'
;                   Data rate mode. Options are {'slow' | 'fast'}
;       INSTR:      in, required, type=string, default='fgm'
;                   FGM strument to use. Options are: {'afg' | 'dfg' | 'fgm'}
;
; :Keywords:
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
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
;       2017/01/17  -   Written by Matthew Argall
;       2018/10/10  -   Added the TEAM_SITE keyword. Fast and slow survey data are
;                           plotted together of MODE='srvy'. - MRA
;-
FUNCTION MrMMS_Plot_DSP, sc, mode, $
FGM_INSTR=fgm_instr, $
NO_LOAD=no_load, $
TEAM_SITE=team_site, $
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
	IF N_Elements(level)     EQ 0 THEN level      = 'l2'
	IF N_Elements(mode)      EQ 0 THEN mode       = 'fast'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr  = 'fgm'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	;DSP
	dsp_instr = 'dsp'
	dsp_mode  = mode EQ 'srvy' ? ['fast', 'slow'] : mode
		
	;FGM
	fgm_mode = 'srvy'
	CASE fgm_instr OF
		'afg': fgm_level = 'l2pre'
		'dfg': fgm_level = 'l2pre'
		'fgm': fgm_level = 'l2'
		ELSE: Message, 'Invalid FGM instrument: "' + instr + '".'
	ENDCASE
	CASE fgm_level OF
		'ql': fgm_coords = 'dmpa'
		ELSE: fgm_coords = 'gse'
	ENDCASE
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------

	;Source names
	b_vname    = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_' )
	bmag_vname = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_' )
	
	bxpsd_vname = StrJoin( [sc, dsp_instr, 'bpsd', 'scm1', dsp_mode[0], level], '_' )
	bypsd_vname = StrJoin( [sc, dsp_instr, 'bpsd', 'scm2', dsp_mode[0], level], '_' )
	bzpsd_vname = StrJoin( [sc, dsp_instr, 'bpsd', 'scm3', dsp_mode[0], level], '_' )
	bpsd_vname  = StrJoin( [sc, dsp_instr, 'bpsd', 'omni', dsp_mode[0], level], '_' )
	
	expsd_vname = StrJoin( [sc, dsp_instr, 'epsd', 'x', dsp_mode[0], level], '_' )
	eypsd_vname = StrJoin( [sc, dsp_instr, 'epsd', 'y', dsp_mode[0], level], '_' )
	ezpsd_vname = StrJoin( [sc, dsp_instr, 'epsd', 'z', dsp_mode[0], level], '_' )
	epsd_vname  = StrJoin( [sc, dsp_instr, 'epsd', 'omni', dsp_mode[0], level], '_' )
	
	IF mode EQ 'srvy' THEN BEGIN
		bxpsd_slow_vname = StrJoin( [sc, dsp_instr, 'bpsd', 'scm1', dsp_mode[1], level], '_' )
		bypsd_slow_vname = StrJoin( [sc, dsp_instr, 'bpsd', 'scm2', dsp_mode[1], level], '_' )
		bzpsd_slow_vname = StrJoin( [sc, dsp_instr, 'bpsd', 'scm3', dsp_mode[1], level], '_' )
		bpsd_slow_vname  = StrJoin( [sc, dsp_instr, 'bpsd', 'omni', dsp_mode[1], level], '_' )
	
		expsd_slow_vname = StrJoin( [sc, dsp_instr, 'epsd', 'x', dsp_mode[1], level], '_' )
		eypsd_slow_vname = StrJoin( [sc, dsp_instr, 'epsd', 'y', dsp_mode[1], level], '_' )
		ezpsd_slow_vname = StrJoin( [sc, dsp_instr, 'epsd', 'z', dsp_mode[1], level], '_' )
		epsd_slow_vname  = StrJoin( [sc, dsp_instr, 'epsd', 'omni', dsp_mode[1], level], '_' )
	ENDIF
	
	;Derived names
	fce_vname   = StrJoin( [sc, dsp_instr, 'fce',      mode, level], '_' )
	fce_2_vname = StrJoin( [sc, dsp_instr, 'halffce',  mode, level], '_' )

;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = fgm_level, $
		                     TEAM_SITE = team_site, $
		                     VARFORMAT = b_vname
		
		;BPSD
		MrMMS_Load_Data, sc, dsp_instr, dsp_mode, level, $
		                 OPTDESC   = 'bpsd', $
		                 TEAM_SITE = team_site, $
		                 VARFORMAT = ['*bpsd_scm*', '*bpsd_omni*']
		
		;ESPD
		MrMMS_Load_Data, sc, dsp_instr, dsp_mode[0], level, $
		                 OPTDESC   = 'epsd', $
		                 SUFFIX    = '_' + dsp_mode[0] + '_' + level, $
		                 TEAM_SITE = team_site, $
		                 VARFORMAT = '*epsd*'
		 
		IF mode EQ 'srvy' THEN BEGIN
			MrMMS_Load_Data, sc, dsp_instr, dsp_mode[1], level, $
			                 OPTDESC   = 'epsd', $
			                 SUFFIX    = '_' + dsp_mode[1] + '_' + level, $
			                 TEAM_SITE = team_site, $
			                 VARFORMAT = '*epsd*'
		ENDIF
	ENDIF
	
;-------------------------------------------
; Gyrofrequency Lines //////////////////////
;-------------------------------------------
	tf_fce = MrVar_IsCached(b_vname)
	IF tf_fce THEN BEGIN
		oBmag  = MrVar_Get(bmag_vname)
		ofce   = MrVar_Freq_Cyclotron(oBmag, 'm_e',  /CACHE, NAME=fce_vname)
		ofce_2 = ofce / 2.0
		ofce_2 -> SetName, fce_2_vname
		ofce_2 -> Cache
		
		;BMAG
		oB = MrVar_Get(b_vname)
		oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, instr, mode, level], ' ' ) )
	
		;FCE
		ofce['COLOR']     = 'White'
		ofce['LINESTYLE'] = '--'
		ofce['NSUM']      = 4
	
		;FCE_2
		ofce_2['COLOR']      = 'Magenta'
		ofce_2['LINESTYLE'] = '.'
		ofce_2['NSUM']      = 4
		ofce_2['UNITS']     = 'Hz'
	ENDIF
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	
	;FREQUENCY
	;   - DELTA_(MINUS|PLUS) does not have the correct number of elements.
	oBx = MrVar_Get(bxpsd_vname)
	bFreq = oBx['DEPEND_1']
	bFreq['TITLE'] = 'Freq!C(Hz)'
	IF bFreq -> HasAttr('DELTA_MINUS') THEN bFreq -> RemoveAttr, 'DELTA_MINUS'
	IF bFreq -> HasAttr('DELTA_PLUS')  THEN bFreq -> RemoveAttr, 'DELTA_PLUS'
	
	;BX PSD
	oBx['LOG']   = 1
	oBx['TITLE'] = 'PSD!CBx'
	
	;BY PSD
	oBy = MrVar_Get(bypsd_vname)
	oBy['LOG']        = 1
	oBy['TITLE'] = 'PSD!CBy'

	;BZ PSD
	oBz = MrVar_Get(bzpsd_vname)
	oBz['LOG']   = 1
	oBz['TITLE'] = 'PSD!CBz!C(nT$\up2$/Hz)'
	
	;B PSD
	oBo = MrVar_Get(bpsd_vname)
	oBo['LOG']   = 1
	oBo['TITLE'] = 'PSD!CB Omni'
	
	;FREQUENCY
	oEx = MrVar_Get(expsd_vname)
	eFreq = oEx['DEPEND_1']
	eFreq['TITLE'] = 'Freq!C(Hz)'
	
	;EX PSD
	oEx['AXIS_RANGE'] = [1e-15, 1e-10]
	oEx['LOG']        = 1
	oEx['TITLE']      = 'PSD!CEx'
	
	;EY PSD
	oEy = MrVar_Get(eypsd_vname)
	oEy['AXIS_RANGE'] = [1e-15, 1e-10]
	oEy['LOG']        = 1
	oEy['TITLE']      = 'PSD!CEy''
	
	;EZ PSD
	oEz = MrVar_Get(ezpsd_vname)
	oEz['AXIS_RANGE'] = [1e-15, 1e-10]
	oEz['LOG']        = 1
	oEz['TITLE']      = 'PSD!CEz!C(mV/m)$\up2$/Hz'
	
	;E PSD
	oEo = MrVar_Get(epsd_vname)
	oEo['AXIS_RANGE'] = [1e-15, 1e-10]
	oEo['LOG']        = 1
	oEo['TITLE']     = 'PSD!CE Omni'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [b_vname, bxpsd_vname, bypsd_vname, bzpsd_vname, bpsd_vname, $
	                     expsd_vname, eypsd_vname, ezpsd_vname, epsd_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 680, $
	                    YSIZE = 700 )
	
	;Slow survey
	IF mode EQ 'srvy' THEN BEGIN
		win = MrVar_OPlotTS( bxpsd_vname, bxpsd_slow_vname )
		win = MrVar_OPlotTS( bypsd_vname, bypsd_slow_vname )
		win = MrVar_OPlotTS( bzpsd_vname, bzpsd_slow_vname )
		win = MrVar_OPlotTS( bpsd_vname,  bpsd_slow_vname )
		
		win = MrVar_OPlotTS( expsd_vname, expsd_slow_vname )
		win = MrVar_OPlotTS( eypsd_vname, eypsd_slow_vname )
		win = MrVar_OPlotTS( ezpsd_vname, ezpsd_slow_vname )
		win = MrVar_OPlotTS( epsd_vname,  epsd_slow_vname )
	ENDIF
	
	;Gyrofrequency lines
	IF tf_fce THEN BEGIN
		win = MrVar_OPlotTS( bxpsd_vname, [fce_vname, fce_2_vname] )
		win = MrVar_OPlotTS( bypsd_vname, [fce_vname, fce_2_vname] )
		win = MrVar_OPlotTS( bzpsd_vname, [fce_vname, fce_2_vname] )
		win = MrVar_OPlotTS( bpsd_vname,  [fce_vname, fce_2_vname] )
		win = MrVar_OPlotTS( expsd_vname, [fce_vname, fce_2_vname] )
		win = MrVar_OPlotTS( eypsd_vname, [fce_vname, fce_2_vname] )
		win = MrVar_OPlotTS( ezpsd_vname, [fce_vname, fce_2_vname] )
		win = MrVar_OPlotTS( epsd_vname,  [fce_vname, fce_2_vname] )
	ENDIF

	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 14]
	win    -> Refresh
	
	RETURN, win
END