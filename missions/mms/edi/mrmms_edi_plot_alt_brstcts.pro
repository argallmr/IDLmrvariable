; docformat = 'rst'
;
; NAME:
;       MrMMS_EDI_Plot_ALT
;
;*****************************************************************************************
;   Copyright (c) 2016, Matthew Argall                                                   ;
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
;   Load EDI alternating-mode data (as opposed to field-aligned and perpendicular).
;   Plot Survey Flux/Counts on top of Burst Flux/Counts to ensure calibration is
;   applied correctly.
;       1.1 Flux 0-PA CH1   2.1 Flux 90-PA GDU1 CH1  3.1 Flux 90-PA GDU2 CH1  4.1 Flux 180-PA CH1
;       1.2 Flux 0-PA CH2   2.2 Flux 90-PA GDU1 CH2  3.2 Flux 90-PA GDU2 CH2  4.2 Flux 180-PA CH2
;       1.3 Flux 0-PA CH3   2.3 Flux 90-PA GDU1 CH3  3.3 Flux 90-PA GDU2 CH3  4.2 Flux 180-PA CH3
;       1.4 Flux 0-PA CH4   2.4 Flux 90-PA GDU1 CH4  3.4 Flux 90-PA GDU2 CH4  4.4 Flux 180-PA CH4
;
; :Categories:
;       MMS, EDI, MrVariable
;
; :Params:
;       SC:                 in, required, type=string/strarr
;                           The MMS spacecraft identifier. Options are:
;                               {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       OPTDESC:            in, optional, type=string, default=''
;                           Optional descriptor of the data. Options are:
;                               {'amb-alt-cc' | 'amb-alt-oc' | 'amb-alt-oom' | 'amb-alt-oob'}
;
; :Keywords:
;       NO_LOAD:            in, optional, type=boolean, default=0
;                           If set, data is not re-read and loaded into the variable cache.
;       TRANGE:             in, optional, type=strarr(2), default=MrVar_GetTRange
;                           The start and end times of the data interval to be loaded.
;                               Formatting is: 'YYYY-MM-DDThh:mm:ss'.
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
;       2016/10/03  -   Written by Matthew Argall
;-
function MrMMS_EDI_Plot_ALT_BrstCts, sc, optdesc, $
NO_LOAD=no_load, $
TRANGE=trange
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if n_elements(win) gt 0 then obj_destroy, win
		MrPrintF, 'LogErr'
		return, obj_new()
	endif
	
	;Defaults
	tf_load = ~keyword_set(no_load)
	if n_elements(trange) gt 0 then MrVar_SetTRange, trange
	if ~MrIsMember(['amb-alt-cc', 'amb-alt-oc', 'amb-alt-oom', 'amb-alt-oob'], optdesc) $
		then message, 'Invalid value for OPTDESC.'
	
	;Variable type
	level = 'l2'
	type  = level eq 'l2' ? 'flux' : 'counts'

;-------------------------------------------
; EDI Data /////////////////////////////////
;-------------------------------------------
	;Get EDI data
	if tf_load then begin
		;Load EDI data
		MrMMS_Load_Data, sc, 'edi', ['brst', 'srvy'], level, $
		                 OPTDESC   = optdesc, $
		                 VARFORMAT = '*' + type + '*'
	endif
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Parts
	prefix      = sc + '_edi_'
	suffix_srvy = '_srvy_' + level
	suffix_brst = '_brst_' + level
	ch          = ['1', '2', '3', '4']
	
	;Names
	cts_srvy_vname = [ prefix + type   + ch[0] + '_0'       + suffix_srvy, $
	                   prefix + type   + ch[0] + '_90_gdu1' + suffix_srvy, $
	                   prefix + type   + ch[0] + '_90_gdu2' + suffix_srvy, $
	                   prefix + type   + ch[0] + '_180'     + suffix_srvy ]
	cts_brst_vname = [ [ prefix + type   + ch  + '_0'       + suffix_brst ], $
	                   [ prefix + type   + ch  + '_90_gdu1' + suffix_brst ], $
	                   [ prefix + type   + ch  + '_90_gdu2' + suffix_brst ], $
	                   [ prefix + type   + ch  + '_180'     + suffix_brst ] ]
	cts_brst_vname = reform(transpose(cts_brst_vname), 16)

;-------------------------------------------
; Set Properties ///////////////////////////
;-------------------------------------------
	if tf_load then begin
		;BRST
		for i = 0, n_elements(cts_brst_vname) - 2 do begin
			oVar = MrVar_Get( cts_brst_vname[i] )
			oVar -> AddAttr, 'AXIS_RANGE', [1, max(oVar['DATA'])]
		endfor
	
		;SRVY
		for i = 0, n_elements(cts_srvy_vname) - 1 do begin
			oVar = MrVar_Get( cts_srvy_vname[i] )
			oVar -> AddAttr, 'COLOR', 'Red'
		endfor
	endif

;-------------------------------------------
; First Row ////////////////////////////////
;-------------------------------------------
	;Plot burst data
	win = MrVar_PlotTS( cts_brst_vname, LAYOUT=[4,4], XSIZE=1000, YSIZE=600 )
	
	;Overplot survey data
	cts_brst_vname = reform( cts_brst_vname, 4, 4 )

	win = MrVar_OPlotTS( cts_brst_vname[*,0], cts_srvy_vname )
;	win = MrVar_OPlotTS( cts_brst_vname[*,1], cts_srvy_vname )
;	win = MrVar_OPlotTS( cts_brst_vname[*,2], cts_srvy_vname )
;	win = MrVar_OPlotTS( cts_brst_vname[*,3], cts_srvy_vname )
	
	;Create a legend
;	lgd = MrLegend( ALIGNMENT    = 'NW', $
;	                LABEL        = ['Brst', 'Srvy'], $
;	                LINESTYLE    = 6, $
;	                POSITION     = [1.0, 1.0], $
;	                /RELATIVE, $
;	                SAMPLE_WIDTH = 0.0, $
;	                TARGET       = win[cts_brst_vname[0,3]], $
;	                TEXT_COLOR   = ['Black', 'Red'] )
	
	win -> SetGlobal, YRANGE=[1,2000]
	win -> Refresh
	return, win
end