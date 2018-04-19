; docformat = 'rst'
;
; NAME:
;       MrVar_Plot
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
;   Plot MrVariable data.
;
;   NOTES
;       For MrVariables with a 'DEPEND_0' attribute that has a MrTimeVar as its value,
;       the x-range is determined from the global time scale set by MrVar_SetTRange.
;       The time range can be obtained using
;
;           trange = MrVar_GetTRange()
;           tr_ssm = MrVar_GetTRange('SSM')
;
;       If the total time range is less than two seconds, time is converted to seconds
;       and displayed as seconds since Floor(tr_ssm[0]). For all other time intervals,
;       time is converted to seconds since midnight on trange[0] and displayed using the
;       XTICKFORMAT='time_labels'.
;
;   Calling Sequence:
;       p = MrVar_Plot(x)
;       p = MrVar_Plot(x, y)
;
; :Params:
;       X:          in, required, type=string/integer/objref
;                   A MrVariable name or index. If `Y` is given, then X is
;                       the independent variable, otherwise it is the dependent variable.
;                       In the latter case, X will be searched for a DEPEND_0 attribute.
;                       If present, it will be used as the dependent variable.
;       Y:          in, optional, type=string/integer/objref
;                   A MrVariable name, index, or object representing the dependent data
;                       to be plotted. 
;
; :Keywords:
;       BUFFER:     in, optional, type=boolean, default=0
;                   If set, graphics will be sent to the Z-buffer.
;       CURRENT:    in, optional, type=boolean, default=0
;                   If set, the plot will be added to the current MrGraphics window.
;       _REF_EXTRA: in, optional, type=any
;                   Any keyword accepted by MrPlot::SetProperty
;
; :Returns:
;       P:          out, required, type=object
;                   A MrPlot object reference.
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
;       2016/06/08  -   Written by Matthew Argall
;       2017/07/26  -   A time interval does not need to be defined (MrVar_Get/SetTRange) - MRA
;       2018/01/24  -   If X is time, then tick labels are formatted with t_ssm_labels.
;                           X-margins are made bigger if legend is created. - MRA
;       2018/01/31  -   Added the BUFFER keyword. - MRA
;-
function MrVar_Plot, x, y, $
BUFFER=buffer, $
CURRENT=current, $
DIMENSION=dimension, $
OVERPLOT=overplot, $
_REF_EXTRA=extra
	compile_opt strictarr

	;Catch errors
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return, !Null
	endif

;-------------------------------------------
; Get MrVariable Objects ///////////////////
;-------------------------------------------
	;Get the x-variable
	oX = size(x, /TNAME) eq 'OBJREF' ? x : MrVar_Get(x)
	
	;If Y was not given, look for DEPEND_0 in X.
	if n_elements(y) eq 0 then begin
		if oX -> HasAttr('DEPEND_0') then begin
			oY = oX
			oX = oY['DEPEND_0']
		endif
	endif else begin
		oY = size(y, /TNAME) eq 'OBJREF' ? y : MrVar_Get(y)
	endelse

;-------------------------------------------
; Setup Graphics ///////////////////////////
;-------------------------------------------
	
	;Get the window
	if n_elements(overplot) gt 0 then begin
		IF Keyword_Set(buffer) THEN MrPrintF, 'LogWarn', 'Using existing window. Ignoring BUFFER keyword.'
		win = overplot.window
		tf_refresh = win -> GetRefresh()
	endif else if keyword_set(current) then begin
		IF Keyword_Set(buffer) THEN MrPrintF, 'LogWarn', 'Using current window. Ignoring BUFFER keyword.'
		win        = GetMrWindows(/CURRENT)
		tf_refresh = win -> GetRefresh()
	endif else begin
		win        = MrWindow(REFRESH=0, BUFFER=buffer)
		tf_refresh = 1
	endelse
	win -> Refresh, /DISABLE
	
	;Keep track of margins
	oxmargin = win.oxmargin
	
	;Plot NaNs instead of fill values
	if oY -> HasAttr('FILLVAL') $
		then y_data = replace_fillval(oY['DATA'], oY['FILLVAL']) $
		else y_data = oY['DATA']

	;Use seconds since midnight, formatted as HH:MM:SS
	if obj_isa(oX, 'MrTimeVar') then begin
		x_data        = oX -> GetData('SSM')
		t_extra       = MrVar_TTicks(LABEL=t_label)
		xticks        = t_extra.xticks
		xminor        = t_extra.xminor
		xtickv        = t_extra.xtickv
		xrange        = t_extra.xrange
		xtickformat   = 't_ssm_ticks' ;time_labels
		xtitle        = '!CTime from ' + oX['DATA', 0, 0:9]
	endif else begin
		x_data = oX['DATA']
	endelse

;-------------------------------------------
; Display Data /////////////////////////////
;-------------------------------------------
	
	;Plot the data
	p1 = MrPlot( temporary(x_data), temporary(y_data), $
	             /CURRENT, $
	             DIMENSION   = oY -> GetAttrValue('DIMENSION',  /NULL), $
	             NAME        = oY.name, $
	             OVERPLOT    = overplot, $
	             COLOR       = oY -> GetAttrValue('COLOR',      /NULL), $
	             FONT        = oY -> GetAttrValue('FONT',       /NULL), $
	             LINESTYLE   = oY -> GetAttrValue('LINESTYLE',  /NULL), $
	             MAX_VALUE   = oY -> GetAttrValue('MAX_VALUE',  /NULL), $
	             MIN_VALUE   = oY -> GetAttrValue('MIN_VALUE',  /NULL), $
	             NOCLIP      = oY -> GetAttrValue('NOCLIP',     /NULL), $
	             NSUM        = oY -> GetAttrValue('NSUM',       /NULL), $
	             POLAR       = oY -> GetAttrValue('POLAR',      /NULL), $
	             PSYM        = oY -> GetAttrValue('SYMBOL',     /NULL), $
	             SYMCOLOR    = oY -> GetAttrValue('SYM_COLOR',  /NULL), $
	             SYMSIZE     = oY -> GetAttrValue('SYM_SIZE',   /NULL), $
	             THICK       = oY -> GetAttrValue('THICK',      /NULL), $
	             TICKLEN     = oY -> GetAttrValue('TICKLEN',    /NULL), $
	             TITLE       = oY -> GetAttrValue('PLOT_TITLE', /NULL) )
	
	;Set graphics properties
	MrVar_SetAxisProperties, p1, oX, /XAXIS
	MrVar_SetAxisProperties, p1, oY, /YAXIS
	
	;Set user-given properties
	p1 -> SetProperty, XTICKFORMAT=xtickformat, XMINOR=xminor, XRANGE=xrange, XTICKS=xticks, XTICKV=xtickv, XTITLE=xtitle
;	p1 -> SetProperty, XRANGE=xrange, XTITLE=xtitle, XTICKFORMAT=xtickformat
	if n_elements(extra) gt 0 then p1 -> SetProperty, _STRICT_EXTRA=extra

;-------------------------------------------
; Create a Legend //////////////////////////
;-------------------------------------------
	if oY -> HasAttr('LABEL') then begin
		;Do not draw the lines
		if ~oY -> HasAttr('SAMPLE_WIDTH') then oY['SAMPLE_WIDTH'] = 0
		
		;Create the legend
		if n_elements(overplot) eq 0 then begin
			l1 = MrVar_Legend( oY, $
			                   COLOR      = '', $
			                   FILL_COLOR = '', $
			                   LINESTYLE  = 'None', $
			                   NAME       = 'Lgd: ' + oY.name, $
			                   TARGET     = p1 )
			                   
			;Make the legend a bit bigger
			oxmargin[1] >= 6
		
		;Add a legend item
		endif else begin
			l1 = MrVar_Legend( oY, $
			                   /ADD, $
			                   COLOR      = '', $
			                   FILL_COLOR = '', $
			                   LINESTYLE  = 'NONE', $
			                   NAME       = 'Lgd: ' + oY.name, $
			                   TARGET     = overplot )
		endelse
	endif

;-------------------------------------------
; Finish ///////////////////////////////////
;-------------------------------------------
	;Update the margins
	IF ~Array_Equal(win.oxmargin, oxmargin) $
		THEN win.oxmargin >= oxmargin
	
	;Return the plot
	if tf_refresh then p1 -> Refresh
	return, p1
end