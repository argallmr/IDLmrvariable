; docformat = 'rst'
;
; NAME:
;       MrVar_PlotTS
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
;   Plot MrVariable data
;
;   Calling Sequence:
;       win = MrVar_PlotTS()
;       win = MrVar_PlotTS(variables)
;
; :Params:
;       VARIABLES:  in, optional, type=strarr/intarr
;                   The names or indices of time series variables to be plotted. If
;                       not provided, the variables given in the previous call to
;                       MrVar_PlotTS will be used.
;
; :Params:
;       BUFFER:     in, optional, type=boolean, default=0
;                   If set, graphics will be sent to the Z-buffer.
;       LAYOUT:     in, optional, type=intarr(2), default=[1\,nVars]
;                   The number of columns and rows in the plot layout: [nCols, nRows].
;                       `VARIABLES` will plotted from right to left then from top to
;                       bottom. 
;       NO_REFRESH: in, optional, type=boolean, default=0
;                   If set, then `WIN` will not be refreshed before being returned. This
;                       means that the GUI window itself will not be realized and the
;                       graphics will not be drawn. This is useful if more graphics are
;                       to be added or graphics properties are to be changed. To refresh
;                       the window, call its ::Refresh method.
;       XSIZE:      in, optional, type=integer, default=640
;                   Width of the graphics window, in pixels.
;       YSIZE:      in, optional, type=integer, default=500
;                   Height of the graphics window, in pixels.
;
; :Returns:
;       WIN:        out, required, type=objref
;                   A MrWindow object reference.
;
; :Common Blocks:
;       MrVar_PlotTS_Common
;           Retains the names of the plotted variable names so that variable properties
;           can be updated and the plot re-invoked.
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
;       2016/08/13  -   Written by Matthew Argall
;       2016/10/03  -   Added the LAYOUT and NO_REFRESH keywords. - MRA
;       2017/03/09  -   Added the CURRENT keyword. - MRA
;       2017/08/02  -   VARIABLES can be MrTimeSeries objects. - MRA
;       2018/01/22  -   Make use of MrVar_TTicks to format x-tick labels. Adjust x-margins
;                           if a time-label is added to the plot. - MRA
;       2018/01/31  -   Added the BUFFER keyword. - MRA
;-
function MrVar_PlotTS, variables, $
BUFFER=buffer, $
CURRENT=current, $
LAYOUT=layout, $
NO_REFRESH=no_refresh, $
XSIZE=xsize, $
YSIZE=ysize
	compile_opt strictarr

	;Catch errors
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if n_elements(win) gt 0 then obj_destroy, win
		MrPrintF, 'LogErr'
		return, !Null
	endif
	
	;Common block to hold previous variables
	common MrVar_PlotTS_Common, vars

;-------------------------------------------
; What to Plot /////////////////////////////
;-------------------------------------------
	;Use what was plotted previously
	if n_elements(variables) eq 0 then begin
		if n_elements(vars) eq 0 $
			then message, 'Specify the variables to be plotted.' $
			else variables = vars
	endif
	
	;Get all of the variables
	allVars = MrVar_Get(variables, COUNT=nVars)

;-------------------------------------------
; Check Inputs /////////////////////////////
;-------------------------------------------
	
	;Number of variables to plot
	tf_current = keyword_set(current)
	if n_elements(no_refresh) gt 0 $
		then tf_refresh = ~keyword_set(no_refresh) $
		else tf_refresh = 1B
	
	;Plot layout
	if n_elements(layout) eq 0 then begin
		layout = [1, nVars]
	endif else if nVars gt layout[0]*layout[1] then begin
		message, 'LAYOUT is not large enough to fit all variables.'
	endif

;-------------------------------------------
; Create the Window ////////////////////////
;-------------------------------------------
	;Get the window
	if keyword_set(current) then begin
		IF Keyword_Set(buffer) THEN MrPrintF, 'LogWarn', 'Using current window. Ignoring BUFFER keyword.'
		win        = GetMrWindows(/CURRENT)
		tf_refresh = win -> GetRefresh()
		win       -> Refresh, /DISABLE
	endif else begin
		win = MrWindow( BUFFER  = buffer, $
		                LAYOUT  = layout, $
		                XSIZE   = xsize, $
		                YGAP    = 0.5, $
		                YSIZE   = ysize, $
		                REFRESH = 0 )
	endelse

;-------------------------------------------
; Determine Tick Marks /////////////////////
;-------------------------------------------
	;Determine where tick marks are placed
	t_extra = MrVar_TTIcks(LABEL=t_label)
	
	;Keep track of margins
	oxmargin = win.oxmargin

;-------------------------------------------
; Step Through Each Variable ///////////////
;-------------------------------------------
	
	;Plot each variable
	for i = 0, nVars - 1 do begin
		iCol = i mod layout[0]
		iRow = i  /  layout[0]

		;Get the variable
		oVar = nVars eq 1 ? allVars : allVars[i]

	;-------------------------------------------
	; Unknown Types ////////////////////////////
	;-------------------------------------------
		if oVar -> HasAttr('DEPEND_3') || oVar -> HasAttr('DEPEND_2') then begin
			MrPrintF, 'LogErr', 'Unknown graphic type for "' + oVar.name + '". Maximum of DEPEND_1 expected.'
			
	;-------------------------------------------
	; Image ////////////////////////////////////
	;-------------------------------------------
		endif else if oVar -> HasAttr('DEPEND_1') then begin
			gfx = MrVar_Image(oVar, /CURRENT)
			
	;-------------------------------------------
	; Plot /////////////////////////////////////
	;-------------------------------------------
		endif else if oVar -> HasAttr('DEPEND_0') then begin
			gfx = MrVar_Plot(oVar, /CURRENT)

	;-------------------------------------------
	; Non-Time-Series //////////////////////////
	;-------------------------------------------
		endif else begin
			MrPrintF, 'LogErr', 'Variable must have a DEPEND_0 attribute: "' + oVar.name + '".'
		endelse
	
	;-------------------------------------------
	; Prettify /////////////////////////////////
	;-------------------------------------------
		if iRow gt 0       then title       = ''     else void = temporary(title)
;		if iRow lt nVars-1 then xtickformat = '(a1)' else void = temporary(xtickformat)
;		if iRow lt nVars-1 then xtitle      = ''     else void = temporary(xtitle)
		gfx -> SetProperty, TITLE       = title, $
		                    XMINOR      = t_extra.xminor, $
		                    XRANGE      = t_extra.xrange, $
		                    XTICKFORMAT = (iRow lt nVars-1 ? '(a1)' : ''), $
		                    XTICKV      = t_extra.xtickv, $
		                    XTICKNAME   = (iRow lt nVars-1 ? '' : t_extra.xtickname), $
		                    XTICKS      = t_extra.xticks, $
		                    XTITLE      = (iRow lt nVars-1 ? '' : '!CTime (UT)')

	;-------------------------------------------
	; Create a Label ///////////////////////////
	;-------------------------------------------
		if iRow eq nVars-1 && iCol eq 0 then begin
			xyo1 = MrText( 0.0, 0.0, t_label, $
			               ALIGNMENT          = 1.0, $
			               /RELATIVE, $
			               TARGET             = gfx, $
			               VERTICAL_ALIGNMENT = 1.0 )
			oxmargin[0] >= 12
		endif
	endfor

;-------------------------------------------
; Finish ///////////////////////////////////
;-------------------------------------------
	;Put the variables into the common block
	vars = allVars

	;Time range
	;   - If the time interval is too short, MrVar_Plot will change units from
	;     in seconds since midnight to seconds since t0. Setting the global XRANGE
	;     here will not account for that.
;	trange  = MrVar_GetTRange('SSM')
;	win    -> SetGlobal, XRANGE=trange

	IF ~Array_Equal(win.oxmargin, oxmargin) $
		THEN win.oxmargin >= oxmargin
	
	;Return the plot
	if tf_refresh then win -> Refresh
	return, win
end