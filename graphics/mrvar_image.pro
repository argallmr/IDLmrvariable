; docformat = 'rst'
;
; NAME:
;       MrVar_Image
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
;   Create an image of MrVariable data
;
;   Calling Sequence:
;       p = MrVar_Image(z)
;       p = MrVar_Image(z, x, y)
;
; :Params:
;       Z:          in, required, type=string/integer/objref
;                   A MrVariable object, name or index of the variable to be displayed.
;       X:          in, optional, type=string/integer/objref
;                   A MrVariable object, name or index of the variable that represents
;                       the independent variable along the abscissa axis. If not given,
;                       `Z` will be searched for a 'DEPEND_0' attribute.
;       Y:          in, optional, type=string/integer/objref
;                   A MrVariable object, name or index of the variable that represents
;                       the independent variable along the ordinate axis. If not given,
;                       `Z` will be searched for a 'DEPEND_1' attribute.
;
; :Keywords:
;       CURRENT:    in, optional, type=boolean, default=0
;                   If set, the plot will be added to the current MrGraphics window.
;       _REF_EXTRA: in, optional, type=any
;                   Any keyword accepted by MrImage::SetProperty
;
; :Returns:
;       IM:         out, required, type=object
;                   A MrImage object reference.
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
;       2016/07/06  -   Written by Matthew Argall
;       2017/07/26  -   A time interval does not need to be defined (MrVar_Get/SetTRange) - MRA
;       2018/02/05  -   If X is time, then tick labels are formatted with t_ssm_labels. - MRA
;-
function MrVar_Image, z, x, y, $
CURRENT=current, $
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
	oZ = IsA(z, 'MrVariable') ? z : MrVar_Get(z)
	
	;X and Y given
	if n_elements(x) gt 0 && n_elements(y) gt 0 then begin
		oX = IsA(x, 'MrVariable') ? x : MrVar_Get(x)
		oY = IsA(y, 'MrVariable') ? y : MrVar_Get(y)
	
	;Use Z to get DEPEND_0 and DEPEND_1
	endif else if n_elements(x) eq 0 && n_elements(y) eq 0 then begin
		if ~oZ -> HasAttr('DEPEND_0') then message, 'Z does not have a DEPEND_0 attribute.'
		if ~oZ -> HasAttr('DEPEND_1') then message, 'Z does not have a DEPEND_1 attribute.'
		
		;Time series variable
		if obj_isa(oZ, 'MrTimeSeries') $
			then oX = oZ['TIMEVAR'] $
			else oX = oZ['DEPEND_0']
		
		oY = oZ['DEPEND_1']
	endif else begin
		message, 'Incorrect number of defined input parameters.'
	endelse

;-------------------------------------------
; Setup Graphics ///////////////////////////
;-------------------------------------------
	;Get the window
	if n_elements(overplot) gt 0 then begin
		win        = overplot.window
		tf_refresh = win -> GetRefresh()
	endif else if keyword_set(current) then begin
		win        = GetMrWindows(/CURRENT)
		tf_refresh = win -> GetRefresh()
	endif else begin
		win        = MrWindow(REFRESH=0)
		tf_refresh = 1
	endelse
	win -> Refresh, /DISABLE
	
	;Make room for the colorbar
	oxmargin    = win.oxmargin
	oxmargin[1] >= 15

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
; Delta Plus and Minus /////////////////////
;-------------------------------------------
	;X DELTA_PLUS
	if oX -> HasAttr('DELTA_PLUS_VAR') then begin
		oDX = MrVar_Get( oX['DELTA_PLUS_VAR'] )
		dxp = oDX['DATA']
	endif else if oX -> HasAttr('DELTA_PLUS') then begin
		dxp = oX['DELTA_PLUS']
	endif
	
	;X DELTA_MINUS
	if oX -> HasAttr('DELTA_MINUS_VAR') then begin
		oDX = MrVar_Get( oX['DELTA_MINUS_VAR'] )
		dxm = oDX['DATA']
	endif else if oX -> HasAttr('DELTA_MINUS') then begin
		dxm= oX['DELTA_MINUS']
	endif
	
	;Y DELTA_PLUS
	if oY -> HasAttr('DELTA_PLUS_VAR') then begin
		oDY = MrVar_Get( oY['DELTA_PLUS_VAR'] )
		dyp = oDY['DATA']
	endif else if oY -> HasAttr('DELTA_PLUS') then begin
		dyp = oY['DELTA_PLUS']
	endif
	
	;Y DELTA_MINUS
	if oY -> HasAttr('DELTA_MINUS_VAR') then begin
		oDY = MrVar_Get( oY['DELTA_MINUS_VAR'] )
		dym = oDY['DATA']
	endif else if oY -> HasAttr('DELTA_MINUS') then begin
		dym = oY['DELTA_MINUS']
	endif
	
	;Make sure DELTA_PLUS and DELTA_MINUS are defined
	if n_elements(dxm) gt 0 && n_elements(dxp) eq 0 then dxp = 0
	if n_elements(dxp) gt 0 && n_elements(dxm) eq 0 then dxm = 0
	if n_elements(dym) gt 0 && n_elements(dyp) eq 0 then dyp = 0
	if n_elements(dyp) gt 0 && n_elements(dym) eq 0 then dym = 0

	;Make sure both x and y are defined
	if n_elements(dxm) eq 0 && n_elements(dym) gt 0 then begin
		MrPrintF, 'LogWarn', 'No DELTA_(PLUS|MINUS) given for X.'
		dxm = 0
		dxp = x_data[1:*] - x_data
		dxp = [dxp, dxp[-1]]
	endif
	if n_elements(dxm) gt 0 && n_elements(dym) eq 0 then begin
		MrPrintF, 'LogWarn', 'No DELTA_(PLUS|MINUS) given for Y.'
		dym = 0
		if size(oY, /N_DIMENSIONS) eq 2 then begin
			dyp = oY['DATA',*,1:*] - oY['DATA']
			dyp = [ [dyp], [dyp[*,-1]] ]
		endif else begin
			dyp = oY['DATA',1:*] - oY['DATA']
			dyp = [dyp, dyp[-1]]
		endelse
	endif

;-------------------------------------------
; Define Pixel Areas ///////////////////////
;-------------------------------------------
	if n_elements(dxm) gt 0 || n_elements(dym) gt 0 then begin
		;Plot the data
		im = MrImage( oZ['DATA'], temporary(x_data), oY['DATA'], dxm, dym, dxp, dyp, $
		              /AXES, $
		              /CURRENT, $, $
		              NAME          = oZ.name, $
		              OVERPLOT      = overplot, $
		              
		              ;MrImage Keywords
		              HIDE        = oZ -> GetAttrValue('HIDE',        /NULL), $
		              KEEP_ASPECT = oZ -> GetAttrValue('KEEP_ASPECT', /NULL), $
		              LAYOUT      = oZ -> GetAttrValue('LAYOUT',      /NULL), $
		              POSITION    = oZ -> GetAttrValue('POSITION',    /NULL), $

		              ;POLAR Plot Options
		              POLAR          = oZ -> GetAttrValue('POLAR',          /NULL), $
		              POL_AXSTYLE    = oZ -> GetAttrValue('POL_AXSTYLE',    /NULL), $
		              POL_RCOLOR     = oZ -> GetAttrValue('POL_RCOLOR',     /NULL), $
		              POL_RLINESTYLE = oZ -> GetAttrValue('POL_RLINESTYLE', /NULL), $
		              POL_TCOLOR     = oZ -> GetAttrValue('POL_TCOLOR',     /NULL), $
		              POL_TLINESTYLE = oZ -> GetAttrValue('POL_TLINESTYLE', /NULL), $
		              POL_THICK      = oZ -> GetAttrValue('POL_THICK',      /NULL), $

		              ;IMAGE_PLOTS Keywords
		              AXISCOLOR     = oZ -> GetAttrValue('AXISCOLOR',     /NULL), $
		              BOTTOM        = oZ -> GetAttrValue('BOTTOM',        /NULL), $
		              CENTER        = oZ -> GetAttrValue('CENTER',        /NULL), $
		              DATA_POS      = oZ -> GetAttrValue('DATA_POS',      /NULL), $
		              LOG           = oZ -> GetAttrValue('LOG',           /NULL), $
		              MISSING_VALUE = oZ -> GetAttrValue('MISSING_VALUE', /NULL), $
		              MISSING_COLOR = oZ -> GetAttrValue('MISSING_COLOR', /NULL), $
		              MISSING_INDEX = oZ -> GetAttrValue('MISSING_INDEX', /NULL), $
		              NAN           = oZ -> GetAttrValue('NAN',           /NULL), $
		              PAINT         = oZ -> GetAttrValue('PAINT',         /NULL), $
		              RGB_TABLE     = oZ -> GetAttrValue('RGB_TABLE',     /NULL), $
		              RANGE         = oZ -> GetAttrValue('AXIS_RANGE',    /NULL), $
		              SCALE         = oZ -> GetAttrValue('SCALE',         /NULL), $
		              TITLE         = oZ -> GetAttrValue('PLOT_TITLE',    /NULL), $
		              TOP           = oZ -> GetAttrValue('TOP',           /NULL) )

;-------------------------------------------
; Display Image ////////////////////////////
;-------------------------------------------
	endif else begin
		;Plot the data
		im = MrImage( oZ['DATA'], temporary(x_data), oY['DATA'], $
		              /AXES, $
		              /CURRENT, $, $
		              NAME          = oZ.name, $
		              OVERPLOT      = overplot, $
		              
		              ;MrImage Keywords
		              HIDE        = oZ -> GetAttrValue('HIDE',        /NULL), $
		              KEEP_ASPECT = oZ -> GetAttrValue('KEEP_ASPECT', /NULL), $
		              LAYOUT      = oZ -> GetAttrValue('LAYOUT',      /NULL), $
		              POSITION    = oZ -> GetAttrValue('POSITION',    /NULL), $

		              ;POLAR Plot Options
		              POLAR          = oZ -> GetAttrValue('POLAR',          /NULL), $
		              POL_AXSTYLE    = oZ -> GetAttrValue('POL_AXSTYLE',    /NULL), $
		              POL_RCOLOR     = oZ -> GetAttrValue('POL_RCOLOR',     /NULL), $
		              POL_RLINESTYLE = oZ -> GetAttrValue('POL_RLINESTYLE', /NULL), $
		              POL_TCOLOR     = oZ -> GetAttrValue('POL_TCOLOR',     /NULL), $
		              POL_TLINESTYLE = oZ -> GetAttrValue('POL_TLINESTYLE', /NULL), $
		              POL_THICK      = oZ -> GetAttrValue('POL_THICK',      /NULL), $

		              ;IMAGE_PLOTS Keywords
		              AXISCOLOR     = oZ -> GetAttrValue('AXISCOLOR',     /NULL), $
		              BOTTOM        = oZ -> GetAttrValue('BOTTOM',        /NULL), $
		              CENTER        = oZ -> GetAttrValue('CENTER',        /NULL), $
		              DATA_POS      = oZ -> GetAttrValue('DATA_POS',      /NULL), $
		              LOG           = oZ -> GetAttrValue('LOG',           /NULL), $
		              MISSING_VALUE = oZ -> GetAttrValue('MISSING_VALUE', /NULL), $
		              MISSING_COLOR = oZ -> GetAttrValue('MISSING_COLOR', /NULL), $
		              MISSING_INDEX = oZ -> GetAttrValue('MISSING_INDEX', /NULL), $
		              NAN           = oZ -> GetAttrValue('NAN',           /NULL), $
		              PAINT         = oZ -> GetAttrValue('PAINT',         /NULL), $
		              RGB_TABLE     = oZ -> GetAttrValue('RGB_TABLE',     /NULL), $
		              RANGE         = oZ -> GetAttrValue('AXIS_RANGE',    /NULL), $
		              SCALE         = oZ -> GetAttrValue('SCALE',         /NULL), $
		              TITLE         = oZ -> GetAttrValue('PLOT_TITLE',    /NULL), $
		              TOP           = oZ -> GetAttrValue('TOP',           /NULL) )
	endelse

;-------------------------------------------
; Set Axis Properties //////////////////////
;-------------------------------------------
	
	;Set graphics properties
;	MrVar_SetGraphicsProperties, im, oZ, /XAXIS
	MrVar_SetAxisProperties, im, oX, /XAXIS
	MrVar_SetAxisProperties, im, oY, /YAXIS
	
	;Set user-given properties
	im -> SetProperty, XTICKFORMAT=xtickformat, XMINOR=xminor, XRANGE=xrange, XTICKS=xticks, XTICKV=xtickv, XTITLE=xtitle
;	im -> SetProperty, XTITLE=xtitle, XTICKFORMAT=xtickformat
	if n_elements(extra) gt 0 then im -> SetProperty, _STRICT_EXTRA=extra

;-------------------------------------------
; Create Colorbar //////////////////////////
;-------------------------------------------
	cb = MrColorbar( NAME        = 'CB: ' + oZ.name, $
	                 ORIENTATION = 1, $
	                 TARGET      = im )
	
	;Set properties
	MrVar_SetAxisProperties, cb, oZ, /CAXIS

;-------------------------------------------
; Finish ///////////////////////////////////
;-------------------------------------------
	;Adjust margins
	IF ~Array_Equal(win.oxmargin, oxmargin) $
		THEN win.oxmargin >= oxmargin
	
	;Return the plot
	if tf_refresh then im -> Refresh
	return, im
end