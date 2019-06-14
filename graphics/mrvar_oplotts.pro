; docformat = 'rst'
;
; NAME:
;       MrVar_OPlotTS
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
;   Overplot time-series MrVariable data
;
; :Params:
;       GFX:        in, optional, type=string/strarr
;                   Names of the graphics objects over which data is plotted.
;       VARNAMES:   in, optional, type=string/strarr
;                   Names of the MrVariable objects to be overplotted.
;
; :Keywords:
;       WIN:        in, optional, type=objref, default=GetMrWindows(/CURRENT)
;                   A MrGraphics window object reference containing `GFXNAMES`.
;
; :Returns:
;       GFXWIN:     out, required, type=objref
;                   A MrGraphics window object reference.
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
;       2016/08/25  -   Graphics names or objects may be given. Variable index, name, or
;                           objrefs may be given. - MRA
;       2017/10/20  -   Iterate the correct number of times if GFX is an array and VARS
;                           is scalar. - MRA
;-
function MrVar_OPlotTS, gfx, vars, $
OVERPLOT=xCrushx, $
WIN=win, $
_REF_EXTRA=extra
	compile_opt strictarr

	;Catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
;-------------------------------------------
; Graphics /////////////////////////////////
;-------------------------------------------
	
	;OBJREFS
	IF Size(gfx, /TNAME) EQ 'OBJREF' THEN BEGIN
		;Extract names from objects
		;   - Objects and ObjArrs are both 'OBJREF' but have different TypeNames
		;   - For single objects, N_Elements can be overloaded
		nGfx   = TypeName(gfx) EQ 'OBJREF' ? N_Elements(gfx) : 1
		theGfx = gfx
	
	;NAMES
	ENDIF ELSE IF Size(gfx, /TNAME) EQ 'STRING' THEN BEGIN
		;Get the graphics window
		IF N_Elements(win) EQ 0 THEN BEGIN
			gfxWin = GetMrWindows(/CURRENT)
		ENDIF ELSE BEGIN
			IF Obj_IsA(win, 'MrWindow') $
				THEN gfxWin = win $
				ELSE Message, 'WIN must be a MrWindow object.'
		ENDELSE
		
		;Get the object reference for the graphics
		count  = 0
		nGfx   = N_Elements(gfx)
		theGfx = ObjArr(nGfx)
		FOR i = 0, nGfx - 1 DO BEGIN
			temp = gfxWin -> FindByName(gfx[i], COUNT=nTemp)
			IF nTemp EQ 1 THEN BEGIN
				theGfx[count] = temp
				count        += 1
			ENDIF ELSE BEGIN
				MrPrintF, 'LogWarn', 'No graphic found for name: "' + gfx[i] + '".'
			ENDELSE
		ENDFOR
		
		nGfx = count
		IF nGfx EQ 1 THEN theGfx = theGfx[0]
		
	;INVALID
	ENDIF ELSE BEGIN
		Message, 'GFX must be the names or objrefs of MrGraphics objects.'
	ENDELSE

;-------------------------------------------
; Variables ////////////////////////////////
;-------------------------------------------
	theVars = MrVar_Get(vars, COUNT=nVars)
	
	;Make same length
	IF nGfx EQ 1 && nVars GT 1 THEN BEGIN
		theGfx = Replicate(theGfx, nVars)
	ENDIF ELSE IF nVars EQ 1 && nGfx GT 1 THEN BEGIN
		theVars = Replicate(theVars, nGfx)
	ENDIF ELSE IF nVars NE nGfx THEN BEGIN
		Message, 'Invalid number of graphics (' + StrTrim(nGfx, 2) + ') and variables (' + StrTrim(nVars, 2) + ') found.'
	ENDIF
	nIter = nGfx > nVars

;-------------------------------------------
; Step Through Variables ///////////////////
;-------------------------------------------
	nWins      = 0
	tf_refresh = BytArr(nIter)
	win        = ObjArr(nIter)
	FOR i = 0, nIter - 1 DO BEGIN
		oVar   = theVars[i]
		oGfx   = theGfx[i]

	;-------------------------------------------
	; Keep Window //////////////////////////////
	;-------------------------------------------
		;Disable refresh
		win[nWins]        = oGfx.window
		tf_refresh[nWins] = win[nWins] -> GetRefresh()
		IF ~tf_refresh[nWins] THEN win[nWins] -> Refresh, /DISABLE
		
		;Keep track of the windows
		;   - Duplicate windows will either be overwritten (next iteration) or thrown out (after loop)
		IF nWins EQ 0 || Array_Equal( Obj_Valid(win[0:nWins-1], /GET_HEAP_IDENTIFIER) EQ Obj_Valid(win[nWins], /GET_HEAP_IDENTIFIER), 0 ) $
			THEN nWins += 1

	;-------------------------------------------
	; Add Graphic //////////////////////////////
	;-------------------------------------------
		
		;UNKNOWN
		IF oVar -> HasAttr('DEPEND_3') || oVar -> HasAttr('DEPEND_2') THEN BEGIN
			MrPrintF, 'LogErr', 'Unknown graphic type for "' + oVar.name + '". Maximum of DEPEND_1 expected.'
		
		;IMAGES
		ENDIF ELSE IF oVar -> HasAttr('DEPEND_1') THEN BEGIN
			tempGfx = MrVar_Image(oVar, OVERPLOT=oGfx, _EXTRA=extra)
		
		;SCATTER PLOTS
		ENDIF ELSE IF oVar -> HasAttr('DEPEND_0') THEN BEGIN
			tempGfx = MrVar_Plot(oVar, OVERPLOT=oGfx, _EXTRA=extra)

		;INVALID
		ENDIF ELSE BEGIN
			MrPrintF, 'LogErr', 'Variable must have a DEPEND_0 attribute: "' + oVar.name + '".'
		ENDELSE
	ENDFOR
	
	;Turn refresh back on
	win = win[0:nWins-1]
	FOR i = 0, nWins - 1 DO IF tf_refresh[i] THEN win[i] -> Refresh
	IF nWins EQ 1 THEN win = win[0]

	;Return the plot
	RETURN, win
END