; docformat = 'rst'
;
; NAME:
;       MrVar_ExportToSPEDAS
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
;   Export MrVariables to TPlot variables
;
;   NOTE:
;       Requires the SPEDAS distribution.
;       http://spedas.org/blog/
;
; :Params:
;       VARIABLES:      in, required, type=string/integer/objref
;                       Name, number, or objref of a MrTimeSeries object to be exported.
;                           Can be an array. Data is loaded into TPlot with the Store_Data
;                           procedure. Variables that do not subclass MrTimeSeries will
;                           silently be skipped. If not provided, or if it is the empty
;                           string, all variables in the cache are exported.
;       FILENAME:       in, optional, type=string
;                       The name of a tplot save file into which the data will be saved.
;                           The extension ".tplot" will be appended to FILENAME.
;
; :RETURN:
;       TNAMES:         out, required, type=string/strarr
;                       Names of the variables that were successfully loaded into tplot.
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
;       2017-01-24  -   Written by Matthew Argall
;       2018-03-06  -   Logic around file creation improved. Catch is cancelled. - MRA
;-
;*****************************************************************************************
FUNCTION MrVar_ExportToSPEDAS, variables, filename
	Compile_Opt idl2

	;Get all of the names
	IF N_Elements(variables) EQ 0 || ( MrIsA(variables, /SCALAR, 'STRING') && variables EQ '' ) $
		THEN MrVar_Names, variables
	
	;Number of names given
	nNames = N_Elements(variables)
	
	;Allocate memory and convert
	tnames = StrArr(nNames)
	count  = 0
	FOREACH var, variables DO BEGIN
		oVar = MrVar_Get(var)
		IF ~Obj_IsA(oVar, 'MrTimeSeries') THEN CONTINUE
		
		Catch, the_error
		IF the_error EQ 0 THEN BEGIN
			temp_name     = oVar -> ExportToSPEDAS()
			tnames[count] = Temporary(temp_name)
			count        += 1
		ENDIF ELSE BEGIN
			Catch, /CANCEL
			MrPrintF, 'LogWarn', 'Could not export "' + oVar.name + '".'
			MrPrintF, 'LogErr'
		ENDELSE
	ENDFOREACH

	;Trim results
	CASE count OF
		0:    tnames = ''
		1:    tnames = tnames[0]
		ELSE: tnames = tnames[0:count-1]
	ENDCASE

	;Save to file
	IF N_Elements(filename) GT 0 THEN BEGIN
		Catch, the_error
		IF the_error EQ 0 THEN BEGIN
			tplot_save, tnames, FILENAME=filename
		ENDIF ELSE BEGIN
			MrPrintF, 'LogText', 'Error creating TPlot save file.'
			MrPrintF, 'LogErr'
		ENDELSE
		Catch, /CANCEL
	ENDIF
	
	RETURN, tnames
END