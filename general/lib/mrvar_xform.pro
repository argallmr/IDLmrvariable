; docformat = 'rst'
;
; NAME:
;       MrVar_xForm
;
;*****************************************************************************************
;   Copyright (c) 2018, Matthew Argall                                                   ;
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
;   Store a 3x3 time-independent transformation matrix.
;
;   Calling Sequence
;       oT = MrVar_xForm_Set(T)
;
; :Categories:
;   Coordinate Systems
;
; :Params:
;       VAR:            in, required, type=string/integer/object
;                       The name, number, or objref of the MrVariable object to be rotated.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the results will be added to the variable cache.
;       NAME:           in, optional, type=string, default='xForm(' + var.name + ')'
;                       Name to be given to the result.
;
; :Returns:
;       OOUT:           out, required, type=objref
;                       The results of rotating `VAR`. The object class is the same as
;                           that of `VAR`.
;
; :Author:
;       Matthew Argall::
;       University of New Hampshire
;       Morse Hall, Room 348
;       8 College Rd.
;       Durham, NH, 03824
;       matthew.argall@unh.edu
;
; :History:
;   Modification History::
;       2018/08/03  -   Written by Matthew Argall
;-
FUNCTION MrVar_xForm, var, $
NAME=name, $
CACHE=cache
	Compile_Opt idl2
	On_Error, 2
	
	Common MrVar_xForm_Comm, xT
	
	IF N_Elements(xT) EQ 0 THEN $
		Message, 'Tranformation matrix must be set with MrVar_xForm_Set.'
	
	;Check inputs
	oVar = MrVar_Get(var)
	IF N_Elements(name) EQ 0 THEN name = 'xForm(' + oVar.name + ')'
	
	;Create a matrix object
	nPts = oVar -> GetNPts()
	oT   = MrMatrixTS( oVar['TIMEVAR'], Rebin(xT, nPts, 3, 3) )
	
	;Transform data
	IF Obj_Class(oVar) EQ 'MRVECTORTS' THEN BEGIN
		oOut = oT -> Rotate_Vector(oVar)
	ENDIF ELSE IF Obj_Class(oVar) EQ 'MRMATRIXTS' THEN BEGIN
		oOut = oT -> Rotate_Matrix(oVAr)
	ENDIF ELSE BEGIN
		Message, 'Cannot rotate object of class ' + Obj_Class(oVar) + '.'
	ENDELSE
	Obj_Destroy, oT
	
	;Name and cache
	oOut -> SetName, name
	IF Keyword_Set(cache) THEN oOut -> Cache
	
	
	RETURN, oOut
END