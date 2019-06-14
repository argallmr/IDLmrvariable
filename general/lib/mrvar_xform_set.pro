; docformat = 'rst'
;
; NAME:
;       MrVar_xForm_Set
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
;       T:              in, required, type=string/integer/object
;                       Transformation matrix from one system to another.
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
PRO MrVar_xForm_Set, T
	Compile_Opt idl2
	On_Error, 2
	
	Common MrVar_xForm_Comm, xT
	
	;Check dimensions
	dims  = Size(T, /DIMENSIONS)
	nDims = N_Elements(dims)
	IF nDims NE 2 || ~Array_Equal(dims, 3) $
		THEN Message, 'T must be a 3x3 transformation matrix.'
	
	;A regular transformation matrix
	IF MrIsA(T, 'OBJREF') THEN BEGIN
		oT = MrVar_Get(T)
		IF Obj_Class(oT) NE 'MRVARIABLE' $
			THEN Message, 'T must be a MrVariable object.'
		xT = oT['DATA']
	ENDIF ELSE BEGIN
		xT = T
	ENDELSE
	
	;Put a time dimension first
	xT = Reform(xT, 1, 3, 3)
END