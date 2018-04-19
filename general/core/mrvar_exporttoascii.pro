; docformat = 'rst'
;
; NAME:
;       MrVar_ExportToCDF
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
;   Export MrVariables to a CDF file
;
;   NOTES:
;       Requires the MrCDF distribution.
;           https://github.com/argallmr/IDLcdf
;
;       If there are two variables by the same name, only the first will be written
;       to the CDF file. This is important when it comes to DEPEND_# variables, which
;       are not incorporated into the MrVariable cache when read by MrVar_ReadCDF.
;       Instead, MrVar_ReadCDF sets them as attribute values without ensuring that they
;       have unique names.
;
; :Params:
;       FILENAME:       in, optional, type=string
;                       The name of the ASCII file into which the data will be saved.
;                           The extension ".txt" will be appended to FILENAME.
;       VARIABLES:      in, required, type=string/integer/objref
;                       Name, number, or objref of a MrTimeSeries object to be exported.
;                           Can be an array. Data is loaded into TPlot with the Store_Data
;                           procedure. Variables that do not subclass MrTimeSeries will
;                           silently be skipped. If not provided, or if it is the empty
;                           string, all variables in the cache are exported.
;       GLOBAL_ATTRS:   in, optional, type=struct/hash
;                       Global attribute name-value pairs. Structure tags or hash keys
;                           are used as attribute names and their values are written as
;                           global attribute values.
;
; :Keywords:
;       CLOBBER:        in, optional, type=boolean, default=0
;                       If set and `FILENAME` already exists, the existing file will be
;                           deleted.
;       COMPRESSION:    in, optional, type=string, default=''
;                       The compression type to be applied to the file.
;       GZIP_LEVEL:     in, optional, type=integer, default=5
;                       The level of gZip compression. This is applied to both `COMPRESSION`
;                           and `VARCOMPRESS` if either are 'GZIP'
;       VARCOMPRESS:    in, optional, type=string, default=''
;                       The compression type to be applied to each variable.
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
;       2017-03-13  -   Written by Matthew Argall
;-
PRO MrVar_ExportToASCII, filename, variables, $
CLOBBER=clobber
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(lun) GT 0 && (FStat(lun)).open THEN Close, lun
		RETURN
	ENDIF

	;Get all of the names
	IF N_Elements(variables) EQ 0 || ( MrIsA(variables, /SCALAR, 'STRING') && variables EQ '' ) THEN BEGIN
		MrVar_Names, variables
		IF variables[0] EQ '' THEN Message, 'No variables in cache.'
	ENDIF
	
	;Create or modify
	tf_clobber = Keyword_Set(clobber)

;-------------------------------------------
; Find DEPEND_0 Variables //////////////////
;-------------------------------------------
	nVars = N_Elements(variables)
	
	;Single variable
	IF nVars EQ 1 THEN BEGIN
		oVar = MrVar_Get(variables[0])
		oVar -> ExportToASCII, filename + '.txt'
	
	;Multiple Variables
	ENDIF ELSE BEGIN
		vars  = ObjArr(nVars)
		sz0   = LonArr(nVars)
		dep0  = ObjArr(nVars)
		fmt   = StrArr(nVars)
		iDep0 = IntArr(nVars) - 1S
		nDep0 = 0
		
		;Find common DEPEND_0 variables so as to not duplicate
		FOR i = 0, nVars - 1 DO BEGIN
			oVar   = MrVar_Get(variables[i])
			sz0[i] = (Size(oVar, /DIMENSIONS))[0]
			
			IF oVar -> HasAttr('DEPEND_0') THEN BEGIN
				tf_member = MrIsMember(dep0, oVar['DEPEND_0'], idx)
				IF tf_member THEN BEGIN
					iDep0[i]    = idx
				ENDIF ELSE BEGIN
					dep0[nDep0] = oVar['DEPEND_0']
					iDep0[i]    = nDep0
					nDep0      += 1
				ENDELSE
			ENDIF
			
			oVar -> ExportToASCII, FORMAT=format
			fmt[i]  = format
			vars[i] = oVar
		ENDFOR
		dep0 = dep0[0:nDep0-1]
	ENDELSE

;-------------------------------------------
; Write Each DEPEND_0 Variable /////////////
;-------------------------------------------
	FOR i = 0, nDep0 - 1 DO BEGIN
		iVars  = Where(iDep0 EQ i, nVars)
		format = fmt[iVars]
		format = '(' + StrMid(format[0], 0, StrPos(format[0], ',')) + ', ' + $
		         StrJoin(Reform(StrMid(Transpose(format), StrPos(transpose(format), ',')+2)), ', ') + ')'
		
		;Open the file for writing
		file = nDep0 EQ 1 ? filename : filename + '_' + String(i, FORMAT='(i0)')
		OpenW, lun, file+'.txt', /GET_LUN
		
		;Write each line individually		
		oDep0 = (vars[iVars[0]])['DEPEND_0']
		FOR j = 0, N_Elements(oDep0) - 1 DO BEGIN
			;Grab the data
			d0 = oDep0['DATA',j]
			SWITCH nVars OF
				12: d12 = (vars[iVars[11]])['DATA',j,*]
				11: d11 = (vars[iVars[10]])['DATA',j,*]
				10: d10 = (vars[iVars[19]])['DATA',j,*]
				 9:  d9 = (vars[iVars[18]])['DATA',j,*]
				 8:  d8 = (vars[iVars[7]])['DATA',j,*]
				 7:  d7 = (vars[iVars[6]])['DATA',j,*]
				 6:  d6 = (vars[iVars[5]])['DATA',j,*]
				 5:  d5 = (vars[iVars[4]])['DATA',j,*]
				 4:  d4 = (vars[iVars[3]])['DATA',j,*]
				 3:  d3 = (vars[iVars[2]])['DATA',j,*]
				 2:  d2 = (vars[iVars[1]])['DATA',j,*]
				 1:  d1 = (vars[iVars[0]])['DATA',j,*]
			ENDSWITCH
			
			;Write the data
			CASE nVars OF
				12: PrintF, lun, d0, d1, d2, d3, d4, d5, d6, d, d8, d9, d10, d11, d12, FORMAT=format
				11: PrintF, lun, d0, d1, d2, d3, d4, d5, d6, d, d8, d9, d10, d11, FORMAT=format
				10: PrintF, lun, d0, d1, d2, d3, d4, d5, d6, d, d8, d9, d10, FORMAT=format
				 9: PrintF, lun, d0, d1, d2, d3, d4, d5, d6, d, d8, d9, FORMAT=format
				 8: PrintF, lun, d0, d1, d2, d3, d4, d5, d6, d, d8, FORMAT=format
				 7: PrintF, lun, d0, d1, d2, d3, d4, d5, d6, d7, FORMAT=format
				 6: PrintF, lun, d0, d1, d2, d3, d4, d5, d6, FORMAT=format
				 5: PrintF, lun, d0, d1, d2, d3, d4, d5, FORMAT=format
				 4: PrintF, lun, d0, d1, d2, d3, d4, FORMAT=format
				 3: PrintF, lun, d0, d1, d2, d3, FORMAT=format
				 2: PrintF, lun, d0, d1, d2, FORMAT=format
				 1: PrintF, lun, d0, d1, FORMAT=format
			ENDCASE
		ENDFOR
		
		;Close the file
		Close, lun
	ENDFOR
END