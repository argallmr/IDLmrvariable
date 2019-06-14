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
;   Import TPlot variables as MrVariables
;
;   NOTE:
;       Requires the SPEDAS distribution.
;       http://spedas.org/blog/
;
; :Params:
;       VARIABLES:      in, optional, type=string/strarr, default='*'
;                       Names of variables to be imported. If the names contain wildcards,
;                           all variable names that match the wildcard pattern will be
;                           imported.
;       FILENAME:       in, optional, type=string
;                       The name of a tplot save file from which data is loaded.
;
; :Keywords:
;       VARNAMES:       out, optional, string/strarr
;                       Names of the variables imported from SPEDAS.
;       VERBOSE:        in, optional, type=boolean, default=false
;                       If set, helpful text will be printed to the console.
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
;       2018-03-06  -   Written by Matthew Argall
;       2019-02-08  -   Added the VARNAMES keyword. Anticipate the type of time
;                           series object to create (e.g. scalar, vector, matrix). - MRA
;-
;*****************************************************************************************
PRO MrVar_ImportFromSPEDAS, variables, filename, $
VARNAMES=out_names, $
VERBOSE=verbose
	Compile_Opt idl2
	
	;Get all of the variables
	tf_verbose = Keyword_Set(verbose)
	IF N_Elements(variables) EQ 0 THEN variables = '*'
	
	;Restore TPlot save file
	IF N_Elements(filename) GT 0 THEN TPlot_Restore, FILENAME=filename, VERBOSE=tf_verbose
	
	;Find the variables to import
	vnames = tnames(variables, nvars)
	out_names = StrArr(nvars)

;-------------------------------------------
; Import Variables /////////////////////////
;-------------------------------------------
	FOR i = 0, nVars - 1 DO BEGIN

	;-------------------------------------------
	; Create Variables /////////////////////////
	;-------------------------------------------
		Get_Data, vnames[i], DATA=d, DLIMITS=dl
		IF tf_verbose THEN MrPrintF, 'LogText', 'Importing "' + vnames[i] + '".'
		
		;Create the time variable
		oDep0 = MrTimeVar( Time_String(d.x, /AUTOPREC), '%Y-%M-%d/%H:%m:%S%f' )
		
		;Guess the type of time series object
		dims = Size(d.y, /DIMENSIONS)
		ndims = N_Elements(dims)
		IF ndims EQ 1 THEN BEGIN
			class = 'MrScalarTS'
		ENDIF ELSE IF ndims EQ 2 && dims[1] EQ 3 THEN BEGIN
			class = 'MrVectorTS'
		ENDIF ELSE IF ndims EQ 3 && Array_Equal(dims[1:2], [3,3]) THEN BEGIN
			class = 'MrMatrixTS'
		ENDIF ELSE BEGIN
			class = 'MrTimeSeries'
		ENDELSE
		
		;Create the MrTimeSeries variable
		oVar = Obj_New( class, oDep0, d.y, $
		                /CACHE, $
		                NAME = vnames[i] )
		out_names[i] = oVar.name
		
		;Add dependency
		tf_dep1 = MrStruct_HasTag(d, 'V')
		IF tf_dep1 THEN BEGIN
			oDep1 = MrVariable(d.v)
			oVar['DEPEND_1'] = oDep1
		ENDIF
		
		;Prepare plotting paramters
		xtitle = StrArr(2)
		ytitle = StrArr(2)
		ztitle = StrArr(2)

	;-------------------------------------------
	; Add Attributes ///////////////////////////
	;-------------------------------------------
		IF Size(dl, /TNAME) EQ 'STRUCT' THEN BEGIN
			tags  = tag_names(dl)
			ntags = n_tags(dl)
			oY    = tf_dep1 ? oDep1 : oVar
			FOR j = 0, ntags - 1 DO BEGIN
			
				;Data Attributes
				IF tags[j] EQ 'DATA_ATT' THEN BEGIN
					dtags = tag_names(dl.data_att)
					ndtags = n_tags(dl.data_att)
					FOR k = 0, ndtags - 1 DO BEGIN
						CASE dtags[k] OF
							'COORD_SYS': oVar['COORDINATES'] = dl.data_att.coord_sys
							ELSE: IF tf_verbose THEN MrPrintF, 'LogText', '  Data Attr "' + dtags[k] + '" ignored.'
						ENDCASE
					ENDFOR
				
				;Plotting Attributes
				ENDIF ELSE BEGIN
					CASE tags[j] OF
						'XSTYLE':    oTime['STYLE']     = dl.xstyle
						'XSUBTITLE': xtitle[1]          = dl.xsubtitle
						'XTITLE':    xtitle[0]          = dl.xtitle
						'YRANGE':    oY['AXIS_RANGE']   = dl.yrange
						'YSTYLE':    oY['STYLE']        = dl.ystyle
						'YSUBTITLE': ytitle[1]          = dl.ysubtitle
						'YTITLE':    ytitle[0]          = dl.ytitle
						'ZRANGE':    oVar['AXIS_RANGE'] = dl.zrange
						'ZSTYLE':    oVar['STYLE']      = dl.zstyle
						'ZSUBTITLE': ztitle[1]          = dl.zsubtitle
						'ZTITLE':    ztitle[0]          = dl.ztitle
						'LOG':       oVar['LOG']        = dl.log
						ELSE: IF tf_verbose THEN MrPrintF, 'LogText', '  Attribute "' + tags[j] + '" ignored.'
					ENDCASE
				ENDELSE
			ENDFOR
			
			;Titles
			IF ~Array_Equal(xtitle, '') $
				THEN oTime['TITLE'] = xtitle[1] EQ '' ? xtitle[0] : StrJoin(xtitle, '!C')
			IF ~Array_Equal(ytitle, '') $
				THEN oY['TITLE']    = ytitle[1] EQ '' ? ytitle[0] : StrJoin(ytitle, '!C')
			IF ~Array_Equal(ztitle, '') $
				THEN oVar['TITLE']  = ztitle[1] EQ '' ? ztitle[0] : StrJoin(ztitle, '!C')
		ENDIF
	ENDFOR
END