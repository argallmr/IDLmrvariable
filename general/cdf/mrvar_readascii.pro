; docformat = 'rst'
;
; NAME:
;       MrVar_ReadAscii
;
;*****************************************************************************************
;   Copyright (c) 2016, University of New Hampshire                                      ;
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
;       * Neither the name of the University of New Hampshire nor the names of its       ;
;         contributors may be used to endorse or promote products derived from this      ;
;         software without specific prior written permission.                            ;
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
;   Read ascii data into the variable cache. If a field in `COLUMN_NAMES`is "Epoch",
;   this variable will be recognized as the time column and time-series data will be
;   created.
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
;       2016-11-19  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Helper function for MrVar_ReadAscii. Retrieve or create a variable.
;
;   1. Check if the variable has been read
;     a. Get the name from the cache
;     b. Look for pre-existing variables with the same name
;         i. Clobber DEPEND_0 attributes
;        ii. Create new variable -- clobber pre-existing var if necessary
;       iii. Keep a record of the variable so we know it has been created
;   2. Return variable
;     
;
;   NOTES:
;       Uses the MrVar_ReadCDF_common common block
;
; :Params:
;       VARNAME:            in, required, type=string
;                           Name of the CDF variable to retrieve or create.
;       VARTYPE:            in, optional, type=string, type='MrVariable'
;                           Type of varible to be created. Should be a subclass
;                               of MrVariable: {'MrScalarTS' | 'MrVectorTS' | 'MrMatrixTS'}
;
; :Keywords:
;       NO_CLOBBER:         in, required, type=boolean, default=0
;                           If set, variables will be renamed if one with the
;                               same name already exists in the cache.
;-
function MrVar_ReadAscii_GetVar, vname, vtype, $
NO_CLOBBER=no_clobber
	compile_opt idl2
	on_error, 2
	
	;Defaults
	;   - Clobber existing variables
	tf_clobber = ~keyword_set(no_clobber)
	
;-----------------------------------------------------
; Clobber Duplicate Name \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Check if there is a variable by the same name already in the cache
	old_var = MrVar_Get(vname, COUNT=count)

	;If there is, and we clobber, then clobber its DEPEND_# attributes
	if count gt 0 && tf_clobber then begin
		;Look for dependent variable
		depVal = old_var -> GetAttrValue('DEPEND_0', /NULL)
		
		;If DEPEND_0 was an object, it may have been made
		;invalid by clobbering a different variable.
		if size(depVal, /TNAME) eq 'OBJREF' then begin
			if obj_valid(depVal) $
				then depName = depVal.name $
				else old_var -> RemoveAttr, 'DEPEND_0'
		endif else begin
			depName = depVal
		endelse

		;Remove dependent variable
		;   - If the variable has already been read, then it has
		;     already been sufficiently clobbered.
		if depName ne !Null && ~MrIsMember(cache_vnames, depName) then begin
			MrVar_Delete, depName
			old_var -> RemoveAttr, depN
		endif
	endif

;-----------------------------------------------------
; Create & Cache Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;We do not want to create new variables every time the same
	;file is read, so default to clobbering any pre-existing
	;variable by the same name.
	oVar = obj_new( vtype, $
	                /CACHE, $
	                NAME       = vname, $
	                NO_CLOBBER = ~tf_clobber )

	;Return the variable
	return, oVar
end


;+
;   Read ascii data into the variable cache. If a field in `COLUMN_NAMES`is "Epoch",
;   this variable will be recognized as the time column and time-series data will be
;   created.
;
; :Params:
;    FILENAME:          in, optional, type=string, default=''
;                       Filename of file to read. If no file is given or the empty string,
;                           a dialog box will appear asking to pick a file. Furthermore,
;                           if `TEMPLATE`, `COLUMN_NAMES`, `COLUMN_TYPES` and `GROUP` are
;                           all undefined, Ascii_Template() will be called to assist in
;                           reading the file.
;
; :Keywords:
;    COLUMN_NAMES:      in, optional, type=strarr
;                       Names for the columns in the data; if there are groups specified,
;                           the column names should be repeated for each column in the
;                           group, so that the number of column names is always equal to
;                           the number of columns. The default is to name the columns
;                           "Field#", where "#" represents the `GROUP` it belongs to.
;    COLUMN_TYPES:      in, optional, type=strarr/lonarr, default=4 (float)
;                       SIZE type codes for the columns in the data; if there are groups
;                           specified, the column types should be repeated for each column
;                           in the group, so that the number of column types is always
;                           equal to the number of columns
;    COLUMN_LOCATIONS:  in, optional, type=lonarr, default=lindgen(nColumsn)
;                       Location of the start of each column within a row of data. A line
;                           begins at location zero.
;    COMMENT_SYMBOL:    in, optional, type=string, default=''
;                       Specifies a comment character for the lines in the file
;    COUNT:             out, optional, type=long
;                       Set to a named variable to get the number of records read
;    DATA_START:        in, optional, TYPE=long, DEFAULT=0L
;                       Number of lines to skip at the beginning of the file. If data
;                           starts on line 10, then `DATA_START` should equal 9. See the
;                           `NHEADER` keyword.
;    DELIMITER:         in, optional, type=string, default=' '
;                       Delimiter between columns of data.
;    DIALOG_PARENT:     in, optional, type=widget_id
;                       The widget ID of that will serve as the parent to the dailog
;                           windows that may appears.
;    GROUPS:            in, optional, type=lonarr
;                       Indices of groups for each column, i.e.::
;
;                               [0, 0, 0, 0, 0, 0, 0]
;
;                           indicates all seven columns are in a single group, where::
;
;                               [0, 1, 2, 3, 4, 5, 6]
;
;                           would put each column in a new group. If columns are grouped
;                           together, they will be put into a 2D array, having the same
;                           number of columns as there are members in the group. The
;                           default is to put each column in its own group.
;    HEADER:            out, optional, type=strarr
;                       Set to a named variable to get the header information skipped by
;                           `DATA_START`
;    MISSING_VALUE:     in, optional, type=scalar, default=!values.f_nan
;                       Value to use for missing items in the data
;    NUM_RECORDS:       in, optional, type=long
;                       Number of records to read; default is to read all available records
;    NHEADER:           in, optional, type=integer, default=0
;                       Number of header lines in the file to skip. A less confusing way
;                           to specify the `DATA_START` keyword.
;    NFOOTER:           in, optional, type=integer, default=`DATA_START-1`
;                       Number of footer lines in the file. Will be removed from the data.
;    RECORD_START:      in, optional, type=long, DEFAULT=0
;                       Set to index of first record to read (after `DATA_START` is taken
;                           into account)
;    SELECT:            in, opitonal, type=boolean, default=0
;                       If set, Ascii_Template() will be called to initiate a GUI from
;                           which the following keywords are specified: `COLUMN_NAMES`,
;                           `COLUMN_TYPES`, `COLUMN_LOCATIONS`, `DATA_START`,
;                           `COMMENT_SYMBOL`, `DELIMITER`, `GROUPS`, and `MISSING_VALUE`.
;    TEMPLATE:          in, out, opitonal, type=structure
;                       A template returned by Ascii_Template(). If provided, all of the
;                           keywords related to Ascii_Template will be ignored. Upon exit,
;                           the template used to read the ascii file will be returned via
;                           this keyword.
;    T_TYPE:            in, opitonal, type=string, default='%Y-%M-%dT%H:%m:%S%f'
;                       A MrTokens pattern describing how epoch variables can be parsed.
;                           If a field in `COLUMN_NAMES` is "Epoch", all groups are
;                           considered time-dependent on this epoch column.
;    T_REF:             in, opitonal, type=string
;                       Set this equal to the time to which `T_TYPE` is referenced.
;    VARFORMAT:         in, optional, type=string/strarr
;                       Only the `COLUMN_NAMES` that match this string via StrMatch will
;                           be loaded into the variable cache.
;    VARNAMES:          out, optional, type=string/strarr
;                       Names of the variables loaded into the cache.
;    VARTYPE:           in, optional, type=string/strarr, default=''
;                       If EPOCH exists in `COLUMN_NAMES`, then VARTYPE describes the
;                           type of time series variable to create (e.g. MrScalarTS). The
;                           default is to create a MrScalarTS object for single column
;                           groups and MrTimeSeries objects for multiple column groups.
;                           There must be one entry per unique group.
;    VERBOSE:           in, optional, type=boolean
;                       Set to print runtime messages
;-
pro MrVar_ReadAscii, file, $
DIALOG_PARENT=dialog_parent, $
NHEADER=nHeader, $
NFOOTER=nFooter, $
SELECT=select, $
;MrVariable
T_TYPE=t_type, $
T_REF=t_ref, $
TRANGE=trange, $
VARNAMES=varnames, $
VARFORMAT=varformat, $
VARTYPE=vartype, $
;Ascii_Template
COLUMN_NAMES=column_names, $
COLUMN_TYPES=column_types, $
COLUMN_LOCATIONS=column_locations, $
COMMENT_SYMBOL=comment_symbol, $
DATA_START=data_start, $
DELIMITER=delimiter, $
GROUPS=groups, $
MISSING_VALUE=missingValue, $
TEMPLATE=template, $
;Read_Ascii
COUNT=count, $
HEADER=header, $
RECORD_START=recordStart, $
VERBOSE=verbose
	compile_opt idl2
	
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return
	endif
	
	;Defaults
	if n_elements(vartype)   eq 0 then vartype   = ''
	if n_elements(varformat) eq 0 then varformat = '*'
	if n_elements(trange)    eq 0 then trange    = MrVar_GetTRange()
	nVarFmt = n_elements(varformat)

;-----------------------------------------------------
; Read Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	data = MrFile_Read_nASCII( file, $
	                           DIALOG_PARENT    = dialog_parent, $
	                           NHEADER          = nHeader, $
	                           NFOOTER          = nFooter, $
	                           ;Ascii_Template
	                           COLUMN_NAMES     = column_names, $
	                           COLUMN_TYPES     = column_types, $
	                           COLUMN_LOCATIONS = column_locations, $
	                           COMMENT_SYMBOL   = comment_symbol, $
	                           DATA_START       = data_start, $
	                           DELIMITER        = delimiter, $
	                           GROUPS           = groups, $
	                           MISSING_VALUE    = missingValue, $
	                           ;Read_Ascii
	                           COUNT            = count, $
	                           HEADER           = header, $
	                           RECORD_START     = recordStart, $
	                           VERBOSE          = verbose )
	
	;Search for an EPOCH variable
	nVars    = n_tags(data)
	varnames = tag_names(data)
	iEpoch   = where( varnames eq 'EPOCH', nEpoch )

;-----------------------------------------------------
; Time Series Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	if nEpoch gt 0 then begin
		if size(data.(iEpoch), /N_DIMENSIONS) eq 2 then begin
			epoch = strjoin(data.epoch, ' ')
			data  = MrStruct_ReplaceValue(data, {EPOCH: epoch})
		endif
		
		;Extract the time
		oTime = MrTimeVar( data.(iEpoch), t_type, $
		                   NAME  = varnames[iEpoch], $
		                   /NO_CLOBBER, $
		                   T_REF = t_ref )
		varnames[iEpoch] = oTime.name
		
		;Create MrTimeSeries variables
		iMatch = intarr(nVars)
		nMatch = 0
		for i = 0, nVars - 1 do begin
			;Skip the time variable that we already read
			if i eq iEpoch then CONTINUE
			
			;Skip if the name does not match the search pattern
			tf_strmatch = 0B
			if nVarFmt eq 1 $
				then tf_strmatch = strmatch(varnames[i], varformat) $
				else for j = 0, nVarFmt-1 do tf_strmatch >= strmatch(varnames[i], varformat[j])
			if ~tf_strmatch then CONTINUE
			
			;Keep track of variables
			iMatch[nMatch]  = i
			nMatch         += 1
			
			;Size of data
			nDims = size(data.(i), /N_DIMENSIONS)
			dims  = size(data.(i), /DIMENSIONS)
			
			;Variable type
			case 1 of
				vartype[0] ne '': oVar = obj_new(vartype[i], oTime, data.(i), NAME=varnames[i], /CACHE )
				nDims      eq  1: oVar = MrScalarTS( oTime, data.(i), NAME=varnames[i], /CACHE )
				else:             oVar = MrTimeSeries( oTime, data.(i), NAME=varnames[i], /CACHE )
			endcase
		endfor
		
		;Destroy the data
		data   = !Null
		iMatch = iMatch[0:nMatch-1]

		;Trim time
		;   - The epoch variable is not in the cache
		if nMatch gt 0 && ~Array_Equal(trange, '') then begin
			varnames = varnames[iMatch]
			MrVar_TLimit, varnames, trange
		endif else begin
			varnames = ''
		endelse

;-----------------------------------------------------
; Other Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	endif else begin
		for i = 0, nVars - 1 do begin
			;Skip if the name does not match the search pattern
			tf_strmatch = 0B
			if nVarFmt eq 1 $
				then tf_strmatch = strmatch(varnames[i], varformat) $
				else for j = 0, nVarFmt-1 do tf_strmatch >= strmatch(varnames[i], varformat[j])
			if ~tf_strmatch then CONTINUE
			
			;Create and cache the variable
			oVar = MrVariable( data.(i), NAME=varnames[i], /CACHE )
		endfor
	endelse
end