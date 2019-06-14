; docformat = 'rst'
;
; NAME:
;       MrVar_ReadCDF
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
;   Read data and its dependencies from multiple CDF data files and cache them
;   as MrVariables.
;
;   MRVARIABLE CACHE NOTES:
;       If a variable being read is already in the cache, and NO_CLOBBER=0, then
;       then it AND its DEPEND_# variables will be clobbered. This assumes that
;       all variables using the DEPEND_# support data come from the data files
;       being read. If other variables have been linked to the same DEPEND_# data,
;       then they will have orphaned variable attributes. Such a case is handled
;       internal to the variable, but may suprise the user when attempting to
;       plot the variable.
;
;       Many CDF files use the variable name "Epoch" for all time-dependencies.
;       In the case of NO_CLOBBER=0 described above, new variables will generate
;       new, unique "DEPEND_0" names (by default, "Epoch" would be renamed to
;       "Epoch_1").
;
;   NO_CLOBBER NOTES:
;       - A CDF variable is read, creating a variable and a DEPEND_0 attribute.
;       - The variable is destroyed, but DEPEND_0 gets carried on to child variables
;       - The next time the original variable is read, its DEPEND_0 attribute creates
;         a NEW cached variable (i.e. Epoch_1 instead of Epoch).
;       - This is to prevent clobbering the child variable's DEPEND_0 attribute.
;
;     Take-away lesson
;       - If you plan to read data more than once, keep the original CDF variable in
;         the cache.
;       - The other option is to rename the DEPEND_0 variables after you read them.
;         This way, the uniq variable name will be over-written each time.
;
;
;   MISC NOTES:
;       If DEPEND_# has no record variance, the values are assumed to be
;       uniform across files and data from the last file is returned.
;
;       `FILES` are assumed to be in chronological order. If not, data
;       will also not be in chronological order.
;
; :Categories:
;       CDF Utilities
;
; :Common Blocks:
;   MrVar_ReadCDF_common::
;       CDF_VNAMES     - CDF variable names of all unique variables from all files
;       CACHE_VNAMES   - MrVariable names of all variables that have been read
;       CDF_VCOUNT     - Total number of variables that have been read
;       FILE_VNAMES    - CDF variable names of variables that have been read from the current file.
;       FILE_VCOUNT    - Number of variables read from the current file
;       SUPPORT_VNAMES - The variables in CACHE_VNAMES that were read as support data.
;       SUPPORT_VCOUNT - The number of variable names in SUPPORT_VNAMES.
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
;       2016/04/30  -   Written by Matthew Argall
;       2016/07/04  -   DEPEND_N was being clobbered for each variable. Now clobber
;                           only for first variable. - MRA
;       2016/07/06  -   VARFORMAT can be an array. - MRA
;       2016/07/13  -   Added the VARNAMES keyword. Implemented TRANGE keyword. - MRA
;       2016/07/20  -   Added MrVar_ReadCDF_VarAttr2GfxKwd. - MRA
;       2016/08/20  -   Prevent variables from being read twice - once as an attribute
;                           and once as a variable - when /SUPPORT_DATA is set. - MRA
;       2016/10/06  -   Added the SUFFIX keyword. - MRA
;       2016/10/06  -   Major rewrite to support MrTimeSeries & simplify program flow. - MRA
;       2016/01/04  -   Support variables that are read via _GetVarAttrs are removed from
;                           the variable cache unless SUPPORT_DATA or VARFORMAT are set.
;                           They still exist as variable attributes, though. - MRA
;-
;*****************************************************************************************
;+
;    Helper function for MrCDF_nRead. Clobber a variable.
;
; :Params:
;       VARNAME:            in, required, type=string/object
;                           Name of the variable to be clobbered.
;-
PRO MrVar_ReadCDF_Clobber, varname
	Compile_Opt idl2
	On_Error, 2

	Common MrVar_ReadCDF_common, cdf_vnames, cdf_vcount, file_vnames, file_vcount, $
	                             cache_vnames, support_vnames, support_vcount, vsuffix
	
	;Get the variable
	iRemove = Where(cdf_vnames EQ varname, nRemove, COMPLEMENT=iKeep, NCOMPLEMENT=nKeep)
	IF nRemove GT 1 THEN Message, 'Duplicate variable names found'
	
	;Remove the variables
	IF nRemove GT 0 THEN BEGIN
		;Remove and delete variable from cache
		MrVar_Delete, cache_vnames[iRemove]
	
		;Remove variable name from common block variables
		cache_vnames[iRemove] = ''
		cdf_vnames[iRemove]   = ''
		cdf_vcount           -= 1
		
		;Reorder
		iSpace = Where(cdf_vnames EQ '', nSpace, COMPLEMENT=iName, NCOMPLEMENT=nName)
		IF nSpace GT 0 && nName GT 0 THEN BEGIN
			cache_vnames = [cache_vnames[iName], cache_vnames[iSpace]]
			cdf_vnames   = [cdf_vnames[iName],   cdf_vnames[iSpace]]
		ENDIF
	ENDIF
	
	;Reorganize the variables from the current file as well.
	iRemove = Where(file_vnames EQ varname or file_vnames EQ '', nRemove, COMPLEMENT=iKeep, NCOMPLEMENT=nKeep)
	IF nRemove GT 0 THEN BEGIN
		file_vnames[iRemove]  = ''
		file_vcount          -= 1
		IF nKeep GT 0 THEN file_vnames = [file_vnames[iKeep], file_vnames[iRemove]]
	ENDIF
END


;+
;    Helper function for MrCDF_nRead. Read and append data.
;
; :Params:
;       CDFID:              in, required, type=string/long
;                           CDF identifier of the CDF file being read.
;       VARNAME:            in, required, type=string/object
;                           Name of the variable whose data will be read.
;
; :Keywords:
;       COUNT:              in, optional, type=intarr
;                           Vector containing the counts to be used when reading each
;                               `VALUE`. The default is to read each record, taking
;                               into account `INTERVAL` and `OFFSET`.
;       NO_CLOBBER:         in, optional, type=boolean, default=0
;                           If set, variables with names that already exist in the
;                               MrVariable cache will be renamed by appending '_#' to
;                               the variable name, where "#" represents a number.
;       VARINQ:             out, optional, type=structure
;                           Result of CDF_VarInq() on the variable for which data is
;                               returned.
;-
FUNCTION MrVar_ReadCDF_GetData, cdfID, varname, $
NO_CLOBBER = no_clobber, $
VARINQ = varinq
	Compile_Opt idl2
	On_Error, 2
	
	Common MrVar_ReadCDF_common, cdf_vnames, cdf_vcount, file_vnames, file_vcount, $
	                             cache_vnames, support_vnames, support_vcount, vsuffix
	
;-----------------------------------------------------
; Read Variable Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Number of records -- read all
	CDF_Control, cdfID, GET_VAR_INFO=var_info, VARIABLE=varname
	rec_count = var_info.maxrec + 1
	
	;
	; If REC_COUNT=0, CDF_VarGet will issue an error. However, IF we undefine
	; REC_COUNT, it will read a PadValue and issue a warning. Here, we choose
	; to read in zero records
	;
	
	IF rec_count GT 0 THEN BEGIN
		;Do not show annoying cdf_varget warnings
		!Quiet = 1

		;Get its data
		CDF_VarGet, cdfID, varname, data, $
		            REC_COUNT = rec_count, $
		            /STRING

		;Turn on normal warnings.
		!Quiet = 0
	ENDIF ELSE BEGIN
		MrPrintF, 'logwarn', 'Variable has zero records: "' + varname + '".'
		data = !Null
	ENDELSE

;-----------------------------------------------------
; Book-keeping \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Store variables from current file
	;   - Collection is necessary once per file when data has been read
	file_vnames[file_vcount]  = varname
	file_vcount              += 1
	
	RETURN, data
END


;+
;   Helper function for MrVar_ReadCDF. Retrieve or create a variable.
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
FUNCTION MrVar_ReadCDF_GetVar, vname, vtype, $
NO_CLOBBER=no_clobber, $
SUFFIX=suffix
	Compile_Opt idl2
	On_Error, 2
	
	Common MrVar_ReadCDF_common, cdf_vnames, cdf_vcount, file_vnames, file_vcount, $
	                             cache_vnames, support_vnames, support_vcount, vsuffix
	
	;Defaults
	;   - Clobber existing variables
	tf_clobber = ~Keyword_Set(no_clobber)
	
	;Has the variable been read yet?
	tf_has = Max(cdf_vnames EQ vname, imax)

;-----------------------------------------------------
; Get Existing Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF tf_has THEN BEGIN
		oVar = MrVar_Get(cache_vnames[imax])

;-----------------------------------------------------
; Create New Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	ENDIF ELSE BEGIN

	;-----------------------------------------------------
	; Clobber Duplicate Name \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Check if there is a variable by the same name already in the cache
		old_var = MrVar_Get(vname+vsuffix, COUNT=count)

		;If there is, and we clobber, then clobber its DEPEND_[0-3] and
		;DELTA_(MINUS|PLUS)_VAR attributes
		IF count GT 0 && tf_clobber THEN BEGIN
			attrClobber = [ 'DEPEND_' + String( [0,1,2,3], FORMAT='(i1)' ), $
			                'DELTA_' + ['MINUS', 'PLUS'] + '_VAR', $
			                'LABL_PTR_' + String( [1,2,3], FORMAT='(i1)' ) ]
		
			;Step through each possible DEPEND_N
			FOR i = 0, N_Elements(attrClobber)-1 DO BEGIN
				;Look for attributes
				attrVal  = old_var -> GetAttrValue(attrClobber[i], /NULL)

				;If ATTRVAL was an object, it may have been made
				;invalid by clobbering a different variable.
				IF Size(attrVal, /TNAME) EQ 'OBJREF' THEN BEGIN
					IF Obj_Valid(attrVal) $
						THEN attrName = attrVal.name $
						ELSE old_var -> RemoveAttr, attrClobber[i]
				ENDIF ELSE BEGIN
					attrName = attrVal
				ENDELSE

				;Remove dependent variable
				;   - If the variable has already been read, then it has
				;     already been sufficiently clobbered.
				IF attrName NE !Null && ~MrIsMember(cache_vnames, attrName) THEN BEGIN
					MrVar_Delete, attrName
					old_var -> RemoveAttr, attrClobber[i]
				ENDIF
			ENDFOR
		ENDIF

	;-----------------------------------------------------
	; Create & Cache Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;We do not want to create new variables every time the same
		;file is read, so default to clobbering any pre-existing
		;variable by the same name.
		oVar = Obj_New( vtype, $
		                /CACHE, $
		                NAME       = vname + vsuffix, $
		                NO_CLOBBER = ~tf_clobber )
		
		;Create an attribute
		oVar -> AddAttr, 'CDF_NAME', vname
		
		;Keep track of variable information
		;   - Collection is necessary only once per call (i.e. not for every file)
		;   - Unique variable names across all files being read.
		cdf_vnames[cdf_vcount]    = vname
		cache_vnames[cdf_vcount]  = oVar.name
		cdf_vcount               += 1
		IF cdf_vcount GE N_Elements(cdf_vnames) THEN BEGIN
			cdf_vnames = [cdf_vnames, StrArr(cdf_vcount)]
			cache_vnames = [cache_vnames, StrArr(cdf_vcount)]
		ENDIF
	ENDELSE

	;Return the variable
	RETURN, oVar
END


;+
;    Helper function for MrCDF_nRead. Read and append data.
;
; :Params:
;       CDFID:              in, required, type=string/long
;                           CDF identifier of the CDF file being read.
;       VARNAME:            in, required, type=string/object
;                           Name of the variable whose data will be read.
;-
PRO MrVar_ReadCDF_GetVarAttrs, cdfID, oVar
	Compile_Opt idl2
	On_Error, 2

	Common MrVar_ReadCDF_common, cdf_vnames, cdf_vcount, file_vnames, file_vcount, $
	                             cache_vnames, support_vnames, support_vcount, vsuffix

	;Number of attributes
	;   - Variable attributes should be saved after global attributes
	;   - I have seen some files where this is not the CASE
	;   - So, loop over all attributes
	cdfinq  = CDF_Inquire(cdfID)
	varname = oVar['CDF_NAME']
	
;-----------------------------------------------------
; Loop Over Attributes \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	FOR i = 0, cdfinq.natts - 1 do BEGIN
		;Check IF the variable has the attribute
		;   - Skip the attribute if it does not exist for this variable
		CDF_AttGet_Entry, cdfID, i, varname, cdf_type, attrValue, tf_exists, ATTRIBUTE_NAME=attrName
		IF tf_exists EQ 0 THEN CONTINUE

	;-----------------------------------------------------
	; DELTA_(PLUS|MINUS)_VAR or DEPEND_[0-3] or LABL_PTR_[1-3] 
	;-----------------------------------------------------
		IF StRegEx(attrName, '(DELTA_(PLUS|MINUS)_VAR|(DEPEND|LABL_PTR)_[0-3])', /BOOLEAN) THEN BEGIN
			
			;
			; Bug in MMS_FPI_DES-DIST files has DELTA_PLUS|MINUS_VAR refer
			; back to the parent variable.
			;
			IF varname EQ attrValue THEN CONTINUE
			
			;
			; The DELTA_PLUS_VAR, DELTA_MINUS_VAR, DEPEND_[0-3], and LABL_PTR_[0-3]
			; attributes are all pointers to other variables. Their attribute value
			; is the name of another variable within the CDF. We will deal with them
			; in slightly different manners:
			;
			; DEPEND_[0-3]
			; Follow the pointer to the variable, read and cache its data, then set
			; the DEPEND_[0-3] attribute value equal to the cached variable name.
			;
			; LABL_PTR_[1-3]
			; Follow the pointer to the variable, read its data, then use it as the
			; variable attribute value.
			;
			; DELTA_(MINUS|PLUS)_VAR
			; Follow the pointer to the variable, read and cache its data, then set
			; the DEPEND_[0-3] attribute value equal to the cached variable name.
			;

			;ATTRVALUE is the name of the variable that serves as DEPEND_[0-3]
			;   - Check if the variable has yet been read.
			iName = Where(file_vnames EQ attrValue, nName)
			IF nName EQ 0 THEN BEGIN
				;Read the data
				;   - Creates and caches the variable
				;   - Do not clobber other DEPEND_0 data. They often
				;     have the same name (e.g. Epoch)
				oAttr = MrVar_ReadCDF_ReadVar( cdfID, attrValue, NO_CLOBBER=(attrName eq 'DEPEND_0') )
				
				;Save the names of the support variables.
				IF ~MrIsMember(support_vnames, attrValue) THEN BEGIN
					support_vnames[support_vcount] = oAttr.name
					support_vcount += 1
					IF support_vcount GE N_Elements(support_vnames) $
						THEN support_vnames = [support_vnames, StrArr(support_vcount)]
				ENDIF
			ENDIF ELSE BEGIN
				oAttr = MrVar_ReadCDF_GetVar(attrValue)
			ENDELSE
			
			oVar[attrName] = oAttr
			
	;-----------------------------------------------------
	; Other Attributes \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		ENDIF ELSE BEGIN
			;Set the attribute-value pair
			oVar[attrName] = attrValue
		ENDELSE
	ENDFOR
END


;+
;    Determine the type of MrVariable object to create.
;
; :Private:
;
; :Params:
;       CDFID:              in, required, type=string/long
;                           CDF identifier of the CDF file being read.
;       VARNAME:            in, required, type=string/object
;                           Name of the variable whose data will be read.
;
; :Returns:
;       VARTYPE:            out, required, type=string
;                           The type of MrVariable object that best represents
;                               data of `VARNAME`.
;-
FUNCTION MrVar_ReadCDF_GetVarType, cdfID, varname
	compile_opt idl2
	on_error, 2
	
	;What type of variable is it?
	;   - For zero-dimensional CDFs, DIMVAR will have one element whose value is zero.
	;   - Variables are still read in as 1xN
	varinq = cdf_varinq(cdfID, varname)
	nDims  = n_elements(varinq.dimvar)
	
	;Is there a DEPEND_0 attribute?
	cdf_attget_entry, cdfID, 'DEPEND_0', varname, attr_type, dep0_vname, tf_exists
	IF tf_exists THEN BEGIN
		dep0_inq      = cdf_varinq(cdfID, dep0_vname)
		tf_timeseries = MrIsMember(['CDF_EPOCH', 'CDF_EPOCH16', 'CDF_TIME_TT2000'], dep0_inq.datatype)
	ENDIF ELSE BEGIN
		tf_timeseries = 0B
	ENDELSE

	;Scalar or time
	IF nDims EQ 1 && varinq.dimvar[0] EQ 0 THEN BEGIN
		;Is the variable an epoch datatype?
		IF MrIsMember(['CDF_EPOCH', 'CDF_EPOCH16', 'CDF_TIME_TT2000'], varinq.datatype) $
			THEN vartype = 'MrTimeVar' $
			ELSE vartype = tf_timeseries ? 'MrScalarTS' : 'MrVariable'
	
	;Vector
	ENDIF ELSE IF nDims EQ 1 && varinq.dim[0] EQ 3 THEN BEGIN
		vartype = tf_timeseries ? 'MrVectorTS' : 'MrVariable'
	
	;Matrix
	;   - TODO: Should ANY 2D variable be a matrix?
	;           Could check FOR a "TENSOR_ORDER" attribute.
	ENDIF ELSE IF nDims EQ 2 THEN BEGIN
		vartype = tf_timeseries ? 'MrMatrixTS' : 'MrVariable'
	
	;Other
	;   - TODO: Transpose time to leading dimension?
	ENDIF ELSE BEGIN
		vartype = tf_timeseries ? 'MrTimeSeries' : 'MrVariable'
	ENDELSE

	RETURN, vartype
END


;+
;   Helper function for MrVar_ReadCDF. Retrieve or create a variable.
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
FUNCTION MrVar_ReadCDF_ReadVar, cdfID, varname, $
NO_CLOBBER=no_clobber, $
SUFFIX=suffix
	compile_opt idl2
	on_error, 2
	
	Common MrVar_ReadCDF_common, cdf_vnames, cdf_vcount, file_vnames, file_vcount, $
	                             cache_vnames, support_vnames, support_vcount, vsuffix

;-----------------------------------------------------
; Read Data & Attributes \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Which MrVariable type?
	vartype = MrVar_ReadCDF_GetVarType(cdfID, varname)
	
	;Get a variable object
	theVar = MrVar_ReadCDF_GetVar( varname, vartype, $
	                               NO_CLOBBER = no_clobber, $
	                               SUFFIX     = suffix )
	
	;
	; Infinite Loop Notice:
	;   Consider variable "A". Its DEPEND_0 attribute points to the
	;   "Epoch" variable, which, in turn, has a time-dependent
	;   DELTA_PLUS_VAR attribute with value "Epoch_Plus". Because
	;   "Epoch_Plus" is time-dependent, its DEPEND_0 attribute
	;   points back to the "Epoch" variable. To prevent an infinite
	;   loop, we must
	;
	;     1) Mark "Epoch" as read before obtaining its attributes
	;
	;   This will allow other variables to pull "Epoch" out of the cache
	;   for use as a DEPEND_0 attribute. It short-circuits the
	;   infinite loop by having "Epoch_Plus" pull "Epoch" out of the
	;   cache instead of cyclically calling MrVar_ReadCDF_ReadVar.
	;
	;     2) Retrieve attributes of non-Epoch variables before marking as read
	;
	;   MrTimeSeries objects require a time array to be provided
	;   with the data array. For that, the DEPEND_0 attribute must
	;   be retrieved and read from the file.
	;
	;   To accomplish this, we
	;     A) Get the attributes of MrTimeVar variables only after
	;        its data has been read.
	;     B) Get the attributes of all other variables types before
	;        trying to create the MrTimeSeries objects.
	;
	
	;Read attributes
	IF varType NE 'MrTimeVar' $
		THEN MrVar_ReadCDF_GetVarAttrs, cdfID, theVar

	;Create IDL graphics keywords attributes
	MrVar_ReadCDF_VarAttr2GfxKwd, theVar
	
	;Read data
	data = MrVar_ReadCDF_GetData(cdfID, varname)

;-----------------------------------------------------
; Add Attributes \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Save the CDF datatype as well
	varinq = cdf_varinq(cdfID, varName)
	theVar -> SetAttrValue, 'CDF_TYPE', varinq.datatype, /CREATE

;-----------------------------------------------------
; Append Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	; Record variance?
	;   - Append to record-varying dimension (always after the last CDF dimension)
	;   - ::Append can handle undefined data.
	; No record variance
	;   - Overwrite data
	;   - Let variables with no data remain as-is
	;
	
	;Variance
	IF varinq.recvar EQ 'VARY' THEN BEGIN
		;Data is read with dimensions of
		;   - [DEPEND_1, DEPEND_2, DEPEND_3, NRECS=DEPEND_0]
		;Transpose so that time is first:
		;   - [NRECS=DEPEND_0, DEPEND_1, DEPEND_2, DEPEND_3]
		ndims = size(data, /N_DIMENSIONS)
		IF ndims GT 1 THEN data = transpose(data, [nDims-1, bindgen(nDims-1)] )
		
		;Get the time object
		IF obj_isa(theVar, 'MrTimeSeries') THEN BEGIN
			IF n_elements(theVar) EQ 0 $
				THEN oTime = MrVar_Get(theVar['DEPEND_0']) $
				ELSE oTime = theVar['TIMEVAR']
		ENDIF
		
		;Append data
		IF vartype EQ 'MrTimeVar' $
			THEN data = [theVar['DATA'], theVar -> toISO(data, varinq.datatype)] $
			ELSE data = [theVar['DATA'], data]
	
	;NoVary
	ENDIF ELSE BEGIN
		IF StRegEx(vartype, '(MrTimeSeries|Mr(Scalar|Vector|Matrix)TS)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
			;Note change
			MrPrintF, 'LogWarn', 'Changing ' + vartype + ' to MrVariable'
			MrPrintF, 'LogWarn', '    VARNAME       = ' + varName
			MrPrintF, 'LogWarn', '    DISPLAY_TYPE  = ' + (theVar -> HasAttr('DISPLAY_TYPE') ? theVar['DISPLAY_TYPE'] : '')
			MrPrintF, 'LogWarn', '    RECVAR        = ' + varinq.recvar
			
			;Swap the variable
			oTemp   = MrVariable( NAME=theVar.name )
			theVar -> CopyAttrTo, oTemp
			MrVar_Replace, theVar, oTemp
			theVar  = oTemp
			vartype = 'MrVariable'
			
			;Remove the DEPEND_0 attribute
			theVar -> RemoveAttr, 'DEPEND_0'
		ENDIF
	ENDELSE

	;Set the data
	CASE vartype of
		'MrVariable':   theVar -> SetData, data, /NO_COPY
		'MrTimeVar':    theVar -> SetData, data, /NO_COPY
		'MrTimeSeries': theVar -> SetData, oTime, data, DIMENSION=1, /NO_COPY
		'MrScalarTS':   theVar -> SetData, oTime, data, DIMENSION=1, /NO_COPY
		'MrVectorTS':   theVar -> SetData, oTime, data, DIMENSION=1, /NO_COPY
		'MrMatrixTS':   theVar -> SetData, oTime, data, DIMENSION=1, /NO_COPY
		ELSE: message, 'Unknown variable type: "' + vartype + '".'
	ENDCASE
	
	;Read attributes
	IF varType EQ 'MrTimeVar' $
		THEN MrVar_ReadCDF_GetVarAttrs, cdfID, theVar
	
	RETURN, theVar
END


;+
;    Helper function for MrCDF_nRead. Read and append data.
;
; :Params:
;       OVAR:               in, required, type=MrVariable objref
;                           Variable for which the graphics keywords are set.
;-
PRO MrVar_ReadCDF_VarAttr2GfxKwd, oVar
	Compile_Opt idl2
	On_Error, 2

	;PLOT_TITLE
	IF ~oVar -> HasAttr('PLOT_TITLE') THEN BEGIN
		CASE 1 of
			oVar -> HasAttr('FIELDNAM'): oVar['PLOT_TITLE'] = oVar['FIELDNAM']
			oVar -> HasAttr('CATDESC'):  oVar['PLOT_TITLE'] = oVar['CATDESC']
			ELSE: ;Do nothing
		ENDCASE
	ENDIF

	;TITLE (Axis)
	IF ~oVar -> HasAttr('TITLE') THEN BEGIN
		CASE 1 of
			oVar -> HasAttr('LABLAXIS'): title = oVar['LABLAXIS']
			oVar -> HasAttr('FIELDNAM'): title = oVar['FIELDNAM']
			ELSE: title = ''
		ENDCASE
		IF oVar -> HasAttr('UNITS') THEN title = (title EQ '') ? oVar['UNITS'] : title + '!C' + oVar['UNITS']
		oVar['TITLE'] = Temporary(title)
	ENDIF
	
	;LABEL (Legend)
	IF ~oVar -> HasAttr('LABEL') THEN BEGIN
		IF oVar -> HasAttr('LABL_PTR_1') THEN oVar['LABEL'] = (oVar['LABL_PTR_1'])['DATA']
	ENDIF
	
	;MIN_VALUE
	IF ~oVar -> HasAttr('MIN_VALUE') THEN BEGIN
		IF oVar -> HasAttr('VALIDMIN')   THEN oVar['MIN_VALUE'] = oVar['VALIDMIN']
	ENDIF
	
	;MAX_VALUE
	IF ~oVar -> HasAttr('MAX_VALUE') THEN BEGIN
		IF oVar -> HasAttr('VALIDMAX')   THEN oVar['MAX_VALUE'] = oVar['VALIDMAX']
	ENDIF

	;SCALETYP
	;   - 'log' or 'linear'
	IF ~oVar -> HasAttr('LOG') THEN BEGIN
		IF oVar -> HasAttr('SCALETYP') THEN BEGIN
			tf_log = oVar['SCALETYP'] EQ 'log' ? 1 : 0
			oVar['LOG'] = tf_log
		ENDIF
	ENDIF
	
	;AXIS_RANGE
	IF ~oVar -> HasAttr('AXIS_RANGE') THEN BEGIN
		IF oVar -> HasAttr('SCALEMIN') && oVar -> HasAttr('SCALEMIN') THEN BEGIN
			oVar['AXIS_RANGE'] = [oVar['SCALEMIN'], oVar['SCALEMAX']]

			;If LOG is set, make sure the min axis range is >= 0
;			IF oVar -> HasAttr('LOG') && oVar['LOG'] THEN BEGIN
;				IF oVar['SCALEMIN'] LE 0 THEN BEGIN
;					iGood = oVar -> Where(0, /GREATER, COUNT=nGood)
;					IF nGood GT 0 $
;						THEN oVar['AXIS_RANGE'] = [ Min(oVar[iGood]), oVar['SCALEMAX'] ] $
;						ELSE oVar -> RemoveAttr, 'AXIS_RANGE'
;				ENDIF
;			ENDIF
		ENDIF
	ENDIF
	
	;DIMENSION
	IF obj_isa(oVar, 'MrVectorTS') THEN BEGIN
		oVar['DIMENSION'] = 1
		oVar['COLOR']     = ['Blue', 'Forest Green', 'Red']
	ENDIF
END


;+
;   A more robust method for obtaining variable data, when compared to the GetVarData
;   method.
;
; :Params:
;       FILE:               in, required, type=string/long
;                           The CDF identifier retured by CDF_Open() of the name of a CDF
;                               file from which data is read.
;
; :Keywords:
;       NO_CLOBBER:         in, optional, type=boolean, default=0
;                           If set, variables with names that already exist in the
;                               MrVariable cache will be renamed by appending '_#' to
;                               the variable name, where "#" represents a number.
;       SUFFIX:             in, optional, type=string, default=''
;                           A suffix to be appended to each variable name.
;       SUPPORT_DATA:       in, optional, type=boolean, default=0
;                           If set, support data will be read as well. Regardless of
;                               the status of this keyword, variable data associated
;                               with DEPEND_# variables will be read. This keyword
;                               is ignored if `VARFORMAT` is set.
;       TRANGE:             in, optional, type=strarr(2), default=MrVar_GetTRange()
;                           The time range over which to read data, formatted as ISO-8601
;                               strings: 'YYYY-MM-DDThh:mm:ss'.
;       STRING:             in, optional, type=boolean, default=1
;                           If set, CDF_CHAR and CDF_UCHAR will be read as strings. This
;                               is the default. STRING is ignored for other datatypes.
;       VALIDATE:           in, optional, type=boolean, default=0
;                           If set, the CDF file will be validated when opened.
;       VARFORMAT:          in, optional, type=string/strarr, default='*'
;                           Variable names that match this search pattern will be read.
;                               Any wildcard accepted by IDL's StrMatch function can be
;                               used.
;       VARNAMES:           out, optional, type=string/strarr
;                           Names of the variables that were read and cached. Names are
;                               as they exist in the cache, not necessarily as they are
;                               in the CDF file.
;       VERBOSE:            out, optional, type=boolean, default=0
;                           If set, status messages will be printed to standard output.
;-
PRO MrVar_ReadCDF, files, $
NO_CLOBBER=no_clobber, $
SUFFIX=suffix, $
SUPPORT_DATA=support_data, $
TRANGE=trange, $
VALIDATE=validate, $
VARFORMAT=varformat, $
VARNAMES=varnames, $
VERBOSE=verbose
	Compile_Opt idl2
	
	Common MrVar_ReadCDF_common, cdf_vnames, cdf_vcount, file_vnames, file_vcount, $
	                             cache_vnames, support_vnames, support_vcount, vsuffix

	;Get the time range (Will throw error if not set.)
	IF N_Elements(trange) EQ 0 THEN BEGIN
		Catch, the_error
		IF the_error EQ 0 $
			THEN trange = MrVar_GetTRange() $
			ELSE trange = ['', '']
		Catch, /CANCEL
	ENDIF

	;catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
	
		;Check that warnings are turned back on
		!Quiet = quiet_in
		IF status EQ 0 THEN status = 1
	
		;Close the file
		IF N_Elements(cdfID) GT 0 && tf_open THEN CDF_Close, cdfID

		;Turn file validation back on
		IF MrCmpVersion('8.0') LE 0 THEN $
			IF validate EQ 0 THEN CDF_Set_Validate, /YES
		
		;Issue error
		IF ~Arg_Present(status) THEN MrPrintF, 'LogErr'
		RETURN
	ENDIF

	;Defaults
	quiet_in   = !Quiet
	status     = 0
	tstart     = trange[0]
	tend       = trange[1]
	validate   = Keyword_Set(validate)
	tf_support = Keyword_Set(support_data)
	tf_verbose = Keyword_Set(verbose)
	nvarfmt    = N_Elements(varformat)
	vsuffix    = N_Elements(suffix) EQ 0 ? '' : suffix
	
;-----------------------------------------------------
; CDF file names or IDs \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;CDF IDs given?
	IF IsA(files, 'LONG') THEN BEGIN
		tf_IDs = 1B
	ENDIF ELSE BEGIN
		tf_IDs = 0B
	
		;Validate the files?
		IF MrCmpVersion('8.0') LE 0 THEN BEGIN
			IF validate $
				THEN CDF_Set_Validate, /YES $
				ELSE CDF_Set_Validate, /NO
		ENDIF
	ENDELSE
	
	;Allocate memory
	;   - Resets common block variables
	;   - We have to keep track of three sets of variable names:
	;       1) All unique CDF variable names from all files
	;       2) All variable names, as they are in the Cache
	;       3) Variable names from the current file
	;       4) Support data that is pointed to by a data variable's attribute value
	;   - #2 is to avoid name conflicts with variable already in the Cache
	;   - #3 is to avoid reading DEPEND_[0-3] variables more than once
	;   - #4 is so that we can remove them from the cache (they are saved as attribute values)
	nalloc         = 100
	cdf_vcount     = 0
	cache_vnames   = StrArr(nalloc)   ;All cache names
	cdf_vnames     = StrArr(nalloc)   ;All CDF names
	support_vcount = 0
	support_vnames = StrArr(nalloc)   ;Support variables
	timevar        = BytArr(nalloc)

;-----------------------------------------------------
; Step Through Each File \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	FOR i = 0, N_Elements(files) - 1 do BEGIN
		;Open the file
		IF tf_IDs THEN BEGIN
			cdfID   = files[i]
			tf_open = 0B
			IF tf_verbose THEN MrPrintF, 'LogText', files[i], FORMAT='(%"Reading CDF file with ID %i.")'
		ENDIF ELSE BEGIN
			cdfID   = CDF_Open(files[i])
			tf_open = 1B
			IF tf_verbose THEN MrPrintF, 'LogText', files[i], FORMAT='(%"Reading file \"%s\".")'
		ENDELSE
		
		;Number of variables in the file
		cdf_info    = CDF_Inquire(cdfID)
		nVars       = cdf_info.nvars + cdf_info.nzvars
		file_vnames = StrArr(nVars)      ;Names in this CDF file
		file_vcount = 0                  ;Counter FOR this CDF

	;-----------------------------------------------------
	; Step Through Each Variable \\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		FOR j = 0, nVars-1 do BEGIN
			;Variable index
			IF j GE cdf_info.nvars THEN BEGIN
				ivar   = j - cdf_info.nvars
				isZVar = 1B
			ENDIF ELSE BEGIN
				ivar   = j
				isZVar = 0B
			ENDELSE

			;Get the variable name and other info
			;   - This is possible IF, e.g., the variable is a DEPEND_0
			;     variable attribute.
			varinq = CDF_VarInq(cdfID, iVar, ZVARIABLE=isZVar)
			
			;Skip IF it has already been read
			IF MrIsMember(file_vnames, varinq.name) $
				THEN CONTINUE

			;Skip IF it does not match VARFORMAT
			IF nvarfmt GE 1 THEN BEGIN
				tf_strmatch = 0B
				IF nvarfmt EQ 1 $
					THEN tf_strmatch = StrMatch(varinq.name, varformat) $
					ELSE FOR k = 0, nvarfmt-1 DO tf_strmatch >= StrMatch(varinq.name, varformat[k])
				IF ~tf_strmatch THEN CONTINUE
			
			;Skip support data
			ENDIF ELSE IF ~tf_support THEN BEGIN
				CDF_AttGet_Entry, cdfID, 'VAR_TYPE', varinq.name, cdf_type, var_type, tf_exists
				IF tf_exists && var_type NE 'data' THEN CONTINUE
			ENDIF
			
			;Read data and attributes
			!Null = MrVar_ReadCDF_ReadVar( cdfID, varinq.name, $
			                               NO_CLOBBER = no_clobber, $
			                               SUFFIX     = suffix )
		
		ENDFOR
		
		;Close the data file
		IF tf_open THEN BEGIN
			CDF_Close, cdfID
			tf_open = 0B
			cdfID   = 0LL
		ENDIF
	ENDFOR
	
	;Trim results
	IF cdf_vcount EQ 0 THEN BEGIN
		cache_vnames = ''
	ENDIF ELSE BEGIN
		cdf_vnames   = cdf_vnames[0:cdf_vcount-1]
		cache_vnames = cache_vnames[0:cdf_vcount-1]
	ENDELSE
	support_vnames = support_vnames[0:support_vcount-1]
	
;-----------------------------------------------------
; Finish Up \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Limit time range
	;   - This must come before support variables are removed
	;     from the cache so that if any support variables are
	;     time-dependent, their time range changes with their
	;     parent variable.
	IF ~Array_Equal(trange, '') && cdf_vcount GT 0 $
		THEN MrVar_TLimit, cache_vnames, trange
	
	;Remove support variables from cache
	;   - They will still exist as variable attributes.
	IF ~tf_support && support_vcount GT 0 THEN BEGIN
		;Keep support variables that were requested explicitly
		tf_keep = BytArr(support_vcount)
		FOR i = 0, nvarfmt-1 DO BEGIN
			tf_keep OR= StrMatch(support_vnames, varformat[i])
		ENDFOR
		
		;Remove support data
		iRemove = Where(~tf_keep, nRemove)
		IF nRemove GT 0 THEN BEGIN
			;From cache
			MrVar_Remove, support_vnames[iRemove]
		
			;From common block variables
			!Null = MrIsMember(support_vnames[iRemove], cache_vnames, COMPLEMENT=iKeep)
			cdf_vnames   = cdf_vnames[iKeep]
			cache_vnames = cache_vnames[iKeep]
		ENDIF
	ENDIF

	;Return variable names
	IF Arg_Present(varnames) THEN varnames = cache_vnames

	;Turn file validation back on
	IF MrCmpVersion('8.0') LE 0 THEN $
		IF validate EQ 0 THEN CDF_Set_Validate, /YES
END
