; docformat = 'rst'
;
; NAME:
;   MrTimeSeries__Define
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
; PURPOSE
;+
;   The purpose of this class is subclass the IDL_Object such that all _Overload*
;   operations behave as they would for a normal array.
;
;   _OVERLOAD:
;       OPERATION:
;       An operator (such as AND, EQ, ^, *, etc.) calls the overload
;       method for the object on the left side first. So, for two MrTimeSeries
;       objects, A and B::
;
;               print, A [operator] B
;
;       will call A's _Overload method first.
;
;       INTERPOLATION:
;       If the operator acts between two MrTimeSeries objects and the size of
;       their time-dependent dimensions is different, the larger will be
;       interpolated via IDL's Interpol(/SPLINE) to the smaller.
;
;       ARRAY TRUNCATION:
;       If the operator acts between a MrTimeSeries object and either a MrVariable
;       object or a numerical expression, IDL's array truncation ruls apply. The
;       results will be truncated to have the same number of elements as the shorter
;       array. If an array and a scalar are compared, each elements of the array is
;       compared against the scalar value.
;
;       RESULT:
;       The result of such an operation will be a subclass of MrVariable::
;
;           MrTimeSeries [op] MrTimeSeries = MrTimeSeries
;           MrTimeSeries [op] MrVariable   = MrVariable
;           MrVariable   [op] MrTimeSeries = MrVariable
;           MrTimeSeries [op] Expression   = MrTimeSeries or MrVariable*
;           Expression   [op] MrTimeSeries = MrTimeSeries or MrVariable*
;
;                                        * MrTimeSeries if expression is
;                                          scalar or the size of its first
;                                          dimension matches that of the
;                                          implicit array. Otherwise,
;                                          MrVariable.
;
;   CACHING
;       The ::Cache method will cache the MrTimeSeries object. All variables in the
;       cache are forced to have unique names. When the object is destroyed, it
;       will be removed automatically from the cache.
;
; :Categories:
;   MrVariable, MrTimeSeries
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
;       2016/10/24  -   Written by Matthew Argall
;       2017/03/31  -   Testing revealed VAR->GetData() is faster than VAR['DATA'],
;                           so the change was made in all ::_Overload methods. Reduced
;                           number of IF ELSE cases in the ::_Overload methods.
;                           _OverloadBracketsRightSide returns and object with DEPEND_# and
;                           DELTA_(PLUS|MINUS)_VAR attributes are propertly reduced. - MRA
;       2017/05/31  -   Added the ::Digital_Filter method. - MRA
;       2018/01/26  -   Added NO_CLOBBER keyword to ::Copy. Default changed to clobber. - MRA
;-
;*****************************************************************************************
;+
;   The initialization method.
;
; :Params:
;       TIME:           in, required, type=NxM array
;                       Name or reference of a MrTimeVar object, or an array
;                           of time stamps. If a name is provided, the assiciated
;                           variable must exist in the variable cache.
;       DATA:           in, required, type=NxM array
;                       Name or reference of a MrVariable object, or the dependent
;                           variable data. If a name is given, the associated variable
;                           must exist in the variable cache. 
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, both `TIME` and `DATA` are added to the variable cache.
;       DIMENSION:      in, optional, type=integer
;                       The time-dependent, 1-based dimension of `DATA`. If not provided,
;                           the dimension of `DATA` that is equal in size to `TIME` is
;                           chose as the default.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;       NO_CLOBBER:     in, optional, type=boolean, default=0
;                       If set, do not clobber variables of the same name. Instead,
;                           rename this variable by appending "_#" to `NAME`, where
;                           "#" represents a unique number. Ignored unless `CACHE` is set.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       If set `DATA` will be copied directly into the object
;                           and will be left undefined (a MrTimeSeries object will not
;                           be destroyed, but its array will be empty).
;       T_NAME:         in, optional, type=integer
;                       Name to be given to the MrTimeVar object. Ignored unless `TIME`
;                           is an array of time stamps.
;       T_REF:              in, optional, type=string, default=''
;                           A reference time from which `TIME` is measured. Applicable
;                               only to certain values of `TYPE`. Must be formatted as
;                               an ISO-8601 string.
;       T_TYPE:         in, optional, type=integer
;                       If `TIME` is an array of time stamps, use this keyword to indicate
;                           the format or time-basis. See MrTimeVar for more details.
;-
function MrTimeSeries::INIT, time, data, $
CACHE=cache, $
DIMENSION=dimension, $
NAME=name, $
NO_CLOBBER=no_clobber, $
NO_COPY=no_copy, $
T_NAME=t_name, $
T_REF=t_ref, $
T_TYPE=t_type
	compile_opt idl2

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return, 0
	endif
	
	;Defaults
	if n_elements(name) eq 0 then name = 'MrTimeSeries'
	
	;Initialize superclass
	success = self -> MrVariable::Init()
	if ~success then message, 'Unable to initialize superclass.'
	
	;Set data
	if n_elements(time) gt 0 then begin
		self -> SetData, time, data, $
		                 DIMENSION = dimension, $
		                 NO_COPY   = no_copy, $
		                 T_NAME    = t_name, $
		                 T_REF     = t_ref, $
		                 T_TYPE    = t_type
	endif
	
	;Set name
	if n_elements(name) gt 0 then self -> SetName, name
	if keyword_set(cache)    then self -> Cache, NO_CLOBBER=no_clobber

	return, 1
end


;+
;   Clean up after the object is destroyed
;-
pro MrTimeSeries::CLEANUP
	compile_opt idl2
	on_error, 2
	
	;
	; Do not destroy the time object.
	;   - It may be used by another MrVariable object
	;   - Automatic garbage collection will destroy it when not used
	;
	
	;Superclass
	self -> MrVariable::Cleanup
end


;+
;   Perform a bit-wise AND comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadAnd, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'AND(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data AND right -> GetData() $
			ELSE result = left -> GetData() AND *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'AND(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) AND right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'AND(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left AND (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   Multiply two expressions together.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadAsterisk, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT MrVariable data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)
	
;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'Multiply(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) * (right -> GetData()) $
			else result = (right -> GetData()) * (*self.data)
		
		;Output object class
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'
	
;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Multiply(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) * right
			
		;RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Multiply(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left * (*self.data)
		endelse
	endelse
	
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
		
	;Create the object
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Allow square-bracket array indexing from the right side of an operator.
;
;   Calling Sequence
;       oTSvar    = oTSvar[0]
;       t0        = oTSvar[[0]]
;       data      = oTSvar['DATA']
;       pData     = oTSvar['POINTER']
;       pData     = oTSvar['PTR']
;       isoTime   = oTSvar['TIME']
;       time      = oTSvar['TIME', T_TYPE]
;       time      = oTSvar['TIME', i2]
;       time      = oTSvar['TIME', i2, <StrMid>]
;       time      = oTSvar['TIME', i2, T_TYPE]
;       oTime     = oTSvar['TIMEVAR']
;       oTime     = oTSvar['TIMEVAR', i1]
;       oTime     = oTSvar['TIMEVAR', i1, <strmid>]
;       oTime     = oTSvar['TIMEVAR', i1, T_TYPE]
;       attrValue = oTSvar[AttrName]
;       data      = oTSvar[i1 [, i2 [, i3 [, i4 [, i5 [, i6 [, i7 [, i8]]]]]]]]
;
; :Params:
;       ISRANGE:            in, required, type=intarr
;                           A vector that has one element for each Subscript argument
;                               supplied by the user; each element contains a zero if the
;                               corresponding input argument was a scalar index value or
;                               array of indices, or a one if the corresponding input
;                               argument was a subscript range.
;       I1:                 in, required, type=integer/intarr(3)
;                           Index subscripts. Either a scalar, an index array, or a 
;                               subscript range in the form [start, stop, step_size]
;       I2:                 in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       I3:                 in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       I4:                 in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       I5:                 in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       I6:                 in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       I7:                 in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       I8:                 in, optional, type=integer/intarr(3)
;                           Index subscripts.
;
; :Returns:
;       RESULT:             in, required, type=numeric array
;                           The subarray accessed by the input parameters.
;-
FUNCTION MrTimeSeries::_OverloadBracketsRightSide, isRange, i1, i2, i3, i4, i5, i6, i7, i8
	Compile_Opt idl2
	On_Error, 2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		Return, !Null
	endif
	
	;Number of subscripts given
	nSubs = N_Elements(isRange)
	
	;String operations
	IF IsA(i1, /SCALAR, 'STRING') THEN BEGIN
		CASE i1 OF
			'DATA': BEGIN
				IF nSubs EQ 1 $
					THEN RETURN, *self.data $
					ELSE RETURN, self -> _OverloadBracketsRightSide_Data(isRange[1:*], i2, i3, i4, i5, i6, i7, i8)
			ENDCASE
			'POINTER': RETURN, self.data
			'PTR':     RETURN, self.data
			'TIMEVAR': BEGIN
				IF nSubs EQ 1 $
					THEN RETURN, self.oTime $
					ELSE RETURN, self.oTime -> _OverloadBracketsRightSide(isRange[1:*], i2, i3)
			ENDCASE
			'TIME': BEGIN
				IF nSubs EQ 1 $
					THEN RETURN, self.oTime['DATA'] $
					ELSE RETURN, self.oTime -> _OverloadBracketsRightSide( [0,isRange[1:*]], 'DATA', i2, i3)
			ENDCASE
			ELSE: RETURN, self -> GetAttrValue(i1)
		ENDCASE

	;Scalar operations
	;   - 0   returns the self object
	;   - [0] returns the first data element
	;   - All other cases RETURN data
	ENDIF ELSE IF nSubs EQ 1 && isRange[0] EQ 0 && IsA(i1, /SCALAR) && i1 EQ 0 THEN BEGIN
		RETURN, self
	ENDIF

;---------------------------------------------------------------------
; Extract the Subarray ///////////////////////////////////////////////
;---------------------------------------------------------------------
	;Data
	data = self -> _OverloadBracketsRightSide_Data(isRange, i1, i2, i3, i4, i5, i6, i7, i8)

	;Attributes
	outDims = Size(data, /DIMENSIONS)
	attrs   = self -> _OverloadBracketsRightSide_Attrs(outDims, isRange, i1, i2, i3, i4, i5, i6, i7, i8)
	
;---------------------------------------------------------------------
; Create Output //////////////////////////////////////////////////////
;---------------------------------------------------------------------
	;MrTimeSeries
	;   - Scalars have a dimension size of 0 with 1 element.
	nDep0 = N_Elements( attrs['DEPEND_0'] )
	IF (outDims[0] EQ nDep0) || (N_Elements(data) EQ 1 && nDep0 EQ 1) THEN BEGIN
		vOut = MrTimeSeries( attrs['DEPEND_0'], data, $
		                     NAME='OverloadBRS(' + self.name + ')', $
		                     /NO_COPY )
	
	;MrVariable
	ENDIF ELSE BEGIN
		;Issue warning
		MrPrintF, 'LogWarn', 'Unable to create MrTimeSeries object. Converting to MrVariable.'
		
		;Create variable
		vOut = MrVariable( data, $
		                   NAME='OverloadBRS(' + self.name + ')', $
		                   /NO_COPY )
	ENDELSE
	
	;Copy all attributes
	self -> CopyAttrTo, vOut
	
	;Update attributes
	vOut -> SetAttrValue, attrs, /OVERWRITE
	
;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------

	RETURN, vOut
END


;+
;   Raise one expression to the power of another.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadCaret, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'Caret(' + left.name + ',' + right.name + ')'

		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) ^ (right -> GetData()) $
			else result = (right -> GetData()) ^ (*self.data)
		
		;Return object
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Caret(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) ^ right
			
		;RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Caret(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left ^ (*self.data)
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
		
	;Create a MrVariable object
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Perform an EQUALS comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadEQ, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'EQ(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data EQ right -> GetData() $
			ELSE result = left -> GetData() EQ *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'EQ(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) EQ right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'EQ(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left EQ (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   The purpose of this method is to provide iterations to the FOREACH operator. Iterate
;   over all of the elements in the array.
;
; :Params:
;       VALUE:              out, required, type=vector
;                           Set to the value of the current object element. If `KEY` is
;                               undefined, then set this to the first object element.
;       KEY:                out, required, type=long
;                           Set to the index (or key) associated with the current element.
;                               In the first iteration, KEY is undefined.
;
; :Returns:
;       NEXT:               Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
function MrTimeSeries::_OverloadForeach, value, key
	compile_opt idl2
	on_error, 2

	nPts = n_elements( (*self.data)[*,0] )
	if n_elements(key) eq 0 then key = 0

	;Get the array element if the index is in range
	if key lt nPts then begin
		next  = 1
		nDims = size(*self.data, /DIMENSIONS)
		
		case nDims of
			1: value = (*self.data)[key]
			2: value = reform((*self.data)[key,*])
			3: value = reform((*self.data)[key,*,*])
			4: value = reform((*self.data)[key,*,*,*])
			5: value = reform((*self.data)[key,*,*,*,*])
			6: value = reform((*self.data)[key,*,*,*,*,*])
			7: value = reform((*self.data)[key,*,*,*,*,*,*])
			8: value = reform((*self.data)[key,*,*,*,*,*,*,*])
			else: message, 'Invalid number of dimensions (' + strtrim(nDims, 2) + ').'
		endcase
	
	;Otherwise, stop iterating 
	endif else next = 0

	;Next element to retrieve
	key += 1

	return, next
end


;+
;   Perform an GREATER THAN OR EQUAL TO comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadGE, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'GE' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data GE right -> GetData() $
			ELSE result = left -> GetData() GE *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) EQ self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'GE(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) GE right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) EQ self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'GE' + dims + ',' + self.name + ')'
			
			;Operation
			result = left GE (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   Add one expression to another.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadGreaterThan, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'GreaterThan(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) > (right -> GetData()) $
			else result = (right -> GetData()) > (*self.data)
		
		;Return object
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;SELF is LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'GreaterThan(' + self.name + ',' + dims + ')'

			;Operation
			result = *self.data > right
			
		;SELF is RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'GreaterThan(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left > *self.data
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	
	;Create a MrVariable object
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Perform an GREATER THAN comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadGT, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'GT(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data GT right -> GetData() $
			ELSE result = left -> GetData() GT *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'GT(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) GT right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'GT(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left GT (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   Perform an LESS THAN OR EQUAL TO comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadLE, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'LE' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data LE right -> GetData() $
			ELSE result = left -> GetData() LE *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) EQ self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'LE(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) LE right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) EQ self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'LE' + dims + ',' + self.name + ')'
			
			;Operation
			result = left LE (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   Add one expression to another.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadLessThan, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT MrVariable data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'LessThan(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) < (right -> GetData()) $
			else result = (right -> GetData()) < (*self.data)
		
		;Output object class
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;SELF is LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'LessThan(' + self.name + ',' + dims + ')'

			;Operation
			result = *self.data < right
			
		;SELF is RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'LessThan(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left < *self.data
		endelse
		
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a MrVariable object
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Perform an LESS THAN comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadLT, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'LT' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data LT right -> GetData() $
			ELSE result = left -> GetData() LT *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) EQ self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'LT(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) LT right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) EQ self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'LT' + dims + ',' + self.name + ')'
			
			;Operation
			result = left LT (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   Raise one expression to the power of another.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadMinus, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'Minus(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) - (right -> GetData()) $
			else result = (right -> GetData()) - (*self.data)
		
		;Output object class
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Minus(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) - right
			
		;RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Minus(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left - (*self.data)
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a MrVariable object
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Take the mod of two expressions.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadMOD, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT MrVariable data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'Mod(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) mod (right -> GetData()) $
			else result = (right -> GetData()) mod (*self.data)
		
		;Output object class
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Mod(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) mod right
			
		;RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Mod(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left mod (*self.data)
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a MrVariable object
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Perform an NOT EQUAL TO comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadNE, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'NE(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data NE right -> GetData() $
			ELSE result = left -> GetData() NE *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'NE(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) NE right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'NE(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left NE (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   Perform a bit-wise OR comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadOr, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'OR(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data OR right -> GetData() $
			ELSE result = left -> GetData() OR *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'OR(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) OR right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'OR(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left OR (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   Add one expression to another.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadPlus, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT MrVariable data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)
	
;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'Plus(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) + (right -> GetData()) $
			else result = (right -> GetData()) + (*self.data)
		
		;Output object class
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;SELF is LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Plus(' + self.name + ',' + dims + ')'

			;Operation
			result = *self.data + right
			
		;SELF is RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Plus(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left + *self.data
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a MrVariable object
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Divide one expression by another.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadPound, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT MrVariable data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'Pound(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) # (right -> GetData()) $
			else result = (right -> GetData()) # (*self.data)
		
		;Output object class
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Pound(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) # right
			
		;RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Pound(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left # (*self.data)
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Take the outer product of two variables.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadPoundPound, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT MrVariable data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'PoundPound(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) ## (right -> GetData()) $
			else result = (right -> GetData()) ## (*self.data)
		
		;Output object class
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'PoundPound(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) ## right
			
		;RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'PoundPound(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left ## (*self.data)
		endelse
		
		;Create a MrVariable object
		if classOut eq 'MrVariable' $
			then result = obj_new( classOut, result, NAME=name, /NO_COPY ) $
			else result = obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   Divide one expression by another.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             The result of multiplying `LEFT` by `RIGHT`.
;-
function MrTimeSeries::_OverloadSlash, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT MrVariable data objects?
	side = self -> _LeftOrRight(left, right, ISMRVARIABLE=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;oVar -> GetData() is faster than oVar['DATA']
	if isa(left, 'MrVariable') && isa(right, 'MrVariable') then begin
		
		;New name
		name = 'Divide(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		if side eq 'LEFT' $
			then result = (*self.data) / (right -> GetData()) $
			else result = (right -> GetData()) / (*self.data)
		
		;Output object class
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	endif else begin
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		if side eq 'LEFT' then begin
			;Object class for output
			if MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(right) eq 1 $
				then dims = strtrim(right[0], 2) $
				else dims = '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Divide(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) / right
			
		;RIGHT
		endif else begin
			;Object class for output
			if MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				then classOut = obj_class(self) $
				else classOut = 'MrVariable'
			
			;Name of output variable
			if n_elements(left) eq 1 $
				then dims = strtrim(left[0], 2) $
				else dims = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Divide(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left / (*self.data)
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a MrVariable object
	if classOut eq 'MrVariable' $
		then return, obj_new( classOut, result, NAME=name, /NO_COPY ) $
		else return, obj_new( classOut, self.oTime, result, NAME=name, /NO_COPY)
end


;+
;   The purpose of this method is to negate the implicit array.
;
; :Returns:
;       RESULT:             Negation of the implicit array.
;-
function MrTimeSeries::_OverloadMinusUnary
	compile_opt idl2
	on_error, 2

	;Negate the array, making positive values negative, and vice versa
	return, self -> New( self.oTime, -(*self.data), NAME='-'+self.name)
end


;+
;   Perform a bit-wise XOR comparison.
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Returns:
;       RESULT:             Returns 1 if there is a current element to retrieve, and 0 if
;                               there are no more elements.
;-
FUNCTION MrTimeSeries::_OverloadXOr, left, right
	Compile_Opt idl2
	On_Error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	IF IsA(left, 'MrVariable') && IsA(right, 'MrVariable') THEN BEGIN
		name = 'XOR(' + left.name + ',' + right.name + ')'
		
		;Perform operation
		IF side EQ 'left' $
			THEN result = *self.data XOR right -> GetData() $
			ELSE result = left -> GetData() XOR *self.data
		
		;Type of output
		if obj_isa(left, 'MrTimeSeries') && obj_isa(right, 'MrTimeSeries') $
			then classOut = 'MrTimeSeries' $
			else classOut = 'MrVariable'

;-------------------------------------------------------
; MrTimeSeries with Expression /////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		;Rules for output
		;  Create time series if:
		;     1. Expression is a scalar
		;     2. Expression's first dimension is same size as that of implicit self
		;  Create MrVariable in all other cases
		
		;LEFT
		IF side EQ 'LEFT' THEN BEGIN
			;Object class for output
			IF MrIsA(right, /SCALAR) || MrDim(right, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(right) eq 1 $
				THEN dims = StrTrim(right[0], 2) $
				ELSE dims = '[' + StrJoin(StrTrim(Size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'XOR(' + self.name + ',' + dims + ')'
			
			;Operation
			result = (*self.data) XOR right
			
		;RIGHT
		endif else begin
			;Object class for output
			IF MrIsA(left, /SCALAR) || MrDim(left, 1) eq self -> GetNPts() $
				THEN classOut = Obj_Class(self) $
				ELSE classOut = 'MrVariable'
			
			;Name of output variable
			IF N_Elements(left) EQ 1 $
				THEN dims = StrTrim(left[0], 2) $
				ELSE dims = Size(left, /TNAME) + '[' + StrJoin(StrTrim(Size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'XOR(' + dims + ',' + self.name + ')'
			
			;Operation
			result = left XOR (*self.data)
		ENDELSE
	ENDELSE
;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	IF classOut EQ 'MrVariable' $
		THEN RETURN, Obj_New(classOut, result, /NO_COPY, NAME=name) $
		ELSE RETURN, Obj_New(classOut, self.oTime, result, /NO_COPY, NAME=name)
END


;+
;   Concatenate an array, or a series of arrays, to the implicit array. If a series of
;   arrays are to be concatenated, as in a FOR loop, memory allocation is handled
;   internally by doubling the current number of elements or by adding the number of
;   elements present in the input array, whichever is greater. Truncation of unfilled
;   elements is also performed internally. See the `START` and `FINISH` methods, as well
;   as the examples below.
;
;   Notes:
;       - Extra elements are not truncated until ::Append is called with the /FINISH
;           keyword set.
;       - If ::Append was initialized with the /START and /BEFORE keywords, elements are
;           appended to the end of the array. When ::Append is called with the /FINISH
;           keyword set, they are shifted to the front.
;
; :Params:
;       TIME:           in, required, type=Nx1 or 1xN array
;                       Array to be concatenated to the implicit array. If undefined or
;                           !Null, then nothing is appended.
;       DATA:           in, required, type=array
;                       Array to be concatenated to the implicit array. If undefined or
;                           !Null, then nothing is appended.
;       TYPE:           in, optional, type=array, default='ISO-8601'
;                       Units of the input time array. See ::toISO for options.
;
; :Keywords:
;       BEFORE:         in, optional, type=boolean, default=0
;                       If set, `DATA` will be appended to the beginning of the implicit
;                           array. The default is to append to the end.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       If set, `DATA` will be copied directly into the object and
;                           will be left undefined.
;-
pro MrTimeSeries::Append, time, data, type, $
BEFORE=before, $
NO_COPY=no_copy
	compile_opt idl2
	on_error, 2
	
	;Allow input data to be empty. In this case, do nothing.
	if n_elements(data) eq 0 then return

	;Defaults
	before  = keyword_set(before)
	no_copy = keyword_set(no_copy)

	;Restrictions
	nPts = n_elements(time)
	dims = size(data, /DIMENSIONS)
	if dims[0] ne nPts then message, 'Size of TIME and first dimension of DATA must agree.'
	
;-----------------------------------------------------
; Append Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Is the object empty?
	if n_elements(*self.data) eq 0 then begin
		self -> SetData, time, data, NO_COPY=no_copy
	
	;Append
	endif else begin
		;Append to time array
		self.oTime -> Append, time, type, BEFORE=before, NO_COPY=no_copy
	
		;Append to data array
		if before then begin
			if no_copy $
				then *self.data = [ temporary(data), *self.data ] $
				else *self.data = [ data, *self.data ]
		endif else begin
			if no_copy $
				then *self.data = [ *self.data, temporary(data) ] $
				else *self.data = [ *self.data, data ]
		endelse
	endelse
end


;+
;   Create a copy of the variable with a new heap identifier. Time, data,
;   and attributes are copied.
;
; :Params:
;       NAME:           in, optional, type=string, default=<name>+'_copy'
;                       Name of the variable copy.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the copy will be added to the variable cache
;       NO_CLOBBER:     in, optional, type=boolean, default=0
;                       If used with `CACHE` and a variable by the same name already
;                           exists in the cache, `NAME` will have "_#" appended to it,
;                           where # is the lowest integer that makes `NAME` unique.
;                           The default is to delete and replace the variable in the cache.
;                          
;-
function MrTimeSeries::Copy, name, $
CACHE=cache, $
NO_CLOBBER=no_clobber
	compile_opt idl2
	on_error, 2

	;Defaults
	tf_cache = keyword_set(cache)
	if n_elements(name) eq 0 then name = 'Copy(' + self.name + ')'

	;Copy the data into a new object
	;   - Use Obj_Class() so subclasses can inherit the method.
	theCopy = obj_new( obj_class(self), self.oTime, *self.data, $
	                   DIMENSION = 1, $
	                   NAME      = name, $
	                   T_NAME    = self.oTime.name )

	;Copy the variable attributes
	self -> CopyAttrTo, theCopy
	
	;Cache the variable
	if tf_cache then theCopy -> Cache, NO_CLOBBER=no_clobber
	return, theCopy
end


;+
;   Detrend data by removing a smoothed version of itself. Acts as a high-pass filter.
;
; :Params:
;       WIDTH:          in, optional, type=integer, default=all
;                       Number of points to use in the boxcar smooth.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the copy will be added to the variable cache
;       NAME:           in, optional, type=string, default=<name>+'_copy'
;                       Name of the variable copy.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by IDL's Smooth function.
;
; :Returns:
;       RESULT:         out, required, type=MrTimeSeries or MrVariable object
;                       The results of the fast fourier transformation. If more than
;                           one FFT is possible, a MrTimeSeries object is returned.
;                           Otherwise, a MrVariable object is returned.
;-
function MrTimeSeries::Detrend, width, $
CACHE=cache, $
NAME=name, $
_REF_EXTRA=extra
	compile_opt idl2
	on_error, 2
	
	;Defaults
	if n_elements(name) eq 0 then name = 'Detrend(' + self.name + ')'
	
	;Smoothing width along each dimension -- only apply to first dimension
	nDims = size(self, /N_DIMENSIONS)
	if nDims le 1 $
		then nSmooth = width $
		else nSmooth = [width, indgen(nDims-1)+1S]
	
	;Remove smoothed data
	data = self - smooth(self['DATA'], nSmooth, _STRICT_EXTRA=extra)
	
	;Set properties
	self -> CopyAttrTo, data
	data['CATDESC'] = self.name + ' detrended with ' + strtrim(width) + ' point boxcar average.'

	;Cache
	data -> SetName, name
	if keyword_set(cache) then data -> Cache
	return, data
end


;+
;   High-pass, low-pass, band-pass, or band-stop filter data. The filter is created
;   using IDL's DIGITAL_FILTER function and applied using IDL's CONVOL function.
;
; :Params:
;       FLOW:           in, required, type=Numeric array
;                       Lower frequency of the filter as a fraction of the nyquist.
;       FHIGH:          in, required, type=Numeric array
;                       Upper frequency of the filter as a fraction of the nyquist.
;       A:              in, required, type=Numeric array
;                       The filter power relative to the Gibbs phenomenon wiggles in
;                           decibels. 50 is a good choice.
;       NTERMS:         in, required, type=Numeric array
;                       The number of terms used to construct the filter.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the result will be added to the cache.
;       DOUBLE:         in, optional, type=boolean, default=0
;                       If set, computations are performed in double-precision.
;       NAME:           in, optional, type=string, default='Digital_Filter(<name>)'
;                       A name to be given to the variable.
;
; :Returns:
;       VAROUT:         out, required, type=object
;                       A MrTimeSeries object containing the filtered data.
;-
FUNCTION MrTimeSeries::Digital_Filter, fLow, fHigh, A, nTerms, $
CACHE=cache, $
DOUBLE=double, $
NAME=name
	compile_opt idl2
	on_error, 2
	
	IF N_Elements(name) EQ 0 THEN name = 'Digital_Filter(' + self.name + ')'
	
	;Create the filter
	coeff = Digital_Filter(fLow, fHigh, A, nTerms, DOUBLE=double)
	
	;Reform to NxM
	dims  = Size(self, /DIMENSIONS)
	nDims = N_Elements(dims)
	IF nDims GT 2 $
		THEN temp = Reform( (*self.data), [dims[0], Product(dims[1:*])] ) $
		ELSE temp = *self.data
	
	;Loop over second dimension
	nLoop = nDims LE 1 ? 1 : Product(dims[1:*])
	FOR i = 0, nLoop - 1 DO temp[0,i] = Convol( temp[*,i], coeff )
	
	;Reform back to orignal size
	IF nDims GT 2 THEN temp = Reform( temp, dims, /OVERWRITE )
	
	;Create the output variable
	theClass = Obj_Class(self)
	varOut = Obj_New( theClass, self.oTime, temp, $
	                  CACHE     = cache, $
	                  NAME      = name, $
	                  /NO_COPY )
	
	;Return the new variable
	RETURN, varOut
END

;+
;   Store variables as TPlot variables.
;
;   NOTE:
;       Requires the SPEDAS distribution.
;       http://spedas.org/blog/
;
; :Returns:
;       TNAME:              out, required, type=string
;                           Name of the TPlot variable created.
;-
FUNCTION MrTimeSeries::ExportToSPEDAS
	Compile_Opt idl2
	On_Error, 2
	
;-----------------------------------------------------
; Create the Data Structure \\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	CASE 1 OF
		;2D Image
		self -> HasAttr('DEPEND_1'): BEGIN
			data = { x: time_double( self.oTime['DATA'] ), $
			         y: *self.data, $
			         v: (self['DEPEND_1'])['DATA'] }
		ENDCASE
		
		;1D Time Series
		self -> HasAttr('DEPEND_0'): BEGIN
			data = { x: time_double( self.oTime['DATA'] ), $
			         y: *self.data }
		ENDCASE
		
		;Unknown
		ELSE: Message, 'I do not know how to load data with > 3 dimensions into TPlot.'
	ENDCASE
	
;-----------------------------------------------------
; Store the Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Store the data
	tname = self.name
	store_data, tname, DATA=temporary(data)
	
;-----------------------------------------------------
; X-Axis Options \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF self -> HasAttr('DEPEND_0') THEN BEGIN
		dep0 = MrVar_Get( self['DEPEND_0'] )
		IF dep0 -> HasAttr('CHARSIZE')     THEN options, tname, 'xcharsize',     dep0['CHARSIZE']
		IF dep0 -> HasAttr('GRIDSTYLE')    THEN options, tname, 'xgridstyle',    dep0['GRIDSTYLE']
		IF dep0 -> HasAttr('LOG')          THEN options, tname, 'xlog',          dep0['LOG']
		IF dep0 -> HasAttr('MINOR')        THEN options, tname, 'xminor',        dep0['MINOR']
		IF dep0 -> HasAttr('AXIS_RANGE')   THEN options, tname, 'xrange',        dep0['AXIS_RANGE']
		IF dep0 -> HasAttr('STYLE')        THEN options, tname, 'xstyle',        dep0['STYLE']
		IF dep0 -> HasAttr('THICK')        THEN options, tname, 'xthick',        dep0['THICK']
		IF dep0 -> HasAttr('TICKFORMAT')   THEN options, tname, 'xtickformat',   dep0['TICKFORMAT']
		IF dep0 -> HasAttr('TICKINTERVAL') THEN options, tname, 'xtickinterval', dep0['TICKINTERVAL']
		IF dep0 -> HasAttr('TICKLAYOUT')   THEN options, tname, 'xticklayout',   dep0['TICKLAYOUT']
		IF dep0 -> HasAttr('TICKLEN')      THEN options, tname, 'xticklen',      dep0['TICKLEN']
		IF dep0 -> HasAttr('TICKNAME')     THEN options, tname, 'xtickname',     dep0['TICKNAME']
		IF dep0 -> HasAttr('MAJOR')        THEN options, tname, 'xticks',        dep0['MAJOR']
		IF dep0 -> HasAttr('TICKUNITS')    THEN options, tname, 'xtickunits',    dep0['TICKUNITS']
		IF dep0 -> HasAttr('TICKV')        THEN options, tname, 'xtickv',        dep0['TICKV']
		IF dep0 -> HasAttr('TITLE')        THEN options, tname, 'xtitle',        dep0['TITLE']
	ENDIF
	
;-----------------------------------------------------
; Y-Axis Options \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Image or Time series?
	IF self -> HasAttr('DEPEND_1') THEN BEGIN
		yvar = MrVar_Get( self['DEPEND_1'] )
		options, tname, 'spec', 1B
	ENDIF ELSE BEGIN
		yvar = self
	ENDELSE

	IF yvar -> HasAttr('CHARSIZE')     THEN options, tname, 'ycharsize',     yvar['CHARSIZE']
	IF yvar -> HasAttr('GRIDSTYLE')    THEN options, tname, 'ygridstyle',    yvar['GRIDSTYLE']
	IF yvar -> HasAttr('LOG')          THEN options, tname, 'ylog',          yvar['LOG']
	IF yvar -> HasAttr('MINOR')        THEN options, tname, 'yminor',        yvar['MINOR']
	IF yvar -> HasAttr('AXIS_RANGE')   THEN options, tname, 'yrange',        yvar['AXIS_RANGE']
	IF yvar -> HasAttr('STYLE')        THEN options, tname, 'ystyle',        yvar['STYLE']
	IF yvar -> HasAttr('THICK')        THEN options, tname, 'ythick',        yvar['THICK']
	IF yvar -> HasAttr('TICKFORMAT')   THEN options, tname, 'ytickformat',   yvar['TICKFORMAT']
	IF yvar -> HasAttr('TICKINTERVAL') THEN options, tname, 'ytickinterval', yvar['TICKINTERVAL']
	IF yvar -> HasAttr('TICKLAYOUT')   THEN options, tname, 'yticklayout',   yvar['TICKLAYOUT']
	IF yvar -> HasAttr('TICKLEN')      THEN options, tname, 'yticklen',      yvar['TICKLEN']
	IF yvar -> HasAttr('TICKNAME')     THEN options, tname, 'ytickname',     yvar['TICKNAME']
	IF yvar -> HasAttr('MAJOR')        THEN options, tname, 'yticks',        yvar['MAJOR']
	IF yvar -> HasAttr('TICKUNITS')    THEN options, tname, 'ytickunits',    yvar['TICKUNITS']
	IF yvar -> HasAttr('TICKV')        THEN options, tname, 'ytickv',        yvar['TICKV']
	IF yvar -> HasAttr('TITLE')        THEN options, tname, 'ytitle',        yvar['TITLE']
	
;-----------------------------------------------------
; General Options \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF self -> HasAttr('COLOR')      THEN options, tname, 'colors', cgColor(self['COLOR'])
	IF self -> HasAttr('PLOT_TITLE') THEN options, tname, 'title', self['PLOT_TITLE']
	IF self -> HasAttr('LABEL') THEN BEGIN
		options, tname, 'labels',  self['LABEL']
		options, tname, 'labflag', 1
	ENDIF
	
;-----------------------------------------------------
; Done \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	RETURN, tname
END


;+
;   Compute the FFT
;
; :Params:
;       NFFT:           in, optional, type=integer, default=all
;                       Number of points per FFT.
;       NSHIFT:         in, optional, type=integer, default=`NFFT`/2.0
;                       Number of points to shift between FFTs.
;
; :Params:
;       AMPLITUDE:      out, optional, type=MrTimeVar/MrVariable object
;                       The amplitude of the signal.
;       PHASE:          out, optional, type=MrTimeVar/MrVariable object
;                       The phase of the signal.
;
; :Returns:
;       RESULT:         out, required, type=MrTimeSeries or MrVariable object
;                       The results of the fast fourier transformation. If more than
;                           one FFT is possible, a MrTimeSeries object is returned.
;                           Otherwise, a MrVariable object is returned.
;-
function MrTimeSeries::FFT_Old, nfft, nshift, $
DOUBLE=double, $
NOISE=noise, $
ONE_SIDED=one_sided, $
POWER=power, $
RMS=rms, $
WINDOW=win
	compile_opt idl2
	on_error, 2
	
	dim          = 1
	tf_double    = keyword_set(double)
	tf_noise     = keyword_set(noise)
	tf_one_sided = keyword_set(one_sided)
	tf_power     = keyword_set(power)
	tf_rms       = keyword_set(rms)
	if n_elements(nfft)   eq 0 then nfft   = self -> GetNPts()
	if n_elements(nshift) eq 0 then nshift = nfft / 2.0
	if n_elements(win)    eq 0 then win    = ''
	
	;Conflicts
	if tf_rms + tf_noise + tf_power gt 1 then message, 'RMS, NOISE, and POWER are mutually exclusive.'

;-------------------------------------------------------
; Tapering Window //////////////////////////////////////
;-------------------------------------------------------
	
	;Name
	if isa(win, 'STRING', /SCALAR) then begin
		case strupcase(win) of
			'BLACKMAN': begin
				xx    = findgen(nfft)/(nfft-1)
				taper = 0.42 - 0.50 * cos(2*!pi*xx) + 0.08 * cos(4*!pi*xx)
			endcase
			'RECTANGULAR': taper = replicate(1.0, nfft)
			'HANNING':     taper = hanning(nfft, ALPHA=0.50)
			'HAMMING':     taper = hanning(nfft, ALPHA=0.54)
			'NONE':        taper = replicate(1.0, nfft)
			'':            taper = replicate(1.0, nfft)
			else: message, 'Window not recognized: "' + win + '".'
		endcase
	
	;User-defined
	endif else begin
		if n_elements(win) ne nfft then message, 'WIN must have NFFT elements.'
		taper = win
	endelse
	
	;Correct dimensions
	dims  = size(self, /DIMENSIONS)
	ndims = size(self, /N_DIMENSIONS)
	if ndims gt 1 then taper = rebin( taper, [nfft, dims[1:*]] )

;-------------------------------------------------------
; Data Intervals ///////////////////////////////////////
;-------------------------------------------------------
	;Convert to seconds
	t  = self.oTime -> GetData('SSM')
	
	;Number of FFTs to perform
	N    = dims[0]
	nMax = floor( (N - nfft) / nshift ) + 1
	
;-----------------------------------------------------
; Step Through Intervals \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Allocate memory
	;   - Only necessary of more than one FFT is being performed
	type = tf_double ? 9 : 6
	if nMax gt 1 then begin
		if nDims gt 7 then message, 'DIMENSIONS must be < 7 for sliding FFTs.'
		outDims   = nDims eq 1 ? [nMax, nfft] : [nMax, nfft, dims[1:*]]
		dft       = make_array( outDims, TYPE=type )
		freqs     = fltarr(nMax,nfft)
		dt_median = fltarr(nMax)
		df_fft    = fltarr(nMax)
		t_fft     = fltarr(nMax)
		dt_fft    = fltarr(nMax)
	endif

	for i = 0, nMax-1 do begin
	;-----------------------------------------------------
	; Define Interval \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;
		; Place at beginning of loop so that CONTINUE
		; advances interval as well as loop index.
		;
	
		;First interval
		if i eq 0 then begin
			istart = 0
			iend   = nfft - 1
		
		;Last interval
		;   - Extend backward from the end of the array
		endif else if i eq nMax-1 then begin
			ilast = iend
			iend  = N - 1
			istart = N - nfft
		
		;Intermediate intervals
		;   - Advance forward by NSHIFT
		endif else begin
			istart += nshift
			iend   += nshift
		endelse
	
	;-----------------------------------------------------
	; Compute Sample Interval \\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Time and sampling interval
		dtTemp = t[istart+1:iend] - t[istart:iend-1]

		;Is the sampling interval consistent throughout?
		;   - A "bad" sampling interval occurs if it is 10% different
		;     from the first sampling interval
		SI = median(dtTemp)
		if total( abs(temporary(dtTemp) - SI) gt 0.1*SI ) ne 0 then begin
			message, 'Sampling interval changes during FFT. Skipping.', /INFORMATIONAL
			continue
		endif
	
	;-----------------------------------------------------
	; Compute FFT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Determine frequencies
		sr    = 1.0 / SI
		fN    = sr / 2.0
		df    = 1.0 / (SI * nfft)
		f     = [0:fN:df]
		f     = nfft MOD 2 EQ 0 ? [ f, -reverse(f[1:-2]) ] : [ f, -reverse(f[1:-1]) ]
	
		;Compute the FFT
		case ndims of
			1: fftTemp = fft( taper * (*self.data)[istart:iend], DIMENSION=1, DOUBLE=tf_double)
			2: fftTemp = fft( taper * (*self.data)[istart:iend,*], DIMENSION=1, DOUBLE=tf_double)
			3: fftTemp = fft( taper * (*self.data)[istart:iend,*,*], DIMENSION=1, DOUBLE=tf_double)
			4: fftTemp = fft( taper * (*self.data)[istart:iend,*,*,*], DIMENSION=1, DOUBLE=tf_double)
			5: fftTemp = fft( taper * (*self.data)[istart:iend,*,*,*,*], DIMENSION=1, DOUBLE=tf_double)
			6: fftTemp = fft( taper * (*self.data)[istart:iend,*,*,*,*,*], DIMENSION=1, DOUBLE=tf_double)
			7: fftTemp = fft( taper * (*self.data)[istart:iend,*,*,*,*,*,*], DIMENSION=1, DOUBLE=tf_double)
			8: fftTemp = fft( taper * (*self.data)[istart:iend,*,*,*,*,*,*,*], DIMENSION=1, DOUBLE=tf_double)
			else: message, 'Incorrect number of dimensions: SELF.'
		endcase

		;Normalize the FFT to account for the window
;		CG = mean(taper)
;		NG = mean(taper^2)
;		case 1 of
;			tf_rms:   fftTemp /= CG
;			tf_noise: ;Do nothing
;			tf_power: ;Do nothing
;			else:     ;Do nothing
;		endcase

	;-----------------------------------------------------
	; Outputs \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Save data
		if nMax eq 1 then begin
			dft   = temporary(fftTemp)
			freqs = temporary(f)
		endif else begin
			case ndims of
				1: dft[i,*] = temporary(fftTemp)
				2: dft[i,*,*] = temporary(fftTemp)
				3: dft[i,*,*,*] = temporary(fftTemp)
				4: dft[i,*,*,*,*] = temporary(fftTemp)
				5: dft[i,*,*,*,*,*] = temporary(fftTemp)
				6: dft[i,*,*,*,*,*,*] = temporary(fftTemp)
				7: dft[i,*,*,*,*,*,*,*] = temporary(fftTemp)
				else: message, 'Incorrect number of dimensions: SELF.'
			endcase
			
			;Save frequencies and time itnerval
			;   - Necessary only for I=0 if TF_CALC_DT is false
			;   - Keep track of all here to simplify logic
			freqs[i,*]   = f
			dt_median[i] = SI
			dt_fft[i]    = nfft * SI
			df_fft[i]    = 1.0 / (nfft * SI)
		
			;Record time stamp
			t_fft[i]  = t[istart] + (nfft * SI) / 2.0
			dt_fft[i] = nfft * SI / 2.0
		endelse
	endfor

;-----------------------------------------------------
; Finishing Touches \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Reduce dimensions if sampling interval never changed
	;   - fN = 1 / T = 1 / (nfft*dt)
	;   - NFFT is constant, but DT is permitted to vary
	;   - If DT changes, then we must keep track of FN as a function of time
	;   - If not, then FN and F do not vary with time, and we can reduce to one copy
	if nMax gt 1 && total( (dt_median - dt_median[0]) gt 0.1*dt_median[0] ) eq 0 then begin
		freqs     = reform(freqs[0,*])
		df_fft    = df_fft[0]
		dt_median = dt_median[0]
		dt_fft    = dt_fft[0]
	endif
	
	;Remove negative frequencies
	;   - TODO: What if DIMENSION ~= 2
	if tf_one_sided then begin
		freqs = freqs[0:nfft/2]
		case ndims of
			1: dft = dft[*,0:nfft/2]
			2: dft = dft[*,0:nfft/2,*]
			3: dft = dft[*,0:nfft/2,*,*]
			4: dft = dft[*,0:nfft/2,*,*,*]
			5: dft = dft[*,0:nfft/2,*,*,*,*]
			6: dft = dft[*,0:nfft/2,*,*,*,*,*]
			7: dft = dft[*,0:nfft/2,*,*,*,*,*,*]
			else: message, 'Incorrect number of dimensions: SELF.'
		endcase
	endif

;-----------------------------------------------------
; Output Object \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;TIME
	if nMax gt 1 then oTime = MrTimeVar( t_fft, 'SSM', /NO_COPY, T_REF=self.oTime['DATA', [0]] )

	;FREQUENCY
	if size(freqs, /DIMENSIONS) eq 2 $
		then oFreq = MrTimeSeries( oTime, freqs, /NO_COPY ) $
		else oFreq = MrVariable( freqs, /NO_COPY )

	;DFT
	if nMax gt 1 $
		then oDFT = MrTimeSeries( oTime, dft ) $
		else oDFT = MrVariable( dft, /NO_COPY )

;-----------------------------------------------------
; Attributes \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;DEPEND_# & DELTA_(MINUS|PLUS)
	if nMax gt 1 then begin
		;dt - Scalar
		if n_elements(dt_fft) eq 1 then begin
			oTime['DELTA_PLUS']  = dt_fft
			oTime['DELTA_MINUS'] = dt_fft
		
		;dt(t)
		endif else begin
			;Create a time-dependent variation
			oDT = MrTimeSeries( oTime, dt_fft, /NO_COPY )
			oDT['TITLE'] = 'Deviation from center time.'
			oDT['UNITS'] = 'seconds'
			
			;Add as attribute to TIME
			oTime['DELTA_PLUS_VAR']  = oDT
			oTime['DELTA_MINUS_VAR'] = oDT
		endelse
		
		;Set DFT dependencies
		oDFT['DEPEND_0']   = oTime
		oDFT['DEPEND_1']   = oFreq
		oDFT['PLOT_TITLE'] = string(nfft, nshift, FORMAT='(%"DFT: NFFT=%i, NSHIFT=%i")')
	
	;Single DFT
	endif else begin
		oDFT['DEPEND_0']   = oFreq
		oDFT['PLOT_TITLE'] = string(nfft, FORMAT='(%"DFT: NFFT=%i")')
	endelse
	
	;TIME
	if nMax gt 1 then begin
		oTime -> SetName, 'fft_tcenter'
		oTime['TITLE'] = 'FFT Center Times'
	endif
	
	;FREQUENCY
	oFreq -> SetName, 'fft_freq'
	oFreq['TITLE'] = 'Freq!C(Hz)'
	oFreq['UNITS'] = 'Hz'
	if tf_one_sided then begin
		oFreq['LOG']        = 1B
		oFreq['AXIS_RANGE'] = [oFreq['DATA',1], oFreq['DATA',-1]]
	endif
	
	;DFT
	oDFT -> SetName, 'FFT'
	oDFT['TITLE'] = 'FFT'
	oDFT['SCALE'] = 1
	
	return, oDFT
end


;+
;   Compute the FFT
;
; :Params:
;       NFFT:           in, optional, type=integer, default=all
;                       Number of points per FFT.
;       NSHIFT:         in, optional, type=integer, default=`NFFT`/2.0
;                       Number of points to shift between FFTs.
;
; :Params:
;       DETREND:        in, optional, type=boolean, default=0
;                       If set, each FFT interval will be fitted with a straight line
;                           which is then used to remove a linear trend in the data.
;       DOUBLE:         in, optional, type=boolean, default=0
;                       If set, the FFT is computed in double precision.
;       ONE_SIDED:      in, optional, type=boolean, default=0
;                       If set, the Fourier coefficients of negative frequencies are discarded.
;       NOISE:          in, optional, type=boolean, default=0
;                       If set, the FFT is normalized based on signal noise.
;       POWER:          in, optional, type=boolean, default=0
;                       If set, the FFT is normalized based on signal power.
;       RANGE:          in, optional, type=strarr(2)/intarr(2), default=[0\,N]
;                       A time or index range outlining a subinterval of data over which
;                           the FFT is to be computed. If a time range is given, it must
;                           be formatted as ISO-8601 strings.
;       RMS:            in, optional, type=boolean, default=0
;                       If set, the FFT is normalized based on RMS values.
;       WINDOW:         in, optional, type=fltarr/string, default='Rectangular'
;                       A tapering window or the name of a tapering window to be applied
;                           before the FFT is performed.
;       ZEROPAD:        in, optional, type=boolean/integer, default=0
;                       If set, the data will be zero-padded with N zeros, where N brings
;                           the total number of points up to the nearest power of 2. ZEROPAD
;                           can also be an integer > 1 indicating the number of zeros with
;                           which to pad the data.
;
; :Returns:
;       RESULT:         out, required, type=MrTimeSeries or MrVariable object
;                       The results of the fast fourier transformation. If more than
;                           one FFT is possible, a MrTimeSeries object is returned.
;                           Otherwise, a MrVariable object is returned.
;-
FUNCTION MrTimeSeries::FFT, nfft, nshift, $
DETREND=detrend, $
DOUBLE=double, $
NOISE=noise, $
ONE_SIDED=one_sided, $
POWER=power, $
RANGE=range, $
RMS=rms, $
WINDOW=win, $
ZEROPAD=zeropad
	Compile_Opt idl2
	On_Error, 2
	
	;Range
	irange = N_Elements(range) EQ 0 ? [0, (self -> GetNPts()) - 1] : range
	IF Size(irange, /TNAME) EQ 'STRING' THEN irange = self -> Value_Locate(range)
	nPts = irange[1] - irange[0] + 1
	
	dim          = 1
	tf_double    = Keyword_Set(double)
	tf_noise     = Keyword_Set(noise)
	tf_one_sided = Keyword_Set(one_sided)
	tf_power     = Keyword_Set(power)
	tf_rms       = Keyword_Set(rms)
	tf_zeropad   = Keyword_Set(zeropad)
	tf_detrend   = Keyword_Set(detrend)
	IF N_Elements(nfft)   EQ 0 THEN nfft   = nPts
	IF N_Elements(nshift) EQ 0 THEN nshift = nfft / 2
	
	;Conflicts
	IF tf_rms + tf_noise + tf_power GT 1 THEN Message, 'RMS, NOISE, and POWER are mutually exclusive.'
	IF nfft GT nPts THEN Message, 'Length of FFT is longer than data interval.'

;-------------------------------------------------------
; Extract Data /////////////////////////////////////////
;-------------------------------------------------------
	;Time
;	time = self.oTime['DATA', irange[0]:irange[1], 'SSM']
	
	;Data
	ndims = Size(*self.data, /N_DIMENSIONS)
;	CASE ndims OF
;		1: data = (*self.data)[irange[0]:irange[1]]
;		2: data = (*self.data)[irange[0]:irange[1],*]
;		3: data = (*self.data)[irange[0]:irange[1],*,*]
;		4: data = (*self.data)[irange[0]:irange[1],*,*,*]
;		5: data = (*self.data)[irange[0]:irange[1],*,*,*,*]
;		6: data = (*self.data)[irange[0]:irange[1],*,*,*,*,*]
;		7: data = (*self.data)[irange[0]:irange[1],*,*,*,*,*,*]
;		8: data = (*self.data)[irange[0]:irange[1],*,*,*,*,*,*,*]
;		ELSE: Message, 'Incorrect number OF dimensions: SELF.'
;	ENDCASE
	
	;New dimension sizes
	dims = Size(*self.data, /DIMENSIONS)

;-------------------------------------------------------
; Tapering Window //////////////////////////////////////
;-------------------------------------------------------
	;Default window
	IF N_Elements(win) EQ 0 THEN BEGIN
		taper = self -> FFT_Window(nfft)
	
	;Named window
	ENDIF ELSE IF MrIsA(win, /SCALAR, 'STRING') THEN BEGIN
		taper = self -> FFT_Window(nfft, win)
		
	;User-defined window
	ENDIF ELSE BEGIN
		;Check length
		IF tf_window && N_Elements(win) NE nfft THEN Message, 'WINDOW must have NFFT elements.'
		taper = win
	ENDELSE
	
	;Correct dimensions
	IF ndims GT 1 THEN taper = Rebin( taper, [nfft, dims[1:*]] )

;-----------------------------------------------------
; Step Through Intervals \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Declare variables
	nMax   = Floor( (nPts - nfft) / nshift ) + 1    ;Number of FFTs to perform
	N_new  = nfft                                   ;Total number of points in FFT after zero-padding
	nZeros = 0                                      ;Number of zeros to pad
	
	;Zero pad
	IF tf_zeropad THEN BEGIN
		;Number of zeros to add
		IF zeropad EQ 1 THEN BEGIN
			N_new  = 2^Ceil(ALog2(nfft))
			nZeros = N_new - nfft
			IF nZeros EQ 0 then tf_zeropad = 0B
		ENDIF ELSE BEGIN
			N_new  = nfft + zeropad
			nZeros = zeropad
		ENDELSE
	ENDIF
	
	;Allocate memory
	;   - Only necessary if more than one FFT is being performed
	type = tf_double ? 9 : 6
	IF nMax GT 1 THEN BEGIN
		IF nDims GT 7 THEN Message, 'DIMENSIONS must be < 7 for sliding FFTs.'
		outDims = nDims EQ 1 ? [nMax, N_new] : [nMax, N_new, dims[1:*]]
		dft     = Make_Array( outDims, TYPE=type )
	ENDIF ELSE BEGIN
		outDims = nDims EQ 1 ? N_new : [N_new, dims[1:*]]
		dft     = Make_Array( outDims, TYPE=type )
	ENDELSE
	freqs     = FltArr(nMax,N_new)
	dt_median = FltArr(nMax)
	df_fft    = FltArr(nMax)
	t_fft     = FltArr(nMax)
	dt_fft    = FltArr(nMax)
	
	;Step through each interval
	FOR i = 0, nMax-1 DO BEGIN
	;-----------------------------------------------------
	; Define Interval \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;
		; Place at beginning of loop so that continue
		; advances interval as well as loop index.
		; Use NFFT, not N_NEW, because we are indexing
		; into the non-zero-padded array.
		;
		
		;First interval
		IF i EQ 0 THEN BEGIN
;			istart = 0
			istart = irange[0]
			iend   = irange[0] + nfft - 1
		
		;Last interval
		;   - Extend backward from the end of the array
		ENDIF ELSE IF i EQ nMax-1 THEN BEGIN
			ilast  = iend
;			iend   = nPts - 1
			iend   = irange[1]
			istart = iend - nfft + 1
		
		;Intermediate intervals
		;   - Advance forward by NSHIFT
		ENDIF ELSE BEGIN
			istart += nshift
			iend   += nshift
		ENDELSE
	
	;-----------------------------------------------------
	; Compute Sample Interval \\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Time and sampling interval
;		ttemp  = time[istart:iend]
		ttemp  = self.oTime['DATA', istart:iend, 'SSM']
		dtTemp = ttemp[1:-1] - ttemp[0:-2]

		;Is the sampling interval consistent throughout?
		;   - A "bad" sampling interval occurs if it is 10% different
		;     from the first sampling interval
		SI = Median(dtTemp)
		IF Total( Abs(Temporary(dtTemp) - SI) GT 0.1*SI ) NE 0 THEN BEGIN
			MrPrintF, 'LogWarn', 'Sampling interval changes during FFT. Skipping.'
			CONTINUE
		ENDIF
	
	;-----------------------------------------------------
	; Select Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Extract the sub-interval
		CASE ndims OF
			1: temp = (*self.data)[istart:iend]
			2: temp = (*self.data)[istart:iend,*]
			3: temp = (*self.data)[istart:iend,*,*]
			4: temp = (*self.data)[istart:iend,*,*,*]
			5: temp = (*self.data)[istart:iend,*,*,*,*]
			6: temp = (*self.data)[istart:iend,*,*,*,*,*]
			7: temp = (*self.data)[istart:iend,*,*,*,*,*,*]
			8: temp = (*self.data)[istart:iend,*,*,*,*,*,*,*]
			ELSE: Message, 'Incorrect number OF dimensions: SELF.'
		ENDCASE

	;-----------------------------------------------------
	; Detrend \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		IF tf_detrend THEN BEGIN
			;Single fit
			IF ndims EQ 1 THEN BEGIN
				params = LADFit(ttemp, temp)
				temp   = temp - (params[0] + params[1]*ttemp)
			
			;Fit each dimension
			ENDIF ELSE BEGIN
				;Reform to 2 dimensions
				dim2 = Product(dims[1:*])
				temp = Reform(temp, [nfft, dim2], /OVERWRITE)
				
				;Detrend each component individually
				FOR j = 0, dim2 - 1 DO BEGIN
					params    = LADFit(ttemp, temp[*,j])
					temp[*,j] = temp - (params[0] + params[1]*ttemp)
				ENDFOR
				
				;Return to original size
				temp = Reform(temp, [nfft, dims[1:*]], /OVERWRITE)
			ENDELSE
		ENDIF
	
	;-----------------------------------------------------
	; Apply Window \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		temp = taper * temp
	
	;-----------------------------------------------------
	; Zero Pad \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		
		IF tf_zeropad THEN BEGIN
			IF ndims GT 1 $
				THEN temp = [ temp, Replicate(0, [nZeros, dims[1:*]]) ] $
				ELSE temp = [ temp, Replicate(0, nZeros) ]
		ENDIF
		
		;Increase number of frequencies
;		IF N_new GT nfft THEN freqs = [ [freqs], [FltArr(nMax,N_new-nfft)] ]
		
	;-----------------------------------------------------
	; Determine Frequencies \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		sr = 1.0 / SI
		fN = sr / 2.0
		df = 1.0 / (SI * N_new)
		f  = [0:fN:df]
		f  = N_new MOD 2 EQ 0 ? [ f, -Reverse(f[1:-2]) ] : [ f, -Reverse(f[1:-1]) ]
		
	;-----------------------------------------------------
	; Compute FFT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------

		;Compute the FFT
		fftTemp = FFT( Temporary(temp), DIMENSION=1, DOUBLE=tf_double )
		
		;Normalize the FFT to account FOR the window
;		CG = mean(taper)
;		NG = mean(taper^2)
;		CASE 1 OF
;			tf_rms:   fftTemp /= CG
;			tf_noise: ;Do nothing
;			tf_power: ;Do nothing
;			ELSE:     ;Do nothing
;		ENDCASE

	;-----------------------------------------------------
	; Outputs \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Save data
		IF nMax EQ 1 THEN BEGIN
			dft   = Temporary(fftTemp)
		ENDIF ELSE BEGIN
			CASE ndims OF
				1: dft[i,*]             = Temporary(fftTemp)
				2: dft[i,*,*]           = Temporary(fftTemp)
				3: dft[i,*,*,*]         = Temporary(fftTemp)
				4: dft[i,*,*,*,*]       = Temporary(fftTemp)
				5: dft[i,*,*,*,*,*]     = Temporary(fftTemp)
				6: dft[i,*,*,*,*,*,*]   = Temporary(fftTemp)
				7: dft[i,*,*,*,*,*,*,*] = Temporary(fftTemp)
				ELSE: Message, 'Incorrect number of dimensions: SELF.'
			ENDCASE
		ENDELSE
			
		;Save frequencies and time itnerval
		;   - Necessary only for I=0 if TF_CALC_DT is false
		;   - Keep track of all here to simplify logic
		freqs[i,*]   = f
		dt_median[i] = SI
		df_fft[i]    = 1.0 / (N_new * SI)
		
		;Record time stamp
;		t_fft[i]  = time[0] + (nfft * SI) / 2.0
		t_fft[i]  = (ttemp[-1] + ttemp[0]) / 2.0
		dt_fft[i] = nfft * SI / 2.0
	ENDFOR
	
	;Scalar values
	IF nMax EQ 1 THEN BEGIN
		freqs     = Reform(freqs)
		dt_median = dt_median[0]
		dt_fft    = dt_fft[0]
		df_fft    = df_fft[0]
		t_fft     = t_fft[0]
		dt_fft    = dt_fft[0]
	ENDIF

;-----------------------------------------------------
; Finishing Touches \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Reduce dimensions IF sampling interval never changed
	;   - fN = 1 / T = 1 / (nfft*dt)
	;   - NFFT is constant, but DT is permitted to vary
	;   - If DT changes, then we must keep track OF FN as a function of time
	;   - If not, then FN and F DO not vary with time, and we can reduce to one copy
	IF nMax GT 1 && Total( (dt_median - dt_median[0]) GT 0.1*dt_median[0] ) EQ 0 THEN BEGIN
		freqs     = Reform(freqs[0,*])
		df_fft    = df_fft[0]
		dt_median = dt_median[0]
		dt_fft    = dt_fft[0]
	ENDIF
	
	;Remove negative frequencies
	;   - TODO: What IF DIMENSION ~= 2
	IF tf_one_sided THEN BEGIN
		freqs = freqs[0:nfft/2]
		CASE ndims OF
			1: dft = dft[*,0:nfft/2]
			2: dft = dft[*,0:nfft/2,*]
			3: dft = dft[*,0:nfft/2,*,*]
			4: dft = dft[*,0:nfft/2,*,*,*]
			5: dft = dft[*,0:nfft/2,*,*,*,*]
			6: dft = dft[*,0:nfft/2,*,*,*,*,*]
			7: dft = dft[*,0:nfft/2,*,*,*,*,*,*]
			ELSE: Message, 'Incorrect number OF dimensions: SELF.'
		ENDCASE
	ENDIF

;-----------------------------------------------------
; Output Object \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;TIME
	IF nMax EQ 1 $
		THEN t_fft = MrCDF_Epoch_Encode( MrCDF_ssm2epoch( t_fft, self.oTime['DATA', 0, 'TT2000'] ) ) $
		ELSE oTime = MrTimeVar( t_fft, 'SSM', /NO_COPY, T_REF=self.oTime['DATA', 0] )

	;FREQUENCY
	IF Size(freqs, /DIMENSIONS) EQ 2 $
		THEN oFreq = MrTimeSeries( oTime, freqs, /NO_COPY ) $
		ELSE oFreq = MrVariable( freqs, /NO_COPY )

	;DFT
	IF nMax GT 1 $
		THEN oDFT = MrTimeSeries( oTime, dft ) $
		ELSE oDFT = MrVariable( dft, /NO_COPY )

;-----------------------------------------------------
; Attributes \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;DEPEND_# & DELTA_(MINUS|PLUS)
	IF nMax GT 1 THEN BEGIN
		;dt - Scalar
		IF N_Elements(dt_fft) EQ 1 THEN BEGIN
			oTime['DELTA_PLUS']  = dt_fft
			oTime['DELTA_MINUS'] = dt_fft
		
		;dt(t)
		ENDIF ELSE BEGIN
			;Create a time-dependent variation
			oDT = MrTimeSeries( oTime, dt_fft, /NO_COPY )
			oDT['TITLE'] = 'Deviation from center time.'
			oDT['UNITS'] = 'seconds'
			
			;Add as attribute to TIME
			oTime['DELTA_PLUS_VAR']  = oDT
			oTime['DELTA_MINUS_VAR'] = oDT
		ENDELSE
		
		;Set DFT dependencies
		oDFT['DEPEND_0']   = oTime
		oDFT['DEPEND_1']   = oFreq
		oDFT['PLOT_TITLE'] = string(nfft, nshift, FORMAT='(%"DFT: NFFT=%i, NSHIFT=%i")')
	
	;Single DFT
	ENDIF ELSE BEGIN
		oDFT['DEPEND_0']   = oFreq
		oDFT['PLOT_TITLE'] = string(nfft, FORMAT='(%"DFT: NFFT=%i")')
		oDFT['TIME_STAMP'] = t_fft
	ENDELSE
	oDFT['SAMPLE_RATE'] = nMax EQ 1 ? 1.0/SI : 1.0/dt_median
	oDFT['NFFT']        = N_new
	oDFT['NZEROPAD']    = nZeros
	oDFT['NSHIFT']      = nshift
	oDFT['WINDOW']      = taper[*,0]
	
	;TIME
	IF nMax GT 1 THEN BEGIN
		oTime -> SetName, 'fft_tcenter'
		oTime['TITLE'] = 'FFT Center Times'
	ENDIF
	
	;FREQUENCY
	oFreq -> SetName, 'fft_freq'
	oFreq['TITLE'] = 'Freq!C(Hz)'
	oFreq['UNITS'] = 'Hz'
	IF tf_one_sided THEN BEGIN
		oFreq['LOG']        = 1B
		oFreq['AXIS_RANGE'] = [oFreq['DATA',1], oFreq['DATA',-1]]
	ENDIF
	
	;DFT
	oDFT -> SetName, 'FFT'
	oDFT['TITLE'] = 'FFT'
	oDFT['SCALE'] = 1
	
	RETURN, oDFT
END


;+
;   Create a tapering window for spectral analysis.
;
; :Params:
;       N:          in, required, type=integer
;                   Number of points for the tapering window.
;       NAME:       in, optional, type=string, default=RECTANGULAR
;                   Name of the tapering window to use. Options are: {"Bartlett" | "Blackman" | 
;                       "Gaussian" | "Triangular" | "Rectangular" | "Hanning" | "Hamming" |
;                       "None" | ""}
;
; :Keywords:
;       WIDTH:      in, optional, type=integer, default=`N`
;                   If `NAME` is "Gaussian", then WIDTH defines the width of the window.
;
; :Returns:
;       TAPER:      out, required, type=fltarr
;                   The tapering window.
;-
FUNCTION MrTimeSeries::FFT_Window, n, name, $
WIDTH=width
	Compile_Opt idl2
	On_Error, 2
	
	IF N_Elements(name)  EQ 0 THEN name = 'RECTANGULAR'
	IF N_Elements(width) EQ 0 THEN width = Float(n)
	IF width GT n THEN Message, 'WIDTH must be <= N.'
	
	;Name
	CASE StrUpCase(name) OF
		'BARTLETT':    taper = 1 - abs(1 - 2*FIndGen(n)/(n-1))
		'BLACKMAN': BEGIN
			xx    = FIndGen(n)/(n-1)
			taper = 0.42 - 0.50 * Cos(2*!pi*xx) + 0.08 * Cos(4*!pi*xx)
		ENDCASE
		'GAUSSIAN':    taper = Exp( -( FIndGen(n) - (n - 1)/2 )^2 / (2*width^2) )
		'TRIANGULAR':  taper = 1 - Abs(1 - 2*FIndGen(n)/(n-1))
		'SQUARE':      taper = Replicate(1.0, n)
		'RECTANGULAR': taper = Replicate(1.0, n)
		'HANNING':     taper = Hanning(n, ALPHA=0.50)
		'HAMMING':     taper = Hanning(n, ALPHA=0.54)
		'NONE':        taper = Replicate(1.0, n)
		'':            taper = Replicate(1.0, n)
		ELSE: Message, 'Window not recognized: "' + name + '".'
	ENDCASE
	
	RETURN, taper
END

;+
;   Return the variable's data.
;
; :Params:
;       TIME:       out, optional, type=strarr
;                   Time tags, formatted as ISO-8601 strings
;       TYPE:       in, optional, type=string, default='ISO-8601'
;                   The time base in which `TIME` should be returned.
;
; :Keywords:
;       TVAR:       in, optional, type=boolean, default=0
;                   If set, `TIME` will be returned as a MrTimeVar object
;                       instead of the individual time tags.
;-
function MrTimeSeries::GetData, time, type, $
TVAR=tvar
	compile_opt idl2
	on_error, 2
	
	;Return time
	if arg_present(time) then begin
		if keyword_set(tvar) $
			then time = self.oTime $
			else time = self.oTime -> GetData(type)
	endif
	
	;Return data
	return, *self.data
end


;+
;   Return the size of the time-dependent dimension.
;
; :Returns:
;       NTIMES:     in, optional, type=intarr
;                   Number of elements along the time-dependent dimension.
;-
function MrTimeSeries::GetNPts
	compile_opt idl2
	on_error, 2
	
	;Number of records
	sz = size(*self.data, /DIMENSIONS)
	return, sz[0]
end


;+
;   Compute the inner product of two time series tensors of rank > 1.
;   This method generalizes the # operator for tensors of rank > 2.
;
; :Params:
;       B:          in, required, type=MrTimeSeries objref
;                   A tensor with > 1 dimension.
;
; :Returns:
;       RESULT:     out, required, type=MrTimeSeries objref
;                   Result of the inner product.
;-
function MrTimeSeries::InnerProduct, B, $
NAME=name
	compile_opt strictarr
	on_error, 2
	
	;B must be a time series
	if ~isa(B, 'MrTimeSeries') $
		then message, 'B must be a MrTimeSeries object.'

;-------------------------------------------------------
; Dimensions ///////////////////////////////////////////
;-------------------------------------------------------
	
	;Number of points
	nPtsA = self -> GetNPts()
	nPtsB = B    -> GetNPts()
	if nPtsA lt nPtsB then begin
		!Null = A -> GetData( oTime, !Null, /TVAR )
		nPts  = nPtsA
	endif else begin
		!Null = B -> GetData( oTime, !Null, /TVAR )
		nPts  = nPtsB
	endelse

	;Dimensionality of inputs
	dimsA  = size(self, /DIMENSIONS)
	dimsB  = size(B, /DIMENSIONS)
	nDimsA = size(self, /N_DIMENSIONS)
	nDimsB = size(B, /N_DIMENSIONS)
	nDims  = nDimsA < nDimsB
	nDims  = (nDimsA > nDimsB) - 1
	if nDimsA eq 1 then message, 'The implicit array must have more than one dimension.'
	if nDimsB eq 1 then message, 'B must have more than one dimension.'
	
	;Output dimensions
	if nDimsA gt nDimsB $
		then outDims = dimsA[0:-2] $
		else outDims = nDimsB eq 2 ? dimsB[0] : [ dimsB[0], dimsB[2:*] ]
	
	;Temporary dimensions
	;   - OUTDIMS, but with non-time-dependent dimensions contracted
	;   - The # operator accepts only two dimensions
	tempDims = nDims eq 1 ? outDims[0] : fix( [outDims[0], product(outDims[1:*])], TYPE=13 )

;-------------------------------------------------------
; Reform Inputs ////////////////////////////////////////
;-------------------------------------------------------
	
	;A is currently D1xD2xD3xD4x ... xDN
	;   - Turn A into D1xDixDN matrix
	if nDimsA le 2 $
		then _A = self['DATA'] $
		else _A = reform( self['DATA'], [ dimsA[0], product(dimsA[1:-2]), dimsA[-1] ] )
	
	;B is currently D1xD2xD3xD4x ... xDN
	;   - Turn B into D1xD2xDi matrix
	if nDimsB le 2 $
		then _B = B['DATA'] $
		else _B = reform( B['DATA'], [ dimsB[0], dimsB[1], product(dimsB[2:*]) ] )

;-------------------------------------------------------
; Inner Product ////////////////////////////////////////
;-------------------------------------------------------
	
	;Allocate memory
	temp = make_array( tempDims, TYPE=size(self[[0]]*B[[0]], /TYPE) )

	;Step through each
	;   - If _A is TxN, we want to keep each time element 1xN so that N is the inner dimension
	;   - If _A is TxNxM, we want each time element to be NxM
	if nDimsA gt 2 $
		then for i = 0, nPts - 1 do temp[i,*] = reform(_A[i,*,*]) # reform(_B[i,*,*]) $
		else for i = 0, nPts - 1 do temp[i,*] = _A[i,*,*]         # reform(_B[i,*,*])
	_A = !Null
	_B = !Null

;-------------------------------------------------------
; Finish Up ////////////////////////////////////////////
;-------------------------------------------------------

	;Reform back to the original shape, minus the contracted dimension.
	temp = reform(temp, outDims, /OVERWRITE)
	
	;Create time series variable
	result = MrTimeSeries( oTime, temp, NAME=name, /NO_COPY )
	return, result
end


;+
;   Perform interpolation on regularly or irregularly vectors.
;
;   Calling Sequence:
;       varOut = MrVar1 -> Interpol(MrVar2)
;       varOut = MrVar1 -> Interpol(X, Xout)
;
; :Params:
;       XOUT:           in, required, type=Numeric array
;                       The new time values.
;       T_TYPE:         in, required, type=string, default='ISO-8601'
;                       The units or time-base of `XOUT`. See MrTimeVar for details.
;
; :Keywords:
;       LSQUADRATIC:    in, optional, type=boolean, default=0
;                       If set, use a least squares quadratic fit to a 4-point neighborhood.
;       NAN:            in, optional, type=boolean, default=0
;                       If set, ignore NaN values.
;       QUADRATIC:      in, optional, type=boolean, default=0
;                       If set, use a quadritic fit to a 3-point neighborhood.
;       SPLINE:         in, optional, type=boolean, default=0
;                       If set, use a cubic spline in a 3-point neighborhood.
;
; :Returns:
;       VAROUT:         out, required, type=object
;                       A MrTimeSeries object containing the interpolated data.
;-
function MrTimeSeries::Interpol, Xout, t_type, $
CACHE=cache, $
NAME=name, $
LSQUADRATIC=lsquadratic, $
NAN=nan, $
QUADRATIC=quadratic, $
SPLINE=spline
	compile_opt idl2
	on_error, 2
	
	if n_elements(name) eq 0 then name = 'Interpol(' + self.name + ')'
	
;-------------------------------------------------------
; MrVariable Object ////////////////////////////////////
;-------------------------------------------------------
	if IsA(Xout, 'MrVariable') then begin
		;Get the time variables
		oT = self.oTime
		if obj_isa(Xout, 'MrTimeVar') $
			then oTout = Xout $
			else oTout = MrVar_Get(Xout['DEPEND_0'])
		
		;
		; Check if the variables are the same
		;
		IF oT -> IsIdentical(oTout) THEN BEGIN
			varOut = self -> Copy()
			varOut -> SetName, name
			IF Keyword_Set(cache) THEN varOut -> Cache
			RETURN, varOut
		ENDIF
		
		;Reference time
		!Null = min( [ oT['DATA', 0, 'JULDAY'], oTout['DATA', 0, 'JULDAY'] ], iMin )
		if iMin eq 0 $
			then t_ref = oT['DATA', 0] $
			else t_ref = oTout['DATA', 0]
		
		;Convert to seconds since midnight
		t     = oT    -> GetData('SSM', t_ref)
		t_out = oTout -> GetData('SSM', t_ref)

;-------------------------------------------------------
; Normal Arrays ////////////////////////////////////////
;-------------------------------------------------------
	endif else begin
		;Convert time to
		t     = self.oTime -> fromISO(t_type)
		t_out = Xout
	endelse

;-------------------------------------------------------
; Interpolate Scalar ///////////////////////////////////
;-------------------------------------------------------
	nDims = size(self, /N_DIMENSIONS)
	dims  = size(self, /DIMENSIONS)
	
	if nDims eq 1 then begin
		y_out = Interpol( *self.data, temporary(t), temporary(t_out), $
		                  LSQUADRATIC = lsquadratic, $
		                  NAN         = nan, $
		                  QUADRATIC   = quadratic, $
		                  SPLINE      = spline)

;-------------------------------------------------------
; Interpolate Array ////////////////////////////////////
;-------------------------------------------------------
	endif else begin
		;Reform Y so that
		;   [ time, everything else ]
		nOut    = n_elements(t_out)
		allDims = product( dims[1:*] )
		y       = reform( *self.data, [dims[0], allDims] )
		y_out   = make_array( nOut, allDims, TYPE=Size(*self.data, /TYPE) )
		
		;Interpolate all
		for i = 0, allDims - 1 do begin
			y_out[0,i] = Interpol( y[*,i], t, t_out, $
			                       LSQUADRATIC = lsquadratic, $
			                       NAN         = nan, $
			                       QUADRATIC   = quadratic, $
			                       SPLINE      = spline)
		endfor
		
		;Free up memory
		Y = !Null
		t = !Null
		
		;Reform back
		y_out = reform(y_out, [nOut, dims[1:*]], /OVERWRITE)
	endelse

;-------------------------------------------------------
; Create New Variable //////////////////////////////////
;-------------------------------------------------------
	
	;Create a new variable with the same attributes
	varOut = self -> Copy(name, CACHE=cache)
	
	;Set the data
	if obj_valid(oTout) $
		then varOut -> SetData, oTout, temporary(y_out), T_TYPE=t_type $
		else varOut -> SetData, t_out, y_out, T_TYPE=t_type, /NO_COPY
	
	;Return the new variable
	return, varOut
end


;+
;   Create a new copy. Properties of the parent object are passed on to the
;   child object.
;-
function MrTimeSeries::New, time, array, $
_REF_EXTRA=extra
	compile_opt idl2
	on_error, 2

	;Return a copy of the array
	;   - Use Obj_Class() so that subclasses can inherit the method.
	return, obj_new(obj_class(self), time, array, _STRICT_EXTRA=extra)
end 


;+
;   Compute the outer product of two time series tensors of rank > 1.
;   This method generalizes the ## operator for tensors of rank > 2.
;
; :Params:
;       B:          in, required, type=MrTimeSeries objref
;                   A tensor with > 1 dimension.
;
; :Returns:
;       RESULT:     out, required, type=MrTimeSeries objref
;                   Result of the outer product.
;-
function MrTimeSeries::OuterProduct, B, $
NAME=name
	compile_opt strictarr
	on_error, 2
	
	;B must be a time series
	if ~isa(B, 'MrTimeSeries') $
		then message, 'B must be a MrTimeSeries object.'

;-------------------------------------------------------
; Dimensions ///////////////////////////////////////////
;-------------------------------------------------------
	;Number of points
	nPtsA = self -> GetNPts()
	nPtsB = B    -> GetNPts()
	if nPtsA lt nPtsB then begin
		!Null = self -> GetData( oTime, /TVAR )
		nPts  = nPtsA
	endif else begin
		!Null = self -> GetData( oTime, /TVAR )
		nPts  = nPtsB
	endelse

	;Dimensionality of inputs
	dimsA  = size(self, /DIMENSIONS)
	dimsB  = size(B, /DIMENSIONS)
	nDimsA = size(self, /N_DIMENSIONS)
	nDimsB = size(B, /N_DIMENSIONS)
	nDims  = (nDimsA > nDimsB) - 1
	if nDimsA eq 1 then message, 'The implicit array must have more than one dimension.'
	if nDimsB eq 1 then message, 'B must have more than one dimension.'
	
	;Output dimensions
	if nDimsA gt nDimsB $
		then outDims = nDimsA eq 2 ? dimsA[0] : [ dimsA[0], dimsA[2:*] ] $
		else outDims = dimsB[0:-2]
	
	;Temporary dimensions
	;   - OUTDIMS, but with non-time-dependent dimensions contracted
	;   - The ## operator accepts only two dimensions
	tempDims = nDims eq 1 ? outDims[0] : fix( [outDims[0], product(outDims[1:*])], TYPE=13 )

;-------------------------------------------------------
; Reform Inputs ////////////////////////////////////////
;-------------------------------------------------------
	;A is currently D1xD2xD3xD4x ... xDN
	;   - Turn A into D1xD2xDi matrix
	if nDimsA le 2 $
		then _A = self['DATA'] $
		else _A = reform( self['DATA'], [ dimsA[0], dimsA[1], product(dimsA[2:*]) ] )
	
	;B is currently D1xD2xD3xD4x ... xDN
	;   - Turn B into D1xDixDN matrix
	if nDimsB le 2 $
		then _B = B['DATA'] $
		else _B = reform( B['DATA'], [ dimsB[0], product(dimsB[1:-2]), dimsB[-1] ] )

;-------------------------------------------------------
; Outer Product ////////////////////////////////////////
;-------------------------------------------------------
	
	;Allocate memory
	temp = make_array( tempDims, TYPE=size(self['DATA',0]*B['DATA',0], /TYPE) )

	;Step through each
	;   - If _B is TxN, we want to keep each time element 1xN so that N is the outer dimension
	;   - If _B is TxNxM, we want B to be NxM
	if nDimsB gt 2 $
		then for i = 0, nPts - 1 do temp[i,*] = reform(_A[i,*,*]) ## reform(_B[i,*,*]) $
		else for i = 0, nPts - 1 do temp[i,*] = reform(_A[i,*,*]) ## _B[i,*,*]
	_A = !Null
	_B = !Null

;-------------------------------------------------------
; Finish Up ////////////////////////////////////////////
;-------------------------------------------------------

	;Reform back to the original shape, minus the contracted dimension.
	temp = reform(temp, outDims, /OVERWRITE)
	
	;Create time series variable
	result = MrTimeSeries( oTime, temp, NAME=name, /NO_COPY )
	return, result
end


;+
;   Compute the PSD
;
; :Params:
;       AMPLITUDE:      out, optional, type=MrTimeVar/MrVariable object
;                       The amplitude of the signal.
;       PHASE:          out, optional, type=MrTimeVar/MrVariable object
;                       The phase of the signal.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, `RESULT` will be added to the variable cache.
;       NDETREND:       in, optional, type=integer, default=0
;                       The number of points in a boxcar average. The average is then 
;                           removed from the data to "detrend" it.
;       NAME:           in, optional, type=string, default='PSD('+<name>+')'
;                       Name to be given to `RESULT`.
;       NOISE:          in, optional, type=boolean, default=0
;                       If set, the PSD will be normalized to the noise value. Cannot be
;                           used with `POWER` or `RMS`.
;       POWER:          in, optional, type=boolean, default=0
;                       If set, the PSD will be normalized to the power. Cannot be
;                           used with `NOISE` or `RMS`.
;       RMS:            in, optional, type=boolean, default=0
;                       If set, the PSD will be normalized to the rms value. Cannot be
;                           used with `NOISE` or `POWER`.
;       WINDOW:         in, optional, type=fltarr, default=fltarr()+1
;                       A tapering window to be applied to the data.
;
; :Returns:
;       RESULT:         out, required, type=MrVariable object
;                       The results of the fast fourier transformation. If more than
;                           one FFT is possible, a MrTimeSeries object is returned.
;                           Otherwise, a MrVariable object is returned.
;-
function MrTimeSeries::PSD_Old, $
CACHE=cache, $
NDETREND=nDetrend, $
NAME=name, $
NOISE=noise, $
POWER=power, $
RMS=rms, $
WINDOW=win
	compile_opt idl2
	on_error, 2
	
	tf_rms   = keyword_set(rms)
	tf_noise = keyword_set(noise)
	tf_power = keyword_set(power)
	nfft = self -> GetNPts()
	if n_elements(nDetrend) eq 0 then nDetrend = 0
	if n_elements(name)     eq 0 then name     = 'PSD(' + self.name + ')'
	
	;Gain from window
;	NG = mean(win^2)
;	CG = mean(win)

;-----------------------------------------------------
; Detrend Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	if nDetrend gt 0 then begin
		case size(self, /N_DIMENSIONS) of
			1: oVar = self - smooth(self['DATA'], nDetrend, /EDGE_TRUNCATE)
			2: oVar = self - smooth(self['DATA'], [nDetrend,1], /EDGE_TRUNCATE)
			3: oVar = self - smooth(self['DATA'], [nDetrend,1,1], /EDGE_TRUNCATE)
			4: oVar = self - smooth(self['DATA'], [nDetrend,1,1,1], /EDGE_TRUNCATE)
			5: oVar = self - smooth(self['DATA'], [nDetrend,1,1,1,1], /EDGE_TRUNCATE)
			6: oVar = self - smooth(self['DATA'], [nDetrend,1,1,1,1,1], /EDGE_TRUNCATE)
			7: oVar = self - smooth(self['DATA'], [nDetrend,1,1,1,1,1,1], /EDGE_TRUNCATE)
			8: oVar = self - smooth(self['DATA'], [nDetrend,1,1,1,1,1,1,1], /EDGE_TRUNCATE)
			else: message, 'Internal array has invalid number of dimensions.'
		endcase
	endif else begin
		oVar = self
	endelse

;-----------------------------------------------------
; Compute FFT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oFFT = oVar -> FFT_Old( nfft, nshift, $
	                        NOISE  = noise, $
	                        POWER  = power, $
	                        RMS    = rms, $
	                        WINDOW = win )
	if ~obj_valid(oFFT) then return, !Null

;-----------------------------------------------------
; One-Sided PSD \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	dims  = size(oFFT, /DIMENSIONS)
	ndims = size(oFFT, /N_DIMENSIONS)

	;One-sided
	case ndims of
		1: pwr = 2.0 * oFFT['DATA',0:nfft/2]               * conj( oFFT['DATA',0:nfft/2] )
		2: pwr = 2.0 * oFFT['DATA',0:nfft/2,*]             * conj( oFFT['DATA',0:nfft/2,*] )
		3: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*]           * conj( oFFT['DATA',0:nfft/2,*,*] )
		4: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*,*]         * conj( oFFT['DATA',0:nfft/2,*,*,*] )
		5: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*,*,*]       * conj( oFFT['DATA',0:nfft/2,*,*,*,*] )
		6: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*,*,*,*]     * conj( oFFT['DATA',0:nfft/2,*,*,*,*,*] )
		7: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*,*,*,*,*]   * conj( oFFT['DATA',0:nfft/2,*,*,*,*,*,*] )
		8: begin
			fft_data = oFFT -> GetData()
			pwr      = 2.0 * fft_data[0:nfft/2,*,*,*,*,*,*,*] * conj( fft_data[0:nfft/2,*,*,*,*,*,*,*] )
		endcase
		else: message, 'FFT must have <= 8 dimensions.'
	endcase

;-----------------------------------------------------
; Set Properties |\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;PWR attributes
	oPWR = MrVariable( real_part(pwr), CACHE=cache, NAME=name, /NO_COPY )
	
	;Frequency attributes
	oFreq = oFFT['DEPEND_0']
	oFreq['AXIS_RANGE'] = [oFreq['DATA',1], oFreq['DATA',nfft/2]]
	oFreq['LOG']        = 1B
	oFreq['TITLE']      = 'Freq (' + oFreq['UNITS'] + ')'
	
	oPWR['DEPEND_0'] = oFreq[0:nfft/2]
	oPWR['TITLE']    = 'PSD'
	oPWR['LOG']      = 1

	;Normalize
;	case 1 of
;		tf_rms:   pwr = NG * pwr / CG^2
;		tf_noise: pwr = pwr / (NG * df)
;		tf_power: pwr = pwr / NG
;		else:     ;Do nothing
;	endcase
	
	;Normalize
;	case 1 of
;		tf_rms:   pwr = NG * pwr / CG^2
;		tf_noise: pwr = NG * df * pwr / CG^2
;		tf_power: pwr = NG * pwr / CG^2
;		else:   ;Nothing
;	endcase
	
	return, oPWR
end


;+
;   Compute the PSD
;
; :Params:
;       NFFT:           in, optional, type=int, default=nPts
;                       Number of points in each FFT. If NFFT is shorter than the
;                           signal length defined by `RANGE`, then the PSD of each NFFT
;                           segment that fits within `RANGE` will be averaged together.
;       NSHIFT:         in, optional, type=integer, default=NFFT/2
;                       Number of points to shift between consecutive FFTs.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, `RESULT` will be added to the variable cache.
;       NAME:           in, optional, type=string, default='PSD('+<name>+')'
;                       Name to be given to `RESULT`.
;       RANGE:          in, optional, type=intarr(2)/strarr(2), default="[0,nPts-1]"
;                       A time or index range outlining a subset of the data for which
;                           the PSD is computed.
;       WINDOW:         in, optional, type=string, default='Rectangular'
;                       The name of a tapering window to be applied to the data. See
;                           the FFT_Window method for valid names.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by the ::FFT2 method is also accepted here.
;
; :Returns:
;       OPSD:           out, required, type=MrVariable object
;                       The power spectral density of the intrinsic time series signal.
;-
FUNCTION MrTimeSeries::PSD, nfft, nshift, $
CACHE=cache, $
NAME=name, $
RANGE=range, $
_REF_EXTRA=extra
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	;Range
	irange = N_Elements(range) EQ 0 ? [0, (self -> GetNPts()) - 1] : range
	IF Size(irange, /TNAME) EQ 'STRING' THEN irange = self -> Value_Locate(range)
	
	;Other defaults
	tf_rms    = Keyword_Set(rms)
	tf_noise  = Keyword_Set(noise)
	tf_power  = Keyword_Set(power)
	IF N_Elements(nfft)   EQ 0 THEN nfft   = irange[1] - irange[0] + 1
	IF N_Elements(nshift) EQ 0 THEN nshift = nfft/2
	IF N_Elements(name)   EQ 0 THEN name   = 'PSD(' + self.name + ')'
	
	;Gain from window
;	NG = mean(win^2)
;	CG = mean(win)

;-----------------------------------------------------
; Compute FFT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	oFFT = self -> FFT( nfft, nshift, $
	                    RANGE         = irange, $
	                    _STRICT_EXTRA = extra )

;-----------------------------------------------------
; Normalizations \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Zero padding
	nNew     = nfft + oFFT['NZEROPAD']
	zeroNorm = nNew^2.0 / Float(nfft)
	
	;Tapering window
	winNorm = Total(oFFT['WINDOW']^2) / nfft
	
;-----------------------------------------------------
; One-Sided PSD \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Time-Series
	tf_ts = IsA(oFFT, 'MrTimeSeries')
	
	;Time-Series
	pFFT = oFFT['PTR']
	IF tf_ts THEN BEGIN
		;If a time series, minimum # dims = 2: [time, freq]
		CASE SiZE(oFFT, /N_DIMENSIONS) OF
			2: temp = (*pFFT)[*,0:nNew/2]             * Conj( (*pFFT)[*,0:nNew/2] )
			3: temp = (*pFFT)[*,0:nNew/2,*]           * Conj( (*pFFT)[*,0:nNew/2,*] )
			4: temp = (*pFFT)[*,0:nNew/2,*,*]         * Conj( (*pFFT)[*,0:nNew/2,*,*] )
			5: temp = (*pFFT)[*,0:nNew/2,*,*,*]       * Conj( (*pFFT)[*,0:nNew/2,*,*,*] )
			6: temp = (*pFFT)[*,0:nNew/2,*,*,*,*]     * Conj( (*pFFT)[*,0:nNew/2,*,*,*,*] )
			7: temp = (*pFFT)[*,0:nNew/2,*,*,*,*,*]   * Conj( (*pFFT)[*,0:nNew/2,*,*,*,*,*] )
			8: temp = (*pFFT)[*,0:nNew/2,*,*,*,*,*,*] * Conj( (*pFFT)[*,0:nNew/2,*,*,*,*,*,*] )
;			8: BEGIN
;				fft_data = oFFT -> GetData()
;				temp     = fft_data[*,0:nNew/2,*,*,*,*,*,*] * Conj( fft_data[*,0:nNew/2,*,*,*,*,*,*] )
;				fft_data = !Null
;			ENDCASE
			ELSE: Message, 'FFT must have <= 8 dimensions.'
		ENDCASE
	
	;Single PSD
	ENDIF ELSE BEGIN
		;One-sided
		CASE SiZE(oFFT, /N_DIMENSIONS) OF
			1: temp = (*pFFT)[0:nNew/2]               * Conj( (*pFFT)[0:nNew/2] )
			2: temp = (*pFFT)[0:nNew/2,*]             * Conj( (*pFFT)[0:nNew/2,*] )
			3: temp = (*pFFT)[0:nNew/2,*,*]           * Conj( (*pFFT)[0:nNew/2,*,*] )
			4: temp = (*pFFT)[0:nNew/2,*,*,*]         * Conj( (*pFFT)[0:nNew/2,*,*,*] )
			5: temp = (*pFFT)[0:nNew/2,*,*,*,*]       * Conj( (*pFFT)[0:nNew/2,*,*,*,*] )
			6: temp = (*pFFT)[0:nNew/2,*,*,*,*,*]     * Conj( (*pFFT)[0:nNew/2,*,*,*,*,*] )
			7: temp = (*pFFT)[0:nNew/2,*,*,*,*,*,*]   * Conj( (*pFFT)[0:nNew/2,*,*,*,*,*,*] )
			8: temp = (*pFFT)[0:nNew/2,*,*,*,*,*,*,*] * Conj( (*pFFT)[0:nNew/2,*,*,*,*,*,*,*] )
;			8: BEGIN
;				fft_data = oFFT -> GetData()
;				temp     = fft_data[0:nNew/2,*,*,*,*,*,*,*] * Conj( fft_data[0:nNew/2,*,*,*,*,*,*,*] )
;				fft_data = !Null
;			ENDCASE
			ELSE: Message, 'FFT must have <= 8 dimensions.'
		ENDCASE
	ENDELSE
	
	;Normalize
	psd = 2.0 * zeroNorm / winNorm / oFFT['SAMPLE_RATE'] * Temporary(temp)
	
	;Average over time
	IF tf_ts THEN psd = Mean(psd, DIMENSION=1)

;-----------------------------------------------------
; Set Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;PWR attributes
	oPSD = MrVariable( Real_Part(psd), CACHE=cache, NAME=name, /NO_COPY )
	oFFT -> CopyAttrTo, oPSD, ['NFFT', 'NSHIFT', 'NZEROPAD', 'SAMPLE_RATE', 'WINDOW']
	IF ~tf_ts THEN oFFT -> CopyAttrTo, oPSD, 'TIME_STAMP'
	
	;Frequency attributes
	oFreq = tf_ts ? oFFT['DEPEND_1'] : oFFT['DEPEND_0']
	oFreq['AXIS_RANGE'] = [oFreq['DATA',1], oFreq['DATA',nNew/2]]
	oFreq['LOG']        = 1B
	oFreq['TITLE']      = 'Freq (' + oFreq['UNITS'] + ')'
	
	oPSD['DEPEND_0']  = oFreq[0:nNew/2]
	oPSD['DIMENSION'] = 1
	oPSD['TITLE']     = 'PSD'
	oPSD['LOG']       = 1

	;Normalize
;	CASE 1 OF
;		tf_rms:   pwr = NG * pwr / CG^2
;		tf_noise: pwr = pwr / (NG * df)
;		tf_power: pwr = pwr / NG
;		ELSE:     ;Do nothing
;	ENDCASE
	
	;Normalize
;	CASE 1 OF
;		tf_rms:   pwr = NG * pwr / CG^2
;		tf_noise: pwr = NG * df * pwr / CG^2
;		tf_power: pwr = NG * pwr / CG^2
;		ELSE:   ;Nothing
;	ENDCASE
	
	RETURN, oPSD
END


;+
;   Determine if time stamps are identical. Time stamps are identical if the
;   MrTimeVar objects have the same heap identifier.
;
; :Params:
;       OVAR:           in, required, type=objref
;                       Time stamps to which the internal time series should be identical.
;                           May be a MrTimeVar or MrTimeSeries class or subclass.
;
; :Returns:
;       TF_SAME:        out, required, type=boolean
;                       If true, time stamps are identical. False otherwise.
;-
function MrTimeSeries::IsTimeIdentical, oVar
	compile_opt idl2
	on_error, 2
	
	tf_same = 0B
	if obj_isa(oVar, 'MrTimeVar') then begin
		tf_same = obj_valid(self.oTime, /GET_HEAP_IDENTIFIER) EQ Obj_Valid(oVar, /GET_HEAP_IDENTIFIER)
	endif else if obj_isa(oVar, 'MrTimeSeries') then begin
		tf_same = obj_valid(self.oTime, /GET_HEAP_IDENTIFIER) EQ Obj_Valid(oVar['TIMEVAR'], /GET_HEAP_IDENTIFIER)
	endif else begin
		message, 'OVAR must be a MrTimeVar or MrTimeSeries object.'
	endelse
	
	return, tf_same
end


;+
;   Set the value of a single attribute.
;
; :Params:
;       ATTRNAME:       in, required, type=string
;                       The name of the attribute for which the value is to be changed.
;       ATTRVALUE:      in, optional, type=any accepted by ::AttrValue_IsValid
;                       The value of the attribute.
;
; :Keyword:
;       CREATE:         in, optional, type=boolean, default=0
;                       If set, the attribute will be created if it does not already exist.
;       OVERWRITE:      in, optional, type=boolean, default=1
;                       If set to zero, the attribute value will not be set if `ATTRNAME`
;                           already exists as an attribute.
;-
PRO MrTimeSeries::SetAttributeValue, attrName, attrValue, $
CREATE=create, $
OVERWRITE=overwrite
	Compile_Opt idl2
	On_Error, 2
	
	;Check inputs
	if ~IsA(attrName, 'STRING', /SCALAR) $
		then message, 'ATTRNAME must be a scalar string.'
	
	;Forbidden attributes
	;   - These are special inputs to ::_OverloadBracketsRightSide
	;   - These are in addition to the forbidden attributes of MrVariable
	IF MrIsMember(['TIME', 'TIMEVAR'], attrName) $
		THEN Message, '"' + attrName + '" cannot be an attribute name.'
	
	;Call superclass
	self -> MrVariable::SetAttributeValue, attrName, attrValue, $
	                                       CREATE    = create, $
	                                       OVERWRITE = overwrite
END


;+
;   Set the array.
;
;   CALLING SEQUENCE:
;       oTS -> SetData, data
;       oTS -> SetData, time, data
;
; :Keywords:
;       TIME:           in, required, type=NxM array
;                       Name or reference of a MrTimeVar object, or an array
;                           of time stamps. If a name is provided, the assiciated
;                           variable must exist in the variable cache.
;       DATA:           in, required, type=NxM array
;                       Name or reference of a MrVariable object, or the dependent
;                           variable data. If a name is given, the associated variable
;                           must exist in the variable cache. 
;
; :Keywords:
;       DIMENSION:      in, optional, type=integer
;                       The time-dependent dimension of `DATA` (1-based). If not
;                           provided, the dimension of `DATA` that is equal in size to
;                           `TIME` is chosen as the default. If zero or multiple
;                           dimensions match in this way, an error will occur.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       If set `DATA` will be copied directly into the object
;                           and will be left undefined (a MrTimeSeries object will not
;                           be destroyed, but its array will be empty).
;       T_NAME:         in, optional, type=integer
;                       Name to be given to the MrTimeVar object. Ignored unless `TIME`
;                           is an array of time stamps.
;       T_REF:              in, optional, type=string, default=''
;                           A reference time from which `TIME` is measured. Applicable
;                               only to certain values of `TYPE`. Must be formatted as
;                               an ISO-8601 string.
;       T_TYPE:         in, optional, type=integer
;                       If `TIME` is an array of time stamps, use this keyword to indicate
;                           the format or time-basis. See MrTimeVar for more details.
;-
PRO MrTimeSeries::SetData, x1, x2, $
DIMENSION=dimension, $
NO_COPY=no_copy, $
T_NAME=t_name, $
T_REF=t_ref, $
T_TYPE=t_type
	Compile_Opt idl2
	On_Error, 2
	
	;oVar -> SetDAta, x1, x2
	IF N_Elements(x2) GT 0 THEN BEGIN
		time   = x1
		ts_var = x2
	
	;oVar -> SetData, x1
	ENDIF ELSE BEGIN
		ts_var = x1
	ENDELSE
	
;-----------------------------------------------------
; Check Time \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Use existing time
	IF N_Elements(time) EQ 0 THEN BEGIN
		IF Obj_Valid(self.oTime) $
			THEN oTime = self.oTime $
			ELSE Message, 'No time variable exists. Please provide.'
	
	;Check given time
	ENDIF ELSE BEGIN
		;Object
		IF Size(time, /TNAME) EQ 'OBJREF' THEN BEGIN
			;MrTimeVar
			IF Obj_IsA(time, 'MrTimeVar') $
				THEN oTime = time $
				ELSE Message, 'X1 must be a MrTimeVar object.'
	
		;Name
		ENDIF ELSE IF MrIsA(time, 'STRING', /SCALAR) && MrVar_IsCached(time) THEN BEGIN
			;Cached?
			oTime = MrVar_Get(time)
			IF ~Obj_IsA(oTime, 'MrTimeVar') $
				THEN Message, 'X1 must be the name of a cached MrTimeVar object.'
	
		;Data
		ENDIF ELSE BEGIN
			oTime = MrTimeVar(time, t_type, NAME=t_name, NO_COPY=no_copy, T_REF=t_ref)
			IF ~Obj_Valid(oTime) THEN Message, 'Could not create X1 variable.'
		ENDELSE
	ENDELSE
	
;-----------------------------------------------------
; Check Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Object
	IF Size(ts_var, /TNAME) EQ 'OBJREF' THEN BEGIN
		;MrVariable
		IF Obj_IsA(ts_var, 'MrVariable') $
			THEN data = ts_var -> GetData() $
			ELSE Message, 'Only "MrVariable" objects can be given.'
	
	;Name
	ENDIF ELSE IF MrIsA(ts_var, 'STRING', /SCALAR) THEN BEGIN
		;Cached?
		IF MrVar_IsCached(ts_var) THEN BEGIN
			oData = MrVar_Get(ts_var)
			IF Obj_IsA(oData, 'MrVariable') $
				THEN data = oData -> GetData() $
				ELSE Message, 'Only names of MrVariable objects may be given.'
		ENDIF ELSE BEGIN
			Message, 'Only names of cached variables may be provided.'
		ENDELSE
	
	;Data
	ENDIF ELSE BEGIN
		data = Temporary(ts_var)
	ENDELSE
	
;-----------------------------------------------------
; Put Time First \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Dimension sizes
	dims  = size(data, /DIMENSIONS)
	nDims = N_Elements(dims)
	nPts  = N_Elements(oTime)
	
	;Find the time-dependent dimension
	IF N_Elements(dimension) EQ 0 THEN BEGIN
		dimension = Where(dims EQ nPts, nt) + 1
		IF nt GT 1 THEN BEGIN
			IF dims[0] EQ nPts THEN BEGIN
				dimension = 1
				MrPrintF, 'LogWarn', 'Multiple dims of size ' + strtrim(nPts, 2) + '. Assuming DIMENSION=1.'
			ENDIF ELSE BEGIN
				Message, 'Cannot determine time-dependent dimension.'
			ENDELSE
		ENDIF
	ENDIF ELSE BEGIN
		IF dims[dimension-1] NE nPts $
			THEN Message, 'The size of TIME and DATA along DIMENSION must be equal.'
	ENDELSE
	
	;Put the time-dependent dimension first
	IF nDims GT 1 && dimension NE 1 THEN BEGIN
		IF dimension EQ nDims $
			THEN p = [dimension-1, BindGen(nDims-1)] $
			ELSE p = [dimension-1, BindGen(dimension-1), BindGen(nDims-dimension)+dimension]
		
		;Transpose the data
		data = Transpose(Temporary(data), Temporary(p))
	ENDIF
	
;-----------------------------------------------------
; Store Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Store the data as properties
	;   - Let automatic garbage collection take care of the OTIME object
	;     and DATA pointer (IDL v8.0+)
	;   - The OTIME property may be shared by another MrTimeSeries object
	;   - The DATA pointer may be saved by a superclass in case additional
	;     requirements on the dimensionality of the data is required.
	;       * e.g. pData = self.data
	;              self  -> MrTimeSeries::SetData, time, newData
	;              IF [dimension check fails] then self.data = pData
	self.oTime = oTime
	self.data  = Ptr_New(data, /NO_COPY)
	
	;Copy attributes
	IF IsA(ts_var, 'MrVariable') $
		THEN ts_var -> CopyAttrTo, self
	
	;Set attributes
	self -> SetAttrValue, 'DEPEND_0',  self.oTime, /CREATE
	self -> SetAttrValue, 'DIMENSION', 1,          /CREATE
	
	;Clear data
	;   - Do not erase a MrTimeVar object. It may be referenced
	;     as a DEPEND_0 attribute for other objects.
	IF Keyword_Set(no_copy) THEN BEGIN
		IF N_Elements(x2) EQ 0 THEN BEGIN
			x1 = !Null
		ENDIF ELSE BEGIN
			IF Size(x1, /TNAME) NE 'OBJREF' THEN x1 = !Null
			x2 = !Null
		ENDELSE
	ENDIF
END


;+
;   Create a spectrogram
;
; :Params:
;       NFFT:           in, optional, type=integer, default=all
;                       Number of points per FFT.
;       NSHIFT:         in, optional, type=integer, default=`NFFT`/2.0
;                       Number of points to shift between FFTs.
;
; :Params:
;       AMPLITUDE:      out, optional, type=MrTimeVar/MrVariable object
;                       The amplitude of the signal.
;       PHASE:          out, optional, type=MrTimeVar/MrVariable object
;                       The phase of the signal.
;
; :Returns:
;       RESULT:         out, required, type=MrTimeSeries or MrVariable object
;                       The results of the fast fourier transformation. If more than
;                           one FFT is possible, a MrTimeSeries object is returned.
;                           Otherwise, a MrVariable object is returned.
;-
function MrTimeSeries::Spectrogram_Old, nfft, nshift, $
CACHE=cache, $
NAME=name, $
NOISE=noise, $
POWER=power, $
RMS=rms, $
WINDOW=win
	compile_opt idl2
	on_error, 2
	
	tf_rms   = keyword_set(rms)
	tf_noise = keyword_set(noise)
	tf_power = keyword_set(power)
	if n_elements(nfft) eq 0 then nfft = self -> GetNPts()
	if n_elements(name) eq 0 then name = 'Spectrogram(' + self.name + ')'
	
	if nfft gt self -> GetNPts() then message, 'NFFT > number of records.'
	
	;Gain from window
;	NG = mean(win^2)
;	CG = mean(win)

;-----------------------------------------------------
; Compute FFT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oFFT = self -> FFT_Old( nfft, nshift, $
	                        NOISE  = noise, $
	                        POWER  = power, $
	                        RMS    = rms, $
	                        WINDOW = win )
	if ~obj_valid(oFFT) then return, !Null

;-----------------------------------------------------
; One-Sided Spectrogram \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	dims  = size(oFFT, /DIMENSIONS)
	ndims = size(oFFT, /N_DIMENSIONS)
	
	;Compute power spectral density
	if obj_isa(oFFT, 'MrTimeSeries') then begin
		;Positive frequencies
		oFreq = oFFT['DEPEND_1']
		if size(oFreq, /N_DIMENSIONS) eq 1 $
			then oFreq -> SetData, oFreq[0:nfft/2] $
			else oFreq -> SetData, oFreq[*,0:nfft/2]
		
		;One-side
		case ndims of
			2: pwr = 2.0 * oFFT['DATA',*,0:nfft/2]             * conj( oFFT['DATA',*,0:nfft/2] )
			3: pwr = 2.0 * oFFT['DATA',*,0:nfft/2,*]           * conj( oFFT['DATA',*,0:nfft/2,*] )
			4: pwr = 2.0 * oFFT['DATA',*,0:nfft/2,*,*]         * conj( oFFT['DATA',*,0:nfft/2,*,*] )
			5: pwr = 2.0 * oFFT['DATA',*,0:nfft/2,*,*,*]       * conj( oFFT['DATA',*,0:nfft/2,*,*,*] )
			6: pwr = 2.0 * oFFT['DATA',*,0:nfft/2,*,*,*,*]     * conj( oFFT['DATA',*,0:nfft/2,*,*,*,*] )
			7: pwr = 2.0 * oFFT['DATA',*,0:nfft/2,*,*,*,*,*]   * conj( oFFT['DATA',*,0:nfft/2,*,*,*,*,*] )
			8: begin
				data = oFFT -> GetData()
				pwr  = 2.0 * data[*,0:nfft/2,*,*,*,*,*,*] * conj( data[*,0:nfft/2,*,*,*,*,*,*] )
				data = !Null
			endcase
			else: message, 'FFT must have <= 8 dimensions.'
		endcase
	
	;Single PSD
	endif else begin
		;Positive frequencies
		oFreq = oFFT['DEPEND_0']
		oFreq -> SetData, oFreq[0:nfft/2]
		
		;One-sided
		case ndims of
			1: pwr = 2.0 * oFFT['DATA',0:nfft/2]               * conj( oFFT['DATA',0:nfft/2] )
			2: pwr = 2.0 * oFFT['DATA',0:nfft/2,*]             * conj( oFFT['DATA',0:nfft/2,*] )
			3: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*]           * conj( oFFT['DATA',0:nfft/2,*,*] )
			4: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*,*]         * conj( oFFT['DATA',0:nfft/2,*,*,*] )
			5: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*,*,*]       * conj( oFFT['DATA',0:nfft/2,*,*,*,*] )
			6: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*,*,*,*]     * conj( oFFT['DATA',0:nfft/2,*,*,*,*,*] )
			7: pwr = 2.0 * oFFT['DATA',0:nfft/2,*,*,*,*,*,*]   * conj( oFFT['DATA',0:nfft/2,*,*,*,*,*,*] )
			8: begin
				data = oFFT -> GetData()
				pwr  = 2.0 * data[0:nfft/2,*,*,*,*,*,*,*] * conj( data[0:nfft/2,*,*,*,*,*,*,*] )
				data = !Null
			endcase
			else: message, 'FFT must have <= 8 dimensions.'
		endcase
	endelse
	
	;PWR attributes
	oPWR  = oFFT -> Copy()
	oPWR -> SetName, name
	oPWR -> SetData, real_part(pwr), /NO_COPY
	oPWR['TITLE'] = 'PSD'
	oPWR['LOG']   = 1
	
	;Frequency attributes
	oFreq['LOG']        = 1B
	oFreq['AXIS_RANGE'] = [oFreq['DATA',1], oFreq['DATA',-1]]

	;Normalize
;	case 1 of
;		tf_rms:   pwr = NG * pwr / CG^2
;		tf_noise: pwr = pwr / (NG * df)
;		tf_power: pwr = pwr / NG
;		else:     ;Do nothing
;	endcase
	
	;Normalize
;	case 1 of
;		tf_rms:   pwr = NG * pwr / CG^2
;		tf_noise: pwr = NG * df * pwr / CG^2
;		tf_power: pwr = NG * pwr / CG^2
;		else:   ;Nothing
;	endcase
	
	if keyword_set(cache) then oPWR -> Cache
	return, oPWR
end


;+
;   Create a spectrogram
;
; :Params:
;       TFFT:           in, optional, type=int, default=TFFT/6
;                       The number of points in each spectral estimate.
;       TSHIFT:         in, optional, type=integer, default=TFFT/2
;                       Number of points to shift between each spectral estimate.
;       NFFT:           in, optional, type=int, default=nPts
;                       Number of points in each FFT. If NFFT is shorter than the
;                           `TFFT`, then the PSD of each NFFT segment that fits within
;                           `TFFT` will be averaged together.
;       NSHIFT:         in, optional, type=integer, default=NFFT/2
;                       Number of points to shift between consecutive FFTs.
;
; :Params:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, `RESULT` will be added to the variable cache.
;       NAME:           in, optional, type=string, default='CrossSpectram('+<name>+','+<name>+')'
;                       Name to be given to `RESULT`.
;       RANGE:          in, optional, type=intarr(2)/strarr(2), default="[0,nPts-1]"
;                       A time or index range outlining a subset of the data for which
;                           the CSD is computed.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by the ::FFT method is also accepted here.
;
; :Returns:
;       RESULT:         out, required, type=objref
;                       A MrTimeSeries variable containing the Spectrogram.
;-
FUNCTION MrTimeSeries::Spectrogram, tfft, tshift, nfft, nshift, $
CACHE=cache, $
NAME=name, $
RANGE=range, $
_REF_EXTRA=extra
	Compile_Opt idl2
	On_Error, 2
	
	N = self -> GetNPts()
	IF N_Elements(tfft)   EQ 0 THEN tfft = N / 6
	IF N_Elements(nfft)   EQ 0 THEN nfft = tfft
	IF N_Elements(tshift) EQ 0 THEN tshift = tfft / 2
	IF N_Elements(nshift) EQ 0 THEN nshift = nfft / 2
	IF N_Elements(name)   EQ 0 THEN name = 'Spectrogram(' + self.name + ')'
	IF N_Elements(range)  GT 0 THEN MrPrintF, 'LogWarn', 'RANGE keyword is ignored.'
	
	;Restrictions
	IF tfft GT N THEN Message, 'NFFT > number of records.'
	
	;Gain from window
;	NG = mean(win^2)
;	CG = mean(win)

;-----------------------------------------------------
; Compute FFT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Total number of FFTs to perform
	nIter = Floor( (N - tfft) / tshift ) + 1
	
	;Initial loop conditions
	iStart = 0
	iEnd   = tfft - 1
	
	;Begin loop
	FOR i = 0, nIter - 1 DO BEGIN
		;Compute the PSD
		oPSD = self -> PSD( nfft, nshift, $
		                    RANGE         = [iStart, iEnd], $
		                    _STRICT_EXTRA = extra )
		IF ~Obj_Valid(oPSD) THEN RETURN, !Null
		
		;Allocate memory
		IF i EQ 0 THEN BEGIN
			;Create variables
			oTime = MrTimeVar( Replicate('0000-00-00T00:00:00', nIter) )
			oSpec = MrTimeSeries( oTime, FltArr( [nIter, Size(oPSD, /DIMENSIONS)] ), CACHE=cache, NAME=name )
			oSpec['DEPEND_1'] = oPSD['DEPEND_0']
			nDims = Size(oPSD, /N_DIMENSIONS)
		ENDIF
		
		;Store the data
		MrTimeParser, oPSD['TIME_STAMP'], '%d-%c-%Y %H:%m:%S%f', '%Y-%M-%dT%H:%m:%S%f%z', tstamp
		oTime[i] = tstamp
		CASE nDims OF
			1: oSpec[i,*]             = oPSD['DATA']
			2: oSpec[i,*,*]           = oPSD['DATA']
			3: oSpec[i,*,*,*]         = oPSD['DATA']
			4: oSpec[i,*,*,*,*]       = oPSD['DATA']
			5: oSpec[i,*,*,*,*,*]     = oPSD['DATA']
			6: oSpec[i,*,*,*,*,*,*]   = oPSD['DATA']
			7: oSpec[i,*,*,*,*,*,*,*] = oPSD['DATA']
			ELSE: Message, 'Invalid number of demensions for SELF.'
		ENDCASE
		
		;Advance to next interval
		iStart += tshift
		iEnd    = iStart + tfft - 1
	ENDFOR
	
	oSpec['SCALE'] = 1B
	oSpec['LOG']   = 1B
	oFreq          = oSpec['DEPEND_1']
	oFreq['TITLE'] = 'Freq!C(Hz)'
	
	RETURN, oSpec
END


;+
;   Transpose the array. Attributes are copied to the results.
;
; :Keywords:
;       P:          in, optional, type=intarr
;                   Vector specifying how the dimensions are to be permuted. Dimensions
;                       start at zero and cannot be repeated. The time dimension must
;                       remain first (i.e. zero must be the first element of P). Zero
;                       may be omitted from P so that P specifies how the non-time
;                       dimensions are permuted.
;-
function MrTimeSeries::Transpose, P
	compile_opt idl2
	on_error, 2

;-----------------------------------------------------
; P Not Given \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;If P is not given, the default is to reverse dimensions
	;   - Time dimension is always first
	if n_elements(p) eq 0 then begin
		;Transpose only does something if there are 3+ dimensions
		nDims = size(self, /N_DIMENSIONS)
		if nDims gt 2 then begin
			permute = [0, reverse( indgen( nDims-1 ) + 1 )]
		
		;Otherwise, keep dimensions the same
		endif else begin
			permute = indgen(nDims)
		endelse
;-----------------------------------------------------
; Permute as Prescribed \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	endif else begin
		permute = P
		
		;If 0 is present, it must be first. If not present, add it first.
		if permute[0] ne 1 then begin
			iZero = where(permute eq 0, nZero)
			if nZero eq 0 $
				then permute = [0, P] $
				else message, 'Zero must be the first element of P.'
		endif
	endelse
	
;-----------------------------------------------------
; Transpose \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Create a time series object
	result  = self -> Copy('Transpose(' + self.name + ')')
	result -> SetData, self.oTime, transpose(*self.data, permute)
	
	;Return
	return, result
end


;+
;   The class definition statement.
;
; :Params:
;       CLASS:          out, optional, type=structure
;-
pro MrTimeSeries__DEFINE
	compile_opt idl2
	
	class = { MrTimeSeries, $
	          inherits MrVariable, $
	          oTime: obj_new() $
	        }
end