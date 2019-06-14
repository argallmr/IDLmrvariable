; docformat = 'rst'
;
; NAME:
;   MrMatrixTS__Define
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
;           MrTimeSeries [op] Expression   = MrVariable
;
;       If the operation is between two time-series objects, the following
;       applies:
;
;           MrMatrixTS [op] MrScalarTS = MrMatrixTS
;           MrMatrixTS [op] MrVectorTS = Error
;           MrMatrixTS [op] MrMatrixTS = MrMatrixTS
;
;   CACHING
;       The ::Cache method will cache the MrMatrixTS object. All variables in the
;       cache are forced to have unique names. When the object is destroyed, it
;       will automatically be removed from the cache.
;
; :Categories:
;   MrVariable, MrTimeSeries, MrMatrixTS
;
; :See Also:
;   MrVariable__Define.pro
;   MrTimeSeries__Define.pro
;   MrScalarTS__Define.pro
;   MrVectorTS__Define.pro
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
;       2016/07/01  -   Written by Matthew Argall
;       2016/10/24  -   Major rewrite to inherit MrTimeSeries object. - MRA
;       2017/03/31  -   Testing revealed VAR->GetData() is faster than VAR['DATA'],
;                           so the change was made in all ::_Overload methods. - MRA
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
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrTimeSeries::Init is accepted here via
;                           keyword inheritance.
;-
function MrMatrixTS::INIT, time, data, $
NAME=name, $
_REF_EXTRA=extra
	compile_opt idl2

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return, 0
	endif
	
	;Defaults
	if n_elements(name) eq 0 then name = 'MrMatrixTS'
	
	;Initialize superclass
	success = self -> MrTimeSeries::Init( time, data, $
	                                      NAME          = name, $
	                                      _STRICT_EXTRA = extra )
	if ~success then message, 'Unable to initialize superclass.'

	return, 1
end


;+
;   Clean up after the object is destroyed
;-
pro MrMatrixTS::CLEANUP
	compile_opt idl2
	on_error, 2
	
	;Superclass
	self -> MrTimeSeries::Cleanup
end


;+
;   Multiply two expressions together.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression or a MrVariable object,
;     normal IDL operations will take effect. This means the output will
;     be the size of the array with the smallest dimensions.
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
FUNCTION MrMatrixTS::_OverloadAsterisk, left, right
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------------------
; Two MrTimeSeries Objects /////////////////////////////
;-------------------------------------------------------
	IF IsA(left, 'MrTimeSeries') && IsA(right, 'MrTimeSeries') THEN BEGIN
		;Pick a name and a side
		name = 'Multiply(' + left.name + ',' + right.name + ')'
		side = self -> _LeftOrRight(left, right)
		
		;Distinguish between the SELF object and the OTHER object
		oOther = side EQ 'LEFT' ? right : left
		type   = Obj_Class(oOther)
		
		;Number of points
		;   - If times are different, output will be truncated to fewer number of points
		;   - TODO: Check if time stamps are equal?
		;           If different: Interpolate? Error?
		nSelf  = self   -> GetNPts()
		nOther = oOther -> GetNPts()
		oTime  = nSelf LT nOther ? self['TIMEVAR'] : oOther['TIMEVAR']
		nPts   = nSelf < nOther
		
		;SCALAR
		IF Obj_IsA(type, 'MrScalarTS') THEN BEGIN
			dims     = Size(*self.data, /DIMENSIONS)
			dims[0]  = nOther
			outClass = 'MrMatrixTS'
		
			;Operate
			IF side EQ 'LEFT' $
				THEN temp = (*self.data) * Rebin( (right -> GetData()), dims ) $
				ELSE temp = Rebin( (left -> GetData()), dims ) * (*self.data)
			
		;VECTOR
		ENDIF ELSE IF Obj_IsA(type, 'MrVectorTS') THEN BEGIN
			Message, 'The * operation is not valid between MrMatrixTS & MrVectorTS objects.'
		
		;MATRIX
		ENDIF ELSE IF Obj_IsA(type, 'MrMatrixTS') THEN BEGIN
			outClass = 'MrMatrixTS'
			temp     = (left -> GetData()) * (right -> GetData())
		
		;OTHER
		ENDIF ELSE BEGIN
			result   = self -> MrTimeSeries::_OverloadAsterisk(left, right)
			outClass = ''
		ENDELSE
		
		;Create a new variable
		IF outClass NE '' THEN BEGIN
			result = Obj_New( outClass, oTime, temp, $
			                  NAME = name, $
			                  /NO_COPY )
		ENDIF

;-------------------------------------------------------
; Expression or MrVariable /////////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		result   = self -> MrTimeSeries::_OverloadAsterisk(left, right)
	ENDELSE

;-------------------------------------------------------
; Done /////////////////////////////////////////////////
;-------------------------------------------------------
	RETURN, result
END


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
FUNCTION MrMatrixTS::_OverloadBracketsRightSide, isRange, i1, i2, i3, i4
	Compile_Opt idl2
	On_Error, 2
	
	;Number of subscripts given
	nSubs = N_Elements(isRange)
	
	;String operations
	IF IsA(i1, /SCALAR, 'STRING') THEN BEGIN
		RETURN, self -> MrTimeSeries::_OverloadBracketsRightSide(isRange, i1, i2, i3, i4, i5, i6, i7, i8)
	
	;Scalar operations
	;   - 0   returns the self object
	;   - [0] returns the first data element
	;   - All other cases return data
	ENDIF ELSE IF nSubs EQ 1 && isRange[0] EQ 0 && IsA(i1, /SCALAR) && i1 EQ 0 THEN BEGIN
		RETURN, self
	ENDIF

;---------------------------------------------------------------------
; Extract the Subarray ///////////////////////////////////////////////
;---------------------------------------------------------------------
	IF nSubs GT 3 THEN Message, 'Too many subscripts given.'
	
	;Data
	data = self -> _OverloadBracketsRightSide_Data(isRange, i1, i2, i3, i4)

	;Attributes
	outDims  = Size(data, /DIMENSIONS)
	nOutDims = Size(data, /N_DIMENSIONS)
	attrs    = self -> _OverloadBracketsRightSide_Attrs(outDims, isRange, i1, i2, i3, i4)
	
;---------------------------------------------------------------------
; Create Output //////////////////////////////////////////////////////
;---------------------------------------------------------------------
	;Which MrTimeSeries object
	;   - The implicit array should start as Nx3
	;   - Scalars have a dimension size of 0 with 1 element.
	nDep0 = N_Elements( attrs['DEPEND_0'] )
	IF (outDims[0] EQ nDep0) || (N_Elements(data) EQ 1 && nDep0 EQ 1) THEN BEGIN
		;Scalar
		IF nOutDims EQ 1 THEN BEGIN
			vOut = MrScalarTS( attrs['DEPEND_0'], data, $
			                   NAME='OverloadBRS(' + self.name + ')', $
			                   /NO_COPY )
		
		;Matrix
		ENDIF ELSE IF nOutDims EQ 3 THEN BEGIN
			vOut = MrVectorTS( attrs['DEPEND_0'], data, $
			                   NAME='OverloadBRS(' + self.name + ')', $
			                   /NO_COPY )
		
		;TimeSeries
		ENDIF ELSE BEGIN
			vOut = MrTimeSeries( attrs['DEPEND_0'], data, $
			                     NAME='OverloadBRS(' + self.name + ')', $
			                     /NO_COPY )
		ENDELSE
	
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
;   NOTE:
;     If one of LEFT or RIGHT is an expression or a MrVariable object,
;     normal IDL operations will take effect. This means the output will
;     be the size of the array with the smallest dimensions.
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
FUNCTION MrMatrixTS::_OverloadCaret, left, right
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------------------
; Two MrTimeSeries Objects /////////////////////////////
;-------------------------------------------------------
	IF IsA(left, 'MrTimeSeries') && IsA(right, 'MrTimeSeries') THEN BEGIN
		;Pick a name and a side
		name = 'Caret(' + left.name + ',' + right.name + ')'
		side = self -> _LeftOrRight(left, right)
		
		;Distinguish between the SELF object and the OTHER object
		oOther = side EQ 'LEFT' ? right : left
		type   = Obj_Class(oOther)
		
		;Number of points
		;   - If times are different, output will be truncated to fewer number of points
		;   - TODO: Check if time stamps are equal?
		;           If different: Interpolate? Error?
		nSelf  = self   -> GetNPts()
		nOther = oOther -> GetNPts()
		oTime  = nSelf LT nOther ? self['TIMEVAR'] : oOther['TIMEVAR']
		nPts   = nSelf < nOther
		
		;SCALAR
		IF Obj_IsA(type, 'MrScalarTS') THEN BEGIN
			dims     = Size(*self.data, /DIMENSIONS)
			dims[0]  = nOther
			outClass = 'MrMatrixTS'
		
			;Operate
			IF side EQ 'LEFT' $
				THEN temp = (*self.data) ^ Rebin( (right -> GetData()), dims ) $
				ELSE temp = Rebin( (left -> GetData()), dims ) ^ (*self.data)
			
		;VECTOR
		ENDIF ELSE IF Obj_IsA(type, 'MrVectorTS') THEN BEGIN
			Message, 'The ^ operation is not valid between MrMatrixTS & MrVectorTS objects.'
		
		;MATRIX
		ENDIF ELSE IF Obj_IsA(type, 'MrMatrixTS') THEN BEGIN
			outClass = 'MrMatrixTS'
			temp     = (left -> GetData()) ^ (right -> GetData())
		
		;OTHER
		ENDIF ELSE BEGIN
			result   = self -> MrTimeSeries::_OverloadCaret(left, right)
			outClass = ''
		ENDELSE
		
		;Create a new variable
		IF outClass NE '' THEN BEGIN
			result = Obj_New( outClass, oTime, temp, $
			                  NAME = name, $
			                  /NO_COPY )
		ENDIF

;-------------------------------------------------------
; Expression or MrVariable /////////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		result = self -> MrTimeSeries::_OverloadCaret(left, right)
	ENDELSE

;-------------------------------------------------------
; Done /////////////////////////////////////////////////
;-------------------------------------------------------
	RETURN, result
END


;+
;   The purpose of this method is to provide iterations to the FOREACH operator. Iterate
;   over all of the elements in the array.
;
; :Params:
;       VALUE:              out, required, type=scalar
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
function MrMatrixTS::_OverloadForeach, value, key
	compile_opt idl2
	on_error, 2

	nPts = n_elements( (*self.data)[*,0,0] )
	if n_elements(key) eq 0 then key = 0 else key += 1

	;Get the array element if the index is in range
	if key lt nPts then begin
		next = 1
		value = reform((*self.data)[key,*,*])
	
	;Otherwise, stop iterating 
	endif else next = 0

	return, next
end
 

;+
;   Subtract one expression from another.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression or a MrVariable object,
;     normal IDL operations will take effect. This means the output will
;     be the size of the array with the smallest dimensions.
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
FUNCTION MrMatrixTS::_OverloadMinus, left, right
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------------------
; Two MrTimeSeries Objects /////////////////////////////
;-------------------------------------------------------
	IF IsA(left, 'MrTimeSeries') && IsA(right, 'MrTimeSeries') THEN BEGIN
		;Pick a name and a side
		name = 'Minus(' + left.name + ',' + right.name + ')'
		side = self -> _LeftOrRight(left, right)
		
		;Distinguish between the SELF object and the OTHER object
		oOther = side EQ 'LEFT' ? right : left
		type   = Obj_Class(oOther)
		
		;Number of points
		;   - If times are different, output will be truncated to fewer number of points
		;   - TODO: Check if time stamps are equal?
		;           If different: Interpolate? Error?
		nSelf  = self   -> GetNPts()
		nOther = oOther -> GetNPts()
		oTime  = nSelf LT nOther ? self['TIMEVAR'] : oOther['TIMEVAR']
		nPts   = nSelf < nOther
		
		;SCALAR
		IF Obj_IsA(type, 'MrScalarTS') THEN BEGIN
			dims     = Size(*self.data, /DIMENSIONS)
			dims[0]  = nOther
			outClass = 'MrMatrixTS'
		
			;Operate
			IF side EQ 'LEFT' $
				THEN temp = (*self.data) - Rebin( (right -> GetData()), dims ) $
				ELSE temp = Rebin( (left -> GetData()), dims ) - (*self.data)
			
		;VECTOR
		ENDIF ELSE IF Obj_IsA(type, 'MrVectorTS') THEN BEGIN
			Message, 'The - operation is not valid between MrMatrixTS & MrVectorTS objects.'
		
		;MATRIX
		ENDIF ELSE IF Obj_IsA(type, 'MrMatrixTS') THEN BEGIN
			outClass = 'MrMatrixTS'
			temp     = (left -> GetData()) - (right -> GetData())
		
		;OTHER
		ENDIF ELSE BEGIN
			result   = self -> MrTimeSeries::_OverloadMinus(left, right)
			outClass = ''
		ENDELSE
		
		;Create a new variable
		IF outClass NE '' THEN BEGIN
			result = Obj_New( outClass, oTime, temp, $
			                  NAME = name, $
			                  /NO_COPY )
		ENDIF

;-------------------------------------------------------
; Expression or MrVariable /////////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		result = self -> MrTimeSeries::_OverloadMinus(left, right)
	ENDELSE

;-------------------------------------------------------
; Done /////////////////////////////////////////////////
;-------------------------------------------------------
	RETURN, result
END


;+
;   Take the modulus of one expression with respect to another.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression or a MrVariable object,
;     normal IDL operations will take effect. This means the output will
;     be the size of the array with the smallest dimensions.
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
FUNCTION MrMatrixTS::_OverloadMOD, left, right
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------------------
; Two MrTimeSeries Objects /////////////////////////////
;-------------------------------------------------------
	IF IsA(left, 'MrTimeSeries') && IsA(right, 'MrTimeSeries') THEN BEGIN
		;Pick a name and a side
		name = 'Mod(' + left.name + ',' + right.name + ')'
		side = self -> _LeftOrRight(left, right)
		
		;Distinguish between the SELF object and the OTHER object
		oOther = side EQ 'LEFT' ? right : left
		type   = Obj_Class(oOther)
		
		;Number of points
		;   - If times are different, output will be truncated to fewer number of points
		;   - TODO: Check if time stamps are equal?
		;           If different: Interpolate? Error?
		nSelf  = self   -> GetNPts()
		nOther = oOther -> GetNPts()
		oTime  = nSelf LT nOther ? self['TIMEVAR'] : oOther['TIMEVAR']
		nPts   = nSelf < nOther
		
		;SCALAR
		IF Obj_IsA(type, 'MrScalarTS') THEN BEGIN
			dims     = Size(*self.data, /DIMENSIONS)
			dims[0]  = nOther
			outClass = 'MrMatrixTS'
		
			;Operate
			IF side EQ 'LEFT' $
				THEN temp = (*self.data) MOD Rebin( (right -> GetData()), dims ) $
				ELSE temp = Rebin( (left -> GetData()), dims ) MOD (*self.data)
			
		;VECTOR
		ENDIF ELSE IF Obj_IsA(type, 'MrVectorTS') THEN BEGIN
			Message, 'The MOD operation is not valid between MrMatrixTS & MrVectorTS objects.'
		
		;MATRIX
		ENDIF ELSE IF Obj_IsA(type, 'MrMatrixTS') THEN BEGIN
			outClass = 'MrMatrixTS'
			temp     = (left -> GetData()) MOD (right -> GetData())
		
		;OTHER
		ENDIF ELSE BEGIN
			result   = self -> MrTimeSeries::_OverloadMOD(left, right)
			outClass = ''
		ENDELSE
		
		;Create a new variable
		IF outClass NE '' THEN BEGIN
			result = Obj_New( outClass, oTime, temp, $
			                  NAME = name, $
			                  /NO_COPY )
		ENDIF

;-------------------------------------------------------
; Expression or MrVariable /////////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		result = self -> MrTimeSeries::_OverloadMOD(left, right)
	ENDELSE

;-------------------------------------------------------
; Done /////////////////////////////////////////////////
;-------------------------------------------------------
	RETURN, result
END


;+
;   Add one expression to another.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression or a MrVariable object,
;     normal IDL operations will take effect. This means the output will
;     be the size of the array with the smallest dimensions.
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
FUNCTION MrMatrixTS::_OverloadPlus, left, right
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------------------
; Two MrTimeSeries Objects /////////////////////////////
;-------------------------------------------------------
	IF IsA(left, 'MrTimeSeries') && IsA(right, 'MrTimeSeries') THEN BEGIN
		;Pick a name and a side
		name = 'Plus(' + left.name + ',' + right.name + ')'
		side = self -> _LeftOrRight(left, right)
		
		;Distinguish between the SELF object and the OTHER object
		oOther = side EQ 'LEFT' ? right : left
		type   = Obj_Class(oOther)
		
		;Number of points
		;   - If times are different, output will be truncated to fewer number of points
		;   - TODO: Check if time stamps are equal?
		;           If different: Interpolate? Error?
		nSelf  = self   -> GetNPts()
		nOther = oOther -> GetNPts()
		oTime  = nSelf LT nOther ? self['TIMEVAR'] : oOther['TIMEVAR']
		nPts   = nSelf < nOther
		
		;SCALAR
		IF Obj_IsA(type, 'MrScalarTS') THEN BEGIN
			dims     = Size(*self.data, /DIMENSIONS)
			dims[0]  = nOther
			outClass = 'MrMatrixTS'
		
			;Operate
			IF side EQ 'LEFT' $
				THEN temp = (*self.data) + Rebin( (right -> GetData()), dims ) $
				ELSE temp = Rebin( (left -> GetData()), dims ) + (*self.data)
			
		;VECTOR
		ENDIF ELSE IF Obj_IsA(type, 'MrVectorTS') THEN BEGIN
			Message, 'The ^ operation is not valid between MrMatrixTS & MrVectorTS objects.'
		
		;MATRIX
		ENDIF ELSE IF Obj_IsA(type, 'MrMatrixTS') THEN BEGIN
			outClass = 'MrMatrixTS'
			temp     = (left -> GetData()) + (right -> GetData())
		
		;OTHER
		ENDIF ELSE BEGIN
			result   = self -> MrTimeSeries::_OverloadPlus(left, right)
			outClass = ''
		ENDELSE
		
		;Create a new variable
		IF outClass NE '' THEN BEGIN
			result = Obj_New( outClass, oTime, temp, $
			                  NAME = name, $
			                  /NO_COPY )
		ENDIF

;-------------------------------------------------------
; Expression or MrVariable /////////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		result = self -> MrTimeSeries::_OverloadPlus(left, right)
	ENDELSE

;-------------------------------------------------------
; Done /////////////////////////////////////////////////
;-------------------------------------------------------
	RETURN, result
END


;+
;   Compute the inner product of two expressions.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression or a MrVariable object,
;     normal IDL operations will take effect. This means the output will
;     be the size of the array with the smallest dimensions.
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
FUNCTION MrMatrixTS::_OverloadPound, left, right
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------------------
; Two MrTimeSeries Objects /////////////////////////////
;-------------------------------------------------------
	IF IsA(left, 'MrTimeSeries') && IsA(right, 'MrTimeSeries') THEN BEGIN
		;Pick a name and a side
		name = 'Pound(' + left.name + ',' + right.name + ')'
		side = self -> _LeftOrRight(left, right)
		
		;Distinguish between the SELF object and the OTHER object
		oOther   = side EQ 'LEFT' ? right : left
		type     = Obj_Class(oOther)
		
		;Number of points
		;   - If times are different, output will be truncated to fewer number of points
		;   - TODO: Check if time stamps are equal?
		;           If different: Interpolate? Error?
		nSelf  = self -> GetNPts()
		nOther = oOther -> GetNPts()
		oTime  = nSelf LT nOther ? self['TIMEVAR'] : oOther['TIMEVAR']
		nPts   = nSelf < nOther
		
		;SCALAR
		IF Obj_IsA(type, 'MrScalarTS') THEN BEGIN
			Message, 'The # operation is not valid between MrMatrixTS & MrScalarTS objects.'
			
		;VECTOR
		ENDIF ELSE IF Obj_IsA(type, 'MrVectorTS') THEN BEGIN
			mDims = Size( *self.data, /DIMENSIONS )
			vDims = Size( oOther, /DIMENSIONS )
			
			;Inner product
			;   - Reform vector from Nx3 to NxMx3
			IF side EQ 'LEFT' THEN BEGIN
				;Make sure the matrix is the correct size
				IF mDims[2] NE 3 THEN Message, 'MrMatrixTS dimensions must be NxMx3.'
				outClass = mDims[1] EQ 3 ? 'MrVectorTS' : 'MrTimeSeries'
				
				;Operate
				temp = Total( *self.data * $
				              Rebin( Reform((right -> GetData()), [vDims[0], 1, vDims[1]]), [vDims[0], mDims[1], vDims[1]] ), $
				              3, /PRESERVE_TYPE )
			
			;NX3XM
			ENDIF ELSE BEGIN
				IF mDims[1] NE 3 THEN Message, 'MrMatrixTS dimensions must be Nx3xM.'
				outClass = mDims[2] EQ 3 ? 'MrVectorTS' : 'MrTimeSeries'
				
				;Operate
				temp  = Total( Rebin( (left -> GetData()), [vDims[0], vDims[1], mDims[2]] ) * $
				               *self.data, $
				               2, /PRESERVE_TYPE )
			ENDELSE
			
		
		;MATRIX
		ENDIF ELSE IF Obj_IsA(type, 'MrMatrixTS') THEN BEGIN
			lDims    = Size( left, /DIMENSIONS )
			rDims    = Size( right, /DIMENSIONS )
			outClass = 'MrMatrixTS'

			;LEFT=Nx3x3, RIGHT=Nx3x3
			IF Array_Equal(lDims[1:2], 3) && Array_Equal(rDims[1:2], 3) THEN BEGIN
				IF side EQ 'LEFT' THEN BEGIN
					temp = [ [ [(*self.data)[*,0,0]*right['DATA',*,0,0] + (*self.data)[*,0,1]*right['DATA',*,1,0] + (*self.data)[*,0,2]*right['DATA',*,2,0]], $
					           [(*self.data)[*,1,0]*right['DATA',*,0,0] + (*self.data)[*,1,1]*right['DATA',*,1,0] + (*self.data)[*,1,2]*right['DATA',*,2,0]], $
					           [(*self.data)[*,2,0]*right['DATA',*,0,0] + (*self.data)[*,2,1]*right['DATA',*,1,0] + (*self.data)[*,2,2]*right['DATA',*,2,0]] ], $
					          
					         [ [(*self.data)[*,0,0]*right['DATA',*,0,1] + (*self.data)[*,0,1]*right['DATA',*,1,1] + (*self.data)[*,0,2]*right['DATA',*,2,1]], $
					           [(*self.data)[*,1,0]*right['DATA',*,0,1] + (*self.data)[*,1,1]*right['DATA',*,1,1] + (*self.data)[*,1,2]*right['DATA',*,2,1]], $
					           [(*self.data)[*,2,0]*right['DATA',*,0,1] + (*self.data)[*,2,1]*right['DATA',*,1,1] + (*self.data)[*,2,2]*right['DATA',*,2,1]] ], $
					           
					         [ [(*self.data)[*,0,0]*right['DATA',*,0,2] + (*self.data)[*,0,1]*right['DATA',*,1,2] + (*self.data)[*,0,2]*right['DATA',*,2,2]], $
					           [(*self.data)[*,1,0]*right['DATA',*,0,2] + (*self.data)[*,1,1]*right['DATA',*,1,2] + (*self.data)[*,1,2]*right['DATA',*,2,2]], $
					           [(*self.data)[*,2,0]*right['DATA',*,0,2] + (*self.data)[*,2,1]*right['DATA',*,1,2] + (*self.data)[*,2,2]*right['DATA',*,2,2]] ] ]
				ENDIF ELSE BEGIN
					temp = [ [ [left['DATA',*,0,0]*(*self.data)[*,0,0] + left['DATA',*,0,1]*(*self.data)[*,1,0] + left['DATA',*,0,2]*(*self.data)[*,2,0]], $
					           [left['DATA',*,1,0]*(*self.data)[*,0,0] + left['DATA',*,1,1]*(*self.data)[*,1,0] + left['DATA',*,1,2]*(*self.data)[*,2,0]], $
					           [left['DATA',*,2,0]*(*self.data)[*,0,0] + left['DATA',*,2,1]*(*self.data)[*,1,0] + left['DATA',*,2,2]*(*self.data)[*,2,0]] ], $
					          
					         [ [left['DATA',*,0,0]*(*self.data)[*,0,1] + left['DATA',*,0,1]*(*self.data)[*,1,1] + left['DATA',*,0,2]*(*self.data)[*,2,1]], $
					           [left['DATA',*,1,0]*(*self.data)[*,0,1] + left['DATA',*,1,1]*(*self.data)[*,1,1] + left['DATA',*,1,2]*(*self.data)[*,2,1]], $
					           [left['DATA',*,2,0]*(*self.data)[*,0,1] + left['DATA',*,2,1]*(*self.data)[*,1,1] + left['DATA',*,2,2]*(*self.data)[*,2,1]] ], $
					           
					         [ [left['DATA',*,0,0]*(*self.data)[*,0,2] + left['DATA',*,0,1]*(*self.data)[*,1,2] + left['DATA',*,0,2]*(*self.data)[*,2,2]], $
					           [left['DATA',*,1,0]*(*self.data)[*,0,2] + left['DATA',*,1,1]*(*self.data)[*,1,2] + left['DATA',*,1,2]*(*self.data)[*,2,2]], $
					           [left['DATA',*,2,0]*(*self.data)[*,0,2] + left['DATA',*,2,1]*(*self.data)[*,1,2] + left['DATA',*,2,2]*(*self.data)[*,2,2]] ] ]
				
				ENDELSE
			
			;LEFT=anything else; RIGHT=anything else
			ENDIF ELSE BEGIN
				;Allocate memory to results
				temp = Make_Array(nPts, lDims[1], rDims[2], $
				                  TYPE=Size( Reform(left['DATA',0,*,*])#Reform(right['DATA',0,*,*]), /TYPE))
				
				;Compute result
				;   - Increment to the least number OF elements
				FOREACH mat, left, idx do BEGIN
					temp[idx,*,*] = mat # Reform(right['DATA',idx,*,*])
					IF idx EQ nPts-1 THEN BREAK
				ENDFOREACH
			ENDELSE
		
		;OTHER
		ENDIF ELSE BEGIN
			result   = self -> MrTimeSeries::_OverloadPound(left, right)
			outClass = ''
		ENDELSE
		
		;Create a new variable
		IF outClass NE '' THEN BEGIN
			result = Obj_New( outClass, oTime, temp, $
			                  NAME = name, $
			                  /NO_COPY )
		ENDIF

;-------------------------------------------------------
; Expression or MrVariable /////////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		result = self -> MrTimeSeries::_OverloadPound(left, right)
	ENDELSE

;-------------------------------------------------------
; Done /////////////////////////////////////////////////
;-------------------------------------------------------
	RETURN, result
END


;+
;   Compute the outer product of two expressions.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression or a MrVariable object,
;     normal IDL operations will take effect. This means the output will
;     be the size of the array with the smallest dimensions.
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
FUNCTION MrMatrixTS::_OverloadPoundPound, left, right
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------------------
; Two MrTimeSeries Objects /////////////////////////////
;-------------------------------------------------------
	IF IsA(left, 'MrTimeSeries') && IsA(right, 'MrTimeSeries') THEN BEGIN
		;Pick a name and a side
		name = 'PoundPound(' + left.name + ',' + right.name + ')'
		side = self -> _LeftOrRight(left, right)
		
		;Distinguish between the SELF object and the OTHER object
		oOther = side EQ 'LEFT' ? right : left
		type   = Obj_Class(oOther)
		
		;Number of points
		;   - If times are different, output will be truncated to fewer number of points
		;   - TODO: Check if time stamps are equal?
		;           If different: Interpolate? Error?
		nSelf  = self -> GetNPts()
		nOther = oOther -> GetNPts()
		oTime  = nSelf LT nOther ? self['TIMEVAR'] : oOther['TIMEVAR']
		nPts   = nSelf < nOther
		
		;SCALAR
		IF Obj_IsA(type, 'MrScalarTS') THEN BEGIN
			Message, 'The ## operation is not valid between MrMatrixTS & MrScalarTS objects.'
			
		;VECTOR
		ENDIF ELSE IF Obj_IsA(type, 'MrVectorTS') THEN BEGIN
			mDims = Size( *self.data, /DIMENSIONS )
			vDims = Size( oOther, /DIMENSIONS )
			
			;Outer product
			;   - Reform vector from Nx3 to NxMx3
			IF side EQ 'LEFT' THEN BEGIN
				IF mDims[1] NE 3 THEN Message, 'MrMatrixTS dimensions must be NxMx3.'
				outClass = mDims[2] EQ 3 ? 'MrVectorTS' : 'MrTimeSeries'
				
				;Operate
				temp  = Total( *self.data * $
				               Rebin( (right -> GetData()), [vDims[0], vDims[1], mDims[2]] ), $
				               2, /PRESERVE_TYPE )
			
			;NX3XM
			ENDIF ELSE BEGIN
				;Make sure the matrix is the correct size
				IF mDims[2] NE 3 THEN Message, 'MrMatrixTS dimensions must be Nx3xM.'
				outClass = mDims[1] EQ 3 ? 'MrVectorTS' : 'MrTimeSeries'
				
				;Operate
				temp = Total( Rebin( Reform((right -> GetData()), [vDims[0], 1, vDims[1]]), [vDims[0], mDims[1], vDims[1]] ) * $
				              *self.data, $
				              3, /PRESERVE_TYPE )
			ENDELSE
			
		
		;MATRIX
		ENDIF ELSE IF Obj_IsA(type, 'MrMatrixTS') THEN BEGIN
			lDims    = Size( left, /DIMENSIONS )
			rDims    = Size( right, /DIMENSIONS )
			outClass = 'MrMatrixTS'
			
			;LEFT=Nx3x3, RIGHT=Nx3x3
			IF Array_Equal(lDims[1:2], 3) && Array_Equal(rDims[1:2], 3) THEN BEGIN
				IF side EQ 'LEFT' THEN BEGIN
					temp = [ [ [(*self.data)[*,0,0]*right['DATA',*,0,0] + (*self.data)[*,1,0]*right['DATA',*,0,1] + (*self.data)[*,2,0]*right['DATA',*,0,2]], $
					           [(*self.data)[*,0,0]*right['DATA',*,1,0] + (*self.data)[*,1,0]*right['DATA',*,1,1] + (*self.data)[*,2,0]*right['DATA',*,1,2]], $
					           [(*self.data)[*,0,0]*right['DATA',*,2,0] + (*self.data)[*,1,0]*right['DATA',*,2,1] + (*self.data)[*,2,0]*right['DATA',*,2,2]] ], $
					         
					         [ [(*self.data)[*,0,1]*right['DATA',*,0,0] + (*self.data)[*,1,1]*right['DATA',*,0,1] + (*self.data)[*,2,1]*right['DATA',*,0,2]], $
					           [(*self.data)[*,0,1]*right['DATA',*,1,0] + (*self.data)[*,1,1]*right['DATA',*,1,1] + (*self.data)[*,2,1]*right['DATA',*,1,2]], $
					           [(*self.data)[*,0,1]*right['DATA',*,2,0] + (*self.data)[*,1,1]*right['DATA',*,2,1] + (*self.data)[*,2,1]*right['DATA',*,2,2]] ], $
					         
					         [ [(*self.data)[*,0,2]*right['DATA',*,0,0] + (*self.data)[*,1,2]*right['DATA',*,0,1] + (*self.data)[*,2,2]*right['DATA',*,0,2]], $
					           [(*self.data)[*,0,2]*right['DATA',*,1,0] + (*self.data)[*,1,2]*right['DATA',*,1,1] + (*self.data)[*,2,2]*right['DATA',*,1,2]], $
					           [(*self.data)[*,0,2]*right['DATA',*,2,0] + (*self.data)[*,1,2]*right['DATA',*,2,1] + (*self.data)[*,2,2]*right['DATA',*,2,2]] ] ]
				ENDIF ELSE BEGIN
					temp = [ [ [left['DATA',*,0,0]*(*self.data)[*,0,0] + left['DATA',*,1,0]*(*self.data)[*,0,1] + left['DATA',*,2,0]*(*self.data)[*,0,2]], $
					           [left['DATA',*,0,0]*(*self.data)[*,1,0] + left['DATA',*,1,0]*(*self.data)[*,1,1] + left['DATA',*,2,0]*(*self.data)[*,1,2]], $
					           [left['DATA',*,0,0]*(*self.data)[*,2,0] + left['DATA',*,1,0]*(*self.data)[*,2,1] + left['DATA',*,2,0]*(*self.data)[*,2,2]] ], $
					         
					         [ [left['DATA',*,0,1]*(*self.data)[*,0,0] + left['DATA',*,1,1]*(*self.data)[*,0,1] + left['DATA',*,2,1]*(*self.data)[*,0,2]], $
					           [left['DATA',*,0,1]*(*self.data)[*,1,0] + left['DATA',*,1,1]*(*self.data)[*,1,1] + left['DATA',*,2,1]*(*self.data)[*,1,2]], $
					           [left['DATA',*,0,1]*(*self.data)[*,2,0] + left['DATA',*,1,1]*(*self.data)[*,2,1] + left['DATA',*,2,1]*(*self.data)[*,2,2]] ], $
					         
					         [ [left['DATA',*,0,2]*(*self.data)[*,0,0] + left['DATA',*,1,2]*(*self.data)[*,0,1] + left['DATA',*,2,2]*(*self.data)[*,0,2]], $
					           [left['DATA',*,0,2]*(*self.data)[*,1,0] + left['DATA',*,1,2]*(*self.data)[*,1,1] + left['DATA',*,2,2]*(*self.data)[*,1,2]], $
					           [left['DATA',*,0,2]*(*self.data)[*,2,0] + left['DATA',*,1,2]*(*self.data)[*,2,1] + left['DATA',*,2,2]*(*self.data)[*,2,2]] ] ]
				ENDELSE
			
			;LEFT=anything else; RIGHT=anything else
			ENDIF ELSE BEGIN
				;Allocate memory to results
				temp = Make_Array(nPts, lDims[2], rDims[1], $
				                  TYPE=Size( Reform(left['DATA',0,*,*])##Reform(right['DATA',0,*,*]), /TYPE))
				
				;Compute result
				;   - Increment to the least number OF elements
				FOREACH mat, left, idx do BEGIN
					temp[idx,*,*] = mat ## Reform(right['DATA',idx,*,*])
					IF idx EQ nPts-1 THEN BREAK
				ENDFOREACH
			ENDELSE
		
		;OTHER
		ENDIF ELSE BEGIN
			result   = self -> MrTimeSeries::_OverloadPoundPound(left, right)
			outClass = ''
		ENDELSE
		
		;Create a new variable
		IF outClass NE '' THEN BEGIN
			result = Obj_New( outClass, oTime, temp, $
			                  NAME = name, $
			                  /NO_COPY )
		ENDIF

;-------------------------------------------------------
; Expression or MrVariable /////////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		result = self -> MrTimeSeries::_OverloadPoundPound(left, right)
	ENDELSE

;-------------------------------------------------------
; Done /////////////////////////////////////////////////
;-------------------------------------------------------
	RETURN, result
END


;+
;   The purpose of this method is to provide information when implied print is used.
;   is called.
;-
function MrMatrixTS::_OverloadPrint
	compile_opt idl2
	on_error, 2

	;Print the array
	return, *self.data
end


;+
;   The purpose of this method is to obtain the results of IDL's Size() function when
;   applied to the applicit array.
;
; :Returns:
;       RESULT:             Results of the Size() function on the implicit array
;-
function MrMatrixTS::_OverloadSize
	compile_opt idl2
	on_error, 2

	;NOTE:
	;   There is a bug that incorrectly reports the number of dimesions of an
	;   undefined variable (returns 1 instead of 0). This will be fixed in
	;   version 8.3.1.
	;
	;   https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/-f8Cxlp5cxQ
	;

	;Return the dimensions of the array. The Size and N_Elements functions
	;will know what to do with them.
	return, size(*self.data, /DIMENSIONS)
end


;+
;   Divide one expression by another.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression or a MrVariable object,
;     normal IDL operations will take effect. This means the output will
;     be the size of the array with the smallest dimensions.
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
FUNCTION MrMatrixTS::_OverloadSlash, left, right
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------------------
; Two MrTimeSeries Objects /////////////////////////////
;-------------------------------------------------------
	IF IsA(left, 'MrTimeSeries') && IsA(right, 'MrTimeSeries') THEN BEGIN
		;Pick a name and a side
		name = 'Divide(' + left.name + ',' + right.name + ')'
		side = self -> _LeftOrRight(left, right)
		
		;Distinguish between the SELF object and the OTHER object
		oOther = side EQ 'LEFT' ? right : left
		type   = Obj_Class(oOther)
		
		;Number of points
		;   - If times are different, output will be truncated to fewer number of points
		;   - TODO: Check if time stamps are equal?
		;           If different: Interpolate? Error?
		nSelf  = self   -> GetNPts()
		nOther = oOther -> GetNPts()
		oTime  = nSelf LT nOther ? self['TIMEVAR'] : oOther['TIMEVAR']
		nPts   = nSelf < nOther
		
		;SCALAR
		IF Obj_IsA(type, 'MrScalarTS') THEN BEGIN
			dims     = Size(*self.data, /DIMENSIONS)
			dims[0]  = nOther
			outClass = 'MrMatrixTS'
		
			;Operate
			IF side EQ 'LEFT' $
				THEN temp = (*self.data) / Rebin( (right -> GetData()), dims ) $
				ELSE temp = Rebin( (left -> GetData()), dims ) / (*self.data)
			
		;VECTOR
		ENDIF ELSE IF Obj_IsA(type, 'MrVectorTS') THEN BEGIN
			Message, 'The / operation is not valid between MrMatrixTS & MrVectorTS objects.'
		
		;MATRIX
		ENDIF ELSE IF Obj_IsA(type, 'MrMatrixTS') THEN BEGIN
			outClass = 'MrMatrixTS'
			temp     = (left -> GetData()) / (right -> GetData())
		
		;OTHER
		ENDIF ELSE BEGIN
			result   = self -> MrTimeSeries::_OverloadSlash(left, right)
			outClass = ''
		ENDELSE
		
		;Create a new variable
		IF outClass NE '' THEN BEGIN
			result = Obj_New( outClass, oTime, temp, $
			                  NAME = name, $
			                  /NO_COPY )
		ENDIF

;-------------------------------------------------------
; Expression or MrVariable /////////////////////////////
;-------------------------------------------------------
	ENDIF ELSE BEGIN
		result = self -> MrTimeSeries::_OverloadSlash(left, right)
	ENDELSE

;-------------------------------------------------------
; Done /////////////////////////////////////////////////
;-------------------------------------------------------
	RETURN, result
END


;+
;   A wrapper for the # operator or the Matrix_Multiply function, generalized to rank N
;   tensors. IDL's Matrix_Multiply truncate dimensions beyond rank 2.
;
;   For two IDL tensors, A and B, the last dimension of A will be contracted with the
;   first dimension of B.
;
; :Examples:
;   See the example program at the end of this document::
;       IDL> .r MrInnerProduct
;
; :Params:
;       A:          in, required, type=tensor rank n
;                   Rank N tensor whose last dimension will be contracted with the
;                       first dimension of `B`.
;       B:          in, required, type=tensor rank n
;                   Rank N tensor whose first dimension will be contracted with the
;                       last dimension of `A`.
;
; :Returns:
;       RESULT:     Result of the inner product.
;-
function MrMatrixTS::InnerProduct, A, B
	compile_opt strictarr
	on_error, 2

	;Dimensionality of inputs
	nDimsA = size(A, /N_DIMENSIONS)
	nDimsB = size(B, /N_DIMENSIONS)
	dimsA = size(A, /DIMENSIONS)
	dimsB = size(B, /DIMENSIONS)

	;Reform _A and B so that they have two dimensions
	_A = (nDimsA le 2) ? A : reform(A, product(dimsA[0:nDimsA-2]), dimsA[nDimsA-1])
	_B = (nDimsB le 2) ? B : reform(B, dimsB[0], product(dimsB[1:*]))

	;Compute the inner product
	result = temporary(_A) # temporary(_B)

	;Reform back to the original shape, minus the contracted dimension.
	if nDimsA gt 2 then dimsOut = dimsA[0:nDimsA-2]
	if nDimsB gt 2 then begin
		if n_elements(dimsOut) gt 0 $
			then dimsOut = [dimsOut, dimsB[1:*]] $
			else dimsOut = dimsB[1:*]
	endif
	if n_elements(dimsOut) gt 0 $
		then result = reform(result, dimsOut, /OVERWRITE)

	return, result
end


;+
;   A wrapper for the ## operator, generalized to rank N tensors. IDL's ## operator
;   truncate dimensions beyond rank 2.
;
;   For two IDL tensors, A and B, the first dimension of A will be contracted with the
;   last dimension of B.
;
; :Examples:
;   See the example program at the end of this document::
;       IDL> .r MrOuterProduct
;
; :Params:
;       A:          in, required, type=tensor rank n
;                   Rank N tensor whose first dimension will be contracted with the
;                       last dimension of `B`.
;       B:          in, required, type=tensor rank n
;                   Rank N tensor whose last dimension will be contracted with the
;                       first dimension of `A`.
;
; :Returns:
;       RESULT:     Result of the outer product.
;-
function MrMatrixTS::OuterProduct, A, B
	compile_opt strictarr
	on_error, 2

	;Dimensionality of inputs
	nDimsA = size(A, /N_DIMENSIONS)
	nDimsB = size(B, /N_DIMENSIONS)
	dimsA = size(A, /DIMENSIONS)
	dimsB = size(B, /DIMENSIONS)

	;Reform _A and B so that they have two dimensions
	_A = (nDimsA le 2) ? A : reform(A, dimsA[0], product(dimsA[1:*]))
	_B = (nDimsB le 2) ? B : reform(B, product(dimsB[0:nDimsB-2]), dimsB[nDimsB-1])

	;Compute the inner product
	result = temporary(_A) ## temporary(_B)

	;Reform back to the original shape, minus the contracted dimension.
	if nDimsA gt 2 then dimsOut = dimsA[1:*]
	if nDimsB gt 2 then begin
		if n_elements(dimsOut) gt 0 $
			then dimsOut = [dimsOut, dimsB[0:nDimsB-2]] $
			else dimsOut = dimsB[0:nDimsB-2]
	endif
	if n_elements(dimsOut) gt 0 $
		then result = reform(result, dimsOut, /OVERWRITE)
	
	return, result
end


;+
;   The purpose of this program is to rotate a 3x3 or 3x3xN matrix via a standard
;   rotation::
;
;           A' = R A R^T
;
;       Where R^T is the transpose of the rotation matrix R.
;
;   The tensor product of two fields is::
;                | HxKx  HxKy  HxKz |
;       P = HK = | HyKx  HyKy  HyKz | = HiKj
;                | HzKx  HzKy  HzKz |
;
;   In a different coordinate system, H and K are given by::
;       H'i = Aij Hj     K'l = Alm Km
;
;   The tensor H'K', then, is::
;       P' = H'i K'j = (Aij Hj) (Alm Km)
;                    = Aij Alm (HK)jm
;                    = Aij (Aml)^T (HK)jm
;                    = Aij (HK)jm (Aml)^T
;                    = A HK A^T
;
; :Keywords:
;       MAT:                in, required, type=NxM array
;                           Matrix to be rotated. Can be a MrVariable or MrMatrixTS object.
;
; :Keywords:
;       CACHE:              in, optional, type=boolean, default=0
;                           If set, the returned MrVectorTS will be cached.
;       NAME:               in, optional, type=string, default=<name>_rot
;                           Name to be given to the resulting MrScalar object. 
;-
function MrMatrixTS::Rotate_Matrix, mat, $
CACHE=cache, $
NAME=name
	compile_opt idl2
	on_error, 2

	if n_elements(name) eq 0 then name = self.name + '_rot'

;-----------------------------------------------------
; MrVariable Object \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	if isa(mat, 'MrMatrixTS') then begin
		;A and B must have the same dimensions.
		Adims = size(self, /DIMENSIONS)
		Bdims = size(mat, /DIMENSIONS)
		if ~array_equal(Adims[1:2], Bdims[1:2]) then message, 'MAT is in correct size.'
		
		;Rotate
		;   - self # mat # transpose(self)
		out = self # mat # self -> Transpose()

;-----------------------------------------------------
; Other Array \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	endif else begin
		;self # mat
		temp = self -> InnerProduct(mat)
		
		;(self # mat) # transpose(self)
		out  = temp -> InnerProduct( self -> Transpose() )
	endelse
	
;-----------------------------------------------------
; Output \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	return, out
end


;+
;   Rotate a vector.
;
; :Keywords:
;       VEC:                in, required, type=Nx3 or 3xN array
;                           Array to be rotated. Can be a MrVariable or MrVectorTS object.
;                               If 3x3, it will be treated as Nx3.
;
; :Keywords:
;       CACHE:              in, optional, type=boolean, default=0
;                           If set, the returned MrVectorTS will be cached.
;       NAME:               in, optional, type=string, default=<name>_rot
;                           Name to be given to the resulting MrScalar object. 
;-
function MrMatrixTS::Rotate_Vector, vec, $
CACHE=cache, $
NAME=name
	compile_opt idl2
	on_error, 2

	dims = size(*self.data, /DIMENSIONS)
	if ~array_equal(dims[1:2], 3) then message, 'MrMatrix must be 3x3.'
	if n_elements(name) eq 0 then name = self.name + '_rot'

;-----------------------------------------------------
; MrVariable Object \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	if isa(vec, 'MrVectorTS') then begin
		;MrVectorTS
		pVec = vec['PTR']
		outVec = [ [(*self.data)[*,0,0]*(*pvec)[*,0] + (*self.data)[*,1,0]*(*pvec)[*,1] + (*self.data)[*,2,0]*(*pvec)[*,2]], $
		           [(*self.data)[*,0,1]*(*pvec)[*,0] + (*self.data)[*,1,1]*(*pvec)[*,1] + (*self.data)[*,2,1]*(*pvec)[*,2]], $
		           [(*self.data)[*,0,2]*(*pvec)[*,0] + (*self.data)[*,1,2]*(*pvec)[*,1] + (*self.data)[*,2,2]*(*pvec)[*,2]] ]

;-----------------------------------------------------
; Normal Array \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	endif else if isa(vec, /NUMBER) || isa(vec, 'MrVariable') then begin
		;Check sizes
		sz = size(vec)

		;VEC is 3xN (Treat 3x3 as Nx3)
		if (sz[0] eq 1 && sz[1] eq 3) || (sz[0] eq 2 && sz[1] eq 3 && sz[2] ne 3) then begin
			outVec = [ [(*self.data)[*,0,0]*vec[0,*] + (*self.data)[*,1,0]*vec[1,*] + (*self.data)[*,2,0]*vec[2,*]], $
			           [(*self.data)[*,0,1]*vec[0,*] + (*self.data)[*,1,1]*vec[1,*] + (*self.data)[*,2,1]*vec[2,*]], $
			           [(*self.data)[*,0,2]*vec[0,*] + (*self.data)[*,1,2]*vec[1,*] + (*self.data)[*,2,2]*vec[2,*]] ]
		
		;VEC is Nx3
		endif else if (sz[0] eq 2 && sz[2] eq 3) then begin
			outVec = [ [(*self.data)[*,0,0]*vec[*,0] + (*self.data)[*,1,0]*vec[*,1] + (*self.data)[*,2,0]*vec[*,1]], $
			           [(*self.data)[*,0,1]*vec[*,0] + (*self.data)[*,1,1]*vec[*,1] + (*self.data)[*,2,1]*vec[*,1]], $
			           [(*self.data)[*,0,2]*vec[*,0] + (*self.data)[*,1,2]*vec[*,1] + (*self.data)[*,2,2]*vec[*,1]] ]
		endif else begin
			message, 'VEC must be 3xN or Nx3.'
		endelse
		
;-----------------------------------------------------
; Incorrect Type \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	endif else begin
		message, 'VEC is the incorrect datatype.'
	endelse
	
;-----------------------------------------------------
; Output \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oVec = MrVectorTS(self.oTime, outVec, NAME=name, CACHE=cache, /NO_COPY)

	return, oVec
end


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
;       _REF_EXTRA:     in, optional, type=integer
;                       Any keyword accepted by MrTimeSeries::SetData is accepted here
;                           via keyword inheritance.
;-
pro MrMatrixTS::SetData, time, data, $
_REF_EXTRA=extra
	compile_opt idl2
	on_error, 2
	
	;Keep the old data
	oTime = self.oTime
	pDAta = self.data

	;Use the superclass
	self -> MrTimeSeries::SetData, time, data, _STRICT_EXTRA=extra

	;Check the results
	;   - Make sure it is NxMxL (3 dimensions of any size)
	sz = size(self)
	if sz[0] ne 3 then begin
		ptr_free, self.data
		self.data  = pData
		self.oTime = oTime
		message, 'Invalid dimensions: Data must be an NxMxL matrix time series.'
	endif else begin
		ptr_free, pData
	endelse
end


;+
;   The class definition statement.
;
; :Params:
;       CLASS:          out, optional, type=structure
;-
pro MrMatrixTS__DEFINE
	compile_opt idl2
	
	class = { MrMatrixTS, $
	          inherits MrTimeSeries $
	        }
end