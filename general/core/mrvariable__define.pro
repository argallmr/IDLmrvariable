; docformat = 'rst'
;
; NAME:
;   MrVariable__Define
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
;       An operator (such as AND, EQ, ^, *, etc.) calls the overload
;       method for the object on the left side first. So, for two MrVariable objects,
;       A and B::
;
;               print, A [operator] B
;
;       will call A's _Overload method first.
;
;   ARRAY TRUNCATION
;       If two arrays are compared, results will be truncated to have the same number
;       of elements as the shorter array. If an array and a scalar are compared, each
;       elements of the array is compared against the scalar value. This is the same
;       as for normal IDL variables.
;
;   CACHING
;       The ::Cache method will cache the MrVariable object. When the object is destroyed,
;       it will automatically be removed from the cache.
;
; :Categories:
;   MrVariable, Graphics
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
;       2014/05/09  -   Written by Matthew Argall
;       2016/10/24  -   Added the InnerProduct and OuterProduct methods. - MRA
;       2016/10/28  -   DEPEND_[0-3] and DELTA_(MINUS|PLUS)_VAR attributes can
;                           have MrVariable objects as values. - MRA
;       2017/03/16  -   If pointers (i.e. variable names) are given to attributes, the
;                           pointer is followed and the object reference is stored as
;                           the attribute value. - MRA
;       2017/03/31  -   Testing revealed VAR->GetData() is faster than VAR['DATA'],
;                           so the change was made in all ::_Overload methods.
;                           _OverloadBracketsRightSide returns an object with DEPEND_# and
;                           DELTA_(PLUS|MINUS)_VAR attributes are propertly reduced.
;                           _OverloadBracketsLeftSide accepts MrVariable objects. - MRA
;       2017/08/03  -   UNITS attribute is now saved as an IDLunit object. Units can
;                           be set initially via bracket overloading, but can only be
;                           changed via the ::SetUnits method. - MRA
;-
;*****************************************************************************************
;+
;   The initialization method.
;
; :Examples:
;   Create a MrVariable object using a pre-existing array::
;       theArray = findgen(24, 36)
;       myArray  = MrVariable(theArray)
;       help, myArray
;           MYARRAY     OBJREF      <ObjHeapVar666(MrVariable)>
;             ARRAY       FLOAT     = Array[24, 36]
;
;   Initialize a MrVariable object via Make_Array::
;       myArray = MrVariable(24, 36, TYPE='ULong')
;       help, myArray
;           MYARRAY     OBJREF      <ObjHeapVar668(MrVariable)>
;             ARRAY       ULONG     = Array[24, 36]
;
;   Initialize a MrVariable object via RandomU::
;       myArray = MrVariable(24, 36, SEED=3)
;       print, myarray[0:3]
;           0.897916     0.558249     0.766930     0.589101
;
;   Initialize a MrVariable with a normal, gaussian distribution::
;       myArray = MrVariable(24, 36, SEED=3, /NORMAL)
;
; :Params:
;       DATA:               in, optional, type=any/integer
;                           If an array, then it is the array to be stored and all other
;                               parameters are ignored. If a scalar value is given,
;                               `MAKE` is set to 1 automatically, unless `RANDOM` or `NORMAL`
;                               are in use.
;       D2:                 in, optional, type=integer
;                           Size of the second dimension. Ignored unless `DATA` is scalar.
;       D3:                 in, optional, type=integer
;                           Size of the third dimension. Ignored unless `DATA` is scalar.
;       D4:                 in, optional, type=integer
;                           Size of the fourth dimension. Ignored unless `DATA` is scalar.
;       D5:                 in, optional, type=integer
;                           Size of the fifth dimension. Ignored unless `DATA` is scalar.
;       D6:                 in, optional, type=integer
;                           Size of the sixth dimension. Ignored unless `DATA` is scalar.
;       D7:                 in, optional, type=integer
;                           Size of the seventh dimension. Ignored unless `DATA` is scalar.
;       D8:                 in, optional, type=integer
;                           Size of the eighth dimension. Ignored unless `DATA` is scalar.
;
; :Keywords:
;       CACHE:              in, optional, type=boolean, default=0
;                           If set, the variable will be placed into the cache.
;       MAKE:               in, optional, type=boolean
;                           If set, all parameters will be passed to the Make_Array
;                               function in order to create the array. This is the default
;                               if `DATA` is a scalar.
;       NORMAL:             in, optional, type=boolean, default=0
;                           If set, then `SEED` will be input to the RandomN function
;                               instead of the RandomU function.
;       NO_CLOBBER:         in, optional, type=boolean, default=0
;                           If set and `NAME` already exists in the cache, then append
;                               '_#' to `NAME` to avoid conflict. "_#" is the smallest
;                               integer that results in a unique name. The default is to
;                               clobber cached variables with the same name.
;       NO_COPY:            in, optional, type=boolean, default=0
;                           If set, `DATA` will be copied directly into the object and
;                               will be left undefined.
;       RANDOM:             in, optional, type=boolean, default=0
;                           If set, the RandomU function will be used to initialize the
;                               array instead of Make_Array. `SEED`, then, is the seed,
;                               and `DATA`, `D2`, etc. specify the dimensions of the
;                               resulting array. See IDL's RandomU function for more
;                               details.
;       SEED:               in, out, optional, type=boolean, default=long
;                           First parameter to the RandomU function. See the `RANDOM`
;                               keyword. If a value is given, `RANDOM` is automatically
;                               set to 1.
;       TO_ARRAY:           in, optional, type=boolean, default=0
;                           If set, `DATA` is a structure, list, or hash whose contents
;                               are to be converted to an array.
;       TYPE:               in, optional, type=int/string, default='FLOAT'
;                           The type of array to be created. Used with `MAKE`.
;       VERBOSE:            in, optional, type=byte, default=1
;                           Level of verboseness of printed messages. Options are::
;                               0 - Quiet
;                               1 - Less verbose
;                               2 - More verbose
;       _REF_EXTRA:         in, optional, type=any
;                           Any keyword accepted by IDL's Make_Array, RandomU, or RandomN
;                               function, or by the ToArray method, depending on the
;                               keywords used.
;-
function MrVariable::INIT, data, D2, D3, D4, D5, D6, D7, D8, $
CACHE=cache, $
MAKE=make, $
NAME=name, $
NO_CLOBBER=no_clobber, $
NORMAL=normal, $
NO_COPY=no_copy, $
RANDOM=random, $
SEED=seed, $
TO_ARRAY=to_array, $
TYPE=type, $
VERBOSE=verbose, $
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
	tf_cache = keyword_set(cache)
	if n_elements(name)    eq 0 then name    = 'MrVariable'
	if n_elements(verbose) eq 0 then verbose = 1B

	;Allocate heap to pointers.
	self.attributes       = hash()
	self.data             = ptr_new(/ALLOCATE_HEAP)
	self._append_unReform = ptr_new(/ALLOCATE_HEAP)
	self._append_unTrans  = ptr_new(/ALLOCATE_HEAP)

	;Create the array
	make     = keyword_set(make)
	normal   = keyword_set(normal)
	random   = keyword_set(random)
	to_array = keyword_set(to_array)
	if n_elements(seed) gt 0 then random = ~normal
	if n_params() gt 1 && ~normal && ~random && ~to_array then make = 1
	if make + normal + random + to_array gt 1 then message, 'MAKE, NORMAL, RANDOM and TO_ARRAY are mutually exclusive.'

	;Were inputs given?
	;   `DATA` can be undefined if RANDOM or NORMAL are set.
	case 1 of
		make:     self -> Make_Array, data, D2, D3, D4, D5, D6, D7, D8, TYPE=type, _STRICT_EXTRA=extra
		normal:   self -> RandomN, seed, data, D2, D3, D4, D5, D6, D7, D8, _STRICT_EXTRA=extra
		random:   self -> RandomU, seed, data, D2, D3, D4, D5, D6, D7, D8, _STRICT_EXTRA=extra
		to_array: self -> ToArray, data, _STRICT_EXTRA=extra
		n_elements(data) gt 0: begin
			if n_elements(extra) gt 0 then message, 'EXTRA keywords not allowed when DATA is an array.'
			self -> SetData, data, NO_COPY=no_copy
		endcase
		else: ;Leave the array empty
	endcase

	;Set properties
	self -> SetName, name
	self -> SetProperty, VERBOSE=verbose
	
	;Cache the array
	if tf_cache then self -> Cache, NO_CLOBBER=no_clobber

	return, 1
end


;+
;   Clean up after the object is destroyed
;-
pro MrVariable::CLEANUP
	compile_opt idl2
	on_error, 2
	
	;Common block
	@mrvar_common
	
	;Remove from the variable container
	if n_elements(MrVarCache) gt 0 then begin
		if MrVarCache -> IsContained(self) then MrVarCache -> Remove, self
	endif
	
	;Destroy objects
	obj_destroy, self.attributes

	;Free pointers
	ptr_free, self.data
	ptr_free, self._append_unReform
	ptr_free, self._append_unTrans
end


;+
;   Helper method to determine if "self" was given on the left or right side of an
;   operator.
;
; :Private:
;
; :Params:
;       LEFT:               out, required, type=any
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;       RIGHT:              out, required, type=long
;                           The argument that appears on the left side of the operator. 
;                               Possibly the implicit "self".
;
; :Keywords:
;       ISMRVARIABLE:       out, optional, type=boolean
;                           Returns true (1) if the non-"self" parameter is also a
;                               MrVariable object and false (0) otherwise.
;
; :Returns:
;       SIDE:               Returns 'RIGHT' if "self" was provided on the right side of
;                               the operator and 'LEFT' if it was provided on the left.
;-
function MrVariable::_LeftOrRight, left, right, $
ISMRVARIABLE=isMrVariable
	compile_opt idl2
	on_error, 2

	;Information about inputs
	lType = size(left, /TNAME)
	rType = size(right, /TNAME)
	heapID = obj_valid(self, /GET_HEAP_IDENTIFIER)

	;Was "self" given on the left of operator?
	side = ''
	if lType eq 'OBJREF' then begin
		lHeapID = obj_valid(left, /GET_HEAP_IDENTIFIER)
		if lHeapID eq heapID then side = 'LEFT'
	endif

	;Was "self" given on the right of the AND operator? Treat the special case
	;where "self" is given on both sides.
	if rType eq 'OBJREF' then begin
		rHeapID = obj_valid(right, /GET_HEAP_IDENTIFIER)
		if rHeapID eq heapID then if side ne 'LEFT' then side = 'RIGHT'
	endif

	;Is the other argument a MrVariable object as well?
	if arg_present(isMrVariable) then begin
		case side of
			'LEFT':  isMrVariable = rtype eq 'OBJREF' && obj_isa(right, 'MrVariable')
			'RIGHT': isMrVariable = ltype eq 'OBJREF' && obj_isa(left,  'MrVariable')
		endcase
	endif

	return, side
end


;+
;   The purpose of this method is to perform a bit-wise comparison between two
;   integer, longword or byte expressions. For other types, `RIGHT` is returned unless
;   `LEFT` is zero or the empty string, in which case 0 (zero) is returned.
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
function MrVariable::_OverloadAnd, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Add
		result = (*self.data) and (right -> GetData())
		
		;Create a new name
		name = 'AND(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) and right['DATA'] $
			else result = left['DATA'] and (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'AND(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'AND(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to multiply two expressions together.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression, normal IDL operations
;     will take effect. This means the output will be the size of the
;     array with the smallest dimensions.
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
function MrVariable::_OverloadAsterisk, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects //////////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Multiply
		result = *self.data * (right -> GetData())

		;Create a new name
		name = 'Multiply(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Multiply the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) * right $
			else result = left * (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Multiply(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Multiply(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to apply exponentiation with the caret operator.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression, normal IDL operations
;     will take effect. This means the output will be the size of the
;     array with the smallest dimensions.
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
;       RESULT:             The result of raising `LEFT` to the power of `RIGHT`.
;-
function MrVariable::_OverloadCaret, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects //////////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Caret
		result = (*self.data) ^ (right -> GetData())
		
		;Create a new name
		name = 'Caret(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression //////////////////////////////
;-------------------------------------------------------
	endif else begin
		;Multiply the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) ^ right $
			else result = left ^ (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Caret(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Caret(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   Allow square-bracket array indexing from the left side of an operator.
;
; Notes:
;   This method works, but was replaced by the other _OverloadBracketsLeftSide method.
;   It is my own attempt at bracket overloading.
;
; :Params:
;       OBJREF:             in, required, type=ObjRef
;                           The object reference variable that is being indexed (i.e. "self")
;                               Use when you want to replace "self" with `VALUE`.
;       VALUE:              in, required, type=numeric array
;                           The value specified on the right-hand side of the equal sign.
;       ISRANGE:            in, required, type=intarr
;                           A vector that has one element for each Subscript argument
;                               supplied by the user; each element contains a zero if the
;                               corresponding input argument was a scalar index value or
;                               array of indices, or a one if the corresponding input
;                               argument was a subscript range.
;       SUBSCRIPT1:         in, required, type=integer/intarr(3)
;                           Index subscripts. Either a scalar, an index array, or a 
;                               subscript range in the form [start, stop, step_size]
;       SUBSCRIPT2:         in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       SUBSCRIPT3:         in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       SUBSCRIPT4:         in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       SUBSCRIPT5:         in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       SUBSCRIPT6:         in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       SUBSCRIPT7:         in, optional, type=integer/intarr(3)
;                           Index subscripts.
;       SUBSCRIPT8:         in, optional, type=integer/intarr(3)
;                           Index subscripts.
;-
pro MrVariable::_OverloadBLS, objRef, value, isRange, subscript1, subscript2, $
                                        subscript3, subscript4, subscript5, $
                                        subscript6, subscript7, subscript8
	compile_opt idl2
	on_error, 2

	;Reform VALUE into a 1D array
	nValues = n_elements(value)
	valueOut = reform(value, nValues)

	;Get the implicit array dimensions and reform the subscripts into a 1D array
	;of indices.
	self -> GetProperty, DIMENSIONS=dims
	indices = MrReformIndices(dims, isRange, subscript1, subscript2, subscript3, $
	                                         subscript4, subscript5, subscript6, $
	                                         subscript7, subscript8)

	;Set the values of the implicit array
	(*self.data)[indices] = valueOut
end


;+
;   Allow square-bracket array indexing from the left side of an operator.
;
; :Examples:
;   Saving a subarray into a 4D array::
;       myArray = MrVariable(4, 2, 5, 3, TYPE='FLOAT')
;       myArray[1:2, 0, 2:4, 1] = RandomU(5, 2,1,3,1)
;
;   Saving a subarray into a 3D array::
;       myArray = MrVariable(fltarr(5,5,5))
;       myArray[*,-3,0:2] = findgen(5,1,3)
;
;
; :Params:
;       OBJREF:             in, required, type=ObjRef
;                           The object reference variable that is being indexed (i.e. "self")
;                               Use when you want to replace "self" with `VALUE`.
;       VALUE:              in, required, type=numeric array
;                           The value specified on the right-hand side of the equal sign.
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
;-
PRO MrVariable::_OverloadBracketsLeftSide, objRef, value, isRange, i1, i2, i3, i4, i5, i6, i7, i8
	Compile_Opt idl2
	On_Error, 2

	;Number OF subscripts given
	nSubscripts = N_Elements(isRange)

;---------------------------------------------------------------------
; Attribute Name Given ///////////////////////////////////////////////
;---------------------------------------------------------------------
	IF nSubscripts EQ 1 && MrIsA(i1, 'STRING', /SCALAR) THEN BEGIN
		IF i1 EQ 'DATA' $
			THEN self -> SetData, value $
			ELSE self -> SetAttributeValue, i1, value, /CREATE
		
		;Done
		RETURN
	ENDIF

;---------------------------------------------------------------------
; Brute Force Subscripting ///////////////////////////////////////////
;---------------------------------------------------------------------
	;
	; Highly optimized code for three dimensions or lower. 
	; Handle all combinations of subscript ranges or indices.
	;
	
	;Was the data given a MrVariable object?
	tf_mrvar = IsA(value, 'MrVariable')

	;<= 3D subscript range
	IF (nSubscripts LE 3) THEN BEGIN
	;---------------------------------------------------------------------
	; 3D Subscripts //////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		IF IsA(i3) THEN BEGIN 
			;[? ,? , min:max:interval]
			IF isRange[2] THEN BEGIN 
				;[? , min:max:interval, min:max:interval]
				IF isRange[1] THEN BEGIN 
					;[min:max:interval , min:max:interval, min:max:interval]
					IF isRange[0] THEN BEGIN
						(*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2],i3[0]:i3[1]:i3[2]] = tf_mrvar ? value['DATA'] : value
					;[index , min:max:interval, min:max:interval]
					ENDIF ELSE BEGIN
						(*self.data)[i1,i2[0]:i2[1]:i2[2],i3[0]:i3[1]:i3[2]] = tf_mrvar ? value['DATA'] : value
					ENDELSE
				;[? , index, min:max:interval]
				ENDIF ELSE BEGIN
					;[min:max:interval , index, min:max:interval]
					IF isRange[0] THEN BEGIN 
						(*self.data)[i1[0]:i1[1]:i1[2],i2,i3[0]:i3[1]:i3[2]] = tf_mrvar ? value['DATA'] : value
					;[index , index, min:max:interval]
					ENDIF ELSE BEGIN 
						(*self.data)[i1,i2,i3[0]:i3[1]:i3[2]] = tf_mrvar ? value['DATA'] : value
					ENDELSE 
				ENDELSE
			
			;[? , ?, index]
			ENDIF ELSE BEGIN
				;[? , min:max:interval, index]
				IF isRange[1] THEN BEGIN
					;[min:max:interval, min:max:interval, index]
					IF isRange[0] THEN BEGIN 
						(*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2],i3] = tf_mrvar ? value['DATA'] : value
					;[index, min:max:interval, index]
					ENDIF ELSE BEGIN 
						(*self.data)[i1,i2[0]:i2[1]:i2[2],i3] = tf_mrvar ? value['DATA'] : value
					ENDELSE 
				;[?, index, index]
				ENDIF ELSE BEGIN
					;[min:max:interval, index, index]
					IF isRange[0] THEN BEGIN 
						(*self.data)[i1[0]:i1[1]:i1[2],i2,i3] = tf_mrvar ? value['DATA'] : value
					;[index, index, index]
					ENDIF ELSE BEGIN 
						(*self.data)[i1,i2,i3] = tf_mrvar ? value['DATA'] : value
					ENDELSE 
				ENDELSE 
			ENDELSE 
	
	;---------------------------------------------------------------------
	; 2D Subscripts //////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE IF IsA(i2) THEN BEGIN
			;[?, min:max:interval]
			IF isRange[1] THEN BEGIN
				;[min:max:interval, min:max:interval]
				IF isRange[0] THEN BEGIN
					(*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2]] = tf_mrvar ? value['DATA'] : value
				;[index, min:max:interval]
				ENDIF ELSE BEGIN
					(*self.data)[i1,i2[0]:i2[1]:i2[2]] = tf_mrvar ? value['DATA'] : value
				ENDELSE
			;[?, index]
			ENDIF ELSE BEGIN
				;[min:max:interval, index]
				IF isRange[0] THEN BEGIN
					(*self.data)[i1[0]:i1[1]:i1[2],i2] = tf_mrvar ? value['DATA'] : value
				;[index, index]
				ENDIF ELSE BEGIN
					(*self.data)[i1,i2] = tf_mrvar ? value['DATA'] : value
				ENDELSE
			ENDELSE
	
	;---------------------------------------------------------------------
	; 1D Subscripts //////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE BEGIN
			;min:max:interval
			IF isRange[0] THEN BEGIN
				(*self.data)[i1[0]:i1[1]:i1[2]] = tf_mrvar ? value['DATA'] : value
			;Index
			ENDIF ELSE BEGIN
				(*self.data)[i1] = tf_mrvar ? value['DATA'] : value
			ENDELSE
		ENDELSE

;---------------------------------------------------------------------
; >= 4D Subscripts ///////////////////////////////////////////////////
;---------------------------------------------------------------------
	ENDIF ELSE BEGIN
		;Works FOR any number of dimensions.
		indices = MrReformIndices( Size(*self.data, /DIMENSIONS), isRange, i1, i2, i3, i4, i5, i6, i7, i8)
		(*self.data)[indices] = tf_mrvar ? value['DATA'] : value
	ENDELSE
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
FUNCTION MrVariable::_OverloadBracketsRightSide_IDLVar, isRange, i1, i2, i3, i4, i5, i6, i7, i8
	Compile_Opt idl2
	On_Error, 2

	;Number of subscripts given
	nSubscripts = N_Elements(isRange)

	;String operations.
;	IF IsA(i1, /SCALAR, 'STRING') THEN BEGIN
;		case i1 of
;			'DATA':    RETURN, *self.data
;			'POINTER': RETURN,  self.data
;			'PTR':     RETURN,  self.data
;			ELSE:      RETURN,  self -> GetAttrValue(i1)
;		endcase

	;Scalar operations
	;   - 0   returns the self object
	;   - [0] returns the first data element
	;   - All other cases RETURN data
;	ENDIF ELSE IF nSubscripts eq 1 && isRange[0] eq 0 && IsA(i1, /SCALAR) && i1 eq 0 THEN BEGIN
;		RETURN, self
;	ENDIF

;---------------------------------------------------------------------
;Optimized Subscripting for <= 3D ////////////////////////////////////
;---------------------------------------------------------------------
	IF (nSubscripts le 3) THEN BEGIN
	;---------------------------------------------------------------------
	;3D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		IF IsA(i3) THEN BEGIN 
			;Subscript range given: [min:max:interval]?
			IF isRange[2] THEN BEGIN 
				;Subscript range for dimensions 2 and 3?
				IF isRange[1] THEN BEGIN 
					;Range: [1,2,3], Index: --
					IF isRange[0] THEN BEGIN
						RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2],i3[0]:i3[1]:i3[2]]
					;Range: [2,3], Index: 1
					ENDIF ELSE BEGIN
						RETURN, (*self.data)[i1,i2[0]:i2[1]:i2[2],i3[0]:i3[1]:i3[2]] 
					ENDELSE
				;Index value for dimension 2?
				ENDIF ELSE BEGIN
					;Range: [3,1], Index: 2
					IF isRange[0] THEN BEGIN 
						RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2,i3[0]:i3[1]:i3[2]]
					;Range: 3, Index: [2,1]
					ENDIF ELSE BEGIN 
						RETURN, (*self.data)[i1,i2,i3[0]:i3[1]:i3[2]] 
					ENDELSE 
				ENDELSE
			;Index for dimension 3?
			ENDIF ELSE BEGIN
				;Range for dimension 2?
				IF isRange[1] THEN BEGIN
					;Range: [2,1]
					IF isRange[0] THEN BEGIN 
						RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2],i3] 
					ENDIF ELSE BEGIN 
						RETURN, (*self.data)[i1,i2[0]:i2[1]:i2[2],i3] 
					ENDELSE 
				ENDIF ELSE BEGIN 
					IF isRange[0] THEN BEGIN 
						RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2,i3] 
					ENDIF ELSE BEGIN 
						RETURN, (*self.data)[i1,i2,i3] 
					ENDELSE 
				ENDELSE 
			ENDELSE 
	;---------------------------------------------------------------------
	;2D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE IF IsA(i2) THEN BEGIN
			IF isRange[1] THEN BEGIN
				;[Range, Range]
				IF isRange[0] THEN BEGIN
					RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2]]
				;[Index, Range]
				ENDIF ELSE BEGIN
					RETURN, (*self.data)[i1,i2[0]:i2[1]:i2[2]]
				ENDELSE
			ENDIF ELSE BEGIN
				;[Range, Index]
				IF isRange[0] THEN BEGIN
					RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2]
				;[Index, Index]
				ENDIF ELSE BEGIN
					RETURN, (*self.data)[i1,i2]
				ENDELSE
			ENDELSE
	;---------------------------------------------------------------------
	;1D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE BEGIN
			;Range?
			IF isRange[0] THEN BEGIN
				RETURN, (*self.data)[i1[0]:i1[1]:i1[2]] 
			;Index
			;   - Compensate for passing in [0] instead of 0 for the first element.
			;   - i.e. RETURN a scalar instead of a 1-element array
			ENDIF ELSE BEGIN
				IF N_Elements(i1) eq 1 && i1 eq 0 $
					THEN RETURN, (*self.data)[0] $
					ELSE RETURN, (*self.data)[i1]
			ENDELSE
		ENDELSE
	ENDIF

;---------------------------------------------------------------------
; Brute Force Code for 4D or Higher Arrays. //////////////////////////
;---------------------------------------------------------------------
	;Works for any number of dimensions.
	dims = Size(*self.data, /DIMENSIONS)
	indices = MrReformIndices(dims, isRange, i1, i2, i3, i4, i5, i6, i7, i8, DIMENSIONS=dimensions);, /IDL_METHOD)

	RETURN, Reform((*self.data)[indices], dimensions, /OVERWRITE)
END


;+
;   Allow square-bracket array indexing from the right side of an operator.
;
;   Calling Sequence
;       oTSvar    = oTSvar[0]
;       t0        = oTSvar[[0]]
;       data      = oTSvar['DATA' [, i2 [, i3 [, i4 [, i5 [, i6 [, i7 [, i8]]]]]]]
;       pData     = oTSvar['POINTER']
;       pData     = oTSvar['PTR']
;       attrValue = oTSvar[AttrName]
;       oVar      = oTSvar[i1 [, i2 [, i3 [, i4 [, i5 [, i6 [, i7 [, i8]]]]]]]]
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
;       RESULT:             out, required, type=MrVariable objref
;                           The subarray accessed by the input parameters.
;-
FUNCTION MrVariable::_OverloadBracketsRightSide_Complete, isRange, i1, i2, i3, i4, i5, i6, i7, i8
	Compile_Opt idl2
	On_Error, 2

	;Number of subscripts given
	nSubs = N_Elements(isRange)

	;String operations.
	IF IsA(i1, /SCALAR, 'STRING') THEN BEGIN
		CASE i1 OF
			'DATA': BEGIN
				;Extract the subset of data
				CASE nSubs OF
					1: RETURN, *self.data
					2: RETURN, self -> _OverloadBracketsRightSide_IDLVar(isRange[1], i2)
					3: RETURN, self -> _OverloadBracketsRightSide_IDLVar(isRange[1:2], i2, i3)
					4: RETURN, self -> _OverloadBracketsRightSide_IDLVar(isRange[1:3], i2, i3, i4)
					5: RETURN, self -> _OverloadBracketsRightSide_IDLVar(isRange[1:4], i2, i3, i4, i5)
					6: RETURN, self -> _OverloadBracketsRightSide_IDLVar(isRange[1:5], i2, i3, i4, i5, i6)
					7: RETURN, self -> _OverloadBracketsRightSide_IDLVar(isRange[1:6], i2, i3, i4, i5, i6, i7)
					8: RETURN, self -> _OverloadBracketsRightSide_IDLVar(isRange[1:7], i2, i3, i4, i5, i6, i7, i8)
					ELSE: Message, 'Incorrect number of subscripts given (' + String(nSubs, FORMAT='(i0)') + ').'
				ENDCASE
			ENDCASE
			'POINTER': RETURN, self.data
			'PTR':     RETURN, self.data
			ELSE:      RETURN, self -> GetAttrValue(i1)
		ENDCASE

	;Scalar operations
	;   - 0   returns the self object
	;   - [0] returns the first data element
	;   - All other cases RETURN data
	ENDIF ELSE IF nSubs EQ 1 && isRange[0] EQ 0 && IsA(i1, /SCALAR) && i1 EQ 0 THEN BEGIN
		RETURN, self
	ENDIF

;---------------------------------------------------------------------
;Optimized Subscripting for <= 3D ////////////////////////////////////
;---------------------------------------------------------------------
	IF (nSubs LE 3) THEN BEGIN
	;---------------------------------------------------------------------
	;3D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		IF IsA(i3) THEN BEGIN 
			;Subscript range given: [min:max:interval]?
			IF isRange[2] THEN BEGIN 
				;Subscript range for dimensions 2 and 3?
				IF isRange[1] THEN BEGIN 
					;Range: [1,2,3], Index: --
					IF isRange[0] THEN BEGIN
						data = (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2],i3[0]:i3[1]:i3[2]]
					;Range: [2,3], Index: 1
					ENDIF ELSE BEGIN
						data = (*self.data)[i1,i2[0]:i2[1]:i2[2],i3[0]:i3[1]:i3[2]] 
					ENDELSE
				;Index value for dimension 2?
				ENDIF ELSE BEGIN
					;Range: [3,1], Index: 2
					IF isRange[0] THEN BEGIN 
						data = (*self.data)[i1[0]:i1[1]:i1[2],i2,i3[0]:i3[1]:i3[2]]
					;Range: 3, Index: [2,1]
					ENDIF ELSE BEGIN 
						data = (*self.data)[i1,i2,i3[0]:i3[1]:i3[2]] 
					ENDELSE 
				ENDELSE
			;Index for dimension 3?
			ENDIF ELSE BEGIN
				;Range for dimension 2?
				IF isRange[1] THEN BEGIN
					;Range: [2,1]
					IF isRange[0] THEN BEGIN 
						data = (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2],i3] 
					ENDIF ELSE BEGIN 
						data = (*self.data)[i1,i2[0]:i2[1]:i2[2],i3] 
					ENDELSE 
				ENDIF ELSE BEGIN 
					IF isRange[0] THEN BEGIN 
						data = (*self.data)[i1[0]:i1[1]:i1[2],i2,i3] 
					ENDIF ELSE BEGIN 
						data = (*self.data)[i1,i2,i3] 
					ENDELSE 
				ENDELSE 
			ENDELSE 
	;---------------------------------------------------------------------
	;2D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE IF IsA(i2) THEN BEGIN
			IF isRange[1] THEN BEGIN
				;[Range, Range]
				IF isRange[0] THEN BEGIN
					data = (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2]]
				;[Index, Range]
				ENDIF ELSE BEGIN
					data = (*self.data)[i1,i2[0]:i2[1]:i2[2]]
				ENDELSE
			ENDIF ELSE BEGIN
				;[Range, Index]
				IF isRange[0] THEN BEGIN
					data = (*self.data)[i1[0]:i1[1]:i1[2],i2]
				;[Index, Index]
				ENDIF ELSE BEGIN
					data = (*self.data)[i1,i2]
				ENDELSE
			ENDELSE
	;---------------------------------------------------------------------
	;1D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE BEGIN
			;Range?
			IF isRange[0] THEN BEGIN
				data = (*self.data)[i1[0]:i1[1]:i1[2]] 
			;Index
			;   - Compensate for passing in [0] instead OF 0 FOR the first element.
			;   - i.e. RETURN a scalar instead OF a 1-element array
			ENDIF ELSE BEGIN
				IF N_Elements(i1) EQ 1 && i1 EQ 0 $
					THEN data = (*self.data)[0] $
					ELSE data = (*self.data)[i1]
			ENDELSE
		ENDELSE
		
;---------------------------------------------------------------------
; Brute Force Code FOR >= 4D /////////////////////////////////////////
;---------------------------------------------------------------------
	ENDIF ELSE BEGIN
		dims    = Size(*self.data, /DIMENSIONS)
		indices = MrReformIndices(dims, isRange, i1, i2, i3, i4, i5, i6, i7, i8, DIMENSIONS=dimensions)
		data    = Reform((*self.data)[indices], dimensions, /OVERWRITE)
	ENDELSE
	
;-------------------------------------------
; Create Output Array //////////////////////
;-------------------------------------------
	;Create the variable
	vOut = MrVariable( data, $
	                   NAME='OverloadBRS(' + self.name + ')', $
	                   /NO_COPY )
	
	;Copy all attributes
	self -> CopyAttrTo, vOut

;-------------------------------------------
; DEPEND_# /////////////////////////////////
;-------------------------------------------
	;
	; Rules:
	;   - NSUBS must be equal to # dependent variables, otherwise subscripting is indeterminate
	;   - Dimensions of implicit array are ordered [DEPEND_0, DEPEND_1, ..., DEPEND_N]
	;   - DEPEND_[1-N] may have record variance, in which case they are 2D
	;   - If a bulk process currently acting on the implicit, and then it may already
	;     have affected the (cached) DEPEND_[1-N] variable. If DEPEND_[1-N] is cached,
	;     allow it to have the same dimensions as the output variable.
	;
	
	;
	; TODO:
	;   Ideally, this should be done with all poitner variables:
	;   https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html
	;       FORM_PTR
	;       LABL_PTR_[1-3]
	;       SCAL_PTR
	;       UNITS_PTR
	;
	
	;Size of implicit & sub- array
	;   - Original size, not the size of the subarray
	inDims   = Size(*self.data, /DIMENSIONS)
	outDims  = Size(vOut, /DIMENSIONS)
	nInDims  = Size(*self.data, /N_DIMENSIONS)
	nOutDims = Size(vOut, /N_DIMENSIONS)
	
	;Number of dependent variables
	;   - IDL allows up to 8 dimensions
	nDeps = ~self -> HasAttr('DEPEND_0')   ? 0 $
	        : ~self -> HasAttr('DEPEND_1') ? 1 $
	        : ~self -> HasAttr('DEPEND_2') ? 2 $
	        : ~self -> HasAttr('DEPEND_3') ? 3 $
	        : ~self -> HasAttr('DEPEND_4') ? 4 $
	        : ~self -> HasAttr('DEPEND_5') ? 5 $
	        : ~self -> HasAttr('DEPEND_6') ? 6 $
	        : ~self -> HasAttr('DEPEND_7') ? 7 $
	        : 8
	
	;RETURN if nSubs LE nDeps
	;   - There does not have to be a DEPEND_# for each subscript
	;   - There DOES     have to be a subscript for each DEPEND_#
	IF nSubs LT nDeps THEN BEGIN
		strNSubs = String(nSubs, FORMAT='(i1)')
		strNDeps = String(nDeps, FORMAT='(i1)')
		IF nDeps GT 0 THEN MrPrintF, 'LogWarn', '# subscripts (' + strNSubs + ') LT # dependent variables (' + strNDeps + ').'
		RETURN, vOut
	ENDIF
	
	;Step through each dependent variable
	FOR i = 0, nDeps-1 DO BEGIN
		depend = 'DEPEND_' + String(i, FORMAT='(i0)')
		IF ~self -> HasAttr(depend) THEN CONTINUE
		
		;Get the dependent variable and its dimensions
		oDep     = self[depend]
		nDepDims = Size(oDep, /N_DIMENSIONS)
		depDims  = Size(oDep, /DIMENSIONS)
		
		;RECORD VARIANCE
		;   - Data     = [N, M, L, P]
		;   - DEPEND_0 = [N]
		;   - DEPEND_1 = [N, M]
		;   - DEPEND_2 = [N, L]
		;   - DEPEND_3 = [N, P]
		IF depend NE 'DEPEND_0' && nDepDims EQ 2 THEN BEGIN
			;Apply to DEPEND_#
			;   - Dimensions of implicit and DEPEND_# arrays match
			IF Array_Equal(depDims, inDims[[0,i]]) THEN BEGIN
				CASE i OF
					1: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i2 )
					2: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i3 )
					3: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i4 )
					4: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i5 )
					5: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i6 )
					6: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i7 )
					7: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i8 )
					ELSE: Message, 'Invalid number of dependent variables.'
				ENDCASE
			
			;Pre-Applied
			;   - It is possible that whatever is happening to the implicit variable
			;     is happening in bulk (e.g. MrVar_TLimit).
			;   - If DEPEND_# is in the cache, allow it to have the same dimensions
			;     as the output variable
			ENDIF ELSE IF MrVar_IsCached(oDep) && Array_Equal(depDims, outDims[[0,i]]) THEN BEGIN
				oDepend = oDep
			
			;UNKNOWN DIMENSION SIZES
			;   - Warning with Expected Dimensions
			ENDIF ELSE BEGIN
				strDims = '[' + StrJoin(String(inDims[[0,i]], FORMAT='(i0)'), ', ') + ']'
			ENDELSE
			
		;NO RECORD VARIANCE
		;   - Data     = [N, M, L, P]
		;   - DEPEND_0 = [N]
		;   - DEPEND_1 = [M]
		;   - DEPEND_2 = [L]
		;   - DEPEND_3 = [P]
		ENDIF ELSE IF nDepDims EQ 1 THEN BEGIN
			;Apply to DEPEND_#
			;   - Dimensions of implicit and DEPEND_# arrays match
			;   - Prevent scalar 0 from returning object by including "[]"
			IF depDims EQ inDims[i] THEN BEGIN
				CASE i OF
					0: oDepend = oDep -> _OverloadBracketsRightSide( isRange[0], [i1] )
					1: oDepend = oDep -> _OverloadBracketsRightSide( isRange[1], [i2] )
					2: oDepend = oDep -> _OverloadBracketsRightSide( isRange[2], [i3] )
					3: oDepend = oDep -> _OverloadBracketsRightSide( isRange[3], [i4] )
					4: oDepend = oDep -> _OverloadBracketsRightSide( isRange[4], [i5] )
					5: oDepend = oDep -> _OverloadBracketsRightSide( isRange[5], [i6] )
					6: oDepend = oDep -> _OverloadBracketsRightSide( isRange[6], [i7] )
					7: oDepend = oDep -> _OverloadBracketsRightSide( isRange[7], [i8] )
					ELSE: Message, 'Invalid number of dependent variables.'
				ENDCASE
			
			;Pre-Applied
			;   - It is possible that whatever is happening to the implicit variable
			;     is happening in bulk (e.g. MrVar_TLimit).
			;   - If DEPEND_# is in the cache, allow it to have the same dimensions
			;     as the output variable
			ENDIF ELSE IF MrVar_IsCached(oDep) && depDims EQ outDims[i] THEN BEGIN
				oDepend = oDep
			
			;UNKNOWN DIMENSION SIZES
			;   - Warning with Expected Dimensions
			ENDIF ELSE BEGIN
				strDims = '[' + StrJoin(String(inDims[i], FORMAT='(i0)'), ', ') + ']'
			ENDELSE
		
		;UNKNOWN DIMENSION SIZES
		;   - Warning with Expected Dimensions
		ENDIF ELSE BEGIN
			strDims = '[' + StrJoin(String(inDims[i], FORMAT='(i0)'), ', ') + ']'
		ENDELSE
		
		;Set the DEPEND_# attribute
		IF Obj_Valid(oDepend) THEN BEGIN
			vOut[depend] = oDepend
		
		;Issue warning
		ENDIF ELSE BEGIN
			;Dimensions
			strInDims  = '[' + StrJoin(String(inDims,    FORMAT='(i0)'), ', ') + ']'
			strOutDims = '[' + StrJoin(String(outDims,   FORMAT='(i0)'), ', ') + ']'
			strDepDims = '[' + StrJoin(String(depDims,   FORMAT='(i0)'), ', ') + ']'
		
			;Warning information
			MrPrintF, 'LogWarn', 'Cannot extract ' + depend + ' data for variable "' + self.name + '".'
			MrPrintF, 'LogWarn', '    ' + depend + ' dims:  ' + strDepDims + '.'
			MrPrintF, 'LogWarn', '    Data dims:      ' + strInDims
			MrPrintF, 'LogWarn', '    Subarray dims:  ' + strOutDims
			MrPrintF, 'LogWarn', '    Expected dims:  ' + strDims
		ENDELSE
		
		;Invalidate the DEPEND_# variable
		oDepend = 0B
	ENDFOR

;-------------------------------------------
; DELTA_(PLUS|MINUS)_VAR ///////////////////
;-------------------------------------------
	deltas  = ['DELTA_MINUS_VAR', 'DELTA_PLUS_VAR']
	nDeltas = N_Elements(deltas)
	
	;Step through each dependent variable
	FOR i = 0, nDeltas-1 DO BEGIN
		IF ~self -> HasAttr(deltas[i]) THEN CONTINUE
		
		;Get the dependent variable and its dimensions
		oDel       = self[deltas[i]]
		nDelDims = Size(oDel, /N_DIMENSIONS)
		delDims  = Size(oDel, /DIMENSIONS)
		
		;SAME
		;   - Both record varying
		;   - Both non-record varying
		IF nDelDims EQ nInDims && Array_Equal(delDims, inDims) THEN BEGIN
			CASE nDelDims OF
				1: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0],   i1)
				2: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:1], i1, i2)
				3: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:2], i1, i2, i3)
				4: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:3], i1, i2, i3, i4)
				5: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:4], i1, i2, i3, i4, i5)
				6: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:5], i1, i2, i3, i4, i5, i6)
				7: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:6], i1, i2, i3, i4, i5, i6, i7)
				8: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:7], i1, i2, i3, i4, i5, i6, i7, i8)
				ELSE: Message, 'Incorrect number of dimensions (' + String(nDelDims, FORMAT='(i0)') + ').'
			ENDCASE
		
		;DIFFERENT
		;   - Parent is record varying
		;   - DELTA is non-record varying
		ENDIF ELSE IF nDelDims EQ nInDims-1 && Array_Equal(delDims, inDims[1:*]) THEN BEGIN
			CASE nDelDims OF
				1: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1],   i2)
				2: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:2], i2, i3)
				3: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:3], i2, i3, i4)
				4: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:4], i2, i3, i4, i5)
				5: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:5], i2, i3, i4, i5, i6)
				6: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:6], i2, i3, i4, i5, i6, i7)
				7: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:7], i2, i3, i4, i5, i6, i7, i8)
				ELSE: Message, 'Incorrect number of dimensions (' + String(nDelDims, FORMAT='(i0)') + ').'
			ENDCASE
		
		;INCOMPATIBLE
		ENDIF ELSE BEGIN
			;Dimensions
			strInDims  = '[' + StrJoin(String(inDims,    FORMAT='(i0)'), ', ') + ']'
			strOutDims = '[' + StrJoin(String(outDims,   FORMAT='(i0)'), ', ') + ']'
			strDims    = '[' + StrJoin(String(inDims[0], FORMAT='(i0)'), ', ') + ']'
			strDelDims = '[' + StrJoin(String(delDims,   FORMAT='(i0)'), ', ') + ']'
		
			;Warning information
			MrPrintF, 'LogWarn', 'Cannot extract ' + deltas[i] + ' data for variable "' + self.name + '".'
			MrPrintF, 'LogWarn', '    ' + deltas[i] + ' dims:  ' + strDelDims + '.'
			MrPrintF, 'LogWarn', '    Data dims:             ' + strInDims
			MrPrintF, 'LogWarn', '    Subarray dims:         ' + strOutDims
			
			;Invalidate the DEPEND_# variable
			oDelta = 0B
		ENDELSE
		
		;Create a variable
		IF Obj_Valid(oDelta) THEN vOut[deltas[i]] = oDelta
	ENDFOR
	
;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------

	RETURN, vOut
END


;+
;   Apply square-bracket array indexing to variable attributes. Relevant attributes are:
;       DEPEND_[0-7]
;       DELTA_(PLUS|MINUS)_VAR
;
; :Params:
;       OUTDIMS:            in, required, type=integer/intarr
;                           The sizes of each dimensions of the subarray after array
;                               indexing has been applied.
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
;       ATTRS:              out, required, type=hash
;                           Attribute names and values after indexing.
;-
FUNCTION MrVariable::_OverloadBracketsRightSide_Attrs, outDims, isRange, i1, i2, i3, i4, i5, i6, i7, i8
	Compile_Opt idl2
	On_Error, 2
	
	;Output variable
	attrs = Hash()
	
	;Number of subscripts given
	nSubs = N_Elements(isRange)
	
	;Size of implicit & sub- array
	;   - Original size, not the size of the subarray
	inDims   = Size(*self.data, /DIMENSIONS)
	nInDims  = Size(*self.data, /N_DIMENSIONS)
	nOutDims = Size(outDims,    /N_DIMENSIONS)

;-------------------------------------------
; DEPEND_# /////////////////////////////////
;-------------------------------------------
	;
	; Rules:
	;   - NSUBS must be equal to # dependent variables, otherwise subscripting is indeterminate
	;   - Dimensions of implicit array are ordered [DEPEND_0, DEPEND_1, ..., DEPEND_N]
	;   - DEPEND_[1-N] may have record variance, in which case they are 2D
	;   - If a bulk process currently acting on the implicit, and then it may already
	;     have affected the (cached) DEPEND_[1-N] variable. If DEPEND_[1-N] is cached,
	;     allow it to have the same dimensions as the output variable.
	;
	
	;
	; TODO:
	;   Ideally, this should be done with all poitner variables:
	;   https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html
	;       FORM_PTR
	;       LABL_PTR_[1-3]
	;       SCAL_PTR
	;       UNITS_PTR
	;
	
	;Number of dependent variables
	;   - IDL allows up to 8 dimensions
	nDeps = ~self -> HasAttr('DEPEND_0')   ? 0 $
	        : ~self -> HasAttr('DEPEND_1') ? 1 $
	        : ~self -> HasAttr('DEPEND_2') ? 2 $
	        : ~self -> HasAttr('DEPEND_3') ? 3 $
	        : ~self -> HasAttr('DEPEND_4') ? 4 $
	        : ~self -> HasAttr('DEPEND_5') ? 5 $
	        : ~self -> HasAttr('DEPEND_6') ? 6 $
	        : ~self -> HasAttr('DEPEND_7') ? 7 $
	        : 8
	
	;RETURN if nSubs LE nDeps
	;   - There does not have to be a DEPEND_# for each subscript
	;   - There DOES     have to be a subscript for each DEPEND_#
	;   - Issure warning from program that called ::_OverloadBracketsRightSide
	IF nSubs LT nDeps THEN BEGIN
		strNSubs = String(nSubs, FORMAT='(i1)')
		strNDeps = String(nDeps, FORMAT='(i1)')
		IF nDeps GT 0 THEN MrPrintF, LEVEL=6, 'LogWarn', '# subscripts (' + strNSubs + ') LT # dependent variables (' + strNDeps + ').'
		RETURN, attrs
	ENDIF
	
	;Step through each dependent variable
	FOR i = 0, nDeps-1 DO BEGIN
		depend = 'DEPEND_' + String(i, FORMAT='(i0)')
		IF ~self -> HasAttr(depend) THEN CONTINUE

		;Get the dependent variable and its dimensions
		oDep     = self[depend]
		nDepDims = Size(oDep, /N_DIMENSIONS)
		depDims  = Size(oDep, /DIMENSIONS)
		
		;RECORD VARIANCE
		;   - Data     = [N, M, L, P]
		;   - DEPEND_0 = [N]
		;   - DEPEND_1 = [N, M]
		;   - DEPEND_2 = [N, L]
		;   - DEPEND_3 = [N, P]
		IF depend NE 'DEPEND_0' && nDepDims EQ 2 THEN BEGIN
			;Apply to DEPEND_#
			;   - Dimensions of implicit and DEPEND_# arrays match
			IF Array_Equal(depDims, inDims[[0,i]]) THEN BEGIN
				CASE i OF
					1: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i2 )
					2: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i3 )
					3: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i4 )
					4: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i5 )
					5: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i6 )
					6: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i7 )
					7: oDepend = oDep -> _OverloadBracketsRightSide( isRange[[0,i]], i1, i8 )
					ELSE: Message, 'Invalid number of dependent variables.'
				ENDCASE
			
			;Pre-Applied
			;   - It is possible that whatever is happening to the implicit variable
			;     is happening in bulk (e.g. MrVar_TLimit).
			;   - If DEPEND_# is in the cache, allow it to have the same dimensions
			;     as the output variable
			ENDIF ELSE IF MrVar_IsCached(oDep) && Array_Equal(depDims, outDims[[0,i]]) THEN BEGIN
				oDepend = oDep
			
			;UNKNOWN DIMENSION SIZES
			;   - Warning with Expected Dimensions
			ENDIF ELSE BEGIN
				strDims = '[' + StrJoin(String(inDims[[0,i]], FORMAT='(i0)'), ', ') + ']'
			ENDELSE
			
		;NO RECORD VARIANCE
		;   - Data     = [N, M, L, P]
		;   - DEPEND_0 = [N]
		;   - DEPEND_1 = [M]
		;   - DEPEND_2 = [L]
		;   - DEPEND_3 = [P]
		ENDIF ELSE IF nDepDims EQ 1 THEN BEGIN
			;Apply to DEPEND_#
			;   - Dimensions of implicit and DEPEND_# arrays match
			;   - Prevent scalar 0 from returning object by including "[]"
			IF depDims EQ inDims[i] THEN BEGIN
				CASE i OF
					0: oDepend = oDep -> _OverloadBracketsRightSide( isRange[0], [i1] )
					1: oDepend = oDep -> _OverloadBracketsRightSide( isRange[1], [i2] )
					2: oDepend = oDep -> _OverloadBracketsRightSide( isRange[2], [i3] )
					3: oDepend = oDep -> _OverloadBracketsRightSide( isRange[3], [i4] )
					4: oDepend = oDep -> _OverloadBracketsRightSide( isRange[4], [i5] )
					5: oDepend = oDep -> _OverloadBracketsRightSide( isRange[5], [i6] )
					6: oDepend = oDep -> _OverloadBracketsRightSide( isRange[6], [i7] )
					7: oDepend = oDep -> _OverloadBracketsRightSide( isRange[7], [i8] )
					ELSE: Message, 'Invalid number of dependent variables.'
				ENDCASE
			
			;Pre-Applied
			;   - It is possible that whatever is happening to the implicit variable
			;     is happening in bulk (e.g. MrVar_TLimit).
			;   - If DEPEND_# is in the cache, allow it to have the same dimensions
			;     as the output variable
			ENDIF ELSE IF MrVar_IsCached(oDep) && depDims EQ outDims[i] THEN BEGIN
				oDepend = oDep
			
			;UNKNOWN DIMENSION SIZES
			;   - Warning with Expected Dimensions
			ENDIF ELSE BEGIN
				strDims = '[' + StrJoin(String(inDims[i], FORMAT='(i0)'), ', ') + ']'
			ENDELSE
		
		;UNKNOWN DIMENSION SIZES
		;   - Warning with Expected Dimensions
		ENDIF ELSE BEGIN
			strDims = '[' + StrJoin(String(inDims[i], FORMAT='(i0)'), ', ') + ']'
		ENDELSE
		
		;Set the DEPEND_# attribute
		IF Obj_Valid(oDepend) THEN BEGIN
			attrs[depend] = oDepend
		
		;Issue warning
		ENDIF ELSE BEGIN
			IF self.verbose EQ 1 THEN BEGIN
				MrPrintF, 'LogWarn', LEVEL=6, 'Cannot exract ' + depend + ' data for variable "' + self.name + '".'
			ENDIF ELSE IF self.verbose EQ 2 THEN BEGIN
				;Dimensions
				strInDims  = '[' + StrJoin(String(inDims,    FORMAT='(i0)'), ', ') + ']'
				strOutDims = '[' + StrJoin(String(outDims,   FORMAT='(i0)'), ', ') + ']'
				strDepDims = '[' + StrJoin(String(depDims,   FORMAT='(i0)'), ', ') + ']'
		
				;Warning information
				;   - Display warning from caller of ::_OverloadBracketsRightSide
				MrPrintF, 'LogWarn', LEVEL=6, 'Cannot extract ' + depend + ' data for variable "' + self.name + '".'
				MrPrintF, 'LogWarn', LEVEL=6, '    ' + depend + ' dims:  ' + strDepDims + '.'
				MrPrintF, 'LogWarn', LEVEL=6, '    Data dims:      ' + strInDims
				MrPrintF, 'LogWarn', LEVEL=6, '    Subarray dims:  ' + strOutDims
				MrPrintF, 'LogWarn', LEVEL=6, '    Expected dims:  ' + strDims
			ENDIF
		ENDELSE
		
		;Invalidate the DEPEND_# variable
		oDepend = 0B
	ENDFOR

;-------------------------------------------
; DELTA_(PLUS|MINUS)_VAR ///////////////////
;-------------------------------------------
	deltas  = ['DELTA_MINUS_VAR', 'DELTA_PLUS_VAR']
	nDeltas = N_Elements(deltas)
	
	;Step through each dependent variable
	FOR i = 0, nDeltas-1 DO BEGIN
		IF ~self -> HasAttr(deltas[i]) THEN CONTINUE
		
		;Get the dependent variable and its dimensions
		oDel       = self[deltas[i]]
		nDelDims = Size(oDel, /N_DIMENSIONS)
		delDims  = Size(oDel, /DIMENSIONS)
		
		;SAME
		;   - Both record varying
		;   - Both non-record varying
		IF nDelDims EQ nInDims && Array_Equal(delDims, inDims) THEN BEGIN
			CASE nDelDims OF
				1: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0],   i1)
				2: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:1], i1, i2)
				3: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:2], i1, i2, i3)
				4: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:3], i1, i2, i3, i4)
				5: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:4], i1, i2, i3, i4, i5)
				6: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:5], i1, i2, i3, i4, i5, i6)
				7: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:6], i1, i2, i3, i4, i5, i6, i7)
				8: oDelta = oDel -> _OverloadBracketsRightSide(isRange[0:7], i1, i2, i3, i4, i5, i6, i7, i8)
				ELSE: Message, 'Incorrect number of dimensions (' + String(nDelDims, FORMAT='(i0)') + ').'
			ENDCASE
		
		;DIFFERENT
		;   - Parent is record varying
		;   - DELTA is non-record varying
		ENDIF ELSE IF nDelDims EQ nInDims-1 && Array_Equal(delDims, inDims[1:*]) THEN BEGIN
			CASE nDelDims OF
				1: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1],   i2)
				2: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:2], i2, i3)
				3: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:3], i2, i3, i4)
				4: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:4], i2, i3, i4, i5)
				5: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:5], i2, i3, i4, i5, i6)
				6: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:6], i2, i3, i4, i5, i6, i7)
				7: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:7], i2, i3, i4, i5, i6, i7, i8)
				ELSE: Message, 'Incorrect number of dimensions (' + String(nDelDims, FORMAT='(i0)') + ').'
			ENDCASE
		
		;INCOMPATIBLE
		ENDIF ELSE BEGIN
			IF self.verbose EQ 1 THEN BEGIN
				MrPrintF, 'LogWarn', LEVEL=6, 'Cannot extract ' + deltas[i] + ' data for variable "' + self.name + '".'
			ENDIF ELSE IF self.verbose EQ 2 THEN BEGIN
				;Dimensions
				strInDims  = '[' + StrJoin(String(inDims,    FORMAT='(i0)'), ', ') + ']'
				strOutDims = '[' + StrJoin(String(outDims,   FORMAT='(i0)'), ', ') + ']'
				strDims    = '[' + StrJoin(String(inDims[0], FORMAT='(i0)'), ', ') + ']'
				strDelDims = '[' + StrJoin(String(delDims,   FORMAT='(i0)'), ', ') + ']'
		
				;Warning information
				;   - Display warning from caller of ::_OverloadBracketsRightSide
				MrPrintF, 'LogWarn', LEVEL=6, 'Cannot extract ' + deltas[i] + ' data for variable "' + self.name + '".'
				MrPrintF, 'LogWarn', LEVEL=6, '    ' + deltas[i] + ' dims:  ' + strDelDims + '.'
				MrPrintF, 'LogWarn', LEVEL=6, '    Data dims:             ' + strInDims
				MrPrintF, 'LogWarn', LEVEL=6, '    Subarray dims:         ' + strOutDims
			ENDIF
			
			;Invalidate the DEPEND_# variable
			oDelta = 0B
		ENDELSE
		
		;Create a variable
		IF Obj_Valid(oDelta) THEN attrs[deltas[i]] = oDelta
	ENDFOR

;-------------------------------------------
; Finished! ////////////////////////////////
;-------------------------------------------
	RETURN, attrs
END


;+
;   Allow square-bracket array indexing from the right side of an operator.
;
;   Calling Sequence
;       data = oTSvar[i1 [, i2 [, i3 [, i4 [, i5 [, i6 [, i7 [, i8]]]]]]]]
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
;       DATA:               out, required, type=numeric array
;                           The subarray accessed by the input parameters.
;-
FUNCTION MrVariable::_OverloadBracketsRightSide_Data, isRange, i1, i2, i3, i4, i5, i6, i7, i8
	Compile_Opt idl2
	On_Error, 2

	;Number of subscripts given
	nSubscripts = N_Elements(isRange)
	
;---------------------------------------------------------------------
;Optimized Subscripting for <= 3D ////////////////////////////////////
;---------------------------------------------------------------------
	IF (nSubscripts le 3) THEN BEGIN
	;---------------------------------------------------------------------
	;3D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		IF IsA(i3) THEN BEGIN 
			;Subscript range given: [min:max:interval]?
			IF isRange[2] THEN BEGIN 
				;Subscript range for dimensions 2 and 3?
				IF isRange[1] THEN BEGIN 
					;Range: [1,2,3], Index: --
					IF isRange[0] THEN BEGIN
						RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2],i3[0]:i3[1]:i3[2]]
					;Range: [2,3], Index: 1
					ENDIF ELSE BEGIN
						RETURN, (*self.data)[i1,i2[0]:i2[1]:i2[2],i3[0]:i3[1]:i3[2]] 
					ENDELSE
				;Index value for dimension 2?
				ENDIF ELSE BEGIN
					;Range: [3,1], Index: 2
					IF isRange[0] THEN BEGIN 
						RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2,i3[0]:i3[1]:i3[2]]
					;Range: 3, Index: [2,1]
					ENDIF ELSE BEGIN 
						RETURN, (*self.data)[i1,i2,i3[0]:i3[1]:i3[2]] 
					ENDELSE 
				ENDELSE
			;Index for dimension 3?
			ENDIF ELSE BEGIN
				;Range for dimension 2?
				IF isRange[1] THEN BEGIN
					;Range: [2,1]
					IF isRange[0] THEN BEGIN 
						RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2],i3] 
					ENDIF ELSE BEGIN 
						RETURN, (*self.data)[i1,i2[0]:i2[1]:i2[2],i3] 
					ENDELSE 
				ENDIF ELSE BEGIN 
					IF isRange[0] THEN BEGIN 
						RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2,i3] 
					ENDIF ELSE BEGIN 
						RETURN, (*self.data)[i1,i2,i3] 
					ENDELSE 
				ENDELSE 
			ENDELSE 
	;---------------------------------------------------------------------
	;2D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE IF IsA(i2) THEN BEGIN
			IF isRange[1] THEN BEGIN
				;[Range, Range]
				IF isRange[0] THEN BEGIN
					RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2[0]:i2[1]:i2[2]]
				;[Index, Range]
				ENDIF ELSE BEGIN
					RETURN, (*self.data)[i1,i2[0]:i2[1]:i2[2]]
				ENDELSE
			ENDIF ELSE BEGIN
				;[Range, Index]
				IF isRange[0] THEN BEGIN
					RETURN, (*self.data)[i1[0]:i1[1]:i1[2],i2]
				;[Index, Index]
				ENDIF ELSE BEGIN
					RETURN, (*self.data)[i1,i2]
				ENDELSE
			ENDELSE
	;---------------------------------------------------------------------
	;1D Subscripts ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE BEGIN
			;Range?
			IF isRange[0] THEN BEGIN
				RETURN, (*self.data)[i1[0]:i1[1]:i1[2]] 
			;Index
			;   - Compensate for passing in [0] instead of 0 for the first element.
			;   - i.e. RETURN a scalar instead of a 1-element array
			ENDIF ELSE BEGIN
				IF N_Elements(i1) eq 1 && i1 eq 0 $
					THEN RETURN, (*self.data)[0] $
					ELSE RETURN, (*self.data)[i1]
			ENDELSE
		ENDELSE
	ENDIF

;---------------------------------------------------------------------
; Brute Force Code for 4D or Higher Arrays. //////////////////////////
;---------------------------------------------------------------------
	;Works for any number of dimensions.
	dims = Size(*self.data, /DIMENSIONS)
	indices = MrReformIndices(dims, isRange, i1, i2, i3, i4, i5, i6, i7, i8, DIMENSIONS=dimensions);, /IDL_METHOD)

	RETURN, Reform((*self.data)[indices], dimensions, /OVERWRITE)
END


;+
;   Allow square-bracket array indexing from the right side of an operator.
;
;   Calling Sequence
;       oVar      = oVar[0]
;       d0        = oVar[[0]]
;       data      = oVar['DATA' [, i2 [, i3 [, i4 [, i5 [, i6 [, i7 [, i8]]]]]]]
;       pData     = oVar['POINTER']
;       pData     = oVar['PTR']
;       attrValue = oVar[AttrName]
;       oVar      = oVar[i1 [, i2 [, i3 [, i4 [, i5 [, i6 [, i7 [, i8]]]]]]]]
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
;       RESULT:             out, required, type=MrVariable objref
;                           The subarray accessed by the input parameters.
;-
FUNCTION MrVariable::_OverloadBracketsRightSide, isRange, i1, i2, i3, i4, i5, i6, i7, i8
	Compile_Opt idl2
	On_Error, 2

	;Number of subscripts given
	nSubs = N_Elements(isRange)

	;String operations.
	IF IsA(i1, /SCALAR, 'STRING') THEN BEGIN
		CASE i1 OF
			'DATA': BEGIN
				IF nSubs EQ 1 $
					THEN RETURN, *self.data $
					ELSE RETURN, self -> _OverloadBracketsRightSide_Data(isRange[1:*], i2, i3, i4, i5, i6, i7, i8)
			ENDCASE
			'POINTER': RETURN, self.data
			'PTR':     RETURN, self.data
			ELSE:      RETURN, self -> GetAttrValue(i1)
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
	vOut = MrVariable( data, $
	                   NAME='OverloadBRS(' + self.name + ')', $
	                   /NO_COPY )
	
	;Copy all attributes
	self -> CopyAttrTo, vOut
	
	;Update attributes
	IF attrs -> Count() GT 0 THEN vOut -> SetAttrValue, attrs, /OVERWRITE
	
;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------

	RETURN, vOut
END




;+
;   The purpose of this method is to check for equality. See IDL's online help for
;   `Relational Operators <http://exelisvis.com/docs/Relational_Operators.html>`
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
;       RESULT:             True (1) if `LEFT` equals `RIGHT`. Otherwise false (0)
;-
function MrVariable::_OverloadEQ, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Add
		result = (*self.data) eq (right -> GetData())
		
		;Create a new name
		name = 'EQ(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) eq right['DATA'] $
			else result = left['DATA'] eq (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'EQ(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'EQ(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


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
function MrVariable::_OverloadForeach, value, key
	compile_opt idl2
	on_error, 2

	nPts = n_elements(*self.data)
	if n_elements(key) eq 0 then key = 0

	;Get the array element if the index is in range
	if key lt nPts then begin
		next = 1
		value = (*self.data)[key]
	
	;Otherwise, stop iterating 
	endif else next = 0

	;Next element to retrieve
	key += 1

	return, next
end


;+
;   The purpose of this method is to check if `LEFT` is greater than or equal to `RIGHT`.
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
;       RESULT:             True (1) if `LEFT` greater than or equal to `RIGHT`.
;                               Otherwise false (0).
;-
function MrVariable::_OverloadGE, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Add
		result = (*self.data) ge (right -> GetData())
		
		;Create a new name
		name = 'GE(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) ge right $
			else result = left ge (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'GE(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'GE(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to check if `LEFT` is greater than `RIGHT`. See IDL's
;   online help for `Relational Operators <http://exelisvis.com/docs/Relational_Operators.html>`
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
;       RESULT:             True (1) if `LEFT` greater than or equal to `RIGHT`.
;                               Otherwise false (0).
;-
function MrVariable::_OverloadGT, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Add
		result = (*self.data) gt (right -> GetData())
		
		;Create a new name
		name = 'GT(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) gt right $
			else result = left gt (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'GT(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'GT(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to return the element-by-element maximum between
;   `LEFT` and `RIGHT`. See IDL's online help for `Relational Operators <http://exelisvis.com/docs/Relational_Operators.html>`
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
;       RESULT:             Returns the element of `LEFT` if it is greater than `RIGHT`,
;                               and vice versa.
;-
function MrVariable::_OverloadGreaterThan, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Pick greater value
		result = (*self.data) > (right -> GetData())
		
		;Create a new name
		name = 'GreaterThan(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) > right $
			else result = left > (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'GreaterThan(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'GreaterThan(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to provide information when the HELP procedure
;   is called.
;
; :Params:
;       VARNAME:        in, required, type=string
;                       Name of the variable supplied to the HELP procedure.
;-
function MrVariable::_OverloadHelp, varname
	compile_opt idl2
	on_error, 2

	;Info about the object
	heapnum  = obj_valid(self, /GET_HEAP_IDENTIFIER)
	type     = size(self, /TNAME)
	class    = obj_class(self)
	dims     = size(*self.data, /DIMENSIONS)
	datatype = size(*self.data, /TNAME)

	;Name, variable name, and object reference info of variable.
	fmt_varname = 'a-' + string(strlen(varname)     > 12, FORMAT='(i0)')
	fmt_name    = 'a-' + string(strlen(self.name)+2 > 12, FORMAT='(i0)')
	str = string(varname, '"' + self.name + '"', '<', heapnum, ' (' + class + ')>', $
	             FORMAT='(' + fmt_varname + ', ' + fmt_name + ', 3x, a1, i0, a0)')

	;Description of datatype and dimension sizes
	;   - If one element, give value.
	if n_elements(*self.data) eq 1 $
		then desc = string(datatype, FORMAT='(a-8)') + ' = ' + strtrim(string(*self.data, /PRINT), 2) $
		else desc = string(datatype, FORMAT='(a-8)') + ' = Array[' + strjoin(string(dims, FORMAT='(i0)'), ',') + ']'

	;Combine the two
	str += '  ' + desc
	return, str
end


;+
;   The purpose of this method is to provide information when implied print is used.
;
; :Params:
;       VARNAME:        in, required, type=string
;                       Name of the variable supplied to Implied Print.
;-
function MrVariable::_OverloadImpliedPrint, varname
	compile_opt idl2
	on_error, 2

	;Print the array
	return, *self.data
end


;+
;   The purpose of this method is to evaluate if the array is true or not. Called when
;   logical operators are used (&&, ||, etc.)
;
; :Returns:
;       RESULT:         True (1) if the implicit array is a scalar or 1 element array
;                           not equal to zero and false if a scalar not equal to zero.
;                           Arrays with more than one element result in an error.
;-
function MrVariable::_OverloadIsTrue
	compile_opt idl2
	on_error, 2

	;Make sure the array is a scalar
	if self.n_elements gt 1 then message, 'Array must be a scalar or 1 element array.'

	;Is it true?
	result = (*self.data)[0] ne 0

	;Print the array
	return, result
end


;+
;   The purpose of this method is to return the element-by-element minimum between
;   `LEFT` and `RIGHT`.
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
;       RESULT:             Returns the element of `LEFT` if it is less than `RIGHT`,
;                               and vice versa.
;-
function MrVariable::_OverloadLessThan, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Pick lesser value
		result = (*self.data) < right['DATA']
		
		;Create a new name
		name = 'LessThan(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) < right $
			else result = left < (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'LessThan(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'LessThan(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to check if `LEFT` is less than or equal to `RIGHT`.
;   See IDL's online help for `Relational Operators <http://exelisvis.com/docs/Relational_Operators.html>`
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
;       RESULT:             True (1) if `LEFT` less than or equal to `RIGHT`.
;                               Otherwise false (0).
;-
function MrVariable::_OverloadLE, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Compare
		result = (*self.data) le (right -> GetData())
		
		;Create a new name
		name = 'LE(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) le right $
			else result = left le (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'LE(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'LE(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to check if `LEFT` is less than `RIGHT`. See IDL's
;   online help for `Relational Operators <http://exelisvis.com/docs/Relational_Operators.html>`
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
;       RESULT:             True (1) if `LEFT` less than to `RIGHT`. Otherwise false (0).
;-
function MrVariable::_OverloadLT, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Compare
		result = (*self.data) lt (right -> GetData())
		
		;Create a new name
		name = 'LT(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) lt right $
			else result = left lt (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'LT(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'LT(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to subract two expressions.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression, normal IDL operations
;     will take effect. This means the output will be the size of the
;     array with the smallest dimensions.
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
;       RESULT:             The result of subtracting `RIGHT` from `LEFT`.
;-
function MrVariable::_OverloadMinus, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects //////////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Subtract
		result = (*self.data) - (right -> GetData())
		
		;Create a new name
		name = 'Minus(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression //////////////////////////////
;-------------------------------------------------------
	endif else begin
		;Subtract the expressions
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) - right $
			else result = left - (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Minus(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Minus(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to negate the implicit array.
;
; :Returns:
;       RESULT:             Negation of the implicit array.
;-
function MrVariable::_OverloadMinusUnary
	compile_opt idl2
	on_error, 2

	;Negate the array, making positive values negative, and vice versa
	return, self -> New(-(*self.data), NAME='-'+self.name)
end


;+
;   The purpose of this method is to take the MOD of `LEFT` and `RIGHT`.
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
;       RESULT:             The remainder of `LEFT` divided by `RIGHT` (i.e. the MOD).
;-
function MrVariable::_OverloadMOD, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Modulus
		result = (*self.data) mod (right -> GetData())
		
		;Create a new name
		name = 'MOD(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) mod right $
			else result = left mod (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'MOD(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'MOD(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to add two expressions.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression, normal IDL operations
;     will take effect. This means the output will be the size of the
;     array with the smallest dimensions.
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
;       RESULT:             The result of adding `RIGHT` to `LEFT`.
;-
function MrVariable::_OverloadPlus, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects //////////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Add
		result = (*self.data) + (right -> GetData())
		
		;Create a new name
		name = 'Plus(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression //////////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) + right $
			else result = left + (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Plus(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Plus(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to check for inequality. See IDL's online help for
;   `Relational Operators <http://exelisvis.com/docs/Relational_Operators.html>`
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
;       RESULT:             True (1) if `LEFT` is not equal to `RIGHT`. Otherwise false (0).
;-
function MrVariable::_OverloadNE, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Not Equal
		result = (*self.data) ne (right -> GetData())
		
		;Create a new name
		name = 'NE(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) ne right $
			else result = left ne (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'NE(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'NE(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to perform a logical NOT of the implicit array. See
;   the IDL Help page for `Bitwise Operators <http://exelisvis.com/docs/Bitwise_Operators.html>`
;   for more information.
;
; :Returns:
;       RESULT:             The logical NOT of the implicit array.
;-
function MrVariable::_OverloadNOT
	compile_opt idl2
	on_error, 2

	;Negate the array, making positive values negative, and vice versa
	return, self -> New(not (*self.data), NAME='NOT('+self.name+')')
end


;+
;   The purpose of this method is to perform an inclusive OR between `LEFT` and `RIGHT`.
;   SEE the IDL Help page for `Bitwise Operators <http://exelisvis.com/docs/Bitwise_Operators.html>`
;   for more information.
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
;       RESULT:             The result of `LEFT` OR `RIGHT`.
;-
function MrVariable::_OverloadOR, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Or
		result = (*self.data) or (right -> GetData())
		
		;Create a new name
		name = 'OR(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) or right $
			else result = left or (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'OR(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'OR(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to multiply the columns of `LEFT` by the rows of `RIGHT`.
;   (i.e., an IDL matrix multiplication, not a mathematical matrix multiplication). See
;   IDL's page for `Matrix Operators <http://exelisvis.com/docs/Matrix_Operators.html>`
;   for more information.
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
;       RESULT:             The result of matrix multiplying `LEFT` by `RIGHT` in the IDL
;                               sense.
;-
function MrVariable::_OverloadPound, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Add
		result = (*self.data) # (right -> GetData())
		
		;Create a new name
		name = 'Pound(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) # right $
			else result = left # (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Pound(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Pound(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to multiply the rows of `LEFT` by the columns of `RIGHT`.
;   (i.e., a mathematical matrix multiplication, not an IDL matrix multiplication). See
;   IDL's page for `Matrix Operators <http://exelisvis.com/docs/Matrix_Operators.html>`
;   for more information.
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
;       RESULT:             The result of matrix multiplying `LEFT` by `RIGHT` in the
;                               mathematical sense.
;-
function MrVariable::_OverloadPoundPound, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Add
		result = (*self.data) ## (right -> GetData())
		
		;Create a new name
		name = 'PoundPound(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) ## right $
			else result = left ## (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'PoundPound(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'PoundPound(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to provide information when implied print is used.
;   is called.
;-
function MrVariable::_OverloadPrint
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
function MrVariable::_OverloadSize
	compile_opt idl2
	on_error, 2

	;NOTE:
	;   There is a bug that incorrectly reports the number of dimesions of an
	;   undefined variable (returns 1 instead of 0). This was fixed in
	;   version 8.3.1.
	;
	;   https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/-f8Cxlp5cxQ
	;
	;   A related bug incorrectly reports the number of elements of a scalar
	;   variable (returns 0 instead of 1).
	;

	;Return the dimensions of the array. The Size and N_Elements functions
	;will know what to do with them.
	return, n_elements(*self.data)          eq 0 ? 0L : $
	        size(*self.data, /n_dimensions) eq 0 ? 1L : $
	        size(*self.data, /dimensions) 
end


;+
;   The purpose of this method is to multiply two expressions together.
;
;   NOTE:
;     If one of LEFT or RIGHT is an expression, normal IDL operations
;     will take effect. This means the output will be the size of the
;     array with the smallest dimensions.
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
function MrVariable::_OverloadSlash, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects //////////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Divide
		result = (*self.data) / (right -> GetData())
		
		;Create a new name
		name = 'Divide(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression //////////////////////////////
;-------------------------------------------------------
	endif else begin
		;Divide the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) / right $
			else result = left / (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'Divide(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'Divide(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   The purpose of this method is to take the logical NOT (~) of the implicit array. See
;   IDL's online help page for `Logical Operators <http://exelisvis.com/docs/Logical_Operators.html>`
;
; :Returns:
;       RESULT:             The logical not of the implicit array.
;-
function MrVariable::_OverloadTilde
	compile_opt idl2
	on_error, 2

	;Negate the array, making positive values negative, and vice versa
	return, self -> New(~(*self.data), NAME='~'+self.name)
end


;+
;   The purpose of this method is to perform an eXlusive OR between `LEFT` and `RIGHT`.
;   SEE the IDL Help page for `Bitwise Operators <http://exelisvis.com/docs/Bitwise_Operators.html>`
;   for more information.
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
;       RESULT:             The result of `LEFT` XOR `RIGHT`.
;-
function MrVariable::_OverloadXOR, left, right
	compile_opt idl2
	on_error, 2

	;Is SELF on the left or right?
	;   - Are both LEFT and RIGHT GDA data objects?
	side = self -> _LeftOrRight(left, right, ISMrVariable=IsMrVariable)

;-------------------------------------------------------
; Two MrVariable Objects ///////////////////////////////
;-------------------------------------------------------
	;Both LEFT and RIGHT are MrVariable objects
	;   - IDL will call _Overload for LEFT first
	;   - oVar -> GetData() is faster than oVar['DATA']
	if IsMrVariable then begin
		;Add
		result = (*self.data) xor (right -> GetData())
		
		;Create a new name
		name = 'XOR(' + self.name + ',' + right.name + ')'

;-------------------------------------------------------
; MrVariable with Expression ///////////////////////////
;-------------------------------------------------------
	endif else begin
		;Add the expressions
		;   - Assume the user knows what they are doing.
		;   - All IDL truncation effects apply (shortest in determines size out).
		if side eq 'LEFT' $
			then result = (*self.data) xor right $
			else result = left xor (*self.data)
		
		;Determine name
		;   - Scalar or TYPE[dims]
		if side eq 'LEFT' then begin
			if n_elements(right) eq 1 $
				then rname = strtrim(right[0], 2) $
				else rname = size(right, /TNAME) + '[' + strjoin(strtrim(size(right, /DIMENSIONS), 2), ',') + ']'
			name = 'XOR(' + self.name + ',' + rname + ')'
		endif else begin
			if n_elements(left) eq 1 $
				then lname = strtrim(left[0], 2) $
				else lname = size(left, /TNAME) + '[' + strjoin(strtrim(size(left, /DIMENSIONS), 2), ',') + ']'
			name = 'XOR(' + lname + ',' + self.name + ')'
		endelse
	endelse

;-------------------------------------------------------
; Output Object ////////////////////////////////////////
;-------------------------------------------------------
	;Create a new object based on the results
	return, self -> New(result, /NO_COPY, NAME=name)
end


;+
;   Set attribute values.
;
; :Params:
;       ATTRNAME:       in, required, type=string/hash/struct
;                       The name of the attribute for which the value is to be changed,
;                           or a hash or structure whos keys/tags are the attribute names.
;                           If a hash, keys must be strings. Values cannot be complex
;                           datatypes.
;       ATTRVALUE:      in, optional, type=any
;                       The value of the attribute(s) to be added. ATTRVALUE must be
;
; :Keyword:
;       OVERWRITE:      in, optional, type=boolean, default=0
;                       If set, attributes that already exist will be over-written. The
;                           default is to issue a warning.
;-
pro MrVariable::AddAttr, attrName, attrValue, $
OVERWRITE=overwrite
	compile_opt idl2
	on_error, 2
	
	;Default to not overwriting
	tf_overwrite = keyword_set(overwrite)
	
	;Create the attribute in addition to setting its value
	self -> SetAttrValue, attrName, attrValue, $
	                      /CREATE, $
	                      OVERWRITE=tf_overwrite
end


;+
;   Get the number of attributes
;
; :Returns:
;       COUNT:          out, required, type=integer
;                       The number of attributes.
;-
FUNCTION MrVariable::AttrCount
	Compile_Opt idl2
	On_Error, 2
	
	RETURN, self.attributes -> Count()
END


;+
;   Get the names of all attributes
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Named variable to recieve the number of attribute names.
;-
function MrVariable::AttrNames, $
COUNT=count
	compile_opt idl2
	on_error, 2
	
	;Hash returns a list that we want to convert to a string array
	names = self.attributes -> Keys()
	names = names -> ToArray()
	
	;Return the count as well?
	if arg_present(count) then count = self.attributes -> Count()
	
	return, names
end


;+
;   Determine if an attribute value is a valid datatype.
;
; :Private:
;
; :Params:
;       VALUE:          in, required, type=any
;                       Value for which the datatype is to be validated.
;
; :Keywords:
;       TYPE_NAME:      in, optional, type=boolean, default=0
;                       If set, `VALUE` is taken to be a datatype name.
;       TYPE_CODE:      in, optional, type=boolean, default=0
;                       If set, `VALUE` is taken to be a datatype code.
;-
function MrVariable::AttrValue_IsValid, value, $
TYPE_NAME=type_name, $
TYPE_CODE=type_code
	compile_opt idl2
	on_error, 2
	
	;Check keywords
	tf_type_name = keyword_set(type_name)
	tf_type_code = keyword_set(type_code)
	if tf_type_name + tf_type_code gt 1 then message, 'TYPE_NAME and TYPE_CODE are mutually exclusive.'
	
	;Get the type code
	if tf_type_name then begin
		tcode = self -> TypeName2Code(value)
	endif else if tf_type_code then begin
		tcode = value
	endif else begin
		tcode = size(value, /TYPE)
	endelse
	
	;Check valid values
	case tcode of
		 0: tf_valid = 0B ;UNDEFINED
		 1: tf_valid = 1B ;BYTE
		 2: tf_valid = 1B ;INT
		 3: tf_valid = 1B ;LONG
		 4: tf_valid = 1B ;FLOAT
		 5: tf_valid = 1B ;DOUBLE
		 6: tf_valid = 0B ;COMPLEX
		 7: tf_valid = 1B ;STRING
		 8: tf_valid = 0B ;STRUCT
		 9: tf_valid = 1B ;DCOMPLEX -- Epoch16
		10: tf_valid = 0B ;POINTER
		11: tf_valid = 0B ;OBJREF
		12: tf_valid = 1B ;UINT
		13: tf_valid = 1B ;ULONG
		14: tf_valid = 1B ;LONG64
		15: tf_valid = 1B ;ULONG64
		else: message, 'Invalid datatype.'
	endcase
	
	return, tf_valid
end


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
; :Examples:
;   Append to the 3rd dimension::
;       myArray  = MrVariable(12, 3, 6, 4, TYPE='FLOAT')
;       myArray -> Append, fltarr(12, 3, 4, 4), 3
;       help, myArray
;            MYARRAY     OBJREF      <ObjHeapVar5(MrVariable)>
;             ARRAY       FLOAT     = Array[12, 3, 10, 4]
;
;  Multiple appends to the third dimension::
;       myArray  = MrVariable(12, 3, 6, 4, TYPE='FLOAT')
;       myArray -> Append, fltarr(12, 3, 4, 4), 3
;       myArray -> Append, fltarr(12, 3, 4, 4), 3
;       myArray -> Append, fltarr(12, 3, 4, 4), 3
;       help, myArray
;           MYARRAY     OBJREF      <ObjHeapVar5(MrVariable)>
;             ARRAY       FLOAT     = Array[12, 3, 18, 4]
;
;   Multiple appends to the beginning of an array::
;       myArray -> MrVariable(lonarr(2, 2))
;       myArray -> Append, 1, /BEFORE
;       myArray -> Append, lonarr(2,2) + 1, 1
;       myArray -> Append, lonarr(2,2) + 2, 1
;       myArray -> Append, lonarr(2,2) + 3, 1
;       help, myArray
;           MYARRAY     OBJREF      <ObjHeapVar17(MrVariable)>
;             ARRAY       LONG      = Array[8, 2]
;       print, myArray
;           1      1     2     2     3     3     0      0
;           1      1     2     2     3     3     0      0
;
; :Params:
;       DATA:           in, required, type=array
;                       Array to be concatenated to the implicit array. If undefined or
;                           !Null, then nothing is appended.
;       DIMENSION:      in, optional, type=integer, default=1
;                       The dimension, starting with 1, along which the two arrays are
;                           to be concatenated. For calls between `START` and `FINISH`,
;                           this parameter is ignored.
;
; :Keywords:
;       BEFORE:         in, optional, type=boolean, default=0
;                       If set, `DATA` will be appended to the beginning of the implicit
;                           array. The default is to append to the end.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       If set, `DATA` will be copied directly into the object and
;                           will be left undefined.
;-
pro MrVariable::Append, data, dimension, $
BEFORE=before, $
NO_COPY=no_copy
	compile_opt idl2
	on_error, 2
	
	;Allow input data to be empty. In this case, do nothing.
	if n_elements(data) eq 0 then return

	;Defaults
	before  = keyword_set(before)
	no_copy = keyword_set(no_copy)
	if n_elements(dimension) eq 0 then dimension = 1

	;Is the object empty?
	if n_elements(*self.data) eq 0 then begin
		self -> SetData, data, NO_COPY=no_copy
	
	;Append
	endif else begin
		if before eq 1 $
			then *self.data = MrConcatenate(data, *self.data, dimension) $
			else *self.data = MrConcatenate(*self.data, data, dimension)
	endelse
end


;+
;   A helper method for cleaning up after the ::Append_Multi method.
;
; :Private:
;-
pro MrVariable::Append_Finalize
	compile_opt idl2, hidden

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		self._append_cont    = 0B
		self._append_before  = 0B
		self._append_no_copy = 0B
		MrPrintF, 'LogErr'
		return
	endif

	;Was data added?
	if n_elements(*self.data) eq 0 then return

	;Put the desired dimension first
	if self._append_putFirst then begin
		*self.data = MrPutDimensionFirst(*self.data, self._dimAppend, $
		                                 UNREFORM=unReform, UNTRANSPOSE=unTranspose)
	endif else begin
		unReform    = *self._append_unReform
		unTranspose = *self._append_unTrans
	endelse

	;Trim extra elements
	;   MrPutDimensionFirst may have removed shallow dimensions
	nDims = size(*self.data, /N_DIMENSIONS)
	case nDims of
		1: *self.data = (*self.data)[0:self._append_iLast]
		2: *self.data = (*self.data)[0:self._append_iLast,*]
		3: *self.data = (*self.data)[0:self._append_iLast,*,*]
		4: *self.data = (*self.data)[0:self._append_iLast,*,*,*]
		5: *self.data = (*self.data)[0:self._append_iLast,*,*,*,*]
		6: *self.data = (*self.data)[0:self._append_iLast,*,*,*,*,*]
		7: *self.data = (*self.data)[0:self._append_iLast,*,*,*,*,*,*]
		8: *self.data = (*self.data)[0:self._append_iLast,*,*,*,*,*,*,*]
	endcase

	;If BEFORE was set, we need to shift the elements to the front of the array
	if self._append_before then begin
		;Shift only the first dimension
		dimShift = lonarr(nDims)
		dimShift[0] = self._append_iLast + 1

		;Shift
		self -> Shift, dimShift
	endif

	;UnTranspose the data.
	;   Number of elements changed, so update the UNREFORM array.
	;   Be careful of leading shallow dimensions that are being re-introduced.
	iUnReform = min(where(unReform ne 1, count))
	if count eq 0 $
		then unReform[0]            = self._append_iLast + 1L $
		else unReform[iUnReform[0]] = self._append_iLast + 1L
	*self.data = transpose(reform(*self.data, unReform, /OVERWRITE), unTranspose)

	self._append_cont   = 0B
	self._append_before = 0B
end


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
; :Examples:
;  Multiple appends to the third dimension::
;       myArray  = MrVariable(12, 3, 6, 4, TYPE='FLOAT')
;       myArray -> Append_Multi, fltarr(12, 3, 4, 4), 3
;       myArray -> Append_Multi, fltarr(12, 3, 4, 4)
;       myArray -> Append_Multi, fltarr(12, 3, 4, 4)
;       myArray -> Append_Multi, /FINISH
;       help, myArray
;           MYARRAY     OBJREF      <ObjHeapVar5(MrVariable)>
;             ARRAY       FLOAT     = Array[12, 3, 18, 4]
;
;   Multiple appends to the beginning of an array::
;       myArray -> MrVariable(lonarr(2, 2))
;       myArray -> Append_Multi, lonarr(2,2) + 1, 1
;       myArray -> Append_Multi, lonarr(2,2) + 2
;       myArray -> Append_Multi, lonarr(2,2) + 3
;       myArray -> Append_Multi, /FINISH
;       help, myArray
;           MYARRAY     OBJREF      <ObjHeapVar17(MrVariable)>
;             ARRAY       LONG      = Array[8, 2]
;       print, myArray
;           1      1     2     2     3     3     0      0
;           1      1     2     2     3     3     0      0
;
; :Params:
;       DATA:          in, required, type=intarr
;                       Array to be concatenated to the current array. If `START` is set,
;                           then omit this parameter.
;       DIMENSION:      in, optional, type=integer, default=1
;                       The dimension, starting with 1, along which the two arrays are
;                           to be concatenated. For calls other than the first, this
;                           parameter is ignored.
;
; :Keywords:
;       BEFORE:         in, optional, type=boolean, default=0
;                       If set, `DATA` will be appended to the beginning of the implicit
;                           array. The default is to append to the end.
;       FINISH:         in, optional, type=boolean, default=0
;                       Set this keyword after the last call to ::Append_Multi to have
;                           extra, unfilled elements truncated and to reset the internal
;                           appending properties.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       If set, `DATA` will be copied directly into the object and
;                           will be left undefined.
;       NO_DEFORM:      in, optional, type=boolean, default=0
;                       During the appending process, the dimension being appended is
;                           made the first dimension, so that the array is deformed.
;                           To save time, the array is left in this state until the
;                           `FINISH` keyword is used. Set NO_DEFORM to transpose the array
;                           back to its original state after each append.
;-
pro MrVariable::Append_Multi, data, dimension, $
BEFORE=before, $
FINISH=finish, $
NO_COPY=no_copy, $
NO_DEFORM=no_deform
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		self._append_cont    = 0B
		self._append_before  = 0B
		self._append_no_copy = 0B
		if ptr_valid(pArray) then ptr_free, pArray
		MrPrintF, 'LogErr'
		return
	endif

;-----------------------------------------------------
; Finalize Append \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	if keyword_set(finish) then begin
		self -> Append_Finalize
		return
	endif

;-----------------------------------------------------
; SetUp to Append \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Dimensions
	if MrCmpVersion('8.3.1') le 0 $
		then nDims = size(self, /N_DIMENSIONS) $
		else nDims = size(*self.data, /N_DIMENSIONS)
	dims   = size(self, /DIMENSIONS)
	dimsIn = size(data, /DIMENSIONS)

	;First time?
	if self._append_cont eq 0 then begin
		;Defaults
		before    = keyword_set(before)
		no_copy   = keyword_set(no_copy)
		no_deform = keyword_set(no_deform)
		if n_elements(dimension) eq 0 then dimension = 1

		;Set up to start appending
		isFirst               = 1B
		self._append_before   = before
		self._append_cont     = 1B
		self._append_dim      = dimension
		self._append_no_copy  = no_copy
		self._append_putFirst = no_deform
		self._append_iLast    = nDims eq 0 ? -1L : (dims[dimension-1]-1 > 0)
	endif else begin
		isFirst               = 0B
	endelse

	;Dimension to which we will append
	dim = self._append_dim
	if isFirst eq 0 && self._append_putFirst eq 0 then dim = 1B

;-----------------------------------------------------
; Append \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Prepare the array being appended
	pArray  = ptr_new(data, NO_COPY=self._append_no_copy)
	*pArray = MrPutDimensionFirst(*pArray, self._append_dim, UNREFORM=unReform, $
	                              UNTRANSPOSE=unTranspose)

	;Ensure we know how to return to the original form
	if isFirst then begin
		*self._append_unReform = unReform
		*self._append_unTrans  = unTranspose
	endif

	;No data yet?
	if nDims eq 0 then begin
		*self.data = *pArray

	;Prep internal array for append
	endif else begin
		;Number of empty elements
		nTot   = dims[dim-1]
		nEmpty = nTot - self._append_iLast + 1
	
		;Extend the current array?
		if nEmpty le dimsIn[self._append_dim-1] then begin
			;Extend to fit or double number of elements, whichever is bigger
			nExtend = (2*dims[dim-1]) > dimsIn[dim-1]
			self -> Extend, nExtend, dim
		endif
	
		;Put the append dimension first
		if isFirst || self._append_putFirst $
			then *self.data = MrPutDimensionFirst(*self.data, self._append_dim)
	endelse

	;Append
	sIndex = self._append_iLast + 1L
	eIndex = sIndex + dimsIn[self._append_dim-1] - 1L
	case nDims of
		0: ;Do nothing
		1: (*self.data)[sIndex:eIndex] = temporary(*pArray)
		2: (*self.data)[sIndex:eIndex,*] = temporary(*pArray)
		3: (*self.data)[sIndex:eIndex,*,*] = temporary(*pArray)
		4: (*self.data)[sIndex:eIndex,*,*,*] = temporary(*pArray)
		5: (*self.data)[sIndex:eIndex,*,*,*,*] = temporary(*pArray)
		6: (*self.data)[sIndex:eIndex,*,*,*,*,*] = temporary(*pArray)
		7: (*self.data)[sIndex:eIndex,*,*,*,*,*,*] = temporary(*pArray)
		8: (*self.data)[sIndex:eIndex,*,*,*,*,*,*,*] = temporary(*pArray)
	endcase
	ptr_free, pArray

	;Update the last element filled and the number of empty elements
	self._append_iLast = eIndex

	;Untranspose the array each time?
	if self._append_putFirst $
		then *self.data = transpose(reform(*self.data, unReform, /OVERWRITE), unTranspose)
end


;+
;   Deterime if the input array is equal to the implicit variable array. This is a
;   wrapper method for IDL's Array_Equal function.
;
; :Params:
;       VALUE:      in, required, type=scalar/array
;                   The variable to be compared.
;
; :Keywords:
;       NO_TYPECONV:    in, optional, type=boolean, default=0
;                       If set, variable types will not be converted. The default is
;                           to convert least precise to most precise before comparing.
;
; :Returns:
;       RESULT:         out, required, type=boolean
;                       Returns true (1) if the variables are equal and false (0) otherwise.
;-
function MrVariable::Array_Equal, value, $
INDEX=index
	compile_opt idl2
	on_error, 2
	
	;Determine array equality
	if size(value, /TNAME) eq 'OBJREF' && obj_isa(value, 'MRVARIABLE') $
		then result = array_equal(*self.data, value['DATA'], NO_TYPECONV=no_typeconv) $
		else result = array_equal(*self.data, value,         NO_TYPECONV=no_typeconv)
	
	return, result
end


;+
;   Convert 1-Dimensional indices to their multi-dimensional counterparts.
;
; :Params:
;       INDEX:          in, required, type=intarr
;                       Array of 1D indices.
;
; :Returns:
;       IMULTI:         Multi-dimensional indices.
;-
function MrVariable::Array_Indices, index
	compile_opt idl2
	on_error, 2

	;Get the dimensions of the implicit array
	self -> GetProperty, DIMENSIONS=dims

	;Convert to multi-dimensional indices
	iMulti = array_indices(dims, index, /DIMENSIONS)

	return, iMulti
end


;+
;   Cache the variable for later use
;
; :Keywords:
;       NO_CLOBBER:     in, optional, type=boolean, default=0
;                       If set, and if the variable's name clashes with that of another
;                           variable, the input variable will be renamed by appending
;                           "_#" to its name, where # represents an integer. The default
;                           is to replace the object already in the container.
;       NAME_OUT:       out, optional, type=string
;                       If `NO_CLOBBER` is set, then this is the new name of `VARIABLE`
;                           after it is renamed. If `VARIABLE` is not renamed, the
;                           original name is returned.
;-
pro MrVariable::Cache, $
NO_CLOBBER=no_clobber, $
NAME_OUT=name_out
	compile_opt idl2
	on_error, 2

	;Setup caching store
	@mrvar_common
	
	;Create a cache in which to store MrVariables
	if ~obj_valid(MrVarCache) then MrVarCache = MrVariable_Cache()
	
	;Add the array to the end of the container
	MrVarCache -> Add, self, NO_CLOBBER=no_clobber, NAME_OUT=name_out
end


;+
;   Calculate the cross-covariance with another variable as a function of lag.
;
; :Params:
;       VAR:            in, required, type=int/string/objref
;                       The name, number, or reference of a MrVariable object.
;       LAG:            in, required, type=intarr
;                       A scalar or n-element integer vector in the interval
;                           [-(n-2), (n-2)], specifying the signed distances between
;                           indexed elements of the implicit array.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the ouput variable is added to the variable cache.
;       COVARIANCE:     in, optional, type=boolean, default=0
;                       If set, the sample covariance will be computed.
;       DOUBLE:         in, optional, type=string
;                       If set, computations are done in double precision.
;       NAME:           in, optional, type=string, default='C_Correlate('+self.name+','+var.name+')'
;                       Name given to the output variable.
;
; :Returns:
;       OCCOR:          out, required, type=objref
;                       MrVariable object with the lagged cross-correlation.
;-
FUNCTION MrVariable::C_Correlate, var, lag, $
CACHE=cache, $
COVARIANCE=covariance, $
DOUBLE=double, $
NAME=name
	Compile_Opt idl2
	On_Error, 2

	;Get the variable
	oVar = MrVar_Get(var)
	
	;Defaults
	tf_covariance = Keyword_Set(covariance)
	IF N_Elements(name) EQ 0 THEN name = (tf_covariance ? 'Covariance(' : 'C_Correlate(') + self.name + ',' + oVar.name + ')'
	
	;Cross-correlation
	result = C_Correlate( *self.data, oVar['DATA'], lag, $
	                      COVARIANCE = covariance, $
	                      DOUBLE     = double )
	
	;Dependent variable
	oDep0 = MrVariable( lag, $
	                    NAME = name + '_lag' )
	
	;Cross-Correlation
	oCCor = MrVariable( result, $
	                    CACHE = cache, $
	                    NAME  = name, $
	                    /NO_COPY )
	
	;Attributes
	oDep0['TITLE']    = 'Lag'
	oCCor['DEPEND_0'] = oDep0
	oCCor['TITLE']    = tf_covariance $
	                        ? 'Covariance' $
	                        : 'Cross-Correlation'
	
	;Done!
	RETURN, oCCor
END


;+
;   The prupose of this method is to erase the data from the array. Alternatively, try
;       myArray.array = !Null
;       myArray[*] = !Null
;-
PRO MrVariable::Clear
	Compile_Opt idl2
	On_Error, 2

	*self.data = !Null
END


;+
;   Create a copy of the variable with a new heap identifier. Both data and
;   attributes are copied.
;
; :Params:
;       NAME:           in, optional, type=string, default=<name>+'_copy'
;                       Name of the variable copy.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the copy will be added to the variable cache
;-
function MrVariable::Copy, name, $
CACHE=cache
	compile_opt idl2
	on_error, 2

	;Defaults
	cache = keyword_set(cache)
	if n_elements(name) eq 0 then name = self.name + '_copy'
	
	;Copy the data into a new object
	;   - Use Obj_Class() so subclasses can inherit the method.
	theCopy = obj_new(obj_class(self), *self.data, $
	                  NAME = name)

	;Copy the variable attributes
	self -> CopyAttrTo, theCopy
	
	;Cache the variable
	if cache then theCopy -> Cache, /NO_CLOBBER
	return, theCopy
end


;+
;   Copy variable attributes to another MrVariable object.
;
; :Params:
;       VAR:            in, required, type=objref
;                       The destination MrVariable object.
;       ATTRNAMES:      in, optional, type=string/strarr, default=all attrs
;                       Names of the attributes to be copied.
;
; :Keywords:
;       OVERWRITE:      in, optional, type=boolean, default=1
;                       If set, then existing attributes with names appearing in
;                           `ATTRNAMES` will have their values overwritten. If set
;                           to zero, issue a warning and skip the attribute.
;-
pro MrVariable::CopyAttrTo, var, attrNames, $
OVERWRITE=overwrite
	compile_opt idl2
	on_error, 2
	
	;VAR must be a MrVariable
	if ~obj_isa(var, 'MRVARIABLE') then message, 'VAR must be a MrVariable object.'
	tf_overwrite = n_elements(overwrite) eq 0 ? 1B : keyword_set(overwrite)
	
	;Default to copying all attributes
	if n_elements(attrNames) eq 0 $
		then attrNames = self -> GetAttrNames(COUNT=nAttrs) $
		else nAttrs    = n_elements(attrNames)

	;Copy the variable attributes
	for i = 0, nAttrs-1 do var -> AddAttr, attrNames[i], self[attrNames[i]], OVERWRITE=tf_overwrite
end


;+
;   Extend a dimension of the implicit array by appending the desired number of
;   elements.
;
; :Params:
;       NELEMENTS:      in, optional, type=intarr, default=2xCurrent number of elements
;                       Number of elements to extend the `DIMENSION` of the implicit
;                           array.
;       DIMENSION:      in, optional, type=integer, default=1
;                       The dimension, starting with 1, to be extended.
;
; :Keywords:
;       BEFORE:         in, optional, type=boolean, default=0
;                       If set, `NELEMENTS` will be appended to the beginning of the
;                           implicit array. The default is to append to the end.
;       VALUE:          in, optional, type=any
;                       The value that the extended elements will take on.
;-
pro MrVariable::Extend, nElements, dimension, $
BEFORE=before, $
VALUE=value
	compile_opt idl2
	on_error, 2

	;Defaults
	before = keyword_set(before)
	if n_elements(nElements) eq 0 then nElements = n_elements(self)
	if n_elements(dimension) eq 0 then dimension = 1

	;Create the 
	self -> GetProperty, DIMENSIONS=ExtendDims, TYPE=type
	ExtendDims[dimension-1] = nElements
	ExtendArray = Make_Array(DIMENSION=ExtendDims, TYPE=type, VALUE=value)

	;Concatenate
	if before eq 1 $
		then *self.data = MrConcatenate(ExtendArray, *self.data, dimension) $
		else *self.data = MrConcatenate(*self.data, ExtendArray, dimension)
end


;+
;   Prepare the variable data to be written to an ASCII file. The output structure
;   can be passed to MrVar_ReadASCII via the _STRICT_EXTRA keyword.
;
; :Params:
;       THECDF:         in, optional, type=long/string/objref
;                       The name or CDF ID of the CDF file to which data is written,
;                           or the MrCDF_File object containing the file information.
;
; :Keywords:
;       CREATE:         in, optional, type=boolean
;                       If set and the variable does not exist in the file, it is created.
;                           The default is to check the file for a variable with the same
;                           name and set the keyword accordingly. If set a variable by the
;                           same name already exists in `THECDF`, an error will occur.
;       CDF_TYPE:       in, optional, type=string
;                       The CDF datatype of the variable. The default determined automatically
;                           from the IDL datatype.
;       TEST:           in, optional, type=boolean, default=0
;                       If set, check if the variable already exists in the file. If it
;                           does, return without doing anything. If not, set `CREATE`=1.
;                           This is useful for, e.g., variables pointed to be the DEPEND_#
;                           variables attribute, which can be shared among several other
;                           variables. Prevents writing multiple times.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrCDF_File::WriteVar is also accepted here.
;
; :Output:
;       ASCII:          out, required, type=struct
;                       A structure of keyword-value pairs that can be passed directly
;                           to MrVar_ReadASCII.
;-
PRO MrVariable::ExportToASCII, file, $
HEADER=header, $
FORMAT=fmt, $
CLOBBER=clobber, $
_REF_EXTRA=extra
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF Size(file, /TNAME) EQ 'STRING' && (file_info(file)).exists $
			THEN IF N_Elements(lun) GT 0 THEN Close, lun
		RETURN
	END
	
	tf_open    = 0
	tf_clobber = Keyword_Set(clobber)
	tf_write   = N_Elements(file) GT 0

;-------------------------------------------
; Open the File ////////////////////////////
;-------------------------------------------
	
	;Open the file
	IF tf_write THEN BEGIN
		;File name was given
		IF MrIsA(file, 'STRING') THEN BEGIN
			;Do not clobber
			IF ~tf_clobber && File_Test(file) $
				THEN Message, 'File exists. Set CLOBBER keyword to overwrite.' $
				ELSE OpenW, lun, file, /GET_LUN
			
		;Logical unit number was given
		ENDIF ELSE IF MrIsA(file, /INTEGER) THEN BEGIN
			stat = FStat(file)
			IF stat.name EQ '' || stat.open EQ 0 $
				THEN Message, 'FILE is not a valid logical unit number.' $
				ELSE lun = file
			
		;Something else was given
		ENDIF ELSE BEGIN
			Message, 'FILE must be a file name or a logical unit number (LUN).'
		ENDELSE
	ENDIF

;-------------------------------------------
; Data Structure ///////////////////////////
;-------------------------------------------
	dims  = Size(self, /DIMENSIONS)
	nDims = Size(self, /N_DIMENSIONS)
	
	;If there are DEPEND_1-3 attributes
	IF ~Array_Equal(self -> HasAttr('DEPEND_' + ['1', '2', '3']), 0) $
		THEN Message, 'I do not know how to export variables with DEPEND_[1,2,3] attributes.'
	
	;Should be equivalent to the DEPEND_1-3 clause
	IF nDims GT 2 $
		THEN Message, 'I do not know how to export variables with > 2 dimensions.'
	
	;Cread a data structure containing information on how the data
	;should be written to the ASCII file.
;	out = { data:   Ptr_New(/ALLOCATE_HEAP), $
;	        header: '', $
;	        fmt:    '' }
	
	;If there are more than one dimensions, expand the structure into
	;and array, with one element per row, etc.
;	IF nDims EQ 2 THEN BEGIN
;		out = Replicate(out, dims[1])
;		FOR i = 0, dims[1]-1 DO *out[i].data = self['DATA',*,i]
;	ENDIF ELSE BEGIN
;		*out.data = self['DATA']
;	ENDELSE
	
;-------------------------------------------
; Format Code //////////////////////////////
;-------------------------------------------
	IF self -> HasAttr('FORMAT') THEN BEGIN
		fmt = self['FORMAT']
	ENDIF ELSE BEGIN
		sample = self['DATA',0]
		CASE 1 OF
			MrIsA(sample, 'STRING'): fmt = 'a' + String(Max(StrLen(self['DATA'])), FORMAT='(i0)')
			MrIsA(sample, /INTEGER): fmt = 'i' + String(Ceil(ALog10(Max((self['DATA']))))+1, FORMAT='(i0)')
			MrIsA(sample, /FLOAT, /REAL): BEGIN
				
				;Range of numbers
				pmax = ALog10(Abs(Max(self['DATA'], MIN=pmin)))
				pmin = ALog10(Abs(pmin))
				pmax = pmax GT 0 ? Ceil(pmax) : Floor(pmax)
				pmin = pmin GT 0 ? Ceil(pmin) : Floor(pmin)
				
				;Width of numbers
				;   - Float:  Print allows 6 digits, 8 with /IMPLIED_PRINT
				;   - Double: Print allows 8 digits, 16 with /IMPLIED_PRINT
				;   - We will use 6 for float and 8 for double
				IF Size(sample, /TNAME) EQ 'DOUBLE' THEN BEGIN
					;Use scientific notation if numbers are too big/small or large dynamic range
					CASE 1 OF
						(pmax - pmin) GT 7: code = 'e'
						pmax GT 8:          code = 'e'
						pmax LT -7:         code = 'e'
						pmin LT -7:         code = 'e'
						ELSE:               code = 'f'
					ENDCASE
					width = code EQ 'e' ? 14 : 10
					dec   = code EQ 'e' ?  7 : (8-pmax) > 0
					fmt   = String( code, width, dec, FORMAT='(%"%1s%0i.%0i")')
				
				;Single precision
				ENDIF ELSE BEGIN
					;Use scientific notation if numbers are too big/small or large dynamic range
					CASE 1 OF
						(pmax - pmin) GT 5: code = 'e'
						pmax GT 6:          code = 'e'
						pmax LT -5:         code = 'e'
						pmin LT -5:         code = 'e'
						ELSE:               code = 'f'
					ENDCASE
					width = code EQ 'e' ? 12 : 8
					dec   = code EQ 'e' ?  5 : (6-pmax) > 0
					fmt   = String( code, width, dec, FORMAT='(%"%1s%0i.%0i")')
				ENDELSE
			ENDCASE
		ENDCASE
	ENDELSE
	
	;Write each column the same as the first
	IF nDims GT 2 THEN fmt = String(dims[1], FORMAT='(i0)') + '('+fmt+',1x)'

;-------------------------------------------
; Header ///////////////////////////////////
;-------------------------------------------
	;Width of 
;	width = StRegEx(fmt, '^[aefi]([0-9]+)', /EXTRACT, /SUBEXP)
;	width = Fix(width[1])
	
	
;	IF self -> HasAttr('LABLAXIS') THEN BEGIN
;		header = self['LABLAXIS']
;	ENDIF ELSE IF self -> HasAttr('LABEL') THEN BEGIN
;		header = self['LABEL']
;		IF N_Elements(label) NE dims[1] THEN header = ''
;	ENDIF ELSE IF self -> HasAttr('LABL_PTR_1') THEN BEGIN
;		header = (self['LABL_PTR_1'])['DATA'] 
;	ENDIF ELSE BEGIN
;		header = ''
;	ENDELSE
	
;	nHeader = N_Elements(header)
;	header  = 'V' + String(IndGen(nHeader), FORMAT='(i0)') + '_' + header
	
;-------------------------------------------
; DEPEND_0 /////////////////////////////////
;-------------------------------------------
	tf_dep0 = self -> HasAttr('DEPEND_0')
	IF tf_dep0 THEN BEGIN
		self['DEPEND_0'] -> ExportToASCII, FORMAT=dep_fmt, HEADER=dep_hdr
		fmt = dep_fmt + ', 1x, ' + fmt
	ENDIF
	
;-------------------------------------------
; Write the Data ///////////////////////////
;-------------------------------------------
	IF tf_write THEN BEGIN
		header    = StrArr(7)
		nHeader   = N_Elements(header)
		header[0] = String(dims[0]+nHead, FORMAT='(%"Lines:         %i")')
		header[1] = String(dims[0]+nHead, FORMAT='(%"Lines:         %i")')
		header[2] = String(6,             FORMAT='(%"Header:        %i")')
		header[3] = String(fmt,           FORMAT='(%"Format Code:   \"%s\"")')
		header[4] = String(self.name,     FORMAT='(%"Variable Name: \"%s\"")')
		header[5] = String('',            FORMAT='(%"Columns:       [%s]")')  ;StrJoin(labels, ', ')
		header[6] = 'DATA:'
		
		oDep0 = self['DEPEND_0']
		IF tf_dep0 $
			THEN FOR i = 0, dims[0]-1 DO PrintF, lun, oDep0['DATA',i], self['DATA',i,*], FORMAT='('+fmt+')' $
			ELSE FOR i = 0, dims[0]-1 DO PrintF, lun, self['DATA',i], FORMAT='('+fmt+')'
		
		Close, lun
	ENDIF
END

;+
;   Write data to a CDF file.
;
;   NOTE:
;       Requires the MrCDF library.
;       https://github.com/argallmr/IDLcdf
;
; :Params:
;       THECDF:         in, optional, type=long/string/objref
;                       The name or CDF ID of the CDF file to which data is written,
;                           or the MrCDF_File object containing the file information.
;
; :Keywords:
;       CREATE:         in, optional, type=boolean
;                       If set and the variable does not exist in the file, it is created.
;                           The default is to check the file for a variable with the same
;                           name and set the keyword accordingly. If set a variable by the
;                           same name already exists in `THECDF`, an error will occur.
;       CDF_TYPE:       in, optional, type=string
;                       The CDF datatype of the variable. The default determined automatically
;                           from the IDL datatype.
;       TEST:           in, optional, type=boolean, default=0
;                       If set, check if the variable already exists in the file. If it
;                           does, return without doing anything. If not, set `CREATE`=1.
;                           This is useful for, e.g., variables pointed to be the DEPEND_#
;                           variables attribute, which can be shared among several other
;                           variables. Prevents writing multiple times.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrCDF_File::WriteVar is also accepted here.
;-
PRO MrVariable::ExportToCDF, theCDF, $
CREATE=create, $
CDF_TYPE=cdf_type, $
TEST=test, $
_REF_EXTRA=extra
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN
	END
	
	;Get a MrCDF_File object
	CASE 1 OF
		MrIsA(theCDF, 'MrCDF_File'): oCDF = theCDF
		MrIsA(theCDF, 'STRING'):     oCDF = MrCDF_File(theCDF)
		MrIsA(theCDF, /NUMBER):      oCDF = MrCDF_File(theCDF)
		ELSE: Message, 'Invalid value for THECDF.'
	ENDCASE

	;Defaults
	IF Keyword_Set(test) THEN BEGIN
		IF oCDF -> HasVar(self.name) $
			THEN RETURN $
			ELSE tf_create = 1B
	ENDIF
	tf_create = N_Elements(create) EQ 0 ? ~oCDF -> HasVar(self.name) : Keyword_Set(create)
	IF self -> HasAttr('CDF_TYPE') THEN cdf_type = self['CDF_TYPE']
	
	;Write the data to the file
	;   - Variable dimensions are [NRECS, DEP1, DEP2, DEP3]
	;   - CDF dimensions are [DEP3, DEP2, DEP1, NRECS]
	oCDF -> WriteVar, self.name, (N_Elements(self) EQ 1 ? self['DATA'] : Transpose(self['DATA']) ), $
	                  CREATE        = tf_create, $
	                  CDF_TYPE      = cdf_type, $
	                  _STRICT_EXTRA = extra

;-------------------------------------------
; Variable Attributes //////////////////////
;-------------------------------------------
	self -> ExportToCDF_Attrs, oCDF
	
;-------------------------------------------
; Clean Up /////////////////////////////////
;-------------------------------------------
	;If the CDF object was created internally, destroy it to close the file.
	IF Size(theCDF, /TNAME) NE 'OBJREF' THEN Obj_Destroy, oCDF
END


;+
;   Write variable attributes to a CDF file.
;
;   NOTE:
;       Requires the MrCDF library.
;       https://github.com/argallmr/IDLcdf
;
; :Private:
;
; :Params:
;       OCDF:           in, optional, type=objref
;                       The MrCDF_File object containing the file information.
;-
PRO MrVariable::ExportToCDF_Attrs, oCDF
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	
	;CDF Attributes
	;   - Includes graphics keywords that serve as substitutes for CDF attributes
	cdfAttrNames = [ 'CATDESC',         'DEPEND_0',   'DEPEND_1',  'DEPEND_2',   'DEPEND_3', 'DISPLAY_TYPE', $
	                 'FIELDNAME',       'FILLVAL',    'FORMAT',    'FMT_PTR',    'LABLAXIS', 'LABL_PTR_1', $
	                 'LABL_PTR_2',      'LABL_PTR_3', 'UNITS',     'UNIT_PTR',   'VALIDMIN', 'VALIDMAX', $
	                 'VARTYPE',         'SCALETYPE',  'SCAL_PTR',  'VAR_NOTES',  'AVG_TYPE', 'DELTA_PLUS_VAR', $
	                 'DELTA_MINUS_VAR', 'DICT_KEY',   'MONOTON',   'SCALEMIN',   'SCALEMAX', 'V_PARENT', $
	                 'DERIVN',          'sig_digits', 'SI_conv',   'PLOT_TITLE', 'TITLE',    'AXIS_RANGE', $
	                 'LOG',             'MIN_VALUE',  'MAX_VALUE', 'LABEL' ]
	
	;Variable attributes
	attrNames = self -> GetAttrNames(COUNT=nAttrs)
	IF nAttrs EQ 0 THEN RETURN

	;Find which variable attributes are NOT CDF attributes
	!Null = MrIsMember(cdfAttrNames, attrNames, COMPLEMENT=iVarAttrs, NCOMPLEMENT=nVarAttrs)

;-------------------------------------------
; Attributes with IDL Alternatives /////////
;-------------------------------------------
	;
	; Set the /CREATE flag in case the attribute does not yet exist in the file.
	;
	
	;CATDESC
	;   - Substitutes: PLOT_TITLE, TITLE
	IF self -> HasAttr('CATDESC') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'CATDESC', self['CATDESC'], /CREATE
	ENDIF ELSE IF self -> HasAttr('PLOT_TITLE') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'CATDESC', self['PLOT_TITLE'], /CREATE
	ENDIF ELSE IF self -> HasAttr('TITLE') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'CATDESC', self['TITLE'], /CREATE
	ENDIF
	
	;FIELDNAM
	;   - Substitutes: TITLE
	IF self -> HasAttr('FIELDNAM') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'FIELDNAM', self['FIELDNAM'], /CREATE
	ENDIF ELSE IF self -> HasAttr('TITLE') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'FIELDNAM', self['TITLE'], /CREATE
	ENDIF
	
	;LABLAXIS
	;   - Substitutes: LABEL (scalar)
	IF self -> HasAttr('LABLAXIS') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'LABLAXIS', self['LABLAXIS'], /CREATE
	ENDIF ELSE IF self -> HasAttr('LABEL') THEN BEGIN
		IF MrIsA(self['LABEL'], /SCALAR) THEN BEGIN
			oCDF -> WriteVarAttr, self.name, 'LABLAXIS', self['LABEL'], /CREATE
		ENDIF
	ENDIF
	
	;LABL_PTR_1
	;   - Substitutes: LABEL (array)
	;   - CDF_AttPut accepts scalar values only
	;   - TODO: Allow LABL_PTR_# attribute to be a MrVariable. This will allow
	;           multiple variables with the same labels to share the same LABL_PTR_#
	;           CDF variable.
	IF self -> HasAttr('LABL_PTR_1') && IsA(self['LABL_PTR_1'], 'OBJREF') THEN BEGIN
		;Write to file
		oCDF -> WriteVarAttr, self.name, 'LABL_PTR_1', (self['LABL_PTR_1']).name, /CREATE
		self['LABL_PTR_1'] -> ExportToCDF, oCDF, /TEST
		
	ENDIF ELSE BEGIN
		attrName = self -> HasAttr('LABL_PTR_1') ? 'LABL_PTR_1' : $
		           self -> HasAttr('LABEL')      ? 'LABEL'      : $
		           ''
		
		;Array of values
		IF attrName NE '' THEN BEGIN
			;Create a variable
			labl_ptr_vname = self.name + '_labl_ptr_1'
			oLabel         = MrVariable( self[attrName], NAME=labl_ptr_vname, /NO_COPY )
			
			;Add attributes
			oLabel['CATDESC']  = 'Labels for ' + self.name
			oLabel['FIELDNAM'] = 'Labels'
			oLabel['VARTYPE']  = 'metadata'
			
			;Write to file
			oCDF   -> WriteVarAttr, self.name, 'LABL_PTR_1', labl_ptr_vname, /CREATE
			oLabel -> ExportToCDF,  oCDF
		ENDIF
	ENDELSE
	
	;SCALEMAX
	;   - Substitutes: AXIS_RANGE
	IF self -> HasAttr('SCALEMAX') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'SCALEMAX', self['SCALEMAX'], /CREATE
	ENDIF ELSE IF self -> HasAttr('AXIS_RANGE') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'SCALEMAX', (self['AXIS_RANGE'])[1], /CREATE
	ENDIF
	
	;SCALEMIN
	;   - Substitutes: AXIS_RANGE
	IF self -> HasAttr('SCALEMIN') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'SCALEMIN', self['SCALEMIN'], /CREATE
	ENDIF ELSE IF self -> HasAttr('AXIS_RANGE') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'SCALEMIN', (self['AXIS_RANGE'])[0], /CREATE
	ENDIF
	
	;SCALETYP
	;   - Substitutes: LOG
	IF self -> HasAttr('SCALETYP') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'SCALETYP', self['SCALETYP'], /CREATE
	ENDIF ELSE IF self -> HasAttr('SCALETYP') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'SCALETYP', (self['LOG'] EQ 1B ? 'log' : 'linear'), /CREATE
	ENDIF
	
	;VALIDMAX
	;   - Substitutes: MAX_VALUE
	IF self -> HasAttr('VALIDMAX') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'VALIDMAX', self['VALIDMAX'], /CREATE
	ENDIF ELSE IF self -> HasAttr('MAX_VALUE') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'VALIDMAX', self['MAX_VALUE'], /CREATE
	ENDIF
	
	;VALIDMIN
	;   - Substitutes: MIN_VALUE
	IF self -> HasAttr('VALIDMIN') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'VALIDMIN', self['VALIDMIN'], /CREATE
	ENDIF ELSE IF self -> HasAttr('MIN_VALUE') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'VALIDMIN', self['MIN_VALUE'], /CREATE
	ENDIF

;-------------------------------------------
; Other CDF Attributes /////////////////////
;-------------------------------------------
	;Missing:
	;   FMT_PTR
	;   UNIT_PTR
	;   SCAL_PTR
	
	;DEPEND_0
	IF self -> HasAttr('DEPEND_0') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'DEPEND_0', (self['DEPEND_0']).name, /CREATE
		self['DEPEND_0'] -> ExportToCDF, oCDF, /TEST
	ENDIF
	
	;DEPEND_1
	IF self -> HasAttr('DEPEND_1') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'DEPEND_1', (self['DEPEND_1']).name, /CREATE
		self['DEPEND_1'] -> ExportToCDF, oCDF, /TEST
	ENDIF
	
	;DEPEND_2
	IF self -> HasAttr('DEPEND_2') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'DEPEND_2', (self['DEPEND_2']).name, /CREATE
		self['DEPEND_2'] -> ExportToCDF, oCDF, /TEST
	ENDIF
	
	;DEPEND_3
	IF self -> HasAttr('DEPEND_3') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'DEPEND_3', (self['DEPEND_3']).name, /CREATE
		self['DEPEND_3'] -> ExportToCDF, oCDF, /TEST
	ENDIF
	
	;DELTA_PLUS_VAR
	IF self -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'DELTA_PLUS_VAR', (self['DELTA_PLUS_VAR']).name, /CREATE
		self['DELTA_PLUS_VAR'] -> ExportToCDF, oCDF, /TEST
	ENDIF
	
	;DELTA_MINUS_VAR
	IF self -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		oCDF -> WriteVarAttr, self.name, 'DELTA_MINUS_VAR', (self['DELTA_MINUS_VAR']).name, /CREATE
		self['DELTA_MINUS_VAR'] -> ExportToCDF, oCDF, /TEST
	ENDIF
	
	;LABL_PTR_2
	IF self -> HasAttr('LABL_PTR_2') THEN BEGIN
		;Create the variable
		labl_ptr_vname = self.name + '_labl_ptr_2'
		oLabel         = MrVariable( self['LABL_PTR_2'], NAME=labl_ptr_vname )
		
		;Add attributes
		oLabel['CATDESC']  = 'Labels for ' + self.name
		oLabel['FIELDNAM'] = 'Labels'
		oLabel['VARTYPE']  = 'metadata'
		
		;Write to file
		oCDF   -> WriteVarAttr, self.name, 'LABL_PTR_2', labl_ptr_vname, /CREATE
		oLabel -> ExportToCDF, oCDF, /TEST
	ENDIF
	
	;LABL_PTR_3
	IF self -> HasAttr('LABL_PTR_3') THEN BEGIN
		;Create the variable
		labl_ptr_vname = self.name + '_labl_ptr_3'
		oLabel         = MrVariable( self['LABL_PTR_3'], NAME=labl_ptr_vname )
		
		;Add attributes
		oLabel['CATDESC']  = 'Labels for ' + self.name
		oLabel['FIELDNAM'] = 'Labels'
		oLabel['VARTYPE']  = 'metadata'
		
		;Write to file
		oCDF   -> WriteVarAttr, self.name, 'LABL_PTR_3', labl_ptr_vname, /CREATE
		oLabel -> ExportToCDF, oCDF, /TEST
	ENDIF
	
	;Non-pointers
	IF self -> HasAttr('DELTA_MINUS')  THEN oCDF -> WriteVarAttr, self.name, 'DELTA_MINUS',  self['DELTA_MINUS'],  /CREATE
	IF self -> HasAttr('DELTA_PLUS')   THEN oCDF -> WriteVarAttr, self.name, 'DELTA_PLUS',   self['DELTA_PLUS'],   /CREATE
	IF self -> HasAttr('DISPLAY_TYPE') THEN oCDF -> WriteVarAttr, self.name, 'DISPLAY_TYPE', self['DISPLAY_TYPE'], /CREATE
	IF self -> HasAttr('FILLVAL')      THEN oCDF -> WriteVarAttr, self.name, 'FILLVAL',      self['FILLVAL'],      /CREATE
	IF self -> HasAttr('FORMAT')       THEN oCDF -> WriteVarAttr, self.name, 'FORMAT',       self['FORMAT'],       /CREATE
	IF self -> HasAttr('UNITS')        THEN oCDF -> WriteVarAttr, self.name, 'UNITS',        self['UNITS'],        /CREATE
	IF self -> HasAttr('VARTYPE')      THEN oCDF -> WriteVarAttr, self.name, 'VARTYPE',      self['VARTYPE'],      /CREATE
	IF self -> HasAttr('VAR_NOTES')    THEN oCDF -> WriteVarAttr, self.name, 'VAR_NOTES',    self['VAR_NOTES'],    /CREATE
	IF self -> HasAttr('AVG_TYPE')     THEN oCDF -> WriteVarAttr, self.name, 'AVG_TYPE',     self['AVG_TYPE'],     /CREATE
	IF self -> HasAttr('DICT_KEY')     THEN oCDF -> WriteVarAttr, self.name, 'DICT_KEY',     self['DICT_KEY'],     /CREATE
	IF self -> HasAttr('MONOTON')      THEN oCDF -> WriteVarAttr, self.name, 'MONOTON',      self['MONOTON'],      /CREATE
	IF self -> HasAttr('V_PARENT')     THEN oCDF -> WriteVarAttr, self.name, 'V_PARENT',     self['V_PARENT'],     /CREATE
	IF self -> HasAttr('DERIVN')       THEN oCDF -> WriteVarAttr, self.name, 'DERIVN',       self['DERIVN'],       /CREATE
	IF self -> HasAttr('sig_digits')   THEN oCDF -> WriteVarAttr, self.name, 'sig_digits',   self['sig_digits'],   /CREATE
	IF self -> HasAttr('SI_conv')      THEN oCDF -> WriteVarAttr, self.name, 'SI_conv',      self['SI_conv'],      /CREATE

;-------------------------------------------
; Variable Attributes //////////////////////
;-------------------------------------------
	FOR i = 0, nVarAttrs - 1 DO BEGIN
		attrName = attrNames[iVarAttrs[i]]
		oCDF -> WriteVarAttr, self.name, attrName, self[attrName], /CREATE
	ENDFOR
END


;+
;   Get attribute names.
;
; :Keywords:
;       COUNT:          in, optional, type=integer
;                       Number of attribute names returned.
;
; :Returns:
;       ATTRNAMES:      Names of all variable attributes.
;-
function MrVariable::GetAttrNames, $
COUNT=count
	compile_opt idl2
	on_error, 2

	;Return null for non-existant attributes
	attrList  = self.attributes -> Keys()
	attrNames = attrList -> ToArray(/NO_COPY)
	count     = n_elements(attrNames)

	return, attrNames
end


;+
;   Get an attribute value.
;
; :Params:
;       ATTRNAME:       in, required, type=string
;                       Name of the attribute for which the value is retreived.
;
; :Keywords:
;       NULL:           in, optional, type=boolean
;                       If set, and `ATTRNAME` is not a valid attribute name,
;                           then silently return !Null. The default is to throw
;                           an error.
;
; :Returns:
;       ATTRVALUE:      The attribute value.
;-
function MrVariable::GetAttrValue, attrName, $
NULL=null
	compile_opt idl2
	on_error, 2

	;Check if the attribute exists
	tf_exists = self.attributes -> HasKey(attrName)
	
	;If it exists
	if tf_exists then begin
		value = self.attributes[attrName]
		
	;Otherwise
	endif else begin
		if keyword_set(null) $
			then value = !Null $
			else message, 'Attribute does not exist: "' + attrName + '".'
	endelse

	return, value
end


;+
;   Get class properties.
;
; :Examples:
;   Get the array::
;       myArray  = MrVariable(11, 91, 6)
;       theArray = myArray -> GetArray()
;       help, theArray
;           THEARRAY        FLOAT     = Array[11, 91, 6]
;
; :Keywords:
;       DESTROY:            in, optional, type=boolean, default=0
;                           If set, the object will be destroyed after returning the array.
;       NO_COPY:            in, optional, type=boolean, default=0
;                           If set, the array will be removed frome the object and return.
;                               This prevents two copies of the array from being in
;                               memory.
;
; :Returns:
;       DATA:               out, optional, type=any
;                           The variable's data.
;-
function MrVariable::GetData, $
DESTROY=destroy, $
NO_COPY=no_copy
	compile_opt idl2
	on_error, 2

	;Destroy the data and the object?
	if keyword_set(destroy) then begin
		;Get the data
		;   - The pointer will be cleaned upon destroy
		data = temporary(*self.data)
		
		;Destroy the object
		obj_destroy, self
	
	;Free memory
	endif else if keyword_set(no_copy) then begin
		data = temporary(*self.data)
		ptr_free, self.data
		self.data = ptr_new(/ALLOCATE_HEAP)
	
	;Return a copy
	endif else begin
		return, *self.data
	endelse
	
	;Return the array
	return, data
end


;+
;   Return the name of the array.
;
; :Returns:
;       NAME:       Name of the object.
;-
function MrVariable::GetName
	return, self.name
end


;+
;   Extract a subarray.
;
; :Parmas:
;       BOUNDS:             in, optional, type=string, default='*'
;                           The array bounds of the data to be returned.
;
; :Keywords:
;       CACHE:              in, optional, type=boolean, default=0
;                           The array bounds of the data to be returned.
;       NAME:               in, optional, type=string, default=self.name + '_subarray'
;                           Name to be given to the output subarray.
;
; :Returns:
;       VOUT:               out, required, type=object
;                           A MrVariable object reference containing the subarray.
;-
FUNCTION MrVariable::GetSubArray, bounds, $
DATA=data, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	On_Error, 2
	
	;Defaults
	tf_data = Keyword_Set(data)
	IF N_Elements(bounds) EQ 0 THEN bounds = ''
	IF N_Elements(name)   EQ 0 THEN name = 'SubArray(' + self.name + ',' + bounds + ')'
	
	;Return the array as is
	IF bounds EQ '' THEN BEGIN
		IF tf_data $
			THEN RETURN, *self.data $
			ELSE RETURN, self -> Copy(name, CACHE=cache)
	ENDIF
	
;-------------------------------------------
; Array Bounds /////////////////////////////
;-------------------------------------------
	
	;Dimensions of the implicit array
	inDims  = Size( *self.data, /DIMENSIONS)
	nInDims = Size( *self.data, /N_DIMENSIONS)
	
	;Get the array indices
	inds = MrArray_Bounds( inDims, bounds, /DIMENSIONS )
	
	;Number of dimensions in the output array
	;   - Indices are [nDimensions, 3]
	nOutDims = Size(inds, /DIMENSIONS)
	nOutDims = nOutDims[0]
	
;-------------------------------------------
; Extract Data /////////////////////////////
;-------------------------------------------
	CASE nOutDims OF
		1: data = (*self.data)[ inds[0,0]:inds[0,1]:inds[0,2] ]
		2: data = (*self.data)[ inds[0,0]:inds[0,1]:inds[0,2], $
		                        inds[1,0]:inds[1,1]:inds[1,2] ]
		3: data = (*self.data)[ inds[0,0]:inds[0,1]:inds[0,2], $
		                        inds[1,0]:inds[1,1]:inds[1,2], $
		                        inds[2,0]:inds[2,1]:inds[2,2] ]
		4: data = (*self.data)[ inds[0,0]:inds[0,1]:inds[0,2], $
		                        inds[1,0]:inds[1,1]:inds[1,2], $
		                        inds[2,0]:inds[2,1]:inds[2,2], $
		                        inds[3,0]:inds[3,1]:inds[3,2] ]
		5: data = (*self.data)[ inds[0,0]:inds[0,1]:inds[0,2], $
		                        inds[1,0]:inds[1,1]:inds[1,2], $
		                        inds[2,0]:inds[2,1]:inds[2,2], $
		                        inds[3,0]:inds[3,1]:inds[3,2], $
		                        inds[4,0]:inds[4,1]:inds[4,2] ]
		6: data = (*self.data)[ inds[0,0]:inds[0,1]:inds[0,2], $
		                        inds[1,0]:inds[1,1]:inds[1,2], $
		                        inds[2,0]:inds[2,1]:inds[2,2], $
		                        inds[3,0]:inds[3,1]:inds[3,2], $
		                        inds[4,0]:inds[4,1]:inds[4,2], $
		                        inds[5,0]:inds[5,1]:inds[5,2] ]
		7: data = (*self.data)[ inds[0,0]:inds[0,1]:inds[0,2], $
		                        inds[1,0]:inds[1,1]:inds[1,2], $
		                        inds[2,0]:inds[2,1]:inds[2,2], $
		                        inds[3,0]:inds[3,1]:inds[3,2], $
		                        inds[4,0]:inds[4,1]:inds[4,2], $
		                        inds[5,0]:inds[5,1]:inds[5,2], $
		                        inds[6,0]:inds[6,1]:inds[6,2] ]
		8: data = (*self.data)[ inds[0,0]:inds[0,1]:inds[0,2], $
		                        inds[1,0]:inds[1,1]:inds[1,2], $
		                        inds[2,0]:inds[2,1]:inds[2,2], $
		                        inds[3,0]:inds[3,1]:inds[3,2], $
		                        inds[4,0]:inds[4,1]:inds[4,2], $
		                        inds[5,0]:inds[5,1]:inds[5,2], $
		                        inds[6,0]:inds[6,1]:inds[6,2], $
		                        inds[7,0]:inds[7,1]:inds[7,2] ]
		ELSE: Message, 'An IDL array cannot have more than 8 dimensions.'
	ENDCASE
	
	;Return just the data
	IF tf_data THEN RETURN, data
	
;-------------------------------------------
; Create Output Array //////////////////////
;-------------------------------------------
	;Create the variable
	vOut = MrVariable(data, NAME=name, /NO_COPY)
	
	;Copy all attributes
	self -> CopyAttrTo, vOut

;-------------------------------------------
; DEPEND_# /////////////////////////////////
;-------------------------------------------
	;
	; Rules:
	;   - NSUBS must be equal to # dependent variables, otherwise subscripting is indeterminate
	;   - Dimensions of implicit array are ordered [DEPEND_0, DEPEND_1, ..., DEPEND_N]
	;   - DEPEND_[1-N] may have record variance, in which case they are 2D
	;   - If a bulk process currently acting on the implicit, and then it may already
	;     have affected the (cached) DEPEND_[1-N] variable. If DEPEND_[1-N] is cached,
	;     allow it to have the same dimensions as the output variable.
	;
	
	;
	; TODO:
	;   Ideally, this should be done with all poitner variables:
	;   https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html
	;       FORM_PTR
	;       LABL_PTR_[1-3]
	;       SCAL_PTR
	;       UNITS_PTR
	;
	
	;Size of implicit & sub- array
	;   - Original size, not the size of the subarray
	outDims  = Size(vOut, /DIMENSIONS)
	nOutDims = Size(vOut, /N_DIMENSIONS)
	
	;Number of dependent variables
	;   - IDL allows up to 8 dimensions
	nDeps = ~self -> HasAttr('DEPEND_0')   ? 0 $
	        : ~self -> HasAttr('DEPEND_1') ? 1 $
	        : ~self -> HasAttr('DEPEND_2') ? 2 $
	        : ~self -> HasAttr('DEPEND_3') ? 3 $
	        : ~self -> HasAttr('DEPEND_4') ? 4 $
	        : ~self -> HasAttr('DEPEND_5') ? 5 $
	        : ~self -> HasAttr('DEPEND_6') ? 6 $
	        : ~self -> HasAttr('DEPEND_7') ? 7 $
	        : 8
	
	;RETURN if nSubs LE nDeps
	;   - There does not have to be a DEPEND_# for each subscript
	;   - There DOES     have to be a subscript for each DEPEND_#
	IF nOutDims LT nDeps THEN BEGIN
		strNSubs = String(nSubs, FORMAT='(i1)')
		strNDeps = String(nDeps, FORMAT='(i1)')
		IF nDeps GT 0 THEN MrPrintF, 'LogWarn', '# subscripts (' + strNSubs + ') LT # dependent variables (' + strNDeps + ').'
		RETURN, vOut
	ENDIF
	
	;Step through each dependent variable
	FOR i = 0, nDeps-1 DO BEGIN
		depend = 'DEPEND_' + String(i, FORMAT='(i0)')
		IF ~self -> HasAttr(depend) THEN CONTINUE
		
		;Get the dependent variable and its dimensions
		oDep     = self[depend]
		nDepDims = Size(oDep, /N_DIMENSIONS)
		depDims  = Size(oDep, /DIMENSIONS)
		
		;RECORD VARIANCE
		;   - Data     = [N, M, L, P]
		;   - DEPEND_0 = [N]
		;   - DEPEND_1 = [N, M]
		;   - DEPEND_2 = [N, L]
		;   - DEPEND_3 = [N, P]
		IF depend NE 'DEPEND_0' && nDepDims EQ 2 THEN BEGIN
			;Apply to DEPEND_#
			;   - Dimensions of implicit and DEPEND_# arrays match
			IF Array_Equal(depDims, inDims[[0,i]]) THEN BEGIN
				CASE i OF
					1: oDepend = oDep[inds[0,0]:inds[0,1]:inds[0,2], inds[1,0]:inds[1,1]:inds[1,2]]
					2: oDepend = oDep[inds[0,0]:inds[0,1]:inds[0,2], inds[2,0]:inds[2,1]:inds[2,2]]
					3: oDepend = oDep[inds[0,0]:inds[0,1]:inds[0,2], inds[3,0]:inds[3,1]:inds[3,2]]
					4: oDepend = oDep[inds[0,0]:inds[0,1]:inds[0,2], inds[4,0]:inds[4,1]:inds[4,2]]
					5: oDepend = oDep[inds[0,0]:inds[0,1]:inds[0,2], inds[5,0]:inds[5,1]:inds[5,2]]
					6: oDepend = oDep[inds[0,0]:inds[0,1]:inds[0,2], inds[6,0]:inds[6,1]:inds[6,2]]
					7: oDepend = oDep[inds[0,0]:inds[0,1]:inds[0,2], inds[7,0]:inds[7,1]:inds[7,2]]
					ELSE: Message, 'Invalid number of dependent variables.'
				ENDCASE
			
			;Pre-Applied
			;   - It is possible that whatever is happening to the implicit variable
			;     is happening in bulk (e.g. MrVar_TLimit).
			;   - If DEPEND_# is in the cache, allow it to have the same dimensions
			;     as the output variable
			ENDIF ELSE IF MrVar_IsCached(oDep) && Array_Equal(depDims, outDims[[0,i]]) THEN BEGIN
				oDepend = oDep
			
			;Incorrect dimension sizes
			ENDIF ELSE BEGIN
				;Dimensions
				strInDims  = '[' + StrJoin(String(inDims,        FORMAT='(i0)'), ', ') + ']'
				strOutDims = '[' + StrJoin(String(outDims,       FORMAT='(i0)'), ', ') + ']'
				strDims    = '[' + StrJoin(String(inDims[[0,i]], FORMAT='(i0)'), ', ') + ']'
				strDepDims = '[' + StrJoin(String(depDims,       FORMAT='(i0)'), ', ') + ']'
			
				;Warning information
				MrPrintF, 'LogWarn', 'Cannot extract ' + depend + ' data for variable "' + self.name + '".'
				MrPrintF, 'LogWarn', '    ' + depend + ' dims:  ' + strDepDims + '.'
				MrPrintF, 'LogWarn', '    Data dims:      ' + strInDims
				MrPrintF, 'LogWarn', '    Subarray dims:  ' + strOutDims
				MrPrintF, 'LogWarn', '    Expected dims:  ' + strDims
				
				;Invalidate the DEPEND_# variable
				oDepend = 0B
			ENDELSE
			
		;NO RECORD VARIANCE
		;   - Data     = [N, M, L, P]
		;   - DEPEND_0 = [N]
		;   - DEPEND_1 = [M]
		;   - DEPEND_2 = [L]
		;   - DEPEND_3 = [P]
		ENDIF ELSE IF nDepDims EQ 1 THEN BEGIN
			;Apply to DEPEND_#
			;   - Dimensions of implicit and DEPEND_# arrays match
			;   - Prevent scalar 0 from returning object by including "[]"
			IF depDims EQ inDims[i] THEN BEGIN
				CASE i OF
					0: oDepend = oDep[inds[0,0]:inds[0,1]:inds[0,2]]
					1: oDepend = oDep[inds[1,0]:inds[1,1]:inds[1,2]]
					2: oDepend = oDep[inds[2,0]:inds[2,1]:inds[2,2]]
					3: oDepend = oDep[inds[3,0]:inds[3,1]:inds[3,2]]
					4: oDepend = oDep[inds[4,0]:inds[4,1]:inds[4,2]]
					5: oDepend = oDep[inds[5,0]:inds[5,1]:inds[5,2]]
					6: oDepend = oDep[inds[6,0]:inds[6,1]:inds[6,2]]
					7: oDepend = oDep[inds[7,0]:inds[7,1]:inds[7,2]]
					ELSE: Message, 'Invalid number of dependent variables.'
				ENDCASE
			
			;Pre-Applied
			;   - It is possible that whatever is happening to the implicit variable
			;     is happening in bulk (e.g. MrVar_TLimit).
			;   - If DEPEND_# is in the cache, allow it to have the same dimensions
			;     as the output variable
			ENDIF ELSE IF MrVar_IsCached(oDep) && depDims EQ outDims[i] THEN BEGIN
				oDepend = oDep
			
			;Incorrect dimension sizes
			ENDIF ELSE BEGIN
				;Dimensions
				strInDims  = '[' + StrJoin(String(inDims,    FORMAT='(i0)'), ', ') + ']'
				strOutDims = '[' + StrJoin(String(outDims,   FORMAT='(i0)'), ', ') + ']'
				strDims    = '[' + StrJoin(String(inDims[i], FORMAT='(i0)'), ', ') + ']'
				strDepDims = '[' + StrJoin(String(depDims,   FORMAT='(i0)'), ', ') + ']'
			
				;Warning information
				MrPrintF, 'LogWarn', 'Cannot extract ' + depend + ' data for variable "' + self.name + '".'
				MrPrintF, 'LogWarn', '    ' + depend + ' dims:  ' + strDepDims + '.'
				MrPrintF, 'LogWarn', '    Data dims:      ' + strInDims
				MrPrintF, 'LogWarn', '    Subarray dims:  ' + strOutDims
				MrPrintF, 'LogWarn', '    Expected dims:  ' + strDims
				
				;Invalidate the DEPEND_# variable
				oDepend = 0B
			ENDELSE
		
		;Unknown dimension sizes
		ENDIF ELSE BEGIN
			;Dimensions
			strInDims  = '[' + StrJoin(String(inDims,    FORMAT='(i0)'), ', ') + ']'
			strOutDims = '[' + StrJoin(String(outDims,   FORMAT='(i0)'), ', ') + ']'
			strDims    = '[' + StrJoin(String(inDims[i], FORMAT='(i0)'), ', ') + ']'
			strDepDims = '[' + StrJoin(String(depDims,   FORMAT='(i0)'), ', ') + ']'
		
			;Warning information
			MrPrintF, 'LogWarn', 'Cannot extract ' + depend + ' data for variable "' + self.name + '".'
			MrPrintF, 'LogWarn', '    ' + depend + ' dims:  ' + strDepDims + '.'
			MrPrintF, 'LogWarn', '    Data dims:      ' + strInDims
			MrPrintF, 'LogWarn', '    Subarray dims:  ' + strOutDims
			MrPrintF, 'LogWarn', '    Expected dims:  ' + strDims
			
			;Invalidate the DEPEND_# variable
			oDepend = 0B
		ENDELSE
		
		;Create a variable
		IF Obj_Valid(oDepend) GT 0 THEN vOut[depend] = oDepend
	ENDFOR

;-------------------------------------------
; DELTA_(PLUS|MINUS)_VAR ///////////////////
;-------------------------------------------
	deltas  = ['DELTA_MINUS_VAR', 'DELTA_PLUS_VAR']
	nDeltas = N_Elements(deltas)
	
	;Step through each dependent variable
	FOR i = 0, nDeltas-1 DO BEGIN
		IF ~self -> HasAttr(deltas[i]) THEN CONTINUE
		
		;Get the dependent variable and its dimensions
		oDel       = self[deltas[i]]
		nDelDims = Size(oDel, /N_DIMENSIONS)
		delDims  = Size(oDel, /DIMENSIONS)
		
		;SAME
		;   - Both record varying
		;   - Both non-record varying
		IF nDelDims EQ nInDims && Array_Equal(delDims, inDims) THEN BEGIN
			CASE nDelDims OF
				1: oDelta = oDel[ inds[0,0]:inds[0,1]:inds[0,2] ]
				2: oDelta = oDel[ inds[0,0]:inds[0,1]:inds[0,2], $
				                  inds[1,0]:inds[1,1]:inds[1,2] ]
				3: oDelta = oDel[ inds[0,0]:inds[0,1]:inds[0,2], $
				                  inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2] ]
				4: oDelta = oDel[ inds[0,0]:inds[0,1]:inds[0,2], $
				                  inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2] ]
				5: oDelta = oDel[ inds[0,0]:inds[0,1]:inds[0,2], $
				                  inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2], $
				                  inds[4,0]:inds[4,1]:inds[4,2] ]
				6: oDelta = oDel[ inds[0,0]:inds[0,1]:inds[0,2], $
				                  inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2], $
				                  inds[4,0]:inds[4,1]:inds[4,2], $
				                  inds[5,0]:inds[5,1]:inds[5,2] ]
				7: oDelta = oDel[ inds[0,0]:inds[0,1]:inds[0,2], $
				                  inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2], $
				                  inds[4,0]:inds[4,1]:inds[4,2], $
				                  inds[5,0]:inds[5,1]:inds[5,2], $
				                  inds[6,0]:inds[6,1]:inds[6,2] ]
				8: oDelta = oDel[ inds[0,0]:inds[0,1]:inds[0,2], $
				                  inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2], $
				                  inds[4,0]:inds[4,1]:inds[4,2], $
				                  inds[5,0]:inds[5,1]:inds[5,2], $
				                  inds[6,0]:inds[6,1]:inds[6,2], $
				                  inds[7,0]:inds[7,1]:inds[7,2] ]
				ELSE: Message, 'Incorrect number of dimensions (' + String(nDelDims, FORMAT='(i0)') + ').'
			ENDCASE
		
		;DIFFERENT
		;   - Parent is record varying
		;   - DELTA is non-record varying
		ENDIF ELSE IF nDelDims EQ nInDims-1 && Array_Equal(delDims, inDims[1:*]) THEN BEGIN
		
			CASE nDelDims OF
				1: oDelta = oDel[ inds[1,0]:inds[1,1]:inds[1,2] ]
				2: oDelta = oDel[ inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2] ]
				3: oDelta = oDel[ inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2] ]
				4: oDelta = oDel[ inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2], $
				                  inds[4,0]:inds[4,1]:inds[4,2] ]
				5: oDelta = oDel[ inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2], $
				                  inds[4,0]:inds[4,1]:inds[4,2], $
				                  inds[5,0]:inds[5,1]:inds[5,2] ]
				6: oDelta = oDel[ inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2], $
				                  inds[4,0]:inds[4,1]:inds[4,2], $
				                  inds[5,0]:inds[5,1]:inds[5,2], $
				                  inds[6,0]:inds[6,1]:inds[6,2] ]
				7: oDelta = oDel[ inds[1,0]:inds[1,1]:inds[1,2], $
				                  inds[2,0]:inds[2,1]:inds[2,2], $
				                  inds[3,0]:inds[3,1]:inds[3,2], $
				                  inds[4,0]:inds[4,1]:inds[4,2], $
				                  inds[5,0]:inds[5,1]:inds[5,2], $
				                  inds[6,0]:inds[6,1]:inds[6,2], $
				                  inds[7,0]:inds[7,1]:inds[7,2] ]
				ELSE: Message, 'Incorrect number of dimensions (' + String(nDelDims, FORMAT='(i0)') + ').'
			ENDCASE
			
			CASE nDelDims OF
				1: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1],   i2)
				2: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:2], i2, i3)
				3: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:3], i2, i3, i4)
				4: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:4], i2, i3, i4, i5)
				5: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:5], i2, i3, i4, i5, i6)
				6: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:6], i2, i3, i4, i5, i6, i7)
				7: oDelta = oDel -> _OverloadBracketsRightSide(isRange[1:7], i2, i3, i4, i5, i6, i7, i8)
				ELSE: Message, 'Incorrect number of dimensions (' + String(nDelDims, FORMAT='(i0)') + ').'
			ENDCASE
		
		ENDIF ELSE BEGIN
			;Dimensions
			strInDims  = '[' + StrJoin(String(inDims,    FORMAT='(i0)'), ', ') + ']'
			strOutDims = '[' + StrJoin(String(outDims,   FORMAT='(i0)'), ', ') + ']'
			strDims    = '[' + StrJoin(String(inDims[0], FORMAT='(i0)'), ', ') + ']'
			strDelDims = '[' + StrJoin(String(delDims,   FORMAT='(i0)'), ', ') + ']'
		
			;Warning information
			MrPrintF, 'LogWarn', 'Cannot extract ' + deltas[i] + ' data for variable "' + self.name + '".'
			MrPrintF, 'LogWarn', '    ' + deltas[i] + ' dims:  ' + strDelDims + '.'
			MrPrintF, 'LogWarn', '    Data dims:             ' + strInDims
			MrPrintF, 'LogWarn', '    Subarray dims:         ' + strOutDims
			
			;Invalidate the DEPEND_# variable
			oDelta = 0B
		ENDELSE
		
		;Create a variable
		IF Obj_Valid(oDelta) GT 0 THEN vOut[deltas[i]] = oDelta
	ENDFOR
	
;-------------------------------------------
; Finish ///////////////////////////////////
;-------------------------------------------
	
	;Return the array
	RETURN, vOut
END


;+
;   Get class properties.
;
; :Keywords:
;       DIMENSIONS:         out, optional, type=lonarr
;                           Sizes of each dimension of `ARRAY`.
;       MAXIMUM:            out, optional, type=any
;                           Maximum value of `ARRAY`.
;       MINIMUM:            out, optional, type=any
;                           Minimum value of `ARRAY`.
;       N_DIMENSIONS:       out, optional, type=int
;                           Number of dimensions of `ARRAY`.
;       N_ELEMENTS:         out, optional, type=long
;                           Number of elements in `ARRAY`.
;       TNAME:              out, optional, type=string
;                           Type-name of `ARRAY`.
;       TYPE:               out, optional, type=int
;                           Type-code of `ARRAY`.
;       VERBOSE:            out, optional, type=byte
;                           Level of verboseness of printed messages. Values are::
;                               0 - Quiet
;                               1 - Less verbose
;                               2 - More verbose
;-
pro MrVariable::GetProperty, $
DIMENSIONS=dimensions, $
N_DIMENSIONS=n_dimensions, $
N_ELEMENTS=n_elements, $
MINIMUM=minimum, $
MAXIMUM=maximum, $
NAME=name, $
TNAME=tname, $
TYPE=type, $
VERBOSE=verbose
	compile_opt idl2
	on_error, 2

	if arg_present(dimensions)   then dimensions   = size(*self.data, /DIMENSIONS)
	if arg_present(n_dimensions) then n_dimensions = size(*self.data, /N_DIMENSIONS)
	if arg_present(n_elements)   then n_elements   = size(*self.data, /N_ELEMENTS)
	if arg_present(name)         then name         = self.name
	if arg_present(tname)        then tname        = size(*self.data, /TNAME)
	if arg_present(type)         then type         = size(*self.data, /TYPE)
	if arg_present(maximum)      then maximum      = max(*self.data)
	if arg_present(minimum)      then minimum      = min(*self.data)
	IF Arg_Present(verbose)      THEN verbose      = self.verbose
end


;+
;   Determine if the variable has an attribute.
;
; :Params:
;       ATTRNAME:       in, required, type=string,hash,struct
;                       The name of the attribute to be check
;-
function MrVariable::HasAttr, attrName
	compile_opt idl2
	on_error, 2

	;Number of names given
	nNames = n_elements(attrName)
	
	;Scalar name
	if nNames eq 1 then begin
		tf_has = self.attributes -> HasKey(attrName)
	
	;Loop over each name
	;   - Hash::HasKey accepts only scalars
	endif else begin
		tf_has = bytarr(nNames)
		foreach key, attrName, idx do tf_has[idx] = self.attributes -> HasKey(key)
	endelse
	
	return, tf_has
end


;+
;   Print information about the variable.
;-
pro MrVariable::Help
	compile_opt idl2
	
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return
	endif

	;Get general help
	help, self, OUTPUT=help_txt
	parts = stregex(help_txt, '[A-Z]+[ ]+(.*)', /SUBEXP, /EXTRACT)
	help_txt = parts[1]
	
	;Get the attribute names
	attrNames = self -> GetAttrNames(COUNT=count)
	if count eq 0 then begin
		print, help_txt
		print, '  Variable does not contain attributes.'
		return
	endif
	
	;Sort attribute names
	attrNames = attrNames[sort(attrNames)]
	
	;Allocate memory
	attr_txt = strarr(count)
	name_fmt = 'a-' + string(max(strlen(attrNames)), FORMAT='(i0)')

	;Step through each variable
	for i = 0, count-1 do begin
		;Get the value
		value = self -> GetAttrValue(attrNames[i])

		;Turn the value into a string
		;   - Show maximum of 10 values
		;   - Added benefit of creating row vector for strjoin
		;   - Add ellipsis if there are more than 10 values
		tname   = size(value, /TNAME)
		nValues = n_elements(value)
		
		;OBJECT
		;   - SetAttributeValue allows only MrVariable objects.
		if tname eq 'OBJREF' then begin
			if obj_valid(value) then begin
				dims = value.dimensions
				dims = 'Dims=[' + strjoin( strtrim(dims, 2), ',' ) + ']'
				str_value = dims + ' ' + value.name + ' <' + obj_class(value) + '>'
			endif else begin
				str_value = '<InvalidObject>'
			endelse
		
		;UNDEFINED
		endif else if tname eq 'UNDEFINED' then begin
			str_value = IsA(value, /NULL) ? '!Null' : '<undefined>'
		
		;SCALAR
		endif else if nValues eq 1 then begin
			str_value = StrTrim(String(value, /PRINT), 2)
		
		;ARRAY
		endif else begin
			nPts = (nValues-1) < 10
			if size(value, /TNAME) eq 'BYTE' $
				then str_value  = '[' + StrJoin(StrTrim(StrSplit(String(value[0:nPts], /PRINT), ' ', /EXTRACT), 2), ',') $
				else str_value  = '[' + StrJoin(StrTrim(String(value[0:nPts]), 2), ',')
			str_value += (nPts lt nValues-1) ? ',...]' : ']'
		endelse
		
		;Store the value
		attr_txt[i] = string(attrNames[i], str_value, FORMAT='(2x, ' + name_fmt + ', 4x, a0)')
	endfor

	;Print info
	print, help_txt
	print, transpose(attr_txt)
end


;+
;   Compute the inner product of two tensors of rank > 1.
;   This method generalizes the # operator for tensors of rank > 2.
;
; :Params:
;       B:          in, required, type=numeric or MrVariable objref
;                   A tensor with > 1 dimension.
;
; :Returns:
;       RESULT:     out, required, type=MrVariable objref
;                   Result of the inner product.
;-
function MrVariable::InnerProduct, B, $
NAME=name
	compile_opt strictarr
	on_error, 2
	
	;Extract the data
	if isa(B, 'MrVariable') $
		then _B = B['DATA'] $
		else _B = B
	
	;Default name
	if n_elements(name) eq 0 then name = 'InnerProduct(' + self.name + ',B)'

;-------------------------------------------------------
; Dimensions ///////////////////////////////////////////
;-------------------------------------------------------
	
	;Dimensionality of inputs
	dimsA  = size(self, /DIMENSIONS)
	dimsB  = size(B, /DIMENSIONS)
	nDimsA = size(self, /N_DIMENSIONS)
	nDimsB = size(B, /N_DIMENSIONS)
	nDims  = nDimsA < nDimsB
	
	;Output dimensions
	if nDims eq 1 $
		then outDims = nDimsA < nDimsB ? dimsA : dimsB $
		else outDims = nDimsA < nDimsB ? dimsA[0:-2] : dimsB[1:*]

;-------------------------------------------------------
; Reform Inputs ////////////////////////////////////////
;-------------------------------------------------------
	
	;A is currently D1xD2xD3xD4x ... xDN
	;   - Turn A into DixDN matrix
	if nDimsA eq 1 $
		then _A = A['DATA'] $
		else _A = reform( A['DATA'], [ product(dimsA[0:-2]), dimsA[-1] ] )
	
	;B is currently D1xD2xD3xD4x ... xDN
	;   - Turn B into D1xDi matrix
	if nDimsB ge 2 $
		then _B = reform( _B, [ dimsB[0], product(dimsB[1:*]) ] )

;-------------------------------------------------------
; Inner Product ////////////////////////////////////////
;-------------------------------------------------------
	
	temp = temporary(_A) # temporary(_B)

;-------------------------------------------------------
; Finish Up ////////////////////////////////////////////
;-------------------------------------------------------

	;Reform back to the original shape, minus the contracted dimension.
	temp = reform(temp, dimsOut, /OVERWRITE)
	
	;Create time series variable
	result = MrVariable( temp, NAME=name, /NO_COPY )
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
;       X:              in, required, type=numeric array
;                       Current abcissa values or a MrVariable object.
;       Xout:           in, required, type=Numeric array
;                       The new abcissa values.
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
;                       A MrVariable object containing the interpolated data.
;-
function MrVariable::Interpol, X, Xout, $
CACHE=cache, $
NAME=name, $
LSQUADRATIC=lsquadratic, $
NAN=nan, $
QUADRATIC=quadratic, $
SPLINE=spline
	compile_opt idl2
	on_error, 2

;-------------------------------------------------------
; MrVariable Object ////////////////////////////////////
;-------------------------------------------------------
	if IsA(Xout, 'MrVariable') then begin
		;Get the current abscissa values
		xx = self['DEPEND_0'] -> GetData()
		
		;Get the new abscissa values
		dep0_var = MrVar_Get(x['DEPEND_0'])
		Xout     = dep0_var -> GetData()
		
		;Interpolate
		Y = Interpol( *self.data, temporary(xx), temporary(Xout), $
		              LSQUADRATIC = lsquadratic, $
		              NAN         = nan, $
		              QUADRATIC   = quadratic, $
		              SPLINE      = spline )
		
		;Create output variable
		;   - Copy this variable
		;   - Set its data to the interpolated values
		;   - Set its DEPEND_0 attribute to the new values
		varOut = self -> Copy(NAME=name, CACHE=cache)
		varOut -> SetData, Y, /NO_COPY
		varOut -> SetAttrValue, 'DEPEND_0', dep0_var.name

;-------------------------------------------------------
; Normal Arrays ////////////////////////////////////////
;-------------------------------------------------------
	endif else begin
		;Interpolate
		Y = Interpol( *self.data, X, Xout, $
		              LSQUADRATIC = lsquadratic, $
		              NAN         = nan, $
		              QUADRATIC   = quadratic, $
		              SPLINE      = spline)
		
		;Create a new variable with the same attributes
		varOut = self -> Copy(name, CACHE=cache)
		varOut -> SetData, Y, /NO_COPY
		
		;Delete the DEPEND_0 attribute, if there is one
		;   - TODO: Also create a DEPEND_0 variable ?????
		if varOut -> HasAttr('DEPEND_0') then varOut -> RemoveAttr, 'DEPEND_0'
	endelse
	
	;Return the new variable
	return, varOut
end


;+
;   Check if a variable is defined, of a specific type/class, etc.
;
; :Params:
;       TYPE_NAME:      in, optional, type=string
;                       IDL data-type name, object class name, or structure name.
;
; :Keywords:
;       ARRAY:          in, optional, type=boolean, default=0
;                       If set, return true if the array is an array, list, or
;                           structure and (optionally) matches `TYPENAME`.
;       NULL:           in, optional, type=boolean, default=0
;                       If set, return true if the array is equal to !NULL and
;                           (optionally) matches `TYPENAME`.
;       NUMBER:         in, optional, type=boolean, default=0
;                       If set, return true if the array a number and (optionally)
;                           matches `TYPENAME`.
;       SCALAR:         in, optional, type=boolean, default=0
;                       If set, return true if the array is a scalar and (optionally)
;                           matches `TYPENAME`.
;-
function MrVariable::IsA, type_name, $
ARRAY=array, $
NULL=null, $
NUMBER=number, $
SCALAR=scalar
	compile_opt idl2
	on_error, 2

	;Return the results.
	if n_elements(type_name) gt 0 $
		then return, IsA(*self.data, ARRAY=array, NULL=null, NUMBER=number, SCALAR=scalar) $
		else return, IsA(*self.data, type_name, ARRAY=array, NULL=null, NUMBER=number, SCALAR=scalar)
end


;+
;   Determine if the given object is identical to the implicit object (i.e. they
;   reference the same heap variable).
;
; :Params:
;       OVAR:           in, required, type=MrVariable objref
;                       Test if this variable is the same as the implicit variable.
;
; :Returns:
;       TF_SAME:        out, required, type=boolean
;                       Returns true if `OVAR` and the implicit variable reference
;                           the same heap variable. Returns false otherwise.
;
;-
function MrVariable::IsIdentical, oVar
	compile_opt idl2
	on_error, 2

	;Return the results.
	return, obj_valid(self, /GET_HEAP_IDENTIFIER) eq obj_valid(oVar, /GET_HEAP_IDENTIFIER)
end


;+
;   Make an array and store it internally.
;
; :Params:
;       D1:                 in, optional, type=any/integer
;                           Size of the second dimension.
;       D2:                 in, optional, type=integer
;                           Size of the second dimension.
;       D3:                 in, optional, type=integer
;                           Size of the third dimension.
;       D4:                 in, optional, type=integer
;                           Size of the fourth dimension.
;       D5:                 in, optional, type=integer
;                           Size of the fifth dimension.
;       D6:                 in, optional, type=integer
;                           Size of the sixth dimension.
;       D7:                 in, optional, type=integer
;                           Size of the seventh dimension.
;       D8:                 in, optional, type=integer
;                           Size of the eight dimension.
;
; :Keywords:
;       TYPE:               in, optional, type=int/float, default='FLOAT'
;                           The type-name or type-code of the array to be created.
;       _REF_EXTRA:         in, optional, type=any
;                           Any keyword accepted by IDL's Make_Array procedure.
;-
pro MrVariable::Make_Array, D1, D2, D3, D4, D5, D6, D7, D8, $
TYPE=type, $
_REF_EXTRA=extra
	compile_opt idl2
	on_error, 2

	;Get a type-code if a name was given.
	if n_elements(type) gt 0 $
		then tcode = size(type, /TNAME) ne 'STRING' ? type : self -> TypeName2Code(type) 

	;Make the array
	case 1 of
		n_elements(D8) gt 0: *self.data = make_array(D1, D2, D3, D4, D5, D6, D7, D8, $
		                                              TYPE=tcode, _STRICT_EXTRA=extra)
		n_elements(D7) gt 0: *self.data = make_array(D1, D2, D3, D4, D5, D6, D7, $
		                                              TYPE=tcode, _STRICT_EXTRA=extra)
		n_elements(D6) gt 0: *self.data = make_array(D1, D2, D3, D4, D5, D6, $
		                                              TYPE=tcode, _STRICT_EXTRA=extra)
		n_elements(D5) gt 0: *self.data = make_array(D1, D2, D3, D4, D5, $
		                                              TYPE=tcode, _STRICT_EXTRA=extra)
		n_elements(D4) gt 0: *self.data = make_array(D1, D2, D3, D4, $
		                                              TYPE=tcode, _STRICT_EXTRA=extra)
		n_elements(D3) gt 0: *self.data = make_array(D1, D2, D3, $
		                                              TYPE=tcode, _STRICT_EXTRA=extra)
		n_elements(D2) gt 0: *self.data = make_array(D1, D2, $
		                                              TYPE=tcode, _STRICT_EXTRA=extra)
		n_elements(D1) gt 0: *self.data = make_array(D1, $
		                                              TYPE=tcode, _STRICT_EXTRA=extra)
		else: message, 'No dimensions specified for result.' 
	endcase
end


;+
;   Find the maximum value of the implicit array.
;
; :Params:
;       SUBSCRIPT_MAX:      in, optional, type=integer
;                           Index of the maximum value.
;
; :Keywords:
;       ABSOLUTE:           in, optional, type=boolean, default=0
;                           Find the minimum of the absolute value of the implicit array.
;       DIMENSION:          in, optional, type=integer
;                           Find the minimum along a particular dimension (begins with 1).
;       MIN:                out, optional, type=any
;                           Name of a variable to recieve the minimum array element.
;       NAN:                in, optional, type=boolean, default=0
;                           If set, NaN's will be excluded from the result.
;       SUBSCRIPT_MIN:      out, optional, type=integer
;                           Name of a variable to recieve the index of the `MIN`.
;-
function MrVariable::Max, subscript_max, $
ABSOLUTE=absolute, $
DIMENSION=dimension, $
MIN=minimum, $
NAN=nan, $
SUBSCRIPT_MIN=iMin
	compile_opt idl2
	on_error, 2

	;Find only the minimum
	if arg_present(minimum) || arg_present(iMin) then begin
		maximum = max(*self.data, subscript_max, ABSOLUTE=absolute, NAN=nan, DIMENSION=dimension, $
		                                          MIN=minimum, SUSBSCRIPT_MIN=iMin)

	;Find the maximum?
	endif else begin
		maximum = max(*self.data, subscript_max, ABSOLUTE=absolute, NAN=nan, DIMENSION=dimension)
	endelse

	return, maximum
end


;+
;   Find the minimum value of the implicit array.
;
; :Params:
;       SUBSCRIPT_MIN:      in, optional, type=integer
;                           Index of the minimum value.
;
; :Keywords:
;       ABSOLUTE:           in, optional, type=boolean, default=0
;                           Find the minimum of the absolute value of the implicit array.
;       DIMENSION:          in, optional, type=integer
;                           Find the minimum along a particular dimension (begins with 1).
;       MAX:                out, optional, type=any
;                           Name of a variable to recieve the maximum array element.
;       NAN:                in, optional, type=boolean, default=0
;                           If set, NaN's will be excluded from the result.
;       SUBSCRIPT_MAX:      out, optional, type=integer
;                           Index of the maximum value.
;-
function MrVariable::Min, iMin, $
ABSOLUTE=absolute, $
DIMENSION=dimension, $
MAX=maximum, $
NAN=nan, $
SUBSCRIPT_MAX=iMax
	compile_opt idl2
	on_error, 2

	;Find the minimum.
	if arg_present(maximum) || arg_present(iMax) then begin
		minimum = min(*self.data, iMin, ABSOLUTE=absolute, NAN=nan, DIMENSION=dimension, $
		                                MAX=maximum, SUBSCRIPT_MAX=iMax)
	endif else begin
		minimum = min(*self.data, iMin, ABSOLUTE=absolute, NAN=nan, DIMENSION=dimension)
	endelse

	return, minimum
end


;+
;   Create a new copy. Properties of the parent object are passed on to the
;   child object.
;-
function MrVariable::New, array, $
NAME=name, $
NO_COPY=no_copy
	compile_opt idl2
	on_error, 2

	;Return a copy of the array
	;   - Use Obj_Class() so that subclasses can inherit the method.
	return, obj_new(obj_class(self), array, NAME=name, NO_COPY=no_copy)
end 


;+
;   Compute the outer product of two tensors of rank > 1.
;   This method generalizes the ## operator for tensors of rank > 2.
;
; :Params:
;       B:          in, required, type=numeric or MrVariable objref
;                   A tensor with > 1 dimension.
;
; :Returns:
;       RESULT:     out, required, type=MrVariable objref
;                   Result of the inner product.
;-
function MrVariable::OuterProduct, B, $
NAME=name
	compile_opt strictarr
	on_error, 2
	
	;Extract the data
	if isa(B, 'MrVariable') $
		then _B = B['DATA'] $
		else _B = B
	
	;Default name
	if n_elements(name) eq 0 then name = 'OuterProduct(' + self.name + ',B)'

;-------------------------------------------------------
; Dimensions ///////////////////////////////////////////
;-------------------------------------------------------
	
	;Dimensionality of inputs
	dimsA  = size(self, /DIMENSIONS)
	dimsB  = size(B, /DIMENSIONS)
	nDimsA = size(self, /N_DIMENSIONS)
	nDimsB = size(B, /N_DIMENSIONS)
	nDims  = nDimsA < nDimsB
	
	;Output dimensions
	if nDims eq 1 $
		then outDims = nDimsA < nDimsB ? dimsA : dimsB $
		else outDims = nDimsA < nDimsB ? dimsA[1:*] : dimsB[0:-2]

;-------------------------------------------------------
; Reform Inputs ////////////////////////////////////////
;-------------------------------------------------------
	
	;A is currently D1xD2xD3xD4x ... xDN
	;   - Turn A into DixDN matrix
	if nDimsA eq 1 $
		then _A = A['DATA'] $
		else _A = reform( A['DATA'], [ dimsA[0], product(dimsA[1:*]) ] )
	
	;B is currently D1xD2xD3xD4x ... xDN
	;   - Turn B into D1xDi matrix
	if nDimsB ge 2 $
		then _B = reform( _B, [ product(dimsB[0:-2]), dimsB[-1] ] )

;-------------------------------------------------------
; Inner Product ////////////////////////////////////////
;-------------------------------------------------------
	
	temp = temporary(_A) ## temporary(_B)

;-------------------------------------------------------
; Finish Up ////////////////////////////////////////////
;-------------------------------------------------------

	;Reform back to the original shape, minus the contracted dimension.
	temp = reform(temp, dimsOut, /OVERWRITE)
	
	;Create time series variable
	result = MrVariable( temp, NAME=name, /NO_COPY )
	return, result
end


;+
;   Make an array and store it internally.
;
; :Params:
;       SEED:               in, out, required
;                           The seed to use. If undefined, the generic state is used
;                               and returned.
;       D1:                 in, optional, type=any/integer
;                           Size of the first dimension. If an array, then its elements
;                               specified the size of each dimension. Other keywords
;                               are ignored.
;       D2:                 in, optional, type=integer
;                           Size of the second dimension.
;       D3:                 in, optional, type=integer
;                           Size of the third dimension.
;       D4:                 in, optional, type=integer
;                           Size of the fourth dimension.
;       D5:                 in, optional, type=integer
;                           Size of the fifth dimension.
;       D6:                 in, optional, type=integer
;                           Size of the sixth dimension.
;       D7:                 in, optional, type=integer
;                           Size of the seventh dimension.
;       D8:                 in, optional, type=integer
;                           Size of the eight dimension.
;
; :Keywords:
;       _REF_EXTRA:         in, optional, type=any
;                           Any keyword accepted by IDL's RandomU function.
;-
pro MrVariable::RandomN, seed, D1, D2, D3, D4, D5, D6, D7, D8, $
_REF_EXTRA=extra
	compile_opt idl2
	on_error, 2

	;Make the array
	case 1 of
		n_elements(D1) gt 1: *self.data = randomn(seed, D1, _STRICT_EXTRA=extra)
		n_elements(D8) gt 0: *self.data = randomn(seed, D1, D2, D3, D4, D5, D6, D7, D8, _STRICT_EXTRA=extra)
		n_elements(D7) gt 0: *self.data = randomn(seed, D1, D2, D3, D4, D5, D6, D7, _STRICT_EXTRA=extra)
		n_elements(D6) gt 0: *self.data = randomn(seed, D1, D2, D3, D4, D5, D6, _STRICT_EXTRA=extra)
		n_elements(D5) gt 0: *self.data = randomn(seed, D1, D2, D3, D4, D5, _STRICT_EXTRA=extra)
		n_elements(D4) gt 0: *self.data = randomn(seed, D1, D2, D3, D4, _STRICT_EXTRA=extra)
		n_elements(D3) gt 0: *self.data = randomn(seed, D1, D2, D3, _STRICT_EXTRA=extra)
		n_elements(D2) gt 0: *self.data = randomn(seed, D1, D2, _STRICT_EXTRA=extra)
		n_elements(D1) gt 0: *self.data = randomn(seed, D1, _STRICT_EXTRA=extra)
		arg_present(seed):   *self.data = [randomn(seed)]
		else: message, 'No dimensions specified for result.' 
	endcase
end


;+
;   Make an array and store it internally.
;
; :Params:
;       SEED:               in, out, required
;                           The seed to use. If undefined, the generic state is used
;                               and returned.
;       D1:                 in, optional, type=any/integer
;                           Size of the first dimension. If an array, then its elements
;                               specified the size of each dimension. Other keywords
;                               are ignored.
;       D2:                 in, optional, type=integer
;                           Size of the second dimension.
;       D3:                 in, optional, type=integer
;                           Size of the third dimension.
;       D4:                 in, optional, type=integer
;                           Size of the fourth dimension.
;       D5:                 in, optional, type=integer
;                           Size of the fifth dimension.
;       D6:                 in, optional, type=integer
;                           Size of the sixth dimension.
;       D7:                 in, optional, type=integer
;                           Size of the seventh dimension.
;       D8:                 in, optional, type=integer
;                           Size of the eight dimension.
;
; :Keywords:
;       _REF_EXTRA:         in, optional, type=any
;                           Any keyword accepted by IDL's RandomU function.
;-
pro MrVariable::RandomU, seed, D1, D2, D3, D4, D5, D6, D7, D8, $
_REF_EXTRA=extra
	compile_opt idl2
	on_error, 2

	;Make the array
	case 1 of
		n_elements(D1) gt 1: *self.data = randomu(seed, D1, _STRICT_EXTRA=extra)
		n_elements(D8) gt 0: *self.data = randomu(seed, D1, D2, D3, D4, D5, D6, D7, D8, _STRICT_EXTRA=extra)
		n_elements(D7) gt 0: *self.data = randomu(seed, D1, D2, D3, D4, D5, D6, D7, _STRICT_EXTRA=extra)
		n_elements(D6) gt 0: *self.data = randomu(seed, D1, D2, D3, D4, D5, D6, _STRICT_EXTRA=extra)
		n_elements(D5) gt 0: *self.data = randomu(seed, D1, D2, D3, D4, D5, _STRICT_EXTRA=extra)
		n_elements(D4) gt 0: *self.data = randomu(seed, D1, D2, D3, D4, _STRICT_EXTRA=extra)
		n_elements(D3) gt 0: *self.data = randomu(seed, D1, D2, D3, _STRICT_EXTRA=extra)
		n_elements(D2) gt 0: *self.data = randomu(seed, D1, D2, _STRICT_EXTRA=extra)
		n_elements(D1) gt 0: *self.data = randomu(seed, D1, _STRICT_EXTRA=extra)
		arg_present(seed):   *self.data = [randomu(seed)]
		else: message, 'No dimensions specified for result.' 
	endcase
end


;+
;   Resize the array to the given dimension sizes. Dimensions must be an integer
;   multiple or factor of the original dimensions.
;
; :Params:
;       D1:         in, required, type=integer/intarr
;                   A scalar or array of up to eight elements denoting the size of
;                       the new array.
;       D2-8:       in, optional, type=integer
;                   The size of dimensions 2-8 of the new array.
;
; :Returns:
;       OVAR:       out, required, type=objref (MrVariable)
;                   A MrVariable with the reformed data name "*_rebin".
;-
function MrVariable::Rebin, d1, d2, d3, d4, d5, d6, d7, d8, $
SAMPLE=sample
	compile_opt idl2
	on_error, 2
	
	;Rebin the data
	case 1 of
		n_elements(d8) gt 0: result = rebin(*self.data, d1, d2, d3, d4, d5, d6, d7, d8, SAMPLE=sample)
		n_elements(d7) gt 0: result = rebin(*self.data, d1, d2, d3, d4, d5, d6, d7, SAMPLE=sample)
		n_elements(d6) gt 0: result = rebin(*self.data, d1, d2, d3, d4, d5, d6, SAMPLE=sample)
		n_elements(d5) gt 0: result = rebin(*self.data, d1, d2, d3, d4, d5, SAMPLE=sample)
		n_elements(d4) gt 0: result = rebin(*self.data, d1, d2, d3, d4, SAMPLE=sample)
		n_elements(d3) gt 0: result = rebin(*self.data, d1, d2, d3, SAMPLE=sample)
		n_elements(d2) gt 0: result = rebin(*self.data, d1, d2, SAMPLE=sample)
		n_elements(d1) gt 0: result = rebin(*self.data, d1, SAMPLE=sample)
		else:                message, 'Usage: oVar -> Rebin, d1[, d2, ..., d8]'
	endcase
	
	;Create result
	oVar = MrVariable(result, NAME=self.name + '_rebin', /NO_COPY)
	self -> CopyAttrTo, oVar
	
	return, oVar
end


;+
;   Resize the array to the given dimension sizes. Dimensions must be an integer
;   multiple or factor of the original dimensions.
;
; :Params:
;       D1:         in, required, type=integer/intarr
;                   A scalar or array of up to eight elements denoting the size of
;                       the new array.
;       D2-8:       in, optional, type=integer
;                   The size of dimensions 2-8 of the new array.
;
; :Returns:
;       OVAR:       out, required, type=objref (MrVariable)
;                   A MrVariable with the reformed data name "*_reform".
;-
function MrVariable::Reform, d1, d2, d3, d4, d5, d6, d7, d8
	compile_opt idl2
	on_error, 2
	
	;Reform the data
	case 1 of
		n_elements(d8) gt 0: result = reform(*self.data, d1, d2, d3, d4, d5, d6, d7, d8)
		n_elements(d7) gt 0: result = reform(*self.data, d1, d2, d3, d4, d5, d6, d7)
		n_elements(d6) gt 0: result = reform(*self.data, d1, d2, d3, d4, d5, d6)
		n_elements(d5) gt 0: result = reform(*self.data, d1, d2, d3, d4, d5)
		n_elements(d4) gt 0: result = reform(*self.data, d1, d2, d3, d4)
		n_elements(d3) gt 0: result = reform(*self.data, d1, d2, d3)
		n_elements(d2) gt 0: result = reform(*self.data, d1, d2)
		n_elements(d1) gt 0: result = reform(*self.data, d1)
		else:                result = reform(*self.data)
	endcase
	
	;Create result
	oVar = MrVariable(result, NAME=self.name + '_reform', /NO_COPY)
	self -> CopyAttrTo, oVar
	
	return, oVar
end


;+
;   Replace a value within the array. Default replacement values follow the ISTP
;   Guidlines for CDF files.
;
;       BYTE    :  255B
;       INT     : -32768S
;       LONG    : -2147483648L
;       FLOAT   :  !values.f_nan
;       DOUBLE  :  !values.d_nan
;       UINT    :  65535US
;       ULONG   :  4294967295UL
;       LONG64  :  -9223372036854775808LL
;
;   http://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#FILLVAL
;
; :Params:
;       VALUE:          in, optional, type=numeric, default=self['FILLVAL']
;                       The value to be replaced.
;       REPLACEMENT:    in, optional, type=numeric, default=CDF Fill value
;                       The value to replace `VALUE`.
;-
PRO MrVariable::ReplaceValue, value, replacement
	Compile_Opt idl2
	On_Error, 2
	
	;Defaults
	IF N_Elements(value) EQ 0 THEN BEGIN
		IF ~self -> HasAttr('FILLVAL') THEN RETURN
		value = self['FILLVAL']
	ENDIF
	IF N_Elements(replacement) EQ 0 THEN BEGIN
		self -> GetProperty, TNAME=tname
		CASE tname OF
			'BYTE':   replacement = 255B
			'INT':    replacement = Fix(-32768, TYPE=2)
			'LONG':   replacement = Fix(-2147483648, TYPE=3) 
			'FLOAT':  replacement = !values.f_nan
			'DOUBLE': replacement = !values.d_nan
			'UINT':   replacement = 65535US
			'ULONG':  replacement = 4294967295UL
			'LONG64': replacement = -9223372036854775808LL
			ELSE: Message, 'Variable of type "' + TNAME  + '" does not have default fill value.'
		ENDCASE
	ENDIF

	;Find and replace
	iValue = self -> Where(value, /EQ, COUNT=count)
	IF count GT 0 THEN self[iValue] = replacement
END


;+
;   Remove attributes.
;
; :Params:
;       ATTRNAMES:      in, required, type=string/hash/struct
;                       Attribute names to be removed. A warning is issued for names
;                           that do not already exist.
;-
pro MrVariable::RemoveAttr, attrNames
	compile_opt idl2
	on_error, 2
	
	;Which attributes exist
	tf_has = self -> HasAttr(attrNames)
	
	;Separate existent from non-existent
	iBad = where(tf_has eq 0, nBad, COMPLEMENT=iGood, NCOMPLEMENT=nGood)
	
	;Issue warning for non-existent
	if nBad gt 0 then begin
		for i = 0, nBad - 1 do MrPrintF, 'LogWarn', 'No such attribute: "' + attrNames[iBad[i]] + '".', $
		                                 LEVEL = 5
	endif
	
	;Remove the attributes
	if nGood gt 0 then self.attributes -> Remove, attrNames[iGood]
end


;+
;   Set attribute values.
;
;   Calling Sequence:
;       var -> SetAttrValue, struct
;       var -> SetAttrValue, hash
;       var -> SetAttrValue, string, value
;       var -> SetAttrValue, strarr, array
;       var -> SetAttrValue, strarr, list
;
; :Params:
;       ATTRNAME:       in, required, type=string/strarr/hash/struct
;                       The name of the attribute for which the value is to be changed.
;                           Or a hash or structure whos keys/tags are the attribute names.
;                           If an array of strings, then each element is an attribute name
;                           and `ATTRVALUE` may be an array or list of values with the same
;                           number of elements as ATTRNAME. If ATTRNAME is a hash, keys
;                           must be strings. Values cannot be complex datatypes.
;       ATTRVALUE:      in, optional, type=any accepted by ::AttrValue_IsValid
;                       The value of the attribute(s) to be added. ATTRVALUE must be
;
; :Keyword:
;       CREATE:         in, optional, type=boolean, default=0
;                       If set, then attributes will be created if they do not exist.
;       OVERWRITE:      in, optional, type=boolean, default=1
;                       If set to zero, the attribute value will not be set if `ATTRNAME`
;                           already exists as an attribute. Set OVERWRITE=0 if you are
;                           creating new attributes but to not want to accidentally over-
;                           write an existing attribute.
;-
pro MrVariable::SetAttrValue, attrName, attrValue, $
CREATE=create, $
OVERWRITE=overwrite
	compile_opt idl2
	on_error, 2
	
;-------------------------------------------------------
; Hash /////////////////////////////////////////////////
;-------------------------------------------------------
	if isa(attrName, 'HASH') then begin
		
		;Step through each key
		foreach val, attrName, key do begin
			;Skip invalid attributes, but issue warning
			catch, the_error
			if the_error eq 0 then begin
				self -> SetAttributeValue, key, val, $
				                           CREATE    = create, $
				                           OVERWRITE = overwrite
			endif else begin
				catch, /CANCEL
				MrPrintF, 'LogWarn', !error_state.msg
			endelse
		endforeach

;-------------------------------------------------------
; Structure ////////////////////////////////////////////
;-------------------------------------------------------
	endif else if isa(attrName, 'STRUCT') then begin
		
		;Number of attributes and their names
		nTags = n_tags(attrName)
		tags  = tag_names(attrName)
		
		;Loop through each value
		for i = 0, nTags - 1 do begin
			
			;Skip invalid attributes, but issue warning
			catch, the_error
			if the_error eq 0 then begin
				self -> SetAttributeValue, tags[i], attrName.(i), $
				                           CREATE    = create, $
				                           OVERWRITE = overwrite
			endif else begin
				catch, /CANCEL
				MrPrintF, 'LogWarn', !error_state.msg
			endelse
		endfor

;-------------------------------------------------------
; String ///////////////////////////////////////////////
;-------------------------------------------------------
	endif else if isa(attrName, 'STRING', /SCALAR) then begin
		
		;Set the attribute
		self -> SetAttributeValue, attrName, attrValue, $
		                           CREATE    = create, $
		                           OVERWRITE = overwrite

;-------------------------------------------------------
; String Array with Array or List of Values ////////////
;-------------------------------------------------------
	endif else if isa(attrName, 'STRING', /ARRAY) then begin
			
		;Restriction of Hash
		if n_elements(attrName) ne n_elements(attrValue) $
			then message, 'ATTRNAME and ATTRVALUE must have the same number of elements.'
		
		;Loop over each value
		foreach val, attrValue, idx do begin
		
			;Skip invalid attributes, but issue warning
			catch, the_error
			if the_error eq 0 then begin
				self -> SetAttributeValue, attrName[idx], val, $
				                           CREATE    = create, $
				                           OVERWRITE = overwrite
			endif else begin
				catch, /CANCEL
				MrPrintF, 'LogWarn', !error_state.msg
			endelse
		endforeach
		
;-------------------------------------------------------
; Other ////////////////////////////////////////////////
;-------------------------------------------------------
	endif else begin
		;
		; Do not accept lists for ATTRNAME because numeric
		; attribute names are not allowed.
		;
	
		message, 'Invalid value for ATTRNAME.'
	endelse
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
pro MrVariable::SetAttributeValue, attrName, attrValue, $
CREATE=create, $
OVERWRITE=overwrite
	compile_opt idl2
	on_error, 2
	
	tf_create    = keyword_set(create)
	tf_overwrite = n_elements(overwrite) eq 0 ? 1B : keyword_set(overwrite)

;-------------------------------------------------------
; Check Attribute Name /////////////////////////////////
;-------------------------------------------------------
	
	;Check inputs
	if ~IsA(attrName, 'STRING', /SCALAR) $
		then message, 'ATTRNAME must be a scalar string.'
	
	;Forbidden attributes
	;   - These are special inputs to ::_OverloadBracketsRightSide
	IF MrIsMember(['DATA', 'PTR', 'POINTER'], attrName) $
		THEN Message, '"' + attrName + '" cannot be an attribute name.'
	
	;Keywords
	tf_has = self -> HasAttr(attrName)
	if tf_has  && ~tf_overwrite then message, 'Cannot set attribute value. Attr already exists: "' + attrName + '".'
	if ~tf_has && ~tf_create    then message, 'Cannot set attribute value. Attr does not exist: "' + attrName + '".'

;-------------------------------------------------------
; Allow Object References for Pointer Attributes ///////
;-------------------------------------------------------
	;Attributes that point to other variables
	ptrAttrs = [ 'DEPEND_0', 'DEPEND_1', 'DEPEND_2', 'DEPEND_3', $
	             'LABL_PTR_1', 'LABL_PTR_2', 'LABL_PTR_3', $
	             'DELTA_MINUS_VAR', 'DELTA_PLUS_VAR' ]

	;Valid value
	if self -> AttrValue_IsValid(attrValue) then begin
		
		;UNITS
		if attrName eq 'UNITS' then begin
			self -> SetUnits, attrValue
			return
			
		;POINTER
		endif else if MrIsMember(ptrAttrs, attrName) then begin
			;Get the variable it points to
			;   - This ensures that the variable with name AttrValue is in the cache
			;   - Save the variable, not its name (i.e. circumvent the pointer)
			if MrIsA(attrValue, /SCALAR, 'STRING') then begin
				theValue = MrVar_Get(attrValue, COUNT=count)
				if count eq 0 then theValue = attrValue
			endif else begin
				theValue = attrValue
			endelse
		
		;VALUE
		endif else begin
			theValue = attrValue
		endelse
	
	;Invalid value
	endif else begin
		
		;POINTER
		if MrIsMember(ptrAttrs, attrName) then begin
			
			;Allow objects, but only if they are of class MrVariable
			if size(attrValue, /TNAME) eq 'OBJREF' then begin
				if obj_isa(attrValue, 'MrVariable') $
					then theValue = attrValue $
					else message, string( obj_class(attrValue), attrName, $
					              FORMAT='(%"Invalid object class (%s) for attribute %s. Must be MrVariable.")')
			
			;Other datatypes not allowed
			endif else begin
				message, string(size(attrValue, /TNAME), attrName, $
				                FORMAT='(%"Invalid attribute datatype (%s) for attribute %s")')
			endelse
		
		;VALUE
		endif else begin
			message, string(size(attrValue, /TNAME), attrName, $
			                FORMAT='(%"Invalid attribute datatype (%s) for attribute %s")')
		endelse
	endelse

;-------------------------------------------------------
; Set Attribute Value //////////////////////////////////
;-------------------------------------------------------
	
	;Set attribute value
	;   - Will simultaneously create the attribute
	self.attributes[attrName] = theValue
end


;+
;   Set the array.
;
; :Examples:
;   Set the array::
;       myArray  = MrVariable(11, 91, 6, TYPE='FLOAT')
;       myArray -> SetArray, indgen(23)
;       help, myArray
;           MYARRAY     OBJREF      <ObjHeapVar680(MrVariable)>
;             ARRAY       INT       = Array[23]
;
; :Keywords:
;       DATA:               in, required, type=array
;                           Array of values to be stored, or a MrVariable object whose
;                               array is to be copied. 
;
; :Keywords:
;       NO_COPY:            in, optional, type=boolean, default=0
;                           If set `DATA` will be copied directly into the object
;                               and will be left undefined (a MrVariable object will not
;                               be destroyed, but its array will be empty).
;-
pro MrVariable::SetData, data, $
NO_COPY=no_copy
	compile_opt idl2
	on_error, 2
	
	;Defaults
	no_copy = keyword_set(no_copy)

	;Take from a MrVariable object?
	if IsA(data, 'OBJREF') then begin
		if obj_isa(data, 'MrVariable') $
			then *self.data = data -> GetData(NO_COPY=no_copy) $
			else message, 'Only "MrVariable" objects can be given.'
	
	;Regular array.
	endif else begin
		if no_copy $
			then *self.data = temporary(data) $
			else *self.data = data
	endelse
end


;+
;   Set class properties.
;
; :Params:
;       NAME:               in, required, type=string
;                           Name to be given to the variable. If the variable is cached,
;                               then NAME must be unique.
;
; :Keywords:
;       NAME_OUT:           out, optional, type=string
;                           In cases where `NAME` must be unique and `NO_CLOBBER` is
;                               set, then this is the modified, unique name. Otherwise,
;                               NAME_OUT=`NAME`.
;       NO_CLOBBER:         in, optional, type=boolean, default=0
;                           If set and `NAME` matches a variable within the cache, then
;                               append "_#" to `NAME`, where "#" represents the smallest
;                               available unique number.
;-
pro MrVariable::SetName, name, $
NAME_OUT=name_out, $
NO_CLOBBER=clobber
	compile_opt idl2
	on_error, 2

;-----------------------------------------------------
; Cached Variables \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Are the current and new names already in the cache?
	;   - Must use SELF.NAME, not SELF because the ::Init method calls
	;     ::SetName before the object is fully defined.
	if MrVar_IsCached(self.name) && MrVar_IsCached(name) then begin
		;Get all of the names
		MrVar_Names, allNames
		
		;Do not include the current name
		;   - Exclude the case where NAME = Self.NAME
		;   - Allows us to determine if NAME is taken by a different variable
		iKeep = where(allNames ne self.name, nKeep)
		if nKeep gt 0 then allNames = allNames[iKeep]
		
		;Is the new name already taken?
		tf_taken = ~array_equal(stregex(allNames, '^'+name+'$', /BOOLEAN), 0)
		if tf_taken then begin
			;Append '_#' to the name
			if keyword_set(no_clobber) then begin
				parts    = stregex(allNames, '^'+name+'(_[0-9]+)?$', /SUBEXP, /EXTRACT)
				nMax     = max(fix(parts[1,*]))
				name_out = name + '_' + string(nMax+1, FORMAT='(i0)')
			
			;Delete the other variable
			endif else begin
				MrVar_Delete, name
				name_out = name
			endelse
		
		;New name is same as current name
		endif else begin
			name_out = name
		endelse
	
;-----------------------------------------------------
; Not in Cache \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	endif else begin
		name_out = name
	endelse
	
;-----------------------------------------------------
; Set Name \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	self.name = name_out
end


;+
;   Set class properties.
;
; :Keywords:
;       VERBOSE:            in, optional, type=byte
;                           Level of verboseness of printed messages. Options are::
;                               0 - Quiet
;                               1 - Less verbose
;                               2 - More verbose
;-
PRO MrVariable::SetProperty, $
VERBOSE=verbose
	Compile_Opt idl2
	On_Error, 2
	
	;VERBOSE
	IF N_Elements(verbose) GT 0 THEN BEGIN
		IF verbose LT 0 || verbose GT 2 THEN Message, 'VERBOSE must be between 0 and 2.'
		self.verbose = verbose
	ENDIF
END


;+
;   Convert the implicit array to the desired type.
;
; :Keywords:
;       TYPE:               in, optional, type=string/integer, default='INT'
;                           The or type-code type-name to which the implicit array will
;                               be converted.
;-
pro MrVariable::SetType, type
	compile_opt idl2
	on_error, 2

	;Default to an integer
	if n_elements(type) eq 0 then type = 'INT'

	;Type-Name?
	tt = size(type, /TNAME)
	typeCode = tt ne 'STRING' ? type : self -> TypeName2Code(tt)

	;Fix the array type
	*self.data = fix(*self.data, TYPE=typeCode)
end


;+
;   Convert to the specified units. Any conversion coefficient will be factored
;   into the implicit data array. If units are not recognized by the IDLUnit object,
;   the input units will replace the current units without conversion.
;
; :Params:
;       UNITS:              in, required, type=string/objref
;                           The name of the physical units to which the implicit array is
;                               to be converted, or and IDLunit object containing the new
;                               units and scaling coefficient.
;-
PRO MrVariable::SetUnits, units
	Compile_Opt idl2
	On_Error, 2
	
	;Get old units
	IF self.attributes -> HasKey('UNITS') $
		THEN oldUnits = self.attributes['UNITS'] $
		ELSE oldUnits = ''
	
	;Unit Name
	newErr = 0
	IF Size(units, /TNAME) EQ 'STRING' THEN BEGIN
		newUnits = units
	
	;Unit Object
	ENDIF ELSE IF IsA(units, 'IDLUnit') THEN BEGIN
		IF units.quantity NE 1.0 THEN  MrPrintF, 'LogWarn', 'New units cannot have quantity. Ignoring.'
		newUnits = units.unit
	
	;Invalid
	ENDIF ELSE BEGIN
		Message, 'UNITS must be a string or IDLUnit object.'
	ENDELSE
	
	;Try to convert units properly with the IDL unit object
	;   - If units are not recognized, error will occur
	;   - If NEWUNITS has a scale factor an error will occur
	Catch, the_error
	IF oldUnits NE '' && the_error EQ 0 THEN BEGIN
		outUnits = IDLunit(oldUnits + '->' + newUnits)
		IF outUnits.quantity NE 1 THEN *self.data *= outUnits.quantity
		self.attributes['UNITS'] = outUnits.unit
	
	;If conversion fails, change units value
	ENDIF ELSE BEGIN
		Catch, /CANCEL
		self.attributes['UNITS'] = newUnits
	ENDELSE
END


;+
;   Shift elements of a an array along any dimension. Positive (negative) shifts move
;   elements right (left).
;
; :Examples:
;   Create a vector and shift its elements::
;       myArray = MrVariable(findgen(10))
;       myArray -> Shift, 5
;       print, myArray
;           5.000  6.000  7.000  8.000  9.000  0.000  1.000  2.000  3.000  4.000
;
;  Create a 2D array and shift one dimension to the right and one to the left.
;       myArray = MrVariable(indgen(10,5))
;       myArray -> Shift, 5, -3
;
; :Params:
;       S1:         in, required, type=integer/intarr, default=0
;                   Either a scalar indicating the number of elements in the first
;                       dimension to shift or an array containing the shift parameters
;                       for each dimension.
;       S2:         in, required, type=integer/intarr, default=0
;                   Number of elements in the second dimension to shift.
;       S3:         in, required, type=integer, default=0
;                   Number of elements in the third dimension to shift.
;       S4:         in, required, type=integer, default=0
;                   Number of elements in the fourth dimension to shift.
;       S5:         in, required, type=integer, default=0
;                   Number of elements in the fifth dimension to shift.
;       S6:         in, required, type=integer, default=0
;                   Number of elements in the sixth dimension to shift.
;       S7:         in, required, type=integer, default=0
;                   Number of elements in the seventh dimension to shift.
;       S8:         in, required, type=integer, default=0
;                   Number of elements in the eighth dimension to shift.
;-
pro MrVariable::Shift, s1, s2, s3, s4, s5, s6, s7, s8
	compile_opt idl2
	on_error, 2

	;Shift *self.data
	if IsA(s1, /ARRAY) then begin
		*self.data = shift(*self.data, s1)

	;Scalars
	endif else begin
		case 1 of
			n_elements(s8) gt 0: *self.data = shift(*self.data, s1, s2, s3, s4, s5, s6, s7, s8)
			n_elements(s7) gt 0: *self.data = shift(*self.data, s1, s2, s3, s4, s5, s6, s7)
			n_elements(s6) gt 0: *self.data = shift(*self.data, s1, s2, s3, s4, s5, s6)
			n_elements(s5) gt 0: *self.data = shift(*self.data, s1, s2, s3, s4, s5)
			n_elements(s4) gt 0: *self.data = shift(*self.data, s1, s2, s3, s4)
			n_elements(s3) gt 0: *self.data = shift(*self.data, s1, s2, s3)
			n_elements(s2) gt 0: *self.data = shift(*self.data, s1, s2)
			n_elements(s1) gt 0: *self.data = shift(*self.data, s1)
			else: message, 'Incorrect number of parameters.'
		endcase
	endelse
end


;+
;   Sort the implicit array.
;
; :Keywords:
;       INDEX:      out, optional, type=intarr
;                   The indices into the original array, in the order needed to sort.
;-
pro MrVariable::Sort, $
INDEX=index
	compile_opt idl2
	on_error, 2
	
	;Return the indices
	if arg_present(index) then begin
		index      = sort(*self.data)
		*self.data = (*self.data)[index]
	
	;Sort in one step
	endif else begin
		*self.data = (*self.data)[sort(*self.data)]
	endelse
end


;+
;   Smooth the implicit array.
;
; :Keywords:
;       INDEX:      out, optional, type=intarr
;                   The indices into the original array, in the order needed to sort.
;-
function MrVariable::Smooth, width, $
CACHE=cache, $
NAME=name, $
_REF_EXTRA=extra
	compile_opt idl2
	on_error, 2
	
	if n_elements(name) eq 0 then name = 'Smooth(' + self.name + ')'
	
	;Smooth
	temp = smooth(self['DATA'], width, _STRICT_EXTRA=extra)
	
	;Create a new variable
	oResult = self -> Copy(name, CACHE=cache)
	oResult -> SetData, temp, /NO_COPY
	oResult['CATDESC'] = 'Smoothed version of "' + self.name + '".'
	
	return, oResult
end


;+
;   Convert a list, hash, or structure to an array. All elements/keys/tags must be
;   of the same type. See MrConcatenate.pro for details.
;
; :Examples:
;   Convert a structure to an array::
;       struct = {A: indgen(13, 4,  8), $
;                 B: indgen(13, 4, 12)}
;       myArray = MrVariable()
;       myArray -> ToArray, struct, DIMENSION=3
;       help, myArray
;           MYARRAY     OBJREF      <ObjHeapVar682(MrVariable)>
;             ARRAY       INT       = Array[13, 4, 20]
;
; :Params:
;       DATA:              in, required, type=any
;                           Structure, List, or Hash whose contents are to be converted
;                               to a single array. All elements must be able to be
;                               concatenated along `DIMENSION`.
;
; :Keywords:
;       DIMENSION:          in, optional, type=integer, default=1
;                           Dimension, starting with 1, to be concatenated.
;       NO_COPY:            in, optional, type=boolean, default=0
;                           If set `DATA` will be copied directly into the object
;                               and will be left undefined.
;-
pro MrVariable::ToArray, data, $
DIMENSION=dimension, $
NO_COPY=no_copy
	compile_opt idl2
	on_error, 2

	no_copy = keyword_set(no_copy)
	if n_elements(dimension) eq 0 then dimension = 0

	type = size(data, /TNAME)
	if type ne 'OBJREF' and type ne 'STRUCT' then $
		message, 'Only array, lists, hashes, and structures may be given.'
	
;---------------------------------------------------------------------
; Structure //////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	if type eq 'STRUCT' then begin
		nTags = n_tags(data)
		for i = 1, nTags-1 do begin
			if i eq 1 $
				then cat_array = MrConcatenate(data.(0), data.(i), dimension) $
				else cat_array = MrConcatenate(temporary(cat_array), data.(i), dimension)
		endfor
	
		*self.data = temporary(cat_array)
	  
;---------------------------------------------------------------------
; List ///////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	endif else if obj_isa(data, 'LIST') then begin
		*self.data = MrConcatenate(data, dimension)
	
;---------------------------------------------------------------------
; Hash ///////////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	endif else if obj_isa(data, 'HASH') then begin
		keys  = data -> Keys()
		nKeys = data -> Count()

		;Step through the keys
		for i = 1, nTags-1 do begin
			if i eq 1 $
				then cat_array = MrConcatenate(data[keys[0]], data[keys[i]], dimension) $
				else cat_array = MrConcatenate(temporary(cat_array), data[keys[i]], dimension)
		endfor
	
		*self.data = temporary(cat_array)

	endif else begin
		message, 'Only Structures, Lists, and Hash objects can be converted to arrays.'
	endelse

	if no_copy then data = !Null
end


;+
;   Transpose the array.
;
; :Keywords:
;       P:          in, optional, type=intarr
;                   Vector specifying how the dimensions are to be permuted. Dimensions
;                       start at zero and cannot be repeated.
;-
function MrVariable::Transpose, P
	compile_opt idl2
	on_error, 2

	;Get the array and transpose it.
	case n_params() of
		0: result = transpose(*self.data)
		1: result = transpose(*self.data, P)
		else: message, 'Usage: myArray -> Transpose([P])'
	endcase
	
	;Create a variable object
	result = MrVariable(result, /NO_COPY, NAME='Transpose(' + self.name + ')')
	self  -> CopyAttrTo, result
	
	;Return
	return, result
end


;+
;   Convert the implicit array to the desired type.
;
; :Keywords:
;       TYPENAME:           in, required, type=string
;                           The type-name to be converted to a type-code
;
; :Returns:
;       TYPECODE:           The type-code belonging to `TYPENAME`
;-
function MrVariable::TypeName2Code, type
	compile_opt idl2
	on_error, 2

	;Type-Name?
	case strupcase(type) of
		'UNDEFINED': typeCode = 0
		'BYTE':      typeCode = 1
		'INT':       typeCode = 2
		'LONG':      typeCode = 3
		'FLOAT':     typeCode = 4
		'DOUBLE':    typeCode = 5
		'COMPLEX':   typeCode = 6
		'STRING':    typeCode = 7
		'STRUCT':    typeCode = 8
		'DCOMPLEX':  typeCode = 9
		'POINTER':   typeCode = 10
		'OBJREF':    typeCode = 11
		'UINT':      typeCode = 12
		'ULONG':     typeCode = 13
		'LONG64':    typeCode = 14
		'ULONG64':   typeCode = 15
		else: message, 'Type name does not exist: "' + type + '".'
	endcase

	return, typeCode
end


;+
;   Indices of the unique elements within the implicit array.
;
; :Keywords:
;       COMPLEMENT:         out, optional, type=intarr
;                           The index values of the non-unique elements of the implicit array.
;       COUNT:              out, optional, type=integer
;                           The number of unique elements in the implicit array.
;       NAN                 in, optional, type=boolean, default=0
;                           If set, all NaN values will be excluded from `IUNIQ`. Instead,
;                               they will be included in `COMPLEMENT` and `NCOMPLEMENT`.
;       NCOMPLEMENT:        out, optional, type=intarr
;                           The number of non-unique values in the implicit array.
;       VALUES:             out, optional, type=any
;                           Unique values of the imipicit array.
;
; :Returns:
;       IUNIQ:              Indices of the unique elements of the implicit array.
;-
function MrVariable::Uniq, $
COMPLEMENT=complement, $
COUNT=count, $
NAN=nan, $
NCOMPLEMENT=nComplement, $
VALUES=values
	compile_opt idl2
	on_error, 2

	;Find uniq indices
	iSort = sort(*self.data)
	iUniq = uniq(*self.data, iSort)

	;Remove NaNs
	if keyword_set(nan) then begin
		iFinite = where(finite((*self.data)[iUniq]) eq 1, nFinite)
		if nFinite gt 0 then begin
			iUniq = iUniq[iFinite]
		endif else begin
			iUniq = !Null
			count = 0
		endelse
	endif

	;Count the number of uniq values
	if arg_present(count) then count = n_elements(unique_inds)

	;Get the index values of the non-unique elements
	if arg_present(complement) or arg_present(ncomplement) then begin
		complement = where(histogram(iUniq, MIN=0, MAX=n_elements(*self.data), BINSIZE=1) eq 0, nComplement, /NULL)
		if nComplement gt 0 then complement = iSort[complement]
	endif

	;Return the uniq values?
	if arg_present(values) then values = (*self.data)[iSort[iUniq]]

	;Return the uniq indices in the correct order.
	return, iSort[iUniq]
end


;+
;   Finds the intervals within a given monotonic vector that brackets a given set of
;   one or more search values. See IDL's `Value_Locate <http://exelisvis.com/docs/VALUE_LOCATE.html>`
;
; :Params:
;       VALUE:              in, required, type=any
;                           Values to be located in the implicit array.
;
; :Keywords:
;       L64:                in, optional, type=boolean, default=0
;                           If set, indices will be returned as type Long64.
;       SORT:               in, optional, type=boolean, default=0
;                           If set, the implicit array will first be sorted into
;                               ascending order, as required by the Value_Locate function.
;
; :Returns:
;       RESULT:             Indices into the implicit array.
;-
function MrVariable::Value_Locate, value, $
L64=l64, $
SORT=sort
	compile_opt idl2
	on_error, 2

	;Sort first?
	if keyword_set(sort) then begin
		iSort = sort(*self.data)
		result = value_locate((*self.data)[iSort], value, L64=l64)
		result = result[iSort]
	endif else begin
		result = value_locate(*self.data, value, L64=l64)
	endelse

	return, result
end


;+
;   Compare the implicit array with a particular value and returns the indices that
;   the make the comparison true.
;
; :Params:
;       VALUE:              in, required, type=any
;                           Value to be the subject of the Where search.
;
; :Keywords:
;       COMPLEMENT:         out, optional, type=intarr
;                           Indices where the comparison failed.
;       COUNT:              out, optional, type=integer
;                           Number of true comparisons.
;       L64:                in, optional, type=boolean, default=0
;                           If set, indices will be returned as type Long64.
;       NCOMPLEMENT:        out, optional, type=integer
;                           Number of false comparisons.
;       EQUAL:              in, optional, type=boolean
;                           If set, use the EQ operator for the comparison. If no other
;                               input keyords are given, this is assumed.
;       GREATER:            in, optional, type=boolean, default=0
;                           If set, use the GT operator for the comparison.
;       GEQ:                in, optional, type=boolean, default=0
;                           If set, use the GE operator for the comparison.
;       NOTEQ:              in, optional, type=boolean, default=0
;                           If set, use the NE operator for the comparison.
;       LESS:               in, optional, type=boolean, default=0
;                           If set, use the LT operator for the comparison.
;       LEQ:                in, optional, type=boolean, default=0
;                           If set, use the LE operator for the comparison.
;       MATCHES:            in, optional, type=boolean
;                           If set, then the elements of the implicit array will be
;                               returned instead of their index locations. Cannot
;                               be used with `MULTID`.
;       MULTID:             in, optional, type=boolean, default=0
;                           If set, `RESULTS` will contain the multi-dimensional array
;                               indices of each match. Normally, 1D array indices are
;                               returned. Cannot be used with `MATCHES`.
;
; :Returns:
;       RESULT:             Indices into the implicit array where the comparison returned
;                               true. If `COUNT`=0, !Null is returned. If `MATCHES` is
;                               set, then the matching elements are returned.
;-
function MrVariable::Where, value, $
COMPLEMENT=complement, $
COUNT=count, $
L64=l64, $
NCOMPLEMENT=ncomplement, $
MATCHES=matches, $
MULTID=multiD, $
;Relational Operators
EQUAL=equal, $
GREATER=greater, $
GEQ=GEQ, $
NOTEQ=notEQ, $
LESS=less, $
LEQ=LEQ
	compile_opt idl2
	on_error, 2

	;Defaults
	equal   = keyword_set(equal)
	greater = keyword_set(greater)
	geq     = keyword_set(geq)
	noteq   = keyword_set(notEQ)
	less    = keyword_set(less)
	leq     = keyword_set(leq)

	;Resolve conflicts
	nKeys = equal + less + leq + greater + geq
	if nKeys gt 1 then message, 'Conflicting keywords. Only one comparison keyword is allowed.'
	if nKeys eq 0 then equal = 1

	;Check where
	case 1 of
		equal:   result = where(*self.data eq value, count, COMPLEMENT=complement, NCOMPLEMENT=ncomplement, L64=l64, /NULL)
		greater: result = where(*self.data gt value, count, COMPLEMENT=complement, NCOMPLEMENT=ncomplement, L64=l64, /NULL)
		geq:     result = where(*self.data ge value, count, COMPLEMENT=complement, NCOMPLEMENT=ncomplement, L64=l64, /NULL)
		noteq:   result = where(*self.data ne value, count, COMPLEMENT=complement, NCOMPLEMENT=ncomplement, L64=l64, /NULL)
		less:    result = where(*self.data lt value, count, COMPLEMENT=complement, NCOMPLEMENT=ncomplement, L64=l64, /NULL)
		leq:     result = where(*self.data le value, count, COMPLEMENT=complement, NCOMPLEMENT=ncomplement, L64=l64, /NULL)
	endcase

	;Return the matches as well.
	if arg_present(matches) then begin
		result = (*self.data)[result]
	
	;Return the multi-dimensional array indices
	endif else if keyword_set(multiD) then begin
		dims   = size(*self.data, /DIMENSIONS)
		result = array_indices(dims, result, /DIMENSIONS)
	endif

	return, result
end


;+
;   The class definition statement.
;
; :Params:
;       CLASS:          out, optional, type=structure
;
; :Fields:
;       DATA:       Data to be accessed via bracket overloading.
;       NAME:       Name of the variable object.
;       VERBOSE:    Level of verboseness for warnings, information, debugging.
;       ATTRIBUTES: Hash of variable attributes and their values.
;-
pro MrVariable__DEFINE
	compile_opt idl2
	
	class = { MrVariable, $
	          inherits IDL_Object, $
	          data:             Ptr_New(), $
	          name:             '', $
	          verbose:          0B, $
	          
	          attributes:       obj_new(), $
	          
	          ;Append
	          _append_before:   0B, $
	          _append_cont:     0B, $
	          _append_putFirst: 0B, $
	          _append_dim:      0B, $
	          _append_iLast:    0L, $
	          _append_no_copy:  0B, $
	          _append_unReform: ptr_new(), $
	          _append_unTrans:  ptr_new() $
	        }
end