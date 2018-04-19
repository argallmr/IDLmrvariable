; docformat = 'rst'
;
; NAME:
;   MrVar_RecipVec__Define
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
;   Compute spatial derivatives via the reciprocal vector method.
;
;   References:
;       Daly, P. W., and G. Paschmann (1998), Analysis Methods for Multi-spacecraft Data,
;           ISSI scientific report SR, book, International Space Science Institute.
;
;       Paschmann, G., and P. W. Daly (Eds.) (2008), Multi-Spacecraft Analysis
;           Methods Revisited.
;
; :Categories:
;   MrVariable
;
; :See Also:
;   MrVariable__Define.pro
;   MrTimeSeries__Define.pro
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
;       2016/11/23  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   The initialization method.
;
;   Calling Sequence:
;       oRV = MrVar_RecipVec()
;       oRV = MrVar_RecipVec(r1, r2, r3, r4)
;
; :Params:
;       R1:         in, optional, type=integer/string/objref
;                   Name, number, or objref of the variable field from spacecraft 1.
;       R2:         in, optional, type=integer/string/objref
;                   Name, number, or objref of the variable field from spacecraft 2.
;       R3:         in, optional, type=integer/string/objref
;                   Name, number, or objref of the variable field from spacecraft 3.
;       R4:         in, optional, type=integer/string/objref
;                   Name, number, or objref of the variable field from spacecraft 4.
;
; :Keywords:
;       TIME:           in, optional, type=MrTimeVar
;                       If provided, `R1`, `R2`, `R3`, `R4` are interpolated these time
;                           stamps. The default is to interpolate to the time stamps of
;                           `R1`.
;-
function MrVar_RecipVec::INIT, r1, r2, r3, r4, $
TIME=time
	compile_opt idl2

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return, 0
	endif
	
;-----------------------------------------------------
; Load Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	if n_params() eq 4 then self -> SetData, r1, r2, r3, r4, TIME=time
	
	return, 1
end


;+
;   Clean up after the object is destroyed
;-
pro MrVar_RecipVec::CLEANUP
	compile_opt idl2
	on_error, 2

	obj_destroy, self.oR1
	obj_destroy, self.oR2
	obj_destroy, self.oR3
	obj_destroy, self.oR4
end


;+
;   Check if the time stamps of the input parameters are identical. Time stamps
;   are considered identical if the time variables have the same heap identifier.
;
; :Params:
;       OP1:        in, required, type=objref
;                   Objref of the variable field from spacecraft 1.
;       OP2:        in, required, type=objref
;                   Objref of the variable field from spacecraft 2.
;       OP3:        in, required, type=objref
;                   Objref of the variable field from spacecraft 3.
;       OP4:        in, required, type=objref
;                   Objref of the variable field from spacecraft 4.
;
; :Returns:
;       SUCCESS:    out, required, type=boolean
;                   Returns true if time stamps are identical and false otherwise.
;-
function MrVar_RecipVec::CheckTimes, oP1, oP2, oP3, oP4
	compile_opt idl2
	on_error, 2

	;Check if times are identical
	oTime = self.oR1['TIMEVAR']
	success = oP1['TIMEVAR'] -> IsIdentical(oTime) && $
	          oP2['TIMEVAR'] -> IsIdentical(oTime) && $
	          oP3['TIMEVAR'] -> IsIdentical(oTime) && $
	          oP4['TIMEVAR'] -> IsIdentical(oTime)
	
	return, success
end


;+
;   Compute the curl using the reciprocal vector method.
;
; :Params:
;       P1:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 1.
;       P2:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 2.
;       P3:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 3.
;       P4:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 4.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the result will be added to the cache.
;       NAME:           in, optional, type=string, default='RecipCurl'
;                       A name to be given to the variable.
;-
function MrVar_RecipVec::Curl, p1, p2, p3, p4, $
NAME=name, $
CACHE=cache
	compile_opt idl2
	on_error, 2

	;Check inputs
	if n_elements(name) eq 0 then name = 'RecipCurl'
	
	;Get the variables
	oP1 = MrVar_Get(p1)
	oP2 = MrVar_Get(p2)
	oP3 = MrVar_Get(p3)
	oP4 = MrVar_Get(p4)
	if ~self -> CheckTimes(oP1, oP2, oP3, oP4) then $
		message, 'Time tags of P1-4 must be identical to those of the position vectors.'

	;Create an array to cycle through quantities.
	oPP = [oP1, oP2, oP3, oP4]
	
	;Get the reciprocal vectors
	;   - Array of MrVectorTS objects
	oRV = self -> ReciprocalVectors()
	
	;Compute the curl
	oCurl = oRV[0] -> Cross( oPP[0] )
	for i = 1, 3 do oCurl += oRV[i] -> Cross( oPP[i] )
	
	;Name and cache
	oCurl -> SetName, name
	if keyword_set(cache) then oCurl -> Cache

	return, oCurl
end


;+
;   Compute the curl using the reciprocal vector method.
;
; :Params:
;       P1:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 1.
;       P2:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 2.
;       P3:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 3.
;       P4:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 4.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the result will be added to the cache.
;       NAME:           in, optional, type=string, default='RecipDiv'
;                       A name to be given to the variable.
;-
function MrVar_RecipVec::Divergence, p1, p2, p3, p4, $
CACHE=cache, $
NAME=name
	compile_opt idl2
	on_error, 2
	
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return, !Null
	endif
	
	;Defaults
	if n_elements(name) eq 0 then name = 'RecipDiv'
	
	;Get the variables
	oP1 = MrVar_Get(p1)
	oP2 = MrVar_Get(p2)
	oP3 = MrVar_Get(p3)
	oP4 = MrVar_Get(p4)
	if ~self -> CheckTimes(oP1, oP2, oP3, oP4) then $
		message, 'Time tags of P1-4 must be identical to those of the position vectors.'

	;Create an array to cycle through quantities.
	oPP = [oP1, oP2, oP3, oP4]
	
	;Get the reciprocal vectors
	;   - Array of MrVectorTS objects
	oRV = self -> ReciprocalVectors()
	
	;
	; Divergence of a vector
	;   - div(V) = ∂_i V_i
	;            = dVx/dx + dVy/dy + dVz/dz
	;
	; Divergence of a tensor
	;   - div(T) = ∂_i T_ij
	;            = (dTxx/dx + dTyx/dy + dTzx/dz)_x
	;              (dTyx/dx + dTyy/dy + dTzy/dz)_y
	;              (dTzx/dx + dTzy/dy + dTzz/dz)_z
	;
	; Sum the contribution at each vertex (Einstein summation notation).
	;
	; Indices in comments refer to
	;   - i : component i
	;   - j : component j
	;   - v : tetrahedron vertex v
	;   - t : time
	; and where
	;   - S : scalar
	;   - V : vector
	;   - T : tensor
	; 

	;Vector
	if obj_isa(oP1, 'MrVectorTS') then begin
		;∂_tiv • V_ti
		oDiv = oRV[0] -> Dot( oPP[0] )
		for i = 1, 3 do oDiv += oRV[i] -> Dot( oPP[i] )
	
	;Matrix
	endif else if obj_isa(oP1, 'MrMatrixTS') then begin
		oDiv = oRV[0] -> Dot( oPP[i] )
		for i = 0, 3 do oDiv += oRV[i] -> Dot( oPP[i] )
	
	;Other
	endif else begin
		message, 'Cannot take divergence of object class "' + obj_class(oP1) + '".'
	endelse

	;Name and cache
	oDiv -> SetName, name
	if keyword_set(cache) then oDiv -> Cache

	return, oDiv
end


;+
;   Get the position vectors of each of the four spacecraft
;
; :Private:
;
; :Params:
;       R1:         out, optional, type=MrVectorTS objref
;                   The position vector for spacecraft 1.
;       R2:         out, optional, type=MrVectorTS objref
;                   The position vector for spacecraft 2.
;       R3:         out, optional, type=MrVectorTS objref
;                   The position vector for spacecraft 3.
;       R4:         out, optional, type=MrVectorTS objref
;                   The position vector for spacecraft 4.
;-
pro MrVar_RecipVec::GetData, r1, r2, r3, r4
	compile_opt idl2
	on_error, 2
	
	if arg_present(r1) then r1 = self.oR1
	if arg_present(r2) then r2 = self.oR2
	if arg_present(r3) then r3 = self.oR3
	if arg_present(r4) then r4 = self.oR4
end


;+
;   Compute the gradient using the reciprocal vector method.
;
; :Params:
;       P1:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 1.
;       P2:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 2.
;       P3:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 3.
;       P4:         in, required, type=objref
;                   Name, number, or objref of the variable field from spacecraft 4.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the result will be added to the cache.
;       NAME:           in, optional, type=string, default='RecipGrad'
;                       A name to be given to the variable.
;-
function MrVar_RecipVec::Gradient, p1, p2, p3, p4, $
CACHE=cache, $
NAME=name
	compile_opt idl2
	on_error, 2
	
	;Defaults
	if n_elements(name) eq 0 then name = 'RecipGrad'
	
	;Get the variables
	oP1 = MrVar_Get(p1)
	oP2 = MrVar_Get(p2)
	oP3 = MrVar_Get(p3)
	oP4 = MrVar_Get(p4)
	if ~self -> CheckTimes(oP1, oP2, oP3, oP4) then $
		message, 'Time tags of P1-4 must be identical to those of the position vectors.'

	;Create arrays to cycle through quantities.
	nPts = oP1 -> GetNPts()
	oPP  = [oP1, oP2, oP3, oP4]
	
	;Get the reciprocal vectors
	;   - Array of MrVectorTS objects
	oRV = self -> ReciprocalVectors()
	
	;
	; Gradient of a scalar
	;   - grad(s) = ∂_i S
	;             = (dS/dx, dS/dy, dS/dz)
	;
	; Gradient of a vector
	;   - grad(v) = ∂_i V_j
	;             = (dvx/dx + dvx/dy + dvx/dz)_x
	;               (dvy/dx + dvy/dy + dvy/dz)_y
	;               (dvz/dx + dvz/dy + dvz/dz)_z
	;
	; Sum the contribution at each vertex.
	;
	; Indices in comments refer to
	;   - i : component i
	;   - j : component j
	;   - v : tetrahedron vertex v
	;   - t : time
	; and where
	;   - S : scalar
	;   - V : vector
	;   - T : tensor
	; 

	;SCALAR Field
	if obj_isa(oP1, 'MrScalarTS') then begin
		oGrad = oRV[0] * oPP[0]
		for i = 1, 3 do oGrad += oRV[i] * oPP[i]
	
	;VECTOR Field
	endif else if obj_isa(oP1, 'MrVectorTS') then begin
		oGrad = oRV[0] ## oPP[0]
		for i = 1, 3 do oGrad += oRV[i] ## oPP[i]
	
	;UNKNOWN field
	endif else begin
		message, 'Cannot take gradient of object class "' + obj_class(oP1) + '".'
	endelse

	;Name and cache
	oGrad -> SetName, name
	if keyword_set(cache) then oGrad -> Cache

	return, oGrad
end


;+
;   Load data from source files.
;-
pro MrVar_RecipVec::Load
	message, 'MrVar_RecipVec::Load must be over-ridden by a subclass.'
end


;+
;   Compute the reciprocal vectors.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the reciprocal vectors will be added to the cache.
;-
function MrVar_RecipVec::ReciprocalVectors, $
CACHE=cache
	compile_opt idl2
	on_error, 2

	;Number of position vectors
	tf_cache = keyword_set(cache)
	nV       = self.oR1 -> GetNPts()

	;Create an object array that we can cycle through
	recvec = objarr(4)
	oRR    = [self.oR1, self.oR2, self.oR3, self.oR4]

	;Initialize some indices
	j = 1
	k = 2
	m = 3
	for i = 0, 3 do begin
		;Calculate the separation between spacecraft
		;r_ij = r_j - r_i   --  the position vector pointing from spacecraft i to
		;                       spacecraft j
		r_ij = (oRR[j] - oRR[i])
		r_ik = (oRR[k] - oRR[i])
		r_im = (oRR[m] - oRR[i])

		;The reciprocal vector for vertex m of the tetrahedron points normal to the area
		;of the face of the tetrahedron opposite to vertex m. The normal vector to this
		;surface can be found by taking the cross product between to vectors that lie in
		;the plane of the surface.
		area = r_ik -> Cross(r_im)
	
		;Calculate the volume of the tetrahedron
		volume = r_ij -> Dot(area)
	
		;The reciprical vector is the area's normal vector normalized to the tetrahedron
		;volume
		recvec[j] = area / volume
		
		;Name and cache
		recvec[j] -> SetName, 'RecipVec_Face' + strjoin(string([j,k,m], FORMAT='(i1)'))
		if tf_cache then recvec[j] -> Cache

		;increase the indices cyclically
		j = (j + 1) mod 4
		k = (k + 1) mod 4
		m = (m + 1) mod 4
	endfor
	
	return, recvec
end


;+
;   Set the position vectors of each of the four spacecraft. All variables will be
;   interpolated to P1 (see ::SetData_P).
;
; :Private:
;
; :Params:
;       R1:         in, required, type=objref
;                   Name, number, or MrVectorTS objref of the position vector for spacecraft 1.
;       R2:         in, required, type=objref
;                   Name, number, or MrVectorTS objref of the position vector for spacecraft 2.
;       R3:         in, required, type=objref
;                   Name, number, or MrVectorTS objref of the position vector for spacecraft 3.
;       R4:         in, required, type=objref
;                   Name, number, or MrVectorTS objref of the position vector for spacecraft 4.
;
; :Keywords:
;       TIME:       in, required, type=integer/string/objref
;                   Name, number, or objref of a time or time-series variable. R1-4 will
;                       be interpolated to these time tags. The default is to interpolate
;                       to the time tags of `R1`.
;-
pro MrVar_RecipVec::SetData, r1, r2, r3, r4, $
TIME=time
	compile_opt idl2
	on_error, 2
	
	;Check object class
	R        = MrVar_Get(r1)
	objClass = obj_class(R)
	if ~obj_isa(R, 'MrVectorTS') $
		then message, 'P1-4 must be a MrVectorTS object.'
	
	;Get the time variable
	;   - Must call ::SetData_P first
	;   - All variables will be interpolated to P1
	if n_elements(time) gt 0 then begin
		if obj_isa(time, 'MrTimeVar') $
			then oTime = time $
			else oTime = MrVar_Get(time['DEPEND_0'])
	endif else begin
		oTime = R['TIMEVAR']
	endelse

	;Loop over all variables
	oRR    = [r1, r2, r3, r4]
	oRRout = objarr(4)
	for i = 0, 3 do begin
		;Get the variable
		R = MrVar_Get(oRR[i])
		
		;Class restriction
		if obj_class(R) ne objClass $
			then message, 'All R1-4 parameters must be the same object class.'
		
		;Interpolate if we need to
		if R['TIMEVAR'] -> IsIdentical(oTime) $
			then oRRout[i] = R -> Copy() $
			else oRRout[i] = R -> Interpol(oTime)
	endfor
	
	;Set data
	self.oR1 = oRRout[0]
	self.oR2 = oRRout[1]
	self.oR3 = oRRout[2]
	self.oR4 = oRRout[3]
end


;+
;   Compute the mean separation between vertices.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the result will be added to the cache.
;       NAME:           in, optional, type=string, default='RecipGrad'
;                       A name to be given to the variable.
;-
function MrVar_RecipVec::Separation, $
CACHE=cache
NAME=name
	compile_opt idl2
	on_error, 2
	
	if n_elements(name) eq 0 then name = 'rmag_mean'
	
	;Compute separation vectors
	oR12 = (self.oR2 - self.oR1) -> Magnitude()
	oR13 = (self.oR3 - self.oR1) -> Magnitude()
	oR14 = (self.oR4 - self.oR1) -> Magnitude()
	oR23 = (self.oR3 - self.oR2) -> Magnitude()
	oR24 = (self.oR4 - self.oR2) -> Magnitude()
	oR34 = (self.oR4 - self.oR3) -> Magnitude()
	
	;Compute mean separation
	oRR = (oR12 + oR13 + oR14 + oR23 + oR24 + oR34) / 6.0
	
	oRR -> SetName, name
	oRR -> AddAttr, 'CATDESC',    'Mean separation length between verticies of the tetrahedron.'
	oRR -> AddAttr, 'PLOT_TITLE', 'Mean separation length'
	oRR -> AddAttr, 'UNITS',      self.oR1['UNITS']
	oRR -> AddAttr, 'TITLE',      'R'
	if keyword_set(cache) then oRR -> Cache
	
	return, oRR
end


;+
;   The class definition statement.
;
; :Params:
;       CLASS:          out, optional, type=structure
;-
pro MrVar_RecipVec__DEFINE
	compile_opt idl2
	
	class = { MrVar_RecipVec, $
	          oR1: obj_new(), $
	          oR2: obj_new(), $
	          oR3: obj_new(), $
	          oR4: obj_new() $
	        }
end