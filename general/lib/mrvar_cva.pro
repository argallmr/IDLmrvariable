; docformat = 'rst'
;
; NAME:
;       MrVar_CVA
;+
;   Compute the velocity of a 1D boundary via the Constant Velocity Approach (CVA).
;   Inputs R1-4 should be time series vectors. These positions are then averaged in
;   time and the mean position of each spacecraft during the interval is used to compute
;   the boundary velocity.
;
; METHOD:
;    The velocity of the boundary can be expressed as a polynomial in time
;
;       V = A0 + A1*t + A2*t^2 + A3*t^3
;
;    If the velocity is constant, then A1 = A2 = A3 = 0. The distance traveled
;    along the normal from spacecraft 0 to spacecraft i can then be expressed
;    as the velocity multiplied by the delay times:
;
;       dR.n = A0 * dt
;
;    To solve, let m = n/A0 and invert the dR matrix. Noting that n is a unit
;    vector, |n| = 1, we find |m| = 1/A0 and
;
;      V = A0 = 1/|m|
;
;
; :Categories:
;   Physics Utilties
;
; :Params:
;       R1:         in, required, type=string/int/objref
;                   Name, number or objref of a MrVectorTS variable containing the
;                       position vectors of the first spacecraft.
;       R2:         in, required, type=string/int/objref
;                   Name, number or objref of a MrVectorTS variable containing the
;                       position vectors of the second spacecraft.
;       R3:         in, required, type=string/int/objref
;                   Name, number or objref of a MrVectorTS variable containing the
;                       position vectors of the third spacecraft.
;       R4:         in, required, type=string/int/objref
;                   Name, number or objref of a MrVectorTS variable containing the
;                       position vectors of the fourth spacecraft.
;       T:          in, required, type=FltArr(4)
;                   Times at which each spacecraft crossed the boundary, in seconds,
;                       ordered as [t1, t2, t3, t4].
;
; :Keywords:
;       N:          out, optional, type=FltArr(3)
;                   Components of the unit vector normal to the boundary.
;
; :Returns:
;       V:          Velocity of the boundary.
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
;       2018/02/12  -   Written by Matthew Argall
;-
FUNCTION MrVar_CVA, R1, R2, R3, R4, t1, t2, t3, t4, $
N=n
	Compile_Opt idl2
	On_Error, 2
	
	;Get the variables
	oR1 = MrVar_Get(R1)
	oR2 = MrVar_Get(R2)
	oR3 = MrVar_Get(R3)
	oR4 = MrVar_Get(R4)
	
	;Check that they are vectors
	IF ~Obj_IsA(oR1, 'MrVectorTS') THEN Message, 'R1 must be a MrVectorTS variable.'
	IF ~Obj_IsA(oR2, 'MrVectorTS') THEN Message, 'R2 must be a MrVectorTS variable.'
	IF ~Obj_IsA(oR3, 'MrVectorTS') THEN Message, 'R3 must be a MrVectorTS variable.'
	IF ~Obj_IsA(oR4, 'MrVectorTS') THEN Message, 'R4 must be a MrVectorTS variable.'
	
	;They must all have the same time stamps
	IF ~oR2 -> IsTimeIdentical(oR1) THEN Message, 'R2 and R1 must have the same time stamps.'
	IF ~oR3 -> IsTimeIdentical(oR1) THEN Message, 'R3 and R1 must have the same time stamps.'
	IF ~oR4 -> IsTimeIdentical(oR1) THEN Message, 'R4 and R1 must have the same time stamps.'
	
	;Separation vectors
	R = [ [Mean(oR1['DATA'], DIMENSION=1)], $
	      [Mean(oR2['DATA'], DIMENSION=1)], $
	      [Mean(oR3['DATA'], DIMENSION=1)], $
	      [Mean(oR4['DATA'], DIMENSION=1)] ]

	;Positions and crossing times relative to the first
	;spacecraft to encounter the boundary.
	dR = R[*,1:*] - rebin(R[*,0], 3, 4)
	dt = [t2, t3, t4] - t1

	;Invert the matrix dR
	;   - LA_INVERT expects the separation vectors to be along the
	;     rows, like a normal math matrix.
;	dRinv = invert(dR)
	dRinv = LA_Invert(dR)

	;Compute the normal vector
	;   - Multiply the rows of DT by the columns of DRINV
	;   - dRinv ## transpose(dt)
	m = Matrix_Multiply(dt, dRinv, /ATRANSPOSE)
	A = 1.0 / Sqrt(Total(m^2))
	n = m * A
	
	;Constant velocity
	V = Temporary(A)
	RETURN, V
END