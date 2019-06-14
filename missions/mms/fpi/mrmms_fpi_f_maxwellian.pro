; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FPI_Moms_Compare
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
;+
;   Construct a Maxwellian distribution with the same density, bulk velocity, and
;   temperature as the given distribution.
;
; :Params:
;       F:          in, required, type=string/int/objref
;                   The name, number, or MrTimeSeries object of the measured distribution
;                       function.
;       SPECIES:    in, required, type=string, default='e'
;                   Particle species. Options are {'e' | 'i'}
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output distribution will be added to the MrVariable cache.
;       DENSITY:        in, optional, type=string/int/objref
;                       The name, number, or MrScalarTS object of the measured density. If
;                           not provided, it will be caluclated by integrating `F`.
;       NAME:           in, optional, type=string, default='_maxwellian'
;                       The name of the output variable. If not provided, "_maxwellian" is
;                           appended to the name of the input distribution function.
;       SCPOT:          in, optional, type=string/int/objref
;                       The name, number, or MrScalarTS object of the measured spacecraft
;                           potential. Used only if `DENSITY`, `VELOCITY`, or `TEMPERATURE`
;                           are being calculated from `F`.
;       TEMPERATURE:    in, optional, type=string/int/objref
;                       The name, number, or MrScalarTS object of the measured scalar
;                           temperature. If not provided, it will be caluclated by
;                           integrating `F`.
;       VELOCITY:       in, optional, type=string/int/objref
;                       The name, number, or MrVectorTS object of the measured velocity. If
;                           not provided, it will be caluclated by integrating `F`.
;
; :Categories:
;    MMS
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
;       2018/06/21  -   Written by Matthew Argall
;-
FUNCTION MrMMS_FPI_F_Maxwellian, f, species, $
CACHE=cache, $
DENSITY=density, $
NAME=name, $
SCPOT=scpot, $
TEMPERATURE=temperature, $
VELOCITY=velocity
	Compile_Opt idl2
	On_Error, 2
	
	;Grab the data
	oDist = MrVar_Get(f)
	
	;Defaults
	IF N_Elements(species) EQ 0 THEN species = 'e'
	IF N_Elements(name)    EQ 0 THEN name    = oDist.name + '_maxwellian'
	
	;Constants
	eV2K = MrConstants('eV2K')
	eV2J = MrConstants('eV2J')
	kB   = MrConstants('k_B')
	mass = species EQ 'i' ? MrConstants('m_H') : MrConstants('m_e')
	
;-------------------------------------------
; Compute Moments //////////////////////////
;-------------------------------------------
	theSpecies = species EQ 'e' ? species : 'H'
	oDF        = MrDist4D(oDist, VSC=scpot, SPECIES=theSpecies)
	
	;Density
	IF N_Elements(n) EQ 0 $
		THEN oN = oDF -> Density() $
		ELSE oN = MrVar_Get(density)
	
	;Velocity
	IF N_Elements(v) EQ 0 $
		THEN oV = oDF -> Velocity() $
		ELSE oV = MrVar_Get(velocity)
	
	;Temperature
	IF N_Elements(t) EQ 0 THEN BEGIN
		oT = oDF -> Temperature()
		oT = (oT[*,0,0] + oT[*,1,1] + oT[*,2,2]) / 3.0
	ENDIF ELSE BEGIN
		oT = MrVar_Get(temperature)
	ENDELSE
	
	Obj_Destroy, oDF
	
;-------------------------------------------
; Compute Maxwellian ///////////////////////
;-------------------------------------------

	;Dimension sizes
	dims = Size(oDist, /DIMENSIONS)
	fmax = FltArr(dims)
	nt   = dims[0]
	np   = dims[1]
	nh   = dims[2]
	nv   = dims[3]

	;Dependent data
	oPhi   = oDist['DEPEND_1']
	oTheta = oDist['DEPEND_2']
	oE     = oDist['DEPEND_3']
	v      = Sqrt(2.0*eV2J/mass*oE['DATA'])  ;m/s

	;Fill the distribution with a Maxwellian with equivalent density and temperature
	FOR i = 0, np - 1 DO BEGIN
		FOR j = 0, nh - 1 DO BEGIN
			FOR k = 0, nv - 1 DO BEGIN
				vxsqr = (-v[*,k] * Sin(oTheta['DATA',j]*!DTOR) * Cos(oPhi['DATA',*,i]*!DTOR) - (1e3*oV['DATA',*,0]))^2
				vysqr = (-v[*,k] * Sin(oTheta['DATA',j]*!DTOR) * Sin(oPhi['DATA',*,i]*!DTOR) - (1e3*oV['DATA',*,1]))^2
				vzsqr = (-v[*,k] * Cos(oTheta['DATA',j]*!DTOR)                               - (1e3*oV['DATA',*,2]))^2
				fmax[*,i,j,k] = oN['DATA'] * 1e-6 * (mass / (2*!pi*kB*eV2K*oT['DATA']))^(3.0/2.0) * Exp(-mass*(vxsqr + vysqr + vzsqr) / (2.0*kB*eV2K*oT['DATA']))
			ENDFOR
		ENDFOR
	ENDFOR
	
	;Create the variable
	oFmax = MrTimeSeries( oDist['TIMEVAR'], fmax, $
	                      CACHE = cache, $
	                      NAME  = name, $
	                      /NO_COPY )
	oDist -> CopyAttrTo, oFmax

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, oFmax
END