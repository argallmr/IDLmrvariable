; docformat = 'rst'
;
; NAME:
;       MrVar_Dist1D_Prep
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
;   Extract 1D cuts along the principle axes of a 2D distribution function. Data in the
;   two bins on either size of +X, -X, +Y, and -Y are averaged together to make the ouput
;   data. The following attributes are used when extracting data.
;
;       Time-Dependent Distribtions:
;           'DEPEND_0' - Time variable
;           'DEPEND_1' - Azimuth angle bins variable
;           'DEPEND_2' - Energy bins variable
;
;       Time-Inependent Distribtions:
;           'DEPEND_0' - Azimuth angle bins variable
;           'DEPEND_1' - Energy bins variable
;
;       Energy:
;           'SI_CONVERSION'   - Conversion factor to SI units
;           'UNITS'           - Energy units. Used when SI_CONVERSION is not present
;           'DELTA_PLUS_VAR'  - Upper-bound of each energy bin. Converted to velocity delta.
;           'DELTA_PLUS'      - Uniform upper-bound of the energy bin (when DELTA_PLUS_VAR not present)
;           'DELTA_MINUS_VAR' - Lower-bound of each energy bin
;           'DELTA_MINUS'     - Uniform lower-bound of the energy bin (when DELTA_MINUS_VAR not present)
;
;       Phi:
;           'UNITS'           - Angle units. Used to convert to radians
;           'DELTA_PLUS_VAR'  - Upper-bound of each energy bin. Used to determin bin centers.
;           'DELTA_PLUS'      - Uniform upper-bound of the energy bin (when DELTA_PLUS_VAR not present)
;           'DELTA_MINUS_VAR' - Lower-bound of each energy bin
;           'DELTA_MINUS'     - Uniform lower-bound of the energy bin (when DELTA_MINUS_VAR not present)
;
;   Calling Sequence:
;       MrVar_Dist_Image(oX, oY, theDist, mass)
;       MrVar_Dist_Image(oX, oY, theDist, mass, idx)
;
; :Params:
;       oX:         out, required, type=objref
;                   A MrVariable containing the 1D cut along the X-axis (PHI=0,180).
;       oY:         out, required, type=objref
;                   A MrVariable containing the 1D cut along the X-axis (PHI=90,270).
;       THEDIST:    in, required, type=string/integer/objref
;                   A MrVariable object, name or index of distribution function
;                       to be displayed. Dimensions should be ordered as
;                       [time, azimuth, energy], with the time dimension being optional.
;                       Also must have corresponding DEPEND_[0-2] attributes.
;       MASS:       in, required, type=string/float
;                   Mass (kg) of the particles represented in the distribution. If a
;                       string, MrConstants('m_'+mass) will be used to obtain the mass.
;       IDX:        in, optional, type=integer, default=0
;                   If `THEDIST` is time-dependent, then this is the index along the
;                       time dimension of the distribution to be displayed.
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
;       2017/02/12  -   Written by Matthew Argall
;       2017/05/10  -   Added the CIRCLE keyword. - MRA
;       2018/05/24  -   Removed the IDX parameter to maintain time dependence of the
;                           distribution function. - MRA
;-
;*****************************************************************************************
;+
;
;-
PRO MrVar_Dist1D_Prep_CenterBins, oPhi
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------
; PHI Bin Centers //////////////////////////
;-------------------------------------------
	;
	; Assumes:
	;    - DELTA_(PLUS|MINUS)_VAR are arrays of the same size
	;    - DELTA_(PLUS|MINUS) are scalars
	;    - DELTA_MINUS(_VAR) is positive
	;    - All deltas have the same units as PHI
	;    - If no deltas are present, bin centers are assumed
	;

	;DELTA_PLUS
	IF oPhi -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		dPlus = MrVar_Get(oPhi['DELTA_PLUS_VAR'])
	ENDIF ELSE IF oPhi -> HasAttr('DELTA_PLUS') THEN BEGIN
		dPlus = oPhi['DELTA_PLUS']
	ENDIF ELSE BEGIN
		dPlus = 0
	ENDELSE
	
	;DELTA_MINUS
	IF oPhi -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		dMinus = MrVar_Get(oPhi['DELTA_MINUS_VAR'])
	ENDIF ELSE IF oPhi -> HasAttr('DELTA_MINUS') THEN BEGIN
		dMinus = oPhi['DELTA_MINUS']
	ENDIF ELSE BEGIN
		dMinus = 0
	ENDELSE

	;Center bins
	delta = (dPlus - dMinus) / 2.0
	oPhi -> SetData, oPhi['DATA'] + delta
	
	;Save deltas
	;   - If DELTA is a MrVariable, then it is a DELTA_(PLUS|MINUS)_VAR
	;   - Otherwise, it can be either a scalar or array
	;   - Deltas are always absolute deviations. Sign is implied with (PLUS|MINUS)
	IF Size(delta, /TNAME) NE 'OBJREF' THEN BEGIN
		IF MrIsA(delta, /SCALAR) THEN BEGIN
			oPhi['DELTA_PLUS']  = Abs(delta)
			oPhi['DELTA_MINUS'] = Abs(delta)
		ENDIF ELSE BEGIN
			delta = MrVariable( Abs(delta) )
			delta['UNITS'] = oPhi['UNITS']
			delta['TITLE'] = 'Half-width of angle bins.'
			oPhi['DELTA_PLUS_VAR']  = delta
			oPhi['DELTA_MINUS_VAR'] = delta
		ENDELSE
	ENDIF ELSE BEGIN
		oPhi['DELTA_PLUS_VAR']  = delta
		oPhi['DELTA_MINUS_VAR'] = delta
	ENDELSE 
END


;+
;
;-
FUNCTION MrVar_Dist1D_Prep_Circle, oDist, energy, $
CACHE=cache
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------
; Find Nearest Energy Bin //////////////////
;-------------------------------------------
	oE  = oDist['DEPEND_1']
	idx = MrNearestNeighbor(oE['DATA'], energy)

;-------------------------------------------
; Create Variable //////////////////////////
;-------------------------------------------
	;Reverse the order of the minus so the velcity maxima are at the edges
	oX = MrTimeSeries( oDist['TIMEVAR'], Reform(oDist['DATA',*,idx,*]), $
	                   CACHE     = cache, $
	                   NAME      = oDist.name + '_EPhi' )
	
;-------------------------------------------
; Set Attributes ///////////////////////////
;-------------------------------------------
	;PHI
	;   - Do not know if it is polar or azimuth angle
	oPhi = oDist['DEPEND_2']
	phi_units = oPhi -> HasAttr('UNITS') ? oPhi['UNITS'] : ''
	IF ~StRegEx(phi_units, 'deg', /BOOLEAN, /FOLD_CASE) THEN BEGIN
		oPhi -> SetData, oPhi['DATA'] * !radeg
		oPhi['UNITS'] = 'degrees'
		IF oPhi -> HasAttr('AXIS_RANGE') THEN oPhi['AXIS_RANGE'] *= !radeg
	ENDIF
	
	;DISTRIBUTION
	oX['DEPEND_1']   = oPhi
	oX['DIMENSION']  = 2B
	oX['LOG']        = 1B
	oX['PLOT_TITLE'] = 'f($\phi$)'
	
	;Legend labels
	IF N_Elements(energy) GT 1 THEN BEGIN
		E_units     = oE -> HasAttr('UNITS') ? oE['UNITS'] : ''
		oX['LABEL'] = String(energy, FORMAT='(i0)') + E_units
		oX['COLOR'] = MrDefaultColor(NCOLORS=N_Elements(energy))
	ENDIF
	
	;Distribution range
	iGT0 = oX -> Where(0, /GREATER, COUNT=nGT0)
	IF nGT0 GT 0 THEN oX['AXIS_RANGE'] = [Min( oX['DATA',*,iGT0], MAX=dmax ), dmax]
	
	;Other attributes
	IF oDist -> HasAttr('UNITS')         THEN oDist -> CopyAttrTo, oX, 'UNITS'
	IF oDist -> HasAttr('SI_CONVERSION') THEN oDist -> CopyAttrTo, oX, 'SI_CONVERSION'

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, oX
END


;+
;
;-
FUNCTION MrVar_Dist1D_Prep_Azimuth, oDist, $
CACHE=cache
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------
; Translate Angles /////////////////////////
;-------------------------------------------
	
	;Move from (-!pi, pi] to [0, 2*!pi)
	oPhi = oDist['DEPEND_2']
	iNeg = oPhi -> Where(0, /LESS, COUNT=nNeg)
	IF nNeg GT 0 THEN oPhi[iNeg] = 2*!pi + oPhi[iNeg]

;-------------------------------------------
; Find +/- X, Y ////////////////////////////
;-------------------------------------------
	;Get 0 +/- Delta
	;   - 0-Delta wraps to 360-Delta
	!Null = oPhi.Min(i0p, SUBSCRIPT_MAX=i0m)
	
	;Get 90 +/- Delta
	;   - Shift values by 270 to get the [0,360] wrapping effect
	!Null = Min( (oPhi['DATA'] + 3.0*!pi/2.0) MOD (2.0*!pi), i90p, SUBSCRIPT_MAX=i90m )
	
	;Get 180 +/- Delta
	;   - Shift values by 180 to get the [0,360] wrapping effect
	!Null = Min( (oPhi['DATA'] + !pi) MOD (2.0*!pi), i180m, SUBSCRIPT_MAX=i180p )
	
	;Get 270 +/- Delta
	;   - Shift values by 90 to get the [0,360] wrapping effect
	!Null = Min( (oPhi['DATA'] + !pi/2.0) MOD (2.0*!pi), i270p, SUBSCRIPT_MAX=i270m )
	
	;
	; Assume dPhi is uniform
	;

	;Average the two channels together
	xp = Mean( oDist[ 'DATA', *, *, [  i0p,   i0m] ], DIMENSION=3 )
	xm = Mean( oDist[ 'DATA', *, *, [i180p, i180m] ], DIMENSION=3 )
	yp = Mean( oDist[ 'DATA', *, *, [ i90p,  i90m] ], DIMENSION=3 )
	ym = Mean( oDist[ 'DATA', *, *, [i270p, i270m] ], DIMENSION=3 )

;-------------------------------------------
; Create Variables /////////////////////////
;-------------------------------------------
	dims = Size(oDist, /DIMENSIONS)
	nTime = dims[0]
	nVel  = dims[1]
	nPhi  = dims[2]
	
	;Reverse the order of the minus so the velcity maxima are at the edges
	oX = MrTimeSeries( oDist['TIMEVAR'], [ [ [Reverse(xm)], [Rebin(!Values.F_NaN, nTime)], [xp] ], $
	                                       [ [Reverse(ym)], [Rebin(!Values.F_NaN, nTime)], [yp] ] ], $
	                   CACHE     = cache, $
	                   NAME      = oDist.name + '_xy' )

;-------------------------------------------
; Set Attributes ///////////////////////////
;-------------------------------------------
	
	;Vx
	;   - Set mirror velocities to -X
	oV_temp             = oDist['DEPEND_1']
	oV                  = MrTimeSeries( oDist['TIMEVAR'], [-Reverse(oV_temp['DATA']), 0, oV_temp['DATA']] )
	oV['TITLE']         = 'Velocity (km/s)'
	oV['UNITS']         = 'km/s'
	oV['SI_CONVERSION'] = '1e3>m/s'
	
	;Distribution
	oX['COLOR']     = ['Blue', 'Red']
	oX['DEPEND_1']  = oV
	oX['DIMENSION'] = 1B
	oX['LABEL']     = ['X', 'Y']
	oX['LOG']       = 1B
	
	;Distribution range
	iGT0 = oX -> Where(0, /GREATER, COUNT=nGT0)
	IF nGT0 GT 0 THEN oX['AXIS_RANGE'] = [Min( oDist['DATA',iGT0], MAX=dmax ), dmax]
	
	;Other attributes
	IF oDist -> HasAttr('UNITS')         THEN oDist -> CopyAttrTo, oX, 'UNITS'
	IF oDist -> HasAttr('SI_CONVERSION') THEN oDist -> CopyAttrTo, oX, 'SI_CONVERSION'

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, oX
END


;+
;
;-
FUNCTION MrVar_Dist1D_Prep_Polar, oDist, $
CACHE=cache
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------
; Translate Angles /////////////////////////
;-------------------------------------------
	
	;Move from (-!pi/2, pi/2] to [0, !pi)
	oPhi = oDist['DEPEND_2']
	iNeg = oPhi -> Where(0, /LESS, COUNT=nNeg)
	IF nNeg GT 0 THEN oPhi[iNeg] = !pi + oPhi[iNeg]

;-------------------------------------------
; Find +/- X, Y ////////////////////////////
;-------------------------------------------
	;Get 0 & 180 Bins
	!Null = oPhi.Min(i0, SUBSCRIPT_MAX=i180)
	
	;Get 90 +/- Delta
	;   - Shift values by 90 to get the [0,360] wrapping effect
	!Null = Min( (oPhi['DATA'] + !pi/2.0) MOD !pi, i90m, SUBSCRIPT_MAX=i90p )
	
	;
	; Assume dPhi is uniform
	;

	;Average the two perpendicular channels together
	x0   = oDist[ 'DATA', *, *, i0 ]
	x180 = oDist[ 'DATA', *, *, i180 ]
	x90  = Mean( oDist[ 'DATA', *, *, [i90m, i90p] ], DIMENSION=3 )

;-------------------------------------------
; Create Variables /////////////////////////
;-------------------------------------------
	
	;Reverse the order of the minus so the velcity maxima are at the edges
	oX = MrVariable( oDist['TIMEVAR'], [ [[x0]], [[x90]], [[x180]] ], $
	                 CACHE     = cache, $
	                 NAME      = oDist.name + '_polarcuts' )

;-------------------------------------------
; Set Attributes ///////////////////////////
;-------------------------------------------
	;Energy
	oE = oDist['DEPEND_1']
	oE['TITLE'] = 'E' + (oE -> HasAttr('UNITS') ? ' (' + oE['UNITS'] + ')' : '')
	
	;DEPEND_0
	oX['COLOR']     = ['Black', 'Red', 'Blue']
	oX['DIMENSION'] = 1B
	oX['LABEL']     = ['0', '90', '180']
	oX['DEPEND_1']  = oDist['DEPEND_1']
	oX['LOG']       = 1B
	
	;AXIS_RANGE
	iGT0 = oX -> Where(0, /GREATER, COUNT=nGT0)
	IF nGT0 GT 0 THEN oX['AXIS_RANGE']   = [Min( oDist['DATA',iGT0], MAX=dmax ), dmax]
	
	IF oDist -> HasAttr('UNITS')         THEN oDist -> CopyAttrTo, oX, 'UNITS'
	IF oDist -> HasAttr('SI_CONVERSION') THEN oDist -> CopyAttrTo, oX, 'SI_CONVERSION'

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, oX
END


;+
;   Extract 1D cuts along the principle axes of a 2D distribution function. Data in the
;   two bins on either size of +X, -X, +Y, and -Y are averaged together to make the ouput
;   data. The following attributes are used when extracting data.
;
;   Calling Sequence:
;       oX = MrVar_Dist1D_Prep( theDist, mass[, idx] )
;       oX = MrVar_Dist1D_Prep( theDist, mass[, idx], /POLAR )
;       oX = MrVar_Dist1D_Prep( theDist, energy, mass[, idx], /CIRCLE )
;
; :Params:
;       DISTFN:     in, required, type=string/integer/objref
;                   A MrVariable object, name or index of distribution function
;                       to be displayed. Dimensions should be ordered as
;                       [time, azimuth, energy], with the time dimension being optional.
;                       Also must have corresponding DEPEND_[0-2] attributes.
;       ENERGY:     in, optional, type=number/numarr
;                   Energies (eV) at which `CIRCULAR` distributions are extracted. Should
;                       be included only if the `CIRCULAR` keyword is set.
;       MASS:       in, required, type=string/float
;                   Mass (kg) of the particles represented in the distribution. If a
;                       string, MrConstants(MASS) will be used to obtain the mass.
;       IDX:        in, optional, type=integer, default=0
;                   If `THEDIST` is time-dependent, then this is the index along the
;                       time dimension of the distribution to be displayed.
;
; :Keywords:
;       CACHE:      in, optional, type=boolean, default=0
;                   If set, returned variables are added to the variable cache.
;       CIRCULAR:   in, optional, type=boolean, default=0
;                   Describes the type of 1D cuts to create. If set, circular cuts around
;                       the core of the distribution are taken. The default is to take
;                       vertical and horizontal cuts through the core distribution.
;       POLAR:      in, optional, type=boolean, default=0
;                   Describes the orientation of `DISTFN`. If set, a polar-energy
;                       distribution is assumed. The default is to assume azimuth-energy.
;
; :Returns:
;       oX:         out, required, type=objref
;                   A MrVariable containing the 1D cut(s) of the distribution.
;-
FUNCTION MrVar_Dist1D_Prep, distfn, energy, mass, $
CACHE=cache, $
CIRCULAR=circular, $
NAME=name, $
POLAR=polar
	Compile_Opt idl2

	;Catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	
	;Defaults
	tf_cache  = Keyword_Set(cache)
	tf_circle = Keyword_Set(circular)
	tf_polar  = Keyword_Set(polar)
	
	IF N_Elements(name) EQ 0 THEN BEGIN
		suffix = '_cuts'
		IF tf_circle THEN suffix = '_circ' + suffix
		IF tf_polar  THEN suffix = '_polar' + suffix
		name = distfn.name + suffix
	ENDIF

;-------------------------------------------
; Order of Inputs //////////////////////////
;-------------------------------------------
	; oX = MrVar_Dist1D_Prep( theDist, energy, mass, /CIRCLE )
	IF tf_circle THEN BEGIN
		theEnergy = energy
		theMass   = mass
	
	; oX = MrVar_Dist1D_Prep( theDist, mass  [, /POLAR] )
	ENDIF ELSE BEGIN
		theMass = energy
		IF N_Elements(mass) GT 0 THEN theIdx = mass
	ENDELSE

;-------------------------------------------
; Create 1D Distributions //////////////////
;-------------------------------------------
	
	;Perform the necessary conversions
	;   - Pick a single distribution
	;   - Energy  --> Velocity
	;   - Degrees --> Radians
	oDist2D = MrVar_Dist2D_Prep(distfn, theMass, ENERGY=(tf_polar || tf_circle))

	;Center angular bins
	MrVar_Dist1D_Prep_CenterBins, oDist2D['DEPEND_2']
	
	;Extract data
	IF tf_circle THEN BEGIN
		oX = MrVar_Dist1D_Prep_Circle( oDist2D, theEnergy, CACHE=cache )
	ENDIF ELSE BEGIN
		IF tf_polar $
			THEN oX = MrVar_Dist1D_Prep_Polar( oDist2D, CACHE=cache ) $
			ELSE oX = MrVar_Dist1D_Prep_Azimuth( oDist2D, CACHE=cache )
	ENDELSE
	
	;Give it a name
	oX -> SetName, name
	
	RETURN, oX
END