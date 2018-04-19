; docformat = 'rst'
;
; NAME:
;   MrDist4D__Define
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
;   Calculate moments and reduced distributions from a 3D distribution FUNCTION.
;
;   REFERENCES:
;       [1] CIS Interface Control Document
;               http://caa.estec.esa.int/documents/ICD/CAA_CIS_ICD_V3.4.2.pdf
;               http://www.cosmos.esa.int/web/csa/documentation
;       [2] FPI Dataset
;               https://lasp.colorado.edu/mms/sdc/public/datasets/fpi/
;
; :Categories:
;   MrVariable, MrTimeSeries, MrDist
;
; :See Also:
;   MrDist3D__Define.PRO
;   MrTimeSeries__Define.PRO
;   MrVariable__Define.PRO
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
;-
;*****************************************************************************************
;+
;   The initialization method.
;
;   CALLING SEQUENCE:
;       oDist = MrDist4D( data )
;       oDist = MrDist4D( time, data )
;       oDist = MrDist4D( time, data, phi, theta, energy )
;
; :Params:
;       TIME:           in, optional, type=NxM array
;                       Name or reference of a MrTimeVar object, or an array
;                           of time stamps. IF a name is provided, the assiciated
;                           variable must exist in the variable cache.
;       DIST4D:         in, required, type=NxM array
;                       Name or reference of a MrVariable object, or the dependent
;                           variable data. IF a name is given, the associated variable
;                           must exist in the variable cache.
;       PHI:            in, optional, type=Nx1 or NxM array
;                       Azimuthal coordinates of the distribution pixels. Can be the name
;                           or reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           It must have dimensions of [phi, theta] or [time, phi, theta]
;       THETA:          in, optional, type=Nx1 or NxM array
;                       Polar coordinates of the distribution pixels. Can be the name or
;                           reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           It must have dimensions of [phi, theta] or [time, phi, theta]
;       ENERGY:         in, optional, type=Nx1 or NxM array
;                       Energy coordinates of the distribution pixels. Can be the name
;                           or reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           IF data has two dimensions, one must be time and the other
;                           must be the same Size as the fourth dimension of the
;                           distribution.
;
; :Keywords:
;       DEGREES:        in, optional, type=boolean, default=0
;                       If set, angles are given in degrees.
;       ELEVATION:      in, optional, type=boolean, default=0
;                       If set, `THETA` is taken to be the elevation angle. By default,
;                           it is interpreted as the polar angle.
;       MASS:           in, optional, type=float
;                       Mass of the particle species represented in the distribution. If
;                           given, `SPECIES` will be determined automatically.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;       RADIANS:        in, optional, type=boolean, default=0
;                       If set, angles are given in radians.
;       SPECIES:        in, optional, type=string, default='e'
;                       Species of particle represented in the distribution. Options are:
;                           {'e', 'H', 'He', 'O'}. If given, `MASS` will be determined
;                           automatically. Takes precedence over `MASS`.
;       UNITS:          in, optional, type=string, default='PSD'
;                       Units of the distribution. Options are: {'PSD', 'EFLUX', 'DIFF FLUX'}
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrTimeSeries::Init is also accepted here.
;-
FUNCTION MrDist4D::INIT, time, dist4D, phi, theta, energy, $
DEGREES=degrees, $
ELEVATION=elevation, $
MASS=mass, $
NAME=name, $
RADIANS=radians, $
SPECIES=species, $
UNITS=units, $
VSC=Vsc, $
_REF_EXTRA=extra
	Compile_Opt idl2

	;Error handling
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, 0
	ENDIF
	
	;Defaults
	IF N_Elements(name) EQ 0 THEN name = 'MrDist4D'
	IF N_Elements(species) EQ 0 && N_Elements(mass) EQ 0 THEN species = 'e'

;-----------------------------------------------------
; Initialize \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	self.oDist = MrTimeSeries( NAME         = name, $
	                          _STRICT_EXTRA = extra )
	IF ~Obj_Valid(self.oDist) THEN Message, 'Unable to initialize distribution function object property.'

;-----------------------------------------------------
; Set Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Distribution function
	IF N_Elements(time) GT 0 THEN BEGIN
		self -> SetData, time, dist4D, phi, theta, energy, $
		                 UNITS   = units, $
		                 DEGREES = degrees, $
		                 RADIANS = radians
	ENDIF
	
	;Spacecraft potential
	IF N_Elements(Vsc) GT 0 $
		THEN self -> SetVsc, Vsc $
		ELSE self.oVsc = MrScalarTS()

;-----------------------------------------------------
; Set Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	self -> SetProperty, ELEVATION = elevation, $
	                     MASS      = mass, $
	                     SPECIES   = species

	RETURN, 1
END


;+
;   Clean up after the object is destroyed
;-
PRO MrDist4D::CLEANUP
	Compile_Opt idl2
	On_Error, 2
	
	;Destroy objects
	Obj_Destroy, self.oDist
	Obj_Destroy, self.oVsc
END


;+
;   Convert the distribution function from one set of units to another. The conversion
;   factors are:
;                                        PARTICLE
;           |  PARTICLE FLUX (F)  |  ENERGY FLUX (FE)  |  PHASE SPACE DENSITY (f)  |
;      -----------------------------------------------------------------------------
;       F   |         1           |      1 / E         |          2*E/m^2          |
;       FE  |         E           |         1          |         2*(E/m)^2         |
;       f   |     m^2/(2*E)       |    m^2/(2*E^2)     |             1             |
;      -----------------------------------------------------------------------------
;
;   Such that f = [ m^2 / (2 * E^2) ] * FE, etc.
;
; :Params:
;       FLUX:           in, required, type=FltArr
;                       Distribution function for which to convert units
;       TO_UNITS:       in, required, type=string
;                       Convert to these units.  Options are::
;                           Energy:                'ENERGY'      - eV
;                           Particle Energy Flux:  'EFLUX'       - eV / cm^2 / s / sr / eV    or    keV / cm^2 / s / sr / keV
;                           Particle Flux:         'DIFF FLUX'   - # / cm^2 / s / sr / keV
;                           Phase Space Density:   'PSD'         - s^2 / cm^6
;                                                  'DF'          - s^2 / cm^6
;       FROM_UNITS:     in, optional, type=string, default='PSD'
;                       Name of the units of `FLUX`.
;
; :Keywords:
;       SPECIES:        in, optional, type=string, default='e'
;                       Particle species measured in the distribution. Particle mass
;                           is taken into account when converting 'PSD' <-> {'EFLUX' | 'DIFF FLUX'}
;       ENERGY:         in, optional, type=float
;                       Energy bins at which the distribution is taken.
;
; :Returns:
;       NEW_FLUX:       Distribution with new units.
;-
PRO MrDist4D::ConvertUnits, to_units
	Compile_Opt idl2
	On_Error, 2

	toUnits = strupcase(to_units)

	;
	; Conversion factor from particle energy flux to phase space density
	;   - (m^2 / 2 * E^2 ) * FE = [ kg^2 / (2 * eV^2) ]                         * (eV / cm^2 / s / sr / eV)
	;                           = ( kg^2 / eV^2 )                               / ( cm^2 * s * sr ) * 0.5
	;                           = [ kg^2 / ( eV * (1.602e-19 J/eV) ) ]          / ( cm^2 * s * sr ) * 0.5
	;                           = [ kg^2 / ( kg m^2 / s^2 )^2 ]                 / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2)
	;                           = [ (N * 1.672e-27 kg)^2 / ( kg^2 m^4 / s^4 ) ] / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2)
	;                           = [ kg^2 / ( kg^2 m^4 / s^4 ) ]                 / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2) * N^2 * 1.672e-27^2
	;                           = [ 1 / (m^4 * 1e8 /cm^4/m^4 / s^4) ]           / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2) * N^2 * 1.672e-27^2
	;                           = [ 1 / (cm^4 / s^4) ]                          / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2) * N^2 * 1.672e-27^2 * 1e-8
	;                           = ( s^4 / cm^4 )                                / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2) * N^2 * 1.672e-27^2 * 1e-8
	;                           = s^3 / cm^6 * N^2 * 1.602e-19^(-2) 1.672e-27^2 * 1e-8
	;                           = s^3 / cm^6 * N^2 * 5.44933e-25
	;                           = s^3 / m^6  * N^2 * 5.44933e-13
	;                           = s^3 / km^6 * N^2 * 5.44933e+5
	;
	
	;Mass number
	CASE self.species of
		'H':  N = 1
		'He': N = 2
		'O':  N = 16
		'e':  N = MrConstants('m_e') / MrConstants('m_p')
		ELSE: Message, 'Species not recognized.'
	ENDCASE
	
	;Must still divide by the values of E^2!!
	eflux_to_psd = N^2 * 5.44933e-25
	
	;Vectorize multiplication
	oEnergy = self.oDist['DEPEND_3']
	dims    = Size(self.oDist, /DIMENSIONS)
	IF Obj_IsA(oEnergy, 'MrTimeSeries') $
		THEN tempE = rebin( reform( oEnergy['DATA'], dims[0], 1, 1, dims[3] ), dims ) $
		ELSE tempE = rebin( reform( oEnergy['DATA'],       1, 1, 1, dims[3] ), dims )
	
	;Energy Flux
	IF self.units EQ 'EFLUX' THEN BEGIN
		
		;Convert to:
		CASE toUnits of
			'DIFF FLUX': new_flux = self['DATA'] / temporary(tempE) * 1e3     ; 1/eV * (1e3 eV/keV) = 1e3/keV
			'PSD':       new_flux = eflux_to_psd * self['DATA'] / temporary(tempE)^2
			'DF':        new_flux = eflux_to_psd * self['DATA'] / temporary(tempE)^2
			'EFLUX':     new_flux = self['DATA']
			ELSE: Message, 'Cannot convert from "' + self.units + '" to "' + to_units + '".'
		ENDCASE
	
	;Differential flux
	ENDIF ELSE IF self.units EQ 'DIFF FLUX' THEN BEGIN
		;Convert from PF to PEF
		eflux = self.oDist['DATA'] * tempE * 1e-3     ; 1/keV * (1e-3 keV/eV) = 1e-3/eV
		
		;Convert to:
		CASE toUnits of
			'EFLUX':     new_flux = temporary(eflux)
			'PSD':       new_flux = eflux_to_psd * temporary(eflux) / temporary(tempE)^2
			'DF':        new_flux = eflux_to_psd * temporary(eflux) / temporary(tempE)^2
			'DIFF FLUX': new_flux = self['DATA']
			ELSE: Message, 'Cannot convert from "' + self.units + '" to "' + to_units + '".'
		ENDCASE
	
	;Phase space density
	ENDIF ELSE IF self.units EQ 'PSD' || self.units EQ 'DF' THEN BEGIN
		;Convert to:
		CASE toUnits of
			'EFLUX':     new_flux = self.oDist['DATA'] / eflux_to_psd * temporary(tempE)^2
			'DIFF FLUX': new_flux = self.oDist['DATA'] / eflux_to_psd * temporary(tempE) * 1e3   ; 1/eV * (1e3 eV/keV) = 1e3/keV
			'PSD':       new_flux = self.oDist['DATA']
			'DF':        new_flux = self.oDist['DATA']
			ELSE: Message, 'Cannot convert from "' + self.units + '" to "' + to_units + '".'
		ENDCASE
	
	;Invalid
	ENDIF ELSE BEGIN
		Message, 'Invalid FROM_UNIT: "' + self.units + '".'
	ENDELSE

;-------------------------------------------
; Set Properties ///////////////////////////
;-------------------------------------------
	
	;Set properties
	self.units = toUnits
	CASE toUnits of
		'EFLUX':     self.oDist['UNITS'] = 'keV / cm^2 / s / sr / keV'
		'ENERGY':    self.oDist['UNITS'] = 'eV'
		'DIFF FLUX': self.oDist['UNITS'] = '# / cm^2 / s / sr / keV'
		'PSD':       self.oDist['UNITS'] = 's^2 / cm^6'
		ELSE: Message, 'Invalid units: "' + to_units + '".'
	ENDCASE
	
	;Set the object properties
	self.oDist -> SetData, Temporary(new_flux)
END


;+
;   Find the spacing of energy bins assuming constant dE/E
;
; :Returns:
;       DE:             Spacing of the energy bins.
;-
FUNCTION MrDist4D::DeltaE, energy
	Compile_Opt idl2
	On_Error, 2

	;
	; dLogE is the spacing of the energy bins in log space
	;
	;   dLogE = log(E1) - log(E0)                      (1)
	;
	; dE/E is just the derivative of a natural log
	;
	;   d( ln(E) ) = dE / E                            (2)
	;
	; If we change bases by log_a(B) = log_x(B) / log_x(A), where x=10, A=exp, B=E
	;
	;   d( ln(E) ) = d [ log(E) / log(exp) ]
	;              = dLogE / log(exp)                  (3)
	;
	; Combining Equations 1 and 2
	;
	;   dE/E = dLogE / log(exp)
	;
	;     dE = E * dLogE / log(exp)
	;
	dims = Size(energy, /DIMENSIONS)
	
	dLogE = ALog10( (energy)[*,1] ) - ALog10( (energy)[*,0] )
	dE    = energy * Rebin(dLogE, dims) / ALog10( Exp(1) )
	
	;Return results
	RETURN, dE
END


;+
;   Compute the 0th moment of the distribution, density.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       ON:             out, required, type=MrScalarTS
;                       Density as a function of time.
;-
FUNCTION MrDist4D::Density, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_Density'
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Constants
	deg2rad = !dpi / 180D
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	dims    = Size(self.oDist, /DIMENSIONS)
	
;-----------------------------------------------------
; Integrate over Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract phi and the bin sizes
	oPhi = self.oDist['DEPEND_1']
	dPhi = (oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR'])['DATA'] * deg2rad

	;Integrate
	ftemp = Total( self.oDist['DATA'] * Rebin(dPhi, dims), 2 )
	
	;Clean up
	dPhi = !Null
	
;-----------------------------------------------------
; Integrate over Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract theta and the bin sizes
	oTheta = self.oDist['DEPEND_2']
	dTheta = (oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR'])['DATA'] * deg2rad

	;Integrate
	ftemp = Total( ftemp * Rebin(Sin(oTheta['DATA']*deg2rad) * dTheta, dims[[0,2,3]]), 2)
	
	;Clean up
	dTheta = !Null
	
;-----------------------------------------------------
; Integrate over Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract energy and the bin sizes
	oEnergy  = self.oDist['DEPEND_3']
	oDEnergy = oEnergy['DELTA_PLUS_VAR'] + oEnergy['DELTA_MINUS_VAR']
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	; Convert energy from eV to J
	v    = Sqrt(2.0 * eV2J * oEnergy['DATA'] / self.mass)
	dv   = eV2J * oDEnergy['DATA'] / (self.mass * v)
	
	;Clean up
	Obj_Destroy, oDEnergy
	
	;SPACECRAFT POTENTIAL
	IF N_Elements(self.oVsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		signQ = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc   = Sqrt( 2.0 * q * (self.oVsc['DATA']) / self.mass )
		VM    = v^2 * (1 + signQ * (vsc^2 / v^2))
	ENDIF ELSE BEGIN
		VM = v^2
	ENDELSE
	
	;Integrate velocity
	;   - V is in m/s while f is in s^3/cm^6
	;   - v^2 * dv --> (m/s)^3  --> (cm/s)^3 * 1e6
	N = Total( Temporary(ftemp) * Temporary(v) * Sqrt(VM) * Temporary(dv), 2 ) * 1e6

;-----------------------------------------------------
; Create Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Energy-time spectrogram
	oN = MrScalarTS( self.oDist['TIMEVAR'], N, $
	                 CACHE = tf_cache, $
	                 NAME  = name, $
	                 /NO_COPY )
	
	;Attributes
	oN['CATDESC']       = 'Number density computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oN['LOG']           = 1B
	oN['UNITS']         = 'cm^-3'
	oN['TITLE']         = 'N!C(cm^-3)'
	oN['SI_CONVERSION'] = '1e-6>m^-3'
	
	RETURN, oN
END


;+
;   Compute entropy.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       OS:             out, required, type=MrScalarTS
;                       Entropy as a function of time.
;-
FUNCTION MrDist4D::Entropy, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oS) GT 0 THEN Obj_Destroy, oS
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_Entropy'
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Convert to radians
	deg2rad = !dpi / 180D
	dims    = Size(self.oDist, /DIMENSIONS)
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	tf_vsc  = N_Elements(self.oVsc) GT 0
	
	;Check for zeros -- mess with ALog
	tf_nan = 0B
	iZero  = self.oDist -> Where(0, /EQUAL, COUNT=nZero)
	IF nZero GT 0 THEN BEGIN
		self.oDist[iZero] = !Values.F_NaN
		tf_nan = 1B
	ENDIF

;-----------------------------------------------------
; Integration Over Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Extract phi and the bin sizes
	oPhi = self.oDist['DEPEND_1']
;	phi  = Rebin(oPhi['DATA'] * deg2rad, dims)
	dPhi = Rebin((oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR'])['DATA'] * deg2rad, dims)
	
	;Integrate
	ftemp = Total( self.oDist['DATA'] * ALog(self.oDist['DATA']) * Temporary(dPhi), 2, NAN=tf_nan)
	
	;Put zeros back
	IF nZero GT 0 THEN self.oDist[iZero] = 0.0
	
;-----------------------------------------------------
; Integration Over Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract theta and the bin sizes
	oTheta = self.oDist['DEPEND_2']
	theta  = Rebin(oTheta['DATA'] * deg2rad, dims[[0,2,3]])
	dTheta = Rebin((oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR'])['DATA'] * deg2rad, dims[[0,2,3]])
	
	;Integrate Theta
	ftemp  = Total( Temporary(ftemp) * Sin(Temporary(theta)) * Temporary(dTheta), 2, NAN=tf_nan)

;-----------------------------------------------------
; Integration Over Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract phi and the bin sizes
	oEnergy = self.oDist['DEPEND_3']
	energy  = oEnergy['DATA']
	dEnergy = (oEnergy['DELTA_PLUS_VAR'] + oEnergy['DELTA_MINUS_VAR'])['DATA']
	
	;VELOCITY
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	; Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * Temporary(energy) / self.mass)
	dv   = eV2J * Temporary(dEnergy) / (self.mass * v)
	
	;SPACECRAFT POTENTIAL
	IF tf_vsc THEN BEGIN
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q / self.mass * (self.oVsc['DATA']) )
	ENDIF
	
	;Integrate velocity
	;   - Sign is charge dependent. E = q * V = 1/2 * m * v^2
	;   - V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6  --> s^3/m^6 * 1e12
	IF tf_Vsc $
		THEN H = Total( ftemp * v * Sqrt(v^2 + sgn*vsc^2) * dv, 2, NAN=tf_nan ) * 1e12 $
		ELSE H = Total( ftemp * v^2 * dv, 2, NAN=tf_nan ) * 1e12
	
	ftemp = !Null
	v     = !Null
	dv    = !Null
	vsc   = !Null

;-----------------------------------------------------
; Create Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------	
	;Energy-time spectrogram
	oS = MrScalarTS( self.oDist['TIMEVAR'], -MrConstants('k_B') * H, $
	                 CACHE = tf_cache, $
	                 NAME  = name, $
	                 /NO_COPY )
	
	;Attributes
	oS['CATDESC']       = 'Entropy computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oS['LOG']           = 1B
	oS['PLOT_TITLE']    = 'Entropy'
	oS['UNITS']         = 'J/K'
	oS['TITLE']         = 'S!C(J/K/m^3)'
	oS['SI_CONVERSION'] = '>'
	oS['VAR_NOTES']     = 'S = -kH, where k is Boltzman constant and H is the Boltzman ' + $
	                      'H-function: \Integral f ln(f) d^3v'
	
	RETURN, oS
END


;+
;   Compute entropy flux.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       OS:             out, required, type=MrScalarTS
;                       Entropy as a function of time.
;-
FUNCTION MrDist4D::EntropyFlux, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oS) GT 0 THEN Obj_Destroy, oS
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_EntropyFlux'
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Compute entropy and velocity
	;   - Convert velocity from km/s to m/s
	oS = self -> Entropy()
	oV = self -> Velocity()
	
	;Multiply each component of V by S
	oSf = oV * oS
	oSf -> SetName, name
	IF Keyword_Set(cache) THEN oSf -> Cache

;-----------------------------------------------------
; Create Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------	
	
	;Attributes
	oSf['CATDESC']       = 'Entropy flux computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oSf['LOG']           = 1B
	oSf['PLOT_TITLE']    = 'Entropy Flux'
	oSf['UNITS']         = 'J/K/s/m^2'
	oSf['TITLE']         = 'S!C(J / K / m^2 / s)'
	oSf['SI_CONVERSION'] = '>'
	oSf['VAR_NOTES']     = 'S = -kH, where k is Boltzman constant and H is the Boltzman ' + $
	                       'H-function: \Integral V f ln(f) d^3v'
	
	RETURN, oSf
END


;+
;   Display the energy bins
;-
PRO MrDist4D::EBins
	Compile_Opt idl2
	On_Error, 2

	;Grab the two energy tables
	oE0 = MrVar_Get(e0_vname)
	oE1 = MrVar_Get(e1_vname)
	
	;Create an index vector
	nEnergy = N_Elements(oE0)
	idx     = indgen(nEnergy)

	;Print a header
	print, 'Index', 'E0', 'E1', FORMAT='(a5, 5x, a2, 9x, a2)'
	
	;Print the energy table
	;   - Scalar integer 0 returns the object itself.
	;   - To RETURN the value at index zero, the index must be an array: index = [0]
	FOR i = 0, nEnergy-1 do print, idx[i], oE0[[i]], oE1[[i]], FORMAT='(2x, i2, 2x, f9.2, 2x, f9.2)'
END


;+
;   Reduce the 3D distribution FUNCTION to a 1D distribution in azimuth angle.
;
; :Keywords:
;       NE_BINS:        in, optional, type=integer
;                       Number of energy bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       PHI_RANGE:      in, optional, type=FltArr(2), default=[0.0\, 360.0]
;                       The range in azimuthal angle (degrees) over which to average.
;       THETA_RANGE:    in, optional, type=FltArr(2), default=[0.0\, 180.0]
;                       The range in polar angle (degrees) over which to average.
;
; :Returns:
;       OESPEC:         out, required, type=MrTimeSeries
;                       A 1D distribution in time, averaged over polar and azimuth angle.
;-
FUNCTION MrDist4D::ESpec, $
CACHE=cache, $
PHI_RANGE=phi_range, $
NAME=name, $
NE_BINS=nE_bins, $
THETA_RANGE=theta_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D) GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oESpec)  GT 0 THEN Obj_Destroy, oESpec
		IF N_Elements(oEBins)  GT 0 THEN Obj_Destroy, oEBins
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_ESpec'

	;Allocate memory
	dims      = Size(self.oDist, /DIMENSIONS)
	nTime     = dims[0]
	nPhi      = dims[1]
	nTheta    = dims[2]
	nEnergy   = dims[3]
	ESpec = FltArr( nTime, nEnergy )

	;Step over each time
	FOR i = 0, nTime - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		ESpec[i,*] = oDist3D -> ESpec( e_bins, dE, $
		                               NE_BINS     = nE_bins, $
		                               PHI_RANGE   = phi_range, $
		                               THETA_RANGE = theta_range )

		;Destroy the object
		Obj_Destroy, oDist3D
	ENDFOR
	
	;Energy-time spectrogram
	oESpec = MrTimeSeries( self.oDist['TIMEVAR'], ESpec, $
	                       CACHE = tf_cache, $
	                       NAME  = name, $
	                       /NO_COPY )
	
	;Ordinate
;	binName = name + '_EBins'
;	oEBins  = Size(e_bins, /N_DIMENSIONS) EQ 2 $
;	               ? MrTimeSeries( self.oTime, e_bins, NAME=binName, /NO_COPY ) $
;	               : MrVariable( e_bins, NAME=binName, /NO_COPY )
	
	;Energy attributes
;	oEBins -> AddAttr, 'UNITS', self.oEnergy['UNITS']
;	oEBins -> AddAttr, 'TITLE', 'Energy'

	;Energy bins have not changed
	;   - MUST ALSO UPDATE MRDIST3D::SPECE
	oEBins = self.oDist['DEPEND_3']

	;Sepctrogram attributes
	oESpec['DEPEND_1'] = oEBins
	oESpec['SCALE']    = 1B
	oESpec['LOG']      = 1B
	oESpec['UNITS']    = self.oDist['UNITS']
	
	RETURN, oESpec
END


;+
;   Extract a single distribution FUNCTION.
;
; :Params:
;       IDX:                in, required, type=integer
;                           Time index FOR which the 3D distribution is returned.
;
; :Returns:
;       DIST3D:             out, required, type=MrVariable object
;                           A 3D distribution FUNCTION.
;-
FUNCTION MrDist4D::GetDist3D, idx
	Compile_Opt idl2
	On_Error, 2

;-----------------------------------------------------
; Extract Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Turn off verboseness
	;   - PHI and THETA are not standard sizes
	;   - Subscripting oDist will quasi-fail and generate warnings
	IF Size(self.oDist['DEPEND_1'], /N_DIMENSIONS) EQ 3 THEN BEGIN
		self.oDist.verbose = 0B
	
		;Get the independent variables
		oDist   = self.oDist[idx,*,*,*]
		oPhi    = (oDist['DEPEND_1'])[idx,*,*]
		oTheta  = (oDist['DEPEND_2'])[idx,*,*]
		oEnergy = oDist['DEPEND_3']
	
		;Turn verboseness back on
		oDist.verbose = 1B
	ENDIF ELSE BEGIN
		oDist   = self.oDist[idx,*,*,*]
		oPhi    = oDist['DEPEND_1']
		oTheta  = oDist['DEPEND_2']
		oEnergy = oDist['DEPEND_3']
	ENDELSE
	IF N_Elements(self.oVsc) GT 0 THEN Vsc = self.oVsc['DATA', idx]

;-----------------------------------------------------
; Phi Deltas \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;PHI DELTA PLUS
	IF oPhi -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		dphi_plus = Reform( ( oPhi['DELTA_PLUS_VAR'] )['DATA'] )
	
	ENDIF ELSE IF oPhi -> HasAttr('DELTA_PLUS') THEN BEGIN
		dphi_plus = oPhi['DELTA_PLUS']
	ENDIF
	
	;PHI DELTA MINUS
	IF oPhi -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		dphi_minus = Reform( ( oPhi['DELTA_MINUS_VAR'] )['DATA'] )
	
	ENDIF ELSE IF oPhi -> HasAttr('DELTA_MINUS') THEN BEGIN
		dphi_minus = oPhi['DELTA_MINUS']
	ENDIF

;-----------------------------------------------------
; Theta Deltas \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;DELTA PLUS
	IF oTheta -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		dtheta_plus = Reform( ( oTheta['DELTA_PLUS_VAR'] )['DATA'] )
	
	ENDIF ELSE IF oTheta -> HasAttr('DELTA_PLUS') THEN BEGIN
		dtheta_plus = oTheta['DELTA_PLUS']
	ENDIF
	
	;DELTA MINUS
	IF oTheta -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		dtheta_minus = Reform( ( oTheta['DELTA_PLUS_VAR'] )['DATA'] )
	
	ENDIF ELSE IF oTheta -> HasAttr('DELTA_MINUS') THEN BEGIN
		dtheta_minus = oTheta['DELTA_MINUS']
	ENDIF

;-----------------------------------------------------
; Energy Deltas \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;DELTA PLUS
	IF oEnergy -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		dE_plus = Reform( ( oEnergy['DELTA_PLUS_VAR'] )['DATA'] )
	
	ENDIF ELSE IF oEnergy -> HasAttr('DELTA_PLUS') THEN BEGIN
		dE_plus = oEnergy['DELTA_PLUS']
	ENDIF
	
	;DELTA MINUS
	IF oEnergy -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		dE_minus = Reform( ( oEnergy['DELTA_PLUS_VAR'] )['DATA'] )
	
	ENDIF ELSE IF oEnergy -> HasAttr('DELTA_MINUS') THEN BEGIN
		dE_minus = oE['DELTA_MINUS']
	ENDIF

;-----------------------------------------------------
; Create Distribution \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Create the 3D distribution
	dist3D = MrDist3D( Reform( oDist['DATA'] ), $
	                   Reform( oPhi['DATA'] ), $
	                   Reform( oTheta['DATA'] ), $
	                   Reform( oEnergy['DATA'] ), $
	                   Temporary( Vsc ), $
	                   /DEGREES, $
	                   DENERGY_MINUS = dE_minus, $
	                   DENERGY_PLUS  = dE_plus, $
	                   DPHI_MINUS    = dphi_minus, $
	                   DPHI_PLUS     = dphi_plus, $
	                   DTHETA_MINUS  = dtheta_minus, $
	                   DTHETA_PLUS   = dtheta_plus, $
	                   ELEVATION     = self.elevation, $
	                   MASS          = self.mass, $
	                   UNITS         = self.units )

;-----------------------------------------------------
; Done \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	RETURN, Dist3D
END


;+
;   Transform a spherical coordinate grid into a cartesian coordinate grid.
;
; :Keywords:
;       ELEVATION:          out, optional, type=boolean
;                           If set, THETA is taken to be the elevation angle.
;       MASS:               out, optional, type=float
;                           Mass (kg) of the species represented in the distribution.
;       SPECIES:            out, optional, type=string
;                           The particle species represented in the distribution.
;       UNITS:              out, optional, type=string
;                           Units of the distribution FUNCTION.
;       _REF_EXTRA:         out, optional, type=any
;                           Any keyword accepted by MrTimeSeries::GetProperty
;-
PRO MrDist4D::GetProperty, $
ELEVATION=elevation, $
MASS=mass, $
SPECIES=species, $
UNITS=units, $
_REF_EXTRA=extra
	Compile_Opt idl2
	On_Error, 2
	
	IF arg_present(mass)      GT 0 THEN mass      = self.mass
	IF arg_present(elevation) GT 0 THEN elevation = self.elevation
	IF arg_present(species)   GT 0 THEN species   = self.species
	IF arg_present(units)     GT 0 THEN units     = self.units
	
	;Distribution properties
	IF N_Elements(extra) GT 0 THEN self.oDist -> GetProperty, _STRICT_EXTRA=extra
END


;+
;   Compute the third moment of the distribution (heat flux).
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       OP:             out, required, type=MrTimeSeries
;                       Heat flux tensor as a function of time.
;-
FUNCTION MrDist4D::HeatFlux, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D) GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oN)      GT 0 THEN Obj_Destroy, oN
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_Pressure'

	;Allocate memory
	nTime = self.oDist -> GetNPts()
	Q     = FltArr( nTime, 3 )

	;Step over each time
	FOR i = 0, nTime - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		Q[i,*] = oDist3D -> HeatFlux2()

		;Destroy the object
		Obj_Destroy, oDist3D
	ENDFOR
	
	;Energy-time spectrogram
	oQ = MrVectorTS( self.oDist['TIMEVAR'], Q, $
	                 CACHE = tf_cache, $
	                 NAME  = name, $
	                 /NO_COPY )
	
	;Attributes
	oQ['CATDESC']       = 'Heat flux computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oQ['LABEL']         = ['Qx', 'Qy', 'Qz']
	oQ['UNITS']         = 'W'
	oQ['PLOT_TITLE']    = 'Heat Flux'
	oQ['TITLE']         = 'Q!C(nW)'
	oQ['SI_CONVERSION'] = '1e-9>W'
	
	RETURN, oQ
END


;+
;   Compute the moments of the distribution function: density, velocity, pressure,
;   temperature, and heatflux. Their computations are inter-dependent.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       OP:             out, required, type=MrTimeSeries
;                       Pressure tensor as a function of time.
;-
PRO MrDist4D::Moments, $
CACHE=cache, $
DENSITY=oN, $
ENTROPY=oS, $
HEATFLUX=oQ, $
PRESSURE=oP, $
TEMPERATURE=oT, $
VELOCITY=oV, $
_REF_EXTRA=extra
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D) GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oN)      GT 0 THEN Obj_Destroy, oN
		RETURN
	ENDIF

;-----------------------------------------------------
; Compute Moments \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Allocate memory
	nTime = self.oDist -> GetNPts()
	N     = FltArr( nTime )
	S     = FltArr( nTime )
	V     = FltArr( nTime, 3 )
	P     = FltArr( nTime, 3, 3 )
	T     = FltArr( nTime, 3, 3 )
	Q     = FltArr( nTime, 3 )

	;Step over each time
	FOR i = 0, nTime - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		oDist3D -> Moments_v2, DENSITY       = n_temp, $
		                       ENTROPY       = s_temp, $
		                       HEATFLUX      = q_temp, $
		                       PRESSURE      = p_temp, $
		                       TEMPERATURE   = t_temp, $
		                       VELOCITY      = v_temp, $
		                       _STRICT_EXTRA = extra

		;Store data
		N[i]     = Temporary(n_temp)
		S[i]     = Temporary(s_temp)
		Q[i,*]   = Temporary(q_temp)
		P[i,*,*] = Temporary(p_temp)
		T[i,*,*] = Temporary(t_temp)
		V[i,*]   = Temporary(v_temp)

		;Destroy the object
		Obj_Destroy, oDist3D
	ENDFOR

;-----------------------------------------------------
; Density \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Energy-time spectrogram
	oN = MrScalarTS( self.oDist['TIMEVAR'], N, $
	                 CACHE = cache, $
	                 NAME  = self.oDist.name + '_density', $
	                 /NO_COPY )

	;Attributes
	oN['CATDESC']       = 'Number density computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oN['DISPLAY_TYPE']  = 'time_series'
	oN['FIELDNAM']      = 'Density'
	oN['FILLVAL']       = -1e31
	oN['FORMAT']        = 'F9.4'
	oN['LOG']           = 1B
	oN['PLOT_TITLE']    = 'Density'
	oN['UNITS']         = 'cm^-3'
	oN['TITLE']         = 'N!C(cm^-3)'
	oN['SCALETYP']      = 'log'
	oN['SI_CONVERSION'] = '1e-6>m^-3'
	oN['VALIDMIN']      = 0.0
	oN['VALIDMAX']      = 1e4
	oN['VAR_TYPE']      = 'data'

;-----------------------------------------------------
; Entropy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Entropy
	oS = MrScalarTS( self.oDist['TIMEVAR'], S, $
	                 CACHE = cache, $
	                 NAME  = self.oDist.name + '_entropy', $
	                 /NO_COPY )
	
	;Attributes
	oS['CATDESC']       = 'Entropy computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oS['DISPLAY_TYPE']  = 'time_series'
	oS['FIELDNAM']      = 'Entropy'
	oS['FILLVAL']       = -1e31
	oS['FORMAT']        = 'E14.5'
	oS['LOG']           = 1B
	oS['PLOT_TITLE']    = 'Entropy'
	oS['UNITS']         = 'J/K'
	oS['TITLE']         = 'S!C(J/K/m^3)'
	oS['SCALETYP']      = 'log'
	oS['SI_CONVERSION'] = '>'
	oS['VAR_NOTES']     = 'S = -kH, where k is Boltzman constant and H is the Boltzman ' + $
	                      'H-function: \Integral f ln(f) d^3v'
	oS['VAR_TYPE']      = 'data'

;-----------------------------------------------------
; Velocity \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Energy-time spectrogram
	oV = MrVectorTS( self.oDist['TIMEVAR'], V, $
	                 CACHE = cache, $
	                 NAME  = self.oDist.name + '_velocity', $
	                 /NO_COPY )
	
	;Attributes
	oV['CATDESC']       = 'Number density computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oV['DISPLAY_TYPE']  = 'time_series'
	oV['FIELDNAM']      = 'Velocity'
	oV['FILLVAL']       = -1e31
	oV['FORMAT']        = 'F14.6'
	oV['LABEL']         = ['Vx', 'Vy', 'Vz']
	oV['PLOT_TITLE']    = 'Velocity'
	oV['UNITS']         = 'km/s'
	oV['TITLE']         = 'V!C(km/s)'
	oV['SCALETYP']      = 'linear'
	oV['SI_CONVERSION'] = '1e3>m/s'
	oV['VAR_TYPE']      = 'data'

;-----------------------------------------------------
; Pressure \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Energy-time spectrogram
	oP = MrMatrixTS( self.oDist['TIMEVAR'], P, $
	                 CACHE = cache, $
	                 NAME  = self.oDist.name + '_pressure', $
	                 /NO_COPY )
	
	;Attributes
	oP['CATDESC']       = 'Pressure computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oP['DISPLAY_TYPE']  = 'time_series'
	oP['FIELDNAM']      = 'Pressure'
	oP['FILLVAL']       = -1e31
	oP['FORMAT']        = 'F11.5'
	oP['LABEL']         = 'P'
	oP['LABL_PTR_1']    = ['x', 'y', 'z']
	oP['LABL_PTR_2']    = ['x', 'y', 'z']
	oP['PLOT_TITLE']    = 'Pressure'
	oP['UNITS']         = 'nPa'
	oP['PLOT_TITLE']    = 'Pressure Tensor'
	oP['TITLE']         = 'P!C(nPa)'
	oP['SCALETYP']      = 'linear'
	oP['SI_CONVERSION'] = '1e-9>Pa'
	oP['VALIDMIN']      = 0.0
	oP['VAR_TYPE']      = 'data'

;-----------------------------------------------------
; Temperature \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Temperature tensor
	oT = MrMatrixTS( self.oDist['TIMEVAR'], T, $
	                 CACHE = cache, $
	                 NAME  = self.oDist.name + '_temperature', $
	                 /NO_COPY )
	
	;Attributes
	oT['CATDESC']       = 'Temperature computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oT['DISPLAY_TYPE']  = 'time_series'
	oT['FIELDNAM']      = 'Temperature'
	oT['FILLVAL']       = -1e31
	oT['FORMAT']        = 'F11.5'
	oT['LABEL']         = 'T'
	oT['LABL_PTR_1']    = ['x', 'y', 'z']
	oT['LABL_PTR_2']    = ['x', 'y', 'z']
	oT['PLOT_TITLE']    = 'Temperature'
	oT['UNITS']         = 'eV'
	oT['PLOT_TITLE']    = 'Temperature Tensor'
	oT['TITLE']         = 'T!C(eV)'
	oT['SCALETYP']      = 'linear'
	oT['SI_CONVERSION'] = '>'
	oT['VALIDMIN']      = 0.0
	oT['VAR_TYPE']      = 'data'

;-----------------------------------------------------
; Heat Flux \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Energy-time spectrogram
	oQ = MrVectorTS( self.oDist['TIMEVAR'], Q, $
	                 CACHE = cache, $
	                 NAME  = self.oDist.name + '_heatflux', $
	                 /NO_COPY )
	
	;Attributes
	oQ['CATDESC']       = 'Heat flux computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oQ['DISPLAY_TYPE']  = 'time_series'
	oQ['FIELDNAM']      = 'Heat Flux'
	oQ['FILLVAL']       = -1e31
	oQ['FORMAT']        = 'F11.5'
	oQ['LABEL']         = ['Qx', 'Qy', 'Qz']
	oQ['PLOT_TITLE']    = 'Heat Flux'
	oQ['UNITS']         = 'W'
	oQ['TITLE']         = 'Q!C(nW)'
	oQ['SI_CONVERSION'] = '1e-9>W'
	oQ['SCALETYP']      = 'linear'
	oQ['VAR_TYPE']      = 'data'
END


;+
;   Reduce the 3D distribution FUNCTION to a 2D distribution in polar angle and energy,
;   averaging over azimuth angle.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output variable will be added to the variable cache.
;       NAME:           in, optional, type=string, default=self.name + '_ThetaE'
;                       Name to be given to the output variable.
;       NE_BINS:        in, optional, type=integer
;                       Number of energy bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       THETA_RANGE:    in, optional, type=FltArr(2), default=[0.0\, 180.0]
;                       The range in polar angle (degrees) over which to average.
;
; :Returns:
;       DIST2D:         out, required, type=MrTimeSeries object
;                       A time-varying 2D distribution in polar angle and energy.
;-
FUNCTION MrDist4D::PhiE, $
CACHE=cache, $
NAME=name, $
NE_BINS=ne_bins, $
NPHI_BINS=nPhi_bins, $
THETA_RANGE=theta_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D)  GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oPhiE)    GT 0 THEN Obj_Destroy, oPhiE
		IF N_Elements(oPhiBins) GT 0 THEN Obj_Destroy, oPhiBins
		IF N_Elements(oEBins)   GT 0 THEN Obj_Destroy, oEBins
		RETURN, !Null
	ENDIF
	
	;Defaults
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_PhiE'

;-------------------------------------------
; Reduce the 4D Distribution ///////////////
;-------------------------------------------

	;Allocate memory
	dims    = Size(self.oDist, /DIMENSIONS)
	nTimes  = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	PhiE  = FltArr( nTimes, nPhi, nEnergy )

	;Step over each time
	FOR i = 0, nTimes - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		PhiE[i,*,*] = oDist3D -> PhiE( phi, energy, dPhi, dE, $
		                               NE_BINS     = ne_bins, $
		                               NPHI_BINS   = nPhi_bins, $
		                               THETA_RANGE = theta_range )

		;Destroy the 3D distribution
		Obj_Destroy, oDist3D
	ENDFOR

;-------------------------------------------
; Datasets /////////////////////////////////
;-------------------------------------------
	;Time variable
	oTime = self.oDist['TIMEVAR']
	
	;Theta-Energy distribution
	oPhiE = MrTimeSeries( oTime, PhiE, $
	                      CACHE = cache, $
	                      NAME  = name, $
	                      /NO_COPY )
	
	;Phi
	binName  = name + '_PhiBins'
	oPhiBins = Size(phi, /N_DIMENSIONS) EQ 2 $
	                 ? MrTimeSeries( oTime, phi, NAME=binName, /NO_COPY ) $
	                 : MrVariable( phi, NAME=binName, /NO_COPY )
	
	;Energy
	binName     = name + '_EnergyBins'
	oEnergyBins = Size(energy, /N_DIMENSIONS) EQ 2 $
	                  ? MrTimeSeries( oTime, energy, NAME=binName, /NO_COPY ) $
	                  : MrVariable( energy, NAME=binName, /NO_COPY )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	
	;Phi attributes
	oPhiBins['DELTA_MINUS'] = dPhi
	oPhiBins['DELTA_PLUS']  = dPhi
	oPhiBins['UNITS']       = 'degrees'
	oPhiBins['TITLE']       = 'Azimuth'
	oPhiBins['PLOT_TITLE']  = 'Azimuthal Bin Centers'
	
	;Energy attributes
;	oEBins['UNITS'] = self.oEnergy['UNITS']
;	oEBins['TITLE'] = 'Energy'

	;Energy bins have not changed
	;   - MUST ALSO UPDATE MRDIST3D::SPECE
	oEBins = self.oDist['DEPEND_3']

	;Distribution attributes
	oPhiE['DEPEND_1'] = oPhiBins
	oPhiE['DEPEND_2'] = oEBins
	oPhiE['SCALE']    = 1B
	oPhiE['LOG']      = 1B
	oPhiE['UNITS']    = self.oDist['UNITS']
	
	;RETURN the 2D distribution
	RETURN, oPhiE
END


;+
;   Reduce the 3D distribution FUNCTION to a 1D distribution in azimuth angle.
;
; :Keywords:
;       E_RANGE:        in, optional, type=FltArr(2), default=[min, max]
;                       The range in energy, in electron volts (eV) over which to average.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       THETA_RANGE:    in, optional, type=FltArr(2), default=[0.0\, 180.0]
;                       The range in polar angle (degrees) over which to average.
;
; :Returns:
;       OPHISPEC:       out, required, type=MrTimeSeries
;                       A 1D distribution in time, averaged over energy and polar angle.
;-
FUNCTION MrDist4D::PhiSpec, $
CACHE=cache, $
E_RANGE=E_range, $
NAME=name, $
NPHI_BINS=nPhi_bins, $
THETA_RANGE=theta_range, $
UNITS=units, $
WEIGHT=weight
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D)  GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oPhiSpec) GT 0 THEN Obj_Destroy, oPhiSpec
		IF N_Elements(oPhiBins) GT 0 THEN Obj_Destroy, oPhiBins
		RETURN, !Null
	ENDIF
	
	;Defaults
	IF N_Elements(units) EQ 0 THEN units = self.units
	IF N_Elements(name)  EQ 0 THEN name  = self.oDist.name + '_PhiSpec'
	
	;Velocity-space weights
	tf_weights = 0
	IF N_Elements(weights) GT 0 THEN BEGIN
		tf_weights = 1B
		oW         = MrVar_Get(weights)
	ENDIF

	;Allocate memory
	dims    = Size(self.oDist, /DIMENSIONS)
	nTime   = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	phiSpec = FltArr( nTime, nPhi )
	
	;Step over each time
	FOR i = 0, nTime - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		IF units NE self.units THEN oDist3D -> SetUnits, units
		
		;Get the weights
		IF tf_weights THEN w = Reform(oW[i,*,*,*])
		
		;Reduce the distribution
		phiSpec[i,*] = oDist3D -> PhiSpec_v2( phi_bins, dPhi, $
		                                      E_RANGE     = e_range, $
		                                      NPHI_BINS   = nPhi_bins, $
		                                      THETA_RANGE = theta_range, $
		                                      WEIGHT      = w )

		;Destroy the object
		Obj_Destroy, oDist3D
	ENDFOR
	
	;Time variable
	oTime = self.oDist['TIMEVAR']
	
	;Phi-time spectrogram
	oPhiSpec = MrTimeSeries( oTime, phiSpec, $
	                         CACHE = cache, $
	                         NAME  = name, $
	                         /NO_COPY )
	
	;Abscissa
	binName  = name + '_PhiBins'
	oPhiBins = Size(phi_bins, /N_DIMENSIONS) EQ 2 $
	                 ? MrTimeSeries( oTime, phi_bins, NAME=binName, /NO_COPY ) $
	                 : MrVariable( phi_bins, NAME=binName, /NO_COPY )
	
	;Phi attributes
	oPhiBins['DELTA_MINUS'] = dPhi
	oPhiBins['DELTA_PLUS']  = dPhi
	oPhiBins['UNITS']       = 'degrees'
	oPhiBins['TITLE']       = 'Azimuth'
	oPhiBins['PLOT_TITLE']  = 'Azimuthal Bin Centers'

	;Sepctrogram attributes
	oPhiSpec['DEPEND_1']   = oPhiBins
	oPhiSpec['SCALE']      = 1B
	oPhiSpec['LOG']        = 1B
	oPhiSpec['UNITS']      = self.oDist['UNITS']
	oPhiSpec['TITLE']      = 'Phi Dist'
	oPhiSpec['PLOT_TITLE'] = 'Distribution in Phi'
	
	RETURN, oPhiSpec
END


;+
;   Reduce the 3D distribution FUNCTION to a 2D distribution in azimuth and polar angles,
;   averaging over energy.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output variable will be added to the variable cache.
;       NAME:           in, optional, type=string, default=self.name + '_ThetaE'
;                       Name to be given to the output variable.
;       NE_BINS:        in, optional, type=integer
;                       Number of energy bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       THETA_RANGE:    in, optional, type=FltArr(2), default=[0.0\, 180.0]
;                       The range in polar angle (degrees) over which to average.
;
; :Returns:
;       DIST2D:         out, required, type=MrTimeSeries object
;                       A time-varying 2D distribution in polar angle and energy.
;-
FUNCTION MrDist4D::PhiTheta, $
CACHE=cache, $
NAME=name, $
E_RANGE=E_Range, $
NPHI_BINS=nPhi_bins, $
NTHETA_BINS=nTheta_bins
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D)    GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oPhiTheta)  GT 0 THEN Obj_Destroy, oPhiTheta
		IF N_Elements(oPhiBins)   GT 0 THEN Obj_Destroy, oPhiBins
		IF N_Elements(oThetaBins) GT 0 THEN Obj_Destroy, oThetaBins
		RETURN, !Null
	ENDIF
	
	;Defaults
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_PhiTheta'

;-------------------------------------------
; Reduce the 4D Distribution ///////////////
;-------------------------------------------

	;Allocate memory
	dims    = Size(self.oDist, /DIMENSIONS)
	nTimes  = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	PhiE    = FltArr( nTimes, nPhi, nTheta )
	
	;Step over each time
	FOR i = 0, N_Elements(nTimes) - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		PhiTheta[i,*,*] = oDist -> ThetaPhi( phi, theta, dPhi, dTheta, $
		                                     E_RANGE     = E_Range, $
		                                     NPHI_BINS   = nPhi_bins, $
		                                     NTHETA_BINS = nTheta_bins )
		
		;Destroy the 3D distribution
		Obj_Destroy, oDist3D
	ENDFOR

;-------------------------------------------
; Datasets /////////////////////////////////
;-------------------------------------------
	;Time variable
	oTime = self.oDist['TIMEVAR']
	
	;Theta-Energy distribution
	oPhiTheta = MrTimeSeries( oT, PhiTheta, $
	                          CACHE = cache, $
	                          NAME  = name, $
	                          /NO_COPY )
	
	;Phi
	binName  = name + '_PhiBins'
	oPhiBins = Size(phi, /N_DIMENSIONS) EQ 2 $
	                 ? MrTimeSeries( oTime, phi, NAME=binName, /NO_COPY ) $
	                 : MrVariable( phi, NAME=binName, /NO_COPY )
	
	;Theta
	binName    = name + '_ThetaBins'
	oThetaBins = Size(theta, /N_DIMENSIONS) EQ 2 $
	                 ? MrTimeSeries( oTime, theta_bins, NAME=binName, /NO_COPY ) $
	                 : MrVariable( theta_bins, NAME=binName, /NO_COPY )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	
	;Phi attributes
	oPhiBins['DELTA_MINUS'] = dPhi
	oPhiBins['DELTA_PLUS']  = dPhi
	oPhiBins['UNITS']       = 'degrees'
	oPhiBins['TITLE']       = 'Azimuth'
	oPhiBins['PLOT_TITLE']  = 'Azimuthal Bin Centers'
	
	;Theta attributes
	oThetaBins['DELTA_MINUS'] = dTheta
	oThetaBins['DELTA_PLUS']  = dTheta
	oThetaBins['UNITS']      = 'degrees'
	oThetaBins['TITLE']      = 'Polar Angle'
	oThetaBins['PLOT_TITLE'] = 'Polar Bin Centers'
	
	;Distribution attributes
	oPhiTheta['DEPEND_1'] = oPhiBins
	oPhiTheta['DEPEND_2'] = oThetaBins
	oPhiTheta['SCALE']    = 1B
	oPhiTheta['LOG']      = 1B
	oPhiTheta['UNITS']    = self.oDist['UNITS']
	
	;RETURN the 2D distribution
	RETURN, oThetaBins
END


;+
;   Compute the second moment of the distribution (pressure).
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       OP:             out, required, type=MrTimeSeries
;                       Pressure tensor as a function of time.
;-
FUNCTION MrDist4D::Pressure, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_Pressure'
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'

	;Convert to radians
	deg2rad = !dpi / 180D
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	dims    = Size(self.oDist, /DIMENSIONS)
	
;-----------------------------------------------------
; Integrate over Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	; Let
	;    VM = Vm^2 - sign(q) Vsc^2
	;
	; Where, here and in what follows,
	;    VM      =  The velocity corrected for spacecraft potential effects
	;    Vm      =  The measured velocity related to the energy bins of the instrument
	;    Vsc     =  The velocity gained by passing through the spacecraft potential sheath
	;    U       =  Plasma bulk velocity
	;    dVm     =  Size of a measured velocity space bin (energy bin width converted to velocity)
	;    dOmega  =  Sin(Theta) dTheta dPhi = Element of solid angle
	;    t       =  Theta, the polar angle
	;    p       =  Phi, the azimuthal angle
	;    m       =  Particle mass
	;    f       =  The distribution function
	;    q       =  Charge of the particle
	;
	; It follows that
	;
	;    Pxx = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 cos(p)^2      - 2 VM Ux sin(t) cos(p)                       + Sqrt(VM) Ux^2  ]
	;    Pxy = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 sin(p) cos(p) -   VM Uy sin(t) cos(p) - VM Ux sin(t) sin(p) + Sqrt(VM) Ux Uy ]
	;    Pxz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)   cos(t) cos(p) -   VM Uz sin(t) cos(p) - VM Ux cos(t)        + Sqrt(VM) Ux Uz ]
	;    Pyy = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 sin(p)^2      - 2 VM Uy sin(t) sin(p)                       + Sqrt(VM) Uy^2  ]
	;    Pyz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)   cos(t) sin(p) -   VM Uz sin(t) sin(p) - VM Uy cos(t)        + Sqrt(VM) Uy Uz ]
	;    Pzz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) cos(t)^2               - 2 VM Uz cos(t)                              + Sqrt(VM) Uz^2  ]
	;
	
	;We need the bulk velocity
	;   - Convert km/s to m/s
	oBulkV = 1e3 * Self -> Velocity()
	vx     = Rebin(oBulkV['DATA',*,0], dims)
	vy     = Rebin(oBulkV['DATA',*,1], dims)
	vz     = Rebin(oBulkV['DATA',*,2], dims)
	Obj_Destroy, oBulkV
	
	;Extract phi and the bin sizes
	oPhi = self.oDist['DEPEND_1']
	phi  = Rebin(oPhi['DATA'] * deg2rad, dims)
	dPhi = Rebin((oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR'])['DATA'] * deg2rad, dims)
	
	;Pxx
	fxx1_temp = Total( self.oDist['DATA'] * Cos(phi)^2 * dPhi,        2 )
	fxx2_temp = Total( self.oDist['DATA'] * Cos(phi)   * dPhi * vx,   2 ) * 2.0
	fxx3_temp = Total( self.oDist['DATA']              * dPhi * vx^2, 2 )
	
	;Pxy
	fxy1_temp = Total( self.oDist['DATA'] * Sin(phi) * Cos(phi) * dPhi,           2 )
	fxy2_temp = Total( self.oDist['DATA'] * Cos(phi)            * dPhi * vy,      2 )
	fxy3_temp = Total( self.oDist['DATA'] * Sin(phi)            * dPhi * vx,      2 )
	fxy4_temp = Total( self.oDist['DATA']                       * dPhi * vx * vy, 2 )
	
	;Pxz
	fxz1_temp = Total( self.oDist['DATA'] * Cos(phi) * dPhi,           2 )
	fxz2_temp = Total( self.oDist['DATA'] * Cos(phi) * dPhi * vz,      2 )
	fxz3_temp = Total( self.oDist['DATA']            * dPhi * vx,      2 )
	fxz4_temp = Total( self.oDist['DATA']            * dPhi * vx * vz, 2 )
	
	;Pyy
	fyy1_temp = Total( self.oDist['DATA'] * Sin(phi)^2 * dPhi,        2 )
	fyy2_temp = Total( self.oDist['DATA'] * Sin(phi)   * dPhi * vy,   2 ) * 2.0
	fyy3_temp = Total( self.oDist['DATA']              * dPhi * vy^2, 2 )
	
	;Pyz
	fyz1_temp = Total( self.oDist['DATA'] * Sin(phi) * dPhi,           2 )
	fyz2_temp = Total( self.oDist['DATA'] * Sin(phi) * dPhi * vz,      2 )
	fyz3_temp = Total( self.oDist['DATA']            * dPhi * vy,      2 )
	fyz4_temp = Total( self.oDist['DATA']            * dPhi * vy * vz, 2 )
	
	;Pzz
	fzz1_temp = Total( self.oDist['DATA'] * dPhi,        2 )
	fzz2_temp = Total( self.oDist['DATA'] * dPhi * vz,   2 ) * 2.0
	fzz3_temp = Total( self.oDist['DATA'] * dPhi * vz^2, 2 )
	
	phi  = !Null
	dPhi = !Null
	vx   = !Null
	vy   = !Null
	vz   = !Null

;-----------------------------------------------------
; Integrate over Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract theta and the bin sizes
	oTheta = self.oDist['DEPEND_2']
	theta  = Rebin(oTheta['DATA'] * deg2rad, dims[[0,2,3]])
	dTheta = Rebin((oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR'])['DATA'] * deg2rad, dims[[0,2,3]])
	
	;Pxx
	fxx1_temp = Total( fxx1_temp * Sin(theta)^3 * dTheta, 2 )
	fxx2_temp = Total( fxx2_temp * Sin(theta)^2 * dTheta, 2 )
	fxx3_temp = Total( fxx3_temp * Sin(theta)   * dTheta, 2 )
	
	;Pxy
	fxy1_temp = Total( fxy1_temp * Sin(theta)^3            * dTheta, 2 )
	fxy2_temp = Total( fxy2_temp * Sin(theta)^2            * dTheta, 2 )
	fxy3_temp = Total( fxy3_temp * Sin(theta) * Cos(theta) * dTheta, 2 )
	fxy4_temp = Total( fxy4_temp * Sin(theta)              * dTheta, 2 )
	
	;Pxz
	fxz1_temp = Total( fxz1_temp * Sin(theta)^2 * Cos(theta) * dTheta, 2 )
	fxz2_temp = Total( fxz2_temp * Sin(theta)^2              * dTheta, 2 )
	fxz3_temp = Total( fxz3_temp * Sin(theta)   * Cos(theta) * dTheta, 2 )
	fxz4_temp = Total( fxz4_temp * Sin(theta)                * dTheta, 2 )
	
	;Pyy
	fyy1_temp = Total( fyy1_temp * Sin(theta)^3 * dTheta, 2 )
	fyy2_temp = Total( fyy2_temp * Sin(theta)^2 * dTheta, 2 )
	fyy3_temp = Total( fyy3_temp * Sin(theta)   * dTheta, 2 )
	
	;Pyz
	fyz1_temp = Total( fyz1_temp * Sin(theta)^2 * Cos(theta) * dTheta, 2 )
	fyz2_temp = Total( fyz2_temp * Sin(theta)^2              * dTheta, 2 )
	fyz3_temp = Total( fyz3_temp * Sin(theta)   * Cos(theta) * dTheta, 2 )
	fyz4_temp = Total( fyz4_temp * Sin(theta)                * dTheta, 2 )
	
	;Pzz
	fzz1_temp = Total( fzz1_temp * Sin(theta) * Cos(theta)^2 * dTheta, 2 )
	fzz2_temp = Total( fzz2_temp * Sin(theta) * Cos(theta)   * dTheta, 2 )
	fzz3_temp = Total( fzz3_temp * Sin(theta)                * dTheta, 2 )

;-----------------------------------------------------
; Integrate over Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract phi and the bin sizes
	oEnergy = self.oDist['DEPEND_3']
	energy  = oEnergy['DATA']
	dEnergy = (oEnergy['DELTA_PLUS_VAR'] +  oEnergy['DELTA_MINUS_VAR'])['DATA']
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	;Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * Temporary(energy) / self.mass)
	dv   = eV2J * Temporary(dEnergy) / (self.mass * v)
	
	;SPACECRAFT POTENTIAL
	IF N_Elements(self.oVsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q / self.mass * self.oVsc['DATA'] )
		VM    = v^2 * (1 + signQ * (vsc^2 / v^2))
	ENDIF ELSE BEGIN
		VM = v^2
	ENDELSE
	
	
	;
	; VELOCITY
	;

	;V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
	
	;Pxx
	fxx1_temp = self.mass * Total( fxx1_temp * v * VM^(3.0/2.0) * dv, 2 ) * 1e12
	fxx2_temp = self.mass * Total( fxx2_temp * v * VM           * dv, 2 ) * 1e12
	fxx3_temp = self.mass * Total( fxx3_temp * v * Sqrt(VM)     * dv, 2 ) * 1e12
	Pxx       = Temporary(fxx1_temp) - Temporary(fxx2_temp) + Temporary(fxx3_temp)
	
	;Pxy
	fxy1_temp = self.mass * Total( fxy1_temp * v * VM^(3.0/2.0) * dv, 2 ) * 1e12
	fxy2_temp = self.mass * Total( fxy2_temp * v * VM           * dv, 2 ) * 1e12
	fxy3_temp = self.mass * Total( fxy3_temp * v * VM           * dv, 2 ) * 1e12
	fxy4_temp = self.mass * Total( fxy4_temp * v * Sqrt(VM)     * dv, 2 ) * 1e12
	Pxy       = Temporary(fxy1_temp) - Temporary(fxy2_temp) - Temporary(fxy3_temp) + Temporary(fxy4_temp)
	
	;Pxz
	fxz1_temp = self.mass * Total( fxz1_temp * v * VM^(3.0/2.0) * dv, 2 ) * 1e12
	fxz2_temp = self.mass * Total( fxz2_temp * v * VM           * dv, 2 ) * 1e12
	fxz3_temp = self.mass * Total( fxz3_temp * v * VM           * dv, 2 ) * 1e12
	fxz4_temp = self.mass * Total( fxz4_temp * v * Sqrt(VM)     * dv, 2 ) * 1e12
	Pxz       = Temporary(fxz1_temp) - Temporary(fxz2_temp) - Temporary(fxz3_temp) + Temporary(fxz4_temp)
	
	;Pyy
	fyy1_temp = self.mass * Total( fyy1_temp * v * VM^(3.0/2.0) * dv, 2 ) * 1e12
	fyy2_temp = self.mass * Total( fyy2_temp * v * VM           * dv, 2 ) * 1e12
	fyy3_temp = self.mass * Total( fyy3_temp * v * Sqrt(VM)     * dv, 2 ) * 1e12
	Pyy       = Temporary(fyy1_temp) - Temporary(fyy2_temp) + Temporary(fyy3_temp)
	
	;Pyz
	fyz1_temp = self.mass * Total( fyz1_temp * v * VM^(3.0/2.0) * dv, 2 ) * 1e12
	fyz2_temp = self.mass * Total( fyz2_temp * v * VM           * dv, 2 ) * 1e12
	fyz3_temp = self.mass * Total( fyz3_temp * v * VM           * dv, 2 ) * 1e12
	fyz4_temp = self.mass * Total( fyz4_temp * v * Sqrt(VM)     * dv, 2 ) * 1e12
	Pyz       = Temporary(fyz1_temp) - Temporary(fyz2_temp) - Temporary(fyz3_temp) + Temporary(fyz4_temp)
	
	;Pzz
	fzz1_temp = self.mass * Total( fzz1_temp * v * VM^(3.0/2.0) * dv, 2 ) * 1e12
	fzz2_temp = self.mass * Total( fzz2_temp * v * VM           * dv, 2 ) * 1e12
	fzz3_temp = self.mass * Total( fzz3_temp * v * Sqrt(VM)     * dv, 2 ) * 1e12
	Pzz       = Temporary(fzz1_temp) - Temporary(fzz2_temp) + Temporary(fzz3_temp)
	
	v  = !Null
	dv = !Null
	VM = !Null
;-----------------------------------------------------
; Create Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Energy-time spectrogram
	;   - Convert Pa to nPa
	oP = MrMatrixTS( self.oDist['TIMEVAR'], [ [[Pxx], [Pxy], [Pxz]], $
	                                          [[Pxy], [Pyy], [Pyz]], $
	                                          [[Pxz], [Pyz], [Pzz]] ] * 1e9, $
	                 CACHE = tf_cache, $
	                 NAME  = name, $
	                 /NO_COPY )
	
	;Attributes
	oP['CATDESC']       = 'Pressure computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oP['LABEL']         = 'P'
	oP['LABEL_PTR_1']   = ['x', 'y', 'z']
	oP['LABEL_PTR_2']   = ['x', 'y', 'z']
	oP['UNITS']         = 'nPa'
	oP['PLOT_TITLE']    = 'Pressure Tensor'
	oP['TITLE']         = 'P!C(nPa)'
	oP['SI_CONVERSION'] = '1e-9>Pa'
	
	RETURN, oP
END


;+
;   Rebin the distribution. This is useful, for example, if the angular bins have
;   been rotated into a new coordinate system. Distinct velocity-space bins from
;   the old distribution may fall into the same velocity-space bin when rebinned
;   after rotation. In this case, data is weighted by the volume of the old
;   velocity-space bin and averaged into the new bin.
;
; :Params:
;       DV:             in, required, type=TxNxMxL fltarr
;                       Volume of each velocity space element of the distribution function
;                           before re-binning.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;       NPHI_BINS:      in, optional, type=integer, default=same as implicit f
;                       Number of azimuthal bins in the output distribution.
;       NTHETA_BINS:    in, optional, type=integer, default=same as implicit f
;                       Number of polar bins in the output distribution.
;       PHI_RANGE:      in, optional, type=fltarr(2), default=[-180.0, 180.0]
;                       Azimuthal range of the output distribution function.
;       THETA_RANGE:    in, optional, type=fltarr(2), default=[0.0, 180.0]
;                       Polar range of the output distribution function.
;
; :Returns:
;       ODIST:          out, required, type=TxNxMxL fltarr
;                       The re-binned distribution function.
;
;-
FUNCTION MrDist4D::RebinAngles, odV, $
CACHE=cache, $
NAME=name, $
NPHI_BINS=nPhi_Bins, $
NTHETA_BINS=nTheta_Bins, $
PHI_RANGE=phi_range, $
THETA_RANGE=theta_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D)    GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oPhiBins)   GT 0 THEN Obj_Destroy, oPhiBins
		IF N_Elements(oThetaBins) GT 0 THEN Obj_Destroy, oThetaBins
		IF N_Elements(oDist)      GT 0 THEN Obj_Destroy, oDist
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_rebinned'

	;Allocate memory
	dims    = Size(self.oDist, /DIMENSIONS)
	nTime   = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	f       = FltArr( nTime, nPhi, nTheta, nEnergy )

	;Step over each time
	FOR i = 0, nTime - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		f[i,*,*,*] = oDist3D -> RebinAngles( Reform(odV[i,*,*,*]), phi, dPhi, theta, dTheta, $
		                                     NPHI_BINS   = nPhi_Bins, $
		                                     NTHETA_BINS = nTheta_Bins, $
		                                     PHI_RANGE   = phi_range, $
		                                     THETA_RANGE = theta_range )

		;Destroy the object
		Obj_Destroy, oDist3D
	ENDFOR

;-------------------------------------------
; Datasets /////////////////////////////////
;-------------------------------------------
	;Time variable
	oTime = self.oDist['TIMEVAR']
	
	;Energy-time spectrogram
	oDist = MrTimeSeries( oTime, f, $
	                      CACHE = tf_cache, $
	                      NAME  = name, $
	                      /NO_COPY )
	
	;Phi
	binName  = name + '_PhiBins'
	oPhiBins = Size(phi, /N_DIMENSIONS) EQ 2 $
	                 ? MrTimeSeries( oTime, phi, NAME=binName, /NO_COPY ) $
	                 : MrVariable( phi, NAME=binName, /NO_COPY )
	
	;Theta
	binName    = name + '_ThetaBins'
	oThetaBins = Size(theta, /N_DIMENSIONS) EQ 2 $
	                 ? MrTimeSeries( oTime, theta, NAME=binName, /NO_COPY ) $
	                 : MrVariable( theta, NAME=binName, /NO_COPY )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	;Phi attributes
	oPhiBins['DELTA_MINUS'] = dPhi/2.0
	oPhiBins['DELTA_PLUS']  = dPhi/2.0
	oPhiBins['UNITS']       = 'degrees'
	oPhiBins['TITLE']       = 'Azimuth'
	oPhiBins['PLOT_TITLE']  = 'Azimuthal Bin Centers'
	
	;Theta attributes
	oThetaBins['DELTA_MINUS'] = dTheta/2.0
	oThetaBins['DELTA_PLUS']  = dTheta/2.0
	oThetaBins['UNITS']      = 'degrees'
	oThetaBins['TITLE']      = 'Polar Angle'
	oThetaBins['PLOT_TITLE'] = 'Polar Bin Centers'
	
	;Distribution attributes
	self.oDist       -> CopyAttrTo, oDist
	oDist['DEPEND_1'] = oPhiBins
	oDist['DEPEND_2'] = oThetaBins
	oDist['SCALE']    = 1B
	oDist['LOG']      = 1B
	
	RETURN, oDist
END


;+
;   Transform a spherical coordinate grid into a cartesian coordinate grid.
;
; :Keywords:
;       ELEVATION:          in, optional, type=boolean
;                           If set, THETA is taken to be the elevation angle.
;       MASS:               in, optional, type=float
;                           Mass (kg) of the species represented in the distribution.
;       SPECIES:            in, optional, type=string
;                           Species of particle represented in the distribution. Options
;                               are: { 'e' | 'p' | 'H' | 'He' | 'O' }
;       _REF_EXTRA:         in, optional, type=any
;                           Any keyword accepted by MrTimeSeries::SetProperty
;-
PRO MrDist4D::SetProperty, $
ELEVATION=elevation, $
MASS=mass, $
SPECIES=species
	Compile_Opt idl2
	On_Error, 2
	
	IF N_Elements(elevation) GT 0 THEN self.elevation = Keyword_Set(elevation)
	
	IF N_Elements(mass) GT 0 THEN BEGIN
		N = round(mass / MrConstants('m_p'))
		CASE N of
			0:    species = 'e'
			1:    species = 'H'
			2:    species = 'He'
			16:   species = 'O'
			ELSE: Message, 'Unable to determine particle species given MASS.'
		ENDCASE
		self.mass    = mass
		self.species = species
	ENDIF
	
	IF N_Elements(species) GT 0 THEN BEGIN
		IF ~MrIsMember(['e', 'i', 'H', 'He', 'O'], species) $
			THEN Message, 'SPECIES must be {"e" | "H" | "He" | "O"}'
		self.mass    = MrConstants('m_' + species)
		self.species = species
	ENDIF
END


;+
;   Set the array.
;
;   CALLING SEQUENCE:
;       oDist -> SetData, data
;       oDist -> SetData, time, data
;       oDist -> SetData, time, data, phi, theta, energy
;
; :Keywords:
;       TIME:           in, optional, type=NxM array
;                       Name or reference of a MrTimeVar object, or an array
;                           of time stamps. IF a name is provided, the assiciated
;                           variable must exist in the variable cache.
;       DATA:           in, required, type=NxM array
;                       Name or reference of a MrVariable object, or the dependent
;                           variable data. IF a name is given, the associated variable
;                           must exist in the variable cache.
;       PHI:            in, optional, type=Nx1 or NxM array
;                       Azimuthal coordinates of the distribution pixels. Can be the name
;                           or reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           It must have dimensions of [phi, theta] or [time, phi, theta]
;       THETA:          in, optional, type=Nx1 or NxM array
;                       Polar coordinates of the distribution pixels. Can be the name or
;                           reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           It must have dimensions of [phi, theta] or [time, phi, theta]
;       ENERGY:         in, optional, type=Nx1 or NxM array
;                       Energy coordinates of the distribution pixels. Can be the name
;                           or reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           IF data has two dimensions, one must be time and the other
;                           must be the same Size as the fourth dimension of the
;                           distribution.
;
; :Keywords:
;       DIMENSION:      in, optional, type=integer
;                       The time-dependent dimension of `DATA` (1-based). IF not
;                           provided, the dimension of `DATA` that is equal in Size to
;                           `TIME` is chosen as the default. IF zero or multiple
;                           dimensions match in this way, an error will occur.
;       T_TYPE:         in, optional, type=integer
;                       IF `TIME` is an array of time stamps, use this keyword to indicate
;                           the format or time-basis. See MrTimeVar FOR more details.
;       T_NAME:         in, optional, type=integer
;                       Name to be given to the MrTimeVar object. Ignored unless `TIME`
;                           is an array of time stamps.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       IF set `DATA` will be copied directly into the object
;                           and will be left undefined (a MrTimeSeries object will not
;                           be destroyed, but its array will be empty).
;-
PRO MrDist4D::SetData, time, data, phi, theta, energy, $
DEGREES=degrees, $
DIMENSION=dimension, $
NO_COPY=no_copy, $
T_NAME=t_name, $
T_TYPE=t_type, $
RADIANS=radians, $
UNITS=units
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN
	ENDIF
	
	;dist -> SetData, data
	IF N_Elements(time) GT 0 && N_Elements(DATA) EQ 0 THEN BEGIN
		data    = MrVar_Get(time)
		theTime = data['DEPEND_0']
		phi     = data['DEPEND_1']
		theta   = data['DEPEND_2']
		energy  = data['DEPEND_3']
		
	;dist -> SetData, time, data, ...
	ENDIF ELSE BEGIN
		theTime = time
	ENDELSE
	
;-----------------------------------------------------
; Check Distribution \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Use the superclass
	self.oDist -> SetData, theTime, data, $
	                       DIMENSION = dimension, $
	                       NO_COPY   = no_copy, $
	                       T_TYPE    = t_type, $
	                       T_NAME    = t_name

	;Check the results
	;   - Make sure it is NxMxLxK (4 dimensions of any Size)
	sz = Size(self.oDist)
	IF sz[0] NE 4 THEN Message, 'Invalid dimensions: Data must be an 4D.'
	
	;Units
	IF N_Elements(units) EQ 0 THEN BEGIN
		IF self.oDist -> HasAttr('UNITS') THEN BEGIN
			units = self.oDist['UNITS']
			IF StRegEx(units, '(psd|phase|space|density|s\^3[ ]*/[ ]*cm\^6)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
				units = 'PSD'
			ENDIF ELSE IF StRegEx(units, '(eflux|energy flux)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
				units = 'EFLUX'
			ENDIF ELSE IF StRegEx(units, '(diff flux)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
				units = 'DIFF FLUX'
			ENDIF ELSE BEGIN
				MrPrintF, 'LogWarn', 'Units "' + units + '" interpreted as "EFLUX".'
				units = 'EFLUX'
			ENDELSE
		ENDIF ELSE BEGIN
			MrPrintF, 'LogWarn', 'No units given. Assuming PSD.'
			units = 'PSD'
		ENDELSE
	ENDIF
	
;-----------------------------------------------------
; Set Dependents \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Set remaining data
	self.units = units
	IF N_Elements(theta)  GT 0 THEN self -> SetTheta, theta, DEGREES=degrees, RADIANS=radians
	IF N_Elements(phi)    GT 0 THEN self -> SetPhi, phi, DEGREES=degrees, RADIANS=radians
	IF N_Elements(energy) GT 0 THEN self -> SetEnergy, energy
END


;+
;   Set the spacecraft potential.
;
; :Params:
;       VSC:            in, required, type=string/integer/objref
;                       Name, number, or MrScalarTS object for the spacecraft potential.
;                           If VSC is a MrScalarTS object, its time property
;                           must be the same as that of the implicit distribution.
;                           If data has two dimensions, one must be time and the other
;                           must be the same Size as the fourth dimension of the
;                           distribution.
;
; :Keywords:
;       NO_COPY:        in, optional, type=boolean, default=0
;                       IF set `PHI` will be copied directly into the object
;                           and will be left undefined (a MrTimeSeries object will not
;                           be destroyed, but its array will be empty).
;-
PRO MrDist4D::SetVsc, Vsc, $
NO_COPY=no_copy
	Compile_Opt idl2
	On_Error, 2
	
	;Dimensions
	theDims = Size(self.oDist, /DIMENSIONS)
	nTimes  = theDims[0]

;-----------------------------------------------------
; Obtain Variable Object \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;DATA
	IF MrIsA(Vsc, /NUMBER, /ARRAY) THEN BEGIN
		;Check size & create MrScalarTS variable
		IF N_Elements(Vsc) EQ nTimes $
			THEN oVsc = MrScalarTS( self.oDist['TIMEVAR'], Vsc, NAME=self.name + '_Vsc'  ) $
			ELSE Message, 'VSC must have the same number of samples as the distribution.'
	
	;VARIABLE
	ENDIF ELSE BEGIN
		oVsc = MrVar_Get(Vsc)
	
		;Compare With Distribution
		IF ~oVsc -> IsTimeIdentical( self.oDist ) $
			THEN oVsc = oVsc -> Interpol(self.oDist)
	ENDELSE

;-----------------------------------------------------
; Set Attribute \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	Obj_Destroy, self.oVsc
	self.oVsc = oVsc
END


;+
;   Set the phi array.
;
; :Params:
;       PHI:            in, required, type=Nx1 or NxM array
;                       Azimuthal coordinates of the distribution pixels. Can be the name
;                           or reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           It must have dimensions of [phi, theta] or [time, phi, theta]
;
; :Keywords:
;       DEGREES:        in, optional, type=boolean, default=1
;                       IF set to zero, `THETA` has units of radians. Degrees are
;                           assumed. Cannot be used with `RADIANS`.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       IF set `THETA` will be copied directly into the object
;                           and will be left undefined (a MrTimeSeries object will not
;                           be destroyed, but its array will be empty).
;       RADIANS:        in, optional, type=boolean, default=1
;                       IF set, `THETA` has units of radians. Otherwise, degrees are
;                           assumed. Cannot be used with `DEGREES`.
;-
PRO MrDist4D::SetEnergy, energy, $
NO_COPY=no_copy
	Compile_Opt idl2
	On_Error, 2
	
	;Dimensions
	oTime   = self.oDist['TIMEVAR']
	theDims = Size(self.oDist, /DIMENSIONS)
	nTimes  = theDims[0]
	nEnergy = theDims[3]
	
	IF N_Elements(degrees) GT 0 && N_Elements(radians) GT 0 $
		THEN Message, 'DEGREES and RADIANS are mutually exclusive.'
	IF N_Elements(degrees) GT 0 THEN tf_degrees = Keyword_Set(degrees)
	IF N_Elements(radians) GT 0 THEN tf_degrees = ~Keyword_Set(radians)
	
	;
	; Steps:
	;   1. Obtain variable object
	;   2. Compare Size to distribution FUNCTION
	;   3. Set DEPEND_2 attribute
	;
	
;-----------------------------------------------------
; Obtain Variable Object \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Use existing energy
	IF N_Elements(energy) EQ 0 THEN BEGIN
		oEnergy = self.oDist['DEPEND_3']
	
	;Check given energy
	ENDIF ELSE BEGIN
		;DATA
		IF MrIsA(energy, /NUMBER, /ARRAY) THEN BEGIN
			dims  = Size(energy, /DIMENSIONS)
			nDims = N_Elements(dims)
			
			CASE nDims OF
				1: oEnergy = MrTimeSeries( oTime, Rebin(Reform(energy, 1, nEnergy), nTimes, nEnergy), $
				                        NAME = self.oDist.name + '_energy'  )
				2: oEnergy = MrTimeSeries( oTime, energy, NAME=self.name + '_energy' )
				ELSE: Message, 'energy must have 1 or 2 dimensions.'
			ENDCASE
		
		;VARIABLE
		ENDIF ELSE BEGIN
			oEnergy = MrVar_Get(energy)
			oEnergy = oEnergy -> Copy(self.oDist.name + '_energy')
			IF ~Obj_IsA(oEnergy, 'MrTimeSeries') THEN BEGIN
				oEnergy = MrTimeSeries( oTime, Rebin(Reform(oEnergy['DATA'], 1, nEnergy), nTimes, nEnergy), $
				                        NAME = self.oDist.name + '_energy' )
			ENDIF
		ENDELSE
	ENDELSE
	
;-----------------------------------------------------
; Delta Plus \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF oEnergy -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		oEnergy_dPlus = oEnergy['DELTA_PLUS_VAR']
		IF ~Obj_IsA(oEnergy_dPlus, 'MrTimeSeries') THEN BEGIN
			oEnergy_dPlus = MrTimeSeries(oTime, Rebin(Reform(oEnergy_dPlus['DATA'], 1, nEnergy), nTimes, nEnergy), $
			                          NAME = oEnergy.name + '_dPlus')
		ENDIF
	
	ENDIF ELSE IF oEnergy -> HasAttr('DELTA_PLUS') THEN BEGIN
		dPlus = oEnergy['DELTA_PLUS']
		oEnergy_dPlus = MrTimeSeries(oTime, Rebin([dPlus], nTimes, nEnergy), $
		                          NAME = oEnergy.name + '_dPlus')
		oEnergy -> RemoveAttr, 'DELTA_PLUS'
	
	;Half of the mean spacing between elements
	ENDIF ELSE BEGIN
		dPlus = self -> DeltaE(oEnergy['DATA']) / 2.0
		oEnergy_dPlus = MrTimeSeries(oTime, dPlus, $
		                             NAME = oEnergy.name + '_dPlus')
	ENDELSE
	
	;Add as attribute
	oEnergy['DELTA_PLUS_VAR'] = oEnergy_dPlus
	
;-----------------------------------------------------
; Delta Minus \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF oEnergy -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		oEnergy_dMinus = oEnergy['DELTA_MINUS_VAR']
		IF ~Obj_IsA(oEnergy_dMinus, 'MrTimeSeries') THEN BEGIN
			oEnergy_dMinus = MrTimeSeries(oTime, Rebin(Reform(oEnergy_dMinus['DATA'], 1, nEnergy), nTimes, nEnergy), $
			                          NAME = oEnergy.name + '_dMinus')
		ENDIF
	
	ENDIF ELSE IF oEnergy -> HasAttr('DELTA_PLUS') THEN BEGIN
		dMinus = oEnergy['DELTA_PLUS']
		oEnergy_dMinus = MrTimeSeries(oTime, Rebin([dMinus], nTimes, nEnergy), $
		                          NAME = oEnergy.name + '_dMinus')
		oEnergy -> RemoveAttr, 'DELTA_PLUS'
	
	;Half of the mean spacing between elements
	ENDIF ELSE BEGIN
		dMinus = self -> DeltaE(oEnergy['DATA']) / 2.0
		oEnergy_dMinus = MrTimeSeries(oTime, dMinus, $
		                              NAME = oEnergy.name + '_dMinus')
	ENDELSE
	
	;Add as attribute
	oEnergy['DELTA_MINUS_VAR'] = oEnergy_dMinus

;-----------------------------------------------------
; Compare With Distribution \\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Make sure DIST and energy use the same time object
	IF ~oEnergy -> IsTimeIdentical( self.oDist ) $
		THEN Message, 'energy and DIST have different time objects.'
	
;-----------------------------------------------------
; SetProperties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	self.oDist['DEPEND_3'] = oEnergy
END


;+
;   Set the phi array.
;
; :Params:
;       PHI:            in, required, type=Nx1 or NxM array
;                       Azimuthal coordinates of the distribution pixels. Can be the name
;                           or reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           It must have dimensions of [phi, theta] or [time, phi, theta]
;
; :Keywords:
;       DEGREES:        in, optional, type=boolean, default=1
;                       IF set to zero, `THETA` has units of radians. Degrees are
;                           assumed. Cannot be used with `RADIANS`.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       IF set `THETA` will be copied directly into the object
;                           and will be left undefined (a MrTimeSeries object will not
;                           be destroyed, but its array will be empty).
;       RADIANS:        in, optional, type=boolean, default=1
;                       IF set, `THETA` has units of radians. Otherwise, degrees are
;                           assumed. Cannot be used with `DEGREES`.
;-
PRO MrDist4D::SetPhi, phi, $
DEGREES=degrees, $
NO_COPY=no_copy, $
RADIANS=radians
	Compile_Opt idl2
	On_Error, 2
	
	;Dimensions
	oTime   = self.oDist['TIMEVAR']
	theDims = Size(self.oDist, /DIMENSIONS)
	nTimes  = theDims[0]
	nPhi    = theDims[1]
	
	IF N_Elements(degrees) GT 0 && N_Elements(radians) GT 0 $
		THEN Message, 'DEGREES and RADIANS are mutually exclusive.'
	
	;
	; Steps:
	;   1. Obtain variable object
	;   2. Compare Size to distribution FUNCTION
	;   3. Set DEPEND_2 attribute
	;
	
;-----------------------------------------------------
; Obtain Variable Object \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Use existing PHI
	IF N_Elements(phi) EQ 0 THEN BEGIN
		oPhi = self.oDist['DEPEND_1']
	
	;Check given PHI
	ENDIF ELSE BEGIN
		;DATA
		IF MrIsA(phi, /NUMBER, /ARRAY) THEN BEGIN
			dims  = Size(phi, /DIMENSIONS)
			nDims = N_Elements(dims)
			
			CASE nDims OF
				1: oPhi = MrTimeSeries( oTime, Rebin(Reform(phi, 1, nPhi), nTimes, nPhi), $
				                        NAME = self.oDist.name + '_phi'  )
				2: oPhi = MrTimeSeries( oTime, phi, NAME=self.name + '_phi' )
				ELSE: Message, 'PHI must have 1 or 2 dimensions.'
			ENDCASE
		
		;VARIABLE
		ENDIF ELSE BEGIN
			oTemp = MrVar_Get(phi)
			oPhi  = oTemp -> Copy(self.oDist.name + '_phi')
			IF ~Obj_IsA(oPhi, 'MrTimeSeries') THEN BEGIN
				oPhi = MrTimeSeries( oTime, Rebin(Reform(oPhi['DATA'], 1, nPhi), nTimes, nPhi), $
				                     NAME = self.oDist.name + '_phi' )
				oTemp -> CopyAttrTo, oPhi
			ENDIF
		ENDELSE
	ENDELSE
	
;-----------------------------------------------------
; Delta Plus \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF oPhi -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		oPhi_dPlus = oPhi['DELTA_PLUS_VAR']
		IF ~Obj_IsA(oPhi_dPlus, 'MrTimeSeries') THEN BEGIN
			oPhi_dPlus = MrTimeSeries(oTime, Rebin(Reform(oPhi_dPlus['DATA'], 1, nPhi), nTimes, nPhi), $
			                          NAME = oPhi.name + '_dPlus')
		ENDIF
	
	ENDIF ELSE IF oPhi -> HasAttr('DELTA_PLUS') THEN BEGIN
		dPlus = oPhi['DELTA_PLUS']
		oPhi_dPlus = MrTimeSeries(oTime, Rebin([dPlus], nTimes, nPhi), $
		                          NAME = oPhi.name + '_dPlus')
		oPhi -> RemoveAttr, 'DELTA_PLUS'
	
	;Half of the mean spacing between elements
	ENDIF ELSE BEGIN
		dPlus = Mean(oPhi['DATA',*,1:*] - oPhi['DATA'])
		oPhi_dPlus = MrTimeSeries(oTime, Rebin([dPlus], nTimes, nPhi), $
		                          NAME = oPhi.name + '_dPlus')
	ENDELSE
	
	;Add as attribute
	oPhi['DELTA_PLUS_VAR'] = oPhi_dPlus
	
;-----------------------------------------------------
; Delta Minus \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF oPhi -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		oPhi_dMinus = oPhi['DELTA_MINUS_VAR']
		IF ~Obj_IsA(oPhi_dMinus, 'MrTimeSeries') THEN BEGIN
			oPhi_dMinus = MrTimeSeries(oTime, Rebin(Reform(oPhi_dMinus['DATA'], 1, nPhi), nTimes, nPhi), $
			                          NAME = oPhi.name + '_dMinus')
		ENDIF
	
	ENDIF ELSE IF oPhi -> HasAttr('DELTA_PLUS') THEN BEGIN
		dMinus = oPhi['DELTA_PLUS']
		oPhi_dMinus = MrTimeSeries(oTime, Rebin([dMinus], nTimes, nPhi), $
		                          NAME = oPhi.name + '_dMinus')
		oPhi -> RemoveAttr, 'DELTA_PLUS'
	
	;Half of the mean spacing between elements
	ENDIF ELSE BEGIN
		dMinus = 0
		oPhi_dMinus = MrTimeSeries(oTime, Rebin([dMinus], nTimes, nPhi), $
		                           NAME = oPhi.name + '_dMinus')
	ENDELSE
	
	;Add as attribute
	oPhi['DELTA_MINUS_VAR'] = oPhi_dMinus

;-----------------------------------------------------
; Compare With Distribution \\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Make sure DIST and PHI use the same time object
	IF ~oPhi -> IsTimeIdentical( self.oDist ) $
		THEN Message, 'PHI and DIST have different time objects.'

;-----------------------------------------------------
; Units \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF N_Elements(degrees) GT 0 && N_Elements(radians) GT 0 THEN BEGIN
		Message, 'DEGREES and RADIANS are mutually exclusive.'
	ENDIF ELSE IF N_Elements(degrees) GT 0 THEN BEGIN
		tf_degrees = Keyword_Set(degrees)
	ENDIF ELSE IF N_Elements(radians) GT 0 THEN BEGIN
		tf_degrees = ~Keyword_Set(radians)
	ENDIF ELSE BEGIN
		IF oPhi -> HasAttr('UNITS') THEN BEGIN
			units = oPhi['UNITS']
			CASE 1 of
				stregex(units, 'deg', /BOOLEAN, /FOLD_CASE): tf_degrees = 1B
				stregex(units, 'rad', /BOOLEAN, /FOLD_CASE): tf_degrees = 0B
				ELSE: BEGIN
					MrPrintF, 'LogWarn', 'Units not recognized "' + units + '". Assuming degrees.'
					tf_degrees = 1B
				ENDELSE
			ENDCASE
		ENDIF ELSE BEGIN
			tf_degrees = N_Elements(radians) EQ 0 ? Keyword_Set(degrees) : ~Keyword_Set(radians)
		ENDELSE
	ENDELSE

	;Convert to degrees
	IF ~tf_degrees THEN BEGIN
		oPhi -> SetData, oPhi['DATA']*!radeg
		oPhi['DELTA_MINUS_VAR'] *= !radeg
		oPhi['DELTA_PLUS_VAR'] *= !radeg
		oPhi['UNITS'] = 'degrees'
	ENDIF

;-----------------------------------------------------
; SetProperties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	self.oDist['DEPEND_1'] = oPhi
END


;+
;   Set the phi array.
;
; :Params:
;       PHI:            in, required, type=Nx1 or NxM array
;                       Azimuthal coordinates of the distribution pixels. Can be the name
;                           or reference of a MrVariable object, or the variable data.
;                           IF the variable is a MrTimeSeries object, its time property
;                           must be the same as that of the implicit distribution.
;                           It must have dimensions of [phi, theta] or [time, phi, theta]
;
; :Keywords:
;       DEGREES:        in, optional, type=boolean, default=1
;                       IF set to zero, `THETA` has units of radians. Degrees are
;                           assumed. Cannot be used with `RADIANS`.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       IF set `THETA` will be copied directly into the object
;                           and will be left undefined (a MrTimeSeries object will not
;                           be destroyed, but its array will be empty).
;       RADIANS:        in, optional, type=boolean, default=1
;                       IF set, `THETA` has units of radians. Otherwise, degrees are
;                           assumed. Cannot be used with `DEGREES`.
;-
PRO MrDist4D::SetTheta, theta, $
DEGREES=degrees, $
NO_COPY=no_copy, $
RADIANS=radians
	Compile_Opt idl2
	On_Error, 2
	
	;Dimensions
	oTime   = self.oDist['TIMEVAR']
	theDims = Size(self.oDist, /DIMENSIONS)
	nTimes  = theDims[0]
	nTheta  = theDims[2]
	
	;
	; Steps:
	;   1. Obtain variable object
	;   2. Compare Size to distribution FUNCTION
	;   3. Set DEPEND_2 attribute
	;
	
;-----------------------------------------------------
; Obtain Variable Object \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Use existing theta
	IF N_Elements(theta) EQ 0 THEN BEGIN
		oTheta = self.oDist['DEPEND_2']
	
	;Check given theta
	ENDIF ELSE BEGIN
		;DATA
		IF MrIsA(theta, /NUMBER, /ARRAY) THEN BEGIN
			dims  = Size(theta, /DIMENSIONS)
			nDims = N_Elements(dims)
			
			CASE nDims OF
				1: oTheta = MrTimeSeries( oTime, Rebin(Reform(theta, 1, nTheta), nTimes, nTheta), $
				                          NAME = self.oDist.name + '_theta'  )
				2: oTheta = MrTimeSeries( oTime, theta, NAME=self.name + '_theta' )
				ELSE: Message, 'theta must have 1 or 2 dimensions.'
			ENDCASE
		
		;VARIABLE
		ENDIF ELSE BEGIN
			oTemp = MrVar_Get(theta)
			oTheta = oTemp -> Copy(self.oDist.name + '_theta')
			IF ~Obj_IsA(oTheta, 'MrTimeSeries') THEN BEGIN
				oTheta = MrTimeSeries( oTime, Rebin(Reform(oTheta['DATA'], 1, nTheta), nTimes, nTheta), $
				                       NAME = self.oDist.name + '_theta' )
				oTemp -> CopyAttrTo, oTheta
			ENDIF
		ENDELSE
	ENDELSE
	
;-----------------------------------------------------
; Delta Plus \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF oTheta -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		oTheta_dPlus = oTheta['DELTA_PLUS_VAR']
		IF ~Obj_IsA(oTheta_dPlus, 'MrTimeSeries') THEN BEGIN
			oTheta_dPlus = MrTimeSeries(oTime, Rebin(Reform(oTheta_dPlus['DATA'], 1, nTheta), nTimes, nTheta), $
			                            NAME = oTheta.name + '_dPlus')
		ENDIF
	
	ENDIF ELSE IF oTheta -> HasAttr('DELTA_PLUS') THEN BEGIN
		dPlus = oTheta['DELTA_PLUS']
		oTheta_dPlus = MrTimeSeries(oTime, Rebin([dPlus], nTimes, nTheta), $
		                            NAME = oTheta.name + '_dPlus')
		oTheta -> RemoveAttr, 'DELTA_PLUS'
	
	;Half of the mean spacing between elements
	ENDIF ELSE BEGIN
		dPlus = Mean(oTheta['DATA',*,1:*] - oTheta['DATA'])
		oTheta_dPlus = MrTimeSeries(oTime, Rebin([dPlus], nTimes, nTheta), $
		                            NAME = oTheta.name + '_dPlus')
	ENDELSE
	
	;Add as attribute
	oTheta['DELTA_PLUS_VAR'] = oTheta_dPlus
	
;-----------------------------------------------------
; Delta Minus \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF oTheta -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		oTheta_dMinus = oTheta['DELTA_MINUS_VAR']
		IF ~Obj_IsA(oTheta_dMinus, 'MrTimeSeries') THEN BEGIN
			oTheta_dMinus = MrTimeSeries(oTime, Rebin(Reform(oTheta_dMinus['DATA'], 1, nTheta), nTimes, nTheta), $
			                             NAME = oTheta.name + '_dMinus')
		ENDIF
	
	ENDIF ELSE IF oTheta -> HasAttr('DELTA_PLUS') THEN BEGIN
		dMinus = oTheta['DELTA_PLUS']
		oTheta_dMinus = MrTimeSeries(oTime, Rebin([dMinus], nTimes, nTheta), $
		                             NAME = oTheta.name + '_dMinus')
		oTheta -> RemoveAttr, 'DELTA_PLUS'
	
	;Half of the mean spacing between elements
	ENDIF ELSE BEGIN
		dMinus = 0
		oTheta_dMinus = MrTimeSeries(oTime, Rebin([dMinus], nTimes, nTheta), $
		                             NAME = oTheta.name + '_dMinus')
	ENDELSE
	
	;Add as attribute
	oTheta['DELTA_MINUS_VAR'] = oTheta_dMinus

;-----------------------------------------------------
; Compare With Distribution \\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Make sure DIST and theta use the same time object
	IF ~oTheta -> IsTimeIdentical( self.oDist ) $
		THEN Message, 'theta and DIST have different time objects.'

;-----------------------------------------------------
; Units \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF N_Elements(degrees) GT 0 && N_Elements(radians) GT 0 THEN BEGIN
		Message, 'DEGREES and RADIANS are mutually exclusive.'
	ENDIF ELSE IF N_Elements(degrees) GT 0 THEN BEGIN
		tf_degrees = Keyword_Set(degrees)
	ENDIF ELSE IF N_Elements(radians) GT 0 THEN BEGIN
		tf_degrees = ~Keyword_Set(radians)
	ENDIF ELSE BEGIN
		IF oTheta -> HasAttr('UNITS') THEN BEGIN
			units = oTheta['UNITS']
			CASE 1 of
				stregex(units, 'deg', /BOOLEAN, /FOLD_CASE): tf_degrees = 1B
				stregex(units, 'rad', /BOOLEAN, /FOLD_CASE): tf_degrees = 0B
				ELSE: BEGIN
					MrPrintF, 'LogWarn', 'Units not recognized "' + units + '". Assuming degrees.'
					tf_degrees = 1B
				ENDELSE
			ENDCASE
		ENDIF ELSE BEGIN
			tf_degrees = N_Elements(radians) EQ 0 ? Keyword_Set(degrees) : ~Keyword_Set(radians)
		ENDELSE
	ENDELSE
	
	;Convert to degrees
	IF ~tf_degrees THEN BEGIN
		oTheta -> SetData, oTheta['DATA']*!radeg
		oTheta['DELTA_MINUS_VAR'] *= !radeg
		oTheta['DELTA_PLUS_VAR'] *= !radeg
		oTheta['UNITS'] = 'degrees'
	ENDIF

;-----------------------------------------------------
; SetProperties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	self.oDist['DEPEND_2'] = oTheta
END


;+
;   Reduce the 3D distribution FUNCTION to a 2D distribution in polar angle and energy,
;   averaging over azimuth angle.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       IF set, the output variable will be added to the variable cache.
;       NAME:           in, optional, type=string, default=self.name + '_ThetaE'
;                       Name to be given to the output variable.
;       NE_BINS:        in, optional, type=integer
;                       Number of energy bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       NTHETA_BINS:    in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       PHI_RANGE:      in, optional, type=FltArr(2), default=[0.0\, 360.0]
;                       The range in azimuthal angle (degrees) over which to average.
;
; :Returns:
;       DIST2D:         out, required, type=MrTimeSeries object
;                       A time-varying 2D distribution in polar angle and energy.
;-
FUNCTION MrDist4D::ThetaE, $
CACHE=cache, $
NAME=name, $
NE_BINS=ne_bins, $
NTHETA_BINS=nTheta_bins, $
PHI_RANGE=phi_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D)    GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oThetaE)    GT 0 THEN Obj_Destroy, oThetaE
		IF N_Elements(oThetaBins) GT 0 THEN Obj_Destroy, oThetaBins
		IF N_Elements(oEBins)     GT 0 THEN Obj_Destroy, oEBins
		RETURN, !Null
	ENDIF
	
	;Defaults
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_ThetaE'

;-------------------------------------------
; Reduce the 4D Distribution ///////////////
;-------------------------------------------

	;Allocate memory
	dims    = Size(self.oDist, /DIMENSIONS)
	nTimes  = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	ThetaE  = FltArr( nTimes, nTheta, nEnergy )
	
	;Step over each time
	FOR i = 0, nTimes - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		ThetaE[i,*,*] = oDist3D -> ThetaE( theta, energy, dTheta, dE, $
		                                   NE_BINS     = ne_bins, $
		                                   NTHETA_BINS = nTheta_bins, $
		                                   PHI_RANGE   = phi_range )
		                                   
		;Destroy the 3D distribution
		Obj_Destroy, oDist3D
	ENDFOR

;-------------------------------------------
; Datasets /////////////////////////////////
;-------------------------------------------
	;Time variable
	oTime = self.oDist['TIMEVAR']
	
	;Theta-Energy distribution
	oThetaE = MrTimeSeries( oTime, ThetaE, $
	                        CACHE = cache, $
	                        NAME  = name, $
	                        /NO_COPY )
	
	;Theta
	binName    = name + '_ThetaBins'
	oThetaBins = Size(theta, /N_DIMENSIONS) EQ 2 $
	                 ? MrTimeSeries( oTime, theta, NAME=binName, /NO_COPY ) $
	                 : MrVariable( theta, NAME=binName, /NO_COPY )
	
	;Energy
	binName     = name + '_EnergyBins'
	oEnergyBins = Size(energy, /N_DIMENSIONS) EQ 2 $
	                  ? MrTimeSeries( oTime, energy, NAME=binName, /NO_COPY ) $
	                  : MrVariable( energy, NAME=binName, /NO_COPY )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	
	;Theta attributes
	oThetaBins['DELTA_MINUS'] = dTheta
	oThetaBins['DELTA_PLUS']  = dTheta
	oThetaBins['UNITS']      = 'degrees'
	oThetaBins['TITLE']      = 'Polar Angle'
	oThetaBins['PLOT_TITLE'] = 'Polar Bin Centers'
	
	;Energy attributes
;	oEBins['UNITS'] = self.oEnergy['UNITS']
;	oEBins['TITLE'] = 'Energy'

	;Energy bins have not changed
	;   - MUST ALSO UPDATE MRDIST3D::SPECE
	oEBins = self.oDist['DEPEND_3']

	;Distribution attributes
	oThetaE['DEPEND_1'] = oThetaBins
	oThetaE['DEPEND_2'] = oEBins
	oThetaE['SCALE']    = 1B
	oThetaE['LOG']      = 1B
	oThetaE['UNITS']    = self.oDist['UNITS']
	
	;RETURN the 2D distribution
	RETURN, oThetaE
END


;+
;   Reduce the 3D distribution FUNCTION to a 1D distribution in polar angle.
;
; :Keywords:
;       E_RANGE:        in, optional, type=FltArr(2), default=[min, max]
;                       The range in energy, in electron volts (eV) over which to average.
;       NTHETA_BINS:    in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       PHI_RANGE:      in, optional, type=FltArr(2), default=[0.0\, 360.0]
;                       The range in azimuthal angle (degrees) over which to average.
;
; :Returns:
;       OTHETASPEC:     out, required, type=MrTimeSeries
;                       A 1D distribution in time, averaged over energy and azimuth.
;-
FUNCTION MrDist4D::ThetaSpec, $
CACHE=cache, $
E_RANGE=e_range, $
NAME=name, $
NTHETA_BINS=nTheta_bins, $
PHI_RANGE=phi_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D)    GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(oThetaSpec) GT 0 THEN Obj_Destroy, oThetaSpec
		IF N_Elements(oThetaBins) GT 0 THEN Obj_Destroy, oThetaBins
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_ThetaSpec'

	;Allocate memory
	dims      = Size(self.oDist, /DIMENSIONS)
	nTime     = dims[0]
	nPhi      = dims[1]
	nTheta    = dims[2]
	nEnergy   = dims[3]
	ThetaSpec = FltArr( nTime, nTheta )
	
	;Step over each time
	FOR i = 0, nTime - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		ThetaSpec[i,*] = oDist3D -> ThetaSpec(theta_bins, dTheta, $
		                                      E_RANGE     = e_range, $
		                                      NTHETA_BINS = nTheta_bins, $
		                                      PHI_RANGE   = phi_range )

		;Destroy the object
		Obj_Destroy, oDist3D
	ENDFOR
	
	;Time variable
	oTime = self.oDist['TIMEVAR']
	
	;Theta-time spectrogram
	oThetaSpec = MrTimeSeries( oTime, thetaSpec, $
	                           CACHE = tf_cache, $
	                           NAME  = name, $
	                           /NO_COPY )
	
	;Abscissa
	binName    = name + '_ThetaBins'
	oThetaBins = Size(theta_bins, /N_DIMENSIONS) EQ 2 $
	                 ? MrTimeSeries( oTime, theta_bins, NAME=binName, /NO_COPY ) $
	                 : MrVariable( theta_bins, NAME=binName, /NO_COPY )
	
	;Theta attributes
	oThetaBins['DELTA_MINUS'] = dTheta
	oThetaBins['DELTA_PLUS']  = dTheta
	oThetaBins['UNITS']      = 'degrees'
	oThetaBins['TITLE']      = 'Polar Angle'
	oThetaBins['PLOT_TITLE'] = 'Polar Bin Centers'

	;Sepctrogram attributes
	oThetaSpec['DEPEND_1']   = oThetaBins
	oThetaSpec['SCALE']      = 1B
	oThetaSpec['LOG']        = 1B
	oThetaSpec['UNITS']      = self['UNITS']
	oThetaSpec['TITLE']      = 'Theta Dist'
	oThetaSpec['PLOT_TITLE'] = 'Distribution in Theta'
	
	RETURN, oThetaSpec
END


;+
;   Compute the temperature from the second moment of the distribution (pressure),
;   using the ideal gas equation of state.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       OP:             out, required, type=MrTimeSeries
;                       Pressure tensor as a function of time.
;-
FUNCTION MrDist4D::Temperature, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_Temperature'

	;Conversion from Kelvin to eV
	eV2K = 11600.0
	
	;Compute the pressure
	oN = self -> Density()
	oP = self -> Pressure()
	
	;Apply the equation of state
	;   - PV = NkT
	;   - T = kP/n
	;   - N = 1/cm^3
	;   - P = nPa
	;   - 1e-15 converts to Kelvin
	oT = oP / ( 1e15 * MrConstants('k_B') * oN )
	oT /= eV2K
	Obj_Destroy, [oN, oP]
	
	oT -> SetName, name
	IF tf_cache THEN oT -> Cache
	
	;Attributes
	oT['CATDESC']       = 'Temperature computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oT['LABEL']         = 'T'
	oT['LABEL_PTR_1']   = ['x', 'y', 'z']
	oT['LABEL_PTR_2']   = ['x', 'y', 'z']
	oT['UNITS']         = 'eV'
	oT['PLOT_TITLE']    = 'Temperature Tensor'
	oT['TITLE']         = 'T!C(eV)'
	oT['SI_CONVERSION'] = '>'
	
	RETURN, oT
END


;+
;   Reduce the 3D distribution FUNCTION to a 1D distribution in azimuth angle.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       ON:             out, required, type=MrScalarTS
;                       Density as a function of time.
;-
FUNCTION MrDist4D::Velocity, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_Velocity'
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'

	;Convert to radians
	deg2rad = !dpi / 180D
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	dims    = Size(self.oDist, /DIMENSIONS)
	
;-----------------------------------------------------
; Integrate over Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract phi and the bin sizes
	oPhi = self.oDist['DEPEND_1']
	phi  = oPhi['DATA'] * deg2rad
	dPhi = (oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR'])['DATA'] * deg2rad
	
	;Integrate phi
	;   - Integral( V * fdist * d^3V )
	;   - V    = (Vx, Vy, Vz)
	;   - Vx   = V * Sin(Theta) * Cos(Phi)
	;   - Vy   = V * Sin(Theta) * Sin(Phi)
	;   - Vz   = V * Cos(Theta)
	;   - (d^3)V = V^2 * Sin(Theta) * dV * dTheta * dPhi
	fx_temp = Total(self.oDist['DATA'] * Rebin(Cos(phi) * dPhi, dims), 2)
	fy_temp = Total(self.oDist['DATA'] * Rebin(Sin(phi) * dPhi, dims), 2)
	fz_temp = Total(self.oDist['DATA'] * Rebin(dPhi, dims), 2)
	
	phi  = !Null
	dPhi = !Null
;-----------------------------------------------------
; Integrate over Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract theta and the bin sizes
	oTheta = self.oDist['DEPEND_2']
	theta  = oTheta['DATA'] * deg2rad
	dTheta = (oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR'])['DATA'] * deg2rad
	
	;Integrate
	fx_temp = Total( Temporary(fx_temp) * Rebin(Sin(theta)^2 * dTheta, dims[[0,2,3]]), 2)
	fy_temp = Total( Temporary(fy_temp) * Rebin(Sin(theta)^2 * dTheta, dims[[0,2,3]]), 2)
	fz_temp = Total( Temporary(fz_temp) * Rebin(Sin(theta)*Cos(theta) * dTheta, dims[[0,2,3]]), 2)
	
	theta  = !Null
	dTheta = !Null
;-----------------------------------------------------
; Integrate over Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract phi and the bin sizes
	oEnergy = self.oDist['DEPEND_3']
	energy  = oEnergy['DATA']
	dEnergy = (oEnergy['DELTA_PLUS_VAR'] + oEnergy['DELTA_MINUS_VAR'])['DATA']
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	;Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * Temporary(energy) / self.mass)
	dv   = eV2J * Temporary(dEnergy) / (self.mass * v)
	
	;Spacecraft potential correction
	IF N_Elements(self.oVsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q * (self.oVsc['DATA']) / self.mass )
		
		Vx = Total( Temporary(fx_temp) * v * (v^2 + sgn*vsc^2) * dv, 2 ) * 1e12
		Vy = Total( Temporary(fy_temp) * v * (v^2 + sgn*vsc^2) * dv, 2 ) * 1e12
		Vz = Total( Temporary(fz_temp) * v * (v^2 + sgn*vsc^2) * dv, 2 ) * 1e12
	
	ENDIF ELSE BEGIN
		;V is in m/s while f is in s^3/cm^6
		;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
		Vx = Total( Temporary(fx_temp) * v^3 * dv, 2 ) * 1e12
		Vy = Total( Temporary(fy_temp) * v^3 * dv, 2 ) * 1e12
		Vz = Total( Temporary(fz_temp) * v^3 * dv, 2 ) * 1e12
	ENDELSE
	
	v  = !Null
	dv = !Null
;-----------------------------------------------------
; Create Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Velocity
	;   - Convert density 1/cm^3 -> 1/m^3
	;   - Convert velocity m/s -> km/s
	;   - Must negate to convert look angles to angles of incidence
	oN = 1e6 * self -> Density()
	oV = MrVectorTS( self.oDist['TIMEVAR'], -[[Temporary(Vx)], [Temporary(Vy)], [Temporary(Vz)]] )
	oV = 1e-3 * oV / oN
	
	;Attributes
	IF tf_cache THEN oV -> Cache
	oV -> SetName, name
	oV['CATDESC']       = 'Bulk velocity computed from the 3D velocity space integral ' + $
	                      'of the distribution function.'
	oV['LABEL']         = ['Vx', 'Vy', 'Vz']
	oV['UNITS']         = 'km/s'
	oV['TITLE']         = 'V!C(km/s)'
	oV['SI_CONVERSION'] = '1e3>m/s'
	
	RETURN, oV
END


;+
;   Compute the size of velocity-space volume elements.
;
;       dV = v^2 * Sin(theta) * dv * dTheta * dPhi
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       DV:             out, required, type=MrScalarTS
;                       Size of each velocity-space volume element.
;-
FUNCTION MrDist4D::VolumeElement, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oDist3D) GT 0 THEN Obj_Destroy, oDist3D
		IF N_Elements(odV)      GT 0 THEN Obj_Destroy, odV
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_dV'

	;Allocate memory
	dims    = Size(self.oDist, /DIMENSIONS)
	nTime   = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	dV    = FltArr( nTime, nPhi, nTheta, nEnergy )

	;Step over each time
	FOR i = 0, nTime - 1 DO BEGIN
		oDist3D = self -> GetDist3D(i)
		
		;Reduce the distribution
		dV[i,*,*,*] = oDist3D -> VolumeElement()

		;Destroy the object
		Obj_Destroy, oDist3D
	ENDFOR
	
	;Energy-time spectrogram
	odV = MrTimeSeries( self.oDist['TIMEVAR'], dV, $
	                    CACHE = tf_cache, $
	                    NAME  = name, $
	                    /NO_COPY )
	
	;Attributes
	self.oDist          -> CopyAttrTo, odV, ['DEPEND_1', 'DEPEND_2', 'DEPEND_3']
	odV['CATDESC']       = 'Size of each velocity space volume element.'
	odV['UNITS']         = 'sr m^3/s^3'
	odV['TITLE']         = 'dV!C(sr m^3/s^3)'
	odV['SI_CONVERSION'] = '1e0>sr m^3/s^3'
	
	RETURN, odV
END


;+
;   The class definition statement.
;
; :Params:
;       CLASS:          out, optional, type=structure
;-
PRO MrDist4D__DEFINE
	Compile_Opt idl2
	
	class = { MrDist4D, $
	          elevation: 0B, $
	          mass:      0.0, $
	          oVsc:      Obj_New(), $
	          oDist:     Obj_New(), $
	          species:   '', $
	          units:     '' $
	        }
END