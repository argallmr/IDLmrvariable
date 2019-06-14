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
;       2018/02/12  -   Written by Matthew Argall
;       2018/10/24  -   Take care when velocity bins are below spacecraft potential. - MRA
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


PRO MrDist4D::Add, oDist
	Compile_Opt idl2
	On_Error, 2

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
;   factors are (start with column):
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
	;                           = [ kg^2 / ( eV * (1.602e-19 J/eV) )^2 ]        / ( cm^2 * s * sr ) * 0.5
	;                           = [ kg^2 / ( kg m^2 / s^2 )^2 ]                 / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2)
	;                           = [ (N * 1.672e-27 kg)^2 / ( kg^2 m^4 / s^4 ) ] / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2)
	;                           = [ kg^2 / ( kg^2 m^4 / s^4 ) ]                 / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2) * N^2 * 1.672e-27^2
	;                           = [ 1 / (m^4 * 1e8 * cm^4/m^4 / s^4) ]          / ( cm^2 * s * sr ) * 0.5 * 1.602e-19^(-2) * N^2 * 1.672e-27^2
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
GROUND=ground, $
NAME=name
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache  = Keyword_Set(cache)
	tf_ground = Keyword_Set(ground)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_Density'
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Constants
	deg2rad = !dpi / 180D
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	dims    = Size(self.oDist, /DIMENSIONS)
	tf_nan  = 0B
	
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
		IF tf_ground THEN BEGIN
			vsc = Rebin( Sqrt( 2.0 * q * (self.oVsc['DATA']) / self.mass ), dims[[0,3]] )
			v  -= vsc
			
			iBad = Where(v LT 0, nBad)
			IF nBad GT 0 THEN BEGIN
				tf_nan  = 1B
				v[iBad] = 0;!Values.D_NaN
			ENDIF
				
			VM  = v^2
		
		ENDIF ELSE BEGIN
			;Sign is charge dependent. E = q * V = 1/2 * m * v^2
			;   - For the case of electrons in a positive spacecraft potential, a further
			;     straightforward effect of the potential is to require, in the discrete
			;     summation, the exclusion of velocity steps, where v < vsc
			signQ = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
			vsc   = Rebin( Sqrt( 2.0 * q * (self.oVsc['DATA']) / self.mass ), dims[[0,3]] )
			VM    = v^2 * (1 + signQ * (vsc^2 / v^2))
		
			;Exclude bins that are below the spacecraft potential
			;   - Total() will not include NaNs
			iBad  = Where(VM LT 0, nBad)
			IF nBad GT 0 THEN BEGIN
				tf_nan   = 1B
				VM[iBad] = !Values.D_NaN
			ENDIF
		ENDELSE
	ENDIF ELSE BEGIN
		VM = v^2
	ENDELSE
	
	;Integrate velocity
	;   - V is in m/s while f is in s^3/cm^6
	;   - v^2 * dv --> (m/s)^3  --> (cm/s)^3 * 1e6
	N = Total( Temporary(ftemp) * Temporary(v) * Sqrt(VM) * Temporary(dv), 2, NAN=tf_nan ) * 1e6

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
	iZero  = self.oDist -> Where(0, /EQUAL, COUNT=nZero)

;-----------------------------------------------------
; Integration Over Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Extract phi and the bin sizes
	oPhi = self.oDist['DEPEND_1']
;	phi  = Rebin(oPhi['DATA'] * deg2rad, dims)
	dPhi = Rebin((oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR'])['DATA'] * deg2rad, dims)
	
	;Integrate
	ftemp = self.oDist['DATA']
	IF nZero GT 0 THEN ftemp[iZero] = 1.0
	ftemp = Total( self.oDist['DATA'] * ALog(Temporary(ftemp)) * Temporary(dPhi), 2)
	
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
	ftemp  = Total( Temporary(ftemp) * Sin(Temporary(theta)) * Temporary(dTheta), 2)

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
		vsc  = Rebin( Sqrt( 2.0 * q / self.mass * (self.oVsc['DATA']) ), dims[[0,3]] )
		VM    = v^2 * (1 + sgn * (vsc^2 / v^2))
		
		;Exclude bins that are below the spacecraft potential
		;   - Total() will not include NaNs
		iBad  = Where(VM LT 0, nBad)
		IF nBad GT 0 THEN BEGIN
			tf_nan   = 1B
			VM[iBad] = !Values.D_NaN
		ENDIF
	ENDIF
	
	;Integrate velocity
	;   - Sign is charge dependent. E = q * V = 1/2 * m * v^2
	;   - V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6  --> s^3/m^6 * 1e12
	IF tf_Vsc $
		THEN H = Total( ftemp * v * Sqrt(VM) * dv, 2, NAN=tf_nan ) * 1e12 $
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
THETA_RANGE=theta_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oESpec)  GT 0 THEN Obj_Destroy, oESpec
		RETURN, !Null
	ENDIF
	
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name)        EQ 0 THEN name        = self.oDist.name + '_ESpec'
	IF N_Elements(phi_range)   EQ 0 THEN phi_range   = [-180.0, 180.0]
	IF N_Elements(theta_range) EQ 0 THEN theta_range = [0.0, 180.0]

	;Allocate memory
	dims      = Size(self.oDist, /DIMENSIONS)
	nTime     = dims[0]
	nPhi      = dims[1]
	nTheta    = dims[2]
	nEnergy   = dims[3]
	Espec     = FltArr(nTime, nEnergy)
	
;-----------------------------------------------------
; Average Over Phi and Theta \\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oPhi    = self.oDist['DEPEND_1']
	oDPhi   = oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR']
	oTheta  = self.oDist['DEPEND_2']
	oDTheta = oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR']
	
	nPts = self.oDist -> GetNPts()
	FOR i = 0, nPts - 1 DO BEGIN
		;Angular bins within desired range
		ip = Where(oPhi['DATA',i,*] GE phi_range[0] AND $
		           oPhi['DATA',i,*] LE phi_range[1], np )
		it = Where(oTheta['DATA',i,*] GE theta_range[0] AND $
		           oTheta['DATA',i,*] LE theta_range[1], nt )
		
		;PHI weight functions
		;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
		;   - w   = v^2 sin(Theta) dTheta dPhi
		;   - v is independent of phi so factors out of the numerator and denominator
		weight = Rebin(oDPhi['DATA',i,ip] * !DToR, 1, np, nTheta, nEnergy, /SAMPLE)
		
		;Average over Phi
		ThetaE = Total( weight * self.oDist['DATA',i,ip,*,*], 2 ) / Total(weight, 2)
		weight = !Null
		
		;THETA weight function
		IF self.elevation $
			THEN weight = Cos( oTheta['DATA',i,it] * !DToR ) * oDTheta['DATA',i,it] * !DToR $
			ELSE weight = Sin( oTheta['DATA',i,it] * !DToR ) * oDTheta['DATA',i,it] * !DToR
		weight = Rebin(weight, 1, nt, nEnergy, /SAMPLE)
		
		;Average over theta
		ESpec[i,*] = Total( weight * ThetaE[0,it,*], 2 ) / Total(weight, 2)
	ENDFOR
	
;-----------------------------------------------------
; Average Over Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Bins over which to average
;	oPhi = self.oDist['DEPEND_1']
;	ip = Where(oPhi['DATA'] GE phi_range[0] AND $
;	           oPhi['DATA'] LE phi_range[1], np )
;	
;	;Compute weights
;	;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
;	;   - w_i = v * dPhi_i
;	;   - v is independent of phi so factors out of the numerator and denominator
;	oDPhi = oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR']
;	weight = Rebin(oDPhi['DATA',*,ip], nTime, np, nTheta, nEnergy, /SAMPLE)
;	
;	;Average over theta
;	ThetaE = Total( weight * self.oDist['DATA',*,ip,*,*], 2 ) / Total(weight, 2)
;	weight = !Null

;-----------------------------------------------------
; Average Over Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
;	;Bins over which to average
;	oTheta = self.oDist['DEPEND_2']
;	it = Where(oTheta['DATA'] GE theta_range[0] AND $
;	           oTheta['DATA'] LE theta_range[1], nt )
;	
;	;Compute weights
;	;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
;	;   - w_i = v * sin(Theta_i) * dTheta_i * dv
;	;   - v and dv are independent of theta so factors out of the numerator and denominator
;	oTheta  = self.oDist['DEPEND_2']
;	oDTheta = oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR']
;	IF self.elevation $
;		THEN weight = Cos( oTheta['DATA',*,it] * !DToR ) * oDTheta['DATA',*,it] * !DToR $
;		ELSE weight = Sin( oTheta['DATA',*,it] * !DToR ) * oDTheta['DATA',*,it] * !DToR
;	weight = Rebin(Rebin(weight, 1, nt, 1, /OVERWRITE), nTime, nt, nEnergy, /SAMPLE)
;	
;	;Average over theta
;	ESpec = Total( weight * ThetaE[*,it,*], 2 ) / Total(weight, 2)
;	weight = !Null

;-------------------------------------------
; Create Variable //////////////////////////
;-------------------------------------------
	
	;Energy-time spectrogram
	oESpec = MrTimeSeries( self.oDist['TIMEVAR'], ESpec, $
	                       CACHE = tf_cache, $
	                       NAME  = name, $
	                       /NO_COPY )

	;Sepctrogram attributes
	oESpec['DEPEND_1'] = self.oDist['DEPEND_3']
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
THETA_RANGE=theta_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oPhiE)    GT 0 THEN Obj_Destroy, oPhiE
		RETURN, !Null
	ENDIF

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Dimension sizes
	dims    = Size(self.oDist, /DIMENSIONS)
	nTime   = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]

	IF N_Elements(name)        EQ 0 THEN name        = self.oDist.name + '_PhiE'
	IF N_Elements(theta_range) EQ 0 THEN theta_range = [0.0, 180.0]

;-----------------------------------------------------
; Average Over Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Bins over which to average
	oTheta  = self.oDist['DEPEND_2']
	it = Where(oTheta['DATA'] GE theta_range[0] AND $
	           oTheta['DATA'] LE theta_range[1], nt )
	
	;Compute weights
	;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
	;   - w_i = v * sin(Theta_i) * dTheta_i
	;   - v is independent of theta so factors out of the numerator and denominator
	oDTheta = oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR']
	IF self.elevation $
		THEN weight = Cos( oTheta['DATA',*,it] * !DToR ) * oDTheta['DATA',*,it] * !DToR $
		ELSE weight = Sin( oTheta['DATA',*,it] * !DToR ) * oDTheta['DATA',*,it] * !DToR
	weight = Rebin(Rebin(weight, nTime, 1, nt, 1, /OVERWRITE), nTime, nPhi, nt, nEnergy, /SAMPLE)
	
	;Average over theta
	PhiE = Total( weight * self.oDist['DATA',*,*,it,*], 3 ) / Total(weight, 3)
	weight = !Null

;-------------------------------------------
; Create Variable //////////////////////////
;-------------------------------------------
	;Theta-Energy distribution
	oPhiE = MrTimeSeries( self.oDist['TIMEVAR'], PhiE, $
	                      CACHE = cache, $
	                      NAME  = name, $
	                      /NO_COPY )
	
	;Attributes
	oPhiE['DEPEND_1'] = self.oDist['DEPEND_1']
	oPhiE['DEPEND_2'] = self.oDist['DEPEND_3']
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
THETA_RANGE=theta_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oPhiSpec) GT 0 THEN Obj_Destroy, oPhiSpec
		RETURN, !Null
	ENDIF

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	eV2J = MrConstants('eV2J')
	IF N_Elements(name)        EQ 0 THEN name        = self.oDist.name + '_PhiSpec'
	IF N_Elements(theta_range) EQ 0 THEN theta_range = [0.0, 180.0]
	IF N_Elements(e_range)     EQ 0 THEN e_range     = [0.0, !values.f_infinity]

	;Allocate memory
	dims    = Size(self.oDist, /DIMENSIONS)
	nTime   = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	phiSpec = FltArr( nTime, nPhi )
	
;-----------------------------------------------------
; Average Over Theta and Energy \\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oTheta = self.oDist['DEPEND_2']
	oDTheta = oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR']
	
	oE  = self.oDist['DEPEND_3']
	oDE = oE['DELTA_PLUS_VAR'] + oE['DELTA_MINUS_VAR']
	v   = Sqrt( 2.0 * eV2J / self.mass * oE['DATA'] )
	dv  = oDE['DATA'] * eV2J / (self.mass * v)
	
	nPts = self.oDist -> GetNPts()
	FOR i = 0, nPts - 1 DO BEGIN
		it = Where(oTheta['DATA',i,*] GE theta_range[0] AND $
		           oTheta['DATA',i,*] LE theta_range[1], nt )
		ie = Where(oE['DATA',i,*] GE e_range[0] AND $
		           oE['DATA',i,*] LE e_range[1], nEn )
		
		;THETA weights
		;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
		;   - w = v * sin(Theta_i) * dTheta * dv
		;   - v and dv are independent of theta so factors out sum over theta
		IF self.elevation $
			THEN weight = Cos( oTheta['DATA',i,it] * !DToR ) * oDTheta['DATA',i,it] * !DToR $
			ELSE weight = Sin( oTheta['DATA',i,it] * !DToR ) * oDTheta['DATA',i,it] * !DToR
		weight = Rebin(Reform(weight, 1, 1, nt, 1, /OVERWRITE), 1, nPhi, nt, nEnergy, /SAMPLE)
		
		;Average over theta
		PhiE = Total( weight * self.oDist['DATA',i,*,it,*], 3 ) / Total(weight, 3)
		weight = !Null
		
		;Compute weights
		weight = Rebin(Reform(v[i,*] * dv[i,*], 1, 1, nEn, /OVERWRITE), 1, nPhi, nEn, /SAMPLE)
		
		;Average over energy
		phiSpec[i,*] = Total( weight * PhiE[0,*,ie], 3 ) / Total(weight, 3)
		weight = !Null
	ENDFOR
	

;-----------------------------------------------------
; Average Over Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
;	;Bins over which to average
;	oTheta = self.oDist['DEPEND_2']
;	it = Where(oTheta['DATA'] GE theta_range[0] AND $
;	           oTheta['DATA'] LE theta_range[1], nt )
;	
;	;Compute weights
;	;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
;	;   - w_i = v * sin(Theta_i) * dTheta_i * dv
;	;   - v and dv are independent of theta so factors out of the numerator and denominator
;	oDTheta = oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR']
;	IF self.elevation $
;		THEN weight = Cos( oTheta['DATA',*,it] * !DToR ) * oDTheta['DATA',*,it] * !DToR $
;		ELSE weight = Sin( oTheta['DATA',*,it] * !DToR ) * oDTheta['DATA',*,it] * !DToR
;	weight = Rebin(Rebin(weight, nTime, 1, nt, 1, /OVERWRITE), nTime, nPhi, nt, nEnergy, /SAMPLE)
;	
;	;Average over theta
;	PhiE = Total( weight * self.oDist['DATA',*,*,it,*], 3 ) / Total(weight, 3)
;	weight = !Null

;-----------------------------------------------------
; Average Over Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
;	oE  = self.oDist['DEPEND_3']
;	oDE = oE['DELTA_PLUS_VAR'] + oE['DELTA_MINUS_VAR']
;	v   = Sqrt( 2.0 * eV2J / self.mass * oE['DATA'] )
;	dv  = oDE['DATA'] * eV2J / (self.mass * v)
;	
;	;Compute weights
;	;   - w = v dv
;	;   - Extra v from theta weights
;	weight = Rebin(Reform(v * dv, 1, 1, nEnergy, /OVERWRITE), nTime, nPhi, nEnergy, /SAMPLE)
;	
;	;Average over energy
;	phiSpec = Total( weight * PhiE, 3 ) / Total(weight, 3)

;-------------------------------------------
; Create Variable //////////////////////////
;-------------------------------------------
	;Phi-time spectrogram
	oPhiSpec = MrTimeSeries( self.oDist['TIMEVAR'], phiSpec, $
	                         CACHE = cache, $
	                         NAME  = name, $
	                         /NO_COPY )

	;Sepctrogram attributes
	oPhiSpec['DEPEND_1']   = self.oDist['DEPEND_1']
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
E_RANGE=E_Range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oPhiTheta)  GT 0 THEN Obj_Destroy, oPhiTheta
		RETURN, !Null
	ENDIF

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_PhiTheta'
	
	;Dimension sizes
	dims    = Size(self.oDist, /DIMENSIONS)
	nTimes  = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]

;-----------------------------------------------------
; Average Over Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oE  = self.oDist['DEPEND_3']
	oDE = oE['DELTA_PLUS_VAR'] + oE['DELTA_MINUS_VAR']
	v   = Sqrt( 2.0 * eV2J / self.mass * oE['DATA'] )
	dv  = oDE['DATA'] * eV2J / (self.mass * v)
	
	;Compute weights
	;   - w = dv
	;   - Extra v from phi weights
	weight = Rebin(Reform(dv, 1, 1, 1, nEnergy, /OVERWRITE), nTime, nPhi, nTheta, nEnergy, /SAMPLE)
	
	;Average over energy
	phiTheta = Total( weight * self.oDist['DATA'], 4 ) / Total(weight, 4)
	
	;Free memory
	weight = !Null
	v      = !Null
	dv     = !Null

;-------------------------------------------
; Create Variable //////////////////////////
;-------------------------------------------
	;Phi-Theta distribution
	oPhiTheta = MrTimeSeries( self.oDist['TIMEVAR'], PhiTheta, $
	                          CACHE = cache, $
	                          NAME  = name, $
	                          /NO_COPY )
	;Attributes
	oPhiTheta['DEPEND_1'] = self.oDist['DEPEND_1']
	oPhiTheta['DEPEND_2'] = self.oDist['DEPEND_2']
	oPhiTheta['SCALE']    = 1B
	oPhiTheta['LOG']      = 1B
	oPhiTheta['UNITS']    = self.oDist['UNITS']
	
	;RETURN the 2D distribution
	RETURN, oPhiTheta
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
GROUND=ground, $
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
	tf_ground = Keyword_Set(ground)
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
	
	;We need density and bulk velocity
	;   - Convert km/s to cm/s
	oN     = 1e6 * self -> Density(GROUND=ground)
	oBulkV = 1e3 * self -> Velocity(GROUND=ground)
	n      = oN['DATA']
	vx     = oBulkV['DATA',*,0]
	vy     = oBulkV['DATA',*,1]
	vz     = oBulkV['DATA',*,2]
	Obj_Destroy, [oN, oBulkV]
	
	;Extract phi and the bin sizes
	oPhi = self.oDist['DEPEND_1']
	phi  = Rebin(oPhi['DATA'] * deg2rad, dims)
	dPhi = Rebin((oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR'])['DATA'] * deg2rad, dims)
	
	;Pxx
	fxx_temp = Total( self.oDist['DATA'] * Cos(phi)^2.0        * dPhi, 2 )
	fxy_temp = Total( self.oDist['DATA'] * Sin(phi) * Cos(phi) * dPhi, 2 )
	fxz_temp = Total( self.oDist['DATA'] * Cos(phi)            * dPhi, 2 )
	fyy_temp = Total( self.oDist['DATA'] * Sin(phi)^2          * dPhi, 2 )
	fyz_temp = Total( self.oDist['DATA'] * Sin(phi)            * dPhi, 2 )
	fzz_temp = Total( self.oDist['DATA']                       * dPhi, 2 )
	
	phi  = !Null
	dPhi = !Null

;-----------------------------------------------------
; Integrate over Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Extract theta and the bin sizes
	oTheta = self.oDist['DEPEND_2']
	theta  = Rebin(oTheta['DATA'] * deg2rad, dims[[0,2,3]])
	dTheta = Rebin((oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR'])['DATA'] * deg2rad, dims[[0,2,3]])
	
	;Pxx
	fxx_temp = Total( fxx_temp * Sin(theta)^3                * dTheta, 2 )
	fxy_temp = Total( fxy_temp * Sin(theta)^3                * dTheta, 2 )
	fxz_temp = Total( fxz_temp * Sin(theta)^2 * Cos(theta)   * dTheta, 2 )
	fyy_temp = Total( fyy_temp * Sin(theta)^3                * dTheta, 2 )
	fyz_temp = Total( fyz_temp * Sin(theta)^2 * Cos(theta)   * dTheta, 2 )
	fzz_temp = Total( fzz_temp * Sin(theta)   * Cos(theta)^2 * dTheta, 2 )

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
	vM  = sqrt(2.0 * eV2J * Temporary(energy) / self.mass)
	dvM = eV2J * Temporary(dEnergy) / (self.mass * vM)
	
	;SPACECRAFT POTENTIAL
	IF N_Elements(self.oVsc) GT 0 THEN BEGIN
		IF tf_ground THEN BEGIN
			vsc = Rebin( Sqrt( 2.0 * q * (self.oVsc['DATA']) / self.mass ), dims[[0,3]] )
			v   = vM - vsc
			
			iBad = Where(v LT 0, nBad)
			IF nBad GT 0 THEN BEGIN
				tf_nan  = 1B
				v[iBad] = 0;!Values.D_NaN
			ENDIF
				
			vMc = v^2
			
		ENDIF ELSE BEGIN
			;Sign is charge dependent. E = q * V = 1/2 * m * v^2
			signQ = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
			vsc   = Rebin( Sqrt( 2.0 * q / self.mass * self.oVsc['DATA'] ), dims[[0,3]] )
			vMc   = vM^2 * (1 + signQ * (vsc^2 / vM^2))
		
			;Exclude bins that are below the spacecraft potential
			;   - Total() will not include NaNs
			tf_nan = 0B
			iBad   = Where(vMc LT 0, nBad)
			IF nBad GT 0 THEN BEGIN
				tf_nan = 1B
				vMc[iBad] = !Values.D_NaN
			ENDIF
		ENDELSE
	ENDIF ELSE BEGIN
		vMc  = vM^2
	ENDELSE
	
	
	;
	; VELOCITY
	;

	;V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
	
	;Pxx
	fxx_temp = Total( fxx_temp * vM * vMc^(3.0/2.0) * dvM, 2, NAN=tf_nan ) * 1e12
	fxy_temp = Total( fxy_temp * vM * VMc^(3.0/2.0) * dvM, 2, NAN=tf_nan ) * 1e12
	fxz_temp = Total( fxz_temp * vM * VMc^(3.0/2.0) * dvM, 2, NAN=tf_nan ) * 1e12
	fyy_temp = Total( fyy_temp * vM * VMc^(3.0/2.0) * dvM, 2, NAN=tf_nan ) * 1e12
	fyz_temp = Total( fyz_temp * vM * VMc^(3.0/2.0) * dvM, 2, NAN=tf_nan ) * 1e12
	fzz_temp = Total( fzz_temp * vM * VMc^(3.0/2.0) * dvM, 2, NAN=tf_nan ) * 1e12
	
	Pxx = self.mass * (Temporary(fxx_temp) - n * vx * vx)
	Pxy = self.mass * (Temporary(fxy_temp) - n * vx * vy)
	Pxz = self.mass * (Temporary(fxz_temp) - n * vx * vz)
	Pyy = self.mass * (Temporary(fyy_temp) - n * vy * vy)
	Pyz = self.mass * (Temporary(fyz_temp) - n * vy * vz)
	Pzz = self.mass * (Temporary(fzz_temp) - n * vz * vz)
	
	vM  = !Null
	dvM = !Null
	vMc = !Null
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
;   Rotate the distribution.
;
; :Params:
;       TMATRIX:        in, required, type=TxNxMxL fltarr
;                       The transformation matrix used to rotate the distribution
;       ORIENTATION:    in, optional, type=boolean, default=1
;                       Orientation of `THETA` and `PHI`. Options are::
;                         1: PHI   - Positive from x-axis
;                            THETA - Polar angle from z-axis
;                         2: PHI   - Positive from y-axis
;                            THETA - Polar angle from z-axis
;                         3: PHI   - Positive from x-axis
;                            THETA - Elevation angle from xy-plane
;                         4: PHI   - Positive from y-axis
;                            THETA - Elevation angle from xy-plane
;                         5: PHI   - Positive from z-axis
;                            THETA - Polar angle from y-axis
;                         6: PHI   - Positive from x-axis
;                            THETA - Polar angle from y-axis
;                         7: PHI   - Positive from z-axis
;                            THETA - Elevation angle from zx-plane
;                         8: PHI   - Positive from x-axis
;                            THETA - Elevation angle from zx-plane
;                         9: PHI   - Positive from y-axis
;                            THETA - Polar angle from x-axis
;                        10: PHI   - Positive from z-axis
;                            THETA - Polar angle from x-axis
;                        11: PHI   - Positive from y-axis
;                            THETA - Elevation angle from yz-plane
;                        12: PHI   - Positive from z-axis
;                            THETA - Elevation angle from yz-plane
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output is added to the variable cache.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;
; :Returns:
;       ODIST:          out, required, type=TxNxMxL fltarr
;                       The re-binned distribution function.
;
;-
FUNCTION MrDist4D::Rotate, tmatrix, orientation
	Compile_Opt idl2
	On_Error, 2
	
	oT = MrVar_Get(tmatrix)
	IF N_Elements(orientation) EQ 0 THEN orientation = 1
	
	;Theta is an elevation angle
	CASE orientation OF
		 3: tf_elevation = 1B
		 4: tf_elevation = 1B
		 7: tf_elevation = 1B
		 8: tf_elevation = 1B
		11: tf_elevation = 1B
		12: tf_elevation = 1B
		ELSE: tf_elevation = 0B
	ENDCASE
	
	;Dimensions
	dims    = Size(self.oDist, /DIMENSIONS)
	nTime   = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	
;-------------------------------------------
; Rotate Angles ////////////////////////////
;-------------------------------------------
	
	;Create a cartesian grid
	MrVar_Grid_MakeCart, self.oDist['DEPEND_1'], self.oDist['DEPEND_2'], oX, oY, oZ, /DEGREES

	;Rotate the cartesian grid
	;   - Negative sign converts from look-dirction to incident trajectory
	MrVar_Grid_Cart2FAC, oT, -temporary(oX), -temporary(oY), -temporary(oZ), oX1, oX2, oX3
	

	;Convert to polar grid
	MrVar_Grid_cart2sphere, oX1, oX2, oX3, oPhi, oTheta, $
	                        /DEGREES, $
	                        ORIENTATION = orientation

;-------------------------------------------
; Rebin Angles /////////////////////////////
;-------------------------------------------
	;Obtain the volume elements of the old configuration to use as weights when
	;rebinning in the new configuration.
	oDV = self -> VolumeElement()
	
	nPhi_bins = 32
	phi_range = [-180.0, 180.0]
	dPhi      = (phi_range[1] - phi_range[0]) / nPhi
	phi_bins  = dPhi * FIndGen(nPhi) + phi_range[0]
	
	nTheta_bins = 16
	theta_range = tf_elevation ? [-90.0, 90.0] : [0.0, 180.0]
	dTheta      = (theta_range[1] - theta_range[0]) / nTheta
	theta_bins  = dTheta * FIndGen(nTheta) + theta_range[0]
	
	rDist = FltArr(nTime, nPhi_bins, nTheta_bins, nEnergy)
	
	;Loop through time
	FOR i = 0, nTime - 1 DO BEGIN
		coords = [ Reform(oPhi['DATA',i,*,*],   1, nPhi*nTheta), $
		           Reform(oTheta['DATA',i,*,*], 1, nPhi*nTheta) ]
		
		;Locate data within new bins
		;   - Keep the number of Theta bins the same
		;   - Average all PHI values within [-DELTA, DELTA] in each theta bin
		cHist  = hist_nd( coords, $
		                 MIN             = [ phi_range[0], theta_range[0] ], $
		                 MAX             = [ phi_range[1], theta_range[1] ], $
		                 NBINS           = [    nPhi_Bins,    nTheta_Bins ], $
		                 REVERSE_INDICES = ri)
		
		;Loop over bins
		FOR j = 0, nEnergy - 1 DO BEGIN
			FOR k = 0, N_Elements(chist) - 1 DO BEGIN
				;Skip empty bins
				IF ri[k] EQ ri[k+1] THEN CONTINUE
			
				;Source indices
				inds = ri[ri[k]:ri[k+1]-1]
				isrc = Array_Indices([nPhi,nTheta], inds, /DIMENSIONS)
			
				;Destination indices
				idest = Array_Indices([nPhi_bins,nTheta_bins], k, /DIMENSIONS)
		
				;Weight
				w = Reform(oDV['DATA', i, isrc[0,*], isrc[1,*], j])
			
				;Re-bin
				rDist[i, idest[0], idest[1], j] = Total( Reform(self.oDist['DATA',i, isrc[0,*], isrc[1,*], j]) * w ) / Total(w)
			ENDFOR
		ENDFOR
	ENDFOR
	

;-------------------------------------------
; Create New MrDist4D Object ///////////////
;-------------------------------------------
	outTime = self.oDist['DEPEND_0']
	outPhi  = MrVariable(phi_bins)
	outPhi['CATDESC']     = 'Azimuthal angle'
	outPhi['DELTA_PLUS']  = dPhi
	outPhi['DELTA_MINUS'] = dPhi
	outPhi['TITLE']         = 'Phi'
	outPhi['UNITS']         = 'degrees'
	outPhi['SI_CONVERSION'] = '0.017453292>radians'
	
	outTheta = MrVariable(theta_bins)
	outTheta['CATDESC']       = (tf_elevation ? 'Polar angle' : 'Elevation angle')
	outTheta['DELTA_PLUS']    = dTheta
	outTheta['DELTA_MINUS']   = dTheta
	outTheta['TITLE']         = 'Theta'
	outTheta['UNITS']         = 'degrees'
	outTheta['SI_CONVERSION'] = '0.017453292>radians'
	
	outDist = MrTimeSeries(outTime, rDist, /NO_COPY)
	outDist['DEPEND_1'] = outPhi
	outDist['DEPEND_2'] = outTheta
	outDist['DEPEND_3'] = self.oDist['DEPEND_3']
	outDist['UNITS']    = self.oDist['UNITS']

	;Create new distribution
	oDist = MrDist4D( outDist, $
	                  ELEVATION = tf_elevation, $
	                  SPECIES   = theSpecies, $
	                  UNITS     = self.units )
	
	
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
;   Set the energy array.
;
;   If the DELTA_PLUS and DELTA_MINUS attributes do not exist, they are created and set
;   such that dE/E is constant and is determined from the first two energy bins. If the
;   attribute values are not time-dependent, they are expanded and turned into
;   MrTimeSeries objects.
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
;   If the DELTA_PLUS and DELTA_MINUS attributes do not exist, they are created and set to
;   half the mean spacing between points. If the attribute values are not time-dependent,
;   they are expanded and turned into MrTimeSeries objects.
;
;   If the UNITS attribute is radians, values are converted to degrees.
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
;   Set the theta array.
;
;   If the DELTA_PLUS and DELTA_MINUS attributes do not exist, they are created and set to
;   half the mean spacing between points. If the attribute values are not time-dependent,
;   they are expanded and turned into MrTimeSeries objects.
;
;   If the UNITS attribute is radians, values are converted to degrees.
;
; :Params:
;       THETA:          in, required, type=Nx1 or NxM array
;                       Polar coordinates of the distribution pixels. Can be the name
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
; Make Delta Plus Time-Dependent \\\\\\\\\\\\\\\\\\\\\
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
; Make Delta Minus Time-Dependent \\\\\\\\\\\\\\\\\\\\
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
; Make Units Degrees \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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
pro MrDist4D::Spectra, $
CACHE=cache, $
E_RANGE=e_range, $
ESPEC=oESpec, $
;NAME=name, $
PHI_RANGE=phi_range, $
PHISPEC=oPhiSpec, $
THETA_RANGE=theta_range, $
THETASPEC=oThetaSpec
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oESpec)     GT 0 THEN Obj_Destroy, oESpec
		IF N_Elements(oPhiSpec)   GT 0 THEN Obj_Destroy, oPhiSpec
		IF N_Elements(oThetaSpec) GT 0 THEN Obj_Destroy, oThetaSpec
		RETURN
	ENDIF

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	eV2J     = MrConstants('eV2J')
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name)        EQ 0 THEN name        = self.oDist.name + '_ThetaSpec'
	IF N_Elements(phi_range)   EQ 0 THEN phi_range   = [-180.0, 180.0]
	IF N_Elements(theta_range) EQ 0 THEN theta_range = [0.0, 180.0]
	IF N_Elements(e_range)     EQ 0 THEN e_range     = [0.0, !values.f_infinity]

	;Dimension sizes
	dims      = Size(self.oDist, /DIMENSIONS)
	nTime     = dims[0]
	nPhi      = dims[1]
	nTheta    = dims[2]
	nEnergy   = dims[3]
	
	phiSpec   = FltArr( nTime, nPhi )
	thetaSpec = FltArr( nTime, nTheta )
	Espec     = FltArr( nTime, nEnergy )
	
;-----------------------------------------------------
; Average Over Phi and Energy \\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oPhi = self.oDist['DEPEND_1']
	oDPhi = oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR']
	
	oTheta  = self.oDist['DEPEND_2']
	oDTheta = oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR']
	
	oE  = self.oDist['DEPEND_3']
	oDE = oE['DELTA_PLUS_VAR'] + oE['DELTA_MINUS_VAR']
	v   = Sqrt( 2.0 * eV2J / self.mass * oE['DATA'] )
	dv  = oDE['DATA'] * eV2J / (self.mass * v)
	
	nPts = self.oDist -> GetNPts()
	FOR i = 0, nPts - 1 DO BEGIN
		;Bins within desired range
		ip = Where(oPhi['DATA',i,*] GE phi_range[0] AND $
		           oPhi['DATA',i,*] LE phi_range[1], np )
		it = Where(oTheta['DATA',i,*] GE theta_range[0] AND $
		           oTheta['DATA',i,*] LE theta_range[1], nt )
		ie = Where(oE['DATA',i,*] GE e_range[0] AND $
		           oE['DATA',i,*] LE e_range[1], nEn )

	;-----------------------------------------------------
	; 3D Spectra \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		
		;THETA-ENERGY
		;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
		;   - w = v * dPhi * dv
		;   - v and dv are independent of phi so factors out of the sum over phi
		w_phi  = Rebin(oDPhi['DATA',i,ip], 1, np, nTheta, nEnergy, /SAMPLE)
		ThetaE = Total( w_phi * self.oDist['DATA',i,ip,*,*], 2 ) / Total(w_phi, 2)
		w_phi  = !Null
	
		;PHI-ENERGY
		IF self.elevation $
			THEN w_theta = Cos( oTheta['DATA',i,it] * !DToR ) * oDTheta['DATA',i,it] * !DToR $
			ELSE w_theta = Sin( oTheta['DATA',i,it] * !DToR ) * oDTheta['DATA',i,it] * !DToR
		w_theta = Rebin(Reform(w_theta, 1, 1, nt, 1, /OVERWRITE), 1, nPhi, nt, nEnergy, /SAMPLE)
		PhiE    = Total( w_theta * self.oDist['DATA',i,*,it,*], 3 ) / Total(w_theta, 3)
		
		;PHI-THETA
		;   - w = dv
		;   - Extra v from phi weights
		w_energy = Rebin(Reform(dv[i,ie], 1, 1, 1, nEn, /OVERWRITE), 1, nPhi, nTheta, nEn, /SAMPLE)
		PhiTheta = Total( w_energy * self.oDist['DATA',i,*,*,ie], 4 ) / Total(w_energy, 4)

	;-----------------------------------------------------
	; 2D Spectra \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		
		;PHI
		w_energy     = Rebin(Reform(v[i,*] * dv[i,*], 1, 1, nEn, /OVERWRITE), 1, nPhi, nEn, /SAMPLE)
		phiSpec[i,*] = Total( w_energy * PhiE[0,*,ie], 3 ) / Total(w_energy, 3)
		
		;THETA
		w_energy       = Rebin(Reform(v[i,ie] * dv[i,ie], 1, 1, nEn, /OVERWRITE), 1, nTheta, nEn, /SAMPLE)
		thetaSpec[i,*] = Total( w_energy * ThetaE[0,*,ie], 3 ) / Total(w_energy, 3)
		
		;ENERGY
		IF self.elevation $
			THEN w_theta = Cos( oTheta['DATA',i,it] * !DToR ) * oDTheta['DATA',i,it] * !DToR $
			ELSE w_theta = Sin( oTheta['DATA',i,it] * !DToR ) * oDTheta['DATA',i,it] * !DToR
		w_theta    = Rebin(w_theta, 1, nt, nEnergy, /SAMPLE)
		ESpec[i,*] = Total( w_theta * ThetaE[0,it,*], 2 ) / Total(w_theta, 2)
	ENDFOR

;-------------------------------------------
; Energy ///////////////////////////////////
;-------------------------------------------
	
	;Energy-time spectrogram
	oESpec = MrTimeSeries( self.oDist['TIMEVAR'], ESpec, $
	                       CACHE = tf_cache, $
	                       NAME  = self.oDist.name + '_espec', $
	                       /NO_COPY )

	;Sepctrogram attributes
	oESpec['DEPEND_1'] = self.oDist['DEPEND_3']
	oESpec['SCALE']    = 1B
	oESpec['LOG']      = 1B
	oESpec['UNITS']    = self.oDist['UNITS']
	
;-------------------------------------------
; Phi //////////////////////////////////////
;-------------------------------------------
	;Phi-time spectrogram
	oPhiSpec = MrTimeSeries( self.oDist['TIMEVAR'], phiSpec, $
	                         CACHE = tf_cache, $
	                         NAME  = self.oDist.name + '_phispec', $
	                         /NO_COPY )

	;Sepctrogram attributes
	oPhiSpec['DEPEND_1']   = self.oDist['DEPEND_1']
	oPhiSpec['SCALE']      = 1B
	oPhiSpec['LOG']        = 1B
	oPhiSpec['UNITS']      = self.oDist['UNITS']
	oPhiSpec['TITLE']      = 'Phi Dist'
	oPhiSpec['PLOT_TITLE'] = 'Distribution in Phi'

;-------------------------------------------
; Theta ////////////////////////////////////
;-------------------------------------------
	
	;Theta-time spectrogram
	oThetaSpec = MrTimeSeries( self.oDist['TIMEVAR'], thetaSpec, $
	                           CACHE = tf_cache, $
	                           NAME  = self.oDist.name + '_thetaspec', $
	                           /NO_COPY )
	
	;Attributes
	oThetaSpec['DEPEND_1']   = self.oDist['DEPEND_2']
	oThetaSpec['SCALE']      = 1B
	oThetaSpec['LOG']        = 1B
	oThetaSpec['UNITS']      = self.oDist['UNITS']
	oThetaSpec['TITLE']      = 'Theta Dist'
	oThetaSpec['PLOT_TITLE'] = 'Distribution in Theta'
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
PHI_RANGE=phi_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oThetaE)    GT 0 THEN Obj_Destroy, oThetaE
		RETURN, !Null
	ENDIF

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF N_Elements(name)      EQ 0 THEN name      = self.oDist.name + '_ThetaE'
	IF N_Elements(phi_range) EQ 0 THEN phi_range = [-180.0, 180.0]
	
	;Dimension sizes
	dims    = Size(self.oDist, /DIMENSIONS)
	nTimes  = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	
;-----------------------------------------------------
; Average Over Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Bins over which to average
	oPhi = self.oDist['DEPEND_11']
	ip = Where(oPhi['DATA'] GE phi_range[0] AND $
	           oPhi['DATA'] LE phi_range[1], np )
	
	;Compute weights
	;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
	;   - w_i = v * dPhi_i
	;   - v is independent of phi so factors out of the numerator and denominator
	oDPhi = oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR']
	weight = Rebin(oDPhi['DATA',*,ip], nTime, np, nTheta, nEnergy, /SAMPLE)
	
	;Average over theta
	ThetaE = Total( weight * self.oDist['DATA',*,ip,*,*], 2 ) / Total(weight, 2)
	weight = !Null

;-------------------------------------------
; Create Variable //////////////////////////
;-------------------------------------------
	;Theta-Energy distribution
	oThetaE = MrTimeSeries( self.oDist['TIMEVAR'], ThetaE, $
	                        CACHE = cache, $
	                        NAME  = name, $
	                        /NO_COPY )

	;Attributes
	oThetaE['DEPEND_1'] = self.oDist['DEPEND_2']
	oThetaE['DEPEND_2'] = self.oDist['DEPEND_3']
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
PHI_RANGE=phi_range
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		IF N_Elements(oThetaSpec) GT 0 THEN Obj_Destroy, oThetaSpec
		RETURN, !Null
	ENDIF

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	eV2J     = MrConstants('eV2J')
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name)      EQ 0 THEN name      = self.oDist.name + '_ThetaSpec'
	IF N_Elements(phi_range) EQ 0 THEN phi_range = [-180.0, 180.0]
	IF N_Elements(e_range)   EQ 0 THEN e_range   = [0.0, !values.f_infinity]

	;Dimension sizes
	dims      = Size(self.oDist, /DIMENSIONS)
	nTime     = dims[0]
	nPhi      = dims[1]
	nTheta    = dims[2]
	nEnergy   = dims[3]
	thetaSpec = FltArr( nTime, nTheta )
	
;-----------------------------------------------------
; Average Over Phi and Energy \\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oPhi = self.oDist['DEPEND_1']
	oDPhi = oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR']
	
	oE  = self.oDist['DEPEND_3']
	oDE = oE['DELTA_PLUS_VAR'] + oE['DELTA_MINUS_VAR']
	v   = Sqrt( 2.0 * eV2J / self.mass * oE['DATA'] )
	dv  = oDE['DATA'] * eV2J / (self.mass * v)
	
	nPts = self.oDist -> GetNPts()
	FOR i = 0, nPts - 1 DO BEGIN
		;Bins within desired range
		ip = Where(oPhi['DATA'] GE phi_range[0] AND $
		           oPhi['DATA'] LE phi_range[1], np )
		ie = Where(oE['DATA',i,*] GE e_range[0] AND $
		           oE['DATA',i,*] LE e_range[1], nEn )
		
		;PHI weights
		;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
		;   - w = v * dPhi * dv
		;   - v and dv are independent of phi so factors out of the sum over phi
		weight = Rebin(oDPhi['DATA',i,ip], 1, np, nTheta, nEnergy, /SAMPLE)
		
		;Average over phi
		ThetaE = Total( weight * self.oDist['DATA',i,ip,*,*], 2 ) / Total(weight, 2)
		weight = !Null
	
		;ENERGY weights
		weight = Rebin(Reform(v[i,ie] * dv[i,ie], 1, 1, nEn, /OVERWRITE), 1, nTheta, nEn, /SAMPLE)

		;Average over energy
		thetaSpec[i,*] = Total( weight * ThetaE[0,*,ie], 3 ) / Total(weight, 3)
	ENDFOR

;-----------------------------------------------------
; Average Over Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
;	;Bins over which to average
;	oPhi = self.oDist['DEPEND_1']
;	ip = Where(oPhi['DATA'] GE phi_range[0] AND $
;	           oPhi['DATA'] LE phi_range[1], np )
;	
;	;Compute weights
;	;   - avg = Sum_i=0^nTheta f_i w_i / Sum w_i, where w_i are the weights
;	;   - w_i = v * dPhi_i * dv
;	;   - v and dv are independent of phi so factors out of the numerator and denominator
;	oDPhi = oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR']
;	weight = Rebin(oDPhi['DATA',*,ip], nTime, np, nTheta, nEnergy, /SAMPLE)
;	
;	;Average over theta
;	ThetaE = Total( weight * self.oDist['DATA',*,ip,*,*], 2 ) / Total(weight, 2)
;	weight = !Null

;-----------------------------------------------------
; Average Over Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
;	oE  = self.oDist['DEPEND_3']
;	oDE = oE['DELTA_PLUS_VAR'] + oE['DELTA_MINUS_VAR']
;	v   = Sqrt( 2.0 * eV2J / self.mass * oE['DATA'] )
;	dv  = oDE['DATA'] * eV2J / (self.mass * v)
;	
;	;Compute weights
;	;   - w = v dv
;	;   - Extra v from phi weights
;	weight = Rebin(Reform(v * dv, 1, 1, nEnergy, /OVERWRITE), nTime, nTheta, nEnergy, /SAMPLE)
;	
;	;Average over energy
;	thetaSpec = Total( weight * ThetaE, 3 ) / Total(weight, 3)

;-------------------------------------------
; Create Variable //////////////////////////
;-------------------------------------------
	
	;Theta-time spectrogram
	oThetaSpec = MrTimeSeries( self.oDist['TIMEVAR'], thetaSpec, $
	                           CACHE = tf_cache, $
	                           NAME  = name, $
	                           /NO_COPY )
	
	;Attributes
	oThetaSpec['DEPEND_1']   = self.oDist['DEPEND_2']
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
	eV2K = MrConstants('eV2K')
	
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
GROUND=ground, $
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
	tf_ground = Keyword_Set(ground)
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
		IF tf_ground THEN BEGIN
			vsc = Rebin( Sqrt( 2.0 * q * (self.oVsc['DATA']) / self.mass ), dims[[0,3]] )
			v  -= vsc
			
			iBad = Where(v LT 0, nBad)
			IF nBad GT 0 THEN BEGIN
				tf_nan  = 1B
				v[iBad] = 0;!Values.D_NaN
			ENDIF
				
			VM  = v^2
		
		ENDIF ELSE BEGIN
			;Sign is charge dependent. E = q * V = 1/2 * m * v^2
			sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
			vsc  = Rebin( Sqrt( 2.0 * q * (self.oVsc['DATA']) / self.mass ), dims[[0,3]] )
			VM   = v^2 * (1 + sgn * (vsc^2 / v^2))
		
			;Exclude bins that are below the spacecraft potential
			;   - Total() will not include NaNs
			tf_nan = 0B
			iBad   = Where(VM LT 0, nBad)
			IF nBad GT 0 THEN BEGIN
				tf_nan = 1B
				VM[iBad] = !Values.D_NaN
			ENDIF
		ENDELSE
		
		;Integrate
		Vx = Total( Temporary(fx_temp) * v * VM * dv, 2, NAN=tf_nan ) * 1e12
		Vy = Total( Temporary(fy_temp) * v * VM * dv, 2, NAN=tf_nan ) * 1e12
		Vz = Total( Temporary(fz_temp) * v * VM * dv, 2, NAN=tf_nan ) * 1e12
;		Vz = Total( Temporary(fz_temp) * v * (v^2 + sgn*vsc^2) * dv, 2 ) * 1e12
	
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
	oN = 1e6 * self -> Density(GROUND=ground)
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
		IF N_Elements(odV) GT 0 THEN Obj_Destroy, odV
		RETURN, !Null
	ENDIF
	
	;Defaults
	eV2J     = MrConstants('eV2J')
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name) EQ 0 THEN name = self.oDist.name + '_dV'

	;Allocate memory
	dims    = Size(self.oDist, /DIMENSIONS)
	nTime   = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	dV    = FltArr( nTime, nPhi, nTheta, nEnergy )

;-----------------------------------------------------
; Create Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	oPhi  = self.oDist['DEPEND_1']
	oDPhi = (oPhi['DELTA_PLUS_VAR'] + oPhi['DELTA_MINUS_VAR'])
	
	oTheta  = self.oDist['DEPEND_2']
	oDTheta = (oTheta['DELTA_PLUS_VAR'] + oTheta['DELTA_MINUS_VAR'])
	
	oE  = self.oDist['DEPEND_3']
	oDE = oE['DELTA_PLUS_VAR'] + oE['DELTA_MINUS_VAR']
	v   = Sqrt( 2.0 * eV2J / self.mass * oE['DATA'] )
	dv  = oDE['DATA'] * eV2J / (self.mass * v)
	
	;Volume
	;   - dV = v^2 Sin(theta) dv dTheta dPhi
	IF self.elevation THEN BEGIN
		dV = Rebin(Reform(v^2,                       nTime,    1,      1, nEnergy), nTime, nPhi, nTheta, nEnergy) $
		   * Rebin(Reform(Cos(oTheta['DATA']*!DToR), nTime,    1, nTheta,       1), nTime, nPhi, nTheta, nEnergy) $
		   * Rebin(Reform(dv,                        nTime,    1,      1, nEnergy), nTIme, nPhi, nTheta, nEnergy) $
		   * Rebin(Reform(oDTheta['DATA']*!DToR,     nTime,    1, nTheta,       1), nTIme, nPhi, nTheta, nEnergy) $
		   * Rebin(Reform(oDPhi['DATA']*!DToR,       nTime, nPhi,      1,       1), nTIme, nPhi, nTheta, nEnergy)
	ENDIF ELSE BEGIN
		dV = Rebin(Reform(v^2,                       nTime,    1,      1, nEnergy), nTime, nPhi, nTheta, nEnergy) $
		   * Rebin(Reform(Sin(oTheta['DATA']*!DToR), nTime,    1, nTheta,       1), nTime, nPhi, nTheta, nEnergy) $
		   * Rebin(Reform(dv,                        nTime,    1,      1, nEnergy), nTIme, nPhi, nTheta, nEnergy) $
		   * Rebin(Reform(oDTheta['DATA']*!DToR,     nTime,    1, nTheta,       1), nTIme, nPhi, nTheta, nEnergy) $
		   * Rebin(Reform(oDPhi['DATA']*!DToR,       nTime, nPhi,      1,       1), nTIme, nPhi, nTheta, nEnergy)
	ENDELSE

;-----------------------------------------------------
; Create Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	odV = MrTimeSeries(self.oDist['TIMEVAR'], dV, /NO_COPY)
	
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
	          Inherits IDL_Container, $
	          elevation: 0B, $
	          mass:      0.0, $
	          oVsc:      Obj_New(), $
	          oDist:     Obj_New(), $
	          species:   '', $
	          units:     '' $
	        }
END