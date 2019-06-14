; docformat = 'rst'
;
; NAME:
;       MrVar_Dist2D_Vel
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
;   Turn a 2D [azimuth,energy] distribution into a 2D [velocity,azimuth] distribution
;   with attributes set such that it can be plotted directly with MrVar_Image. The
;   following attributes of the input distribution are used:
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
;           'DELTA_PLUS_VAR'  - Upper-bound of each energy bin
;           'DELTA_PLUS'      - Uniform upper-bound of the energy bin (when DELTA_PLUS_VAR not present)
;           'DELTA_MINUS_VAR' - Lower-bound of each energy bin
;           'DELTA_MINUS'     - Uniform lower-bound of the energy bin (when DELTA_MINUS_VAR not present)
;
;       Phi:
;           'UNITS'           - Angle units. Used to convert to radians
;           'DELTA_PLUS_VAR'  - Upper-bound of each energy bin
;           'DELTA_PLUS'      - Uniform upper-bound of the energy bin (when DELTA_PLUS_VAR not present)
;           'DELTA_MINUS_VAR' - Lower-bound of each energy bin
;           'DELTA_MINUS'     - Uniform lower-bound of the energy bin (when DELTA_MINUS_VAR not present)
;
;   Calling Sequence:
;       img = MrVar_Dist_Image(theDist, mass)
;       img = MrVar_Dist_Image(theDist, mass, idx)
;
; :Params:
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
; :Keywords:
;       _REF_EXTRA: in, optional, type=any
;                   Any keyword accepted by MrVar_Image.pro
;
; :Returns:
;       Dist2D:     out, required, type=object
;                   A MrVariable object containing the 2D velocity-space distribution.
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
;       2017/02/11  -   Written by Matthew Argall
;       2018/05/24  -   Removed the IDX parameter to maintain time dependence of the
;                           distribution function. - MRA
;-
;*****************************************************************************************
;+
;   Convert a distribution from [phi,energy] to [energy,phi]. If the input distribution
;   is also time-dependent, a single distribution is extracted.
;
; :Params:
;       THEDIST:    in, required, type=string/integer/objref
;                   A MrVariable object, name or index of distribution function
;                       to be displayed. Dimensions should be ordered as
;                       [time, azimuth, energy], with the time dimension being optional.
;                       Also must have corresponding DEPEND_[0-2] attributes.
;       IDX:        in, optional, type=integer, default=0
;                   If `THEDIST` is time-dependent, then this is the index along the
;                       time dimension of the distribution to be displayed.
;
; :Returns:
;       oDist:      out, required, type=object
;                   A MrVariable object containing the 2D velocity-space distribution.
;-
FUNCTION MrVar_Dist2D_Prep_Pick, oDistIn, idx
	Compile_Opt idl2
	On_Error, 2
	
;-------------------------------------------
; Time-Dependent Distributions /////////////
;-------------------------------------------
	;Time-Dependent distribution
	;   - Select one
	IF Obj_IsA(oDistIn, 'MrTimeSeries') THEN BEGIN
		
		;Extract one time
		;   - oDist = oDistIn[idx,*,*] will do all of this work automatically,
		;     but will leave MrTimeSeries variables as time-dependent.
		;   - Our goal is to get time-independent MrVariable objects.
		IF N_Elements(idx) EQ 0 THEN idx = 0
		oDist = MrVariable( Reform(oDistIn['DATA',idx,*,*]) )

	;-------------------------------------------
	; PHI: Single Element //////////////////////
	;-------------------------------------------
		oDep1 = MrVar_Get(oDistIn['DEPEND_1'])
		IF Obj_IsA(oDep1, 'MrTimeSeries') THEN BEGIN
			;Extract single time
			oPhi   = MrVariable( Reform(oDep1['DATA',idx,*]) )
			oDep1 -> CopyAttrTo, oPhi
			oPhi  -> RemoveAttr, 'DEPEND_0'
			
			;DELTA_PLUS_VAR
			IF oPhi -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
				;Extract single time
				oDelta     = MrVar_Get(oPhi['DELTA_PLUS_VAR'])
				oPhi_dPlus = MrVariable( Reform(oDelta['DATA',idx,*]) )
				
				;Copy attributes
				oDelta     -> CopyAttrTo, oPhi_dPlus
				oPhi_dPlus -> RemoveAttr, 'DEPEND_0'
				
				;Set as delta variable
				oPhi['DELTA_PLUS_VAR'] = oPhi_dPlus
			ENDIF
			
			;DELTA_MINUS_VAR
			IF oPhi -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
				;Extract single time
				oDelta      = MrVar_Get(oPhi['DELTA_MINUS_VAR'])
				oPhi_dMinus = MrVariable( Reform(oDelta['DATA',idx,*]) )
				
				;Copy attributes
				oDelta      -> CopyAttrTo, oPhi_dMinus
				oPhi_dMinus -> RemoveAttr, 'DEPEND_0'
				
				;Set as delta variable
				oPhi['DELTA_MINUS_VAR'] = oPhi_dMinus
			ENDIF
		ENDIF ELSE BEGIN
			oPhi = oDep1
		ENDELSE

	;-------------------------------------------
	; ENERGY: Single Element ///////////////////
	;-------------------------------------------
		oDep2 = MrVar_Get(oDistIn['DEPEND_2'])
		IF Obj_IsA(oDep2, 'MrTimeSeries') THEN BEGIN
			;Extract single time
			oEnergy  = MrVariable( Reform(oDep2['DATA',idx,*]) )
			oDep2   -> CopyAttrTo, oEnergy
			oEnergy -> RemoveAttr, 'DEPEND_0'
			
			;DELTA_PLUS_VAR
			IF oEnergy -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
				;Extract single time
				oDelta    = MrVar_Get(oEnergy['DELTA_PLUS_VAR'])
				oE_dPlus  = MrVariable( Reform(oDelta['DATA',idx,*]) )
				
				;Copy attributes
				oDelta   -> CopyAttrTo, oE_dPlus
				oE_dPlus -> RemoveAttr, 'DEPEND_0'
				
				;Set as delta variable
				oEnergy['DELTA_PLUS_VAR'] = oE_dPlus
			ENDIF
			
			;DELTA_MINUS_VAR
			IF oEnergy -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
				;Extract single time
				oDelta    = MrVar_Get(oEnergy['DELTA_MINUS_VAR'])
				oE_dMinus = MrVariable( Reform(oDelta['DATA',idx,*]) )
				
				;Copy attributes
				oDelta   -> CopyAttrTo, oE_dMinus
				oE_dMinus -> RemoveAttr, 'DEPEND_0'
				
				;Set as delta variable
				oEnergy['DELTA_MINUS_VAR'] = oE_dMinus
			ENDIF
		ENDIF ELSE BEGIN
			oEnergy = oDep2
		ENDELSE
		
		;Add attributes
		;   - Want to make polar plot, so order as (v,phi)
		;Add attributes
		oDistIn -> CopyAttrTo, oDist
		oDist   -> RemoveAttr, 'DEPEND_' + ['0', '1', '2']
		oDist   -> SetData, Transpose(oDist['DATA'])
		oDist['DEPEND_0'] = oEnergy
		oDist['DEPEND_1'] = oPhi

;-------------------------------------------
; Time-Independent Distribution ////////////
;-------------------------------------------
	ENDIF ELSE BEGIN
		;Rearrange dimensions from [phi, energy] to [energy, phi]
		oDist  = oDistIn -> Copy()
		temp   = oDist['DEPEND_0']
		oDist -> SetData, Transpose(oDist['DATA'])
		oDist['DEPEND_0'] = oDist['DEPEND_1']
		oDist['DEPEND_1'] = Temporary(temp)
	ENDELSE

	RETURN, oDist
END


;+
;   Convert energy to velocity.
;
; :Params:
;       ODIST:      in, required, type=objref
;                   Convert the energy variable to velocity.
;-
PRO MrVar_Dist2D_Prep_E2V, oDist, mass
	Compile_Opt idl2
	On_Error, 2
	
	m = Size(mass, /TNAME) EQ 'STRING' ? MrConstants(mass) : mass

;-------------------------------------------
; Convert Energy to Velocity ///////////////
;-------------------------------------------
	;
	; Calculation
	;   - E = 1/2 m v^2
	;   - v = Sqrt( 2 E / m )
	;

	;Get energy
	oEnergy = MrVar_Get(oDist['DEPEND_1'])
	
	;Conversion factor of energy units to SI (Joules)
	IF oEnergy -> HasAttr('SI_CONVERSION') THEN BEGIN
		parts = StrSplit(oEnergy['SI_CONVERSION'], '>', /EXTRACT)
		EtoSI = Float(parts[0])
	
	;Guess conversion factor based on units
	ENDIF ELSE BEGIN
		units =  oEnergy -> HasAttr('UNITS') ? oEnergy['UNITS'] : ''
		
		CASE oEnergy['UNITS'] OF
			'eV': EtoSI = MrConstants('eV2J')
			'J':  EtoSI = 1.0
			ELSE: BEGIN
				MrPrintF, 'LogWarn', 'Unknown energy units: "' + units + '". Assuming eV'
				EtoSI = MrConstants('eV2J')
			ENDCASE
		ENDCASE
	ENDELSE

	;Convert to velocity
	velocity = 1e-3 * Sqrt( 2.0 * EtoSI * oEnergy['DATA'] / m )
	oV       = MrVariable(velocity)
	oV['UNITS']         = 'km/s'
	oV['SI_CONVERSION'] = '1e3>m/s'
	
	;Store as dependent variable
	oDist['DEPEND_1'] = oV

;-------------------------------------------
; Velocity Deltas //////////////////////////
;-------------------------------------------
	;
	; Delta computation:
	;   - E = 1/2 m v^2
	;   - dE = m v dv
	;   - dv = dE / (m v)
	;
	; Assume DELTA_(PLUS|MINUS) have the same units as Energy.
	;
	
	;DELTA_PLUS_VAR
	IF oEnergy -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
		;Extract deltas and compute velocity
		;   - Convert from meters to kilometers (J = kg m^2 / s^2)
		oDelta  = MrVar_Get(oEnergy['DELTA_PLUS_VAR'])
		oV_Plus = 1e-6 * EtoSI * oDelta / (m * velocity)
		
		;Set attributes
		oV_Plus['TITLE']         = 'dV+!C(km/s)'
		oV_Plus['UNITS']         = 'km/s'
		oV_Plus['SI_CONVERSION'] = '1e3>m/s'
		oV['DELTA_PLUS_VAR']     = oV_Plus
	
	;DELTA_PLUS
	ENDIF ELSE IF oEnergy -> HasAttr('DELTA_PLUS') THEN BEGIN
		;Calculate delta & set attribute
		;   - Convert meters to kilometers (J = kg m^2 / s^2)
		oV_Plus = 1e-6 * oEnergy['DELTA_PLUS'] * EtoSI / (m * velocity)
		oV['DELTA_PLUS'] = oV_Plus
	
	;Compute deltas
	ENDIF ELSE BEGIN
		MrPrintF, 'LogWarn', 'More consideration of dV/V = constant required.'
	ENDELSE
	
	;DELTA_MINUS_VAR
	IF oEnergy -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
		;Extract deltas and compute velocity
		;   - Convert from meters to kilometers
		oDelta  = MrVar_Get(oEnergy['DELTA_PLUS_VAR'])
		oV_Plus = 1e-6 * EtoSI * oDelta / (m * velocity)
		
		;Set attributes
		oV_Plus['TITLE']         = 'dV-!C(km/s)'
		oV_Plus['UNITS']         = 'km/s'
		oV_Plus['SI_CONVERSION'] = '1e3>m/s'
		oV['DELTA_MINUS_VAR']    = oV_Plus
	
	;DELTA_MINUS
	ENDIF ELSE IF oEnergy -> HasAttr('DELTA_PLUS') THEN BEGIN
		;Calculate delta & set attribute
		;   - Convert meters to kilometers (J = kg m^2 / s^2)
		oV_Plus = 1e-6 * EtoSI * oEnergy['DELTA_PLUS'] / (m * velocity)
		oV['DELTA_PLUS'] = oV_Plus
	
	;Compute deltas
	ENDIF ELSE BEGIN
		MrPrintF, 'LogWarn', 'More consideration of dV/V = constant required.'
	ENDELSE
END


;+
;   Convert degrees to radians.
;
; :Params:
;       ODIST:      in, required, type=string/integer/objref
;                   Convert the angle variable from degrees to radians.
;-
PRO MrVar_Dist2D_Prep_Deg2Rad, oDist
	Compile_Opt idl2
	On_Error, 2
	
	oDep1   = MrVar_Get(oDist['DEPEND_2'])
	deg2rad = !pi / 180.0
	
;-------------------------------------------
; Convert Angle to Radians /////////////////
;-------------------------------------------
	
	;UNITS
	tf_deg2rad = 0B
	IF oDep1 -> HasAttr('UNITS') THEN BEGIN
		units = oDep1['UNITS']
		
		;DEGREES
		IF StRegEx(units, 'deg', /FOLD_CASE, /BOOLEAN) THEN BEGIN
			;Convert to radians
			tf_deg2rad = 1B
			oPhi       = oDep1 * deg2rad
			
			;Attributes
			oDep1    -> CopyAttrTo, oPhi
			oPhi['UNITS'] = 'radians'
			
			;Replace DEPEND_1
			oDist['DEPEND_1'] = oPhi
		
		;Not RADIANS
		ENDIF ELSE IF ~StRegEx(units, 'rad', /FOLD_CASE, /BOOLEAN) THEN BEGIN
			MrPrintF, 'LogWarn', 'Unknown angle units: "' + units + '". Assuming radians.'
			oPhi = oDep1
		ENDIF
	
	;No units
	ENDIF ELSE BEGIN
		MrPrintF, 'LogWarn', 'No angle units provided. Assuming degrees.'
		tf_deg2rad    = 1B
		oPhi          = oDep1 * deg2rad
		oDep1        -> CopyAttrTo, oPhi
		oPhi['UNITS'] = 'radians'
	ENDELSE

;-------------------------------------------
; Phi Deltas ///////////////////////////////
;-------------------------------------------
	;Convert to radians
	IF tf_deg2rad THEN BEGIN
		
		;DELTA_PLUS_VAR
		IF oPhi -> HasAttr('DELTA_PLUS_VAR') THEN BEGIN
			;Convert to degrees
			oDelta = MrVar_Get(oPhi['DELTA_PLUS_VAR'])
			oPlus  = oDelta * deg2rad
			
			;Set attributes
			oDelta -> CopyAttrTo, oPlus
			oPlus['UNITS']         = 'radians'
			oPhi['DELTA_PLUS_VAR'] = oPlus
		
		;DELTA_PLUS
		ENDIF ELSE IF oPhi -> HasAttr('DELTA_PLUS') THEN BEGIN
			oPhi['DELTA_PLUS'] = oPhi['DELTA_PLUS'] * deg2rad
		ENDIF
		
		;DELTA_MINUS_VAR
		IF oPhi -> HasAttr('DELTA_MINUS_VAR') THEN BEGIN
			;Convert to degrees
			oDelta = MrVar_Get(oPhi['DELTA_MINUS_VAR'])
			oMinus = oDelta * deg2rad
			
			;Set attributes
			oDelta -> CopyAttrTo, oMinus
			oMinus['UNITS']         = 'radians'
			oPhi['DELTA_MINUS_VAR'] = oMinus
		
		;DELTA_MINUS
		ENDIF ELSE IF oPhi -> HasAttr('DELTA_MINUS') THEN BEGIN
			oPhi['DELTA_MINUS'] = oPhi['DELTA_MINUS'] * deg2rad
		ENDIF
	ENDIF
END


;+
;   Convert a distribution from [phi,energy] to [energy,phi]. If the input distribution
;   is also time-dependent, a single distribution is extracted.
;
; :Params:
;       ODISTIN:    in, required, type=objref
;                   A MrTimeSeries object containing a distribution function
;                       to be displayed. Dimensions should be ordered as
;                       [time, azimuth, energy], with the time dimension being optional.
;                       Also must have corresponding DEPEND_[0-2] attributes.
;
; :Returns:
;       oDist:      out, required, type=object
;                   A MrTimeSeries object containing a distribution function
;                       to be displayed. Dimensions should be ordered as
;                       [time, energy, azimuth]
;-
FUNCTION MrVar_Dist2D_Prep_Transpose, oDistIn
	Compile_Opt idl2
	On_Error, 2
	
	;Transpose Data
	oDist = oDistIn -> Copy()
	oDist -> SetData, Transpose(oDist['DATA'], [0,2,1])
	
	;Swap dependent variables
	oPhi    = oDist['DEPEND_1']
	oEnergy = oDist['DEPEND_2']
	oDist['DEPEND_1'] = oEnergy
	oDist['DEPEND_2'] = oPhi
	
	;Done
	RETURN, oDist
END


;+
;   Turn a 2D [azimuth,energy] distribution into a 2D [velocity,azimuth] distribution
;   with attributes set such that it can be plotted directly with MrVar_Image.
;
;   Calling Sequence:
;       img = MrVar_Dist_Image(theDist, mass)
;       img = MrVar_Dist_Image(theDist, mass, idx)
;
; :Params:
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
; :Keywords:
;       CACHE:      in, optional, type=boolean, default=0
;                   If set, output distribution will be added to the variable cache.
;       ENERGY:     in, optional, type=boolean, default=0
;                   If set, the output distribution will be a function of energy.
;                       The default is to convert energy to velocity.
;       NAME:       in, optional, type=boolean, default=`THEDIST`.Name + '_vspace'
;                   Name to be given to the output variable.
;       POLAR:      in, optional, type=boolean, default=1
;                   If set, the distribution will be configure to be displayed in polar
;                       coordinates. This is the default. Set explicitly to 0 to display
;                       an energy-angle spectrogram.
;       _REF_EXTRA: in, optional, type=any
;                   Any keyword accepted by MrVar_Image.pro
;
; :Returns:
;       oDist:      out, required, type=object
;                   A MrVariable object containing the 2D velocity-space distribution.
;-
FUNCTION MrVar_Dist2D_Prep, theDist, mass, idx, $
CACHE=cache, $
ENERGY=energy, $
NAME=name, $
POLAR=polar
	Compile_Opt strictarr

	;Catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	;Constants & Defaults
	tf_cache    = Keyword_Set(cache)
	tf_velocity = ~Keyword_Set(energy)
	tf_polar    = N_Elements(polar) EQ 0 ? 1B : Keyword_Set(polar)
	
	;Grab the distribution. Pick a name
	oDistIn = MrVar_Get(theDist)
	IF N_Elements(name) EQ 0 THEN name = oDistIn.name + '_vspace'

;-------------------------------------------
; Make Necessary Conversions ///////////////
;-------------------------------------------
	
	;Pick the desired distribution
;	oDist = MrVar_Dist2D_Prep_Pick(oDistIn, idx)
	oDist = MrVar_Dist2D_Prep_Transpose(oDistIn)
	
	;Convert E to V
	IF tf_velocity THEN MrVar_Dist2D_Prep_E2V, oDist, mass

	;Convert Degrees to Radians
	MrVar_Dist2D_Prep_Deg2Rad, oDist

;-------------------------------------------
; Create the Image /////////////////////////
;-------------------------------------------
	
	oDist -> SetName, name
	IF tf_cache THEN oDist -> Cache
	
	;Set Attributes
	oDist['POLAR'] = tf_polar
	
	;Done
	RETURN, oDist
END