; docformat = 'rst'
;
; NAME:
;       MrMMS_FEEPS_Load_Data
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
;   Generate a plot of EPD-EIS ElectronEnergy quantities:
;       1. FGM Bxyz
;       2. Integraged flux for Telescopes 0-5
;       3. Energy-time spectrogram for T0
;       4. Energy-time spectrogram for T1
;       5. Energy-time spectrogram for T2
;       6. Energy-time spectrogram for T3
;       7. Energy-time spectrogram for T4
;       8. Energy-time spectrogram for T5
;
; :Categories:
;   MMS
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
;       2017/03/25  -   Written by Matthew Argall
;       2018/02/02  -   Sun contamination CSV files are time-dependent. Updated read
;                           routine. Set of bad sectors has also been updated. - MRA
;-
;*****************************************************************************************
;+
;
;-
FUNCTION MrMMS_FEEPS_Load_Data_Active_Eyes, mode, species
	Compile_Opt idl2
	On_Error, 2
	
	;Create a hash table
	active = hash()
	
	;BRST mode
	IF mode EQ 'brst' THEN BEGIN
		IF species EQ 'electron' THEN BEGIN
			active['mms1-top']    = [1, 2, 3, 4, 5, 9, 10, 11, 12]
			active['mms1-bottom'] = [1, 2, 3, 4, 5, 9, 10, 11, 12]
			active['mms2-top']    = [1, 2, 3, 4, 5, 9, 10, 11, 12]
			active['mms2-bottom'] = [1, 2, 3, 4, 5, 9, 10, 11, 12]
			active['mms3-top']    = [1, 2, 3, 4, 5, 9, 10, 11, 12]
			active['mms3-bottom'] = [1, 2, 3, 4, 5, 9, 10, 11, 12]
			active['mms4-top']    = [1, 2, 3, 4, 5, 9, 10, 11, 12]
			active['mms4-bottom'] = [1, 2, 3, 4, 5, 9, 10, 11, 12]
		ENDIF ELSE BEGIN
			active['mms1-top']    = [6, 7, 8]
			active['mms1-bottom'] = [6, 7, 8]
			active['mms2-top']    = [6, 7, 8]
			active['mms2-bottom'] = [6, 7, 8]
			active['mms3-top']    = [6, 7, 8]
			active['mms3-bottom'] = [6, 7, 8]
			active['mms4-top']    = [6, 7, 8]
			active['mms4-bottom'] = [6, 7, 8]
		ENDELSE
	
	;SRVY mode
	;   - Which eyes are active switched during the mission
	ENDIF ELSE BEGIN
		;Beginning of data interval & switch dates
		tswitch = MrCDF_Compute_TT2000(2017, 08, 16)
		t0      = (MrVar_GetTRange('TT2000'))[0]
		
		;BEFORE 2017-08-16
		IF t0 LT tswitch THEN BEGIN
			IF species EQ 'electron' THEN BEGIN
				active['mms1-top']    = [3, 4, 5, 11, 12]
				active['mms1-bottom'] = [3, 4, 5, 11, 12]
				active['mms2-top']    = [3, 4, 5, 11, 12]
				active['mms2-bottom'] = [3, 4, 5, 11, 12]
				active['mms3-top']    = [3, 4, 5, 11, 12]
				active['mms3-bottom'] = [3, 4, 5, 11, 12]
				active['mms4-top']    = [3, 4, 5, 11, 12]
				active['mms4-bottom'] = [3, 4, 5, 11, 12]
			ENDIF ELSE BEGIN
				active['mms1-top']    = [6, 7, 8]
				active['mms1-bottom'] = [6, 7, 8]
				active['mms2-top']    = [6, 7, 8]
				active['mms2-bottom'] = [6, 7, 8]
				active['mms3-top']    = [6, 7, 8]
				active['mms3-bottom'] = [6, 7, 8]
				active['mms4-top']    = [6, 7, 8]
				active['mms4-bottom'] = [6, 7, 8]
			ENDELSE
			
		;AFTER 2017-08-16
		ENDIF ELSE BEGIN
			IF species EQ 'electron' THEN BEGIN
				active['mms1-top']    = [3, 5, 9, 10, 12]
				active['mms1-bottom'] = [2, 4, 5, 9, 10]
				active['mms2-top']    = [1, 2, 3, 5, 10, 11]
				active['mms2-bottom'] = [1, 4, 5, 9, 11]
				active['mms3-top']    = [3, 5, 9, 10, 12]
				active['mms3-bottom'] = [1, 2, 3, 9, 10]
				active['mms4-top']    = [3, 4, 5, 9, 10, 11]
				active['mms4-bottom'] = [3, 5, 9, 10, 12]
			ENDIF ELSE BEGIN
				active['mms1-top']    = [6, 7, 8]
				active['mms1-bottom'] = [6, 7, 8]
				active['mms2-top']    = [6, 8]
				active['mms2-bottom'] = [6, 7, 8]
				active['mms3-top']    = [6, 7, 8]
				active['mms3-bottom'] = [6, 7, 8]
				active['mms4-top']    = [6, 8]
				active['mms4-bottom'] = [6, 7, 8]
			ENDELSE
		ENDELSE
	ENDELSE
	
	RETURN, active
END


;+
;       Apply flat field correction factors to FEEPS ion/electron data;
;       correct factors are from the gain factor found in:
;       
;           FlatFieldResults_V3.xlsx
;           
;       from Drew Turner, 1/19/2017
;
; NOTES:
; 
;   From Drew Turner, 1/18/17:
;       Here are the correction factors that we need to apply to the current 
;       ION counts/rates/fluxes in the CDF files.  
;       NOTE, THIS IS A DIFFERENT TYPE OF CORRECTION THAN THAT FOR THE ELECTRONS!  
;       These shifts should be applied to the counts/rates/fluxes data EYE-BY-EYE on each spacecraft.  
;       These are multiplication factors (i.e., Jnew = Jold * Gcorr). 
;       For those equations, Jold is the original count/rate/flux array and
;       Jnew is the corrected version of the arrays using the factors listed below.
;
; NOTES:
;     BAD EYES are replaced by NaNs
;
;     See the spedas distribution
;         idl/projects/mms/feeps/mms_feeps_flat_field_corrections
;-
PRO MrMMS_FEEPS_Load_Data_Bad_Sectors, theVar
	Compile_Opt idl2
	On_Error, 2
	
	;Which spacecraft and eye?
	parts   = StrSplit(theVar.name, '_', /EXTRACT)
	sc      = parts[0]
	instr   = parts[1]
	instrid = parts[2]
	mode    = parts[3]
	level   = parts[4]
	species = parts[5]
	eye     = parts[6]
	units   = parts[7]
	sensor  = parts[8]
	sid     = parts[9]
	
	;Get the sector mask for this variable
	mask_vname = StrJoin([sc, instr, instrid, mode, level, species, eye, 'sector', 'mask', sensor, sid], '_')
	oMask      = MrVar_Get(mask_vname)
	
	;Mask the data
	iBad = oMask -> Where(1, /EQUAL, COUNT=nBad)
	IF nBad GT 0 THEN theVar[iBad] = !Values.F_NaN
END


;+
;       Apply flat field correction factors to FEEPS ion/electron data;
;       correct factors are from the gain factor found in:
;       
;           FlatFieldResults_V3.xlsx
;           
;       from Drew Turner, 1/19/2017
;
; NOTES:
; 
;   From Drew Turner, 1/18/17:
;       Here are the correction factors that we need to apply to the current 
;       ION counts/rates/fluxes in the CDF files.  
;       NOTE, THIS IS A DIFFERENT TYPE OF CORRECTION THAN THAT FOR THE ELECTRONS!  
;       These shifts should be applied to the counts/rates/fluxes data EYE-BY-EYE on each spacecraft.  
;       These are multiplication factors (i.e., Jnew = Jold * Gcorr). 
;       For those equations, Jold is the original count/rate/flux array and
;       Jnew is the corrected version of the arrays using the factors listed below.
;
; NOTES:
;     BAD EYES are replaced by NaNs
;
;     See the spedas distribution
;         idl/projects/mms/feeps/mms_feeps_flat_field_corrections
;-
PRO MrMMS_FEEPS_Load_Data_Bad_Eyes, theVar
	Compile_Opt idl2
	On_Error, 2
	
	;Which spacecraft and eye?
	parts   = StRegEx(theVar.name, '^(mms[1-4]).*(electron|ion).*(top|bottom).*sensorid_([1-9]|1[0-2])', /SUBEXP, /EXTRACT)
	sc      = parts[1]
	species = parts[2]
	eye     = parts[3]
	sensor  = Fix(parts[4])
	
	;
	; Eyes that are 100% bad
	;
	
	ebad                = hash()
	ebad['mms1-top']    = [1]
	ebad['mms1-bottom'] = [1, 11]
	ebad['mms2-top']    = [5, 12]
;	ebad['mms2-bottom'] = []
	ebad['mms3-top']    = [2, 12]
	ebad['mms3-bottom'] = [2, 5, 11]
	ebad['mms4-top']    = [1, 2]
	ebad['mms4-bottom'] = [2, 4, 5, 10, 11]
	
	ibad                = hash()
;	ibad['mms1-top']    = []
;	ibad['mms1-bottom'] = []
	ibad['mms2-top']    = [7]
	ibad['mms2-bottom'] = [7]
;	ibad['mms3-top']    = []
;	ibad['mms3-bottom'] = []
	ibad['mms4-top']    = [7]
;	ibad['mms4-bottom'] = []
	
	;Get the energy table
	key = sc + '-' + eye
	
	;Ions
	IF species EQ 'ion' THEN BEGIN
		IF ibad -> HasKey(key) THEN BEGIN
			IF ~Array_Equal(ibad[key] EQ sensor, 0) $
				THEN theVar[*] = !Values.F_NaN
		ENDIF
	
	;Electrons
	ENDIF ELSE BEGIN
		IF ebad -> HasKey(key) THEN BEGIN
			IF ~Array_Equal(ebad[key] EQ sensor, 0) $
				THEN theVar[*] = !Values.F_NaN
		ENDIF
	ENDELSE
	
	
	;
	; Eyes for which the first energy channel is bad
	;
	
	ebad                = hash()
	ebad['mms1-top']    = [2, 5]
	ebad['mms1-bottom'] = [2, 3, 4, 5, 9, 11, 12] ;11
	ebad['mms2-top']    = [1, 2, 3, 4, 9, 10, 11, 12]
	ebad['mms2-bottom'] = [1, 2, 3, 4, 5, 9, 10, 11, 12]
	ebad['mms3-top']    = [4, 5, 9, 10, 11]
	ebad['mms3-bottom'] = [1, 3, 4, 9, 10, 11, 12]
	ebad['mms4-top']    = [3, 4, 5, 9, 10, 11, 12]
	ebad['mms4-bottom'] = [1, 3, 9, 12]
	
	ibad                = hash()
	ibad['mms1-top']    = [6]
	ibad['mms1-bottom'] = [7, 8]
	ibad['mms2-top']    = [8] ;8
	ibad['mms2-bottom'] = [6, 8, 12] ;12
	ibad['mms3-top']    = [2, 6, 7] ;2
	ibad['mms3-bottom'] = [6, 7]
;	ibad['mms4-top']    = []
	ibad['mms4-bottom'] = [6, 7] ;6
	
	;Ions
	IF species EQ 'ion' THEN BEGIN
		IF ibad -> HasKey(key) THEN BEGIN
			IF ~Array_Equal(ibad[key] EQ sensor, 0) $
				THEN theVar[*,0] = !Values.F_NaN
		ENDIF
	
	;Electrons
	ENDIF ELSE BEGIN
		IF ebad -> HasKey(key) THEN BEGIN
			IF ~Array_Equal(ebad[key] EQ sensor, 0) $
				THEN theVar[*,0] = !Values.F_NaN
		ENDIF
	ENDELSE
	
	
	;
	; Eyes for which the second energy channel is bad
	;
	
	ebad                = hash()
;	ebad['mms1-top']    = []
;	ebad['mms1-bottom'] = []
;	ebad['mms2-top']    = []
	ebad['mms2-bottom'] = [12]
	ebad['mms3-top']    = [1]
;	ebad['mms3-bottom'] = []
;	ebad['mms4-top']    = []
;	ebad['mms4-bottom'] = []
	
	ibad                = hash()
	ibad['mms1-top']    = [6, 7]
	ibad['mms1-bottom'] = [6, 7, 8]
	ibad['mms2-top']    = [8]
	ibad['mms2-bottom'] = [6, 8]
	ibad['mms3-top']    = [6, 7]
	ibad['mms3-bottom'] = [6, 7]
	ibad['mms4-top']    = [6]
	ibad['mms4-bot']    = [7]
	
	;Ions
	IF species EQ 'ion' THEN BEGIN
		IF ibad -> HasKey(key) THEN BEGIN
			IF ~Array_Equal(ibad[key] EQ sensor, 0) $
				THEN theVar[*,1] = !Values.F_NaN
		ENDIF
	
	;Electrons
	ENDIF ELSE BEGIN
		IF ebad -> HasKey(key) THEN BEGIN
			IF ~Array_Equal(ebad[key] EQ sensor, 0) $
				THEN theVar[*,1] = !Values.F_NaN
		ENDIF
	ENDELSE
END


;+
;       Apply flat field correction factors to FEEPS ion/electron data;
;       correct factors are from the gain factor found in:
;       
;           FlatFieldResults_V3.xlsx
;           
;       from Drew Turner, 1/19/2017
;
; NOTES:
; 
;   From Drew Turner, 1/18/17:
;       Here are the correction factors that we need to apply to the current 
;       ION counts/rates/fluxes in the CDF files.  
;       NOTE, THIS IS A DIFFERENT TYPE OF CORRECTION THAN THAT FOR THE ELECTRONS!  
;       These shifts should be applied to the counts/rates/fluxes data EYE-BY-EYE on each spacecraft.  
;       These are multiplication factors (i.e., Jnew = Jold * Gcorr). 
;       For those equations, Jold is the original count/rate/flux array and
;       Jnew is the corrected version of the arrays using the factors listed below.
;
; NOTES:
;
;     See the spedas distribution
;         idl/projects/mms/feeps/mms_feeps_flat_field_corrections
;-
FUNCTION MrMMS_FEEPS_Load_Data_Flat_Field, theVar
	Compile_Opt idl2
	On_Error, 2
	
	G_corr                  = hash()
	G_corr['mms1-top-6']    = 0.7
	G_corr['mms1-top-7']    = 2.5
	G_corr['mms1-top-8']    = 1.5
	G_corr['mms1-bottom-5'] = 1.2    ; updated 1/24
	G_corr['mms1-bottom-6'] = 0.9
	G_corr['mms1-bottom-7'] = 2.2    ; updated 1/24
	G_corr['mms1-bottom-8'] = 1.0

	G_corr['mms2-top-4']    = 1.2    ; added 1/24
	G_corr['mms2-top-6']    = 1.3
	G_corr['mms2-top-7']    = 0      ; bad eye
	G_corr['mms2-top-8']    = 0.8
	G_corr['mms2-bottom-6'] = 1.4
	G_corr['mms2-bottom-7'] = 0      ; bad eye
	G_corr['mms2-bottom-8'] = 1.5

	G_corr['mms3-top-6']    = 0.7
	G_corr['mms3-top-7']    = 0.8
	G_corr['mms3-top-8']    = 1.0
	G_corr['mms3-bottom-6'] = 0.9
	G_corr['mms3-bottom-7'] = 0.9
	G_corr['mms3-bottom-8'] = 1.3

	G_corr['mms4-top-6']    = 0.8
	G_corr['mms4-top-7']    = 0      ; bad eye
	G_corr['mms4-top-8']    = 1.0
	G_corr['mms4-bottom-6'] = 0.8
	G_corr['mms4-bottom-7'] = 0.6
	G_corr['mms4-bottom-8'] = 0.9
	G_corr['mms4-bottom-9'] = 1.5    ; added 1/24
	
	;Which spacecraft and eye?
	parts  = StRegEx(theVar.name, '^(mms[1-4]).*(top|bottom).*sensorid_([1-9]|1[0-2])', /SUBEXP, /EXTRACT)
	sc     = parts[1]
	eye    = parts[2]
	sensor = parts[3]
	
	;Correct the counts
	tag = sc + '-' + eye + '-' + sensor
	IF G_corr -> HasTag(tag) $
		THEN theVar -> SetData, theVar['DATA'] * G_corr[tag]
END


;+
;
;       This function returns the energy table based on
;       each spacecraft and eye; based on the table from:
;       
;               FlatFieldResults_V3.xlsx
;               
;       from Drew Turner, 1/19/2017
;
; NOTES:
;     BAD EYES are replaced by NaNs
;
;     See the spedas distribution
;         idl/projects/mms/feeps/mms_feeps_correct_energies.pro
;         idl/projects/mms/feeps/mms_feeps_energy_table.pro
;-
PRO MrMMS_FEEPS_Load_Data_Fix_Energies, theVar
	Compile_Opt idl2
	On_Error, 2
	
	;For easy referencing
	NaN = !values.f_nan
	
	;Create the correction table
	table                = hash()
	table['mms1-top']    = [14.0,  7.0, 16.0, 14.0, 14.0,  0.0, 0.0,  0.0, 14.0, 14.0,  17.0, 15.0]
	table['mms1-bottom'] = [ NaN, 14.0, 14.0, 13.0, 14.0,  0.0, 0.0,  0.0, 14.0, 14.0, -25.0, 14.0]
	table['mms2-top']    = [-1.0,  6.0, -2.0, -1.0,  NaN,  0.0, NaN,  0.0,  4.0, -1.0,  -1.0,  0.0]
	table['mms2-bottom'] = [-2.0, -1.0, -2.0,  0.0, -2.0, 15.0, NaN, 15.0, -1.0, -2.0,  -1.0, -3.0]
	table['mms3-top']    = [-3.0,  NaN,  2.0, -1.0, -5.0,  0.0, 0.0,  0.0, -3.0, -1.0,  -3.0,  NaN]
	table['mms3-bottom'] = [-7.0,  NaN, -5.0, -6.0,  NaN,  0.0, 0.0, 12.0,  0.0, -2.0,  -3.0, -3.0]
	table['mms4-top']    = [ NaN,  NaN, -2.0, -5.0, -5.0,  0.0, NaN,  0.0, -1.0, -3.0,  -6.0, -6.0]
	table['mms4-bottom'] = [-8.0,  NaN, -2.0,  NaN,  NaN, -8.0, 0.0,  0.0, -2.0,  NaN,   NaN, -4.0]
	
	;Which spacecraft and eye?
	parts  = StRegEx(theVar.name, '^(mms[1-4]).*(top|bottom).*sensorid_([1-9]|1[0-2])', /SUBEXP, /EXTRACT)
	sc     = parts[1]
	eye    = parts[2]
	sensor = Fix(parts[3])
	
	;Correct the energy table
	;   - Sensors start at 1, IDL index starts at 0
	oE_old = theVar['DEPEND_1']
	oE_new = oE_old + table[sc+'-'+eye, sensor-1]
	oE_old -> CopyAttrTo, oE_new
	theVar['DEPEND_1'] = oE_new
END


;+
;    this procedure splits the last integral channel from the FEEPS spectra, 
;    creating 2 new tplot variables:
;    
;       [original variable]_clean - spectra with the integral channel removed
;       [original variable]_500keV_int - the integral channel that was removed
;-
FUNCTION MrMMS_FEEPS_Load_Data_Integral_Channel, theVar
	Compile_Opt idl2
	On_Error, 2
	
	;Name of integral channel
	;   - Put "500keV" after eye and before data product
	parts     = StrSplit(theVar.name, '_', /EXTRACT)
	int_vname = StrJoin( [parts[0:7], '500keV', parts[8:*]], '_' )
	
	;Create variables
	o500keV  = theVar[*, -1]
	
	;Remove the 500keV channel
	;   - Make sure to output the new variable via THEVAR
	newVar  = theVar[*,0:-2]
	newVar -> SetName, theVar.name
	MrVar_Replace, theVar, newVar
	theVar  = newVar
	
	;Add 500keV to the cache
	o500keV -> SetName, int_vname
	o500keV -> Cache
	
	;Attributes
	o500keV -> RemoveAttr, 'DEPEND_1'
	
	;Return its variable name
	RETURN, o500keV
END


;+
;    this procedure splits the last integral channel from the FEEPS spectra, 
;    creating 2 new tplot variables:
;    
;       [original variable]_clean - spectra with the integral channel removed
;       [original variable]_500keV_int - the integral channel that was removed
;-
FUNCTION MrMMS_FEEPS_Load_Data_Omni, varname
	Compile_Opt idl2
	On_Error, 2
	
	;Energy and gain correction factors
	;   - inter-spacecraft calibrations
	;   - The element corresponds to the spacecraft
	eEcorr = [14.0,  -1.0, -3.0, -3.0]
	iEcorr = [ 0.0,   0.0,  0.0,  0.0]
	eGfact = [ 1.0,   1.0,  1.0,  1.0]
	iGfact = [ 0.84,  1.0,  1.0,  1.0]
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	
	;Dissect the variable names
	;   - The eye and sensor have been replaced with the "#" character
	;   - The sc, instr, mode, and level occur before the first "#"
	segs    = StrSplit(varname[0], '*', /EXTRACT, COUNT=nSegs)
	parts   = [StrSplit(segs[0], '_', /EXTRACT), StrSplit(segs[1], '_', /EXTRACT)]
	sc      = parts[0]
	instr   = parts[1] + '_' + parts[2]
	mode    = parts[3]
	level   = parts[4]
	species = parts[5]
	units   = parts[6]
	sensor  = parts[7]
	
	;
	; Available sensors
	;
	
	sensors = MrMMS_FEEPS_Load_Data_Active_Eyes(mode, species)
	
	;Variable names of each sensor
	top_vnames    = StrJoin( [parts[0:5], 'top',    units,           sensor], '_') + '_' + String(sensors[sc+'-top'],    FORMAT='(i0)')
	top500_vnames = StrJoin( [parts[0:5], 'top',    units, '500keV', sensor], '_') + '_' + String(sensors[sc+'-top'],    FORMAT='(i0)')
	bot_vnames    = StrJoin( [parts[0:5], 'bottom', units,           sensor], '_') + '_' + String(sensors[sc+'-bottom'], FORMAT='(i0)')
	bot500_vnames = StrJoin( [parts[0:5], 'bottom', units, '500keV', sensor], '_') + '_' + String(sensors[sc+'-bottom'], FORMAT='(i0)')
	IF nSegs GT 2 THEN BEGIN
		top_vnames += segs[2:*]
		bot_vnames += segs[2:*]
	ENDIF
	
	;Ouput name
	omni_vname    = StrJoin( [parts[0:5], units, 'omni'          ], '_')
	omni500_vname = StrJoin( [parts[0:5], units, 'omni', '500keV'], '_')
	
	;Are the names present?
	IF ~Array_Equal( MrVar_IsCached( [top_vnames, bot_vnames] ), 1) THEN BEGIN
		MrPrintF, 'LogText', 'Not all sensors found: "' + varname + '".'
		RETURN, !Null
	ENDIF

;-------------------------------------------
; Energy Bins //////////////////////////////
;-------------------------------------------
	
	;Energy table
	IF species EQ 'electron' THEN BEGIN
		energies = [ 33.200000D,  51.900000D,  70.600000d,  89.400000d, 107.10000d, $
		            125.20000D,  146.50000D,  171.30000d,  200.20000d,  234.00000d, $
		            273.40000D,  319.40000D,  373.20000d,  436.00000d,  509.20000d ]
		
		CASE sc OF
			'mms1': energies += eECorr[0]
			'mms2': energies += eECorr[1]
			'mms3': energies += eECorr[2]
			'mms4': energies += eECorr[3]
			ELSE: Message, 'Invalid spacecraft: "' + sc + '".'
		ENDCASE
		
	;Ions
	ENDIF ELSE BEGIN
		energies = [ 57.900000d,  76.800000d,  95.400000d, 114.10000d, 133.00000d, $
		            153.70000d,  177.60000d,  205.10000d,  236.70000d, 273.20000d, $
		            315.40000d,  363.80000d,  419.70000d,  484.20000d, 558.60000d]
		
		CASE sc OF
			'mms1': energies += iECorr[0]
			'mms2': energies += iECorr[1]
			'mms3': energies += iECorr[2]
			'mms4': energies += iECorr[3]
			ELSE: Message, 'Invalid spacecraft: "' + sc + '".'
		ENDCASE
	ENDELSE

;-------------------------------------------
; Omni-Directional Flux ////////////////////
;-------------------------------------------
	
	;SITL
	IF level EQ 'sitl' THEN BEGIN
		MrPrintF, 'LogErr', 'SITL data products are not yet supported.'
		RETURN, !Null
	ENDIF
	
	;Allocate memory
	oVar     = MrVar_Get(top_vnames[0])
	dims     = oVar.dimensions
	nTime    = dims[0]
	nEnergy  = dims[1]
	nTop     = N_Elements(sensors[sc+'-top'])
	nBottom  = N_Elements(sensors[sc+'-bottom'])
	nSensors = N_Elements(sensors)
	omni     = FltArr(nTime, nEnergy, nTop + nBottom) + !Values.F_NaN
	omni500  = FltArr(nTime, nTop + nBottom) + !Values.F_NaN
	
	;TOP
	FOR iSensor = 0, nTop - 1 DO BEGIN
		IF MrVar_IsCached(top_vnames[iSensor]) THEN BEGIN
			;Collect the data
			oVar               = MrVar_Get(top_vnames[iSensor])
			oVar500            = MrVar_Get(top500_vnames[iSensor])
			omni[0,0,iSensor]  = oVar['DATA']
			omni500[0,iSensor] = oVar500['DATA']
			
			;Energies
			oE   = oVar['DEPEND_1']
			iBad = Where( Abs( oE['DATA'] - Energies ) GT 0.1*Energies, nBad)
			IF nBad GT 0 THEN omni[*,iBad,iSensor] = !Values.F_NaN
		ENDIF ELSE BEGIN
			MrPrintF, 'LogWarn', 'Sensor not found: "' + top_vnames[iSensor] + '".'
		ENDELSE
	ENDFOR
	
	;BOTTOM
	FOR iSensor = 0, nBottom-1 DO BEGIN
		IF MrVar_IsCached(bot_vnames[iSensor]) THEN BEGIN
			;Collect the data
			oVar                    = MrVar_Get(bot_vnames[iSensor])
			oVar500                 = MrVar_Get(top500_vnames[iSensor])
			omni[0,0,iSensor+nTop]  = oVar['DATA']
			omni500[0,iSensor+nTop] = oVar500['DATA']
			
			;Energies
			oE   = oVar['DEPEND_1']
			iBad = Where( Abs( oE['DATA'] - Energies ) GT 0.1*Energies, nBad)
			IF nBad GT 0 THEN omni[*,iBad,iSensor+nTop] = !Values.F_NaN
		ENDIF ELSE BEGIN
			MrPrintF, 'LogWarn', 'Sensor not found: "' + bot_vnames[iSensor] + '".'
		ENDELSE
	ENDFOR
	
	;Average the data
	flux_omni    = Mean( omni,    DIMENSION=3, /NAN )
	flux_omni500 = Mean( omni500, DIMENSION=2, /NAN )
	
;-------------------------------------------
; Gain Correction //////////////////////////
;-------------------------------------------
	
	IF StRegEx(species, 'electron', /BOOLEAN) THEN BEGIN
		omni    *= eGfact[Fix(StrMid(sc, 3, 1)) - 1]
		omni500 *= eGfact[Fix(StrMid(sc, 3, 1)) - 1]
	ENDIF ELSE BEGIN
		omni    *= iGfact[Fix(StrMid(sc, 3, 1)) - 1]
		omni500 *= iGfact[Fix(StrMid(sc, 3, 1)) - 1]
	ENDELSE

;-------------------------------------------
; Output Variable //////////////////////////
;-------------------------------------------
	
	;Energy variable
	oEnergy = MrVariable(energies, NAME=e_vname, /NO_COPY)
	oEnergy['LOG']   = 1B
	oEnergy['TITLE'] = 'Energy!C(keV)'
	oEnergy['UNITS'] = 'keV'
	
	;Omni
	oOmni = MrTimeSeries( oVar['TIMEVAR'], flux_omni, $
	                      /CACHE, $
	                      NAME = omni_vname )
	
	;Attributes
	oVar -> CopyAttrTo, oOmni, [ 'FILLVAL', 'FORMAT', 'LOG', 'SCALETYP', 'SCALEMAX', $
	                             'SCALEMIN', 'SI_CONVERSION', 'UNITS', 'VALIDMIN', 'VALIDMAX']
	
;	range = [oOmni.min, oOmni.max]
;	IF oOmni['LOG'] && range[0] LE 0 THEN range[0] = 1e-1
;	oOmni['AXIS_RANGE'] = range
	oOmni['DEPEND_1']   = oEnergy
	oOmni['NAN']        = 1B
	oOmni['TITLE']      = 'Omni!C' + oOmni['UNITS']
	oOmni['PLOT_TITLE'] = 'FEEPS Omni-Directional Flux'
	
	;Omni-500
	oOmni500 = MrScalarTS( oVar['TIMEVAR'], flux_omni500, $
	                       /CACHE, $
	                       NAME = omni500_vname )
	
	;Attributes
	oVar -> CopyAttrTo, oOmni500, [ 'FILLVAL', 'FORMAT', 'LOG', 'SCALETYP', 'SCALEMAX', $
	                                'SCALEMIN', 'SI_CONVERSION', 'UNITS', 'VALIDMIN', 'VALIDMAX']
	oOmni500['TITLE']      = 'Omni >500keV!C' + oOmni500['UNITS']
	oOmni500['PLOT_TITLE'] = 'FEEPS Integrated Omni-Directional Flux >500keV'
END



;+
;    Create pitch angle distributions of FEEPS ion and electron data.
;
; :Params:
;       VARNAME:        in, required, type=string
;                       Name of a MrTimeSeries variable for which a pitch angle
;                           distribution is made.
;       OPA:            in, required, type=objref
;                       A MrVariable containing the pitch angles of each sector.
;
; :Keywords:
;       DPA:            in, optional, type=float
;                       Width of the pitch angle bins in the output distribution. If
;                           neither DPA nor `NPA_BINS` are given, the default is 60.0/22.0.
;                           If `NPA_BINS` is given, DPA is calculated from it and `PA_RANGE`.
;       E_RANGE:        in, optional, type=fltarr(2), default=[70, 600]
;                       Energy range (keV) over which the PAD is created.
;       NPA_BINS:       in, optional, type=float
;                       Number of pitch angle bins to create. If `DPA` is given, the
;                           default is calculated from it and `PA_RANGE`.
;       PA_RANGE:       in, optional, type=float, default=[0.0, 180.0]
;                       The pitch angle range over which the PAD is created.
;
; :Returns:
;       OPAD:           out, required, type=objref
;                       A MrTimeSeries variable containing the pitch angle distribution.
;-
FUNCTION MrMMS_FEEPS_Load_Data_PAD, varname, oPA, $
DPA=dPA, $
E_RANGE=e_range, $
NPA_BINS=nPA_bins, $
PA_RANGE=pa_range
	Compile_Opt idl2
	On_Error, 2
	
	;Defaults
	IF N_Elements(pa_range) EQ 0 THEN pa_range = [0.0, 180.0]   ;degrees
	IF N_Elements(e_range)  EQ 0 THEN e_range  = [70, 600]       ;keV
	IF N_Elements(dPA) EQ 0 && N_Elements(nPA_bins) EQ 0 THEN dPA = 360.0/22.0
	
	;Conflicts
	IF N_Elements(dPA) GT 1 && N_Elements(nPA_bins) GT 1 $
		THEN Message, 'DPA and NPA_BINS are mutually exclusive.'
	
;-------------------------------------------
; Bin Information //////////////////////////
;-------------------------------------------
	IF N_Elements(dPA) GT 0 $
		THEN nPA_bins = Long( (pa_range[1] - pa_range[0]) / dPA ) + 1 $
		ELSE dPA      = (pa_range[1] - pa_range[0]) / nPA_bins
	pa_bins = linspace(pa_range[0], pa_range[1], nPA_bins)
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	
	;Dissect the variable names
	;   - The eye and sensor have been replaced with the "#" character
	;   - The sc, instr, mode, and level occur before the first "#"
	segs    = StrSplit(varname[0], '*', /EXTRACT, COUNT=nSegs)
	parts   = [StrSplit(segs[0], '_', /EXTRACT), StrSplit(segs[1], '_', /EXTRACT)]
	sc      = parts[0]
	instr   = parts[1] + '_' + parts[2]
	mode    = parts[3]
	level   = parts[4]
	species = parts[5]
	units   = parts[6]
	sensor  = parts[7]
	
	;Available sensors
	;   - DANGRESP = Angular response (finite field of view) of instruments
	dAngResp = species EQ 'electron' ? 21.4 : 10
	sensors  = MrMMS_FEEPS_Load_Data_Active_Eyes(mode, species)
	
	;Variable names of each sensor
	top_vnames    = StrJoin( [parts[0:5], 'top',    units,           sensor], '_') + '_' + String(sensors[sc+'-top'],    FORMAT='(i0)')
	top500_vnames = StrJoin( [parts[0:5], 'top',    units, '500keV', sensor], '_') + '_' + String(sensors[sc+'-top'],    FORMAT='(i0)')
	bot_vnames    = StrJoin( [parts[0:5], 'bottom', units,           sensor], '_') + '_' + String(sensors[sc+'-bottom'], FORMAT='(i0)')
	bot500_vnames = StrJoin( [parts[0:5], 'bottom', units, '500keV', sensor], '_') + '_' + String(sensors[sc+'-bottom'], FORMAT='(i0)')
	IF nSegs GT 2 THEN BEGIN
		top_vnames += segs[2:*]
		bot_vnames += segs[2:*]
	ENDIF
	
	;Ouput name
	pa_bin_vname = StrJoin( [parts[0:5], units, 'pa',  'bins'  ], '_')
	pad_vname    = StrJoin( [parts[0:5], units, 'pad'          ], '_')
	pad500_vname = StrJoin( [parts[0:5], units, 'pad', '500keV'], '_')
	
	;Are the names present?
	IF ~Array_Equal( MrVar_IsCached( [top_vnames, bot_vnames] ), 1) THEN BEGIN
		MrPrintF, 'LogText', 'Not all sensors found: "' + varname + '".'
		RETURN, !Null
	ENDIF
	
;-------------------------------------------
; Sector Indices ///////////////////////////
;-------------------------------------------
	;Get the channel IDs
	theMode    = mode EQ 'brst' ? 'burst' : 'survey'
	chan_vname = StrJoin( [parts[0:5], theMode, species, 'channel', 'ids'], '_')
	oChan      = MrVar_Get(chan_vname)
	
	;Find the indices of the top and bottom channels
	top_sensors = 'TOP_'    + String(sensors[sc+'-top'],    FORMAT='(i0)')
	bot_sensors = 'BOTTOM_' + String(sensors[sc+'-bottom'], FORMAT='(i0)')
	tf_top = MrIsMember( StrTrim(oChan['DATA'], 2), top_sensors, A_INDICES=iTop )
	tf_bot = MrIsMember( StrTrim(oChan['DATA'], 2), bot_sensors, A_INDICES=iBot )
	
;-------------------------------------------
; Average Energy & Combine /////////////////
;-------------------------------------------
	
	;Allocate memory
	oVar     = MrVar_Get(top_vnames[0])
	dims     = oVar.dimensions
	nTime    = dims[0]
	nEnergy  = dims[1]
	nTop     = N_Elements(iTop)
	nBot     = N_Elements(iBot)
	nSensors = nTop + nBot
	
	;All indices
	iEyes = [iTop, iBot]
	eye_vnames = [top_vnames, bot_vnames]
	
	;Combine sensors into a single array
	alleyes = FltArr(nTime, nSensors); + !Values.F_NaN
	allpas  = FltArr(nTime, nSensors) + !Values.F_NaN
	FOR i = 0, nSensors - 1 DO BEGIN
		;Average over the energy range
		oVar = MrVar_Get(eye_vnames[iEyes[i]])
		oE   = oVar['DEPEND_1']
		iE   = Where(oE['DATA'] GE e_range[0] AND oE['DATA'] LE e_range[1], cE)
		
		IF cE EQ 0 THEN BEGIN
			MrPrintF, 'LogWarn', 'No energies in range.'
			CONTINUE
		ENDIF
		
		;Set zeros equal to NaN
		iZero = oVar -> Where(0.0, /EQUAL, COUNT=nZero)
		IF nZero GT 0 THEN oVar[iZero] = !Values.F_NaN
		
		;Store into single variable
		allpas[*,i] = oPA['DATA',*,iEyes[i]]
		IF cE EQ 1 $
			THEN alleyes[0,i] = oVar['DATA',*,iE] $
			ELSE alleyes[0,i] = Mean(oVar['DATA',*,iE], DIMENSION=2, /NAN)
	ENDFOR
	
;-------------------------------------------
; Bin By Pitch Angle ///////////////////////
;-------------------------------------------
	;Allocate memory
	pad      = FltArr(nTime, nPA_bins); + !Values.F_NaN
;	pad500   = FltArr(nTime, nPA_bins) + !Values.F_NaN
	
	;Select only the sensors that are active
	oPA = oPA[*,iEyes]
	
	;Bin the data
	FOR i = 0, nTime - 1 DO BEGIN
		FOR j = 0, nPA_Bins - 1 DO BEGIN
			idx = Where( oPA['DATA',i,*]-dAngResp LE pa_bins[j]+dPA AND $
			             oPA['DATA',i,*]+dAngResp GE pa_bins[j], count )

			IF count GT 0 THEN BEGIN
				IF count EQ 1 $
					THEN pad[i,j] = alleyes[i,idx] $
					ELSE pad[i,j] = Mean(alleyes[i,idx], DIMENSION=2, /NAN)
			ENDIF
		ENDFOR
	ENDFOR
	
	;Loop over all points
;	FOR i = 0, nTime - 1 DO BEGIN
;		;Locate the data within the PA bins
;		h = Histogram( oPA['DATA',i,*], $
;		               MAX             = pa_range[1], $
;		               MIN             = pa_range[0], $
;		               NBINS           = nPA_bins, $
;		               REVERSE_INDICES = ri )
;		
;		;Loop over each bin
;		FOR j = 0, N_Elements(h) - 1 DO BEGIN
;			nPA = ri[j+1] - ri[j]
;			IF nPA EQ 0 THEN CONTINUE
;			
;			;Indices within the source
;			iPA  = ri[ri[j]:ri[j+1]-1]
;			
;			;Form the PAD
;			IF nPA EQ 1 $
;				THEN pad[i,j] = alleyes[i,iPA] $
;				ELSE pad[i,j] = Mean(alleyes[i,iPA], DIMENSION=2, /NAN)
;		ENDFOR
;	ENDFOR
	
;-------------------------------------------
; Output Variable //////////////////////////
;-------------------------------------------
	;Pitch Angle Variable
	oPitchAngle = MrVariable( pa_bins, $
	                          NAME = pa_bin_vname )
	oPitchAngle['CATDESC']    = 'Pitch angle bins.'
	oPitchAngle['DELTA_PLUS'] = dPA
	oPitchAngle['LABLAXIS']   = 'PA!C(deg)
	oPitchAngle['LOG']        = 0
	oPitchAngle['UNITS']      = 'deg'
	oPitchAngle['VALIDMIN']   = 0.0
	oPitchAngle['VALIDMAX']   = 180.0
	
	;Energy variable
;	oEnergy = MrVariable(energies, NAME=e_vname, /NO_COPY)
;	oEnergy['LOG']   = 1B
;	oEnergy['TITLE'] = 'Energy!C(eV)'
;	oEnergy['UNITS'] = 'keV'
	
	;Omni
	oPAD = MrTimeSeries( oVar['TIMEVAR'], pad, $
	                     /CACHE, $
	                     NAME = pad_vname )
	
	;Attributes
	oVar -> CopyAttrTo, oPAD, [ 'FILLVAL', 'FORMAT', 'LOG', 'SCALETYP', 'SCALEMAX', $
	                            'SCALEMIN', 'SI_CONVERSION', 'UNITS', 'VALIDMIN', 'VALIDMAX']
	oPAD['DEPEND_1']   = oPitchAngle
	oPAD['NAN']        = 1B
	oPAD['SCALE']      = 1B
	oPAD['TITLE']      = 'PAD!C' + oVar['UNITS']
	oPAD['PLOT_TITLE'] = 'FEEPS Pitch Angle Distribution'
END


;+
;    Read the CSV files containing the sun masked sectors.
;
; :Returns:
;    MASK:          out, required, type=hash
;                   A mask indicating whether a sector is contaminated by sunlight
;                       or not. Keys are 'sc-eye-sensor', where SC is the spacecraft ID,
;                       EYE is either "top" or "bottom", and SENSOR is the sensor ID.
;                       A value of 1 (0) indicates that a sensor is (not) contaminated
;                       by sunlight.
;-
FUNCTION MrMMS_FEEPS_Load_Data_Read_CSV
	Compile_Opt idl2
	On_Error, 2
	
	;File to read
	path = FilePath( '', $
	                 ROOT_DIR     = File_DirName( File_Which('mrmms_feeps_load_data.pro') ), $
	                 SUBDIRECTORY = 'sun_tables' )
	
	;Find the CSV file that is closest to (without going over) the data range
	files  = File_Search(path, StrJoin(['MMS1', 'FEEPS', 'ContaminatedSectors', '*.csv'], '_'))
	times  = StRegEx(files, '([0-9]{4})([0-9]{2})([0-9]{2})', /EXTRACT, /SUBEXP)
	tt2000 = MrCDF_Epoch_Compute(Fix(times[1,*]), Fix(times[2,*]), Fix(times[3,*]))
	t0     = (MrVar_GetTRange('TT2000'))[0]
	iFile  = Value_Locate(tt2000, t0) > 0
	
	;Convert the structure to a hash
	mask = hash()
	
	;Each spacecraft
	FOR i = 0, 3 DO BEGIN
		sc   = 'mms' + String(i+1, FORMAT='(i1)')
		file = FilePath( StrJoin( [StrUpCase(sc), 'FEEPS', 'ContaminatedSectors', times[iFile]+'.csv'], '_' ), $
		                 ROOT_DIR = path )
	
		;Read the file
		data = Read_CSV(file)
		
		;TOP
		FOR j = 0, 11 DO BEGIN
			iBad = Where(data.(j) EQ 1, nBad)
			IF nBad GT 0 THEN mask[sc + '-top-' + String(j+1, FORMAT='(i0)')] = iBad
		ENDFOR
		
		;BOTTOM
		FOR j = 0, 11 DO BEGIN
			iBad = Where(data.(j+12) EQ 1, nBad)
			IF nBad GT 0 THEN mask[sc + '-bottom-' + String(j+1, FORMAT='(i0)')] = iBad
		ENDFOR
	ENDFOR
	
	RETURN, mask
END


;+
;    Remove data from sectors that have been contaminated by sunlight.
;
; :Params:
;       THEVAR:     in, required, type=objref
;                   A MrTimeSeries variable containing the sensor data.
;       SUNSECTOR:  in, required, type=objref
;                   A MrTimeSeries variable containing the sun sector data.
;       MASK:       in, required, type=hash
;                   A table of bad sectors.
;-
PRO MrMMS_FEEPS_Load_Data_Remove_Sun, theVar, sunSector, mask
	Compile_Opt idl2
	On_Error, 2
	
	;Which spacecraft, eye, and sensor?
	parts  = StRegEx(theVar.name, '^(mms[1-4]).*(top|bottom).*sensorid_([1-9]|1[0-2])', /SUBEXP, /EXTRACT)
	sc     = parts[1]
	eye    = parts[2]
	sensor = parts[3]
	
	;Remove bad sectors
	key = sc + '-' + eye + '-' + sensor
	IF mask -> HasKey(key) THEN BEGIN
		tf_bad = MrIsMember( mask[key], sunSector['DATA'], iBad, COUNT=nBad )
		IF nBad GT 0 THEN theVar[iBad,*] = !Values.F_NaN
	ENDIF
END


;+
;   Find FEEPS data then load it into the MrVariable cache.
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' }
;
; :Keywords:
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
;       TRANGE:     in, optional, type=string/strarr(2), default=MrVar_GetTRange()
;                   The start and end times of the data interval to be plotted, formatted
;                       as 'YYYY-MM-DDThh:mm:ss'
;-
PRO MrMMS_FEEPS_Load_Data, sc, mode, $
LEVEL=level, $
OPTDESC=optdesc, $
SUFFIX=suffix, $
TRANGE=trange, $
VARFORMAT=varformat, $
VARNAMES=varnames
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN
	ENDIF
	
	instr = 'epd_feeps'
	IF N_Elements(varformat) EQ 0 THEN varformat = '*intensity*'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(mode)      EQ 0 THEN mode      = 'srvy'
	IF N_Elements(optdesc)   EQ 0 THEN optdesc   = 'electron'
	IF N_Elements(suffix)    EQ 0 THEN suffix    = ''
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
	;Also get the spin sector number to make the sun correction
	varformat = [varformat, '*mask*', '*spinsectnum*', '*pitch_angle', '*channel_ids']

;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	;FEEPS
	MrMMS_Load_Data, sc, 'feeps', mode, level, $
	                 OPTDESC   = optdesc, $
	                 SUFFIX    = suffiX, $
	                 VARFORMAT = varformat, $
	                 VARNAMES  = varnames
	
	;Get the sun-sensor mask
	mask = MrMMS_FEEPS_Load_Data_Read_CSV()
	
	;Get the sun-sector variable
	iSun = Where( StrMatch(varnames, '*spinsectnum*'), nSun )
	oSunSector = MrVar_Get(varnames[iSun])
	
	;Get the pitch angles
	iPA = Where( StrMatch(varnames, '*pitch_angle'), nPA )
	oPA = MrVar_Get(varnames[iPA])
	
;-------------------------------------------
; Fix Data /////////////////////////////////
;-------------------------------------------
	;Fix each variable
	FOR i = 0, N_Elements(varnames) - 1 DO BEGIN
		tf_ion  = StRegEx(varnames[i], 'ion', /BOOLEAN)
		tf_top  = StRegEx(varnames[i], 'top', /BOOLEAN)
		tf_data = StRegEx(varnames[i], '(count_rate|intensity)', /BOOLEAN)
		
		;SKIP
		IF ~tf_data THEN CONTINUE
		
		;Extract the variable
		oVar = MrVar_Get(varnames[i])
		
		;Energy table corrections
		MrMMS_FEEPS_Load_Data_Fix_Energies, oVar
		
		;Flat-field correction
		;   - Only for ions
		IF tf_ion THEN MrMMS_FEEPS_Load_Data_Flat_Field, oVar
	
		;Remove data from bad eyes
;		MrMMS_FEEPS_Load_Data_Bad_Eyes, oVar
		MrMMS_FEEPS_Load_Data_Bad_Sectors, oVar
	
		;Separate 500keV channel
		;   - It integrates from 500keV and above
		;   - Energy bin is disproportionate
		MrMMS_FEEPS_Load_Data_Remove_Sun, oVar, oSunSector, mask
		
		;Separate 500keV channel
		;   - It integrates from 500keV and above
		;   - Energy bin is disproportionate
		o500keV = MrMMS_FEEPS_Load_Data_Integral_Channel(oVar)
		varnames = [varnames, o500keV.name]
	ENDFOR
	
;-------------------------------------------
; Omni-Directional Flux ////////////////////
;-------------------------------------------
	dataset = ['count_rate', 'intensity']
	
	IF N_Elements(b_field) GT 0 THEN BEGIN
		oPAD = MrMMS_FEEPS_Load_Data_PAD(NPA_BINS=nPA_Bins, ERANGE=erange)
	ENDIF
	
	;Create the omin-directional flux
	FOR i = 0, N_Elements(sc)      - 1 DO $
	FOR j = 0, N_Elements(mode)    - 1 DO $
	FOR k = 0, N_Elements(level)   - 1 DO $
	FOR l = 0, N_Elements(optdesc) - 1 DO $
	FOR m = 0, N_Elements(dataset) - 1 DO BEGIN

		;Look for relevant data
		vtest = StrJoin( [sc[i], instr, mode[j], level[k], optdesc[l], '*', dataset[m], 'sensorid', '*'], '_') + suffix
		MrVar_Names, varnames, vtest
		IF varnames[0] NE '' THEN BEGIN
			oOmni = MrMMS_FEEPS_Load_Data_Omni(vtest)
			oPAD  = MrMMS_FEEPS_Load_Data_PAD(vtest, oPA)
		ENDIF
		
		;Add name
		IF Obj_Valid(oOmni) THEN varnames = [varnames, oOmni.name, oPAD.name]
	ENDFOR

;-------------------------------------------
; Clean Up /////////////////////////////////
;-------------------------------------------
	;Delete data from the individual sensors
	MrVar_Names, names, '.*_(electron|ion)_(bottom|top)_(count_rate|intensity)(_500keV)?_sensorid_.*'+suffix, /REGEX
	MrVar_Delete, names
	
	;Delete spin sector
	MrVar_Names, names, '*spinsectnum'
	MrVar_Delete, names
	
	;Delete sector masks
	MrVar_Names, names, '*sector_mask*'
	MrVar_Delete, names
	
	;Delete sector masks
	MrVar_Names, names, '*channel_ids'
	MrVar_Delete, names
	
	;Delete sector masks
	MrVar_Names, names, '*pitch_angle'
	MrVar_Delete, names
END