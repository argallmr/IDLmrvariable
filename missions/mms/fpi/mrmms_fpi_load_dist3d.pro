; docformat = 'rst'
;
; NAME:
;       MrMMS_FPI_Load_Dist3D
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
; PURPOSE:
;+
;   Load FPI distribution function data and rotate the instrument look directions
;   into a field-aligned coordinate system.
;
; :Categories:
;       CDF Utilities
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
;       2016/08/20  -   Written by Matthew Argall
;       2017/01/04  -   Remove internally generated photoelectrons. - MRA
;       2017/10/16  -   All variables associated with the photoelectron model are removed
;                           from the cache. - MRA
;       2017/10/23  -   Added the RAGER keyword to load 7.5ms moments data. Pass LEVEL
;                           on to helper functions because of change to MrMMS_Get_*. - MRA
;-
;*****************************************************************************************
;+
;   Create new energy tables and set the distribution function's DEPEND_1-3 variable
;   attributes to be variables with physical units.
;
; :Params:
;       NAMES:          in, required, type=string
;                       Name of the distribution function variable.
;-
PRO MrMMS_FPI_Load_Dist3D_Meta, name
	Compile_Opt idl2
	On_Error, 2
	
	;Variable name parts
	;   - sc_instr_energy_mode[_suffix]
	;      0   1     2     3      4+
	;    | prefix |      |   suffix   |
	parts  = StrSplit(name, '_', /EXTRACT)
	prefix = StrJoin(parts[0:1], '_') + '_'
	suffix = '_' + StrJoin(parts[3:*], '_')
	
	;Names of phi, theta, and new energy tables
	dist_vname   = prefix + 'dist'             + suffix
	phi_vname    = prefix + 'phi'              + suffix
	theta_vname  = prefix + 'theta'            + suffix
	e0_vname     = prefix + 'energy0'          + suffix
	e1_vname     = prefix + 'energy1'          + suffix
	parity_vname = prefix + 'steptable_parity' + suffix

;-----------------------------------------------------
; Energy Tables \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;BRST mode
	;   - Two energy tables alternate based on STEPTABLE_PARITY
	;   - Energies are combined into a single array as of v3.0.0
	IF MrVar_IsCached(e0_vname) THEN BEGIN
		;Get the sector, pixel, and energy tables
		oParity  = MrVar_Get( parity_vname )
		oEnergy0 = MrVar_Get( e0_vname )
		oEnergy1 = MrVar_Get( e1_vname )
		oDist    = MrVar_get( dist_vname )
	
		;Names of new energy tables
		eTable_vname = prefix + 'energy_table'   + suffix
		eMean_vname  = prefix + 'energy_geomean' + suffix
	
		;Create new energy table
		;   - One time-dependent energy table
		;   - One combined (average) energy table
		energy      = transpose( [ [oEnergy0['DATA']], [oEnergy1['DATA']] ] )
		energy_full = MrVariable( energy[oParity['DATA'], *],                NAME=eTable_vname, /CACHE)
		energy_mean = MrVariable( reform( sqrt(energy[0,*] * energy[1,*]) ), NAME=eMean_vname,  /CACHE)
	
	;SRVY mode
	ENDIF ELSE BEGIN
		;There is only one energy table
		eTable_vname = prefix + 'energy' + suffix
	ENDELSE

;-----------------------------------------------------
; Attributes \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Set the distribution function dependencies
	oDist = MrVar_Get(dist_vname)
	oDist['DEPEND_1'] = MrVar_Get(phi_vname)
	oDist['DEPEND_2'] = MrVar_Get(theta_vname)
	oDist['DEPEND_3'] = MrVar_Get(eTable_vname)

;-----------------------------------------------------
; Adjust DELTA_(PLUS|MINUS)_VAR \\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	; The DELTA_(PLUS|MINUS)_VAR for Phi are twice what they
	; should be
	;
	oPhi   = oDist['DEPEND_1']
	oTheta = oDist['DEPEND_2']
	
	;PHI: Check if the correction is needed
	;   - DELTA_MINUS_VAR & DELTA_PLUS_VAR are the same variable
	odMinus = oPhi['DELTA_MINUS_VAR']
	IF Abs(odMinus['DATA', 0] - 11.25/2.0) GT 1.0 THEN BEGIN
		MrPrintF, 'LogWarn', 'Correcting FPI Dist PHI DELTA+/-.'
		odMinus -> SetData, odMinus['DATA'] / 2.0
	ENDIF ELSE BEGIN
		MrPrintF, 'LogWarn', 'FPI Dist PHI DELTA+/- Fixed.'
	ENDELSE

;-----------------------------------------------------
; Clean Up Variables \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	; The *energyspectr* variables in the FPI moments files
	; also have and energy table as a support data variable.
	; So that the one created here is not overwritten by
	; moments files, it is removed from the cache and left
	; in-tact as a variable attribute.
	;
	; Do not delete other energy tables or parity state. They
	; can be useful elsewhere (e.g. calculating the photo-
	; electron model distribution) and are not stored anywhere
	; else (i.e. as a variable attribute).
	;
	MrVar_Remove, [phi_vname, theta_vname, eTable_vname]
END


;+
;   Load the internal photoelectron model and correct the distribution.
;
; :Params:
;       NAMES:          in, required, type=string
;                       Name of the distribution function variable.
;-
PRO MrMMS_FPI_Load_Dist3D_ePhoto, name
	Compile_Opt idl2
	On_Error, 2
	
	;Variable name parts
	;   - sc_instr_energy_mode[_suffix]
	;      0   1     2     3      4+
	;    | prefix |      |   suffix   |
	parts  = strsplit(name, '_', /EXTRACT)
	sc     = parts[0]
	instr  = parts[1]
	mode   = parts[3]
	suffix = N_Elements(parts) EQ 4 ? '' : StrJoin(parts[4:*], '_')
	
	;Source Distribution
	f_vname      = StrJoin( [sc, instr, 'dist',                  mode], '_' ) + suffix
	dphi_vname   = StrJoin( [sc, instr, 'startdelphi', 'count',  mode], '_' ) + suffix
	parity_vname = StrJoin( [sc, instr, 'steptable',   'parity', mode], '_' ) + suffix

	;Source Photoelectron model
	model_suffix   = '_model'
	ph_scl_vname   = StrJoin( [sc,    instr, 'scl',         'model',  mode], '_' ) + model_suffix
	ph_dphi_vname  = StrJoin( ['mms', instr, 'startdelphi', 'counts', mode], '_' ) + model_suffix
	ph_phi_vname   = StrJoin( ['mms', instr, 'phi',                   mode], '_' ) + model_suffix
	ph_theta_vname = StrJoin( ['mms', instr, 'theta',                 mode], '_' ) + model_suffix
	IF mode EQ 'brst' THEN BEGIN
		ph_f0_vname = StrJoin( ['mms', instr, 'bgdist',  'p0', mode], '_' ) + model_suffix
		ph_f1_vname = StrJoin( ['mms', instr, 'bgdist',  'p1', mode], '_' ) + model_suffix
		ph_e0_vname = StrJoin( ['mms', instr, 'energy0',       mode], '_' ) + model_suffix
		ph_e1_vname = StrJoin( ['mms', instr, 'energy1',       mode], '_' ) + model_suffix
	ENDIF ELSE BEGIN
		ph_f0_vname = StrJoin( ['mms', instr, 'bgdist', mode], '_' ) + model_suffix
		ph_e0_vname = StrJoin( ['mms', instr, 'energy', mode], '_' ) + model_suffix
	ENDELSE

	;Derived names
	ph_f_vname = StrJoin([sc, instr, 'dist', 'photo', mode], '_')

	;Photoelectron model
	;   - Only exist for electrons
	MrMMS_FPI_Load_Models, sc, mode, 'e', $
	                       LEVEL  = level, $
	                       SUFFIX = model_suffix

	;Create f_photo
	IF mode EQ 'brst' THEN BEGIN
		ofPhoto = MrMMS_FPI_Dist_Photo( f_vname, dphi_vname, ph_scl_vname, ph_dphi_vname, $
		                                ph_f0_vname, ph_e0_vname, ph_f1_vname, ph_e1_vname, $
		                                parity_vname, $
		                                /CACHE, $
		                                NAME = ph_f_vname )
	ENDIF ELSE BEGIN
		ofPhoto = MrMMS_FPI_Dist_Photo( f_vname, dphi_vname, ph_scl_vname, ph_dphi_vname, $
		                                ph_f0_vname, ph_e0_vname, $
		                                /CACHE, $
		                                NAME = ph_f_vname )
	ENDELSE
	
	;Subtract the model from the distribution
	;   - Set all negative contributions to zero
	oDist   = MrVar_Get(f_vname)
	oDistTS = oDist - ofPhoto > 0
	
	;Replace the old variable with the new
	oDist -> CopyAttrTo, oDistTS
	MrVar_Replace, oDist, oDistTS
	oDistTS -> SetName, f_vname
	
	;Remove variables
	MrVar_Delete, [ph_f_vname, ph_dphi_vname, ph_f0_vname, ph_e0_vname, ph_scl_vname, ph_phi_vname, ph_theta_vname]
	IF mode EQ 'brst' THEN MrVar_Delete, [ph_f1_vname, ph_e1_vname]
END


;+
;   Rotate the distribution into a field-aligned coordinate system.
;
; :Params:
;       NAMES:          in, required, type=string
;                       Name of the distribution function variable.
;-
PRO MrMMS_FPI_Load_Dist3D_Rotate, name, par, perp, fac, $
ORIENTATION=orientation, $
VARNAMES=varnames
	Compile_Opt idl2
	On_Error, 2
	
	;Variable name parts
	;   - sc_instr_energy_mode[_suffix]
	;      0   1     2     3      4+
	;    | prefix |      |   suffix   |
	parts   = strsplit(name, '_', /EXTRACT)
	sc      = parts[0]
	instr   = parts[1]
	mode    = parts[3]
	species = StrMid(instr, 1, 1)
	suffix  = N_Elements(parts) EQ 4 ? '' : StrJoin(parts[4:*], '_')
	
	;Source Distribution
	f_vname = StrJoin( [sc, instr, 'dist', mode], '_' ) + suffix
	
	;Derived names
	dist_fac_vname  = StrJoin( [sc, instr, 'dist',  'fac', mode], '_' ) + suffix
	phi_fac_vname   = StrJoin( [sc, instr, 'phi',   'fac', mode], '_' ) + suffix
	theta_fac_vname = StrJoin( [sc, instr, 'theta', 'fac', mode], '_' ) + suffix

;-----------------------------------------------------
; Rotate Coords to FAC \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Grab the distribution
	oDist = MrVar_Get(name)

	;Convert from instrument coordinates to field-aligned coordinates
	MrVar_Grid_sphere2fac, par, oDist['DEPEND_1'], oDist['DEPEND_2'], oPhi_FAC, oTheta_FAC, $
	                       ORIENTATION = orientation, $
	                       /SPHERE, $
	                       FAC         = fac, $
	                       VEC         = perp

;-----------------------------------------------------
; Store Variables \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Set names
	oPhi_FAC   -> SetName, phi_fac_vname
	oTheta_FAC -> SetName, theta_fac_vname
	
	;Add to cache
	oPhi_FAC   -> Cache
	oTheta_FAC -> Cache
	
	;Set time dependence
	IF mode EQ 'brst' THEN BEGIN
MrPrintF, 'LogText', 'Stopping to see if DEPEND_0 was set correctly.'
stop
;		oPhi_FAC['DEPEND_0']   = oDist['DEPEND_0']
;		oTheta_FAC['DEPEND_0'] = oDist['DEPEND_0']
	ENDIF

	;Distribution
	oDist_FAC  = oDist -> Copy()
	oDist_FAC -> SetName, dist_fac_vname
	oDist_FAC -> Cache
	oDist_FAC['DEPEND_1'] = oPhi_FAC
	oDist_FAC['DEPEND_2'] = oTheta_FAC

	;Variable names
	varnames = [dist_fac_vname, phi_fac_vname, theta_fac_vname, oDist['DEPEND_3']]
END


;+
;   Rotate the distribution into a field-aligned coordinate system.
;
; :Params:
;       NAMES:          in, required, type=string
;                       Name of the distribution function variable.
;-
PRO MrMMS_FPI_Load_Dist3D_Orient, name, par, perp, fac, $
ORIENTATION=orientation, $
VARNAMES=varnames
	Compile_Opt idl2
	On_Error, 2
	
	;Variable name parts
	;   - sc_instr_energy_mode[_suffix]
	;      0   1     2     3      4+
	;    | prefix |      |   suffix   |
	parts   = strsplit(name, '_', /EXTRACT)
	sc      = parts[0]
	instr   = parts[1]
	mode    = parts[3]
	species = StrMid(instr, 1, 1)
	suffix  = N_Elements(parts) EQ 4 ? '' : StrJoin(parts[4:*], '_')
	
	;Source Distribution
	f_vname     = StrJoin( [sc, instr, 'dist',  mode], '_' ) + suffix
	phi_vname   = StrJoin( [sc, instr, 'phi',   mode], '_' ) + suffix
	theta_vname = StrJoin( [sc, instr, 'theta', mode], '_' ) + suffix
	
	;Derived names
	dist_fac_vname  = StrJoin( [sc, instr, 'dist',  'grid', mode], '_' ) + suffix
	phi_fac_vname   = StrJoin( [sc, instr, 'phi',   'grid', mode], '_' ) + suffix
	theta_fac_vname = StrJoin( [sc, instr, 'theta', 'grid', mode], '_' ) + suffix

;-----------------------------------------------------
; Change Orientations \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Create the cartesian grid
	oDist = MrVar_Get(f_name)
	MrVar_Grid_MakeCart, oDist['DEPEND_1'], oDist['DEPEND_2'], oX, oY, oZ, /DEGREES
	
	;Convert back to spherical grid
	MrVar_Grid_Cart2Sphere, oX, oY, oZ, oPhiGrid, oThetaGrid, $
	                        /DEGREES, $
	                        ORIENTATION = orientation

MrPrintF, 'LogText', 'Stopping to see if DEPEND_0 was set correctly.'
stop
	
	;Save PHI
	oPhiGrid -> SetName, phi_grid_vname
	oPhiGrid -> Cache
	oPhiGrid['DEPEND_0'] = oDistTS['DEPEND_0']
	oPhiGrid['DEPEND_1'] = phi_vname
	oPhiGrid['DEPEND_2'] = theta_vname
	
	;Save THETA
	oThetaGrid -> SetName, theta_grid_vname
	oThetaGrid -> Cache
	oThetaGrid['DEPEND_0'] = oDist['DEPEND_0']
	oThetaGrid['DEPEND_1'] = phi_vname
	oThetaGrid['DEPEND_2'] = theta_vname
	
	;Update varnames
	varnames = [f_vname, phi_grid_vname, theta_grid_vname]
END


;+
;   Find and read MMS FPI data.
;
; :Params:
;       SC:                 in, required, type=string/strarr
;                           The MMS spacecraft identifier. Options are:
;                               {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:               in, required, type=string/strarr
;                           Telemetry mode of the data. Options include: {'fast' | 'brst'}
;       SPECIES:            in, optional, type=string, default='e'
;                           Particle species for which the distribution function is
;                               to be loaded. Options are: {'e' | 'i'} for electrons
;                               and ions, respectively.
;
; :Keywords:
;       COORD_SYS:          in, optional, type=string, default='gse'
;                           Coordinate system of the original distribution function.
;                               Options are: {'dbcs' | 'gse'}.
;       FAC:                in, optional, type=string
;                           The field-aligned coordinate system into which the instrument
;                               look directions should be rotated. Used only if `PAR` and
;                               `PERP` are given. See MrVar_FAC for options.
;       LEVEL:              in, optional, type=string/strarr, default='l2'
;                           Data quality level. Options include:
;                               {'sitl' | 'l1b' | 'l2'}
;       ORIENTATION:        in, optional, type=integer, default=3
;                           Orientation of the sphereical grid. See
;                               MrVar_Grid_Cart2Sphere.pro for options.
;       PAR:                in, optional, type=objref/string/integer
;                           The name, number or MrVectorTS object of the vector that
;                               defines the parallel- (z-) axis of the field-aligned
;                               coordinate system into which the distribution is rotated.
;       PERP:               in, optional, type=objref/string/integer
;                           The name, number or MrVectorTS object of the vector that
;                               defines the perpendicular axis of the field-aligned
;                               coordinate system into which the distribution is rotated.
;                               Used only if `PAR` is given.
;       PHOTO_MODEL:        in, optional, type=boolean, default=0
;                           If set, the photoelectron model that corrects this distribution
;                               will be loaded into the cache.
;       TEAM_SITE:          in, optional, type=boolean, default=0
;                           If set, requests will be sent to the team site (password
;                               required, L1A and above). Automatically set if `LEVEL`
;                               is below level 2. This option is sticky.
;       TRANGE:             out, optional, type=string, default=MrVar_GetTRange()
;                           Start and end times over which to read data.
;       VARNAMES:           out, optional, type=strarr
;                           Names of all variables that have been loaded.
;-
pro MrMMS_FPI_Load_Dist3D, sc, mode, species, $
APPLY_MODEL=apply_model, $
COORD_SYS=coord_sys, $
FAC=fac, $
LEVEL=level, $
ORIENTATION=orientation, $
PAR=par, $
PERP=perp, $
RAGER=rager, $
SUFFIX=suffix, $
TEAM_SITE=team_site, $
TRANGE=trange, $
VARNAMES=varnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return
	endif
	
	;Check inputs
	cs = n_elements(coord_sys) eq 0 ? 'gse' : strlowcase(coord_sys)
	tf_apply_model = keyword_set(apply_model)
	tf_rager       = keyword_set(rager)
	if n_elements(level)   eq 0 then level   = 'l2'
	if n_elements(species) eq 0 then species = 'e'
	if n_elements(suffix)  eq 0 then suffix  = ''
	
	;Conflicts
	if n_elements(sc)    ne 1 then message, 'SC must be scalar.'
	if n_elements(mode)  ne 1 then message, 'MODE must be scalar.'
	
	;Fast and Brst are organized differently
	if mode eq 'fast' then begin
		varformat = [ StrJoin( [sc, 'd'+species+'s', 'dist',                  mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'phi',                   mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'theta',                 mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'energy',                mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'startdelphi', 'count',  mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'startdelphi', 'angle',  mode], '_' ) ]
	ENDIF ELSE BEGIN
		varformat = [ StrJoin( [sc, 'd'+species+'s', 'dist',                  mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'phi',                   mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'theta',                 mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'energy',                mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'energy0',               mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'energy1',               mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'startdelphi', 'count',  mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'startdelphi', 'angle',  mode], '_' ), $
		              StrJoin( [sc, 'd'+species+'s', 'steptable',   'parity', mode], '_' ) ]
	ENDELSE

;-----------------------------------------------------
; Load the FPI \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	dist_level = level
	IF tf_rager THEN BEGIN
		!MrMMS -> SetProperty, OFFLINE=offline, DROPBOX_ROOT=dropbox_root
		!MrMMS -> SetProperty, OFFLINE=1B,      DROPBOX_ROOT='/home/argall/data/rager/'
		dist_level = 'acr'
	ENDIF
	
	;Get the FPI distribution function
	MrMMS_Load_Data, sc, 'fpi', mode, dist_level, $
	                 OPTDESC   = 'd' + species + 's-dist', $
	                 SUFFIX    = suffix, $
	                 TEAM_SITE = team_site, $
	                 VARFORMAT = varformat
	
	IF tf_rager THEN !MrMMS -> SetProperty, OFFLINE=offline, DROPBOX_ROOT=dropbox_root

;-----------------------------------------------------
; Add Metadata \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Associate variable attributes with DEPEND_[0-3]
	dist_vname = sc + '_d' + species + 's_dist_' + mode
	MrMMS_FPI_Load_Dist3D_Meta, dist_vname
	varnames = dist_vname

;-----------------------------------------------------
; Correct for Internal Photoelectrons \\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Correct the distribution for internally generated photoelectrons
	IF tf_apply_model && species EQ 'e' $
		THEN MrMMS_FPI_Load_Dist3D_ePhoto, dist_vname
	
;-----------------------------------------------------
; FAC \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Rotate to FAC
	IF N_Elements(par) GT 0 THEN BEGIN
		MrMMS_FPI_Load_Dist3D_Rotate, dist_vname, par, perp, fac, $
		                              ORIENTATION = orientation, $
		                              VARNAMES    = varnames

;-----------------------------------------------------
; Create Grid \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	ENDIF ELSE IF N_Elements(orientation) GT 0 THEN BEGIN
		MrMMS_FPI_Load_Dist3D_Orient, dist_vname, par, perp, fac, $
		                              ORIENTATION = orientation, $
		                              VARNAMES    = varnames
	ENDIF
	
	;Remove extra variables
	MrVar_Remove, varformat[1:*] + suffix
END