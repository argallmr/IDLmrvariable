; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FPI_Spectra_Compare
;
;*****************************************************************************************
;   Copyright (c) 2018, Matthew Argall                                                   ;
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
;   Compute 2D spectral plots from the FPI distribution function files and compare them
;   with the spectra provided in the moments files. The calculations take into account
;   the FPI internal photo-electron model.
;
;       1. Bxyz, |B|
;       2. Energy-time spectrogram
;       2. Computed E-t spectrogram
;       3. PA-time spectrogram
;       3. Computed PA-t spectrogram
;       4. Computed Gyrophase-t spectrogram
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;       SPECIES:    in, required, type=string
;                   Particle species. Options are {'e' | 'i'}
;
; :Keywords:
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source files.
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
;       2017/02/07  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_FPI_Spectra_Compare, sc, mode, species, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
NO_LOAD=no_load
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, Obj_New()
	ENDIF
	
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(species)   EQ 0 THEN species   = 'e'
	instr   = 'd' + species + 's'

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source B-Field
	b_vname    = StrJoin( [sc, fgm_instr, 'b',    coords, mode, level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec', coords, mode, level], '_' )
	bmag_vname = StrJoin( [sc, fgm_instr, 'bmag', coords, mode, level], '_' )
	
	;Source Moments
	scpot_vname = StrJoin( [sc, 'edp', 'scpot', 'fast', level], '_' )
	espec_vname = StrJoin( [sc, instr, 'energyspectr', 'omni', mode], '_')
	pad_vname   = StrJoin( [sc, instr, 'pitchangdist',         mode], '_')
	V_vname     = StrJoin( [sc, instr, 'bulkv',        coords, mode], '_')
	
	;Source Distribution
	f_vname = StrJoin( [sc, instr, 'dist', mode], '_')
	
	;Derived names
	espec_calc_vname = StrJoin( [sc, instr, 'energyspectr', 'omni', 'calc', mode], '_')
	pad_calc_vname   = StrJoin( [sc, instr, 'pitchangdist',         'calc', mode], '_')
	gpd_calc_vname   = StrJoin( [sc, instr, 'gyrophasedist',        'calc', mode], '_')

;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, mode, $
		                     INSTR  = fgm_instr, $
		                     LEVEL  = level, $
		                     VARFORMAT = '*b_gse*', $
		                     SUFFIX = suffix
		
		;FPI
		MrMMS_FPI_Load_Dist3D, sc, mode, species, $
		                       /APPLY_MODE
		
		;Load FPI Moments
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = instr + '-moms', $
		                     TEAM_SITE = team_site, $
		                     VARFORMAT = [ '*energyspectr_omni*', '*pitchangdist*', $
		                                   '*bulkv_'+coords+'*' ]
		
		;Spacecraft potential
		MrMMS_Load_Data, sc, 'edp', 'fast', 'l2', $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
	ENDIF

;-------------------------------------------
; Compute Moments //////////////////////////
;-------------------------------------------
	;Energy-time distribution
	theSpecies = species EQ 'i' ? 'H' : species
	oDist      = MrDist4D(f_vname, VSC=scpot_vname, SPECIES=theSpecies)
	oDist     -> ConvertUnits, 'EFLUX'
;	oESpec     = oDist -> ESpec(/CACHE, NAME=espec_calc_vname)
	
	oT         = MrVar_FAC(bvec_vname, v_vname, 'VXB', TIME=v_vname)
	oDist_fac  = oDist -> Rotate(oT)
	
	;Spectra
	oDist_fac  -> Spectra, ESPEC=oESpec, PHISPEC=oPhiSpec, THETASPEC=oThetaSpec, /CACHE
	oESpec     -> SetName, espec_calc_vname
	oPhiSpec   -> SetName, gpd_calc_vname
	oThetaSpec -> SetName, pad_calc_vname
	;Free memory
	Obj_Destroy, [oDist, oT, oDist_fac]

;-------------------------------------------
; Compute Moments //////////////////////////
;-------------------------------------------
	;Energy-time distribution
;	theSpecies = species EQ 'i' ? 'H' : species
;	oDist  = MrDist4D(f_vname, VSC=scpot_vname, SPECIES=theSpecies)
;	oEspec = oDist -> ESpec(/CACHE, NAME=espec_calc_vname)
	
	;PAD
;	oT        = MrVar_FAC(bvec_vname, v_vname, 'VXB', TIME=v_vname)
;	oDist_fac = oDist -> Rotate(oT)
;	oPAD      = oDist_fac -> ThetaSpec(/CACHE, NAME=pad_calc_vname)
	
	;GPD
;	oGPD = oDist_fac -> PhiSpec(/CACHE, NAME=gpd_calc_vname)
	
	;Free memory
;	Obj_Destroy, [oDist, oDist_fac]

;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [b_vname, espec_vname, espec_calc_vname, pad_vname, pad_calc_vname, gpd_calc_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )

	win[0] -> SetLayout, [1,1]
	win -> TrimLayout

	win -> Refresh
	RETURN, win
END