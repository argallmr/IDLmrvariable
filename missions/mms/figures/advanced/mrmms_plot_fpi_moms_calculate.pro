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
;   Calculate moments of the distribution function and plot them against the official
;   FPI L2 dataset. The moments calculation takes into account the FPI internal photo-
;   electron model, but the method of integration is different.
;
;       1. Bxyz, |B|
;       2. Energy-time spectrogram
;       3. PA-time spectrogram
;       4. Density
;       5. Entropy
;       6. Velocity
;       7. Scalar Pressure
;       8. Scalar Temperature
;       9. Heat Flux
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
FUNCTION MrMMS_Plot_FPI_Moms_Compare, sc, mode, species, $
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
	n_vname     = StrJoin( [sc, instr, 'numberdensity',        mode], '_')
	v_vname     = StrJoin( [sc, instr, 'bulkv',      coords, mode], '_')
	p_vname     = StrJoin( [sc, instr, 'prestensor', coords, mode], '_')
	t_vname     = StrJoin( [sc, instr, 'temptensor', coords, mode], '_')
	q_vname     = StrJoin( [sc, instr, 'heatq',      coords, mode], '_')
	
	;Source Distribution
	f_vname      = StrJoin( [sc, instr, 'dist',                  mode], '_')
	dphi_vname   = StrJoin( [sc, instr, 'startdelphi', 'count',  mode], '_' )
	parity_vname = StrJoin( [sc, instr, 'steptable',   'parity', mode], '_' )
	
	;Source Photoelectron model
	ph_dphi_vname = StrJoin( ['mms', instr, 'startdelphi', 'counts', mode], '_' ) + '_model'
	ph_f0_vname   = StrJoin( ['mms', instr, 'bgdist',      'p0',     mode], '_' ) + '_model'
	ph_f1_vname   = StrJoin( ['mms', instr, 'bgdist',      'p1',     mode], '_' ) + '_model'
	ph_e0_vname   = StrJoin( ['mms', instr, 'energy0',               mode], '_' ) + '_model'
	ph_e1_vname   = StrJoin( ['mms', instr, 'energy1',               mode], '_' ) + '_model'
	ph_scl_vname  = StrJoin( [sc,    instr, 'scl',         'model',  mode], '_' ) + '_model'
	
	;Derived names
	dist_vname   = StrJoin([sc, instr, 'dist4d',        mode], '_')
	fph_vname    = StrJoin([sc, instr, 'dist', 'photo', mode], '_')
	p_scal_vname = StrJoin([sc, instr, 'presscalar',    mode], '_')
	t_scal_vname = StrJoin([sc, instr, 'tempscalar',    mode], '_')
	n_calc_vname       = StrJoin([sc, instr, 'numberdensity', 'calc', mode], '_')
	s_calc_vname       = StrJoin([sc, instr, 'entropy',       'calc', mode], '_')
	v_calc_vname       = StrJoin([sc, instr, 'bulkv',         'calc', mode], '_')
	p_calc_vname       = StrJoin([sc, instr, 'prestensor',    'calc', mode], '_')
	p_scal_calc_vname  = StrJoin([sc, instr, 'presscalar',    'calc', mode], '_')
	t_calc_vname       = StrJoin([sc, instr, 'temptensor',    'calc', mode], '_')
	t_scal_calc_vname  = StrJoin([sc, instr, 'tempscalar',    'calc', mode], '_')
	q_calc_vname       = StrJoin([sc, instr, 'heatflux',      'calc', mode], '_')


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
;		                       /APPLY_MODEL, $
		                       COORD_SYS   = coord_sys, $
		                       LEVEL       = level, $
		                       ORIENTATION = orientation
		
		;Load FPI Moments
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = instr + '-moms', $
		                     TEAM_SITE = team_site, $
		                     VARFORMAT = [ '*energyspectr_omni*', '*pitchangdist*', $
		                                   '*numberdensity*', '*bulkv_'+coords+'*', $
		                                   '*pres*'+coords+'*', '*temp*'+coords+'*', '*heatq_'+coords+'*' ]
		
		;Spacecraft potential
		MrMMS_Load_Data, sc, 'edp', 'fast', 'l2', $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
	ENDIF

;-------------------------------------------
; Compute Moments //////////////////////////
;-------------------------------------------
	;Distribution function
	theSpecies = species EQ 'i' ? 'H' : species
	ofcorr = MrVar_Get(f_vname)
	oDist  = MrDist4D(ofcorr, VSC=scpot_vname, SPECIES=theSpecies)
	oDist -> Moments, /CACHE, $
	                  DENSITY     = oDensity, $
	                  ENTROPY     = oEntropy, $
	                  HEATFLUX    = oHeatFlux, $
	                  PRESSURE    = oPres, $
	                  TEMPERATURE = oTemp, $
	                  VELOCITY    = oVelocity
	
	;Set names
	oDensity  -> SetName, n_calc_vname
	oEntropy  -> SetName, s_calc_vname
	oHeatFlux -> SetName, q_calc_vname
	oPres     -> SetName, p_calc_vname
	oTemp     -> SetName, t_calc_vname
	oVelocity -> SetName, v_calc_vname

;-------------------------------------------
; Scalar Pressure & Temperature ////////////
;-------------------------------------------
	
	;FPI Pressure
	op_fpi = MrVar_Get(T_vname)
	op_fpi = MrScalarTS( oPres['TIMEVAR'], (oPres[*,0,0] + oPres[*,1,1] + oPres[*,2,2])/3.0, $
	                     /CACHE, $
	                     NAME = p_scal_vname )
	
	;FPI Temperature
	oT_fpi = MrVar_Get(T_vname)
	oT_fpi = MrScalarTS( oTemp['TIMEVAR'], (oTemp[*,0,0] + oTemp[*,1,1] + oTemp[*,2,2])/3.0, $
	                     /CACHE, $
	                     NAME = t_scal_vname )

	;Calculated Pressure
	oP = MrScalarTS( oPres['TIMEVAR'], (oPres[*,0,0] + oPres[*,1,1] + oPres[*,2,2])/3.0, $
	                 /CACHE, $
	                 NAME = p_scal_calc_vname )
	
	;Calculated Temperature
	oT = MrScalarTS( oTemp['TIMEVAR'], (oTemp[*,0,0] + oTemp[*,1,1] + oTemp[*,2,2])/3.0, $
	                 /CACHE, $
	                 NAME = t_scal_calc_vname )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	;B
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, coords], ' ' ) )
	
	;FPI Density
	oN = MrVar_Get(n_vname)
	oN['LABEL'] = 'moms'
	
	;FPI Velocity
	oV = MrVar_Get(v_vname)
	oV['LABEL']     = ['Vx', 'Vy', 'Vz']
	oV['COLOR']     = ['Black', 'Black', 'Black']
;	oV['LINESTYLE'] = ['-', '--', '-..']
	oV['TITLE']     = 'V!C(km/s)'
	
	;FPI Scalar Pressure
	op_fpi['TITLE'] = 'P!C(nPa)'
	op_fpi['LABEL'] = 'moms'
	
	;FPI Scalar Temperature
	ot_fpi['TITLE'] = 'T!C(eV)'
	ot_fpi['LABEL'] = 'moms'
	
	;FPI Heat Flux
	oQ_fpi = MrVar_Get(q_vname)
	oQ_fpi['LABEL']     = ['Qx', 'Qy', 'Qz']
	oQ_fpi['COLOR']     = ['Black', 'Black', 'Black']
;	oQ_fpi['LINESTYLE'] = ['-', '--', '-..']
	oQ_fpi['TITLE']     = 'Q!CmW/m^2'
	
	;DIST Density
	oDensity['COLOR'] = 'Blue'
	oDensity['LABEL'] = 'Dist'
	
	;DIST Velocity
	oVelocity -> RemoveAttr, 'LABEL'
	
	;DIST Scalar Pressure
	oP['COLOR'] = 'Blue'
	oP['LABEL'] = 'Dist'
	
	;DIST Scalar Temperature
	oT['COLOR'] = 'Blue'
	oT['LABEL'] = 'Dist'
	
	;DIST Heat Flux
	oHeatFlux -> RemoveAttr, 'LABEL'
;	oHeatFlux['LABEL'] = ['Qx', 'Qy', 'Qz']

;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS(  [b_vname, espec_vname, pad_vname, n_vname, s_calc_vname, v_vname, $
	                      p_scal_vname, t_scal_vname, q_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )

	win = MrVar_OPlotTS( [n_vname, v_vname, p_scal_vname, t_scal_vname, q_vname], $
	                     [n_calc_vname, v_calc_vname, p_scal_calc_vname, t_scal_calc_vname, q_calc_vname] )

	win[0] -> SetLayout, [1,1]
	win -> TrimLayout

	win -> Refresh
	RETURN, win
END