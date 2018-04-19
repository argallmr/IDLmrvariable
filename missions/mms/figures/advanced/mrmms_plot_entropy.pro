; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_Entropy
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
;       2. density
;       3. Entropy: S = Integral{ f ln(f) } / n
;       4. Entropy per particle: S / n
;       5. Entropy: Sb = P/n^(5/3)
;       6. Scalar Pressure
;       7. Temperature Pressure
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
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
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
;       2017/02/14  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_Entropy, sc, mode, species, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext
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
	
	;Derived names
	dist_vname   = StrJoin([sc, instr, 'dist4d',        mode], '_')
	fph_vname    = StrJoin([sc, instr, 'dist', 'photo', mode], '_')
	p_scal_vname = StrJoin([sc, instr, 'presscalar',    mode], '_')
	t_scal_vname = StrJoin([sc, instr, 'tempscalar',    mode], '_')
	n_calc_vname       = StrJoin([sc, instr, 'numberdensity',  'calc', mode], '_')
	s_calc_vname       = StrJoin([sc, instr, 'entropydensity', 'calc', mode], '_')
	v_calc_vname       = StrJoin([sc, instr, 'bulkv',          'calc', mode], '_')
	p_calc_vname       = StrJoin([sc, instr, 'prestensor',     'calc', mode], '_')
	p_scal_calc_vname  = StrJoin([sc, instr, 'presscalar',     'calc', mode], '_')
	t_calc_vname       = StrJoin([sc, instr, 'temptensor',     'calc', mode], '_')
	t_scal_calc_vname  = StrJoin([sc, instr, 'tempscalar',     'calc', mode], '_')
	q_calc_vname       = StrJoin([sc, instr, 'heatflux',       'calc', mode], '_')
	sn_calc_vname      = StrJoin([sc, instr, 'entropy',        'calc', mode], '_')
	sb_calc_vname      = StrJoin([sc, instr, 'entropymaxwell', 'calc', mode], '_')
	
	ppar_calc_vname   = StrJoin([sc, instr, 'prespar',   'calc', mode], '_')
	pperp1_calc_vname = StrJoin([sc, instr, 'presperp1', 'calc', mode], '_')
	pperp2_calc_vname = StrJoin([sc, instr, 'presperp2', 'calc', mode], '_')
	
	tpar_calc_vname   = StrJoin([sc, instr, 'temppar',   'calc', mode], '_')
	tperp1_calc_vname = StrJoin([sc, instr, 'tempperp1', 'calc', mode], '_')
	tperp2_calc_vname = StrJoin([sc, instr, 'tempperp2', 'calc', mode], '_')


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
		                       /APPLY_MODEL, $
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
	oDist  = MrDist4D(f_vname, VSC=scpot_vname, SPECIES=theSpecies)
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
; Entropy Per Particle /////////////////////
;-------------------------------------------

	;Entropy per particle
	oSn = oEntropy / oDensity
	oSn -> SetName, sn_calc_vname
	oSn -> Cache

;-------------------------------------------
; Parallel & Perpendicular Temperature /////
;-------------------------------------------
	
	;Field-aligned coordinates
	oTx = MrVar_FAC( bvec_vname, v_vname, 'VXB', TIME=oSn['TIMEVAR'] )
	
	;Rotate variables
	oP     = MrVar_Get(p_vname)
	oT     = MrVar_Get(t_vname)
	oP_fac = oTx ## oP
	oT_fac = oTx ## oT
	
	Obj_Destroy, oTx
	
	;Pressures
	p_par   = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,0,0], /CACHE, NAME=ppar_calc_vname )
	p_perp1 = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,1,1], /CACHE, NAME=pperp1_calc_vname )
	p_perp2 = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,2,2], /CACHE, NAME=pperp2_calc_vname )
	p_scl   = (p_par + p_perp1 + p_perp2) / 3.0
	
	;Temperatures
	t_par   = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,0,0], /CACHE, NAME=ppar_calc_vname)
	t_perp1 = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,1,1], /CACHE, NAME=pperp1_calc_vname)
	t_perp2 = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,2,2], /CACHE, NAME=pperp2_calc_vname)
	t_scl   = (t_par + t_perp1 + t_perp2) / 3.0

;-------------------------------------------
; Maxwellian Entropy ///////////////////////
;-------------------------------------------
	oSb = (1e-9 * p_scl) / ( MrConstants('m_e') * 1e6 * oDensity )^(5.0/3.0)
	oSb -> SetName, sb_calc_vname
	oSb -> Cache

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------

	;B
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, coords], ' ' ) )
	
	;Entropy per particle
	oSn['TITLE'] = 'S/N!C(J/K)'
	
	;Entropy of Equivalent Maxwellian Distribution
	oSb['TITLE'] = 'p/$\rho$$\up5/3$'
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [b_vname, n_vname, s_calc_vname, sn_calc_vname, sb_calc_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )

	win[0] -> SetLayout, [1,1]
	win -> TrimLayout
	win.oxmargin = [15,4]

;-------------------------------------------
; Save the File ////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
			MrPrintF, 'LogText', 'Saving file to: "' + output_dir + '".'
		ENDIF
		
		;Save the figure
		fname = StrJoin([sc, instr, mode, level, 'entropy'], '_')
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	win -> Refresh
	RETURN, win
END