; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FPI_CalcMoms
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
;       COORDS:     in, optional, type=string, default='gse'
;                   Coordinate system in which to load the data. Options are: {'dbcs' | 'gse' | 'gsm'}
;       FGM_INSTR:  in, optional, type=string, default='fgm'
;                   The FGM instrument to use. Options are: {'afg' | 'dfg' | 'fgm'}
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'ql' | 'l2'}
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source files.
;       OUTPUT_DIR: in, optional, type=string, default='~/figures/'
;                   Directory in which to save the figure. If neither `OUTPUT_DIR` or
;                       `OUTPUT_EXT` are given, no file is made.
;       OUTPUT_EXT: in, optional, type=string/strarr, default='png'
;                   Extension (and file type) of the figure. Options include
;                       'png', 'jpeg', 'tiff', 'ps', 'eps', 'pdf'. If neither
;                       `OUTPUT_DIR` or `OUTPUT_EXT` are given, no file is made.
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
FUNCTION MrMMS_Plot_FPI_Entropy, sc, mode, $
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
	des_instr = 'des'
	dis_instr = 'dis'

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;
	; SOURCE
	;
	
	
	;FGM
	b_vname    = StrJoin( [sc, fgm_instr, 'b',             coords, mode, level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec',          coords, mode, level], '_' )
	bmag_vname = StrJoin( [sc, fgm_instr, 'bmag',          coords, mode, level], '_' )
	
	;SCPOT
	scpot_vname = StrJoin( [sc, 'edp', 'scpot', 'fast', level], '_' )
	
	;DES-MOMS
;	espec_vname  = StrJoin( [sc, des_instr, 'energyspectr', 'omni', mode], '_')
;	pad_vname    = StrJoin( [sc, des_instr, 'pitchangdist',         mode], '_')
	ne_vname     = StrJoin( [sc, des_instr, 'numberdensity',        mode], '_')
	ve_vname     = StrJoin( [sc, des_instr, 'bulkv',      coords, mode], '_')
	pe_vname     = StrJoin( [sc, des_instr, 'prestensor', coords, mode], '_')
	te_vname     = StrJoin( [sc, des_instr, 'temptensor', coords, mode], '_')
	qe_vname     = StrJoin( [sc, des_instr, 'heatq',      coords, mode], '_')
	
	;DIS-MOMS
;	espec_vname  = StrJoin( [sc, dis_instr, 'energyspectr', 'omni', mode], '_')
;	pad_vname    = StrJoin( [sc, dis_instr, 'pitchangdist',         mode], '_')
	ni_vname     = StrJoin( [sc, dis_instr, 'numberdensity',        mode], '_')
	vi_vname     = StrJoin( [sc, dis_instr, 'bulkv',      coords, mode], '_')
	pi_vname     = StrJoin( [sc, dis_instr, 'prestensor', coords, mode], '_')
	ti_vname     = StrJoin( [sc, dis_instr, 'temptensor', coords, mode], '_')
	qi_vname     = StrJoin( [sc, dis_instr, 'heatq',      coords, mode], '_')

	;DES-DIST
	f_des_vname      = StrJoin( [sc, des_instr, 'dist',                  mode], '_')
	dphi_des_vname   = StrJoin( [sc, des_instr, 'startdelphi', 'count',  mode], '_' )
	parity_des_vname = StrJoin( [sc, des_instr, 'steptable',   'parity', mode], '_' )
	
	;DES-MODEL
	ph_dphi_vname = StrJoin( ['mms', des_instr, 'startdelphi', 'counts', mode], '_' ) + '_model'
	ph_f0_vname   = StrJoin( ['mms', des_instr, 'bgdist',      'p0',     mode], '_' ) + '_model'
	ph_f1_vname   = StrJoin( ['mms', des_instr, 'bgdist',      'p1',     mode], '_' ) + '_model'
	ph_e0_vname   = StrJoin( ['mms', des_instr, 'energy0',               mode], '_' ) + '_model'
	ph_e1_vname   = StrJoin( ['mms', des_instr, 'energy1',               mode], '_' ) + '_model'
	ph_scl_vname  = StrJoin( [sc,    des_instr, 'scl',         'model',  mode], '_' ) + '_model'
	
	;DIS-DIST
	f_dis_vname   = StrJoin( [sc, dis_instr, 'dist',                  mode], '_')
	
	;
	; DERIVED
	;
	
	;Pressures
	b_press_vname  = StrJoin( [sc, fgm_instr, 'b', 'press', mode, level], '_' )
	ptot_fpi_vname = StrJoin( [sc, 'fpi',     'ptot',       mode, level], '_' )
	ptot_vname     = StrJoin( [sc,            'ptot',       mode, level], '_' )
	
	;DES-DIST
	dist_vname   = StrJoin([sc, des_instr, 'dist4d',        mode], '_')
	fph_vname    = StrJoin([sc, des_instr, 'dist', 'photo', mode], '_')
	
	;DES-MOMS
	pe_par_vname   = StrJoin([sc, des_instr, 'prespar',    mode], '_')
	pe_perp1_vname = StrJoin([sc, des_instr, 'presperp1',  mode], '_')
	pe_perp2_vname = StrJoin([sc, des_instr, 'presperp2',  mode], '_')
	pe_perp_vname  = StrJoin([sc, des_instr, 'presperp',   mode], '_')
	pe_scal_vname  = StrJoin([sc, des_instr, 'presscalar', mode], '_')
	te_par_vname   = StrJoin([sc, des_instr, 'temppar',    mode], '_')
	te_perp1_vname = StrJoin([sc, des_instr, 'tempperp1',  mode], '_')
	te_perp2_vname = StrJoin([sc, des_instr, 'tempperp2',  mode], '_')
	te_perp_vname  = StrJoin([sc, des_instr, 'tempperp',   mode], '_')
	te_scal_vname  = StrJoin([sc, des_instr, 'tempscalar', mode], '_')
	sne_vname      = StrJoin([sc, des_instr, 'entropy',    mode], '_')
	sbe_vname      = StrJoin([sc, des_instr, 'entropymaxwell', mode], '_')
	
	;DES-MOMS
	pi_par_vname   = StrJoin([sc, dis_instr, 'prespar',    mode], '_')
	pi_perp1_vname = StrJoin([sc, dis_instr, 'presperp1',  mode], '_')
	pi_perp2_vname = StrJoin([sc, dis_instr, 'presperp2',  mode], '_')
	pi_perp_vname  = StrJoin([sc, dis_instr, 'presperp',   mode], '_')
	pi_scal_vname  = StrJoin([sc, dis_instr, 'presscalar', mode], '_')
	ti_par_vname   = StrJoin([sc, dis_instr, 'temppar',    mode], '_')
	ti_perp1_vname = StrJoin([sc, dis_instr, 'tempperp1',  mode], '_')
	ti_perp2_vname = StrJoin([sc, dis_instr, 'tempperp2',  mode], '_')
	ti_perp_vname  = StrJoin([sc, dis_instr, 'tempperp',   mode], '_')
	ti_scal_vname  = StrJoin([sc, dis_instr, 'tempscalar', mode], '_')
	sni_vname      = StrJoin([sc, dis_instr, 'entropy',        mode], '_')
	sbi_vname      = StrJoin([sc, dis_instr, 'entropymaxwell', mode], '_')
	
	;DES Integrated Moments
	n_des_vname      = StrJoin([sc, des_instr, 'numberdensity',  'calc', mode], '_')
	s_des_vname      = StrJoin([sc, des_instr, 'entropydensity', 'calc', mode], '_')
	v_des_vname      = StrJoin([sc, des_instr, 'bulkv',          'calc', mode], '_')
	p_des_vname      = StrJoin([sc, des_instr, 'prestensor',     'calc', mode], '_')
	t_des_vname      = StrJoin([sc, des_instr, 'temptensor',     'calc', mode], '_')
	q_des_vname      = StrJoin([sc, des_instr, 'heatflux',       'calc', mode], '_')
	sn_des_vname     = StrJoin([sc, des_instr, 'entropy',        'calc', mode], '_')
	sb_des_vname     = StrJoin([sc, des_instr, 'entropymaxwell', 'calc', mode], '_')
	ppar_des_vname   = StrJoin([sc, des_instr, 'prespar',    'calc', mode], '_')
	pperp1_des_vname = StrJoin([sc, des_instr, 'presperp1',  'calc', mode], '_')
	pperp2_des_vname = StrJoin([sc, des_instr, 'presperp2',  'calc', mode], '_')
	p_scal_des_vname = StrJoin([sc, des_instr, 'presscalar', 'calc', mode], '_')
	tpar_des_vname   = StrJoin([sc, des_instr, 'temppar',    'calc', mode], '_')
	tperp1_des_vname = StrJoin([sc, des_instr, 'tempperp1',  'calc', mode], '_')
	tperp2_des_vname = StrJoin([sc, des_instr, 'tempperp2',  'calc', mode], '_')
	t_scal_des_vname = StrJoin([sc, des_instr, 'tempscalar', 'calc', mode], '_')
	
	;DIS Integrated Moments
	n_dis_vname      = StrJoin([sc, dis_instr, 'numberdensity',  'calc', mode], '_')
	s_dis_vname      = StrJoin([sc, dis_instr, 'entropydensity', 'calc', mode], '_')
	v_dis_vname      = StrJoin([sc, dis_instr, 'bulkv',          'calc', mode], '_')
	p_dis_vname      = StrJoin([sc, dis_instr, 'prestensor',     'calc', mode], '_')
	t_dis_vname      = StrJoin([sc, dis_instr, 'temptensor',     'calc', mode], '_')
	q_dis_vname      = StrJoin([sc, dis_instr, 'heatflux',       'calc', mode], '_')
	sn_dis_vname     = StrJoin([sc, dis_instr, 'entropy',        'calc', mode], '_')
	sb_dis_vname     = StrJoin([sc, dis_instr, 'entropymaxwell', 'calc', mode], '_')
	ppar_dis_vname   = StrJoin([sc, dis_instr, 'prespar',    'calc', mode], '_')
	pperp1_dis_vname = StrJoin([sc, dis_instr, 'presperp1',  'calc', mode], '_')
	pperp2_dis_vname = StrJoin([sc, dis_instr, 'presperp2',  'calc', mode], '_')
	p_scal_dis_vname = StrJoin([sc, dis_instr, 'presscalar', 'calc', mode], '_')
	tpar_dis_vname   = StrJoin([sc, dis_instr, 'temppar',    'calc', mode], '_')
	tperp1_dis_vname = StrJoin([sc, dis_instr, 'tempperp1',  'calc', mode], '_')
	tperp2_dis_vname = StrJoin([sc, dis_instr, 'tempperp2',  'calc', mode], '_')
	t_scal_dis_vname = StrJoin([sc, dis_instr, 'tempscalar', 'calc', mode], '_')

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
		
		;DIS-MOMS
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = 'dis-moms', $
		                     TEAM_SITE = team_site, $
		                     VARFORMAT = [ '*energyspectr_omni*', '*pitchangdist*', $
		                                   '*numberdensity*', '*bulkv_'+coords+'*', $
		                                   '*pres*'+coords+'*', '*temp*'+coords+'*', '*heatq_'+coords+'*' ]

		;DES-MOMS
		;   - Must come before MrMMS_FPI_Load_Dist3D
		;   - Will cause 'mms1_dis_energy_delta_brst' to be destroyed
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = 'des-moms', $
		                     TEAM_SITE = team_site, $
		                     VARFORMAT = [ '*energyspectr_omni*', '*pitchangdist*', $
		                                   '*numberdensity*', '*bulkv_'+coords+'*', $
		                                   '*pres*'+coords+'*', '*temp*'+coords+'*', '*heatq_'+coords+'*' ]
		
		;DIS
		MrMMS_FPI_Load_Dist3D, sc, mode, 'i', $
		                       COORD_SYS   = coord_sys, $
		                       LEVEL       = level, $
		                       ORIENTATION = orientation
		
		;DES
		MrMMS_FPI_Load_Dist3D, sc, mode, 'e', $
		                       /APPLY_MODEL, $
		                       COORD_SYS   = coord_sys, $
		                       LEVEL       = level, $
		                       ORIENTATION = orientation
		
		;SCPOT
		MrMMS_Load_Data, sc, 'edp', 'fast', 'l2', $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
	ENDIF

;-------------------------------------------
; DES: Compute Moments /////////////////////
;-------------------------------------------
	;Distribution function
	species = 'e'
	oDist   = MrDist4D_v1(f_des_vname, VSC=scpot_vname, SPECIES=species)
	oDist  -> Moments, /CACHE, $
	                   DENSITY     = oN_des, $
	                   ENTROPY     = oS_des, $
	                   HEATFLUX    = oQ_des, $
	                   PRESSURE    = oP_des, $
	                   TEMPERATURE = oT_des, $
	                   VELOCITY    = oV_des

	;Set names
	oN_des -> SetName, n_des_vname
	oS_des -> SetName, s_des_vname
	oQ_des -> SetName, q_des_vname
	oP_des -> SetName, p_des_vname
	oT_des -> SetName, t_des_vname
	oV_des -> SetName, v_des_vname

;-------------------------------------------
; DIS: Compute Moments /////////////////////
;-------------------------------------------
	;Distribution function
	species = 'H'
	ofDist  = MrVar_Get(f_dis_vname)
	oDist   = MrDist4D_v1(ofDist, VSC=scpot_vname, SPECIES=species)
	oDist  -> Moments, /CACHE, $
	                   DENSITY     = oN_dis, $
	                   ENTROPY     = oS_dis, $
	                   HEATFLUX    = oQ_dis, $
	                   PRESSURE    = oP_dis, $
	                   TEMPERATURE = oT_dis, $
	                   VELOCITY    = oV_dis
	
	;Set names
	oN_dis -> SetName, n_dis_vname
	oS_dis -> SetName, s_dis_vname
	oQ_dis -> SetName, q_dis_vname
	oP_dis -> SetName, p_dis_vname
	oT_dis -> SetName, t_dis_vname
	oV_dis -> SetName, v_dis_vname

;-------------------------------------------
; Entropy Per Particle /////////////////////
;-------------------------------------------

	;DES
	oSn_des = oS_des / oN_des
	oSn_des -> SetName, sn_des_vname
	oSn_des -> Cache

	;DIS
	oSn_dis = oS_dis / oN_dis
	oSn_dis = oSn_dis + (oSn_des.min - oSn_dis.min)
	oSn_dis -> SetName, sn_dis_vname
	oSn_dis -> Cache

;-------------------------------------------
; DES: T & P -- Par & Perp /////////////////
;-------------------------------------------
	;Field-aligned coordinates
	oTx = MrVar_FAC( bvec_vname, v_des_vname, 'VXB', TIME=oSn_des['TIMEVAR'] )
	
	;Rotate variables
	oP     = MrVar_Get(p_des_vname)
	oT     = MrVar_Get(t_des_vname)
	oP_fac = oTx # oP # oTx -> Transpose()
	oT_fac = oTx # oT # oTx -> Transpose()

	Obj_Destroy, oTx
	
	;Pressures
	oP_par_des   = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,0,0], /CACHE, NAME=ppar_des_vname )
	oP_perp1_des = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,1,1], /CACHE, NAME=pperp1_des_vname )
	oP_perp2_des = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,2,2], /CACHE, NAME=pperp2_des_vname )
	oP_scl_des   = (oP_par_des + oP_perp1_des + oP_perp2_des) / 3.0
	oP_scl_des  -> SetName, p_scal_des_vname
	oP_scl_des  -> Cache
	
	;Temperatures
	oT_par_des   = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,0,0], /CACHE, NAME=tpar_des_vname)
	oT_perp1_des = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,1,1], /CACHE, NAME=tperp1_des_vname)
	oT_perp2_des = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,2,2], /CACHE, NAME=tperp2_des_vname)
	oT_scl_des   = (oT_par_des + oT_perp1_des + oT_perp2_des) / 3.0
	oT_scl_des  -> SetName, t_scal_des_vname
	oT_scl_des  -> Cache
	
	Obj_Destroy, [oT_fac, oP_fac]

;-------------------------------------------
; DIS: T & P -- Par & Perp /////////////////
;-------------------------------------------
	
	;Field-aligned coordinates
	oTx = MrVar_FAC( bvec_vname, v_des_vname, 'VXB', TIME=oSn_dis['TIMEVAR'] )
	
	;Rotate variables
	oP     = MrVar_Get(p_dis_vname)
	oT     = MrVar_Get(t_dis_vname)
	oP_fac = oTx # oP # oTx -> Transpose()
	oT_fac = oTx # oT # oTx -> Transpose()

;	Obj_Destroy, oTx
	
	;Pressures
	oP_par_dis   = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,0,0], /CACHE, NAME=ppar_dis_vname )
	oP_perp1_dis = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,1,1], /CACHE, NAME=pperp1_dis_vname )
	oP_perp2_dis = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,2,2], /CACHE, NAME=pperp2_dis_vname )
	
	;Scalar T
	oP_scl_dis   = (oP_par_dis + oP_perp1_dis + oP_perp2_dis) / 3.0
	oP_scl_dis  -> SetName, p_scal_dis_vname
	oP_scl_dis  -> Cache

	;Temperatures
	oT_par_dis   = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,0,0], /CACHE, NAME=tpar_dis_vname)
	oT_perp1_dis = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,1,1], /CACHE, NAME=tperp1_dis_vname)
	oT_perp2_dis = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,2,2], /CACHE, NAME=tperp2_dis_vname)
	oT_scl_dis   = (oT_par_dis + oT_perp1_dis + oT_perp2_dis) / 3.0
	oT_scl_dis  -> SetName, t_scal_dis_vname
	oT_scl_dis  -> Cache

	Obj_Destroy, [oT_fac, oP_fac]

;-------------------------------------------
; FPI-DES: T & P -- Par & Perp /////////////
;-------------------------------------------
	;Pressure
	oP = MrVar_Get(pe_vname)
	
	;Field-aligned coordinates
	oTx = MrVar_FAC( bvec_vname, ve_vname, 'VXB', TIME=oP['TIMEVAR'] )
	
	;Rotate variables
	oP     = MrVar_Get(pe_vname)
	oT     = MrVar_Get(te_vname)
	oP_fac = oTx # oP # oTx -> Transpose()
	oT_fac = oTx # oT # oTx -> Transpose()

	Obj_Destroy, oTx
	
	;Pressures
	oPe_par   = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,0,0], /CACHE, NAME=pe_par_vname )
	oPe_perp1 = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,1,1], /CACHE, NAME=pe_perp1_vname )
	oPe_perp2 = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,2,2], /CACHE, NAME=pe_perp2_vname )
	
	;P perp
	oPe_perp = (oPe_perp1 + oPe_perp2) / 2.0
	oPe_perp  -> SetName, pe_perp_vname
	oPe_perp  -> Cache
	
	;P Scalar
	oPe_scl   = (oPe_par + oPe_perp1 + oPe_perp2) / 3.0
	oPe_scl  -> SetName, pe_scal_vname
	oPe_scl  -> Cache
	
	;Temperatures
	oTe_par   = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,0,0], /CACHE, NAME=te_par_vname)
	oTe_perp1 = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,1,1], /CACHE, NAME=te_perp1_vname)
	oTe_perp2 = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,2,2], /CACHE, NAME=te_perp2_vname)
	
	;T Perp
	oTe_perp   = (oTe_perp1 + oTe_perp2) / 2.0
	oTe_perp  -> SetName, te_perp_vname
	oTe_perp  -> Cache
	
	;T Scalar
	oTe_scl   = (oTe_par + oTe_perp1 + oTe_perp2) / 3.0
	oTe_scl  -> SetName, te_scal_vname
	oTe_scl  -> Cache
	
	Obj_Destroy, [oT_fac, oP_fac]

;-------------------------------------------
; FPI-DIS: T & P -- Par & Perp /////////////
;-------------------------------------------
	;Pressure
	oP = MrVar_Get(pi_vname)
	
	;Field-aligned coordinates
	oTx = MrVar_FAC( bvec_vname, ve_vname, 'VXB', TIME=oP['TIMEVAR'] )
	
	;Rotate variables
	oP     = MrVar_Get(pi_vname)
	oT     = MrVar_Get(ti_vname)
	oP_fac = oTx # oP # oTx -> Transpose()
	oT_fac = oTx # oT # oTx -> Transpose()

	Obj_Destroy, oTx
	
	;Pressures
	oPi_par   = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,0,0], /CACHE, NAME=pi_par_vname )
	oPi_perp1 = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,1,1], /CACHE, NAME=pi_perp1_vname )
	oPi_perp2 = MrScalarTS( oP_fac['TIMEVAR'], oP_fac[*,2,2], /CACHE, NAME=pi_perp2_vname )
	
	;P Perp
	oPi_perp  = (oPi_perp1 + oPi_perp2) / 2.0
	oPi_perp -> SetName, pi_perp_vname
	oPi_perp -> Cache
	
	;P Scalar
	oPi_scl  = (oPi_par + oPi_perp1 + oPi_perp2) / 3.0
	oPi_scl -> SetName, pi_scal_vname
	oPi_scl -> Cache
	
	;Temperatures
	oTi_par   = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,0,0], /CACHE, NAME=ti_par_vname)
	oTi_perp1 = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,1,1], /CACHE, NAME=ti_perp1_vname)
	oTi_perp2 = MrScalarTS( oT_fac['TIMEVAR'], oT_fac[*,2,2], /CACHE, NAME=ti_perp2_vname)
	
	;T Perp
	oTi_perp  = (oTi_perp1 + oTi_perp2) / 2.0
	oTi_perp -> SetName, ti_perp_vname
	oTi_perp -> Cache
	
	;T Scalar
	oTi_scl   = (oTi_par + oTi_perp1 + oTi_perp2) / 3.0
	oTi_scl -> SetName, ti_scal_vname
	oTi_scl -> Cache
	
	Obj_Destroy, [oT_fac, oP_fac]

;-------------------------------------------
; Maxwellian Entropy ///////////////////////
;-------------------------------------------
	;DES
	oSb_des = (1e-9 * oP_scl_des) / ( MrConstants('m_e') * 1e6 * oN_des )^(5.0/3.0)
	oSb_des -> SetData, alog10(oSb_des['DATA'])
	oSb_des -> SetName, sb_des_vname
	oSb_des -> Cache
	
	;DIS
	oSb_dis = (1e-9 * oP_scl_dis) / ( MrConstants('m_H') * 1e6 * oN_dis )^(5.0/3.0)
	oSb_dis -> SetData, alog10(oSb_dis['DATA'])
	oSb_dis -> SetData, oSb_dis['DATA'] + (oSb_des.min - oSb_dis.min)
	oSb_dis -> SetName, sb_dis_vname
	oSb_dis -> Cache

;-------------------------------------------
; FPI: Maxwellian Entropy //////////////////
;-------------------------------------------
	;DES
	oNe  = MrVar_Get(ne_vname)
	oSbe = (1e-9 * oPe_scl) / ( MrConstants('m_e') * 1e6 * oNe )^(5.0/3.0)
	oSbe -> SetData, alog10(oSbe['DATA'])
	oSbe -> SetName, sbe_vname
	oSbe -> Cache
	
	;DIS
	oNi  = MrVar_Get(ni_vname)
	oSbi = (1e-9 * oPi_scl) / ( MrConstants('m_H') * 1e6 * oNi )^(5.0/3.0)
	oSbi -> SetData, alog10(oSbi['DATA'])
	oSbi -> SetData, oSbi['DATA'] + (oSbe.min - oSbi.min)
	oSbi -> SetName, sbi_vname
	oSbi -> Cache

;-------------------------------------------
; Pressure /////////////////////////////////
;-------------------------------------------
	;Magnetic Pressure
	;   - 1e-9 converst to nPa
	oB  = MrVar_Get(bmag_vname)
	oPb = oB^2 / (1e9 * 2 * MrConstants('mu_0'))
	oPb -> SetName, b_press_vname
	oPb -> Cache
	
	;CALC: Total Pressure
	oPb_temp    = oPb        -> Interpol(oP_scl_des)
	oP_dis_temp = oP_scl_dis -> Interpol(oP_scl_des)
	
	oPt = oPb_temp + oP_scl_des + oP_dis_temp
	oPt -> SetName, ptot_vname
	oPt -> Cache
	
	Obj_Destroy, [oPb_temp, oP_dis_temp]
	
	;FPI Total pressure
	oPb_temp = oPb     -> Interpol(oPe_scl)
	oPi_temp = oPi_scl -> Interpol(oPe_scl)

	oPt_fpi = oPb_temp + oPe_scl + oPi_temp
	oPt_fpi -> SetName, ptot_fpi_vname
	oPt_fpi -> Cache
	
	Obj_Destroy, [oPb_temp, oPi_temp]

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------

	;B
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, coords], ' ' ) )
	
	;
	; Density
	;
	
	;DIS
	oN_dis['AXIS_RANGE'] = [ Min([oN_dis.min, oN_des.min]), Max([oN_dis.max, oN_des.max]) ]
	oN_dis['COLOR']      = 'Blue'
	oN_dis['LABEL']      = 'DIS'
	
	;DES
	oN_des['AXIS_RANGE'] = oN_dis['AXIS_RANGE']
	oN_des['COLOR']      = 'Red'
	oN_des['LABEL']      = 'DES'
	
	;
	; Entropy per particle
	;
	
	;DIS
	oSn_dis['AXIS_RANGE'] = [ Min([oSn_dis.min, oSn_dis.min]), Max([oSn_dis.max, oSn_dis.max]) ]
	oSn_dis['COLOR'] = 'Blue'
	oSn_dis['LABEL'] = 'DIS'
	oSn_dis['TITLE'] = 'S/N!C(J/K)'
	
	;DES
	oSn_des['AXIS_RANGE'] = oSn_dis['AXIS_RANGE']
	oSn_des['COLOR'] = 'Red'
	oSn_des['LABEL'] = 'DES'
	oSn_des['TITLE'] = 'S/N!C(J/K)'
	
	;
	; Entropy of Equivalent Maxwellian Distribution
	;
	
	;DIS
	oSb_dis['AXIS_RANGE'] = [ Min([oSb_dis.min, oSb_des.min]), Max([oSb_dis.max, oSb_des.max]) ]
	oSb_dis['COLOR']      = 'Blue'
	oSb_dis['LABEL']      = 'DIS'
	oSb_dis['TITLE']      = 'p/$\rho$$\up5/3$'
	
	;DES
	oSb_des['AXIS_RANGE'] = oSb_dis['AXIS_RANGE']
	oSb_des['COLOR']      = 'Red'
	oSb_des['LABEL']      = 'DES'
	oSb_des['TITLE']      = 'Log(p/$\rho$$\up5/3$)'
	
	;
	; FPI: Entropy of Equivalent Maxwellian Distribution
	;
	
	;DIS
	oSbi['AXIS_RANGE'] = [ Min([oSbi.min, oSbe.min]), Max([oSbi.max, oSbe.max]) ]
	oSbi['COLOR']      = 'Blue'
	oSbi['LABEL']      = 'DIS'
	oSbi['TITLE']      = 'p/$\rho$$\up5/3$'
	
	;DES
	oSbe['AXIS_RANGE'] = oSbi['AXIS_RANGE']
	oSbe['COLOR']      = 'Red'
	oSbe['LABEL']      = 'DES'
	oSbe['TITLE']      = 'Log(p/$\rho$$\up5/3$)'
	
	;
	; Pressure
	;
	
	;Total
	oPt['AXIS_RANGE'] = [0, oPt.max*1.1]
	oPt['COLOR']      = 'Black'
	oPt['LABEL']      = 'Total'
	oPt['TITLE']      = 'P!C(nT)'
	
	;FGM
	oPb['COLOR'] = 'Forest Green'
	oPb['LABEL'] = 'FGM'
	oPb['TITLE'] = 'P!C(nT)'
	
	;DIS
	oP_scl_dis['COLOR'] = 'Blue'
	oP_scl_dis['LABEL'] = 'DIS'
	oP_scl_dis['TITLE'] = 'P!C(nT)'
	
	;DES
	oP_scl_des['COLOR'] = 'Red'
	oP_scl_des['LABEL'] = 'DES'
	oP_scl_des['TITLE'] = 'P!C(nT)'
	
	;
	
	;
	; FPI Pressure
	;
	pirange = [ Min( [oPi_par.min, oPi_perp.min] ), $
	            Max( [oPi_par.max, oPi_perp.max] ) ]
	tirange = [ Min( [oTi_par.min, oTi_perp.min] ), $
	            Max( [oTi_par.max, oTi_perp.max] ) ]
	perange = [ Min( [oPe_par.min, oPe_perp.min] ), $
	            Max( [oPe_par.max, oPe_perp.max] ) ]
	terange = [ Min( [oTe_par.min, oTe_perp.min] ), $
	            Max( [oTe_par.max, oTe_perp.max] ) ]
	
	;Total
	oPt_fpi['AXIS_RANGE'] = [0, oPt_fpi.max*1.1]
	oPt_fpi['COLOR']      = 'Black'
	oPt_fpi['LABEL']      = 'Total'
	oPt_fpi['TITLE']      = 'P!C(nPa)'
	
	;DIS - PAR
	oPi_par['AXIS_RANGE'] = pirange
	oPi_par['COLOR']      = 'Blue'
	oPi_par['LABEL']      = 'Par'
	oPi_par['TITLE']      = 'Pi!C(nPa)'
	
	;DIS - PERP
	oPi_perp['AXIS_RANGE'] = pirange
	oPi_perp['COLOR']      = 'Red'
	oPi_perp['LABEL']      = 'Perp'
	oPi_perp['TITLE']      = 'Pi!C(nPa)'
	
	;DIS - SCALAR
	oPi_scl['COLOR'] = 'Blue'
	oPi_scl['LABEL'] = 'DIS'
	oPi_scl['TITLE'] = 'P!C(nPa)'
	
	;DES - PAR
	oPe_par['AXIS_RANGE'] = perange
	oPe_par['COLOR']      = 'Blue'
	oPe_par['LABEL']      = 'Par'
	oPe_par['TITLE']      = 'Pe!C(nPa)'
	
	;DES - PERP
	oPe_perp['AXIS_RANGE'] = perange
	oPe_perp['COLOR']      = 'Red'
	oPe_perp['LABEL']      = 'Perp'
	oPe_perp['TITLE']      = 'Pe!C(nPa)'
	
	;DES - SCALAR
	oPe_scl['COLOR'] = 'Red'
	oPe_scl['LABEL'] = 'DES'
	oPe_scl['TITLE'] = 'P!C(nPa)'
	
	;
	; FPI Temperature
	;
	
	;DIS - PAR
	oTi_par['AXIS_RANGE'] = tirange
	oTi_par['COLOR']      = 'Blue'
	oTi_par['LABEL']      = 'Par'
	oTi_par['TITLE']      = 'Ti!C(eV)'
	
	;DIS - PERP
	oTi_perp['AXIS_RANGE'] = tirange
	oTi_perp['COLOR']      = 'Red'
	oTi_perp['LABEL']      = 'Perp'
	oTi_perp['TITLE']      = 'Ti!C(eV)'
	
	;DES - PAR
	oTe_par['AXIS_RANGE'] = terange
	oTe_par['COLOR']      = 'Blue'
	oTe_par['LABEL']      = 'Par'
	oTe_par['TITLE']      = 'Te!C(eV)'
	
	;DES - PERP
	oTe_perp['AXIS_RANGE'] = terange
	oTe_perp['COLOR']      = 'Red'
	oTe_perp['LABEL']      = 'Perp'
	oTe_perp['TITLE']      = 'Te!C(eV)'
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------

	win = MrVar_PlotTS(  [ b_vname, n_des_vname, sn_des_vname, sbe_vname, ptot_fpi_vname, $
	                       pi_par_vname, pe_par_vname, ti_par_vname, te_par_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )
	
	win = MrVar_OPlotTS( [n_des_vname, sn_des_vname, sbe_vname, ptot_fpi_vname], $
	                     [n_dis_vname, sn_dis_vname, sbi_vname, b_press_vname] )
	win = MrVar_OPlotTS( ptot_fpi_vname, pi_scal_vname )
	win = MrVar_OPlotTS( ptot_fpi_vname, pe_scal_vname )
	win = MrVar_OPlotTS( pi_par_vname, pi_perp_vname )
	win = MrVar_OPlotTS( pe_par_vname, pe_perp_vname )
	win = MrVar_OPlotTS( ti_par_vname, ti_perp_vname )
	win = MrVar_OPlotTS( te_par_vname, te_perp_vname )

	win[0] -> SetLayout, [1,1]
	win -> TrimLayout
	win.oxmargin = [15,10]

;-------------------------------------------
; Save the Figure //////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
			MrPrintF, 'LogText', 'Saving file to: "' + output_dir + '".'
		ENDIF
		
		;Save the figure
		fname = StrJoin([sc, 'fpi', mode, level, 'entropy'], '_')
		fname = FilePath(fname, ROOT_DIR=output_dir)
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF
	
;-------------------------------------------
; Finish ///////////////////////////////////
;-------------------------------------------

	win -> Refresh
	RETURN, win
END