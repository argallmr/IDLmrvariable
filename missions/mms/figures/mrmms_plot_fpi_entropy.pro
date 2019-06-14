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
	
	m_e  = MrConstants('m_e')
	m_i  = MrConstants('m_H')
	q    = MrConstants('q')
	kB   = MrConstants('k_b')
	J2eV = MrConstants('J2eV')

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
	sb_des_vname     = StrJoin([sc, des_instr, 'entropydensity', 'maxwell',    'calc', mode], '_')
	sbn_des_vname    = StrJoin([sc, des_instr, 'entropy',        'maxwell',    'calc', mode], '_')
	mbar_des_vname   = StrJoin([sc, des_instr, 'entropydensity', 'nonmaxwell', 'calc', mode], '_')
	mnbar_des_vname  = StrJoin([sc, des_instr, 'entropy',        'nonmaxwell', 'calc', mode], '_')
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
	sb_dis_vname     = StrJoin([sc, dis_instr, 'entropydensity', 'maxwell',   'calc', mode], '_')
	sbn_dis_vname    = StrJoin([sc, dis_instr, 'entropy',        'maxwell',   'calc', mode], '_')
	mbar_dis_vname   = StrJoin([sc, dis_instr, 'entropydensity', 'nonmaxwell', 'calc', mode], '_')
	mnbar_dis_vname  = StrJoin([sc, dis_instr, 'entropy',        'nonmaxwell', 'calc', mode], '_')
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
; Compute Moments //////////////////////////
;-------------------------------------------
	;DES
	oDist = MrDist4D(f_des_vname, VSC=scpot_vname, SPECIES='e')
	oN_des = oDist -> Density(/CACHE, NAME=n_des_vname)
	oS_des = oDist -> Entropy()
	oP_des = oDist -> Pressure(/CACHE, NAME=p_des_vname)
	oT_des = oDist -> Temperature(/CACHE, NAME=t_des_vname)
	oV_des = oDist -> Velocity(/CACHE, NAME=v_des_vname)
	
	;DIS
	oDist = MrDist4D(f_dis_vname, VSC=scpot_vname, SPECIES='H')
	oN_dis = oDist -> Density(/CACHE, NAME=n_dis_vname)
	oS_dis = oDist -> Entropy()
	oP_dis = oDist -> Pressure(/CACHE, NAME=p_dis_vname)
	oT_dis = oDist -> Temperature(/CACHE, NAME=t_dis_vname)
	oV_dis = oDist -> Velocity(/CACHE, NAME=v_dis_vname)
	
	;Density
	oN_dis['AXIS_RANGE'] = [ Min([oN_dis.min, oN_des.min]), Max([oN_dis.max, oN_des.max]) ]
	oN_dis['COLOR']      = 'Blue'
	oN_dis['LABEL']      = 'DIS'
	
	oN_des['AXIS_RANGE'] = oN_dis['AXIS_RANGE']
	oN_des['COLOR']      = 'Red'
	oN_des['LABEL']      = 'DES'

;-------------------------------------------
; Entropy //////////////////////////////////
;-------------------------------------------
	
	;DIS
	oS_dis *= J2eV
	oS_dis -> SetName, s_dis_vname
	oS_dis -> Cache
	
	;DES
	oS_des *= J2eV
	oS_des -> SetName, s_des_vname
	oS_des -> Cache
	
	oS_dis['COLOR']      = 'Red'
	oS_dis['LABEL']      = 's'
	oS_dis['TITLE']      = 'si!C(eV/K/m$\up3$)'
	oS_dis['UNITS']      = 'eV/K/m^3 ln(s^3/m^6)'
	
	oS_des['COLOR']      = 'Red'
	oS_des['LABEL']      = 's'
	oS_des['TITLE']      = 'se!C(eV/K/m$\up3$)'
	oS_des['UNITS']      = 'eV/K/m^3 ln(s^3/m^6)'
	
;-------------------------------------------
; Entropy Per Particle /////////////////////
;-------------------------------------------

	;Calculations
	;   - Match minimum values so they are on the same scale
	oSn_des = oS_des / (1e6 * oN_des)
	oSn_dis = oS_dis / (1e6 * oN_dis)
	
	;DIS
	oSn_dis -> SetName, sn_dis_vname
	oSn_dis -> Cache
	oSn_dis['COLOR'] = 'Red'
	oSn_dis['LABEL'] = 's/N'
	oSn_dis['TITLE'] = 'si/N!C(eV/K)'
	
	;DES
	oSn_des -> SetName, sn_des_vname
	oSn_des -> Cache
	oSn_des['COLOR'] = 'Red'
	oSn_des['LABEL'] = 's/N'
	oSn_des['TITLE'] = 'se/N!C(eV/K)'

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
	oP_scl_des['COLOR'] = 'Red'
	oP_scl_des['LABEL'] = 'DES'
	oP_scl_des['TITLE'] = 'P!C(nT)'
	
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
	oP_scl_dis['COLOR'] = 'Blue'
	oP_scl_dis['LABEL'] = 'DIS'
	oP_scl_dis['TITLE'] = 'P!C(nT)'

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
	
	
	;
	; ATTRIBUTES
	;
	perange = [ Min( [oPe_par.min, oPe_perp.min] ), $
	            Max( [oPe_par.max, oPe_perp.max] ) ]
	terange = [ Min( [oTe_par.min, oTe_perp.min] ), $
	            Max( [oTe_par.max, oTe_perp.max] ) ]
	
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
	
	;
	; ATTRIBUTES
	;
	pirange = [ Min( [oPi_par.min, oPi_perp.min] ), $
	            Max( [oPi_par.max, oPi_perp.max] ) ]
	tirange = [ Min( [oTi_par.min, oTi_perp.min] ), $
	            Max( [oTi_par.max, oTi_perp.max] ) ]
	
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

;-------------------------------------------
; Maxwellian Entropy ///////////////////////
;-------------------------------------------
	;Integration of the distribution function is performed in MKS units. As such,
	;the logarithm should be also.

	;DES
	;   - Entroy per unit volume
	;   - Constant factors convert to MKS
	oSb_des = 3.0/2.0 * kB * J2eV * 1e6 * oN_des $
	          * ALog( 2 * !pi * 1e-11 / m_e * oP_scl_des['DATA'] / oN_des['DATA']^(5.0/3.0) + 1 )
	
	;DIS
	oSb_dis = 3.0/2.0 * kB * J2eV * 1e6 * oN_dis $
	          * ALog( 2 * !pi * 1e-11 / m_i * oP_scl_dis['DATA'] / oN_dis['DATA']^(5.0/3.0) + 1 )
	
	;DIS
	oSb_dis -> SetName, sb_dis_vname
	oSb_dis -> Cache
	oSb_dis['COLOR']      = 'Black'
	oSb_dis['LABEL']      = 's$\downB$'
	oSb_dis['TITLE']      = 'si$\downB$!C(eV/K/m$\up3$)'
	oSb_dis['UNITS']      = 'eV/K/m^3 ln(J m^2/kg)'
	
	;DES
	oSb_des -> SetName, sb_des_vname
	oSb_des -> Cache
	oSb_des['COLOR']      = 'Black'
	oSb_des['LABEL']      = 's$\downB$'
	oSb_des['TITLE']      = 'se$\downB$!C(eV/K/m$\up3$)'
	oSb_des['UNITS']      = 'eV/K/m^3 ln(J m^2/kg)'

;-------------------------------------------
; Maxwellian Entropy per Particle //////////
;-------------------------------------------
	;DIS
	oSbn_dis = oSb_dis / (1e6 * oN_dis)
	oSbn_dis -> SetName, sbn_dis_vname
	oSbn_dis -> Cache
	oSbn_dis['COLOR']      = 'Black'
	oSbn_dis['LABEL']      = 's$\downB$/N'
	oSbn_dis['TITLE']      = 'si$\downB$/N!C(eV/K/m$\up3$)'
	oSbn_dis['UNITS']      = 'eV/K/m^3 ln(J m^2/kg)'
	
	;DES
	oSbn_des = oSb_des / (1e6 * oN_des)
	oSbn_des -> SetName, sbn_des_vname
	oSbn_des -> Cache
	oSbn_des['COLOR']      = 'Black'
	oSbn_des['LABEL']      = 's$\downB$/N'
	oSbn_des['TITLE']      = 'si$\downB$/N!C(eV/K)'
	oSbn_des['UNITS']      = 'eV/K ln(J m^2/kg)'

;-------------------------------------------
; Non-Maxwellian Entropy ///////////////////
;-------------------------------------------
	
	;MBAR
	oMbar_des  = oS_des  - oS_des
	oMnbar_des = oSn_des - oSbn_des
	oMbar_dis  = oS_dis  - oSb_dis
	oMnbar_dis = oSn_dis - oSbn_dis
	
	oMbar_des -> SetName, mbar_des_vname
	oMbar_des -> Cache
	oMbar_des['COLOR'] = 'Red'
	oMbar_des['LABEL'] = 'DES'
	oMbar_des['TITLE'] = 'M!C(eV/K/m$\up3$)'
	oMbar_des['UNITS'] = 'eV/K/m^3 ln(J m^2/kg)'
	
	oMnbar_des -> SetName, mnbar_des_vname
	oMnbar_des -> Cache
	oMnbar_des['COLOR'] = 'Red'
	oMnbar_des['LABEL'] = 'DES'
	oMnbar_des['TITLE'] = 'M/N!C(eV/K)'
	oMnbar_des['UNITS'] = 'eV/K ln(J m^2/kg)'
	
	oMbar_des -> SetName, mbar_dis_vname
	oMbar_dis -> Cache
	oMbar_dis['COLOR'] = 'Blue'
	oMbar_dis['LABEL'] = 'DIS'
	oMbar_dis['TITLE'] = 'M!C(eV/K/m$\up3$)'
	oMbar_dis['UNITS'] = 'eV/K/m^3 ln(J m^2/kg)'
	
	oMnbar_dis -> SetName, mnbar_dis_vname
	oMnbar_dis -> Cache
	oMnbar_dis['COLOR'] = 'Blue'
	oMnbar_dis['LABEL'] = 'DIS'
	oMnbar_dis['TITLE'] = 'M/N!C(eV/K)'
	oMnbar_dis['UNITS'] = 'eV/K ln(J m^2/kg)'
	

;-------------------------------------------
; FPI: Maxwellian Entropy //////////////////
;-------------------------------------------
	;DES
	oNe  = MrVar_Get(ne_vname)
	oSbe = 3.0/2.0 * kB * J2eV * 1e6 * oN_des $
	       * ALog( 2 * !pi * 1e-11 / m_e * oPe_scl['DATA'] / oNe['DATA']^(5.0/3.0) + 1 )
	
	;DIS
	oNi  = MrVar_Get(ni_vname)
	oSbi = 3.0/2.0 * kB * J2eV * 1e6 * oNi $
	       * ALog( 2 * !pi * 1e-11 / m_i * oPi_scl['DATA'] / oNi['DATA']^(5.0/3.0) + 1 )
	
	;DIS
	oSbi -> SetName, sbi_vname
	oSbi -> Cache
	oSbi['AXIS_RANGE'] = [ Min([oSbi.min, oSbe.min]), Max([oSbi.max, oSbe.max]) ]
	oSbi['COLOR']      = 'Blue'
	oSbi['LABEL']      = 'DIS'
	oSbi['TITLE']      = 's$\downB$!C(eV/K/m$\up3$)'
	oSbi['UNITS']      = 'eV/K/m^3 ln(J m^2/kg)'
	
	;DES
	oSbe -> SetName, sbe_vname
	oSbe -> Cache
	oSbe['AXIS_RANGE'] = oSbi['AXIS_RANGE']
	oSbe['COLOR']      = 'Red'
	oSbe['LABEL']      = 'DES'
	oSbe['TITLE']      = 's$\downB$!C(eV/K/m$\up3$)'
	oSbe['UNITS']      = 'eV/K/m^3 ln(J m^2/kg)'

;-------------------------------------------
; Pressure /////////////////////////////////
;-------------------------------------------
	;Magnetic Pressure
	;   - 1e-9 converst to nPa
	oB  = MrVar_Get(bmag_vname)
	oPb = oB^2 / (1e9 * 2 * MrConstants('mu_0'))
	oPb -> SetName, b_press_vname
	oPb -> Cache
	oPb['COLOR'] = 'Forest Green'
	oPb['LABEL'] = 'FGM'
	oPb['TITLE'] = 'P!C(nT)'
	
	;CALC: Total Pressure
	oPb_temp    = oPb        -> Interpol(oP_scl_des)
	oP_dis_temp = oP_scl_dis -> Interpol(oP_scl_des)
	
	oPt = oPb_temp + oP_scl_des + oP_dis_temp
	oPt -> SetName, ptot_vname
	oPt -> Cache
	oPt['AXIS_RANGE'] = [0, oPt.max*1.1]
	oPt['COLOR']      = 'Black'
	oPt['LABEL']      = 'Total'
	oPt['TITLE']      = 'P!C(nT)'
	
	Obj_Destroy, [oPb_temp, oP_dis_temp]
	
	;FPI Total pressure
	oPb_temp = oPb     -> Interpol(oPe_scl)
	oPi_temp = oPi_scl -> Interpol(oPe_scl)

	oPt_fpi = oPb_temp + oPe_scl + oPi_temp
	oPt_fpi -> SetName, ptot_fpi_vname
	oPt_fpi -> Cache
	oPt_fpi['AXIS_RANGE'] = [0, oPt_fpi.max*1.1]
	oPt_fpi['COLOR']      = 'Black'
	oPt_fpi['LABEL']      = 'Total'
	oPt_fpi['TITLE']      = 'P!C(nPa)'
	
	Obj_Destroy, [oPb_temp, oPi_temp]

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------

	;B
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, coords], ' ' ) )
	
	oSn_dis['AXIS_RANGE']    = [ Min([oSn_dis.min, oSbn_dis.min]), Max([oSn_dis.max, oSbn_dis.max]) ]
	oSbn_dis['AXIS_RANGE']   = oSn_dis['AXIS_RANGE']
	oSn_des['AXIS_RANGE']    = [ Min([oSn_des.min, oSbn_des.min]), Max([oSn_des.max, oSbn_des.max]) ]
	oSbn_des['AXIS_RANGE']   = oSn_des['AXIS_RANGE']
	oMnbar_des['AXIS_RANGE'] = [ Min([oMnbar_des.min, oMnbar_dis.min]), Max([oMnbar_des.max, oMnbar_dis.max]) ]
	oMnbar_dis['AXIS_RANGE'] = oMnbar_des['AXIS_RANGE']
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS(  [ b_vname, n_des_vname, sn_des_vname, sn_dis_vname, mnbar_des_vname, ptot_fpi_vname, $
	                       pi_par_vname, pe_par_vname, ti_par_vname, te_par_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )
	
	win = MrVar_OPlotTS( [n_des_vname, sn_des_vname, sn_dis_vname, mnbar_des_vname, ptot_fpi_vname], $
	                     [n_dis_vname, sbn_des_vname, sbn_dis_vname, mnbar_dis_vname, b_press_vname] )
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