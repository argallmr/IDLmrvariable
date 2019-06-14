; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_KEntropy_Grad
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
;   Calculate moments of the distribution function and plot them against the official
;   FPI L2 dataset. The moments calculation takes into account the FPI internal photo-
;   electron model, but the method of integration is different.
;
;       1. Bxyz, |B|
;       2. n: Density
;       3. N: Total number of particles
;       4. n: Maximum number of particles in any bin
;       5. sn, sBn: Entropy per particle
;       6. Mnbar = sn - sBn: Non-maxwellian entropy per particle
;       7. Q
;       8. J.E'
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
;       2018/05/01  -   Written by Matthew Argall
;       2018/06/07  -   Use calculated moments for sB so that the comparison with
;                           s is fair. Fix units error in calculation of sB - MRA
;-
FUNCTION MrMMS_Plot_FPI_KEntropy_Grad, sc, mode, species, $
COORDS=coords, $
EPHDESC=ephdesc, $
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
		IF Obj_Valid(win) THEN Obj_Destroy, win
		RETURN, Obj_New()
	ENDIF
	
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(ephdesc)   EQ 0 THEN ephdesc   = 'epht89d'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	
	q = MrConstants('q')

;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	fgm_coords  = (coords EQ 'dbcs' || coords EQ 'dsl') ? 'dmpa' : coords
	
	fpi_instr  = 'd' + species + 's'
	fpi_mode   = mode EQ 'brst' ? mode : 'fast'
	fpi_coords = (coords EQ 'dmpa' || coords EQ 'dsl') ? 'dbcs' : coords
	
	edp_mode   = mode EQ 'brst' ? mode : 'fast'
	edp_coords = (coords EQ 'dbcs' || coords EQ 'dmpa') ? 'dsl' : coords
	
	eph_mode = mode

;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, mode, $
		                     INSTR  = fgm_instr, $
		                     LEVEL  = level, $
		                     VARFORMAT = '*b_' + fgm_coords + '*', $
		                     SUFFIX = suffix
		
		;EDP
		MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = '*dce_' + edp_coords + '*'
		
		;DIS
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = 'dis-moms', $
		                     VARFORMAT = '*' + ['numberdensity', 'bulkv_'+fpi_coords, 'prestensor_'+fpi_coords] $
		                                 + '_' + fpi_mode
		
		;DES
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = 'des-moms', $
		                     VARFORMAT = '*' + ['numberdensity', 'bulkv_'+fpi_coords, 'prestensor_'+fpi_coords] $
		                                 + '_' + fpi_mode
	ENDIF
	
	MrMMS_FPI_Load_Entropy, sc, mode, species, $
	                        NO_LOAD  = no_load, $
	                        VARNAMES = s_varnames

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source Names
	b_vname     = StrJoin( [sc, fgm_instr, 'b',    coords, mode, level], '_' )
	bvec_vname  = StrJoin( [sc, fgm_instr, 'bvec', coords, mode, level], '_' )
	bmag_vname  = StrJoin( [sc, fgm_instr, 'bmag', coords, mode, level], '_' )
	e_vname     = StrJoin( [sc, 'edp', 'dce', edp_coords, edp_mode, level], '_')
	ni_vname    = StrJoin( [sc, 'dis', 'numberdensity',             fpi_mode], '_' )
	vi_vname    = StrJoin( [sc, 'dis', 'bulkv',         fpi_coords, fpi_mode], '_' )
	pi_vname    = StrJoin( [sc, 'dis', 'prestensor',    fpi_coords, fpi_mode], '_' )
	ne_vname    = StrJoin( [sc, 'des', 'numberdensity',             fpi_mode], '_' )
	ve_vname    = StrJoin( [sc, 'des', 'bulkv',         fpi_coords, fpi_mode], '_' )
	pe_vname    = StrJoin( [sc, 'des', 'prestensor',    fpi_coords, fpi_mode], '_' )
	
	IF species EQ 'i' THEN BEGIN
		n_vname  = ni_vname
		v_vname  = vi_vname
		p_vname  = pi_vname
	ENDIF ELSE BEGIN
		n_vname  = ne_vname
		v_vname  = ve_vname
		p_vname  = pe_vname
	ENDELSE
	
	;Derived Names
;	n_vname       = s_varnames[0]
;	pt_vname      = s_varnames[4]
	cts_max_vname   = s_varnames[5]
	cts_tot_vname   = s_varnames[6]
	mcts_max_vname  = s_varnames[7]
	mcts_tot_vname  = s_varnames[8]
	sn_vname        = s_varnames[11]
	sbn_vname       = s_varnames[15]
	mfrac_vname     = s_varnames[17]
	sn_tot_vname    = s_varnames[24]
	sbn_tot_vname   = s_varnames[28]
	mfrac_tot_vname = s_varnames[30]
	
	q_vname      = StrJoin([sc, fpi_instr, 'agyrotropy-factor',                fpi_mode], '_')
	j_vname      = StrJoin([sc, 'fpi', 'j',                                    fpi_mode], '_')
	eprime_vname = StrJoin([sc, 'edp', 'eprime', edp_coords, edp_mode, level], '_')
	jdote_vname  = StrJoin([sc, 'edp', 'jdote',  edp_coords, edp_mode, level], '_')

;-------------------------------------------
; Agyrotropy Parameter /////////////////////
;-------------------------------------------
	oQ = MrVar_Pres_QFactor(bvec_vname, p_vname, /CACHE, NAME=q_vname)

;-------------------------------------------
; Current Density //////////////////////////
;-------------------------------------------
	;Pull data
	oVi = MrVar_Get(vi_vname)
	oVe = MrVar_Get(ve_vname)
	oNe = MrVar_Get(ne_vname)
	
	;Interpolate to DES times
	oVi_des = oVi -> Interpol(oVe)
	
	;Calculate current density
	;   - 1e15 converts to uA/m^2
	oJ = 1e15 * q * oNe * (oVi_des - oVe)
	Obj_Destroy, oVi_des
	
	;Attributes
	oJ -> SetName, j_vname
	oJ -> Cache
	oJ['CATDESC']       = 'Current density derived from moments of the distribution function.'
	oJ['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oJ['LABEL']         = 'J$\down' + ['X', 'Y', 'Z'] + '$'
	oJ['TITLE']         = 'J!C($\mu$A/m^2)'
	oJ['UNITS']         = 'uA/m^2'
	oJ['SI_CONVERSION'] = '1e-6>A/m^2'

;-------------------------------------------
; Electric Field in e- Frame ///////////////
;-------------------------------------------
	oE      = MrVar_Get(e_vname)
	oB      = MrVar_Get(bvec_vname)
	oB_des  = oB -> Interpol(oVe)
	oE_des  = oE -> Interpol(oVe)
	oVexB   = MrVar_E_VxB(oVe, oB_des)
	oEprime = oE_des - oVexB
	Obj_Destroy, [oE_des, oB_des, oVexB]
	
	oEprime -> SetName, eprime_vname
	oEprime -> Cache
	oEprime['CATDESC']       = 'Electric field in electron bulk velocity frame.'
	oEprime['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oEprime['LABEL']         = "E'$\down" + ['X', 'Y', 'Z'] + '$'
	oEprime['TITLE']         = 'E!C(mV/m)'
	oEprime['UNITS']         = 'mV/m'
	oEprime['SI_CONVERSION'] = '1e-3>V/m'

;-------------------------------------------
; Electron Frame Dissipation Measure ///////
;-------------------------------------------
	;J.E = mV/m * uA/m^2 * 1e-3 V/mV * 1e-6 uA/A
	;    = V/m * A/m^2 * 1e-9
	;    = kg m^2 / (C s^2 m) * C / (s m^2) * 1e-9
	;    = kg * C/C * m^2/(m m^2) * 1/(s^2 s) * 1e-9
	;    = kg / (m s^3) * 1e-9
	;    ---> W/m^3 = Nm/s = kg m^2 / (s^3 m^3) = kg / (m s^3)
	;    = W / m^3 * 1e-9 * 1e9 nW/W
	;    = nW / m^3
	oJdotE = oJ -> Dot(oEprime)
	oJdotE -> SetName, jdote_vname
	oJdotE -> Cache
	oJdotE['CATDESC']       = 'Electron frame dissipation measure.'
	oJdotE['TITLE']         = "J.E'!C(nW/m^3)"
	oJdotE['UNITS']         = 'nW/m^3'
	oJdotE['SI_CONVERSION'] = '1e9>W/m^3'

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------

	;B
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, fpi_instr, mode, level, coords], ' ' ) )
	
	;Entropy per particle
	oSn  = MrVar_Get(sn_vname)
	oSbn = MrVar_Get(sbn_vname)
	
	oSn['AXIS_RANGE'] = [ Min([oSn.min, oSbn.min]), Max([oSn.max, oSbn.max]) ]
	oSn['COLOR']      = 'Black'
	
	oSbn['AXIS_RANGE'] = oSn['AXIS_RANGE']
	oSbn['COLOR']      = 'Blue'
	
	;Maximum Counts
	oNmax = MrVar_Get(cts_max_vname)
	oMNmax = MrVar_Get(mcts_max_vname)
	
	oNmax['AXIS_RANGE'] = [ Min([oNmax.min, oMNmax.min]), Max([oNmax.max, oMNmax.max]) ]
	oNmax['COLOR']      = 'Black'
	
	oMNmax['AXIS_RANGE'] = oNmax['AXIS_RANGE']
	oMNmax['COLOR']      = 'Blue'
	
	;Total Counts
	oNtot  = MrVar_Get(cts_tot_vname)
	oMNtot = MrVar_Get(mcts_tot_vname)
	
	oNtot['AXIS_RANGE'] = [ Min([oNtot.min, oMNtot.min]), Max([oNtot.max, oMNtot.max]) ]
	oNtot['COLOR']      = 'Black'
	
	oMNtot['AXIS_RANGE'] = oNtot['AXIS_RANGE']
	oMNtot['COLOR']      = 'Blue'
	
	;Non-Maxwellian Entropy
	oMfrac    = MrVar_Get(mfrac_vname)
	oMfrac_tot = MrVar_Get(mfrac_tot_vname)
	
	oMfrac['AXIS_RANGE'] = [ Min([oMfrac.min, oMfrac_tot.min]), Max([oMfrac.max, oMfrac_tot.max]) ]
	oMfrac['COLOR']      = 'Blue'
	
	oMfrac_tot['AXIS_RANGE'] = oMfrac['AXIS_RANGE']
	oMfrac_tot['COLOR']      = 'Black'
	
	;Boltzmann Entropy per Particle
	oSn_tot  = MrVar_Get(sn_tot_vname)
	oSbn_tot = MrVar_Get(sbn_tot_vname)
	
	oSn_tot['AXIS_RANGE'] = [ Min([oSn_tot.min, oSbn_tot.min]), Max([oSn_tot.max, oSbn_tot.max]) ]
	oSn_tot['COLOR']      = 'Black'
	
	oSbn_tot['AXIS_RANGE'] = oSbn['AXIS_RANGE']
	oSbn_tot['COLOR']      = 'Blue'
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [b_vname, n_vname, cts_max_vname, cts_tot_vname, sn_vname, sn_tot_vname, mfrac_vname, q_vname, jdote_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )
	win = MrVar_OPlotTS( cts_max_vname, mcts_max_vname )
	win = MrVar_OPlotTS( cts_tot_vname, mcts_tot_vname )
	win = MrVar_OPlotTS( sn_vname, sbn_vname )
	win = MrVar_OPlotTS( sn_tot_vname, sbn_tot_vname )
	win = MrVar_OPlotTS( mfrac_vname, mfrac_tot_vname )

	win[0] -> SetLayout, [1,1]
	win -> TrimLayout
	win.oxmargin = [15,9]

;-------------------------------------------
; Save the File ////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
			MrPrintF, 'LogText', 'Saving file to: "' + output_dir + '".'
		ENDIF
		
		;File name
		fname = StrJoin([sc, fpi_instr, fpi_mode, level, 'kentropy'], '_')
		fname = FilePath(fname, ROOT_DIR=output_dir)
		
		;Save
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	win -> Refresh
	RETURN, win
END