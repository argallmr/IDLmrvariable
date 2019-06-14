; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_FPI_KEntropy
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
;       1. Bxyz, |B| at barycenter
;       2. Density                    (4sc)
;       3. Entropy density            (4sc)
;       4. Maxwellian entropy density (4sc)
;       5. M-bar = s - sB             (4sc)
;       6. Grad(M-bar)
;       7. Q                          (4sc)
;       8. J.E'                       (4sc)
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
;       2018/05/20  -   Written by Matthew Argall
;       2018/06/07  -   Use calculated moments for sB so that the comparison with
;                           s is fair. Fix units error in calculation of sB - MRA
;-
FUNCTION MrMMS_Plot_4sc_FPI_KEntropy, mode, species, $
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
	
	q     = MrConstants('q')
	kB    = MrConstants('k_B')
	m_e   = MrConstants('m_e')
	m_i   = MrConstants('m_H')
	nabla = '!9' + String(71B) + '!X'
	J2eV  = MrConstants('J2eV')
	colors = ['Black', 'Red', 'Forest Green', 'Blue']

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
; Variable Names ///////////////////////////
;-------------------------------------------
	sc  = 'mms' + ['1', '2', '3', '4']

	;Source Names
	b_vnames     = sc + '_' + StrJoin( [fgm_instr, 'b',    coords, mode, level], '_' )
	bvec_vnames  = sc + '_' + StrJoin( [fgm_instr, 'bvec', coords, mode, level], '_' )
	bmag_vnames  = sc + '_' + StrJoin( [fgm_instr, 'bmag', coords, mode, level], '_' )
	e_vnames     = sc + '_' + StrJoin( ['edp', 'dce', edp_coords, edp_mode, level], '_')
	scpot_vnames = sc + '_' + StrJoin( ['edp', 'scpot',           edp_mode, level], '_')
	fi_vnames    = sc + '_' + StrJoin( ['dis', 'dist',           fpi_mode], '_' )
	ni_vnames    = sc + '_' + StrJoin( ['dis', 'numberdensity',             fpi_mode], '_' )
	vi_vnames    = sc + '_' + StrJoin( ['dis', 'bulkv',         fpi_coords, fpi_mode], '_' )
	pi_vnames    = sc + '_' + StrJoin( ['dis', 'prestensor',    fpi_coords, fpi_mode], '_' )
	fe_vnames    = sc + '_' + StrJoin( ['des', 'dist',           fpi_mode], '_' )
	ne_vnames    = sc + '_' + StrJoin( ['des', 'numberdensity',             fpi_mode], '_' )
	ve_vnames    = sc + '_' + StrJoin( ['des', 'bulkv',         fpi_coords, fpi_mode], '_' )
	pe_vnames    = sc + '_' + StrJoin( ['des', 'prestensor',    fpi_coords, fpi_mode], '_' )
	r_vnames     = sc + '_' + StrJoin( ['mec', 'r', coords], '_' )
	
	IF species EQ 'i' THEN BEGIN
		f_vnames = fi_vnames
		n_vnames = ni_vnames
		v_vnames = vi_vnames
		p_vnames = pi_vnames
		mass     = m_i
		type     = 'H'
	ENDIF ELSE BEGIN
		f_vnames = fe_vnames
		n_vnames = ne_vnames
		v_vnames = ve_vnames
		p_vnames = pe_vnames
		mass     = m_e
		type     = 'e'
	ENDELSE
	
	;Derived names
	b_bc_vname    = StrJoin( ['mms', fgm_instr, 'bvec', coords, 'barycenter', mode, level], '_' )
	e_bc_vname    = StrJoin( ['mms', 'edp', 'dce', edp_coords, 'barycenter', edp_mode, level], '_')

	sn_vnames     = sc + '_' + StrJoin([fpi_instr, 'sn',    fpi_mode], '_')
	sbn_vnames    = sc + '_' + StrJoin([fpi_instr, 'sbn', fpi_mode], '_')
	mnbar_vnames  = sc + '_' + StrJoin([fpi_instr, 'mnbar',          fpi_mode], '_')
	gradmn_vname  = StrJoin(['mms', fpi_instr, 'gradmn',             fpi_mode], '_')
	q_vnames      = sc + '_' + StrJoin([fpi_instr, 'agyrotropy-factor',          fpi_mode], '_')
	j_vnames      = sc + '_' + StrJoin(['fpi', 'j',                              fpi_mode], '_')
	jcurl_vname   = StrJoin( ['mms', fgm_instr, 'j', coords, 'curlometer', mode, level], '_' )
	eprime_vnames = sc + '_' + StrJoin(['edp', 'eprime', edp_coords, edp_mode, level], '_')
	jdote_vnames  = sc + '_' + StrJoin(['edp', 'jdote',  edp_coords, edp_mode, level], '_')
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, '', mode, $
		                     INSTR  = fgm_instr, $
		                     LEVEL  = level, $
		                     VARFORMAT = '*b_' + fgm_coords + '_' + mode + '*', $
		                     SUFFIX = suffix
		
		;EDP
		MrMMS_Load_Data, '', 'edp', edp_mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = '*dce_' + edp_coords + '*'
		
		;FPI-dist
		MrMMS_FPI_Load_Dist3D, 'mms1', mode, species, /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms2', mode, species, /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms3', mode, species, /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms4', mode, species, /APPLY_MODEL
		
		;DIS-Moms
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = 'dis-moms', $
		                     VARFORMAT = [ '*numberdensity*', '*bulkv_'+coords+'*', $
		                                   '*pres*'+coords+'*' ]
		
		;DES-Moms
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = 'des-moms', $
		                     VARFORMAT = [ '*numberdensity*', '*bulkv_'+coords+'*', $
		                                   '*pres*'+coords+'*' ]
		
		;Spacecraft potential
		MrMMS_Load_Data, '', 'edp', edp_mode, level, $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
		
		;Ephemeris
		IF mode EQ 'brst' THEN BEGIN
			;BRST data is not always available
			IF (MrMMS_Get_FileNames(sc, 'mec', mode, level, OPTDESC=ephdesc))[0] EQ '' $
				THEN eph_mode = 'srvy'
		ENDIF
		MrMMS_Load_Data, '', 'mec', eph_mode, level, $
		                 OPTDESC   = ephdesc, $
		                 VARFORMAT = '*r_'+coords
	ENDIF
	
;-------------------------------------------
; Barycenter ///////////////////////////////
;-------------------------------------------
	;B-Field
	oB1 = MrVar_Get(bvec_vnames[0])
	oB2 = MrVar_Get(bvec_vnames[1])
	oB3 = MrVar_Get(bvec_vnames[2])
	oB4 = MrVar_Get(bvec_vnames[3])
	oB2 = oB2 -> Interpol(oB1)
	oB3 = oB3 -> Interpol(oB1)
	oB4 = oB4 -> Interpol(oB1)
	
	oBbc = (oB1 + oB2 + oB3 + oB4) / 4.0
	oBbc -> SetName, b_bc_vname
	oBbc -> Cache
	oBbc['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oBbc['LABEL']         = 'B$\down' + ['X', 'Y', 'Z'] + '$'
	oBbc['SI_CONVERSION'] = '1e-9>T'
	oBbc['TITLE']         = 'B$\downBC$!C(nT)'
	oBbc['UNITS']         = 'nT'
	
	;E-Field
	oE1 = MrVar_Get(e_vnames[0])
	oE2 = MrVar_Get(e_vnames[1])
	oE3 = MrVar_Get(e_vnames[2])
	oE4 = MrVar_Get(e_vnames[3])
	oE2 = oE2 -> Interpol(oE1)
	oE3 = oE3 -> Interpol(oE1)
	oE4 = oE4 -> Interpol(oE1)
	
	oEbc = (oE1 + oE2 + oE3 + oE4) / 4.0
	oEbc -> SetName, e_bc_vname
	oEbc -> Cache
	oEbc['TITLE']         = 'E$\downBC$!C(mV/m)'
	oEbc['UNITS']         = 'mV/m'
	oEbc['SI_CONVERSION'] = '1e-3>V/m'

;-------------------------------------------
; Loop Over Each Spacecraft ////////////////
;-------------------------------------------
	FOR i = 0, N_Elements(sc) - 1 DO BEGIN
	;-------------------------------------------
	; Density //////////////////////////////////
	;-------------------------------------------
		oNe = MrVar_Get(ne_vnames[i])
		oNe['COLOR'] = colors[i]
		oNe['LABEL'] = sc[i]
	
	;-------------------------------------------
	; Moments //////////////////////////////////
	;-------------------------------------------
		oDF       = MrDist4D(f_vnames[i], VSC=scpot_vnames[i], SPECIES=type)
		oN        = oDF -> Density()
		oP_tensor = oDF -> Pressure()
		oP        = (oP_tensor[*,0,0] + oP_tensor[*,1,1] + oP_tensor[*,2,2]) / 3.0
	
	;-------------------------------------------
	; Kinetic Entropy //////////////////////////
	;-------------------------------------------
		;   - Convert 1/cm^3 to 1/m^3
		;   - Convert J/K to eV/K
		oSn  = (J2eV * oDF -> Entropy()) / (1e6 * oN)
		
		;Sn
		oSn -> SetName, sn_vnames[i]
		oSn -> Cache
		oSn['CATDESC']       = 'Entropy per particle computed via integration of the full velocity ' + $
		                       'space distribution function.'
		oSn['COLOR']         = colors[i]
		oSn['PLOT_TITLE']    = 'Entropy Per Particle'
;		oSn['LABEL']         = sc[i]
		oSn['TITLE']         = 's/N!C(eV/K)'
		oSn['UNITS']         = 'eV/K ln(s^3/m^6)'
		oSn['SI_CONVERSION'] = '1.602e-19>J/K'
	
	;-------------------------------------------
	; Maxwellian Entropy ///////////////////////
	;-------------------------------------------
;		oN        = MrVar_Get(n_vnames[i])
;		oP_tensor = MrVar_Get(p_vnames[i])
		oSbn = 3.0/2.0 * kB * J2eV * 1e6 * oN $
		       * ALog( 2 * !pi * 1e-11 / mass * oP['DATA'] / oN['DATA']^(5.0/3.0) + 1 )
		oSbn = oSbn / (1e6 * oN)
		Obj_Destroy, [oN, oP_tensor, oP]
	
		;Set Attributes
		oSbn -> SetName, sbn_vnames[i]
		oSbn -> Cache
		oSbn['CATDESC']       = 'Maxwellian entropy per particle.'
		oSbn['COLOR']         = colors[i]
		oSbn['PLOT_TITLE']    = 'Maxwellian Entropy per Particle'
;		oSbn['LABEL']         = 's$\downB$'
		oSbn['TITLE']         = 's$\downB$/N!C(eV/K)'
		oSbn['UNITS']         = 'eV/K ln(J m^2/kg)'
		oSbn['SI_CONVERSION'] = '1.602e-19>J/K ln(J m^2/kg)'
	
	;-------------------------------------------
	; M-Bar ////////////////////////////////////
	;-------------------------------------------
		oMnbar = oSn - oSbn
		oMnbar -> SetName, mnbar_vnames[i]
		oMnbar -> Cache
		
		;Set Attributes
		oMnbar['CATDESC']       = 'Difference between the entropy per particle of the ' + $
		                          'full distribution and that of an equivalent Maxwellian.'
		oMnbar['COLOR']         = colors[i]
		oMnbar['PLOT_TITLE']    = 'Maxwellian Entropy Density'
;		oMnbar['LABEL']         = sc[i]
		oMnbar['TITLE']         = 'M/N!C(eV/K)'
		oMnbar['UNITS']         = 'eV/K [ln(s^3/m^6) - ln(J m^2/kg)]'
		oMnbar['SI_CONVERSION'] = '1.602e-19>J/K [ln(s^3/m^6) - ln(J m^2/kg)]'
	
	;-------------------------------------------
	; Agyrotropy ///////////////////////////////
	;-------------------------------------------
		oQ = MrVar_Pres_QFactor(bvec_vnames[i], p_vnames[i], /CACHE, NAME=q_vnames[i])
		oQ['COLOR'] = colors[i]
	
	;-------------------------------------------
	; Current Density //////////////////////////
	;-------------------------------------------
		;Pull data
		oVi = MrVar_Get(vi_vnames[i])
		oVe = MrVar_Get(ve_vnames[i])
		oNe = MrVar_Get(ne_vnames[i])
	
		;Interpolate to DES times
		oVi_des = oVi -> Interpol(oVe)
	
		;Calculate current density
		;   - 1e15 converts to uA/m^2
		oJ = 1e15 * q * oNe * (oVi_des - oVe)
		Obj_Destroy, oVi_des
	
		;Attributes
		oJ -> SetName, j_vnames[i]
		oJ -> Cache
		oJ['CATDESC']       = 'Current density derived from moments of the distribution function.'
		oJ['COLOR']         = ['Blue', 'Forest Green', 'Red']
		oJ['LABEL']         = 'J$\down' + ['X', 'Y', 'Z'] + '$'
		oJ['TITLE']         = 'J!C($\mu$A/m^2)'
		oJ['UNITS']         = 'uA/m^2'
		oJ['SI_CONVERSION'] = '1e-6>A/m^2'
	
	;-------------------------------------------
	; E-Prime //////////////////////////////////
	;-------------------------------------------
		oE      = MrVar_Get(e_vnames[i])
		oB      = MrVar_Get(bvec_vnames[i])
		oB_des  = oB -> Interpol(oVe)
		oE_des  = oE -> Interpol(oVe)
		oVexB   = MrVar_E_VxB(oVe, oB_des)
		oEprime = oE_des - oVexB
		Obj_Destroy, [oE_des, oB_des, oVexB]
	
		oEprime -> SetName, eprime_vnames[i]
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
		oJdotE -> SetName, jdote_vnames[i]
		oJdotE -> Cache
		oJdotE['CATDESC']       = 'Electron frame dissipation measure.'
		oJdotE['COLOR']         = colors[i]
		oJdotE['TITLE']         = "J.E'!C(nW/m^3)"
		oJdotE['UNITS']         = 'nW/m^3'
		oJdotE['SI_CONVERSION'] = '1e9>W/m^3'
	ENDFOR

;-------------------------------------------
; Gradient and Divergence //////////////////
;-------------------------------------------
	oR1 = MrVar_Get(r_vnames[0])
	oR2 = MrVar_Get(r_vnames[1])
	oR3 = MrVar_Get(r_vnames[2])
	oR4 = MrVar_Get(r_vnames[3])
	oMn1 = MrVar_Get(mnbar_vnames[0])
	oMn2 = MrVar_Get(mnbar_vnames[1])
	oMn3 = MrVar_Get(mnbar_vnames[2])
	oMn4 = MrVar_Get(mnbar_vnames[3])
	
	;Interpolate to MMS1
	oR1   = oR1  -> Interpol(oMn1)
	oR2   = oR2  -> Interpol(oMn1)
	oR3   = oR3  -> Interpol(oMn1)
	oR4   = oR4  -> Interpol(oMn1)
	oMn2   = oMn2  -> Interpol(oMn1)
	oMn3   = oMn3  -> Interpol(oMn1)
	oMn4   = oMn4  -> Interpol(oMn1)
	
	;Compute the gradient
	;   - Convert 1/km to 1/m
	oRecip  = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	oGradMn = 1e-3 * oRecip -> Gradient(oMn1, oMn2, oMn3, oMn4)
	
	;Attributes
	oGradMn -> SetName, gradmn_vname
	oGradMn -> Cache
	oGradMn['CATDESC']       = 'Gradient of the non-Maxwellian entropy per particle.'
	oGradMn['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oGradMn['LABEL']         = nabla + '(M/N)$\down' + ['X', 'Y', 'Z'] + '$'
;	oGradMn['LINESTYLE']     = '--'
	oGradMn['TITLE']         = nabla + '(M/N)!C(eV/K/m)'
	oGradMn['UNITS']         = 'eV/K/m [ln(s^3/m^6) - ln(J m^2/kg)]'
	oGradMn['SI_CONVERSION'] = '1.602e-19>J/K/m [ln(s^3/m^6) - ln(J m^2/kg)]'
	
	
	;CURLOMETER
	oR1 = MrVar_Get(r_vnames[0])
	oR2 = MrVar_Get(r_vnames[1])
	oR3 = MrVar_Get(r_vnames[2])
	oR4 = MrVar_Get(r_vnames[3])
	oB1 = MrVar_Get(bvec_vnames[0])
	oB2 = MrVar_Get(bvec_vnames[1])
	oB3 = MrVar_Get(bvec_vnames[2])
	oB4 = MrVar_Get(bvec_vnames[3])

	;Interpolate to MMS1
	oR1   = oR1  -> Interpol(oB1)
	oR2   = oR2  -> Interpol(oB1)
	oR3   = oR3  -> Interpol(oB1)
	oR4   = oR4  -> Interpol(oB1)
	oB2   = oB2  -> Interpol(oB1)
	oB3   = oB3  -> Interpol(oB1)
	oB4   = oB4  -> Interpol(oB1)
	
	;Compute the curl
	;   - Converts to uA/m^2
	oRecip = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	oCurlB = 1e-6 * oRecip -> Curl(oB1, oB2, oB3, oB4)
	
	;Attributes
	oCurlB -> SetName, jcurl_vname
	oCurlB -> Cache
	oCurlB['CATDESC']       = 'Current density calculated via the curlometer technique.'
	oCurlB['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oCurlB['LABEL']         = 'J$\down' + ['X', 'Y', 'Z'] + '$'
	oCurlB['TITLE']         = 'J$\down' + nabla + 'xB$!C(uA/m$\up2$)'
	oCurlB['UNITS']         = 'uA/m^2'
	oCurlB['SI_CONVERSION'] = '1.602e-6>A/m^2'

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------

	;B
	oB = MrVar_Get(b_bc_vname)
	oB['PLOT_TITLE'] = 'MMS1-4 ' + StrUpCase(fpi_instr)
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [b_bc_vname, ne_vnames[0], sn_vnames[0], sbn_vnames[0], mnbar_vnames[0], gradmn_vname[0], q_vnames[0], jdote_vnames[0]], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )
	win = MrVar_OPlotTS( ne_vnames[0], ne_vnames[1:*] )
	win = MrVar_OPlotTS( sn_vnames[0], sn_vnames[1:*] )
	win = MrVar_OPlotTS( sbn_vnames[0], sbn_vnames[1:*] )
	win = MrVar_OPlotTS( mnbar_vnames[0], mnbar_vnames[1:*] )
	win = MrVar_OPlotTS( q_vnames[0], q_vnames[1:*] )
	win = MrVar_OPlotTS( jdote_vnames[0], jdote_vnames[1:*] )

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
		fname = StrJoin(['mms', fpi_instr, fpi_mode, level, '4sc-kentropy'], '_')
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