; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_PoyntingTheorem
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
;       2. S: Entropy Density
;       3. Grad(S)
;       4. Sf: Entropy Density Flux
;       5. Div(Sf)
;       6. Mbar: S - Sb
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
;       2018/05/22  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_PoyntingTheorem, $
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
	
	mode = 'brst'
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(ephdesc)   EQ 0 THEN ephdesc   = 'epht89d'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	
	q     = MrConstants('q')
	mu0   = MrConstants('mu_0')
	e0    = MrConstants('epsilon_0')
	m_e   = MrConstants('m_e')
	m_i   = MrConstants('m_H')
	nabla = '!9' + String(71B) + '!X'
	J2eV  = MrConstants('J2eV')
	colors = ['Black', 'Red', 'Forest Green', 'Blue']

;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	fgm_mode    = mode
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	fgm_coords  = (coords EQ 'dbcs' || coords EQ 'dsl') ? 'dmpa' : coords
	
;	fpi_instr  = 'd' + species + 's'
	fpi_mode   = mode EQ 'brst' ? mode : 'fast'
	fpi_coords = (coords EQ 'dmpa' || coords EQ 'dsl') ? 'dbcs' : coords
	
	edp_mode   = 'fast'
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
	ni_vnames    = sc + '_' + StrJoin( ['dis', 'numberdensity',             fpi_mode], '_' )
	vi_vnames    = sc + '_' + StrJoin( ['dis', 'bulkv',         fpi_coords, fpi_mode], '_' )
	ne_vnames    = sc + '_' + StrJoin( ['des', 'numberdensity',             fpi_mode], '_' )
	ve_vnames    = sc + '_' + StrJoin( ['des', 'bulkv',         fpi_coords, fpi_mode], '_' )
	r_vnames     = sc + '_' + StrJoin( ['mec', 'r', coords], '_' )
	
	;Derived names
	b_des_vnames  = sc + '_' + StrJoin( [fgm_instr, 'bvec', coords, 'des', mode, level], '_' )
	e_des_vnames  = sc + '_' + StrJoin( ['edp', 'dce', edp_coords, 'des', edp_mode, level], '_')
	ni_des_vnames = sc + '_' + StrJoin( ['dis', 'numberdensity', 'des', fpi_mode], '_' )
	vi_des_vnames = sc + '_' + StrJoin( ['dis', 'bulkv', fpi_coords, 'des', fpi_mode], '_' )
	ne_des_vnames = sc + '_' + StrJoin( ['des', 'numberdensity', 'des', fpi_mode], '_' )
	ve_des_vnames = sc + '_' + StrJoin( ['des', 'bulkv', fpi_coords, 'des', fpi_mode], '_' )
	r_des_vnames  = sc + '_' + StrJoin( ['mec', 'r', coords, 'des'], '_' )
	
	b_bc_vname     = StrJoin( ['mms', fgm_instr, 'bvec', coords, 'barycenter', mode, level], '_' )
	e_bc_vname     = StrJoin( ['mms', 'edp', 'dce', edp_coords, 'barycenter', edp_mode, level], '_')
	jdote_bc_vname = StrJoin( ['mms', 'edp', 'jdote', edp_coords, 'barycenter', edp_mode, level], '_')
	dudt_bc_vname  = StrJoin( ['mms', 'fields', 'dudt', edp_coords, 'barycenter', edp_mode, level], '_')

	j_vnames      = sc + '_' + StrJoin(['fpi', 'j',                              fpi_mode], '_')
	jcurl_vname   = StrJoin( ['mms', fgm_instr, 'j', coords, 'curlometer', mode, level], '_' )
	eprime_vnames = sc + '_' + StrJoin(['edp',    'eprime', edp_coords, edp_mode, level], '_')
	jdote_vnames  = sc + '_' + StrJoin(['edp',    'jdote',  edp_coords, edp_mode, level], '_')
	s_vnames      = sc + '_' + StrJoin(['fields', 's',      edp_coords, edp_mode, level], '_')
	u_vnames      = sc + '_' + StrJoin(['fields', 'u',      edp_coords, edp_mode, level], '_')
	dudt_vnames   = sc + '_' + StrJoin(['fields', 'dudt',   edp_coords, edp_mode, level], '_')
	divs_vname    = StrJoin( ['mms', 'fields', 'divs', edp_coords, edp_mode, level], '_')
	
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
; Interpolate to DES on MMS1 ///////////////
;-------------------------------------------
	oNe = MrVar_Get(ne_vnames[0])
	oT  = oNe['TIMEVAR']
	FOR i = 0, N_Elements(sc) - 1 DO BEGIN
		oPos = MrVar_Get(r_vnames[i])
		oB   = MrVar_Get(bvec_vnames[i])
		oE   = MrVar_Get(e_vnames[i])
		oNe  = MrVar_Get(ne_vnames[i])
		oVe  = MrVar_Get(ve_vnames[i])
		oNi  = MrVar_Get(ni_vnames[i])
		oVi  = MrVar_Get(vi_vnames[i])
		
		oPos = oPos -> Interpol(oT, /CACHE, NAME=r_des_vnames[i])
		oB   = oB   -> Interpol(oT, /CACHE, NAME=b_des_vnames[i])
		oE   = oE   -> Interpol(oT, /CACHE, NAME=e_des_vnames[i])
		oNe  = oNe  -> Interpol(oT, /CACHE, NAME=ne_des_vnames[i])
		oVe  = oVe  -> Interpol(oT, /CACHE, NAME=ve_des_vnames[i])
		oNi  = oNi  -> Interpol(oT, /CACHE, NAME=ni_des_vnames[i])
		oVi  = oVi  -> Interpol(oT, /CACHE, NAME=vi_des_vnames[i])
	ENDFOR

;-------------------------------------------
; Loop Over Each Spacecraft ////////////////
;-------------------------------------------
	FOR i = 0, N_Elements(sc) - 1 DO BEGIN
	;-------------------------------------------
	; Current Density //////////////////////////
	;-------------------------------------------
		;Pull data
		oVi = MrVar_Get(vi_des_vnames[i])
		oVe = MrVar_Get(ve_des_vnames[i])
		oNe = MrVar_Get(ne_des_vnames[i])
	
		;Calculate current density
		;   - 1e15 converts to uA/m^2
		oJ = 1e15 * q * oNe * (oVi - oVe)
	
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
		oE_des  = MrVar_Get(e_des_vnames[i])
		oB_des  = MrVar_Get(b_des_vnames[i])
		oVexB   = MrVar_E_VxB(oVe, oB_des)
		oEprime = oE_des - oVexB
		Obj_Destroy, oVexB
	
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

	;-------------------------------------------
	; Poynting Flux ////////////////////////////
	;-------------------------------------------
		oS = 1e-12/mu0 * J2eV * oE_des -> Cross(oB_des)
		oS -> SetName, s_vnames[i]
		oS -> Cache
		oS['CATDESC']       = 'Poynting vector.'
		oS['COLOR']         = ['Blue', 'Forest Green', 'Red']
		oS['TITLE']         = 'S!CeV/s/m$\up2$'
		oS['UNITS']         = 'eV/s/m^2'
		oS['SI_CONVERSION'] = '1.602e-19>J/s/m^2'
		
	;-------------------------------------------
	; Electromagnetic Energy Density ///////////
	;-------------------------------------------
		oE = oE_des -> Magnitude()
		oB = oB_des -> Magnitude()
		ou = (1e-6 * J2eV * e0 * oE^2) + (1e-18 * J2eV / mu0 * oB^2)
		ou -> SetName, u_vnames[i]
		ou -> Cache
		ou['CATDESC']       = 'Electromagnetic energy density.'
		ou['COLOR']         = 'Black'
		ou['TITLE']         = 'u!CeV/m$\up3$'
		ou['UNITS']         = 'eV/m^3'
		ou['SI_CONVERSION'] = '1.602e-19>J/m^3'
		
	;-------------------------------------------
	; du/dt ////////////////////////////////////
	;-------------------------------------------
		dt           = ou['TIMEVAR'] -> GetSI()
		odudt        = ou -> Copy()
		odudt[1:*] = 1e9/J2eV * (ou[1:*] - ou[0:-2]) / dt
		odudt[[0]] = !values.f_nan
		odudt -> SetName, dudt_vnames[i]
		odudt -> Cache
		odudt['CATDESC']       = 'Time derivative of the electromagnetic energy density.'
		odudt['COLOR']         = 'Black'
		odudt['TITLE']         = 'u!CnW/m$\up3$'
		odudt['UNITS']         = 'nW/m^3'
		odudt['SI_CONVERSION'] = '1.602e-19>J/s/m^3'
	ENDFOR

;-------------------------------------------
; Gradient and Divergence //////////////////
;-------------------------------------------
	oR1 = MrVar_Get(r_des_vnames[0])
	oR2 = MrVar_Get(r_des_vnames[1])
	oR3 = MrVar_Get(r_des_vnames[2])
	oR4 = MrVar_Get(r_des_vnames[3])
	oB1 = MrVar_Get(b_des_vnames[0])
	oB2 = MrVar_Get(b_des_vnames[1])
	oB3 = MrVar_Get(b_des_vnames[2])
	oB4 = MrVar_Get(b_des_vnames[3])
	oS1 = MrVar_Get(s_vnames[0])
	oS2 = MrVar_Get(s_vnames[1])
	oS3 = MrVar_Get(s_vnames[2])
	oS4 = MrVar_Get(s_vnames[3])
	
	;Compute the spatial derivatives
	;   - Converts to uA/m^2
	;   - Converts to nW/m^3
	oRecip = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	oCurlB = 1e-6 * oRecip -> Curl(oB1, oB2, oB3, oB4)
	oDivS  = 1e6/J2eV * oRecip -> Divergence(oS1, oS2, oS3, oS4)
	
	;Attributes
	oCurlB -> SetName, jcurl_vname
	oCurlB -> Cache
	oCurlB['CATDESC']       = 'Current density calculated via the curlometer technique.'
	oCurlB['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oCurlB['LABEL']         = 'J$\down' + ['X', 'Y', 'Z'] + '$'
	oCurlB['TITLE']         = 'J$\down' + nabla + 'xB$!C(uA/m$\up2$)'
	oCurlB['UNITS']         = 'uA/m^2'
	oCurlB['SI_CONVERSION'] = '1.602e-6>A/m^2'
	
	oDivS -> SetName, divs_vname
	oDivS -> Cache
	oDivS['CATDESC']       = 'Divergence of the Poynting flux.'
	oDivS['TITLE']         = nabla + '.S!CnW/m^3'
	oDivS['UNITS']         = 'nW/m^3'
	oDivS['SI_CONVERSION'] = '1.602e-19>nW/m^3'
	
;-------------------------------------------
; Barycenter ///////////////////////////////
;-------------------------------------------
	;B-Field
	oB1 = MrVar_Get(b_des_vnames[0])
	oB2 = MrVar_Get(b_des_vnames[1])
	oB3 = MrVar_Get(b_des_vnames[2])
	oB4 = MrVar_Get(b_des_vnames[3])
	
	oBbc = (oB1 + oB2 + oB3 + oB4) / 4.0
	oBbc -> SetName, b_bc_vname
	oBbc -> Cache
	oBbc['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oBbc['LABEL']         = 'B$\down' + ['X', 'Y', 'Z'] + '$'
	oBbc['SI_CONVERSION'] = '1e-9>T'
	oBbc['TITLE']         = 'B$\downBC$!C(nT)'
	oBbc['UNITS']         = 'nT'
	
	;E-Field
	oE1 = MrVar_Get(e_des_vnames[0])
	oE2 = MrVar_Get(e_des_vnames[1])
	oE3 = MrVar_Get(e_des_vnames[2])
	oE4 = MrVar_Get(e_des_vnames[3])
	
	oEbc = (oE1 + oE2 + oE3 + oE4) / 4.0
	oEbc -> SetName, e_bc_vname
	oEbc -> Cache
	oEbc['TITLE']         = 'E$\downBC$!C(mV/m)'
	oEbc['UNITS']         = 'mV/m'
	oEbc['SI_CONVERSION'] = '1e-3>V/m'
	
	;J.E
	oJdotE1 = MrVar_Get(jdote_vnames[0])
	oJdotE2 = MrVar_Get(jdote_vnames[1])
	oJdotE3 = MrVar_Get(jdote_vnames[2])
	oJdotE4 = MrVar_Get(jdote_vnames[3])
	
	oJdotEbc = (oJdotE1 + oJdotE2 + oJdotE3 + oJdotE4) / 4.0
	oJdotEbc -> SetName, jdote_bc_vname
	oJdotEbc -> Cache
	oJdotEbc['CATDESC']       = 'Electron frame dissipation measure at the barycenter.'
	oJdotEbc['COLOR']         = 'Black';['Blue', 'Forest Green', 'Red']
	oJdotEbc['TITLE']         = "J.E!C(nW/m^3)"
	oJdotEbc['UNITS']         = 'nW/m^3'
	oJdotEbc['SI_CONVERSION'] = '1e9>W/m^3'
	
	;du/dt
	odudt1 = MrVar_Get(dudt_vnames[0])
	odudt2 = MrVar_Get(dudt_vnames[1])
	odudt3 = MrVar_Get(dudt_vnames[2])
	odudt4 = MrVar_Get(dudt_vnames[3])
	
	odudtbc = (odudt1 + odudt2 + odudt3 + odudt4) / 4.0
	odudtbc -> SetName, dudt_bc_vname
	odudtbc -> Cache
	odudtbc['CATDESC']       = 'Time derivative of the electromagnetic energy density at the barycenter.'
	odudtbc['COLOR']         = 'Black';['Blue', 'Forest Green', 'Red']
	odudtbc['TITLE']         = 'du/dt!CnW/m$\up3$'
	odudtbc['UNITS']         = 'nW/m^3'
	odudtbc['SI_CONVERSION'] = '1.602e-19>J/s/m^3'

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------

	;B
	oB = MrVar_Get(b_bc_vname)
	oB['PLOT_TITLE'] = 'MMS 1-4'
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [b_bc_vname, jdote_bc_vname, divs_vname, dudt_bc_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )

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
		fname = StrJoin(['mms', fpi_instr, fpi_mode, level, '4sc-poynting'], '_')
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