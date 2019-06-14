; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_KEntropy
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
;       2. sn: Entropy per particle             (i,e)
;       3. sbn: Maxwellian entropy per particle (i,e)
;       4. Grad(sbn)                            (i,e)
;       5. sfn: Entropy flux per particle       (i,e)
;       6. Div(sfn)                             (i,e)
;       7. Q: Agyrotropy parameter
;       8: J.E'
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
;       2018/04/03  -   Written by Matthew Argall
;       2018/06/07  -   Use calculated moments for sB so that the comparison with
;                           s is fair. Fix units error in calculation of sB - MRA
;-
FUNCTION MrMMS_Plot_FPI_KEntropy, sc, mode, $
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
	nabla = '!9' + String(71B) + '!X'
	J2eV  = MrConstants('J2eV')
	kB    = MrConstants('k_B')
	m_e   = MrConstants('m_e')
	m_i   = MrConstants('m_H')

;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	fgm_coords  = (coords EQ 'dbcs' || coords EQ 'dsl') ? 'dmpa' : coords
	
	fpi_mode   = mode EQ 'brst' ? mode : 'fast'
	fpi_coords = (coords EQ 'dmpa' || coords EQ 'dsl') ? 'dbcs' : coords
	
	edp_mode   = mode EQ 'brst' ? mode : 'fast'
	edp_coords = (coords EQ 'dbcs' || coords EQ 'dmpa') ? 'dsl' : coords
	
	eph_mode = mode

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	scraft = 'mms' + ['1', '2', '3', '4']
	isc    = where(scraft eq sc, n)
	IF n EQ 0 THEN Message, 'Invalid spacecraft ID "' + sc + '".'

	;Source Names
	b_vname     = StrJoin( [sc, fgm_instr, 'b',    coords, mode, level], '_' )
	bvec_vname  = StrJoin( [sc, fgm_instr, 'bvec', coords, mode, level], '_' )
	bmag_vname  = StrJoin( [sc, fgm_instr, 'bmag', coords, mode, level], '_' )
	e_vname     = StrJoin( [sc, 'edp', 'dce', edp_coords, edp_mode, level], '_')
	scpot_vname = StrJoin( [sc, 'edp', 'scpot',           edp_mode, level], '_')
	fi_vnames   = scraft + '_' + StrJoin( ['dis', 'dist',           fpi_mode], '_' )
	ni_vname    = StrJoin( [sc, 'dis', 'numberdensity',             fpi_mode], '_' )
	vi_vname    = StrJoin( [sc, 'dis', 'bulkv',         fpi_coords, fpi_mode], '_' )
	pi_vname    = StrJoin( [sc, 'dis', 'prestensor',    fpi_coords, fpi_mode], '_' )
	fe_vnames   = scraft + '_' + StrJoin( ['des', 'dist',           fpi_mode], '_' )
	ne_vname    = StrJoin( [sc, 'des', 'numberdensity',             fpi_mode], '_' )
	ve_vname    = StrJoin( [sc, 'des', 'bulkv',         fpi_coords, fpi_mode], '_' )
	pe_vname    = StrJoin( [sc, 'des', 'prestensor',    fpi_coords, fpi_mode], '_' )
	r_vnames    = scraft + '_' + StrJoin( ['mec', 'r', coords], '_' )
	
	;Derived names
	ni_calc_vname = StrJoin( [sc, 'dis', 'numberdensity', 'calc',             fpi_mode], '_' )
	pi_calc_vname = StrJoin( [sc, 'dis', 'prestensor',    'calc', fpi_coords, fpi_mode], '_' )
	ne_calc_vname = StrJoin( [sc, 'des', 'numberdensity', 'calc',             fpi_mode], '_' )
	pe_calc_vname = StrJoin( [sc, 'des', 'prestensor',    'calc', fpi_coords, fpi_mode], '_' )
	si_vnames     = scraft + '_' + StrJoin(['dis', 'kinetic-entropy-density',           fpi_mode], '_')
	sni_vnames    = scraft + '_' + StrJoin(['dis', 'kinetic-entropy-per-particle',      fpi_mode], '_')
	sfni_vnames   = scraft + '_' + StrJoin(['dis', 'kinetic-entropy-flux-per-particle', fpi_mode], '_')
	sbi_vname     = StrJoin([sc,    'dis', 'maxwellian-entropy-density',            fpi_mode], '_')
	sbni_vname    = StrJoin([sc,    'dis', 'maxwellian-entropy-per-particle',       fpi_mode], '_')
	gradsni_vname = StrJoin(['mms', 'dis', 'grad-kinetic-entropy-per-particle',     fpi_mode], '_')
	divsfni_vname = StrJoin(['mms', 'dis', 'div-kinetic-entropy-flux-per-particle', fpi_mode], '_')
	qi_vname      = StrJoin([sc,    'dis', 'agyrotropy-factor',                     fpi_mode], '_')
	se_vnames     = scraft + '_' + StrJoin(['des', 'kinetic-entropy-density',           fpi_mode], '_')
	sne_vnames    = scraft + '_' + StrJoin(['des', 'kinetic-entropy-per-particle',      fpi_mode], '_')
	sfne_vnames   = scraft + '_' + StrJoin(['des', 'kinetic-entropy-flux-per-particle', fpi_mode], '_')
	sbe_vname     = StrJoin([sc,    'des', 'maxwellian-entropy',                    fpi_mode], '_')
	sbne_vname    = StrJoin([sc,    'des', 'maxwellian-entropy-per-particle',       fpi_mode], '_')
	gradsne_vname = StrJoin(['mms', 'des', 'grad-kinetic-entropy-per-particle',     fpi_mode], '_')
	divsfne_vname = StrJoin(['mms', 'des', 'div-kinetic-entropy-flux-per-particle', fpi_mode], '_')
	qe_vname      = StrJoin([sc,    'des', 'sqrtq',                        fpi_mode], '_')
	j_vname       = StrJoin([sc,    'fpi', 'j',                            fpi_mode], '_')
	eprime_vname  = StrJoin([sc,    'edp', 'eprime', edp_coords, edp_mode, level], '_')
	jdote_vname   = StrJoin([sc,    'edp', 'jdote',  edp_coords, edp_mode, level], '_')
	
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
		
		;DIS-Dist
		MrMMS_FPI_Load_Dist3D, 'mms1', mode, 'i', /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms2', mode, 'i', /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms3', mode, 'i', /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms4', mode, 'i', /APPLY_MODEL
		
		;DIS-Moms
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = 'dis-moms', $
		                     VARFORMAT = [ '*numberdensity*', '*bulkv_'+coords+'*', $
		                                   '*pres*'+coords+'*' ]
		
		;DES-Dist
		MrMMS_FPI_Load_Dist3D, 'mms1', mode, 'e', /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms2', mode, 'e', /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms3', mode, 'e', /APPLY_MODEL
		MrMMS_FPI_Load_Dist3D, 'mms4', mode, 'e', /APPLY_MODEL
		
		;DES-Moms
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = 'des-moms', $
		                     VARFORMAT = [ '*numberdensity*', '*bulkv_'+coords+'*', $
		                                   '*pres*'+coords+'*' ]
		
		;Spacecraft potential
		MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
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
; Loop Over Each Spacecraft ////////////////
;-------------------------------------------
	
	FOR i = 0, N_Elements(scraft) - 1 DO BEGIN
	
	;-------------------------------------------
	; Entropy Per Particle /////////////////////
	;-------------------------------------------
		;DIS
		;   - Convert 1/cm^3 to 1/m^3
		;   - Convert J/K to eV/K
		oDFi  = MrDist4D(fi_vnames[i], VSC=scpot_vname, SPECIES='H')
		
		;Density and pressure
		oNi = oDFi -> Density(CACHE=(i EQ iSC), NAME=ni_calc_vname)
		oPi = oDFi -> Pressure(CACHE=(i EQ iSC), NAME=pi_calc_vname)
		
		;Entropy per particle
		oSi  = J2eV * oDFi -> Entropy(/CACHE, NAME=si_vnames[i])
		oSni = oSi / (1e6 * oNi)
;		Obj_Destroy, oDFi
	
		;DES
		;   - Convert 1/cm^3 to 1/m^3
		;   - Convert J/K to eV/K
		oDFe  = MrDist4D(fe_vnames[i], VSC=scpot_vname, SPECIES='e')
		
		;Density and pressure
		oNe = oDFe -> Density(CACHE=(i EQ iSC), NAME=ne_calc_vname)
		oPe = oDFe -> Pressure(CACHE=(i EQ iSC), NAME=pe_calc_vname)
		
		;Entropy per particle
		oSe  = J2eV * oDFe -> Entropy(/CACHE, NAME=se_vnames[i])
		oSne = oSe / (1e6 * oNe)
;		Obj_Destroy, oDFe
	
		;Set attributes
		oSni -> SetName, sni_vnames[i]
		oSni -> Cache
		oSni['CATDESC']       = 'Entropy per particle computed via integration of the full ' + $
		                        'velocity space distribution function.'
		oSni['COLOR']         = 'Blue'
		oSni['PLOT_TITLE']    = 'Entropy per Particle'
		oSni['LABEL']         = 'dis'
		oSni['TITLE']         = 's/N!C(eV/K)'
		oSni['UNITS']         = 'eV/K'
		oSni['SI_CONVERSION'] = '1.602e-19>J/K'
	
		oSne -> SetName, sne_vnames[i]
		oSne -> Cache
		oSne['CATDESC']       = 'Entropy per particle computed via integration of the full ' + $
		                        'velocity space distribution function.'
		oSne['COLOR']         = 'Red'
		oSne['PLOT_TITLE']    = 'Entropy per Particle'
		oSne['LABEL']         = 'des'
		oSne['TITLE']         = 's/N!C(eV/K)'
		oSne['UNITS']         = 'eV/K'
		oSne['SI_CONVERSION'] = '1.602e-19>J/K'
	
	;-------------------------------------------
	; Entropy Flux Per Particle ////////////////
	;-------------------------------------------
		;DIS
		;   - Convert J/K to eV/K
		;   - Convert 1/cm^3 to 1/m^3
		oSfi  = J2eV * oDFi -> EntropyFlux()
		oSfni = oSfi / (1e6 * oNi)
		Obj_Destroy, [oDFi, oSfi]
	
		;DES
		;   - Convert J m / K s to eV m / K s
		;   - Convert 1/cm^3 to 1/m^3
		oSfe  = J2eV * oDFe -> EntropyFlux()
		oSfne = oSfe / (1e6 * oNe)
		Obj_Destroy, [oDFe, oSfe]
	
		;Set attributes
		oSfni -> SetName, sfni_vnames[i]
		oSfni -> Cache
		oSfni['CATDESC']       = 'Entropy flux per particle computed via integration of the full ' + $
		                         'velocity space distribution function.'
		oSfni['COLOR']         = ['Blue', 'Forest Green', 'Red']
		oSfni['PLOT_TITLE']    = 'Entropy Flux per Particle'
		oSfni['LABEL']         = 'Sf$\down' + ['X', 'Y', 'Z'] + '$'
		oSfni['LINESTYLE']     = '--'
		oSfni['TITLE']         = 'sf/N!C(eV m K$\up-1$ s$\up-1$)'
		oSfni['UNITS']         = 'eV m K^-1 s^-1'
		oSfni['SI_CONVERSION'] = '1.602e-19>J m K^-1 s^-1'
	
		oSfne -> SetName, sfne_vnames[i]
		oSfne -> Cache
		oSfne['CATDESC']       = 'Entropy flux per particle computed via integration of the full ' + $
		                         'velocity space distribution function.'
		oSfne['COLOR']         = ['Blue', 'Forest Green', 'Red']
		oSfne['PLOT_TITLE']    = 'Entropy Flux per Particle'
;		oSfne['LABEL']         = 'des'
		oSfne['TITLE']         = 'sf/N!C(eV m K$\up-1$ s$\up-1$)'
		oSfne['UNITS']         = 'eV m K^-1 s^-1'
		oSfne['SI_CONVERSION'] = '1.602e-19>J m K^-1 s^-1'
	ENDFOR
	
;-------------------------------------------
; Gradient and Divergence //////////////////
;-------------------------------------------
	;DIS
	oR1    = MrVar_Get(r_vnames[0])
	oR2    = MrVar_Get(r_vnames[1])
	oR3    = MrVar_Get(r_vnames[2])
	oR4    = MrVar_Get(r_vnames[3])
	oSni1  = MrVar_Get(sni_vnames[0])
	oSni2  = MrVar_Get(sni_vnames[1])
	oSni3  = MrVar_Get(sni_vnames[2])
	oSni4  = MrVar_Get(sni_vnames[3])
	oSfni1 = MrVar_Get(sfni_vnames[0])
	oSfni2 = MrVar_Get(sfni_vnames[1])
	oSfni3 = MrVar_Get(sfni_vnames[2])
	oSfni4 = MrVar_Get(sfni_vnames[3])
	
	;Interpolate to MMS1
	oR1    = oR1    -> Interpol(oSni1)
	oR2    = oR2    -> Interpol(oSni1)
	oR3    = oR3    -> Interpol(oSni1)
	oR4    = oR4    -> Interpol(oSni1)
	oSni2  = oSni2  -> Interpol(oSni1)
	oSni3  = oSni3  -> Interpol(oSni1)
	oSni4  = oSni4  -> Interpol(oSni1)
	oSfni2 = oSfni2 -> Interpol(oSni1)
	oSfni3 = oSfni3 -> Interpol(oSni1)
	oSfni4 = oSfni4 -> Interpol(oSni1)
	
	;Compute the gradient
	;   - Convert 1/km to 1/m
	oRecip   = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	oGradSni = 1e-3 * oRecip -> Gradient(oSni1, oSni2, oSni3, oSni4)
	oDivSfni = 1e-3 * oRecip -> Divergence(oSfni1, oSfni2, oSfni3, oSfni4)
	
	;Attributes
	oGradSni -> SetName, gradsni_vname
	oGradSni -> Cache
	oGradSni['CATDESC']       = 'Gradient of the kinetic energy per particle'
	oGradSni['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oGradSni['LABEL']         = nabla + 'Sn$\down' + ['X', 'Y', 'Z'] + '$'
	oGradSni['LINESTYLE']     = '--'
	oGradSni['TITLE']         = nabla + 'Sni!C(eV K$\up-1$ s$\up-1$)'
	oGradSni['UNITS']         = 'eV K^-1 m^-1'
	oGradSni['SI_CONVERSION'] = '1.602e-19>J K^-1 m^-1'
	
	oDivSfni -> SetName, divsfni_vname
	oDivSfni -> Cache
	oDivSfni['CATDESC']       = 'Divergence of the kinetic energy flux per particle'
	oDivSfni['COLOR']         = 'Blue'
	oDivSfni['LABEL']         = nabla + '.Sfni'
	oDivSfni['TITLE']         = nabla + '.Sfni!C(eV K$\up-1$ m$\up-1$)'
	oDivSfni['UNITS']         = 'eV K^-1 m^-1'
	oDivSfni['SI_CONVERSION'] = '1.602e-19>J K-1 m^-1'
	
	
	;DES
	oR1    = MrVar_Get(r_vnames[0])
	oR2    = MrVar_Get(r_vnames[1])
	oR3    = MrVar_Get(r_vnames[2])
	oR4    = MrVar_Get(r_vnames[3])
	oSne1  = MrVar_Get(sne_vnames[0])
	oSne2  = MrVar_Get(sne_vnames[1])
	oSne3  = MrVar_Get(sne_vnames[2])
	oSne4  = MrVar_Get(sne_vnames[3])
	oSfne1 = MrVar_Get(sfne_vnames[0])
	oSfne2 = MrVar_Get(sfne_vnames[1])
	oSfne3 = MrVar_Get(sfne_vnames[2])
	oSfne4 = MrVar_Get(sfne_vnames[3])
	
	;Interpolate to MMS1
	oR1    = oR1    -> Interpol(oSne1)
	oR2    = oR2    -> Interpol(oSne1)
	oR3    = oR3    -> Interpol(oSne1)
	oR4    = oR4    -> Interpol(oSne1)
	oSne2  = oSni2  -> Interpol(oSne1)
	oSne3  = oSni3  -> Interpol(oSne1)
	oSne4  = oSni4  -> Interpol(oSne1)
	oSfne2 = oSfni2 -> Interpol(oSne1)
	oSfne3 = oSfni3 -> Interpol(oSne1)
	oSfne4 = oSfni4 -> Interpol(oSne1)
	
	;Compute the gradient
	;   - Convert 1/km to 1/m
	oRecip   = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	oGradSne = 1e-3 * oRecip -> Gradient(oSne1, oSne2, oSne3, oSne4)
	oDivSfne = 1e-3 * oRecip -> Divergence(oSfne1, oSfne2, oSfne3, oSfne4)
	
	;Attributes
	oGradSne -> SetName, gradsne_vname
	oGradSne -> Cache
	oGradSne['CATDESC']       = 'Gradient of the kinetic energy per particle'
	oGradSne['COLOR']         = ['Blue', 'Forest Green', 'Red']
;	oGradSne['LABEL']         = nabla + 'Sne$\down' + ['X', 'Y', 'Z'] + '$'
	oGradSne['TITLE']         = nabla + 'Sne!C(eV K$\up-1$ s$\up-1$)'
	oGradSne['UNITS']         = 'eV K^-1 m^-1'
	oGradSne['SI_CONVERSION'] = '1.602e-19>J K^-1 m^-1'
	
	oDivSfne -> SetName, divsfne_vname
	oDivSfne -> Cache
	oDivSfne['CATDESC']       = 'Gradient of the kinetic energy per particle'
	oDivSfne['COLOR']         = 'Red'
	oDivSfne['LABEL']         = nabla + '.Sfne'
	oDivSfne['TITLE']         = nabla + '.Sfne!C(eV K$\up-1$ m$\up-1$)'
	oDivSfne['UNITS']         = 'eV K^-1 m^-1'
	oDivSfne['SI_CONVERSION'] = '1.602e-19>J K-1 m^-1'

;-------------------------------------------
; Maxwellian Entropy Density ///////////////
;-------------------------------------------
	;DIS
	;   1. Calculate scalar pressure
	;   2. Calculate P/rho^gamma
	oNi        = MrVar_Get(ni_calc_vname)
	oPi_tensor = MrVar_Get(pi_calc_vname)
	oPi        = (oPi_tensor[*,0,0] + oPi_tensor[*,1,1] + oPi_tensor[*,2,2]) / 3.0
	oSbi       = 3.0/2.0 * kB * J2eV * 1e6 * oNi $
	             * ALog( 2 * !pi * 1e-11 / m_i * oPi['DATA'] / oNi['DATA']^(5.0/3.0) + 1 )
	Obj_Destroy, oPi
	
	;DES
	;   1. Calculate scalar pressure
	;   2. Calculate P/rho^gamma
	oNe        = MrVar_Get(ne_calc_vname)
	oPe_tensor = MrVar_Get(pe_calc_vname)
	oPe        = (oPe_tensor[*,0,0] + oPe_tensor[*,1,1] + oPe_tensor[*,2,2]) / 3.0
	oSbe       = 3.0/2.0 * kB * J2eV * 1e6 * oNe $
	             * ALog( 2 * !pi * 1e-11 / m_e * oPe['DATA'] / oNe['DATA']^(5.0/3.0) + 1 )
	Obj_Destroy, oPe
	
	;Set Attributes
	oSbi -> SetName, sbi_vname
	oSbi -> Cache
	oSbi['CATDESC']       = 'Maxwellian entropy density.'
	oSbi['COLOR']         = 'Blue'
	oSbi['PLOT_TITLE']    = 'Entropy Density'
	oSbi['LABEL']         = 'dis'
	oSbi['TITLE']         = 's$\downB$!C(eV/K/m$\up3$)'
	oSbi['UNITS']         = 'eV/K/m^3 ln(J m^2/kg)'
	oSbi['SI_CONVERSION'] = '1.602e-19>J/K/m^3 ln(J m^2/kg)'
	
	oSbe -> SetName, sbe_vname
	oSbe -> Cache
	oSbe['CATDESC']       = 'Maxwellian entropy density.'
	oSbe['COLOR']         = 'Red'
	oSbe['PLOT_TITLE']    = 'Entropy Density'
	oSbe['LABEL']         = 'des'
	oSbe['TITLE']         = 's$\downB$!C(eV/K/m$\up3$)'
	oSbe['UNITS']         = 'eV/K/m^3 ln(J m^2/kg)'
	oSbe['SI_CONVERSION'] = '1.602e-19>J/K/m^3 ln(J m^2/kg)'

;-------------------------------------------
; Maxwellian Entropy per Particle //////////
;-------------------------------------------
	;1e-6 converts density from cm^3 to m^3
	oSbni = 1e-6 * oSbi / oNi
	oSbne = 1e-6 * oSbe / oNe
	
	oSbni -> SetName, sbni_vname
	oSbni -> Cache
	oSbni['CATDESC']       = 'Maxwellian entropy per particle.'
	oSbni['COLOR']         = 'Blue'
	oSbni['PLOT_TITLE']    = 'Entropy per Particle'
	oSbni['LABEL']         = 'dis'
	oSbni['TITLE']         = 's$\downB$/N!C(eV/K)'
	oSbni['UNITS']         = 'eV/K ln(J m^2/kg)'
	oSbni['SI_CONVERSION'] = '1.602e-19>J/K ln(J m^2/kg)'
	
	oSbne -> SetName, sbne_vname
	oSbne -> Cache
	oSbne['CATDESC']       = 'Maxwellian entropy per particle.'
	oSbne['COLOR']         = 'Red'
	oSbne['PLOT_TITLE']    = 'Entropy per Particle'
	oSbne['LABEL']         = 'des'
	oSbne['TITLE']         = 's$\downB$/N!C(eV/K)'
	oSbne['UNITS']         = 'eV/K ln(J m^2/kg)'
	oSbne['SI_CONVERSION'] = '1.602e-19>J/K ln(J m^2/kg)'
	
;-------------------------------------------
; Agyrotropy Parameter /////////////////////
;-------------------------------------------
	oQi = MrVar_Pres_QFactor(bvec_vname, pi_vname, /CACHE, NAME=qi_vname)
	oQe = MrVar_Pres_QFactor(bvec_vname, pe_vname, /CACHE, NAME=qe_vname)

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
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, coords], ' ' ) )
	
	;DES Entropy Flux
;	oSfne = MrVar_Get(sfne_vnames[isc])
;	oSfne['DATA'] = oSfne['DATA']*10
	
	;Kinetic Entropy per Particle
	oSni = MrVar_Get(sni_vnames[isc])
	oSne = MrVar_Get(sne_vnames[isc])
	oSni['AXIS_RANGE']  = [Min([oSni.min,  oSne.min]),  Max([oSni.max,  oSne.max])]
	oSbni['AXIS_RANGE'] = [Min([oSbni.min, oSbne.min]), Max([oSbni.max, oSbne.max])]
	
	;Maxwellian Entropy Density
	oSi = MrVar_Get(si_vnames[isc])
	oSe = MrVar_Get(se_vnames[isc])
	oSi['AXIS_RANGE']  = [Min([oSbi.min, oSbe.min]), Max([oSbi.max, oSbe.max])]
	oSbi['AXIS_RANGE'] = [Min([oSbi.min, oSbe.min]), Max([oSbi.max, oSbe.max])]
	
	;Entropy Flux
	oSfni = MrVar_Get(sfni_vnames[isc])
	oSfne = MrVar_Get(sfne_vnames[isc])
	oSfni['AXIS_RANGE'] = [Min([oSfni.min, oSfne.min]), Max([oSfni.max, oSfne.max])]
	
	;Gradient
	oGradSni['AXIS_RANGE'] = [Min([oGradSni.min, oGradSne.min]), Max([oGradSni.max, oGradSne.max])]
	
	;Divergence
	oDivSfni['AXIS_RANGE'] = [Min([oDivSfni.min, oDivSfne.min]), Max([oDivSfni.max, oDivSfne.max])]
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [b_vname, sni_vnames[isc], sbni_vname, gradsni_vname, sfni_vnames[isc], divsfni_vname, qe_vname, jdote_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )
	win = MrVar_OPlotTS( sni_vnames[isc], sne_vnames[isc] )
	win = MrVar_OPlotTS( sbni_vname, sbne_vname )
	win = MrVar_OPlotTS( gradsni_vname, gradsne_vname )
	win = MrVar_OPlotTS( sfni_vnames[isc], sfne_vnames[isc] )
	win = MrVar_OPlotTS( divsfni_vname, divsfne_vname )

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
		fname = StrJoin([sc, 'fpi', fpi_mode, level, 'kinetic-entropy'], '_')
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