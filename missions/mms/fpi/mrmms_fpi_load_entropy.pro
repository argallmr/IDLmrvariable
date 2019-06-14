; docformat = 'rst'
;
; NAME:
;       MrMMS_FPI_Load_Entropy
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
;   Compute entropy by integrating the distribution function. The measured entropy and
;   the equivalent Maxwellian entropy are computed.
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
;       EPHDESC:    in, optional, type=string, default='epht89d'
;                   Optional descriptor of the ephemeris file.
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data level
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, no data will be loaded, but VARNAMES will be filled
;       VARNAMES:   out, optional, type=strarr
;                   Names of the variables produced by the program. Order is:
;                        0: n            - density
;                        1: v            - bulk velocity
;                        2. ptensor      - pressure tensor
;                        3: p            - scalar pressure
;                        4: nf           - density flux
;                        5: s            - entropy density
;                        6: sf           - entropy density flux
;                        7: sn           - entropy per particle
;                        8: sfn          - entropy flux per particle
;                        9: sb           - equivalent Maxwellian entropy density
;                       10: sbf          - equivalent Maxwellian entropy density flux
;                       11: sbn          - equivalent Maxwellian entropy per particle
;                       12: sbfn         - equivalent Maxwellian entropy flux per particle
;                       13: mbar         - M = s - sb
;                       14: mfbar        - M = sf - sbf
;                       15: mnbar        - Mn = (s - sb) / n
;                       16: mfnbar       - Mn = (sf - sbf) / n
;                       17: cts_max      - Maximum number of counts in any one bin
;                       18: cts_tot      - Total counts in all bins
;                       
;                       iSC*19+iParam    - Parameter IPARAM for spacecraft ISC
;
;                       -14: Grad(n)     - Density gradient
;                       -13: Grad(s)     - Entropy density gradient
;                       -12: Grad(sn)    - Entropy par particle gradient
;                       -11: Grad(sb)    - Maxwellian entropy density gradient
;                       -10: Grad(sbn)   - Maxwellian entropy per particle gradient
;                        -9: Grad(mbar)  - Non-Maxwellian entropy density gradient
;                        -8: Grad(mnbar) - Non-Maxwellian entropy per particle gradient
;                        -7: Div(nf)     - Density flux divergence
;                        -6: Div(sf)     - Entropy density flux divergence
;                        -5: Div(sfn)    - Entropy flux per particle divergence
;                        -4: Div(sbf)    - Maxwellian entropy flux density divergence
;                        -3: Div(sbfn)   - Maxwellian entropy flux per particle divergence
;                        -2: Div(mfbar)  - Non-Maxwellian entropy flux density divergence
;                        -1: Div(mfnbar) - Non-Maxwellian entropy flux per particle divergence
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
;       2018/06/05  -   Written by Matthew Argall
;       2018/06/19  -   Compute full Boltzmann entropy. - MRA
;-
PRO MrMMS_FPI_Load_Entropy, sc, mode, species, $
EPHDESC=ephdesc, $
LEVEL=level, $
NO_LOAD=no_load, $
VARNAMES=varnames
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN
	ENDIF
	
	IF N_Elements(ephdesc)   EQ 0 THEN ephdesc   = 'epht89d'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	
	;Allow one or four spacecraft
	nSC = N_Elements(sc)
	IF nSC NE 1 && nSC NE 4 THEN Message, 'SC must be a single or all four spacecraft.'
	
	q     = MrConstants('q')
	kB    = MrConstants('k_B')
	m_e   = MrConstants('m_e')
	m_i   = MrConstants('m_H')
	nabla = '!9' + String(71B) + '!X'
	J2eV  = MrConstants('J2eV')
	eV2K  = MrConstants('eV2K')
	eV2J  = MrConstants('eV2J')
	sc_colors   = ['Black', 'Red', 'Forest Green', 'Blue']
	comp_colors = ['Blue', 'Forest Green', 'Red']

;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	fpi_instr  = 'd' + species + 's'
	fpi_mode   = mode EQ 'brst' ? mode : 'fast'
	fpi_coords = 'dbcs'
	
	eph_mode   = mode
	eph_coords = 'gse'
	
	edp_mode = 'fast'

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source Names
	suffix        = '_raw'
	scpot_vnames  = sc + '_' + StrJoin( ['edp', 'scpot', edp_mode, level], '_')
	f_vnames      = sc + '_' + StrJoin( [fpi_instr, 'dist',                      fpi_mode], '_' )
	f_raw_vnames  = sc + '_' + StrJoin( [fpi_instr, 'dist',                      fpi_mode], '_' ) + suffix
	df_raw_vnames = sc + '_' + StrJoin( [fpi_instr, 'disterr',                   fpi_mode], '_' ) + suffix
	r_vnames      = sc + '_' + StrJoin( ['mec', 'r', eph_coords], '_' )
	
	IF species EQ 'i' THEN BEGIN
		mass          = m_i
		type          = 'H'
	ENDIF ELSE BEGIN
		mass          = m_e
		type          = 'e'
	ENDELSE
	
	;Derived names
	cts_tot_vnames  = sc + '_' + StrJoin([fpi_instr, 'counts', 'total',   fpi_mode], '_')
	cts_max_vnames  = sc + '_' + StrJoin([fpi_instr, 'counts', 'max',     fpi_mode], '_')
	mcts_tot_vnames = sc + '_' + StrJoin([fpi_instr, 'counts', 'total', 'maxwellian', fpi_mode], '_')
	mcts_max_vnames = sc + '_' + StrJoin([fpi_instr, 'counts', 'max',   'maxwellian', fpi_mode], '_')
	n_calc_vnames   = sc + '_' + StrJoin([fpi_instr, 'n',         'calc', fpi_mode], '_')
	nf_calc_vnames  = sc + '_' + StrJoin([fpi_instr, 'nf',        'calc', fpi_mode], '_')
	v_calc_vnames   = sc + '_' + StrJoin([fpi_instr, 'v', 'dbcs', 'calc', fpi_mode], '_')
	pt_calc_vnames  = sc + '_' + StrJoin([fpi_instr, 'ptensor',   'calc', fpi_mode], '_')
	p_calc_vnames   = sc + '_' + StrJoin([fpi_instr, 'p',         'calc', fpi_mode], '_')
	tt_calc_vnames  = sc + '_' + StrJoin([fpi_instr, 'ttensor',   'calc', fpi_mode], '_')
	t_calc_vnames   = sc + '_' + StrJoin([fpi_instr, 't',         'calc', fpi_mode], '_')
	
	s_vnames       = sc + '_' + StrJoin([fpi_instr, 's',        fpi_mode], '_')
	sf_vnames      = sc + '_' + StrJoin([fpi_instr, 'sf',       fpi_mode], '_')
	sn_vnames      = sc + '_' + StrJoin([fpi_instr, 'sn',       fpi_mode], '_')
	sfn_vnames     = sc + '_' + StrJoin([fpi_instr, 'sfn',      fpi_mode], '_')
	sb_vnames      = sc + '_' + StrJoin([fpi_instr, 'sb',       fpi_mode], '_')
	sbf_vnames     = sc + '_' + StrJoin([fpi_instr, 'sbf',      fpi_mode], '_')
	sbn_vnames     = sc + '_' + StrJoin([fpi_instr, 'sbn',      fpi_mode], '_')
	sbfn_vnames    = sc + '_' + StrJoin([fpi_instr, 'sbfn',     fpi_mode], '_')
	
	mbar_vnames   = sc + '_' + StrJoin([fpi_instr, 'mbar',   fpi_mode], '_')
	mfbar_vnames  = sc + '_' + StrJoin([fpi_instr, 'mfbar',  fpi_mode], '_')
	mnbar_vnames  = sc + '_' + StrJoin([fpi_instr, 'mnbar',  fpi_mode], '_')
	mfnbar_vnames = sc + '_' + StrJoin([fpi_instr, 'mfnbar', fpi_mode], '_')
	mfrac_vnames  = sc + '_' + StrJoin([fpi_instr, 'mfrac',  fpi_mode], '_')
	
	gradn_vname     = StrJoin(['mms', fpi_instr, 'grad', 'n',      fpi_mode], '_')
	grads_vname     = StrJoin(['mms', fpi_instr, 'grad', 's',      fpi_mode], '_')
	gradsn_vname    = StrJoin(['mms', fpi_instr, 'grad', 'sn',     fpi_mode], '_')
	gradsb_vname    = StrJoin(['mms', fpi_instr, 'grad', 'sb',     fpi_mode], '_')
	gradsbn_vname   = StrJoin(['mms', fpi_instr, 'grad', 'sbn',    fpi_mode], '_')
	gradm_vname     = StrJoin(['mms', fpi_instr, 'grad', 'mbar',   fpi_mode], '_')
	gradmn_vname    = StrJoin(['mms', fpi_instr, 'grad', 'mnbar',  fpi_mode], '_')
	gradmfrac_vname = StrJoin(['mms', fpi_instr, 'grad', 'mfrac',  fpi_mode], '_')
	divnf_vname     = StrJoin(['mms', fpi_instr, 'div',  'nf',     fpi_mode], '_')
	divsf_vname     = StrJoin(['mms', fpi_instr, 'div',  'sf',     fpi_mode], '_')
	divsfn_vname    = StrJoin(['mms', fpi_instr, 'div',  'sfn',    fpi_mode], '_')
	divsbf_vname    = StrJoin(['mms', fpi_instr, 'div',  'sbf',    fpi_mode], '_')
	divsbfn_vname   = StrJoin(['mms', fpi_instr, 'div',  'sbfn',   fpi_mode], '_')
	divmf_vname     = StrJoin(['mms', fpi_instr, 'div',  'mfbar',  fpi_mode], '_')
	divmfn_vname    = StrJoin(['mms', fpi_instr, 'div',  'mfnbar', fpi_mode], '_')
	
	s_tot_vnames    = sc + '_' + StrJoin([fpi_instr, 's',    'tot', fpi_mode], '_')
	sf_tot_vnames   = sc + '_' + StrJoin([fpi_instr, 'sf',   'tot', fpi_mode], '_')
	sn_tot_vnames   = sc + '_' + StrJoin([fpi_instr, 'sn',   'tot', fpi_mode], '_')
	sfn_tot_vnames  = sc + '_' + StrJoin([fpi_instr, 'sfn',  'tot', fpi_mode], '_')
	sb_tot_vnames   = sc + '_' + StrJoin([fpi_instr, 'sb',   'tot', fpi_mode], '_')
	sbf_tot_vnames  = sc + '_' + StrJoin([fpi_instr, 'sbf',  'tot', fpi_mode], '_')
	sbn_tot_vnames  = sc + '_' + StrJoin([fpi_instr, 'sbn',  'tot', fpi_mode], '_')
	sbfn_tot_vnames = sc + '_' + StrJoin([fpi_instr, 'sbfn', 'tot', fpi_mode], '_')
	
	mbar_tot_vnames   = sc + '_' + StrJoin([fpi_instr, 'mbar',   'tot', fpi_mode], '_')
	mfbar_tot_vnames  = sc + '_' + StrJoin([fpi_instr, 'mfbar',  'tot', fpi_mode], '_')
	mnbar_tot_vnames  = sc + '_' + StrJoin([fpi_instr, 'mnbar',  'tot', fpi_mode], '_')
	mfnbar_tot_vnames = sc + '_' + StrJoin([fpi_instr, 'mfnbar', 'tot', fpi_mode], '_')
	mfrac_tot_vnames  = sc + '_' + StrJoin([fpi_instr, 'mfrac',  'tot', fpi_mode], '_')
	
	gradn_tot_vname     = StrJoin(['mms', fpi_instr, 'grad', 'n',      'tot', fpi_mode], '_')
	grads_tot_vname     = StrJoin(['mms', fpi_instr, 'grad', 's',      'tot', fpi_mode], '_')
	gradsn_tot_vname    = StrJoin(['mms', fpi_instr, 'grad', 'sn',     'tot', fpi_mode], '_')
	gradsb_tot_vname    = StrJoin(['mms', fpi_instr, 'grad', 'sb',     'tot', fpi_mode], '_')
	gradsbn_tot_vname   = StrJoin(['mms', fpi_instr, 'grad', 'sbn',    'tot', fpi_mode], '_')
	gradm_tot_vname     = StrJoin(['mms', fpi_instr, 'grad', 'mbar',   'tot', fpi_mode], '_')
	gradmn_tot_vname    = StrJoin(['mms', fpi_instr, 'grad', 'mnbar',  'tot', fpi_mode], '_')
	gradmfrac_tot_vname = StrJoin(['mms', fpi_instr, 'grad', 'mfrac',  'tot', fpi_mode], '_')
	divnf_tot_vname     = StrJoin(['mms', fpi_instr, 'div',  'nf',     'tot', fpi_mode], '_')
	divsf_tot_vname     = StrJoin(['mms', fpi_instr, 'div',  'sf',     'tot', fpi_mode], '_')
	divsfn_tot_vname    = StrJoin(['mms', fpi_instr, 'div',  'sfn',    'tot', fpi_mode], '_')
	divsbf_tot_vname    = StrJoin(['mms', fpi_instr, 'div',  'sbf',    'tot', fpi_mode], '_')
	divsbfn_tot_vname   = StrJoin(['mms', fpi_instr, 'div',  'sbfn',   'tot', fpi_mode], '_')
	divmf_tot_vname     = StrJoin(['mms', fpi_instr, 'div',  'mfbar',  'tot', fpi_mode], '_')
	divmfn_tot_vname    = StrJoin(['mms', fpi_instr, 'div',  'mfnbar', 'tot', fpi_mode], '_')

;-------------------------------------------
; Output Names /////////////////////////////
;-------------------------------------------
	
	varnames = [ [n_calc_vnames],    [v_calc_vnames],   [pt_calc_vnames],   [p_calc_vnames], [nf_calc_vnames], $
	             [cts_max_vnames],   [cts_tot_vnames],  [mcts_max_vnames],  [mcts_tot_vnames], $
	             [s_vnames],         [sf_vnames],       [sn_vnames],        [sfn_vnames], $
	             [sb_vnames],        [sbf_vnames],      [sbn_vnames],       [sbfn_vnames], $
	             [mfrac_vnames],     [mbar_vnames],     [mfbar_vnames],     [mnbar_vnames],      [mfnbar_vnames], $
	             [s_tot_vnames],     [sf_tot_vnames],   [sn_tot_vnames],    [sfn_tot_vnames], $
	             [sb_tot_vnames],    [sbf_tot_vnames],  [sbn_tot_vnames],   [sbfn_tot_vnames], $
	             [mfrac_tot_vnames], [mbar_tot_vnames], [mfbar_tot_vnames], [mnbar_tot_vnames],  [mfnbar_tot_vnames] ]
	varnames = Reform(Transpose(varnames), N_Elements(varnames))
	IF nSC EQ 4 THEN BEGIN
		varnames = [varnames, $
		            gradn_vname, $
		            grads_vname, gradsn_vname, gradsb_vname, gradsbn_vname, $
		            gradm_vname, gradmn_vname, $
		            divnf_vname, $
		            divsf_vname, divsfn_vname, divsbf_vname, divsbfn_vname, $
		            divmf_vname, divmfn_vname]
	ENDIF
	
	;RETURN if we only want the variable names
	IF Keyword_Set(no_load) THEN RETURN
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	;FPI-dist
	FOR i = 0, nSC - 1 $
		DO MrMMS_FPI_Load_Dist3D, sc[i], mode, species, /APPLY_MODEL
	
	;Distribution
	MrMMS_Load_Data, sc, 'fpi', mode, level, $
	                 OPTDESC   = fpi_instr + '-dist', $
	                 SUFFIX    = suffix, $
	                 VARFORMAT = '*_dist*'
	
	;Spacecraft potential
	MrMMS_Load_Data, sc, 'edp', 'fast', level, $
	                 OPTDESC   = 'scpot', $
	                 VARFORMAT = '*scpot*'
	
	;Ephemeris
	IF nSC EQ 4 THEN BEGIN
		IF mode EQ 'brst' THEN BEGIN
			;BRST data is not always available
			IF (MrMMS_Get_FileNames(sc, 'mec', mode, level, OPTDESC=ephdesc))[0] EQ '' $
				THEN eph_mode = 'srvy'
		ENDIF
		MrMMS_Load_Data, sc, 'mec', eph_mode, level, $
		                 OPTDESC   = ephdesc, $
		                 VARFORMAT = '*r_'+eph_coords
	ENDIF
	
;-------------------------------------------
; Loop Over Each Spacecraft ////////////////
;-------------------------------------------
	FOR i = 0, nSC - 1 DO BEGIN
		oDist = MrVar_Get(f_vnames[i])
		dims = Size(oDist, /DIMENSIONS)
		
	;-------------------------------------------
	; Moments //////////////////////////////////
	;-------------------------------------------
		oDF = MrDist4D(f_vnames[i], VSC=scpot_vnames[i], SPECIES=type)
		oN  = oDF -> Density(/CACHE, NAME=n_calc_vnames[i])
		oV  = oDF -> Velocity(/CACHE, NAME=v_calc_vnames[i])
		oPt = oDF -> Pressure(/CACHE, NAME=pt_calc_vnames[i])
		oP  = (oPt[*,0,0] + oPt[*,1,1] + oPt[*,2,2]) / 3.0
		oP -> SetName, p_calc_vnames[i]
		oP -> Cache
		oTt = oDF -> Temperature(/CACHE, NAME=tt_calc_vnames[i])
		oT  = (oTt[*,0,0] + oTt[*,1,1] + oTt[*,2,2]) / 3.0
		oT -> SetName, t_calc_vnames[i]
		oT -> Cache
		
		;Density
		oN['LABEL'] = nSC EQ 4 ? sc[i] : 'N' + species
		IF nSC EQ 4 THEN oN['COLOR'] = sc_colors[i]
		
		;Bulk Velocity
		oV['COLOR'] = comp_colors
		oV['LABEL'] = nSC EQ 4 ? sc[i] : 'N' + species
		
		;Scalar Pressure
		oP['LABEL'] = nSC EQ 4 ? sc[i] : 'P' + species
		IF nSC EQ 4 THEN oP['COLOR'] = sc_colors[i]
		
	;-------------------------------------------
	; Density Flux /////////////////////////////
	;-------------------------------------------
		;   - Convert km/s to m/s
		oNf = 1e3 * oV * oN
		
		;S
		oNf -> SetName, nf_calc_vnames[i]
		oNf -> Cache
		oNf['CATDESC']       = 'Density flux computed via integration of the full velocity ' + $
		                       'space distribution function.'
		oNf['COLOR']         = comp_colors
		oNf['PLOT_TITLE']    = 'Density flux'
		oNf['LABEL']         = nSC EQ 4 ? sc[i] : 'Nf'
		oNf['TITLE']         = 'N!C(1/cm$\up3$'
		oNf['UNITS']         = '1/cm^3'
		oNf['SI_CONVERSION'] = '1e6>m^-3'
	
	;-------------------------------------------
	; Entropy Density //////////////////////////
	;-------------------------------------------
		;   - Convert J/K to eV/K
		oS  = J2eV * oDF -> Entropy()
		
		;S
		oS -> SetName, s_vnames[i]
		oS -> Cache
		oS['CATDESC']       = 'Entropy density computed via integration of the full velocity ' + $
		                      'space distribution function.'
		oS['LABEL']         = nSC EQ 4 ? sc[i] : 's'
		oS['PLOT_TITLE']    = 'Entropy Density'
		oS['TITLE']         = 's!C(eV/K/m$\up3$)'
		oS['UNITS']         = 'eV/K/m^3 ln(s^3/m^6)'
		oS['SI_CONVERSION'] = '1.602e-19>J/K/m^3'
		IF nSC EQ 4 THEN oS['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Entropy Density Flux /////////////////////
	;-------------------------------------------
		;   - Convert km/s to m/s
		oSf = 1e3 * oV * oS
		
		;S
		oSf -> SetName, sf_vnames[i]
		oSf -> Cache
		oSf['CATDESC']       = 'Entropy density flux computed via integration of the full velocity ' + $
		                      'space distribution function.'
		oSf['COLOR']         = comp_colors
		oSf['LABEL']         = nSC EQ 4 ? sc[i] : 'sf'
		oSf['PLOT_TITLE']    = 'Entropy Density Flux'
		oSf['TITLE']         = 'sf!C(eV/K/m$\up2$/s)'
		oSf['UNITS']         = 'eV/K/m^2/s ln(s^3/m^6)'
		oSf['SI_CONVERSION'] = '1.602e-19>J/K/m^2/s'
	
	;-------------------------------------------
	; Entropy per Particle /////////////////////
	;-------------------------------------------
		;Convert 1/cm^3 to 1/m^3
		oSn = oS / (1e6 * oN)
		
		;Sn
		oSn -> SetName, sn_vnames[i]
		oSn -> Cache
		oSn['CATDESC']       = 'Entropy per particle computed via integration of the full velocity ' + $
		                       'space distribution function.'
		oSn['LABEL']         = nSC EQ 4 ? sc[i] : 's/N'
		oSn['PLOT_TITLE']    = 'Entropy per Particle'
		oSn['TITLE']         = 's/N!C(eV/K)'
		oSn['UNITS']         = 'eV/K ln(s^3/m^6)'
		oSn['SI_CONVERSION'] = '1.602e-19>J/K'
		IF nSC EQ 4 THEN oSn['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Entropy Flux per Particle ////////////////
	;-------------------------------------------
		;Convert km/s to m/s
		oSfn = 1e3 * oV * oSn
		
		;Sn
		oSfn -> SetName, sfn_vnames[i]
		oSfn -> Cache
		oSfn['CATDESC']       = 'Entropy flux per particle computed via integration of the full velocity ' + $
		                        'space distribution function.'
		oSfn['COLOR']         = comp_colors
		oSfn['LABEL']         = nSC EQ 4 ? sc[i] : 'sf/N'
		oSfn['PLOT_TITLE']    = 'Entropy Flux per Particle'
		oSfn['TITLE']         = 'sf/N!C(eV/K)'
		oSfn['UNITS']         = 'eV m / K s ln(s^3/m^6)'
		oSfn['SI_CONVERSION'] = '1.602e-19>J m / K s'
	
	;-------------------------------------------
	; Maxwellian Entropy Density ///////////////
	;-------------------------------------------
		oSb = 3.0/2.0 * kB * J2eV * 1e6 * oN $
		      * ALog( 2 * !pi * 1e-11 / mass * oP['DATA'] / oN['DATA']^(5.0/3.0) + 1 )
	
		;Set Attributes
		oSb -> SetName, sb_vnames[i]
		oSb -> Cache
		oSb['CATDESC']       = 'Maxwellian entropy density calculated with the equivalent ' + $
		                       'density and temperature of the measured distribution.'
		oSb['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downB$'
		oSb['PLOT_TITLE']    = 'Maxwellian Entropy Density'
		oSb['TITLE']         = 's$\downB$!C(eV/K/m$\up3$)'
		oSb['UNITS']         = 'eV/K/m^3 ln(J m^2/kg)'
		oSb['SI_CONVERSION'] = '1.602e-19>J/K/m^3 ln(J m^2/kg)'
		IF nSC EQ 4 THEN oSb['COLOR'] = sc_colors[i]
		
	;-------------------------------------------
	; Maxwellian Entropy Density Flux //////////
	;-------------------------------------------
		;Convert km/s to m/s
		oSbf = 1e3 * oV * oSb
	
		;Set Attributes
		oSbf -> SetName, sbf_vnames[i]
		oSbf -> Cache
		oSbf['CATDESC']       = 'Maxwellian entropy density flux calculated with the equivalent ' + $
		                        'density and temperature of the measured distribution.'
		oSbf['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downBf$'
		oSbf['COLOR']         = comp_colors
		oSbf['PLOT_TITLE']    = 'Entropy Density Flux'
		oSbf['TITLE']         = 's$\downB,f$!C(eV/K/m$\up2$/s)'
		oSbf['UNITS']         = 'eV/K/m^2/s ln(J m^2/kg)'
		oSbf['SI_CONVERSION'] = '1.602e-19>J/K/m^2/s ln(J m^2/kg)'
	
	;-------------------------------------------
	; Maxwellian Entropy per Particle //////////
	;-------------------------------------------
		oSbn = oSb / (1e6 * oN)
		
		;Set Attributes
		oSbn -> SetName, sbn_vnames[i]
		oSbn -> Cache
		oSbn['CATDESC']       = 'Maxwellian entropy per particle.'
		oSbn['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downB$/N'
		oSbn['PLOT_TITLE']    = 'Maxwellian Entropy per Particle'
		oSbn['TITLE']         = 's$\downB$/N!C(eV/K)'
		oSbn['UNITS']         = 'eV/K ln(J m^2/kg)'
		oSbn['SI_CONVERSION'] = '1.602e-19>J/K ln(J m^2/kg)'
		IF nSC EQ 4 THEN oSbn['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Maxwellian Entropy Flux per Particle /////
	;-------------------------------------------
		;Convert km/s to m/s
		oSbfn = 1e3 * oV * oSbn
		
		;Set Attributes
		oSbfn -> SetName, sbfn_vnames[i]
		oSbfn -> Cache
		oSbfn['CATDESC']       = 'Maxwellian entropy flux per particle.'
		oSbfn['COLOR']         = comp_colors
		oSbfn['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downBf$/N'
		oSbfn['PLOT_TITLE']    = 'Maxwellian Entropy Flux per Particle'
		oSbfn['TITLE']         = 's$\downB,f$/N!C(eV/K m/s)'
		oSbfn['UNITS']         = 'eV/K m/s ln(J m^2/kg)'
		oSbfn['SI_CONVERSION'] = '1.602e-19>J/K m/s ln(J m^2/kg)'
	
	;-------------------------------------------
	; M-Bar Fraction ///////////////////////////
	;-------------------------------------------
		oMfrac = (oSb - oS) / oSb
		oMfrac -> SetName, mfrac_vnames[i]
		oMfrac -> Cache
		
		;Set Attributes
		oMfrac['AXIS_RANGE']    = [0.0, 1.0]
		oMfrac['CATDESC']       = 'Fractional difference between the entropy density of the ' + $
		                          'full distribution and that of an equivalent Maxwellian.'
		oMfrac['PLOT_TITLE']    = 'Fractional Non-Maxwellian Entropy Density'
		oMfrac['LABEL']         = nSC EQ 4 ? sc[i] : 'M'
		oMfrac['TITLE']         = 'M'
		oMfrac['UNITS']         = ''
		oMfrac['SI_CONVERSION'] = '>'
		IF nSC EQ 4 THEN oMfrac['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; M-Bar ////////////////////////////////////
	;-------------------------------------------
		oMbar = oS - oSb
		oMbar -> SetName, mbar_vnames[i]
		oMbar -> Cache
		
		;Set Attributes
		oMbar['CATDESC']       = 'Difference between the entropy density of the ' + $
		                         'full distribution and that of an equivalent Maxwellian.'
		oMbar['PLOT_TITLE']    = 'Non-Maxwellian Entropy Density'
		oMbar['LABEL']         = nSC EQ 4 ? sc[i] : 'M'
		oMbar['TITLE']         = 's-s$\downB$!C(eV/K/m$\up3$)'
		oMbar['UNITS']         = 'eV/K/m^3 [ln(s^3/m^6) - ln(J m^2/kg)]'
		oMbar['SI_CONVERSION'] = '1.602e-19>J/K/m^3 [ln(s^3/m^6) - ln(J m^2/kg)]'
		IF nSC EQ 4 THEN oMbar['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; M-Bar Flux ///////////////////////////////
	;-------------------------------------------
		oMfbar = 1e3 * oV * oMbar
		oMfbar -> SetName, mfbar_vnames[i]
		oMfbar -> Cache
		
		;Set Attributes
		oMfbar['CATDESC']       = 'Difference between the entropy density flux of the ' + $
		                          'full distribution and that of an equivalent Maxwellian.'
		oMfbar['COLOR']         = comp_colors
		oMfbar['PLOT_TITLE']    = 'Non-Maxwellian Entropy Density Flux'
		oMfbar['LABEL']         = nSC EQ 4 ? sc[i] : 'Mf'
		oMfbar['TITLE']         = 's-s$\downB$!C(eV/K/m$\up2$/s)'
		oMfbar['UNITS']         = 'eV/K/m^2/s [ln(s^3/m^6) - ln(J m^2/kg)]'
		oMfbar['SI_CONVERSION'] = '1.602e-19>J/K/m^2/s [ln(s^3/m^6) - ln(J m^2/kg)]'
	
	;-------------------------------------------
	; M-Bar Per Particle ///////////////////////
	;-------------------------------------------
		oMnbar = oMbar / (1e6 * oN)
		oMnbar -> SetName, mnbar_vnames[i]
		oMnbar -> Cache
		
		;Set Attributes
		oMnbar['CATDESC']       = 'Difference between the entropy per particle of the ' + $
		                          'full distribution and that of an equivalent Maxwellian.'
		oMnbar['PLOT_TITLE']    = 'Non-Maxwellian Entropy per Particle'
		oMnbar['LABEL']         = nSC EQ 4 ? sc[i] : 'M/N'
		oMnbar['TITLE']         = '(s-s$\downB$)/N!C(eV/K)'
		oMnbar['UNITS']         = 'eV/K [ln(s^3/m^6) - ln(J m^2/kg)]'
		oMnbar['SI_CONVERSION'] = '1.602e-19>J/K [ln(s^3/m^6) - ln(J m^2/kg)]'
		IF nSC EQ 4 THEN oMnbar['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; M-Bar Flux Per Particle //////////////////
	;-------------------------------------------
		oMfnbar = 1e3 * oV * oMnBar
		oMfnbar -> SetName, mfnbar_vnames[i]
		oMfnbar -> Cache
		
		;Set Attributes
		oMfnbar['CATDESC']       = 'Difference between the entropy flux per particle of the ' + $
		                           'full distribution and that of an equivalent Maxwellian.'
		oMfnbar['COLOR']         = comp_colors
		oMfnbar['PLOT_TITLE']    = 'Non-Maxwellian Entropy Flux per Particle'
		oMfnbar['LABEL']         = nSC EQ 4 ? sc[i] : 'Mf/N'
		oMfnbar['TITLE']         = 'Mf/N!C(eV/K m/s)'
		oMfnbar['UNITS']         = 'eV/K m/s [ln(s^3/m^6) - ln(J m^2/kg)]'
		oMfnbar['SI_CONVERSION'] = '1.602e-19>J/K m/s [ln(s^3/m^6) - ln(J m^2/kg)]'
	
	;-------------------------------------------
	; Boltzmann Entropy Per Particle ///////////
	;-------------------------------------------
		oF_dist  = MrVar_Get(f_raw_vnames[i])
		odF_dist = MrVar_Get(df_raw_vnames[i])
		oSCPot   = MrVar_Get(scpot_vnames[i])
		
		;Compute counts. Careful of zeros!
		iZero  = odF_dist -> Where(0.0, /EQUAL, COUNT=nZero)
		counts = Round((oF_dist['DATA'] / odF_dist['DATA'])^2.0)
		IF nZero GT 0 THEN counts[iZero] = 0
		
		;Set counts in bins with energy less than the spacecraft potential equal to zero
		oE = oF_dist['DEPEND_3']
		FOR iE = 0, dims[3] - 1 DO BEGIN
			iBad = Where(oE['DATA',*,iE] LT oSCPot['DATA'], nBad)
			IF nBad GT 0 THEN counts[iBad,*,*,iE] = 0
		ENDFOR
		
		;Boltzmann Entropy
		counts = Reform(counts, [dims[0], Product(dims[1:3])], /OVERWRITE)
		sn_tot = J2eV * kB * (LNGamma(Total(counts, 2)+1) - Total(LNGamma(counts+1), 2))
		
		;Create variable
		oSn_tot = MrScalarTS( oF_dist['TIMEVAR'], sn_tot, $
		                      /CACHE, $
		                      NAME = sn_tot_vnames[i], $
		                      /NO_COPY )
		oSn_tot['CATDESC']       = 'Full Boltzmann entropy per particle.'
		oSn_tot['PLOT_TITLE']    = 'Boltzmann Entropy per Particle'
		oSn_tot['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downT$/N'
		oSn_tot['TITLE']         = 's$\downT$/N!C(eV/K)'
		oSn_tot['UNITS']         = 'eV/K'
		oSn_tot['SI_CONVERSION'] = '1.602e-19>J/K'
	
	;-------------------------------------------
	; Boltzmann Entropy Flux Per Particle //////
	;-------------------------------------------
		;Convert km/s to m/s
		oSnf_tot = 1e3 * oV * oSn_tot
		oSnf_tot -> SetName, sfn_tot_vnames[i]
		oSnf_tot -> Cache
		
		;Create variable
		oSnf_tot['CATDESC']       = 'Boltzmann entropy flux per particle.'
		oSnf_tot['PLOT_TITLE']    = 'Boltzmann Entropy Flux per Particle'
		oSnf_tot['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downTf$/N'
		oSnf_tot['TITLE']         = 's$\downTf$/N!C(eV/K m/s)'
		oSnf_tot['UNITS']         = 'eV/K m/s'
		oSnf_tot['SI_CONVERSION'] = '1.602e-19>J/K m/s'
	
	;-------------------------------------------
	; Boltzmann Entropy Density ////////////////
	;-------------------------------------------
		;Convert cm^-3 to m^-3
		oS_tot = 1e6 * oN * oSn_tot
		oS_tot -> SetName, s_tot_vnames[i]
		oS_tot -> Cache
		
		;Create variable
		oS_tot['CATDESC']       = 'Full Boltzmann entropy density.'
		oS_tot['PLOT_TITLE']    = 'Boltzmann Entropy Density'
		oS_tot['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downT$'
		oS_tot['TITLE']         = 's!C(eV/K/m$\up3$)'
		oS_tot['UNITS']         = 'eV/K/m^3'
		oS_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^3'
	
	;-------------------------------------------
	; Boltzmann Entropy Density Flux ///////////
	;-------------------------------------------
		;Convert km/s to m/s
		oSf_tot = 1e3 * oV * oS_tot
		oSf_tot -> SetName, sf_tot_vnames[i]
		oSf_tot -> Cache
		
		;Create variable
		oSf_tot['CATDESC']       = 'Boltzmann entropy density flux.'
		oSf_tot['PLOT_TITLE']    = 'Boltzmann Entropy Density Flux'
		oSf_tot['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downTf$'
		oSf_tot['TITLE']         = 's$\downTf$!C(eV/K/m^2/s)'
		oSf_tot['UNITS']         = 'eV/K/m^2/s'
		oSf_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^2/s'
	
	;-------------------------------------------
	; Maxwell-Boltzmann Entropy Per Particle ///
	;-------------------------------------------
		;Compute the Maxwellian distribution
		oF_max  = MrMMS_FPI_F_Maxwellian(oF_dist, type, DENSITY=oN, VELOCITY=oV, TEMPERATURE=oT)
		
		;Compute counts. Careful of zeros!
		;   - Convert from s^3/cm^6 to s^3/km^6 to prevent underflow errors
		iZero      = odF_dist -> Where(0.0, /EQUAL, COUNT=nZero)
		counts_max = Round( ( (1e30*oF_dist['DATA']) / (1e30*odF_dist['DATA'])^2 ) * (1e30*oF_max['DATA']) )
		IF nZero GT 0 THEN counts_max[iZero] = 0
		
		;Total counts
		counts_max = Reform(counts_max, [dims[0], Product(dims[1:3])], /OVERWRITE)
		stot_max   = J2eV * kB * (LNGamma(Total(counts_max, 2)+1) - Total(LNGamma(counts_max+1), 2))
		
		;Create variable
		oSbn_tot = MrScalarTS( oF_dist['TIMEVAR'], stot_max, $
		                       /CACHE, $
		                       NAME = sbn_tot_vnames[i], $
		                       /NO_COPY )
		oSbn_tot['CATDESC']       = 'Full Boltzmann entropy per particle of a Maxwellian ' + $
		                            'distribution with the same density and temperature of ' + $
		                            'the measured distribution.'
		oSbn_tot['PLOT_TITLE']    = 'Boltzmann Entropy per Particle of a Maxwellian'
		oSbn_tot['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downBT$/N'
		oSbn_tot['TITLE']         = 's$\downBT$/N!C(eV/K)'
		oSbn_tot['UNITS']         = 'eV/K'
		oSbn_tot['SI_CONVERSION'] = '1.602e-19>J/K'
	
	;--------------------------------------------
	; Maxwell-Boltzmann Entropy Flux Per Particle
	;--------------------------------------------
		;Convert km/s to m/s
		oSbfn_tot = 1e3 * oV * oSbn_tot

		;Attributes
		oSbfn_tot['CATDESC']       = 'Boltzmann entropy flux per particle of a Maxwellian ' + $
		                             'distribution with the same density and temperature of ' + $
		                             'the measured distribution.'
		oSbfn_tot['PLOT_TITLE']    = 'Boltzmann Entropy Flux per Particle of a Maxwellian'
		oSbfn_tot['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downBTf$/N'
		oSbfn_tot['TITLE']         = 's$\downBTf$/N!C(eV/K m/s)'
		oSbfn_tot['UNITS']         = 'eV/K m/s'
		oSbfn_tot['SI_CONVERSION'] = '1.602e-19>J/K m/s'
	
	;-------------------------------------------
	; Maxwell-Boltzmann Entropy Density ////////
	;-------------------------------------------
		;Convert cm^-3 to m^-3
		oSb_tot = 1e6 * oN * oSbn_tot
		oSb_tot -> SetName, sb_tot_vnames[i]
		oSb_tot -> Cache
		
		;Attributes
		oSb_tot['CATDESC']       = 'Full Boltzmann entropy density of a Maxwellian ' + $
		                           'distribution with the same density and temperature of ' + $
		                           'the measured distribution.'
		oSb_tot['PLOT_TITLE']    = 'Boltzmann Entropy Density of a Maxwellian'
		oSb_tot['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downBT$'
		oSb_tot['TITLE']         = 's$\downBT$!C(eV/K/m^3)'
		oSb_tot['UNITS']         = 'eV/K/m^3'
		oSb_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^3'
	
	;--------------------------------------------
	; Maxwell-Boltzmann Entropy Density Flux ////
	;--------------------------------------------
		;Convert km/s to m/s
		oSbf_tot = 1e3 * oV * oSb_tot

		;Attributes
		oSbf_tot['CATDESC']       = 'Boltzmann entropy density flux of a Maxwellian ' + $
		                            'distribution with the same density and temperature of ' + $
		                            'the measured distribution.'
		oSbf_tot['PLOT_TITLE']    = 'Boltzmann Entropy Density Flux of a Maxwellian'
		oSbf_tot['LABEL']         = nSC EQ 4 ? sc[i] : 's$\downBTf$'
		oSbf_tot['TITLE']         = 's$\downBTf$!C(eV/K/m^2/s)'
		oSbf_tot['UNITS']         = 'eV/K/m^2/s'
		oSbf_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^2/s'
	
	;-------------------------------------------
	; Boltzmann M-Bar Fraction /////////////////
	;-------------------------------------------
		oMfrac_tot = (oSb_tot - oS_tot) / oSb_tot
		oMfrac_tot -> SetName, mfrac_tot_vnames[i]
		oMfrac_tot -> Cache
		
		;Set Attributes
		oMfrac_tot['AXIS_RANGE']    = [0.0, 1.0]
		oMfrac_tot['CATDESC']       = 'Fractional difference between the entropy density of the ' + $
		                              'full distribution and that of an equivalent Maxwellian.'
		oMfrac_tot['PLOT_TITLE']    = 'Fractional Non-Maxwellian Entropy Density'
		oMfrac_tot['LABEL']         = nSC EQ 4 ? sc[i] : 'M$\downT$'
		oMfrac_tot['TITLE']         = 'M$\downT$'
		oMfrac_tot['UNITS']         = ''
		oMfrac_tot['SI_CONVERSION'] = '>'
		IF nSC EQ 4 THEN oMfrac_tot['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Boltzmann M-Bar //////////////////////////
	;-------------------------------------------
		oMbar_tot = oS_tot - oSb_tot
		oMbar_tot -> SetName, mbar_tot_vnames[i]
		oMbar_tot -> Cache
		
		;Set Attributes
		oMbar_tot['CATDESC']       = 'Difference between the Boltzmann entropy density of the ' + $
		                             'full distribution and that of an equivalent Maxwellian.'
		oMbar_tot['PLOT_TITLE']    = 'Non-Maxwellian Entropy Density'
		oMbar_tot['LABEL']         = nSC EQ 4 ? sc[i] : 'M'
		oMbar_tot['TITLE']         = 'M$\downT$$!C(eV/K/m$\up3$)'
		oMbar_tot['UNITS']         = 'eV/K/m^3'
		oMbar_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^3'
		IF nSC EQ 4 THEN oMbar_tot['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Boltzmann M-Bar Flux /////////////////////
	;-------------------------------------------
		oMfbar_tot = 1e3 * oV * oMbar_tot
		oMfbar_tot -> SetName, mfbar_tot_vnames[i]
		oMfbar_tot -> Cache
		
		;Set Attributes
		oMfbar_tot['CATDESC']       = 'Difference between the Boltzmann entropy density flux of the ' + $
		                              'full distribution and that of an equivalent Maxwellian.'
		oMfbar_tot['COLOR']         = comp_colors
		oMfbar_tot['PLOT_TITLE']    = 'Non-Maxwellian Entropy Density Flux'
		oMfbar_tot['LABEL']         = nSC EQ 4 ? sc[i] : 'Mf'
		oMfbar_tot['TITLE']         = 'M$\downTf$!C(eV/K/m$\up2$/s)'
		oMfbar_tot['UNITS']         = 'eV/K/m^2/s'
		oMfbar_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^2/s'
	
	;-------------------------------------------
	; Boltzmann M-Bar Per Particle /////////////
	;-------------------------------------------
		oMnbar_tot = oMbar_tot / (1e6 * oN)
		oMnbar_tot -> SetName, mnbar_tot_vnames[i]
		oMnbar_tot -> Cache
		
		;Set Attributes
		oMnbar_tot['CATDESC']       = 'Difference between the Boltzmann entropy per particle of the ' + $
		                              'full distribution and that of an equivalent Maxwellian.'
		oMnbar_tot['PLOT_TITLE']    = 'Non-Maxwellian Entropy per Particle'
		oMnbar_tot['LABEL']         = nSC EQ 4 ? sc[i] : 'M/N'
		oMnbar_tot['TITLE']         = 'M$\downT$/N!C(eV/K)'
		oMnbar_tot['UNITS']         = 'eV/K'
		oMnbar_tot['SI_CONVERSION'] = '1.602e-19>J/K'
		IF nSC EQ 4 THEN oMnbar_tot['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Boltzmann M-Bar Flux Per Particle ////////
	;-------------------------------------------
		oMfnbar_tot = 1e3 * oV * oMnBar_tot
		oMfnbar_tot -> SetName, mfnbar_tot_vnames[i]
		oMfnbar_tot -> Cache
		
		;Set Attributes
		oMfnbar_tot['CATDESC']       = 'Difference between the Boltzmann entropy flux per particle of the ' + $
		                               'full distribution and that of an equivalent Maxwellian.'
		oMfnbar_tot['COLOR']         = comp_colors
		oMfnbar_tot['PLOT_TITLE']    = 'Non-Maxwellian Entropy Flux per Particle'
		oMfnbar_tot['LABEL']         = nSC EQ 4 ? sc[i] : 'M$\downTf$/N'
		oMfnbar_tot['TITLE']         = 'M$\downTf$/N!C(eV/K m/s)'
		oMfnbar_tot['UNITS']         = 'eV/K m/s'
		oMfnbar_tot['SI_CONVERSION'] = '1.602e-19>J/K m/s'
	
	;-------------------------------------------
	; Maximum Counts ///////////////////////////
	;-------------------------------------------
		ctemp = Float(counts)
		IF nZero GT 0 THEN ctemp[iZero] = !Values.F_NaN
		nmax = Fix(Max(Temporary(ctemp), DIMENSION=2, NAN=(nZero GT 0)), TYPE=3)
		oNmax = MrScalarTS( oF_dist['TIMEVAR'], nmax, $
		                    /CACHE, $
		                    NAME = cts_max_vnames[i], $
		                    /NO_COPY )
		
		;Set Attributes
		oNmax['CATDESC']       = 'Maximum number of counts in any single bin of the distribution function.'
		oNmax['PLOT_TITLE']    = StrUpCase(fpi_instr) + ' Maximum Counts'
		oNmax['LABEL']         = nSC EQ 4 ? sc[i] : 'N'
		oNmax['TITLE']         = 'N-max'
		oNmax['UNITS']         = ''
		oNmax['SI_CONVERSION'] = '>'
		IF nSC EQ 4 THEN oNmax['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Total Counts /////////////////////////////
	;-------------------------------------------
		oNtot = MrScalarTS( oF_dist['TIMEVAR'], Total(counts, 2), $
		                    /CACHE, $
		                    NAME = cts_tot_vnames[i] )
		
		;Set Attributes
		oNtot['CATDESC']       = 'Total number of counts in the distribution function.'
		oNtot['PLOT_TITLE']    = StrUpCase(fpi_instr) + ' Total Counts'
		oNtot['LABEL']         = nSC EQ 4 ? sc[i] : 'N'
		oNtot['TITLE']         = 'N-tot'
		oNtot['UNITS']         = ''
		oNtot['SI_CONVERSION'] = '>'
		IF nSC EQ 4 THEN oNtot['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Boltzmann Maximum Counts /////////////////
	;-------------------------------------------
		ctemp = Float(counts_max)
		IF nZero GT 0 THEN ctemp[iZero] = !Values.F_NaN
		nmax = Fix(Max(Temporary(ctemp), DIMENSION=2, NAN=(nZero GT 0)), TYPE=3)
		oNmax = MrScalarTS( oF_dist['TIMEVAR'], nmax, $
		                    /CACHE, $
		                    NAME = mcts_max_vnames[i], $
		                    /NO_COPY )
		
		;Set Attributes
		oNmax['CATDESC']       = 'Maximum number of counts in any single bin of the Maxwellian distribution function.'
		oNmax['PLOT_TITLE']    = StrUpCase(fpi_instr) + ' Maximum Counts'
		oNmax['LABEL']         = nSC EQ 4 ? sc[i] : 'N$\downB$'
		oNmax['TITLE']         = 'N-max'
		oNmax['UNITS']         = ''
		oNmax['SI_CONVERSION'] = '>'
		IF nSC EQ 4 THEN oNmax['COLOR'] = sc_colors[i]
	
	;-------------------------------------------
	; Boltzmann Total Counts ///////////////////
	;-------------------------------------------
		oNtot = MrScalarTS( oF_dist['TIMEVAR'], Total(counts_max, 2), $
		                    /CACHE, $
		                    NAME = mcts_tot_vnames[i] )
		
		;Set Attributes
		oNtot['CATDESC']       = 'Total number of counts in the Maxwellian distribution function.'
		oNtot['PLOT_TITLE']    = StrUpCase(fpi_instr) + ' Total Counts'
		oNtot['LABEL']         = nSC EQ 4 ? sc[i] : 'N$\downB$'
		oNtot['TITLE']         = 'N-tot'
		oNtot['UNITS']         = ''
		oNtot['SI_CONVERSION'] = '>'
		IF nSC EQ 4 THEN oNtot['COLOR'] = sc_colors[i]
	ENDFOR

;-------------------------------------------
; 4-SC Gradients ///////////////////////////
;-------------------------------------------
	IF nSC EQ 4 THEN BEGIN
		oPos  = ObjArr(nSC)
		oN    = ObjArr(nSC)
		oNf   = ObjArr(nSC)
		oS    = ObjArr(nSC)
		oSf   = ObjArr(nSC)
		oSn   = ObjArr(nSC)
		oSfn  = ObjArr(nSC)
		oSb   = ObjArr(nSC)
		oSbf  = ObjArr(nSC)
		oSbn  = ObjArr(nSC)
		oSbfn = ObjArr(nSC)
		oMfc  = ObjArr(nSC)
		oM    = ObjArr(nSC)
		oMf   = ObjArr(nSC)
		oMn   = ObjArr(nSC)
		oMfn  = ObjArr(nSC)
		
		oSt    = ObjArr(nSC)
		oSft   = ObjArr(nSC)
		oSnt   = ObjArr(nSC)
		oSfnt  = ObjArr(nSC)
		oSbt   = ObjArr(nSC)
		oSbft  = ObjArr(nSC)
		oSbnt  = ObjArr(nSC)
		oSbfnt = ObjArr(nSC)
		oMfct  = ObjArr(nSC)
		oMt    = ObjArr(nSC)
		oMft   = ObjArr(nSC)
		oMnt   = ObjArr(nSC)
		oMfnt  = ObjArr(nSC)
		
		oRef  = MrVar_Get(s_vnames[0])
		oTref = oRef['TIMEVAR']
		
		FOR i = 0, nSC - 1 DO BEGIN
			;Get the variables
			oNN    = MrVar_Get(n_calc_vnames[i])
			oNNf   = MrVar_Get(nf_calc_vnames[i])
			oRR    = MrVar_Get(r_vnames[i])
			oSS    = MrVar_Get(s_vnames[i])
			oSSf   = MrVar_Get(sf_vnames[i])
			oSSn   = MrVar_Get(sn_vnames[i])
			oSSfn  = MrVar_Get(sfn_vnames[i])
			oSSb   = MrVar_Get(sb_vnames[i])
			oSSbf  = MrVar_Get(sbf_vnames[i])
			oSSbn  = MrVar_Get(sbn_vnames[i])
			oSSbfn = MrVar_Get(sbfn_vnames[i])
			oMMfc  = MrVar_Get(mfrac_vnames[i])
			oMM    = MrVar_Get(mbar_vnames[i])
			oMMf   = MrVar_Get(mfbar_vnames[i])
			oMMn   = MrVar_Get(mnbar_vnames[i])
			oMMfn  = MrVar_Get(mfnbar_vnames[i])
			
			oSSt    = MrVar_Get(s_tot_vnames[i])
			oSSft   = MrVar_Get(sf_tot_vnames[i])
			oSSnt   = MrVar_Get(sn_tot_vnames[i])
			oSSfnt  = MrVar_Get(sfn_tot_vnames[i])
			oSSbt   = MrVar_Get(sb_tot_vnames[i])
			oSSbft  = MrVar_Get(sbf_tot_vnames[i])
			oSSbnt  = MrVar_Get(sbn_tot_vnames[i])
			oSSbfnt = MrVar_Get(sbfn_tot_vnames[i])
			oMMfct  = MrVar_Get(mfrac_tot_vnames[i])
			oMMt    = MrVar_Get(mbar_tot_vnames[i])
			oMMft   = MrVar_Get(mfbar_tot_vnames[i])
			oMMnt   = MrVar_Get(mnbar_tot_vnames[i])
			oMMfnt  = MrVar_Get(mfnbar_tot_vnames[i])
			
			;Interpolate to common time
			oPos[i]  = oRR    -> Interpol(oTref)
			oN[i]    = oNN    -> Interpol(oTref)
			oNf[i]   = oNNf   -> Interpol(oTref)
			
			oS[i]    = oSS    -> Interpol(oTref)
			oSf[i]   = oSSf   -> Interpol(oTref)
			oSn[i]   = oSSn   -> Interpol(oTref)
			oSfn[i]  = oSSfn  -> Interpol(oTref)
			oSb[i]   = oSSb   -> Interpol(oTref)
			oSbf[i]  = oSSbf  -> Interpol(oTref)
			oSbn[i]  = oSSbn  -> Interpol(oTref)
			oSbfn[i] = oSSbfn -> Interpol(oTref)
			oM[i]    = oMM    -> Interpol(oTref)
			oMfc[i]  = oMMfc  -> Interpol(oTref)
			oMf[i]   = oMMf   -> Interpol(oTref)
			oMn[i]   = oMMn   -> Interpol(oTref)
			oMfn[i]  = oMMfn  -> Interpol(oTref)
			
			oSt[i]    = oSSt    -> Interpol(oTref)
			oSft[i]   = oSSft   -> Interpol(oTref)
			oSnt[i]   = oSSnt   -> Interpol(oTref)
			oSfnt[i]  = oSSfnt  -> Interpol(oTref)
			oSbt[i]   = oSSbt   -> Interpol(oTref)
			oSbft[i]  = oSSbft  -> Interpol(oTref)
			oSbnt[i]  = oSSbnt  -> Interpol(oTref)
			oSbfnt[i] = oSSbfnt -> Interpol(oTref)
			oMfct[i]  = oMMfct  -> Interpol(oTref)
			oMt[i]    = oMMt    -> Interpol(oTref)
			oMft[i]   = oMMft   -> Interpol(oTref)
			oMnt[i]   = oMMnt   -> Interpol(oTref)
			oMfnt[i]  = oMMfnt  -> Interpol(oTref)
		ENDFOR
		
		;Compute gradients
		;   - Convert 1/km to 1/m
		oRecip   = MrVar_RecipVec(oPos[0], oPos[1], oPos[2], oPos[3])
		oGradN   = 1e-3 * oRecip -> Gradient(oN[0],   oN[1],   oN[2],   oN[3])
		oGradS   = 1e-3 * oRecip -> Gradient(oS[0],   oS[1],   oS[2],   oS[3])
		oGradSn  = 1e-3 * oRecip -> Gradient(oSn[0],  oSn[1],  oSn[2],  oSn[3])
		oGradSb  = 1e-3 * oRecip -> Gradient(oSb[0],  oSb[1],  oSb[2],  oSb[3])
		oGradSbn = 1e-3 * oRecip -> Gradient(oSbn[0], oSbn[1], oSbn[2], oSbn[3])
		oGradMfc = 1e-3 * oRecip -> Gradient(oMfc[0], oMfc[1], oMfc[2], oMfc[3])
		oGradM   = 1e-3 * oRecip -> Gradient(oM[0],   oM[1],   oM[2],   oM[3])
		oGradMn  = 1e-3 * oRecip -> Gradient(oMn[0],  oMn[1],  oMn[2],  oMn[3])
		
		oGradS_tot   = 1e-3 * oRecip -> Gradient(oSt[0],   oSt[1],   oSt[2],   oSt[3])
		oGradSn_tot  = 1e-3 * oRecip -> Gradient(oSnt[0],  oSnt[1],  oSnt[2],  oSnt[3])
		oGradSb_tot  = 1e-3 * oRecip -> Gradient(oSbt[0],  oSbt[1],  oSbt[2],  oSbt[3])
		oGradSbn_tot = 1e-3 * oRecip -> Gradient(oSbnt[0], oSbnt[1], oSbnt[2], oSbnt[3])
		oGradMfc_tot = 1e-3 * oRecip -> Gradient(oMfct[0], oMfct[1], oMfct[2], oMfct[3])
		oGradM_tot   = 1e-3 * oRecip -> Gradient(oMt[0],   oMt[1],   oMt[2],   oMt[3])
		oGradMn_tot  = 1e-3 * oRecip -> Gradient(oMnt[0],  oMnt[1],  oMnt[2],  oMnt[3])
		
		oDivNf   = 1e-3 * oRecip -> Divergence(oNf[0],   oNf[1],   oNf[2],   oNf[3])
		oDivSf   = 1e-3 * oRecip -> Divergence(oSf[0],   oSf[1],   oSf[2],   oSf[3])
		oDivSfn  = 1e-3 * oRecip -> Divergence(oSfn[0],  oSfn[1],  oSfn[2],  oSfn[3])
		oDivSbf  = 1e-3 * oRecip -> Divergence(oSbf[0],  oSbf[1],  oSbf[2],  oSbf[3])
		oDivSbfn = 1e-3 * oRecip -> Divergence(oSbfn[0], oSbfn[1], oSbfn[2], oSbfn[3])
		oDivMf   = 1e-3 * oRecip -> Divergence(oMf[0],   oMf[1],   oMf[2],   oMf[3])
		oDivMfn  = 1e-3 * oRecip -> Divergence(oMfn[0],  oMfn[1],  oMfn[2],  oMfn[3])
		
		oDivSf_tot   = 1e-3 * oRecip -> Divergence(oSft[0],   oSft[1],   oSft[2],   oSft[3])
		oDivSfn_tot  = 1e-3 * oRecip -> Divergence(oSfnt[0],  oSfnt[1],  oSfnt[2],  oSfnt[3])
		oDivSbf_tot  = 1e-3 * oRecip -> Divergence(oSbft[0],  oSbft[1],  oSbft[2],  oSbft[3])
		oDivSbfn_tot = 1e-3 * oRecip -> Divergence(oSbfnt[0], oSbfnt[1], oSbfnt[2], oSbfnt[3])
		oDivMf_tot   = 1e-3 * oRecip -> Divergence(oMft[0],   oMft[1],   oMft[2],   oMft[3])
		oDivMfn_tot  = 1e-3 * oRecip -> Divergence(oMfnt[0],  oMfnt[1],  oMfnt[2],  oMfnt[3])
		
		;
		;Attributes
		;
		
		;N
		oGradN -> SetName, gradn_vname
		oGradN -> Cache
		oGradN['CATDESC']       = 'Gradient of the density.'
		oGradN['COLOR']         = comp_colors
		oGradN['LABEL']         = nabla + 'N$\down' + ['X', 'Y', 'Z'] + '$'
		oGradN['TITLE']         = nabla + 'N!C(1/cm$\up4$)'
		oGradN['UNITS']         = '1/cm^4'
		oGradN['SI_CONVERSION'] = '1e8>1/m^4'
		
		;
		; Gradients: Sterling Approx
		;
		
		;Grad(S)
		oGradS -> SetName, grads_vname
		oGradS -> Cache
		oGradS['CATDESC']       = 'Gradient of the entropy density computed via integration ' + $
		                          'of the full velocity space distribution function.'
		oGradS['COLOR']         = comp_colors
		oGradS['PLOT_TITLE']    = 'Gradient of the Entropy Density'
		oGradS['LABEL']         = nabla + 's$\down' + ['X', 'Y', 'Z'] + '$'
		oGradS['TITLE']         = nabla + 's!C(eV/K/m$\up4$)'
		oGradS['UNITS']         = 'eV/K/m^4 ln(s^3/m^6)'
		oGradS['SI_CONVERSION'] = '1.602e-19>J/K/m^4'
		
		;Grad(Sn)
		oGradSn -> SetName, gradsn_vname
		oGradSn -> Cache
		oGradSn['CATDESC']       = 'Gradient of the entropy per particle computed via ' + $
		                           'integration of the full velocity space distribution function.'
		oGradSn['COLOR']         = comp_colors
		oGradSn['PLOT_TITLE']    = 'Gradient of the Entropy per Particle'
		oGradSn['LABEL']         = nabla + '(s/N)$\down' + ['X', 'Y', 'Z'] + '$'
		oGradSn['TITLE']         = nabla + 's/N!C(eV/K/m)'
		oGradSn['UNITS']         = 'eV/K/m ln(s^3/m^6)'
		oGradSn['SI_CONVERSION'] = '1.602e-19>J/K/m ln(s^3/m^6)'
		
		;Grad(Sb)
		oGradSb -> SetName, gradsb_vname
		oGradSb -> Cache
		oGradSb['CATDESC']       = 'Gradient of the Maxwellian entropy density.'
		oGradSb['COLOR']         = comp_colors
		oGradSb['PLOT_TITLE']    = 'Gradient of the Maxwellian Entropy per Particle'
		oGradSb['LABEL']         = nabla + 's$\downB,' + ['X', 'Y', 'Z'] + '$'
		oGradSb['TITLE']         = nabla + 's$\downB$!C(eV/K/m$\up4$)'
		oGradSb['UNITS']         = 'eV/K/m^4 ln(J m^2/kg)'
		oGradSb['SI_CONVERSION'] = '1.602e-19>J/K/m^4 ln(J m^2/kg)'
		
		;Grad(Sbn)
		oGradSbn -> SetName, gradsbn_vname
		oGradSbn -> Cache
		oGradSbn['CATDESC']       = 'Gradient of the Maxwellian Entropy per particle.'
		oGradSbn['COLOR']         = comp_colors
		oGradSbn['PLOT_TITLE']    = 'Gradient of the Maxwellian Entropy per Particle'
		oGradSbn['LABEL']         = nabla + 's$\downB,N,' + ['X', 'Y', 'Z'] + '$'
		oGradSbn['TITLE']         = nabla + 's$\downB$/N!C(eV/K/m)'
		oGradSbn['UNITS']         = 'eV/K/m ln(J m^2/kg)'
		oGradSbn['SI_CONVERSION'] = '1.602e-19>J/K/m ln(J m^2/kg)'
		
		;Grad(Mfrac)
		oGradMfrac -> SetName, gradmfrac_vname
		oGradMfrac -> Cache
		oGradMfrac['CATDESC']       = 'Gradient of the fractional non-Maxwellian entropy density.'
		oGradMfrac['COLOR']         = comp_colors
		oGradMfrac['LABEL']         = nabla + 'M$\down' + ['X', 'Y', 'Z'] + '$'
		oGradMfrac['TITLE']         = nabla + 'M'
		oGradMfrac['UNITS']         = ''
		oGradMfrac['SI_CONVERSION'] = '>'
		
		;Grad(M)
		oGradM -> SetName, gradm_vname
		oGradM -> Cache
		oGradM['CATDESC']       = 'Gradient of the non-Maxwellian entropy density.'
		oGradM['COLOR']         = comp_colors
		oGradM['LABEL']         = nabla + 'M$\down' + ['X', 'Y', 'Z'] + '$'
		oGradM['TITLE']         = nabla + 'M!C(eV/K/m$\up4$)'
		oGradM['UNITS']         = 'eV/K/m^4 [ln(s^3/m^6) - ln(J m^2/kg)]'
		oGradM['SI_CONVERSION'] = '1.602e-19>J/K/m^4 [ln(s^3/m^6) - ln(J m^2/kg)]'
		
		;Grad(Mn)
		oGradMn -> SetName, gradmn_vname
		oGradMn -> Cache
		oGradMn['CATDESC']       = 'Gradient of the non-Maxwellian entropy per particle.'
		oGradMn['COLOR']         = comp_colors
		oGradMn['LABEL']         = nabla + 'M$\down' + ['X', 'Y', 'Z'] + '$'
		oGradMn['TITLE']         = nabla + 'M/N!C(eV/K/m)'
		oGradMn['UNITS']         = 'eV/K/m [ln(s^3/m^6) - ln(J m^2/kg)]'
		oGradMn['SI_CONVERSION'] = '1.602e-19>J/K/m [ln(s^3/m^6) - ln(J m^2/kg)]'
		
		;
		; Gradients: Boltzmann
		;
		
		;Grad(St)
		oGradS_tot -> SetName, grads_tot_vname
		oGradS_tot -> Cache
		oGradS_tot['CATDESC']       = 'Gradient of the Boltzmann entropy density computed via integration ' + $
		                              'of the full velocity space distribution function.'
		oGradS_tot['COLOR']         = comp_colors
		oGradS_tot['PLOT_TITLE']    = 'Gradient of the Boltzmann Entropy Density'
		oGradS_tot['LABEL']         = nabla + 's$\downT' + ['X', 'Y', 'Z'] + '$'
		oGradS_tot['TITLE']         = nabla + 's$\downT$!C(eV/K/m$\up4$)'
		oGradS_tot['UNITS']         = 'eV/K/m^4'
		oGradS_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^4'
		
		;Grad(Snt)
		oGradSn_tot -> SetName, gradsn_tot_vname
		oGradSn_tot -> Cache
		oGradSn_tot['CATDESC']       = 'Gradient of the Boltzmann entropy per particle computed via ' + $
		                               'integration of the full velocity space distribution function.'
		oGradSn_tot['COLOR']         = comp_colors
		oGradSn_tot['PLOT_TITLE']    = 'Gradient of the Boltzmann Entropy per Particle'
		oGradSn_tot['LABEL']         = nabla + '(s/N)$\downT' + ['X', 'Y', 'Z'] + '$'
		oGradSn_tot['TITLE']         = nabla + 's$\downT$/N!C(eV/K/m)'
		oGradSn_tot['UNITS']         = 'eV/K/m'
		oGradSn_tot['SI_CONVERSION'] = '1.602e-19>J/K/m'
		
		;Grad(Sbt)
		oGradSb_tot -> SetName, gradsb_tot_vname
		oGradSb_tot -> Cache
		oGradSb_tot['CATDESC']       = 'Gradient of the Boltzmann entropy density of a Maxwellian.'
		oGradSb_tot['COLOR']         = comp_colors
		oGradSb_tot['PLOT_TITLE']    = 'Gradient of the Boltzmann Entropy per Particle of a Maxwellian'
		oGradSb_tot['LABEL']         = nabla + 's$\downBT,' + ['X', 'Y', 'Z'] + '$'
		oGradSb_tot['TITLE']         = nabla + 's$\downBT$!C(eV/K/m$\up4$)'
		oGradSb_tot['UNITS']         = 'eV/K/m^4'
		oGradSb_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^4'
		
		;Grad(Sbnt)
		oGradSbn_tot -> SetName, gradsbn_tot_vname
		oGradSbn_tot -> Cache
		oGradSbn_tot['CATDESC']       = 'Gradient of the Boltzmann Entropy per particle of a Maxwellian.'
		oGradSbn_tot['COLOR']         = comp_colors
		oGradSbn_tot['PLOT_TITLE']    = 'Gradient of the Boltzmann Entropy per Particle of a Maxwellian'
		oGradSbn_tot['LABEL']         = nabla + 's$\downB,N,T,' + ['X', 'Y', 'Z'] + '$'
		oGradSbn_tot['TITLE']         = nabla + 's$\downBT$/N!C(eV/K/m)'
		oGradSbn_tot['UNITS']         = 'eV/K/m'
		oGradSbn_tot['SI_CONVERSION'] = '1.602e-19>J/K/m'
		
		;Grad(Mfrac)
		oGradMfrac_tot -> SetName, gradmfrac_tot_vname
		oGradM_fractot -> Cache
		oGradMfrac_tot['CATDESC']       = 'Gradient of the fractional non-Maxwellian entropy density.'
		oGradMfrac_tot['COLOR']         = comp_colors
		oGradMfrac_tot['LABEL']         = nabla + 'M$\downT' + ['X', 'Y', 'Z'] + '$'
		oGradMfrac_tot['TITLE']         = nabla + 'M$\downT$'
		oGradMfrac_tot['UNITS']         = ''
		oGradMfrac_tot['SI_CONVERSION'] = '>'
		
		;Grad(Mt)
		oGradM_tot -> SetName, gradm_tot_vname
		oGradM_tot -> Cache
		oGradM_tot['CATDESC']       = 'Gradient of the non-Maxwellian entropy density.'
		oGradM_tot['COLOR']         = comp_colors
		oGradM_tot['LABEL']         = nabla + 'M$\downT' + ['X', 'Y', 'Z'] + '$'
		oGradM_tot['TITLE']         = nabla + 'M$\downT$!C(eV/K/m$\up4$)'
		oGradM_tot['UNITS']         = 'eV/K/m^4'
		oGradM_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^4'
		
		;Grad(Mnt)
		oGradMn_tot -> SetName, gradmn_tot_vname
		oGradMn_tot -> Cache
		oGradMn_tot['CATDESC']       = 'Gradient of the non-Maxwellian entropy per particle.'
		oGradMn_tot['COLOR']         = comp_colors
		oGradMn_tot['LABEL']         = nabla + 'M$\downT' + ['X', 'Y', 'Z'] + '$'
		oGradMn_tot['TITLE']         = nabla + 'M$\downT$/N!C(eV/K/m)'
		oGradMn_tot['UNITS']         = 'eV/K/m'
		oGradMn_tot['SI_CONVERSION'] = '1.602e-19>J/K/m'
		
		;
		; Divergence: Sterling Approx
		;
		
		;Div(Nf)
		oDivNf -> SetName, divnf_vname
		oDivNf -> Cache
		oDivNf['CATDESC']       = 'Divergence of the density flux.'
		oDivNf['LABEL']         = nabla + '.Nf'
		oDivNf['TITLE']         = nabla + '.Nf!C(1/cm$\up3/s$)'
		oDivNf['UNITS']         = '1/cm^3/s'
		oDivNf['SI_CONVERSION'] = '1e6>1/m^3/s'
		
		;Div(Sf)
		oDivSf -> SetName, divsf_vname
		oDivSf -> Cache
		oDivSf['CATDESC']       = 'Divergence of the entropy density flux computed via integration ' + $
		                          'of the full velocity space distribution function.'
		oDivSf['PLOT_TITLE']    = 'Divergence of the Entropy Density Flux'
		oDivSf['LABEL']         = nabla + '.sf'
		oDivSf['TITLE']         = nabla + '.sf!C(eV/K/m$\up3$/s)'
		oDivSf['UNITS']         = 'eV/K/m^3/s ln(s^3/m^6)'
		oDivSf['SI_CONVERSION'] = '1.602e-19>J/K/m^3/s'
		
		;Div(Sfn)
		oDivSfn -> SetName, divsfn_vname
		oDivSfn -> Cache
		oDivSfn['CATDESC']       = 'Divergence of the entropy flux per particle computed via ' + $
		                           'integration of the full velocity space distribution function.'
		oDivSfn['PLOT_TITLE']    = 'Divergence of the Entropy Flux per Particle'
		oDivSfn['LABEL']         = nabla + '.sf/N'
		oDivSfn['TITLE']         = nabla + '.sf/N!C(eV/K/s)'
		oDivSfn['UNITS']         = 'eV/K/s ln(s^3/m^6)'
		oDivSfn['SI_CONVERSION'] = '1.602e-19>J/K/s ln(s^3/m^6)'
		
		;Div(Sbf)
		oDivSbf -> SetName, divsbf_vname
		oDivSbf -> Cache
		oDivSbf['CATDESC']       = 'Divergence of the Maxwellian entropy density flux.'
		oDivSbf['PLOT_TITLE']    = 'Divergence of the Maxwellian Entropy Density Flux'
		oDivSbf['LABEL']         = nabla + '.s$\downBf$'
		oDivSbf['TITLE']         = nabla + '.s$\downB,f$!C(eV/K/m$\up3$/s)'
		oDivSbf['UNITS']         = 'eV/K/m^3/s ln(J m^2/kg)'
		oDivSbf['SI_CONVERSION'] = '1.602e-19>J/K/m^3/s ln(J m^2/kg)'
		
		;Div(Sbfn)
		oDivSbfn -> SetName, divsbfn_vname
		oDivSbfn -> Cache
		oDivSbfn['CATDESC']       = 'Divergence of the Maxwellian Entropy Flux per particle.'
		oDivSbfn['PLOT_TITLE']    = 'Gradient of the Maxwellian Entropy per Particle'
		oDivSbfn['LABEL']         = nabla + 's$\downBf$/N'
		oDivSbfn['TITLE']         = nabla + 's$\downB,f$/N!C(eV/K/s)'
		oDivSbfn['UNITS']         = 'eV/K/s ln(J m^2/kg)'
		oDivSbfn['SI_CONVERSION'] = '1.602e-19>J/K/s ln(J m^2/kg)'
		
		;Div(Mf)
		oDivMf -> SetName, divmf_vname
		oDivMf -> Cache
		oDivMf['CATDESC']       = 'Divergence of the non-Maxwellian entropy density flux.'
		oDivMf['LABEL']         = nabla + '.Mf'
		oDivMf['TITLE']         = nabla + '.Mf!C(eV/K/m$\up3$/s)'
		oDivMf['UNITS']         = 'eV/K/m^3/s [ln(s^3/m^6) - ln(J m^2/kg)]'
		oDivMf['SI_CONVERSION'] = '1.602e-19>J/K/m^3/s [ln(s^3/m^6) - ln(J m^2/kg)]'
		
		;Div(Mfn)
		oDivMfn -> SetName, divmfn_vname
		oDivMfn -> Cache
		oDivMfn['CATDESC']       = 'Divergence of the non-Maxwellian entropy flux per particle.'
		oDivMfn['LABEL']         = nabla + '.M$\downf$/N'
		oDivMfn['TITLE']         = nabla + '.Mf/N!C(eV/K/s)'
		oDivMfn['UNITS']         = 'eV/K/s [ln(s^3/m^6) - ln(J m^2/kg)]'
		oDivMfn['SI_CONVERSION'] = '1.602e-19>J/K/s [ln(s^3/m^6) - ln(J m^2/kg)]'
		
		
		;
		; Divergence: Boltzmann
		;
		
		;Div(Nf)
		oDivNf_tot -> SetName, divnf_tot_vname
		oDivNf_tot -> Cache
		oDivNf_tot['CATDESC']       = 'Divergence of the density flux.'
		oDivNf_tot['LABEL']         = nabla + '.Nf'
		oDivNf_tot['TITLE']         = nabla + '.Nf!C(1/cm$\up3/s$)'
		oDivNf_tot['UNITS']         = '1/cm^3/s'
		oDivNf_tot['SI_CONVERSION'] = '1e6>1/m^3/s'
		
		;Div(Sf)
		oDivSf_tot -> SetName, divsf_tot_vname
		oDivSf_tot -> Cache
		oDivSf_tot['CATDESC']       = 'Divergence of the Boltzmann entropy density flux computed via integration ' + $
		                              'of the full velocity space distribution function.'
		oDivSf_tot['PLOT_TITLE']    = 'Divergence of the Boltzmann Entropy Density Flux'
		oDivSf_tot['LABEL']         = nabla + '.s$\downTf$'
		oDivSf_tot['TITLE']         = nabla + '.s$\downTf$!C(eV/K/m$\up3$/s)'
		oDivSf_tot['UNITS']         = 'eV/K/m^3/s'
		oDivSf_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^3/s'
		
		;Div(Sfn)
		oDivSfn_tot -> SetName, divsfn_tot_vname
		oDivSfn_tot -> Cache
		oDivSfn_tot['CATDESC']       = 'Divergence of the Boltzmann entropy flux per particle computed via ' + $
		                               'integration of the full velocity space distribution function.'
		oDivSfn_tot['PLOT_TITLE']    = 'Divergence of the Boltzmann Entropy Flux per Particle'
		oDivSfn_tot['LABEL']         = nabla + '.s$\downTf$/N'
		oDivSfn_tot['TITLE']         = nabla + '.s$\downTf$/N!C(eV/K/s)'
		oDivSfn_tot['UNITS']         = 'eV/K/s'
		oDivSfn_tot['SI_CONVERSION'] = '1.602e-19>J/K/s'
		
		;Div(Sbf)
		oDivSbf_tot -> SetName, divsbf_tot_vname
		oDivSbf_tot -> Cache
		oDivSbf_tot['CATDESC']       = 'Divergence of the Boltzmann entropy density flux of a Maxwellian.'
		oDivSbf_tot['PLOT_TITLE']    = 'Divergence of the Boltzmann Entropy Density Flux of a Maxwellian'
		oDivSbf_tot['LABEL']         = nabla + '.s$\downBTf$'
		oDivSbf_tot['TITLE']         = nabla + '.s$\downBTf$!C(eV/K/m$\up3$/s)'
		oDivSbf_tot['UNITS']         = 'eV/K/m^3/s'
		oDivSbf_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^3/s'
		
		;Div(Sbfn)
		oDivSbfn_tot -> SetName, divsbfn_tot_vname
		oDivSbfn_tot -> Cache
		oDivSbfn_tot['CATDESC']       = 'Divergence of the Boltzmann Entropy Flux per particle of a Maxwellian.'
		oDivSbfn_tot['PLOT_TITLE']    = 'Gradient of the Boltzmann Entropy per Particle of a Maxwellian'
		oDivSbfn_tot['LABEL']         = nabla + 's$\downBTf$/N'
		oDivSbfn_tot['TITLE']         = nabla + 's$\downBTf$/N!C(eV/K/s)'
		oDivSbfn_tot['UNITS']         = 'eV/K/s'
		oDivSbfn_tot['SI_CONVERSION'] = '1.602e-19>J/K/s'
		
		;Div(Mf)
		oDivMf_tot -> SetName, divmf_tot_vname
		oDivMf_tot -> Cache
		oDivMf_tot['CATDESC']       = 'Divergence of the non-Maxwellian entropy density flux.'
		oDivMf_tot['LABEL']         = nabla + '.M$\downTf$'
		oDivMf_tot['TITLE']         = nabla + '.M$\downTf$!C(eV/K/m$\up3$/s)'
		oDivMf_tot['UNITS']         = 'eV/K/m^3/s'
		oDivMf_tot['SI_CONVERSION'] = '1.602e-19>J/K/m^3/s'
		
		;Div(Mfn)
		oDivMfn_tot -> SetName, divmfn_tot_vname
		oDivMfn_tot -> Cache
		oDivMfn_tot['CATDESC']       = 'Divergence of the non-Maxwellian entropy flux per particle.'
		oDivMfn_tot['LABEL']         = nabla + '.M$\downTf$/N'
		oDivMfn_tot['TITLE']         = nabla + '.M$\downTf$/N!C(eV/K/s)'
		oDivMfn_tot['UNITS']         = 'eV/K/s'
		oDivMfn_tot['SI_CONVERSION'] = '1.602e-19>J/K/s'
	ENDIF

;-------------------------------------------
; Finish Up ////////////////////////////////
;-------------------------------------------
	MrVar_Delete, [r_vnames, f_vnames, f_raw_vnames, df_raw_vnames, scpot_vnames]
END