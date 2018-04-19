; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_Whistler_Omni
;
;*****************************************************************************************
;   Copyright (c) 2016, Matthew Argall                                                   ;
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
; PURPOSE:
;+
;   Generate a plot to provide an overview of reconnection quantities:
;       1. B & |B|
;       2. Ve
;       3. PAD + loss-cone
;       4. Curvature scattering parameters
;       5. Ee Spectra + E Resonant
;       6. B PSD Omni
;       7. E PSD Omni
;
;   Figure 3 of:
;       Lavraud, B. et al. (2016), Currents and associated electron scattering and
;           bouncing near the diffusion region at Earth’s magnetopause, Geophys. Res.
;           Lett., 3042–3050, doi:10.1002/2016GL068359.
;
; :Categories:
;   MMS
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;
; :Keywords:
;       TRANGE:     in, optional, type=string/strarr(2), default=MrVar_GetTRange()
;                   The start and end times of the data interval to be plotted, formatted
;                       as 'YYYY-MM-DDThh:mm:ss'
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
;       2016/12/08  -   Written by Matthew Argall
;-
function MrMMS_Plot_Whistler_Omni, sc, mode, species, $
ENERGIES=energies, $
NO_LOAD=no_load, $
TRANGE=trange
	compile_opt idl2
	
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if n_elements(win) gt 0 then obj_destroy, win
		MrPrintF, 'LogErr'
		return, !Null
	endif
	
	tf_load = ~keyword_set(no_load)
	if n_elements(nfft)      eq 0 then nfft      = 1024
	if n_elements(fgm_instr) eq 0 then fgm_instr = 'fgm'
	if n_elements(b_mirror)  eq 0 then b_mirror  = 42                     ;nT
	if n_elements(f_wave)    eq 0 then f_wave    = 100                    ;Hz
	if n_elements(energies)  eq 0 then energies  = [20, 100, 200, 500]    ;eV
	if n_elements(species)   eq 0 then species   = 'e'
	if n_elements(trange)    gt 0 then MrVar_SetTRange, trange
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	case fgm_instr of
		'afg': fgm_level = 'l2pre'
		'dfg': fgm_level = 'l2pre'
		'fgm': fgm_level = 'l2'
		else:  message, 'Invalid FGM instrument: "' + fgm_instr + '".'
	endcase
	fpi_instr = 'd' + species + 's'
	case mode of
		'slow': scm_optdesc = 'scs'
		'fast': scm_optdesc = 'scf'
		'srvy': scm_optdesc = 'scsrvy'
		'brst': scm_optdesc = 'scb'
		else: message, 'Invalid MODE: "' + mode + '".'
	endcase
	
	level    = 'l2'
	edp_mode = mode eq 'srvy' ? 'fast' : mode
	fgm_mode = mode
	dsp_mode = 'fast'
	fpi_mode = mode eq 'brst' ? mode : 'fast'
	coords   = 'gse'

	;Source names
	b_vname       = sc + '_' + fgm_instr + '_b_'    + coords + '_' + fgm_mode + '_' + fgm_level
	bvec_vname    = sc + '_' + fgm_instr + '_bvec_' + coords + '_' + fgm_mode + '_' + fgm_level
	bmag_vname    = sc + '_' + fgm_instr + '_bmag_' + coords + '_' + fgm_mode + '_' + fgm_level
	b1_vec_vname  = 'mms1_'  + fgm_instr + '_bvec_' + coords + '_' + fgm_mode + '_' + fgm_level
	b2_vec_vname  = 'mms2_'  + fgm_instr + '_bvec_' + coords + '_' + fgm_mode + '_' + fgm_level
	b3_vec_vname  = 'mms3_'  + fgm_instr + '_bvec_' + coords + '_' + fgm_mode + '_' + fgm_level
	b4_vec_vname  = 'mms4_'  + fgm_instr + '_bvec_' + coords + '_' + fgm_mode + '_' + fgm_level
	r1_vec_vname  = 'MMS1_DEFEPH_R'
	r2_vec_vname  = 'MMS2_DEFEPH_R'
	r3_vec_vname  = 'MMS3_DEFEPH_R'
	r4_vec_vname  = 'MMS4_DEFEPH_R'
	n_vname       = strjoin([sc, fpi_instr, 'numberdensity',     fpi_mode], '_')
	v_vname       = strjoin([sc, fpi_instr, 'bulkv', coords,     fpi_mode], '_')
	pad_vname     = strjoin([sc, fpi_instr, 'pitchangdist',      fpi_mode], '_')
	espec_vname   = strjoin([sc, fpi_instr, 'energyspectr_omni', fpi_mode], '_')
	acb_vname     = strjoin([sc, 'scm', 'acb', coords, scm_optdesc, mode, level], '_')
	dce_vname     = strjoin([sc, 'edp', 'dce', coords, mode, level], '_')
	
	;Derived names
	mir1_hi_vname = sc + '_' + fgm_instr + '_mirror1_low_'           + fgm_mode + '_' + fgm_level
	mir1_lo_vname = sc + '_' + fgm_instr + '_mirror1_high_'          + fgm_mode + '_' + fgm_level
	mir2_hi_vname = sc + '_' + fgm_instr + '_mirror2_low_'           + fgm_mode + '_' + fgm_level
	mir2_lo_vname = sc + '_' + fgm_instr + '_mirror2_high_'          + fgm_mode + '_' + fgm_level
	k1_vname      = sc + '_' + fgm_instr + '_k1_'                    + fgm_mode + '_' + fgm_level
	k2_vname      = sc + '_' + fgm_instr + '_k2_'                    + fgm_mode + '_' + fgm_level
	k3_vname      = sc + '_' + fgm_instr + '_k3_'                    + fgm_mode + '_' + fgm_level
	k4_vname      = sc + '_' + fgm_instr + '_k4_'                    + fgm_mode + '_' + fgm_level
	scm_psd_vname = strjoin([sc, 'scm', 'acb', 'psd', 'omni', mode, level], '_')
	edp_psd_vname = strjoin([sc, 'edp', 'dce', 'psd', 'omni', mode, level], '_')
	eres_vname    = strjoin([sc, 'resonant', 'energy', mode, level], '_')
	fce_vname     = strjoin([sc, 'fce', mode, level], '_')
	fpe_vname     = strjoin([sc, 'fpe', mode, level], '_')
	flh_vname     = strjoin([sc, 'flh', mode, level], '_')
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	if tf_load then begin
	;-------------------------------------------
	; FGM //////////////////////////////////////
	;-------------------------------------------
		MrMMS_FGM_Load_Data, '', mode, $
		                     INSTR     = fgm_instr, $
		                     VARFORMAT = ['*' + fgm_instr + '_b_' + coords + '*', $
		                                  '*' + fgm_instr + '_*_'   + coords]
	
	;-------------------------------------------
	; EPHEMERIS ////////////////////////////////
	;-------------------------------------------
		MrMMS_Load_Data, ANC_PRODUCT ='defeph', $
		                 VARFORMAT   = '*R'
	
	;-------------------------------------------
	; FPI //////////////////////////////////////
	;-------------------------------------------
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = 'des-moms', $
		                     VARFORMAT = ['*energyspectr_omni*', $
		                                  '*pitchangdist*', $
		                                  '*density*', $
		                                  '*bulk?_' + coords + '*']
	
	;-------------------------------------------
	; B PSD ////////////////////////////////////
	;-------------------------------------------
		;DSP for SRVY
		if mode eq 'srvy' then begin
			MrMMS_Load_Data, sc, 'dsp', dsp_mode, level, $
			                 OPTDESC   = 'bpsd', $
			                 VARFORMAT = '*bpsd_omni*'
		
		;SCM for BRST
		endif else begin
			MrMMS_Load_Data, sc, 'scm', mode, level, $
			                 OPTDESC   = 'scb', $
			                 VARFORMAT = '*_scb_brst_l2'
		endelse
	
	;-------------------------------------------
	; E PSD ////////////////////////////////////
	;-------------------------------------------
		;DSP for SRVY
		if mode eq 'srvy' then begin
			MrMMS_Load_Data, sc, 'dsp', dsp_mode, level, $
			                 OPTDESC   = 'epsd', $
			                 VARFORMAT = '*epsd_omni*'
		
		;SCM for BRST
		endif else begin
			MrMMS_Load_Data, sc, 'edp', mode, level, $
			                 OPTDESC   = 'dce', $
			                 VARFORMAT = '*_dce_' + coords + '_*'
		endelse
	endif

	if ~MrVar_IsCached(bvec_vname) then begin
		b_vname      = sc + '_' + fgm_instr + '_'     + fgm_mode + '_' + fgm_level + '_' + coords
		bvec_vname   = sc + '_' + fgm_instr + '_vec_' + fgm_mode + '_' + fgm_level + '_' + coords
		bmag_vname   = sc + '_' + fgm_instr + '_mag_' + fgm_mode + '_' + fgm_level + '_' + coords
		b1_vec_vname = 'mms1_'  + fgm_instr + '_vec_' + fgm_mode + '_' + fgm_level + '_' + coords
		b2_vec_vname = 'mms2_'  + fgm_instr + '_vec_' + fgm_mode + '_' + fgm_level + '_' + coords
		b3_vec_vname = 'mms3_'  + fgm_instr + '_vec_' + fgm_mode + '_' + fgm_level + '_' + coords
		b4_vec_vname = 'mms4_'  + fgm_instr + '_vec_' + fgm_mode + '_' + fgm_level + '_' + coords

		if ~MrVar_IsCached(bvec_vname) $
			then message, 'Problem with FGM variable names.'
	endif

;-------------------------------------------
; Mirror Angle /////////////////////////////
;-------------------------------------------
	oBmag = MrVar_Get(bmag_vname)
	
	;Mirror angle and its complement
	anglo = ASin( Sqrt(oBmag['DATA'] / b_mirror) ) * 180.0/!pi
	anghi = 180.0 - anglo
	
	;Create time series variables
	oAngLo = MrScalarTS(oBmag['TIMEVAR'], anglo, /CACHE, NAME=mir1_lo_vname)
	oAngHi = MrScalarTS(oBmag['TIMEVAR'], anghi, /CACHE, NAME=mir1_hi_vname)
	
;-------------------------------------------
; Scattering ///////////////////////////////
;-------------------------------------------
	;Particle mass
	case species of
		'e': mass = MrConstants('m_e')
		'i': mass = MrConstants('m_i')
		else: message, 'SPECIES must be {"e" | "i"}.'
	endcase
	
	;Create multi-spacecraft methods object
	oFGM = MrMMS_FGM_4sc( b1_vec_vname, b2_vec_vname, b3_vec_vname, b4_vec_vname, $
	                      r1_vec_vname, r2_vec_vname, r3_vec_vname, r4_vec_vname )
	
	;Adiabatic parameter
	oK1 = oFGM -> Kappa(energies[0], mass, NAME=k1_vname, /CACHE)
	oK2 = oFGM -> Kappa(energies[1], mass, NAME=k2_vname, /CACHE)
	oK3 = oFGM -> Kappa(energies[2], mass, NAME=k3_vname, /CACHE)
	oK4 = oFGM -> Kappa(energies[3], mass, NAME=k4_vname, /CACHE)
	
	;Cleanup
	obj_destroy, oFGM
	
;-------------------------------------------
; B PSD ////////////////////////////////////
;-------------------------------------------
	;Detrend data
	oB  = MrVar_Get(acb_vname)
	oE  = MrVar_Get(dce_vname)
	odB = oB -> Detrend(nfft)
	odE = oE -> Detrend(nfft)

	;B PSD Omni
	odB_psd = odB -> Spectrogram(nfft, nshift, NAME=scm_psd_vname, /CACHE, WINDOW='hanning')
	odB_psd = MrTimeSeries( odB_psd['TIMEVAR'], total(odB_psd['DATA'], 3) )
	
	;E PSD Omni
	odE_psd = odE -> Spectrogram(nfft, nshift, NAME=edp_psd_vname, /CACHE, WINDOW='hanning')
	odE_psd = MrTimeSeries( odE_psd['TIMEVAR'], total(odE_psd['DATA'], 3) )
	
	;Frequency lines
	f_ce = MrVar_Freq_Cyclotron(bmag_vname, 'm_e', NAME=fce_vname, /CACHE)
	f_lh = MrVar_Freq_LowerHybrid(bmag_vname, n_vname, NAME=flh_vname, /CACHE)
	f_pe = MrVar_Freq_Plasma(n_vname, 'm_e', NAME=fpe_vname, /CACHE)
	
	;Cleanup
	obj_destroy, [odB, odE]

;-----------------------------------------------------
; Calculate Resonant Energy \\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	mass  = species eq 'i' ? MrConstants('m_p') : MrConstants('m_e')
	oEr   = MrVar_Whistler_Eres(bmag_vname, n_vname, mass, f_wave, NAME=eres_vname, /CACHE)

;-------------------------------------------
; Set Properties ///////////////////////////
;-------------------------------------------
	title = strupcase(sc)

	;B
	oB = MrVar_Get(b_vname)
	oB['AXIS_RANGE'] = [oB.min*1.1, oB.max*1.1]
	oB['PLOT_TITLE'] = title
	oB['LABEL']      = ['Bx', 'By', 'Bz', '|B|']
	
	;V
	oV = MrVar_Get(v_vname)
	oV['TITLE'] = 'Ve!C(km/s)'
	oV['LABEL'] = ['Vx', 'Vy', 'Vz']
	
	;PAD
	oPAD = MrVar_Get(pad_vname)
	oPA  = oPAD['DEPEND_1']
	oPAD['TITLE'] = 'EFlux!C' + oPAD['UNITS']
	oPA['TITLE']  = 'PA!C(deg)'
	
	;K1
	oK1['AXIS_RANGE'] = [oK4.min < 1e-1, oK1.max*1.1]
	oK1['COLOR']      = 'Blue'
	oK1['LABEL']      = string(energies[0], FORMAT='(i0)') + 'eV'
	oK1['LOG']        = 1
	oK1['TITLE']      = 'K$\up2$'
	
	;K2
	oK2['COLOR'] = 'Green'
	oK2['LABEL'] = string(energies[1], FORMAT='(i0)') + 'eV'
	oK2['LOG']   = 1
	
	;K3
	oK3 = MrVar_Get( k3_vname )
	oK3['COLOR'] = 'Red'
	oK3['LABEL'] = string(energies[2], FORMAT='(i0)') + 'eV'
	oK3['LOG']   = 1
	
	;K4
	oK4 = MrVar_Get( k4_vname )
	oK4['COLOR'] = 'Magenta'
	oK4['LABEL'] = string(energies[3], FORMAT='(i0)') + 'eV'
	oK4['LOG']   = 1
	
	;ESpec
	oEe = MrVar_Get(espec_vname)
	oEn = oEe['DEPEND_1']
	oEn['TITLE'] = 'Energy!C(eV)'
	
	;B PSD
	odB_psd['TITLE'] = 'B Omni!C(nT$\up2$/Hz)'
	
	;E PSD
	odE_psd['TITLE'] = 'E Omni!C(mV/m)$\up2$/Hz'
	
	;F_LH
	f_lh['THICK'] = 2
	f_ce['THICK'] = 2
	oEr['THICK']  = 2

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Create plot
	win = MrVar_PlotTS( [b_vname, v_vname, pad_vname, k1_vname, espec_vname, scm_psd_vname, edp_psd_vname], $
	                    XSIZE = 800, $
	                    YSIZE = 700)
	win    -> Refresh, /DISABLE
	
	;Overplot loss cone and scattering parameters
	win = MrVar_OPlotTS( [pad_vname,  pad_vname], [mir1_lo_vname, mir1_hi_vname] )
	win = MrVar_OPlotTS( [k1_vname, k1_vname, k1_vname], $
	                     [k2_vname, k3_vname, k4_vname] )
	win = MrVar_OPlotTS( espec_vname, eres_vname )
	win = MrVar_OPlotTS( [scm_psd_vname, scm_psd_vname, scm_psd_vname, edp_psd_vname, edp_psd_vname, edp_psd_vname], $
	                     [fce_vname,     fpe_vname,     flh_vname,     fce_vname,     fpe_vname,     fpe_vname] )
;	if nMirror eq 2 then win = MrVar_OPlotTS( [pad_vname, pad_vname], [mir2_lo_vname, mir2_hi_vname] )
	
	;Draw horizontal lines
	trange = MrVar_GetTRange('SSM')
	oHorz  = MrPlotS( trange, [10, 10], LINESTYLE='--', NAME='k=10', TARGET=win[k1_vname])
	oHorz  = MrPlotS( trange, [25, 25], LINESTYLE='--', NAME='k=25', TARGET=win[k1_vname])
	
	win.oxmargin = [12,13]
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> Refresh
	return, win
end