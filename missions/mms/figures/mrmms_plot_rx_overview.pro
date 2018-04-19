; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_Rx_Overview
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
;       1. B
;       2. E
;       3. ni & ne
;       4. Vi
;       5. E spectra (e)
;       6. E spectra (i)
;
; :Categories:
;   MMS, EDI, MrVariable
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;
; :Keywords:
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
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
;       2016/09/21  -   Written by Matthew Argall
;-
function MrMMS_Plot_Rx_Overview, sc, mode, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
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
	if n_elements(trange) gt 0 then MrVar_SetTRange, trange

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	level    = 'l2'
	edp_mode = mode eq 'srvy' ? 'fast' : mode
	fgm_mode = mode
	fpi_mode = mode eq 'brst' ? mode : 'fast'
	
	coords     = 'gse'
	fpi_coords = 'dbcs'
	
	;Source names
	b_vname       = sc + '_fgm_b_'             +     coords + '_' + fgm_mode + '_' + level
	bvec_vname    = sc + '_fgm_bvec_'          +     coords + '_' + fgm_mode + '_' + level
	E_vname       = sc + '_edp_dce_'           +     coords + '_' + edp_mode + '_' + level
	ni_vname      = sc + '_dis_numberdensity_'                    + fpi_mode
	ne_vname      = sc + '_des_numberdensity_'                    + fpi_mode
	vi_vname      = sc + '_dis_bulkv_'         + fpi_coords + '_' + fpi_mode
	ve_vname      = sc + '_des_bulkv_'         + fpi_coords + '_' + fpi_mode
	ti_para_vname = sc + '_dis_temppara_'                         + fpi_mode
	te_para_vname = sc + '_des_temppara_'                         + fpi_mode
	ti_perp_vname = sc + '_dis_tempperp_'                         + fpi_mode
	te_perp_vname = sc + '_des_tempperp_'                         + fpi_mode
	Ee_vname      = sc + '_des_energyspectr_'                     + fpi_mode
	Ei_vname      = sc + '_dis_energyspectr_'                     + fpi_mode
	
	;Derived names
	j_vname = sc + '_fpi_currentdensity_' + fpi_coords + '_' + fpi_mode
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	
	;FGM
	if tf_load then begin
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     VARFORMAT = '*b_' + coords + '*'
		
		;EDP E-Field
		MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
		                 OPTDESC   = 'dce'
		                 VARFORMAT = '*dce_' + coords + '*'
	
		;DES
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = 'des-moms', $
		                     VARFORMAT = ['*energyspectr*', $
		                                  '*density*', $
		                                  '*bulk?_' + fpi_coords + '*', $
		                                  '*tempp*']
	
		;DIS
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = 'dis-moms', $
		                     VARFORMAT = ['*energyspectr*', $
		                                  '*density*', $
		                                  '*bulk?_'   + fpi_coords + '*', $
		                                  '*tempp*']
	endif
	
;-------------------------------------------
; Current Density //////////////////////////
;-------------------------------------------
	;Grab data
	n_e = MrVar_Get(ne_vname)
	vi  = MrVar_Get(vi_vname)
	ve  = MrVar_Get(ve_vname)

	;Interpolate ions
	vi_des = vi -> Interpol(ve)
	
	;Compute current density
	j = 1e-3 * n_e * (vi_des - ve)
	j -> SetName, j_vname
	j -> Cache
	
	;Free data
	obj_destroy, vi_des
	
;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	title = strupcase(sc)
	
	;B
	oB = MrVar_Get(b_vname)
	oB['LABEL'] = ['Bx', 'By', 'Bz', '|B|']
	oB['PLOT_TITLE'] = title
	
	;e- Energy Spectrogram
	oEe = MrVar_Get(Ee_vname)
	oEnergy = oEe['DEPEND_1']
	oEnergy['TITLE'] = 'E e-!C(eV)'
	
	;Ion Energy Spectrogram
	oEi = MrVar_Get(Ei_vname)
	oEnergy = oEI['DEPEND_1']
	oEnergy['TITLE'] = 'E ion!C(eV)'
	
	;NI
	ni     = MrVar_Get(ni_vname)
	nrange = [ Min([ni.min, n_e.min]), Max([ni.max, n_e.max]) ]
	ni['AXIS_RANGE'] = nrange
	ni['COLOR']      = 'Blue'
	ni['LABEL']      = 'ni'
	ni['TITLE']      = 'N!Ccm$\up-3$'
	
	;NE
	n_e['AXIS_RANGE'] = nrange
	n_e['COLOR']      = 'Red'
	n_e['LABEL']      = 'ne'
	
	;VI
	vi['TITLE'] = 'Vi!C(km/s)'
	vi['LABEL'] = ['Vx', 'Vy', 'Vz']
	
	;VE
	ve['LABEL'] = ['Vx', 'Vy', 'Vz']
	ve['TITLE'] = 'Ve!C(km/s)'
	
	;J
	j['COLOR'] = ['Blue', 'Forest Green', 'Red']
	j['LABEL'] = ['Jx', 'Jy', 'Jz']
	j['TITLE'] = 'J!C$\mu$A/m$\up2$'
	
	;TI_PARA
	ti_perp = MrVar_Get(ti_perp_vname)
	ti_para = MrVar_Get(ti_para_vname)
	ti_para['TITLE'] = 'Ti!C(eV)'
	ti_para['COLOR'] = 'Blue'
	ti_para['LABEL'] = 'Para'
	
	;TI_PERP
	ti_para['AXIS_RANGE'] = [0, max( [ti_para.max, ti_perp.max] ) * 1.1]
	ti_perp['TITLE']      = 'Ti!C(eV)'
	ti_perp['COLOR']      = 'Red'
	ti_perp['LABEL']      = 'Perp'
	
	;TE_PARA
	te_perp = MrVar_Get(te_perp_vname)
	te_para = MrVar_Get(te_para_vname)
	te_para['AXIS_RANGE'] = [0, max( [te_para.max, te_perp.max] ) * 1.1]
	te_para['COLOR']      = 'Blue'
	te_para['LABEL']      = 'Para'
	te_para['TITLE']      = 'Te!C(eV)'
	
	;TE_PERP
	te_perp['TITLE'] = 'Te!C(eV)'
	te_perp['COLOR'] = 'Red'
	te_perp['LABEL'] = 'Perp'
	
	;E
	E = MrVar_Get(e_vname)
	E['TITLE'] = 'E!C(mV/m)'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	tf_save = N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0
	
	;Plot
	win = MrVar_PlotTS( [B_vname, Ee_vname, Ei_vname, ni_vname, vi_vname, ve_vname, j_vname, ti_para_vname, te_para_vname, E_vname], $
;	                    BUFFER = tf_save, $
	                    /NO_REFRESH, $
	                    XSIZE  = 800, $
	                    YSIZE  = 800)
	
	;Overplot
	win = MrVar_OPlotTS( [ni_vname, ti_para_vname, te_para_vname], $
	                     [ne_vname, ti_perp_vname, te_perp_vname] )
	
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[12,10]

	win -> Refresh

;-------------------------------------------
; Save Figure //////////////////////////////
;-------------------------------------------
	IF tf_save THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
		ENDIF ELSE IF ~File_Test(output_dir, /DIRECTORY) THEN BEGIN
			MrPrintF, 'LogText', 'Creating directory: "' + output_dir + '".'
			File_MKDir, output_dir
		ENDIF
		
		;File name
		fname = StrJoin( [sc, 'fpi', mode, level, 'rx-overview'], '_' )
		fname = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Finished! ////////////////////////////////
;-------------------------------------------
	return, win
end