; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_Wave_PowerSpectra
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
;   Plot quantities related to the Poynting flux:
;       1. Bxyz, |B|
;       2. Field
;       3. Detrended Field
;       4. PSD X
;       5. PSD Y
;       6. PSD Z
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
;       NFFT:       in, optional, type=integer, default=512
;                   Number of points to use per FFT interval.
;       NSHIFT:     in, optional, type=integer, default=`NFFT`/2
;                   Number of points to shift after each FFT.
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
;       OPTDESC:    in, optional, type=string
;                   SCM optional descriptor to use. Choices are: 'scs', 'scf', 'scsrvy',
;                       'scb', 'sch'. The default is to use 'scsrvy' when `MODE`='srvy'
;                       and 'scb' when `MODE`='brst'.
;       OUTPUT_DIR: in, optional, type=string, default='~/figures/'
;                   Directory in which to save the figure. If neither `OUTPUT_DIR` or
;                       `OUTPUT_EXT` are given, no file is made.
;       OUTPUT_EXT: in, optional, type=string/strarr, default='png'
;                   Extension (and file type) of the figure. Options include
;                       'png', 'jpeg', 'tiff', 'ps', 'eps', 'pdf'. If neither
;                       `OUTPUT_DIR` or `OUTPUT_EXT` are given, no file is made.
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
;       2016/11/27  -   Written by Matthew Argall
;       2017/02/01  -   Simplify plot to 1-column focusing of Poynting Flux - MRA
;       2017/06/06  -   Window the FFT. - MRA
;-
FUNCTION MrMMS_Plot_Wave_PowerSpectra, sc, instr, mode, nfft, nshift, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
NO_LOAD=no_load, $
OPTDESC=optdesc, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
TRANGE=trange
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF N_Elements(win) GT 0 THEN Obj_Destroy, win
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords)   EQ 0 THEN coords   = 'gse'
	IF N_Elements(level)    EQ 0 THEN level    = 'l2'
	IF N_Elements(nfft)     EQ 0 THEN nfft     = 512
	IF N_Elements(nshift)   EQ 0 THEN nshift   = nfft / 2
	IF N_Elements(ndetrend) EQ 0 THEN ndetrend = nfft
	IF N_Elements(trange)   GT 0 THEN MrVar_SetTRange, trange
	
	;Restrictions
	IF Array_Equal(instr EQ ['afg', 'dfg', 'fgm', 'scm', 'edp'], 0) $
		THEN Message, 'INSTR invalid: "' + instr + '". Must be {"afg" | "dfg" | "fgm" | "scm" | "edp"}.'
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;FSM
	IF instr EQ 'fsm' && mode NE 'brst' THEN Message, 'FSM data is only available in burst mode.'
	fsm_optdesc = '8khz'
	fsm_level   = 'l3'
	fsm_mode    = 'brst'
	fsm_coords  = coords
	
	;FGM parameters
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = (instr NE 'scm' && instr NE 'edp') ? instr : 'fgm'
	fgm_coords = (coords EQ 'dsl'  || coords EQ 'dbcs') ? 'dmpa' : coords
	fgm_mode   = (mode   EQ 'fast' || mode   EQ 'slow') ? 'srvy' : mode
	
	;SCM
	scm_coords = coords EQ 'gse' ? coords : 'gse'
	scm_mode   = (mode  EQ 'fast' || mode EQ 'slow') ? 'srvy' : mode
	IF coords NE 'gse' THEN MrPrintF, 'SCM data is only available in GSE.'
	IF n_elements(optdesc) EQ 0 THEN BEGIN
		CASE mode OF
			'slow': scm_optdesc = 'scs'
			'fast': scm_optdesc = 'scf'
			'srvy': scm_optdesc = 'scsrvy'
			'brst': scm_optdesc = 'scb'       ;also schb
			ELSE: Message, 'Invalid MODE: "' + mode + '".'
		ENDCASE
	ENDIF
	
	;EDP
	edp_coords = (coords EQ 'dmpa' || coords EQ 'dbcs') ? 'dsl' : coords
	edp_mode   = mode EQ 'srvy' ? 'fast' : mode
	IF mode EQ 'srvy' THEN MrPrintF, 'LogWarn', 'EDP does not have srvy data. Using fast.'
	
	;FPI
	fpi_coords = (coords EQ 'dmpa' || coords EQ 'dsl') ? 'dbcs' : coords
	fpi_mode   = mode EQ 'brst' ? mode : 'fast'

	;Source names
	b_vname    = StrJoin( [sc, fgm_instr, 'b',    fgm_coords,              fgm_mode, level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords,              fgm_mode, level], '_' )
	bmag_vname = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords,              fgm_mode, level], '_' )
	bfsm_vname = StrJoin( [sc, 'fsm',     'b',    fsm_coords,              fsm_mode, fsm_level], '_') 
	dcb_vname  = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords,              fgm_mode, level], '_' )
	acb_vname  = StrJoin( [sc, 'scm',     'acb',  scm_coords, scm_optdesc, scm_mode, level], '_' )
	dce_vname  = StrJoin( [sc, 'edp',     'dce',  edp_coords,              edp_mode, level], '_' )
	ni_vname   = StrJoin( [sc, 'dis',     'numberdensity',                 fpi_mode], '_' )
	ne_vname   = StrJoin( [sc, 'des',     'numberdensity',                 fpi_mode], '_' )
	
	CASE instr OF
		'edp': field_vname = dce_vname
		'scm': field_vname = acb_vname
		'fsm': field_vname = fsm_vname
		ELSE:  field_vname = dcb_vname
	ENDCASE
	
	;Output names
	dv_vname    = field_vname + '_delta'
	f_ceh_vname = StrJoin( [sc, 'plasma', 'fce', 'half', mode, level], '_' )
	f_ce_vname  = StrJoin( [sc, 'plasma', 'fce',          mode, level], '_' )
	f_pe_vname  = StrJoin( [sc, 'plasma', 'fpe',          mode, level], '_' )
	f_ci_vname  = StrJoin( [sc, 'plasma', 'fci',          mode, level], '_' )
	f_pi_vname  = StrJoin( [sc, 'plasma', 'fpi',          mode, level], '_' )
	f_lh_vname  = StrJoin( [sc, 'plasma', 'flh',          mode, level], '_' )
	psd_x_vname = field_vname + '_x_psd'
	psd_y_vname = field_vname + '_y_psd'
	psd_z_vname = field_vname + '_z_psd'
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;SCM
		IF instr EQ 'scm' THEN BEGIN
			MrMMS_Load_Data, sc, 'scm', mode, level, $
			                 OPTDESC = optdesc
	
		;EDP
		ENDIF ELSE IF instr EQ 'edp' THEN BEGIN
			MrMMS_Load_Data, sc, 'edp', edp_mode, level, $
			                 OPTDESC   = 'dce', $
			                 VARFORMAT = dce_vname
		
		;FSM
		ENDIF ELSE IF instr EQ 'fsm' THEN BEGIN
			MrMMS_Load_Data, sc, 'fsm', fsm_mode, fsm_level, $
			                 OPTDESC   = fsm_optdesc, $
			                 VARFORMAT = '*_b_'+fsm_coords+'_'+fsm_mode+'*'
		
		;FGM
		ENDIF ELSE IF instr NE fgm_instr THEN BEGIN
			MrMMS_FGM_Load_Data, sc, fgm_mode, $
			                     VARFORMAT = '*b_'+fgm_coords+'*'+level
		ENDIF
		
		;|B|
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     VARFORMAT = '*b_'+fgm_coords+'*'+level
		
		;DENSITY
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = ['des-moms', 'dis-moms'], $
		                     VARFORMAT = '*numberdensity_'+fpi_mode
	ENDIF
	
;-------------------------------------------
; Poynting Flux ////////////////////////////
;-------------------------------------------
	oVec = MrVar_Get(field_vname)
	
	;Detrend data
	odV = oVec -> Detrend(nfft, /CACHE, NAME=dv_vname)
	
	;Spectrogram
	oSpec   = odV -> Spectrogram( nfft, nshift, WINDOW='hanning' )
	oSpec_X = oSpec[*,*,0]
	oSpec_Y = oSpec[*,*,1]
	oSpec_Z = oSpec[*,*,2]
	
;-------------------------------------------
; Frequency Lines //////////////////////////
;-------------------------------------------
	oBmag = MrVar_Get(bmag_vname)
	oNi   = MrVar_Get(ni_vname)
	
	;"Upsample"
	oN_fgm = oNi -> Interpol(oBmag)

	;Calculate key frequencies
	oFce = MrVar_Freq_Cyclotron(oBmag, 'm_e', /CACHE, NAME=f_ce_vname)
	oFci = MrVar_Freq_Cyclotron(oBmag, 'm_p', /CACHE, NAME=f_ci_vname)
	oFpe = MrVar_Freq_Plasma(ne_vname, 'm_e', /CACHE, NAME=f_pe_vname)
	oFpi = MrVar_Freq_Plasma(oN_fgm,   'm_p', /CACHE, NAME=f_pi_vname)
	oFlh = MrVar_Freq_LowerHybrid(oBmag, ne_vname, /CACHE, NAME=f_lh_vname)
	
	;Half cyclotron
	oFce2  = oFce / 2.0
	oFce2 -> SetName, f_ceh_vname
	oFce2 -> Cache
	
	;Determine if the frequenies should be plotted
	f_vnames = []
	oFreq    = oSpec['DEPEND_1']
	f_range  = [ Min(oFreq['DATA',1:*], MAX=fmax), fmax ]
	IF ~(oFce.min  GT f_range[1] || oFce.max  LT f_range[0]) THEN f_vnames = [f_vnames, f_ce_vname]
	IF ~(oFce2.min GT f_range[1] || oFce2.max LT f_range[0]) THEN f_vnames = [f_vnames, f_ceh_vname]
	IF ~(oFci.min  GT f_range[1] || oFci.max  LT f_range[0]) THEN f_vnames = [f_vnames, f_ci_vname]
	IF ~(oFpe.min  GT f_range[1] || oFpe.max  LT f_range[0]) THEN f_vnames = [f_vnames, f_pe_vname]
	IF ~(oFpi.min  GT f_range[1] || oFpi.max  LT f_range[0]) THEN f_vnames = [f_vnames, f_pi_vname]
	IF ~(oFlh.min  GT f_range[1] || oFlh.max  LT f_range[0]) THEN f_vnames = [f_vnames, f_lh_vname]

;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	title = StrJoin( StrUpCase([sc, instr]), ' ' )
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = title
	
	;EDP
	IF instr EQ 'edp' THEN BEGIN
		oVec['TITLE'] = 'E!C(mV/m)'
		
		;DELTA
		odV['TITLE']  = 'dE!C(mV/m)'
		odV['LABEL']  = 'dE$\down' + ['X', 'Y', 'Z'] + '$'
	
	;SCM
	ENDIF ELSE IF instr EQ 'scm' THEN BEGIN
		oVec['TITLE'] = 'B!C(nT)'
		
		;DELTA
		odV['TITLE'] = 'dB!C(nT)'
		odV['LABEL'] = 'dB$\down' + ['X', 'Y', 'Z'] + '$'
	
	;FGM
	ENDIF ELSE BEGIN
		odV['TITLE'] = 'dB!C(nT)'
		odV['LABEL'] = 'dB$\down' + ['X', 'Y', 'Z'] + '$'
	ENDELSE
	
	;X PSD
	oSpec_X -> SetName, psd_x_vname
	oSpec_X -> Cache
	oSpec_X['TITLE'] = 'PSD!C(' + oVec['UNITS'] + ')^2/Hz'
	
	;Y PSD
	oSpec_Y -> SetName, psd_y_vname
	oSpec_Y -> Cache
	oSpec_Z['TITLE'] = 'PSD!C(' + oVec['UNITS'] + ')^2/Hz'
	
	;Z PSD
	oSpec_Z -> SetName, psd_z_vname
	oSpec_Z -> Cache
	oSpec_Z['TITLE'] = 'PSD!C(' + oVec['UNITS'] + ')^2/Hz'
	
	;FCE
	oFce['COLOR']     = 'Black'
	oFce['LABEL']     = 'fce'
	oFce['LINESTYLE'] = '-'
	
	;FCE/2
	oFce2['COLOR']     = 'Black'
	oFce2['LABEL']     = 'fce/2'
	oFce2['LINESTYLE'] = '--'
	
	;FCI
	oFci['COLOR']     = 'White'
	oFci['LABEL']     = 'fci'
	oFci['LINESTYLE'] = '--'
	
	;FPE
	oFpe['COLOR']     = 'Black'
	oFpe['LABEL']     = 'fpe'
	oFpe['LINESTYLE'] = '-:'
	
	;FPI
	oFpi['COLOR']     = 'White'
	oFpi['LABEL']     = 'fpi'
	oFpi['LINESTYLE'] = '-:'
	
	;FLH
	oFlh['COLOR'] = 'Red'
	oFlh['LABEL'] = 'flh'
	
;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Variables to plot
	IF field_vname NE bvec_vname $
		THEN vars = [b_vname, field_vname, dv_vname, psd_x_vname, psd_y_vname, psd_z_vname] $
		ELSE vars = [b_vname, dv_vname, psd_x_vname, psd_y_vname, psd_z_vname]
	
	;Create the plot
	win = MrVar_PlotTS( vars, $
	                    /NO_REFRESH, $
	                    XSIZE  = 680, $
	                    YSIZE  = 750 )
	
	;Overplot frequencies
	win = MrVar_OPlotTS( psd_x_vname, f_vnames )
	FOR i = 0, N_Elements(f_vnames) - 1 DO BEGIN
		oFreq  = MrVar_Get(f_vnames[i])
		oFreq -> RemoveAttr, 'LABEL'
	ENDFOR
	
	win = MrVar_OPlotTS( psd_y_vname, f_vnames )
	win = MrVar_OPlotTS( psd_z_vname, f_vnames )

	;Set legend properties
	lgds = win -> Get(/ALL, COUNT=nLgds, ISA='MrLegend')
	FOR i = 0, nLgds - 1 DO BEGIN
		IF StRegEx(lgds[i].name, 'f(ce|ci|pe|pi|lh)', /BOOLEAN) THEN BEGIN
			lgds[i] -> SetProperty, ALIGNMENT  ='NE', $
			                        COLOR      = 'Black', $
			                        FILL_COLOR = 'Light Grey', $
			                        LINESTYLE  = '-'
		ENDIF
	ENDFOR
	
	;Make pretty
	win.oxmargin = [12, 15]
	win[0]       -> SetLayout, [1,1]
	win          -> TrimLayout

;-------------------------------------------
; Save Figure //////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 $
			THEN output_dir = FilePath( '', ROOT_DIR=File_Search('~', /TEST_DIRECTORY), SUBDIRECTORY='figures' )
		
		;File name
		optdesc = 'n' + String(nfft, FORMAT='(i0)') + '-ns' + String(nshift, FORMAT='(i0)')
		fname = StrJoin( [sc, instr, mode, level, 'psd-'+optdesc], '_' )
		fname = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	win -> Refresh
	RETURN, win
END