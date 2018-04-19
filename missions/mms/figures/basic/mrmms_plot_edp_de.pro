; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_EDP_dE
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
; PURPOSE:
;+
;   Generate a plot of EDP quantities:
;        1. FGM Bxyz
;        2. Ex
;        3. dEx -- Ex filtered between F_PASS[0] and F_PASS[1]
;        4. Ex PSD + Gyrofrequency lines
;        5. Ey
;        6. dEy -- Ez filtered between F_PASS[0] and F_PASS[1]
;        7. Ey PSD + Gyrofrequency lines
;        8. Ez
;        9. dEz -- Ez filtered between F_PASS[0] and F_PASS[1]
;       10. Ez PSD + Gyrofrequency lines
;
; :Categories:
;   MMS
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' }
;
; :Keywords:
;       COORDS:     in, optional, type=string, default='gse'
;                   Coordinate system in which data is plotted. Options are: {'dsl' | 'gse' | 'fac'}
;       FGM_INSTR:  in, optional, type=string, default='fgm'
;                   FGM instrument to use. Options are: { 'afg' | 'dfg' | 'fgm' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       OPTDESC:    in, optional, type=string, default='dce'
;                   Optional filename descriptor. Options are: {'dce' | 'hmfe'}
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
;       F_PASS:     in, optional, type=fltarr(2), default=[20.0, 0.0]
;                   Frequency range in which the components of E are filtered. A value
;                       of 0.0 in the second element is equivalent to the Nyquist frequency.
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
;       2017/01/13  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_EDP_dE, sc, mode, nfft, nshift, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
OPTDESC=optdesc, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
F_PASS=f_pass, $
NO_LOAD=no_load, $
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
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(mode)      EQ 0 THEN mode      = 'fast'
	IF N_Elements(nfft)      EQ 0 THEN nfft      = 2048
	IF N_Elements(nshift)    EQ 0 THEN nshift    = nfft/2
	IF N_Elements(f_pass)    EQ 0 THEN f_pass    = [20.0, 0.0]
	IF N_Elements(optdesc)   EQ 0 THEN optdesc   = 'dce'
	IF N_Elements(coords)    EQ 0 THEN coords    = optdesc EQ 'hmfe' ? 'dsl' : 'gse'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange 
	
;-------------------------------------------
; Dataset Parameters ///////////////////////
;-------------------------------------------
	
	
	;EDP
	;   - OPDESC of 'scpot' is not allowed
	edp_instr  = 'edp'
	edp_coords = MrIsMember(['dmpa', 'dbcs'], coords) ? 'dsl' : coords
	IF coords EQ 'fac' THEN edp_coords = 'gse'
	edp_mode   = mode EQ 'brst' ? mode : 'fast'
	IF ~MrIsMember(['dce', 'hmfe'], optdesc) $
		THEN Message, 'Optional descriptor "' + optdesc + '" not allowed. Must be {"dce" | "hmfe"}' $
		ELSE edp_optdesc = optdesc
	
	;FGM
	fgm_coords = MrIsMember(['dbcs', 'dsl'], coords) ? 'dmpa' : coords
	IF coords EQ 'fac' THEN fgm_coords = 'gse'
	fgm_mode   = mode EQ 'brst' ? mode : 'srvy'
	CASE fgm_instr OF
		'fsm': fgm_level = 'l3'
		'fgm': fgm_level = 'l2'
		'dfg': fgm_level = 'l2pre'
		'afg': fgm_level = 'l2pre'
		ELSE: Message, 'FGM_INSTR not recognized: "' + fgm_instr + '".'
	ENDCASE
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------

	;Source names
	b_vname    = StrJoin( [sc, fgm_instr, 'b',     fgm_coords, fgm_mode, fgm_level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec',  fgm_coords, fgm_mode, fgm_level], '_' )
	bmag_vname = StrJoin( [sc, fgm_instr, 'bmag',  fgm_coords, fgm_mode, fgm_level], '_' )
	e_vname    = StrJoin( [sc, edp_instr, edp_optdesc, edp_coords, mode, level], '_' )
	epar_vname = StrJoin( [sc, edp_instr, edp_optdesc, 'par', 'epar', mode, level], '_' )
	
	;Output names
	ex_vname    = e_vname + '_x'
	ey_vname    = e_vname + '_y'
	ez_vname    = e_vname + '_z'
	dex_vname   = StrJoin( [sc, edp_instr, 'dex',       edp_coords, edp_mode, level], '_' )
	dey_vname   = StrJoin( [sc, edp_instr, 'dey',       edp_coords, edp_mode, level], '_' )
	dez_vname   = StrJoin( [sc, edp_instr, 'dez',       edp_coords, edp_mode, level], '_' )
	expsd_vname = StrJoin( [sc, edp_instr, 'ex', 'psd', edp_coords, edp_mode, level], '_' )
	eypsd_vname = StrJoin( [sc, edp_instr, 'ey', 'psd', edp_coords, edp_mode, level], '_' )
	ezpsd_vname = StrJoin( [sc, edp_instr, 'ez', 'psd', edp_coords, edp_mode, level], '_' )
	
	;Gyrofrequencies
	IF mode EQ 'brst' THEN BEGIN
		f1_vname = StrJoin( [sc, fgm_instr, 'fcE',     fgm_mode, fgm_level], '_' )
		f2_vname = StrJoin( [sc, fgm_instr, 'halffce', fgm_mode, fgm_level], '_' )
	ENDIF ELSE BEGIN
		f1_vname = StrJoin( [sc, fgm_instr, 'fcH',  fgm_mode, fgm_level], '_' )
		f2_vname = StrJoin( [sc, fgm_instr, 'fcHe', fgm_mode, fgm_level], '_' )
		f3_vname = StrJoin( [sc, fgm_instr, 'fcO',  fgm_mode, fgm_level], '_' )
	ENDELSE

;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = fgm_level, $
		                     VARFORMAT = b_vname

		;EDP
		MrMMS_Load_Data, sc, edp_instr, edp_mode, level, $
		                 OPTDESC   = edp_optdesc, $
		                 VARFORMAT = '*'+edp_optdesc+'_'+[edp_coords, 'par_epar']+'_'+edp_mode+'*'
	ENDIF
	
;-------------------------------------------
; Rotate to FAC ////////////////////////////
;-------------------------------------------
	IF coords EQ 'fac' THEN BEGIN
		;Get the data
		oE    = MrVar_Get(e_vname)
		oEpar = MrVar_Get(epar_vname)
		oB    = MrVar_Get(bvec_vname)
		
		;Rotate to FAC
		oT = MrVar_FAC(oB, '', 'CROSSX', TIME=oE)
		oE = oT ## oE
		oE -> Split, oEx, oEy, oEz
		
		;Replace parallel component
		;   - DCE data joins the Epar error with the data
		IF optdesc EQ 'hmfe' THEN BEGIN
			oEz = oEpar
		ENDIF ELSE BEGIN
			oEz = oEpar[*,1]
			oEz -> RemoveAttr, 'LABEL'
		ENDELSE
		
		;Ex
		oEx -> SetName, ex_vname
		oEx -> Cache
		oEx['CATDESC'] = 'Electric field perp1-component in field-aligned coordinates.'
		oEx['COORDS']  = 'FAC: PAR=B, PERP2=(PAR x [1,0,0]), PERP1=(PERP2 x PAR)'
		oEx['TITLE']   = 'E!D!9' + String(120B) + '!X1!N!C(mV/m)'
		oEx['UNITS']   = 'mV/m'
		
		;Ey
		oEy -> SetName, ey_vname
		oEy -> Cache
		oEy['CATDESC'] = 'Electric field perp2-component in field-aligned coordinates.'
		oEy['COORDS']  = 'FAC: PAR=B, PERP2=(PAR x [1,0,0]), PERP1=(PERP2 x PAR)'
		oEy['TITLE']   = 'E!D!9' + String(120B) + '!X2!C(mV/m)'
		oEy['UNITS']   = 'mV/m'
		
		;Ez
		oEz -> SetName, ez_vname
		oEz -> Cache
		oEz['TITLE']   = 'E$\down||$!C(mV/m)'
		oEz['UNITS']   = 'mV/m'
	
	;Split into components without rotating
	ENDIF ELSE BEGIN
		oE = MrVar_Get(e_vname)
		oE -> Split, oEx, oEy, oEz, /CACHE

		;Ex
		oEx['TITLE'] = 'Ex!C(mV/m)'
	
		;Ey
		oEy['TITLE'] = 'Ey!C(mV/m)'
	
		;Ez
		oEz['TITLE'] = 'Ez!C(mV/m)'
	ENDELSE

;-------------------------------------------
; Detrend Data /////////////////////////////
;-------------------------------------------
	sfrange = String(f_pass, FORMAT='(%"[%0.1f,, %0.1f]")')
	dt = oE['TIMEVAR'] -> GetSI(RATE=fs)
	fN = fs / 2.0
	f0 = (f_pass[0] EQ 0.0 ? 0.0 : f_pass[0]/fN) > 0.0
	f1 = (f_pass[1] EQ 0.0 ? 1.0 : f_pass[1]/fN) < 1.0
	A  = 75
	m  = 1024
	odEx = oEx -> Digital_Filter(f0, f1, A, m, /CACHE, NAME=dex_vname)
	odEy = oEy -> Digital_Filter(f0, f1, A, m, /CACHE, NAME=dey_vname)
	odEz = oEz -> Digital_Filter(f0, f1, A, m, /CACHE, NAME=dez_vname)
	
	;Ex
	odEx['CATDESC'] = 'Detrended electric field filtered between ' + sfrange + ' Hz.'
	odEx['TITLE']   = 'd' + oEx['TITLE']
	odEx['UNITS']   = oEx['UNITS']
	
	;Ey
	odEy['CATDESC'] = 'Detrended electric field filtered between ' + sfrange + ' Hz.'
	odEy['TITLE']   = 'd' + oEy['TITLE']
	odEy['UNITS']   = oEy['UNITS']
	
	;Ez
	odEz['CATDESC'] = 'Detrended electric field filtered between ' + sfrange + ' Hz.'
	odEz['TITLE']   = 'd' + oEz['TITLE']
	odEz['UNITS']   = oEz['UNITS']

;-------------------------------------------
; Power Spectral Density ///////////////////
;-------------------------------------------
	
	;Ex
	oEx_psd = oEx -> Spectrogram(nfft, nshift, NAME=expsd_vname, /CACHE, WINDOW='hanning')
	oEx_psd['AXIS_RANGE'] = [1e-7, 1e-0]
	oEx_psd['TITLE']      = 'Ex PSD'
	
	;Ey
	oEy_psd = oEy -> Spectrogram(nfft, nshift, NAME=eypsd_vname, /CACHE, WINDOW='hanning')
	oEy_psd['AXIS_RANGE'] = [1e-7, 1e-0]
	oEy_psd['TITLE']      = 'Ey PSD!C(mV/m)$\up2$/Hz'
	
	;Ez
	oEz_psd = oEz -> Spectrogram(nfft, nshift, NAME=ezpsd_vname, /CACHE, WINDOW='hanning')
	oEz_psd['AXIS_RANGE'] = [1e-7, 1e-0]
	oEz_psd['TITLE']      = 'Ez PSD'

;-------------------------------------------
; Gyrofrequency Lines //////////////////////
;-------------------------------------------
	oBmag  = MrVar_Get(bmag_vname)
	
	;Electron cyclotron frequency
	IF mode EQ 'brst' THEN BEGIN
		;fce
		of1 = MrVar_Freq_Cyclotron(oBmag, 'm_e', /CACHE, NAME=f1_vname)
		
		;0.5  * fce
		of2 = of1 / 2.0
		of2 -> SetName, f2_vname
		of2 -> Cache
		
		;FC_E
		of1['COLOR'] = 'White'
		of1['NSUM']  = 4
		
		;0.5*FC_$
		of2['COLOR']     = 'White'
		of2['LINESTYLE'] = '--'
		of2['NSUM']      = 4

	;Ion cyclotron frequencies
	ENDIF ELSE BEGIN
		of1 = MrVar_Freq_Cyclotron(oBmag, 'm_H',  /CACHE, NAME=f1_vname)
		of2 = MrVar_Freq_Cyclotron(oBmag, 'm_He', /CACHE, NAME=f2_vname)
		of3 = MrVar_Freq_Cyclotron(oBmag, 'm_O',  /CACHE, NAME=f3_vname)
		
		;FC_H
		of1['COLOR'] = 'Blue'
		of1['NSUM']  = 4
	
		;FC_HE
		of2['COLOR'] = 'Magenta'
		of2['NSUM']  = 4
	
		;FC_HE
		of3['COLOR'] = 'Orange'
		of3['NSUM']  = 4
	ENDELSE
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	;BMAG
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, optdesc[0]], ' ' ) )

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [b_vname, $
	                     ex_vname, dex_vname, expsd_vname, $
	                     ey_vname, dey_vname, eypsd_vname, $
	                     ez_vname, dez_vname, ezpsd_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 680, $
	                    YSIZE = 700 )
	
	;Cyclotron frequencies
	IF mode EQ 'brst' THEN BEGIN
		win = MrVar_OPlotTS( expsd_vname, [f1_vname, f2_vname] )
		win = MrVar_OPlotTS( eypsd_vname, [f1_vname, f2_vname] )
		win = MrVar_OPlotTS( ezpsd_vname, [f1_vname, f2_vname] )
	ENDIF ELSE BEGIN
		win = MrVar_OPlotTS( expsd_vname, [f1_vname, f2_vname, f3_vname] )
		win = MrVar_OPlotTS( eypsd_vname, [f1_vname, f2_vname, f3_vname] )
		win = MrVar_OPlotTS( ezpsd_vname, [f1_vname, f2_vname, f3_vname] )
	ENDELSE

	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 14]
	win    -> Refresh

;-------------------------------------------
; Save Figure //////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
		ENDIF ELSE IF ~File_Test(output_dir, /DIRECTORY) THEN BEGIN
			MrPrintF, 'LogText', 'Creating directory: "' + output_dir + '".'
			File_MKDir, output_dir
		ENDIF
		
		outdesc = String(optdesc, coords, nfft, nshift, f_pass, FORMAT='(%"%s-de-%s-n%i-s%i-f%i-%i")')
		
		;File name
		fname = StrJoin( [sc, 'edp', mode, level, outdesc], '_' )
		fname = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Finished! ////////////////////////////////
;-------------------------------------------
	RETURN, win
END