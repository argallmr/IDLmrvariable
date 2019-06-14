; docformat = 'rst'
;
; NAME:
;       MrMMS_EDI_Plot_ALT_PSD
;
;*****************************************************************************************
;   Copyright (c) 2019, Matthew Argall                                                   ;
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
;   Load EDI alternating-mode data (as opposed to field-aligned and perpendicular).
;   Plot Survey Flux/Counts on top of Burst Flux/Counts to ensure calibration is
;   applied correctly.
;       1.1 Flux 0-PA CH1   2.1 Flux 90-PA GDU1 CH1  3.1 Flux 90-PA GDU2 CH1  4.1 Flux 180-PA CH1
;       1.2 Flux 0-PA CH2   2.2 Flux 90-PA GDU1 CH2  3.2 Flux 90-PA GDU2 CH2  4.2 Flux 180-PA CH2
;       1.3 Flux 0-PA CH3   2.3 Flux 90-PA GDU1 CH3  3.3 Flux 90-PA GDU2 CH3  4.2 Flux 180-PA CH3
;       1.4 Flux 0-PA CH4   2.4 Flux 90-PA GDU1 CH4  3.4 Flux 90-PA GDU2 CH4  4.4 Flux 180-PA CH4
;
; :Categories:
;       MMS, EDI, MrVariable
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
;       2019/04/02  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Splice Field-Aligned and Perpendicular mode data together to create time series
;   products that are continuous in time (i.e. they do not have the 32ms data gap that
;   exists as EDI alternates between the two modes).
;
; :Params:
;       CTS_0_VNAMES:       in, required, type=string/strarr
;                           Names of the 0-degree pitch angle particle flux variables.
;       CTS_90_GDU1_VNAMES: in, required, type=string/strarr
;                           Names of the 90-degree pitch angle GDU 1 particle flux variables.
;       CTS_90_GDU2_VNAMES: in, required, type=string/strarr
;                           Names of the 90-degree pitch angle GDU 2 particle flux variables.
;       CTS_180_VNAMES:     in, required, type=string/strarr
;                           Names of the 180-degree pitch angle particle flux variables.
;
; :Returns:
;       OVARS:              out, required, type=ObjArr
;                           MrScalarTS object references for the concatenated inputs, each
;                           with the 32ms data gaps filled with zeroes (0s).
;-
FUNCTION edi_alt_splice, cts_0_vnames, cts_90_gdu1_vnames, cts_90_gdu2_vnames, cts_180_vnames
	Compile_Opt idl2
	On_Error, 2
	
	;EDI is in either field-aligned or perpendicular mode at any one time
	c_fa   = MrVar_Get(cts_0_vnames[0])
	c_perp = MrVar_Get(cts_90_gdu1_vnames[0])
	t_fa   = c_fa['TIME', 'SSM']
	t_perp = c_perp['TIME', 'SSM']
	n_fa   = N_Elements(t_fa)
	n_perp = N_Elements(t_perp)
	
	;Combine times
	time   = [t_fa, t_perp]
	isort  = Sort(time)
	time   = time[isort]
	i_fa = Where(isort LT n_fa, COMPLEMENT=i_perp)
	
	;Splice together the data
	cts_vnames = [cts_0_vnames, cts_90_gdu1_vnames, cts_90_gdu2_vnames, cts_180_vnames]
	nvars = N_Elements(cts_vnames)
	oVars = ObjArr(nvars)
	FOR i = 0, nvars - 1 DO BEGIN
		;Splice data together
		cts = FltArr(n_fa + n_perp)
		oCts = MrVar_Get(cts_vnames[i])
		IF StRegEx(oCts.name, '90', /BOOLEAN) $
			THEN cts[i_perp] = oCts['DATA'] $
			ELSE cts[i_fa]   = oCts['DATA']
		
		;Create a new variable
		oVars[i] = MrScalarTS( time, cts, $
;		                       /CACHE, $
;		                       NAME   = vnames[i], $
		                       T_REF  = c_fa['TIME', 0], $
		                       T_TYPE = 'SSM' )
		oCts -> CopyAttrTo, oVars[i]
	ENDFOR
	
	RETURN, oVars
END


;+
;   Plot the power spectral density of EDI alternating mode particle fluxes. Both survey
;   and burst data are plotted, as well as the noise due to statistical uncertainty. Plots
;   are sorted by pitch angle, each with spectra from the four channels.
;
; :Params:
;       SC:                 in, required, type=string/strarr
;                           The MMS spacecraft identifier. Options are:
;                               {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       FFT_INTERVAL:       in, optional, type=float, default=16.0
;                           Seconds worth of data to be included in each spectra.
;                               Consecutive spectral windows have 50% overlap and are
;                               averaged together.
;
; :Keywords:
;       NO_LOAD:            in, optional, type=boolean, default=0
;                           If set, data is not re-read and loaded into the variable cache.
;       OUTPUT_DIR:         in, optional, type=string, default=pwd
;                           A directory in which to save the figure. If neither `OUTPUT_DIR`
;                               nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT:         in, optional, type=string, default=pwd
;                           File extensions for the output figure. Options include: 'eps', 'gif',
;                               'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                               `OUTPUT_EXT` are defined, no file is generated.
;       TRANGE:             in, optional, type=strarr(2), default=MrVar_GetTRange
;                           The start and end times of the data interval to be loaded.
;                               Formatting is: 'YYYY-MM-DDThh:mm:ss'.
;
;-
FUNCTION MrMMS_EDI_Plot_ALT_PSD, sc, fft_interval, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
TRANGE=trange
	Compile_Opt idl2

	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		catch, /CANCEL
		IF N_Elements(wins) GT 0 THEN Obj_Destroy, wins
		MrPrintF, 'LogErr'
		RETURN, Obj_New()
	ENDIF
	
	;Defaults
	tf_load = ~Keyword_Set(no_load)
	T = N_Elements(fft_interval) EQ 0 ? 16.0 : fft_interval
	IF N_Elements(trange) GT 0 THEN MrVar_SetTRange, trange
	
	;Variable type
	level = 'l2'
	type  = level EQ 'l2' ? 'flux' : 'counts'

;-------------------------------------------
; Get File Types ///////////////////////////
;-------------------------------------------
	files = MrMMS_Get_Filenames(sc, 'edi', ['srvy','brst'], level, OPTDESC=['amb-alt-oob', 'amb-alt-oc', 'amb-alt-cc'], COUNT=count)
	IF count EQ 0 THEN Message, 'No EDI files found.'
		
	;Which files
	MrMMS_Parse_Filename, files, OPTDESC=optdesc
	iUniq    = Uniq(optdesc, Sort(optdesc))
	optdesc  = optdesc[iUniq]
	IF N_Elements(optdesc) GT 1 THEN BEGIN
		MrPrintF, 'LogText', 'More than one EDI optional descriptor found. Choosing first: "' + optdesc[0] + '".'
		optdesc = optdesc[0]
	ENDIF
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	
	ch    = ['1', '2', '3', '4']
	nChan = N_Elements(ch)
	type  = level EQ 'l2' ? 'flux' : 'counts'
	
	;SRVY
	cts_0_srvy_vname        = StrJoin( [sc, 'edi', type+ch[0],   '0', 'srvy', level], '_' )
	cts_90_gdu1_srvy_vname  = StrJoin( [sc, 'edi', type+ch[0],  '90', 'gdu1', 'srvy', level], '_' )
	cts_90_gdu2_srvy_vname  = StrJoin( [sc, 'edi', type+ch[0],  '90', 'gdu2', 'srvy', level], '_' )
	cts_180_srvy_vname      = StrJoin( [sc, 'edi', type+ch[0], '180', 'srvy', level], '_' )
	
	delta_0_srvy_vname        = StrJoin( [sc, 'edi', type+ch[0],   '0', 'delta', 'srvy', level], '_' )
	delta_90_gdu1_srvy_vname  = StrJoin( [sc, 'edi', type+ch[0],  '90', 'delta', 'gdu1', 'srvy', level], '_' )
	delta_90_gdu2_srvy_vname  = StrJoin( [sc, 'edi', type+ch[0],  '90', 'delta', 'gdu2', 'srvy', level], '_' )
	delta_180_srvy_vname      = StrJoin( [sc, 'edi', type+ch[0], '180', 'delta', 'srvy', level], '_' )

	;BRST
	cts_0_vnames        = StrJoin( [sc, 'edi', type,   '0', 'brst', level], '_' )
	cts_90_gdu1_vnames  = StrJoin( [sc, 'edi', type,  '90', 'gdu1', 'brst', level], '_' )
	cts_90_gdu2_vnames  = StrJoin( [sc, 'edi', type,  '90', 'gdu2', 'brst', level], '_' )
	cts_180_names       = StrJoin( [sc, 'edi', type, '180', 'brst', level], '_' )
	
	delta_0_vnames        = StrJoin( [sc, 'edi', type,   '0', 'delta', 'brst', level], '_' )
	delta_90_gdu1_vnames  = StrJoin( [sc, 'edi', type,  '90', 'delta', 'gdu1', 'brst', level], '_' )
	delta_90_gdu2_vnames  = StrJoin( [sc, 'edi', type,  '90', 'delta', 'gdu2', 'brst', level], '_' )
	delta_180_names       = StrJoin( [sc, 'edi', type, '180', 'delta', 'brst', level], '_' )
	
	;Add channel index to burst variable names
	cts_0_vnames        = StrMid(cts_0_vnames, 0, 13)       + ch + StrMid(cts_0_vnames, 13)
	cts_90_gdu1_vnames  = StrMid(cts_90_gdu1_vnames, 0, 13) + ch + StrMid(cts_90_gdu1_vnames, 13)
	cts_90_gdu2_vnames  = StrMid(cts_90_gdu2_vnames, 0, 13) + ch + StrMid(cts_90_gdu2_vnames, 13)
	cts_180_vnames      = StrMid(cts_180_names, 0, 13)      + ch + StrMid(cts_180_names, 13)
	
	delta_0_vnames        = StrMid(delta_0_vnames, 0, 13)       + ch + StrMid(delta_0_vnames, 13)
	delta_90_gdu1_vnames  = StrMid(delta_90_gdu1_vnames, 0, 13) + ch + StrMid(delta_90_gdu1_vnames, 13)
	delta_90_gdu2_vnames  = StrMid(delta_90_gdu2_vnames, 0, 13) + ch + StrMid(delta_90_gdu2_vnames, 13)
	delta_180_vnames      = StrMid(delta_180_names, 0, 13)      + ch + StrMid(delta_180_names, 13)
	
	;Output PSD
	psd_0_srvy_vname       = StrMid(cts_0_srvy_vname, 0, 9)       + 'psd' + StrMid(cts_0_srvy_vname, 13, 1)       + StrMid(cts_0_srvy_vname, 14)
	psd_90_gdu1_srvy_vname = StrMid(cts_90_gdu1_srvy_vname, 0, 9) + 'psd' + StrMid(cts_90_gdu1_srvy_vname, 13, 1) + StrMid(cts_90_gdu1_srvy_vname, 14)
	psd_90_gdu2_srvy_vname = StrMid(cts_90_gdu2_srvy_vname, 0, 9) + 'psd' + StrMid(cts_90_gdu2_srvy_vname, 13, 1) + StrMid(cts_90_gdu2_srvy_vname, 14)
	psd_180_srvy_vname     = StrMid(cts_180_srvy_vname, 0, 9)     + 'psd' + StrMid(cts_180_srvy_vname, 13, 1)     + StrMid(cts_180_srvy_vname, 14)
	psd_0_vnames       = StrMid(cts_0_vnames, 0, 9)       + 'psd' + StrMid(cts_0_vnames, 13, 1)       + StrMid(cts_0_vnames, 14)
	psd_90_gdu1_vnames = StrMid(cts_90_gdu1_vnames, 0, 9) + 'psd' + StrMid(cts_90_gdu1_vnames, 13, 1) + StrMid(cts_90_gdu1_vnames, 14)
	psd_90_gdu2_vnames = StrMid(cts_90_gdu2_vnames, 0, 9) + 'psd' + StrMid(cts_90_gdu2_vnames, 13, 1) + StrMid(cts_90_gdu2_vnames, 14)
	psd_180_vnames     = StrMid(cts_180_vnames, 0, 9)     + 'psd' + StrMid(cts_180_vnames, 13, 1)     + StrMid(cts_180_vnames, 14)

	psd_delta_0_srvy_vname       = StrMid(delta_0_srvy_vname, 0, 9)       + 'psd' + StrMid(delta_0_srvy_vname, 13, 1)       + StrMid(delta_0_srvy_vname, 14)
	psd_delta_90_gdu1_srvy_vname = StrMid(delta_90_gdu1_srvy_vname, 0, 9) + 'psd' + StrMid(delta_90_gdu1_srvy_vname, 13, 1) + StrMid(delta_90_gdu1_srvy_vname, 14)
	psd_delta_90_gdu2_srvy_vname = StrMid(delta_90_gdu2_srvy_vname, 0, 9) + 'psd' + StrMid(delta_90_gdu2_srvy_vname, 13, 1) + StrMid(delta_90_gdu2_srvy_vname, 14)
	psd_delta_180_srvy_vname     = StrMid(delta_180_srvy_vname, 0, 9)     + 'psd' + StrMid(delta_180_srvy_vname, 13, 1)     + StrMid(delta_180_srvy_vname, 14)
	psd_delta_0_vnames       = StrMid(delta_0_vnames, 0, 9)       + 'psd' + StrMid(delta_0_vnames, 13, 1)       + StrMid(delta_0_vnames, 14)
	psd_delta_90_gdu1_vnames = StrMid(delta_90_gdu1_vnames, 0, 9) + 'psd' + StrMid(delta_90_gdu1_vnames, 13, 1) + StrMid(delta_90_gdu1_vnames, 14)
	psd_delta_90_gdu2_vnames = StrMid(delta_90_gdu2_vnames, 0, 9) + 'psd' + StrMid(delta_90_gdu2_vnames, 13, 1) + StrMid(delta_90_gdu2_vnames, 14)
	psd_delta_180_vnames     = StrMid(delta_180_vnames, 0, 9)     + 'psd' + StrMid(delta_180_vnames, 13, 1)     + StrMid(delta_180_vnames, 14)
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		MrMMS_Load_Data, sc, 'edi', ['srvy', 'brst'], level, $
		                 OPTDESC   = optdesc, $
		                 VARFORMAT = ['*_'+type+'?_0_*', $
		                              '*_'+type+'?_90_*', $
		                              '*_'+type+'?_180_*']
	ENDIF
	
;-------------------------------------------
; Power Spectral Density ///////////////////
;-------------------------------------------
if tf_load then begin
	;FFT Parameters
	oSrvy = MrVar_Get(cts_0_srvy_vname)
	dt_srvy = oSrvy['TIMEVAR'] -> GetSI()
	fs_srvy = Round(1.0/dt_srvy)
	dt_srvy = 1.0/fs_srvy
	
	oBrst = MrVar_Get(cts_0_vnames[0])
	dt_brst = oBrst['TIMEVAR'] -> GetSI()
	fs_brst = Round(1.0/dt_brst)
	dt_brst = 1.0/fs_brst
	
	;Compute PSD
	in_vnames = [ [cts_0_srvy_vname,       cts_0_vnames], $
	              [cts_90_gdu1_srvy_vname, cts_90_gdu1_vnames], $
	              [cts_90_gdu2_srvy_vname, cts_90_gdu2_vnames], $
	              [cts_180_srvy_vname,     cts_180_vnames] ]
	out_vnames = [ [psd_0_srvy_vname,       psd_0_vnames], $
	               [psd_90_gdu1_srvy_vname, psd_90_gdu1_vnames], $
	               [psd_90_gdu2_srvy_vname, psd_90_gdu2_vnames], $
	               [psd_180_srvy_vname,     psd_180_vnames] ]
	colors = ['Black', 'Red', 'Forest Green', 'Blue']
	labels = ['Ch1', 'Ch2', 'Ch3', 'Ch4']
	titles = ['PSD!CPA0', 'PSD!CPA90 GDU1', 'PSD!CPA90 GDU2', 'PSD!CPA180']
	yrange = [!values.f_nan, -!values.f_nan]
	FOR iCh = 0, N_Elements(in_vnames[*,0]) - 1 DO BEGIN
		oCts = edi_alt_splice(in_vnames[iCh,0], in_vnames[iCh,1], in_vnames[iCh,2], in_vnames[iCh,3])
		
		;PSD
		fs = iCh EQ 0 ? fs_srvy : fs_brst
		FOR iPA = 0, N_Elements(oCts) - 1 DO BEGIN
			oPSD = oCts[iPA] -> PSD(fs*T, WINDOW='hamming', /ZEROPAD)
			oPSD -> SetName, out_vnames[iCh,iPA]
			oPSD -> Cache
			
			;Set plot properties
			oPSD['TITLE'] = titles[iPA]
			IF iPA EQ 0 THEN BEGIN
				IF iCh EQ 0 $
					THEN oPSD['LABEL'] = 'srvy' $
					ELSE oPSD['LABEL'] = labels[iCh-1]
			ENDIF
			IF iCh EQ 0 $
				THEN oPSD['COLOR'] = 'Magenta' $
				ELSE oPSD['COLOR'] = colors[iCh-1]
			IF iPA LT 3 THEN BEGIN
				oFreq = oPSD['DEPEND_0']
				oFreq['TICKFORMAT'] = '(a1)'
				oFreq['TITLE']      = ''
			ENDIF
			yrange[0] <= oPSD.min
			yrange[1] >= oPSD.max
			
			;Clean up data
			Obj_Destroy, oCts[iPA]
		ENDFOR
	ENDFOR

;-------------------------------------------
; Errors ///////////////////////////////////
;-------------------------------------------
	;https://www.mathworks.com/help/matlab/math/random-numbers-with-specific-mean-and-variance.html
	in_vnames  = [delta_0_vnames[0], delta_90_gdu1_vnames[0], delta_90_gdu2_vnames[0], delta_180_vnames[0]]
	out_vnames = [psd_delta_0_vnames[0], psd_delta_90_gdu1_vnames[0], psd_delta_90_gdu2_vnames[0], psd_delta_180_vnames[0]]
	oDelta     = edi_alt_splice(in_vnames[0], in_vnames[1], in_vnames[2], in_vnames[3])
	FOR i = 0, N_Elements(in_vnames) - 1 DO BEGIN
		npts  = N_Elements(oDelta[i])
		oRand = oDelta[i] * RandomN(seed, npts)
		oPSD  = oRand -> PSD(fs_brst*T, WINDOW='hamming', /ZEROPAD)
		oPSD -> SetName, out_vnames[i]
		oPSD -> Cache
		oPSD['COLOR'] = 'Grey'
		IF i EQ 0 THEN oPSD['LABEL'] = 'noise'
		yrange[0] <= oPSD.min
		yrange[1] <= oPSD.max
		
		Obj_Destroy, oDelta[i]
		Obj_Destroy, oRand
	ENDFOR
	
	yrange[0] = 10.0^Floor(ALog10(yrange[0]))
	yrange[1] = 10.0^Ceil(ALog10(yrange[1]))
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	title = MrVar_GetTRange()
	title = 'EDI ' + StrUpCase(optdesc) + '!C' + StrMid(title[0], 0, 10) + ' ' + StrJoin(StrMid(title, 11, 8), ' - ') 
	win = MrWindow( LAYOUT  = [1, 4], $
	                NAME    = 'amb-alt-psd', $
	                REFRESH = 0, $
	                XSIZE   = 500, $
	                YGAP    = 0.5, $
	                YSIZE   = 700 )
	
	p1 = MrVar_Plot( psd_0_vnames[0],       /CURRENT, TITLE=title )
	p2 = MrVar_Plot( psd_90_gdu1_vnames[0], /CURRENT )
	p3 = MrVar_Plot( psd_90_gdu2_vnames[0], /CURRENT )
	p4 = MrVar_Plot( psd_180_vnames[0],     /CURRENT )
	FOR i = 1, nChan - 1 DO BEGIN
		op1 = MrVar_Plot( psd_0_vnames[i],       OVERPLOT=p1 )
		op2 = MrVar_Plot( psd_90_gdu1_vnames[i], OVERPLOT=p2 )
		op3 = MrVar_Plot( psd_90_gdu2_vnames[i], OVERPLOT=p3 )
		op4 = MrVar_Plot( psd_180_vnames[i],     OVERPLOT=p4 )
	ENDFOR
	
	;Srvy
	op1 = MrVar_Plot( psd_0_srvy_vname,        OVERPLOT=p1 )
	op2 = MrVar_Plot( psd_90_gdu1_srvy_vname,  OVERPLOT=p2 )
	op3 = MrVar_Plot( psd_90_gdu2_srvy_vname,  OVERPLOT=p3 )
	op4 = MrVar_Plot( psd_180_srvy_vname,      OVERPLOT=p4 )
	
	;Delta
	op1 = MrVar_Plot( psd_delta_0_vnames[0],        OVERPLOT=p1 )
	op2 = MrVar_Plot( psd_delta_90_gdu1_vnames[0],  OVERPLOT=p2 )
	op3 = MrVar_Plot( psd_delta_90_gdu2_vnames[0],  OVERPLOT=p3 )
	op4 = MrVar_Plot( psd_delta_180_vnames[0],      OVERPLOT=p4 )
	
;-------------------------------------------
; Save /////////////////////////////////////
;-------------------------------------------
	
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
		ENDIF ELSE IF ~File_Test(output_dir, /DIRECTORY) THEN BEGIN
			MrPrintF, 'LogText', 'Creating directory: "' + output_dir + '".'
			File_MKDir, output_dir
		ENDIF
		
		;File name
		fname1 = StrJoin( [sc, 'edi', 'brst', level, optdesc+'-psd'], '_' )
		fname1 = FilePath( fname1, ROOT_DIR=output_dir )
		
		;Save the figure
		fout1 = MrVar_PlotTS_Save( win, fname1, output_ext )
	ENDIF
	
;-------------------------------------------
; Finish! //////////////////////////////////
;-------------------------------------------
	win[0] -> SetLayout, [1,1]
	win.oymargin = [4,4]
	win.oxmargin = [10,7]
	win -> SetGlobal, YRANGE=yrange
	win -> TrimLayout
	win -> Refresh
	
	RETURN, win
END