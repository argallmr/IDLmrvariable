; docformat = 'rst'
;
; NAME:
;       MrMMS_EDI_Plot_ALT
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
; :Params:
;       SC:                 in, required, type=string/strarr
;                           The MMS spacecraft identifier. Options are:
;                               {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       OUTPUT_DIR:         in, optional, type=string, default=pwd
;                           A directory in which to save the figure. If neither `OUTPUT_DIR`
;                               nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT:         in, optional, type=string, default=pwd
;                           File extensions for the output figure. Options include: 'eps', 'gif',
;                               'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                               `OUTPUT_EXT` are defined, no file is generated.
;       OPTDESC:            in, optional, type=string, default=''
;                           Optional descriptor of the data. Options are:
;                               {'amb' | 'amb-pm2' | 'amb-alt-oom' | 'amb-alt-oob'}
;
; :Keywords:
;       NO_LOAD:            in, optional, type=boolean, default=0
;                           If set, data is not re-read and loaded into the variable cache.
;       TRANGE:             in, optional, type=strarr(2), default=MrVar_GetTRange
;                           The start and end times of the data interval to be loaded.
;                               Formatting is: 'YYYY-MM-DDThh:mm:ss'.
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
;       2016/10/03  -   Written by Matthew Argall
;-
function MrMMS_EDI_Plot_AMB_BrstCts, sc, $
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
		return, obj_new()
	endif
	
	;Defaults
	tf_load = ~keyword_set(no_load)
	if n_elements(trange) gt 0 then MrVar_SetTRange, trange
	
	;Variable type
	level = 'l2'
	type  = level eq 'l2' ? 'flux' : 'counts'

;-------------------------------------------
; EDI Data /////////////////////////////////
;-------------------------------------------
	;Get EDI data
	if tf_load then begin
		;Load EDI data
		MrMMS_Load_Data, sc, 'edi', 'brst', level, $
		                 OPTDESC   = ['amb', 'amb-pm2'], $
		                 VARFORMAT = '*' + type + '*'
		
		;Load EDI data
		MrMMS_Load_Data, sc, 'edi', 'srvy', level, $
		                 OPTDESC   = ['amb', 'amb-pm2'], $
		                 VARFORMAT = '*' + type + '*'
		
		;Load FPI data
		MrMMS_FPI_Load_Data, sc, 'brst', $
		                     OPTDESC   = 'des-moms', $
		                     VARFORMAT = '*energyspectr_' + ['par', 'perp', 'anti'] + '*'
	endif
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Parts
	prefix      = sc + '_edi_'
	suffix_srvy = '_srvy_' + level
	suffix_brst = '_brst_' + level
	ch          = ['1', '2', '3', '4']
	
	des_par_vname  = StrJoin( [sc, 'des', 'energyspectr', 'par',  'brst'], '_' )
	des_perp_vname = StrJoin( [sc, 'des', 'energyspectr', 'perp', 'brst'], '_' )
	des_anti_vname = StrJoin( [sc, 'des', 'energyspectr', 'anti', 'brst'], '_' )
	
	des_par_500eV_vname  = StrJoin( [sc, 'des', 'flux', 'par',  '500eV', 'brst'], '_' )
	des_perp_500eV_vname = StrJoin( [sc, 'des', 'flux', 'perp', '500eV', 'brst'], '_' )
	des_anti_500eV_vname = StrJoin( [sc, 'des', 'flux', 'anti', '500eV', 'brst'], '_' )
	
	;Names
	cts_srvy_vname = [ prefix   + type + ch[0] + '_0'   + suffix_srvy, $
	                   prefix   + type + ch[0] + '_180' + suffix_srvy ]
	cts_brst_vname = [ [ prefix + type + ch    + '_0'   + suffix_brst ], $
	                   [ prefix + type + ch    + '_180' + suffix_brst ] ]
	cts_brst_vname = reform(transpose(cts_brst_vname), 8)

;-------------------------------------------
; Set Properties ///////////////////////////
;-------------------------------------------
	;BRST
	yrange = [!values.f_infinity, -!values.f_infinity]
	for i = 0, n_elements(cts_brst_vname) - 2 do begin
		oVar = MrVar_Get( cts_brst_vname[i] )
		iCts = oVar -> Where(0, /GREATER)
		cmin = Min(oVar['DATA',iCts], MAX=cmax)
		yrange[0] <= cmin
		yrange[1] >= cmax
		if i eq 0 then oVar['LABEL'] = 'brst'
	endfor

	;SRVY
	for i = 0, n_elements(cts_srvy_vname) - 1 do begin
		oVar = MrVar_Get( cts_srvy_vname[i] )
		iCts = oVar -> Where(0, /GREATER)
		cmin = Min(oVar['DATA',iCts], MAX=cmax)
		yrange[0] <= cmin
		yrange[1] >= cmax
		oVar['COLOR'] = 'Red'
		oVar['SYMBOL'] = 1
		if i eq 0 then oVar['LABEL'] = 'srvy'
	endfor
	
	yrange[0] = 10.0^Floor(ALog10(yrange[0]))
	yrange[1] = 10.0^Ceil(ALog10(yrange[1]))

;-------------------------------------------
; Extract 500 eV Channel from FPI //////////
;-------------------------------------------
	oPar  = MrVar_Get(des_par_vname)
	oPerp = MrVar_Get(des_perp_vname)
	oAnti = MrVar_Get(des_anti_vname)
	
	;Indices of closest bin
	oE   = oPar['DEPEND_1']
	iE0  = MrNearestNeighbor(oE['DATA', 0, *], 500.0)
	iE1  = MrNearestNeighbor(oE['DATA', 1, *], 500.0)
	nPts = oE -> GetNPts()
	iE   = LonArr(nPts)
	iE[0:*:2] = iE0
	iE[1:*:2] = iE1
	
	;Energy scale factor
	oPar      = oPar[LIndGen(nPts), iE]
	odE_plus  = (oPar['DEPEND_1'])['DELTA_MINUS_VAR']
	odE_minus = (oPar['DEPEND_1'])['DELTA_MINUS_VAR']
	odE       = odE_plus + odE_minus
	e_scale   = 50.0 / odE['DATA']

	;Select the energy bin
	oPar  = e_scale * oPar
	oPar -> SetName, des_par_500eV_vname
	oPar -> Cache
	oPar['COLOR']     = 'Blue'
	oPar['LABEL']     = 'DES'
	oPar['SYMBOL']    = 9
	oPar['SYMTHICK']  = 2
	oPar -> RemoveAttr, 'DEPEND_1'
	
	oPerp  = e_scale * oPerp[LIndGen(nPts), iE]
	oPerp -> SetName, des_perp_500eV_vname
	oPerp -> Cache
	oPerp['COLOR']     = 'Blue'
	oPerp['LABEL']     = 'DES'
	oPerp['SYMBOL']    = 9
	oPerp['SYMTHICK']  = 2
	oPerp -> RemoveAttr, 'DEPEND_1'
	
	
	oAnti  = e_scale * oAnti[LIndGen(nPts), iE]
	oAnti -> SetName, des_anti_500eV_vname
	oAnti -> Cache
	oAnti['COLOR']     = 'Blue'
;	oAnti['LABEL']     = 'DES'
	oAnti['SYMCOLOR']  = 'Blue'
	oAnti['SYMBOL']    = 9
	oAnti['SYMTHICK']  = 2
	oAnti -> RemoveAttr, 'DEPEND_1'

;-------------------------------------------
; First Row ////////////////////////////////
;-------------------------------------------
	;Plot burst data
	win = MrVar_PlotTS( cts_brst_vname, $
	                    LAYOUT = [2,4], $
	                    /NO_REFRESH, $
	                    XSIZE  = 1000, $
;	                    /YERR, $
	                    YSIZE  = 600 )
	
	;YLOG -> Do not let error be negative
;	oPoly = win -> Get(ISA='MrPolygon', /ALL, COUNT=nPoly)
;	FOR i = 0, nPoly - 1 DO BEGIN
;		oPoly[i] -> GetData, x, y
;		oPoly[i] -> SetData, x, y>1
;		oPoly[i].noclip = 0
;	ENDFOR
	
	;Overplot survey data
	cts_brst_vname = Reform( cts_brst_vname, 2, 4 )
	win = MrVar_OPlotTS( cts_brst_vname[*,0], cts_srvy_vname )
	win = MrVar_OPlotTS( cts_brst_vname[*,1], cts_srvy_vname )
	win = MrVar_OPlotTS( cts_brst_vname[*,2], cts_srvy_vname )
	win = MrVar_OPlotTS( cts_brst_vname[*,3], cts_srvy_vname )
	win = MrVar_OPlotTS( cts_brst_vname[*,0], [des_par_500eV_vname, des_anti_500eV_vname] )
	win = MrVar_OPlotTS( cts_brst_vname[*,1], [des_par_500eV_vname, des_anti_500eV_vname] )
	win = MrVar_OPlotTS( cts_brst_vname[*,2], [des_par_500eV_vname, des_anti_500eV_vname] )
	win = MrVar_OPlotTS( cts_brst_vname[*,3], [des_par_500eV_vname, des_anti_500eV_vname] )
	
	;Create a legend
;	lgd = MrLegend( ALIGNMENT    = 'NW', $
;	                FILL_COLOR   = '', $
;	                LABEL        = ['Brst', 'Srvy'], $
;	                LINESTYLE    = 6, $
;	                POSITION     = [1.0, 1.0], $
;	                /RELATIVE, $
;	                SAMPLE_WIDTH = 0.0, $
;	                TARGET       = win[cts_brst_vname[0,3]], $
;	                TEXT_COLOR   = ['Black', 'Red'] )
	
	win -> SetGlobal, YRANGE=yrange
	win[cts_brst_vname[0]] -> SetLayout, [1,1]
	win -> TrimLayout
	win -> Refresh
	win -> Refresh

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
		
		;File name
		fname   = StrJoin( [sc, 'edi', 'brst', level, 'amb-srvy-brst-des-compare'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------
	return, win
end