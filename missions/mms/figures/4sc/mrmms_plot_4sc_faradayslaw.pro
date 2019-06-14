; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_FaradaysLaw
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
; PURPOSE:
;+
;   Generate a plot to provide an overview of Faraday's Law:
;       1. B barycentric average
;       2. -dBx/dt
;       3. -dBy/dt
;       4. -dBz/dt
;       5. Ex
;       6. Ey
;       7. Ez
;       8. Curl(E)
;
; :Categories:
;   MMS
;
; :Params:
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;
; :Keywords:
;       EPHDESC:    in, optional, type=string, default='ephts04d'
;                   Optional descriptor of the definitive ephemeris datatype to use.
;                       Options are: { 'epht89d' | 'epht89q' | 'ephts04d' | 'defeph' | 'predeph' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
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
;       2018/09/17  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_FaradaysLaw, mode, $
FC = fc, $
FGM_INSTR = fgm_instr, $
EPHDESC=ephdesc, $
LEVEL=level, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
NO_LOAD=no_load, $
TRANGE=trange
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF N_Elements(win) GT 0 THEN Obj_Destroy, win
		IF N_Elements(w2)  GT 0 THEN Obj_Destroy, w2
		IF N_Elements(w3)  GT 0 THEN Obj_Destroy, w3
		IF N_Elements(w4)  GT 0 THEN Obj_Destroy, w4
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	tf_load   = ~Keyword_Set(no_load)
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(coords) EQ 0 THEN coords = 'gse'
	IF N_Elements(fc)     EQ 0 THEN fc     = [0.0, 1.0]
	IF N_Elements(level)  EQ 0 THEN level  = 'l2'
	IF N_Elements(trange) GT 0 THEN MrVar_SetTRange, trange
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	
	e0        = MrConstants('epsilon_0')
	q         = MrConstants('q')
	sc_colors = ['Black', 'Red', 'Green', 'Blue']
	dot       = '!9'+String(46B)+'!X'
	nabla     = '!9'+String(71B)+'!X'
	partial   = '!9'+String(68B)+'!X'
	fancyf    = '!9'+String(70B)+'!X'
	
;-------------------------------------------
; Variable Parameters //////////////////////
;-------------------------------------------
	;FGM
	fgm_mode    = mode EQ 'brst' ? mode : 'srvy'
	fgm_coords  = MrIsMember(['dsl', 'dbcs'], coords) ? 'dmpa' : coords
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	CASE fgm_instr OF
		'fsm': fgm_level = 'l3'
		'fgm': fgm_level = 'l2'
		'dfg': fgm_level = 'l2pre'
		'afg': fgm_level = 'l2pre'
		ELSE: Message, 'Invalid FGM instrument: "' + fgm_instr + '".'
	ENDCASE
	
	;EDP
	edp_instr  = 'edp'
	edp_coords = MrIsMember(['dmpa', 'dbcs'], coords) ? 'dsl' : coords
	edp_mode   = mode
	IF mode EQ 'srvy' THEN BEGIN
		MrPrintF, 'LogWarn', 'EDP does not have "srvy" data. Using "fast".'
		edp_mode = 'fast'
	ENDIF
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	sc  = 'mms' + ['1', '2', '3', '4']
	nSC = N_Elements(sc)
	
	;Source names
	b_vnames    = sc + '_' + StrJoin([fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_')
	IF fgm_instr EQ 'fsm' THEN BEGIN
		bvec_vnames = sc + '_' + StrJoin([fgm_instr, 'b', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([fgm_instr, 'b', 'mag',      fgm_mode, fgm_level], '_')
	ENDIF ELSE BEGIN
		bvec_vnames = sc + '_' + StrJoin([fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_') 
		bmag_vnames = sc + '_' + StrJoin([fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_')
	ENDELSE
	e_vnames    = sc + '_' + StrJoin([edp_instr, 'dce',  edp_coords, edp_mode, level], '_' )
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	bf_vnames     = sc + '_' + StrJoin([fgm_instr, 'bvec', fgm_coords, 'lpfilt', fgm_mode, fgm_level], '_')
	bx_vnames     = bf_vnames + '_x'
	by_vnames     = bf_vnames + '_y'
	bz_vnames     = bf_vnames + '_z'
	ef_vnames     = sc + '_' + StrJoin([edp_instr, 'dce',  edp_coords, 'lpfilt', edp_mode, level], '_' )
	ex_vnames     = ef_vnames + '_x'
	ey_vnames     = ef_vnames + '_y'
	ez_vnames     = ef_vnames + '_z'
	dB_dt_vnames  = sc + '_' + StrJoin( ['fgm', 'bdot', fgm_coords, 'bary', fgm_mode, fgm_level], '_' )
	dBx_dt_vnames = dB_dt_vnames + '_x'
	dBy_dt_vnames = dB_dt_vnames + '_y'
	dBz_dt_vnames = dB_dt_vnames + '_z'
	b_bary_vname  = StrJoin( ['mms', 'fgm', 'b', fgm_coords, 'bary', fgm_mode, fgm_level], '_' )
	bx_bary_vname = b_bary_vname + '_x'
	by_bary_vname = b_bary_vname + '_y'
	bz_bary_vname = b_bary_vname + '_z'
	e_bary_vname  = StrJoin( ['mms', edp_instr, 'e', edp_coords, 'bary', edp_mode, level], '_' )
	ex_bary_vname = e_bary_vname + '_x'
	ey_bary_vname = e_bary_vname + '_y'
	ez_bary_vname = e_bary_vname + '_z'
	db_dt_bary_vname = StrJoin( ['mms', 'fgm', 'bdot', fgm_coords, 'bary', fgm_mode, fgm_level], '_' )
	dbx_dt_bary_vname = db_dt_bary_vname + '_x'
	dby_dt_bary_vname = db_dt_bary_vname + '_y'
	dbz_dt_bary_vname = db_dt_bary_vname + '_z'
	db_dt_bary_mag_vname = StrJoin( ['mms', 'fgm', 'bdot', 'mag', 'bary', fgm_mode, fgm_level], '_' )
	db_dt_bary_corr_vname = StrJoin( ['mms', 'fgm', 'bdot', fgm_coords, 'bary', 'corrected', fgm_mode, fgm_level], '_' )
	dbx_dt_bary_corr_vname = db_dt_bary_corr_vname + '_x'
	dby_dt_bary_corr_vname = db_dt_bary_corr_vname + '_y'
	dbz_dt_bary_corr_vname = db_dt_bary_corr_vname + '_z'
	divE_vname     = StrJoin( ['mms', 'dce', 'dive', edp_mode, level], '_' )
	divE_mag_vname = StrJoin( ['mms', 'dce', 'dive', 'mag', edp_mode, level], '_' )
	divE_err_vname = StrJoin( ['mms', 'dce', 'dive', 'err', edp_mode, level], '_' )
	curlE_vname    = StrJoin( ['mms', 'dce', 'curle', edp_mode, level], '_' )
	curlEx_vname   = curlE_vname + '_x'
	curlEy_vname   = curlE_vname + '_y'
	curlEz_vname   = curlE_vname + '_z'
	curlE_mag_vname    = StrJoin( ['mms', 'dce', 'curle', 'mag', edp_mode, level], '_' )
	curlE_corr_vname   = StrJoin( ['mms', 'dce', 'curle', 'corrected', edp_mode, level], '_' )
	curlEx_corr_vname  = curlE_corr_vname + '_x'
	curlEy_corr_vname  = curlE_corr_vname + '_y'
	curlEz_corr_vname  = curlE_corr_vname + '_z'
		
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;B-Field
		IF fgm_instr EQ 'fsm' THEN BEGIN
			!MrMMS.dropbox_root = '/nfs/fsm/temp/'
			MrMMS_Load_Data, '', fgm_instr, fgm_mode, fgm_level, $
			                 /TEAM_SITE, $
			                 OPTDESC  = fgm_optdesc, $
			                 VARFORMAT = '*_b_'+fgm_coords+'_'+fgm_mode+'*'
		ENDIF ELSE BEGIN
			MrMMS_FGM_Load_Data, '', fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = fgm_level, $
			                     OPTDESC   = fgm_optdesc, $
			                     VARFORMAT = '*_b_'+fgm_coords+'_'+fgm_mode+'*'
		ENDELSE
		
		;E-Field
		MrMMS_Load_Data, '', instr, mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = '*dce_'+coords+'*'
		
		;MEC
		IF mode EQ 'brst' THEN BEGIN
			IF (mrmms_get_filenames('mms1', 'mec', mode, level, OPTDESC=ephdesc))[0] EQ '' $
				THEN mec_mode = 'srvy' $
				ELSE mec_mode = mode
		ENDIF
		
		;Ephemeris
		MrMMS_Load_Data, '', 'mec', mec_mode, level, $
		                 OPTDESC   = ephdesc, $
		                 VARFORMAT = r_vnames
	ENDIF
	
;-------------------------------------------
; Interpolate //////////////////////////////
;-------------------------------------------
	;Interpolate to MMS1 FGM
	oB1   = MrVar_Get(b_vnames[0])
	oTref = oB1['TIMEVAR']
	dt    = oTref -> GetSI(RATE=fs)
	
	;Allocate space for inteprolated arrays
	aB      = ObjArr(nSC)
	aBx     = ObjArr(nSC)
	aBy     = ObjArr(nSC)
	aBz     = ObjArr(nSC)
	adB_dt  = ObjArr(nSC)
	adBx_dt = ObjArr(nSC)
	adBy_dt = ObjArr(nSC)
	adBz_dt = ObjArr(nSC)
	aE      = ObjArr(nSC)
	aEx     = ObjArr(nSC)
	aEy     = ObjArr(nSC)
	aEz     = ObjArr(nSC)
	aR      = ObjArr(nSC)
	
	bxrange = [!Values.f_infinity, -!Values.f_infinity]
	byrange = [!Values.f_infinity, -!Values.f_infinity]
	bzrange = [!Values.f_infinity, -!Values.f_infinity]
	dbxdtrange = [!Values.f_infinity, -!Values.f_infinity]
	dbydtrange = [!Values.f_infinity, -!Values.f_infinity]
	dbzdtrange = [!Values.f_infinity, -!Values.f_infinity]
	exrange = [!Values.f_infinity, -!Values.f_infinity]
	eyrange = [!Values.f_infinity, -!Values.f_infinity]
	ezrange = [!Values.f_infinity, -!Values.f_infinity]
	
	;Interpolate all
	FOR i = 0, nSC - 1 DO BEGIN
		sSC = 'MMS' + String(i+1, FORMAT='(i1)')
		
		;Grab the vectors
		oB   = MrVar_Get(bvec_vnames[i])
		oE   = MrVar_Get(e_vnames[i])
		oPos = MrVar_Get(r_vnames[i])

		;Interpolate
		aB[i] = MrVar_Resample(oB, oTref, /CACHE, NAME=bf_vnames[i])
		aE[i] = MrVar_Resample(oE, oTref, /CACHE, NAME=ef_vnames[i])
		aR[i] = MrVar_Resample(oPos, oTref)
		
		;Filter Data
		IF ~Array_Equal(fc, [0.0, 0.0]) THEN BEGIN
			;Filter Parameters
			fN = fs / 2.0
			f0 = fc[0] / fN
			f1 = fc[1] / fN
			A  = 75
			N  = Round(fs)*8.0
			
			aB[i] = (aB[i] -> Digital_Filter(f0, f1, A, N, /CACHE))[N:-N-1,*]
			aE[i] = (aE[i] -> Digital_Filter(f0, f1, A, N, /CACHE))[N:-N-1,*]
			
			IF i EQ 0 THEN oT_filt = oTref[N:-N-1]
			(aB[i]) -> SetData, oT_filt, aB[i]
			(aE[i]) -> SetData, oT_filt, aE[i]
			(aR[i]) -> SetData, oT_filt, (aR[i])[N:-N-1,*]
			
			aB[i] -> SetName, bf_vnames[i]
			aE[i] -> SetName, ef_vnames[i]
		ENDIF ELSE BEGIN
			oT_filt = oTref
		ENDELSE
	
	;-------------------------------------------
	; dB/dt ////////////////////////////////////
	;-------------------------------------------
		;dB/dt
		dB_dt  = ( Double((aB[i])['DATA', 1:-1, *]) - Double((aB[i])['DATA', 0:-2, *]) ) / dt
		odB_dt = MrVectorTS( oT_filt[1:-1], dB_dt )
	
		;Set properties
		odB_dt -> SetName, dB_dt_vnames[i]
		odB_dt -> Cache
		odB_dt['CATDESC'] = 'Time derivative of the magnetic field.'
		odB_dt['COLOR']   = ['Blue', 'Forest Green', 'Red']
		odB_dt['LABEL']   = '-' + partial + 'B$\down' + ['X', 'Y', 'Z'] + '$/' + partial + 't'
		odB_dt['TITLE']   = '-' + partial + 'B/' + partial + 't!C(nT/s)'
	
	;-------------------------------------------
	; Split into Components ////////////////////
	;-------------------------------------------
		;B-Field
		aB[i] -> Split, oBx, oBy, oBz, /CACHE
		oBx['COLOR']      = sc_colors[i]
		oBx['LABEL']      = sSC
		oBx['PLOT_TITLE'] = "Faraday's Law"
		oBx['TITLE']      = 'B$\downX$!C(nT)'
		oBx['UNITS']      = 'nT'
		
		oBy['COLOR'] = sc_colors[i]
;		oBy['LABEL'] = sSC
		oBy['TITLE'] = 'B$\downY$!C(nT)'
		oBy['UNITS'] = 'nT'
		
		oBz['COLOR'] = sc_colors[i]
;		oBz['LABEL'] = sSC
		oBz['TITLE'] = 'B$\downZ$!C(nT)'
		oBz['UNITS'] = 'nT'
		
		;E-Field
		aE[i] -> Split, oEx, oEy, oEz, /CACHE
		oEx['COLOR'] = sc_colors[i]
;		oEx['LABEL'] = sSC
		oEx['TITLE'] = 'E$\downX$!C(mV/m)'
		oEx['UNITS'] = 'mV/m'
		
		oEy['COLOR'] = sc_colors[i]
;		oEy['LABEL'] = sSC
		oEy['TITLE'] = 'E$\downY$!C(mV/m)'
		oEy['UNITS'] = 'mV/m'
		
		oEz['COLOR'] = sc_colors[i]
;		oEz['LABEL'] = sSC
		oEz['TITLE'] = 'E$\downZ$!C(mV/m)'
		oEz['UNITS'] = 'mV/m'
		
		;dB/dt
		odB_dt -> Split, odBx_dt, odBy_dt, odBz_dt, /CACHE
		odBx_dt['COLOR'] = sc_colors[i]
;		odBx_dt['LABEL'] = sSC
		odBx_dt['TITLE'] = '-' + partial + 'B$\downX$/' + partial + 't!C(nT/s)'
		odBx_dt['UNITS'] = 'nT/s'
		
		odBy_dt['COLOR'] = sc_colors[i]
;		odBy_dt['LABEL'] = sSC
		odBy_dt['TITLE'] = '-' + partial + 'B$\downY$/' + partial + 't!C(nT/s)'
		odBy_dt['UNITS'] = 'nT/s'
		
		odBz_dt['COLOR'] = sc_colors[i]
;		odBz_dt['LABEL'] = sSC
		odBz_dt['TITLE'] = '-' + partial + 'B$\downZ$/' + partial + 't!C(nT/s)'
		odBz_dt['UNITS'] = 'nT/s'
	
	;-------------------------------------------
	; Ranges ///////////////////////////////////
	;-------------------------------------------
		bxrange[0] <= oBx.min
		bxrange[1] >= oBx.max
		byrange[0] <= oBy.min
		byrange[1] >= oBy.max
		bzrange[0] <= oBz.min
		bzrange[1] >= oBz.max
		
		dbxdtrange[0] <= odBx_dt.min
		dbxdtrange[1] >= odBx_dt.max
		dbydtrange[0] <= odBy_dt.min
		dbydtrange[1] >= odBy_dt.max
		dbzdtrange[0] <= odBz_dt.min
		dbzdtrange[1] >= odBz_dt.max
		
		exrange[0] <= oEx.min
		exrange[1] >= oEx.max
		eyrange[0] <= oEy.min
		eyrange[1] >= oEy.max
		ezrange[0] <= oEz.min
		ezrange[1] >= oEz.max
	
	;-------------------------------------------
	; Store ////////////////////////////////////
	;-------------------------------------------
		aBx[i] = oBx
		aBy[i] = oBy
		aBz[i] = oBz
		aEx[i] = oEx
		aEy[i] = oEy
		aEz[i] = oEz
		adBx_dt[i] = odBx_dt
		adBy_dt[i] = odBy_dt
		adBz_dt[i] = odBz_dt
	ENDFOR
	
;-------------------------------------------
; Barycentric Averages /////////////////////
;-------------------------------------------
	;Bx
	oBx_bary = (aBx[0] + aBx[1] + aBx[2] + aBx[3]) / 4.0
	oBx_bary -> SetName, bx_bary_vname
	oBx_bary -> Cache
	oBx_bary['CATDESC'] = 'Barycentric average of the x-component of the magnetic field.'
	oBx_bary['COLOR']   = 'Magenta'
;	oBx_bary['LABEL']   = 'B$\downX,BC$'
	oBx_bary['LABEL']   = 'Barycenter'
	oBx_bary['TITLE']   = 'B$\downX,BC$!C(nT)'
	oBx_bary['UNITS']   = 'nT'
	
	;By
	oBy_bary = (aBy[0] + aBy[1] + aBy[2] + aBy[3]) / 4.0
	oBy_bary -> SetName, by_bary_vname
	oBy_bary -> Cache
	oBy_bary['CATDESC'] = 'Barycentric average of the y-component of the magnetic field.'
	oBy_bary['COLOR']   = 'Magenta'
;	oBy_bary['LABEL']   = 'B$\downY,BC$'
	oBy_bary['TITLE']   = 'B$\downY,BC$!C(nT)'
	oBy_bary['UNITS']   = 'nT'
	
	;Bz
	oBz_bary = (aBz[0] + aBz[1] + aBz[2] + aBz[3]) / 4.0
	oBz_bary -> SetName, bz_bary_vname
	oBz_bary -> Cache
	oBz_bary['CATDESC'] = 'Barycentric average of the z-component of the magnetic field.'
	oBz_bary['COLOR']   = 'Magenta'
;	oBz_bary['LABEL']   = 'B$\downZ,BC$'
	oBz_bary['TITLE']   = 'B$\downZ,BC$!C(nT)'
	oBz_bary['UNITS']   = 'nT'
	
	;Ex
	oEx_bary = (aEx[0] + aEx[1] + aEx[2] + aEx[3]) / 4.0
	oEx_bary -> SetName, ex_bary_vname
	oEx_bary -> Cache
	oEx_bary['CATDESC'] = 'Barycentric average of the x-component of the electric field.'
	oEx_bary['COLOR']   = 'Magenta'
;	oEx_bary['LABEL']   = 'E$\downX,BC$'
	oEx_bary['TITLE']   = 'E$\downX,BC$!C(mV/m)'
	oEx_bary['UNITS']   = 'mV/m'
	
	;Ey
	oEy_bary = (aEy[0] + aEy[1] + aEy[2] + aEy[3]) / 4.0
	oEy_bary -> SetName, ey_bary_vname
	oEy_bary -> Cache
	oEy_bary['CATDESC'] = 'Barycentric average of the y-component of the electric field.'
	oEy_bary['COLOR']   = 'Magenta'
;	oEy_bary['LABEL']   = 'E$\downY,BC$'
	oEy_bary['TITLE']   = 'E$\downY,BC$!C(mV/m)'
	oEy_bary['UNITS']   = 'mV/m'
	
	;Ez
	oEz_bary = (aEz[0] + aEz[1] + aEz[2] + aEz[3]) / 4.0
	oEz_bary -> SetName, ez_bary_vname
	oEz_bary -> Cache
	oEz_bary['CATDESC'] = 'Barycentric average of the time derivative of the  x-component of the electric field.'
	oEz_bary['COLOR']   = 'Magenta'
;	oEz_bary['LABEL']   = 'E$\downZ,BC$'
	oEz_bary['TITLE']   = 'E$\downZ,BC$!C(mV/m)'
	oEz_bary['UNITS']   = 'mV/m'
	
	;dBx/dt
	odBx_dt_bary = (adBx_dt[0] + adBx_dt[1] + adBx_dt[2] + adBx_dt[3]) / 4.0
	odBx_dt_bary -> SetName, dBx_dt_bary_vname
	odBx_dt_bary -> Cache
	odBx_dt_bary['CATDESC'] = 'Barycentric average of the time derivative of the  x-component of the magnetic field.'
	odBx_dt_bary['COLOR']   = 'Magenta'
;	odBx_dt_bary['LABEL']   = '(-' + partial + 'B$\downX$/' + partial + 't)$\downBC$'
	odBx_dt_bary['TITLE']   = '(-' + partial + 'B$\downX$/' + partial + 't)$\downBC$!C(nT/s)'
	odBx_dt_bary['UNITS']   = 'nT/s'
	
	;dBy/dt
	odBy_dt_bary = (adBy_dt[0] + adBy_dt[1] + adBy_dt[2] + adBy_dt[3]) / 4.0
	odBy_dt_bary -> SetName, dBy_dt_bary_vname
	odBy_dt_bary -> Cache
	odBy_dt_bary['CATDESC'] = 'Barycentric average of the time derivative of the  y-component of the magnetic field.'
	odBy_dt_bary['COLOR']   = 'Magenta'
;	odBy_dt_bary['LABEL']   = '(-' + partial + 'B$\downY$/' + partial + 't)$\downBC$'
	odBy_dt_bary['TITLE']   = '(-' + partial + 'B$\downY$/' + partial + 't)$\downBC$!C(nT/s)'
	odBy_dt_bary['UNITS']   = 'nT/s'
	
	;dBz/dt
	odBz_dt_bary = (adBz_dt[0] + adBz_dt[1] + adBz_dt[2] + adBz_dt[3]) / 4.0
	odBz_dt_bary -> SetName, dBz_dt_bary_vname
	odBz_dt_bary -> Cache
	odBz_dt_bary['CATDESC'] = 'Barycentric average of the time derivative of the z-component of the magnetic field.'
	odBz_dt_bary['COLOR']   = 'Magenta'
;	odBz_dt_bary['LABEL']   = '(-' + partial + 'B$\downZ$/' + partial + 't)$\downBC$'
	odBz_dt_bary['TITLE']   = '(-' + partial + 'B$\downZ$/' + partial + 't)$\downBC$!C(nT/s)'
	odBz_dt_bary['UNITS']   = 'nT/s'
	
	;|dB/dt|
	odB_dt_bary_mag = odBx_dt_bary^2.0 + odBy_dt_bary^2.0 + odBz_dt_bary^2.0
	odB_dt_bary_mag -> SetData, Sqrt(odB_dt_bary_mag['DATA'])
	odB_dt_bary_mag -> SetName, db_dt_bary_mag_vname
	odB_dt_bary_mag -> Cache
	odB_dt_bary_mag['CATDESC'] = 'Magnitude of dB/dt'
	odB_dt_bary_mag['COLOR']   = 'Red'
	odB_dt_bary_mag['LABEL']   = '|'+partial+'B/'+partial+'t|'
	odB_dt_bary_mag['TITLE']   = '|'+partial+'B/'+partial+'t|!C(nT/s)'
	odB_dt_bary_mag['UNITS']   = 'nT/s'
	
;-------------------------------------------
; Curl(E) //////////////////////////////////
;-------------------------------------------
	oRecipVec = MrVar_RecipVec(aR[0], aR[1], aR[2], aR[3])
	
	;Curl of E
	;   - 1e3 converts to nV/m^2 = kg/Cs^2 * 1e6
	oCurlE =  1e3 * oRecipVec -> Curl( aE[0], aE[1], aE[2], aE[3] )
	oCurlE -> SetName, curlE_vname
	oCurlE -> Cache
	oCurlE['COLOR'] = ['Blue', 'Forest Green', 'Red']
	oCurlE['LABEL'] = '(' + nabla + 'xE)$\down' + ['X', 'Y', 'Z'] + '$'
	oCurlE['TITLE'] = nabla + 'xE!C(nV/m^2)'
	
	;Attributes
	oCurlE['CATDESC'] = 'Curl(E)'
	oCurlE['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oCurlE['LABEL']   = '(' + nabla + 'xE)$\down' + ['X', 'Y', 'Z'] + '$'
	oCurlE['TITLE']   = nabla + 'xE!C(nT/s)'
	oCurlE['UNITS']   = ''
	
	;Split Curl(E) into components
	oCurlE -> Split, oCurlEx, oCurlEy, oCurlEz
	
	;Curl(E)x
	oCurlEx -> SetName, curlEx_vname
	oCurlEx -> Cache
	oCurlEx['AXIS_RANGE'] = [oCurlE.min, oCurlE.max]
	oCurlEx['COLOR']     = 'Blue'
	oCurlEx['LABEL']     = '(' + nabla + 'xE)$\downX$'
	oCurlEx['TITLE']     = '(' + nabla + 'xE)!C(nT/s)'
	
	;Curl(E)y
	oCurlEy -> SetName, curlEy_vname
	oCurlEy -> Cache
	oCurlEy['AXIS_RANGE'] = oCurlEx['AXIS_RANGE']
	oCurlEy['COLOR']      = 'Forest Green'
	oCurlEy['LABEL']      = '(' + nabla + 'xE)$\downY$'
	oCurlEy['TITLE']      = '(' + nabla + 'xE)!C(nT/s)'
	
	;Curl(E)z
	oCurlEz -> SetName, curlEz_vname
	oCurlEz -> Cache
	oCurlEy['AXIS_RANGE'] = oCurlEx['AXIS_RANGE']
	oCurlEz['COLOR']      = 'Red'
	oCurlEz['LABEL']      = '(' + nabla + 'xE)$\downZ$'
	oCurlEz['TITLE']      = '(' + nabla + 'xE)!C(nV/m$\up2$)'
	
	;Magnitude
	oCurlEmag = oCurlE -> Magnitude(/CACHE, NAME=curle_mag_vname)
	oCurlEmag['CATDESC'] = 'Magnitude of the curl of the electric field.'
	oCurlEmag['COLOR']   = 'Black'
	oCurlEmag['LABEL']   = '|'+nabla+'xE|'
	oCurlEmag['TITLE']   = '$\Delta$E!C(nV/m$\up2$)'
	oCurlEmag['UNITS']   = 'nV/m^2'
	
;-------------------------------------------
; Div(E) ///////////////////////////////////
;-------------------------------------------
	;Charge Density
	;   - 1e3 converts to nV/m^2
	oDivE = 1e3 * oRecipVec -> Divergence(aE[0], aE[1], aE[2], aE[3])
	oDivE -> SetName, divE_vname
	oDivE -> Cache
	oDivE['CATDESC']    = 'Divergence of the electric field.'
	oDivE['COLOR']      = 'Blue'
	oDivE['LABEL']      = nabla+dot+'E'
	oDivE['TITLE']      = nabla+dot+'E!C(nV/m$\up2$)'
	oDivE['UNITS']      = 'nV/m^2'
	
	;Absolute Value
	;   - For error computations later
	oAbsDivE = oDivE -> Copy(dive_mag_vname, /CACHE)
	oAbsDivE -> SetData, Abs(oAbsDivE['DATA'])
	oAbsDivE['AXIS_RANGE'] = [oAbsDivE.min < oCurlEmag.min < odB_dt_bary_mag.min, $
	                          oAbsDivE.max > oCurlEmag.max > odB_dt_bary_mag.max]
	oAbsDivE['CATDESC']    = 'Absolute value of the electric field divergence.'
	oAbsDivE['COLOR']      = 'Blue'
	oAbsDivE['LABEL']      = '|'+nabla+dot+'E|'
	oAbsDivE['TITLE']      = '|'+nabla+dot+'E|!C(nV/m$\up2$)'
	oAbsDivE['UNITS']      = 'nV/m^2'
	
	;Standardize axis ranges
	odB_dt_bary_mag['AXIS_RANGE'] = oAbsDivE['AXIS_RANGE']
	oCurlEmag['AXIS_RANGE']       = oAbsDivE['AXIS_RANGE']
	
;-------------------------------------------
; Errors in Curl(E) ////////////////////////
;-------------------------------------------
	;Mean error among four spacecraft
	IF 0B THEN BEGIN
		dEx = ( (aE[0])['DELTA_PLUS_VAR'] + $
		        (aE[1])['DELTA_PLUS_VAR'] + $
		        (aE[2])['DELTA_PLUS_VAR'] + $
		        (aE[3])['DELTA_PLUS_VAR'] ) / 4.0
		
		;Mean separation among spacecraft
		aR   = MrVar_Get(r_vnames)
		oR12 = (aR[1] - aR[0]) -> Magnitude()
		oR13 = (aR[2] - aR[0]) -> Magnitude()
		oR14 = (aR[3] - aR[0]) -> Magnitude()
		oR23 = (aR[2] - aR[1]) -> Magnitude()
		oR24 = (aR[3] - aR[1]) -> Magnitude()
		oR34 = (aR[3] - aR[2]) -> Magnitude()
		
		;Mean position
		r12 = Reform( Mean(oR12['DATA']) )
		r13 = Reform( Mean(oR13['DATA']) )
		r14 = Reform( Mean(oR14['DATA']) )
		r23 = Reform( Mean(oR23['DATA']) )
		r24 = Reform( Mean(oR24['DATA']) )
		r34 = Reform( Mean(oR34['DATA']) )
		r   = Mean( [r12, r13, r14, r23, r24, r34] )
		
		;Error of Curl(E)
		odCurlE = Sqrt(4.0/3.0) / r * dEx
		odCurlE['CATDESC'] = 'Standard deviation of the Curl of E.'
		odCurlE['TITLE']   = '$\sigma$$\downE$!C(nV/m$\up2$)'
		odCurlE['UNITS']   = 'nV/m^2'
		
		;Set as error bars
		oCurlEx['DELTA_PLUS_VAR']  = odCurlE
		oCurlEx['DELTA_MINUS_VAR'] = odCurlE
	ENDIF
	
;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	(aBx[0])['AXIS_RANGE'] = bxrange
	(aBy[0])['AXIS_RANGE'] = byrange
	(aBz[0])['AXIS_RANGE'] = bzrange
	
	(adBx_dt[0])['AXIS_RANGE'] = dbxdtrange
	(adBy_dt[0])['AXIS_RANGE'] = dbydtrange
	(adBz_dt[0])['AXIS_RANGE'] = dbzdtrange
	
	(aEx[0])['AXIS_RANGE'] = exrange
	(aEy[0])['AXIS_RANGE'] = eyrange
	(aEz[0])['AXIS_RANGE'] = ezrange

	;Plot data
	win = MrVar_PlotTS( [ bx_vnames[0], by_vnames[0], bz_vnames[0], $
	                      dbx_dt_vnames[0], dby_dt_vnames[0], dbz_dt_vnames[0], $
	                      ex_vnames[0], ey_vnames[0], ez_vnames[0], $
	                      curlex_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	win = MrVar_OPlotTS( bx_vnames[0], [bx_vnames[1:*], bx_bary_vname] )
	win = MrVar_OPlotTS( by_vnames[0], [by_vnames[1:*], by_bary_vname] )
	win = MrVar_OPlotTS( bz_vnames[0], [bz_vnames[1:*], bz_bary_vname] )
	win = MrVar_OPlotTS( ex_vnames[0], [ex_vnames[1:*], ex_bary_vname] )
	win = MrVar_OPlotTS( ey_vnames[0], [ey_vnames[1:*], ey_bary_vname] )
	win = MrVar_OPlotTS( ez_vnames[0], [ez_vnames[1:*], ez_bary_vname] )
	win = MrVar_OPlotTS( dbx_dt_vnames[0], [dbx_dt_vnames[1:*], dbx_dt_bary_vname] )
	win = MrVar_OPlotTS( dby_dt_vnames[0], [dby_dt_vnames[1:*], dby_dt_bary_vname] )
	win = MrVar_OPlotTS( dbz_dt_vnames[0], [dbz_dt_vnames[1:*], dbz_dt_bary_vname] )
	win = MrVar_OPlotTS( curlex_vname, curley_vname )
	win = MrVar_OPlotTS( curlex_vname, curlez_vname )
	
	;Pretty-up the window
	win.name = 'FaradaysLaw'
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 13]
	win    -> Refresh

;-------------------------------------------
; Curl(E) vs. -dB/dt ///////////////////////
;-------------------------------------------
	;Scatter plots of -dB/dt and Curl(E) to view gain and offset parameters
	odBxdt = MrVar_Get(dbx_dt_bary_vname)
	odBydt = MrVar_Get(dby_dt_bary_vname)
	odBzdt = MrVar_Get(dbz_dt_bary_vname)
	oCurlEx = MrVar_Get(curlEx_vname)
	oCurlEy = MrVar_Get(curlEy_vname)
	oCurlEz = MrVar_Get(curlEz_vname)
	
	;Data range
;	xrange = [ Min( [odBxdt[i0:i1].min, oCurlEx[i0+1:i1].min] ), Max( [odBxdt[i0:i1].max, oCurlEx[i0+1:i1].max] ) ]
;	yrange = [ Min( [odBydt[i0:i1].min, oCurlEy[i0+1:i1].min] ), Max( [odBydt[i0:i1].max, oCurlEy[i0+1:i1].max] ) ]
;	zrange = [ Min( [odBzdt[i0:i1].min, oCurlEz[i0+1:i1].min] ), Max( [odBzdt[i0:i1].max, oCurlEz[i0+1:i1].max] ) ]
	
	w2 = MrWindow( ASPECT  = 1.0, $
	               LAYOUT  = [3,1], $
	               NAME    = 'CurlEvdBdt', $
	               REFRESH = 0, $
	               XSIZE   = 1000 )
	
	;Scatter plots
	p4 = MrVar_Plot( odBxdt, oCurlEx[1:-1], $
	                 COLOR         = 'Black', $
	                 TITLE         = '', $
;	                 XRANGE        = xrange, $
	                 XTITLE        = '-'+partial+'B$\downX$/'+partial+'t (nT/s)', $
;	                 YRANGE        = xrange, $
	                 YTITLE        = '('+nabla+'xE)$\downX$ (nV/m$\up2$)', $
	                 /CURRENT )
	p5 = MrVar_Plot( odBydt, oCurlEy[1:-1], $
	                 COLOR  = 'Black', $
;	                 XRANGE = yrange, $
	                 XTITLE = '-'+partial+'B$\downY$/'+partial+'t (nT/s)', $
;	                 YRANGE = yrange, $
	                 YTITLE = '('+nabla+'xE)$\downY$ (nV/m$\up2$)', $
	                 /CURRENT )
	p6 = MrVar_Plot( odBzdt, oCurlEz[1:-1], $
	                 COLOR  = 'Black', $
	                 TITLE  = '', $
;	                 XRANGE = zrange, $
	                 XTITLE = '-'+partial+'B$\downZ$/'+partial+'t (nT/s)', $
;	                 YRANGE = zrange, $
	                 YTITLE = '('+nabla+'xE)$\downZ$ (nV/m$\up2$)', $
	                 /CURRENT )
	
	;Fits to data
	p4 -> GetData, x, y
	ifinite = where(finite(x) and finite(y))
	isort   = sort(x[ifinite])
	fit4    = LADFit(x[ifinite[isort]], y[ifinite[isort]])
	r4      = Correlate(x[ifinite[isort]], y[ifinite[isort]])
	xrange  = [ Min(x, MAX=xmax) < Min(y, MAX=ymax), xmax > ymax ]
	txt4    = String(fit4[1], fit4[0], r4, FORMAT='(%"m=%0.4e!Cb=%0.4e!CR=%0.4f")')
	lin4    = fit4[0] + fit4[1]*xrange
	
	p5 -> GetData, x, y
	ifinite = where(finite(x) and finite(y))
	isort   = sort(x[ifinite])
	fit5    = LADFit(x[ifinite[isort]], y[ifinite[isort]])
	r5      = Correlate(x[ifinite[isort]], y[ifinite[isort]])
	yrange  = [ Min(x, MAX=xmax) < Min(y, MAX=ymax), xmax > ymax ]
	txt5    = String(fit5[1], fit5[0], r5, FORMAT='(%"m=%0.4e!Cb=%0.4e!CR=%0.4f")')
	lin5    = fit5[0] + fit5[1]*yrange
	
	p6 -> GetData, x, y
	ifinite = where(finite(x) and finite(y))
	isort   = sort(x[ifinite])
	fit6    = LADFit(x[ifinite[isort]], y[ifinite[isort]])
	r6      = Correlate(x[ifinite[isort]], y[ifinite[isort]])
	zrange  = [ Min(x, MAX=xmax) < Min(y, MAX=ymax), xmax > ymax ]
	txt6    = String(fit6[1], fit6[0], r6, FORMAT='(%"m=%0.4e!Cb=%0.4e!CR=%0.4f")')
	lin6    = fit6[0] + fit6[1]*zrange
	
	;Text
	txt4 = MrText( 0.05, 0.9, txt4, $
	               CHARSIZE = 1.0, $
	               NAME     = 'Txt: Eqn-X', $
	               /RELATIVE, $
	               TARGET   = p4 )
	
	txt5 = MrText( 0.05, 0.9, txt5, $
	               CHARSIZE = 1.0, $
	               NAME     = 'Txt: Eqn-Y', $
	               /RELATIVE, $
	               TARGET   = p5 )
	
	txt6 = MrText( 0.05, 0.9, txt6, $
	               CHARSIZE = 1.0, $
	               NAME     = 'Txt: Eqn-Z', $
	               /RELATIVE, $
	               TARGET   = p6 )
	
	;Draw Line Fits
	ps4 = MrPlotS( xrange, lin4, $
	               COLOR     = 'Red', $
	               LINESTYLE = '--', $
	               NAME      = 'PS: X-Fit', $
	               TARGET    = p4 )
	
	ps5 = MrPlotS( yrange, lin5, $
	               COLOR     = 'Red', $
	               LINESTYLE = '--', $
	               NAME      = 'PS: Y-Fit', $
	               TARGET    = p5 )
	
	ps6 = MrPlotS( zrange, lin6, $
	               COLOR     = 'Red', $
	               LINESTYLE = '--', $
	               NAME      = 'PS: Z-Fit', $
	               TARGET    = p6 )
	
	;Make Pretty
	p4 -> SetLayout, [1,1]
;	p4 -> SetProperty, XRANGE=xrange, YRANGE=xrange
;	p5 -> SetProperty, XRANGE=yrange, YRANGE=yrange
;	p6 -> SetProperty, XRANGE=zrange, YRANGE=zrange
	trange   = MrVar_GetTRange()
	p5.title = StrJoin(StrSplit(StrMid(trange[0], 0, 19), 'T', /EXTRACT), ' ') + ' -- ' + StrMid(trange[1], 11, 8)
	w2 -> TrimLayout
	w2 -> Remove, w2 -> Get(/ALL, ISA='MrLegend')
	w2 -> Refresh
;	w2 -> Save, '/home/argall/figures/20151206/mms_edp_brst_l2_4sc-e-maxwell-curlE-scatter_20151206_233827_233836.png'

;-------------------------------------------
; Corrected Curl(E) with -dB/dt ////////////
;-------------------------------------------
	;Correct dB/dt by dividing by the slope of the best fit line.
	;   - uV/m^2 = uT/s
	oCurlEx_cor = oCurlEx / Abs(fit4[1])
	oCurlEy_cor = oCurlEy / Abs(fit5[1])
	oCurlEz_cor = oCurlEz / Abs(fit6[1])
;	odBxdt_cor = odBxdt * fit4[1]
;	odBydt_cor = odBydt * fit5[1]
;	odBzdt_cor = odBzdt * fit6[1]
	
	xrange = [oCurlEx_cor.min < odBx_dt_bary.min, oCurlEx_cor.Max > odBx_dt_bary.max]
	yrange = [oCurlEy_cor.min < odBy_dt_bary.min, oCurlEy_cor.Max > odBy_dt_bary.max]
	zrange = [oCurlEz_cor.min < odBz_dt_bary.min, oCurlEz_cor.Max > odBz_dt_bary.max]
;	xrange = [oCurlEx.min < odBxdt_cor.min, oCurlEx.Max > odBxdt_cor.max]
;	yrange = [oCurlEy.min < odBydt_cor.min, oCurlEy.Max > odBydt_cor.max]
;	zrange = [oCurlEz.min < odBzdt_cor.min, oCurlEz.Max > odBzdt_cor.max]
	
	oCurlEx_cor -> SetName, curlex_corr_vname
	oCurlEx_cor -> Cache
	oCurlEx_cor['AXIS_RANGE'] = xrange
	oCurlEx_cor['COLOR']      = 'Black'
	oCurlEx_cor['LABEL']      = nabla + 'xE'
	oCurlEx_cor['PLOT_TITLE'] = 'Corrected '+nabla+'xE vs. -' + partial+'B/'+partial+'t'
	oCurlEx_cor['TITLE']      = fancyf+'$\downX$!C(nV/m$\up2$)'
	oCurlEz_cor['UNITS']      = 'nV/m^2'
	
	oCurlEy_cor -> SetName, curley_corr_vname
	oCurlEy_cor -> Cache
	oCurlEy_cor['AXIS_RANGE'] = yrange
	oCurlEy_cor['COLOR']      = 'Black'
	oCurlEy_cor['LABEL']      = nabla + 'xE'
	oCurlEy_cor['TITLE']      = fancyf+'$\downY$!C(nV/m$\up2$)'
	oCurlEy_cor['UNITS']      = 'nV/m^2'
	
	oCurlEz_cor -> SetName, curlez_corr_vname
	oCurlEz_cor -> Cache
	oCurlEz_cor['AXIS_RANGE'] = zrange
	oCurlEz_cor['COLOR']      = 'Black'
	oCurlEz_cor['LABEL']      = nabla + 'xE'
	oCurlEz_cor['TITLE']      = fancyf+'$\downZ$!C(nV/m$\up2$)'
	oCurlEz_cor['UNITS']      = 'nV/m^2'
	
	odBx_dt_bary['COLOR'] = 'Blue'
	odBx_dt_bary['LABEL'] = '-'+partial+'B$\downX$/'+partial+'t'
	
	odBy_dt_bary['COLOR'] = 'Blue'
	odBy_dt_bary['LABEL'] = '-'+partial+'B/'+partial+'t'
	
	odBz_dt_bary['COLOR'] = 'Blue'
	odBz_dt_bary['LABEL'] = '-'+partial+'B/'+partial+'t'
	
;	oCurlEx['AXIS_RANGE'] = xrange
;	oCurlEx['COLOR']      = 'Black'
;	oCurlEx['LABEL']      = nabla + 'xE'
;	oCurlEx['TITLE']      = fancyf+'!C(nV/m$\up2$)'
	
;	oCurlEy['AXIS_RANGE'] = yrange
;	oCurlEy['COLOR']      = 'Black'
;	oCurlEy['LABEL']      = nabla + 'xE'
;	oCurlEy['TITLE']      = fancyf+'!C(nV/m$\up2$)'
	
;	oCurlEz['AXIS_RANGE'] = zrange
;	oCurlEz['COLOR']      = 'Black'
;	oCurlEz['LABEL']      = nabla + 'xE'
;	oCurlEz['TITLE']      = fancyf+'!C(nV/m$\up2$)'
	
;	odBxdt_cor -> SetName, dbx_dt_bary_corr_vname
;	odBxdt_cor -> Cache
;	odBxdt_cor['CATDESC'] = 'dBx/dt at the barycenter, corrected by slope of correlation with Curl(E).'
;	odBxdt_cor['COLOR']   = 'Blue'
;	odBxdt_cor['LABEL']   = '-'+partial+'B$\downX$/'+partial+'t'
;	odBxdt_cor['TITLE']   = '-'+partial+'B$\downX$/'+partial+'t!CnT/s'
;	odBxdt_cor['UNITS']   = 'nT/s'
	
;	odBydt_cor -> SetName, dby_dt_bary_corr_vname
;	odBydt_cor -> Cache
;	odBydt_cor['CATDESC'] = 'dBy/dt at the barycenter, corrected by slope of correlation with Curl(E).'
;	odBydt_cor['COLOR']   = 'Blue'
;	odBydt_cor['LABEL']   = '-'+partial+'B$\downY$/'+partial+'t'
;	odBydt_cor['TITLE']   = '-'+partial+'B$\downY$/'+partial+'t!CnT/s'
;	odBydt_cor['UNITS']   = 'nT/s'
	
;	odBzdt_cor -> SetName, dbz_dt_bary_corr_vname
;	odBzdt_cor -> Cache
;	odBzdt_cor['CATDESC'] = 'dBz/dt at the barycenter, corrected by slope of correlation with Curl(E).'
;	odBzdt_cor['COLOR']   = 'Blue'
;	odBzdt_cor['LABEL']   = '-'+partial+'B$\downZ$/'+partial+'t'
;	odBzdt_cor['TITLE']   = '-'+partial+'B$\downZ$/'+partial+'t!CnT/s'
;	odBzdt_cor['UNITS']   = 'nT/s'
	
	;Plot the data
	w3 = MrVar_PlotTS( [curlEx_corr_vname, curlEy_corr_vname, curlEz_corr_vname], $
	                   /NO_REFRESH )
	
	w3 = MrVar_OPlotTS( curlEx_corr_vname, dbx_dt_bary_vname )
	w3 = MrVar_OPlotTS( curlEy_corr_vname, dby_dt_bary_vname )
	w3 = MrVar_OPlotTS( curlEz_corr_vname, dbz_dt_bary_vname )
	
	w3.name = 'Faraday-Corrected'
	w3[0] -> SetLayout, [1,1]
	w3    -> TrimLayout
	w3.oxmargin = [13, 10]
	w3    -> Refresh
;-------------------------------------------
; dB/dt, CurlE, DivE ///////////////////////
;-------------------------------------------
	oErr = oAbsDivE[1:-1] / (oCurlE[1:-1,*] - [[odBx_dt_bary['DATA']], [odBy_dt_bary['DATA']], [odBz_dt_bary['DATA']]]) -> Magnitude()
	oErr -> SetName, dive_err_vname
	oErr -> Cache
	oErr['CATDESC'] = 'Quality estimate for Div(E) computed as |Div(E)| / |Curl(E) - (dB/dt)|'
	oErr['TITLE']   = '$\sigma$$\down'+nabla+'.E$'
	
	w4 = MrVar_PlotTS( [curle_mag_vname, dive_err_vname], /NO_REFRESH )
	w4 = MrVar_OPlotTS( curle_mag_vname, [dive_mag_vname, db_dt_bary_mag_vname] )
	w4.name = 'Faraday-Curl-Div-E'
	w4[0] -> SetLayout, [1,1]
	w4.oxmargin = [12, 9]
	w4    -> TrimLayout
	w4    -> Refresh

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
		
		suffix = String(fc, FORMAT='(%"-fc=%0.1f-%0.1f")')
		
		;File name
		fname1 = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-faraday'], '_' ) + suffix
		fname2 = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-faraday-scatter'], '_' ) + suffix
		fname3 = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-faraday-corrected'], '_' ) + suffix
		fname4 = StrJoin( ['mms', edp_instr, edp_mode, level, '4sc-faraday-curl-div-e'], '_' ) + suffix
		
		fname1 = FilePath( fname1, ROOT_DIR=output_dir )
		fname2 = FilePath( fname2, ROOT_DIR=output_dir )
		fname3 = FilePath( fname3, ROOT_DIR=output_dir )
		fname4 = FilePath( fname4, ROOT_DIR=output_dir )
		
		;Save the figure
		fout1 = MrVar_PlotTS_Save( win, fname1, output_ext )
		fout2 = MrVar_PlotTS_Save( w2,  fname2, output_ext )
		fout3 = MrVar_PlotTS_Save( w3,  fname3, output_ext )
		fout4 = MrVar_PlotTS_Save( w4,  fname3, output_ext )
	ENDIF


;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END