; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_EGradV
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
;   Generate a plot to provide an overview of reconnection quantities:
;       1. B barycentric average
;       2. Curl(E)x & -dBx/dt
;       3. Curl(E)y & -dBy/dt
;       4. Curl(E)z & -dBz/dt
;       5. Charge density: e0 * Div(E)
;       6. Ex barycenter & -Grad(V)x
;       7. Ey barycenter & -Grad(V)y
;       8. Ez barycenter & -Grad(V)z
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
;       2018/02/04  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_EGradV, mode, $
FC = fc, $
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
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	tf_load   = ~Keyword_Set(no_load)
	fgm_instr = mode EQ 'brst' ? 'fsm' : 'fgm'
	IF N_Elements(coords) EQ 0 THEN coords = 'gse'
	IF N_Elements(fc)     EQ 0 THEN fc     = 0.0
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
	nabla     = '!9'+String(71B)+'!X'
	partial   = '!9'+String(68B)+'!X'
	
;-------------------------------------------
; Variable Parameters //////////////////////
;-------------------------------------------
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
	e_vnames  = sc + '_' + StrJoin([edp_instr, 'dce',  edp_coords, edp_mode, level], '_' )
	ex_vnames = e_vnames + '_x'
	ey_vnames = e_vnames + '_y'
	ez_vnames = e_vnames + '_z'
	v_vnames  = sc + '_' + StrJoin([edp_instr, 'scpot', edp_mode, level], '_' )
	v_srvy_vnames = sc + '_' + StrJoin([edp_instr, 'scpot', 'fast', level], '_' )
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames    = sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames    = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;Output names
	ef_vnames     = sc + '_' + StrJoin([edp_instr, 'dce',  edp_coords, 'lpfilt', edp_mode, level], '_' )
	v_pca_vnames  = sc + '_' + StrJoin([edp_instr, 'scpot', 'pca',               edp_mode, level], '_' )
	v_out_vnames  = sc + '_' + StrJoin([edp_instr, 'scpot', 'out',               edp_mode, level], '_' )
	e_bary_vname  = StrJoin( ['mms', edp_instr, 'e', edp_coords, 'bary', edp_mode, level], '_' )
	ex_bary_vname = e_bary_vname + '_x'
	ey_bary_vname = e_bary_vname + '_y'
	ez_bary_vname = e_bary_vname + '_z'
	gradv_vname   = StrJoin( ['mms', 'dce', 'gradv', edp_mode, level], '_' )
	gradvx_vname  = gradv_vname + '_x'
	gradvy_vname  = gradv_vname + '_y'
	gradvz_vname  = gradv_vname + '_z'
	gradv_corr_vname  = StrJoin( ['mms', 'dce', 'gradv', 'corrected', edp_mode, level], '_' )
	gradvx_corr_vname = gradv_corr_vname + '_x'
	gradvy_corr_vname = gradv_corr_vname + '_y'
	gradvz_corr_vname = gradv_corr_vname + '_z'
		
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;E-Field
		MrMMS_Load_Data, '', edp_instr, edp_mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = '*dce_'+coords+'*'
		
		;Spacecraft Potential
		IF 1 THEN BEGIN
			MrMMS_Load_Data, '', edp_instr, edp_mode, level, $
			                 OPTDESC   = 'scpot', $
			                 VARFORMAT = '*scpot*'
		ENDIF ELSE BEGIN
			file = '/home/argall/data/mms_Vns_fast_20170711.tplot'
			MrVar_ImportFromSpedas, '*', file, VARNAMES=vnames, /VERBOSE
			MrVar_TLimit, vnames
			FOR i = 0, N_Elements(vnames) - 1 DO BEGIN
				oVar = MrVar_Get(vnames[i])
				oVar -> SetName, v_vnames[i]
			ENDFOR
		ENDELSE
		
		;Survey data
;		tr_in  = MrVar_GetTRange()
;		tr_ssm = MrVar_GetTRange('SSM')
;		dt_ssm = tr_ssm[1] - tr_ssm[0]
;		IF dt_ssm LT 20.0*18.0 THEN BEGIN
;			;Adjust time interval
;			dt      = (20.0*18.0 - dt_ssm) / 2.0
;			tr_srvy = tr_ssm + [-dt, dt]
;			MrVar_SetTRange, tr_srvy, 'SSM', tr_in[0]
;			
;			;Load data
;			MrMMS_Load_Data, '', edp_instr, 'fast', level, $
;			                 OPTDESC   = 'scpot', $
;			                 VARFORMAT = '*scpot*'
;			
;			;Restore original time interval
;			MrVar_SetTRange, tr_in
;		ENDIF
		
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
; Spin Epoch ///////////////////////////////
;-------------------------------------------
;	oPCA = ObjArr(nSC)
;	FOR i = 0, nSC - 1 DO BEGIN
;		sc = 'mms' + String(i+1, FORMAT='(i1)')
;		oPCA[i] = MrMMS_DSS_PCA(sc, v_srvy_vnames[i], NO_LOAD=~tf_load)
;	ENDFOR
	
;-------------------------------------------
; Interpolate //////////////////////////////
;-------------------------------------------
	;Interpolate everything to E
	oE1   = MrVar_Get(e_vnames[0])
	oTref = oE1['TIMEVAR']
	dt    = oTref -> GetSI(RATE=fs)
	
	;ranges
	exrange = [!Values.f_infinity, -!Values.f_infinity]
	eyrange = [!Values.f_infinity, -!Values.f_infinity]
	ezrange = [!Values.f_infinity, -!Values.f_infinity]
	vrange  = [!Values.f_infinity, -!Values.f_infinity]
	
	oE   = ObjArr(nSC)
	oEx  = ObjArr(nSC)
	oEy  = ObjArr(nSC)
	oEz  = ObjArr(nSC)
	oV   = ObjArr(nSC)
	oPos = ObjArr(nSC)
	FOR i = 0, nSC - 1 DO BEGIN
		oE[i]   = MrVar_Resample(e_vnames[i], oTref)
		oV[i]   = MrVar_Resample(v_vnames[i], oTref)
		oPos[i] = MrVar_Resample(r_vnames[i], oTref)
;		oPCA[i] = MrVar_Resample(oPCA[i], oTref)
		
		;Remove principal component
;		oV[i] = oV[i] - oPCA[i]
	
	;-------------------------------------------
	; Filter Data //////////////////////////////
	;-------------------------------------------
		IF fc GT 0 THEN BEGIN
			;Filter
			fN = fs / 2.0
			f0 = 0.0
			f1 = fc / fN
			A  = 50
			N  = Round(fs)
			oE[i] = oE[i] -> Digital_Filter(f0, f1, A, N)
			oV[i] = oV[i] -> Digital_Filter(f0, f1, A, N)
		ENDIF
	
	;-------------------------------------------
	; Break into Components ////////////////////
	;-------------------------------------------
		oE[i] -> Split, Ex, Ey, Ez
		
		;Ex
		Ex -> SetName, ex_vnames[i]
		Ex -> Cache
		Ex['PLOT_TITLE'] = 'E = -'+nabla+'V'
		Ex['COLOR']      = sc_colors[i]
		Ex['LABEL']      = 'MMS' + String(i+1, FORMAT='(i1)')
		Ex['TITLE']      = 'E$\downX$!C(mV/m)'
		
		;Ey
		Ey -> SetName, ey_vnames[i]
		Ey -> Cache
		Ey['COLOR'] = sc_colors[i]
		Ey['TITLE'] = 'E$\downY$!C(mV/m)'
		
		;Ez
		Ez -> SetName, ez_vnames[i]
		Ez -> Cache
		Ez['COLOR'] = sc_colors[i]
		Ez['TITLE'] = 'E$\downZ$!C(mV/m)'
		
		;V
		(oV[i]) -> SetName, v_out_vnames[i]
		(oV[i]) -> Cache
		(oV[i])['COLOR'] = sc_colors[i]
		(oV[i])['TITLE'] = 'V$\downSC$!C(V)'
		
		exrange[0] <= Ex.min
		exrange[1] >= Ex.max
		eyrange[0] <= Ey.min
		eyrange[1] >= Ey.max
		ezrange[0] <= Ez.min
		ezrange[1] >= Ez.max
		vrange[0]  <= (oV[i]).min
		vrange[1]  >= (oV[i]).max
	
	;-------------------------------------------
	; Store Components /////////////////////////
	;-------------------------------------------
		oEx[i] = Ex
		oEy[i] = Ey
		oEz[i] = Ez
	ENDFOR
	
	(oEx[0])['AXIS_RANGE'] = exrange
	(oEy[0])['AXIS_RANGE'] = eyrange
	(oEz[0])['AXIS_RANGE'] = ezrange
	(oV[0])['AXIS_RANGE']  = vrange

	
;-------------------------------------------
; Barycentric Averages /////////////////////
;-------------------------------------------
	;E
	oE_bary = (oE[0] + oE[1] + oE[2] + oE[3]) / 4.0
	oE_bary -> SetName, e_bary_vname
	oE_bary -> Cache
	
	;Split into components
	oE_bary -> Split, oEx_bary, oEy_bary, oEz_bary, /CACHE
	oE_bary['COLOR']      = ['Blue', 'Forest Green', 'Red']
	oE_bary['LABEL']      = 'E$\down' + ['X', 'Y', 'Z'] + '$'
	oE_bary['TITLE']      = 'E$\downBC$!C(mV/m^2)'
	
	;Ex
	oEx_bary['COLOR'] = 'Magenta'
	oEx_bary['LABEL'] = 'E$\downBC$'
	oEx_bary['TITLE'] = 'E!C(mV/m)'
	
	;Ey
	oEy_bary['COLOR'] = 'Magenta'
;	oEy_bary['LABEL'] = 'E$\downBC,Y$'
	oEy_bary['TITLE'] = 'E!C(mV/m)'
	
	;Ez
	oEz_bary['COLOR'] = 'Magenta'
;	oEz_bary['LABEL'] = 'E$\downBC,Z$'
	oEz_bary['TITLE'] = 'E!C(mV/m)'
	
;-------------------------------------------
; -Grad(V) /////////////////////////////////
;-------------------------------------------
	oRecipVec = MrVar_RecipVec(oPos[0], oPos[1], oPos[2], oPos[3])
	
	
	;E = -Grad(V)
	;   - 1e0 converts V/km to mV/m
	oGradV = -oRecipVec -> Gradient( oV[0], oV[1], oV[2], oV[3] )
	oGradV -> SetName, gradv_vname
	oGradV -> Cache
	oGradV['CATDESC'] = '-Grad(V)'
	oGradV['COLOR']   = ['Blue', 'Green', 'Red']
	oGradV['LABEL']   = ['X', 'Y', 'Z']
	oGradV['TITLE']   = '-Grad(V)!C(mV/m)'
	oGradV['UNITS']   = 'mV/m'
	
	;Split into components
	oGradV -> Split, oGradVx, oGradVy, oGradVz, /CACHE
	
;	oGradVx['COLOR']     = 'Magenta'
;	oGradVx['LABEL']     = '(' + nabla + 'V)$\downX$'
;	oGradVx['LINESTYLE'] = '--'
	oGradVx['TITLE']     = '(-'+nabla+'V)$\downX$!C(mV/m)'
	
;	oGradVy['COLOR']     = 'Magenta'
;	oGradVy['LABEL']     = '(' + nabla + 'V)$\downY$'
;	oGradVy['LINESTYLE'] = '--'
	oGradVy['TITLE']     = '(-'+nabla+'V)$\downY$!C(mV/m)'
	
;	oGradVz['COLOR']     = 'Magenta'
;	oGradVz['LABEL']     = '(' + nabla + 'V)$\downZ$'
;	oGradVz['LINESTYLE'] = '--'
	oGradVz['TITLE']     = '(-'+nabla+'V)$\downZ$!C(mV/m)'
	
	Obj_Destroy, oRecipVec
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
;	title = StrUpCase(StrJoin(['mms1234', instr, mode, level], ' '))
;
;	oEx1 = MrVar_Get(ex_vnames[0])
;	oEx2 = MrVar_Get(ex_vnames[1])
;	oEx3 = MrVar_Get(ex_vnames[2])
;	oEx4 = MrVar_Get(ex_vnames[3])
;	odVx = MrVar_Get(gradvx_vname)
;	oEx1['PLOT_TITLE'] = title
;	oEx1['TITLE'] = 'Ex!C(mV/m)'
;	oEx1['LABEL'] = 'mms1'
;	oEx2['LABEL'] = 'mms2'
;	oEx3['LABEL'] = 'mms3'
;	oEx4['LABEL'] = 'mms4'
;	odVx['LABEL'] = '-Grad(V)'
;	
;	oEy = MrVar_Get(ey_vnames[0])
;	oEy['TITLE'] = 'Ey!C(mV/m)'
;	
;	oEz = MrVar_Get(ez_vnames[0])
;	oEz['TITLE'] = 'Ez!C(mV/m)'
;	
;	oEpar = MrVar_Get(epar_vnames[0])
;	oEpar['TITLE'] = 'E$\down||$!C(mV/m)'

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [ ex_vnames[0], ey_vnames[0], ez_vnames[0], v_out_vnames[0], $
	                      gradvx_vname, gradvy_vname, gradvz_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	win = MrVar_OPlotTS( ex_vnames[0], [ex_vnames[1:nSC-1], ex_bary_vname] )
	win = MrVar_OPlotTS( ey_vnames[0], [ey_vnames[1:nSC-1], ey_bary_vname] )
	win = MrVar_OPlotTS( ez_vnames[0], [ez_vnames[1:nSC-1], ez_bary_vname] )
	win = MrVar_OPlotTS( v_out_vnames[0], v_out_vnames[1:nSC-1] )
;	win = MrVar_OPlotTS( v_out_vnames[0], v_pca_vnames[0] )
;	win = MrVar_OPlotTS( v_out_vnames[1], v_pca_vnames[1] )
;	win = MrVar_OPlotTS( v_out_vnames[2], v_pca_vnames[2] )
;	win = MrVar_OPlotTS( v_out_vnames[3], v_pca_vnames[3] )
;	win = MrVar_OPlotTS( gradvx_vname, [gradvy_vname, gradvz_vname] )
	
	;Pretty-up the window
	win.name = 'E-GradV'
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 11]
	win    -> Refresh

;-------------------------------------------
; E vs. -Grad(V) Scatter Plot //////////////
;-------------------------------------------
	oEx_bary['COLOR'] = 'Black'
	oEy_bary['COLOR'] = 'Black'
	oEz_bary['COLOR'] = 'Black'
	
	;Scatter plots of Grad(V) and E_bary to view gain and offset parameters
	w2 = MrWindow( ASPECT  = 1.0, $
	               LAYOUT  = [3,1], $
	               NAME    = 'EvGradV', $
	               REFRESH = 0, $
	               XSIZE   = 1000 )
	
	;Scatter Plots
	p1 = MrVar_Plot( gradvx_vname, ex_bary_vname, $
	                 XTITLE        = '(' + nabla + 'V)$\downX$ (mV/m)', $
	                 XTICKINTERVAL = 0.04, $
	                 YTITLE        = 'E$\downX$ (mV/m)', $
	                 /CURRENT )
	p2 = MrVar_Plot( gradvy_vname, ey_bary_vname, $
	                 XTITLE = '(' + nabla + 'V)$\downY$ (mV/m)', $
	                 YTITLE = 'E$\downY$ (mV/m)', $
	                 /CURRENT )
	p3 = MrVar_Plot( gradvz_vname, ez_bary_vname, $
	                 XTITLE = '(' + nabla + 'V)$\downZ$ (mV/m)', $
	                 YTITLE = 'E$\downZ$ (mV/m)', $
	                 /CURRENT )
	
	;Fit data
	p1 -> GetData, x, y
	ifinite = where(finite(x) and finite(y))
	isort   = sort(x[ifinite])
	fit1    = LADFit(x[ifinite[isort]], y[ifinite[isort]])
	r1      = Correlate(x[ifinite[isort]], y[ifinite[isort]])
	xrange  = [ Min(x, MAX=xmax) < Min(y, MAX=ymax), xmax > ymax ]
	txt1    = String(fit1[1], fit1[0], r1, FORMAT='(%"m=%0.4e!Cb=%0.4e!CR=%0.4f")')
	lin1    = fit1[0] + fit1[1]*xrange
	
	p2 -> GetData, x, y
	ifinite = where(finite(x) and finite(y))
	isort   = sort(x[ifinite])
	fit2    = LADFit(x[ifinite[isort]], y[ifinite[isort]])
	r2      = Correlate(x[ifinite[isort]], y[ifinite[isort]])
	yrange  = [ Min(x, MAX=xmax) < Min(y, MAX=ymax), xmax > ymax ]
	txt2    = String(fit2[1], fit2[0], r2, FORMAT='(%"m=%0.4e!Cb=%0.4e!CR=%0.4f")')
	lin2    = fit2[0] + fit2[1]*yrange
	
	p3 -> GetData, x, y
	ifinite = where(finite(x) and finite(y))
	isort   = sort(x[ifinite])
	fit3    = LADFit(x[ifinite[isort]], y[ifinite[isort]])
	r3      = Correlate(x[ifinite[isort]], y[ifinite[isort]])
	zrange  = [ Min(x, MAX=xmax) < Min(y, MAX=ymax), xmax > ymax ]
	txt3    = String(fit3[1], fit3[0], r3, FORMAT='(%"m=%0.4e!Cb=%0.4e!CR=%0.4f")')
	lin3    = fit3[0] + fit3[1]*zrange
	
	;Text
	txt1 = MrText( 0.05, 0.9, txt1, $
	               CHARSIZE = 1.0, $
	               NAME     = 'Txt: Eqn-X', $
	               /RELATIVE, $
	               TARGET   = p1 )
	
	txt2 = MrText( 0.05, 0.9, txt2, $
	               CHARSIZE = 1.0, $
	               NAME     = 'Txt: Eqn-Y', $
	               /RELATIVE, $
	               TARGET   = p2 )
	
	txt3 = MrText( 0.05, 0.9, txt3, $
	               CHARSIZE = 1.0, $
	               NAME     = 'Txt: Eqn-Z', $
	               /RELATIVE, $
	               TARGET   = p3 )
	
	;Draw Line Fits
	ps1 = MrPlotS( xrange, lin1, $
	               COLOR     = 'Red', $
	               LINESTYLE = '--', $
	               NAME      = 'PS: X-Fit', $
	               TARGET    = p1 )
	
	ps2 = MrPlotS( yrange, lin2, $
	               COLOR     = 'Red', $
	               LINESTYLE = '--', $
	               NAME      = 'PS: Y-Fit', $
	               TARGET    = p2 )
	
	ps3 = MrPlotS( zrange, lin3, $
	               COLOR     = 'Red', $
	               LINESTYLE = '--', $
	               NAME      = 'PS: Z-Fit', $
	               TARGET    = p3 )
	
	;Make pretty
	p1 -> SetLayout, [1,1]
	trange   = MrVar_GetTRange()
	p2.title = StrJoin(StrSplit(StrMid(trange[0], 0, 19), 'T', /EXTRACT), ' ') + ' -- ' + StrMid(trange[1], 11, 8)
	w2.name = 'E-GradV-Scatter'
	w2 -> TrimLayout
	w2 -> Remove, w2 -> Get(/ALL, ISA='MrLegend')
	w2 -> Refresh

;-------------------------------------------
; E vs. -Grad(V) Corrected /////////////////
;-------------------------------------------

	;Correct Grad(V) by multiplying by the slope of the best fit line.
	oGradVx_cor = oGradVx * Abs(fit1[1])
	oGradVy_cor = oGradVy * Abs(fit2[1])
	oGradVz_cor = oGradVz * Abs(fit3[1])
	
	xrange = [oGradVx_cor.min < oEx_bary.min, oGradVx_cor.Max > oEx_bary.max]
	yrange = [oGradVy_cor.min < oEy_bary.min, oGradVy_cor.Max > oEy_bary.max]
	zrange = [oGradVz_cor.min < oEz_bary.min, oGradVz_cor.Max > oEz_bary.max]
	
	oEx_bary['AXIS_RANGE'] = xrange
	oEx_bary['LABEL']      = 'E$\downBC$'
	oEx_bary['PLOT_TITLE'] = 'Corrected E & -'+nabla+'V'
	oEy_bary['AXIS_RANGE'] = yrange
	oEy_bary['LABEL']      = 'E$\downBC$'
	oEz_bary['AXIS_RANGE'] = zrange
	oEz_bary['LABEL']      = 'E$\downBC$'
	
	oGradVx_cor -> SetName, gradvx_corr_vname
	oGradVx_cor -> Cache
	oGradVx_cor['AXIS_RANGE'] = xrange
	oGradVx_cor['COLOR']      = 'Red'
	oGradVx_cor['LABEL']      = '-' + nabla + 'V'
	oGradVx_cor['PLOT_TITLE'] = 'E vs. Corrected -'+nabla+'V'
	oGradVx_cor['TITLE']      = '(-'+nabla+'V)$\downX$!C(mV/m)'
	oGradVx_cor['UNITS']      = 'mV/m'
	
	oGradVy_cor -> SetName, gradvy_corr_vname
	oGradVy_cor -> Cache
	oGradVy_cor['AXIS_RANGE'] = yrange
	oGradVy_cor['COLOR']      = 'Red'
	oGradVy_cor['LABEL']      = '-' + nabla + 'V'
	oGradVy_cor['PLOT_TITLE'] = 'E vs. Corrected -'+nabla+'V'
	oGradVy_cor['TITLE']      = '(-'+nabla+'V)$\downY$!C(mV/m)'
	oGradVy_cor['UNITS']      = 'mV/m'
	
	oGradVz_cor -> SetName, gradvz_corr_vname
	oGradVz_cor -> Cache
	oGradVz_cor['AXIS_RANGE'] = zrange
	oGradVz_cor['COLOR']      = 'Red'
	oGradVz_cor['LABEL']      = '-' + nabla + 'V'
	oGradVz_cor['PLOT_TITLE'] = 'E vs. Corrected -'+nabla+'V'
	oGradVz_cor['TITLE']      = '(-'+nabla+'V)$\downZ$!C(mV/m)'
	oGradVz_cor['UNITS']      = 'mV/m'
	
	;Plot the data
	w3 = MrVar_PlotTS( [ex_bary_vname, ey_bary_vname, ez_bary_vname], $
	                   /NO_REFRESH )
	
	w3 = MrVar_OPlotTS( ex_bary_vname, gradvx_corr_vname )
	w3 = MrVar_OPlotTS( ey_bary_vname, gradvy_corr_vname )
	w3 = MrVar_OPlotTS( ez_bary_vname, gradvz_corr_vname )
	
	w3.name = 'E-GradV-Corrected'
	w3[0] -> SetLayout, [1,1]
	w3    -> TrimLayout
	w3.oxmargin = [13, 10]
	w3    -> Refresh

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
		fname1   = StrJoin( ['mms', edp_instr, edp_mode, level, 'egradv'], '_' )
		fname2   = StrJoin( ['mms', edp_instr, edp_mode, level, 'egradv-scatter'], '_' )
		fname3   = StrJoin( ['mms', edp_instr, edp_mode, level, 'egradv-corrected'], '_' )
		
		fname1   = FilePath( fname1, ROOT_DIR=output_dir )
		fname2   = FilePath( fname2, ROOT_DIR=output_dir )
		fname3   = FilePath( fname3, ROOT_DIR=output_dir )
		
		;Save the figure
		fout1 = MrVar_PlotTS_Save( win, fname1, output_ext )
		fout2 = MrVar_PlotTS_Save( w2,  fname2, output_ext )
		fout3 = MrVar_PlotTS_Save( w3,  fname3, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END