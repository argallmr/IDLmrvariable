; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_FPI_fMap
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
;   Generate a plot of FPI distribution functions:
;       1. Par-Perp1
;       2. Par-Perp2
;       3. Perp1-Perp2
;       4. Cuts of Par, Perp1, Perp2
;       5. Cuts of Par, Perp, Anti-Par
;
; :Categories:
;   MMS
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' }
;       SPECIES:    in, required, type=string, default='e'
;                   Particle species. Options are: { 'e' | 'i' }
;
; :Keywords:
;       FGM_INSTR:  in, optional, type=string, default='dfg'
;                   FGM instrument to use. Options are: { 'afg' | 'dfg' | 'fgm' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       OPTDESC:    in, optional, type=string, default=''
;                   Optional filename descriptor.
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
;       2018/05/25  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_FPI_fMap_SetUp, data
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------
; Setup Window /////////////////////////////
;-------------------------------------------
	xsize = 900
	ysize = 650
	
	;Positions of time series plots
	pos = MrLayout( [4,1], $
	                CHARSIZE = 2.0, $
	                OXMARGIN = [10,5], $
	                OYMARGIN = [25,2], $
	                WDIMS    = [xsize, ysize], $
	                XGAP     = 0, $
	                YGAP     = 0 )
	
	;Window
	win = MrWindow( ASPECT   = 1.0, $
	                LAYOUT   = [data.nCols,4], $
	                OXMARGIN = [12,5], $
	                OYMARGIN = [4,15], $
	                XGAP     = 0, $
	                XSIZE    = xsize, $
	                YGAP     = 0, $
	                YSIZE    = ysize, $
	                REFRESH  = 0 )
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	nSC   = N_Elements(data.sc)     ;Number of spacecraft (rows)
	nCols = 8                       ;Number of distributions (columns)
	t_iso = MrVar_GetTrange()       ;Total time interval as string
	t_ssm = MrVar_GetTRange('SSM')  ;Total time interval in SSM
	f_ttot = nCols * Mean(data.dt)  ;Total time spanned by distributions
		
	f_trange = DblArr(2,nSC)                               ;Time range of the NCOLS distributions
	f_xrange = [!Values.F_Infinity, -!Values.F_Infinity]   ;Common Energy range of the distributions
	f_yrange = [!Values.F_Infinity, -!Values.F_Infinity]   ;Common PSD range of the distributions
	

;-------------------------------------------
; Row of Distributions Per Spacecraft //////
;-------------------------------------------
	FOR i = 0, nSC - 1 DO BEGIN
		strsc  = String(i+1, FORMAT='(i1)')
		oEPara = data.oEPara[i]
		oEPerp = data.oEPerp[i]
		oEAnti = data.oEAnti[i]
		oEOmni = data.oEOmni[i]
		oEMax  = data.oEMax[i]
		oE     = oEPara['DEPEND_1']
		
		;Mark distribution times from the start of the previous minute
		f_t0     = (oEPara['TIMEVAR'])['DATA', data.istart, 'SSM']
		f_t0     = f_t0 - (f_t0 MOD 60.0)
		f_trange[*,i] = (oEPara['TIMEVAR'])['DATA', [data.istart, data.iend], 'SSM']
		
	;-------------------------------------------
	; Distributions ////////////////////////////
	;-------------------------------------------
		
		FOR j = data.istart, data.iend DO BEGIN
			jstr  = String(j, FORMAT='(i1)')
			f_xrange[0] <= Min(oE['DATA',oE -> Where(0, /GREATER)])
			f_xrange[1] >= oE.max
			f_yrange[0] <= Min(oEPara['DATA',oEPara -> Where(0, /GREATER)])
			f_yrange[1] >= (oEPara.max > oEPerp.max > oEAnti.max > oEMax.max)
			
			pomni = MrPlot( oE['DATA',j,*], oEOmni['DATA',j,*], $
			                /CURRENT, $
			                COLOR       = 'Magenta', $
			                NAME        = 'MMS' + strsc + ' Col ' + jstr + ' EOmni', $
			                /XLOG, $
			                XTICKFORMAT = '(a1)', $
			                /YLOG, $
			                YTICKFORMAT = '(a1)' )
			ppara = MrPlot( oE['DATA',j,*], oEPara['DATA',j,*], $
			                COLOR       = 'Blue', $
			                NAME        = 'MMS' + strsc + ' Col ' + jstr + ' EPara', $
			                OVERPLOT    = pomni, $
			                XTICKFORMAT = '(a1)', $
			                YTICKFORMAT = '(a1)' )
			panti = MrPlot( oE['DATA',j,*], oEAnti['DATA',j,*], $
			                COLOR       = 'Forest Green', $
			                NAME        = 'MMS' + strsc + ' Col ' + jstr + ' EAnti', $
			                OVERPLOT    = pomni, $
			                XTICKFORMAT = '(a1)', $
			                YTICKFORMAT = '(a1)' )
			pperp = MrPlot( oE['DATA',j,*], oEPerp['DATA',j,*], $
			                COLOR       = 'Red', $
			                NAME        = 'MMS' + strsc + ' Col ' + jstr + ' EPerp', $
			                OVERPLOT    = pomni, $
			                XTICKFORMAT = '(a1)', $
			                YTICKFORMAT = '(a1)' )
			pmaxb = MrPlot( oE['DATA',j,*], oEMax['DATA',j,*], $
			                COLOR       = 'Black', $
			                NAME        = 'MMS' + strsc + ' Col ' + jstr + ' EMax', $
			                OVERPLOT    = pomni, $
			                XTICKFORMAT = '(a1)', $
			                YTICKFORMAT = '(a1)' )
			
			;Legend with time
			tstr = String( (oEPara['TIMEVAR'])['DATA', j, 'SSM'] - f_t0, FORMAT='(f6.3)' )
			oTxt = MrText( 0.0, 0.05, tstr, $
			               ALIGNMENT = 0.0, $
			               CHARSIZE  = 1.0, $
			               NAME      = 'Txt: MMS' + strsc + ' Col ' + jstr, $
			               /RELATIVE, $
			               TARGET    = pomni )
			
			IF j EQ data.istart THEN pomni -> SetProperty, YTITLE='MMS' + strsc
			IF i EQ nSC - 1 && j EQ data.istart THEN pomni -> SetProperty, YTICKFORMAT='', YTITLE=oEPara['TITLE']
			IF i EQ nSC - 1 && j EQ data.istart THEN pomni -> SetProperty, XTICKFORMAT='logtickformat', XTITLE='E (eV)'
		ENDFOR
	ENDFOR
	
	;Create a legend
	oL = MrLegend( ALIGNMENT    = 'NE', $
	               FILL_COLOR   = '', $
	               LABEL        = ['Par', 'Perp', 'Anti', 'Omni', 'Max'], $
	               LINESTYLE    = 'None', $
	               NAME         = 'Lgd: f', $
	               /NORMAL, $
	               POSITION     = [1.0, 0.65], $
	               SAMPLE_WIDTH = 0.0, $
	               TARGET       = pomni, $
	               TEXT_COLOR   = ['Blue', 'Red', 'Forest Green', 'Magenta', 'Black'] )
	
	win -> SetGlobal, XRANGE=f_xrange, YRANGE=f_yrange

;-------------------------------------------
; Plot Time Series /////////////////////////
;-------------------------------------------
	;Display more data in the time series plots
	;   - Make the time range symmetric about the center if the distribution functions
	ts_tcenter = Min(f_trange[0,*]) + (Max(f_trange[1,*]) - Min(f_trange[0,*])) / 2.0
	ts_trange  = ts_tcenter + 4.0*f_ttot*[-1.0, 1.0]
	ts_t0      = ts_trange[0] - (ts_trange[0] MOD 60.0)
	ts_trange  = ts_trange - ts_t0
	
	FOR i = 0, nSC - 1 DO BEGIN
		scstr = String(i+1, FORMAT='(i1)')
		p = MrVar_Plot( data.varnames[i], $
		                /CURRENT, $
		                NAME      = 'MMS' + scstr + ' Time Series', $
		                NO_LEGEND = (i LT nSC-1), $
		                POSITION  = pos[*,i], $
		                TITLE     = 'MMS' + scstr, $
		                XTITLE    = '' )
		
		;Convert time to seconds past t0
		p -> GetData, x, y
		p -> SetData, x-ts_t0, y
		p -> SetProperty, XTICKV=0, XTICKS=0, XTICKFORMAT='', XRANGE=ts_trange
		IF i GT 0 THEN p -> SetProperty, YTITLE='', YTICKFORMAT='(a1)'
		
		;Draw lines marking distribution functions
		oS = MrPlotS( f_trange[[0,0],i]-ts_t0, p.yrange, $
		              NAME   = 'Line: MMS' + scstr + ' t0', $
		              TARGET = p )
		oS = MrPlotS( f_trange[[1,1],i]-ts_t0, p.yrange, $
		              NAME   = 'Line: MMS' + scstr + ' t1', $
		              TARGET = p )
	ENDFOR
	
	oText = MrText( 0.5, 0.675, 'Seconds Past ' + StrMid(ssm_to_hms(Floor(ts_t0)), 0, 12), $
	                ALIGNMENT = 0.5, $
	                NAME      = 'Txt: TS XTitle' )

;-------------------------------------------
; Next Iteration ///////////////////////////
;-------------------------------------------
	;Clean up this window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	RETURN, win
END


;+
;
;-
FUNCTION MrMMS_Plot_4sc_FPI_fMap_Update, win, data
	Compile_Opt idl2
	On_Error, 2
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	nSC   = N_Elements(data.sc)     ;Number of spacecraft (rows)
	nCols = 8                       ;Number of distributions (columns)
	t_iso = MrVar_GetTrange()       ;Total time interval as string
	t_ssm = MrVar_GetTRange('SSM')  ;Total time interval in SSM
	f_ttot = nCols * Mean(data.dt)  ;Total time spanned by distributions
		
	f_trange = DblArr(2,nSC)                               ;Time range of the NCOLS distributions
	f_xrange = [!Values.F_Infinity, -!Values.F_Infinity]   ;Common Energy range of the distributions
	f_yrange = [!Values.F_Infinity, -!Values.F_Infinity]   ;Common PSD range of the distributions

;-------------------------------------------
; Row of Distributions Per Spacecraft //////
;-------------------------------------------
	FOR i = 0, nSC - 1 DO BEGIN
		strsc  = String(i+1, FORMAT='(i1)')
		oEOmni = data.oEOmni[i]
		oEPara = data.oEPara[i]
		oEAnti = data.oEAnti[i]
		oEPerp = data.oEPerp[i]
		oEMax  = data.oEMax[i]
		oE     = oEPara['DEPEND_1']
		
		;Mark distribution times from the start of the previous minute
		f_t0     = (oEOmni['TIMEVAR'])['DATA', data.istart, 'SSM']
		f_t0     = f_t0 - (f_t0 MOD 60.0)
		f_trange[*,i] = (oEOmni['TIMEVAR'])['DATA', [data.istart, data.iend], 'SSM']
		
	;-------------------------------------------
	; Distributions ////////////////////////////
	;-------------------------------------------
		
		FOR j = data.istart, data.iend DO BEGIN
			jstr  = String(j-data.istart, FORMAT='(i1)')
			f_xrange[0] <= Min(oE['DATA',oE -> Where(0, /GREATER)])
			f_xrange[1] >= oE.max
			f_yrange[0] <= Min(oEOmni['DATA',oEOmni -> Where(0, /GREATER)])
			f_yrange[1] >= (oEOmni.max > oEPara.max > oEPerp.max > oEAnti.max > oEMax.max)
			
			
			pomni = win['MMS' + strsc + ' Col ' + jstr + ' EOmni']
			ppara = win['MMS' + strsc + ' Col ' + jstr + ' EPara']
			panti = win['MMS' + strsc + ' Col ' + jstr + ' EAnti']
			pperp = win['MMS' + strsc + ' Col ' + jstr + ' EPerp']
			pmaxb = win['MMS' + strsc + ' Col ' + jstr + ' EMax']
			
			pomni -> SetData, oE['DATA',j,*], oEOmni['DATA',j,*]
			panti -> SetData, oE['DATA',j,*], oEAnti['DATA',j,*]
			ppara -> SetData, oE['DATA',j,*], oEPara['DATA',j,*]
			pperp -> SetData, oE['DATA',j,*], oEPerp['DATA',j,*]
			pmaxb -> SetData, oE['DATA',j,*], oEMax['DATA',j,*]
			
			oTxt = win['Txt: MMS' + strsc + ' Col ' + jstr]
			oTxt.String = String( (oEOmni['TIMEVAR'])['DATA', j, 'SSM'] - f_t0, FORMAT='(f6.3)' )
		ENDFOR
	ENDFOR

;-------------------------------------------
; Plot Time Series /////////////////////////
;-------------------------------------------
	;Display more data in the time series plots
	;   - Make the time range symmetric about the center if the distribution functions
	ts_tcenter = Min(f_trange[0,*]) + (Max(f_trange[1,*]) - Min(f_trange[0,*])) / 2.0
	ts_trange  = ts_tcenter + 4.0*f_ttot*[-1.0, 1.0]
	ts_t0      = ts_trange[0] - (ts_trange[0] MOD 60.0)
	ts_trange  = ts_trange - ts_t0
	
	ts_yrange = [!Values.F_Infinity, -!Values.F_Infinity]
	
	FOR i = 0, nSC - 1 DO BEGIN
		scstr    = String(i+1, FORMAT='(i1)')
		p        = win['MMS' + scstr + ' Time Series']
		p -> GetData, x, y
		ts_yrange[0] <= Min(y)
		ts_yrange[1] >= Max(y)
		
		;Draw lines marking distribution functions
		oS1 = win['Line: MMS' + scstr + ' t0']
		oS2 = win['Line: MMS' + scstr + ' t1']
		oS1 -> SetData, f_trange[[0,0],i]-ts_t0, p.yrange
		oS2 -> SetData, f_trange[[1,1],i]-ts_t0, p.yrange
	ENDFOR
	
	oText = win['Txt: TS XTitle']
	oText.string = 'Seconds Past ' + StrMid(ssm_to_hms(Floor(ts_t0)), 0, 12)
	
;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------
	win -> SetGlobal, XRANGE=f_xrange, YRANGE=f_yrange
	win['MMS1 Time Series'] -> SetProperty, XRANGE=ts_trange, YRANGE=ts_yrange
	win['MMS2 Time Series'] -> SetProperty, XRANGE=ts_trange, YRANGE=ts_yrange
	win['MMS3 Time Series'] -> SetProperty, XRANGE=ts_trange, YRANGE=ts_yrange
	win['MMS4 Time Series'] -> SetProperty, XRANGE=ts_trange, YRANGE=ts_yrange
	RETURN, win
END


;+
;
;-
FUNCTION MrMMS_Plot_4sc_FPI_fMap, mode, species, time, $
FGM_INSTR=fgm_instr, $
FRANGE=frange, $
LEVEL=level, $
NO_LOAD=no_load, $
OPTDESC=optdesc, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
RAGER=rager, $
TAIL=tail, $
TRANGE=trange, $
VARNAMES=varnames, $
VRANGE=vrange
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF Obj_Valid(oVid) THEN Obj_Destroy, oVid
		IF N_Elements(win) GT 0 THEN Obj_Destroy, win
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	kB       = MrConstants('k_B')
	m_e      = MrConstants('m_e')
	m_i      = MrConstants('m_H')
	eV2J     = MrConstants('eV2J')
	eV2K     = MrConstants('eV2K')
	sc       = 'mms' + ['1', '2', '3', '4']
	tf_movie = N_Elements(time) EQ 0
	tf_load  = ~Keyword_Set(no_load)
	IF N_Elements(coords)  EQ 0 THEN coords  = 'gse'
	IF N_Elements(level)   EQ 0 THEN level   = 'l2'
	IF N_Elements(mode)    EQ 0 THEN mode    = 'fast'
	IF N_Elements(species) EQ 0 THEN species = 'e'
	IF N_Elements(trange)  GT 0 THEN MrVar_SetTRange, trange

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	fpi_instr   = 'd' + species + 's'
	fpi_optdesc = fpi_instr + '-dist'
	fpi_mode    = mode EQ 'brst' ? mode : 'fast'
	IF N_Elements(coords) EQ 0 THEN BEGIN
		CASE level OF
			'ql': fpi_coords = 'dbcs'
			ELSE: fpi_coords = 'gse'
		ENDCASE
	ENDIF
	
	;FGM
	IF N_Elements(fgm_instr) EQ 0 THEN BEGIN
		CASE level OF
			'l2': fgm_instr = 'fgm'
			ELSE: fgm_instr = 'dfg'
		ENDCASE
	ENDIF
	fgm_level  = fgm_instr EQ 'fgm'  ? 'l2'   : 'l2pre'
	fgm_coords = coords    EQ 'dbcs' ? 'dmpa' : coords
	fgm_mode   = mode      EQ 'brst' ? mode   : 'srvy'

	;Source names
	b_vnames      = sc + '_' + StrJoin( [fgm_instr, 'b',     fgm_coords, fgm_mode, level], '_' )
	bvec_vnames   = sc + '_' + StrJoin( [fgm_instr, 'bvec',  fgm_coords, fgm_mode, level], '_' )
	bmag_vnames   = sc + '_' + StrJoin( [fgm_instr, 'bmag',  fgm_coords, fgm_mode, level], '_' )
	f_vnames      = sc + '_' + StrJoin( [fpi_instr, 'dist',          mode], '_' )
	n_vnames      = sc + '_' + StrJoin( [fpi_instr, 'numberdensity', mode], '_' )
	v_vnames      = sc + '_' + StrJoin( [fpi_instr, 'bulkv', coords, mode], '_' )
	t_para_vnames = sc + '_' + StrJoin( [fpi_instr, 'temppara', mode], '_' )
	t_perp_vnames = sc + '_' + StrJoin( [fpi_instr, 'tempperp', mode], '_' )
	e_para_vnames = sc + '_' + StrJoin( [fpi_instr, 'energyspectr', 'par', mode], '_' )
	e_perp_vnames = sc + '_' + StrJoin( [fpi_instr, 'energyspectr', 'perp', mode], '_' )
	e_anti_vnames = sc + '_' + StrJoin( [fpi_instr, 'energyspectr', 'anti', mode], '_' )
	e_omni_vnames = sc + '_' + StrJoin( [fpi_instr, 'energyspectr', 'omni', mode], '_' )

	;Derived names
	t_vnames          = sc + '_' + StrJoin( [fpi_instr, 'temp', mode], '_' )
	f_max_vnames      = sc + '_' + StrJoin( [fpi_instr, 'dist', 'maxwell', mode], '_' )
	e_para_psd_vnames = sc + '_' + StrJoin( [fpi_instr, 'energyspectr', 'par',  'psd', mode], '_' )
	e_perp_psd_vnames = sc + '_' + StrJoin( [fpi_instr, 'energyspectr', 'perp', 'psd', mode], '_' )
	e_anti_psd_vnames = sc + '_' + StrJoin( [fpi_instr, 'energyspectr', 'anti', 'psd', mode], '_' )
	e_omni_psd_vnames = sc + '_' + StrJoin( [fpi_instr, 'energyspectr', 'omni', 'psd', mode], '_' )
	
	;Variable to plot
	IF N_Elements(varnames) EQ 0 THEN varnames = b_vnames
	
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, '', fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = fgm_level, $
		                     VARFORMAT = b_vname

		;DIST
		;   - The distribution is corrected for internally
		;     generated photoelectrons.
		MrMMS_FPI_Load_Dist3d, 'mms1', mode, species, $
		                       /APPLY_MODEL, $
		                       COORD_SYS = coords, $
		                       LEVEL     = level, $
		                       RAGER     = rager
		MrMMS_FPI_Load_Dist3d, 'mms2', mode, species, $
		                       /APPLY_MODEL, $
		                       COORD_SYS = coords, $
		                       LEVEL     = level, $
		                       RAGER     = rager
		MrMMS_FPI_Load_Dist3d, 'mms3', mode, species, $
		                       /APPLY_MODEL, $
		                       COORD_SYS = coords, $
		                       LEVEL     = level, $
		                       RAGER     = rager
		MrMMS_FPI_Load_Dist3d, 'mms4', mode, species, $
		                       /APPLY_MODEL, $
		                       COORD_SYS = coords, $
		                       LEVEL     = level, $
		                       RAGER     = rager
		
		;MOMENTS
		MrMMS_FPI_Load_Data, '', mode, $
		                     RAGER     = rager, $
		                     OPTDESC   = fpi_instr+'-moms', $
		                     VARFORMAT = ['*numberdensity_'+mode, '*bulkv_'+coords+'*', $
		                                  '*energyspectr_' + ['par', 'anti', 'perp', 'omni'] + '*', $
		                                  '*tempperp*', '*temppara*']
	ENDIF
	
;-------------------------------------------
; Prep Data ////////////////////////////////
;-------------------------------------------
	nSC  = N_Elements(sc)
	dt   = FltArr(nSC)
	npts = LIndGen(nSC)
	
	FOR i = 0, nSC - 1 DO BEGIN
		thesc = 'MMS' + String(i+1, FORMAT='(i1)')
		
	;-------------------------------------------
	; Convert EFLUX to PSD /////////////////////
	;-------------------------------------------
		;Grab data
		oEPara = MrVar_Get(e_para_vnames[i])
		oEPerp = MrVar_Get(e_perp_vnames[i])
		oEAnti = MrVar_Get(e_anti_vnames[i])
		oEOmni = MrVar_Get(e_omni_vnames[i])
		oE     = oEPara['DEPEND_1']
		dims   = Size(oEPara, /DIMENSIONS)
		nTime  = dims[0]
		nEnergy = dims[1]
		
		;Convert from EFlux to PSD
		N         = species EQ 'i' ? 1 : MrConstants('m_e', /DOUBLE) / MrConstants('m_H', /DOUBLE)
		eflux2psd = N^2.0 * 5.44933e-25
		oEOmni = eflux2psd * oEOmni / oE^2
		oEPara = eflux2psd * oEPara / oE^2
		oEAnti = eflux2psd * oEAnti / oE^2
		oEPerp = eflux2psd * oEPerp / oE^2
		
		;Cache
		oEOmni -> SetName, e_omni_psd_vnames[i]
		oEOmni -> Cache
		oEOmni['DEPEND_1']      = oE
		oEOmni['TITLE']         = thesc + '!C(s^3/cm^6)'
		oEOmni['UNITS']         = 's^3/cm^6'
		oEOmni['SI_CONVERSION'] = '1e12>s^3/m^6'
		
		oEPara -> SetName, e_para_psd_vnames[i]
		oEPara -> Cache
		oEPara['DEPEND_1']      = oE
		oEPara['TITLE']         = thesc + '!C(s^3/cm^6)'
		oEPara['UNITS']         = 's^3/cm^6'
		oEPara['SI_CONVERSION'] = '1e12>s^3/m^6'
		
		oEAnti -> SetName, e_anti_psd_vnames[i]
		oEAnti -> Cache
		oEAnti['DEPEND_1']      = oE
		oEAnti['TITLE'] = thesc + '!C(s^3/cm^6)'
		oEAnti['UNITS'] = 's^3/cm^6'
		oEAnti['SI_CONVERSION'] = '1e12>s^3/m^6'
		
		oEPerp -> SetName, e_perp_psd_vnames[i]
		oEPerp -> Cache
		oEPerp['DEPEND_1']      = oE
		oEPerp['TITLE']         = thesc + '!C(s^3/cm^6)'
		oEPerp['UNITS']         = 's^3/cm^6'
		oEPerp['SI_CONVERSION'] = '1e12>s^3/m^6'
		
		
	;-------------------------------------------
	; Scalar Temperature ///////////////////////
	;-------------------------------------------
		;Grab the data
		oTPara = MrVar_Get(t_para_vnames[i])
		oTPerp = MrVar_Get(t_perp_Vnames[i])
		
		;Total temperature
		oT = (oTPara + 2.0*oTPerp) / 3.0
		oT -> SetName, t_vnames[i]
		oT -> Cache
		
	;-------------------------------------------
	; Maxwellian Distribution //////////////////
	;-------------------------------------------
		oEOmni = MrVar_Get(e_omni_vnames[i])
		oN     = MrVar_Get(n_vnames[i])
		oE     = oEOmni['DEPEND_1']
		vsqr   = 2.0*eV2J*oE['DATA']/m_e
		ttemp  = Rebin(oT['DATA'], nTime, nEnergy)
		ntemp  = Rebin(oN['DATA'], nTime, nEnergy)
;		f_max  = 2.0 * Sqrt(oE['DATA']/!pi) / (eV2K*kB*temp)^(3.0/2.0) * Exp(-eV2J*oE['DATA']/(eV2K*kB*temp))
		f_max  = ntemp * 1e-6 * (m_e / (2*!pi*eV2K*kB*ttemp))^(3.0/2.0) * Exp(-m_e*vsqr / (2.0*eV2K*kB*ttemp))
		ntemp  = !Null
		ttemp  = !Null

		ofMax  = MrTimeSeries(oEOmni['TIMEVAR'], f_max, /NO_COPY)
		ofMax -> SetName, f_max_vnames[i]
		ofMax -> Cache
		ofMax['DEPEND_1'] = oE
		
	;-------------------------------------------
	; Useful Info //////////////////////////////
	;-------------------------------------------
		dt[i]   = oEOmni['TIMEVAR'] -> GetSI()
		nPts[i] = oEOmni -> GetNPts()
	ENDFOR

;-------------------------------------------
; Output File //////////////////////////////
;-------------------------------------------
	IF ~tf_movie THEN BEGIN
		oEOmni = MrVar_Get(e_omni_psd_vnames[0])
		i0     = oEPara['TIMEVAR'] -> Value_Locate(time)
	ENDIF ELSE BEGIN
		i0 = 0
	ENDELSE
	
	;Create a video object for saving
	tf_save = 0B
	IF N_Elements(output_ext) GT 0 || N_Elements(output_dir) GT 0 THEN BEGIN
		tf_save = 1B
		IF N_Elements(output_dir) EQ 0 THEN output_dir = File_Search('~', /TEST_DIRECTORY)
		IF N_Elements(output_ext) EQ 0 THEN output_ext = tf_movie ? 'avi' : 'png'
		
		;Time stamp of file
		ftime = MrVar_GetTRange()
		ftime = StrMid(ftime[0], 0, 10) + '_' + StrMid(ftime[0], 11, 8) + '_' + StrMid(ftime[1], 11, 8)
		ftime = StrJoin(StrSplit(ftime, '-:T', /EXTRACT), '')
		
		;Save the file
		itime = tf_movie ? '' : '-i' + String(i0, FORMAT='(i0)')
		fname = StrJoin(['mms', fpi_instr, fpi_mode, level, 'fmap'+itime, ftime + '.' + output_ext], '_')
		fname = FilePath(fname, ROOT_DIR=output_dir)
		
		;Video object
		IF tf_movie THEN BEGIN
			xsize = 900
			ysize = 650
			oVid = IDLffVideoWrite(fname, FORMAT='avi')
			iVid = oVid -> AddVideoStream(xsize, ysize, 3)
		ENDIF
	ENDIF
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	nCols = 8                       ;Number of distributions (columns)
	data = { sc:     sc, $
	         istart: i0, $
	         iend:   i0 + nCols - 1, $
	         nCols:  nCols, $
	         dt:     dt, $
	         oEOmni: MrVar_Get(e_omni_psd_vnames), $
	         oEPara: MrVar_Get(e_para_psd_vnames), $
	         oEAnti: MrVar_Get(e_anti_psd_vnames), $
	         oEPerp: MrVar_Get(e_perp_psd_vnames), $
	         oEMax:  MrVar_Get(f_max_vnames), $
	         varnames: varnames $
	       }
	
	;Loop over each time interval
	nMax = tf_movie ? Min(nPts) : i0+nCols - 1
	WHILE data.istart LT nMax DO BEGIN
		Print, data.istart, data.iend, nMax, FORMAT='(%"Points %i-%i of %i")'
		
		;Create the first
		IF data.istart EQ i0 $
			THEN win = MrMMS_Plot_4sc_FPI_fMap_SetUp(data) $
			ELSE win = MrMMS_Plot_4sc_FPI_fMap_Update(win, data)
		
	;-------------------------------------------
	; Next Iteration ///////////////////////////
	;-------------------------------------------
		;Clean up this window
		win -> Refresh
		IF tf_save && tf_movie THEN BEGIN
			win -> GetProperty, XSIZE=xsize, YSIZE=ysize
			img  = TVRD(0, 0, xsize, ysize, TRUE=1)
			time = oVid -> Put(iVid, img)
		ENDIF
		
		;Next
		win         -> Refresh, /DISABLE
		data.istart += nCols
		data.iend   += nCols
	ENDWHILE

;-------------------------------------------
; Finish ///////////////////////////////////
;-------------------------------------------
	IF tf_save THEN BEGIN
		IF tf_movie $
			THEN Obj_Destroy, oVid $
			ELSE win -> Save, fname
	ENDIF

	RETURN, win
END