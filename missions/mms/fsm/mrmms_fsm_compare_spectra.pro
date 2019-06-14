

;http://www.vibrationresearch.com/university/lesson/what-is-the-psd/
;https://www.dsprelated.com/freebooks/mdft/
;	https://www.dsprelated.com/freebooks/mdft/DFT_Derived.html
;http://holometer.fnal.gov/GH_FFT.pdf
;https://dsp.stackexchange.com/questions/28399/dft-normalization-for-amplitude-estimation
;http://www.ws.binghamton.edu/fowler/fowler%20personal%20page/EE521_files/VI-02%20Practical%20Classical%20Methods_2007.pdf


;mrvar_settrange, '2018-01-08T' + ['06:40:50', '06:41:30']
mrvar_settrange, '2017-11-24T' + ['01:10:03', '02:10:03']
;mrvar_settrange, '2017-11-24T' + ['01:10:03', '01:12:03']
sc     = 'mms1'
mode   = 'brst'
level  = 'l2'
coords = 'gse'

output_dir = '/home/argall/figures/20151206/'
output_ext = 'png'
tf_load   = 1B
tf_hipass = 0B
tf_edp    = 0B

;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------

IF tf_load THEN BEGIN
	;FSM Data
	MrMMS_Load_Data, sc, 'fsm', 'brst', 'l3', $
	                 /TEAM_SITE, $
	                 OPTDESC  = '8khz', $
	                 VARFORMAT = '*_b_gse_brst_l3'
	
	;FGM Data
	MrMMS_FGM_Load_Data, sc, mode, $
	                     VARFORMAT = '*_b_gse_brst_l2'
	
	;SCM Data
	MrMMS_Load_Data, sc, 'scm', mode, 'l2', $
	                 OPTDESC   = 'scb', $
	                 VARFORMAT = '*_acb_gse_scb_'+mode+'_'+level
	
	;EDP Data
	IF tf_edp THEN BEGIN
		MrMMS_Load_Data, sc, 'edp', mode, 'l2', $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = '*_dce_gse_'+mode+'_'+level
	ENDIF
ENDIF

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------

;FSM
fsm_level = 'l3'
fsm_mode = 'brst'
fsm_optdesc = '8khz'

;SCM
scm_optdesc = 'scb'

;Variable names
b_fsm_vname = StrJoin([sc, 'fsm', 'b', coords, fsm_mode, fsm_level], '_')
b_fgm_vname = StrJoin([sc, 'fgm', 'bvec', coords, mode, level], '_')
b_scm_vname = StrJoin([sc, 'scm', 'acb', coords, scm_optdesc, mode, level], '_')
e_edp_vname = StrJoin([sc, 'edp', 'dce', coords, mode, level], '_')
psd_fsm_vname = StrJoin([sc, 'fsm', 'psd', coords, fsm_mode, fsm_level], '_')
psd_fgm_vname = StrJoin([sc, 'fgm', 'psd', coords, mode, level], '_')
psd_scm_vname = StrJoin([sc, 'scm', 'psd', coords, scm_optdesc, mode, level], '_')
psd_edp_vname = StrJoin([sc, 'edp', 'psd', coords, mode, level], '_')

;Derived names
b_fsm_clip_vname = StrJoin([sc, 'fsm', 'b', coords, 'clipped', fsm_mode, fsm_level], '_')
b_fgm_hipass_vname = StrJoin([sc, 'fgm', 'bvec', coords, 'hipass', mode, level], '_')
nf_fgm_vname = StrJoin([sc, 'fgm', 'inflight', 'noise', 'floor'], '_')
nf_scm_vname = StrJoin([sc, 'scm', 'inflight', 'noise', 'floor'], '_')

;-------------------------------------------
; FGM: High-Pass Filter ////////////////////
;-------------------------------------------
IF tf_hipass THEN BEGIN
	fs = 128.0       ;Sample frequency (S/s)
	T  = 1000.0      ;Signal duration (s)
	N  = T*fs        ;Total number of samples (S)
	fn = fs / 2.0    ;Nyquist frequency (Hz)
	f0 = 1.0         ;Low-pass filter cutoff (Hz)
	f1 = 0.0         ;High-pass filter cutoff (Hz)
	A  = 50          ;Power
	m  = 128         ;Number of filter coefficients
	df = fs / N      ;Frequency resolution
	fh = [0:fn:df]   ;Filter frequencies
	it = [-m:m]      ;Sample number of the time-domain filter

	;Impulse response function
	oB_fgm = MrVar_Get(b_fgm_vname)
	oB_fgm = oB_fgm -> Digital_Filter(f0/fn, f1/fn, A, m, $
	                                  /CACHE, $
	                                  NAME=b_fgm_hipass_vname)
	b_fgm_vname = b_fgm_hipass_vname
ENDIF

;-------------------------------------------
; Remove End Points ////////////////////////
;-------------------------------------------
oB_fsm = MrVar_Get(b_fsm_vname)
iGood  = Where( Finite(oB_fsm['DATA',*,0]), nGood)
oB_fsm = oB_fsm[iGood,*]
oB_fsm -> SetName, b_fsm_clip_vname
oB_fsm -> Cache
b_fsm_vname = b_fsm_clip_vname

;-------------------------------------------
; Compute PSD //////////////////////////////
;-------------------------------------------

;Grab the data
;oB_fsm = MrVar_Get(b_fsm_vname)
oB_fgm = MrVar_Get(b_fgm_vname)
oB_scm = MrVar_Get(b_scm_vname)

;Spectral analysis
t = 60.0 ;seconds
overlap = 0.5 ;Number of points to shift = (1.0 - overlap)*NFFT
taper = 'hanning'
tf_zeropad = 1B
dt_fsm = oB_fsm['TIME',1,'SSM'] - oB_fsm['TIME',0,'SSM']
dt_fgm = oB_fgm['TIME',1,'SSM'] - oB_fgm['TIME',0,'SSM']
dt_scm = oB_scm['TIME',1,'SSM'] - oB_scm['TIME',0,'SSM']
;dt_fsm = oB_fsm['TIMEVAR'] -> GetSI()
;dt_fgm = oB_fgm['TIMEVAR'] -> GetSI()
;dt_scm = oB_scm['TIMEVAR'] -> GetSI()
;dt_edp = oE_edp['TIMEVAR'] -> GetSI()
n_fsm = Round(t/dt_fsm)
n_fgm = Round(t/dt_fgm)
n_scm = Round(t/dt_scm)

;trange  = '2018-01-08T' + ['06:41:10.65', '06:41:10.8']
;itrange_fgm = oB_fgm['TIMEVAR'] -> Value_Locate(trange)
;itrange_scm = oB_scm['TIMEVAR'] -> Value_Locate(trange)
;itrange_fsm = oB_fsm['TIMEVAR'] -> Value_Locate(trange)
;itrange_edp = oE_edp['TIMEVAR'] -> Value_Locate(trange)


;NOTE: If NSHIFT is a floating point number, FFT freezes on the 35th iteration.

;FGM
oB_fgm_psd = oB_fgm -> PSD(n_fgm, Round(n_fgm*(1.0 - overlap)), $
                           /CACHE, $
                           NAME    = psd_fgm_vname, $
                           RANGE   = itrange_fgm, $
                           WINDOW  = taper, $
                           ZEROPAD = tf_zeropad)

;SCM
oB_scm_psd = oB_scm -> PSD(n_scm, Round(n_scm*(1.0 - overlap)), $
                           /CACHE, $
                           NAME    = psd_scm_vname, $
                           RANGE   = itrange_scm, $
                           WINDOW  = taper, $
                           ZEROPAD = tf_zeropad)

;FSM
oB_fsm_psd = oB_fsm -> PSD(n_fsm, n_fsm/2, $ ;n_fsm*(1.0 - overlap), $
                           /CACHE, $
                           NAME    = psd_fsm_vname, $
                           RANGE   = itrange_fsm, $
                           WINDOW  = taper, $
                           ZEROPAD = tf_zeropad)


IF tf_edp THEN BEGIN
	oE_edp = MrVar_Get(e_edp_vname)
	dt_edp = oE_edp['TIME',1,'SSM'] - oE_edp['TIME',0,'SSM']
	n_edp = Round(t/dt_edp)
	
	;EDP
	oE_edp_psd = oE_edp -> PSD(n_edp, n_edp*(1.0 - overlap), $
	                           /CACHE, $
	                           NAME    = psd_edp_vname, $
	                           RANGE   = itrange_edp, $
	                           WINDOW  = taper, $
	                           ZEROPAD = tf_zeropad)
ENDIF

;-------------------------------------------
; IDL PSD //////////////////////////////////
;-------------------------------------------
;Catch, the_error
;IF the_error EQ 0 THEN BEGIN
;
;pB_fgm = oB_fgm['PTR']
;nPts   = oB_fgm -> GetNPts()
;istart = 0
;iend   = n_fgm - 1
;nShift = n_fgm / 2
;nMax   = Floor( (nPts - n_fgm) / nShift ) + 1
;
;psd_fgm = FltArr(n_fgm/2+1, nMax, 3)
;FOR i = 0, nMax - 1 DO BEGIN
;	FOR j = 0, 2 DO BEGIN
;		psd_fgm[0,i,j] = FFT_PowerSpectrum((*pB_fgm)[istart:iend,j], dt_fgm, FREQ=f_fgm)
;	ENDFOR
;	istart += nshift
;	iend   += nshift
;ENDFOR
;psd_fgm = Mean(psd_fgm, DIMENSION=2) / dt_fgm
;
;of_fgm   = MrVariable(f_fgm, /NO_COPY)
;opsd_fgm = MrVariable(psd_fgm, /NO_COPY)
;of_fgm['TITLE']      = 'f (Hz)'
;opsd_fgm['DEPEND_0'] = of_fgm
;
;ENDIF ELSE STOP
;Catch, /CANCEL

;-------------------------------------------
; In-Flight Noise Floors ///////////////////
;-------------------------------------------

of_fgm = oB_fgm_psd['DEPEND_0']
of_scm = oB_scm_psd['DEPEND_0']

;FGM
oNoise_fgm = 5.81e-6 * of_fgm^(-0.198)
oNoise_fgm -> SetName, nf_fgm_vname
oNoise_fgm -> Cache
oNoise_fgm['DEPEND_0'] = of_fgm

;SCM
oNoise_scm = 2.01e-4 * of_scm^(-2.19)
oNoise_scm -> SetName, nf_scm_vname
oNoise_scm -> Cache
oNoise_scm['DEPEND_0'] = of_scm

;-------------------------------------------
; Plot B PSD ///////////////////////////////
;-------------------------------------------

;RANGES
frange = [ Min( [oB_fgm_psd['DEPEND_0'].min, oB_scm_psd['DEPEND_0'].min, oB_fsm_psd.min] ), $
           Max( [oB_fgm_psd['DEPEND_0'].max, oB_scm_psd['DEPEND_0'].max, oB_fsm_psd.max] ) ]
frange[0] = 1e-2
bpsdrange = [ Min( [oB_fgm_psd.min, oB_scm_psd.min, oB_fsm_psd.min] ), $
              Max( [oB_fgm_psd.max, oB_scm_psd.max, oB_fsm_psd.max] ) ]

;TITLE
trange = MrVar_GetTRange()
ttitle = StrMid(StrJoin(StrSplit(trange[0], 'T', /EXTRACT), ' '), 0, 19) + '-' + StrMid(trange[1], 11, 8)
title = String(ttitle, t, 0.5*100.0, taper, tf_zeropad, FORMAT='(%"%s!CT=%0.1fs Overlap=%i\% Win=%s Zeropad=%i")')
;title  = ttitle+'!CT='+String(t, FORMAT='(f0.1)')+'s Overlap=50%, Win='+taper + ' Zeropad='+(tf_zeropad ? 'True' : 'False')

;Set up the window
w1 = MrWindow( LAYOUT   = [3,1], $
               NAME     = 'fsm-compare-b', $
               OYMARGIN = [4,4], $
               REFRESH  = 0, $
               XGAP     = 12, $
               XSIZE    = 900, $
               YSIZE    = 350 )

;Plot B-Field
p1a = MrVar_Plot( oB_scm_psd[*,0], COLOR='Red', /CURRENT, XRANGE=frange, XTICKFORMAT='logtickformat', YRANGE=bpsdrange, YTITLE='Bx PSD!C(nT$\up2$/Hz)' )
p1b = MrVar_Plot( oB_fgm_psd[*,0], COLOR='Blue', OVERPLOT=p1a )
p1c = MrVar_Plot( oB_fsm_psd[*,0], COLOR='Black', OVERPLOT=p1a )
;p1d = MrVar_Plot( opsd_fgm[*,0], COLOR='Magenta', OVERPLOT=p1a)
p1d = MrVar_Plot( oNoise_fgm, COLOR='Blue', LINESTYLE='--', OVERPLOT=p1a )
p1e = MrVar_Plot( oNoise_scm, COLOR='Red',  LINESTYLE='--', OVERPLOT=p1a )
l1 = MrLegend( ALIGNMENT    = 'NE', $
               FILL_COLOR   = '', $
               LABEL        = ['FGM', 'SCM', 'FSM'], $
               LINESTYLE    = 'None', $
               POSITION     = [1.0,1.0], $
               /RELATIVE, $
               SAMPLE_WIDTH = 0.0, $
               TARGET       = p1a, $
               TEXT_COLOR   = ['Blue', 'Red', 'Black'] )

p2a = MrVar_Plot( oB_fgm_psd[*,1], COLOR='Blue', /CURRENT, TITLE=title, XRANGE=frange, XTICKFORMAT='logtickformat', YRANGE=bpsdrange, YTITLE='By PSD!C(nT$\up2$/Hz)' )
p2b = MrVar_Plot( oB_scm_psd[*,1], COLOR='Red', OVERPLOT=p2a )
p2c = MrVar_Plot( oB_fsm_psd[*,1], COLOR='Black', OVERPLOT=p2a )
p2d = MrVar_Plot( oNoise_fgm, COLOR='Blue', LINESTYLE='--', OVERPLOT=p2a )
p2e = MrVar_Plot( oNoise_scm, COLOR='Red',  LINESTYLE='--', OVERPLOT=p2a )
l2 = MrLegend( ALIGNMENT    = 'NE', $
               FILL_COLOR   = '', $
               LABEL        = ['FGM', 'SCM', 'FSM'], $
               LINESTYLE    = 'None', $
               POSITION     = [1.0,1.0], $
               /RELATIVE, $
               SAMPLE_WIDTH = 0.0, $
               TARGET       = p2a, $
               TEXT_COLOR   = ['Blue', 'Red', 'Black'] )

p3a = MrVar_Plot( oB_fgm_psd[*,2], COLOR='Blue', /CURRENT, XRANGE=frange, XTICKFORMAT='logtickformat', YRANGE=bpsdrange, YTITLE='Bz PSD!C(nT$\up2$/Hz)' )
p3b = MrVar_Plot( oB_scm_psd[*,2], COLOR='Red', OVERPLOT=p3a )
p3c = MrVar_Plot( oB_fsm_psd[*,2], COLOR='Black', OVERPLOT=p3a )
p3d = MrVar_Plot( oNoise_fgm, COLOR='Blue', LINESTYLE='--', OVERPLOT=p3a )
p3e = MrVar_Plot( oNoise_scm, COLOR='Red',  LINESTYLE='--', OVERPLOT=p3a )
l3 = MrLegend( ALIGNMENT    = 'NE', $
               FILL_COLOR   = '', $
               LABEL        = ['FGM', 'SCM', 'FSM'], $
               LINESTYLE    = 'None', $
               POSITION     = [1.0,1.0], $
               /RELATIVE, $
               SAMPLE_WIDTH = 0.0, $
               TARGET       = p3a, $
               TEXT_COLOR   = ['Blue', 'Red', 'Black'] )


w1[0] -> SetLayout, [1,1]
w1 -> TrimLayout
w1 -> Refresh

;-------------------------------------------
; Plot E PSD ///////////////////////////////
;-------------------------------------------

IF tf_edp THEN BEGIN
	epsdrange = [ oE_edp_psd.min, oE_edp_psd.max ]
	
	;Set up the window
	w2 = MrWindow( LAYOUT   = [3,1], $
	               NAME     = 'fsm-compare-e', $
	               OYMARGIN = [4,4], $
	               REFRESH  = 0, $
	               XGAP     = 12, $
	               XSIZE    = 900, $
	               YSIZE    = 350 )
	
	;Plot E-Field
	p4 = MrVar_Plot( oE_edp_psd[*,0], /CURRENT, XRANGE=frange, XTICKFORMAT='logtickformat', YRANGE=epsdrange, YTITLE='Ex PSD!C(mV$\up2$/Hz)' )
	p5 = MrVar_Plot( oE_edp_psd[*,1], /CURRENT, XRANGE=frange, XTICKFORMAT='logtickformat', YRANGE=epsdrange, YTITLE='Ey PSD!C(mV$\up2$/Hz)', TITLE=title )
	p6 = MrVar_Plot( oE_edp_psd[*,2], /CURRENT, XRANGE=frange, XTICKFORMAT='logtickformat', YRANGE=epsdrange, YTITLE='Ez PSD!C(mV$\up2$/Hz)' )
	
	w2[0] -> SetLayout, [1,1]
	w2 -> TrimLayout
	w2 -> Refresh
ENDIF

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
	
	foptdesc = String(taper, t, overlap*100, tf_zeropad, FORMAT='(%"%s-%is-%i\%Overlap-%iZero")')
	
	;File name
	fname1 = StrJoin( [sc, 'fsm', mode, level, 'fgm-scm-spectra-nwss-'+foptdesc], '_' )
	fname2 = StrJoin( [sc, 'fsm', mode, level, 'edp-spectra-'+foptdesc], '_' )
	fname1 = FilePath( fname1, ROOT_DIR=output_dir )
	fname2 = FilePath( fname2, ROOT_DIR=output_dir )
	
	;Save the figure
	fout1 = MrVar_PlotTS_Save( w1, fname1, output_ext )
	IF tf_edp THEN fout2 = MrVar_PlotTS_Save( w2, fname2, output_ext )
ENDIF

END