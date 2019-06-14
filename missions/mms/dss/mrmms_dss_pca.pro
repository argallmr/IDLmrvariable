; docformat = 'rst'
;
; NAME:
;       mrmms_dss_pca
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
;   Perform a spin-epoch principal component analysis.
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
;       2018/08/12  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Plot the results of the principal component analysis.
;-
FUNCTION MrMMS_DSS_PCA_Test, oVar, oVar_pca, oPCA
	Compile_Opt idl2
	On_Error, 2
	
	oDelta = oVar - oVar_pca
	oDelta['COLOR'] = 'Blue'
	oDelta['LABEL'] = '$\Delta$'
	
	;Plot the variable and its reconstruction
	win = MrVar_PlotTS( [oVar, oVar_pca], /NO_REFRESH)
	win = MrVar_OPlotTS( win[oVar.name], oDelta)
	
	;Mark the start of each spin
	oGfx     = win[oVar_pca.name]
	oT_Epoch = oPCA['DEPEND_0']
	oT_Spin  = oT_Epoch['DEPEND_0']
	FOR i = 0, N_Elements(oT_Spin) - 1 DO BEGIN
		ps = MrPlotS( Replicate(oT_Spin['DATA', i, 'SSM'], 2), oGfx.yrange, $
		              COLOR     = 'Red', $
		              LINESTYLE = '--', $
		              TARGET    = oGfx )
	ENDFOR
	
	;Realize the graph
	win[oVar_pca.name] -> SetProperty, YTITLE='PCA'
	win[0] -> SetLayout, [1,1]
	win -> TrimLayout
	win.oxmargin = [10, 8]
	win -> Refresh

;-------------------------------------------
; Plot PCA Results /////////////////////////
;-------------------------------------------
	dims     = Size(oPCA, /DIMENSIONS)
	npspin   = dims[0]
	nSpins   = dims[1]
	nSpinFit = dims[2]
	
	;Plot
	w1 = MrWindow( LAYOUT   = [1, nSpinFit], $
	               OXMARGIN = [10,8], $
		           REFRESH  = 0, $
		           XSIZE    = 500, $
		           YGAP     = 0.5, $
		           YSIZE    = 700 )
	
	colors = MrDefaultColor(NCOLORS=nSpins)
	FOR i = 0, nSpinFit-1 DO BEGIN
		;Least-Squares fit
	;	lad = LADFit(Reform(v_temp[i,*]), Reform(v_out[i,*]))
		
		;Principal Components
		p1 = MrPlot( oT_Epoch['DATA', 'SSM'], oPCA['DATA',*,*,i], $
		             COLOR       = colors, $
		             /CURRENT, $
		             DIMENSION   = 1, $
		             TITLE       = (i EQ 0 ? 'Principal Components' : ''), $
		             XTICKFORMAT = (i EQ nSpinFit-1 ? '' : '(a1)'), $
		             XTITLE      = (i EQ nSpinFit-1 ? 'Spin Epoch (s)' : ''), $
		             YTITLE      = 'PC' + String(i+1, FORMAT='(i0)') )
	ENDFOR
	
	l1 = MrLegend( ALIGNMENT    = 'NW', $
	               COLOR        = '', $
	               FILL_COLOR   = '', $
	               LABEL        = 'Spin' + String(LIndGen(nSpins)+1, FORMAT='(i0)'), $
	               POSITION     = [1.0, 1.0], $
	               /RELATIVE, $
	               SAMPLE_WIDTH = 0, $
	               TARGET       = w1[0], $
	               TEXT_COLOR   = colors, $
	               TEXT_SIZE    = 1.0 )
	
;;-------------------------------------------
;; Plot Spin Epoch Results //////////////////
;;-------------------------------------------
;	
;	w2 = MrWindow( LAYOUT   = [1,3], $
;	               OXMARGIN = [10,8], $
;	               REFRESH  = 0, $
;	               XSIZE    = 500, $
;	               YGAP     = 0.5, $
;	               YSIZE    = 600 )
;	
;	;Original data
;	p2 = MrPlot( t_epoch, v_out, $
;	             COLOR       = colors, $
;	             /CURRENT, $
;	             DIMENSION   = 2, $
;	             TITLE       = '', $
;	             XTICKFORMAT = '(a1)', $
;	             YTITLE      = 'SCPot!C(V)' )
;	
;	;Reconstructed from first N principal components
;	p3 = MrPlot( t_epoch, v_pca, $
;	             COLOR       = colors, $
;	             /CURRENT, $
;	             DIMENSION   = 2, $
;	             XTICKFORMAT = '(a1)', $
;	             XTITLE      = '', $
;	             YTITLE      = 'V$\downPCA$!C(V)' )
;	
;	;Reconstruction from frist N principal components + correction for /STANDARDIZE
;	p4 = MrPlot( t_epoch, v_out - v_pca, $
;	             COLOR     = colors, $
;	             /CURRENT, $
;	             DIMENSION = 2, $
;	             XTITLE    = 'Spin Epoch (s)', $
;	             YTITLE    = 'V-V$\downPCA$!C(V)' )
;	
;	
;	l2 = MrLegend( ALIGNMENT    = 'NW', $
;	               COLOR        = '', $
;	               FILL_COLOR   = '', $
;	               LABEL        = 'Spin' + String(LIndGen(nSpins)+1, FORMAT='(i0)'), $
;	               POSITION     = [1.0, 1.0], $
;	               /RELATIVE, $
;	               SAMPLE_WIDTH = 0, $
;	               TARGET       = p2, $
;	               TEXT_COLOR   = colors, $
;	               TEXT_SIZE    = 1.0 )
	
;-------------------------------------------
; Finish ///////////////////////////////////
;-------------------------------------------
	
	w1[0] -> SetLayout, [1,1]
	w1 -> TrimLayout
	w1 -> Refresh

	RETURN, [win, w1]
END


;+
;   Perform a spin-epoch principal component analysis. The start of each spin
;   is identified by the sun pulse from the digital sun sensor. Data is broken
;   into spins and a principal component analysis is performed over a set of
;   spins. The highest ranked principal components are transformed back into
;   data coordinates and are considered the spin-related noise component.
;
; :Params:
;    SC:            in, required, type=string
;                   MMS spacecraft identifier.
;    VARIABLE:      in, required, type=int/string/objref
;                   The name, number, or objref of a MrScalarTS variable for which
;                       the PCA analysis is performed.
;    NFIT:          in, optional, type=int, default=9
;                   Number of spins per PCA.
;
; :Keywords:
;    CACHE:         in, optional, type=boolean, default=0
;                   If set, the output will be stored in the variable cache.
;    MINVAR:        in, optional, type=float, default=2.0
;                   The minimum percentage of the total variance that still classifies
;                       a principal component as spin-related noise. By default, any
;                       principal component that accounts for 2% or more of the total
;                       variance is considered to be caused by spin-related noise.
;    NAME:          in, optional, type=string, default='SpinEpochPCA(' + variable.name + ')'
;                   Name to be given to the output variable.
;    NO_LOAD:       in, optional, type=boolean, default=0
;                   If zero, DSS sun pulse data will be loaded from file into the variable
;                       cache. If set, assume the data has already been loaded.
;    PCTVAR:        in, optional, type=float, default=85.0
;                   The cumulative percent variance of the principal components that
;                       is caused by spin-related noise. Both `MINVAR` and `PCTVAR`
;                       are considered when determining which principal components are
;                       noise related. When summing the variance of each principal
;                       component in order from maximum to minumum variance, only the
;                       first N components for which the percent variance is less than
;                       85% of total variance (and which have a percent variance greater
;                       than `MINVAR`) are selected as contributions to spin-related
;                       noise.
;    TEST:          in, optional, type=boolean, default=0
;                   If set, the results will be plotted and a MrGraphics window will be
;                       returned.
;    TRANGE:        in, optional, type=strarr(2)
;                   The start and end times of the data interval, in ISO_8601 format.
;
; :Returns:
;    oPCA:          out, required, type=objref
;                   A MrScalarTS object with the spin-related principal components
;                       transformed into data space.
;-
FUNCTION MrMMS_DSS_PCA, sc, variable, nFit, $
CACHE=cache, $
MINVAR=minvar, $
NAME=name, $
NO_LOAD=no_load, $
PCTVAR=pctvar, $
TEST=test, $
TRANGE=trange, $
VERBOSE=verbose
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		Return, !Null
	ENDIF
	
	oVar    = MrVar_Get(variable)
	tf_load = ~Keyword_Set(no_load)
	tf_test = Keyword_Set(test)
	tf_verbose = Keyword_Set(verbose)
	IF N_Elements(nSpinFit) EQ 0 THEN nSpinFit = 9
	IF N_Elements(minvar)   EQ 0 THEN minvar   = 2.0
	IF N_Elements(pctvar)   EQ 0 THEN pctvar   = 85.0
	IF N_Elements(trange)   GT 0 THEN MrVar_SetTRange, trange
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	
	sp_vname = StrJoin( [sc, '101', 'sunpulse'], '_' )
	
	pca_vname      = 'SpinEpochPCA('          + oVar.name + ')'
	eigval_vname   = 'SpinEpochEigenvalues('  + oVar.name + ')'
	eigvec_vname   = 'SpinEpochEigenvectors(' + oVar.name + ')'
	variance_vname = 'SpinEpochVariances('    + oVar.name + ')'
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	
	IF tf_load THEN BEGIN
		MrMMS_Load_Data, sc, 'fields', 'hk', 'l1b', $
		                 /HK, $
		                 OPTDESC   = '101', $
		                 VARFORMAT = '*sunpulse'
	ENDIF
	
;-------------------------------------------
; Filter Sunpulse Times ////////////////////
;-------------------------------------------
	
	;Check for triple pulse
	;oSP   = MrVar_Get(sp_vname)
	;t_ssm = oSP['DATA', 'SSM']
	;dt    = t_ssm[1:-1] - t_ssm[0:-2]
	;count = 0
	;FOR i = 0, N_Elements(dt)-3, 3 DO IF Array_Equal(dt[i:i+3], [0,0,0]) THEN count += 1
	
	;A sunpulse is indicated by a double pulse. A 1PPS pulse is indicated by a triple pulse.
	;   - Assume all double pulses. Remove duplicates
	oSP    = MrVar_Get(sp_vname)
	iUniq  = oSP -> Uniq()
	oSP    = oSP[iUniq]
	
	;Spin period
	t_ssm_spin = oSP['DATA','SSM']
	dt_spin    = t_ssm_spin[1:-1] - t_ssm_spin[0:-2]
	dt_spin    = [dt_spin, dt_spin[-1]]
	T_spin     = Mean(dt_spin)
	
;-------------------------------------------
; Spin-Epoch Determination /////////////////
;-------------------------------------------
	
	;Sample interval
	;   - Ensure that the sampling rate, 1/dt, is an integer
	dt   = oVar['TIMEVAR'] -> GetSI()
	dt   = 1.0/Round(1.0/dt)
	
	;Locate data within each spin
	t_ssm = (oVar['TIMEVAR'])['DATA','SSM']
	iSpin = Value_Locate(t_ssm_spin, t_ssm)
	iGood = Where( (iSpin GE 0) AND (t_ssm LT t_ssm_spin[-1]+T_spin), nGood )
	IF nGood GT 0 THEN BEGIN
		iSpin = iSpin[iGood]
		t_ssm = t_ssm[iGood]
		oVar  = oVar[iGood]
	ENDIf
	
	;Find the indices associated with each spin
	nHist  = Histogram(iSpin, MIN=iSpin[0], MAX=iSpin[-1], REVERSE_INDICES=ri)
	nSpins = N_Elements(nHist)
	
	;Loop through each spin
	npspin  = Floor(T_spin/dt)+1     ;Number of samples per spin
	it_out  = LIndGen(npspin)        ;Common time grid for each spin
	pVar    = oVar['PTR']            ;Data
	v_out   = FltArr(nSpins,npspin)
	t_out   = FltArr(nSpins,npspin)
	t_epoch = LIndGen(npspin) * dt
	FOR i = 0, nSpins - 1 DO BEGIN
		IF ri[i] EQ ri[i+1] THEN BEGIN
			MrPrintF, 'LogWarn', 'No points in spin ', i+1, FORMAT='(a0, i0)'
			CONTINUE
		ENDIF
		IF (dt_spin[i]-T_spin)/T_spin GT 0.01 THEN MrPrintF, 'LogWarn', 'Spin period changed by >1%.'
		
		i0 = ri[ri[i]]
		i1 = ri[ri[i+1]-1]
		
		;Interpolate onto common time grid
		t_out[i,*] = t_epoch + t_ssm_spin[i]
		v_out[i,*] = Interpol((*pVar)[i0:i1], t_ssm[i0:i1], t_out[i,*])
	ENDFOR
	
;-------------------------------------------
; Principal Component Analysis /////////////
;-------------------------------------------
	pca     = FltArr(nSpins, nSpinFit, npspin)
	v_pca   = FltArr(nSpins, npspin)
	eigvals = FltArr(nSpins, nSpinFit)
	eigvecs = FltArr(nSpins, nSpinFit, nSpinFit)
	sigma   = FltArr(nSpins, nSpinFit)
	FOR i = 0, nSpins - nSpinFit DO BEGIN
		n0  = Floor(nSpinFit/2)
		avg = Total(v_out[i:i+nSpinFit-1,*], 2)/npspin
		x    = v_out[i:i+nSpinFit-1,*] - Rebin(avg, nSpinFit, npspin)
		temp = PComp( x, $
		              COEFFICIENTS = coeffs, $
		              /COVARIANCE, $
		              /DOUBLE, $
		              EIGENVALUES  = eigenvalues, $
		              /STANDARDIZE, $
		              VARIANCES    = variances )
		
		;Determine which principal components are to be removed from the data
		iMin  = Where(variances GT minvar/100.0, nMin)
		iPct  = Where( Total(variances[iMin], /CUMULATIVE) GE pctvar/100.0, nPct)
		iKeep = LIndGen(iMin[iPct[0]]+1)
		
		;Eigenvectors
		;   - The eigen vectors are scaled such that the sum of the squares
		;     is equal to the corresponding eigenvalue
		;   - eigvecs[*,0] corresponds to the first principal component
		eigenvectors = coeffs / Rebin(eigenvalues, nSpinFit, nSpinFit)
		
		;Reconstructed result
		v_temp = temp[iKeep,*] ## eigenvectors[*,iKeep]
		
		;Reconstructed error
		err = Total( (temp##eigenvectors  - x)^2.0 )
		
		;Percent variance
		IF tf_verbose THEN BEGIN
			PRINT, nSpins - nSpinFit, FORMAT='(%"Fit #%i")'
			PRINT, 'Mode     EigVal    %Var      CumVar'
			FOR j = 0, nSpinFit - 1 $
				DO Print, j+1, eigenvalues[j], variances[j]*100.0, Total(variances[0:j])*100, $
				          FORMAT='(i2, 3x, f10.4, 3x, f7.2, 3x, f7.2)'
			PRINT, ''
			PRINT, 'Coefficients:'
			FOR j = 0, nSpinFit-1 $
				DO Print, j+1, coeffs[*,j], FORMAT='("MODE#", i02, '+StrTrim(nSpinFit, 2)+'(f10.4))'
			PRINT, ''
			PRINT, 'Reconstruction Error: ', err
			PRINT, '------------------------------'
			PRINT, ''
		ENDIF
	
		IF i EQ 0 THEN BEGIN
			pca[i:i+n0,*,*]     = Rebin(Reform(Temporary(temp), 1, nSpinFit, npspin), n0+1, nSpinFit, npspin)
			eigvals[i:i+n0,*]   = Rebin(Reform(Temporary(eigenvalues), 1, nSpinFit), n0+1, nSpinFit)
			eigvecs[i:i+n0,*,*] = Rebin(Reform(Temporary(eigenvectors), 1, nSpinFit, nSpinFit), n0+1, nSpinFit, nSpinFit)
			sigma[i:i+n0,*]     = Rebin(Reform(Temporary(variances), 1, nSpinFit), n0+1, nSpinFit)
			v_pca[i:i+n0, *]    = (Temporary(v_temp))[0:n0,*]
		ENDIF ELSE IF i EQ nSpins - nSpinFit THEN BEGIN
			pca[i+n0:*, *, *]   = Rebin(Reform(Temporary(temp), 1, nSpinFit, npspin), n0+1, nSpinFit, npspin)
			eigvals[i+n0:*, *]  = Rebin(Reform(Temporary(eigenvalues), 1, nSpinFit), n0+1, nSpinFit)
			eigvecs[i+n0:*,*,*] = Rebin(Reform(Temporary(eigenvectors), 1, nSpinFit, nSpinFit), n0+1, nSpinFit, nSpinFit)
			sigma[i+n0:*, *]    = Rebin(Reform(Temporary(variances), 1, nSpinFit), n0+1, nSpinFit)
			v_pca[i+n0:*, *]    = (Temporary(v_temp))[n0:*, *]
		ENDIF ELSE BEGIN
			pca[i+n0,*,*]     = Temporary(temp)
			eigvals[i+n0,*]   = Temporary(eigenvalues)
			eigvecs[i+n0,*,*] = Temporary(eigenvectors)
			sigma[i+n0,*]     = Temporary(variances)
			v_pca[i+n0, *]    = (Temporary(v_temp))[n0, *]
		ENDELSE
	ENDFOR
	
	strPC   = String(LIndGen(nSpinFit), FORMAT='(i0)')
	strSpin = String(LIndGen(nSpins), FORMAT='(i0)')
	
	;Spin Epoch
	trange   = MrVar_GetTRange()
	oT_epoch = MrTimeVar(t_epoch, 'SSM', T_REF=trange[0])
	oT_epoch['CATDESC']    = 'Spin epoch times'
	oT_epoch['DEPEND_0']   = MrTimeVar(t_ssm_spin, 'SSM', T_REF=trange[0])
	oT_epoch['LABL_PTR_1'] = strPC
	
	;Principal Components
	oPCA = MrTimeSeries( oT_epoch, Transpose(pca, [2, 0, 1]), $
	                     CACHE = cache, $
	                     NAME  = pca_vname )
	oPCA['CATDESC']    = 'Principal components'
	oPCA['LABL_PTR_1'] = 'Spin' + strSpin
	oPCA['LABL_PTR_2'] = 'PC' + strPC
	
	;Eigenvalues
	oEigVals = MrVariable( eigvals, $
	                       CACHE = cache, $
	                       NAME  = eigval_vname )
	oEigVals['CATDESC']    = 'Eigenvalues of the spin-epoch covariance matrix.'
	oEigVals['DIMENSION']  = 1
	oEigVals['LABL_PTR_1'] = StrSpin
	oEigVals['LABL_PTR_2'] = StrPC
	
	;Eigenvectors
	oEigVecs = MrVariable( eigvecs, $
	                       CACHE = cache, $
	                       NAME  = eigvec_vname )
	oEigVecs['CATDESC']    = 'Eigenvectors of the spin-epoch covariance matrix. Matrix ' + $
	                         'is organized as [Spin, Component, Eigenvector].'
	oEigVecs['DIMENSION']  = 1
	oEigVecs['LABL_PTR_1'] = StrSpin
	oEigVecs['LABL_PTR_2'] = StrPC
	oEigVecs['LABL_PTR_3'] = StrPC
	
	;Variances
	oSigma = MrVariable( sigma, $
	                     CACHE = cache, $
	                     NAME  = variance_vname )
	oSigma['CATDESC']    = 'Variances associated with each principal component.'
	oSigma['LABL_PTR_1'] = StrSpin
	oSigma['LABL_PTR_2'] = StrPC
	
;-------------------------------------------
; Apply Spin Epoch Analysis ////////////////
;-------------------------------------------
	
	;Spin epoch variable
	oV_pca = MrScalarTS( Reform(Transpose(t_out), nSpins*npspin), Reform(Transpose(v_pca), nSpins*npspin), $
	                     T_TYPE = 'SSM', $
	                     T_REF = trange[0], $
	                     CACHE = cache, $
	                     NAME  = pca_vname )
	oV_pca['CATDESC'] = 'Principal components transformed into data coordinates. Only ' + $
	                    'those components that meet the minimum percentage of the total ' + $
	                    'variance and the cumulative percent variance conditions are kept.'
	
;-------------------------------------------
; Test? ////////////////////////////////////
;-------------------------------------------
	IF tf_test $
		THEN RETURN, MrMMS_DSS_PCA_Test( oVar, oV_pca, oPCA ) $
		ELSE RETURN, oV_pca
END

