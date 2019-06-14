; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_4sc_FPI_KEntropy_NDep
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
;+
;   Calculate moments of the distribution function and plot them against the official
;   FPI L2 dataset. The moments calculation takes into account the FPI internal photo-
;   electron model, but the method of integration is different.
;
;       1. Bxyz, |B|
;       2. density
;       3. Entropy: S = Integral{ f ln(f) } / n
;       4. Entropy per particle: S / n
;       5. Entropy: Sb = P/n^(5/3)
;       6. Scalar Pressure
;       7. Temperature Pressure
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;       SPECIES:    in, required, type=string
;                   Particle species. Options are {'e' | 'i'}
;
; :Keywords:
;       COORDS:     in, optional, type=string, default='gse'
;                   Coordinate system in which to load the data. Options are: {'dbcs' | 'gse' | 'gsm'}
;       FGM_INSTR:  in, optional, type=string, default='fgm'
;                   The FGM instrument to use. Options are: {'afg' | 'dfg' | 'fgm'}
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'ql' | 'l2'}
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source files.
;       OUTPUT_DIR: in, optional, type=string, default='~/figures/'
;                   Directory in which to save the figure. If neither `OUTPUT_DIR` or
;                       `OUTPUT_EXT` are given, no file is made.
;       OUTPUT_EXT: in, optional, type=string/strarr, default='png'
;                   Extension (and file type) of the figure. Options include
;                       'png', 'jpeg', 'tiff', 'ps', 'eps', 'pdf'. If neither
;                       `OUTPUT_DIR` or `OUTPUT_EXT` are given, no file is made.
;
; :Categories:
;    MMS
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
;       2018/06/09  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_4sc_FPI_KEntropy_NDep, mode, species, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
FLUX=flux, $
LEVEL=level, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
PER_PARTICLE=per_particle
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, Obj_New()
	ENDIF
	
	tf_load         = ~Keyword_Set(no_load)
	tf_per_particle = Keyword_Set(per_particle)
	tf_flux         = Keyword_Set(flux)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	
	m_e  = MrConstants('m_e')
	m_i  = MrConstants('m_H')
	q    = MrConstants('q')
	kB   = MrConstants('k_b')
	J2eV = MrConstants('J2eV')
	
;-------------------------------------------
; Variable Parameters //////////////////////
;-------------------------------------------

	sc  = 'mms' + ['1', '2', '3', '4']
	isc = [0, 1, 2, 3]
	nSC = N_Elements(sc)
	
	fpi_instr = 'd' + species + 's'

;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, mode, $
		                     INSTR  = fgm_instr, $
		                     LEVEL  = level, $
		                     VARFORMAT = '*b_' + coords + '_*', $
		                     SUFFIX = suffix
		
	ENDIF

	;Entropy
	MrMMS_FPI_Load_Entropy, sc, mode, species, $
	                        NO_LOAD  = no_load, $
	                        VARNAMES = s_varnames

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	
	;FGM
	b_vnames    = sc + '_' + StrJoin( [fgm_instr, 'b', coords, mode, level], '_' )
	bvec_vnames = sc + '_' + StrJoin( [fgm_instr, 'bvec', coords, mode, level], '_' )
	bmag_vnames = sc + '_' + StrJoin( [fgm_instr, 'bmag', coords, mode, level], '_' )
	bbc_vname   = StrJoin( ['mms', fgm_instr, 'b', coords, 'barycenter', mode, level], '_' )
	btbc_vname  = StrJoin( ['mms', fgm_instr, 'b', 'mag',  'barycenter', mode, level], '_' )
	
	CASE 1 OF
		tf_flux && tf_per_particle: BEGIN
			n_vnames     = s_varnames[isc*18+3]
			s_vnames     = s_varnames[isc*18+4+3]
			sb_vnames    = s_varnames[isc*18+8+3]
			mbar_vnames  = s_varnames[isc*18+12+3]
			gradn_vname  = s_varnames[-7]
			grads_vname  = s_varnames[-5]
			gradsb_vname = s_varnames[-3]
			gradm_vname  = s_varnames[-1]
		ENDCASE
		tf_flux: BEGIN
			n_vnames     = s_varnames[isc*18+3]
			s_vnames     = s_varnames[isc*18+4+1]
			sb_vnames    = s_varnames[isc*18+8+1]
			mbar_vnames  = s_varnames[isc*18+12+1]
			gradn_vname  = s_varnames[-7]
			grads_vname  = s_varnames[-6]
			gradsb_vname = s_varnames[-4]
			gradm_vname  = s_varnames[-2]
		ENDCASE
		tf_per_particle: BEGIN
			n_vnames     = s_varnames[isc*18]
			s_vnames     = s_varnames[isc*18+4+2]
			sb_vnames    = s_varnames[isc*18+8+2]
			mbar_vnames  = s_varnames[isc*18+12+2]
			gradn_vname  = s_varnames[-14]
			grads_vname  = s_varnames[-12]
			gradsb_vname = s_varnames[-10]
			gradm_vname  = s_varnames[-8]
		ENDCASE
		ELSE: BEGIN
			n_vnames     = s_varnames[isc*18]
			s_vnames     = s_varnames[isc*18+4]
			sb_vnames    = s_varnames[isc*18+8]
			mbar_vnames  = s_varnames[isc*18+12]
			gradn_vname  = s_varnames[-14]
			grads_vname  = s_varnames[-13]
			gradsb_vname = s_varnames[-11]
			gradm_vname  = s_varnames[-9]
		ENDCASE
	ENDCASE
	
;-------------------------------------------
; Barycenter ///////////////////////////////
;-------------------------------------------
	;B
	oB    = ObjArr(4)
	oB[0] = MrVar_Get(bvec_vnames[0])
	oT    = (oB[0])['TIMEVAR']
	FOR i = 1, nSC - 1 DO BEGIN
		oBB   = MrVar_Get(bvec_vnames[i])
		oB[i] = oBB -> Interpol(oT)
	ENDFOR
	
	oBbc = (oB[0] + oB[1] + oB[2] + oB[3]) / 4.0
	oBbc -> SetName, bbc_vname
	oBbc -> Cache
	
	oBbc['COLOR']         = ['Blue', 'Forest Green', 'Red']
	oBbc['LABEL']         = 'B$\down' + ['X', 'Y', 'Z'] + '$'
	oBbc['PLOT_TITLE']    = StrUpCase( StrJoin( ['mms', fpi_instr, mode, level], ' ' ) )
	oBbc['TITLE']         = 'B$\downBC$!C(nT)'
	oBbc['UNITS']         = 'nT'
	oBbc['SI_CONVERSION'] = '1e-9>T'
	
	oBtbc = oBbc -> Magnitude(/CACHE, NAME=btbc_vname)
	oBtbc['COLOR']         = 'Black'
	oBtbc['LABEL']         = '|B|'
	oBtbc['TITLE']         = 'B$\downBC$!C(nT)'
	oBtbc['UNITS']         = 'nT'
	oBtbc['SI_CONVERSION'] = '1e-9>T'
	
;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	FOR i = 0, nSC - 1 DO BEGIN
		oN  = MrVar_Get(n_vnames[i])
		oSb = MrVar_Get(sb_vnames[i])
		oM  = MrVar_Get(mbar_vnames[i])
		
		IF oN  -> HasAttr('LABEL') THEN oN  -> RemoveAttr, 'LABEL'
		IF oSb -> HasAttr('LABEL') THEN oSb -> RemoveAttr, 'LABEL'
		IF oM  -> HasAttr('LABEL') THEN oM  -> RemoveAttr, 'LABEL'
	ENDFOR
	
;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------

	win = MrVar_PlotTS( [ bbc_vname, n_vnames[0], s_vnames[0], sb_vnames[0], mbar_vnames[0], $
	                      gradn_vname, grads_vname, gradsb_vname, gradm_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )
	
	win = MrVar_OPlotTS( bbc_vname,      btbc_vname )
	win = MrVar_OPlotTS( n_vnames[0],    n_vnames[1:3] )
	win = MrVar_OPlotTS( s_vnames[0],    s_vnames[1:3] )
	win = MrVar_OPlotTS( sb_vnames[0],   sb_vnames[1:3] )
	win = MrVar_OPlotTS( mbar_vnames[0], mbar_vnames[1:3] )

	win[0] -> SetLayout, [1,1]
	win -> TrimLayout
	win.oxmargin = [15,10]

;-------------------------------------------
; Save the Figure //////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
			MrPrintF, 'LogText', 'Saving file to: "' + output_dir + '".'
		ENDIF
		
		CASE 1 OF
			tf_per_particle && tf_flux: optdesc = 'sfn'
			tf_per_particle:            optdesc = 'sn'
			tf_flux:                    optdesc = 'sf'
			ELSE:                       optdesc = 's'
		ENDCASE
		
		;Save the figure
		fname = StrJoin(['mms', fpi_instr, mode, level, 'entropy-'+optdesc+'-ndep'], '_')
		fname = FilePath(fname, ROOT_DIR=output_dir)
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF
	
;-------------------------------------------
; Finish ///////////////////////////////////
;-------------------------------------------

	win -> Refresh
	RETURN, win
END