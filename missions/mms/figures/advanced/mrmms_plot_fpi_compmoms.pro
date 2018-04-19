; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FPI_CompMoms
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
;+
;   Calculate moments of the distribution function and plot them against the official
;   FPI L2 dataset. The moments calculation takes into account the FPI internal photo-
;   electron model, but the method of integration is different.
;
;       1.  Bxyz, |B|             1. Bxyz, |B|
;       2.  Density               2. Pxx
;       3.  Vx                    3. Pxy
;       4.  Vy                    4. Pxz
;       5.  Vz                    5. Pyy
;       6.  |V|                   6. Pyz
;       7.  Qx                    7. Pzz
;       8.  Qy                    8. p
;       9.  Qz
;       10. |Q|
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
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source files.
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
;       2017/02/07  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_FPI_CompMoms, sc, mode, species, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
NO_LOAD=no_load
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, Obj_New()
	ENDIF
	
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(species)   EQ 0 THEN species   = 'e'
	instr   = 'd' + species + 's'

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source B-Field
	b_vname    = StrJoin( [sc, fgm_instr, 'b',    coords, mode, level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec', coords, mode, level], '_' )
	bmag_vname = StrJoin( [sc, fgm_instr, 'bmag', coords, mode, level], '_' )
	
	;Source Moments
	scpot_vname = StrJoin( [sc, 'edp', 'scpot', 'fast', level], '_' )
	n_vname     = StrJoin( [sc, instr, 'numberdensity',        mode], '_')
	v_vname     = StrJoin( [sc, instr, 'bulkv',      coords, mode], '_')
	p_vname     = StrJoin( [sc, instr, 'prestensor', coords, mode], '_')
	t_vname     = StrJoin( [sc, instr, 'temptensor', coords, mode], '_')
	q_vname     = StrJoin( [sc, instr, 'heatq',      coords, mode], '_')
	
	;Source Distribution
	f_vname      = StrJoin( [sc, instr, 'dist',                  mode], '_')
	dphi_vname   = StrJoin( [sc, instr, 'startdelphi', 'count',  mode], '_' )
	parity_vname = StrJoin( [sc, instr, 'steptable',   'parity', mode], '_' )
	
	;Source Photoelectron model
	ph_dphi_vname = StrJoin( ['mms', instr, 'startdelphi', 'counts', mode], '_' ) + '_model'
	ph_f0_vname   = StrJoin( ['mms', instr, 'bgdist',      'p0',     mode], '_' ) + '_model'
	ph_f1_vname   = StrJoin( ['mms', instr, 'bgdist',      'p1',     mode], '_' ) + '_model'
	ph_e0_vname   = StrJoin( ['mms', instr, 'energy0',               mode], '_' ) + '_model'
	ph_e1_vname   = StrJoin( ['mms', instr, 'energy1',               mode], '_' ) + '_model'
	ph_scl_vname  = StrJoin( [sc,    instr, 'scl',         'model',  mode], '_' ) + '_model'
	
	;Derived names
	dist_vname   = StrJoin([sc, instr, 'dist4d',        mode], '_')
	fph_vname    = StrJoin([sc, instr, 'dist', 'photo', mode], '_')
	vmag_vname   = StrJoin([sc, instr, 'bulkvmag',      mode], '_')
	
	;FPI Moments
	vx_vname    = v_vname + '_x'
	vy_vname    = v_vname + '_y'
	vz_vname    = v_vname + '_z'
	vmag_vname  = v_vname + '_mag'
	pxx_vname   = p_vname + '_xx'
	pxy_vname   = p_vname + '_xy'
	pxz_vname   = p_vname + '_xz'
	pyy_vname   = p_vname + '_yy'
	pyz_vname   = p_vname + '_yz'
	pzz_vname   = p_vname + '_zz'
	p_scl_vname = p_vname + '_scl'
	qx_vname    = q_vname + '_x'
	qy_vname    = q_vname + '_y'
	qz_vname    = q_vname + '_z'
	qmag_vname  = q_vname + '_mag'
	
	;Calculated moments
	n_calc_vname     = StrJoin([sc, instr, 'numberdensity', 'calc', mode], '_')
	v_calc_vname     = StrJoin([sc, instr, 'bulkvx',        'calc', mode], '_')
	vx_calc_vname    = v_calc_vname + '_x'
	vy_calc_vname    = v_calc_vname + '_y'
	vz_calc_vname    = v_calc_vname + '_z'
	vmag_calc_vname  = v_calc_vname + '_mag'
	p_calc_vname     = StrJoin([sc, instr, 'pres',          'calc', mode], '_')
	pxx_calc_vname   = p_calc_vname + '_xx'
	pxy_calc_vname   = p_calc_vname + '_xy'
	pxz_calc_vname   = p_calc_vname + '_xz'
	pyy_calc_vname   = p_calc_vname + '_yy'
	pyz_calc_vname   = p_calc_vname + '_yz'
	pzz_calc_vname   = p_calc_vname + '_zz'
	p_scl_calc_vname = p_calc_vname + '_scl'
	q_calc_vname     = StrJoin([sc, instr, 'heatflux',      'calc', mode], '_')
	qx_calc_vname    = q_calc_vname + '_x'
	qy_calc_vname    = q_calc_vname + '_y'
	qz_calc_vname    = q_calc_vname + '_z'
	qmag_calc_vname  = q_calc_vname + '_mag'

;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, mode, $
		                     INSTR  = fgm_instr, $
		                     LEVEL  = level, $
		                     VARFORMAT = '*b_gse*', $
		                     SUFFIX = suffix
		
		;FPI
		MrMMS_FPI_Load_Dist3D, sc, mode, species, $
		                       /APPLY_MODEL, $
		                       COORD_SYS   = coord_sys, $
		                       LEVEL       = level, $
		                       ORIENTATION = orientation
		
		;Load FPI Moments
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = instr + '-moms', $
		                     TEAM_SITE = team_site, $
		                     VARFORMAT = [ '*numberdensity*',   '*bulkv_'+coords+'*', $
		                                   '*pres*'+coords+'*', '*heatq_'+coords+'*' ]
		
		;Spacecraft potential
		MrMMS_Load_Data, sc, 'edp', 'fast', 'l2', $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
	ENDIF

;-------------------------------------------
; Compute Moments //////////////////////////
;-------------------------------------------
	;Distribution function
	theSpecies = species EQ 'i' ? 'H' : species
	ofcorr = MrVar_Get(f_vname)
	oDist  = MrDist4D(ofcorr, VSC=scpot_vname, SPECIES=theSpecies)
	oDist -> Moments, /CACHE, $
	                  DENSITY     = oDensity, $
;	                  HEATFLUX    = oHeatFlux, $
	                  PRESSURE    = oPres, $
;	                  TEMPERATURE = oTemp, $
	                  VELOCITY    = oVelocity
	
	;Set names
	oDensity  -> SetName, n_calc_vname
;	oHeatFlux -> SetName, q_calc_vname
	oPres     -> SetName, p_calc_vname
;	oTemp     -> SetName, t_calc_vname
	oVelocity -> SetName, v_calc_vname
	
	oHeatFlux = oDist -> HeatFlux( NAME=q_calc_vname, /CACHE )

;-------------------------------------------
; Extract Components ///////////////////////
;-------------------------------------------
	;FPI Velocity
	oV    = MrVar_Get(v_vname)
	oV   -> Split, oVx, oVy, oVz, /CACHE
	oVmag = oV -> Magnitude( NAME=vmag_vname, /CACHE )
	
	;FPI Pressure
	oP   = MrVar_Get(p_vname)
	ot   = oP['TIMEVAR']
	oPxx = MrScalarTS( ot, oP[*,0,0], /CACHE, NAME=pxx_vname )
	oPxy = MrScalarTS( ot, oP[*,1,0], /CACHE, NAME=pxy_vname )
	oPxz = MrScalarTS( ot, oP[*,2,0], /CACHE, NAME=pxz_vname )
	oPyy = MrScalarTS( ot, oP[*,1,1], /CACHE, NAME=pyy_vname )
	oPyz = MrScalarTS( ot, oP[*,2,1], /CACHE, NAME=pyz_vname )
	oPzz = MrScalarTS( ot, oP[*,2,2], /CACHE, NAME=pzz_vname )
	op   = MrScalarTS( ot, oP[*,0,0] + oP[*,1,1] + oP[*,2,2], /CACHE, NAME=p_scl_vname )
	
	;FPI Heat Flux
	oQ    = MrVar_Get(q_vname)
	oQ   -> Split, oQx, oQy, oQz, /CACHE
	oQmag = oQ -> Magnitude( NAME=qmag_vname, /CACHE )
	
	
	;Calculated Velocity
	oV    = MrVar_Get(v_calc_vname)
	oV   -> Split, oVx, oVy, oVz, /CACHE
	oVmag = oV -> Magnitude( NAME=vmag_calc_vname, /CACHE )
	
	;Calculated Pressure
	oP   = MrVar_Get(p_calc_vname)
	ot   = oP['TIMEVAR']
	oPxx = MrScalarTS( ot, oP[*,0,0], /CACHE, NAME=pxx_calc_vname )
	oPxy = MrScalarTS( ot, oP[*,1,0], /CACHE, NAME=pxy_calc_vname )
	oPxz = MrScalarTS( ot, oP[*,2,0], /CACHE, NAME=pxz_calc_vname )
	oPyy = MrScalarTS( ot, oP[*,1,1], /CACHE, NAME=pyy_calc_vname )
	oPyz = MrScalarTS( ot, oP[*,2,1], /CACHE, NAME=pyz_calc_vname )
	oPzz = MrScalarTS( ot, oP[*,2,2], /CACHE, NAME=pzz_calc_vname )
	op   = MrScalarTS( ot, oP[*,0,0] + oP[*,1,1] + oP[*,2,2], /CACHE, NAME=p_scl_calc_vname )
	
	;Calculated Heat Flux
	oQ    = MrVar_Get(q_calc_vname)
	oQ   -> Split, oQx, oQy, oQz, /CACHE
	oQmag = oQ -> Magnitude( NAME=qmag_calc_vname, /CACHE )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	;B
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, coords], ' ' ) )
	
	;
	; FPI Moments
	;
	
	;Density
	oN = MrVar_Get(n_vname)
	oN['LABEL'] = 'moms'
	
	;Velocity
	oVx = MrVar_Get(vx_vname)
	oVx['LABEL'] = 'moms'
	oVx['TITLE'] = 'Vx!C(km/s)'
	
	oVy = MrVar_Get(vy_vname)
	oVy['LABEL'] = 'moms'
	oVy['TITLE'] = 'Vy!C(km/s)'
	
	oVz = MrVar_Get(vz_vname)
	oVz['LABEL'] = 'moms'
	oVz['TITLE'] = 'Vz!C(km/s)'
	
	oVmag = MrVar_Get(vmag_vname)
	oVmag['LABEL'] = 'moms'
	oVmag['TITLE'] = '|V|!C(km/s)'
	
	;Pressure
	oPxx = MrVar_Get(pxx_vname)
	oPxx['LABEL'] = 'moms'
	oPxx['TITLE'] = 'Pxx!C(nPa)'
	
	oPxy = MrVar_Get(pxy_vname)
	oPxy['LABEL'] = 'moms'
	oPxy['TITLE'] = 'Pxy!C(nPa)'
	
	oPxz = MrVar_Get(pxz_vname)
	oPxz['LABEL'] = 'moms'
	oPxz['TITLE'] = 'Pxz!C(nPa)'
	
	oPyy = MrVar_Get(pyy_vname)
	oPyy['LABEL'] = 'moms'
	oPyy['TITLE'] = 'Pyy!C(nPa)'
	
	oPyz = MrVar_Get(pyz_vname)
	oPyz['LABEL'] = 'moms'
	oPyz['TITLE'] = 'Pyz!C(nPa)'
	
	oPzz = MrVar_Get(pzz_vname)
	oPzz['LABEL'] = 'moms'
	oPzz['TITLE'] = 'Pzz!C(nPa)'

	op = MrVar_Get(p_scl_vname)
	op['LABEL'] = 'moms'
	op['TITLE'] = 'p!C(nPa)'
	
	;Heat Flux
	oQx = MrVar_Get(qx_vname)
	oQx['LABEL'] = 'moms'
	oQx['TITLE'] = 'Qx!C(mW/m^2)'
	
	oQy = MrVar_Get(qy_vname)
	oQy['LABEL'] = 'moms'
	oQy['TITLE'] = 'Qy!C(mW/m^2)'
	
	oQz = MrVar_Get(qz_vname)
	oQz['LABEL'] = 'moms'
	oQz['TITLE'] = 'Qz!C(mW/m^2)'
	
	oQmag = MrVar_Get(qmag_vname)
	oQmag['AXIS_RANGE'] = [0, oQmag.max]
	oQmag['LABEL']      = 'moms'
	oQmag['TITLE']      = '|Q|!C(mW/m^2)'
	
	;
	; Calculated Moments
	;
	
	;Density
	oN = MrVar_Get(n_calc_vname)
	oN['COLOR'] = 'Blue'
	oN['LABEL'] = 'dist'
	
	;Velocity
	oVx = MrVar_Get(vx_calc_vname)
	oVx['COLOR'] = 'Blue'
	oVx['LABEL'] = 'dist'
	
	oVy = MrVar_Get(vy_calc_vname)
	oVy['COLOR'] = 'Blue'
	oVy['LABEL'] = 'dist'
	
	oVz = MrVar_Get(vz_calc_vname)
	oVz['COLOR'] = 'Blue'
	oVz['LABEL'] = 'dist'
	
	oVmag = MrVar_Get(vmag_calc_vname)
	oVmag['COLOR'] = 'Blue'
	oVmag['LABEL'] = 'dist'
	
	;Pressure
	oPxx = MrVar_Get(pxx_calc_vname)
	oPxx['COLOR'] = 'Blue'
	oPxx['LABEL'] = 'dist'
	
	oPxy = MrVar_Get(pxy_calc_vname)
	oPxy['COLOR'] = 'Blue'
	oPxy['LABEL'] = 'dist'
	
	oPxz = MrVar_Get(pxz_calc_vname)
	oPxz['COLOR'] = 'Blue'
	oPxz['LABEL'] = 'dist'
	
	oPyy = MrVar_Get(pyy_calc_vname)
	oPyy['COLOR'] = 'Blue'
	oPyy['LABEL'] = 'dist'
	
	oPyz = MrVar_Get(pyz_calc_vname)
	oPyz['COLOR'] = 'Blue'
	oPyz['LABEL'] = 'dist'
	
	oPzz = MrVar_Get(pzz_calc_vname)
	oPzz['COLOR'] = 'Blue'
	oPzz['LABEL'] = 'dist'
	
	op = MrVar_Get(p_scl_calc_vname)
	op['COLOR'] = 'Blue'
	op['LABEL'] = 'dist'
	
	;Heat Flux
	oQx = MrVar_Get(qx_calc_vname)
	oQx['COLOR'] = 'Blue'
	oQx['LABEL'] = 'dist'
	oQx['TITLE'] = 'Qx!C(mW/m^2)'
	
	oQy = MrVar_Get(qy_calc_vname)
	oQy['COLOR'] = 'Blue'
	oQy['LABEL'] = 'dist'
	oQy['TITLE'] = 'Qy!C(mW/m^2)'
	
	oQz = MrVar_Get(qz_calc_vname)
	oQz['COLOR'] = 'Blue'
	oQz['LABEL'] = 'dist'
	oQz['TITLE'] = 'Qz!C(mW/m^2)'
	
	oQmag = MrVar_Get(qmag_calc_vname)
	oQmag['COLOR'] = 'Blue'
	oQmag['LABEL'] = 'dist'
	oQmag['TITLE'] = '|Q|!C(mW/m^2)'
	

;-------------------------------------------
; Window1 //////////////////////////////////
;-------------------------------------------
	;Window1: N, V, Q
	win1 = MrVar_PlotTS( [ b_vname, n_vname, vx_vname, vy_vname, vz_vname, vmag_vname, $
	                       qx_vname, qy_vname, qz_vname, qmag_vname ], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )

	;Overplot: Calculated N, V, Q
	win1 = MrVar_OPlotTS( [ n_vname, vx_vname, vy_vname, vz_vname, vmag_vname ], $
	                      [ n_calc_vname, vx_calc_vname, vy_calc_vname, vz_calc_vname, vmag_calc_vname ] )
	win1 = MrVar_OPlotTS( [ qx_vname, qy_vname, qz_vname, qmag_vname ], $
	                      [ qx_calc_vname, qy_calc_vname, qz_calc_vname, qmag_calc_vname ] )

	win1[0] -> SetLayout, [1,1]
	win1 -> TrimLayout
	win1 -> SetProperty, OXMARGIN=[10,6]

;-------------------------------------------
; Window2 //////////////////////////////////
;-------------------------------------------
	;Window2: P
	win2 = MrVar_PlotTS( [ b_vname, p_scl_vname, pxx_vname, pyy_vname, pzz_vname, pxy_vname, pxz_vname, pyz_vname ], $
	                     /NO_REFRESH, $
	                     YSIZE = 750 )

	;Overplot: Calculated P
	win2 = MrVar_OPlotTS( [ p_scl_vname, pxx_vname, pyy_vname, pzz_vname, pxy_vname, pxz_vname, pyz_vname ], $
	                      [ p_scl_calc_vname, pxx_calc_vname, pyy_calc_vname, pzz_calc_vname, pxy_calc_vname, pxz_calc_vname, pyz_calc_vname ] )

	win2[0] -> SetLayout, [1,1]
	win2 -> TrimLayout
	win2 -> SetProperty, OXMARGIN=[12,6]

	win1 -> Refresh
	win2 -> Refresh
	RETURN, win2
END