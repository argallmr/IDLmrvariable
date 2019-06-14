; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FPI_Moms_Compare
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
;       1.  Bxyz, |B|             1. Bxyz, |B|       1. Bxyz, |B|
;       2.  Density               2. Pxx             2. Txx
;       3.  Vx                    3. Pxy             3. Txy
;       4.  Vy                    4. Pxz             4. Txz
;       5.  Vz                    5. Pyy             5. Tyy
;       6.  |V|                   6. Pyz             6. Tyz
;       7.  Qx                    7. Pzz             7. Tzz
;       8.  Qy                    8. p               8. T
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
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
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
;       2018/06/21  -   Include temperature in the comparison. - MRA
;-
FUNCTION MrMMS_Plot_FPI_Moms_Compare, sc, mode, species, $
COORDS=coords, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
MAXWELLIAN=maxwellian, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, Obj_New()
	ENDIF
	
	tf_load = ~Keyword_Set(no_load)
	tf_maxwellian = Keyword_Set(maxwellian)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'dbcs'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(species)   EQ 0 THEN species   = 'e'
	instr   = 'd' + species + 's'

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	fgm_coords = coords EQ 'dbcs' || coords EQ 'dsl' ? 'dmpa' : coords

	;Source B-Field
	b_vname    = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, mode, level], '_' )
	bvec_vname = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, mode, level], '_' )
	bmag_vname = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, mode, level], '_' )
	
	;Source Moments
	scpot_vname = StrJoin( [sc, 'edp', 'scpot', 'fast', level], '_' )
	n_vname     = StrJoin( [sc, instr, 'numberdensity',        mode], '_')
	v_vname     = StrJoin( [sc, instr, 'bulkv',      coords, mode], '_')
	p_vname     = StrJoin( [sc, instr, 'prestensor', coords, mode], '_')
	t_vname     = StrJoin( [sc, instr, 'temptensor', coords, mode], '_')
	q_vname     = StrJoin( [sc, instr, 'heatq',      coords, mode], '_')
	
	;Source Distribution
	f_vname      = StrJoin( [sc, instr, 'dist',                  mode], '_')
	
	;Derived names
	dist_vname   = StrJoin([sc, instr, 'dist4d',        mode], '_')
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
	txx_vname   = t_vname + '_xx'
	txy_vname   = t_vname + '_xy'
	txz_vname   = t_vname + '_xz'
	tyy_vname   = t_vname + '_yy'
	tyz_vname   = t_vname + '_yz'
	tzz_vname   = t_vname + '_zz'
	t_scl_vname = t_vname + '_scl'
	qx_vname    = q_vname + '_x'
	qy_vname    = q_vname + '_y'
	qz_vname    = q_vname + '_z'
	qmag_vname  = q_vname + '_mag'
	
	;Calculated moments
	optdesc          = tf_maxwellian ? 'max' : 'calc'
	n_calc_vname     = StrJoin([sc, instr, 'numberdensity', optdesc, mode], '_')
	v_calc_vname     = StrJoin([sc, instr, 'bulkvx',        optdesc, mode], '_')
	vx_calc_vname    = v_calc_vname + '_x'
	vy_calc_vname    = v_calc_vname + '_y'
	vz_calc_vname    = v_calc_vname + '_z'
	vmag_calc_vname  = v_calc_vname + '_mag'
	p_calc_vname     = StrJoin([sc, instr, 'pres', optdesc, mode], '_')
	pxx_calc_vname   = p_calc_vname + '_xx'
	pxy_calc_vname   = p_calc_vname + '_xy'
	pxz_calc_vname   = p_calc_vname + '_xz'
	pyy_calc_vname   = p_calc_vname + '_yy'
	pyz_calc_vname   = p_calc_vname + '_yz'
	pzz_calc_vname   = p_calc_vname + '_zz'
	p_scl_calc_vname = p_calc_vname + '_scl'
	t_calc_vname     = StrJoin([sc, instr, 'temp', optdesc, mode], '_')
	txx_calc_vname   = t_calc_vname + '_xx'
	txy_calc_vname   = t_calc_vname + '_xy'
	txz_calc_vname   = t_calc_vname + '_xz'
	tyy_calc_vname   = t_calc_vname + '_yy'
	tyz_calc_vname   = t_calc_vname + '_yz'
	tzz_calc_vname   = t_calc_vname + '_zz'
	t_scl_calc_vname = t_calc_vname + '_scl'
	q_calc_vname     = StrJoin([sc, instr, 'heatflux', optdesc, mode], '_')
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
		                     VARFORMAT = '*b_'+fgm_coords+'*', $
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
		                                   '*pres*'+coords+'*', '*temptensor*'+coords+'*', $
		                                   '*heatq_'+coords+'*' ]
		
		;Spacecraft potential
		MrMMS_Load_Data, sc, 'edp', 'fast', 'l2', $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
	ENDIF

;-------------------------------------------
; FPI: Extract Components //////////////////
;-------------------------------------------
	;FPI Velocity
	oV    = MrVar_Get(v_vname)
	oV   -> Split, oVx, oVy, oVz, /CACHE
	oVmag = oV -> Magnitude( NAME=vmag_vname, /CACHE )
	
	;FPI Pressure
	oP    = MrVar_Get(p_vname)
	oTime = oP['TIMEVAR']
	oPxx  = MrScalarTS( oTime, oP[*,0,0], /CACHE, NAME=pxx_vname )
	oPxy  = MrScalarTS( oTime, oP[*,1,0], /CACHE, NAME=pxy_vname )
	oPxz  = MrScalarTS( oTime, oP[*,2,0], /CACHE, NAME=pxz_vname )
	oPyy  = MrScalarTS( oTime, oP[*,1,1], /CACHE, NAME=pyy_vname )
	oPyz  = MrScalarTS( oTime, oP[*,2,1], /CACHE, NAME=pyz_vname )
	oPzz  = MrScalarTS( oTime, oP[*,2,2], /CACHE, NAME=pzz_vname )
	op    = MrScalarTS( oTime, (oP[*,0,0] + oP[*,1,1] + oP[*,2,2]) / 3.0, /CACHE, NAME=p_scl_vname )
	
	;FPI Temperature
	oT = MrVar_Get(t_vname)
	oTxx = MrScalarTS( oTime, oT[*,0,0], /CACHE, NAME=txx_vname )
	oTxy = MrScalarTS( oTime, oT[*,1,0], /CACHE, NAME=txy_vname )
	oTxz = MrScalarTS( oTime, oT[*,2,0], /CACHE, NAME=txz_vname )
	oTyy = MrScalarTS( oTime, oT[*,1,1], /CACHE, NAME=tyy_vname )
	oTyz = MrScalarTS( oTime, oT[*,2,1], /CACHE, NAME=tyz_vname )
	oTzz = MrScalarTS( oTime, oT[*,2,2], /CACHE, NAME=tzz_vname )
	oT   = MrScalarTS( oTime, (oT[*,0,0] + oT[*,1,1] + oT[*,2,2]) / 3.0, /CACHE, NAME=t_scl_vname )
	
	;FPI Heat Flux
	oQ    = MrVar_Get(q_vname)
	oQ   -> Split, oQx, oQy, oQz, /CACHE
	oQmag = oQ -> Magnitude( NAME=qmag_vname, /CACHE )
	
;-------------------------------------------
; Compute Moments //////////////////////////
;-------------------------------------------
	;Distribution function
	theSpecies = species EQ 'i' ? 'H' : species
	IF tf_maxwellian $
		THEN oDist = MrMMS_FPI_F_Maxwellian(f_vname, theSpecies, DENSITY=n_vname, VELOCITY=v_vname, TEMPERATURE=t_scl_vname) $
		ELSE oDist = MrVar_Get(f_vname)
	
	;Calculate moments
	oDist4D    = MrDist4D(oDist, VSC=(tf_maxwellian ? !Null : scpot_vname), SPECIES=theSpecies)
	oDist4D    = MrDist4D(oDist, SPECIES=theSpecies)
	oDensity   = oDist4D -> Density(/CACHE, NAME=n_calc_vname)
	oPres      = oDist4D -> Pressure(/CACHE, NAME=p_calc_vname)
	oTemp      = oDist4D -> Temperature(/CACHE, NAME=t_calc_vname)
	oVelocity  = oDist4D -> Velocity(/CACHE, NAME=v_calc_vname)
	oHeatFlux  = oDist4D -> HeatFlux(/CACHE, NAME=q_calc_vname)

;-------------------------------------------
; Calc: Extract Components /////////////////
;-------------------------------------------
	
	;Calculated Velocity
	oV    = MrVar_Get(v_calc_vname)
	oV   -> Split, oVx, oVy, oVz, /CACHE
	oVmag = oV -> Magnitude( NAME=vmag_calc_vname, /CACHE )
	
	;Calculated Pressure
	oP    = MrVar_Get(p_calc_vname)
	oTime = oP['TIMEVAR']
	oPxx  = MrScalarTS( oTime, oP[*,0,0], /CACHE, NAME=pxx_calc_vname )
	oPxy  = MrScalarTS( oTime, oP[*,1,0], /CACHE, NAME=pxy_calc_vname )
	oPxz  = MrScalarTS( oTime, oP[*,2,0], /CACHE, NAME=pxz_calc_vname )
	oPyy  = MrScalarTS( oTime, oP[*,1,1], /CACHE, NAME=pyy_calc_vname )
	oPyz  = MrScalarTS( oTime, oP[*,2,1], /CACHE, NAME=pyz_calc_vname )
	oPzz  = MrScalarTS( oTime, oP[*,2,2], /CACHE, NAME=pzz_calc_vname )
	op    = MrScalarTS( oTime, (oP[*,0,0] + oP[*,1,1] + oP[*,2,2]) / 3.0, /CACHE, NAME=p_scl_calc_vname )
	
	;Calculated Temperature
	oT   = MrVar_Get(t_calc_vname)
	oTxx = MrScalarTS( oTime, oT[*,0,0], /CACHE, NAME=txx_calc_vname )
	oTxy = MrScalarTS( oTime, oT[*,1,0], /CACHE, NAME=txy_calc_vname )
	oTxz = MrScalarTS( oTime, oT[*,2,0], /CACHE, NAME=txz_calc_vname )
	oTyy = MrScalarTS( oTime, oT[*,1,1], /CACHE, NAME=tyy_calc_vname )
	oTyz = MrScalarTS( oTime, oT[*,2,1], /CACHE, NAME=tyz_calc_vname )
	oTzz = MrScalarTS( oTime, oT[*,2,2], /CACHE, NAME=tzz_calc_vname )
	oT   = MrScalarTS( oTime, (oT[*,0,0] + oT[*,1,1] + oT[*,2,2]) / 3.0, /CACHE, NAME=t_scl_calc_vname )
	
	;Calculated Heat Flux
	oQ    = MrVar_Get(q_calc_vname)
	oQ   -> Split, oQx, oQy, oQz, /CACHE
	oQmag = oQ -> Magnitude( NAME=qmag_calc_vname, /CACHE )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	;B
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, instr, mode, level, coords], ' ' ) )
	
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
	
	;Temperature
	oTxx = MrVar_Get(txx_vname)
	oTxx['LABEL'] = 'moms'
	oTxx['TITLE'] = 'Txx!C(eV)'
	
	oTxy = MrVar_Get(txy_vname)
	oTxy['LABEL'] = 'moms'
	oTxy['TITLE'] = 'Txy!C(eV)'
	
	oTxz = MrVar_Get(txz_vname)
	oTxz['LABEL'] = 'moms'
	oTxz['TITLE'] = 'Txz!C(eV)'
	
	oTyy = MrVar_Get(tyy_vname)
	oTyy['LABEL'] = 'moms'
	oTyy['TITLE'] = 'Tyy!C(eV)'
	
	oTyz = MrVar_Get(tyz_vname)
	oTyz['LABEL'] = 'moms'
	oTyz['TITLE'] = 'Tyz!C(eV)'
	
	oTzz = MrVar_Get(tzz_vname)
	oTzz['LABEL'] = 'moms'
	oTzz['TITLE'] = 'Tzz!C(eV)'

	oT = MrVar_Get(t_scl_vname)
	oT['LABEL'] = 'moms'
	oT['TITLE'] = 'T!C(eV)'
	
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
	
	;Temperature
	oTxx = MrVar_Get(txx_calc_vname)
	oTxx['COLOR'] = 'Blue'
	oTxx['LABEL'] = 'dist'
	
	oTxy = MrVar_Get(txy_calc_vname)
	oTxy['COLOR'] = 'Blue'
	oTxy['LABEL'] = 'dist'
	
	oTxz = MrVar_Get(txz_calc_vname)
	oTxz['COLOR'] = 'Blue'
	oTxz['LABEL'] = 'dist'
	
	oTyy = MrVar_Get(tyy_calc_vname)
	oTyy['COLOR'] = 'Blue'
	oTyy['LABEL'] = 'dist'
	
	oTyz = MrVar_Get(tyz_calc_vname)
	oTyz['COLOR'] = 'Blue'
	oTyz['LABEL'] = 'dist'
	
	oTzz = MrVar_Get(tzz_calc_vname)
	oTzz['COLOR'] = 'Blue'
	oTzz['LABEL'] = 'dist'
	
	oT = MrVar_Get(t_scl_calc_vname)
	oT['COLOR'] = 'Blue'
	oT['LABEL'] = 'dist'
	
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
	win1 -> SetProperty, OXMARGIN=[14,8]
	win1.name = 'FPI MOMS NVQ'

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
	win2 -> SetProperty, OXMARGIN=[12,8]
	win2.name = 'FPI MOMS P'

;-------------------------------------------
; Window3 //////////////////////////////////
;-------------------------------------------
	;Window2: P
	win3 = MrVar_PlotTS( [ b_vname, t_scl_vname, txx_vname, tyy_vname, tzz_vname, txy_vname, txz_vname, tyz_vname ], $
	                     /NO_REFRESH, $
	                     YSIZE = 750 )

	;Overplot: Calculated P
	win3 = MrVar_OPlotTS( [ t_scl_vname, txx_vname, tyy_vname, tzz_vname, txy_vname, txz_vname, tyz_vname ], $
	                      [ t_scl_calc_vname, txx_calc_vname, tyy_calc_vname, tzz_calc_vname, txy_calc_vname, txz_calc_vname, tyz_calc_vname ] )

	win3[0] -> SetLayout, [1,1]
	win3 -> TrimLayout
	win3 -> SetProperty, OXMARGIN=[12,8]
	win3.name = 'FPI MOMS T'

;-------------------------------------------
; Save the File ////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
			MrPrintF, 'LogText', 'Saving file to: "' + output_dir + '".'
		ENDIF
		
		;File name
		fname1 = StrJoin([sc, instr, mode, level, 'moms-comp-nvq'], '_')
		fname2 = StrJoin([sc, instr, mode, level, 'moms-comp-p'], '_')
		fname3 = StrJoin([sc, instr, mode, level, 'moms-comp-t'], '_')
		fname1 = FilePath(fname1, ROOT_DIR=output_dir)
		fname2 = FilePath(fname2, ROOT_DIR=output_dir)
		fname3 = FilePath(fname3, ROOT_DIR=output_dir)
		
		;Save
		fout = MrVar_PlotTS_Save( win1, fname1, output_ext )
		fout = MrVar_PlotTS_Save( win2, fname2, output_ext )
		fout = MrVar_PlotTS_Save( win3, fname3, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------

	win1 -> Refresh
	win2 -> Refresh
	win3 -> Refresh
	RETURN, win2
END