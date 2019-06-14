; docformat = 'rst'
;
; NAME:
;    MrMMS_FPI_Test_Moments
;
; PURPOSE:
;+
;   Compare UNH-made particle moments with those of FPI.
;
; :Categories:
;    MMS, SPEDAS
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID: {'mms1' | 'mms2' | 'mms3' | 'mm4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;       SPECIES:    in, required, type=string
;                   Particle species: {'i' | 'e'}
;
;
; :Keywords:
;       GROUND:     in, optional, type=boolean, default=0
;                   If set, moments are calculated by shifting the energy targets by
;                       the spacecraft potential (as is typically done in ground
;                       processing). Otherwise, phase space is shifted via Liouville
;                       mapping.
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;
; :Author:
;    Matthew Argall::
;    University of New Hampshire
;    Morse Hall Room 348
;    8 College Road
;    Durham, NH 03824
;    matthew.argall@unh.edu
;
; :History:
;    Modification History::
;       2016/02/10  -   Written by Matthew Argall
;       2018/09/29  -   Photoelectron models now applied within FPI_LOAD_DATA, moments
;                           calculated with updated methods, plot results in two columns. - MRA
;       2018/09/30  -   Save plot to file. - MRA
;       2018/10/24  -   Compare to equivalent Maxwellian moments. - MRA
;-
FUNCTION MrMMS_FPI_Test_Moments, sc, mode, species, $
GROUND=ground, $
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
	
	instr         = 'd' + species + 's'
	level         = 'l2'
	tf_ground     = Keyword_Set(ground)
	tf_load       = ~Keyword_Set(no_load)
	tf_maxwellian = Keyword_Set(maxwellian)

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;Source names
	scpot_vname = StrJoin( [sc, 'edp', 'scpot', 'fast', level], '_' )
	f_vname     = StrJoin( [sc, instr, 'dist',              mode], '_')
	n_vname     = StrJoin( [sc, instr, 'numberdensity',      mode], '_')
	v_vname     = StrJoin( [sc, instr, 'bulkv',      'dbcs', mode], '_')
	vmag_vname  = StrJoin( [sc, instr, 'bulkv',      'mag',  mode], '_')
	p_vname     = StrJoin( [sc, instr, 'prestensor', 'dbcs', mode], '_')
	p_scl_vname = StrJoin( [sc, instr, 'pressure',           mode], '_')
	t_vname     = StrJoin( [sc, instr, 'temptensor', 'dbcs', mode], '_')
	q_vname     = StrJoin( [sc, instr, 'heatq',      'dbcs', mode], '_')
	
	;Derived names
	f_max_vname = StrJoin([sc, instr, 'f', 'maxwellian', mode], '_')
	ncorr_vname = StrJoin([sc, instr, 'numberdensity', 'corr', mode], '_')
	vcorr_vname = StrJoin([sc, instr, 'bulkv',         'corr', mode], '_')
	vc_mag_vname = StrJoin([sc, instr, 'bulkv', 'mag', 'corr', mode], '_')
	pcorr_vname = StrJoin([sc, instr, 'pres',          'corr', mode], '_')
	pc_scl_vname = StrJoin([sc, instr, 'presscl',   'corr', mode], '_')
	qcorr_vname = StrJoin([sc, instr, 'heatflux',      'corr', mode], '_')
	dist_vname  = StrJoin([sc, instr, 'dist4d',        mode], '_')
	
	vx_vname    = v_vname     + '_x'
	vy_vname    = v_vname     + '_y'
	vz_vname    = v_vname     + '_z'
	pxx_vname   = p_vname     + '_xx'
	pxy_vname   = p_vname     + '_xy'
	pxz_vname   = p_vname     + '_xz'
	pyy_vname   = p_vname     + '_yy'
	pyz_vname   = p_vname     + '_yz'
	pzz_vname   = p_vname     + '_zz'
	vc_x_vname  = vcorr_vname + '_x'
	vc_y_vname  = vcorr_vname + '_y'
	vc_z_vname  = vcorr_vname + '_z'
	pc_xx_vname = pcorr_vname + '_xx'
	pc_xy_vname = pcorr_vname + '_xy'
	pc_xz_vname = pcorr_vname + '_xz'
	pc_yy_vname = pcorr_vname + '_yy'
	pc_yz_vname = pcorr_vname + '_yz'
	pc_zz_vname = pcorr_vname + '_zz'


;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;Load FPI Distribution
		MrMMS_FPI_Load_Dist3D, sc, mode, species, $;, fac, $
		                       /APPLY_MODEL, $ 
		                       COORD_SYS   = coord_sys, $
		                       LEVEL       = level, $
		                       ORIENTATION = orientation, $
		                       SUFFIX      = suffix, $
		                       TEAM_SITE   = team_site, $
		                       TRANGE      = trange, $
		                       VARNAMES    = varnames
		
		;Load FPI Moments
		MrMMS_FPI_Load_Data, sc, mode, $
		                     OPTDESC   = instr + '-moms', $
		                     TEAM_SITE = team_site, $
		                     VARFORMAT = ['*numberdensity*', '*bulkv_dbcs*', $
		                                  '*pres*dbcs*', '*temp*dbcs*', '*heat*']
		
		;Spacecraft potential
		MrMMS_Load_Data, sc, 'edp', 'fast', 'l2', $
		                 OPTDESC  = 'scpot', $
		                 VARNAMES = '*scpot*'
	ENDIF

;-------------------------------------------
; Maxwellian Distribution //////////////////
;-------------------------------------------
	IF tf_maxwellian THEN BEGIN
		;Scalar temperature
		oT = MrVar_Get(t_vname)
		oT = (oT[*,0,0] + oT[*,1,1] + oT[*,2,2]) / 3.0
		
		;Maxwellian distribution
		ofb = MrMMS_FPI_f_Maxwellian( f_vname, species, $
		                              /CACHE, $
		                              DENSITY     = n_vname, $
		                              NAME        = f_max_vname, $
		                              TEMPERATURE = oT, $
		                              VELOCITY    = v_vname )
		
		;Integrate the Maxwellian, not the FPI distribution
		f_vname = f_max_vname
	ENDIF

;-------------------------------------------
; Compute Moments //////////////////////////
;-------------------------------------------
	;Distribution function
	oDist = MrDist4D(f_vname, VSC=scpot_vname, SPECIES=species)
	oN    = oDist -> Density(/CACHE,  GROUND=ground, NAME=ncorr_vname)
	oV    = oDist -> Velocity(/CACHE, GROUND=ground, NAME=vcorr_vname)
	oP    = oDist -> Pressure(/CACHE, GROUND=ground, NAME=pcorr_vname)
	
;-------------------------------------------
; FPI: Split Into Components ///////////////
;-------------------------------------------
	;Velocity
	oV_fpi = MrVar_Get(v_vname)
	oV_fpi -> Split, oVx, oVy, oVz, /CACHE
	
	;Magnitude
	oVmag = oV_fpi -> Magnitude(/CACHE, NAME=vmag_vname)
	
	;Pressure
	oP   = MrVar_Get(p_vname)
	oPxx = MrScalarTS(oP['TIMEVAR'], oP[*,0,0], NAME=oP.name+'_xx', /CACHE)
	oPxy = MrScalarTS(oP['TIMEVAR'], oP[*,0,1], NAME=oP.name+'_xy', /CACHE)
	oPxz = MrScalarTS(oP['TIMEVAR'], oP[*,0,2], NAME=oP.name+'_xz', /CACHE)
	oPyy = MrScalarTS(oP['TIMEVAR'], oP[*,1,1], NAME=oP.name+'_yy', /CACHE)
	oPyz = MrScalarTS(oP['TIMEVAR'], oP[*,1,2], NAME=oP.name+'_yz', /CACHE)
	oPzz = MrScalarTS(oP['TIMEVAR'], oP[*,2,2], NAME=oP.name+'_zz', /CACHE)
	
	;Scalar Pressure
	oPscl = (oPxx + oPyy + oPzz) / 3.0
	oPscl -> SetName, p_scl_vname
	oPscl -> Cache
	
;-------------------------------------------
; Calc: Split Into Components //////////////
;-------------------------------------------
	
	;Velocity
	oVcorr = MrVar_Get(vcorr_vname)
	oVcorr -> Split, oVc_x, oVc_y, oVc_z, /CACHE
	
	;Magnitude
	oVc_mag = oVcorr -> Magnitude(/CACHE, NAME=vc_mag_vname) 
	
	;Pressure
	oPcorr = MrVar_Get(pcorr_vname)
	oPc_xx = MrScalarTS(oPcorr['TIMEVAR'], oPcorr[*,0,0], NAME=oPcorr.name+'_xx', /CACHE)
	oPc_xy = MrScalarTS(oPcorr['TIMEVAR'], oPcorr[*,0,1], NAME=oPcorr.name+'_xy', /CACHE)
	oPc_xz = MrScalarTS(oPcorr['TIMEVAR'], oPcorr[*,0,2], NAME=oPcorr.name+'_xz', /CACHE)
	oPc_yy = MrScalarTS(oPcorr['TIMEVAR'], oPcorr[*,1,1], NAME=oPcorr.name+'_yy', /CACHE)
	oPc_yz = MrScalarTS(oPcorr['TIMEVAR'], oPcorr[*,1,2], NAME=oPcorr.name+'_yz', /CACHE)
	oPc_zz = MrScalarTS(oPcorr['TIMEVAR'], oPcorr[*,2,2], NAME=oPcorr.name+'_zz', /CACHE)

	;Scalar Pressure
	oPc_scl = (oPxx + oPyy + oPzz) / 3.0
	oPc_scl -> SetName, pc_scl_vname
	oPc_scl -> Cache
	
;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	clabel = tf_maxwellian ? 'max' : 'dist'
	
	;Density FPI
	oN_fpi = MrVar_Get(n_vname)
	oN_fpi['PLOT_TITLE'] = StrUpCase(StrJoin([sc, mode], ' ')) + ' Moments'
	oN_fpi['COLOR']  = 'Black'
;	oN_fpi['LABEL']  = 'moms'
	
	;Density Corrected
	oNcorr = MrVar_Get(ncorr_vname)
	oNcorr['COLOR']     = 'Blue'
;	oNcorr['LABEL']     = 'dist'
	
	;FPI Velocity
;	oVmag['LABEL'] = 'moms'
	oVmag['TITLE'] = '|V|!C(km/s)'
	
;	oVx['LABEL'] = 'moms'
	oVx['TITLE'] = 'Vx!C(km/s)'
	
;	oVy['LABEL'] = 'moms'
	oVy['TITLE'] = 'Vy!C(km/s)'
	
;	oVz['LABEL'] = 'moms'
	oVz['TITLE'] = 'Vz!C(km/s)'
	
	;Velocity
	oVc_mag['COLOR'] = 'Blue'
;	oVc_mag['LABEL'] = 'dist'
	oVc_mag['TITLE'] = '|V|!C(km/s)'
	
	oVc_x['COLOR'] = 'Blue'
;	oVc_x['LABEL'] = 'dist'
	oVc_x['TITLE'] = 'Vx!C(km/s)'
	
	oVc_y['COLOR'] = 'Blue'
;	oVc_y['LABEL'] = 'dist'
	oVc_z['TITLE'] = 'Vx!C(km/s)'
	
	oVc_z['COLOR'] = 'Blue'
;	oVc_z['LABEL'] = 'dist'
	oVc_z['TITLE'] = 'Vx!C(km/s)'
	
	;FPI Pressure
;	oPscl['LABEL'] = 'moms'
	oPscl['TITLE'] = 'p!C(nPa)'
	
	oPxx['LABEL'] = 'moms'
	oPxx['TITLE'] = 'Pxx!C(nPa)'
	
	oPxy['LABEL'] = 'moms'
	oPxy['TITLE'] = 'Pxy!C(nPa)'
	
	oPxz['LABEL'] = 'moms'
	oPxz['TITLE'] = 'Pxz!C(nPa)'
	
	oPyy['LABEL'] = 'moms'
	oPyy['TITLE'] = 'Pyy!C(nPa)'
	
	oPyz['LABEL'] = 'moms'
	oPyz['TITLE'] = 'Pyz!C(nPa)'
	
	oPzz['LABEL'] = 'moms'
	oPzz['TITLE'] = 'Pzz!C(nPa)'
	
	;Pressure
	oPc_scl['COLOR'] = 'Blue'
;	oPc_scl['LABEL'] = 'dist'
	oPc_scl['TITLE'] = 'p!C(nPa)'
	
	oPc_xx['COLOR'] = 'Blue'
	oPc_xx['LABEL'] = clabel
	oPc_xx['TITLE'] = 'Pxx!C(nPa)'
	
	oPc_xy['COLOR'] = 'Blue'
	oPc_xy['LABEL'] = clabel
	oPc_xy['TITLE'] = 'Pxy!C(nPa)'
	
	oPc_xz['COLOR'] = 'Blue'
	oPc_xz['LABEL'] = clabel
	oPc_xz['TITLE'] = 'Pxz!C(nPa)'
	
	oPc_yy['COLOR'] = 'Blue'
	oPc_yy['LABEL'] = clabel
	oPc_yy['TITLE'] = 'Pyy!C(nPa)'
	
	oPc_yz['COLOR'] = 'Blue'
	oPc_yz['LABEL'] = clabel
	oPc_yz['TITLE'] = 'Pyz!C(nPa)'
	
	oPc_zz['COLOR'] = 'Blue'
	oPc_zz['LABEL'] = clabel
	oPc_zz['TITLE'] = 'Pzz!C(nPa)'

;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	;Create window
	win = MrWindow( LAYOUT  = [2,6], $
	                REFRESH = 0, $
	                XSIZE   = 800, $
	                YGAP    = 0.5, $
	                YSIZE   = 650 )
	
	;Compare moments
	win = MrVar_PlotTS( [ n_vname,     pxx_vname, $
	                      vmag_vname,  pxy_vname, $
	                      vx_vname,    pxz_vname, $
	                      vy_vname,    pyy_vname, $
	                      vz_vname,    pyz_vname, $
	                      p_scl_vname, pzz_vname ], $
	                    /CURRENT )
	win = MrVar_OPlotTS( [n_vname, vmag_vname, vx_vname, vy_vname, vz_vname, p_scl_vname, $
	                      pxx_vname, pxy_vname, pxz_vname, pyy_vname, pyz_vname, pzz_vname], $
	                     [ncorr_vname, vc_mag_vname, vc_x_vname, vc_y_vname, vc_z_vname, pc_scl_vname, $
	                      pc_xx_vname, pc_xy_vname, pc_xz_vname, pc_yy_vname, pc_yz_vname, pc_zz_vname] )
	
	;Fix layout
	win[0] -> SetLayout, [1,1]
	win -> TrimLayout
	win -> SetProperty, OXMARGIN=[12,8]
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
		optdesc = 'test-moms' + (tf_ground ? '-gnd' : '') + (tf_maxwellian ? '-maxwell' : '')
		fname   = StrJoin( [sc, instr, mode, level, optdesc], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END