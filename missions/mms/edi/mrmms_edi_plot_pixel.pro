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
FUNCTION MrMMS_EDI_Plot_Pixel, sc, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
TRANGE=trange
	compile_opt idl2

	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF N_Elements(win) GT 0 THEN Obj_Destroy, win
		MrPrintF, 'LogErr'
		RETURN, Obj_New()
	ENDIF
	
	sc     = 'mms1'
	mode   = 'brst'
	level  = 'l2'
	coords = 'dbcs'
	tf_load = ~Keyword_Set(no_load)
	mrvar_settrange, '2017-07-06T' + ['13:53:35', '13:54:40']
	
	me   = MrConstants('m_e', /DOUBLE)
	eV2J = MrConstants('eV2J', /DOUBLE)
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	
	chan = mode EQ 'brst' ? ['1', '2', '3', '4'] : '1'
	
	;Variable names
	f_vname             = StrJoin( [sc, 'des', 'dist', mode], '_' )
	edi_flux_0_vnames   = StrJoin( [sc, 'edi', 'flux',   '0', mode, level], '_' )
	edi_flux_180_vnames = StrJoin( [sc, 'edi', 'flux', '180', mode, level], '_' )
	edi_traj_0_vnames   = StrJoin( [sc, 'edi', 'traj', coords,   '0', mode, level], '_' )
	edi_traj_180_vnames = StrJoin( [sc, 'edi', 'traj', coords, '180', mode, level], '_' )
	
	edi_flux_0_vnames   = StrMid(edi_flux_0_vnames, 0, 13)   + chan + StrMid(edi_flux_0_vnames, 13)
	edi_flux_180_vnames = StrMid(edi_flux_180_vnames, 0, 13) + chan + StrMid(edi_flux_180_vnames, 13)
	edi_traj_0_vnames   = StrMid(edi_traj_0_vnames, 0, 13)   + chan + StrMid(edi_traj_0_vnames, 13)
	edi_traj_180_vnames = StrMid(edi_traj_180_vnames, 0, 13) + chan + StrMid(edi_traj_180_vnames, 13)
	edi_flux_vnames = [edi_flux_0_vnames, edi_flux_180_vnames]
	edi_traj_vnames = [edi_traj_0_vnames, edi_traj_180_vnames]
	
	;Output
	flux_edi_vnames = [StrMid(edi_flux_0_vnames, 0, 16) + '_resample' + StrMid(edi_flux_0_vnames, 16), $
	                   StrMid(edi_flux_180_vnames, 0, 18) + '_resample' + StrMid(edi_flux_180_vnames, 18)]
	traj_edi_vnames = [StrMid(edi_traj_0_vnames, 0, 21) + '_resample' + StrMid(edi_traj_0_vnames, 21), $
	                   StrMid(edi_traj_180_vnames, 0, 23) + '_resample' + StrMid(edi_traj_180_vnames, 23)]
	flux_des_vnames = StrMid(flux_edi_vnames, 0, 5) + 'des' + StrMid(flux_edi_vnames, 8)
	traj_des_vnames = StrMid(traj_edi_vnames, 0, 5) + 'des' + StrMid(traj_edi_vnames, 8)
	traj_dot_vnames = [StrMid(edi_traj_0_vnames, 0, 21) + '_dot_des' + StrMid(edi_traj_0_vnames, 21), $
	                   StrMid(edi_traj_180_vnames, 0, 23) + '_dot_des' + StrMid(edi_traj_180_vnames, 23)]
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	
	IF tf_load THEN BEGIN
		;Load EDI data
		MrMMS_Load_Data, sc, 'edi', mode, level, $
		                 OPTDESC   = ['amb', 'amb-pm2'], $
		                 VARFORMAT = ['*flux?_0', '*flux?_180', $
		                              '*traj?_dbcs_0', '*traj?_dbcs_180'] + '_'+mode+'*'
		
		;Load DES-Dist data
		MrMMS_FPI_Load_Dist3D, sc, 'brst', 'e', $
		                       /APPLY_MODEL
	
	;-------------------------------------------
	; Center FPI Time Stamps ///////////////////
	;-------------------------------------------
	
		;Center FPI timestamps
		oDist   = MrVar_Get(f_vname)
		oT_des = oDist['DEPEND_0']
		dt     = (oT_des['DELTA_PLUS_VAR'] - oT_des['DELTA_MINUS_VAR'])/2.0
		oT     = MrTimeVar( oT_des['DATA', 'SSM'] + dt['DATA'], 'SSM', $
							NAME  = 'des_center_times', $
							T_REF = oT_des['DATA', 0] )
		oT_des -> CopyAttrTo, oT
		oT['DELTA_PLUS'] = dt['DATA']
		oT['DELTA_MINUS'] = dt['DATA']
		oT -> RemoveAttr, ['DELTA_PLUS_VAR', 'DELTA_MINUS_VAR']
		oDist['DEPEND_0'] = oT
	ENDIF
	
;-------------------------------------------
; Adjust DES Grid //////////////////////////
;-------------------------------------------
	oDist = MrVar_Get(f_vname)
	dims    = Size(oDist, /DIMENSIONS)
	nTime   = dims[0]
	nPhi    = dims[1]
	nTheta  = dims[2]
	nEnergy = dims[3]
	
	oT = oDist['DEPEND_0']
	
	;Convert PHI to trajectories and put in the range [-180, 180]
	oPhi = oDist['DEPEND_1']
	phi  = oPhi['DATA'] + 180.0
	iPhi = Where(phi GT 360.0)
	phi[iPhi] -= 360.0
	iPhi = Where(phi GT 180)
	phi[iPhi] = phi[iPhi] - 360
	
	;Convert THETA to trajectory
	oTheta = oDist['DEPEND_2']
	theta  = 180.0 - oTheta['DATA']
	
	;Convert ENERGY to velocity (cm/s)
	oEnergy  = oDist['DEPEND_3']
	energy   = oEnergy['DATA']
	velocity = Sqrt(0.5 * energy * eV2J / me) * 100
	
	;DELTAS (half-width of bin)
	pdelta = (oPhi['DELTA_PLUS_VAR'])['DATA'] + (oPhi['DELTA_MINUS_VAR'])['DATA']
	IF pdelta[0] EQ (oPhi['DELTA_MINUS_VAR'])['DATA', 0] THEN pdelta /= 2.0
	tdelta = (oTheta['DELTA_PLUS_VAR'])['DATA'] + (oPhi['DELTA_MINUS_VAR'])['DATA']
	IF tdelta[0] EQ (oTheta['DELTA_MINUS_VAR'])['DATA',0] THEN tdelta /= 2.0
	edelta = (oEnergy['DELTA_PLUS_VAR'])['DATA'] + (oEnergy['DELTA_MINUS_VAR'])['DATA']
	;These two are equivalent to within 0.1%
	vdelta = 1.0/Sqrt(2.0*me*oEnergy['DATA']*eV2J) * edelta*eV2J * 100/2.0
	;dv = (Sqrt(2.0*(oEnergy['DATA']+(oEnergy['DELTA_PLUS_VAR'])['DATA'])*eV2J/me) - $
	;      Sqrt(2.0*(oEnergy['DATA']-(oEnergy['DELTA_MINUS_VAR'])['DATA'])*eV2J/me)) * 100/2.0
	
;-------------------------------------------
; DES Flux of Pixel Nearest to EDI /////////
;-------------------------------------------
	frange = [!Values.F_Infinity, -!Values.F_Infinity]
	ch_color = ['Black', 'Red', 'Green', 'Blue']
	
	;Convert trajectories to unit vectors
	nTraj = N_Elements(edi_traj_vnames)
	oTraj_edi = ObjArr(nTraj)
	oFlux_edi = ObjArr(nTraj)
	oTraj_des = ObjArr(nTraj)
	oFlux_des = ObjArr(nTraj)
	FOR i = 0, nTraj - 1 DO BEGIN
		;Convert trajectories to unit vectors
		oTraj_cart = MrMMS_EDI_sphr2cart(edi_traj_vnames[i])
		
		;Interpolate to FPI time tags
		oFlux = MrVar_Resample(edi_flux_vnames[i], oT, /CACHE, NAME=flux_edi_vnames[i])
		oTraj_cart = MrVar_Resample(oTraj_cart, oT, /CACHE, NAME=traj_edi_vnames[i])
		
		;Convert back to spherical coordinates
		oTraj = MrMMS_EDI_cart2sphr(oTraj_cart)
		
		;Find bin closest to EDI incident trajectories
		iTheta = MrNearestNeighbor(theta, oTraj['DATA',*,1])
		!Null  = Min( Abs(phi - Rebin(oTraj['DATA',*,0], nTime, nPhi)), iPhi, DIMENSION=2 )
		!Null  = Min( Abs(energy - 500), iE, DIMENSION=2 )
		iPhi   = Array_Indices([nTime, nPhi], iPhi, /DIMENSIONS)
		iE     = Array_Indices([nTime, nEnergy], iE,   /DIMENSIONS)
		f_des  = oDist['DATA']
		f_des  = f_des[ Reform(iPhi[0,*]), Reform(iPhi[1,*]), iTheta, Reform(iE[1,*]) ]
		
		;Integrate over bin to get flux
		;   - Flux = v*N = Integral( V * f * v^2 * dv * dPhi * dTheta )
		;   - Integrating over a single bin, so no integral/summation.
		p  = oPhi['DATA', Reform(iPhi[0,*]), Reform(iPhi[1,*])] * !dtor
		t  = oTheta['DATA', iTheta] * !dtor
		v  = velocity[Reform(iE[0,*]), Reform(iE[1,*])]
;		dp = pdelta[Reform(iPhi[1,*]) * !dtor
;		dt = tdelta(iTheta] * !dtor
		dv = vdelta[Reform(iE[0,*]), Reform(iE[1,*])]
		fx = v * Sin(t) * Cos(p) * f_des * v^2.0 * dv
		fy = v * Sin(t) * Sin(p) * f_des * v^2.0 * dv
		fz = v * Cos(t)          * f_des * v^2.0 * dv
		f  = Sqrt(fx^2.0 + fy^2.0 + fz^2.0)

		frange[0] <= Min( [f[Where(f GT 0)], oFlux['DATA', oFlux -> Where(0, /GREATER)]] )
		frange[1] >= Max( [f[Where(f GT 0)], oFlux['DATA', oFlux -> Where(0, /GREATER)]] )
		
		;DES Trajectories
		xdes = Sin(t) * Cos(p)
		ydes = Sin(t) * Sin(p)
		zdes = Cos(t)
		
		;Save variables
		oFlux_edi[i] = oFlux
		oTraj_edi[i] = oTraj_cart
		oFlux_des[i] = MrScalarTS(oT, f, /CACHE, NAME=flux_des_vnames[i])
		oTraj_des[i] = MrVectorTS(oT, -[ [xdes], [ydes], [zdes] ], /CACHE, NAME=traj_des_vnames[i])
		
		;Attributes
		name_parts = StrSplit(flux_edi_vnames[i], '_', /EXTRACT)
		IF i EQ 0 THEN oFlux['LABEL'] = 'EDI'
		IF i EQ 0 THEN oFlux['PLOT_TITLE'] = 'Pitch Angle 0'
		IF i EQ 4 THEN oFlux['PLOT_TITLE'] = 'Pitch Angle 180'
		oFlux['COLOR'] = 'Blue'
		oFlux['TITLE'] = 'Flux!CCh' + StrMid(name_parts[2], 4, 1) + ' PA' + name_parts[3] + '!C(s$\up-1$ cm$\up-2$)'
		oFlux['LOG']   = 1B
		
		IF i EQ 0 THEN (oFlux_des[i])['LABEL'] = 'DES'
		(oFlux_des[i])['COLOR'] = 'Red'
		
		;Angle between EDI and DES trajectories
		oDot = oTraj_edi[i] -> Dot(oTraj_des[i], /CACHE, NAME=traj_dot_vnames[i])
		oDot -> SetData, ACos(oDot['DATA']) * !radeg
		IF i EQ 0  THEN oDot['PLOT_TITLE'] = 'Angle between EDI and DES Nearest Pixel'
		oDot['COLOR'] = ch_color[i MOD 4]
		oDot['LABEL'] = 'Ch' + String( (i MOD 4)+1, FORMAT='(i1)')
		oDot['TITLE'] = '$\theta$!C(Deg)'
	ENDFOR
	
	;Set axis range
	FOR i = 0, N_Elements(oFlux_edi) - 1 DO BEGIN
		(oFlux_edi[i])['AXIS_RANGE'] = frange
		(oFlux_des[i])['AXIS_RANGE'] = frange
	ENDFOR
	
;-------------------------------------------
; Plot Fluxes //////////////////////////////
;-------------------------------------------
	w1 = MrVar_PlotTS( Transpose(Reform(flux_edi_vnames, 4, 2)), $
	                   LAYOUT = [2, 4], $
	                   /NO_REFRESH, $
	                   XSIZE = 800, $
	                   YSIZE = 700 )
	w2 = MrVar_OPlotTS( Transpose(Reform(flux_edi_vnames, 4, 2)), $
	                    Transpose(Reform(flux_des_vnames, 4, 2)) )
	
	w1.name = 'EDI-DES Pixel Flux'
	w1[0] -> SetLayout, [1,1]
	w1 -> TrimLayout
	w1 -> Refresh
	
;-------------------------------------------
; Plot Trajectories ////////////////////////
;-------------------------------------------
	
	w2 = MrVar_PlotTS( traj_dot_vnames[[0,4]], $
	                   /NO_REFRESH, $
	                   XSIZE = 600, $
	                   YSIZE = 400 )
	w2 = MrVar_OPlotTS( traj_dot_vnames[0], traj_dot_vnames[1:3] )
	w2 = MrVar_OPlotTS( traj_dot_vnames[4], traj_dot_vnames[5:7] )
	
	w2.name = 'EDI-DES Pixel Traj'
	w2[0] -> SetLayout, [1,1]
	w2 -> TrimLayout
	w2 -> Refresh
	
;-------------------------------------------
; Scatter Plot /////////////////////////////
;-------------------------------------------
	w3 = MrWindow( ASPECT  = 1.0, $
	               LAYOUT  = [2,4], $
	               REFRESH = 0B, $
	               XSIZE   = 600, $
	               YSIZE   = 700 )
	
	FOR i = 0, N_Elements(flux_edi_vnames) - 1 DO BEGIN
		chan = String( (i MOD 4) + 1, FORMAT='(i1)')
		CASE i OF
			0: title = 'PA0 Ch' + chan
			1: title = 'PA180 Ch' + chan
			ELSE: title = 'Ch' + chan
		ENDCASE
		
		;Plot the variables
		p = MrVar_Plot( flux_edi_vnames[i], flux_des_vnames[i], $
		                /CURRENT, $
		                TITLE  = title, $
		                XTITLE = 'EDI', $
		                YTITLE = 'DES' )
		
		;Fit a line
		p -> GetData, x, y
		idx = where(Finite(x) AND x NE 0 AND y NE 0, count)
		fit = LADFit(x[idx], y[idx], ABSDEV=err)
		eqn = String(Reverse(fit), FORMAT='(%"m=%10.3e!Cb=%10.4e")')
		txt = MrText( 0.1, 0.8, eqn, $
		              CHARSIZE = 1.0, $
		              NAME     = 'Fit: ' + flux_edi_vnames[i] + '-' + flux_des_vnames[i], $
		              /RELATIVE, $
		              TARGET   = p )
	ENDFOR
	
	w3.name = 'EDI-DES Pixel Scatter'
	w3[0] -> SetLayout, [1,1]
	w3 -> TrimLayout
	w3 -> Refresh
	
;-------------------------------------------
; Save Figures /////////////////////////////
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
		fname1 = StrJoin( [sc, 'edi', mode, level, 'flux-des-pixel'], '_' )
		fname2 = StrJoin( [sc, 'edi', mode, level, 'flux-des-traj'], '_' )
		fname3 = StrJoin( [sc, 'edi', mode, level, 'flux-des-scatter'], '_' )
		fname1 = FilePath( fname1, ROOT_DIR=output_dir )
		fname2 = FilePath( fname2, ROOT_DIR=output_dir )
		fname3 = FilePath( fname3, ROOT_DIR=output_dir )
		
		;Save the figure
		f1_out = MrVar_PlotTS_Save( w1, fname1, output_ext )
		f2_out = MrVar_PlotTS_Save( w2, fname2, output_ext )
		f3_out = MrVar_PlotTS_Save( w3, fname3, output_ext )
	ENDIF
	
	RETURN, w3
END