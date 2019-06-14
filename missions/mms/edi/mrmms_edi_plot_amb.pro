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
;       OPTDESC:            in, optional, type=string, default=''
;                           Optional descriptor of the data. Options are:
;                               {'amb-alt-cc' | 'amb-alt-oc' | 'amb-alt-oom' | 'amb-alt-oob'}
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
FUNCTION edi_amb_skip, optdesc, names
	Compile_Opt idl2
	On_Error, 2
	
	tf_skip = 0B
	
	;Sample variable
	CASE 1 OF
		StRegEx(optdesc, '(^amb$|amb-pm2)', /BOOLEAN, /FOLD_CASE): theName = names.flux.fa_0[0]
		StRegEx(optdesc, '(amb-alt)',       /BOOLEAN, /FOLD_CASE): theName = names.flux.fa_0[0]
		StRegEx(optdesc, '(amb-perp)',      /BOOLEAN, /FOLD_CASE): theName = names.flux.perp_gdu1[0]
		StRegEx(optdesc, 'efield',          /BOOLEAN, /FOLD_CASE): RETURN, 1B
		ELSE: Message, 'Unknown optional descriptor: "' + optdesc + '".'
	ENDCASE
	
	;Is the variable in the cache (i.e. has it been skipped previously)?
	oVar = MrVar_Get(thename, COUNT=nVar)
	IF nVar EQ 0 THEN RETURN, 1B
	
	;Are there data in the interval?
	oTime  = oVar['TIMEVAR']
	trange = MrVar_GetTRange('TT2000')
	
	IF N_Elements(oTime) EQ 2 THEN BEGIN
		IF oTime['DATA', 0, 'TT2000'] GT trange[1] || oTime['DATA', 1, 'TT2000'] LT trange[0] $
			THEN tf_skip = 1B
	ENDIF
	
	;Remove variables
	IF tf_skip THEN BEGIN
		;Field-Aligned
		IF StRegEx(optdesc, '(^amb$|amb-pm2|amb-alt)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
			vars = [names.flux.fa_0, names.flux.fa_180, names.traj.fa_0, names.traj.fa_180]
			tf_cached = MrVar_IsCached(vars)
			FOR i = 0, N_Elements(vars) - 1 DO BEGIN
				IF tf_cached[i] THEN MrVar_Delete, vars[i]
			ENDFOR
		ENDIF
		
		;Perpendicular
		IF StRegEx(optdesc, '(amb-perp|amb-alt)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
			vars = [names.flux.perp_gdu1, names.flux.perp_gdu2, names.traj.perp_gdu1, names.traj.perp_gdu2]
			tf_cached = MrVar_IsCached(vars)
			FOR i = 0, N_Elements(vars) - 1 DO BEGIN
				IF tf_cached[i] THEN MrVar_Delete, vars[i]
			ENDFOR
		ENDIF
	ENDIF
	
	RETURN, tf_skip
END


FUNCTION edi_alt_splice, cts_0_vnames, cts_90_gdu1_vnames, cts_90_gdu2_vnames, cts_180_vnames
	Compile_Opt idl2
	On_Error, 2
	
	;EDI is in either field-aligned or perpendicular mode at any one time
	c_fa   = MrVar_Get(cts_0_vnames[0])
	c_perp = MrVar_Get(cts_90_gdu1_vnames[0])
	t_fa   = c_fa['TIME', 'SSM']
	t_perp = c_perp['TIME', 'SSM']
	n_fa   = N_Elements(t_fa)
	n_perp = N_Elements(t_perp)
	
	;Combine times
	time   = [t_fa, t_perp]
	isort  = Sort(time)
	time   = time[isort]
	i_fa = Where(isort LT n_fa, COMPLEMENT=i_perp)
	
	;Splice together the data
	cts_vnames = [cts_0_vnames, cts_90_gdu1_vnames, cts_90_gdu2_vnames, cts_180_vnames]
	nvars = N_Elements(cts_vnames)
	oVars = ObjArr(nvars)
	FOR i = 0, nvars - 1 DO BEGIN
		;Splice data together
		cts = FltArr(n_fa + n_perp)
		oCts = MrVar_Get(cts_vnames[i])
		IF StRegEx(oCts.name, '90', /BOOLEAN) $
			THEN cts[i_perp] = oCts['DATA'] $
			ELSE cts[i_fa]   = oCts['DATA']
		
		;Create a new variable
		oVars[i] = MrScalarTS( time, cts, $
;		                       /CACHE, $
;		                       NAME   = vnames[i], $
		                       T_REF  = c_fa['TIME', 0], $
		                       T_TYPE = 'SSM' )
		oCts -> CopyAttrTo, oVars[i]
	ENDFOR
	
	RETURN, oVars
END


;+
;
;-
FUNCTION edi_cts_srvy, optdesc, names
	Compile_Opt idl2
	On_Error, 2
	
	;FIELD-ALIGNED & CENTERED
	IF StRegEx(optdesc, '(^amb$|amb-alt-c.*)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
		;New Names
		fa_0_vname   = StrMid(names.flux.fa_0[3], 0, 13)   + StrMid(names.flux.fa_0[3], 14)
		fa_180_vname = StrMid(names.flux.fa_180[3], 0, 13) + StrMid(names.flux.fa_180[3], 14)
		
		;Parallel
		oFlux2_0 = MrVar_Get(names.flux.fa_0[3])
		oFlux3_0 = MrVar_Get(names.flux.fa_0[4])
		oFlux_0  = (oFlux2_0 + oFlux3_0) / 2.0
		oFlux_0 -> SetName, fa_0_vname
		oFlux_0 -> Cache
		oFlux2_0 -> CopyAttrTo, oFlux_0
		oFlux_0['CATDESC'] = 'EDI field-aligned ambient data averaged over the two ' + $
		                     'channels centered on the magnetic field direction.'
		
		;Anti-parallel
		oFlux2_180 = MrVar_Get(names.flux.fa_180[3])
		oFlux3_180 = MrVar_Get(names.flux.fa_180[4])
		oFlux_180  = (oFlux2_180 + oFlux3_180) / 2.0
		oFlux_180 -> SetName, fa_180_vname
		oFlux_180 -> Cache
		oFlux2_180 -> CopyAttrTo, oFlux_180
		oFlux_180['CATDESC'] = 'EDI anti-field-aligned ambient data averaged over the two ' + $
		                       'channels centered on the magnetic field direction.'
		
		;Store names
		names.flux.fa_0[-1]   = fa_0_vname
		names.flux.fa_180[-1] = fa_180_vname
		
	;FIELD-ALIGNED & ONE-SIDED
	;   - amb-pm2, amb-alt-o*
	ENDIF ELSE BEGIN
		names.flux.fa_0[-1]   = names.flux.fa_0[1]
		names.flux.fa_180[-1] = names.flux.fa_180[1]
	ENDELSE
	
	;PERPENDICULAR & CENTERED
	IF StRegEx(optdesc, '(amb-alt-.c|amb-perp-c).*', /BOOLEAN, /FOLD_CASE) THEN BEGIN
		;New Names
		perp_gdu1_vname = StrMid(names.flux.perp_gdu1[3], 0, 13) + StrMid(names.flux.perp_gdu1[3], 14)
		perp_gdu2_vname = StrMid(names.flux.perp_gdu2[3], 0, 13) + StrMid(names.flux.perp_gdu2[3], 14)
		
		;GDU1
		oFlux2_gdu1 = MrVar_Get(names.flux.perp_gdu1[3])
		oFlux3_gdu1 = MrVar_Get(names.flux.perp_gdu1[4])
		oFlux_gdu1  = (oFlux2_gdu1 + oFlux3_gdu1) / 2.0
		oFlux_gdu1 -> SetName, perp_gdu1_vname
		oFlux_gdu1 -> Cache
		oFlux2_gdu1 -> CopyAttrTo, oFlux_gdu1
		oFlux_gdu1['CATDESC'] = 'EDI perpendicular ambient data averaged over the two ' + $
		                        'channels centered on the direction perpendicular to ' + $
		                        'the magnetic field.'
		
		;GDU2
		oFlux2_gdu2 = MrVar_Get(names.flux.perp_gdu2[3])
		oFlux3_gdu2 = MrVar_Get(names.flux.perp_gdu2[4])
		oFlux_gdu2  = (oFlux2_gdu2 + oFlux3_gdu2) / 2.0
		oFlux_gdu2 -> SetName, perp_gdu2_vname
		oFlux_gdu2 -> Cache
		oFlux2_gdu2 -> CopyAttrTo, oFlux_gdu2
		oFlux_gdu2['CATDESC'] = 'EDI perpendicular ambient data averaged over the two ' + $
		                        'channels centered on the direction perpendicular to ' + $
		                        'the magnetic field.'
		
		names.flux.perp_gdu1[-1] = perp_gdu1_vname
		names.flux.perp_gdu2[-1] = perp_gdu2_vname
	
	;PERPENDICULAR & ONE-SIDED
	;   - amb-alt-.o., amb-perp-o.
	ENDIF ELSE BEGIN
		names.flux.perp_gdu1[-1] = names.flux.perp_gdu1[1]
		names.flux.perp_gdu2[-1] = names.flux.perp_gdu2[1]
	ENDELSE
END


;+
;
;-
FUNCTION edi_traj_to_pa, oT, traj
	Compile_Opt idl2
	On_Error, 2
	
	;Rotate to field-aligned system
	oTraj = MrVar_Get(traj)
	oVec  = MrMMS_EDI_sphr2cart(oTraj)
	oVec  = oT ## oVec
	oTemp = MrMMS_EDI_cart2sphr(oVec)
	
	;Extract gyrophase angle
	oGP = oTemp[*,0]
	oGP['AXIS_RANGE']   = [-180.0, 180.0]
	oGP['COLOR']        = 'Black'
	oGP['LABEL']        = '$\phi$'
	oGP['TITLE']        = '!7' + String(117B) + '!X$\downB$!C(deg)' ;phi
	oGP['UNITS']        = 'deg'
	oGP['TICKINTERVAL'] = 90.0
	
	;Extract pitch angle
	oPA = oTemp[*,1]
	oPA['AXIS_RANGE']   = [0, 180.0]
	oPA['COLOR']        = 'Black'
	oPA['LABEL']        = '$\theta$'
	oPA['TITLE']        = '!7' + String(104B) + '!X$\downB$!C(deg)' ;theta
	oPA['UNITS']        = 'deg'
	oPA['TICKINTERVAL'] = 90.0
	
	Obj_Destroy, oTemp
	
	;New names
	;   - amb-pm2:     oTraj.name = mms1_edi_traj1_gse_0_srvy_l2_amb_pm2
	;   - amb-perp-ob: oTraj.name = mms1_edi_traj1_gse_90_gdu1_srvy_l2_amb_perp_ob
	CASE 1 OF
		StRegEx(oTraj.name,   '_0_', /BOOLEAN): pos = 20
		StRegEx(oTraj.name,  '_90_', /BOOLEAN): pos = 21
		StRegEx(oTraj.name, '_180_', /BOOLEAN): pos = 22
		ELSE: Message, 'Trajectory name not recognized: "' + oTraj.name + '".'
	ENDCASE
	pa_vname = StrMid(oTraj.name, 0, 9) + 'pa' + StrMid(oTraj.name, 13, 1) + StrMid(oTraj.name, 18)
	ga_vname = StrMid(oTraj.name, 0, 9) + 'ga' + StrMid(oTraj.name, 13, 1) + StrMid(oTraj.name, 18)
	
	;Name and cache
	oPA -> SetName, pa_vname
	oGP -> SetName, ga_vname
	oPA -> Cache
	oGP -> Cache
	
	RETURN, [pa_vname, ga_vname]
END


;+
;
;-
FUNCTION edi_traj_fac, b, traj
	Compile_Opt idl2
	On_Error, 2
	
	;Field-Aligned coordinate transformation
	oB = MrVar_Get(b)
	oB = oB -> Normalize()
	oB = MrVar_Resample(oB, traj[0])
	oT = MrVar_FAC(oB, '', 'CROSSX')
	
	;Rotate
	nTraj = N_Elements(traj)
	pa_vnames = StrArr(nTraj)
	ga_vnames = StrArr(nTraj)
	FOR i = 0, N_Elements(traj) - 1 DO BEGIN
		names = edi_traj_to_pa(oT, traj[i])
		pa_vnames[i] = names[0]
		ga_vnames[i] = names[1]
	ENDFOR
	
	RETURN, [[pa_vnames], [ga_vnames]]
END


;+
;
;-
PRO edi_amb_rename, innames, outnames
	Compile_Opt idl2
	On_Error, 2
	
	;Input names
	in = [innames.flux.fa_0, innames.flux.perp_gdu1, innames.flux.perp_gdu2, innames.flux.fa_180, $
	      innames.traj.fa_0, innames.traj.perp_gdu1, innames.traj.perp_gdu2, innames.traj.fa_180]
	out = [outnames.flux.fa_0, outnames.flux.perp_gdu1, outnames.flux.perp_gdu2, outnames.flux.fa_180, $
	       outnames.traj.fa_0, outnames.traj.perp_gdu1, outnames.traj.perp_gdu2, outnames.traj.fa_180]

;-------------------------------------------
; Rename ///////////////////////////////////
;-------------------------------------------
	FOR i = 0, N_Elements(in) - 1 DO BEGIN
		oVar = MrVar_Get(in[i], COUNT=nVar)
		IF nVar GT 0 THEN oVar -> SetName, out[i]
	ENDFOR
END


;+
;
;-
FUNCTION edi_amb_setup
	Compile_Opt idl2
	On_Error, 2
	
;-------------------------------------------
; Create False Variables ///////////////////
;-------------------------------------------
	oFlux = MrScalarTS( MrVar_GetTRange(), [0,0] )
	oFlux -> SetName, 'Flux'
	oFlux['AXIS_RANGE'] = [1e4, 1e7]
	oFlux['LOG']        = 1B
	oFlux['NODATA']     = 1B
	oFlux['TITLE']      = 'Flux!C(cm$\up-2$ s$\up-1$ sr$\up-1$)'
	
	oPA1 = oFlux -> Copy('PA1')
	oPA1['AXIS_RANGE']    = [0, 180]
	oPA1['LOG']           = 0B
	oPA1['NODATA']        = 1B
	oPA1['TITLE']         = 'PA!C(deg)'
	oPA1['YTICKINTERVAL'] = 90
	
	oGA1 = oPA1 -> Copy('GA1')
	oGA1['AXIS_RANGE']    = [-180, 180]
	oPA1['LOG']           = 0B
	oGA1['NODATA']        = 1B
	oGA1['TITLE']         = 'GA!C(deg)'
	oGA1['YTICKINTERVAL'] = 90
	
;-------------------------------------------
; Survey Data //////////////////////////////
;-------------------------------------------
	
	
	win_srvy = MrWindow( LAYOUT  = [1,3], $
	                     NAME    = 'edi-amb-srvy', $
	                     REFRESH = 0, $
	                     XSIZE   = 550, $
	                     YGAP    = 0.5, $
	                     YSIZE   = 600 )
	
	w_srvy = MrVar_PlotTS( [oFlux, oPA1, oGA1], /CURRENT )
	w_srvy[0] -> SetLayout, [1,1]
	w_srvy -> TrimLayout
	
	
;-------------------------------------------
; Fluxes ///////////////////////////////////
;-------------------------------------------
	
	win_flux = MrWindow( LAYOUT  = [1,5], $
	                     NAME    = 'edi-amb-flux', $
	                     REFRESH = 0, $
	                     XSIZE   = 550, $
	                     YGAP    = 0.5, $
	                     YSIZE   = 700 )
	
	oFluxN = ObjArr(4)
	FOR i = 0, 3 DO BEGIN
		ch = String(i+1, FORMAT='(i1)')
		oFluxN[i]            = oFlux -> Copy('Flux'+ch)
		(oFluxN[i])['TITLE'] = 'Flux Ch'+ch+'!C(cm$\up-2$ s$\up-1$ sr$\up-1$)'
	ENDFOR
;	oFlux2 = oFlux -> Copy('Flux2')
;	oFlux3 = oFlux -> Copy('Flux3')
;	oFlux4 = oFlux -> Copy('Flux4')
	w_flux = MrVar_PlotTS( [oFlux, oFluxN], /CURRENT )
	Obj_Destroy, [oFlux, oFluxN]
	
	w_flux[0] -> SetLayout, [1,1]
	w_flux -> TrimLayout
	
;-------------------------------------------
; Pitch Angles /////////////////////////////
;-------------------------------------------
	
	win_pa = MrWindow( LAYOUT  = [1,4], $
	                   NAME    = 'edi-amb-pa', $
	                   REFRESH = 0, $
	                   XSIZE   = 550, $
	                   YGAP    = 0.5, $
	                   YSIZE   = 700 )
	
	oPA2 = oPA1 -> Copy('PA2')
	oPA3 = oPA1 -> Copy('PA3')
	oPA4 = oPA1 -> Copy('PA4')
	w_pa = MrVar_PlotTS( [oPA1, oPA2, oPA3, oPA4], /CURRENT )
	Obj_Destroy, [oPA1, oPA2, oPA3, oPA4]
	
	w_pa[0] -> SetLayout, [1,1]
	w_pa -> TrimLayout
	
;-------------------------------------------
; Gyrophase Angles /////////////////////////
;-------------------------------------------
	
	win_ga = MrWindow( LAYOUT  = [1,4], $
	                   NAME    = 'edi-amb-ga', $
	                   REFRESH = 0, $
	                   XSIZE   = 550, $
	                   YGAP    = 0.5, $
	                   YSIZE   = 700 )
	
	oGA2 = oGA1 -> Copy('GA2')
	oGA3 = oGA1 -> Copy('GA3')
	oGA4 = oGA1 -> Copy('GA4')
	w_ga = MrVar_PlotTS( [oGA1, oGA2, oGA3, oGA4], /CURRENT )
	Obj_Destroy, [oGA1, oGA2, oGA3, oGA4]
	
	w_ga[0] -> SetLayout, [1,1]
	w_ga -> TrimLayout
	
;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, [w_srvy, w_flux, w_pa, w_ga]
END


;+
;
;-
FUNCTION edi_amb_varnames, sc, level, coords, optdesc, $
OPTNAMES=optnames
	Compile_Opt idl2
	On_Error, 2
	
	ch = ['1', '2', '3', '4']
	type = level EQ 'l2' ? 'flux' : 'counts'
	
	flux = { fa_0:      StrArr(6), $
	         perp_gdu1: StrArr(6), $
	         perp_gdu2: StrArr(6), $
	         fa_180:    StrArr(6) }
	         
	traj = { fa_0:      StrArr(5), $
	         perp_gdu1: StrArr(5), $
	         perp_gdu2: StrArr(5), $
	         fa_180:    StrArr(5) }
	         
	pa = { fa_0:      StrArr(5), $
	       perp_gdu1: StrArr(5), $
	       perp_gdu2: StrArr(5), $
	       fa_180:    StrArr(5) }
	         
	ga = { fa_0:      StrArr(5), $
	       perp_gdu1: StrArr(5), $
	       perp_gdu2: StrArr(5), $
	       fa_180:    StrArr(5) }
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;SRVY
	cts_0_srvy_vname        = StrJoin( [sc, 'edi', type+ch[0],   '0', 'srvy', level], '_' )
	cts_90_gdu1_srvy_vname  = StrJoin( [sc, 'edi', type+ch[0],  '90', 'gdu1', 'srvy', level], '_' )
	cts_90_gdu2_srvy_vname  = StrJoin( [sc, 'edi', type+ch[0],  '90', 'gdu2', 'srvy', level], '_' )
	cts_180_srvy_vname      = StrJoin( [sc, 'edi', type+ch[0], '180', 'srvy', level], '_' )
	traj_0_srvy_vname       = StrJoin( [sc, 'edi', 'traj'+ch[0], coords,   '0', 'srvy', level], '_' )
	traj_90_gdu1_srvy_vname = StrJoin( [sc, 'edi', 'traj'+ch[0], coords,  '90', 'gdu1', 'srvy', level], '_' )
	traj_90_gdu2_srvy_vname = StrJoin( [sc, 'edi', 'traj'+ch[0], coords,  '90', 'gdu2', 'srvy', level], '_' )
	traj_180_srvy_vname     = StrJoin( [sc, 'edi', 'traj'+ch[0], coords, '180', 'srvy', level], '_' )
	
	;BRST
	cts_0_vnames        = StrJoin( [sc, 'edi', type,   '0', 'brst', level], '_' )
	cts_90_gdu1_vnames  = StrJoin( [sc, 'edi', type,  '90', 'gdu1', 'brst', level], '_' )
	cts_90_gdu2_vnames  = StrJoin( [sc, 'edi', type,  '90', 'gdu2', 'brst', level], '_' )
	cts_180_names       = StrJoin( [sc, 'edi', type, '180', 'brst', level], '_' )
	traj_0_vnames       = StrJoin( [sc, 'edi', 'traj', coords,   '0', 'brst', level], '_' )
	traj_90_gdu1_vnames = StrJoin( [sc, 'edi', 'traj', coords,  '90', 'gdu1', 'brst', level], '_' )
	traj_90_gdu2_vnames = StrJoin( [sc, 'edi', 'traj', coords,  '90', 'gdu2', 'brst', level], '_' )
	traj_180_vnames     = StrJoin( [sc, 'edi', 'traj', coords, '180', 'brst', level], '_' )
	
	cts_0_vnames        = StrMid(cts_0_vnames, 0, 13)        + ch + StrMid(cts_0_vnames, 13)
	cts_90_gdu1_vnames  = StrMid(cts_90_gdu1_vnames, 0, 13)  + ch + StrMid(cts_90_gdu1_vnames, 13)
	cts_90_gdu2_vnames  = StrMid(cts_90_gdu2_vnames, 0, 13)  + ch + StrMid(cts_90_gdu2_vnames, 13)
	cts_180_vnames      = StrMid(cts_180_names, 0, 13)       + ch + StrMid(cts_180_names, 13)
	traj_0_vnames       = StrMid(traj_0_vnames, 0, 13)       + ch + StrMid(traj_0_vnames, 13)
	traj_90_gdu1_vnames = StrMid(traj_90_gdu1_vnames, 0, 13) + ch + StrMid(traj_90_gdu1_vnames, 13)
	traj_90_gdu2_vnames = StrMid(traj_90_gdu2_vnames, 0, 13) + ch + StrMid(traj_90_gdu2_vnames, 13)
	traj_180_vnames     = StrMid(traj_180_vnames, 0, 13)     + ch + StrMid(traj_180_vnames, 13)
	
;-------------------------------------------
; Output Variable Names ////////////////////
;-------------------------------------------
	;FLUXES
	flux.fa_0      = [cts_0_srvy_vname, cts_0_vnames]
	flux.perp_gdu1 = [cts_90_gdu1_srvy_vname, cts_90_gdu1_vnames]
	flux.perp_gdu2 = [cts_90_gdu2_srvy_vname, cts_90_gdu2_vnames]
	flux.fa_180    = [cts_180_srvy_vname, cts_180_vnames]
	
	;TRAJECTORIES
	traj.fa_0      = [traj_0_srvy_vname, traj_0_vnames]
	traj.perp_gdu1 = [traj_90_gdu1_srvy_vname, traj_90_gdu1_vnames]
	traj.perp_gdu2 = [traj_90_gdu2_srvy_vname, traj_90_gdu2_vnames]
	traj.fa_180    = [traj_180_srvy_vname, traj_180_vnames]
	
	;PITCH ANGLES
	pa.fa_0      = StrMid(traj.fa_0, 0, 9)      + 'pa' + StrMid(traj.fa_0, 13, 1)      + StrMid(traj.fa_0, 18)
	pa.perp_gdu1 = StrMid(traj.perp_gdu1, 0, 9) + 'pa' + StrMid(traj.perp_gdu1, 13, 1) + StrMid(traj.perp_gdu1, 18)
	pa.perp_gdu2 = StrMid(traj.perp_gdu2, 0, 9) + 'pa' + StrMid(traj.perp_gdu2, 13, 1) + StrMid(traj.perp_gdu2, 18)
	pa.fa_180    = StrMid(traj.fa_180, 0, 9)    + 'pa' + StrMid(traj.fa_180, 13, 1)    + StrMid(traj.fa_180, 18)
	
	;GYROPHASE ANGLES
	ga.fa_0      = StrMid(traj.fa_0, 0, 9)      + 'ga' + StrMid(traj.fa_0, 13, 1)      + StrMid(traj.fa_0, 18)
	ga.perp_gdu1 = StrMid(traj.perp_gdu1, 0, 9) + 'ga' + StrMid(traj.perp_gdu1, 13, 1) + StrMid(traj.perp_gdu1, 18)
	ga.perp_gdu2 = StrMid(traj.perp_gdu2, 0, 9) + 'ga' + StrMid(traj.perp_gdu2, 13, 1) + StrMid(traj.perp_gdu2, 18)
	ga.fa_180    = StrMid(traj.fa_180, 0, 9)    + 'ga' + StrMid(traj.fa_180, 13, 1)    + StrMid(traj.fa_180, 18)

;-------------------------------------------
; Finished /////////////////////////////////
;-------------------------------------------
	
	;Variables
	varnames = { flux: Temporary(flux), $
	             traj: Temporary(traj), $
	             pa:   Temporary(pa), $
	             ga:   Temporary(ga) }
	
	;Append optional descriptor
	IF Arg_Present(optnames) THEN BEGIN
		optnames = varnames
		desc = '_' + StrJoin(StrSplit(optdesc, '-', /EXTRACT), '_')
		FOR i = 0, N_Tags(optnames) - 1 DO BEGIN
			FOR j = 0, N_Tags(optnames.(i)) - 1 DO BEGIN
				optnames.(i).(j) += desc
			ENDFOR
		ENDFOR
	ENDIF
	
	RETURN, varnames
END


;+
;
;-
FUNCTION MrMMS_EDI_Plot_Amb, sc, $
COORDS=coords, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
TRANGE=trange
	Compile_Opt idl2

	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		catch, /CANCEL
		IF N_Elements(wins) GT 0 THEN Obj_Destroy, wins
		MrPrintF, 'LogErr'
		RETURN, Obj_New()
	ENDIF
	
	;Defaults
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords) EQ 0 THEN coords = 'gse'
	IF N_Elements(trange) GT 0 THEN MrVar_SetTRange, trange
	
	;Variable type
	level = 'l2'
	type  = level EQ 'l2' ? 'flux' : 'counts'

;-------------------------------------------
; Get File Types ///////////////////////////
;-------------------------------------------
	files = MrMMS_Get_Filenames(sc, 'edi', ['srvy','brst'], level, OPTDESC='amb-alt-oob', COUNT=count)
	IF count EQ 0 THEN Message, 'No EDI files found.'
		
	;Which files
	MrMMS_Parse_Filename, files, OPTDESC=optdesc
	iUniq    = Uniq(optdesc, Sort(optdesc))
	optdesc  = optdesc[iUniq]

;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	
	fgm_instr  = 'dfg'
	fgm_level  = 'l2pre'
	fgm_coords = coords EQ 'dbcs' ? 'dmpa' : coords
;	fgm_mode   = mode EQ 'brst' ? mode : 'srvy'

;-------------------------------------------
; FGM Variable Names ///////////////////////
;-------------------------------------------
	b_srvy_vname    = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, 'srvy', fgm_level], '_' )
	bvec_srvy_vname = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, 'srvy', fgm_level], '_' )
	bmag_srvy_vname = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, 'srvy', fgm_level], '_' )
	b_brst_vname    = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, 'brst', fgm_level], '_' )
	bvec_brst_vname = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, 'brst', fgm_level], '_' )
	bmag_brst_vname = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, 'brst', fgm_level], '_' )
	
;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		
		;FGM
		MrMMS_FGM_Load_Data, sc, ['brst', 'srvy'], $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = fgm_level, $
		                     VARFORMAT = '*_b_'+fgm_coords+'*'+'_'+fgm_level
		
		IF fgm_level NE 'l2' && ~MrVar_IsCached(b_srvy_vname) THEN BEGIN
			b_srvy_vname     = StrJoin( [sc, fgm_instr,        'srvy', fgm_level, fgm_coords], '_' )
			bvec_srvy_vname  = StrJoin( [sc, fgm_instr, 'vec', 'srvy', fgm_level, fgm_coords], '_' )
			bmag_srvy_vname  = StrJoin( [sc, fgm_instr, 'mag', 'srvy', fgm_level, fgm_coords], '_' )
			b_brst_vname     = StrJoin( [sc, fgm_instr,        'brst', fgm_level, fgm_coords], '_' )
			bvec_brst_vname  = StrJoin( [sc, fgm_instr, 'vec', 'brst', fgm_level, fgm_coords], '_' )
			bmag_brst_vname  = StrJoin( [sc, fgm_instr, 'mag', 'brst', fgm_level, fgm_coords], '_' )
			
			MrMMS_FGM_Load_Data, sc, fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = fgm_level, $
			                     VARFORMAT = '*'+'_'+fgm_level+'_'+fgm_coords
			
			IF ~MrVar_IsCached(b_srvy_vname) $
				THEN Message, 'Unknown FGM variable naming convention.'
		ENDIF
	ENDIF
	
;-------------------------------------------
; Set Up Graphics Windows //////////////////
;-------------------------------------------
	wins = edi_amb_setup()
	
;-------------------------------------------
; Process Each File Type ///////////////////
;-------------------------------------------
	
	;
	;TODO:
	;   X. Load each file type individually (prevent variables from being overwritten)
	;   X. Rename variables to contain optional descriptor
	;   X. Compute equivalent survey counts ("amb" averages channels 2 & 3)
	;   X. Compute pitch and gyrophase angles
	;   X. Splice alternating counts and trajectories
	;
	
	FOR i = 0, N_Elements(optdesc) - 1 DO BEGIN
		;Variable names
		tempnames = edi_amb_varnames(sc, level, coords, optdesc[i], OPTNAMES=vnames)
		
		;Load EDI data
		IF tf_load THEN BEGIN
			MrMMS_Load_Data, sc, 'edi', ['brst', 'srvy'], level, $
			                 OPTDESC   = optdesc[i], $
			                 VARFORMAT = ['*_flux?_0_[sb]r[vs][yt]_*', $
			                              '*_flux?_90_gdu?_[sb]r[vs][yt]_*', $
			                              '*_flux?_180_[sb]r[vs][yt]_*', $
			                              '*_traj?_'+coords+'_*']
			
			;Skip this variable?
			tf_skip = edi_amb_skip(optdesc[i], tempnames)
			IF tf_skip THEN CONTINUE
		
			;Rename variables
			edi_amb_rename, Temporary(tempnames), vnames
		ENDIF
		
		
		;
		; Counts
		;
		oVars = edi_cts_srvy(optdesc[i], vnames)
		
		;
		; Trajectories
		;
		
		;PERPENDICULAR
		IF StRegEx(optdesc[i], 'amb-perp', /BOOLEAN, /FOLD_CASE) THEN BEGIN
			ang_names = edi_traj_fac( bvec_srvy_vname, [vnames.traj.perp_gdu1[0], vnames.traj.perp_gdu2[0]])
			ang_names = edi_traj_fac( bvec_brst_vname, [vnames.traj.perp_gdu1[1:*], vnames.traj.perp_gdu2[1:*]])
		
		;ALTERNATING
		ENDIF ELSE IF StRegEx(optdesc[i], 'amb-alt', /BOOLEAN, /FOLD_CASE) THEN BEGIN
			ang_names = edi_traj_fac(bvec_srvy_vname, [vnames.traj.fa_0[0], vnames.traj.fa_180[0], $
			                                          vnames.traj.perp_gdu1[0], vnames.traj.perp_gdu2[0]])
			ang_names = edi_traj_fac(bvec_brst_vname, [vnames.traj.fa_0[1:*], vnames.traj.fa_180[1:*], $
			                                          vnames.traj.perp_gdu1[1:*], vnames.traj.perp_gdu2[1:*]])
		
		;FIELD-ALIGNED
		ENDIF ELSE BEGIN
			ang_names = edi_traj_fac( bvec_srvy_vname, [vnames.traj.fa_0[0], vnames.traj.fa_180[0]])
			ang_names = edi_traj_fac( bvec_brst_vname, [vnames.traj.fa_0[1:*], vnames.traj.fa_180[1:*]])
		ENDELSE
		
		;
		; Plot
		;
		IF StRegEx(optdesc[i], '(^amb$|amb-pm2|amb-alt)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
			;
			;Update label and color
			;
			
			;BRST
			names_0   = [vnames.flux.fa_0[1:-2], vnames.pa.fa_0[1:-1], vnames.ga.fa_0[1:-1]]
			names_180 = [vnames.flux.fa_180[1:-2], vnames.pa.fa_180[1:-1], vnames.ga.fa_180[1:-1]]
			FOR j = 0, N_Elements(names_0) - 1 DO BEGIN
				oVar = MrVar_Get(names_0[j])
				oVar['COLOR'] = 'Blue'
				oVar['LABEL'] = '0'
				
				oVar = MrVar_Get(names_180[j])
				oVar['COLOR'] = 'Red'
				oVar['LABEL'] = '180'
			ENDFOR
			
			;SRVY
			names_0   = [vnames.flux.fa_0[0], vnames.pa.fa_0[0], vnames.ga.fa_0[0]]
			names_180 = [vnames.flux.fa_180[0], vnames.pa.fa_180[0], vnames.ga.fa_180[0]]
			FOR j = 0, N_Elements(names_0) - 1 DO BEGIN
				oVar = MrVar_Get(names_0[j])
				oVar['COLOR'] = 'Medium Grey'
				oVar['LABEL'] = '0'
				
				oVar = MrVar_Get(names_180[j])
				oVar['COLOR'] = 'Sky Blue'
				oVar['LABEL'] = '180'
			ENDFOR
			
			;SRVY
			wins[0] -> SetCurrent
			p = MrVar_OPlotTS( 'Flux', [vnames.flux.fa_0[0], vnames.flux.fa_180[0]] )
			p = MrVar_OPlotTS( 'PA1', [vnames.pa.fa_0[0], vnames.pa.fa_180[0]] )
			p = MrVar_OPlotTS( 'GA1', [vnames.ga.fa_0[0], vnames.ga.fa_180[0]] )
			
			;FLUXES
			wins[1] -> SetCurrent
			p = MrVar_OPlotTS( 'Flux', vnames.flux.fa_0[[0,5]] )
			p = MrVar_OPlotTS( 'Flux1', [vnames.flux.fa_0[1], vnames.flux.fa_180[1]] )
			p = MrVar_OPlotTS( 'Flux2', [vnames.flux.fa_0[2], vnames.flux.fa_180[2]] )
			p = MrVar_OPlotTS( 'Flux3', [vnames.flux.fa_0[3], vnames.flux.fa_180[3]] )
			p = MrVar_OPlotTS( 'Flux4', [vnames.flux.fa_0[4], vnames.flux.fa_180[4]] )
			
			;PITCH ANGLES
			wins[2] -> SetCurrent
			p = MrVar_OPlotTS( 'PA1', [vnames.pa.fa_0[1], vnames.pa.fa_180[1]] )
			p = MrVar_OPlotTS( 'PA2', [vnames.pa.fa_0[2], vnames.pa.fa_180[2]] )
			p = MrVar_OPlotTS( 'PA3', [vnames.pa.fa_0[3], vnames.pa.fa_180[3]] )
			p = MrVar_OPlotTS( 'PA4', [vnames.pa.fa_0[4], vnames.pa.fa_180[4]] )
			p = MrVar_OPlotTS( 'PA1', [vnames.pa.fa_0[0], vnames.pa.fa_180[0]] )
			
			;GYROPHASE ANGLES
			wins[3] -> SetCurrent
			p = MrVar_OPlotTS( 'GA1', [vnames.ga.fa_0[1], vnames.ga.fa_180[1]] )
			p = MrVar_OPlotTS( 'GA2', [vnames.ga.fa_0[2], vnames.ga.fa_180[2]] )
			p = MrVar_OPlotTS( 'GA3', [vnames.ga.fa_0[3], vnames.ga.fa_180[3]] )
			p = MrVar_OPlotTS( 'GA4', [vnames.ga.fa_0[4], vnames.ga.fa_180[4]] )
			p = MrVar_OPlotTS( 'GA1', [vnames.ga.fa_0[0], vnames.ga.fa_180[0]] )
		ENDIF
		
		IF StRegEx(optdesc[i], '(amb-perp|amb-alt)', /BOOLEAN, /FOLD_CASE) THEN BEGIN
			;
			;Update label and color
			;
			
			;BRST
			names_gdu1 = [vnames.flux.perp_gdu1[1:-2], vnames.pa.perp_gdu1[1:-1], vnames.ga.perp_gdu1[1:-1]]
			names_gdu2 = [vnames.flux.perp_gdu2[1:-2], vnames.pa.perp_gdu2[1:-1], vnames.ga.perp_gdu2[1:-1]]
			FOR j = 0, N_Elements(names_gdu1) - 1 DO BEGIN
				oVar = MrVar_Get(names_gdu1[j])
				oVar['COLOR'] = 'Forest Green'
				oVar['LABEL'] = '90 GDU1'
				
				oVar = MrVar_Get(names_gdu2[j])
				oVar['COLOR'] = 'Magenta'
				oVar['LABEL'] = '90 GDU2'
			ENDFOR
			
			;SRVY
			names_gdu1 = [vnames.flux.perp_gdu1[0], vnames.pa.perp_gdu1[0], vnames.ga.perp_gdu1[0]]
			names_gdu2 = [vnames.flux.perp_gdu2[0], vnames.pa.perp_gdu2[0], vnames.ga.perp_gdu2[0]]
			FOR j = 0, N_Elements(names_gdu1) - 1 DO BEGIN
				oVar = MrVar_Get(names_gdu1[j])
				oVar['COLOR'] = 'Pur4'
				oVar['LABEL'] = '90 GDU1'
				
				oVar = MrVar_Get(names_gdu2[j])
				oVar['COLOR'] = 'Tan7'
				oVar['LABEL'] = '90 GDU2'
			ENDFOR
			
			;SRVY
			wins[0] -> SetCurrent
			p = MrVar_OPlotTS( 'Flux', [vnames.flux.perp_gdu1[0], vnames.flux.perp_gdu2[0]] )
			p = MrVar_OPlotTS( 'PA1', [vnames.pa.perp_gdu1[0], vnames.pa.perp_gdu2[0]] )
			p = MrVar_OPlotTS( 'GA1', [vnames.ga.perp_gdu1[0], vnames.ga.perp_gdu2[0]] )
			
			;FLUXES
			wins[1] -> SetCurrent
			p = MrVar_OPlotTS( 'Flux',  [vnames.flux.perp_gdu1[5], vnames.flux.perp_gdu2[5], $
			                             vnames.flux.perp_gdu1[0], vnames.flux.perp_gdu2[0]] )
			p = MrVar_OPlotTS( 'Flux1', [vnames.flux.perp_gdu1[1], vnames.flux.perp_gdu2[1]] )
			p = MrVar_OPlotTS( 'Flux2', [vnames.flux.perp_gdu1[2], vnames.flux.perp_gdu2[2]] )
			p = MrVar_OPlotTS( 'Flux3', [vnames.flux.perp_gdu1[3], vnames.flux.perp_gdu2[3]] )
			p = MrVar_OPlotTS( 'Flux4', [vnames.flux.perp_gdu1[4], vnames.flux.perp_gdu2[4]] )
			
			;PITCH ANGLES
			wins[2] -> SetCurrent
			p = MrVar_OPlotTS( 'PA1', [vnames.pa.perp_gdu1[1], vnames.pa.perp_gdu2[1]] )
			p = MrVar_OPlotTS( 'PA2', [vnames.pa.perp_gdu1[2], vnames.pa.perp_gdu2[2]] )
			p = MrVar_OPlotTS( 'PA3', [vnames.pa.perp_gdu1[3], vnames.pa.perp_gdu2[3]] )
			p = MrVar_OPlotTS( 'PA4', [vnames.pa.perp_gdu1[4], vnames.pa.perp_gdu2[4]] )
			p = MrVar_OPlotTS( 'PA1', [vnames.pa.perp_gdu1[0], vnames.pa.perp_gdu2[0]] )
			
			;GYROPHASE ANGLES
			wins[3] -> SetCurrent
			p = MrVar_OPlotTS( 'GA1', [vnames.ga.perp_gdu1[1], vnames.ga.perp_gdu2[1]] )
			p = MrVar_OPlotTS( 'GA2', [vnames.ga.perp_gdu1[2], vnames.ga.perp_gdu2[2]] )
			p = MrVar_OPlotTS( 'GA3', [vnames.ga.perp_gdu1[3], vnames.ga.perp_gdu2[3]] )
			p = MrVar_OPlotTS( 'GA4', [vnames.ga.perp_gdu1[4], vnames.ga.perp_gdu2[4]] )
			p = MrVar_OPlotTS( 'GA1', [vnames.ga.perp_gdu1[0], vnames.ga.perp_gdu2[0]] )
		ENDIF
	ENDFOR
	
;-------------------------------------------
; Save /////////////////////////////////////
;-------------------------------------------
	
	FOR i = 0, N_Elements(wins) - 1 DO BEGIN
		wins[i].oxmargin = [13,10]
		wins[i] -> Refresh
	ENDFOR
	
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
		ENDIF ELSE IF ~File_Test(output_dir, /DIRECTORY) THEN BEGIN
			MrPrintF, 'LogText', 'Creating directory: "' + output_dir + '".'
			File_MKDir, output_dir
		ENDIF
		
		;File name
		fname1   = StrJoin( [sc, 'edi', 'all', level, 'amb-srvy'], '_' )
		fname2   = StrJoin( [sc, 'edi', 'all', level, 'amb-flux'], '_' )
		fname3   = StrJoin( [sc, 'edi', 'all', level, 'amb-pa'], '_' )
		fname4   = StrJoin( [sc, 'edi', 'all', level, 'amb-ga'], '_' )
		
		fname1   = FilePath( fname1, ROOT_DIR=output_dir )
		fname2   = FilePath( fname2, ROOT_DIR=output_dir )
		fname3   = FilePath( fname3, ROOT_DIR=output_dir )
		fname4   = FilePath( fname4, ROOT_DIR=output_dir )
		
		;Save the figure
		fout1 = MrVar_PlotTS_Save( wins[0], fname1, output_ext )
		fout2 = MrVar_PlotTS_Save( wins[1], fname2, output_ext )
		fout3 = MrVar_PlotTS_Save( wins[2], fname3, output_ext )
		fout4 = MrVar_PlotTS_Save( wins[3], fname4, output_ext )
	ENDIF
	
	RETURN, wins
END