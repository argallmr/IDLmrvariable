; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_EDI_Amb
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
; PURPOSE:
;+
;   Load EDI data from source files.
;
; :Categories:
;   MMS
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
;       2017/03/09  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Update variables for AMB data.
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' 
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'ql' | 'l2'}
;       COORDS:     in, optional, type=string, default='gse'
;                   The coordinate system in which trajectories are to be loaded. Options
;                       are: {'gse' | 'gsm'}
;
; :Keywords:
;       SUFFIX:     in, optional, type=string, default=''
;                   A suffix to be appendend to the end of EDI variable names.
;-
PRO MrMMS_EDI_Load_Data_Amb, sc, mode, level, coords, $
SUFFIX=suffix
	Compile_Opt idl2
	On_Error, 2
		
	;EDI
	instr  = 'edi'
	
	;Variable names
	channels      = ['1', '2', '3', '4']
	flux_0_vname   = sc + '_' + instr + '_' + 'flux' + channels + '_'                +   '0' + '_' + mode + '_' + level + suffix
	flux_180_vname = sc + '_' + instr + '_' + 'flux' + channels + '_'                + '180' + '_' + mode + '_' + level + suffix
	traj_0_vname   = sc + '_' + instr + '_' + 'traj' + channels + '_' + coords + '_' +   '0' + '_' + mode + '_' + level + suffix
	traj_180_vname = sc + '_' + instr + '_' + 'traj' + channels + '_' + coords + '_' + '180' + '_' + mode + '_' + level + suffix
	
	;Output names
	f_0_vname   = StrJoin( [sc, instr, 'flux',   '0',         mode, level], '_' ) + suffix
	f_180_vname = StrJoin( [sc, instr, 'flux', '180',         mode, level], '_' ) + suffix
	t_0_vname   = StrJoin( [sc, instr, 'traj',   '0', coords, mode, level], '_' ) + suffix
	t_180_vname = StrJoin( [sc, instr, 'traj', '180', coords, mode, level], '_' ) + suffix
	
;-------------------------------------------
; 0-Degree Flux ////////////////////////////
;-------------------------------------------

	;Mean of middle two channels
	;   - B-field is centered between channels 2 & 3
	IF Array_Equal(MrVar_IsCached(flux_0_vname[1:2]), 1) THEN BEGIN
		;Get the data
		oFlux2_0 = MrVar_Get(flux_0_vname[1])
		oFlux3_0 = MrVar_Get(flux_0_vname[2])
		
		;Average the channels
		oFlux_0  = (oFlux2_0 + oFlux3_0) / 2.0
		oFlux_0 -> SetName, f_0_vname
		oFlux_0 -> Cache
		
		;Attributes
		oFlux2_0 -> CopyAttrTo, oFlux_0
		oFlux_0['CATDESC'] = 'Flux for electrons with trajectories given by traj_0. ' + $
		                     'This is the average flux from channels 2 and 3.'
		oFlux_0 -> RemoveAttr, 'CDF_NAME'
	ENDIF
	
;-------------------------------------------
; 0-Degree Trajectories ////////////////////
;-------------------------------------------
	IF Array_Equal(MrVar_IsCached(traj_0_vname[1:2]), 1) THEN BEGIN
		oTraj2_0 = MrVar_Get(traj_0_vname[1])
		oTraj3_0 = MrVar_Get(traj_0_vname[2])
		
		;Average the angles
		;   - Average the x- and y-components, then find the resulting angles.
		;   - The factor of two is ommitted because it would otherwise cancel anyway
		deg2rad = !dpi / 180D
		rad2deg = 180D / !dpi
		sinA    = Sin( oTraj2_0['DATA'] * deg2rad ) + Sin( oTraj3_0['DATA'] * deg2rad )
		cosA    = Cos( oTraj2_0['DATA'] * deg2rad ) + Cos( oTraj3_0['DATA'] * deg2rad )
		traj    = Float( ATan( Temporary(sinA), Temporary(cosA) ) * rad2deg )
		
		;Create a new variable
		oTraj_0 = MrTimeSeries( oTraj2_0['TIMEVAR'], traj, /CACHE, NAME=t_0_vname, /NO_COPY )
		
		;Attributes
		oTraj2_0 -> CopyAttrTo, oTraj_0
		oTraj_0['CATDESC'] = 'Trajectory of 0-degree pitch angle electrons in GSE coordinates. ' + $
		                     'This is an average of the trajectories from channels 2 and 3.'
		oTraj_0 -> RemoveAttr, 'CDF_NAME'
	ENDIF
	
;-------------------------------------------
; 180-Degree Flux //////////////////////////
;-------------------------------------------
	
	;Mean of middle two channels
	;   - B-field is centered between channels 2 & 3
	IF Array_Equal(MrVar_IsCached(flux_180_vname[1:2]), 1) THEN BEGIN
		;Get the data
		oFlux2_180 = MrVar_Get(flux_180_vname[1])
		oFlux3_180 = MrVar_Get(flux_180_vname[2])
		
		;Average the channels
		oFlux_180  = (oFlux2_180 + oFlux3_180) / 2.0
		oFlux_180 -> SetName, f_180_vname
		oFlux_180 -> Cache
		
		;Attributes
		oFlux2_180 -> CopyAttrTo, oFlux_180
		oFlux_180['CATDESC'] = 'Flux for electrons with trajectories given by traj_180. ' + $
		                       'This is the average flux from channels 2 and 3.'
		oFlux_180 -> RemoveAttr, 'CDF_NAME'
	ENDIF
	
;-------------------------------------------
; 180-Degree Trajectories //////////////////
;-------------------------------------------
	IF Array_Equal(MrVar_IsCached(traj_180_vname[1:2]), 1) THEN BEGIN
		oTraj2_180 = MrVar_Get(traj_180_vname[1])
		oTraj3_180 = MrVar_Get(traj_180_vname[2])
		
		;Average the angles
		;   - Average the x- and y-components, then find the resulting angles.
		;   - The factor of two is ommitted because it would otherwise cancel anyway
		deg2rad = !dpi / 180D
		rad2deg = 180D / !dpi
		sinA    = Sin( oTraj2_180['DATA'] * deg2rad ) + Sin( oTraj3_180['DATA'] * deg2rad )
		cosA    = Cos( oTraj2_180['DATA'] * deg2rad ) + Cos( oTraj3_180['DATA'] * deg2rad )
		traj    = Float( ATan( Temporary(sinA), Temporary(cosA) ) * rad2deg )
		
		;Create a new variable
		oTraj_180 = MrTimeSeries( oTraj2_180['TIMEVAR'], traj, /CACHE, NAME=t_180_vname, /NO_COPY )
		
		;Attributes
		oTraj2_180 -> CopyAttrTo, oTraj_180
		oTraj_180['CATDESC'] = 'Trajectory of 180-degree pitch angle electrons in GSE coordinates. ' + $
		                       'This is an average of the trajectories from channels 2 and 3.'
		oTraj_180 -> RemoveAttr, 'CDF_NAME'
	ENDIF
END


;+
;   Update variables for AMB-PM2 data.
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' 
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'ql' | 'l2'}
;       COORDS:     in, optional, type=string, default='gse'
;                   The coordinate system in which trajectories are to be loaded. Options
;                       are: {'gse' | 'gsm'}
;
; :Keywords:
;       SUFFIX:     in, optional, type=string, default=''
;                   A suffix to be appendend to the end of EDI variable names.
;-
PRO MrMMS_EDI_Load_Data_AmbPM2, sc, mode, level, coords, $
SUFFIX=suffix
	Compile_Opt idl2
	On_Error, 2
	
	;EDI
	instr  = 'edi'
	
	;Variable names
	channels      = ['1', '2', '3', '4']
	flux0_vname   = sc + '_' + instr + '_' + 'flux' + channels + '_'                +   '0' + '_' + mode + '_' + level
	flux180_vname = sc + '_' + instr + '_' + 'flux' + channels + '_'                + '180' + '_' + mode + '_' + level
	traj0_vname   = sc + '_' + instr + '_' + 'traj' + channels + '_' + coords + '_' +   '0' + '_' + mode + '_' + level
	traj180_vname = sc + '_' + instr + '_' + 'traj' + channels + '_' + coords + '_' + '180' + '_' + mode + '_' + level
	
	;Output names
	f0_vname   = StrJoin( [sc, instr, 'flux',   '0',         mode, level], '_' ) + suffix
	f180_vname = StrJoin( [sc, instr, 'flux', '180',         mode, level], '_' ) + suffix
	t0_vname   = StrJoin( [sc, instr, 'traj',   '0', coords, mode, level], '_' ) + suffix
	t180_vname = StrJoin( [sc, instr, 'traj', '180', coords, mode, level], '_' ) + suffix
	
;-------------------------------------------
; 0-Degree Flux ////////////////////////////
;-------------------------------------------
	
	;Center channel
	;   - B-field is centered on channel 1
	IF Array_Equal(MrVar_IsCached(flux0_vname[0]), 1) THEN BEGIN
		;Get and copy the data
		oFlux1_0 = MrVar_Get(flux0_vname[0])
		oFlux_0  = oFlux1_0 -> Copy(f0_vname, /CACHE)
	ENDIF
	
;-------------------------------------------
; 0-Degree Trajectories ////////////////////
;-------------------------------------------
	
	;Center channel
	;   - B-field is centered on channel 1
	IF Array_Equal(MrVar_IsCached(traj0_vname[0]), 1) THEN BEGIN
		;Get and copy the data
		oTraj1_0 = MrVar_Get(traj0_vname[0])
		oTraj_0  = oTraj1_0 -> Copy(t0_vname, /CACHE)
	ENDIF
	
;-------------------------------------------
; 180-Degree Flux //////////////////////////
;-------------------------------------------
	
	;Center channel
	;   - B-field is centered on channel 1
	IF Array_Equal(MrVar_IsCached(flux180_vname[0]), 1) THEN BEGIN
		;Get and copy the data
		oFlux1_180 = MrVar_Get(flux180_vname[0])
		oFlux_180 = oFlux1_180 -> Copy(f180_vname, /CACHE)
	ENDIF
	
;-------------------------------------------
; 180-Degree Trajectories //////////////////
;-------------------------------------------
	
	;Center channel
	;   - B-field is centered on channel 1
	IF Array_Equal(MrVar_IsCached(traj180_vname[0]), 1) THEN BEGIN
		;Get and copy the data
		oTraj1_180 = MrVar_Get(traj180_vname[0])
		oTraj_180  = oTraj1_180 -> Copy(t180_vname, /CACHE)
	ENDIF
END


;+
;   Update variables for AMB-PM2 data.
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' 
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'ql' | 'l2'}
;       COORDS:     in, optional, type=string, default='gse'
;                   The coordinate system in which trajectories are to be loaded. Options
;                       are: {'gse' | 'gsm'}
;
; :Keywords:
;       SUFFIX:     in, optional, type=string, default=''
;                   A suffix to be appendend to the end of EDI variable names.
;-
PRO MrMMS_EDI_Load_Data_Q0, sc, mode, level, coords, $
SUFFIX=suffix
	Compile_Opt idl2
	On_Error, 2
	
	;EDI
	instr  = 'edi'
	
	cts_gdu1_vname       = StrJoin( [sc, instr, 'counts',         'gdu1', mode, level], '_' ) + suffix
	cts_gdu2_vname       = StrJoin( [sc, instr, 'counts',         'gdu2', mode, level], '_' ) + suffix
	traj_bcs_gdu1_vname  = StrJoin( [sc, instr, 'traj',   'bcs',  'gdu1', mode, level], '_' ) + suffix
	traj_bcs_gdu2_vname  = StrJoin( [sc, instr, 'traj',   'bcs',  'gdu2', mode, level], '_' ) + suffix
	traj_dbcs_gdu1_vname = StrJoin( [sc, instr, 'traj',   'dbcs', 'gdu1', mode, level], '_' ) + suffix
	traj_dbcs_gdu2_vname = StrJoin( [sc, instr, 'traj',   'dbcs', 'gdu2', mode, level], '_' ) + suffix
	traj_gse_gdu1_vname  = StrJoin( [sc, instr, 'traj',   'gse',  'gdu1', mode, level], '_' ) + suffix
	traj_gse_gdu2_vname  = StrJoin( [sc, instr, 'traj',   'gse',  'gdu2', mode, level], '_' ) + suffix
	abscal_gdu1_vname    = StrJoin( [sc, instr, 'abs',            'gdu1', 'cal', 'l2'], '_' ) + suffix
	abscal_gdu2_vname    = StrJoin( [sc, instr, 'abs',            'gdu2', 'cal', 'l2'], '_' ) + suffix
	
	;Output names
	f_gdu1_vname     = StrJoin( [sc, instr, 'flux',           'gdu1', 'q0', mode, level], '_' ) + suffix
	f_gdu2_vname     = StrJoin( [sc, instr, 'flux',           'gdu2', 'q0', mode, level], '_' ) + suffix
	traj1_bcs_vname  = StrJoin( [sc, instr, 'traj',   'bcs',  'gdu1', 'q0', mode, level], '_' ) + suffix
	traj2_bcs_vname  = StrJoin( [sc, instr, 'traj',   'bcs',  'gdu2', 'q0', mode, level], '_' ) + suffix
	traj1_dbcs_vname = StrJoin( [sc, instr, 'traj',   'dbcs', 'gdu1', 'q0', mode, level], '_' ) + suffix
	traj2_dbcs_vname = StrJoin( [sc, instr, 'traj',   'dbcs', 'gdu2', 'q0', mode, level], '_' ) + suffix
	traj1_gse_vname  = StrJoin( [sc, instr, 'traj',   'gse',  'gdu1', 'q0', mode, level], '_' ) + suffix
	traj2_gse_vname  = StrJoin( [sc, instr, 'traj',   'gse',  'gdu2', 'q0', mode, level], '_' ) + suffix
	ang_diff_vname   = StrJoin( [sc, instr, 'dlook',                  'q0', mode, level], '_' ) + suffix
	
;-------------------------------------------
; Load Calibrations ////////////////////////
;-------------------------------------------
	
	;EDI Calibrations
	MrMMS_EDI_Load_Cals, sc, $
	                     VARFORMAT='*abs*gdu*cal*'
	
;-------------------------------------------
; Apply Absolute Calibrations //////////////
;-------------------------------------------
	;GDU1
	IF MrVar_IsCached( cts_gdu1_vname ) THEN BEGIN
		;Grab the data
		oQ0_gdu1  = MrVar_Get(cts_gdu1_vname)
		oAbs_gdu1 = MrVar_Get(abscal_gdu1_vname)

		;Load calibrations
		ot_gdu1   = oQ0_gdu1['TIMEVAR']
		iAbs_gdu1 = oAbs_gdu1['TIMEVAR'] -> Value_Locate( ot_gdu1['DATA', 0, 'TT2000'], 'TT2000' )
	
		;Apply absolute calibrations
		oQ0_gdu1 = oQ0_gdu1 * oAbs_GDU1[iAbs_gdu1]
	
		;Name the variables
		oQ0_gdu1 -> SetName, f_gdu1_vname
	
		;Stash them away
		oQ0_gdu1 -> Cache
	ENDIF
	
	;GDU2
	IF MrVar_IsCached( cts_gdu2_vname ) THEN BEGIN
		;Grab the data
		oQ0_gdu2  = MrVar_Get(cts_gdu2_vname)
		oAbs_gdu2 = MrVar_Get(abscal_gdu2_vname)
	
		;Load calibrations
		ot_gdu2   = oQ0_gdu2['TIMEVAR']
		iAbs_gdu2 = oAbs_gdu2['TIMEVAR'] -> Value_Locate( ot_gdu2['DATA', 0, 'TT2000'], 'TT2000' )
	
		;Apply absolute calibrations
		oQ0_gdu2 = oQ0_gdu2 * oAbs_GDU2[iAbs_gdu2]
	
		;Name the variable
		oQ0_gdu2 -> SetName, f_gdu2_vname
	
		;Stash it away
		oQ0_gdu2 -> Cache
	ENDIF

;-------------------------------------------
; Fix DBCS Data ////////////////////////////
;-------------------------------------------
	;
	; In v2.1.z, the trajectories in DBCS do not have
	; any metadata. Current version is v4.0.0
	;

	;GDU1
	IF MrVar_IsCached(traj_dbcs_gdu1_vname) THEN BEGIN
		oTraj = MrVar_Get(traj_dbcs_gdu1_vname)
		
		;Metadata is missing if OTRAJ is not a MrTimeSeries
		;   - Cannot fix because we do not know bounds of MrVar_TLimit
		IF Obj_Class(oTraj) EQ 'MRVARIABLE' $
			THEN MrPrintF, 'LogWarn', 'Trajectories in DBCS are missing metadata.'
	ENDIF
	
	;GDU2
	IF MrVar_IsCached(traj_dbcs_gdu2_vname) THEN BEGIN
		oTraj = MrVar_Get(traj_dbcs_gdu2_vname)
		
		;Metadata is missing if OTRAJ is not a MrTimeSeries
		;   - Cannot fix because we do not know bounds of MrVar_TLimit
		IF Obj_Class(oTraj) EQ 'MRVARIABLE' $
			THEN MrPrintF, 'LogWarn', 'Trajectories in DBCS are missing metadata.'
	ENDIF

;-------------------------------------------
; Difference in Look Directions ////////////
;-------------------------------------------
	;Trajectories
	IF Array_Equal( MrVar_IsCached( [traj_gse_gdu1_vname, traj_gse_gdu2_vname] ), 1 ) THEN BEGIN
		oTraj_GDU1 = MrVar_Get(traj_gse_gdu1_vname)
		oTraj_GDU2 = MrVar_Get(traj_gse_gdu2_vname)
	ENDIF ELSE IF Array_Equal( MrVar_IsCached( [traj_dbcs_gdu1_vname, traj_dbcs_gdu2_vname] ), 1 ) THEN BEGIN
		oTraj_GDU1 = MrVar_Get(traj_dbcs_gdu1_vname)
		oTraj_GDU2 = MrVar_Get(traj_dbcs_gdu2_vname)
	ENDIF ELSE IF Array_Equal( MrVar_IsCached( [traj_bcs_gdu1_vname, traj_bcs_gdu2_vname] ), 1 ) THEN BEGIN
		oTraj_GDU1 = MrVar_Get(traj_bcs_gdu1_vname)
		oTraj_GDU2 = MrVar_Get(traj_bcs_gdu2_vname)
	ENDIF

	;DELTA LOOK
	IF Obj_Valid(oTraj_GDU1) && Obj_Valid(oTraj_GDU2) THEN BEGIN
		;Find all unique points
		t1   = oTraj_GDU1['TIME', 'TT2000']
		t2   = oTraj_GDU2['TIME', 'TT2000']
		tout = [ t1, t2 ]
		it   = Uniq(tout, Sort(tout))
		tout = tout[it]
		
		;Extract angles
		;   - Interpolate data onto unique points
		phi_gdu1   = Interpol( oTraj_GDU1['DATA', *, 0] * !dtor, t1, tout )
		theta_gdu1 = Interpol( oTraj_GDU1['DATA', *, 1] * !dtor, t1, tout )
		phi_gdu2   = Interpol( oTraj_GDU2['DATA', *, 0] * !dtor, t2, tout )
		theta_gdu2 = Interpol( oTraj_GDU2['DATA', *, 1] * !dtor, t2, tout )
	
		;Convert to cartesian coords
		x1 = Sin( theta_gdu1 ) * Cos( phi_gdu1 )
		y1 = Sin( theta_gdu1 ) * Sin( phi_gdu1 )
		z1 = Cos( theta_gdu1 )
		
		;GDU2
		x2 = Sin( theta_gdu2 ) * Cos( phi_gdu2 )
		y2 = Sin( theta_gdu2 ) * Sin( phi_gdu2 )
		z2 = Cos( theta_gdu2 )
		
		;Angular difference between GDU1 & GDU2
		angdiff = ACos( x1*x2 + y1*y2 + z1*z2 ) * !radeg
		oDAng   = MrScalarTS( tout, angdiff, $
		                      /CACHE, $
		                      NAME   = ang_diff_vname, $
		                      T_TYPE = 'TT2000' )
		
		t1 = !Null
		x1 = !Null
		y1 = !Null
		z1 = !Null
		t2 = !Null
		x2 = !Null
		y2 = !Null
		z2 = !Null
	ENDIF

;-------------------------------------------
; Rename Trajectories //////////////////////
;-------------------------------------------
	IF MrVar_IsCached(traj_bcs_gdu1_vname)  THEN MrVar_Rename, traj_bcs_gdu1_vname,  traj1_bcs_vname
	IF MrVar_IsCached(traj_bcs_gdu2_vname)  THEN MrVar_Rename, traj_bcs_gdu2_vname,  traj2_bcs_vname
	IF MrVar_IsCached(traj_dbcs_gdu1_vname) THEN MrVar_Rename, traj_dbcs_gdu1_vname, traj1_dbcs_vname
	IF MrVar_IsCached(traj_dbcs_gdu2_vname) THEN MrVar_Rename, traj_dbcs_gdu2_vname, traj2_dbcs_vname
	IF MrVar_IsCached(traj_gse_gdu1_vname)  THEN MrVar_Rename, traj_gse_gdu1_vname,  traj1_gse_vname
	IF MrVar_IsCached(traj_gse_gdu2_vname)  THEN MrVar_Rename, traj_gse_gdu2_vname,  traj2_gse_vname
	
;-------------------------------------------
; Clean Up /////////////////////////////////
;-------------------------------------------
	MrVar_Delete, [abscal_gdu1_vname, abscal_gdu2_vname]
END


;+
;   Load EDI data from source files.
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' }
;       OPTDESC:    in, out, optional, type=string/strarr
;                   Either the optional descriptor of EDI file names or a named variable
;                       in which all of the available optional descriptors is returned.
;                       Optional descriptors include: {'amb' | 'amb-pm2' | 'amb-alt-cc' |
;                       'amb-alt-oc' | 'amb-alt-ooc' | 'amb-alt-oob'}
;
; :Keywords:
;       COORDS:     in, optional, type=string, default='gse'
;                   The coordinate system in which trajectories are to be loaded. Options
;                       are: {'gse' | 'gsm'}
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'ql' | 'l2'}
;       SUFFIX:     in, optional, type=string, default=''
;                   A suffix to be appendend to the end of EDI variable names.
;       TRANGE:     in, optional, type=string/strarr(2), default=MrVar_GetTRange()
;                   The start and end times of the data interval to be plotted, formatted
;                       as 'YYYY-MM-DDThh:mm:ss'
;       VARFORMAT:  in, optional, type=string/strarr, default="['*flux*', '*counts*', '*traj*'+`COORDS`+'*']"
;                   A variable name filter.
;-
PRO MrMMS_EDI_Load_Data, sc, mode, optdesc, $
COORDS=coords, $
LEVEL=level, $
SUFFIX=suffix, $
TRANGE=trange, $
VARFORMAT=varformat
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN
	ENDIF
	
	;Defaults
	instr   = 'edi'
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(mode)      EQ 0 THEN mode      = 'srvy'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(suffix)    EQ 0 THEN suffix    = ''
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	IF N_Elements(varformat) EQ 0 THEN varformat = [ '*flux?_0_'+mode+'*', $
	                                                 '*flux?_180_'+mode+'*', $
	                                                 '*flux_gdu?_'+mode+'*', $
	                                                 '*counts*', $
	                                                 '*traj?_'+coords+'_0_'+mode+'*', $
	                                                 '*traj?_'+coords+'_180_'+mode+'*', $
	                                                 '*traj_'+coords+'_gdu?_'+mode+'*']
	
	;
	; TODO: If VARFORMAT is supplied by the user, COORDS can be broken.
	;

;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	;Determine which EDI file to load
	IF N_Elements(optdesc) EQ 0 THEN BEGIN
		fnames = MrMMS_Get_FileNames(sc, 'edi', mode, level, OPTDESC=optdesc)
		MrMMS_Parse_Filename, fnames, OPTDESC=optdesc
		iUniq = Uniq(optdesc, Sort(optdesc))
		IF N_Elements(iUniq) NE 1 THEN BEGIN
			MrPrintF, 'LogWarn', 'More than one EDI file type found.'
			MrPrintF, 'LogWarn', '   ' + '[' + StrJoin(optdesc, ', ') + ']'
		ENDIF
	ENDIF
	
	;EDI
	MrMMS_Load_Data, sc, instr, mode, level, $
	                 OPTDESC   = optdesc, $
	                 VARFORMAT = varformat

;-------------------------------------------
; Alter Data ///////////////////////////////
;-------------------------------------------
	IF mode EQ 'brst' THEN BEGIN
		FOR i = 0, N_Elements(optdesc) - 1 DO BEGIN
			CASE optdesc[i] OF
				'amb':     MrMMS_EDI_Load_Data_Amb,    sc, mode, level, coords, SUFFIX=suffix
				'amb-pm2': MrMMS_EDI_Load_Data_AmbPM2, sc, mode, level, coords, SUFFIX=suffix
				'q0':      MrMMS_EDI_Load_Data_Q0,     sc, mode, level, coords, SUFFIX=suffix
				ELSE: Message, 'Optdesc not implemented: "' + optdesc[i] + '".'
			ENDCASE
		ENDFOR
	ENDIF
END