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
FUNCTION MrMMS_EDI_Plot_PitchAngles, sc, mode, $
COORDS=coords, $
NO_LOAD=no_load, $
TRANGE=trange
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		if n_elements(win) gt 0 then obj_destroy, win
		MrPrintF, 'LogErr'
		return, obj_new()
	endif
	
	;Defaults
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(coords) EQ 0 THEN coords = 'gse'
	IF N_Elements(trange) GT 0 THEN MrVar_SetTRange, trange
	
	;Variable type
	level = 'l2'
	type  = level EQ 'l2' ? 'flux' : 'counts'
	
;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	;Parts
	edi_ch    = mode EQ 'brst' ? ['1', '2', '3', '4'] : '1'
	edi_instr = 'edi'
	edi_param = level EQ 'l2' ? 'flux' : 'counts'
	
	;Find which optional descriptor is available
	fnames = MrMMS_Get_Filenames(sc, edi_instr, mode, level)
	MrMMS_Parse_Filename, fnames, OPTDESC=optdesc, VERSION=version
	iUniq = Uniq(optdesc, Sort(optdesc))
	IF N_Elements(iUniq) NE 1 THEN BEGIN
		MrPrintF, 'LogWarn', 'More than one EDI file type found.'
		MrPrintF, 'LogWarn', '   ' + '[' + StrJoin(optdesc[iUniq], ', ') + ']'
		MrPrintF, 'LogWarn', '   Choosing "' + optdesc[0] + '".'
	ENDIF
	edi_optdesc = optdesc[0]
	
	CASE edi_optdesc OF
		'amb': 
		'amb-pm2': 
		'amb-alt-cc':
		'amb-alt-oc':
		'amb-alt-ooc':
		'amb-alt-oob':
		ELSE: Message, 'Invalid optional descriptor: "' + optdesc + '".'
	ENDCASE
	
	fgm_instr  = 'dfg'
	fgm_level  = 'l2pre'
	fgm_coords = coords EQ 'dbcs' ? 'dmpa' : coords
	fgm_mode   = mode EQ 'brst' ? mode : 'srvy'
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	
	;EDI
	IF StrMatch(edi_optdesc, 'amb-alt-oob') THEN BEGIN
		edi_flux_0       = StrJoin([sc, edi_instr, edi_param], '_') + edi_ch + '_' + StrJoin([  '0',         mode, level], '_')
		edi_flux_90_gdu1 = StrJoin([sc, edi_instr, edi_param], '_') + edi_ch + '_' + StrJoin([ '90', 'gdu1', mode, level], '_')
		edi_flux_90_gdu2 = StrJoin([sc, edi_instr, edi_param], '_') + edi_ch + '_' + StrJoin([ '90', 'gdu2', mode, level], '_')
		edi_flux_180     = StrJoin([sc, edi_instr, edi_param], '_') + edi_ch + '_' + StrJoin(['180',         mode, level], '_')
		
		edi_traj_0       = StrJoin([sc, edi_instr, 'traj'], '_') + edi_ch + '_' + StrJoin([coords,  '0',         mode, level], '_')
		edi_traj_90_gdu1 = StrJoin([sc, edi_instr, 'traj'], '_') + edi_ch + '_' + StrJoin([coords, '90', 'gdu1', mode, level], '_')
		edi_traj_90_gdu2 = StrJoin([sc, edi_instr, 'traj'], '_') + edi_ch + '_' + StrJoin([coords, '90', 'gdu2', mode, level], '_')
		edi_traj_180     = StrJoin([sc, edi_instr, 'traj'], '_') + edi_ch + '_' + StrJoin([coords,'180',         mode, level], '_')
		traj_vnames      = [edi_traj_0, edi_traj_90_gdu1, edi_traj_90_gdu2, edi_traj_180]
		
		;Derived names
		edi_fac_0       = StrJoin([sc, edi_instr, 'traj'], '_') + edi_ch + '_' + StrJoin(['fac',  '0',         mode, level], '_')
		edi_fac_90_gdu1 = StrJoin([sc, edi_instr, 'traj'], '_') + edi_ch + '_' + StrJoin(['fac', '90', 'gdu1', mode, level], '_')
		edi_fac_90_gdu2 = StrJoin([sc, edi_instr, 'traj'], '_') + edi_ch + '_' + StrJoin(['fac', '90', 'gdu2', mode, level], '_')
		edi_fac_180     = StrJoin([sc, edi_instr, 'traj'], '_') + edi_ch + '_' + StrJoin(['fac','180',         mode, level], '_')
		fac_vnames      = [edi_fac_0, edi_fac_90_gdu1, edi_fac_90_gdu2, edi_fac_180]
	ENDIF ELSE BEGIN
		Message, 'EDI optional descriptor not recognized: "' + edi_optdesc + '".'
	ENDELSE
	
	;FGM
	b_vname     = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_' )
	bvec_vname  = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_' )
	bmag_vname  = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_' )

;-------------------------------------------
; EDI Data /////////////////////////////////
;-------------------------------------------
	;Get EDI data
	IF tf_load THEN BEGIN
		;Load EDI data
		MrMMS_Load_Data, sc, edi_instr, mode, level, $
		                 OPTDESC   = edi_optdesc, $
		                 VARFORMAT = ['*' + edi_param + '?_0_'       + mode + '*', $
		                              '*' + edi_param + '?_90_gdu1_' + mode + '*', $
		                              '*' + edi_param + '?_90_gdu2_' + mode + '*', $
		                              '*' + edi_param + '?_180_'     + mode + '*', $
		                              '*traj?_' + coords + '_0_'       + mode + '*', $
		                              '*traj?_' + coords + '_90_gdu1_' + mode + '*', $
		                              '*traj?_' + coords + '_90_gdu2_' + mode + '*', $
		                              '*traj?_' + coords + '_180_'     + mode + '*']
		
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = fgm_level, $
		                     VARFORMAT = b_vname
		
		IF fgm_level NE 'l2' && ~MrVar_IsCached(b_vname) THEN BEGIN
			b_vname     = StrJoin( [sc, fgm_instr,        fgm_mode, fgm_level, fgm_coords], '_' )
			bvec_vname  = StrJoin( [sc, fgm_instr, 'vec', fgm_mode, fgm_level, fgm_coords], '_' )
			bmag_vname  = StrJoin( [sc, fgm_instr, 'mag', fgm_mode, fgm_level, fgm_coords], '_' )
			
			MrMMS_FGM_Load_Data, sc, fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = fgm_level, $
			                     VARFORMAT = b_vname
			
			IF ~MrVar_IsCached(fgm_b_vname) $
				THEN Message, 'Unknown FGM variable naming convention.'
		ENDIF
	ENDIF
	
;-------------------------------------------
; Rotate to FAC ////////////////////////////
;-------------------------------------------
	oB = MrVar_Get(bvec_vname)
	oB_hat = oB -> Normalize()
	FOR i = 0, N_Elements(traj_vnames) - 1 DO BEGIN
		oTraj = MrVar_Get(traj_vnames[i])
		oT    = MrVar_FAC(oB, '', 'CROSSX', TIME=oTraj)
		oVec  = MrMMS_EDI_sphr2cart(oTraj)
		oVec  = oT ## oVec
		
;		theta = oVec -> Normalize()
;;		theta = theta -> Dot(oB_hat)
;		theta = acos(theta['DATA'])
;		oTheta = MrScalarTS(oVec['TIMEVAR']

		oTraj = MrMMS_EDI_cart2sphr(oVec, /CACHE, NAME=fac_vnames[i])
		oTraj['TITLE']         = '!9' + String(97B) + '!X$\downB$!C(deg)'
		oTraj['UNITS']         = 'deg'
		oTraj['YTICKINTERVAL'] = 90.0
		oTraj['YRANGE']        = [-180.0, 180.0]
		Obj_Destroy, oVec
	ENDFOR
	
;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	oCh1 = MrVar_Get(fac_vnames[0])
	oCh2 = MrVar_Get(fac_vnames[1])
	oCh3 = MrVar_Get(fac_vnames[2])
	oCh4 = MrVar_Get(fac_vnames[3])
	oCh1['PLOT_TITLE'] = 'Channel 1'
	oCh2['PLOT_TITLE'] = 'Channel 2'
	oCh3['PLOT_TITLE'] = 'Channel 3'
	oCh4['PLOT_TITLE'] = 'Channel 4'

;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	;Plot burst data
	win = MrVar_PlotTS( fac_vnames, LAYOUT=[4,4], XSIZE=1000, YSIZE=600, /NO_REFRESH )

	win[0].SetLayout, [1,1]
	win.TrimLayout
	win -> Refresh
	RETURN, win
END