; docformat = 'rst'
;
; NAME:
;       MrMMS_Get_Data
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
;   Read CDF varialbes into MrVariables and cache the data.
;
; :Categories:
;       MMS, FGM
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
;       2016/08/16  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   FGM CDF files combine B and |B| into a 4D array. Split up the vector and magnitude.
;
; :Params:
;       B_NAME:         in, required, type=string
;                       Name of the magnetic field variable.
;-
pro MrMMS_FGM_Load_Data_SplitB, b_name, bvec_name, bmag_name
	compile_opt idl2
	on_error, 2

	;Find a place to split the variable name
	;   - Add "vec" and "mag"
	pos = strpos(b_name, '_b_')
	if pos eq -1 then begin
		pos       = stregex(b_name, '_(afg|dfg|fgm)_', LEN=len)
		bvec_name = strmid(b_name, 0, pos+len) + 'vec_' + strmid(b_name, pos+len)
		bmag_name = strmid(b_name, 0, pos+len) + 'mag_' + strmid(b_name, pos+len)
	endif else begin
		bvec_name = strmid(b_name, 0, pos+2) + 'vec' + strmid(b_name, pos+2)
		bmag_name = strmid(b_name, 0, pos+2) + 'mag' + strmid(b_name, pos+2)
	endelse

	;Grab the variable
	oB = MrVar_Get(b_name)
	oT = oB['TIMEVAR']

	;Separate magnitude from vector
	Bxyz  = MrVectorTS( oT, oB['DATA',*,0:2], NAME=bvec_name, /CACHE )
	Bmag  = MrScalarTS( oT, oB['DATA',*,3],   NAME=bmag_name, /CACHE )
	
	;Copy over all attributes
	oB -> CopyAttrTo, Bxyz
	oB -> CopyAttrTo, Bmag
	
	min_value = min(Bxyz['MIN_VALUE'])
	max_value = max(Bxyz['MAX_VALUE'])
	
	;B
	oB['COLOR']     = ['Blue', 'Forest Green', 'Red', 'Black']
	oB['LABEL']     = ['Bx', 'By', 'Bz', '|B|']
	oB['MIN_VALUE'] = min_value
	oB['MAX_VALUE'] = max_value
	oB['TITLE']     = 'B!C(nT)'

	;BVEC - Set new attributes
	Bxyz['COLOR']     = ['blue', 'forest green', 'red']
	Bxyz['DIMENSION'] = 1
	Bxyz['LABEL']     = ['Bx', 'By', 'Bz']
	Bxyz['MIN_VALUE'] = min_value
	Bxyz['MAX_VALUE'] = max_value
	Bxyz['TITLE']     = 'B!C(nT)'
	
	;BMAG - Set new attributes
	Bmag['COLOR']     = 'black'
	Bmag['MIN_VALUE'] = 0
	Bmag['MAX_VALUE'] = max(Bmag['MAX_VALUE'])
	Bmag['TITLE']     = '|B|!C(nT)'
	Bmag             -> RemoveAttr, 'LABEL'
end


;+
;   Find and read MMS FGM data.
;
;   The following quantities are split. When using VARFORMAT, the original
;   names must be provided.
;       b_bcs    ->    bvec_bcs  & bmag_bcs
;       b_dmpa   ->    bvec_dmap & bmag_dmpa
;       b_gse    ->    bvec_gse  & bmag_gse
;       b_gsm    ->    bvec_gsm  & bmag_gsm
;
;
; :Params:
;       SC:                 in, required, type=string/strarr
;                           The MMS spacecraft identifier. Options are:
;                               {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:               in, required, type=string/strarr
;                           Telemetry mode of the data. Options include:
;                               {'slow' | 'fast' | 'srvy' | 'brst'}
;
; :Keywords:
;       INSTR:              in, optional, type=string/strarr, default='dfg'
;                           Instrument ID for which data is read. Options are 
;                               {'fgm', 'dfg', 'afg'}
;       LEVEL:              in, optional, type=string/strarr, default='l2'
;                           Data quality level. Options include:
;                               {'l1a' | 'l1b' | 'l2pre' | 'l2'}
;       OPTDESC:            in, optional, type=string, default=''
;                           Optional descriptor of the data.
;       TEAM_SITE:          in, optional, type=boolean, default=0
;                           If set, requests will be sent to the team site (password
;                               required, L1A and above). Automatically set if `LEVEL`
;                               is below level 2. This option is sticky.
;       SUFFIX:             in, optional, type=string, default=''
;                           A suffix to be appended to variable names.
;       TRANGE:             out, optional, type=string, default=MrVar_GetTRange()
;                           Start and end times over which to read data.
;       VARFORMAT:          out, optional, type=string, default='*'
;                           Variables that match this search pattern will be read,
;                               others are ignored.
;-
pro MrMMS_FGM_Load_Data, sc, mode, $
INSTR=instr, $
LEVEL=level, $
OPTDESC=optdesc, $
SUFFIX=suffix, $
TEAM_SITE=team_site, $
TRANGE=trange, $
VARFORMAT=varformat, $
VARNAMES=varnames
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return
	endif
	
	;Defaults
	if n_elements(suffix) eq 0 then suffix = ''
	if n_elements(instr)  eq 0 then instr  = 'fgm'
	if n_elements(level)  eq 0 then begin
		case instr of
			'afg': level = 'l2pre'
			'dfg': level = 'l2pre'
			'fgm': level = 'l2'
			else:  message, 'Invalid FGM instrument: "' + instr + '".'
		endcase
	endif
	
	;Check spacecraft
	tf_instr = MrIsMember(['afg', 'dfg', 'fgm'], instr)
	if ~array_equal(tf_instr, 1) then message, 'Invalid value for INSTR: "' + instr + '".'
	
	;DFG || AFG
	if (instr eq 'dfg' || instr eq 'afg') && level eq 'l2' then begin
		MrPrintF, 'LogWarn', 'L2 data is not avalable for ' + instr + '. Switching to L2Pre.'
		fgm_level = 'l2pre'
	
	;FGM
	endif else if instr eq 'fgm' && level ne 'l2' then begin
		MrPrintF, 'LogWarn', StrUpCase(level) + ' data is not available for FGM. Switching to L2.'
		fgm_level = 'l2'
	endif else begin
		fgm_level = level
	endelse

	;Get the data
	MrMMS_Load_Data, sc, instr, mode, fgm_level, $
	                 OPTDESC   = optdesc, $
	                 SUFFIX    = suffix, $
	                 TEAM_SITE = team_site, $
	                 TRANGE    = trange, $
	                 VARFORMAT = varformat, $
	                 VARNAMES  = varnames
	
;-----------------------------------------------------
; Prettify Things \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Build variable names
	nVars = n_elements(varnames)
	for i = 0, nVars - 1 do begin
		;Separate B from |B|
		;   - L2:   mms#_afg_b_gse_srvy_l2
		;   - QL:   mms#_afg_srvy_dmpa has different names from L1A, L1B, L2Pre
		;   - ELSE: mms#_afg_srvy_l1a_dmpa
		if stregex(varnames[i], '(afg|dfg|fgm)_b_(gse|gsm|dmpa|bcs)_', /BOOLEAN) || $
		   stregex(varnames[i], '(afg|dfg|fgm)_(slow|fast|srvy|brst)_(l1a|l1b|l2pre_)?(dmpa|gse|gsm|gsm|gse)_.*' + suffix + '$', /BOOLEAN) $
		then begin
			MrMMS_FGM_Load_Data_SplitB, varnames[i], bvec_name, bmag_name
			varnames = [varnames, bvec_name, bmag_name]
		endif
	endfor
end
