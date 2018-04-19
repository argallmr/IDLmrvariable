; docformat = 'rst'
;
; NAME:
;       MrMMS_Load_Data
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
;       MrVariable, MMS, CDF
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
;       2016/09/02  -   Rename "Epoch" variables to follow the MMS variable naming
;                           convention. - MRA
;       2016/10/06  -   Added the SUFFIX keyword. Renamed Epoch variables are now
;                           reflected in VARNAMES. - MRA
;       2016/11/19  -   Load ancillary data products. - MRA
;       2017/01/24  -   Fix indexing error when multiple files are found for more than
;                           one file type. - MRA
;-
;*****************************************************************************************
;+
;   Provide unique variable names to CDF Epoch variables.
;
;   CDF data variables are saved to the cache with their original
;   names. Many have a DEPEND_0 variable that is support data, e.g.
;   Epoch. If the original data variable is removed from the cache,
;   there will no longer be a connection between the CDF and the
;   cache. The next time the CDF is read, Epoch will be read as a
;   new variable "Epoch_1".
;
;   To prevent duplicate data from being stored in the cache, we
;   must rename the Epoch variable before deleting data variables.
;
; :Private:
;
; :Params:
;       VARNAMES:       in, required, type=string/strarr
;                       Names of CDF variables for which the DEPEND_0 attribute
;                           is to be renamed.
;       SUFFIX:         in, optional, type=string, default=''
;                       A suffix to be appended to variable names.
;-
pro MrMMS_Load_Data_RenameEpoch, varnames, suffix
	compile_opt idl2
	on_error, 2
	
	;Number of segments in SUFFIX
	if n_elements(suffix) gt 0 then begin
		nSeg = n_elements( strsplit(suffix, '_') )
	endif else begin
		nSeg   = 0
		suffix = ''
	endelse
	
;-----------------------------------------------------
; Step Through Variables \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	for i = 0, n_elements(varnames) - 1 do begin
		;Get the variable
		;   - Epoch variables are being renamed, so it may no longer exist
		oVar = MrVar_Get(varnames[i], COUNT=count)
		if count eq 0 then continue
		
		;Skip if no DEPEND_0
		if ~oVar -> HasAttr('DEPEND_0') then continue
		
		;Get the name of the DEPEND_0 variable
		dep0      = oVar['DEPEND_0']
		dep0_type = size(dep0, /TNAME)
		if dep0_type eq 'STRING' $
			then old_name = dep0 $
			else old_name = dep0.name
		
		;Skip if DEPEND_0 is not an "Epoch  variable
		;   - DEPEND_0 can be a variable or variable name
		;   - Allow "_#" in case /NO_CLOBBER caused Epoch to be renamed.
		if ~stregex(old_name, '^Epoch', /BOOLEAN) then continue
	
	;-----------------------------------------------------
	; Create New Name \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Dissect the variable name
		parts = strsplit(varnames[i], '_', /EXTRACT)

		;Remove a suffix from OLD_NAME
		;   - TODO: If /NO_CLOBBER was set, "_#" might be appended to the
		;           end of the variable name. Change the regular expression
		;           to '^(.*)' + suffix + '(_[0-9]+)?' and follow through.
		epoch_base = stregex(old_name, '^(.*)' + suffix + '$', /SUBEXP, /EXTRACT)
		epoch_base = epoch_base[1]

		;Create a new name
		;   - Variable names are suppose to follow the format
		;       sc_instr_param[_coordsys][_optdesc]_mode_level[_suffix]
		;   - We want to make "Epoch" look like
		;       sc_instr_epoch_mode_level[_suffix]
		;   - But some instruments do not follow the veriable naming
		;     convention exactly.
		;   - Here, we follow the naming convention of the instrument
		case 1 of
			;DIS/DES
			stregex(varnames[i], '(des|dis)', /BOOLEAN): begin
				new_name = strjoin( [parts[0:1], strlowcase(epoch_base), parts[-1-nSeg:-1]], '_' )
			endcase
			
			;MEC
			stregex(varnames[i], 'mec', /BOOLEAN): begin
				if nSeg eq 0 $
					then new_name = strjoin( [parts[0:1], strlowcase(epoch_base)], '_' ) $
					else new_name = strjoin( [parts[0:1], strlowcase(epoch_base), parts[-nSeg:-1]], '_' )
			endcase
			
			;ANYTHING ELSE
			else: new_name = strjoin( [parts[0:1], strlowcase(epoch_base), parts[-2-nSeg:-1]], '_' )
		endcase
			
	;-----------------------------------------------------
	; Update DEPEND_0 Attribute \\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		oVar  -> SetAttrValue, 'DEPEND_0', new_name

	;-----------------------------------------------------
	; Rename Epoch Variable \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;DEPEND_0 is a name
		if dep0_type eq 'STRING' then begin
			
			;DEPEND_0's name needs to be updated if it is still in the cache
			;   - Do only the first time DEPEND_0 is being renamed.
			if MrVar_IsCached(old_name) then begin
				;Set the variable name
				oEpoch = MrVar_Get(old_name)
				oDep0 -> SetName, new_name
			
				;Update the variable name in VARNAMES
				idx = where(varnames eq old_name)
				varnames[idx] = new_name
			endif
		
		;DEPEND_0 is an object
		;   - Do only the first time DEPEND_0 is being renamed.
		endif else if dep0.name ne new_name then begin
			;Update DEPEND_0's name
			dep0 -> SetName, new_name
			
			;Update the variable name in VARNAMES
			idx           = where(varnames eq old_name)
			varnames[idx] = new_name
		endif
	endfor
end


;+
; :Params:
;       SC:                 in, required, type=string/strarr
;                           The MMS spacecraft identifier. Options are:
;                               {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       INSTR:              in, required, type=string/strarr
;                           Instrument ID for which data is read.
;       MODE:               in, required, type=string/strarr
;                           Telemetry mode of the data. Options include:
;                               {'slow' | 'fast' | 'srvy' | 'brst'}
;       LEVEL:              in, required, type=string/strarr
;                           Data quality level. Options include:
;                               {'l1a' | 'l1b' | 'l2pre' | 'l2'}
;
; :Keywords:
;       OPTDESC:            in, optional, type=string, default=''
;                           Optional descriptor of the data.
;       SUFFIX:             in, optional, type=string, default=''
;                           A suffix to be appended to variable names.
;       SUPPORT_DATA:       in, optional, type=boolean, default=0
;                           If set, support data will be read as well. Regardless of
;                               the status of this keyword, variable data associated
;                               with DEPEND_# variables will be read. This keyword
;                               is ignored if `VARFORMAT` is set.
;       TEAM_SITE:          in, optional, type=boolean, default=0
;                           If set, requests will be sent to the team site (password
;                               required, L1A and above). Automatically set if `LEVEL`
;                               is below level 2. This option is sticky.
;       TRANGE:             out, optional, type=string, default=MrVar_GetTRange()
;                           Start and end times over which to read data.
;       VARFORMAT:          in, optional, type=string, default=''
;                           Variables that match this search pattern will be read,
;                               others are ignored.
;       VARNAMES:           out, optional, type=string/strarr
;                           Names of the variables that have been loaded into the cache.
;       _REF_EXTRA:         in, optional, type=any
;                           Any keyword accepted by MrMMS_SDC_Query::CD is also accepted
;                               via keyword inheritance.
;-
pro MrMMS_Load_Data, sc, instr, mode, level, $
OPTDESC=optdesc, $
SUPPORT_DATA=support_data, $
SUFFIX=suffix, $
TEAM_SITE=team_site, $
TRANGE=trange, $
VARFORMAT=varformat, $
VARNAMES=varnames, $
_REF_EXTRA=extra
	compile_opt idl2

	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return
	endif
	
	;Defaults
	if n_elements(sc)     eq 0 then sc     = ''
	if n_elements(instr)  eq 0 then instr  = ''
	if n_elements(mode)   eq 0 then mode   = ''
	if n_elements(level)  eq 0 then level  = ''
	if n_elements(trange) eq 0 then trange = MrVar_GetTRange()
	
	;Do not let variable names get passed in
	if n_elements(varnames) gt 0 then varnames = !Null
	
	;Initialize MMS
	MrMMS_Init

;-----------------------------------------------------
; Web \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Change direstories to the proper request
	;   - Always start with a new request.
	!MrMMS -> SetProperty, /RESET_PATH, $
	                       /RESET_QUERY, $
	                       /DOWNLOAD, $
	                       SC_ID         = sc, $
	                       INSTR         = instr, $
	                       MODE          = mode, $
	                       LEVEL         = level, $
	                       OPTDESC       = optdesc, $
	                       PUBLIC_SITE   = public_site, $
	                       TEAM_SITE     = team_site, $
	                       DATE_START    = trange[0], $
	                       DATE_END      = trange[1], $
	                       _STRICT_EXTRA = extra
	
	;Attempt to get the data
	files = !MrMMS -> Get(COUNT=count)
	if count eq 0 then return

;-----------------------------------------------------
; Read Ancillary Files \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Find ancillary files
	iAnc = where( !MrMMS -> IsAncillary(files), nAnc, $
	              COMPLEMENT=iData, NCOMPLEMENT=nData )
	
	;Read ancillary files
	if nAnc gt 0 then begin
		MrMMS_Anc_Load, files[iAnc], $
		                VARFORMAT = varformat, $
		                VARNAMES  = varnames
	endif
	
;-----------------------------------------------------
; Read Data Files \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	if nData gt 0 then begin
		;Parse the file names
		files = files[iData]
		MrMMS_Parse_Filename, files, SC=fsc, INSTR=finstr, MODE=fmode, LEVEL=flevel, OPTDESC=foptdesc
		
		;Determine file types
		fType = fsc + '_' + finstr + '_' + fmode + '_' + flevel
		iOptDesc = where(foptdesc ne '', nOptDesc)
		if nOptDesc gt 0 then fType[iOptDesc] += '_' + foptdesc[iOptDesc]
		
		;Pick unique file types
		fbase = file_basename(files)
		iUniq = uniq(fType, sort(fType))
		nUniq = n_elements(iUniq)
	
		;Step through each file type
		for i = 0, nUniq - 1 do begin
			;Find similar file types
			iFiles = where(strmid( fbase, 0, strlen(fType[iUniq[i]])) eq fType[iUniq[i]], nFiles)
			if nFiles eq 0 then continue
			
			;
			; TODO: Sort files by time (and slow/fast srvy)
			;
			
			;Read files
			MrVar_ReadCDF, files[iFiles], $
			               SUFFIX       = suffix, $
			               SUPPORT_DATA = support_data, $
			               TRANGE       = trange, $
			               VARFORMAT    = varformat, $
			               VARNAMES     = temp_names, $
			               VERBOSE      = !MrMMS.verbose
			
			;Save all of the variables
			if n_elements(varnames) gt 0 $
				then varnames = [varnames, temporary(temp_names)] $
				else varnames = temporary(temp_names)
		endfor
	endif

;-----------------------------------------------------
; Rename the Epoch Variables \\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	MrMMS_Load_Data_RenameEpoch, varnames, suffix
end
