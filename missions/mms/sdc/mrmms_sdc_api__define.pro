; docformat = 'rst'
;
; NAME:
;       MrMMS_SDC_API__DEFINE
;
;*****************************************************************************************
;   Copyright (c) 2016, University of New Hampshire                                      ;
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
;       * Neither the name of the University of New Hampshire nor the names of its       ;
;         contributors may be used to endorse or promote products derived from this      ;
;         software without specific prior written permission.                            ;
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
;   An interface into the MMS Science Data Center
;
;   Class Structure:
;       MrMMS_SDC_API
;           MrWebURI
;               IDLnetURL
;               MrURI
;                   IDL_Object
;
;   OFFLINE MODE:
;      ::SetURI  -- Parses the URI and sets object properties. Does not check if the URI
;                   is valid. This is because ::CD will pass it a URI with a query string,
;                   which means there is not a 1-to-1 correspondence between the remote
;                   URI and the local URI.
;      ::Get     -- Will call ::GetLocalFiles instead of ::GetFileNames
;      ::uri2dir -- Since ::GetLocalFiles already returns local file paths, this method
;                   will return the URI path instead of trying to convert a remote path to
;                   a local path.
;
;   NO_DOWNLOAD:
;      ::Get     -- Search for files remotely, then find their local counter-parts. If files
;                   are not found locally, they are ignored.
;
; :Categories:
;       MMS
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
;       2016/06/30  -   Written by Matthew Argall
;       2016/08/30  -   Bump "End_Date" field up to next day if time is not 00:00:00. - MRA
;       2016/11/19  -   Incorporate ancillary data products. - MRA
;       2018/01/25  -   In ::SetProperty, RESET_PATH was overwriting input arguments. - MRA
;       2018/03/08  -   When picking auto defaults for TEAM, PUBLIC in ::SetProperty, do
;                           the same for SITE. - MRA
;-
;*****************************************************************************************
;+
;   The purpose of this method is to provide information when implied print is used.
;   is called.
;-
FUNCTION MrMMS_SDC_API::_OverloadPrint
	Compile_Opt idl2
	On_Error, 2
	
	;Get info from the superclass
	web_print = self.oWeb -> _OverloadPrint()

	;Create strings
	site      = String('  Site',      '=', self.site,      FORMAT='(a-26, a-2, a0)')
	info_type = String('  Info_Type', '=', self.info_type, FORMAT='(a-26, a-2, a0)')
	data_type = String('  Data_Type', '=', self.data_type, FORMAT='(a-26, a-2, a0)')
	sc        = String('  SC',        '=', self.sc,        FORMAT='(a-26, a-2, a0)')
	instr     = String('  Instr',     '=', self.instr,     FORMAT='(a-26, a-2, a0)')
	mode      = String('  Mode',      '=', self.mode,      FORMAT='(a-26, a-2, a0)')
	level     = String('  Level',     '=', self.level,     FORMAT='(a-26, a-2, a0)')
	optdesc   = String('  OptDesc',   '=', self.optdesc,   FORMAT='(a-26, a-2, a0)')
	version   = String('  Version',   '=', self.version,   FORMAT='(a-26, a-2, a0)')
	ancprod   = String('  AncProd',   '=', self.ancprod,   FORMAT='(a-26, a-2, a0)')
	files     = String('  Files',     '=', self.files,     FORMAT='(a-26, a-2, a0)')

	;Output array
	outStr = [ [ web_print ], $
	           [ site      ], $
	           [ info_type ], $
	           [ data_type ], $
	           [ sc        ], $
	           [ instr     ], $
	           [ mode      ], $
	           [ level     ], $
	           [ optdesc   ], $
	           [ version   ], $
	           [ ancprod   ], $
	           [ files     ] $
	         ]
	
	;Sort alphabetically
	outStr = outStr[0,Sort(outStr)]

	;Print the array
	RETURN, outStr
END


;+
;   Filter files by time. Burst and Srvy files should be filtered separately.
;
; :Private:
;
; :Params:
;       FILENAMES:      in, required, type=string/strarr
;                       File names to be filtered.
;       TSTART:         in, required, type=string
;                       Start time of filter.
;       TEND:           in, required, type=string
;                       End time of filter.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of files that pass the filter.
;
; :Returns:
;       FILES_OUT:      Files that pass the filter.
;-
function MrMMS_SDC_API::Anc_FilterTime, filenames, $
COUNT=count
	compile_opt idl2
	on_error, 2
	
	;Get current file parameters
	self.oWeb -> GetProperty, FPATTERN  = fPattern, $
	                          TPATTERN  = tpattern

	;Set ancillary file parameters
	self.oWeb -> SetProperty, FPATTERN  = 'MMS*_%Y%D_%Y%D.V*'
	
	;Filter using superclass
	files_out = self.oWeb -> FilterTime( filenames, $
	                                     COUNT = count )
	
	;Return to normal
	self.oWeb -> SetProperty, FPATTERN  = fPattern, $
	                          TPATTERN  = tpattern
	
	return, files_out
end


;+
;   Filter ancillary files by version number.
;
; :Private:
;
; :Params:
;       FILENAMES:      in, required, type=string/strarr
;                       File names to be filtered.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of files that pass the filter.
;       LATEST_VERSION: in, optional, type=boolean, default=0
;                       If set, only the newest version of the files is returned. This
;                           is the default if `MIN_VERSION` and `VERSION` are not set.
;       MIN_VERSION:    in, optional, type=string, default=''
;                       A version number formatted as 'X.Y.Z'. Files with
;                           versions smaller this version will be filtered out.
;       VERSION:        in, optional, type=string, default=''
;                       A version number formatted as 'X.Y.Z'. Only files of
;                           this version will be returned.
;
; :Returns:
;       FILES_OUT:      Files that pass the filter.
;-
function MrMMS_SDC_API::Anc_FilterVersion, filenames, $
COUNT=count, $
LATEST_VERSION=latest_version, $
MIN_VERSION=min_version, $
VERSION=version
	compile_opt idl2
	on_error, 2
	
	;Defaults
	tf_latest = keyword_set(latest_version)
	tf_checkv = n_elements(version)     gt 0
	tf_minv   = n_elements(min_version) gt 0
	if tf_checkv + tf_minv + tf_latest  eq 0 then tf_latest = 1B
	if tf_checkv + tf_minv + tf_latest  gt 1 $
		then message, 'VERSION, MIN_VERSION and LATEST_VERSION are mutually exclusive.'

	;Results
	files_out = filenames
	count     = n_elements(filenames)

	;Extract version number from file
	fversion = stregex(files_out, 'V([0-9]+)$', /SUBEXP, /EXTRACT)
	fv = fix( reform(fversion[1,*]) )

;------------------------------------;
; Minimum Version                    ;
;------------------------------------;
	if tf_minv then begin
		;Select file versions
		iv = where( fv ge fix(min_version), count )
		if count gt 0 then files_out = files_out[iv]

;------------------------------------;
; Exact Version                      ;
;------------------------------------;
	endif else if tf_checkv then begin
		;Select secific versions
		iv = where( fv eq fix(version), count )
		if count gt 0 then files_out = files_out[iv]

;------------------------------------;
; Latest Version                     ;
;------------------------------------;
	endif else begin
		;Select latest version
		iv = where( fv eq max(fv), count )
		if count gt 0 then files_out = files_out[iv]
	endelse
	
	;Return
	if count eq 0 then files_out = ''
	if count eq 1 then files_out = files_out[0]
	return, files_out
end


;+
;   Build an MMS URI using properties of the API provided by the science data center.
;
; :Keywords:
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrMMS_SDC_Build_URI is accepted here.
;                           By default, object property values are used.
;-
FUNCTION MrMMS_SDC_API::BuildURI, $
_REF_EXTRA=extra
	Compile_Opt idl2
	On_Error, 2
	
	;Use the independent function.
	;   - Keywords passed in explicitly are over-ridden by those in EXTRA
	uri = MrMMS_SDC_Build_URI( ANC_PRODUCT   = self.ancprod, $
	                           ANCILLARY     = self.ancillary, $
	                           DOWNLOAD      = self.download, $
	                           HK            = self.hk, $
	                           FILES         = self.files, $
	                           FILE_INFO     = self.info, $
	                           FILE_NAMES    = self.names, $
	                           INSTR         = self.instr, $
	                           LEVEL         = self.level, $
	                           MODE          = self.mode, $
	                           OPTDESC       = self.optdesc, $
	                           PUBLIC_SITE   = self.public, $
	                           SCIENCE       = self.science, $
	                           SC_ID         = self.sc_id, $
	                           TEAM_SITE     = self.team, $
	                           TEND          = self.tend, $
	                           TSTART        = self.tstart, $
	                           V_INFO        = self.v_info, $
	                           VERSION       = self.version, $
	                           _STRICT_EXTRA = extra )
	
	RETURN, uri
END


;+
;   Build a URI
;
; :Keywords:
;       DOWNLOAD:       in, optional, type=boolean, default=0
;                       If set, the path will be returned to the default path. This is
;                           equivalent to setting `PUBLIC`, `DOWNLOAD`, and `SCIENCE`.
;-
function MrMMS_SDC_API::BuildLocalName, $
DROPBOX=dropbox
	compile_opt idl2
	on_error, 2
	
	;Defaults
	tf_dropbox = keyword_set(dropbox)
	
	;Define the directory structure to use
	if tf_dropbox $
		then dropbox_root = self.dropbox_root $
		else sdc_root     = self.local_root

;-----------------------------------------------------
; Use File Names \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	if (*self.field_values).file[0] ne '' then begin
		;Separate ancillary files from the rest
		ftemp = strsplit((*self.field_values).file, ',', /EXTRACT)
		iAnc  = where( self -> IsAncillary(ftemp), nAnc, $
		               COMPLEMENT=iData, NCOMPLEMENT=nData)
		
		;DROPBOX is a flat directory structure
		if tf_dropbox then begin
			uri = filepath(temporary(ftemp), ROOT_DIR=dropbox_root)
		
		;SDC directory structure
		endif else begin
			uri_root = strarr(nAnc + nData)
			if nAnc  gt 0 then uri_root[iAnc]  = MrMMS_Anc_Build_Path( ftemp[iAnc], SDC_ROOT=sdc_root )
			if nData gt 0 then uri_root[iData] = MrMMS_Build_Path( ftemp[iData], SDC_ROOT=sdc_root )
			uri = temporary(uri_root) + temporary(ftemp)
		endelse
		
;-----------------------------------------------------
; Use Query Field Values \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Directories from field-value pairs
	;   - TSTART (& TEND) will be made from MrTokens
	;   - VERSION will be the asterisk
	endif else begin
		;Ancillary
		if (*self.field_values).product ne '' then begin
			uri = MrMMS_Anc_Build_Filename( strsplit( (*self.field_values).sc_id,   ',', /EXTRACT ), $
			                                strsplit( (*self.field_values).product, ',', /EXTRACT ), $
			                                DIRECTORY = dropbox_root, $
			                                SDC_ROOT  = sdc_root )
		
		;Data
		endif else begin
			uri = MrMMS_Build_Filename( strsplit( (*self.field_values).sc_id,          ',', /EXTRACT ), $
			                            strsplit( (*self.field_values).instrument_id,  ',', /EXTRACT ), $
			                            strsplit( (*self.field_values).data_rate_mode, ',', /EXTRACT ), $
			                            strsplit( (*self.field_values).data_level,     ',', /EXTRACT ), $
			                            DIRECTORY = dropbox_root, $
			                            OPTDESC   = strsplit( (*self.field_values).descriptor, ',', /EXTRACT ), $
			                            SDC_ROOT  = sdc_root )
		endelse
	endelse
	
	return, uri
end


;+
;   Download a single file from the SDC.
;
; :Params:
;       FILENAME:       in, optional, type=string
;                       Name given to the file being downloaded. By default, the
;                           server-provided name will be used.
;
; :Returns:
;       FILEOUT:         Names of the downloaded files.
;-
FUNCTION MrMMS_SDC_API::Download, filename, $
STRING=str
	Compile_Opt idl2

	;Catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF n_elements(fileOut) GT 0 && file_test(fileOut) THEN file_delete, fileOut
		MrPrintF, 'LogErr'
		RETURN, ''
	ENDIF

;---------------------------------------------------------------------
; Download ///////////////////////////////////////////////////////////
;---------------------------------------------------------------------

	;Get data.
	;   - Cookies (credentials) are stored in HEADER, which is reset by GET.
	;   - 
	fileOut = self.oWeb -> IDLnetURL::Get( FILENAME = filename, $
	                                       STRING   = str )
	
;---------------------------------------------------------------------
; Clean Up ///////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	
	;Return the file
	RETURN, fileOut
END


;+
;   Get information about files.
;
; :Params:
;       FILENAMES:      out, optional, type=strarr
;                       Names of the files that match search criteria.
;       FILESIZE:       out, optional, type=ulonarr
;                       Names of the files that match search criteria.
;       MODIFIED:       out, optional, type=strarr
;                       Dates files were modified, formatted as 'YYYY-MM-DDThh:mm:ss'.
;       START_DATE:     out, optional, type=strarr
;                       Time tag of files, formatted as 'YYYY-MM-DDThh:mm:ss'.
;       END_DATE:       out, optional, type=strarr
;                       End time tag of files, formatted as 'YYYY-MM-DDThh:mm:ss'. Values
;                           are equal to the empty string for non-ancillary files.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of file names returned.
;       SORT:           in, optional, type=boolean, default=0
;                       If set, results will be sorted alphabetically by `FILENAMES`.
;       TT2000:         in, optional, type=boolean, default=0
;                       If set, `TIMETAG` and `MODIFIED` will be returned as tt2000 values.
;-
PRO MrMMS_SDC_API::FileInfo, filenames, filesize, modified, start_date, end_date, $
COUNT=nFiles, $
SORT=tf_sort, $
TT2000=tt2000
	Compile_Opt idl2

	;Catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		
		;Delete the tar file
		IF N_Elements(file) GT 0 && File_Test(file, /REGULAR) $
			THEN File_Delete, file
		
		;Delete the JSON file
		IF N_Elements(jsonfile) GT 0 && File_Test(jsonfile, /REGULAR) $
			THEN File_Delete, jsonfile
			
		MrPrintF, 'LogErr'
		
		;Return to the previous setting
		IF N_Elements(info_type) GT 0 && info_type NE 'file_info' $
			THEN self -> SetProperty, INFO_TYPE=info_type
		
		RETURN
	ENDIF
	
	;Offline mode
	IF self.offline THEN Message, 'OFFLINE=1: Cannot obtain file info.'
	
;-----------------------------------------------------
; Download File Info \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;CD to the Info service
	info_type = self.info_type
	IF info_type NE 'file_info' THEN self -> SetProperty, /FILE_INFO
	
	;Download the JSON file
	file = self -> Download()
	
	;Return to the previous setting
	IF info_type NE 'file_info' THEN self -> SetProperty, INFO_TYPE=info_type

	;Replace the ".tar.gz" extension with ".json"
	jsonfile  = StRegEx(file, '(.*)\.tar\.gz', /SUBEXP, /EXTRACT)
	jsonfile  = jsonfile[1]
	jsonfile += '.json'
	
	;Rename the file
	File_Move, file, jsonfile
	
;-----------------------------------------------------
; Read Info from File \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Read the file
	jsonstr = ''
	OpenR, lun, jsonfile, /GET_LUN
	ReadF, lun, jsonstr, FORMAT='(a)'
	Free_LUN, lun
	
	;Delete the file
	File_Delete, jsonfile
	
;-----------------------------------------------------
; Parse Info \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Parse string
	json = JSON_Parse(jsonstr)

	;Parse further
	nFiles     = N_Elements(json['files'])
	filenames  = StrArr(nFiles)
	start_date = StrArr(nFiles)
	end_date   = StrArr(nFiles)
	modified   = StrArr(nFiles)
	filesize   = ULonArr(nFiles)
	foreach value, json['files'], idx DO BEGIN
		;Get the file info
		file_hash = json['files', idx]
		
		;Extract the file info
		filenames[idx] = file_hash['file_name']
		modified[idx]  = file_hash['modified_date']
		filesize[idx]  = file_hash['file_size']
		
		;Extract file times
		;   - Normal data files have TIMETAG key
		;   - Ancillary data files have START_DATE and END_DATE keys
		IF file_hash -> HasKey('timetag') THEN BEGIN
			start_date[idx] = file_hash['timetag']
		ENDIF ELSE BEGIN
			start_date[idx] = file_hash['start_date']
			end_date[idx]   = file_hash['end_date']
		ENDELSE
	endforeach
	
;-----------------------------------------------------
; Organize Results \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Sort
	IF Keyword_Set(tf_sort) THEN BEGIN
		isort      = Sort(filenames)
		filenames  = filenames[isort]
		modified   = modified[isort]
		start_date = start_date[isort]
		end_date   = end_date[isort]
		filesize   = filesize[Temporary(isort)]
	ENDIF
	
	;Scalar
	IF nFiles EQ 1 THEN BEGIN
		filenames  = filenames[0]
		modified   = modified[0]
		start_date = start_date[0]
		end_date   = end_date[0]
		filesize   = filesize[0]
	ENDIF
	
	;TT2000 times
	IF Keyword_Set(tt2000) THEN BEGIN
		;Convert START_DATE to TT2000
		CDF_TT2000, start_date, Fix(StrMid(start_date,  0, 4)), Fix(strmid(start_date,  5, 2)), Fix(StrMid(start_date,  8, 2)), $
		                        Fix(strmid(start_date, 11, 2)), Fix(strmid(start_date, 14, 2)), Fix(StrMid(start_date, 17, 2)), $
		                        /COMPUTE_EPOCH
		
		;Convert END_DATE to TT2000
		cdf_tt2000, end_date, Fix(StrMid(end_date,  0, 4)), Fix(StrMid(end_date,  5, 2)), Fix(StrMid(end_date,  8, 2)), $
		                      Fix(StrMid(end_date, 11, 2)), Fix(StrMid(end_date, 14, 2)), Fix(StrMid(end_date, 17, 2)), $
		                      /COMPUTE_EPOCH
		
		;Convert MODIFIED to TT2000
		cdf_tt2000, modified, Fix(StrMid(modified,  0, 4)), Fix(StrMid(modified,  5, 2)), Fix(StrMid(modified,  8, 2)), $
		                      Fix(StrMid(modified, 11, 2)), Fix(StrMid(modified, 14, 2)), Fix(StrMid(modified, 17, 2)), $
		                      /COMPUTE_EPOCH
	ENDIF
END


;+
;   Get names of files
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of file names returned.
;
; :Returns:
;       FILENAMES:      Names of files that match search criteria
;-
FUNCTION MrMMS_SDC_API::FileNames, $
COUNT=count
	Compile_Opt idl2
	On_Error, 2
	
	IF self.oWeb.offline THEN Message, 'OFFLINE=1: Cannot search for files.'
	
;-----------------------------------------------------
; Download Filenames \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;CD to the File Names service
	count     = 0
	info_type = self.info_type
	IF info_type NE 'file_names' THEN self -> SetProperty, /FILE_NAMES

	;Download the JSON file
	fnames = self -> Download(/STRING)
	
	;Return to the original URI
	IF info_type NE 'file_names' THEN self -> SetProperty, INFO_TYPE=info_type
	
;-----------------------------------------------------
; Finish \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Create string array
	fnames = StrSplit(fnames, ',', /EXTRACT, COUNT=count)
	
	;Filter the results
	IF count GT 0 THEN fnames = self -> FilterFiles(fnames, COUNT=count, /TIME)
	
	;Return
	RETURN, fnames
END


;+
;   Filter files by time and version.
;
; :Private:
;
; :Params:
;       FILENAMES:      in, required, type=string/strarr
;                       File names to be filtered.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of files that pass the filter.
;       TIME:           in, optional, type=boolean
;                       If set, files will be filtered by start (and end) times. If
;                           neither TIME nor `VERSION` are set, then the default is to
;                           filter by both time and version.
;       VERSION:        in, optional, type=boolean
;                       If set, files will be filtered by file version number. If
;                           neither `TIME` nor VERSION are set, then the default is to
;                           filter by both time and version.
;
; :Returns:
;       FILES:          Files that pass the filter.
;-
FUNCTION MrMMS_SDC_API::FilterFiles, filenames, $
COUNT=count, $
TIME=time, $
VERSION=version
	Compile_Opt idl2
	On_Error, 2
	
	tf_time    = Keyword_Set(time)
	tf_version = Keyword_Set(version)
	IF tf_time EQ 0 && tf_version EQ 0 THEN BEGIN
		tf_time    = 1B
		tf_version = 1B
	ENDIF

;---------------------------------------------------------------------
; Find Unique File Types /////////////////////////////////////////////
;---------------------------------------------------------------------
	;Start and end time
	self.oWeb -> GetProperty, DATE_START=tstart, DATE_END=tend
	
	;File's base name
	MrMMS_Parse_Filename, filenames, SC=sc, INSTR=instr, MODE=mode, LEVEL=level, OPTDESC=optdesc
	fbase = sc + '_' + instr + '_' + mode + '_' + level
	iOpt = Where(optdesc NE '', nOpt)
	IF nOpt GT 0 THEN fbase[iOpt] += '_' + optdesc[iOpt]

	;Unique file types
	iUniq = Uniq(fbase, Sort(fbase))
	nUniq = N_Elements(iUniq)
	
;---------------------------------------------------------------------
; Filter Each File Type //////////////////////////////////////////////
;---------------------------------------------------------------------
	count = 0
	FOR i = 0, nUniq - 1 DO BEGIN
		;Get files of similar type
		fType = fbase[iUniq[i]]
		idx   = Where( fbase eq fType, nIdx )
	
	;---------------------------------------------------------------------
	; Ancillary Files ////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		nTemp = 0
		fTemp = ''
		IF StRegEx(fType, '(DEF|PRED|ATTVAL|NAV|TIMEBIAS)', /BOOLEAN) THEN BEGIN
			;Time filter
			IF tf_time THEN BEGIN
				fTemp = self -> Anc_FilterTime( filenames[idx], $
				                                COUNT = nTemp)
			ENDIF ELSE BEGIN
				fTemp = filenames[idx]
				nTemp = nIdx
			ENDELSE
			
			;Version filter
			IF nTemp GT 0 && tf_version THEN BEGIN
				fTemp = self -> Anc_FilterVersion( fTemp, COUNT=nTemp )
			ENDIF
	
	;---------------------------------------------------------------------
	; Normal Files ///////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		ENDIF ELSE BEGIN
			;Time filter
			IF tf_time THEN BEGIN
				fTemp = self -> FilterTime( filenames[idx], $
				                            COUNT = nTemp)
			ENDIF ELSE BEGIN
				fTemp = filenames[idx]
				nTemp = nIdx
			ENDELSE
			
			;Version filter
			IF nTemp GT 0 && tf_version THEN BEGIN
				fTemp = self -> FilterVersion(fTemp, COUNT=nTemp)
			ENDIF
		ENDELSE
	
	;---------------------------------------------------------------------
	; Record Results /////////////////////////////////////////////////////
	;---------------------------------------------------------------------
		
		;Skip if no files pass filter
		IF nTemp EQ 0 THEN CONTINUE
		
		;Otherwise keep the files
		files  = count EQ 0 ? fTemp : [files, temporary(fTemp)]
		count += nTemp
	ENDFOR
	
	;Return
	IF COUNT EQ 0 THEN files = ''
	RETURN, files
END


;+
;   Filter files by time. Burst and Srvy files should be filtered separately.
;
; :Private:
;
; :Params:
;       FILENAMES:      in, required, type=string/strarr
;                       File names to be filtered.
;       TSTART:         in, required, type=string
;                       Start time of filter.
;       TEND:           in, required, type=string
;                       End time of filter.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of files that pass the filter.
;
; :Returns:
;       FILES_OUT:      Files that pass the filter.
;-
FUNCTION MrMMS_SDC_API::FilterTime, filenames, $
COUNT=count
	Compile_Opt idl2
	On_Error, 2

	self.oWeb -> GetProperty, DATE_START=tstart, DATE_END=tend

;------------------------------------;
; Sort Time                          ;
;------------------------------------;
	;Parse the file names
	MrMMS_Parse_Filename, filenames, TSTART=fstart

	;Parse the start times
	;   - Convert to TT2000
	MrMMS_Parse_Time, Temporary(fstart), TT2000=tt2000
;	CDF_TT2000, tt2000, year, month, day, hour, minute, second, /COMPUTE_EPOCH

	;Sort by time
	nfiles    = N_Elements(filenames)
	isort     = Sort(tt2000)
	files_out = filenames[isort]
	tt2000    = tt2000[isort]

;------------------------------------;
; Filter Time                        ;
;------------------------------------;

	;TT2000 values for the start and end times
	MrMMS_Parse_Time, StrJoin(StrSplit(tstart, '-T:', /EXTRACT)), TT2000=tt2000_start
	MrMMS_Parse_Time, StrJoin(StrSplit(tend, '-T:', /EXTRACT)),   TT2000=tt2000_end

	;Filter files by end time
	;   - Any files that start after TEND can be discarded
	iend = where( tt2000 LE tt2000_end[0], count )
	IF count GT 0 THEN BEGIN
		files_out = files_out[iend]
		tt2000    = tt2000[iend]
	ENDIF

	;Filter files by BEGIN time
	;   - Any file with TSTART < TRANGE[0] can potentially have data
	;     in our time interval of interest.
	;   - Assume the start time of one file marks the end time of the previous file.
	;   - With this, we look for the file that begins just prior to TRANGE[0] and
	;     throw away any files that start before it.
	istart = Where( tt2000 LE tt2000_start[0], nstart )
	IF nstart GT 0 THEN BEGIN
		;Select the file time that starts closest to the given time without
		;going over.
		istart = istart[nstart-1]
	
		;Find all files with start time on or after the selected time
		ifiles = Where(tt2000 GE tt2000[istart], count)
		IF nstart GT 0 THEN BEGIN
			files_out = files_out[ifiles]
			tt2000    = tt2000[ifiles]
		ENDIF
	ENDIF

	;Number of files kept
	;   - If TRANGE[0] < TSTART of first file, COUNT will be zero
	;   - However, there still may be files with TRANGE[0] < TSTART < TRANGE[1]
	IF count EQ 0 THEN RETURN, ''

	;The last caveat:
	;   - Our filter may be too lenient. The first file may or may not contain
	;     data within our interval.
	;   - Check if it starts on the same day. If not, toss it
	;   - There may be many files with the same start time, but different
	;     version numbers. Make sure we get all copies OF the first start
	;     time.
	CDF_TT2000, tt2000_start, syr, smo, sday, /BREAKDOWN, /TOINTEGER
	CDF_TT2000, tt2000[0], year, month, day, hour, /BREAKDOWN, /TOINTEGER
	IF year NE syr || month NE smo || day NE sday THEN BEGIN
		;Filter out all files matching the first starting date
		ibad = Where( tt2000 EQ tt2000[0], nbad, COMPLEMENT=igood, NCOMPLEMENT=count)
	
		;Extract files
		IF count GT 0 $
			THEN files_out = files_out[igood] $
			ELSE files_out = ''
	ENDIF
	
	RETURN, files_out
END


;+
;   Filter files by version number.
;
; :Private:
;
; :Params:
;       FILENAMES:      in, required, type=string/strarr
;                       File names to be filtered.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of files that pass the filter.
;       LATEST_VERSION: in, optional, type=boolean, default=0
;                       If set, only the newest version of the files is returned. This
;                           is the default if `MIN_VERSION` and `VERSION` are not set.
;       MIN_VERSION:    in, optional, type=string, default=''
;                       A version number formatted as 'X.Y.Z'. Files with
;                           versions smaller this version will be filtered out.
;       VERSION:        in, optional, type=string, default=''
;                       A version number formatted as 'X.Y.Z'. Only files of
;                           this version will be returned.
;
; :Returns:
;       FILES_OUT:      Files that pass the filter.
;-
FUNCTION MrMMS_SDC_API::FilterVersion, filenames, $
COUNT=count, $
LATEST_VERSION=latest_version, $
MIN_VERSION=min_version, $
VERSION=version
	Compile_Opt idl2
	On_Error, 2

	;Defaults
	tf_latest = Keyword_Set(latest_version)
	tf_checkv = N_Elements(version)     GT 0
	tf_minv   = N_Elements(min_version) GT 0
	IF tf_checkv + tf_minv + tf_latest  EQ 0 THEN tf_latest = 1B
	IF tf_checkv + tf_minv + tf_latest  GT 1 $
		THEN Message, 'VERSION, MIN_VERSION and LATEST_VERSION are mutually exclusive.'

	;Allocate memory to results
	files_out = StrArr(N_Elements(filenames))

;------------------------------------;
; Filter Unique File Types           ;
;------------------------------------;
	
	;Find unique file names
	;   - Remake the files name without version numbers
	MrMMS_Parse_Filename, filenames, SC=sc, INSTR=instr, MODE=mode, LEVEL=level, TSTART=tstart
	f_no_version = sc + '_' + instr + '_' + mode + '_' + level + '_' + tstart
	iUniq = Uniq(f_no_version, Sort(f_no_version))
	nUniq = N_Elements(iUniq)

	;Step through each unique file type
	count = 0
	FOR i = 0, nUniq - 1 DO BEGIN
		;Files that match the current file type
		iFiles = Where(f_no_version EQ f_no_version[iUniq[i]], nFiles)
		temp   = filenames[iFiles]

		;Extract X, Y, Z version numbers from file
		fversion = StRegEx(temp, 'v([0-9]+)\.([0-9]+)\.([0-9]+)\.cdf$', /SUBEXP, /EXTRACT)
		fv = Fix(fversion[1:3,*])
	
	;------------------------------------;
	; Minimum Version                    ;
	;------------------------------------;
		IF tf_minv THEN BEGIN
			;Version numbers to match
			min_vXYZ = Fix(StrSplit(min_version, '.', /EXTRACT))
		
			;Select file versions
			iv = Where(   (fv[0,*] GT min_vXYZ[0]) or $
			            ( (fv[0,*] EQ min_vXYZ[0]) and (fv[1,*] GT min_vXYZ[1]) ) or $
			            ( (fv[0,*] EQ min_vXYZ[0]) and (fv[1,*] EQ min_vXYZ[1]) and (fv[2,*] GE min_vXYZ[2]) ), nMatch )
			
			;Keep matches
			IF nMatch GT 0 THEN BEGIN
				files_out[count:nMatch-1] = temp[iv]
				count += nMatch
			ENDIF
	
	;------------------------------------;
	; Exact Version                      ;
	;------------------------------------;
		ENDIF ELSE IF tf_checkv THEN BEGIN
			;Version numbers to match
			vXYZ  = Fix(StrSplit(version, '.', /EXTRACT))
			
			;Select file versions
			iv = Where( (fv[0,*] EQ vXYZ[0]) and (fv[1,*] EQ vXYZ[1]) and (fv[2,*] EQ vXYZ[2]), nMatch )
			
			;Keep matches
			IF count GT 0 THEN BEGIN
				files_out[count:nMatch-1] = T[iv]
				count += nMatch
			ENDIF
	
	;------------------------------------;
	; Latest Version                     ;
	;------------------------------------;
		ENDIF ELSE BEGIN
			;Step through X, Y, Z version numbers
			FOR j = 0, 2 DO BEGIN
				;Find latest X, Y, Z version
				iv = Where( fv[j,*] EQ max(fv[j,*]), nMatch)
				
				;Select version
				fv   = fv[*,iv]
				temp = temp[iv]
			ENDFOR
			
			;Keep matches
			IF nMatch GT 0 THEN BEGIN
				files_out[count:count+nMatch-1] = temp
				count += nMatch
			ENDIF
		ENDELSE
	ENDFOR

	;Return
	CASE count OF
		0:    files_out = ''
		1:    files_out = files_out[0]
		ELSE: files_out = files_out[0:count-1]
	ENDCASE

	RETURN, files_out
end


;+
;   Download files one-by-one.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of file names returned.
;
; :Returns:
;       FILENAMES:      Names of files that match search criteria
;-
FUNCTION MrMMS_SDC_API::Get, $
COUNT=count
	Compile_Opt idl2

	;Catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		
		;Delete progress bar
		IF N_Elements(progress_bar) GT 0 && Obj_Valid(progress_bar) $
			THEN Obj_Destroy, progress_bar
		
		;Report error
		MrPrintF, 'LogErr'
		
		;CD to original directory
		;   - ::CD calls ::GetResponseHeader, which purposely generates an
		;     error. Report this error before calling ::CD
		IF N_Elements(state) GT 0 THEN self -> SetProperty, _STRICT_EXTRA=state
		self.oWeb -> SetProperty, CALLBACK_FUNCTION='MrWebURI_Callback'
		
		count = 0
		RETURN, ''
	ENDIF
	
	;Get the current state so that we can return to it.
	state = self -> GetState()
	
	;Search FOR local & remote files
	self -> Search, flocal, fremote, $
	                NLOCAL  = nLocal, $
	                NREMOTE = nRemote
	IF nLocal + nRemote EQ 0 THEN Message, 'No local or remote files found.'

	;Get directories
	self.oWeb -> GetProperty, LOCAL_ROOT  = local_root, $
	                          NO_DOWNLOAD = no_download, $
	                          VERBOSE     = verbose
	
	;Allocate memory
	count = nLocal + nRemote
	fout  = StrArr(count)

;-----------------------------------------------------
; Local Files \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF nLocal GT 0 THEN BEGIN
		self.oWeb       -> ParseURI, Temporary(flocal), PATH=fpath
		fout[0:nLocal-1] = Temporary(fpath)
	ENDIF

;-----------------------------------------------------
; No Download \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Do not download
	IF nRemote EQ 0 || no_download THEN BEGIN
		
		;Information
		IF verbose THEN BEGIN
			nSkip = count-nLocal
			IF nSkip GT 0 THEN MrPrintF, 'LogText', nSkip, FORMAT='(%"%i remote files not downloaded.")'
			MrPrintF, 'LogText', 'Local files:'
			MrPrintF, 'LogText', '    ' + fout[0:nLocal-1]
		ENDIF
		
		;Trim results
		count = nLocal
		fout  = count EQ 0 ? '' : fout[0:count-1]

;-----------------------------------------------------
; Download \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	ENDIF ELSE BEGIN 
		IF verbose THEN BEGIN
			MrPrintF, 'LogText', nLocal, count, FORMAT='(%"%i of %i files found locally.")'
			MrPrintF, 'LogText', '    ' + fout[0:nLocal-1]
			MrPrintF, 'LogText', nRemote, FORMAT='(%"Downloading %i files.")'
		ENDIF
	
	;-----------------------------------------------------
	; Progress Bar \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------

		;theText = string(FORMAT='(%"File %i OF %i.")', i+1, nFiles)    ;Another choice FOR the text string.
		progressBar = Obj_New()
;		progressBar = Obj_New('cgProgressBar', TEXT=fileOut, TITLE='Downloading...', $
;		                      /CANCEL, /START)
	
		;Create and set the callback structure
		callbackData = { bar:          progressBar, $
		                 last_update:  0.0, $
		                 last_size:    0L, $
		                 oNet:         self.oWeb, $
		                 percent:      0B, $
		                 total_size:   0LL $
		               }
		
		self.oWeb -> SetProperty, CALLBACK_FUNCTION='MrWebURI_Callback_QuietDL' ;'MrWebURI_Callback_Download'
	
	;-----------------------------------------------------
	; Download \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Name of downloaded file
		temp_local = MrMMS_Build_Filename(fremote, SDC_ROOT=local_root)
		local_dest = File_Dirname(temp_local)
		temp_fname = File_Basename(temp_local)
		
		;Switch to Download api
		info_type = self.info_type
		IF info_type NE 'download' THEN self -> SetProperty, /DOWNLOAD

		;Download remote files
		;   - Only the latest version are downloaded. If an older version
		;     is specified, no file will be returned.
		FOR i = 0, nRemote-1 DO BEGIN
			IF verbose THEN MrPrintF, 'LogText', 'Downloading "' + temp_fname[i] + '".'
			
			;Reset the progress bar and callback data
;			progressBar -> Update, 0.0
			self.oWeb -> SetProperty, CALLBACK_DATA=callbackData
			
			;Prep for download
			;   - Set the file for download
			;   - Make sure the destination directory exists
			self -> SetProperty, FILES=temp_fname[i]
			IF ~File_Test(local_dest[i], /DIRECTORY) THEN File_MKDir, local_dest[i]

			;Download the file.
			;   - TODO: If ::Download fails, remove element from FILES_OUT.
			;           Also, if operation is cancelled, let MrMMS_Load_Data know
			fout[nLocal+i] = self -> Download(temp_local[i])
			print, ''
			
			IF verbose THEN MrPrintF, 'LogText', temp_local[i], FORMAT='(%"File downloaded to \"%s\".")'
		ENDFOR
	
	;-----------------------------------------------------
	; Clean Up \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Destroy progress bar
		IF nRemote GT 0 THEN BEGIN
			;Destroy the progress bar
			obj_destroy, progressBar
			
			;Reset the callback FUNCTION
			self.oWeb -> SetProperty, CALLBACK_DATA     = 0B, $
			                          CALLBACK_FUNCTION = 'MrWebURI_Callback'
		ENDIF
	ENDELSE

;-----------------------------------------------------
; Return \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Return to the previous configuration
	self -> SetProperty, _STRICT_EXTRA=state
	
	IF count EQ 1 THEN fout = fout[0]
	RETURN, fout
END


;+
;   Download files in bulk. A zip file is retrieved then unzipped.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of file names returned.
;
; :Returns:
;       FILENAMES:      Names of files that match search criteria
;-
function MrMMS_SDC_API::GetBulk, $
UNPACK=unpack, $
COUNT=count
	compile_opt idl2

	;Catch errors
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		
		;Report error
		MrPrintF, 'LogErr'
		
		;CD to original directory
		;   - ::CD calls ::GetResponseHeader, which purposely generates an
		;     error. Report this error before calling ::CD
		if n_elements(uri) gt 0 then self -> CD, uri
		
		count = 0
		return, ''
	endif
	
	tf_unpack = keyword_set(unpack)
	
	;Set URI to the Download service
	uri = self -> GetURI()
	self -> CD, /DOWNLOAD
	
	;Get the zip file
	file = self -> Download()
	
	;Unpack
	if tf_unpack then file = self -> GetBulk_Unpack(file, COUNT=count)
	
	;Return files
	return, file
end


;+
;   Unpackage a zip file of data files that have been downloaded in bulk.
;
; :Private:
;
; :Params:
;       ZIPFILE:        in, required, type=string
;                       Name of the zip file to be unpacked.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of file names returned.
;
; :Returns:
;       FILENAMES:      Names of files that match search criteria
;-
function MrMMS_SDC_API::GetBulk_Unpack, zipfile
	compile_opt idl2
	on_error, 2

	;Note to self.
	if MrCmpVersion('8.2.3') le 0 $
		then message, 'FILE_ZIP routines were introduced in IDL 8.2.3. Update this program.', /INFORMATIONAL

	;
	; The zip file directory structure is:
	;
	;     mms/public/[sc]/[instr]/[mode]/[level]/[year]/[month]
	;
	; Wheras the SDC directory structure is
	;
	;     mms/public/data/[sc]/[instr]/[mode]/[level]/[year]/[month]
	;

	;Check the LOCAL_DATA_ROOT system variable
	;   - Default to creating an SDC-like directory structure
	mms_data_root = getenv('MMS_DATA_ROOT')
	if mms_data_root eq '' then begin
		;Set it equal to $HOME/MrWebData
		mms_data_root = filepath('', ROOT_DIR=getenv('HOME'), SUBDIRECTORY='MrWebData')
		if ~file_test(mms_data_root, /DIRECTORY) then file_mkdir, mms_data_root
	
		;Set environment variable
		MrPrintF, 'LogWarn', 'LOCAL_DATA_ROOT environment variable not defined.'
		MrPrintF, 'LogWarn', 'Setting it equal to "' + mms_data_root + '".'
		setenv, 'MMS_DATA_ROOT=' + mms_data_root
	endif

	;
	;Download structure
	;   self.directory/file.zip                          [zipPath/zipfile]
	;       ->  unzip
	;                 /CAA_Download_[date]_[time]               [/zipRoot]
	;                 /CAA_Download_[date]_[time]/Dataset
	;                 /CAA_Download_[date]_[time]/Dataset/datafile.cdf
	;
	;Variable names
	;   zipPath     -   self.directory
	;   zipfile     -   zipPath/filename.zip            -   CAA_Download_[date]_[time].zip
	;   zipRoot     -   CAA_Download_[date]_[time]      -   [directory]
	;   unzipDir    -   zipPath/zipRoot
	;


	;Change to the directory into which the zip file was downloaded
	;   - pwd -> self.directory [= zipPath]
	cd, CURRENT=pwd
	download_dir = file_dirname(zipfile)
	cd, download_dir

	;Unzip & Delete file
	;   - mms/public/[sc]/[instr]/[mode]/[level]/[year]/[month]/
	spawn, 'unzip -qo ' + zipfile
	file_delete, zipfile
	
	;Get all files just downloaded
	;   - Use only the relative directory
	unzip_dir = filepath('', ROOT_DIR=download_dir, SUBDIRECTORY='mms')
	filenames = file_search(unzip_dir, 'mms*.cdf', COUNT=nFiles)
	files_out = strarr(nFiles)
	
	;Move to destination
	for i = 0, nFiles - 1 do begin
		;Build the output directory
		dir_parts = strsplit(file_dirname(filenames[i]), path_sep(), /EXTRACT)
		dir_out   = filepath('', ROOT_DIR=mms_data_root, SUBDIRECTORY=dir_parts[2:*])
		if ~file_test(dir_out, /DIRECTORY) then file_mkdir, dir_out
		
		;Filename
		files_out[i] = filepath(file_basename(filenames[i]), ROOT_DIR=dir_out)
		
		;Move the file
		file_move, filenames[i], files_out[i], /OVERWRITE
	endfor
	
	;Delete the unzipped directory tree
	file_delete, unzip_dir, /RECURSIVE
	
	;Change back to the current directory
	cd, pwd
	return, files_out
end


;+
;   Get object properties. Properties are grouped by function.
;
; :Keywords:
;       ANCILLARY:      out, optional, type=boolean
;                       If set, the request will be for anciallary data products.
;                           Cannot be used with the `HK`, `SCIENCE`, or `VERSION`.
;       DATA_TYPE:      out, optional, type=string
;                       Identifies which type of data is being downloaded. Possible values
;                           are: {'ancillary' | 'hk' | 'science'}.
;       HK:             out, optional, type=boolean
;                       If set, the request will be for house keeping data products.
;                           Cannot be used with the `ANCILLARY`, `SCIENCE`, or `VERSION`.
;       SCIENCE:        out, optional, type=boolean
;                       If set, the request will be for science data products.
;                           Cannot be used with the `ANCILLARY`, or `HK`. If
;                           None of the three keywords are set, the current path will
;                           be checked for one of those options. If still nothing is found
;                           then SCIENCE will be the default.
;
;       DOWNLOAD:       out, optional, type=boolean
;                       If set, the URL will request that files be downloaded. Cannot
;                           be used with `INFO`, `NAMES`, or `V_INFO`. If none of the
;                           four keywords are set, the current URI will be examined for
;                           a choice. If still nothing has been selected, DOWNLOAD will
;                           be the default.
;       FILE_INFO:      out, optional, type=boolean
;                       If set, the URL will request file information. Cannot be used
;                           with `DOWNLOAD`, `NAMES`, or `V_INFO`
;       FILE_NAMES:     out, optional, type=boolean
;                       If set, the URL will requets file names. Cannot be used with
;                           `DOWNLOAD`, `INFO`, or `V_INFO`
;       INFO_TYPE:      out, optional, type=string
;                       Identifies which type of info is being requested. Possible values
;                           are: {'download' | 'file_info' | 'file_names' | 'version_info'}.
;       V_INFO:         out, optional, type=boolean
;                       If set, the request will be for file version information.
;                           Cannot be used with the `ANCILLARY` or `HK`, or with
;                           `DOWNLOAD`, `INFO`, or `NAMES`.
;
;       PUBLIC_SITE:    out, optional, type=boolean
;                       If set, the public site will be queried (no password required,
;                           contains data of quality L2 and above). Cannot be used
;                           with `TEAM`. If neither keywords are set, the current URL
;                           is examined. If still neither are chosen, PUBLIC is the default.
;       SITE:           out, optional, type=string
;                       Identifies the site from which info is obtained. Possible values
;                           are {'public' | 'team'}
;       TEAM_SITE:      out, optional, type=string/strarr
;                       If set, the public site will be queried (password required,
;                           contains data of quality L1A and above). Cannot be used
;                           with `TEAM`.
;
;       ANC_PRODUCT:    out, optional, type=string/strarr
;                       Ancillary data product names. Setting this keyword will
;                           automatically set `ANCILLARY`.
;       FILES:          out, optional, type=string/strarr
;                       Names of the files to be downloaded. If set, then `RESET_QUERY`, 
;                           will be set.
;       INSTR:          out, optional, type=string/strarr
;                       Instrument ID. If INSTR='hk', the `HK` keyword will be auto-set.
;       LEVEL:          out, optional, type=string/strarr
;                       Data quality level.
;       MODE:           out, optional, type=string/strarr
;                       Data rate mode.
;       OPTDESC:        out, optional, type=string/strarr
;                       Optional file name descriptor.
;       SC_ID:          out, optional, type=string/strarr
;                       Spacecraft ID. Options are: {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       VERSION:        out, optional, type=string/strarr
;                       File name version number, formatted as X.Y.Z, where X, Y, and Z
;                           are integers.
;
;       _REF_EXTRA:     out, optional, type=any
;                       Any keyword accepted by MrWebURI::GetProperty is also accepted here.
;-
PRO MrMMS_SDC_API::GetProperty, $
ANC_PRODUCT=ancprod, $
ANCILLARY=ancillary, $
DATA_TYPE=data_type, $
DOWNLOAD=download, $
FILES=files, $
FILE_INFO=info, $
FILE_NAMES=names, $
HK=hk, $
INFO_TYPE=info_type, $
INSTR=instr, $
LEVEL=level, $
MODE=mode, $
OPTDESC=optdesc, $
PUBLIC_SITE=public, $
SC_ID=sc_id, $
SCIENCE=science, $
TEAM_SITE=team, $
V_INFO=v_info, $
VERSION=version, $
_REF_EXTRA=extra
	Compile_Opt idl2
	On_Error, 2
	
	;Get properties
	IF Arg_Present(ancprod)   THEN ancprod   = StrSplit(self.ancprod, ',', /EXTRACT)
	IF Arg_Present(ancillary) THEN ancillary = self.data_type EQ 'ancillary'
	IF Arg_Present(data_type) THEN data_type = self.data_type
	IF Arg_Present(download)  THEN download  = self.info_type EQ 'download'
	IF Arg_Present(files)     THEN files     = StrSplit(self.files, ',', /EXTRACT)
	IF Arg_Present(info)      THEN info      = self.info_type EQ 'file_info'
	IF Arg_Present(names)     THEN names     = self.info_type EQ 'file_names'
	IF Arg_Present(hk)        THEN hk        = self.instr EQ 'hk'
	IF Arg_Present(info_type) THEN info_type = self.info_type
	IF Arg_Present(instr)     THEN instr     = StrSplit(self.instr, ',', /EXTRACT)
	IF Arg_Present(level)     THEN level     = StrSplit(self.level, ',', /EXTRACT)
	IF Arg_Present(mode)      THEN mode      = StrSplit(self.mode, ',', /EXTRACT)
	IF Arg_Present(optdesc)   THEN optdesc   = StrSplit(self.optdesc, ',', /EXTRACT)
	IF Arg_Present(public)    THEN public    = self.site EQ 'public'
	IF Arg_Present(sc_id)     THEN sc_id     = StrSplit(self.sc, ',', /EXTRACT)
	IF Arg_Present(science)   THEN science   = self.data_type EQ 'science'
	IF Arg_Present(team)      THEN team      = self.site EQ 'team'
	IF Arg_Present(site)      THEN site      = self.site
	IF Arg_Present(v_info)    THEN v_info    = self.info_type EQ 'version_info'
	IF Arg_Present(version)   THEN version   = self.version

	;MrWebURI properties
	IF N_Elements(extra) GT 0 THEN self.oWeb -> GetProperty, _STRICT_EXTRA=extra
END


;+
;   Create a structure containing the current object properties. This is useful
;   for the ::Get method, since the FILES keyword to ::SetProperty resets all
;   other properties. The output structure can be passed to the SetProperty
;   method via _STRICT_EXTRA.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of file names returned.
;
; :Returns:
;       STATE:          out, required, type=struct
;                       Current property values of the object.
;-
FUNCTION MrMMS_SDC_API::GetState
	Compile_Opt idl2
	On_Error, 2
	
	;Get the date/time range
	self.oWeb -> GetProperty, DATE_START=date_start, DATE_END=date_end
	
	;Current state of the object.
	state = { site:        self.site, $
	          info_type:   self.info_type, $
	          data_type:   self.data_type, $
	          sc_id:       self.sc, $
	          instr:       self.instr, $
	          mode:        self.mode, $
	          level:       self.level, $
	          optdesc:     self.optdesc, $
	          version:     self.version, $
	          date_start:       date_start, $
	          date_end:         date_end, $
	          anc_product: self.ancprod, $
	          files:       self.files $
	        }

	RETURN, state
END


;+
;   Determine if file names correspond to ancillary files.
;
; :Private:
;
; :Params:
;       FILENAMES:      in, required, type=string/strarr
;                       Names of files.
;
; :Returns:
;       TF_ANCE:        out, required, type=byte/bytarr
;                       Returns true for each file that is an ancillary file and fase for
;                           all other files.
;-
function MrMMS_SDC_API::IsAncillary, filenames
	compile_opt idl2
	on_error, 2
	
	;Return
	return, stregex(filenames, 'MMS[1-4]_[A-Z]+_[0-9]{4}[0-9]{3}_[0-9]{4}[0-9]{3}\.V[0-9]+', /BOOLEAN)
end


;+
;    Print the list of files in the current directory.
;
; :Params:
;       SEARCHSTR:      in, optional, type=string, default=''
;                       A pattern to match against files in the current directory.
;                           Files and directoires matching this search string will
;                           be returned. If the search string is preceeded by a
;                           directory path, the search will take place in that
;                           directory instead of the current directory.
;
; :Keywords:
;       COUNT:          in, optional, type=integer
;                       Number of files found.
;       INFO:           in, optional, type=boolean, default=0
;                       In addition to displaying file names, show file size, date
;                           modified, and file time tag.
;       OUTPUT:         out, optional, type=strarr/struct
;                       A named variable into which the directory contents are returned.
;                           If present, the contents are not printed to the display. If
;                           `INFO` is set, a structure is returned. Otherwise a string
;                           array of of file names is returned.
;-
pro MrMMS_SDC_API::LS, srchstr, $
COUNT=count, $
INFO=info, $
OUTPUT=output
	compile_opt idl2
	on_error, 2
	
	tf_info = keyword_set(info)

;---------------------------------------------------------------------
; File Info //////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	if tf_info then begin
		;Get Info
		self -> GetFileInfo, filenames, filesize, modified, timetag, COUNT=count

		;Return the information
		if arg_present(output) then begin
			output = { filenames:     temporary(filename), $
			           file_size:     temporary(filesize), $
			           date_modified: temporary(modified), $
			           timetag:       temporary(timetag) }
		
		;Display the information
		endif else begin
			;Format the output
			size_fmt = strtrim(max(strlen(string(filesize, FORMAT='(i0)'))), 2)
			mod_fmt  = 'a-' + strtrim(max(strlen(modified)), 2)
			tag_fmt  = 'a-' + strtrim(max(strlen(timetag)), 2)
			name_fmt = 'a-' + strtrim(max(strlen(filenames)), 2)
			
			;Header
			print, 'SIZE', 'MODIFIED', 'TIMETAG', 'FILENAME', $
			       FORMAT='(a-' + size_fmt + ', 2x, ' + mod_fmt + ', 2x, ' + tag_fmt + ', 2x, ' + name_fmt + ')'
			
			;Info
			fmt = '(i' + size_fmt + ', 2x, ' + mod_fmt + ', 2x, ' + tag_fmt + ', 2x, ' + name_fmt + ')'
			for i = 0, count-1 do print, filesize[i], modified[i], timetag[i], filenames[i], FORMAT=fmt
		endelse

;---------------------------------------------------------------------
; File Names /////////////////////////////////////////////////////////
;---------------------------------------------------------------------
	endif else begin
		files = self -> GetFileNames(COUNT=count)

		;Print or return?
		if arg_present(output) $
			then output = temporary(files) $
			else print, transpose(files)
	endelse
	
end


;+
;   Search in REMOTE_ROOT, LOCAL_ROOT, DROPBOX_ROOT, and MIRROR_ROOT for files
;   that match the input URI. Results are filtered by time and version.
;
; :Params:
;       FLOCAL:         out, optional, type=string/strarr
;                       Names of local files that match search criteria.
;       FREMOTE:        out, optional, type=string/strarr
;                       Names of remote files that match search criteria.
;
; :Keywords:
;       NLOCAL:         out, optional, type=integer
;                       Number of local files found.
;       NREMOTE:        out, optional, type=integer
;                       Number of remote files found.
;-
FUNCTION MrMMS_SDC_API::Local_FileNames, $
MIRROR=mirror, $
DROPBOX=dropbox, $
COUNT=count
	Compile_Opt idl2
	On_Error, 2
	
	;Where to look
	count      = 0
	tf_dropbox = Keyword_Set(dropbox)
	tf_mirror  = Keyword_Set(mirror)
	IF tf_dropbox && tf_mirror THEN Message, 'DROPBOX and MIRROR are mutually exclusive.'

	;To mimic the SDC, default to all possible values
	;   - To exactly mimic the SDC, we would have to do the same with OPTDESC
	;   - There are too many values for OPTDESC
	;   - Using wildcards would cause many of the searches to fail because of how MrTokens are used for dates
	;   - It is better for the user to define OPTDESC
	sc     = self.sc    NE '' ? StrSplit(self.sc,    ',', /EXTRACT) : ['mms1', 'mms2', 'mms3', 'mms4']
	instr  = self.instr NE '' ? StrSplit(self.instr, ',', /EXTRACT) : ['aspoc', 'dsp', 'edi', 'edp', 'epd-eis', 'feeps', 'fgm', 'fpi', 'hpca', 'hk', 'mec', 'scm']
	mode   = self.mode  NE '' ? StrSplit(self.mode,  ',', /EXTRACT) : ['slow', 'fast', 'srvy', 'brst']
	level  = self.level NE '' ? StrSplit(self.level, ',', /EXTRACT) : ['l1a', 'l1b', 'ql', 'sitl', 'l2pre', 'l2']
	
	;Root directory
	;   - DROPBOX is flat. The rest are SDC-like directory structures.
	dir  = ''
	root = ''
	CASE 1 OF
		tf_dropbox: dir  = self.oWeb.dropbox_root
		tf_mirror:  root = self.oWeb.mirror_root
		ELSE:       root = self.oWeb.local_root
	ENDCASE
	
	;Pre-emptive return
	IF dir EQ '' && ROOT EQ '' THEN RETURN, ''
	
	;Build the file names
	;   - Search for all versions.
	self -> GetProperty, SC_ID=sc, INSTR=instr, MODE=mode, LEVEL=level, OPTDESC=optdesc
	fpattern = MrMMS_Build_Filename( sc, instr, mode, level, $
	                                 COUNT     = nPatterns, $
	                                 OPTDESC   = optdesc, $
	                                 DIRECTORY = dir, $
	                                 SDC_ROOT  = root )
	
	;File names
	uri_root = 'file://' + (root EQ '' ? dir : root)
	oFile = MrFileURI( uri_root )
	
	;Loop over each name
	nAlloc = 100
	count  = 0
	files  = StrArr(nAlloc)
	FOR i = 0, nPatterns - 1 DO BEGIN
		;Find files
		temp = oFile -> Search( fpattern[i], $
		                        COUNT = nFiles )

		;Make room
		IF count + nFiles GE nAlloc THEN BEGIN
			files  = [files, StrArr(nAlloc > nFiles)]
			nAlloc = (2*nAlloc) > (nAlloc + nFiles)
		ENDIF
		
		;Store results
		IF nFiles GT 0 THEN BEGIN
			files[count:count+nFiles-1] = Temporary(temp)
			count += nFiles
		ENDIF
	ENDFOR
	
	;Clean up the file object
	Obj_Destroy, oFile
	
	;Trim results
	CASE count OF
		0:    files = ''
		1:    files = files[0]
		ELSE: files = files[0:count-1]
	ENDCASE
	
	RETURN, files
END


;+
;   Get information about files.
;
; :Params:
;       FILENAMES:      out, optional, type=strarr
;                       Names of the files that match search criteria.
;       FILESIZE:       out, optional, type=ulonarr
;                       Names of the files that match search criteria.
;       MODIFIED:       out, optional, type=strarr
;                       Dates files were modified, formatted as 'YYYY-MM-DDThh:mm:ss'.
;       START_DATE:     out, optional, type=strarr
;                       Time tag of files, formatted as 'YYYY-MM-DDThh:mm:ss'.
;       END_DATE:       out, optional, type=strarr
;                       End time tag of files, formatted as 'YYYY-MM-DDThh:mm:ss'. Values
;                           are equal to the empty string for non-ancillary files.
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of file names returned.
;       SORT:           in, optional, type=boolean, default=0
;                       If set, results will be sorted alphabetically by `FILENAME`.
;       TT2000:         in, optional, type=boolean, default=0
;                       If set, `TIMETAG` and `MODIFIED` will be returned as tt2000 values.
;-
pro MrMMS_SDC_API::Local_FileInfo, filenames, filesize, modified, start_date, end_date, $
COUNT=count, $
SORT=tf_sort, $
TT2000=tt2000
	compile_opt idl2
	on_error, 2

	;Defaults
	tf_tt2000 = keyword_set(tt2000)

	;Get files
	file_uris = self -> GetLocalNames(COUNT=count, /FILTER)
	iAnc      = where( self -> IsAncillary(file_uris), nAnc, COMPLEMENT=iData, NCOMPLEMENT=nData)

	;Parse the files for their paths
	self -> ParseURI, file_uris, PATH=file_path 

	;FILENAMES
	filenames = self -> Path_Basename(file_path)
	
	;START_DATE & END_DATE
	start_date = tf_tt2000 ? lon64arr(count) : strarr(count)
	end_date   = tf_tt2000 ? lon64arr(count) : strarr(count)
	if nData gt 0 then begin
		;Extract and parse the timetag
		MrMMS_Parse_Filename, filenames[iData], TSTART=temp
		MrMMS_Parse_Time, temp, year, month, day, hour, minute, second, INTEGER=tf_tt2000
		
		;Create return value
		if tf_tt2000 $
			then cdf_tt2000, temp, year, month, day, hour, minute, second, /COMPUTE_TT2000 $
			else temp = year + '-' + month + '-' + day + 'T' + hour + ':' + minute + ':' + second
		
		;Store data
		start_date[iData] = temporary(temp)
	endif
	if nAnc gt 0 then begin
		MrMMS_Anc_Parse_Filename, filenames[iAnc], TSTART=tstart, TEND=tend
		
		;Date_Start
		MrMMS_Anc_Parse_Time, tstart, year, doy, month, day, INTEGER=tf_tt2000
		if tf_tt2000 $
			then cdf_tt2000, temp_start, year, month, day, 0, 0, 0, /COMPUTE_TT2000 $
			else temp_start = year + '-' + month + '-' + day + 'T00:00:00'
		
		;Date_End
		MrMMS_Anc_Parse_Time, tend, year, doy, month, day, INTEGER=tf_tt2000
		if tf_tt2000 $
			then cdf_tt2000, temp_end, year, month, day, 0, 0, 0, /COMPUTE_TT2000 $
			else temp_end = year + '-' + month + '-' + day + 'T00:00:00'
		
		;Store data
		start_date[iAnc] = temporary(temp_start)
		end_date[iAnc]   = temporary(temp_end)
	endif
	
	;FILESIZE
	finfo    = file_info(file_path)
	filesize = finfo.size
	
	;MODIFIED
	;   - SysTime() does not accept arrays
	modified = strarr(count)
	for i = 0, count-1 do begin
		;Breakdown modified time
		caldat, systime(finfo[i].mtime, /JULIAN), month, day, year, hour, minute, second
		
		;Create TT2000 time
		if tf_tt2000 then begin
			cdf_tt2000, tt2000, year, month, day, hour, minute, second, /COMPUTE_TT2000
			modified[i] = tt2000
		
		;Create ISO time
		endif else begin
			modified[i] = string(year,   FORMAT='(i04)') + '-' + $
			              string(month,  FORMAT='(i02)') + '-' + $
			              string(day,    FORMAT='(i02)') + 'T' + $
			              string(hour,   FORMAT='(i02)') + ':' + $
			              string(minute, FORMAT='(i02)') + ':' + $
			              string(second, FORMAT='(i02)')
		endelse
	endfor
	
	;SORT
	if keyword_set(tf_sort) then begin
		isort      = sort(filenames)
		filenames  = filenames[isort]
		filesize   = filesize[isort]
		modified   = modified[isort]
		start_date = start_date[isort]
		end_date   = start_date[temporary(isort)]
	endif
end


;+
;   Get local file versions.
;
; :Params:
;       MAX_VERSION:    out, optional, type=strarr
;                       Maximum version number for each MMS data file type.
;       FILE_TYPE:      out, optional, type=strarr
;                       Types of files associated with `MAX_VERSION`
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of elements returned.
;-
pro MrMMS_SDC_API::Local_Versions, max_version, file_type, $
COUNT=count
	compile_opt idl2
	on_error, 2

	;Get files
	file_uris = self -> GetLocalNames(COUNT=count, /FILTER)
	
	;Parse the uris for their paths
	self -> ParseURI, file_uris, PATH=file_path 
	
	;Separate file name from its path
	filenames = self -> Path_Basename(file_path)
	
	;Parse file names
	MrMMS_Parse_Filename, filenames, $
	                      SC      = sc, $
	                      INSTR   = instr, $
	                      MODE    = mode, $
	                      LEVEL   = level, $
	                      OPTDESC = optdesc, $
	                      VERSION = max_version
	
	;Create file types
	file_type = sc + '_' + instr + '_' + mode + '_' + level
	iOptDesc = where(optdesc ne '', nOptDesc)
	if nOptDesc ne 0 then file_type[iOptDesc] += '_' + optdesc[iOptDesc]
	
	;Search for unique file types and versions
	;   - THIS SHOULD NOT BE NECESSARY AFTER THE FILTER!
;	iuniq       = uniq(file_type, sort(file_type))
;	count       = n_elements(iuniq)
;	file_type   = file_type[iuniq]
;	max_version = max_version[iuniq]
end


;+
;   Search in REMOTE_ROOT, LOCAL_ROOT, DROPBOX_ROOT, and MIRROR_ROOT for files
;   that match the input URI. Results are filtered by time and version.
;
; :Params:
;       FLOCAL:         out, optional, type=string/strarr
;                       Names of local files that match search criteria.
;       FREMOTE:        out, optional, type=string/strarr
;                       Names of remote files that match search criteria.
;
; :Keywords:
;       NLOCAL:         out, optional, type=integer
;                       Number of local files found.
;       NREMOTE:        out, optional, type=integer
;                       Number of remote files found.
;-
PRO MrMMS_SDC_API::Search, flocal, fremote, $
NLOCAL=nLocal, $
NREMOTE=nRemote
	Compile_Opt idl2
	On_Error, 2

	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN
	ENDIF
	
	;Initialize variables
	count   = 0
	nMirror = 0
	nLocal  = 0
	nRemote = 0
	
	self.oWeb -> GetProperty, COPY_TO_LOCAL = copy_to_local, $
	                          DROPBOX_ROOT  = dropbox_root, $
	                          LOCAL_ROOT    = local_root, $
	                          MIRROR_ROOT   = mirror_root, $
	                          OFFLINE       = offline, $
	                          REMOTE_ROOT   = remote_root

;-----------------------------------------------------
; Search Offline \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF offline THEN BEGIN
		;Search the local or mirror directory
		files = self -> Local_FileNames( COUNT=count, MIRROR=(mirror_root NE '') )
		
		;Filter by time before any additional searches
		IF count GT 0 $
			THEN files = self -> FilterFiles(files, COUNT=count, /TIME)

;-----------------------------------------------------
; Search Online \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	ENDIF ELSE BEGIN
		;Search remote
		files = self -> FileNames(COUNT=count)
		
		;REMOTE is the authority on files
		;   - Assume its files are the most recent, correct version
		IF count GT 0 THEN BEGIN
			;Search locally for equivalent files
			;   - Remove remote files that are found locally
			local_files = MrMMS_Build_Filename(files, SDC_ROOT=local_root)
			ilocal      = Where( File_Test( local_files ), nLocal, COMPLEMENT=iRemote, NCOMPLEMENT=count )
			IF nLocal GT 0 THEN local_files = local_files[iLocal]
			IF count  GT 0 THEN files       = files[iRemote]
		
			;Search mirror for equivalent files
			;   - Remove remote files that are found on the mirror
			IF count GT 0 && mirror_root NE '' THEN BEGIN
				mirror_files = MrMMS_Build_Filename(files, SDC_ROOT=mirror_root)
				iMirror      = Where( File_Test( mirror_files ), nMirror, COMPLEMENT=iRemote, NCOMPLEMENT=count )
				IF nMirror GT 0 THEN mirror_files = mirror_files[iMirror]
				IF count   GT 0 THEN files        = files[iRemote]
			ENDIF
		ENDIF
	ENDELSE

;-----------------------------------------------------
; Search Dropbox \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Full search of dropbox files
	dropbox_files = self -> Local_FileNames(COUNT=nDropbox, /DROPBOX)
	IF nDropbox GT 0 THEN BEGIN
		;Filter by time
		dropbox_files = self -> FilterFiles(dropbox_files, COUNT=nDropbox, /TIME)
		
		;We want to compare the file names, not paths
		dropboxBase = File_Basename(dropbox_files)
		fileBase    = self.oWeb -> Path_Basename(files)
		
		;OFFLINE
		;   - LOCAL & MIRROR take precedence over DROPBOX
		IF offline && count GT 0 THEN BEGIN
			tf_member = MrIsMember( fileBase, dropboxBase, COMPLEMENT=iKeep )
			IF N_Elements(iKeep) GT 0 THEN dropbox_files = dropbox_files[iKeep]
		
		;ONLINE
		;   - DROPBOX takes precedence over REMOTE
		;   - LOCAL and MIRROR files have already been removed from REMOTE
		ENDIF ELSE IF ~offline && count GT 0 THEN BEGIN
			tf_member = MrIsMember( dropboxBase, fileBase, COMPLEMENT=iKeep )
			count     = N_Elements(iKeep)
			IF count GT 0 THEN files = files[iKeep]
		ENDIF
	ENDIF

;-----------------------------------------------------
; Filter Version \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	; At this point, we have files from REMOTE, MIRROR, LOCAL, DROPBOX. The files
	; are all unique in version number, if not also in name. The next step is to
	; combine all files and filter by version.
	;
	
	;Mirror
	IF nMirror GT 0 THEN BEGIN
		files  = count EQ 0 ? Temporary(mirror_files) : [files, Temporary(mirror_files)]
		count += Temporary(nMirror)
	ENDIF
	
	;Local
	IF nLocal GT 0 THEN BEGIN
		files  = count EQ 0 ? Temporary(local_files) : [files, Temporary(local_files)]
		count += Temporary(nLocal)
	ENDIF
	
	;Dropbox
	IF nDropbox GT 0 THEN BEGIN
		files  = count EQ 0 ? Temporary(dropbox_files) : [files, Temporary(dropbox_files)]
		count += Temporary(nDropbox)
	ENDIF
	
	;Filter by version
	IF count GT 0 THEN files = self -> FilterFiles(files, COUNT=count, /VERSION)
	
	;Filter by time
	;   - Individual searches return files with FTIME <= TSTART
	;   - Try to remove as many files with FTIME < TSTART as possible
	IF count GT 0 THEN files = self -> FilterFiles(files, COUNT=count, /TIME)

;-----------------------------------------------------
; Copy DROPBOX and MIRROR Files to LOCAL \\\\\\\\\\\\\
;-----------------------------------------------------
	;Local Files
	;   - DROPBOX_ROOT and MIRROR_ROOT can share a base path. However, MIRROR_ROOT
	;     must be followed immediately by "mms[1-4]" in order to match the SDC whereas
	;     DROPBOX_ROOT should be the end of the directory structure.
	iLocal = Where( StRegEx(files, local_root,  /BOOLEAN), nLocal )
	IF mirror_root  NE '' THEN iMirror  = Where( StRegEx(files, mirror_root+path_sep()+'?mms[1-4]', /BOOLEAN), nMirror )
	IF dropbox_root NE '' THEN iDropbox = Where( StRegEx(files, dropbox_root, /BOOLEAN), nDropbox )
	
	;Remote files
	;   - ::FileNames returns files with a partial directory structure.
	;   - The easiest way (for me) to separate them is to find files that are not local
	tf_remote = count EQ 0 ? 0B : BytArr(count) + 1B
	IF nLocal   GT 0 THEN tf_remote[iLocal]   = 0B
	IF nMirror  GT 0 THEN tf_remote[iMirror]  = 0B
	IF nDropbox GT 0 THEN tf_remote[iDropbox] = 0B
	iRemote = Where(tf_remote, nRemote)
	
	;Copy
	IF copy_to_local THEN BEGIN
		;DROPBOX
		IF nDropbox GT 0 THEN BEGIN
			;Build equivalent local file names
			temp_files = MrMMS_Build_Filename(files[iDropbox], SDC_ROOT=local_root)
			
			;Copy from DROPBOX to LOCAL
			File_Copy, files[iDropbox], temp_files
			
			;Swap to local
			files[iDropbox] = Temporary(temp_files)
			iLocal          = nLocal EQ 0 ? Temporary(iDropbox) : [iLocal, Temporary(iDropbox)]
			nDropbox        = 0
		ENDIF
		
		;MIRROR
		IF nMirror GT 0 THEN BEGIN
			;Build equivalent local file names
			temp_files = FilePath( StrMid(files[iMirror], StrLen(mirror_root)), ROOT_DIR=local_dir )
			
			;Copy from DROPBOX to LOCAL
			File_Copy, files[iMirror], temp_files
			
			;Swap to local
			files[iMirror] = Temporary(temp_files)
			iLocal         = nLocal EQ 0 ? Temporary(iMirror) : [iLocal, Temporary(iMirror)]
			nMirror        = 0
		ENDIF
	ENDIF

;-----------------------------------------------------
; Separate REMOTE from LOCAL \\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Remote files
	fremote = nRemote EQ 0 ? '' : files[iRemote]
	
	;Local files
	IF nLocal GT 0 THEN flocal = files[iLocal]
	
	;Dropbox files
	IF nDropbox GT 0 THEN BEGIN
		flocal  = nLocal EQ 0 ? files[iDropbox] : [flocal, files[iDropbox]]
		nLocal += nDropbox
	ENDIF
	
	;Mirror files
	IF nMirror GT 0 THEN BEGIN
		flocal  = nLocal EQ 0 ? files[iMirror] : [flocal, files[iMirror]]
		nLocal += nMirror
	ENDIF
	
	IF nLocal  EQ 1 THEN flocal  = flocal[0]
	IF nLocal  EQ 0 THEN flocal  = ''
	IF nRemote EQ 1 THEN fremote = fremote[0]
	IF nRemote EQ 0 THEN fremote = ''
END


;+
;   Set object properties.
;
; :Keywords:
;       RESET_PATH:     in, optional, type=boolean, default=0
;                       If set, the path will be returned to the default path. This is
;                           equivalent to setting `PUBLIC`, `DOWNLOAD`, and `SCIENCE`.
;       RESET_QUERY:    in, optional, type=boolean, default=0
;                       If set, no information from previous queries will be retained.
;                           The default is to use values from the previous query unless
;                           otherwise provided via the `SC`, `INSTR`, `MODE`, `LEVEL`,
;                           `OPTDESC`, `TSTART`, `TEND`, or `ANC_PRODUCT` keywords.
;                           This keyword takes precedence over the aforementioned
;                           keywords.
;
;       ANCILLARY:      in, optional, type=boolean, default=0
;                       If set, the request will be for anciallary data products.
;                           Cannot be used with the `HK`, `SCIENCE`, or `VERSION`.
;       HK:             in, optional, type=boolean, default=0
;                       If set, the request will be for house keeping data products.
;                           Cannot be used with the `ANCILLARY`, `SCIENCE`, or `VERSION`.
;       SCIENCE:        in, optional, type=boolean, default=0
;                       If set, the request will be for science data products.
;                           Cannot be used with the `ANCILLARY`, or `HK`. If
;                           None of the three keywords are set, the current path will
;                           be checked for one of those options. If still nothing is found
;                           then SCIENCE will be the default.
;
;       DOWNLOAD:       in, optional, type=boolean, default=0
;                       If set, the URL will request that files be downloaded. Cannot
;                           be used with `INFO`, `NAMES`, or `V_INFO`. If none of the
;                           four keywords are set, the current URI will be examined for
;                           a choice. If still nothing has been selected, DOWNLOAD will
;                           be the default.
;       FILE_INFO:      in, optional, type=boolean, default=0
;                       If set, the URL will request file information. Cannot be used
;                           with `DOWNLOAD`, `NAMES`, or `V_INFO`
;       FILE_NAMES:     in, optional, type=boolean, default=0
;                       If set, the URL will requets file names. Cannot be used with
;                           `DOWNLOAD`, `INFO`, or `V_INFO`
;       V_INFO:         in, optional, type=boolean, default=0
;                       If set, the request will be for file version information.
;                           Cannot be used with the `ANCILLARY` or `HK`, or with
;                           `DOWNLOAD`, `INFO`, or `NAMES`.
;
;       PUBLIC:         in, optional, type=boolean, default=0
;                       If set, the public site will be queried (no password required,
;                           contains data of quality L2 and above). Cannot be used
;                           with `TEAM`. If neither keywords are set, the current URL
;                           is examined. If still neither are chosen, PUBLIC is the default.
;       TEAM:           in, optional, type=string/strarr
;                       If set, the public site will be queried (password required,
;                           contains data of quality L1A and above). Cannot be used
;                           with `TEAM`.
;
;       ANC_PRODUCT:    in, optional, type=string/strarr
;                       Ancillary data product names. Setting this keyword will
;                           automatically set `ANCILLARY`.
;       FILES:          in, optional, type=string/strarr, default=''
;                       Names of the files to be downloaded. If set, then `RESET_QUERY`, 
;                           will be set.
;       INSTR:          in, optional, type=string/strarr
;                       Instrument ID. If INSTR='hk', the `HK` keyword will be auto-set.
;       LEVEL:          in, optional, type=string/strarr
;                       Data quality level.
;       MODE:           in, optional, type=string/strarr
;                       Data rate mode.
;       OPTDESC:        in, optional, type=string/strarr
;                       Optional file name descriptor.
;       SC_ID:          in, optional, type=string/strarr
;                       Spacecraft ID. Options are: {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       TEND:           in, optional, type=string
;                       End time of data interval. Must be formatted as YYYY-MM-DD or
;                           YYYY-MM-DDThh:mm:ss
;       TSTART:         in, optional, type=string/strarr
;                       Start time of data interval. Must be formatted as YYYY-MM-DD or
;                           YYYY-MM-DDThh:mm:ss
;       VERSION:        in, optional, type=string/strarr
;                       File name version number, formatted as X.Y.Z, where X, Y, and Z
;                           are integers.
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrWebURI::SetProperty is also accepted here.
;-
PRO MrMMS_SDC_API::SetProperty, $
ANC_PRODUCT=ancprod, $
ANCILLARY=ancillary, $
DATA_TYPE=data_type, $
DATE_END=date_end, $
DATE_START=date_start, $
DOWNLOAD=download, $
HK=hk, $
FILES=files, $
FILE_INFO=info, $
FILE_NAMES=names, $
INFO_TYPE=info_type, $
INSTR=instr, $
LEVEL=level, $
MODE=mode, $
OPTDESC=optdesc, $
PUBLIC_SITE=public, $
RESET_PATH=reset_path, $
RESET_QUERY=reset_query, $
SCIENCE=science, $
SC_ID=sc_id, $
SITE=site, $
TEAM_SITE=team, $
V_INFO=v_info, $
VERSION=version, $
_REF_EXTRA=extra
	Compile_Opt idl2
	On_Error, 2
	
	;Web properties
	IF N_Elements(extra) GT 0 THEN self.oWeb -> SetProperty, _STRICT_EXTRA=extra
	
;---------------------------------------------------------------------
; Default Path & Query ///////////////////////////////////////////////
;---------------------------------------------------------------------
	tf_reset_path = Keyword_Set(reset_path)
	IF N_Elements(files) GT 0 && files[0] NE '' THEN reset_query = 1B
	IF Keyword_Set(reset_query) THEN BEGIN
		new_sc      = ''
		new_instr   = ''
		new_mode    = ''
		new_level   = ''
		new_optdesc = ''
		new_version = ''
		new_ancprod = ''
		new_files   = ''
		new_tstart  = ''
		new_tend    = ''
	ENDIF ELSE BEGIN
		new_sc      = self.sc
		new_instr   = self.instr
		new_mode    = self.mode
		new_level   = self.level
		new_optdesc = self.optdesc
		new_version = self.version
		new_ancprod = self.ancprod
		new_files   = self.files
		self.oWeb  -> GetProperty, DATE_START=new_tstart, DATE_END=new_tend
	ENDELSE
	
;---------------------------------------------------------------------
; Dependencies ///////////////////////////////////////////////////////
;---------------------------------------------------------------------
	
	;Dependencies
	;   - Set ANCILLARY and HK automatically
	;   - If FILES is provided, other query parameters are ignored.
	IF N_Elements(ancprod) GT 0 && ancprod[0] NE '' THEN ancillary = 1B
	
	IF N_Elements(instr) GT 0 THEN BEGIN
		IF instr[0] EQ 'hk' THEN BEGIN
			hk        = 1B
			ancillary = 0B
			science   = 0B
		ENDIF ELSE IF Max(MrIsMember(['afg', 'dfg'], instr)) THEN BEGIN
			team   = 1B
			public = 0B
			site   = ''
		ENDIF
	ENDIF
	
	IF N_Elements(level) GT 0 && Max(MrIsMember(['l1a', 'l1b', 'sitl', 'ql', 'l2pre'], level) ) THEN BEGIN
		team   = 1B
		public = 0B
		site   = ''
	ENDIF

;---------------------------------------------------------------------
; Team or Public /////////////////////////////////////////////////////
;---------------------------------------------------------------------
	new_site = self.site
	IF N_Elements(site) GT 0 || N_Elements(team) GT 0 || N_Elements(public) GT 0 THEN BEGIN
		tf_site   = Keyword_set(site)
		tf_team   = Keyword_Set(team)
		tf_public = Keyword_Set(public)
		
		IF tf_site + tf_team + tf_public GT 1 THEN Message, 'SITE, TEAM, and PUBLIC are mutually exclusive.'
		
		;Default to public site
		CASE 1 OF
			tf_site: BEGIN
				IF ~StRegEx(site, '^(sitl|public|team)$', /BOOLEAN, /FOLD_CASE) $
					THEN Message, 'Invalid value for SITE.'
				new_site = site EQ 'public' ? site : 'sitl'
			ENDCASE
			tf_team:   new_site = 'sitl'
			tf_public: new_site = 'public'
			ELSE:      new_site = 'public'
		ENDCASE
	ENDIF ELSE IF tf_reset_path THEN BEGIN
		new_site = 'public'
	ENDIF

;---------------------------------------------------------------------
; Names | Info | Version | Download //////////////////////////////////
;---------------------------------------------------------------------
	new_info_type = self.info_type
	IF N_Elements(info_type) GT 0 || N_Elements(info)     GT 0 || N_Elements(names) GT 0 || $
	   N_Elements(v_info)    GT 0 || N_Elements(download) GT 0 $
	THEN BEGIN
		tf_type     = Keyword_Set(type)
		tf_info     = Keyword_Set(info)
		tf_names    = Keyword_Set(names)
		tf_v_info   = Keyword_Set(v_info)
		tf_download = Keyword_Set(download)
		IF tf_type + tf_names + tf_info + tf_v_info + tf_download GT 1 $
			THEN Message, 'Only one of DOWNLOAD, FILE_NAMES, FILE_INFO, INFO_TYPE, and V_INFO may be selected.'
	
		;Default to downloading information
		CASE 1 OF
			tf_type: BEGIN
				If ~StRegEx(info_type, '^(file_names|file_info|version_info|download)$', /BOOLEAN, /FOLD_CASE) $
					THEN Message, 'Invalid value for INFO_TYPE.'
				self.info_type = StrLowCase(info_type)
			ENDCASE
			tf_names:    new_info_type = 'file_names'
			tf_info:     new_info_type = 'file_info'
			tf_v_info:   new_info_type = 'version_info'
			tf_download: new_info_type = 'download'
			ELSE:        new_info_type = 'download'
		ENDCASE
	ENDIF ELSE IF tf_reset_path THEN BEGIN
		new_info_type = 'download'
	ENDIF

;---------------------------------------------------------------------
; Ancillary | HK | Science ///////////////////////////////////////////
;---------------------------------------------------------------------
	new_data_type = self.data_type
	IF N_Elements(date_type) GT 0 || N_Elements(ancillary) GT 0 || $
	   N_Elements(science)   GT 0 || N_Elements(hk)        GT 0 $
	THEN BEGIN
		tf_data_type = Keyword_Set(data_type)
		tf_science   = Keyword_Set(science)
		tf_anc       = Keyword_Set(ancillary)
		tf_hk        = Keyword_Set(hk)
		IF tf_data_type && tf_science && tf_anc && tf_hk THEN Message, 'ANCILLARY, DATA_TYPE, HK, and SCIENCE are mutually exclusive.'
		
		;Default to science
		CASE 1 OF
			tf_data_type: BEGIN
				IF ~StRegEx(data_type, '^(ancillary|hk|science)$', /BOOLEAN, /FOLD_CASE) $
					THEN Message, 'Invalid value for DATA_TYPE'
				self.data_type = StrLowCase(data_type)
			ENDCASE
			tf_anc:     new_data_type = 'ancillary'
			tf_hk:      new_data_type = 'hk'
			tf_science: new_data_type = 'science'
			ELSE:       new_data_type = 'science'
		ENDCASE
	ENDIF ELSE IF tf_reset_path THEN BEGIN
		new_data_type = 'science'
	ENDIF

;---------------------------------------------------------------------
; Query Parameters ///////////////////////////////////////////////////
;---------------------------------------------------------------------
	;FILES is not used with anything else
	IF N_Elements(files) GT 0 THEN BEGIn
		new_files = files
	
	;QUERY
	ENDIF ELSE BEGIN
		;SC
		IF N_Elements(sc_id) GT 0 THEN BEGIN
			IF ~Array_Equal( MrIsMember(['', 'mms1', 'mms2', 'mms3', 'mms4'], sc_id), 1 ) $
				THEN Message, 'Invalid value(s) for SC_ID.'
			new_sc = StrJoin(sc_id, ',')
		ENDIF
		
		;INSTR
		IF N_Elements(instr) GT 0 THEN BEGIN
			IF ~Array_Equal( MrIsMember(['', 'afg', 'aspoc', 'dfg', 'dsp', 'edi', 'edp', $
			                             'epd-eis', 'feeps', 'fgm', 'fpi', 'fsm', 'hpca', $
			                             'hk', 'mec', 'scm'], instr), 1 ) $
				THEN Message, 'Invalid value(s) for INSTR.'
			new_instr = StrJoin(instr, ',')
		ENDIF
		
		;MODE
		IF N_Elements(mode) GT 0 THEN BEGIN
			IF ~Array_Equal( MrIsMember(['', 'slow', 'fast', 'srvy', 'brst'], mode), 1 ) $
				THEN Message, 'Invalid value(s) for MODE.'
			new_mode = StrJoin(mode, ',')
		ENDIF
		
		;LEVEL
		IF N_Elements(level) GT 0 THEN BEGIN
			IF ~Array_Equal( MrIsMember(['', 'l1a', 'l1b', 'ql', 'sitl', 'l2pre', 'l2', 'l2plus', 'l3', 'acr'], level), 1 ) $
				THEN Message, 'Invalid value(s) for LEVEL.'
			new_level = StrJoin(level, ',')
		ENDIF
		
		;OTHER
		IF N_Elements(optdesc)    GT 0 THEN new_optdesc = StrJoin(optdesc, ',')
		IF N_Elements(version)    GT 0 THEN new_version = version
		IF N_Elements(ancprod)    GT 0 THEN new_ancprod = StrJoin(ancprod, ',')
		IF N_Elements(files)      GT 0 THEN new_files   = StrJoin(files, ',')
		IF N_Elements(date_start) GT 0 THEN new_tstart  = date_start
		IF N_Elements(date_end)   GT 0 THEN new_tend    = date_end
	ENDELSE

;---------------------------------------------------------------------
; Scheme, Host, Path & Query Strings /////////////////////////////////
;---------------------------------------------------------------------
	
	;Query
	uri_query = ''
	IF new_sc      NE '' THEN uri_query += 'sc_id='          + new_sc      + '&'
	IF new_instr   NE '' THEN uri_query += 'instrument_id='  + new_instr   + '&'
	IF new_mode    NE '' THEN uri_query += 'data_rate_mode=' + new_mode    + '&'
	IF new_level   NE '' THEN uri_query += 'data_level='     + new_level   + '&'
	IF new_optdesc NE '' THEN uri_query += 'descriptor='     + new_optdesc + '&'
	IF new_version NE '' THEN uri_query += 'version='        + new_version + '&'
	IF new_ancprod NE '' THEN uri_query += 'product='        + new_ancprod + '&'
	IF new_files   NE '' THEN uri_query += 'file='           + new_files   + '&'
	
	;The query requires two dates from two different days
	;   - If we want data from a single day, THEN some fudging needs to be done
	;   - Also, IF TSTART=2015-288 and TEND=2015-289, THEN the ancillary file
	;     returned will be 2015-288 to 2015-289. This file, however, contains
	;     only 2 minutes OF data from 288. There is another file from 2015-287
	;     to 2015-288. In order to find it, TSTART must be 2015-287. Thus, here,
	;     we bump TSTART and TEND by one day FOR ancillary files
	date_start = new_tstart
	date_end   = new_tend
	IF date_start NE '' && date_end NE '' THEN BEGIN
		;Extract the date part
		date_start = strmid(new_tstart, 0, 10)
		date_end   = strmid(new_tend,   0, 10)
		time_end   = strmid(new_tend,  11,  8)
		IF time_end EQ '' THEN time_end = '00:00:00'

		;Bump down DATE_START IF:
		;   - ANCILLARY data is being downloaded
		IF new_data_type EQ 'ancillary' THEN BEGIN
			;Parse the date
			month = fix(strmid(date_start, 5, 2))
			day   = fix(strmid(date_start, 8, 2))
			year  = fix(strmid(date_start, 0, 4))
			
			;Bump the END date up by one day
			cdf_tt2000, tt_end, year, month, day-1, /COMPUTE_EPOCH
			date_start = strmid(cdf_encode_tt2000(tt_end), 0, 10)
		ENDIF

		;Bump up DATE_END IF:
		;   - DATE_START = DATE_END
		;   - ANCILLARY data is being downloaded
		IF date_start EQ date_end || time_end NE '00:00:00' || new_data_type EQ 'ancillary' THEN BEGIN
			;Parse the date
			month = fix(strmid(date_end, 5, 2))
			day   = fix(strmid(date_end, 8, 2))
			year  = fix(strmid(date_end, 0, 4))
			
			;Bump the END date up by one day
			cdf_tt2000, tt_end, year, month, day, 24, 0, 0, /COMPUTE_EPOCH
			date_end = strmid(cdf_encode_tt2000(tt_end), 0, 10)
		ENDIF
	ENDIF

	;TSTART
	IF new_tstart NE '' THEN BEGIN
		;2015-05-14T12:00:00
		IF MrTokens_IsMatch(new_tstart, '%Y-%M-%dT%H:%m:%S') THEN BEGIN
			uri_query  += 'start_date=' + date_start + '&'

		;2015-05-14
		ENDIF ELSE IF MrTokens_IsMatch(new_tstart, '%Y-%M-%d') THEN BEGIN
			uri_query  += 'start_date=' + strmid(date_start, 0, 10) + '&'
			new_tstart += 'T00:00:00'
		;Unrecognized format
		ENDIF ELSE BEGIN
			message, 'Invalid format for TSTART: "' + tstart + '".'
		ENDELSE
	ENDIF

	;TEND
	IF new_tend NE '' THEN BEGIN
		;2015-05-14T12:00:00
		IF MrTokens_IsMatch(new_tend, '%Y-%M-%dT%H:%m:%S') THEN BEGIN
			uri_query += 'end_date=' + strmid(date_end, 0, 10) + '&'

		;2015-05-14
		ENDIF ELSE IF MrTokens_IsMatch(new_tend, '%Y-%M-%d') THEN BEGIN
			uri_query += 'end_date=' + date_end + '&'
			new_tend  += 'T24:00:00'
		;Unrecognized format
		ENDIF ELSE BEGIN
			message, 'Invalid format for TEND: "' + tend + '".'
		ENDELSE
	ENDIF

	;Trim the trailing "&"
	uri_query = StrMid(uri_query, 0, StrLen(uri_query)-1)

;---------------------------------------------------------------------
; Set URI Properties /////////////////////////////////////////////////
;---------------------------------------------------------------------
	;Path
	;   - Prevent the query string from being empty. Tries to download everything.
	;   - Stop at "site". Further does not exist.
	uri_scheme = 'https'
	uri_host   = 'lasp.colorado.edu'
	IF uri_query EQ '' $
		THEN uri_path = '/mms/sdc/' + new_site $
		ELSE uri_path = '/mms/sdc/' + new_site + '/files/api/v1/' + new_info_type + '/' + new_data_type
	
	;Previous URI
	self.oWeb -> GetProperty, SCHEME = old_scheme, $
	                          HOST   = old_host, $
	                          PATH   = old_path, $
	                          QUERY  = old_query
	
	;Set the URI only if it has changed
	;   - It is possible for files to be found locally in DROPBOX_ROOT,
	;     MIRROR_ROOT, or LOCAL_ROOT that do not exist yet at the official
	;     SDC. This will prevent some repeated warnings and errors concerning
	;     attempts to access a link with Code 204: "NO CONTENT" (and the like).
	IF uri_scheme NE old_scheme || uri_host  NE old_host || $
	   uri_path   NE old_path   || uri_query NE old_query $
	THEN BEGIN
		self.oWeb -> SetURI, '', success, $
		                     HOST   = uri_host, $
		                     PATH   = uri_path, $
		                     QUERY  = uri_query, $
		                     SCHEME = uri_scheme
	ENDIF

	;Set properties
	self.sc        = Temporary(new_sc)
	self.instr     = Temporary(new_instr)
	self.mode      = Temporary(new_mode)
	self.level     = Temporary(new_level)
	self.optdesc   = Temporary(new_optdesc)
	self.version   = Temporary(new_version)
	self.ancprod   = Temporary(new_ancprod)
	self.files     = Temporary(new_files)
	self.info_type = Temporary(new_info_type)
	self.data_type = Temporary(new_data_type)
	self.site      = Temporary(new_site)
	IF new_tstart NE '' || new_tend NE '' $
		THEN self.oWeb -> SetProperty, DATE_START=new_tstart, DATE_END=new_tend
END


;+
;   Get file versions.
;
; :Params:
;       MAX_VERSION:    out, optional, type=strarr
;                       Maximum version number for each MMS data file type.
;       FILE_TYPE:      out, optional, type=strarr
;                       Types of files associated with `MAX_VERSION`
;
; :Keywords:
;       COUNT:          out, optional, type=integer
;                       Number of elements returned.
;-
PRO MrMMS_SDC_API::VersionInfo, max_version, file_type, $
COUNT=nVersions
	Compile_Opt idl2

	;Catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		
		;Delete the tar file
		IF N_Elements(file) GT 0 && File_Test(file, /REGULAR) $
			THEN File_Delete, file
		
		;Delete the JSON file
		IF N_Elements(jsonfile) GT 0 && File_Test(jsonfile, /REGULAR) $
			THEN File_Delete, jsonfile
		
		;Print the error
		;   - ::CD will call ::GetResponseHeader which will intentionally
		;     throw an error. Must report this error first.
		MrPrintF, 'LogErr'
		
		;Return to the previous setting
		IF N_Elements(info_type) && info_type NE 'version_info' $
			THEN self -> SetProperty, INFO_TYPE=info_type
		
		RETURN
	ENDIF
	
	;Offline mode
	IF self.offline THEN Message, 'OFFLINE=1: Cannot obtain version info.'
	
	;CD to the Versions service
	info_type = self.info_type
	IF info_type NE 'version_info' THEN self -> SetProperty, /V_INFO
	
	;Download the JSON file
	file = self -> Download()
	
	;Return to the previous setting
	IF info_type NE 'version_info' THEN self -> SetProperty, INFO_TYPE=info_type

	;Replace the ".tar.gz" extension with ".json"
	jsonfile  = StRegEx(file, '(.*)\.tar\.gz', /SUBEXP, /EXTRACT)
	jsonfile  = jsonfile[1]
	jsonfile += '.json'
	
	;Rename the file
	File_Move, file, jsonfile
	
	;Read the file
	jsonstr = ''
	OpenR, lun, jsonfile, /GET_LUN
	ReadF, lun, jsonstr, FORMAT='(a)'
	Free_LUN, lun
	
	;Delte the file
	File_Delete, jsonfile
	
	;Parse string
	json = JSON_Parse(jsonstr)
	
	;Parse further
	nVersions   = N_Elements(json['versions'])
	file_type   = StrArr(nVersions)
	max_version = StrArr(nVersions)
	FOREACH value, json['versions'], idx DO BEGIN
		;Get the file info
		file_hash = json['versions', idx]
		
		;Extract the file info
		file_type[idx]   = file_hash['file_type']
		max_version[idx] = file_hash['max_version']
	ENDFOREACH
END


;+
;   Cleanup after the object is destroyed. This will destroy the widget, if it exists.
;-
PRO MrMMS_SDC_API::Cleanup
	Compile_Opt idl2
	On_Error, 2

	;Objects
	Obj_Destroy, self.oWeb
END


;+
;   Initialization method.
;
; :Keywords:
;       RESET_PATH:     in, optional, type=boolean, default=0
;                       If set, the path will be returned to the default path. This is
;                           equivalent to setting `PUBLIC`, `DOWNLOAD`, and `SCIENCE`.
;       RESET_QUERY:    in, optional, type=boolean, default=0
;                       If set, no information from previous queries will be retained.
;                           The default is to use values from the previous query unless
;                           otherwise provided via the `SC`, `INSTR`, `MODE`, `LEVEL`,
;                           `OPTDESC`, `TSTART`, `TEND`, or `ANC_PRODUCT` keywords.
;                           This keyword takes precedence over the aforementioned
;                           keywords.
;
;       ANCILLARY:      in, optional, type=boolean, default=0
;                       If set, the request will be for anciallary data products.
;                           Cannot be used with the `HK`, `SCIENCE`, or `VERSION`.
;       HK:             in, optional, type=boolean, default=0
;                       If set, the request will be for house keeping data products.
;                           Cannot be used with the `ANCILLARY`, `SCIENCE`, or `VERSION`.
;       SCIENCE:        in, optional, type=boolean, default=0
;                       If set, the request will be for science data products.
;                           Cannot be used with the `ANCILLARY`, or `HK`. If
;                           None of the three keywords are set, the current path will
;                           be checked for one of those options. If still nothing is found
;                           then SCIENCE will be the default.
;
;       DOWNLOAD:       in, optional, type=boolean, default=0
;                       If set, the URL will request that files be downloaded. Cannot
;                           be used with `INFO`, `NAMES`, or `V_INFO`. If none of the
;                           four keywords are set, the current URI will be examined for
;                           a choice. If still nothing has been selected, DOWNLOAD will
;                           be the default.
;       FILE_INFO:      in, optional, type=boolean, default=0
;                       If set, the URL will request file information. Cannot be used
;                           with `DOWNLOAD`, `NAMES`, or `V_INFO`
;       FILE_NAMES:     in, optional, type=boolean, default=0
;                       If set, the URL will requets file names. Cannot be used with
;                           `DOWNLOAD`, `INFO`, or `V_INFO`
;       V_INFO:         in, optional, type=boolean, default=0
;                       If set, the request will be for file version information.
;                           Cannot be used with the `ANCILLARY` or `HK`, or with
;                           `DOWNLOAD`, `INFO`, or `NAMES`.
;       FILES:          in, optional, type=string/strarr, default=''
;                       Names of the files to be downloaded. If set, then `SCIENCE`, 
;                           `DOWNLOAD`, and `RESET_QUERY` will be set.
;
;       PUBLIC:         in, optional, type=boolean, default=0
;                       If set, the public site will be queried (no password required,
;                           contains data of quality L2 and above). Cannot be used
;                           with `TEAM`. If neither keywords are set, the current URL
;                           is examined. If still neither are chosen, PUBLIC is the default.
;       TEAM:           in, optional, type=string/strarr
;                       If set, the public site will be queried (password required,
;                           contains data of quality L1A and above). Cannot be used
;                           with `TEAM`.
;
;       ANC_PRODUCT:    in, optional, type=string/strarr
;                       Ancillary data product names. Should only be used with `ANCILLARY`.
;       INSTR:          in, optional, type=string/strarr
;                       Instrument ID.
;       LEVEL:          in, optional, type=string/strarr
;                       Data quality level.
;       MODE:           in, optional, type=string/strarr
;                       Data rate mode.
;       OPTDESC:        in, optional, type=string/strarr
;                       Optional file name descriptor.
;       SC_ID:          in, optional, type=string/strarr
;                       Spacecraft ID. Options are: {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       TEND:           in, optional, type=string
;                       End time of data interval. Must be formatted as YYYY-MM-DD or
;                           YYYY-MM-DDThh:mm:ss
;       TSTART:         in, optional, type=string/strarr
;                       Start time of data interval. Must be formatted as YYYY-MM-DD or
;                           YYYY-MM-DDThh:mm:ss
;       VERSION:        in, optional, type=string/strarr
;                       File name version number, formatted as X.Y.Z, where X, Y, and Z
;                           are integers.
;
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by MrWebURI::Init.
;-
FUNCTION MrMMS_SDC_API::init, uri, $
ANC_PRODUCT=ancprod, $
ANCILLARY=ancillary, $
DATE_END=date_end, $
DATE_START=date_start, $
DOWNLOAD=download, $
HK=hk, $
FILES=files, $
FILE_INFO=info, $
FILE_NAMES=names, $
INSTR=instr, $
LEVEL=level, $
LOCAL_ROOT=local_root, $
MODE=mode, $
OPTDESC=optdesc, $
PUBLIC_SITE=public, $
SCIENCE=science, $
SC_ID=sc_id, $
TEAM_SITE=team, $
V_INFO=v_info, $
VERSION=version, $
_REF_EXTRA = extra
	Compile_Opt idl2

	;Catch errors
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, 0
	ENDIF
	
	;Set the local data root
	IF n_elements(local_root) EQ 0 THEN BEGIN
		local_root = filepath('', ROOT_DIR=file_search('~', /TEST_DIRECTORY, /EXPAND_TILDE), $
		                          SUBDIRECTORY=['MrWebData', 'mms'] )
		MrPrintF, 'LogWarn', 'Setting local data root to: "' + local_root + '".'
	ENDIF

	;Superclass INIT.
	self.oWeb = MrWebURI( DATE_START   = date_start, $
	                      DATE_END     = date_end, $
	                      LOCAL_ROOT   = local_root, $
	                      TPATTERN     = '%Y-%M-%dT%H:%m:%S', $
	                     _STRICT_EXTRA = extra )
	IF ~Obj_Valid(self.oWeb) THEN Message, 'Unable to initialize oWeb property.'
	
	;Defaults
	tf_info     = Keyword_Set(info)
	tf_names    = Keyword_Set(names)
	tf_v_info   = Keyword_Set(v_info)
	tf_download = Keyword_Set(download)
	IF ~tf_info && ~tf_names && ~tf_v_info && ~tf_download THEN tf_download = 1B
	
	tf_team   = Keyword_Set(team)
	tf_public = Keyword_Set(public)
	IF ~tf_team && ~tf_public THEN tf_public = 1B
	
	tf_ancillary = Keyword_Set(ancillary)
	tf_science   = Keyword_Set(science)
	tf_hk        = Keyword_Set(hk)
	IF tf_ancillary + tf_science + tf_hk EQ 0 THEN tf_science = 1B
	
	;Allowed query values
;	*self.field_values = { sc_id:          '', $
;	                       instrument_id:  '', $
;	                       data_rate_mode: '', $
;	                       data_level:     '', $
;	                       descriptor:     '', $
;	                       version:        '', $
;	                       product:        '', $
;	                       start_date:     '', $
;	                       end_date:       '', $
;	                       file:           '' }

	;Create the SDC GUI
	IF N_Elements(uri) GT 0 THEN self.oWeb -> CD, uri
	
	;Set properties via the BuildURI method
	self -> SetProperty, ANC_PRODUCT = ancprod, $
	                     ANCILLARY   = tf_ancillary, $
	                     DOWNLOAD    = tf_download, $
	                     HK          = hk, $
	                     FILES       = files, $
	                     FILE_INFO   = tf_info, $
	                     FILE_NAMES  = tf_names, $
	                     INSTR       = instr, $
	                     LEVEL       = level, $
	                     MODE        = mode, $
	                     OPTDESC     = optdesc, $
	                     PUBLIC_SITE = tf_public, $
	                     SCIENCE     = tf_science, $
	                     SC_ID       = sc_id, $
	                     TEAM_SITE   = tf_team, $
	                     V_INFO      = tf_v_info, $
	                     VERSION     = version

	RETURN, 1
END


;+
;   Object definition.
;
; :Params:
;       CLASS:          out, optional, type=structure
;                       The class definition structure.
;
; :Fields:
;       TEND:               End time of data interval to be queried.
;       TSTART:             Start time of data interval to be queried.
;-
PRO MrMMS_SDC_API__define, class

	class = { MrMMS_SDC_API, $
	          Inherits IDL_Object, $
	          oWeb:      Obj_New(), $
	          site:      '', $
	          info_type: '', $
	          data_type: '', $
	          sc:        '', $
	          instr:     '', $
	          mode:      '', $
	          level:     '', $
	          optdesc:   '', $
	          version:   '', $
	          ancprod:   '', $
	          files:     '' $
	        }
END