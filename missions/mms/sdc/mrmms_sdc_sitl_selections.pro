; docformat = 'rst'
;
; NAME:
;       MrMMS_SDC_SITL_Selections
;
;*****************************************************************************************
;   Copyright (c) 2019, Matthew Argall                                                   ;
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
;   Return information on the SITL selections within a given interval.
;
; :Categories:
;   MMS
;
; :Params:
;       TSTART:         in, optional, type=string, default=MrVar_Get_TRange()
;                       Start time of interval during which SITL selections are desired,
;                           formatted as YYYY-MM-DDTHH:MM:SS
;       TEND:           in, optional, type=string, default=MrVar_Get_TRange()
;                       End time of interval during which SITL selections are desired,
;                           formatted as YYYY-MM-DDTHH:MM:SS
;
; :Keywords:
;       COMBINE:        in, optional, type=boolean, default=0
;                       If set, contiguous burst intervals are merged into a single entry.
;       COMPLETE:       in, optional, type=boolean, default=0
;                       If set, only completed, downlinked selections are returned.
;       COUNT:          out, optional, type=int
;                       Number of SITL selections returned.
;       FILENAME:       in, optional, type=string, default=''
;                       Name of a file to which the selections are saved.
;       PRINT:          in, optional, type=boolean, default=0
;                       If set, results are printed to the terminal window.
;       SEARCH_STRING:  in, optional, type=string, default=''
;                       Only those entries that match this string are returned. Matching
;                           is done via StRegEx.
;       TEAM:           in, optional, type=boolean, default=0
;                       If set, the query will be made to the team site. User is prompted
;                           for username and password.
;
; :Returns:
;       DATA:           out, required, structure
;                       Data structure containing sitl selections. If `COUNT`=0,
;                           !Null is returned
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
;       2019/01/31  -   Written by Matthew Argall
;-
FUNCTION MrMMS_SDC_SITL_Selections, tstart, tend, $
COMBINE=combine, $
COMPLETE=complete, $
COUNT=count, $
FILENAME=filename, $
PRINT=print, $
SEARCH_STRING=search_string, $
TEAM=team
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF Obj_Valid(oURL) THEN Obj_Destroy, oURL
		IF N_Elements(file) GT 0 THEN File_Delete, file, /ALLOW_NONEXISTENT
		IF N_Elements(lun) GT 0 && (FStat(lun)).open THEN Close, lun
		MrPrintF, 'LogErr'
		count = 0
		RETURN, !Null
	ENDIF
	
	tf_complete = Keyword_Set(complete)
	tf_combine  = Keyword_Set(combine)
	tf_print    = Keyword_Set(print)
	tf_team     = Keyword_Set(team)
	tf_events   = Keyword_Set(team)
	IF N_Elements(filename)      EQ 0 THEN filename      = ''
	IF N_Elements(search_string) EQ 0 THEN search_string = ''
	
	CASE N_Elements(tstart) OF
		2: BEGIN
			t0 = tstart[0]
			t1 = tstart[1]
		ENDCASE
		1: BEGIN
			t0 = tstart
			t1 = tend
		ENDCASE
		0: BEGIN
			trange = MrVar_GetTRange()
			t0     = trange[0]
			t1     = trange[1]
		ENDCASE
		ELSE: Message, 'Incorrect inputs.'
	ENDCASE

	;Burst Segment Status
;	basic  = 'https://lasp.colorado.edu/mms/sdc/sitl/latis/dap/mms_burst_data_segment.csv?CREATETIME>=01-Jun-2017'
;	public = 'https://lasp.colorado.edu/mms/sdc/public/service/latis/mms_burst_data_segment.csv?CREATETIME>=01-Jun-2017'

	;Mission Events
;	public = 'https://lasp.colorado.edu/mms/sdc/public/service/latis/mms_events_view.csv?start_time_utc>2017-06-09'
	
	;Convert times to TT2000 times
	t_1958 = MrCDF_Epoch_Compute(1958, 01, 01, 0, 0, 0, /TT2000)
	t0     = Long64((MrCDF_Epoch_Parse(t0, /TO_TT2000, PATTERN='%Y-%M-%dT%H:%m:%S') - t_1958) / 1d9)
	t1     = Long64((MrCDF_Epoch_Parse(t1, /TO_TT2000, PATTERN='%Y-%M-%dT%H:%m:%S') - t_1958) / 1d9)
	
;-------------------------------------------
; Download Data ////////////////////////////
;-------------------------------------------
	IF tf_team $
		THEN oURL = MrMMS_SDC_LogIn() $
		ELSE oURL = IDLnetURL()
	url_path = '/mms/sdc/' + (tf_team ? 'sitl/latis/dap/' : 'public/service/latis/')
	url_file = 'mms_burst_data_segment.csv'
	oURL -> SetProperty, URL_HOSTNAME = 'lasp.colorado.edu', $
	                     URL_PATH     = url_path + url_file, $
	                     URL_QUERY    = 'TAISTARTTIME>='+String(t0, FORMAT='(i0)')+'&TAIENDTIME<='+String(t1, FORMAT='(i0)'), $
	                     URL_SCHEME   = 'https'
	file = oURL -> Get(FILENAME=url_file)
	Obj_Destroy, oURL
	
;-------------------------------------------
; Read Downloaded File /////////////////////
;-------------------------------------------
	
	;No selections in time interval
	IF File_Lines(file) LE 2 THEN BEGIN
		File_Delete, file
		count = 0
		RETURN, {}
	ENDIF
	
	;Read header
	header = ''
	OpenR, lun, file, /GET_LUN
	ReadF, lun, header, FORMAT='(a)' 
	Close, lun
	
	;Read data
	temp = Read_CSV( file, $
	                 COUNT          = count, $
	                 N_TABLE_HEADER = 0, $
	                 RECORD_START   = 1, $
	                 TYPE           = ['LONG', 'LONG64', 'LONG64', 'STRING', 'FLOAT', 'BYTE', 'BYTE', 'STRING', $
	                                   'LONG', 'STRING', 'DATETIME', 'DATETIME', 'LONG', 'LONG', 'LONG', 'LONG', $
	                                   'LONG', 'LONG', 'LONG', 'LONG', 'LONG', 'LONG', 'LONG', 'LONG', 'STRING'] )
	
	;Delete the file
	File_Delete, file, /ALLOW_NONEXISTENT
	
	;Update structure fields
	data = {}
	header = StrSplit(header, ',', /EXTRACT, COUNT=nHeader)
	FOR i = 0, nHeader - 1 DO data = create_struct(data, (StrSplit(header[i], ' ', /EXTRACT))[0], temp.(i))
	temp = !Null
	
	;Encode times as strings
	data = Create_Struct( data, $
	                      'DT',     Float(data.taiendtime - data.taistarttime), $
	                      'TSTART', MrCDF_Epoch_Encode(Long64(data.taistarttime*1d9) + t_1958, EPOCH=3), $
	                      'TEND',   MrCDF_Epoch_Encode(Long64(data.taiendtime*1d9) + t_1958, EPOCH=3) )
	

;-------------------------------------------
; Select Downlinked Intervals //////////////
;-------------------------------------------
	IF count GT 0 && tf_complete THEN BEGIn
		idx  = Where( StRegEx(data.status, '(^COMPLETE|[^IN]+COMPLETE)', /BOOLEAN, /FOLD_CASE), count )
		temp = {}
		IF count GT 0 THEN BEGIN
			tags = Tag_Names(data)
			FOR i = 0, N_Elements(tags) - 1 $
				DO temp = Create_Struct(temp, tags[i], data.(i)[idx])
		ENDIF
		data = Temporary(temp)
	ENDIF

;-------------------------------------------
; Filter Discussion ////////////////////////
;-------------------------------------------
	IF search_string NE '' THEN BEGIN
		;Find matches
		idx = Where( StRegEx(data.discussion, search_string, /BOOLEAN, /FOLD_CASE), count )
		
		;Select elements
		temp = {}
		IF count GT 0 THEN BEGIN
			tags = Tag_Names(data)
			FOR i = 0, N_Elements(tags) - 1 $
				DO temp = Create_Struct(temp, tags[i], data.(i)[idx])
		ENDIF
		data = Temporary(temp)
	ENDIF

;-------------------------------------------
; Combine Contiguous Intervals /////////////
;-------------------------------------------
	IF count GT 0 && tf_combine THEN BEGIN
		;Contiguous intervals are separated by 10s
		dt    = [1000, data.taistarttime[1:-1] - data.taiendtime[0:-2]]
		temp  = data
		ntags = N_Tags(data)
		idx   = 0
		FOR i = 0, count - 1 DO BEGIN
			;Append the discussion
			IF dt[i] EQ 10 THEN BEGIN
				;Only unique discussions
				IF ~StrMatch(temp.discussion[idx-1], '*'+data.discussion[i]+'*') $
					THEN temp.discussion[idx-1] += ('; ' + data.discussion[i])
				
				;The last interval
				IF i EQ count - 1 THEN BEGIN
					temp.taiendtime[idx-1] = data.taiendtime[i]
					temp.tend[idx-1]       = data.tend[i]
					temp.dt[idx-1]         = temp.taiendtime[idx-1] - temp.taistarttime[idx-1]
				ENDIF
				
				;Jump to next
				CONTINUE
			ENDIF
			
			;Update end time of previous entry
			IF idx GT 0 THEN BEGIN
				temp.taiendtime[idx-1] = data.taiendtime[i-1]
				temp.tend[idx-1]       = data.tend[i-1]
				temp.dt[idx-1]         = temp.taiendtime[idx-1] - temp.taistarttime[idx-1]
			ENDIF
			
			;Copy data from new interval
			FOR j = 0, ntags - 1 DO temp.(j)[idx] = data.(j)[i]
			
			;Next interval
			idx += 1
		ENDFOR
		
		data  = Temporary(temp)
		count = idx
		tags  = Tag_Names(data)
		temp  = {}
		FOR i = 0, ntags - 1 DO temp = Create_Struct(temp, tags[i], data.(i)[0:count-1])
		data = Temporary(temp)
	ENDIF

;-------------------------------------------
; Print to Terminal ////////////////////////
;-------------------------------------------
	IF tf_print THEN BEGIN
		Print, 'TSTART', 'TEND', 'DT', 'FOM', 'STATUS', 'DISCUSSION', $
		       FORMAT='(a6, 16x, a4, 18x, a2, 5x, a3, 2x, a6, 15x, a10)'
		FOR i = 0, count - 1 DO BEGIN
			Print, StrMid(data.tstart[i], 0, 10), StrMid(data.tstart[i], 11, 8), $
			       StrMid(data.tend[i], 0, 10),   StrMid(data.tend[i], 11, 8), $
			       data.dt[i], data.fom[i], data.status[i], data.discussion[i], $
			       FORMAT='(2(a10, 2x, a8, 2x), i5, 2x, i3, 2x, a-19, 2x, a0)'
		ENDFOR
	ENDIF

;-------------------------------------------
; Write to File ////////////////////////////
;-------------------------------------------
	IF filename NE '' THEN BEGIN
		OpenW, lun, filename, /GET_LUN
		PrintF, lun, 'TSTART', 'TEND', 'DT', 'FOM', 'STATUS', 'DISCUSSION', $
		             FORMAT='(a6, 16x, a4, 18x, a2, 5x, a3, 2x, a6, 15x, a10)'
		FOR i = 0, count - 1 DO BEGIN
			PrintF, lun, StrMid(data.tstart[i], 0, 10), StrMid(data.tstart[i], 11, 8), $
			       StrMid(data.tend[i], 0, 10),   StrMid(data.tend[i], 11, 8), $
			       data.dt[i], data.fom[i], data.status[i], data.discussion[i], $
			       FORMAT='(2(a10, 2x, a8, 2x), i5, 2x, i3, 2x, a-19, 2x, a0)'
		ENDFOR
		Close, lun
	ENDIF
	
	RETURN, data
END