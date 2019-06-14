; docformat = 'rst'
;
; NAME:
;   MrVar_TTicks
;
;*****************************************************************************************
;   Copyright (c) 2018, Matthew Argall                                                   ;
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
; PURPOSE
;+
;   .
;
; :Categories:
;   Plot Utilities
;
; :Params:
;       TRANGE:         in, required, type=strarr(2)
;                       The begin and end time, formatted as ISO-8601 strings
;
; :Keywords:
;       LABEL:          out, optional, type=string
;                       Upper and lower axis labels that indicate the time units used.
;       XTICKS:         in, optional, type=integer, default=3
;                       Number of tick intervals to use (one less than the number of tick marks).
;       XTICKINTERVAL:  in, optional, type=double
;                       Interval between tick marks to use.
;
; :Returns:
;       OUT:            out, required, type=struct
;                       A structure that can be passed to the PLOT procedure via
;                           use of the _STRICT_EXTRA keyword. Structure tags include:
;                               XTICKNAME - Labels for each major tick mark
;                               XTICKV    - Values at which major tick marks are placed
;                               XTICKS    - Number of tick intervals used
;                               XMINOR    - Number of minor tick marks between major ticks
;                               XRANGE    - Axis range
;                               XSTYLE    - Axis style
;                               XTITLE    - Axis title
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
;       2018-01-22  -   Created by Matthew Argall
;       2018-09-22  -   Lower tick label was being duplicated when units were minutes.
;                           Lower tick labels were being placed on last of common
;                           labels instead of first. Fixed. - MRA
;-
FUNCTION mrvar_tticks, trange, $
LABEL=label, $
XTICKS=xticks, $
XTICKINTERVAL=xtickinterval
	Compile_Opt idl2
	On_Error, 2

;	COMMON MrVar_TTick_Common, ct0, ct1
	
	;Initialize
	IF N_Elements(trange) EQ 0 THEN BEGIN
		trange     = MrVar_GetTRange()
		trange_ssm = MrVar_GetTRange('SSM')
		IF Array_Equal(trange, ['', '']) THEN Message, 'No time range defined.'
	ENDIF ELSE BEGIN
		oT = MrTimeVar(trange)
		trange_ssm = oT['DATA', 'SSM']
	ENDELSE
	
	;Defaults
	IF N_Elements(xticks) EQ 0 THEN xticks = 3
	
	;Base time units in seconds
	day    = 86400D
	year   = day*365.25
	month  = day*30
	hour   = 3600.0
	minute = 60.0
	second = 1.0
	dtimes = [year, month, day, hour, minute, second, 1e-3, 1e-6, 1e-9]
	
	;2017 Jan 22 10:08:52.134
	units  = ['', 'year', 'date', 'month', 'day', 'hour', 'minute',  'second',   'milli', 'micro', 'nano']
	labels = ['', 'year', 'date', 'month', 'day', 'hhmm',   'hhmm', 'seconds', 'seconds',    'us',   'ns']
	pos    = [0,      0,      5,       5,     9,     12,       15,        17,        19,      23,     26]
	width  = [0,      4,      6,       3,     2,      4,        4,         2,         4,       3,      3]
	
	; 
	;   - DT        = time interval used for XTICKINTERVAL
	;   - UNIT      = index into DTIMES
	;   - UNIT2     = index into UNITS
	;   - INCREMENT = number of DTIMES between tickmarks
	;   - NMINOR    = number of minor tick marks between major tick marks
	;   - FLD1      = index into LABELS, POS, WIDTH; used for upper axis label
	;   - FLD2      = index into LABELS, POS, WIDTH; used for lower axis label
	setup = [ $
		{ dt: 0.0, unit: 0, unit2: 0, increment: 1000, nminor: 10, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:  500, nminor:  5, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:  200, nminor:  4, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:  100, nminor: 10, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:   50, nminor:  5, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:   20, nminor:  4, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:   10, nminor: 10, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:    5, nminor:  5, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:    2, nminor:  4, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 0, unit2: 0, increment:    1, nminor: 12, fld1: 1, fld2: 0 }, $    ;years
		{ dt: 0.0, unit: 1, unit2: 1, increment:    6, nminor:  6, fld1: 3, fld2: 1 }, $    ;months
		{ dt: 0.0, unit: 1, unit2: 1, increment:    4, nminor:  4, fld1: 3, fld2: 1 }, $    ;months
		{ dt: 0.0, unit: 1, unit2: 1, increment:    3, nminor:  3, fld1: 3, fld2: 1 }, $    ;months
		{ dt: 0.0, unit: 1, unit2: 1, increment:    2, nminor:  2, fld1: 3, fld2: 1 }, $    ;months
		{ dt: 0.0, unit: 1, unit2: 1, increment:    1, nminor:  2, fld1: 3, fld2: 1 }, $    ;months
		{ dt: 0.0, unit: 2, unit2: 1, increment:   20, nminor:  4, fld1: 2, fld2: 1 }, $    ;days
		{ dt: 0.0, unit: 2, unit2: 1, increment:   10, nminor: 10, fld1: 2, fld2: 1 }, $    ;days
		{ dt: 0.0, unit: 2, unit2: 3, increment:    5, nminor:  5, fld1: 4, fld2: 3 }, $    ;days
		{ dt: 0.0, unit: 2, unit2: 3, increment:    2, nminor:  4, fld1: 4, fld2: 3 }, $    ;days
		{ dt: 0.0, unit: 2, unit2: 3, increment:    1, nminor:  6, fld1: 4, fld2: 3 }, $    ;days
		{ dt: 0.0, unit: 3, unit2: 2, increment:   12, nminor:  4, fld1: 5, fld2: 2 }, $    ;hours
		{ dt: 0.0, unit: 3, unit2: 2, increment:    8, nminor:  4, fld1: 5, fld2: 2 }, $    ;hours
		{ dt: 0.0, unit: 3, unit2: 2, increment:    6, nminor:  6, fld1: 5, fld2: 2 }, $    ;hours
		{ dt: 0.0, unit: 3, unit2: 2, increment:    4, nminor:  4, fld1: 5, fld2: 2 }, $    ;hours
		{ dt: 0.0, unit: 3, unit2: 2, increment:    3, nminor:  6, fld1: 5, fld2: 2 }, $    ;hours
		{ dt: 0.0, unit: 3, unit2: 2, increment:    2, nminor:  4, fld1: 5, fld2: 2 }, $    ;hours
		{ dt: 0.0, unit: 3, unit2: 2, increment:    1, nminor:  6, fld1: 5, fld2: 2 }, $    ;hours
		{ dt: 0.0, unit: 4, unit2: 2, increment:   30, nminor:  6, fld1: 5, fld2: 2 }, $    ;minutes
		{ dt: 0.0, unit: 4, unit2: 2, increment:   20, nminor:  4, fld1: 5, fld2: 2 }, $    ;minutes
		{ dt: 0.0, unit: 4, unit2: 2, increment:   10, nminor: 10, fld1: 5, fld2: 2 }, $    ;minutes
		{ dt: 0.0, unit: 4, unit2: 2, increment:    5, nminor:  5, fld1: 5, fld2: 2 }, $    ;minutes
		{ dt: 0.0, unit: 4, unit2: 2, increment:    2, nminor:  4, fld1: 5, fld2: 2 }, $    ;minutes
		{ dt: 0.0, unit: 4, unit2: 2, increment:    1, nminor:  6, fld1: 5, fld2: 2 }, $    ;minutes
		{ dt: 0.0, unit: 5, unit2: 6, increment:   30, nminor:  6, fld1: 7, fld2: 5 }, $    ;seconds
		{ dt: 0.0, unit: 5, unit2: 6, increment:   20, nminor:  4, fld1: 7, fld2: 5 }, $    ;seconds
		{ dt: 0.0, unit: 5, unit2: 6, increment:   10, nminor: 10, fld1: 7, fld2: 5 }, $    ;seconds
		{ dt: 0.0, unit: 5, unit2: 6, increment:    5, nminor:  5, fld1: 7, fld2: 5 }, $    ;seconds
		{ dt: 0.0, unit: 5, unit2: 6, increment:    2, nminor:  4, fld1: 7, fld2: 5 }, $    ;seconds
		{ dt: 0.0, unit: 5, unit2: 6, increment:    1, nminor: 10, fld1: 7, fld2: 5 }, $    ;seconds
		{ dt: 0.0, unit: 6, unit2: 7, increment:  500, nminor:  5, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 6, unit2: 7, increment:  200, nminor:  4, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 6, unit2: 7, increment:  100, nminor: 10, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 6, unit2: 7, increment:   50, nminor:  5, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 6, unit2: 7, increment:   20, nminor:  4, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 6, unit2: 7, increment:   10, nminor: 10, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 6, unit2: 7, increment:    5, nminor:  5, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 6, unit2: 7, increment:    2, nminor:  4, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 6, unit2: 7, increment:    1, nminor: 10, fld1: 8, fld2: 7 }, $    ;milli
		{ dt: 0.0, unit: 7, unit2: 7, increment:  500, nminor:  5, fld1: 9, fld2: 7 }, $    ;micro
		{ dt: 0.0, unit: 7, unit2: 7, increment:  200, nminor:  4, fld1: 9, fld2: 7 }, $    ;micro
		{ dt: 0.0, unit: 7, unit2: 7, increment:  100, nminor: 10, fld1: 9, fld2: 7 }, $    ;micro
		{ dt: 0.0, unit: 7, unit2: 7, increment:   50, nminor:  5, fld1: 9, fld2: 7 }, $    ;micro
		{ dt: 0.0, unit: 7, unit2: 7, increment:   20, nminor:  4, fld1: 9, fld2: 7 }, $    ;micro
		{ dt: 0.0, unit: 7, unit2: 7, increment:   10, nminor: 10, fld1: 9, fld2: 7 }, $    ;micro
		{ dt: 0.0, unit: 7, unit2: 7, increment:    5, nminor:  4, fld1: 9, fld2: 7 }, $    ;micro
		{ dt: 0.0, unit: 7, unit2: 7, increment:    2, nminor:  4, fld1: 9, fld2: 7 }, $    ;micro
		{ dt: 0.0, unit: 7, unit2: 7, increment:    1, nminor: 10, fld1: 9, fld2: 7 } ]     ;micro
	
	;Base increment for each level
	setup.dt = dtimes[setup.unit] * setup.increment
	
	;Pick the time interval closest to the desired tick spacing
	;   - Intervals in SETUP are ordered largest to smallest
	tickinterval = (trange_ssm[1] - trange_ssm[0]) / xticks
	iTick        = Where(setup.dt LE tickinterval, nTick)
	IF nTick GT 0 $
		THEN tickinterval = setup[iTick[0]].dt $
		ELSE Message, 'Time range too small.'
	
	;Create a series of tickmark labels
	;   - Round the initial time down to the nearest multiple of the interval
	;   - IDL allows a maximum of 60 tickmarks
	t0    = trange_ssm[0] - (trange_ssm[0] MOD tickinterval)
	ticks = t0 + FIndGen(60.0) * tickinterval
	iKeep = Where(ticks GE trange_ssm[0] AND ticks LE trange_ssm[1], nKeep)
	ticks = ticks[iKeep]
	oT    = MrTimeVar(ticks, 'SSM', T_REF=trange[0])
	times = oT -> GetData('CUSTOM', TOKEN_FMT='%Y %c %d %H%m:%S.%1%2%3')
	
	MrTimeParser_Breakdown, times, '%Y %c %d %H%m:%S.%1%2%3', $
	                        YEAR   = year, $
	                        MONTH  = month, $
	                        DAY    = day, $
	                        HOUR   = hour, $
	                        MINUTE = minute, $
	                        SECOND = second, $
	                        MILLI  = milli, $
	                        MICRO  = micro, $
	                        NANO   = nano
	
	;Side label
	s1 = StrMid(times, pos[setup[iTick[0]].fld1], width[setup[iTick[0]].fld1])
	s2 = StrMid(times, pos[setup[iTick[0]].fld2], width[setup[iTick[0]].fld2])
	
	;Do not duplicate lower label on adjacent tickmarks
	IF setup[iTick[0]].fld2 NE 0 THEN BEGIN
		CASE setup[iTick[0]].unit2 OF
			1: u2 = Fix(year)   - Fix(Shift(year,   1))
			2: u2 = Fix(month)  - Fix(Shift(month,  1))
			3: u2 = Fix(day)    - Fix(Shift(day,    1))
			4: u2 = Fix(hour)   - Fix(Shift(hour,   1))
			5: u2 = Fix(minute) - Fix(Shift(minute, 1))
			6: u2 = Fix(minute) - Fix(Shift(minute, 1))
			7: u2 = Fix(second) - Fix(Shift(second, 1))
			8: u2 = Fix(micro)  - Fix(Shift(micro,  1))
			ELSE: Message, 'Unexpected unit for lower string label.'
		ENDCASE
		iZero = Where(u2 EQ 0, nZero)
		IF nZero GT 0 THEN s2[iZero] = ''
		
		;Does the upper level change?
		iu2 = Where(u2 NE 0, nu2)
		IF nu2 EQ 0 $
			THEN side_string = StrMid(times[0], 0, pos[setup[iTick[0]].fld1]) $
			ELSE side_string = StrMid(times[0], 0, pos[setup[iTick[0]].fld2])
	ENDIF ELSE BEGIN
		side_string = ''
	ENDELSE
	
	str   = s1 + '!C' + s2
	label = labels[setup[iTick[0]].fld1] + '!C' + side_string
	
	out = { XTICKNAME: str, $ 
	        XTICKV:    ticks, $
	        XTICKS:    N_Elements(ticks) - 1, $
	        XMINOR:    setup[iTick[0]].nminor, $
	        XRANGE:    trange_ssm, $
	        XSTYLE:    1, $
	        XTITLE:    '' }
	
	RETURN, out
END