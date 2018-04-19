; docformat = 'rst'
;
; NAME:
;       MrVar_TLimit
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
;   Restrict time range for certain variables.
;
;   NOTE:
;       The same operation can be achieved via bracket overloading. Namely, if MrVAR
;       is an NxM MrVariable object that depends on time T, an Nx1 MrTimeVar object,
;       and i0 and i1 are indices of the time interval to be kept, then
;
;           RESULT = MrVar[i0:i1,*]
;
;       Is similar to 
;
;           MrVar_TLimit, MrVar, TRANGE
;
;       with the following differences. In the former case, a new (uncached) variable
;       is created with a new name and new dependent variables. In the latter case,
;       MrVar is altered directly, it remains in the cache (if it is there to begin with),
;       and any variable passed in with MrVar that shares the same time object will
;       continue to share the same time object.
;
;   Calling Sequence::
;      MrVar_TLimit, varnames
;      MrVar_TLimit, varnames, trange
;
; :Categories:
;   MrVariable
;
; :Params:
;       VARIABLES:      in, required, type=int/str/obj
;                       Names, numbers, or objrefs of the variables for which to limit
;                           the time range. If undefined, the names of all variables
;                           present in the cache are used.
;       TRANGE:         in, optional, type=strarr(2), default=MrVar_GetTRange()
;                       Time interval by which to restrict data, formatted as
;                           ISO-8601 date-time strings. Sets the global time interval via
;                           MrVar_SetTRange. If TRANGE is not set and no time range has
;                           been establish, the program will exit without altering the
;                           data.
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
;       2016-11-22  -   Written by Matthew Argall
;       2016-12-09  -   Variable numbers or objrefs can be passed in. Variable objrefs
;                           do not have to exist in the cache. - MRA
;       2017-05-03  -   If no time range has been established via the `TRANGE` parameter
;                           or MrVar_SetTRange, silently return. - MRA
;-
PRO MrVar_TLimit, variables, trange
	Compile_Opt idl2
	On_Error, 2

	;Defaults
	IF N_Elements(variables) EQ 0 THEN MrVar_Names, variables
	IF N_Elements(trange) GT 0 && trange[0] NE '' THEN MrVar_SetTRange, trange

	;Number OF variables
	;   - Scalar objects with ::_overloadSize can return any number of elements
	;   - Size(/TNAME) does not distinguish between objects and object arrays
	IF Size(variables, /TNAME) EQ 'OBJREF' $
		THEN nVars = TypeName(variables) EQ 'OBJREF' ? N_Elements(variables) : 1 $
		ELSE nVars = N_Elements(variables)
	
	;Get the time range
	;   - If no time range has been established, return without doing anything.
	Catch, the_error
	IF the_error EQ 0 THEN BEGIN
		trange     = MrVar_GetTRange()
		trange_ssm = MrVar_GetTRange('SSM')
	ENDIF ELSE BEGIN
		Catch, /CANCEL
		RETURN
	ENDELSE
	Catch, /CANCEL
	On_Error, 2

;-----------------------------------------------------
; Pick out the Epoch Variables \\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Find epoch variables
	nTime = 0
	oTime = ObjArr(nVars)
	tName = StrArr(nVars)
	tID   = ULonArr(nVars)
	theVars = MrVar_Get(variables)

	;Loop over each variable
	FOR i = 0, nVars - 1 DO BEGIN
		;Ensire the variable is valid
		theVar = theVars[i]
		IF ~Obj_Valid(theVar) THEN continue
		
		;TIME
		IF Obj_IsA(theVar, 'MrTimeVar') THEN BEGIN
			tempTime = theVar
		
		;DATA - Get DEPEND_0 object
		ENDIF ELSE IF theVar -> HasAttr('DEPEND_0') THEN BEGIN
			tempTime = MrVar_Get(theVar['DEPEND_0'])
			IF ~Obj_IsA(tempTime, 'MrTimeVar') THEN CONTINUE
		
		;NEITHER
		ENDIF ELSE BEGIN
			tempTime = Obj_New()
		ENDELSE
		
		;Store the epoch variable
		heapID = Obj_Valid( tempTime, /GET_HEAP_IDENTIFIER )
		IF heapID GT 0 THEN BEGIN
			IF ~MrIsMember( tID, heapID ) THEN BEGIN
				oTime[nTime] = tempTime
				tName[nTime] = tempTime.name
				tID[nTime]   = heapID
				nTime       += 1
			ENDIF
		ENDIF
	ENDFOR

	;Return IF no time variables
	IF nTime EQ 0 THEN RETURN
	
	;Trim down results
	tName = tName[0:nTime-1]
	oTime = oTime[0:nTime-1]

;-----------------------------------------------------
; Loop Over Epoch Variables \\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	FOR i = 0, nTime - 1 DO BEGIN

	;-----------------------------------------------------
	; Set Epoch Range \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		;Get the variable
		theTime = oTime[i]
		
		;Convert to SSM and find relevant data
		t_ssm = theTime -> GetData('SSM', trange[0])
		iKeep = Where( (t_ssm GE trange_ssm[0]) and (t_ssm LE trange_ssm[1]), nKeep)
		IF nKeep EQ 0 THEN BEGIN
			iKeep = Value_Locate(t_ssm, trange_ssm[0]) > 0
			CASE 1 OF
				iKeep EQ N_Elements(t_ssm)-1: iKeep = [iKeep-1, iKeep]
				iKeep EQ 0:                   iKeep = [iKeep, iKeep+1]
				ELSE: BEGIN
					IF t_ssm[iKeep]-trange_ssm[0] LT trange_ssm[0]-t_ssm[iKeep] $
						THEN iKeep = [iKeep, iKeep+1] $
						ELSE iKeep = [iKeep-1, iKeep]
				ENDCASE
			ENDCASE
			MrPrintF, 'LogText', 'No data found in time interval for variable "' + theTime.name + '".'
			MrPrintF, 'LogText', '  Interval: [' + StrJoin(trange, ', ') + ']'
			MrPrintF, 'LogText', '  Closest:  [' + StrJoin(theTime['DATA', iKeep], ', ') + ']'
			MrPrintF, 'LogWarn', 'Choosing closest two points.'
		
		ENDIF ELSE IF nKeep EQ 1 THEN BEGIN
			IF t_ssm[iKeep]-trange_ssm[0] LT trange_ssm[0]-t_ssm[iKeep] $
				THEN iKeep = [iKeep, iKeep+1] $
				ELSE iKeep = [iKeep-1, iKeep]
			MrPrintF, 'LogText', 'Only one point found in time interval for variable "' + theTime.name + '".'
			MrPrintF, 'LogText', '  Interval: [' + StrJoin(trange, ', ') + ']'
			MrPrintF, 'LogText', '  Closest:  [' + StrJoin(theTime['DATA', iKeep], ', ') + ']'
			MrPrintF, 'LogWarn', 'Including next closest point.'
		ENDIF

		;Trim the time data
		newEpoch  = theTime[iKeep]
		newEpoch -> SetName, 'TLimit(' + theTime.name + ')'

	;-----------------------------------------------------
	; Set Data Range \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	;-----------------------------------------------------
		FOR j = 0, nVars - 1 DO BEGIN
			;Skip IF the variable does not have an DEPEND_0 attribute
			theVar = theVars[j]
			IF ~theVar -> HasAttr('DEPEND_0') THEN continue
			
			;Skip IF the DEPEND_0 attribute is not the one we are looking at
			;   - DEPEND_0 may be an object or a variable name
			IF ~theVar -> IsTimeIdentical( theTime ) THEN continue
			
			;Trim the variable data
			CASE obj_class(theVar) OF
				'MRSCALARTS': theVar -> SetData, newEpoch, theVar['DATA',iKeep], DIMENSION=1
				'MRVECTORTS': theVar -> SetData, newEpoch, theVar['DATA',iKeep,*], DIMENSION=1
				'MRMATRIXTS': theVar -> SetData, newEpoch, theVar['DATA',iKeep,*,*], DIMENSION=1
				'MRTIMESERIES': BEGIN
					CASE theVar.n_dimensions OF
						1: theVar -> SetData, newEpoch, theVar['DATA',iKeep], DIMENSION=1
						2: theVar -> SetData, newEpoch, theVar['DATA',iKeep,*], DIMENSION=1
						3: theVar -> SetData, newEpoch, theVar['DATA',iKeep,*,*], DIMENSION=1
						4: theVar -> SetData, newEpoch, theVar['DATA',iKeep,*,*,*], DIMENSION=1
						5: theVar -> SetData, newEpoch, theVar['DATA',iKeep,*,*,*,*], DIMENSION=1
						6: theVar -> SetData, newEpoch, theVar['DATA',iKeep,*,*,*,*,*], DIMENSION=1
						7: theVar -> SetData, newEpoch, theVar['DATA',iKeep,*,*,*,*,*,*], DIMENSION=1
						8: BEGIN
							data = theVar -> GetData()
							theVar -> SetData, newEpoch, (Temporary(data))[iKeep,*,*,*,*,*,*,*], DIMENSION=1
						ENDCASE
						ELSE: Message, 'Variable has unexpected number of dimensions: "' + theVar.name + '".'
					ENDCASE
				ENDCASE
				
				;TODO: Determine which is the time-varying dimension
				'MrVariable': MrPrintF, 'LogWarn', 'TLimit not possible for object class MrVariable.'
				ELSE: MrPrintF, 'LogWarn', 'Unknown object of class "' + Obj_Class(theVar) + '".'
			ENDCASE
		ENDFOR
	ENDFOR
END