; docformat = 'rst'
;
; NAME:
;       MrMMS_FPI_Load_Models
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
;   Find and read MMS FPI internal photoelectron models. Involves multiple steps.
;
;       1. Find moments file
;       2. Extract stepper-id from metadata
;       3. Extract scale factor from metadata
;       4. Find photoelectron model files
;       5. Load model
;
; :Categories:
;       MMS, FPI
;
; :Params:
;       MODE:               in, required, type=string
;                           Data telemetry mode. Options are: {'fast' | 'brst'}
;       STEPPER_ID:         in, required, type=string
;                           The energy table tag.
;
; :Keywords:
;       SUFFIX:             in, optional, type=string, default=''
;                           A suffix to be appended to the variable names.
;       VARFORMAT:          in, optional, type=strarr/strarr
;                           A variable name filter. Those that do not match VARFORMAT
;                               will not be read into the cache.
;       VARNAMES:           out, optional, type=strarr
;                           A named variable to receive the names of each variable read
;                               into the cache.
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
;       2016/01/18  -   Written by Matthew Argall
;       2016/01/24  -   Extract information from moments files. - MRA
;       2016/01/25  -   Substitute angle and energy table indices for their targets. - MRA
;       2016/05/01  -   Do not assume moments files are already saved locally. - MRA
;       2017/09/13  -   Updated to work with survey data. - MRA
;-
PRO MrMMS_FPI_Load_Models, sc, mode, species, $
LEVEL=level, $
SUFFIX=suffix, $
VARFORMAT=varformat, $
VARNAMES=varnames
	Compile_Opt idl2

	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF Obj_Valid(oURI) THEN Obj_Destroy, oURI
		MrPrintF, 'LogErr'
		RETURN
	ENDIF
	
	IF N_Elements(suffix) EQ 0 THEN suffix = ''

;-----------------------------------------------------
; Extract Metadata from Moments File \\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Get the local file name
	fmoms = MrMMS_Get_Data( sc, 'fpi', mode, level, $
	                        COUNT   = nMoms, $
	                        OPTDESC = 'd' + species + 's-moms' )
	
	;Read the photoelectron model global attributes
	IF N_Elements(nMoms) GT 0 THEN MrPrintF, 'LogWarn', 'More than one moments file. Taking model from first.'
	oCDF = MrCDF_File(fmoms[0])
	scl  = oCDF -> GetGlobalAttrValue('Photoelectron_model_scaling_factor')
	fphe = oCDF -> GetGlobalAttrValue('Photoelectron_model_filenames')
	Obj_Destroy, oCDF

;-----------------------------------------------------
; Find Photo-Electron Model \\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Get the photo electron model
	stepper_id = StRegEx( fphe, 'p([0-9]+-[0-9]+)', /SUBEXP, /EXTRACT )
	stepper_id = stepper_id[1]
	
	;Create the photoelectron model file name
	instr = 'd' + species + 's'
	url   = 'https://lasp.colorado.edu/mms/sdc/public/data/models/fpi/'
	fname = StrJoin( ['mms', 'fpi', mode, 'l2', 'd'+species+'s-bgdist', 'v*', 'p'+stepper_id], '_' ) + '.cdf'
	
	;URI Object
	oURI  = MrWebURI( url, $
	                  LOCAL_ROOT = !MrMMS.local_root, $
	                  NSEGMENT   = 7 )
	
	;Get the files
	fname = oURI -> Path_Append(fname, ROOT=url)
	files = oURI -> Get(fname, COUNT=count)
	Obj_Destroy, oURI
	
	IF count EQ 0 THEN Message, 'Unable to find model files: "' + fname + '".'

;-----------------------------------------------------
; Scale Factor \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Create a variable for the scale factor
	scl_vname = StrJoin( [sc, instr, 'scl', 'model', mode], '_' ) + suffix
	oScl = MrVariable( [Float(scl)], $
	                   /CACHE, $
	                   NAME = scl_vname )
	oScl['NOTES'] = 'After reconstructing the photo electron distribution function, ' + $
	                'multiply by this scaling factor.'

;-----------------------------------------------------
; Load the Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Load the data
	MrVar_ReadCDF, files, $
	               SUFFIX    = suffix, $
	               VARFORMAT = ['*bgdist*', '*phi*', '*theta*', '*energy?_'+mode+'*'], $
	               VARNAMES  = varnames
	
	;Attach the proper tables
	;   - Energy tables will have to be adjusted when fphoto is built.
	;   - See MrMMS_FPI_Dist_Photo.pro
	phi_vname   = StrJoin( ['mms', 'd'+species+'s', 'phi',            mode], '_' ) + suffix
	theta_vname = StrJoin( ['mms', 'd'+species+'s', 'theta',          mode], '_' ) + suffix
	
	IF mode EQ 'brst' THEN BEGIN
		;One distribution per polarity is given
		f0_vname    = StrJoin( ['mms', 'd'+species+'s', 'bgdist', 'p0',   mode], '_' ) + suffix
		f1_vname    = StrJoin( ['mms', 'd'+species+'s', 'bgdist', 'p1',   mode], '_' ) + suffix
		
		;Set the dependent variables
		of0             = MrVar_Get(f0_vname)
		of0['DEPEND_1'] = phi_vname
		of0['DEPEND_2'] = theta_vname
		of1             = MrVar_Get(f1_vname)
		of1['DEPEND_1'] = phi_vname
		of1['DEPEND_2'] = theta_vname
	ENDIF ELSE BEGIN
		f_vname = StrJoin( ['mms', 'd'+species+'s', 'bgdist', mode], '_' ) + suffix
		of0 = MrVar_Get(f_vname)
		of0['DEPEND_1'] = phi_vname
		of0['DEPEND_2'] = theta_vname
	ENDELSE
	
	;Add the scaling factor to the variable names
	varnames = [varnames, scl_vname]
END
