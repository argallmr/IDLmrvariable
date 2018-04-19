; docformat = 'rst'
;
; NAME:
;       MrMMS_EDI_traj2cart
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
;   Transform trajectory vectors from spherical to cartesian coordinates.
;
; :Categories:
;   MrVariable, MMS, EDI
;
; :Params:
;       TRAJ:       in, required, type=int/string/objref
;                   The index, name, or MrTimeSeries object of EDI incident trajectory
;                       vectors in spherical coordinates [azimuth, polar].
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the output will be added to the variable cache.
;       NAME:           in, optional, type=string, default='Cyclotron_Frequency'
;                       Name to be given to the output variable.
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
;       2017/09/30  -   Written by Matthew Argall
;       2018/03/08  -   Renamed from MrMMS_EDI_traj2cart to MrMMS_EDI_sphr2cart. - MRA
;-
FUNCTION MrMMS_EDI_sphr2cart, traj, $
CACHE=cache, $
NAME=name
	Compile_Opt idl2
	On_Error, 2
	
	oTraj = MrVar_Get(traj)
	
	;Convert to radians
	azimuth = oTraj['DATA',*,0] * !dpi / 180D
	polar   = oTraj['DATA',*,1] * !dpi / 180D
	
	;Compute cartesian coordinates
	x = Sin(polar) * Cos(azimuth)
	y = Sin(polar) * Sin(azimuth)
	z = Cos(polar)
	
	;Create vector
	oV = MrVectorTS( oTraj['TIMEVAR'], Float([ [x], [y], [z] ]), $
	                 CACHE = cache, $
	                 NAME  = name )
	
	;Attributes
	oTraj -> CopyAttrTo, oV, ['CATDESC', 'FIELDNAM']
	oV['AXIS_RANGE'] = [-1.0, 1.0]
	oV['LABEL']      = ['X', 'Y', 'Z']
	
	RETURN, oV
END