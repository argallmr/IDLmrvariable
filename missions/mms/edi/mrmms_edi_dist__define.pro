; docformat = 'rst'
;
; NAME:
;   MrMMS_EDI_Dist__Define
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
; PURPOSE
;+
;   Return distribution function of EDI counts.
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
;       2017/01/22  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   Get the angles of indicence of the particle trajectories.
;
; :Params:
;       OPHI:           out, optional, type=objref
;                       Azimuth angle of particle trajectories (MrScalarTS variable).
;       OTHETA:         out, optional, type=objref
;                       Polar angle of particle trajectories (MrScalarTS variable).
;
; :Keywords:
;       TMATRIX:        in, optional, type=objref
;                       A MrMatrixTS object that rotates the angles into a new coordinate
;                           system.
;-
PRO MrMMS_EDI_Channel::GetAngles, oPhi, oTheta, $
TMATRIX=TMatrix
	Compile_Opt idl2
	On_Error, 2
	
	;Grab the data
	IF N_Elements(TMatrix) EQ 0 THEN BEGIN
		oPhi   = self.oTraj[*,0]
		oTheta = self.oTraj[*,1]
	
	;Rotate the data
	ENDIF ELSE BEGIN
		;Rotate
		self -> Rotate, TMatrix, phi, theta
		
		;PHI
		oPhi  = MrScalarTS( self.oTraj['TIMEVAR'], phi, /NO_COPY )
		oPhi['AXIS_RANGE'] = [0, 360]
		oPhi['PLOT_TITLE'] = 'Azimuthal Angle of Particle Trajectories'
		oPhi['TITLE']      = 'Phi'
		oPhi['UNITS']      = 'degrees'

		;THETA
		oTheta  = MrScalarTS( self.oTraj['TIMEVAR'], theta, /NO_COPY )
		oTheta['AXIS_RANGE'] = [0, 180]
		oTheta['PLOT_TITLE'] = 'Polar Angle of Particle Trajectories'
		oTheta['TITLE']      = 'Theta'
		oTheta['UNITS']      = 'degrees'
	ENDELSE
END


;+
;   Get object properties.
;
; :Params:
;       COUNTS:         out, optional, type=objref
;                       Counts or flux measurements
;       CHANNEL:        out, optional, type=int
;                       Channel from which meaurments were made.
;       GDU:            out, optional, type=int
;                       GDU that took the measurement.
;       PA_STATE:       out, optional, type=int
;                       Pitch-angle state.
;       TRAJ:           out, optional, type=objref
;                       Trajectory of incoming electrons that registar as `COUNTS`
;-
PRO MrMMS_EDI_Channel::GetProperty, $
COUNTS=counts, $
CHANNEL=channel, $
GDU=gdu, $
PA_STATE=pa_state, $
TRAJ=traj
	Compile_Opt idl2
	On_Error, 2
	
	;Get properties
	IF Arg_Present(counts)   THEN counts   = self.oCounts
	IF Arg_Present(traj)     THEN traj     = self.oTraj
	IF Arg_Present(channel)  THEN channel  = self.channel
	IF Arg_Present(gdu)      THEN gdu      = self.gdu
	IF Arg_Present(pa_state) THEN pa_state = self.pa_state
END


;+
;   Rotate the look direction.
;
; :Params:
;       T:              in, required, type=objref
;                       A MrMatrixTS variable of the coordinate system transformation matrix
;       PHI:            out, optional, type=fltarr
;                       The rotated azimuthal look direction.
;       THETA:          out, optional, type=fltarr
;                       The rotated polar look direction.
;-
PRO MrMMS_EDI_Channel::Rotate, T, phi, theta
	Compile_Opt idl2
	On_Error, 2
	
	;Get the matrix
	oT = MrVar_Get(T)
	IF ~Obj_IsA(oT, 'MRMATRIXTS') THEN Message, 'T must be a MrMatrixTS object.'
	IF ~oT -> IsTimeIdentical(self.oTraj['TIMEVAR']) THEN Message, 'T must have the same time stamps as TRAJ.'
	
	;Size of input
	nDims = size(self.oTraj, /N_DIMENSIONS)
	dims  = size(self.oTraj, /DIMENSIONS)
	nTime = dims[0]
	nChan = nDims eq 2 ? 1 : dims[2]
	
;-----------------------------------------------------
; Rotate \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	deg2rad = !dpi / 180D
	phi     = self.oTraj['DATA',*,0] * deg2rad
	theta   = self.oTraj['DATA',*,1] * deg2rad
	
	;Compute the vector components
	traj_vec      = fltarr(nTime, 3)
	traj_vec[*,0] = Sin(theta) * Cos(phi)
	traj_vec[*,1] = Sin(theta) * Sin(phi)
	traj_vec[*,2] = Cos(theta)
	
	;Rotate vector components
	traj_fac      = FltArr(nTime, 3)
	traj_fac[*,0] = oT['DATA',*,0,0] * traj_vec[*,0] + oT['DATA',*,1,0] * traj_vec[*,1] + oT['DATA',*,2,0] * traj_vec[*,2]
	traj_fac[*,1] = oT['DATA',*,0,1] * traj_vec[*,0] + oT['DATA',*,1,1] * traj_vec[*,1] + oT['DATA',*,2,1] * traj_vec[*,2]
	traj_fac[*,2] = oT['DATA',*,0,2] * traj_vec[*,0] + oT['DATA',*,1,2] * traj_vec[*,1] + oT['DATA',*,2,2] * traj_vec[*,2]

	;Compute polar coordinates
	rad2deg = 180D / !dpi
	phi     = ATan( Reform(traj_fac[*,1]), Reform(traj_fac[*,0]) ) * rad2deg
	theta   = ACos( Reform(traj_fac[*,2]) ) * rad2deg

	;Clear data
	traj_fac = !Null
	traj_vec = !Null
END


;+
;   Object destructor method.
;-
PRO MrMMS_EDI_Channel::CleanUp
	;
	; Do not destroy the oCounts or oTraj properties.
	;   - They contain data that is stored in the variable
	;     cache and will be cleaned up independently.
	;
END


;+
;   The initialization method.
;
; :Params:
;       COUNTS:         in, required, type=string/int/objref
;                       The name, number, or object reference of a MrScalarTS variable
;                          containing the counts or flux measurements
;       TRAJ:           in, required, type=Nx2 fltarr
;                       The name, number, or object reference of a MrScalarTS variable
;                          containing the trajectory of incoming electrons that registar
;                          as `COUNTS`
;       CHANNEL:        in, optional, type=int
;                       Channel from which meaurments were made. Options are: { 1 | 2 | 3 | 4 }
;       GDU:            in, optional, type=int
;                       GDU that took the measurement. Options are: { 1 | 2 }
;       PA_STATE:       in, optional, type=int
;                       Pitch-angle state. Options are: { 'q0' | '0' | '90' | '180' }
;
; :Returns:
;       If instantiation is successful, a valid MrMMS_EDI_Channel object is returned.
;-
FUNCTION MrMMS_EDI_Channel::INIT, counts, traj, channel, gdu, pa_state
	Compile_Opt idl2

	;Error handling
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, 0
	ENDIF
	
	;Check inputs
	IF N_Elements(counts) EQ 0 THEN Message, 'COUNTS must be defined.'
	IF N_Elements(traj)   EQ 0 THEN Message, 'TRAJ must be defined.'
	IF ~MrIsMember( [1, 2, 3], gdu)          THEN Message, 'GDU must be {1 | 2 | 3}.'
	IF ~MrIsMember( [1, 2, 3, 4], channel)   THEN Message, 'CHANNEL must be {1 | 2 | 3 | 4}.'
	IF ~MrIsMember( [0, 90, 180], pa_state ) THEN Message, 'PA_STATE must be {0 |  90 | 180}.'

	;Set properties
	self.oCounts  = MrVar_Get(counts)
	self.oTraj    = MrVar_Get(traj)
	self.channel  = channel
	self.gdu      = gdu
	self.pa_state = pa_state

	RETURN, 1
END


;+
;   The class definition statement.
;
; :Params:
;       CLASS:          out, optional, type=structure
;-
PRO MrMMS_EDI_Channel__DEFINE
	Compile_Opt idl2
	
	class = { MrMMS_EDI_Channel, $
	          Inherits IDL_Object, $
	          oCounts:  Obj_New(), $
	          oTraj:    Obj_New(), $
	          gdu:      0B, $
	          channel:  0B, $
	          pa_state: 0B $
	        }
END


;+
;   The initialization method.
;
;   Calling Sequence:
;       oEDI = MrMMS_EDI_Dist()
;       oEDI = MrMMS_EDI_Dist(sc, mode)
;       oEDI = MrMMS_EDI_Dist(sc, mode, optdesc)
;
; :Params:
;       SC:             in, required, type=string
;                       MMS spacecraft identifier. Options are: {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:           in, required, type=string
;                       Telemetry mode. Options are: {'slow' | 'fast' | 'srvy' | 'brst'}
;       OPTDESC:        in, optional, type=string, defualt='q0'
;                       EDI-specific optional descriptor. Options are: {'amb' | 'amb-pm2' |
;                           'alt-cc' | 'alt-oc' | 'alt-oob' | 'alt-ooc', 'q0'}
;
; :Keywords:
;       _REF_EXTRA:     in, optional, type=boolean, default=0
;                       Any keyword accepted by the ::Load method is also accepted here.
;-
FUNCTION MrMMS_EDI_Dist::INIT, sc, mode, optdesc, $
_REF_EXTRA=extra
	Compile_Opt idl2

	;Error handling
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, 0
	ENDIF
	
	self.fac = 'NONE'
	
	;Load data
	IF N_Elements(sc) GT 0 THEN BEGIN
		self -> Load, sc, mode, optdesc, $
		              _STRICT_EXTRA = extra
	ENDIF

	RETURN, 1
END


;+
;   Clean up after the object is destroyed
;-
PRO MrMMS_EDI_Dist::CLEANUP
	Compile_Opt idl2
	On_Error, 2
	
	;Cleanup superclasses
	self -> IDL_Container::Cleanup
END


;+
;   Add a channel to the container.
;
; :Params:
;       OCHANNELS:      in, required, type=MrMMS_EDI_Channel object
;                       The channel(s) to be added.
;
; :Keywords:
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by IDL_Container::Add.
;-
PRO MrMMS_EDI_Dist::Add, oChannel, $
_REF_EXTRA=extra
	Compile_Opt idl2
	On_Error, 2
	
	IF Obj_Class(oChannel) NE 'MRMMS_EDI_CHANNEL' $
		THEN Message, 'OCHANNEL must be an MrMMS_EDI_Channel object.'
		
	;Add to the container
	self -> IDL_Container::Add, oChannel, $
	                            _STRICT_EXTRA = extra
END


;+
;   Bin EDI data by pitch or gyrophase angle. Resulting data is added to the
;   MrVariable cache.
;
; :Params:
;       OCOUNTS:        in, required, type=MrScalarTS object
;                       Counts data to be binned.
;       OPHI:           in, required, type=MrScalarTS object
;                       Azimuthal angle at which `OCOUNTS` was measured.
;       OTHETA:         in, required, type=MrScalarTS object
;                       Polar angle at which `OCOUNTS` was measured.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the distribution is added to the variable cache.
;       DPHI:           in, optional, type=float, default=1.5
;                       Width of each azimuthal bin (degrees). Cannot be used with `NPHI_BINS`.
;       DTHETA:         in, optional, type=float, default=0.5
;                       Width of each polar angle bin (degrees). Cannot be used with `NTHETA_BINS`.
;       DTIME:          in, optional, type=float, default=5.0
;                       Duration, in seconds, of each time bin.
;       NAME:           in, optional, type=string, default=sc_instr_phidist_optdesc_mode_level
;                       Name to be given to the distribution.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of bins spanning `PHI_RANGE`. Cannot be used with `DPHI`.
;       NTHETA_BINS:    in, optional, type=integer
;                       Number of bins spanning `THETA_RANGE`. Cannot be used with `DTHETA`.
;       PHI_RANGE:      in, optional, type=fltarr(2)
;                       Minimum and maximum azimuthal angles over which to bin (degrees).
;       THETA_RANGE:    in, optional, type=fltarr(2)
;                       Minimum and maximum polar angles over which to bin (degrees).
;-
FUNCTION MrMMS_EDI_Dist::Bin, oCounts, oPhi, oTheta, $
CACHE=cache, $
DPHI=dphi, $
DTIME=dt, $
DTHETA=dtheta, $
NAME=name, $
NPHI_BINS=nPhi_bins, $
NTHETA_BINS=nTheta_bins, $
PHI_RANGE=phi_range, $
THETA_RANGE=theta_range
	Compile_Opt idl2
	On_Error, 2

;-------------------------------------------
; Defaults /////////////////////////////////
;-------------------------------------------
	;Defaults
	tf_cache = Keyword_Set(cache)
	IF N_Elements(name)        EQ 0 THEN name        = StrJoin([self.sc, 'edi', 'dist', self.optdesc, self.mode, self.level], '_')
	IF N_Elements(dt)          EQ 0 THEN dt          = 1.0
	IF N_Elements(phi_range)   EQ 0 THEN phi_range   = [-180,180]
	IF N_Elements(theta_range) EQ 0 THEN theta_range = [0, 180]
	IF N_Elements(dPhi)   EQ 0 && N_Elements(nPhi_bins)   EQ 0 THEN dPhi   = 5.0
	IF N_Elements(dTheta) EQ 0 && N_Elements(nTheta_bins) EQ 0 THEN dTheta = 0.5

	;Restrictions
	IF N_Elements(nPhi_bins) GT 0 && N_Elements(dPhi) GT 0 $
		THEN Message, 'NPHI_BINS and DPHI are mutually exclusive.'
	IF N_Elements(nTheta_bins) GT 0 && N_Elements(dTheta) GT 0 $
		THEN Message, 'NTHETA_BINS and DTHETA are mutually exclusive.'

;-------------------------------------------
; Bin Information //////////////////////////
;-------------------------------------------
	
	;Time range
	oTime         = oCounts['TIMEVAR']
	trange        = [oTime['DATA', 0, 'TT2000'], oTime['DATA', -1, 'TT2000']]
	trange_ssm    = MrCDF_epoch2ssm(trange)
	trange_ssm[0] = trange_ssm[0] - (trange_ssm[0] mod dt)
	trange_ssm[1] = trange_ssm[1] - (trange_ssm[1] mod dt) + dt

	;Time bins
	;   - BINSIZE = (MAX – MIN) / (NBINS – 1)
	t_ssm   = oTime -> GetData('SSM')
	nt_bins = Long( (trange_ssm[1] - trange_ssm[0]) / dt ) + 1
	t_bins  = linspace(trange_ssm[0], trange_ssm[1], nt_bins)
	
	;Phi bins
	IF N_Elements(dPhi) GT 0 $
		THEN nPhi_bins = Long( (phi_range[1] - phi_range[0]) / dPhi ) + 1 $
		ELSE dPhi      = (phi_range[1] - phi_range[0]) / nPhi_bins
	phi_bins = linspace(phi_range[0], phi_range[1], nPhi_bins)
	
	;Theta bins
	IF N_Elements(dTheta) GT 0 $
		THEN nTheta_bins = Long( (theta_range[1] - theta_range[0]) / dTheta ) + 1 $
		ELSE dTheta      = (theta_range[1] - theta_range[0]) / nTheta_bins
	theta_bins = linspace(theta_range[0], theta_range[1], nTheta_bins)

;-------------------------------------------
; Bin by Time //////////////////////////////
;-------------------------------------------
	;Allocate memory
	theDist = ULonArr(nt_bins, nPhi_bins, nTheta_bins)
	
	;Histogram in time
	tHist  = Histogram( Temporary(t_ssm), $
	                    MIN             = trange_ssm[0], $
	                    MAX             = trange_ssm[1], $
	                    NBINS           = nt_bins, $
	                    REVERSE_INDICES = rit )
	
	;Loop over each tiem bin
	FOR i = 0, N_Elements(tHist) - 1 DO BEGIN
		nt = rit[i+1] - rit[i]
		IF nt EQ 0 THEN continue
		
		;Indices within the source
		it = rit[rit[i]:rit[i+1]-1]

	;-------------------------------------------
	; Bin by Phi ///////////////////////////////
	;-------------------------------------------
		;Histogram in Phi
		phiHist = Histogram( oPhi['DATA', it], $
		                     MIN             = phi_range[0], $
		                     MAX             = phi_range[1], $
		                     NBINS           = nPhi_bins, $
		                     REVERSE_INDICES = riPhi )
		
		;Loop over each Phi bin
		FOR j = 0, N_Elements(phiHist) - 1 DO BEGIN
			nPhi = riPhi[j+1] - riPhi[j]
			IF nPhi EQ 0 THEN continue

			;Indices within the source
			iPhi  = riPhi[riPhi[j]:riPhi[j+1]-1]
			itPhi = it[iPhi]

		;-------------------------------------------
		; Bin by Theta /////////////////////////////
		;-------------------------------------------
			;Histogram in theta
			tempTheta = nPhi EQ 1 ? [oTheta['DATA', itPhi]] : oTheta['DATA', itPhi]
			thetaHist = Histogram( Temporary(tempTheta), $
			                       MIN             = theta_range[0], $
			                       MAX             = theta_range[1], $
			                       NBINS           = nTheta_bins, $
			                       REVERSE_INDICES = riTheta )

			;Loop over each Theta bin
			FOR k = 0, N_Elements(thetaHist) - 1 DO BEGIN
				nTheta = riTheta[k+1] - riTheta[k]
				IF nTheta EQ 0 THEN continue
				
				;Inidices within the source
				iTheta     = riTheta[riTheta[k]:riTheta[k+1]-1]
				itPhiTheta = itPhi[iTheta]

				;Save the data
				;   - Sum all counts
				;   - Normalize by time spent in each bin
				theDist[i,j,k] = Total(oCounts['DATA', itPhiTheta]) / N_Elements(itPhiTheta)
			ENDFOR
		ENDFOR
	ENDFOR
	theDist = Round(theDist)

;-------------------------------------------
; Output Data //////////////////////////////
;-------------------------------------------

	;Time Bins
	oBinnedT  = MrTimeVar(t_bins+dt/2.0, 'SSM', /NO_COPY, T_REF=oTime['DATA', [0]])
	oBinnedT -> AddAttr, 'DELTA_PLUS',  dt/2.0
	oBinnedT -> AddAttr, 'DELTA_MINUS', dt/2.0
	oBinnedT -> AddAttr, 'TIME_BASE',   'J2000'
	oBinnedT -> AddAttr, 'TITLE',       'Center time tags FOR EDI Q0 Phi Spectrogram'
	oBinnedT -> AddAttr, 'UNITS',       'ns'

	;Phi Bins
	IF nPhi_bins GT 1 THEN BEGIN
		oPhiBins  = MrVariable(phi_bins + dPhi/2.0, /NO_COPY)
		oPhiBins -> AddAttr, 'AXIS_RANGE',  phi_range
		oPhiBins -> AddAttr, 'DELTA_PLUS',  dPhi/2.0
		oPhiBins -> AddAttr, 'DELTA_MINUS', dPhi/2.0
		oPhiBins -> AddAttr, 'PLOT_TITLE',  'Azimuthal Angle Bins'
		oPhiBins -> AddAttr, 'TITLE',       'Phi!C(deg)'
		oPhiBins -> AddAttr, 'UNITS',       'degrees'
	ENDIF

	;Theta Bins
	IF nTheta_bins GT 1 THEN BEGIN
		oThetaBins  = MrVariable(theta_bins + dTheta/2.0, /NO_COPY)
		oThetaBins -> AddAttr, 'AXIS_RANGE',  theta_range
		oThetaBins -> AddAttr, 'DELTA_PLUS',  dTheta/2.0
		oThetaBins -> AddAttr, 'DELTA_MINUS', dTheta/2.0
		oThetaBins -> AddAttr, 'PLOT_TITLE',  'Polar Angle Bins'
		oThetaBins -> AddAttr, 'TITLE',       'Theta!C(deg)'
		oThetaBins -> AddAttr, 'UNITS',       'degrees'
	ENDIF

	;Binned Counts
	oBinnedCts  = MrTimeSeries(oBinnedT, Reform(theDist), /NO_COPY, DIMENSION=1)
	oBinnedCts -> AddAttr, 'PLOT_TITLE',    'EDI Quality 0 Binned Electron Counts'
	oBinnedCts -> AddAttr, 'MISSING_VALUE', 0
	oBinnedCts -> AddAttr, 'SCALE',         1
	oBinnedCts -> AddAttr, 'TITLE',         'Normalized Counts'
	oBinnedCts -> AddAttr, 'UNITS',         'Counts'
	
	IF nPhi_bins GT 1 && nTheta_bins GT 1 THEN BEGIN
		oBinnedCts -> AddAttr, 'DEPEND_1', oPhiBins
		oBinnedCts -> AddAttr, 'DEPEND_2', oThetaBins
	ENDIF ELSE IF nPhi_bins GT 1 THEN BEGIN
		oBinnedCts -> AddAttr, 'DEPEND_1', oPhiBins
	ENDIF ELSE IF nTheta_bins GT 1 THEN BEGIN
		oBinnedCts -> AddAttr, 'DEPEND_1', oThetaBins
	ENDIF

;-------------------------------------------
; Name, Cache & Output /////////////////////
;-------------------------------------------
	;Name & Cache
	oBinnedCts -> SetName, name
	IF tf_cache THEN oBinnedCts -> Cache
	
	;Return
	RETURN, oBinnedCts
end


;+
;   Find data from a particular channel, gdu, and pitch-angle state.
;
; :Keywords:
;       COUNT:          out, optional, type=int, default=[1,2,3,4]
;                       Number of channels found.
;       CHANNEL:        in, optional, type=int, default=[1,2,3,4]
;                       Channel from which meaurments were made.
;       GDU:            in, optional, type=int, default=[1,2]
;                       GDU that took the measurement.
;       PA_STATE:       in, optional, type=int, default=['q0', '0', '90', '180']
;                       Pitch-angle state.
;
; :Returns:
;       CHANNELS:       out, required, type=objref/objarr
;                       MrMMS_EDI_Channel objects that match the search criteria.
;-
FUNCTION MrMMS_EDI_Dist::Find, $
COUNT=count, $
CHANNEL=channel, $
GDU=gdu, $
PA_STATE=pa_state
	Compile_Opt idl2
	On_Error, 2
	
	;Short-cut
	nGDU     = N_Elements(gdu)
	nChannel = N_Elements(channel)
	nPA      = N_Elements(pa_state)
	IF nGDU + nChannel + nPA EQ 0 THEN BEGIN
		channels = self -> Get(/ALL, COUNT=count)
		RETURN, channels
	ENDIF
	
	;Defaults
	IF nGDU     EQ 0 THEN gdu      = [1,2,3]
	IF nChannel EQ 0 THEN channel  = [1,2,3,4]
	IF nPA      EQ 0 THEN pa_state = [0, 90, 180]
	
	;Get all objects
	allChannels = self -> Get(/ALL, COUNT=nChannels)
	
	;Allocate memory for the found channels
	outChannels = ObjArr(nChannels)
	
	;Loop over each channel
	count = 0
	FOREACH chan, allChannels DO BEGIN
		IF MrIsMember( gdu,      chan.gdu      ) && $
		   MrIsMember( channel,  chan.channel  ) && $
		   MrIsMember( pa_state, chan.pa_state ) $
		THEN BEGIN
			outChannels[count] = chan
			count += 1
		ENDIF
	ENDFOREACH
	
	;Trim results
	CASE count OF
		0:    outChannels = Obj_New()
		1:    outChannels = outChannels[0]
		ELSE: outChannels = outChannesl[0:count-1]
	ENDCASE
	
	RETURN, outChannels
END


;+
;   Retrieve data.
;
; :Params:
;       OCOUNTS:        out, required, type=MrScalarTS object
;                       Electron counts.
;       OPHI:           out, required, type=MrScalarTS object
;                       Azimuthal angle of electron incident trajectories.
;       OTHETA:         out, required, type=MrScalarTS object
;                       Polar angle of electron incident trajectories.
;
; :Keywords:
;       CHANNEL:        in, optional, type=int, default=[1,2,3,4]
;                       Channel from which meaurments were made.
;       GDU:            in, optional, type=int, default=[1,2]
;                       GDU that took the measurement.
;       PA_STATE:       in, optional, type=int, default=['q0', '0', '90', '180']
;                       Pitch-angle state.
;-
PRO MrMMS_EDI_Dist::GetData, oCounts, oPhi, oTheta, $
CHANNEL=channel, $
GDU=gdu, $
PA_STATE=pa_state
	Compile_Opt idl2
	On_Error, 2
	
	;Get the desired channels
	channels = self -> Find( COUNT    = count, $
	                         CHANNEL  = channel, $
	                         GDU      = gdu, $
	                         PA_STATE = pa_state )

;-------------------------------------------
; Concatenate Data /////////////////////////
;-------------------------------------------
	FOREACH chan, channels, idx DO BEGIN
		
		;FAC
		IF self.fac NE 'NONE' THEN BEGIN
			;Define the time variable
			IF idx EQ 0 THEN oTime = MrTimeVar()
			
			;Update the FAC if the times have changed
			IF ~chan.counts -> IsTimeIdentical(oTime) THEN BEGIN
				oT = MrVar_FAC( self.oPar, self.oPerp, self.fac, $
				                TIME = chan.counts['TIMEVAR'] )
			ENDIF
		ENDIF

		;Grab the angles
		chan -> GetAngles, oPhi_temp, oTheta_temp, TMATRIX=oT
		
		;Concatenate variables
		IF idx EQ 0 THEN BEGIN
			time   = chan.counts['TIME', 'TT2000']
			counts = chan.counts['DATA']
			phi    = oPhi_temp['DATA']
			theta  = oTheta_temp['DATA']
		ENDIF ELSE BEGIN
			time   = [ time,   chan.counts['TIME', 'TT2000'] ]
			counts = [ counts, chan.counts['DATA']           ]
			phi    = [ phi,    oPhi_temp['DATA']             ]
			theta  = [ theta,  oTheta_temp['DATA']           ]
		ENDELSE
	ENDFOREACH

;-------------------------------------------
; Create MrVariables ///////////////////////
;-------------------------------------------
	;First, sort the data
	iSort  = Sort(time)
	time   = time[iSort]
	counts = counts[iSort]
	phi    = phi[iSort]
	theta  = theta[iSort]
	
	;Create variables
	oCounts = MrScalarTS( time, Temporary(counts), T_TYPE='TT2000', NAME='mms_edi_counts_agregate')
	oPhi    = MrScalarTS( time, Temporary(phi),    T_TYPE='TT2000', NAME='mms_edi_phi_agregate')
	oTheta  = MrScalarTS( time, theta, /NO_COPY,   T_TYPE='TT2000', NAME='mms_edi_theta_agregate')

	;Copy attributes over
	chan.counts -> CopyAttrTo, oCounts
	oPhi_temp   -> CopyAttrTo, oPhi
	oTheta_temp -> CopyAttrTo, oTheta
END


;+
;   Create a distribution as a function of time, azimuth angle, and polar angle.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the distribution is added to the variable cache.
;       DPHI:           in, optional, type=float, default=1.5
;                       Width of each azimuthal bin (degrees). Cannot be used with `NPHI_BINS`.
;       DTHETA:         in, optional, type=float, default=0.1
;                       Width of each polar angle bin (degrees). Cannot be used with `NTHETA_BINS`.
;       DTIME:          in, optional, type=float, default=1.0
;                       Duration, in seconds, of each time bin.
;       GDU1:           in, optional, type=boolean
;                       If set, only data from GDU1 is returned. The default is to
;                           combine data from both GDU1 and GDU2 into a single array.
;                           Cannot be used with `GDU2`.
;       GDU2:           in, optional, type=boolean
;                       If set, only data from GDU2 is returned. The default is to
;                           combine data from both GDU1 and GDU2 into a single array.
;                           Cannot be used with `GDU1`.
;       NAME:           in, optional, type=string, default=sc_instr_dist_optdesc_mode_level
;                       Name to be given to the distribution.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of bins spanning `PHI_RANGE`. Cannot be used with `DPHI`.
;       NTHETA_BINS:    in, optional, type=integer
;                       Number of bins spanning `THETA_RANGE`. Cannot be used with `DTHETA`.
;       PHI_RANGE:      in, optional, type=fltarr(2)
;                       Minimum and maximum azimuthal angles over which to bin (degrees).
;       THETA_RANGE:    in, optional, type=fltarr(2)
;                       Minimum and maximum polar angles over which to bin (degrees).
;       TYPE_FAC:       in, optional, type=string
;                       Name of the field-aligned coordinate system into which data
;                           is rotated. Options are: {'None' | '' | 'ExB' | 'VxB'}
;
; :Returns:
;       ODIST:          out, required, type=MrTimeSeries
;                       The distribution of counts as a function of time, azimuth angle,
;                           and polar angle.
;-
function MrMMS_EDI_Dist::Dist, $
CACHE=cache, $
DPHI=dphi, $
DTIME=dt, $
DTHETA=dtheta, $
GDU1=gdu1, $
GDU2=gdu2, $
NAME=name, $
NPHI_BINS=nPhi_bins, $
NTHETA_BINS=nTheta_bins, $
PHI_RANGE=phi_range, $
THETA_RANGE=theta_range, $
TYPE_FAC=type_fac
	compile_opt idl2
	
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return, !Null
	endif
	
	;Default name
	if n_elements(name) eq 0 then name = strjoin([self.sc, 'edi', 'dist', self.optdesc, self.mode, self.level], '_')
	
	;Get the data
	self -> GetData, oCounts, oPhi, oTheta, $
	                 GDU1     = gdu1, $
	                 GDU2     = gdu2, $
	                 TYPE_FAC = type_fac
	
	;Bin the data
	oDist = self -> Bin( oCounts, oPhi, oTheta, $
	                     DPHI        = dphi, $
	                     DTIME       = dt, $
	                     DTHETA      = dtheta, $
	                     NPHI_BINS   = nPhi_bins, $
	                     NTHETA_BINS = nTheta_bins, $
	                     NTIME_BINS  = nt_bins, $
	                     PHI_RANGE   = phi_range, $
	                     THETA_RANGE = theta_range )
	
	;Set the name and cache variable
	oDist -> SetName, name
	if keyword_set(cache) then oDist -> Cache
	
	;Return the distribution
	return, oDist
end


;+
;   Create a distribution as a function of time and azimuth angle. Data is averaged
;   over all polar angles.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the distribution is added to the variable cache.
;       DPHI:           in, optional, type=float, default=1.5
;                       Width of each azimuthal bin (degrees). Cannot be used with `NPHI_BINS`.
;       DTIME:          in, optional, type=float, default=1.0
;                       Duration, in seconds, of each time bin.
;       NAME:           in, optional, type=string, default=sc_instr_phidist_optdesc_mode_level
;                       Name to be given to the distribution.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of bins spanning `PHI_RANGE`. Cannot be used with `DPHI`.
;       PHI_RANGE:      in, optional, type=fltarr(2)
;                       Minimum and maximum azimuthal angles over which to bin (degrees).
;       THETA_RANGE:    in, optional, type=fltarr(2)
;                       Minimum and maximum polar angles over which to bin (degrees).
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by the ::GetData method is also accepted here.
;
; :Returns:
;       ODIST:          out, required, type=MrTimeSeries
;                       The distribution of counts as a function of time and azimuth angle,
;                           averaged over polar angle.
;-
function MrMMS_EDI_Dist::PhiSpec, $
CACHE=cache, $
DPHI=dPhi, $
DTIME=dt, $
NAME=name, $
NPHI_BINS=nPhi_bins, $
PHI_RANGE=phi_range, $
THETA_RANGE=theta_range, $
_REF_EXTRA=extra
	compile_opt idl2
	on_error, 2
	
	;Defaults
	nTheta_bins = 1
	if n_elements(name) eq 0 then name = strjoin([self.sc, 'edi', 'phispec', self.optdesc, self.mode, self.level], '_')
	
	;Get the data
	self -> GetData, oCounts, oPhi, oTheta, $
	                 _STRICT_EXTRA = extra
	
	;Create the distribution
	oPhiDist = self -> Bin( oCounts, oPhi, oTheta, $
	                        CACHE       = cache, $
	                        DPHI        = dphi, $
	                        DTIME       = dt, $
	                        NAME        = name, $
	                        NPHI_BINS   = nPhi_bins, $
	                        NTHETA_BINS = nTheta_bins, $
	                        PHI_RANGE   = phi_range, $
	                        THETA_RANGE = theta_range )
	
	return, oPhiDist
end


;+
;   Create a distribution as a function of time and polar angle. Data is averaged
;   over all azimuthal angles.
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, the distribution is added to the variable cache.
;       DTHETA:         in, optional, type=float, default=0.1
;                       Width of each polar angle bin (degrees). Cannot be used with `NTHETA_BINS`.
;       DTIME:          in, optional, type=float, default=1.0
;                       Duration, in seconds, of each time bin.
;       NAME:           in, optional, type=string, default=sc_instr_thetadist_optdesc_mode_level
;                       Name to be given to the distribution.
;       NTHETA_BINS:    in, optional, type=integer
;                       Number of bins spanning `THETA_RANGE`. Cannot be used with `DTHETA`.
;       PHI_RANGE:      in, optional, type=fltarr(2)
;                       Minimum and maximum azimuthal angles over which to bin (degrees).
;       THETA_RANGE:    in, optional, type=fltarr(2)
;                       Minimum and maximum polar angles over which to bin (degrees).
;       _REF_EXTRA:     in, optional, type=any
;                       Any keyword accepted by the ::GetData method is also accepted here.
;
; :Returns:
;       ODIST:          out, required, type=MrTimeSeries
;                       The distribution of counts as a function of time, polar angle,
;                           averaged over azimuthal angle.
;-
FUNCTION MrMMS_EDI_Dist::ThetaSpec, $
CACHE=cache, $
DTIME=dt, $
DTHETA=dTheta, $
NAME=name, $
NTHETA_BINS=nTheta_bins, $
PHI_RANGE=phi_range, $
THETA_RANGE=theta_range, $
_REF_EXTRA=extra
	Compile_Opt idl2

	;Error handling
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, 0
	ENDIF
	
	;Defaults
	nPhi_bins = 1
	IF n_elements(name) EQ 0 THEN name = StrJoin([self.sc, 'edi', 'thetaspec', self.optdesc, self.mode, self.level], '_')
	
	;Get the data
	self -> GetData, oCounts, oPhi, oTheta, $
	                 _STRICT_EXTRA=extra

	;Create the distribution
	oThetaDist = self -> Bin( oCounts, oPhi, oTheta, $
	                          CACHE       = cache, $
	                          DTHETA      = dTheta, $
	                          DTIME       = dt, $
	                          NAME        = name, $
	                          NPHI_BINS   = nPhi_bins, $
	                          NTHETA_BINS = nTheta_bins, $
	                          PHI_RANGE   = phi_range, $
	                          THETA_RANGE = theta_range )
	
	RETURN, oThetaDist
END


;+
;   Load EDI data from source.
;
; :Params:
;       SC:             in, required, type=string
;                       MMS spacecraft identifier. Options are: {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:           in, required, type=string
;                       Telemetry mode. Options are: {'slow' | 'fast' | 'srvy' | 'brst'}
;       OPTDESC:        in, optional, type=string, defualt='q0'
;                       EDI-specific optional descriptor. Options are: {'amb' | 'amb-pm2' |
;                           'alt-cc' | 'alt-oc' | 'alt-oob' | 'alt-ooc'}
;
; :Keywords:
;       COORD_SYS:      in, required, type=string
;                       Coordinate system in which data is loaded. Options are: { 'dbcs' | 'gse' }.
;       LEVEL:          in, required, type=string, default='l2'
;                       Data quality level. Options are: { 'l1a' | 'ql' | 'l2' }
;-
PRO MrMMS_EDI_Dist::Load, sc, mode, optdesc, $
COORD_SYS=coord_sys, $
LEVEL=level
	Compile_Opt idl2
	On_Error, 2
	
	instr = 'edi'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(coord_sys) EQ 0 THEN coord_sys = 'gse'
	IF N_Elements(optdesc)   EQ 0 THEN optdesc   = ['q0', 'amb', 'amb-pm2', 'amb-alt-cc', 'amb-alt-oc', 'amb-alt-oob', 'amb-alt-ooc']
	
	;Search for available files
	filenames = MrMMS_Get_Filenames(sc, 'edi', mode, level, OPTDESC=optdesc)
	
	;Extract the optional descriptors
	MrMMS_Parse_Filename, filenames, OPTDESC=optdesc
	
	;If multiple, pick the first (unless it is efield mode)
	IF ~Array_Equal(optdesc, optdesc[0]) THEN BEGIN
		iopt = Where(optdesc NE 'efield', nopt)
		IF nopt EQ 0 THEN BEGIN
			Message, 'No ambient or Q0 data found.'
		ENDIF ELSE BEGIN
			optdesc = optdesc[iopt[0]]
			MrPrintF, 'LogWarn', 'Multiple optional descriptors found. Choosing "' + optdesc + '".'
		ENDELSE
	ENDIF ELSE BEGIN
		optdesc = optdesc[0]
	ENDELSE
	
	;Load the data
	MrMMS_Load_Data, sc, instr, mode, level, $
	                 OPTDESC   = optdesc, $
;	                 VARFORMAT = ['*flux*', '*counts*', '*traj_' + coord_sys + '_gdu?_' + mode + '*'], $
	                 VARNAMES  = varnames

	;Add variables to the cache
	FOR i = 0, N_ELEMENTS(varnames) - 1 DO BEGIN
		;Search for counts
		parts = StRegEx(varnames[i], '(counts|flux)([1-4])?_(gdu1|gdu2|0|90|180)?_'+mode, /EXTRACT, /SUBEXP)
		IF parts[0] EQ '' THEN CONTINUE
		
		;Find corresponding trajectory
		tf_traj = StRegEx(varnames, 'traj' + parts[2] + '_' + coord_sys + '_' + parts[3], /BOOLEAN)
		iTraj   = Where(tf_traj, nTraj)
		IF nTraj EQ 0 THEN BEGIN
			MrPrintF, 'LogWarn', 'No matching trajectories found for: "' + varnames[i] + '".'
			CONTINUE
		ENDIF
		
		;Need:
		;  - Channel
		;  - GDU
		;  - PA State
		channel  = parts[2]
		IF channel EQ '' THEN BEGIN
			channel  = 1
			gdu      = Fix(StrMid(parts[3], 3, 1))
			pa_state = 90
		ENDIF ELSE BEGIN
			gdu      = 3B
			pa_state = Fix(parts[3])
		ENDELSE
		
		;Add the variable
		oChannel = MrMMS_EDI_Channel(varnames[i], varnames[iTraj], channel, gdu, pa_state)
		self -> Add, oChannel
	ENDFOR
	
	;Set properties
	self.sc        = sc
	self.instr     = instr
	self.mode      = mode
	self.level     = level
	self.optdesc   = N_Elements(optdesc) EQ 0 ? '' : optdesc
	self.coord_sys = coord_sys
END


;+
;   Set the field-aligned coordinate system.
;
;   Calling Sequence:
;       oEDI -> SetFAC, par
;       oEDI -> SetFAC, par, perp, fac
;
; :Params:
;       PAR:            in, required, type=string/integer/objref
;                       The name, number, or object reference of a MrVectorTS variable
;                           that defines the parallel direction of a field-aligned
;                           coordinate system..
;       PERP:           in, optional, type=string
;                       The name, number, or object reference of a MrVectorTS variable
;                           that defines the perpendicular direction of a field-aligned
;                           coordinate system.
;       FAC:            in, optional, type=string
;                       Type of field-aligned coordinate system to use.
;-
PRO MrMMS_EDI_Dist::SetFAC, par, perp, fac
	Compile_Opt idl2
	On_Error, 2
	
	IF N_Elements(par)  GT 0 THEN self.oPar  = MrVar_Get(par)
	IF N_Elements(perp) GT 0 THEN self.oPerp = MrVar_Get(perp)
	IF N_Elements(fac)  GT 0 THEN self.fac   = StrUpCase(fac)
END


;+
;   Set the counts and trajectories data (also, see the ::Load method).
;
; :Params:
;       COUNTS:         in, required, type=string/int/objref
;                       The name, number, or object reference of a MrScalarTS variable
;                          containing the counts or flux measurements
;       TRAJ:           in, required, type=Nx2 fltarr
;                       The name, number, or object reference of a MrScalarTS variable
;                          containing the trajectory of incoming electrons that registar
;                          as `COUNTS`
;       CHANNEL:        in, optional, type=int
;                       Channel from which meaurments were made. Options are: { 1 | 2 | 3 | 4 }
;       GDU:            in, optional, type=int
;                       GDU that took the measurement. Options are: { 1 | 2 }
;       PA_STATE:       in, optional, type=int
;                       Pitch-angle state. Options are: { 'q0' | '0' | '90' | '180' }
;-
PRO MrMMS_EDI_Dist::SetData, counts, traj, channel, gdu, pa_state
	Compile_Opt idl2
	On_Error, 2
	
	;Number of variables to add
	nVars = N_Elements(counts)
	FOR i = 0, nVars - 1 DO BEGIN
		;Create the channel variable
		theChan = MrMMS_EDI_Channel(counts[i], traj[i], channel[i], gdu[i], pa_state[i])
		IF ~Obj_Valid(theChan) THEN CONTINUE
		
		;Add the channel to the container
		self -> Add, theChan
	ENDFOR
END


;+
;   The class definition statement.
;
; :Params:
;       CLASS:          out, optional, type=structure
;-
PRO MrMMS_EDI_Dist__DEFINE
	Compile_Opt idl2
	
	class = { MrMMS_EDI_Dist, $
	          Inherits IDL_Container, $
	          sc:        '', $
	          instr:     '', $
	          mode:      '', $
	          level:     '', $
	          optdesc:   '', $
	          coord_sys: '', $
	          oPar:      Obj_New(), $
	          oPerp:     Obj_New(), $
	          fac:       '' $
	        }
END