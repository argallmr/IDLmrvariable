; docformat = 'rst'
;
; NAME:
;   MrDist3D__Define
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
;   Calculate moments and reduced distributions from a 3D distribution function.
;
; :Categories:
;   MrDist
;
; :See Also:
;   MrDist4D__Define.pro
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
;       2016/10/24  -   Written by Matthew Argall
;-
;*****************************************************************************************
;+
;   The initialization method.
;
; :Params:
;       TIME:           in, required, type=NxM array
;                       Name or reference of a MrTimeVar object, or an array
;                           of time stamps. If a name is provided, the assiciated
;                           variable must exist in the variable cache.
;       DATA:           in, required, type=NxM array
;                       Name or reference of a MrVariable object, or the dependent
;                           variable data. If a name is given, the associated variable
;                           must exist in the variable cache. 
;
; :Keywords:
;       CACHE:          in, optional, type=boolean, default=0
;                       If set, both `TIME` and `DATA` are added to the variable cache.
;       ELEVATION:      in, optional, type=boolean, default=0
;                       If set, `THETA` is taken to be the elevation angle. By default,
;                           it is interpreted as the polar angle.
;       DIMENSION:      in, optional, type=integer
;                       The time-dependent, 1-based dimension of `DATA`. If not provided,
;                           the dimension of `DATA` that is equal in size to `TIME` is
;                           chose as the default.
;       NAME:           in, optional, type=integer
;                       Name to be given to the variable object.
;       NO_CLOBBER:     in, optional, type=boolean, default=0
;                       If set, do not clobber variables of the same name. Instead,
;                           rename this variable by appending "_#" to `NAME`, where
;                           "#" represents a unique number. Ignored unless `CACHE` is set.
;       NO_COPY:        in, optional, type=boolean, default=0
;                       If set `DATA` will be copied directly into the object
;                           and will be left undefined (a MrTimeSeries object will not
;                           be destroyed, but its array will be empty).
;       T_TYPE:         in, optional, type=integer
;                       If `TIME` is an array of time stamps, use this keyword to indicate
;                           the format or time-basis. See MrTimeVar for more details.
;       T_NAME:         in, optional, type=integer
;                       Name to be given to the MrTimeVar object. Ignored unless `TIME`
;                           is an array of time stamps.
;-
function MrDist3D::INIT, dist3D, phi, theta, energy, Vsc, $
ELEVATION=elevation, $
DEGREES=degrees, $
MASS=mass, $
SPECIES=species, $
UNITS=units, $
_REF_EXTRA=extra
	compile_opt idl2

	;Error handling
	catch, the_error
	if the_error ne 0 then begin
		catch, /CANCEL
		MrPrintF, 'LogErr'
		return, 0
	endif
	
	;Allocate heap
	self.distfn  = ptr_new(/ALLOCATE_HEAP)
	self.dphi    = ptr_new(/ALLOCATE_HEAP)
	self.dtheta  = ptr_new(/ALLOCATE_HEAP)
	self.denergy = ptr_new(/ALLOCATE_HEAP)
	self.phi     = ptr_new(/ALLOCATE_HEAP)
	self.theta   = ptr_new(/ALLOCATE_HEAP)
	self.energy  = ptr_new(/ALLOCATE_HEAP)
	self.Vsc     = ptr_new(/ALLOCATE_HEAP)
	
	;Store data
	if n_elements(dist3d) gt 0 then begin
		self -> SetData2, dist3D, phi, theta, energy, Vsc, $
		                  DEGREES       = degrees, $
		                  UNITS         = units, $
		                  RADIANS       = radians, $
		                  _STRICT_EXTRA = extra
	endif
	
	;Set object properties
	self -> SetProperty, ELEVATION = elevation, $
	                     MASS      = mass, $
	                     SPECIES   = species

	return, 1
end


;+
;   Clean up after the object is destroyed
;-
pro MrDist3D::CLEANUP
	compile_opt idl2
	on_error, 2
	
	ptr_free, self.distfn
	ptr_free, self.denergy
	ptr_free, self.dphi
	ptr_free, self.dtheta
	ptr_free, self.energy
	ptr_free, self.phi
	ptr_free, self.theta
	ptr_free, self.Vsc
end


;+
;   Find the spacing of energy bins assuming constant dE/E
;
; :Returns:
;       DE:             Spacing of the energy bins.
;-
function MrDist3D::DeltaE
	compile_opt idl2
	on_error, 2

	;
	; dLogE is the spacing of the energy bins in log space
	;
	;   dLogE = log(E1) - log(E0)                      (1)
	;
	; dE/E is just the derivative of a natural log
	;
	;   d( ln(E) ) = dE / E                            (2)
	;
	; If we change bases by log_a(B) = log_x(B) / log_x(A), where x=10, A=exp, B=E
	;
	;   d( ln(E) ) = d [ log(E) / log(exp) ]
	;              = dLogE / log(exp)                  (3)
	;
	; Combining Equations 1 and 2
	;
	;   dE/E = dLogE / log(exp)
	;
	;     dE = E * dLogE / log(exp)
	;
	
	dLogE = alog10( (*self.energy)[1] ) - alog10( (*self.energy)[0] )
	dE    = *self.energy * dLogE / alog10( exp(1) )
	
	;Return results
	return, dE
end


;+
;   Compute the 0th order moment, denstiy.
;
;       n = \sum_{i} \sum_{j} \sum_{k} f v_{i}^{2} \sin{\theta_{k}} \Delta v_{i} \Delta \phi_{j} \Delta \theta_{k}
;
; :Returns:
;       N:          out, required, type=float
;                   Density (1/cm^3) of particles calculated from the distribution.
;-
FUNCTION MrDist3D::Density
	Compile_Opt idl2
	On_Error, 2
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Constants
	deg2rad = !dpi / 180D
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	dims    = Size(*self.distfn, /DIMENSIONS)

;-----------------------------------------------------
; Integration Bins and their Sizes \\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Phi
	dphi  = Median( (*self.phi)[1:*,*] - (*self.phi)[0:-2,*]) * deg2rad
	phi   = Rebin( *self.phi, dims ) * deg2rad
	
	;Theta
	theta  = Rebin(Reform((*self.theta)[0,*]), dims[1], dims[2]) * deg2rad
	dTheta = Median( theta[1:-1,*] - theta[0:-2,*] )
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	; Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * *self.energy / self.mass)
	dE   = self -> DeltaE()
	dv   = eV2J * dE / (self.mass * v)
	
	;SPACECRAFT POTENTIAL
	IF N_Elements(*self.Vsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		signQ = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc   = Sqrt( 2.0 * q * (*self.Vsc) / self.mass )
		VM    = v^2 * (1 + signQ * (vsc^2 / v^2))
	ENDIF ELSE BEGIN
		VM = v^2
	ENDELSE

;-----------------------------------------------------
; Density \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Integrate phi
	ftemp = Total( *self.distfn, 1) * dPhi

	;Integrate theta
	ftemp  = Total( ftemp * Sin(theta), 1) * dTheta
	
	;Integrate velocity
	;   - V is in m/s while f is in s^3/cm^6
	;   - v^2 * dv --> (m/s)^3  --> (cm/s)^3 * 1e6
	N = Total( ftemp * v * Sqrt(VM) * dv ) * 1e6
	
	RETURN, N
END


;+
;   Compute entropy via Boltzmann's H-theorem.
;
;       H = \integral f \ln(f) d^{3}v
;
;   References
;       Kaufmann, R. L., and W. R. Paterson (2009), Boltzmann H function and entropy in
;           the plasma sheet, J. Geophys. Res. Sp. Phys., 114(A9), doi:10.1029/2008JA014030.
;
; :Returns:
;       H:          out, required, type=float
;                   Density of particles calculated from the distribution.
;-
FUNCTION MrDist3D::Entropy
	Compile_Opt idl2
	On_Error, 2
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Convert to radians
	deg2rad = !dpi / 180D
	dims    = Size(*self.distfn, /DIMENSIONS)
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	tf_vsc  = N_Elements(*self.Vsc) GT 0
	
	;Check for zero values -- bad for ALog()
	tf_nan = 0B
	ftemp  = *self.distfn
	iZero  = Where(*self.distfn EQ 0, nZero)
	IF nZero GT 0 THEN BEGIN
		ftemp[iZero] = !Values.F_NaN
		tf_nan = 1B
	ENDIF

;-----------------------------------------------------
; Integration Bins and their Sizes \\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;SPACECRAFT POTENTIAL
	IF tf_vsc THEN BEGIN
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q * (*self.Vsc) / self.mass )
	ENDIF
	
	;PHI
	dphi  = Median( (*self.phi)[1:*,*] - (*self.phi)[0:-2,*]) * deg2rad
	phi   = Rebin( *self.phi, dims ) * deg2rad

	;THETA
	theta  = Rebin(Reform((*self.theta)[0,*]), dims[1], dims[2]) * deg2rad
	dTheta = Median( theta[1:-1,*] - theta[0:-2,*] )
	
	;VELOCITY
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	; Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * *self.energy / self.mass)
	dE   = self -> DeltaE()
	dv   = eV2J * dE / (self.mass * v)

;-----------------------------------------------------
; Entropy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Integrate Phi
	ftemp = Total( ftemp * ALog(ftemp), 1, NAN=tf_nan) * dPhi
	
	;Integrate Theta
	ftemp  = Total( ftemp * Sin(theta), 1, NAN=tf_nan) * dTheta
	
	;Integrate velocity
	;   - Sign is charge dependent. E = q * V = 1/2 * m * v^2
	;   - V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6  --> s^3/m^6 * 1e12
	IF tf_Vsc $
		THEN H = Total( ftemp * v * Sqrt(v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e12 $
		ELSE H = Total( ftemp * v^2 * dv, NAN=tf_nan ) * 1e12
	
	;Entropy:
	;   - S = -kH
	S = -MrConstants('k_B') * H
	
	RETURN, S
END


;+
;   Reduce a single 3D distribution function to 1D by averaging over
;   energy and pitch angle.
;
; :Params:
;       ENERGY:         out, optional, type=fltarr
;                       Energy bin centers of the distribution.
;       DE:             out, optional, type=fltarr
;                       Half-width of the energy bins.
;
; :Keywords:
;       NE_BINS:        in, optional, type=integer
;                       Number of energy bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       PHI_RANGE:      in, optional, type=fltarr(2), default=[0.0\, 360.0]
;                       The range in azimuthal angle (degrees) over which to average.
;       THETA_RANGE:    in, optional, type=fltarr(2), default=[0.0\, 180.0]
;                       The range in polar angle (degrees) over which to average.
;
; :Returns:
;       DIST1D:         out, required, type=fltarr(L)
;                       The 1D distribution with size nEnergy.
;-
function MrDist3D::ESpec, energy, dE, $
PHI_RANGE=phi_range, $
NE_BINS=nE_bins, $
THETA_RANGE=theta_range
	compile_opt idl2
	on_error, 2

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Dimension sizes
	dims    = size(*self.distFn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

	;PA range over which to average
	if n_elements(phi_range)   eq 0 then phi_range   = [0.0, 360.0]
	if n_elements(theta_range) eq 0 then theta_range = [0.0, 180.0]
	if n_elements(nE_bins)     gt 0 then MrPrintF, 'LogWarn', 'TODO: Allow NE_BINS to vary.'

	;Use original number of energy bins
	;   - Do not use NE_BINS to prevent cyclic messages there ---^
	nBins = nEnergy

;-----------------------------------------------------
; Coordinate Space \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Combine data into 3-vectors of [phi, theta, energy]
	coords = transpose( [ [reform( *self.phi,   nPhi*nTheta)], $
	                      [reform( *self.theta, nPhi*nTheta)] ] )

;-----------------------------------------------------
; Weight \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	;   - w = v^2 * sin(theta) * dv * dTheta * dPhi
	;   - Applied as total(w*psd)/total(w)
	;   - dTheta and dPhi are constants, so will factor out
	;   - v^2 and dv are constants, because we are not avaraging over different energies.
	;   - Elements from two different old theta bins can fall into the same
	;     new theta bin, so the sin(theta) term does not cancel.
	;   => w = sin(theta)
	;
	
	;Weight
	if self.elevation $
		then weight = cos( *self.theta * !dtor ) $
		else weight = sin( *self.theta * !dtor )

;-----------------------------------------------------
; Re-bin Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Locate data within new bins
	;   - Keep the number of Theta bins the same
	;   - Average all relevant PHI and ENERGY values.
	cHist  = hist_nd( coords, $
	                  MIN             = [phi_range[0], theta_range[0]], $
	                  MAX             = [phi_range[1], theta_range[1]], $
	                  NBINS           = [           1,              1], $
	                  REVERSE_INDICES = ri )

;-----------------------------------------------------
; Reduce to 1D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Allocate memory to reduced 2D distribution
	dist1D = fltarr(nBins)

	;Loop over bins
	k = 0
	for i = 0, nEnergy - 1 do begin
		;Skip empty bins
		if ri[k] eq ri[k+1] then continue

		;Source indices
		inds = ri[ri[k]:ri[k+1]-1]
;		isrc = array_indices([nPhi,nTheta], inds, /DIMENSIONS)
	
		;Weight
;		w = weight[isrc[0,*], isrc[1,*]]

		;Re-bin
		temp      = (*self.distFn)[*, *, i]
;		dist1D[i] = total( temp[isrc[0,*], isrc[1,*]] * w ) / total(w)
		dist1D[i] = total( temp[inds] * weight[inds]) / total( weight[inds] )
	endfor

;-----------------------------------------------------
; Output Bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Create bins
	energy = *self.energy

	;Return the 1D distribution
	return, dist1D
end


;+
;   Compute the 3rd order moment, heat flux.
;
; :Returns:
;       Q:          out, required, type=float
;                   Heat flux (mW/m^2) of particles calculated from the distribution.
;-
FUNCTION MrDist3D::HeatFlux2
	Compile_Opt idl2
	On_Error, 2

	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Constants
	deg2rad = !dpi / 180D
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	dims    = Size(*self.distfn, /DIMENSIONS)

;-----------------------------------------------------
; Integration Bins and their Sizes \\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Phi
	dphi = Median( (*self.phi)[1:*,*] - (*self.phi)[0:-2,*]) * deg2rad
	phi  = Rebin( *self.phi, dims ) * deg2rad
	
	;Theta
	theta  = Rebin(Reform((*self.theta)[0,*]), dims[1], dims[2]) * deg2rad
	dTheta = Median( theta[1:-1,*] - theta[0:-2,*] )
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	;Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * *self.energy / self.mass)
	dE   = self -> DeltaE()
	dv   = eV2J * dE / (self.mass * v)
	
	;BULK VELOCITY
	;   - We need to be in the rest frame of the plasma
	;   - Convert BULKV from km/s to m/s
	bulkv    = 1e3 * self -> Velocity()
	bulkvmag = Sqrt( Total(bulkv^2) )
	
	;SPACECRAFT POTENTIAL
	IF N_Elements(*self.Vsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		signQ = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc   = Sqrt( 2.0 * q * (*self.Vsc) / self.mass )
		VM    = v^2 * (1 + signQ * (vsc^2 / v^2))
	ENDIF ELSE BEGIN
		VM = v^2
	ENDELSE

;-----------------------------------------------------
; Integration Bins and their Sizes \\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;
	; Let
	;    VM = Vm^2 - sign(q) Vsc^2
	;    v  = Sqrt( VM )
	;    dv = Vm / Sqrt( VM ) * dVm
	;
	; Where, here and in what follows,
	;    VM      =  The velocity corrected for spacecraft potential effects
	;    Vm      =  The measured velocity related to the energy bins of the instrument
	;    Vsc     =  The velocity gained by passing through the spacecraft potential sheath
	;    U       =  Plasma bulk velocity
	;    dVm     =  Size of a measured velocity space bin (energy bin width converted to velocity)
	;    dOmega  =  Sin(Theta) dTheta dPhi = Element of solid angle
	;    t       =  Theta, the polar angle
	;    p       =  Phi, the azimuthal angle
	;    m       =  Particle mass
	;    f       =  The distribution function
	;    q       =  Charge of the particle
	;
	; It follows that
	;
	;    H  = m/2 * \Integral{ f v^2 dv dOmega [ (\vec{v} - \vec{V}) * |\vec{v} - \vec{V}|^2 ] }
	;    H  = m/2 * \Integral{ f v^2 dv dOmega (\vec{v} - \vec{V}) (v^2 - 2 * v \cdot V + V^2) }
	;    Hx = m/2 * \Integral{ f v^2 dv dOmega (vx - Vx) (v^2 - 2 * vx * Vx + V^2) }
	;    Hy = m/2 * \Integral{ f v^2 dv dOmega (vy - Vy) (v^2 - 2 * vy * Vy + V^2) }
	;    Hz = m/2 * \Integral{ f v^2 dv dOmega (vz - Vz) (v^2 - 2 * vz * Vz + V^2) }
	;
	
	
	;
	; PHI
	;
	
	;Hx
	hx0 = Total( *self.distfn * Cos(phi),            1 ) * dPhi
	hx1 = Total( *self.distfn,                       1 ) * dPhi * bulkv[0]
	hx2 = Total( *self.distfn * Cos(phi),            1 ) * dPhi * bulkvmag^2
	hx3 = Total( *self.distfn,                       1 ) * dPhi * bulkv[0] * bulkvmag^2
	hx4 = Total( *self.distfn * Cos(phi)^2,          1 ) * dPhi * bulkv[0]
	hx5 = Total( *self.distfn * Sin(phi) * Cos(phi), 1 ) * dPhi * bulkv[1]
	hx6 = Total( *self.distfn * Cos(phi),            1 ) * dPhi * bulkv[2]
	hx7 = Total( *self.distfn * Cos(phi),            1 ) * dPhi * bulkv[0]^2
	hx8 = Total( *self.distfn * Sin(phi),            1 ) * dPhi * bulkv[0] * bulkv[1]
	hx9 = Total( *self.distfn,                       1 ) * dPhi * bulkv[0] * bulkv[2]
	
	;Hy
	hy0 = Total( *self.distfn * Sin(phi),            1 ) * dPhi
	hy1 = Total( *self.distfn,                       1 ) * dPhi * bulkv[1]
	hy2 = Total( *self.distfn * Sin(phi),            1 ) * dPhi * bulkvmag^2
	hy3 = Total( *self.distfn,                       1 ) * dPhi * bulkv[1] * bulkvmag^2
	hy4 = Total( *self.distfn * Sin(phi) * Cos(phi), 1 ) * dPhi * bulkv[0]
	hy5 = Total( *self.distfn * Sin(phi)^2,          1 ) * dPhi * bulkv[1]
	hy6 = Total( *self.distfn * Sin(phi),            1 ) * dPhi * bulkv[2]
	hy7 = Total( *self.distfn * Cos(phi),            1 ) * dPhi * bulkv[0] * bulkv[1]
	hy8 = Total( *self.distfn * Sin(phi),            1 ) * dPhi * bulkv[1]^2
	hy9 = Total( *self.distfn,                       1 ) * dPhi * bulkv[1] * bulkv[2]
	
	;Hz
	hz0 = Total( *self.distfn,                       1 ) * dPhi
	hz1 = Total( *self.distfn,                       1 ) * dPhi * bulkv[2]
	hz2 = Total( *self.distfn,                       1 ) * dPhi * bulkvmag^2
	hz3 = Total( *self.distfn,                       1 ) * dPhi * bulkv[2] * bulkvmag^2
	hz4 = Total( *self.distfn * Cos(phi),            1 ) * dPhi * bulkv[0]
	hz5 = Total( *self.distfn * Sin(phi),            1 ) * dPhi * bulkv[1]
	hz6 = Total( *self.distfn,                       1 ) * dPhi * bulkv[2]
	hz7 = Total( *self.distfn * Cos(phi),            1 ) * dPhi * bulkv[0] * bulkv[2]
	hz8 = Total( *self.distfn * Sin(phi),            1 ) * dPhi * bulkv[1] * bulkv[2]
	hz9 = Total( *self.distfn,                       1 ) * dPhi * bulkv[2]^2
	
	;
	; THETA
	;
	
	;Hx
	hx0 = Total( hx0 * Sin(theta)^2,              1 ) * dTheta
	hx1 = Total( hx1 * Sin(theta),                1 ) * dTheta
	hx2 = Total( hx2 * Sin(theta)^2,              1 ) * dTheta
	hx3 = Total( hx3 * Sin(theta),                1 ) * dTheta
	hx4 = Total( hx4 * Sin(theta)^3,              1 ) * dTheta
	hx5 = Total( hx5 * Sin(theta)^3,              1 ) * dTheta
	hx6 = Total( hx6 * Sin(theta)^2 * Cos(theta), 1 ) * dTheta
	hx7 = Total( hx7 * Sin(theta)^2,              1 ) * dTheta
	hx8 = Total( hx8 * Sin(theta)^2,              1 ) * dTheta
	hx9 = Total( hx9 * Sin(theta)   * Cos(theta), 1 ) * dTheta
	
	;Hy
	hy0 = Total( hy0 * Sin(theta)^2,              1 ) * dTheta
	hy1 = Total( hy1 * Sin(theta),                1 ) * dTheta
	hy2 = Total( hy2 * Sin(theta)^2,              1 ) * dTheta
	hy3 = Total( hy3 * Sin(theta),                1 ) * dTheta
	hy4 = Total( hy4 * Sin(theta)^3,              1 ) * dTheta
	hy5 = Total( hy5 * Sin(theta)^3,              1 ) * dTheta
	hy6 = Total( hy6 * Sin(theta)^2 * Cos(theta), 1 ) * dTheta
	hy7 = Total( hy7 * Sin(theta)^2,              1 ) * dTheta
	hy8 = Total( hy8 * Sin(theta)^2,              1 ) * dTheta
	hy9 = Total( hy9 * Sin(theta)   * Cos(theta), 1 ) * dTheta
	
	;Hz
	hz0 = Total( hz0 * Sin(theta)   * Cos(theta), 1 ) * dTheta
	hz1 = Total( hz1 * Sin(theta),                1 ) * dTheta
	hz2 = Total( hz2 * Sin(theta)   * Cos(theta), 1 ) * dTheta
	hz3 = Total( hz3 * Sin(theta),                1 ) * dTheta
	hz4 = Total( hz4 * Sin(theta)^2 * Cos(theta), 1 ) * dTheta
	hz5 = Total( hz5 * Sin(theta)^2 * Cos(theta), 1 ) * dTheta
	hz6 = Total( hz6 * Cos(theta)^3,              1 ) * dTheta
	hz7 = Total( hz7 * Sin(theta)^2,              1 ) * dTheta
	hz8 = Total( hz8 * Sin(theta)^2,              1 ) * dTheta
	hz9 = Total( hz9 * Sin(theta)   * Cos(theta), 1 ) * dTheta
	
	;
	; VELOCITY
	;
	
	
	;Hx
	;   - V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
	hx0 = 0.5 * self.mass * Total( hx0 * v * VM^2.0       * dv ) * 1e12
	hx1 = 0.5 * self.mass * Total( hx1 * v * VM^(3.0/2.0) * dv ) * 1e12
	hx2 = 0.5 * self.mass * Total( hx2 * v * VM           * dv ) * 1e12
	hx3 = 0.5 * self.mass * Total( hx3 * v * Sqrt(VM)     * dv ) * 1e12
	hx4 =       self.mass * Total( hx4 * v * VM^(3.0/2.0) * dv ) * 1e12
	hx5 =       self.mass * Total( hx5 * v * VM^(3.0/2.0) * dv ) * 1e12
	hx6 =       self.mass * Total( hx6 * v * VM^(3.0/2.0) * dv ) * 1e12
	hx7 =       self.mass * Total( hx7 * v * VM           * dv ) * 1e12
	hx8 =       self.mass * Total( hx8 * v * VM           * dv ) * 1e12
	hx9 =       self.mass * Total( hx9 * v * VM           * dv ) * 1e12
	Hx  = hx0 - hx1 + hx2 - hx3 - hx4 - hx5 - hx6 + hx7 + hx8 + hx9
	
	;Hy
	hy0 = 0.5 * self.mass * Total( hy0 * v * VM^2.0       * dv ) * 1e12
	hy1 = 0.5 * self.mass * Total( hy1 * v * VM^(3.0/2.0) * dv ) * 1e12
	hy2 = 0.5 * self.mass * Total( hy2 * v * VM           * dv ) * 1e12
	hy3 = 0.5 * self.mass * Total( hy3 * v * Sqrt(VM)     * dv ) * 1e12
	hy4 =       self.mass * Total( hy4 * v * VM^(3.0/2.0) * dv ) * 1e12
	hy5 =       self.mass * Total( hy5 * v * VM^(3.0/2.0) * dv ) * 1e12
	hy6 =       self.mass * Total( hy6 * v * VM^(3.0/2.0) * dv ) * 1e12
	hy7 =       self.mass * Total( hy7 * v * VM           * dv ) * 1e12
	hy8 =       self.mass * Total( hy8 * v * VM           * dv ) * 1e12
	hy9 =       self.mass * Total( hy9 * v * VM           * dv ) * 1e12
	Hy  = hy0 - hy1 + hy2 - hy3 - hy4 - hy5 - hy6 + hy7 + hy8 + hy9
	
	;Hz
	hz0 = 0.5 * self.mass * Total( hz0 * v * VM^2.0       * dv ) * 1e12
	hz1 = 0.5 * self.mass * Total( hz1 * v * VM^(3.0/2.0) * dv ) * 1e12
	hz2 = 0.5 * self.mass * Total( hz2 * v * VM           * dv ) * 1e12
	hz3 = 0.5 * self.mass * Total( hz3 * v * Sqrt(VM)     * dv ) * 1e12
	hz4 =       self.mass * Total( hz4 * v * VM^(3.0/2.0) * dv ) * 1e12
	hz5 =       self.mass * Total( hz5 * v * VM^(3.0/2.0) * dv ) * 1e12
	hz6 =       self.mass * Total( hz6 * v * VM^(3.0/2.0) * dv ) * 1e12
	hz7 =       self.mass * Total( hz7 * v * VM           * dv ) * 1e12
	hz8 =       self.mass * Total( hz8 * v * VM           * dv ) * 1e12
	hz9 =       self.mass * Total( hz9 * v * VM           * dv ) * 1e12
	Hz  = hz0 - hz1 + hz2 - hz3 - hz4 - hz5 - hz6 + hz7 + hz8 + hz9
	
	;Return mW/m^2
	RETURN, -1e3 * [Hx, Hy, Hz]
END


;+
;   Compute the 3rd order moment, heat flux.
;
; :Returns:
;       Q:          out, required, type=float
;                   Heat flux of particles calculated from the distribution.
;-
FUNCTION MrDist3D::HeatFlux
	Compile_Opt idl2
	On_Error, 2

	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Convert to radians
	deg2rad = !dpi / 180D
	dims    = Size(*self.distfn, /DIMENSIONS)
	
	;Integrate phi
	dphi    = Median( (*self.phi)[1:*,*] - (*self.phi)[0:-2,*]) * deg2rad
	phi     = Rebin( *self.phi, dims ) * deg2rad
	fx_temp = Total( *self.distfn * Cos(phi), 1) * dPhi
	fy_temp = Total( *self.distfn * Sin(phi), 1) * dPhi
	fz_temp = Total( *self.distfn,            1) * dPhi

	;Integrate theta
	theta   = Rebin(Reform((*self.theta)[0,*]), dims[1], dims[2]) * deg2rad
	dTheta  = Median( theta[1:-1,*] - theta[0:-2,*] )
	fx_temp = Total( fx_temp * Sin(theta)^2, 1) * dTheta
	fy_temp = Total( fy_temp * Sin(theta)^2, 1) * dTheta
	fz_temp = Total( fz_temp * Sin(theta)*Cos(theta), 1) * dTheta
	
	;Integrate velocity
	q    = MrConstants('q')
	eV2J = MrConstants('eV2J')
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	;Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * *self.energy / self.mass)
	dE   = self -> DeltaE()
	dv   = eV2J * dE / (self.mass * v)
	
	;We need to be in the rest frame of the plasma
	;   - Convert BULKV from km/s to m/s
	bulkv = self -> Velocity()
	v     = v - Sqrt( Total(bulkv^2) ) * 1e3
	
	;Spacecraft potential correction
	IF N_Elements(*self.Vsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q * (*self.Vsc) / self.mass )
	
		;V is in m/s while f is in s^3/cm^6
		;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
		Qx = Total( fx_temp * v * (v^2 + sgn*vsc^2)^2 * dv ) * 1e12
		Qy = Total( fy_temp * v * (v^2 + sgn*vsc^2)^2 * dv ) * 1e12
		Qz = Total( fz_temp * v * (v^2 + sgn*vsc^2)^2 * dv ) * 1e12
	ENDIF ELSE BEGIN
		;V is in m/s while f is in s^3/cm^6
		;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
		Qx = Total( fx_temp * v^5 * dv ) * 1e12
		Qy = Total( fy_temp * v^5 * dv ) * 1e12
		Qz = Total( fz_temp * v^5 * dv ) * 1e12
	ENDELSE
	
	; Unit check:
	;   kg * s^3/m^6 * (m/s)^6   --> kg / s^3
	;   W / A = kg m^2 / s^3 / m^2 = kg / s^3
	
	;Form vector
	;   - V --> W/m^2 --> mW/m^2 * 1e3
	;   - Negative sign because the angles represent look directions,
	;     not incident directions
	Q = -1e3 * self.mass * [Qx, Qy, Qz] / 2.0

	RETURN, Q
END


;+
;   Set object properties.
;
; :Params:
;       DISTFN:             out, optional, type=NxMxL float
;                           Three-dimensional distribution function.
;       THETA:              out, optional, type=NxM float
;                           Polar location of each pixel in `DISTFN`.
;       PHI:                out, optional, type=NxM float
;                           Azimuthal location of each pixel in `DISTFN`.
;       ENERGY:             out, optional, type=Lx1 float
;                           Energy location of each pixel in `DISTFN`.
;-
pro MrDist3D::GetData, distfn, phi, theta, energy
	compile_opt idl2
	on_error, 2
	
	;Get the data
	if arg_present(distfn) then distfn = *self.distfn
	if arg_present(theta)  then theta  = *self.theta
	if arg_present(phi)    then phi    = *self.phi
	if arg_present(energy) then energy = *self.energy
end


;+
;   Get object properties.
;
; :Keywords:
;       DEGREES:            out, optional, type=boolean
;                           If set, the THETA and PHI properties have units of degrees.
;       ELEVATION:          out, optional, type=boolean
;                           If set, THETA is taken to be the elevation angle.
;       RADIANS:            out, optional, type=boolean
;                           If set, the THETA and PHI properties have units of radians.
;       UNITS:              out, optional, type=string
;                           Units of the distribution function.
;-
pro MrDist3D::GetProperty, $
ELEVATION=elevation, $
DEGREES=degrees, $
MASS=mass, $
RADIANS=radians, $
UNITS=units
	compile_opt idl2
	on_error, 2
	
	if arg_present(degrees)   then degrees   = self.degrees
	if arg_present(elevation) then elevation = self.elevation
	if arg_present(mass)      then mass      = self.mass
	if arg_present(radians)   then radians   = ~self.degrees
	if arg_present(units)     then units     = self.units
end


;+
;   Compute moments of the distribution: density, velocity, pressure, temperature,
;   and heat flux. This routine is faster than computing the moments individually,
;   and the calculations are inter-dependent.
;
; :Keywords:
;       DENSITY:        out, optional, type=float
;                       Plasma density.
;       ENTROPY:        out, optional, type=float
;                       Entropy.
;       VELOCITY:       out, optional, type=float
;                       Plasma bulk velocity.
;       PRESSURE:       out, optional, type=float
;                       Plasma pressure tensor.
;       TEMPERATURE:    out, optional, type=float
;                       Plasma temperature tensor.
;       HEATFLUX:       out, optional, type=float
;                       Heat flux.
;-
PRO MrDist3D::Moments_v2, $
DENSITY=density, $
ENTROPY=entropy, $
HEATFLUX=heatflux, $
PRESSURE=pressure, $
TEMPERATURE=temperature, $
VELOCITY=velocity, $
ENERGY_RANGE=energy_range, $
PHI_RANGE=phi_range, $
THETA_RANGE=theta_range
	Compile_Opt idl2
	On_Error, 2
	
	dims    = Size(*self.distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]
	
	;Must integrate over phase space
	IF N_Elements(phi_range)    EQ 0 THEN phi_range    = [0, nPhi    - 1]
	IF N_Elements(theta_range)  EQ 0 THEN theta_range  = [0, nTheta  - 1]
	IF N_Elements(energy_range) EQ 0 THEN energy_range = [0, nEnergy - 1]
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Convert to radians
	deg2rad = !dpi / 180D
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	eV2K    = 11600.0
	kB      = MrConstants('k_B')
	
	tf_nan = 0B
	f      = *self.distfn
	iZero  = Where(*self.distfn EQ 0, nZero)
	IF nZero GT 0 THEN BEGIN
		f[iZero] = !Values.F_NaN
		tf_nan = 1B
	ENDIF
	
	;Which products to return
	tf_n = Arg_Present(density)
	tf_s = Arg_Present(entropy)
	tf_q = Arg_Present(heatflux)
	tf_p = Arg_Present(pressure)
	tf_t = Arg_Present(temperature)
	tf_v = Arg_Present(velocity)
	tf_p = tf_p || tf_t
	tf_v = tf_v || tf_p || tf_q
	tf_n = tf_n || tf_v || tf_p || tf_t || tf_q
	
;-----------------------------------------------------
; Integration Bins and their Sizes \\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF MrIsA(phi_range, /INTEGER) THEN BEGIN
		iPhi = 0 > phi_range < (nPhi-1)
	ENDIF ELSE BEGIN
		iPhi = Value_Locate((*self.phi)[*,0], phi_range)
		IF iPhi[1] LT nPhi-1 THEN iPhi[1] += 1
	ENDELSE
	
	IF MrIsA(theta_range, /INTEGER) THEN BEGIN
		iTheta = 0 > theta_range < (nTheta-1)
	ENDIF ELSE BEGIN
		iTheta = Value_Locate((*self.theta)[0,*], theta_range)
		IF iTheta[1] LT nTheta-1 THEN iTheta[1] += 1
	ENDELSE
	
	IF MrIsA(energy_range, /INTEGER) THEN BEGIN
		iEnergy = 0 > energy_range < (nEnergy - 1)
	ENDIF ELSE BEGIN
		iEnergy = Value_Locate(*self.energy, energy_range)
	ENDELSE
	
	nPhi    = iPhi[1]    - iPhi[0]    + 1
	nTheta  = iTheta[1]  - iTheta[0]  + 1
	nEnergy = iEnergy[1] - iEnergy[0] + 1
	
	;Phi
	dPhi = Median( (*self.dphi)[iPhi[0]:iPhi[1]] ) * deg2rad
	phi  = (*self.phi)[iPhi[0]:iPhi[1], iTheta[0]:iTheta[1]] * deg2rad
	phi  = Rebin( phi, nPhi, nTheta, nEnergy )
	
	;Theta
	dTheta = Median( (*self.dtheta)[iTheta[0]:iTheta[1]] ) * deg2rad
	theta  = (*self.theta)[0, iTheta[0]:iTheta[1]] * deg2rad
	theta = Rebin( Reform(theta), nTheta, nEnergy )
	
	;Energy
	E  = (*self.energy)[iEnergy[0]:iEnergy[1]]
	dE = (*self.dEnergy)[iEnergy[0]:iEnergy[1]]
	v  = Sqrt(2.0 * eV2J * E / self.mass)
	dv = eV2J * dE / (self.mass * v)
	
	;Distribution function
	f = f[iPhi[0]:iPhi[1], iTheta[0]:iTheta[1], iEnergy[0]:iEnergy[1]]

	;PHI
;	dphi  = Median( (*self.phi)[1:*,*] - (*self.phi)[0:-2,*]) * deg2rad
;	phi   = Rebin( *self.phi, dims ) * deg2rad
	
	;THETA
;	theta  = Rebin(Reform((*self.theta)[0,*]), dims[1], dims[2]) * deg2rad
;	dTheta = Median( theta[1:-1,*] - theta[0:-2,*] )
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	;Convert energy from eV to J
;	v  = Sqrt(2.0 * eV2J * *self.energy / self.mass)
;	dE = self -> DeltaE()
;	dv = eV2J * dE / (self.mass * v)
	
	;SPACECRAFT POTENTIAL
	tf_Vsc = 0B
	IF N_Elements(*self.Vsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q * (*self.Vsc) / self.mass )
		VM   = v^2 * ( 1 + sgn * (vsc^2 / v^2) )
		tf_Vsc = 1B
	ENDIF ELSE BEGIN
		VM   = v^2
	ENDELSE
	
	;If the particle velocity is below the spacecraft potential
	;set its value to NaN to prevent Floating illegal operands.
	IF tf_Vsc THEN BEGIN
		iNaN = Where(v LT vsc, nNaN)
		IF nNaN GT 0 THEN BEGIN
			tf_nan  = 1B
			v[iNaN] = !values.f_nan
			VM[iNaN] = !values.f_nan
		ENDIF
	ENDIF

;-----------------------------------------------------
; Density \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF tf_n THEN BEGIN
		;Phi integral
		ftemp = Total( f, 1, NAN=tf_nan) * dPhi

		;Theta integral
		ftemp  = Total( ftemp * Sin(theta), 1, NAN=tf_nan) * dTheta

		;Integrate velocity
		;   - V is in m/s while f is in s^3/cm^6
		;   - v^2 * dv --> (m/s)^3  --> (cm/s)^3 * 1e6
		IF tf_Vsc $
			THEN density = Total( Temporary(ftemp) * v * Sqrt(v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e6 $
			ELSE density = Total( Temporary(ftemp) * v^2 * dv, NAN=tf_nan ) * 1e6
	ENDIF

;-----------------------------------------------------
; Entropy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF tf_s THEN BEGIN
		;Integrate Phi
		tf_s_nan = tf_nan
		iZero    = Where( f LE 0, nZero )
		ftemp    = f
		IF nZero GT 0 THEN BEGIN
			ftemp[iZero] = !Values.F_NaN
			tf_s_nan     = 1B
		ENDIF
		ftemp = Total( ftemp * alog(ftemp), 1, NAN=tf_s_nan) * dPhi

		;Integrate Theta
		ftemp  = Total( ftemp * Sin(theta), 1, NAN=tf_s_nan) * dTheta
	
		;Integrate velocity
		;   - Sign is charge dependent. E = q * V = 1/2 * m * v^2
		;   - V is in m/s while f is in s^3/cm^6
		;   - f --> s^3/cm^6  --> s^3/m^6 * 1e12
		IF tf_Vsc $
			THEN H = Total( ftemp * v * Sqrt(v^2 + sgn*vsc^2) * dv, NAN=tf_s_nan ) * 1e12 $
			ELSE H = Total( ftemp * v^2 * dv,                       NAN=tf_s_nan ) * 1e12

		;Entropy:
		;   - S = -kH
		entropy = -kB * Temporary(H)
	ENDIF

;-----------------------------------------------------
; Velocity \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF tf_v THEN BEGIN
		;Phi integral
		fx_temp = Total( f * Cos(phi), 1, NAN=tf_nan) * dPhi
		fy_temp = Total( f * Sin(phi), 1, NAN=tf_nan) * dPhi
		fz_temp = Total( f,            1, NAN=tf_nan) * dPhi
	
		;Theta integral
		fx_temp = Total( fx_temp * Sin(theta)^2,          1, NAN=tf_nan) * dTheta
		fy_temp = Total( fy_temp * Sin(theta)^2,          1, NAN=tf_nan) * dTheta
		fz_temp = Total( fz_temp * Sin(theta)*Cos(Theta), 1, NAN=tf_nan) * dTheta

		;Velocity integral
		;   - V is in m/s while f is in s^3/cm^6
		;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
		IF tf_Vsc THEN BEGIN
			Vx = Total( fx_temp * v * (v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e12
			Vy = Total( fy_temp * v * (v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e12
			Vz = Total( fz_temp * v * (v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e12
		ENDIF ELSE BEGIN
			Vx = Total( fx_temp * v^3 * dv, NAN=tf_nan ) * 1e12
			Vy = Total( fy_temp * v^3 * dv, NAN=tf_nan ) * 1e12
			Vz = Total( fz_temp * v^3 * dv, NAN=tf_nan ) * 1e12
		ENDELSE
	
		;Create vectory
		;   - V --> m/s --> km/s * 1e-3
		;   - N --> 1/cm^3 --> 1/m^3 * 1e6 ------> 1/N --> 1e-6
		;   - Negative sign because the angles represent look directions,
		;     not incident directions
		velocity = -1e-9 * [Vx, Vy, Vz] / density
		bulkv    = -1e-6 * [Vx, Vy, Vz] / density
		bulkvmag = Sqrt(Total(velocity^2, NAN=tf_nan))
	ENDIF

;-----------------------------------------------------
; Pressure \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	; Let
	;    VM = Vm^2 - sign(q) Vsc^2
	;
	; Where, here and in what follows,
	;    VM      =  The velocity corrected for spacecraft potential effects
	;    Vm      =  The measured velocity related to the energy bins of the instrument
	;    Vsc     =  The velocity gained by passing through the spacecraft potential sheath
	;    U       =  Plasma bulk velocity
	;    dVm     =  Size of a measured velocity space bin (energy bin width converted to velocity)
	;    dOmega  =  Sin(Theta) dTheta dPhi = Element of solid angle
	;    t       =  Theta, the polar angle
	;    p       =  Phi, the azimuthal angle
	;    m       =  Particle mass
	;    f       =  The distribution function
	;    q       =  Charge of the particle
	;
	; It follows that
	;
	;    Pxx = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 cos(p)^2      - 2 VM Ux sin(t) cos(p)                       + Sqrt(VM) Ux^2  ]
	;    Pxy = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 sin(p) cos(p) -   VM Uy sin(t) cos(p) - VM Ux sin(t) sin(p) + Sqrt(VM) Ux Uy ]
	;    Pxz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)   cos(t) cos(p) -   VM Uz sin(t) cos(p) - VM Ux cos(t)        + Sqrt(VM) Ux Uz ]
	;    Pyy = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 sin(p)^2      - 2 VM Uy sin(t) sin(p)                       + Sqrt(VM) Uy^2  ]
	;    Pyz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)   cos(t) sin(p) -   VM Uz sin(t) sin(p) - VM Uy cos(t)        + Sqrt(VM) Uy Uz ]
	;    Pzz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) cos(t)^2               - 2 VM Uz cos(t)                              + Sqrt(VM) Uz^2  ]
	;
	
	IF tf_p THEN BEGIN
		;
		;PHI
		;
		
		;Pxx
		fxx1_temp = Total( f * Cos(phi)^2, 1, NAN=tf_nan ) * dPhi
		fxx2_temp = Total( f * Cos(phi),   1, NAN=tf_nan ) * dPhi * 2.0 * bulkv[0]
		fxx3_temp = Total( f,              1, NAN=tf_nan ) * dPhi       * bulkv[0]^2
		
		;Pxy
		fxy1_temp = Total( f * Sin(phi) * Cos(phi), 1, NAN=tf_nan ) * dPhi
		fxy2_temp = Total( f * Cos(phi),            1, NAN=tf_nan ) * dPhi * bulkv[1]
		fxy3_temp = Total( f * Sin(phi),            1, NAN=tf_nan ) * dPhi * bulkv[0]
		fxy4_temp = Total( f,                       1, NAN=tf_nan ) * dPhi * bulkv[0] * bulkv[1]
		
		;Pxz
		fxz1_temp = Total( f * Cos(phi), 1, NAN=tf_nan ) * dPhi
		fxz2_temp = Total( f * Cos(phi), 1, NAN=tf_nan ) * dPhi * bulkv[2]
		fxz3_temp = Total( f,            1, NAN=tf_nan ) * dPhi * bulkv[0]
		fxz4_temp = Total( f,            1, NAN=tf_nan ) * dPhi * bulkv[0] * bulkv[2]
		
		;Pyy
		fyy1_temp = Total( f * Sin(phi)^2, 1, NAN=tf_nan ) * dPhi
		fyy2_temp = Total( f * Sin(phi),   1, NAN=tf_nan ) * dPhi * 2.0 * bulkv[1]
		fyy3_temp = Total( f,              1, NAN=tf_nan ) * dPhi       * bulkv[1]^2
		
		;Pyz
		fyz1_temp = Total( f * Sin(phi), 1, NAN=tf_nan ) * dPhi
		fyz2_temp = Total( f * Sin(phi), 1, NAN=tf_nan ) * dPhi * bulkv[2]
		fyz3_temp = Total( f,            1, NAN=tf_nan ) * dPhi * bulkv[1]
		fyz4_temp = Total( f,            1, NAN=tf_nan ) * dPhi * bulkv[1] * bulkv[2]
		
		;Pzz
		fzz1_temp = Total( f, 1, NAN=tf_nan ) * dPhi
		fzz2_temp = Total( f, 1, NAN=tf_nan ) * dPhi * 2 * bulkv[2]
		fzz3_temp = Total( f, 1, NAN=tf_nan ) * dPhi     * bulkv[2]^2
	
		;
		; THETA
		;
		
		;Pxx
		fxx1_temp = Total( fxx1_temp * Sin(theta)^3, 1, NAN=tf_nan ) * dTheta
		fxx2_temp = Total( fxx2_temp * Sin(theta)^2, 1, NAN=tf_nan ) * dTheta
		fxx3_temp = Total( fxx3_temp * Sin(theta),   1, NAN=tf_nan ) * dTheta
		
		;Pxy
		fxy1_temp = Total( fxy1_temp * Sin(theta)^3,             1, NAN=tf_nan ) * dTheta
		fxy2_temp = Total( fxy2_temp * Sin(theta)^2,             1, NAN=tf_nan ) * dTheta
		fxy3_temp = Total( fxy3_temp * Sin(theta) * Cos(theta),  1, NAN=tf_nan ) * dTheta
		fxy4_temp = Total( fxy4_temp * Sin(theta),               1, NAN=tf_nan ) * dTheta
		
		;Pxz
		fxz1_temp = Total( fxz1_temp * Sin(theta)^2 * Cos(theta), 1, NAN=tf_nan ) * dTheta
		fxz2_temp = Total( fxz2_temp * Sin(theta)^2,              1, NAN=tf_nan ) * dTheta
		fxz3_temp = Total( fxz3_temp * Sin(theta)   * Cos(theta), 1, NAN=tf_nan ) * dTheta
		fxz4_temp = Total( fxz4_temp * Sin(theta),                1, NAN=tf_nan ) * dTheta
		
		;Pyy
		fyy1_temp = Total( fyy1_temp * Sin(theta)^3, 1, NAN=tf_nan ) * dTheta
		fyy2_temp = Total( fyy2_temp * Sin(theta)^2, 1, NAN=tf_nan ) * dTheta
		fyy3_temp = Total( fyy3_temp * Sin(theta),   1, NAN=tf_nan ) * dTheta
		
		;Pyz
		fyz1_temp = Total( fyz1_temp * Sin(theta)^2 * Cos(theta), 1, NAN=tf_nan ) * dTheta
		fyz2_temp = Total( fyz2_temp * Sin(theta)^2,              1, NAN=tf_nan ) * dTheta
		fyz3_temp = Total( fyz3_temp * Sin(theta)   * Cos(theta), 1, NAN=tf_nan ) * dTheta
		fyz4_temp = Total( fyz4_temp * Sin(theta),                1, NAN=tf_nan ) * dTheta
		
		;Pzz
		fzz1_temp = Total( fzz1_temp * Sin(theta) * Cos(theta)^2, 1, NAN=tf_nan ) * dTheta
		fzz2_temp = Total( fzz2_temp * Sin(theta) * Cos(theta),   1, NAN=tf_nan ) * dTheta
		fzz3_temp = Total( fzz3_temp * Sin(theta),                1, NAN=tf_nan ) * dTheta
		
		
		;
		; ENERGY
		;
	
		;V is in m/s while f is in s^3/cm^6
		;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
		
		;Pxx
		fxx1_temp = self.mass * Total( fxx1_temp * v * VM^(3.0/2.0) * dv, NAN=tf_nan ) * 1e12
		fxx2_temp = self.mass * Total( fxx2_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fxx3_temp = self.mass * Total( fxx3_temp * v * Sqrt(VM)     * dv, NAN=tf_nan ) * 1e12
		Pxx       = Temporary(fxx1_temp) - Temporary(fxx2_temp) + Temporary(fxx3_temp)
		
		;Pxy
		fxy1_temp = self.mass * Total( fxy1_temp * v * VM^(3.0/2.0) * dv, NAN=tf_nan ) * 1e12
		fxy2_temp = self.mass * Total( fxy2_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fxy3_temp = self.mass * Total( fxy3_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fxy4_temp = self.mass * Total( fxy4_temp * v * Sqrt(VM)     * dv, NAN=tf_nan ) * 1e12
		Pxy       = Temporary(fxy1_temp) - Temporary(fxy2_temp) - Temporary(fxy3_temp) + Temporary(fxy4_temp)
		
		;Pxz
		fxz1_temp = self.mass * Total( fxz1_temp * v * VM^(3.0/2.0) * dv, NAN=tf_nan ) * 1e12
		fxz2_temp = self.mass * Total( fxz2_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fxz3_temp = self.mass * Total( fxz3_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fxz4_temp = self.mass * Total( fxz4_temp * v * Sqrt(VM)     * dv, NAN=tf_nan ) * 1e12
		Pxz       = Temporary(fxz1_temp) - Temporary(fxz2_temp) - Temporary(fxz3_temp) + Temporary(fxz4_temp)
		
		;Pyy
		fyy1_temp = self.mass * Total( fyy1_temp * v * VM^(3.0/2.0) * dv, NAN=tf_nan ) * 1e12
		fyy2_temp = self.mass * Total( fyy2_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fyy3_temp = self.mass * Total( fyy3_temp * v * Sqrt(VM)     * dv, NAN=tf_nan ) * 1e12
		Pyy       = Temporary(fyy1_temp) - Temporary(fyy2_temp) + Temporary(fyy3_temp)
		
		;Pyz
		fyz1_temp = self.mass * Total( fyz1_temp * v * VM^(3.0/2.0) * dv, NAN=tf_nan ) * 1e12
		fyz2_temp = self.mass * Total( fyz2_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fyz3_temp = self.mass * Total( fyz3_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fyz4_temp = self.mass * Total( fyz4_temp * v * Sqrt(VM)     * dv, NAN=tf_nan ) * 1e12
		Pyz       = Temporary(fyz1_temp) - Temporary(fyz2_temp) - Temporary(fyz3_temp) + Temporary(fyz4_temp)
		
		;Pzz
		fzz1_temp = self.mass * Total( fzz1_temp * v * VM^(3.0/2.0) * dv, NAN=tf_nan ) * 1e12
		fzz2_temp = self.mass * Total( fzz2_temp * v * VM           * dv, NAN=tf_nan ) * 1e12
		fzz3_temp = self.mass * Total( fzz3_temp * v * Sqrt(VM)     * dv, NAN=tf_nan ) * 1e12
		Pzz       = Temporary(fzz1_temp) - Temporary(fzz2_temp) + Temporary(fzz3_temp)
		
		; Unit check:
		;   kg * s^3/m^6 * (m/s)^5     --> kg / s^2 / m
		;   P = F / A = kg m / s^2 / m^2 = kg / s^2 / m
		
		;Form vector
		;   - Pa --> nPa * 1e9
		;   - Negative sign because the angles represent look directions,
		;     not incident directions
		Pressure = [ [Pxx, Pxy, Pxz], $
		             [Pxy, Pyy, Pyz], $
		             [Pxz, Pyz, Pzz] ] * 1e9
	ENDIF

;-----------------------------------------------------
; Temperature \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF tf_t THEN BEGIN
		;T = P/nk
		;   - 1/(1e15 * eV2K) converts units to eV
		temperature = pressure / (1e15 * eV2K * kB * density)
	ENDIF

;-----------------------------------------------------
; Heat Flux \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	IF tf_q THEN BEGIN
		vrest = v - bulkvmag
	
		;Phi integral
		fx_temp = Total( f * Cos(phi), 1, NAN=tf_nan) * dPhi
		fy_temp = Total( f * Sin(phi), 1, NAN=tf_nan) * dPhi
		fz_temp = Total( f,            1, NAN=tf_nan) * dPhi
		
		;Theta integral
		fx_temp = Total( fx_temp * Sin(theta)^2,          1, NAN=tf_nan) * dTheta
		fy_temp = Total( fy_temp * Sin(theta)^2,          1, NAN=tf_nan) * dTheta
		fz_temp = Total( fz_temp * Sin(theta)*Cos(theta), 1, NAN=tf_nan) * dTheta
		
		;Velocity integral
		IF tf_Vsc THEN BEGIN
			Qx = Total( fx_temp * vrest * (vrest^2 + sgn*vsc^2)^2 * dv, NAN=tf_nan ) * 1e12
			Qy = Total( fy_temp * vrest * (vrest^2 + sgn*vsc^2)^2 * dv, NAN=tf_nan ) * 1e12
			Qz = Total( fz_temp * vrest * (vrest^2 + sgn*vsc^2)^2 * dv, NAN=tf_nan ) * 1e12
		ENDIF ELSE BEGIN
			Qx = Total( fx_temp * vrest^5 * dv, NAN=tf_nan ) * 1e12
			Qy = Total( fy_temp * vrest^5 * dv, NAN=tf_nan ) * 1e12
			Qz = Total( fz_temp * vrest^5 * dv, NAN=tf_nan ) * 1e12
		ENDELSE
		
		; Unit check:
		;   kg * s^3/m^6 * (m/s)^6   --> kg / s^3
		;   W / A = kg m^2 / s^3 / m^2 = kg / s^3
	
		;Form vector
		;   - V --> W/m^2 --> mW/m^2 * 1e3
		;   - Negative sign because the angles represent look directions,
		;     not incident directions
		heatflux = -1e3 * self.mass * [Qx, Qy, Qz] / 2.0
	ENDIF
END


;+
;   Compute moments of the distribution: density, velocity, pressure, temperature,
;   and heat flux. This routine is faster than computing the moments individually,
;   and the calculations are inter-dependent.
;
; :Keywords:
;       DENSITY:        out, optional, type=float
;                       Plasma density.
;       ENTROPY:        out, optional, type=float
;                       Entropy.
;       VELOCITY:       out, optional, type=float
;                       Plasma bulk velocity.
;       PRESSURE:       out, optional, type=float
;                       Plasma pressure tensor.
;       TEMPERATURE:    out, optional, type=float
;                       Plasma temperature tensor.
;       HEATFLUX:       out, optional, type=float
;                       Heat flux.
;-
PRO MrDist3D::Moments, $
DENSITY=density, $
ENTROPY=entropy, $
HEATFLUX=heatflux, $
PRESSURE=pressure, $
TEMPERATURE=temperature, $
VELOCITY=velocity
	Compile_Opt idl2
	On_Error, 2
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Convert to radians
	deg2rad = !dpi / 180D
	q       = MrConstants('q')
	eV2J    = MrConstants('eV2J')
	eV2K    = 11600.0
	kB      = MrConstants('k_B')
	dims    = Size(*self.distfn, /DIMENSIONS)
	tf_nan  = ~Array_Equal( finite(*self.distfn), 1 )

;-----------------------------------------------------
; Integration Bins and their Sizes \\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;PHI
	dphi  = Median( (*self.phi)[1:*,*] - (*self.phi)[0:-2,*]) * deg2rad
	phi   = Rebin( *self.phi, dims ) * deg2rad
	
	;THETA
	theta  = Rebin(Reform((*self.theta)[0,*]), dims[1], dims[2]) * deg2rad
	dTheta = Median( theta[1:-1,*] - theta[0:-2,*] )
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	;Convert energy from eV to J
	v  = Sqrt(2.0 * eV2J * *self.energy / self.mass)
	dE = self -> DeltaE()
	dv = eV2J * dE / (self.mass * v)
	
	;SPACECRAFT POTENTIAL
	tf_Vsc = 0B
	IF N_Elements(*self.Vsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q * (*self.Vsc) / self.mass )
		VM   = v^2 * ( 1 + sgn * (vsc^2 / v^2) )
		tf_Vsc = 1B
	ENDIF ELSE BEGIN
		VM   = v^2
	ENDELSE

;-----------------------------------------------------
; Density \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Phi integral
	ftemp = Total( *self.distfn, 1, NAN=tf_nan) * dPhi

	;Theta integral
	ftemp  = Total( ftemp * Sin(theta), 1, NAN=tf_nan) * dTheta

	;Integrate velocity
	;   - V is in m/s while f is in s^3/cm^6
	;   - v^2 * dv --> (m/s)^3  --> (cm/s)^3 * 1e6
	IF tf_Vsc $
		THEN density = Total( Temporary(ftemp) * v * Sqrt(v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e6 $
		ELSE density = Total( Temporary(ftemp) * v^2 * dv, NAN=tf_nan ) * 1e6

;-----------------------------------------------------
; Entropy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Integrate Phi
	tf_s_nan = tf_nan
	iZero    = Where( *self.distfn LE 0, nZero )
	ftemp    = *self.distfn
	IF nZero GT 0 THEN BEGIN
		ftemp[iZero] = !Values.F_NaN
		tf_s_nan     = 1B
	ENDIF
	ftemp = Total( ftemp * alog(ftemp), 1, NAN=tf_s_nan) * dPhi

	;Integrate Theta
	ftemp  = Total( ftemp * Sin(theta), 1, NAN=tf_s_nan) * dTheta
	
	;Integrate velocity
	;   - Sign is charge dependent. E = q * V = 1/2 * m * v^2
	;   - V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6  --> s^3/m^6 * 1e12
	IF tf_Vsc $
		THEN H = Total( ftemp * v * Sqrt(v^2 + sgn*vsc^2) * dv, NAN=tf_s_nan ) * 1e12 $
		ELSE H = Total( ftemp * v^2 * dv,                       NAN=tf_s_nan ) * 1e12

	;Entropy:
	;   - S = -kH
	entropy = -kB * Temporary(H)

;-----------------------------------------------------
; Velocity \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Phi integral
	fx_temp = Total( *self.distfn * Cos(phi), 1, NAN=tf_nan) * dPhi
	fy_temp = Total( *self.distfn * Sin(phi), 1, NAN=tf_nan) * dPhi
	fz_temp = Total( *self.distfn,            1, NAN=tf_nan) * dPhi
	
	;Theta integral
	fx_temp = Total( fx_temp * Sin(theta)^2,          1, NAN=tf_nan) * dTheta
	fy_temp = Total( fy_temp * Sin(theta)^2,          1, NAN=tf_nan) * dTheta
	fz_temp = Total( fz_temp * Sin(theta)*Cos(Theta), 1, NAN=tf_nan) * dTheta

	;Velocity integral
	;   - V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
	IF tf_Vsc THEN BEGIN
		Vx = Total( fx_temp * v * (v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e12
		Vy = Total( fy_temp * v * (v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e12
		Vz = Total( fz_temp * v * (v^2 + sgn*vsc^2) * dv, NAN=tf_nan ) * 1e12
	ENDIF ELSE BEGIN
		Vx = Total( fx_temp * v^3 * dv, NAN=tf_nan ) * 1e12
		Vy = Total( fy_temp * v^3 * dv, NAN=tf_nan ) * 1e12
		Vz = Total( fz_temp * v^3 * dv, NAN=tf_nan ) * 1e12
	ENDELSE
	
	;Create vectory
	;   - V --> m/s --> km/s * 1e-3
	;   - N --> 1/cm^3 --> 1/m^3 * 1e6 ------> 1/N --> 1e-6
	;   - Negative sign because the angles represent look directions,
	;     not incident directions
	velocity = -1e-9 * [Vx, Vy, Vz] / density
	bulkvmag = Sqrt(Total(velocity^2, NAN=tf_nan))

;-----------------------------------------------------
; Pressure \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	
	;
	; Let
	;    VM = Vm^2 - sign(q) Vsc^2
	;
	; Where, here and in what follows,
	;    VM      =  The velocity corrected for spacecraft potential effects
	;    Vm      =  The measured velocity related to the energy bins of the instrument
	;    Vsc     =  The velocity gained by passing through the spacecraft potential sheath
	;    U       =  Plasma bulk velocity
	;    dVm     =  Size of a measured velocity space bin (energy bin width converted to velocity)
	;    dOmega  =  Sin(Theta) dTheta dPhi = Element of solid angle
	;    t       =  Theta, the polar angle
	;    p       =  Phi, the azimuthal angle
	;    m       =  Particle mass
	;    f       =  The distribution function
	;    q       =  Charge of the particle
	;
	; It follows that
	;
	;    Pxx = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 cos(p)^2      - 2 VM Ux sin(t) cos(p)                       + Sqrt(VM) Ux^2  ]
	;    Pxy = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 sin(p) cos(p) -   VM Uy sin(t) cos(p) - VM Ux sin(t) sin(p) + Sqrt(VM) Ux Uy ]
	;    Pxz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)   cos(t) cos(p) -   VM Uz sin(t) cos(p) - VM Ux cos(t)        + Sqrt(VM) Ux Uz ]
	;    Pyy = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 sin(p)^2      - 2 VM Uy sin(t) sin(p)                       + Sqrt(VM) Uy^2  ]
	;    Pyz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)   cos(t) sin(p) -   VM Uz sin(t) sin(p) - VM Uy cos(t)        + Sqrt(VM) Uy Uz ]
	;    Pzz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) cos(t)^2               - 2 VM Uz cos(t)                              + Sqrt(VM) Uz^2  ]
	;
	
	;
	;PHI
	;
	
	;Pxx
	fxx1_temp = Total( *self.distfn * Cos(phi)^2, 1 ) * dPhi
	fxx2_temp = Total( *self.distfn * Cos(phi),   1 ) * dPhi * 2.0 * velocity[0]
	fxx3_temp = Total( *self.distfn,              1 ) * dPhi       * velocity[0]^2
	
	;Pxy
	fxy1_temp = Total( *self.distfn * Sin(phi) * Cos(phi), 1 ) * dPhi
	fxy2_temp = Total( *self.distfn * Cos(phi),            1 ) * dPhi * velocity[1]
	fxy3_temp = Total( *self.distfn * Sin(phi),            1 ) * dPhi * velocity[0]
	fxy4_temp = Total( *self.distfn,                       1 ) * dPhi * velocity[0] * velocity[1]
	
	;Pxz
	fxz1_temp = Total( *self.distfn * Cos(phi), 1 ) * dPhi
	fxz2_temp = Total( *self.distfn * Cos(phi), 1 ) * dPhi * velocity[2]
	fxz3_temp = Total( *self.distfn,            1 ) * dPhi * velocity[0]
	fxz4_temp = Total( *self.distfn,            1 ) * dPhi * velocity[0] * velocity[2]
	
	;Pyy
	fyy1_temp = Total( *self.distfn * Sin(phi)^2, 1 ) * dPhi
	fyy2_temp = Total( *self.distfn * Sin(phi),   1 ) * dPhi * 2.0 * velocity[1]
	fyy3_temp = Total( *self.distfn,              1 ) * dPhi       * velocity[1]^2
	
	;Pyz
	fyz1_temp = Total( *self.distfn * Sin(phi), 1 ) * dPhi
	fyz2_temp = Total( *self.distfn * Sin(phi), 1 ) * dPhi * velocity[2]
	fyz3_temp = Total( *self.distfn,            1 ) * dPhi * velocity[1]
	fyz4_temp = Total( *self.distfn,            1 ) * dPhi * velocity[1] * velocity[2]
	
	;Pzz
	fzz1_temp = Total( *self.distfn, 1 ) * dPhi
	fzz2_temp = Total( *self.distfn, 1 ) * dPhi * 2 * velocity[2]
	fzz3_temp = Total( *self.distfn, 1 ) * dPhi     * velocity[2]^2

	;
	; THETA
	;
	
	;Pxx
	fxx1_temp = Total( fxx1_temp * Sin(theta)^3, 1 ) * dTheta
	fxx2_temp = Total( fxx2_temp * Sin(theta)^2, 1 ) * dTheta
	fxx3_temp = Total( fxx3_temp * Sin(theta),   1 ) * dTheta
	
	;Pxy
	fxy1_temp = Total( fxy1_temp * Sin(theta)^3,             1 ) * dTheta
	fxy2_temp = Total( fxy2_temp * Sin(theta)^2,             1 ) * dTheta
	fxy3_temp = Total( fxy3_temp * Sin(theta) * Cos(theta),  1 ) * dTheta
	fxy4_temp = Total( fxy4_temp * Sin(theta),               1 ) * dTheta
	
	;Pxz
	fxz1_temp = Total( fxz1_temp * Sin(theta)^2 * Cos(theta), 1 ) * dTheta
	fxz2_temp = Total( fxz2_temp * Sin(theta)^2,              1 ) * dTheta
	fxz3_temp = Total( fxz3_temp * Sin(theta)   * Cos(theta), 1 ) * dTheta
	fxz4_temp = Total( fxz4_temp * Sin(theta),                1 ) * dTheta
	
	;Pyy
	fyy1_temp = Total( fyy1_temp * Sin(theta)^3, 1 ) * dTheta
	fyy2_temp = Total( fyy2_temp * Sin(theta)^2, 1 ) * dTheta
	fyy3_temp = Total( fyy3_temp * Sin(theta),   1 ) * dTheta
	
	;Pyz
	fyz1_temp = Total( fyz1_temp * Sin(theta)^2 * Cos(theta), 1 ) * dTheta
	fyz2_temp = Total( fyz2_temp * Sin(theta)^2,              1 ) * dTheta
	fyz3_temp = Total( fyz3_temp * Sin(theta)   * Cos(theta), 1 ) * dTheta
	fyz4_temp = Total( fyz4_temp * Sin(theta),                1 ) * dTheta
	
	;Pzz
	fzz1_temp = Total( fzz1_temp * Sin(theta) * Cos(theta)^2, 1 ) * dTheta
	fzz2_temp = Total( fzz2_temp * Sin(theta) * Cos(theta),   1 ) * dTheta
	fzz3_temp = Total( fzz3_temp * Sin(theta),                1 ) * dTheta
	
	
	;
	; ENERGY
	;

	;V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
	
	;Pxx
	fxx1_temp = self.mass * Total( fxx1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fxx2_temp = self.mass * Total( fxx2_temp * v * VM           * dv ) * 1e12
	fxx3_temp = self.mass * Total( fxx3_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pxx       = Temporary(fxx1_temp) - Temporary(fxx2_temp) + Temporary(fxx3_temp)
	
	;Pxy
	fxy1_temp = self.mass * Total( fxy1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fxy2_temp = self.mass * Total( fxy2_temp * v * VM           * dv ) * 1e12
	fxy3_temp = self.mass * Total( fxy3_temp * v * VM           * dv ) * 1e12
	fxy4_temp = self.mass * Total( fxy4_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pxy       = Temporary(fxy1_temp) - Temporary(fxy2_temp) - Temporary(fxy3_temp) + Temporary(fxy4_temp)
	
	;Pxz
	fxz1_temp = self.mass * Total( fxz1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fxz2_temp = self.mass * Total( fxz2_temp * v * VM           * dv ) * 1e12
	fxz3_temp = self.mass * Total( fxz3_temp * v * VM           * dv ) * 1e12
	fxz4_temp = self.mass * Total( fxz4_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pxz       = Temporary(fxz1_temp) - Temporary(fxz2_temp) - Temporary(fxz3_temp) + Temporary(fxz4_temp)
	
	;Pyy
	fyy1_temp = self.mass * Total( fyy1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fyy2_temp = self.mass * Total( fyy2_temp * v * VM           * dv ) * 1e12
	fyy3_temp = self.mass * Total( fyy3_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pyy       = Temporary(fyy1_temp) - Temporary(fyy2_temp) + Temporary(fyy3_temp)
	
	;Pyz
	fyz1_temp = self.mass * Total( fyz1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fyz2_temp = self.mass * Total( fyz2_temp * v * VM           * dv ) * 1e12
	fyz3_temp = self.mass * Total( fyz3_temp * v * VM           * dv ) * 1e12
	fyz4_temp = self.mass * Total( fyz4_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pyz       = Temporary(fyz1_temp) - Temporary(fyz2_temp) - Temporary(fyz3_temp) + Temporary(fyz4_temp)
	
	;Pzz
	fzz1_temp = self.mass * Total( fzz1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fzz2_temp = self.mass * Total( fzz2_temp * v * VM           * dv ) * 1e12
	fzz3_temp = self.mass * Total( fzz3_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pzz       = Temporary(fzz1_temp) - Temporary(fzz2_temp) + Temporary(fzz3_temp)
	
	; Unit check:
	;   kg * s^3/m^6 * (m/s)^5     --> kg / s^2 / m
	;   P = F / A = kg m / s^2 / m^2 = kg / s^2 / m
	
	;Form vector
	;   - Pa --> nPa * 1e9
	;   - Negative sign because the angles represent look directions,
	;     not incident directions
	Pressure = [ [Pxx, Pxy, Pxz], $
	             [Pxy, Pyy, Pyz], $
	             [Pxz, Pyz, Pzz] ] * 1e9

;-----------------------------------------------------
; Temperature \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;T = P/nk
	;   - 1/(1e15 * eV2K) converts units to eV
	temperature = pressure / (1e15 * eV2K * kB * density)

;-----------------------------------------------------
; Heat Flux \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	vrest = v - bulkvmag

	;Phi integral
	fx_temp = Total( *self.distfn * Cos(phi), 1, NAN=tf_nan) * dPhi
	fy_temp = Total( *self.distfn * Sin(phi), 1, NAN=tf_nan) * dPhi
	fz_temp = Total( *self.distfn,            1, NAN=tf_nan) * dPhi
	
	;Theta integral
	fx_temp = Total( fx_temp * Sin(theta)^2,          1, NAN=tf_nan) * dTheta
	fy_temp = Total( fy_temp * Sin(theta)^2,          1, NAN=tf_nan) * dTheta
	fz_temp = Total( fz_temp * Sin(theta)*Cos(theta), 1, NAN=tf_nan) * dTheta
	
	;Velocity integral
	IF tf_Vsc THEN BEGIN
		Qx = Total( fx_temp * vrest * (vrest^2 + sgn*vsc^2)^2 * dv, NAN=tf_nan ) * 1e12
		Qy = Total( fy_temp * vrest * (vrest^2 + sgn*vsc^2)^2 * dv, NAN=tf_nan ) * 1e12
		Qz = Total( fz_temp * vrest * (vrest^2 + sgn*vsc^2)^2 * dv, NAN=tf_nan ) * 1e12
	ENDIF ELSE BEGIN
		Qx = Total( fx_temp * vrest^5 * dv, NAN=tf_nan ) * 1e12
		Qy = Total( fy_temp * vrest^5 * dv, NAN=tf_nan ) * 1e12
		Qz = Total( fz_temp * vrest^5 * dv, NAN=tf_nan ) * 1e12
	ENDELSE
	
	; Unit check:
	;   kg * s^3/m^6 * (m/s)^6   --> kg / s^3
	;   W / A = kg m^2 / s^3 / m^2 = kg / s^3

	;Form vector
	;   - V --> W/m^2 --> mW/m^2 * 1e3
	;   - Negative sign because the angles represent look directions,
	;     not incident directions
	heatflux = -1e3 * self.mass * [Qx, Qy, Qz] / 2.0
END


;+
;   Reduce a single 3D distribution function to 2D by averaging over polar angle.
;
; :Params:
;       PHI:            out, optional, type=fltarr
;                       Azimuthal angle (degrees) bins of the distribution.
;       ENERGY:         out, optional, type=fltarr
;                       Energy (eV) bins of the distribution.
;
; :Keywords:
;       NE_BINS:        in, optional, type=integer
;                       Number of energy bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of phi bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       THETA_RANGE:    in, optional, type=fltarr(2), default=[0.0\, 180.0]
;                       The range in polar angle (degrees) over which to average.
;
; :Returns:
;       DIST2D:         out, required, type=fltarr(L)
;                       The 2D distribution with size nEnergy.
;-
function MrDist3D::PhiE, phi, energy, dPhi, dE, $
NE_BINS=nE_bins, $
NPHI_BINS=nphi_bins, $
THETA_RANGE=theta_range
	compile_opt idl2
	on_error, 2

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Dimension sizes
	dims    = size(*self.distFn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

	;Amount to average over THETA
	phi_range = [-180.0, 180.0]
	if n_elements(theta_range) eq 0 then theta_range = [0.0, 180.0]
	if n_elements(nPhi_bins)   eq 0 then nPhi_bins   = nPhi
	if n_elements(nE_bins)     gt 0 then MrPrintF, 'logwarn', 'TODO: Allow nE_bins to change.'
	
	;Keep the number of energy bins the same.
	;   - Do not use NE_BINS to prevent cyclic messages there ---^
	nBins = nEnergy

;-----------------------------------------------------
; Coordinate Space \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Locate data within new bins
	;   - Keep the number of Theta bins the same
	;   - Average all PHI values within [-DELTA, DELTA] in each theta bin
	coords = transpose( [ [reform(*self.phi,   nPhi*nTheta)], $
	                      [reform(*self.theta, nPhi*nTheta)] ] )

;-----------------------------------------------------
; Weight \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	;Weight
	;   - Weight by size of old element
	;   - theta element: v*dtheta.
	;   - Weight is applied as total(w*psd)/total(w).
	;   - v and dtheta are constant, so the weight is 1.
	;
	
	;Weight function
	if self.elevation $
		then weight = rebin( cos( *self.theta * !dtor ),  nPhi, nTheta, nEnergy ) $
		else weight = rebin( sin( *self.theta * !dtor ),  nPhi, nTheta, nEnergy )

;-----------------------------------------------------
; Re-bin Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	cHist  = hist_nd(coords, $
	                 MIN             = [ phi_range[0], theta_range[0]], $
	                 MAX             = [ phi_range[1], theta_range[1]], $
	                 NBINS           = [    nPhi_bins,              1], $
	                 REVERSE_INDICES = ri)

;-----------------------------------------------------
; Reduce to 2D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Allocate memory to reduced 2D distribution
	dist2D = fltarr(nPhi_bins, nBins)

	;Loop over each PHI bin
	for j = 0, nEnergy - 1 do begin
		for k = 0, n_elements(chist) - 1 do begin
			;Skip empty bins
			if ri[k] eq ri[k+1] then continue
			
			;Source indices
			count = ri[k+1] - ri[k]
			inds  = ri[ri[k]:ri[k+1]-1]
			isrc  = array_indices([nPhi,nTheta], inds, /DIMENSIONS)
			
			;Weight
			w = weight[isrc[0,*], isrc[1,*], replicate(j, count)]
			
			;Re-bin
			dist2D[k,j] = total( (*self.distFn)[isrc[0,*], isrc[1,*], replicate(j, count)] * w ) / total(w)
		endfor
	endfor

;-----------------------------------------------------
; Output Bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Create bins
	dPhi  = (phi_range[1] - phi_range[0]) / nPhi_bins
	phi   = phi_range[0] + dPhi*findgen(nPhi_bins) + dPhi/2.0
	dPhi /= 2.0

	energy = *self.energy

	;Return the 2D distribution
	return, dist2D
end


;+
;   Reduce a single 3D distribution function to 1D by averaging over
;   energy and polar angle.
;
; :Params:
;       PHI:            out, optional, type=fltarr
;                       Azimuth angles (degrees) of the bin centers of the distribution.
;       DPHI:           out, optional, type=float
;                       Half-width of the `PHI` bins.
;
; :Keywords:
;       E_RANGE:        in, optional, type=fltarr(2), default=[min, max]
;                       The range in energy, in electron volts (eV) over which to average.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       THETA_RANGE:    in, optional, type=fltarr(2), default=[0.0\, 180.0]
;                       The range in polar angle (degrees) over which to average.
;
; :Returns:
;       DIST1D:         out, required, type=fltarr(L)
;                       The 1D distribution with size nEnergy.
;-
function MrDist3D::PhiSpec, phi, dPhi, $
E_RANGE=e_range, $
NPHI_BINS=nPhi_bins, $
THETA_RANGE=theta_range
	compile_opt idl2
	on_error, 2

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Dimension sizes
	dims    = size(*self.distFn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

	;PA range over which to average
	phi_range = [-180.0, 180.0]
	if n_elements(theta_range) eq 0 then theta_range = [0.0, 180.0]
	if n_elements(nPhi_bins)   eq 0 then nPhi_bins   = nPhi

	;Energy range over which to average
	dE = self -> DeltaE()
	if n_elements(e_range) eq 0 then begin
		e_range = [ (*self.energy)[0]  - (*self.energy)[0]  * dE[0]  / exp(1), $
		            (*self.energy)[-1] + (*self.energy)[-1] * dE[-1] / exp(1) ]
	endif


;-----------------------------------------------------
; Coordinate Space \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Combine data into 3-vectors of [phi, theta, energy]
	coords = transpose( [ [reform( rebin( *self.phi,   nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ], $
	                      [reform( rebin( *self.theta, nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ], $
	                      [reform( rebin( reform( alog10( *self.energy ), 1, 1, nEnergy ), nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ] ] )

;-----------------------------------------------------
; Weight \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;
	;   - w = v^2 * sin(theta) * dv * dTheta * dPhi
	;   - Applied as total(w*psd)/total(w)
	;   - Units of w cancel, so does not need same units as psd.
	;   - dTheta and dPhi are constants, so will factor out
	;   - Elements from two different old theta bins can fall into the same
	;     new theta bin, so the sin(theta) term does not cancel.
	;   => w = v^2 * sin(theta) * dv 
	;
	
	;Velocity element
	vsqr     = 2.0 * *self.energy / self.mass
	dv       = sqrt( 1.0 / (2.0 * self.mass * *self.energy ) ) * dE
	
	
	if self.elevation $
		then wTheta = cos(*self.theta * !dtor)$
		else wTheta = sin(*self.theta * !dtor)
	
	;Weight function
	weight = rebin( reform(vsqr, 1, 1, nEnergy), nPhi, nTheta, nEnergy ) * $
	         rebin( wTheta,                      nPhi, nTheta, nEnergy ) * $
	         rebin( reform(dv, 1, 1, nEnergy),   nPhi, nTheta, nEnergy )

;-----------------------------------------------------
; Re-bin Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Locate data within new bins
	;   - Keep the number of Theta bins the same
	;   - Average all PHI values within [-DELTA, DELTA] in each theta bin
	cHist  = hist_nd(coords, $
	                 MIN             = [ phi_range[0], theta_range[0], e_range[0]], $
	                 MAX             = [ phi_range[1], theta_range[1], e_range[1]], $
	                 NBINS           = [    nPhi_bins,              1,          1], $
	                 REVERSE_INDICES = ri)

;-----------------------------------------------------
; Reduce to 1D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Allocate memory to reduced 2D distribution
	dist1D = fltarr(nPhi)

	;Loop over bins
	for k = 0, n_elements(chist) - 1 do begin
		;Skip empty bins
		if ri[k] eq ri[k+1] then continue

		;Source indices
		inds = ri[ri[k]:ri[k+1]-1]
;		isrc = array_indices([nPhi,nTheta,nEnergy], inds, /DIMENSIONS)
		
		;Weight
;		w = weight[isrc[0,*], isrc[1,*], isrc[2,*]]

		;Re-bin
;		dist1D[k] = total( reform(distFn[isrc[0,*], isrc[1,*], isrc[2,*]]) * w ) / total(w)
		dist1D[k] = total( (*self.distFn)[inds] * weight[inds] ) / total(weight[inds])
	endfor

;-----------------------------------------------------
; Output Bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Create bins
	dPhi  = (phi_range[1] - phi_range[0]) / nPhi_bins
	phi   = phi_range[0] + dPhi*findgen(nPhi_bins) + dPhi/2.0
	dPhi /= 2.0

	;Return the 1D distribution
	return, dist1D
end


;+
;   Reduce a single 3D distribution function to 1D by averaging over
;   energy and polar angle.
;
; :Params:
;       PHI:            out, optional, type=fltarr
;                       Azimuth angles (degrees) of the bin centers of the distribution.
;       DPHI:           out, optional, type=float
;                       Half-width of the `PHI` bins.
;
; :Keywords:
;       E_RANGE:        in, optional, type=fltarr(2), default=[min, max]
;                       The range in energy, in electron volts (eV) over which to average.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       THETA_RANGE:    in, optional, type=fltarr(2), default=[0.0\, 180.0]
;                       The range in polar angle (degrees) over which to average.
;       WEIGHT:         in, optional, type=NxMxN fltarr, default=v^2 * sin(theta) * dv * dTheta * dPhi
;                       Weights for averaging bins on the velocity-space sphere. This is
;                           necessary only if the look directions have been rotated into
;                           a new coordinate system. They represent the velocity-space
;                           volume of each element in the distribution prior to rotation.
;
; :Returns:
;       DIST1D:         out, required, type=fltarr(L)
;                       The 1D distribution with size nEnergy.
;-
FUNCTION MrDist3D::PhiSpec_v2, phi, dPhi, $
E_RANGE=e_range, $
NPHI_BINS=nPhi_bins, $
THETA_RANGE=theta_range, $
WEIGHT=weight
	Compile_Opt idl2
	On_Error, 2

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Dimension sizes
	dims    = Size(*self.distFn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]
	IF N_Elements(units) EQ 0 THEN units = self.units

	;PA range over which to average
	phi_range = [-180.0, 180.0]
	IF N_Elements(nPhi_bins) EQ 0 THEN nPhi_bins = nPhi

	;Default range
	IF N_Elements(e_range) EQ 0 THEN BEGIN
		erange = [ (*self.energy)[0]  - (*self.dEnergy)[0]/2.0, $
		           (*self.energy)[-1] + (*self.dEnergy)[-1]/2.0 ]
	
	;Index range given
	ENDIF ELSE IF MrIsA(e_range, /INTEGER) THEN BEGIN
		erange = [ (*self.energy)[e_range[0]] - (*self.dEnergy)[e_range[0]]/2.0, $
		           (*self.energy)[e_range[1]] + (*self.dEnergy)[e_range[1]]/2.0 ]
	
	;Energy range
	ENDIF ELSE BEGIN
		erange = e_range
	ENDELSE
		
	
	;Default range
	IF N_Elements(theta_range) EQ 0 THEN BEGIN
		trange = [ (*self.theta)[0,0]  - (*self.dTheta)[0]/2.0, $
		           (*self.theta)[0,-1] + (*self.dTheta)[-1]/2.0 ]
	
	;Index range given
	ENDIF ELSE IF MrIsA(theta_range, /INTEGER) THEN BEGIN
		trange = [ (*self.theta)[0,theta_range[0]] - (*self.dTheta)[theta_range[0]]/2.0, $
		           (*self.theta)[0,theta_range[1]] + (*self.dTheta)[theta_range[1]]/2.0 ]
	
	;Angular range given
	ENDIF ELSE BEGIN
		trange = theta_range
	ENDELSE

;-----------------------------------------------------
; Coordinate Space \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	lnE    = ALog( *self.energy )
	erange = ALog(erange)
	
	;Combine data into 3-vectors of [phi, theta, energy]
	coords = Transpose( [ [Reform( Rebin( *self.phi,   nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ], $
	                      [Reform( Rebin( *self.theta, nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ], $
	                      [Reform( Rebin( Reform( lnE, 1, 1, nEnergy ), nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ] ] )

;-----------------------------------------------------
; Weight \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	;   - w = v^2 * sin(theta) * dv * dTheta * dPhi
	;   - Applied as total(w*psd)/total(w)
	;   - Units of w cancel, so does not need same units as psd.
	;   - dTheta and dPhi are constants, so will factor out
	;   - Elements from two different old theta bins can fall into the same
	;     new theta bin, so the sin(theta) term does not cancel.
	;   => w = v^2 * sin(theta) * dv 
	;
	
	IF N_Elements(weight) EQ 0 THEN BEGIN
		;Velocity element
		vsqr     = 2.0 * *self.energy / self.mass
		dv       = Sqrt( 1.0 / (2.0 * self.mass * *self.energy ) ) * *self.dEnergy
		
		
		IF self.elevation $
			THEN wTheta = Cos(*self.theta * !dtor) $
			ELSE wTheta = Sin(*self.theta * !dtor)
		
		;Weight function
		weight = Rebin( Reform(vsqr, 1, 1, nEnergy), nPhi, nTheta, nEnergy ) * $
		         Rebin( Temporary(wTheta),           nPhi, nTheta, nEnergy ) * $
		         Rebin( Reform(dv, 1, 1, nEnergy),   nPhi, nTheta, nEnergy )
	ENDIF

;-----------------------------------------------------
; Re-bin Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Locate data within new bins
	;   - Keep the number of Theta bins the same
	;   - Average all PHI values within [-DELTA, DELTA] in each theta bin
	cHist  = hist_nd(coords, $
	                 MIN             = [ phi_range[0], trange[0], erange[0]], $
	                 MAX             = [ phi_range[1], trange[1], erange[1]], $
	                 NBINS           = [    nPhi_bins,         1,         1], $
	                 REVERSE_INDICES = ri)

;-----------------------------------------------------
; Reduce to 1D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Allocate memory to reduced 2D distribution
	dist1D = FltArr(nPhi)

	;Loop over bins
	FOR k = 0, N_Elements(chist) - 1 DO BEGIN
		;Skip empty bins
		IF ri[k] eq ri[k+1] THEN CONTINUE

		;Source indices
		inds = ri[ri[k]:ri[k+1]-1]
;		isrc = array_indices([nPhi,nTheta,nEnergy], inds, /DIMENSIONS)
		
		;Weight
;		w = weight[isrc[0,*], isrc[1,*], isrc[2,*]]

		;Re-bin
;		dist1D[k] = total( reform(distFn[isrc[0,*], isrc[1,*], isrc[2,*]]) * w ) / total(w)
		dist1D[k] = Total( (*self.distFn)[inds] * weight[inds] ) / Total(weight[inds])
	ENDFOR

;-----------------------------------------------------
; Output Bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Create bins
	dPhi  = (phi_range[1] - phi_range[0]) / nPhi_bins
	phi   = phi_range[0] + dPhi*FIndGen(nPhi_bins) + dPhi/2.0
	dPhi /= 2.0

	;Return the 1D distribution
	RETURN, dist1D
END


;+
;   Compute the 0th order moment, denstiy.
;
;       vx = \sum_{i} \sum_{j} \sum_{k} f v_{i}^{2} \sin{\theta_{k}} \Delta v_{i} \Delta \phi_{j} \Delta \theta_{k}
;
; :Returns:
;       P:          out, required, type=float
;                   Pressure tensor calculated from the distribution.
;-
FUNCTION MrDist3D::Pressure
	Compile_Opt idl2
	On_Error, 2
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'
	
	;Convert to radians
	deg2rad = !dpi / 180D
	dims    = Size(*self.distfn, /DIMENSIONS)
	
	;Constants
	q     = MrConstants('q')
	eV2J  = MrConstants('eV2J')
	signQ = Round(self.mass / MrConstants('m_H')) EQ 0 ? -1 : 1

;-----------------------------------------------------
; Integration Bins and their Sizes \\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;PHI
	dphi = Median( (*self.phi)[1:*,*] - (*self.phi)[0:-2,*]) * deg2rad
	phi  = Rebin( *self.phi, dims ) * deg2rad
	
	;THETA
	theta    = Rebin(Reform((*self.theta)[0,*]), dims[1], dims[2]) * deg2rad
	dTheta   = Median( theta[1:-1,*] - theta[0:-2,*] )
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	;Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * *self.energy / self.mass)
	dE   = self -> DeltaE()
	dv   = eV2J * dE / (self.mass * v)
	
	;BULK VELOCITY
	;   - We need to be in the rest frame of the plasma
	;   - Convert BULKV from km/s to m/s
	bulkv = 1e3 * self -> Velocity()
	
	;SPACECRAFT POTENTIAL
	IF N_Elements(*self.Vsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q * (*self.Vsc) / self.mass )
		VM    = v^2 * (1 + signQ * (vsc^2 / v^2))
	ENDIF ELSE BEGIN
		VM = v^2
	ENDELSE

;-----------------------------------------------------
; Integrate \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	
	;
	; Let
	;    VM = Vm^2 - sign(q) Vsc^2
	;
	; Where, here and in what follows,
	;    VM      =  The velocity corrected for spacecraft potential effects
	;    Vm      =  The measured velocity related to the energy bins of the instrument
	;    Vsc     =  The velocity gained by passing through the spacecraft potential sheath
	;    U       =  Plasma bulk velocity
	;    dVm     =  Size of a measured velocity space bin (energy bin width converted to velocity)
	;    dOmega  =  Sin(Theta) dTheta dPhi = Element of solid angle
	;    t       =  Theta, the polar angle
	;    p       =  Phi, the azimuthal angle
	;    m       =  Particle mass
	;    f       =  The distribution function
	;    q       =  Charge of the particle
	;
	; It follows that
	;
	;    Pxx = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 cos(p)^2      - 2 VM Ux sin(t) cos(p)                       + Sqrt(VM) Ux^2  ]
	;    Pxy = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 sin(p) cos(p) -   VM Uy sin(t) cos(p) - VM Ux sin(t) sin(p) + Sqrt(VM) Ux Uy ]
	;    Pxz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)   cos(t) cos(p) -   VM Uz sin(t) cos(p) - VM Ux cos(t)        + Sqrt(VM) Ux Uz ]
	;    Pyy = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)^2 sin(p)^2      - 2 VM Uy sin(t) sin(p)                       + Sqrt(VM) Uy^2  ]
	;    Pyz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) sin(t)   cos(t) sin(p) -   VM Uz sin(t) sin(p) - VM Uy cos(t)        + Sqrt(VM) Uy Uz ]
	;    Pzz = m * \Integral{ f dVm dOmega Vm [ VM^(3/2) cos(t)^2               - 2 VM Uz cos(t)                              + Sqrt(VM) Uz^2  ]
	;
	
	;
	;PHI
	;
	
	;Pxx
	fxx1_temp = Total( *self.distfn * Cos(phi)^2, 1 ) * dPhi
	fxx2_temp = Total( *self.distfn * Cos(phi),   1 ) * dPhi * 2.0 * bulkv[0]
	fxx3_temp = Total( *self.distfn,              1 ) * dPhi       * bulkv[0]^2
	
	;Pxy
	fxy1_temp = Total( *self.distfn * Sin(phi) * Cos(phi), 1 ) * dPhi
	fxy2_temp = Total( *self.distfn * Cos(phi),            1 ) * dPhi * bulkv[1]
	fxy3_temp = Total( *self.distfn * Sin(phi),            1 ) * dPhi * bulkv[0]
	fxy4_temp = Total( *self.distfn,                       1 ) * dPhi * bulkv[0] * bulkv[1]
	
	;Pxz
	fxz1_temp = Total( *self.distfn * Cos(phi), 1 ) * dPhi
	fxz2_temp = Total( *self.distfn * Cos(phi), 1 ) * dPhi * bulkv[2]
	fxz3_temp = Total( *self.distfn,            1 ) * dPhi * bulkv[0]
	fxz4_temp = Total( *self.distfn,            1 ) * dPhi * bulkv[0] * bulkv[2]
	
	;Pyy
	fyy1_temp = Total( *self.distfn * Sin(phi)^2, 1 ) * dPhi
	fyy2_temp = Total( *self.distfn * Sin(phi),   1 ) * dPhi * 2.0 * bulkv[1]
	fyy3_temp = Total( *self.distfn,              1 ) * dPhi       * bulkv[1]^2
	
	;Pyz
	fyz1_temp = Total( *self.distfn * Sin(phi), 1 ) * dPhi
	fyz2_temp = Total( *self.distfn * Sin(phi), 1 ) * dPhi * bulkv[2]
	fyz3_temp = Total( *self.distfn,            1 ) * dPhi * bulkv[1]
	fyz4_temp = Total( *self.distfn,            1 ) * dPhi * bulkv[1] * bulkv[2]
	
	;Pzz
	fzz1_temp = Total( *self.distfn, 1 ) * dPhi
	fzz2_temp = Total( *self.distfn, 1 ) * dPhi * 2 * bulkv[2]
	fzz3_temp = Total( *self.distfn, 1 ) * dPhi     * bulkv[2]^2

	;
	; THETA
	;
	
	;Pxx
	fxx1_temp = Total( fxx1_temp * Sin(theta)^3, 1 ) * dTheta
	fxx2_temp = Total( fxx2_temp * Sin(theta)^2, 1 ) * dTheta
	fxx3_temp = Total( fxx3_temp * Sin(theta),   1 ) * dTheta
	
	;Pxy
	fxy1_temp = Total( fxy1_temp * Sin(theta)^3,             1 ) * dTheta
	fxy2_temp = Total( fxy2_temp * Sin(theta)^2,             1 ) * dTheta
	fxy3_temp = Total( fxy3_temp * Sin(theta) * Cos(theta),  1 ) * dTheta
	fxy4_temp = Total( fxy4_temp * Sin(theta),               1 ) * dTheta
	
	;Pxz
	fxz1_temp = Total( fxz1_temp * Sin(theta)^2 * Cos(theta), 1 ) * dTheta
	fxz2_temp = Total( fxz2_temp * Sin(theta)^2,              1 ) * dTheta
	fxz3_temp = Total( fxz3_temp * Sin(theta)   * Cos(theta), 1 ) * dTheta
	fxz4_temp = Total( fxz4_temp * Sin(theta),                1 ) * dTheta
	
	;Pyy
	fyy1_temp = Total( fyy1_temp * Sin(theta)^3, 1 ) * dTheta
	fyy2_temp = Total( fyy2_temp * Sin(theta)^2, 1 ) * dTheta
	fyy3_temp = Total( fyy3_temp * Sin(theta),   1 ) * dTheta
	
	;Pyz
	fyz1_temp = Total( fyz1_temp * Sin(theta)^2 * Cos(theta), 1 ) * dTheta
	fyz2_temp = Total( fyz2_temp * Sin(theta)^2,              1 ) * dTheta
	fyz3_temp = Total( fyz3_temp * Sin(theta)   * Cos(theta), 1 ) * dTheta
	fyz4_temp = Total( fyz4_temp * Sin(theta),                1 ) * dTheta
	
	;Pzz
	fzz1_temp = Total( fzz1_temp * Sin(theta) * Cos(theta)^2, 1 ) * dTheta
	fzz2_temp = Total( fzz2_temp * Sin(theta) * Cos(theta),   1 ) * dTheta
	fzz3_temp = Total( fzz3_temp * Sin(theta),                1 ) * dTheta
	
	
	;
	; VELOCITY
	;

	;V is in m/s while f is in s^3/cm^6
	;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
	
	;Pxx
	fxx1_temp = self.mass * Total( fxx1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fxx2_temp = self.mass * Total( fxx2_temp * v * VM           * dv ) * 1e12
	fxx3_temp = self.mass * Total( fxx3_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pxx       = Temporary(fxx1_temp) - Temporary(fxx2_temp) + Temporary(fxx3_temp)
	
	;Pxy
	fxy1_temp = self.mass * Total( fxy1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fxy2_temp = self.mass * Total( fxy2_temp * v * VM           * dv ) * 1e12
	fxy3_temp = self.mass * Total( fxy3_temp * v * VM           * dv ) * 1e12
	fxy4_temp = self.mass * Total( fxy4_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pxy       = Temporary(fxy1_temp) - Temporary(fxy2_temp) - Temporary(fxy3_temp) + Temporary(fxy4_temp)
	
	;Pxz
	fxz1_temp = self.mass * Total( fxz1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fxz2_temp = self.mass * Total( fxz2_temp * v * VM           * dv ) * 1e12
	fxz3_temp = self.mass * Total( fxz3_temp * v * VM           * dv ) * 1e12
	fxz4_temp = self.mass * Total( fxz4_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pxz       = Temporary(fxz1_temp) - Temporary(fxz2_temp) - Temporary(fxz3_temp) + Temporary(fxz4_temp)
	
	;Pyy
	fyy1_temp = self.mass * Total( fyy1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fyy2_temp = self.mass * Total( fyy2_temp * v * VM           * dv ) * 1e12
	fyy3_temp = self.mass * Total( fyy3_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pyy       = Temporary(fyy1_temp) - Temporary(fyy2_temp) + Temporary(fyy3_temp)
	
	;Pyz
	fyz1_temp = self.mass * Total( fyz1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fyz2_temp = self.mass * Total( fyz2_temp * v * VM           * dv ) * 1e12
	fyz3_temp = self.mass * Total( fyz3_temp * v * VM           * dv ) * 1e12
	fyz4_temp = self.mass * Total( fyz4_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pyz       = Temporary(fyz1_temp) - Temporary(fyz2_temp) - Temporary(fyz3_temp) + Temporary(fyz4_temp)
	
	;Pzz
	fzz1_temp = self.mass * Total( fzz1_temp * v * VM^(3.0/2.0) * dv ) * 1e12
	fzz2_temp = self.mass * Total( fzz2_temp * v * VM           * dv ) * 1e12
	fzz3_temp = self.mass * Total( fzz3_temp * v * Sqrt(VM)     * dv ) * 1e12
	Pzz       = Temporary(fzz1_temp) - Temporary(fzz2_temp) + Temporary(fzz3_temp)
	
	; Unit check:
	;   kg * s^3/m^6 * (m/s)^5     --> kg / s^2 / m
	;   P = F / A = kg m / s^2 / m^2 = kg / s^2 / m
	
	;Form vector
	;   - Pa --> nPa * 1e9
	;   - Negative sign because the angles represent look directions,
	;     not incident directions
	P = [ [Pxx, Pxy, Pxz], $
	      [Pxy, Pyy, Pyz], $
	      [Pxz, Pyz, Pzz] ] * 1e9

	RETURN, P
END


;+
;   Rebin the distribution.
;
; :Params:
;       DV:             in, required, type=NxMxL fltarr, default=self -> VolumeElement()
;                       Volume of each velocity space element of the distribution function
;                           before re-binning.
;
; :Keywords:
;       NPHI_BINS:      in, optional, type=integer, default=same as implicit f
;                       Number of azimuthal bins in the output distribution.
;       NTHETA_BINS:    in, optional, type=integer, default=same as implicit f
;                       Number of polar bins in the output distribution.
;       PHI_RANGE:      in, optional, type=fltarr(2), default=[-180.0, 180.0]
;                       Azimuthal range of the output distribution function.
;       THETA_RANGE:    in, optional, type=fltarr(2), default=[0.0, 180.0]
;                       Polar range of the output distribution function.
;
; :Returns:
;       DIST3D:         out, required, type=NxMxL fltarr
;                       The re-binned distribution function.
;-
FUNCTION MrDist3D::RebinAngles, dV, phi, dphi, theta, dtheta, $
NPHI_BINS=nPhi_bins, $
NTHETA_BINS=nTheta_bins, $
PHI_RANGE=phi_range, $
THETA_RANGE=theta_range
	Compile_Opt idl2
	On_Error, 2
	
;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	dims    = Size(*self.distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]
	IF N_Elements(nPhi_Bins)   EQ 0 THEN nPhi_Bins   = nPhi
	IF N_Elements(nTheta_Bins) EQ 0 THEN nTheta_Bins = nTheta
	IF N_Elements(phi_range)   EQ 0 THEN phi_range   = [-180.0, 180.0]
	IF N_Elements(theta_range) EQ 0 THEN theta_range = [   0.0, 180.0]

;-----------------------------------------------------
; Re-bin Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Combine data into 3-vectors of [phi, theta, energy]
	coords = Transpose( [ [ Reform( Rebin( *self.phi,   [nPhi, nTheta] ), nPhi*nTheta ) ], $
	                      [ Reform( Rebin( *self.theta, [nPhi, nTheta] ), nPhi*nTheta ) ] ] )
	
	;Locate data within new bins
	;   - Keep the number of Theta bins the same
	;   - Average all PHI values within [-DELTA, DELTA] in each theta bin
	cHist  = hist_nd(coords, $
	                 MIN             = [ phi_range[0], theta_range[0] ], $
	                 MAX             = [ phi_range[1], theta_range[1] ], $
	                 NBINS           = [    nPhi_Bins,    nTheta_Bins ], $
	                 REVERSE_INDICES = ri)

;-----------------------------------------------------
; Reduce to 2D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Allocate memory to reduced 2D distribution
	dist3D = FltArr(nPhi_Bins, nTheta_Bins, nEnergy)

	;Loop over bins
	FOR j = 0, nEnergy - 1 DO BEGIN
		FOR k = 0, N_Elements(chist) - 1 DO BEGIN
			;Skip empty bins
			IF ri[k] EQ ri[k+1] THEN CONTINUE
			
			;Source indices
			inds = ri[ri[k]:ri[k+1]-1]
			isrc = Array_Indices([nPhi,nTheta], inds, /DIMENSIONS)
			
			;Destination indices
			idest = Array_Indices([nPhi,nTheta], k, /DIMENSIONS)
		
			;Weight
			w = dV[isrc[0,*], isrc[1,*], j]
			
			;Re-bin
			dist3D[idest[0], idest[1], j] = Total( Reform((*self.distFn)[isrc[0,*], isrc[1,*], j]) * w ) / Total(w)
		ENDFOR
	ENDFOR

;-----------------------------------------------------
; Output Bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Phi bins
	dPhi  = (phi_range[1] - phi_range[0]) / nPhi_bins
	phi   = phi_range[0] + dPhi*findgen(nPhi_bins) + dPhi/2.0
	dPhi /= 2.0
	
	;Theta bins
	dTheta  = (theta_range[1] - theta_range[0]) / nTheta_bins
	theta   = theta_range[0] + dTheta*findgen(nTheta_bins) + dTheta/2.0
	dTheta /= 2.0

	;Return the 2D distribution
	RETURN, dist3D
END


;+
;   Set object properties.
;
; :Params:
;       DISTFN:             in, required, type=NxMxL float
;                           Three-dimensional distribution function.
;       PHI:                in, required, type=NxM float
;                           Azimuthal location of each pixel in `DISTFN`.
;       THETA:              in, required, type=NxM float
;                           Polar location of each pixel in `DISTFN`.
;       ENERGY:             in, required, type=Lx1 float
;                           Energy location of each pixel in `DISTFN`.
;-
PRO MrDist3D::SetData2, distfn, phi, theta, energy, vsc, $
DEGREES=degrees, $
DENERGY_MINUS=dEnergy_minus, $
DENERGY_PLUS=dEnergy_plus, $
DPHI_MINUS=dPhi_minus, $
DPHI_PLUS=dPhi_plus, $
DTHETA_MINUS=dTheta_minus, $
DTHETA_PLUS=dTheta_plus, $
RADIANS=radians, $
UNITS=units
	Compile_Opt idl2
	On_Error, 2
	
	;Unis
	IF N_Elements(units) EQ 0 THEN units = 'E FLUX'
	
	nDims   = size(distfn, /N_DIMENSIONS)
	dims    = size(distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

;-----------------------------------------------------
; Distribution Function \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Distribution function must be 3D: [nPhi, nTheta, nEnergy]
	IF nDims NE 3 THEN Message, 'DISTFN must have 3 dimensions.'

;-----------------------------------------------------
; Spacecraft Potential \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	nVsc = n_elements(Vsc)
	if nVsc gt 0 && nVsc ne 1 then message, 'Vsc must be scalar.'

;-----------------------------------------------------
; Set Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Set data
	*self.distfn = distfn
	self -> SetEnergy, energy, dEnergy_minus, dEnergy_plus
	self -> SetPhi,    phi,    dPhi_minus,    dPhi_plus,   DEGREES=degrees
	self -> SetTheta,  theta,  dTheta_minus,  dTheta_plus, DEGREES=degrees
	IF nVsc GT 0 THEN *self.Vsc = Vsc[0]
	self.units   = units
END


;+
;   Set object properties.
;
; :Params:
;       DISTFN:             in, required, type=NxMxL float
;                           Three-dimensional distribution function.
;       PHI:                in, required, type=NxM float
;                           Azimuthal location of each pixel in `DISTFN`.
;       THETA:              in, required, type=NxM float
;                           Polar location of each pixel in `DISTFN`.
;       ENERGY:             in, required, type=Lx1 float
;                           Energy location of each pixel in `DISTFN`.
;-
pro MrDist3D::SetData, distfn, phi, theta, energy, vsc, $
DEGREES=degrees, $
RADIANS=radians, $
UNITS=units
	compile_opt idl2
	on_error, 2
	
	;Unis
	if n_elements(units) eq 0 then units = 'E FLUX'
	if n_elements(degrees) eq 0 && n_elements(radians) eq 0 then degrees = 1B
	tf_degrees = n_elements(radians) gt 0 ? ~keyword_set(radians) : keyword_set(degrees)
	
	nDims   = size(distfn, /N_DIMENSIONS)
	dims    = size(distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

;-----------------------------------------------------
; Distribution Function \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Distribution function must be 3D: [nPhi, nTheta, nEnergy]
	if nDims ne 3 then message, 'DISTFN must have 3 dimensions.'

;-----------------------------------------------------
; Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Phi must be 2D: [nPhi, nTheta]
	nDims = size(phi, /N_DIMENSIONS)
	dims  = size(phi, /DIMENSIONS)
	
	;1D
	if nDims eq 1 then begin
		;Reform from [nPhi] to [nPhi, nTheta]
		if dims[0] ne nPhi $
			then message, 'PHI must be the size of the first dimension of DISTFN.' $
			else phOut = rebin( phi, nPhi, nTheta )
		
	;2D
	endif else if nDims eq 2 then begin
		if (dims[0] ne nPhi) || (dims[1] ne nTheta) $
			then message, 'PHI must have dimensions of [nPhi, nTheta].' $
			else phOut = phi
	
	;3+D
	endif else begin
		message, 'PHI must have one or two dimensions.'
	endelse

;-----------------------------------------------------
; Theta \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Theta must be 2D: [nPhi, nTheta]
	nDims = size(theta, /N_DIMENSIONS)
	dims  = size(theta, /DIMENSIONS)
	
	;1D
	if nDims eq 1 then begin
		;Reform from [nTheta] to [nPhi, nTheta]
		if dims[0] ne nTheta $
			then message, 'THETA must be the size of the second dimension of DISTFN.' $
			else thOut = rebin( reform(theta, 1, nTheta), nPhi, nTheta )
		
	;2D
	endif else if nDims eq 2 then begin
		if (dims[0] ne nPhi) || (dims[1] ne nTheta) $
			then message, 'THETA must have dimensions of [nPhi, nTheta].' $
			else thOut = theta
	
	;3+D
	endif else begin
		message, 'THETA must have one or two dimensions.'
	endelse

;-----------------------------------------------------
; Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Energy must be a row or column vector
	nDims = size(energy, /N_DIMENSIONS)
	dims  = size(energy, /DIMENSIONS)
	if (nDims eq 2 && dims[0] ne 1) || nDims ge 3 then message, 'ENERGY must be a column or row vector.'
	if n_elements(energy) ne nEnergy $
		then message, 'ENERGY must be the size of third dimension of DISTFN.'

;-----------------------------------------------------
; Spacecraft Potential \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	nVsc = n_elements(Vsc)
	if nVsc gt 0 && nVsc ne 1 then message, 'Vsc must be scalar.'

;-----------------------------------------------------
; Set Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Set data
	*self.distfn = distfn
	*self.phi    = tf_degrees ? temporary(phOut) : temporary(phOut) * !radeg
	*self.theta  = tf_degrees ? temporary(thOut) : temporary(thOut) * !radeg
	*self.energy = energy
	if nVsc gt 0 then *self.Vsc = Vsc[0]
	self.degrees = 1B
	self.units   = units
end


;+
;   Set object properties.
;
; :Params:
;       DISTFN:             in, required, type=NxMxL float
;                           Three-dimensional distribution function.
;       PHI:                in, required, type=NxM float
;                           Azimuthal location of each pixel in `DISTFN`.
;       THETA:              in, required, type=NxM float
;                           Polar location of each pixel in `DISTFN`.
;       ENERGY:             in, required, type=Lx1 float
;                           Energy location of each pixel in `DISTFN`.
;-
PRO MrDist3D::SetEnergy, energy, dE_Minus, dE_Plus, $
DEGREES=degrees, $
RADIANS=radians
	Compile_Opt idl2
	On_Error, 2
	
	dims    = Size(*self.distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

;-----------------------------------------------------
; Energy \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Energy must be a row or column vector
	nDims = Size(energy, /N_DIMENSIONS)
	dims  = Size(energy, /DIMENSIONS)
	IF (nDims EQ 2 && dims[0] NE 1) || nDims GE 3 THEN Message, 'ENERGY must be a column or row vector.'
	IF N_Elements(energy) ne nEnergy $
		THEN Message, 'ENERGY must be the size of third dimension of DISTFN.'
	
	;Set the energy
	*self.energy = energy

;-----------------------------------------------------
; Deltas \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Assume values given are the bin centers
	IF N_Elements(dE_Plus) EQ 0 && N_Elements(dE_Minus) EQ 0 THEN BEGIN
		dE_Plus  = 0.5 * self -> DeltaE()
		dE_Minus = dE_Plus
	ENDIF
	
	;DPLUS
	;   - Get a vector of values
	CASE N_Elements(dE_Plus) OF
		0:       dPlus = Replicate(0.0, nEnergy)
		1:       dPlus = Replicate(dE_Plus, nEnergy)
		nEnergy: dPlus = dE_Plus
		ELSE: Message, 'DPLUS had incorrect number of elements.'
	ENDCASE
	
	;DPLUS
	;   - Get a vector of values
	CASE N_Elements(dE_Minus) OF
		0:       dMinus = Replicate(0.0, nEnergy)
		1:       dMinus = Replicate(dE_Minus, nEnergy)
		nEnergy: dMinus = dE_Minus
		ELSE: Message, 'DMINUS had incorrect number of elements.'
	ENDCASE
	
	;Move PHI to bottom of bin
	E      = energy + (dPlus - dMinus)/2.0
	delta  = dPlus + dMinus

;-----------------------------------------------------
; Set Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Set data
	*self.energy  = Temporary(E)
	*self.dEnergy = Temporary(delta)
END


;+
;   Set object properties.
;
; :Params:
;       DISTFN:             in, required, type=NxMxL float
;                           Three-dimensional distribution function.
;       PHI:                in, required, type=NxM float
;                           Azimuthal location of each pixel in `DISTFN`.
;       THETA:              in, required, type=NxM float
;                           Polar location of each pixel in `DISTFN`.
;       ENERGY:             in, required, type=Lx1 float
;                           Energy location of each pixel in `DISTFN`.
;-
PRO MrDist3D::SetPhi, phi, dPhi_Minus, dPhi_Plus, $
DEGREES=degrees, $
RADIANS=radians
	Compile_Opt idl2
	On_Error, 2
	
	;Unis
	IF N_Elements(degrees) eq 0 && N_Elements(radians) eq 0 then degrees = 1B
	tf_degrees = N_Elements(radians) GT 0 ? ~Keyword_Set(radians) : Keyword_Set(degrees)
	
	dims    = Size(*self.distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

;-----------------------------------------------------
; Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Phi must be 2D: [nPhi, nTheta]
	nDims = Size(phi, /N_DIMENSIONS)
	dims  = Size(phi, /DIMENSIONS)
	
	;1D
	IF nDims EQ 1 THEN BEGIN
		;Reform from [nPhi] to [nPhi, nTheta]
		IF dims[0] NE nPhi $
			THEN Message, 'PHI must be the size of the first dimension of DISTFN.' $
			ELSE phOut = Rebin( phi, nPhi, nTheta )
		
	;2D
	ENDIF ELSE IF nDims EQ 2 THEN BEGIN
		IF (dims[0] NE nPhi) || (dims[1] NE nTheta) $
			THEN Message, 'PHI must have dimensions of [nPhi, nTheta].' $
			ELSE phOut = phi
	
	;3+D
	ENDIF ELSE BEGIN
		Message, 'PHI must have one or two dimensions.'
	ENDELSE

;-----------------------------------------------------
; Deltas \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Assume values given are the bin centers
	IF N_Elements(dPhi_Plus) EQ 0 && N_Elements(dPhi_Minus) EQ 0 THEN BEGIN
		dPhi_Plus  = Median(phOut[1:*,0] - phOut[*,0]) / 2.0
		dPhi_Minus = dPhi_Plus
	ENDIF

	;DPLUS
	;   - Get a vector of values
	CASE N_Elements(dPhi_Plus) OF
		0:    dPlus = Replicate(0.0, nPhi)
		1:    dPlus = Replicate(dPhi_Plus, nPhi)
		nPhi: dPlus = dPhi_Plus
		ELSE: Message, 'DPLUS had incorrect number of elements.'
	ENDCASE
	
	;DPLUS
	;   - Get a vector of values
	CASE N_Elements(dPhi_Minus) OF
		0:    dMinus = Replicate(0.0, nPhi)
		1:    dMinus = Replicate(dPhi_Minus, nPhi)
		nPhi: dMinus = dPhi_Minus
		ELSE: Message, 'DMINUS had incorrect number of elements.'
	ENDCASE

	;Move PHI to center of bin
	phOut  = phOut + Rebin( (dPlus - dMinus) / 2.0, nPhi, nTheta)
	delta  = dPlus + dMinus

;-----------------------------------------------------
; Set Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Set data
	*self.phi    = tf_degrees ? Temporary(phOut) : Temporary(phOut) * !radeg
	*self.dPhi   = tf_degrees ? Temporary(delta) : Temporary(delta) * !radeg
	self.degrees = 1B
end


;+
;   Set object properties.
;
; :Params:
;       DISTFN:             in, required, type=NxMxL float
;                           Three-dimensional distribution function.
;       PHI:                in, required, type=NxM float
;                           Azimuthal location of each pixel in `DISTFN`.
;       THETA:              in, required, type=NxM float
;                           Polar location of each pixel in `DISTFN`.
;       ENERGY:             in, required, type=Lx1 float
;                           Energy location of each pixel in `DISTFN`.
;-
PRO MrDist3D::SetTheta, theta, dTheta_Minus, dTheta_Plus, $
DEGREES=degrees, $
RADIANS=radians
	Compile_Opt idl2
	On_Error, 2
	
	;Unis
	IF N_Elements(degrees) eq 0 && N_Elements(radians) eq 0 then degrees = 1B
	tf_degrees = N_Elements(radians) GT 0 ? ~Keyword_Set(radians) : Keyword_Set(degrees)
	
	dims    = Size(*self.distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

;-----------------------------------------------------
; Phi \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Phi must be 2D: [nPhi, nTheta]
	nDims = Size(theta, /N_DIMENSIONS)
	dims  = Size(theta, /DIMENSIONS)
	
	;1D
	IF nDims EQ 1 THEN BEGIN
		;Reform from [nPhi] to [nPhi, nTheta]
		IF dims[0] NE nTheta $
			THEN Message, 'THETA must be the size of the second dimension of DISTFN.' $
			ELSE thOut = Rebin( Reform(theta, 1, nTheta), nPhi, nTheta )
		
	;2D
	ENDIF ELSE IF nDims EQ 2 THEN BEGIN
		IF (dims[0] NE nPhi) || (dims[1] NE nTheta) $
			THEN Message, 'THETA must have dimensions of [nPhi, nTheta].' $
			ELSE thOut = theta
	
	;3+D
	ENDIF ELSE BEGIN
		Message, 'THETA must have one or two dimensions.'
	ENDELSE

;-----------------------------------------------------
; Deltas \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Assume values given are the bin centers
	IF N_Elements(dTheta_Plus) EQ 0 && N_Elements(dTheta_Minus) EQ 0 THEN BEGIN
		dTheta_Plus  = Median(thOut[0,1:*] - thOut[0,*]) / 2.0
		dTheta_Minus = dTheta_Plus
	ENDIF
	
	;DPLUS
	;   - Get a vector of values
	CASE N_Elements(dTheta_Plus) OF
		0:      dPlus = Replicate(0.0, nTheta)
		1:      dPlus = Replicate(dTheta_Plus, nTheta)
		nTheta: dPlus = dTheta_Plus
		ELSE: Message, 'DPLUS had incorrect number of elements.'
	ENDCASE
	
	;DMINUS
	;   - Get a vector of values
	CASE N_Elements(dTheta_Minus) OF
		0:      dMinus = Replicate(0.0, nTheta)
		1:      dMinus = Replicate(dTheta_Minus, nTheta)
		nTheta: dMinus = dTheta_Minus
		ELSE: Message, 'DMINUS had incorrect number of elements.'
	ENDCASE
	
	;Move PHI to bottom of bin
	thOut  = thOut + Rebin( Reform( (dPlus - dMinus) / 2.0, 1, nTheta), nPhi, nTheta )
	delta  = dPlus + dMinus

;-----------------------------------------------------
; Set Properties \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Set data
	*self.theta  = tf_degrees ? Temporary(thOut) : Temporary(thOut) * !radeg
	*self.dTheta = tf_degrees ? Temporary(delta) : Temporary(delta) * !radeg
	self.degrees = 1B
end


;+
;   Convert the distribution function from one set of units to another.
;
; :Params:
;       FLUX:           in, required, type=fltarr
;                       Distribution function for which to convert units
;       TO_UNITS:       in, required, type=string
;                       Convert to these units.  Options are::
;                           'ENERGY'      - eV
;                           'EFLUX'       - eV / cm^2 / s / sr / eV    or    keV / cm^2 / s / sr / keV
;                           'DIFF FLUX'   - # / cm^2 / s / sr / keV
;                           'PSD'         - s^2 / km^6
;                           'DF'          - s^2 / km^6
;       FROM_UNITS:     in, optional, type=string, default='PSD'
;                       Name of the units of `FLUX`.
;
; :Keywords:
;       SPECIES:        in, optional, type=string, default='e'
;                       Particle species measured in the distribution. Particle mass
;                           is taken into account when converting 'PSD' <-> {'EFLUX' | 'DIFF FLUX'}
;       ENERGY:         in, optional, type=float
;                       Energy bins at which the distribution is taken.
;
; :Returns:
;       NEW_FLUX:       Distribution with new units.
;-
pro MrDist3D::SetUnits, to_units
	compile_opt idl2
	on_error, 2

	to = strupcase(to_units)

	;Mass relative to the proton mass.
	A = self.mass / MrConstants('m_H')

	;Conversion factor from energy flux to phase space density
	;   - mass^2 * (sr / eV^2)
	eflux_to_psd = A^2 * 0.5447 * 1.0e6
	
	dims    = Size(*self.distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]
	energy  = Rebin( Reform( *self.energy, 1, 1, nEnergy ), nPhi, nTheta, nEnergy )
	
	;Energy Flux
	if self.units eq 'EFLUX' then begin
		
		;Convert to:
		case to of
			'DIFF FLUX': new_flux = *self.distfn / energy*1e3
			'PSD':       new_flux = eflux_to_psd * *self.distfn / energy^2
			'DF':        new_flux = eflux_to_psd * *self.distfn / energy^2
			'EFLUX':     new_flux = *self.distfn
			else: message, 'Cannot convert from "' + from + '" to "' + to + '".'
		endcase
	
	;Differential flux
	endif else if self.units eq 'DIFF FLUX' then begin
		eflux = *self.distfn * energy * 1e-3
		
		;Convert to:
		case to of
			'EFLUX':     new_flux = *self.distfn
			'PSD':       new_flux = eflux_to_psd * eflux / energy^2
			'DF':        new_flux = eflux_to_psd * eflux / energy^2
			'DIFF FLUX': new_flux = *self.distfn
			else: message, 'Cannot convert from "' + from + '" to "' + to + '".'
		endcase
	
	;Phase space density
	endif else if self.units eq 'PSD' || self.units eq 'DF' then begin
		eflux = *self.distfn * energy * 1e-3
		
		;Convert to:
		case to of
			'EFLUX':     new_flux = *self.distfn / eflux_to_psd * energy^2
			'DIFF FLUX': new_flux = *self.distfn / eflux_to_psd * energy*1e3
			'PSD':       new_flux = *self.distfn
			'DF':        new_flux = *self.distfn
			else: message, 'Cannot convert from "' + from + '" to "' + to + '".'
		endcase
	
	;Invalid
	endif else begin
		message, 'Invalid units: "' + self.units + '".'
	endelse
	
	;Return results
	self.units  = to
	*self.distfn = temporary(new_flux)
end


;+
;   Set object properties.
;
; :Keywords:
;       ELEVATION:          in, optional, type=boolean
;                           If set, THETA is taken to be the elevation angle.
;       MASS:               in, optional, type=float
;                           Mass (kg) of the species represented in the distribution.
;       SPECIES:            in, optional, type=string
;                           Species of particle represented in the distribution. Options
;                               are: { 'e' | 'p' | 'H' | 'He' | 'O' }. This keyword
;                               sets the MASS property and takes precendence over the
;                               `MASS` keyword.
;-
pro MrDist3D::SetProperty, $
ELEVATION=elevation, $
MASS=mass, $
SPECIES=species
	compile_opt idl2
	on_error, 2
	
	if n_elements(elevation) gt 0 then self.elevation = keyword_set(elevation)
	if n_elements(mass)      gt 0 then self.mass      = mass
	if n_elements(species)   gt 0 then self.mass      = MrConstants('m_' + species)
end


;+
;   Compute the temperature from the 2nd moment of the distribution (pressure)
;
; :Returns:
;       T:          out, required, type=float
;                   Temperature tensor calculated from the distribution.
;-
FUNCTION MrDist3D::Temperature
	Compile_Opt idl2
	On_Error, 2
	
	;Conversion from Kelvin to eV
	eV2K = 11600.0
	
	;Compute the pressure
	N = self -> Density()
	P = self -> Pressure()
	
	;Apply the equation of state
	;   - PV = NkT
	;   - T = kP/n
	;   - N = 1/cm^3
	;   - P = nPa
	;   - 1e-15 converts to Kelvin
	T = P / ( 1e15 * MrConstants('k_B') * N )
	T /= eV2K

	RETURN, T
END


;+
;   Reduce a single 3D distribution function to 2D by averaging over energy.
;
; :Categories:
;   Distribution Function
;
;
; :Params:
;       THETA:          out, optional, type=fltarr
;                       Polar angle bins (degrees) of the distribution.
;       PHI:            out, optional, type=fltarr
;                       Azimuthal angle bins (degrees) of the distribution.
;
; :Keywords:
;       E_RANGE:        in, optional, type=fltarr(2), default=[min, max]
;                       The range in energy, in electron volts (eV) over which to average.
;       NPHI_BINS:      in, optional, type=integer
;                       Number of phi bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       NTHETA_BINS:    in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;
; :Returns:
;       DIST1D:         out, required, type=fltarr(L)
;                       The 1D distribution with size nEnergy.
;-
function MrDist3D::ThetaPhi, phi, theta, dPhi, dTheta, $
E_RANGE=E_Range, $
NPHI_BINS=nPhi_bins, $
NTHETA_BINS=nTheta_bins
	compile_opt idl2
	on_error, 2

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Dimension sizes
	dims    = size(*self.distFn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

	;Energy range over which to average
	dE = self -> DeltaE()
	if n_elements(ERange) eq 0 then begin
		ERange = [ (*self.energy)[0]  - (*self.energy)[0]  * dE[0]  / exp(1), $
		           (*self.energy)[-1] + (*self.energy)[-1] * dE[-1] / exp(1) ]
	endif
	
	;Bins & mass
	phi_range   = [-180, 180]
	theta_range = [0, 180]
	if n_elements(nTheta_bins) eq 0 then nTheta_bins = nTheta_bins
	if n_elements(nPhi_bins)   eq 0 then nPhi_bins   = nPhi_bins

;-----------------------------------------------------
; Coordinate Space \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Combine data into 3-vectors of [phi, theta, energy]
	coords = transpose( [ [reform( rebin( *self.phi,   nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ], $
	                      [reform( rebin( *self.theta, nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ], $
	                      [reform( rebin( reform( alog10( *self.energy ), 1, 1, nEnergy ), nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ] ] )

;-----------------------------------------------------
; Weight \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	;Weight
	;   - w = v^2 * sin(theta) * dv * dTheta * dPhi
	;   - Applied as total(w*psd)/total(w)
	;   - Units of w cancel, so does not need same units as psd.
	;   - dTheta and dPhi are constants, so will factor out
	;   - Elements from two different old theta bins can fall into a new theta bin,
	;     so the sin(theta) term does not cancel.
	;   => w = v^2 * sin(theta) * dv 
	;
	
	;Velocity element
	vsqr     = 2.0 * *self.energy / self.mass
	dv       = sqrt( 1.0 / (2.0 * self.mass * *self.energy ) ) * dE
	
	
	if self.elevation $
		then wTheta = cos( *self.theta * !dtor) $
		else wTheta = sin( *self.theta * !dtor)
	
	;Weight function
	weight = rebin( reform(vsqr, 1, 1, nEnergy), nPhi, nTheta, nEnergy ) * $
	         rebin( wTheta,                      nPhi, nTheta, nEnergy ) * $
	         rebin( reform(dv, 1, 1, nEnergy),   nPhi, nTheta, nEnergy )

;-----------------------------------------------------
; Re-bin Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Locate data within new bins
	;   - Keep the number of Theta bins the same
	;   - Average all PHI values within [-DELTA, DELTA] in each theta bin
	cHist  = hist_nd(coords, $
	                 MIN             = [ phi_range[0], theta_range[0], E_Range[0]], $
	                 MAX             = [ phi_range[1], theta_range[1], E_Range[1]], $
	                 NBINS           = [   nPhi_bins,     nTheta_bins,          1], $
	                 REVERSE_INDICES = ri)

;-----------------------------------------------------
; Reduce to 2D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Allocate memory to reduced 2D distribution
	dist2D = fltarr(nPhi_bins, nTheta_bins)

	;Loop over bins
	for k = 0, n_elements(chist) - 1 do begin
		;Skip empty bins
		if ri[k] eq ri[k+1] then continue

		;Source indices
		inds = ri[ri[k]:ri[k+1]-1]
		isrc = array_indices([nPhi,nTheta,nEnergy], inds, /DIMENSIONS)
		
		;Destination indices
		idest = array_indices([nPhi,nTheta], k, /DIMENSIONS)
		
		;Weight
		w = weight[isrc[0,*], isrc[1,*], isrc[2,*]]

		;Re-bin
		dist2D[idest[0], idest[1]] = total( reform((*self.distFn)[isrc[0,*], isrc[1,*], isrc[2,*]]) * w ) / total(w)
	endfor

;-----------------------------------------------------
; Output Bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Phi bins
	dPhi  = (phi_range[1] - phi_range[0]) / nPhi_bins
	phi   = phi_range[0] + dPhi*findgen(nPhi_bins) + dPhi/2.0
	dPhi /= 2.0
	
	;Theta bins
	dTheta  = (theta_range[1] - theta_range[0]) / nTheta_bins
	theta   = theta_range[0] + dTheta*findgen(nTheta_bins) + dTheta/2.0
	dTheta /= 2.0

	;Return the 2D distribution
	return, dist2D
end


;+
;   Reduce the 3D distribution function to a 1D distribution in polar angle.
;
; :Params:
;       THETA:          out, optional, type=fltarr
;                       Polar angle (degrees) bins of the distribution.
;       ENERGY:         out, optional, type=fltarr
;                       Energy (eV) bins of the distribution.
;
; :Keywords:
;       NE_BINS:        in, optional, type=integer
;                       Number of energy bins in the reduced distribution. The default
;                           is to use the same bins and the original distribution.
;       NTHETA_BINS:    in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       PHI_RANGE:      in, optional, type=fltarr(2), default=[0.0\, 360.0]
;                       The range in azimuthal angle (degrees) over which to average.
;-
function MrDist3D::ThetaE, theta, energy, dTheta, dE, $
NE_BINS=nE_bins, $
NTHETA_BINS=nTheta_bins, $
PHI_RANGE=phi_range
	compile_opt idl2
	on_error, 2

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Dimension sizes
	dims    = size(*self.distFn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

	;Amount to average over PHI
	theta_range = [0.0, 180.0]
	if n_elements(phi_range)   eq 0 then phi_range   = [-180.0, 180.0]
	if n_elements(nTheta_bins) eq 0 then nTheta_bins = nTheta
	if n_elements(nE_bins)     gt 0 then MrPrintF, 'logwarn', 'TODO: Allow nE_bins to change.'
	
	;Keep the number of energy bins the same.
	;   - Do not use NE_BINS to prevent cyclic messages there ---^
	nBins = nEnergy

;-----------------------------------------------------
; Coordinate Space \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Combine data into 3-vectors of [phi, theta, energy]
	coords = transpose( [ [reform(*self.phi,   nPhi*nTheta)], $
	                      [reform(*self.theta, nPhi*nTheta)] ] )

;-----------------------------------------------------
; Weight \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Weight
	;   - Weight by size of old element
	;   - phi element: v*sin(theta)*dphi.
	;   - Weight is applied as total(w*psd)/total(w).
	;   - v is constant, so cancels out of w.
	;   - dphi is constant, so concels out of w
	;
	
	;Weight function
	if self.elevation $
		then weight = cos( *self.theta * !dtor ) $
		else weight = sin( *self.theta * !dtor )

;-----------------------------------------------------
; Re-bin Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;Locate data within new bins
	;   - Keep the number of Theta bins the same
	;   - Average all PHI values within [-DELTA, DELTA] in each theta bin
	cHist  = hist_nd( coords, $
	                  MIN             = [ phi_range[0], theta_range[0] ], $
	                  MAX             = [ phi_range[1], theta_range[1] ], $
	                  NBINS           = [            1,    nTheta_bins ], $
	                  REVERSE_INDICES = ri)

;-----------------------------------------------------
; Reduce to 2D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Allocate memory to reduced 2D distribution
	dist2D = fltarr(nTheta_bins, nBins)
	
	;Loop over bins
	for j = 0, nEnergy - 1 do begin
		for k = 0, n_elements(chist) - 1 do begin
			;Skip empty bins
			if ri[k] eq ri[k+1] then continue

			;Source indices
			count = ri[k+1] - ri[k]
			inds  = ri[ri[k]:ri[k+1]-1]
			isrc  = array_indices([nPhi,nTheta], inds, /DIMENSIONS)
			
			;Weight
			w = weight[isrc[0,*], isrc[1,*], replicate(j, count)]

			;Re-bin
			dist2D[k,j] = total( (*self.distFn)[isrc[0,*], isrc[1,*], replicate(j, count)] * w ) / total(w)
		endfor
	endfor

;-----------------------------------------------------
; Output Bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Create bins
	dTheta  = (theta_range[1] - theta_range[0]) / nTheta_bins
	theta   = theta_range[0] + dTheta*findgen(nTheta_bins) + dTheta/2.0
	dTheta /= 2.0
	
	energy = *self.energy

	;Return the 2D distribution
	return, dist2D
end


;+
;   Reduce a single 3D distribution function to 1D by averaging over
;   energy and azimuth.
;
; :Params:
;       THETA:          out, optional, type=fltarr
;                       Angles (degrees) of the polar bin centers of the distribution.
;       DTHETA:         out, optional, type=float
;                       Half-width of the `THETA` bins.
;
; :Keywords:
;       E_RANGE:        in, optional, type=fltarr(2), default=[min, max]
;                       The range in energy, in electron volts (eV) over which to average.
;       NTHETA_BINS:    in, optional, type=integer
;                       Number of polar angle bins in the reduced distribution. The
;                           default is to use the same bins and the original distribution.
;       PHI_RANGE:      in, optional, type=fltarr(2), default=[0.0\, 360.0]
;                       The range in azimuthal angle (degrees) over which to average.
;
; :Returns:
;       DIST1D:         out, required, type=fltarr(L)
;                       The 1D distribution in polar angle.
;-
function MrDist3D::ThetaSpec, theta, dTheta, $
E_RANGE=e_range, $
NTHETA_BINS=nTheta_bins, $
PHI_RANGE=phi_range
	compile_opt idl2
	on_error, 2

;-----------------------------------------------------
; Defaults \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------

	;Constants
	eV_to_J = 1.602e-19
	m_to_km = 1e-3

	;Dimension sizes
	dims    = size(*self.distFn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]

	;PA range over which to average
	theta_range = [0, 180]
	if n_elements(phi_range)   eq 0 then phi_range   = [0.0, 360.0]
	if n_elements(nTheta_bins) eq 0 then nTheta_bins = nTheta

	;Energy range over which to average
	dE = self -> DeltaE()
	if n_elements(e_range) eq 0 then begin
		e_range = [ (*self.energy)[0]  - (*self.energy)[0]  * dE[0]  / exp(1), $
		            (*self.energy)[-1] + (*self.energy)[-1] * dE[-1] / exp(1) ]
	endif

;-----------------------------------------------------
; Coordinate Space \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Combine data into 3-vectors of [phi, theta, energy]
	coords = transpose( [ [reform( rebin( *self.phi,   nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ], $
	                      [reform( rebin( *self.theta, nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ], $
	                      [reform( rebin( reform( *self.energy, 1, 1, nEnergy ), nPhi, nTheta, nEnergy ), nPhi*nTheta*nEnergy ) ] ] )

;-----------------------------------------------------
; Weight \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	;
	;   - w = v^2 * sin(theta) * dv * dTheta * dPhi
	;   - Applied as total(w*psd)/total(w)
	;   - Units of w cancel, so does not need same units as psd.
	;   - dTheta and dPhi are constants, so will factor out
	;   - Elements from two different old theta bins can fall into the same
	;     new theta bin, so the sin(theta) term does not cancel.
	;   => w = v^2 * sin(theta) * dv 
	;
	
	;Velocity element
	v_sqr    = *self.energy * ( (2*eV_to_J*m_to_km^2)/self.mass )
	dv       = sqrt(1.0 / (2.0 * self.mass * *self.energy) ) * dE * sqrt(eV_to_J)*m_to_km
	sinTheta = sin(*self.theta * !dtor)
	
	;Weight function
	weight = rebin( reform(v_sqr, 1, 1, nEnergy), nPhi, nTheta, nEnergy ) * $
	         rebin( sinTheta,                     nPhi, nTheta, nEnergy ) * $
	         rebin( reform(dv, 1, 1, nEnergy),    nPhi, nTheta, nEnergy )
;	weight = replicate(1.0, nPhi, nTheta, nEnergy)

;-----------------------------------------------------
; Re-bin Data \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Locate data within new bins
	;   - Keep the number of Theta bins the same
	;   - Average all relevant PHI and ENERGY values.
	cHist  = hist_nd(coords, $
	                 MIN             = [phi_range[0], theta_range[0], e_range[0]], $
	                 MAX             = [phi_range[1], theta_range[1], e_range[1]], $
	                 NBINS           = [           1,    nTheta_bins,          1], $
	                 REVERSE_INDICES = ri)

;-----------------------------------------------------
; Reduce to 1D \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Allocate memory to reduced 2D distribution
	dist1D = fltarr(nTheta_bins)

	;Loop over bins
	for k = 0, n_elements(chist) - 1 do begin
		;Skip empty bins
		if ri[k] eq ri[k+1] then continue

		;Source indices
		inds = ri[ri[k]:ri[k+1]-1]

		;Re-bin
		dist1D[k] = total( (*self.distFn)[inds] * weight[inds] ) / total(weight[inds])
	endfor

;-----------------------------------------------------
; Output Bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
;-----------------------------------------------------
	
	;Create bins
	dTheta  = (theta_range[1] - theta_range[0]) / nTheta_bins
	theta   = theta_range[0] + dTheta*findgen(nTheta_bins) + dTheta/2.0
	dTheta /= 2.0

	;Return the 1D distribution
	return, dist1D
end


;+
;   Compute the 1th order moment, velocity.
;
; :Returns:
;       V:          out, required, type=float
;                   Velocity (km/s) of particles calculated from the distribution.
;-
FUNCTION MrDist3D::Velocity
	Compile_Opt idl2
	On_Error, 2
	
	;Must integrate over phase space
	IF self.units NE 'PSD' THEN Message, 'Units must be "PSD". They are "' + self.units + '".'

	;Convert to radians
	deg2rad = !dpi / 180D
	dims    = Size(*self.distfn, /DIMENSIONS)
	
	;Integrate phi
	;   - Integral( V * fdist * d^3V )
	;   - V    = (Vx, Vy, Vz)
	;   - Vx   = V * Sin(Theta) * Cos(Phi)
	;   - Vy   = V * Sin(Theta) * Sin(Phi)
	;   - Vz   = V * Cos(Theta)
	;   - (d^3)V = V^2 * Sin(Theta) * dV * dTheta * dPhi
	dphi    = Median( (*self.phi)[1:*,*] - (*self.phi)[0:-2,*]) * deg2rad
	phi     = Rebin( *self.phi, dims ) * deg2rad
	fx_temp = Total( *self.distfn * Cos(phi), 1) * dPhi
	fy_temp = Total( *self.distfn * Sin(phi), 1) * dPhi
	fz_temp = Total( *self.distfn, 1) * dPhi

	;Integrate theta
	theta   = Rebin(Reform((*self.theta)[0,*]), dims[1], dims[2]) * deg2rad
	dTheta  = Median( theta[1:-1,*] - theta[0:-2,*] )
	fx_temp = Total( fx_temp * Sin(theta)^2, 1) * dTheta
	fy_temp = Total( fy_temp * Sin(theta)^2, 1) * dTheta
	fz_temp = Total( fz_temp * Sin(theta)*Cos(Theta), 1) * dTheta
	
	;Integrate velocity
	q    = MrConstants('q')
	eV2J = MrConstants('eV2J')
	
	;
	; Convert energy to velocity
	;   - E  = 1/2 * m * v^2
	;   - dE = m * v * dv
	;   - v  = Sqrt( 2 * E / m )
	;   - dv = dE / (m v)
	;
	
	;Convert energy from eV to J
	v    = sqrt(2.0 * eV2J * *self.energy / self.mass)
	dE   = self -> DeltaE()
	dv   = eV2J * dE / (self.mass * v)
	
	;Spacecraft potential correction
	IF N_Elements(*self.Vsc) GT 0 THEN BEGIN
		;Sign is charge dependent. E = q * V = 1/2 * m * v^2
		sgn  = Round(self.mass / MrConstants('m_p')) EQ 0 ? -1.0 : 1.0
		vsc  = Sqrt( 2.0 * q * (*self.Vsc) / self.mass )
		
		Vx = Total( fx_temp * v * (v^2 + sgn*vsc^2) * dv ) * 1e12
		Vy = Total( fy_temp * v * (v^2 + sgn*vsc^2) * dv ) * 1e12
		Vz = Total( fz_temp * v * (v^2 + sgn*vsc^2) * dv ) * 1e12
	
	ENDIF ELSE BEGIN
		;V is in m/s while f is in s^3/cm^6
		;   - f --> s^3/cm^6 --> s^3/m^6 * 1e12
		Vx = Total( fx_temp * v^3 * dv ) * 1e12
		Vy = Total( fy_temp * v^3 * dv ) * 1e12
		Vz = Total( fz_temp * v^3 * dv ) * 1e12
	ENDELSE
	
	;Compute the number density (1/cm^3)
	;   - N --> 1/cm^3 --> 1/m^3 * 1e6
	N = self -> Density() * 1e6
	
	;Form vector
	;   - V --> m/s --> km/s * 1e-3
	;   - Negative sign because the angles represent look directions,
	;     not incident directions
	V = -1e-3 * [Vx, Vy, Vz] / N

	RETURN, V
END


;+
;   Compute the size of velocity-space volume elements.
;
;       dV = v^2 * Sin(theta) * dv * dTheta * dPhi
;
; :Returns:
;       DVOL:       out, required, type=NxMxL FltArr
;                   Size of each volume element (m^3 / s^3 * sr) of the distribution.
;-
FUNCTION MrDist3D::VolumeElement
	Compile_Opt idl2
	On_Error, 2
	
	dims    = Size(*self.distfn, /DIMENSIONS)
	nPhi    = dims[0]
	nTheta  = dims[1]
	nEnergy = dims[2]
	
	;Constants
	deg2rad = !pi / 180.0
	eV2J    = MrConstants('eV2J')

	;Center angular bins
	phi   = *self.phi   + Rebin(*self.dphi, nPhi, nTheta) / 2.0
	theta = *self.theta + Rebin( Reform(*self.dtheta, 1, nTheta), nPhi, nTheta ) / 2.0
	
	;Compute velocity
	v  = Sqrt( 2.0 * eV2J * (*self.energy + *self.denergy/2.0) / self.mass )
	dv = Sqrt( 2.0 * eV2J * *self.denergy / self.mass )
	
	;Rebin to same size
	;   - [nPhi, nTheta, nEnergy]
	phi    = Rebin( phi   * deg2rad,   nPhi, nTheta, nEnergy )
	theta  = Rebin( theta * deg2rad, nPhi, nTheta, nEnergy ) * deg2rad
	v      = Rebin( Reform( v,  1, 1, nEnergy ), nPhi, nTheta, nEnergy )
	dv     = Rebin( Reform( dv, 1, 1, nEnergy ), nPhi, nTheta, nEnergy )
	dPhi   = Rebin( *self.dphi * deg2rad, nPhi, nTheta, nEnergy )
	dTheta = Rebin( Reform( *self.dTheta * deg2rad, 1, nTheta ), nPhi, nTheta, nEnergy )
	
	;Compute volume element
	IF self.elevation $
		THEN dVol = v^2 * Cos(theta) * dv * dTheta * dPhi $
		ELSE dVol = v^2 * Sin(theta) * dv * dTheta * dPhi
	
	RETURN, dVol
END


;+
;   The class definition statement.
;
; :Params:
;       CLASS:          out, optional, type=structure
;-
pro MrDist3D__DEFINE
	compile_opt idl2
	
	class = { MrDist3D, $
	          distfn:      ptr_new(), $
	          denergy:     ptr_new(), $
	          dphi:        ptr_new(), $
	          dtheta:      ptr_new(), $
	          theta:       ptr_new(), $
	          phi:         ptr_new(), $
	          energy:      ptr_new(), $
	          Vsc:         ptr_new(), $
	          elevation:   0B, $
	          degrees:     0B, $
	          units:       '', $
	          mass:        0.0 $
	        }
end