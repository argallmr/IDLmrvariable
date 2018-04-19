; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_EDI_Amb
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
;   Generate a plot of EDI Ambient quantities:
;        1. FGM Bxyz, |B|
;        2. EDP Exyz
;        3. Flux 0-degrees
;        4. PSD Flux 0-degrees
;        5. Flux 180-degrees
;        6. PSD Flux 180-degrees
;        7. Phi: All channels & pitch angles
;        8. GPD
;        9. Theta: All channels & pitch angles
;       10. PAD
;
; :Categories:
;   MMS
;
; :Params:
;       SC:         in, required, type=string
;                   Spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4' }
;       MODE:       in, required, type=string, default='srvy'
;                   Data telemetry rate of the data. Options are: { 'slow' | 'fast' | 'srvy' | 'brst' }
;
; :Keywords:
;       FGM_INSTR:  in, optional, type=string, default='fgm'
;                   FGM instrument to use. Options are: { 'afg' | 'dfg' | 'fgm' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'ql' | 'l2'}
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;       TRANGE:     in, optional, type=string/strarr(2), default=MrVar_GetTRange()
;                   The start and end times of the data interval to be plotted, formatted
;                       as 'YYYY-MM-DDThh:mm:ss'
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
FUNCTION MrMMS_Plot_EDI_Amb, sc, mode, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
TRANGE=trange
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		IF N_Elements(win) GT 0 THEN Obj_Destroy, win
		MrPrintF, 'LogErr'
		RETURN, !Null
	ENDIF
	
	tf_load = ~Keyword_Set(no_load)
	IF N_Elements(mode)   EQ 0 THEN mode       = 'srvy'
	IF N_Elements(level)  EQ 0 THEN level      = 'l2'
	IF N_Elements(nfft)   EQ 0 THEN nfft       = 512
	IF N_Elements(nshift) EQ 0 THEN nshift     = nfft/2
	IF N_Elements(trange) GT 0 THEN MrVar_SetTRange, trange
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	;EDI
	instr   = 'edi'
	coords  = 'gse'
	optdesc = ['amb', 'amb-pm2', 'amb-alt-cc', 'amb-alt-oc', 'amb-alt-ooc', 'amb-alt-oob']
	
	;FGM
	fgm_instr  = mode EQ 'brst' ? 'fgm' : 'dfg'
	fgm_level  = mode EQ 'brst' ? 'l2'  : 'l2pre'
	fgm_mode   = mode EQ 'brst' ? mode : 'srvy'
	fgm_coords = coords EQ 'dbcs'  ? 'dmpa' : coords
	
	;EDP
	edp_instr  = 'edp'
	edp_mode   = mode eq 'brst' ? mode : 'fast'
	edp_coords = coords EQ 'dbcs'  ? 'dsl'  : coords

	;Source names
	;   - Different version of L2Pre use different naming conventions
	fgm_b_vname     = StrJoin( [sc, fgm_instr, 'b',    fgm_coords, fgm_mode, fgm_level], '_' )
	fgm_bvec_vname  = StrJoin( [sc, fgm_instr, 'bvec', fgm_coords, fgm_mode, fgm_level], '_' )
	fgm_bmag_vname  = StrJoin( [sc, fgm_instr, 'bmag', fgm_coords, fgm_mode, fgm_level], '_' )
	e_vname       = StrJoin( [sc, edp_instr, 'dce',  edp_coords, edp_mode, level], '_' )
	
	channels      = ['1', '2', '3', '4']
	flux0_vname   = sc + '_' + instr + '_' + 'flux' + channels + '_'                + '0'   + '_' + mode + '_' + level
	flux180_vname = sc + '_' + instr + '_' + 'flux' + channels + '_'                + '180' + '_' + mode + '_' + level
	traj0_vname   = sc + '_' + instr + '_' + 'traj' + channels + '_' + coords + '_' + '0'   + '_' + mode + '_' + level
	traj180_vname = sc + '_' + instr + '_' + 'traj' + channels + '_' + coords + '_' + '180' + '_' + mode + '_' + level

	;Derived names
	f0_vname       = StrJoin( [sc, instr, 'flux',   '0', mode, level], '_' )
	f180_vname     = StrJoin( [sc, instr, 'flux', '180', mode, level], '_' )
	psd0_vname     = StrJoin( [sc, instr, 'psd',    '0', mode, level], '_' )
	psd180_vname   = StrJoin( [sc, instr, 'psd',  '180', mode, level], '_' )
	phi0_vname     = sc + '_' + instr + '_' + 'phi'   + channels + '_' + 'fac' + '_' +   '0' + '_' + mode + '_' + level
	theta0_vname   = sc + '_' + instr + '_' + 'theta' + channels + '_' + 'fac' + '_' +   '0' + '_' + mode + '_' + level
	phi180_vname   = sc + '_' + instr + '_' + 'phi'   + channels + '_' + 'fac' + '_' + '180' + '_' + mode + '_' + level
	theta180_vname = sc + '_' + instr + '_' + 'theta' + channels + '_' + 'fac' + '_' + '180' + '_' + mode + '_' + level
	pad_vname      = StrJoin( [sc, instr, 'pad', mode, level], '_' )
	gpd_vname      = StrJoin( [sc, instr, 'gpd', mode, level], '_' )

;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, fgm_mode, $
		                     INSTR     = fgm_instr, $
		                     LEVEL     = fgm_level, $
		                     VARFORMAT = fgm_b_vname
		
		IF fgm_level NE 'l2' && ~MrVar_IsCached(fgm_b_vname) THEN BEGIN
			fgm_b_vname     = StrJoin( [sc, fgm_instr,        fgm_mode, fgm_level, fgm_coords], '_' )
			fgm_bvec_vname  = StrJoin( [sc, fgm_instr, 'vec', fgm_mode, fgm_level, fgm_coords], '_' )
			fgm_bmag_vname  = StrJoin( [sc, fgm_instr, 'mag', fgm_mode, fgm_level, fgm_coords], '_' )
			
			MrMMS_FGM_Load_Data, sc, fgm_mode, $
			                     INSTR     = fgm_instr, $
			                     LEVEL     = fgm_level, $
			                     VARFORMAT = fgm_b_vname
			
			IF ~MrVar_IsCached(fgm_b_vname) $
				THEN Message, 'Unknown FGM variable naming convention.'
		ENDIF

		;EDP
		MrMMS_Load_Data, sc, edp_instr, edp_mode, level, $
		                 OPTDESC   = 'dce', $
		                 VARFORMAT = e_vname

		;EDI
		oEDI = MrMMS_EDI_Dist(sc, mode, optdesc, $
		                      COORD_SYS=coords)
	ENDIF
stop
	;Determine which EDI file was loaded
	fnames = MrMMS_Get_FileNames(sc, 'edi', mode, level, OPTDESC=optdesc)
	MrMMS_Parse_Filename, fnames, OPTDESC=optdesc
	iUniq = Uniq(optdesc, Sort(optdesc))
	IF N_Elements(iUniq) NE 1 THEN BEGIN
		MrPrintF, 'LogWarn', 'More than one EDI file type found.'
		MrPrintF, 'LogWarn', '   ' + '[' + StrJoin(optdesc[iUniq], ', ') + ']'
		MrPrintF, 'LogWarn', '   Choosing "' + optdesc[0] + '".'
	ENDIF
	optdesc = optdesc[0]

;-------------------------------------------
; EDI: PAD & GPD ///////////////////////////
;-------------------------------------------
	IF ~Obj_Valid(oEDI) THEN BEGIN
		oEDI  = MrMMS_EDI_Dist()
		oEDI -> SetData, flux0_vname,   traj0_vname,   channels, Replicate(3, 4), Replicate(0, 4)
		oEDI -> SetData, flux180_vname, traj180_vname, channels, Replicate(3, 4), Replicate(180, 4)
	ENDIF
	
	;Field-Aligned Coordinates
	oEDI -> SetFAC, fgm_bvec_vname, e_vname, 'EXB'

	;Distributions
	oPAD = oEDI -> ThetaSpec(/CACHE, NAME=pad_vname, THETA_RANGE=[0,180], DTHETA=2.5)
	oGPD = oEDI -> PhiSpec(/CACHE, NAME=gpd_vname)

;-------------------------------------------
; EDI: Flux ////////////////////////////////
;-------------------------------------------
	;AMB
	IF optdesc EQ 'amb' THEN BEGIN
		;Mean of middle two channels
		;   - B-field is centered between channels 2 & 3
		oFlux2_0 = MrVar_Get(flux0_vname[1])
		oFlux3_0 = MrVar_Get(flux0_vname[2])

		oFlux_0  = (oFlux2_0 + oFlux3_0) / 2.0
		oFlux_0 -> SetName, f0_vname
		oFlux_0 -> Cache
		
		;Mean of middle two channels
		;   - B-field is centered between channels 2 & 3
		oFlux2_180 = MrVar_Get(flux180_vname[1])
		oFlux3_180 = MrVar_Get(flux180_vname[2])
		oFlux_180  = (oFlux2_180 + oFlux3_180) / 2.0
		oFlux_180 -> SetName, f180_vname
		oFlux_180 -> Cache
	
	;AMB-PM2
	ENDIF ELSE IF optdesc EQ 'amb-pm2' THEN BEGIN
		;B-field is centered on channel 1
		oFlux1_0 = MrVar_Get(flux0_vname[1])
		oFlux_0  = oFlux1_0 -> Copy(f0_vname, /CACHE)
		
		;B-field is centered on channel 1
		oFlux1_180 = MrVar_Get(flux180_vname[1])
		oFlux_180  = oFlux1_0 -> Copy(f180_vname, /CACHE)
	ENDIF ELSE BEGIN
		Message, 'Optdesc not implemented: "' + optdesc + '".'
	ENDELSE

;-------------------------------------------
; EDI: Trajectories ////////////////////////
;-------------------------------------------
	colors = MrDefaultColor(NCOLORS=8)
	
	;Get the data
	FOR i = 0, 3 DO BEGIN
		oEDI -> GetData, !Null,   oPhi0,   oTheta0, CHANNEL=i+1, PA_STATE=0
		oEDI -> GetData, !Null, oPhi180, oTheta180, CHANNEL=i+1, PA_STATE=180
		
		;Set names
		oPhi0     -> SetName, phi0_vname[i]
		oPhi180   -> SetName, phi180_vname[i]
		oTheta0   -> SetName, theta0_vname[i]
		oTheta180 -> SetName, theta180_vname[i]
		
		;Cache
		oPhi0     -> Cache
		oPhi180   -> Cache
		oTheta0   -> Cache
		oTheta180 -> Cache
		
		;Attributes
		oPhi0['COLOR']     = colors[i]
		oPhi0['SYMBOL']    = 3
		oPhi0['LINESTYLE'] = 'None'
		
		oPhi180['COLOR']     = colors[i+4]
		oPhi180['SYMBOL']    = 3
		oPhi180['LINESTYLE'] = 'None'
		
		oTheta0['COLOR']     = colors[i]
		oTheta0['SYMBOL']    = 3
		oTheta0['LINESTYLE'] = 'None'
		
		oTheta180['COLOR']     = colors[i+4]
		oTheta180['SYMBOL']    = 3
		oTheta180['LINESTYLE'] = 'None'
	ENDFOR
	
	;Clean Up
	Obj_Destroy, oEDI

;-------------------------------------------
; EDI: PSD /////////////////////////////////
;-------------------------------------------
	CASE optdesc OF
		'amb':     idx = 1
		'amb-pm2': idx = 0
		ELSE:      Message, 'OPTDESC not implemented yet: "' + optdesc + '".'
	ENDCASE
	
	;PSD-0
	oFlux0 = MrVar_Get(flux0_vname[idx])
	oPSD0  = oFlux0 -> Spectrogram( nfft, nshift, $
	                                /CACHE, $
	                                NAME   = psd0_vname, $
	                                WINDOW = 'hanning' )
	
	;PSD-180
	oFlux_180 = MrVar_Get(flux180_vname[idx])
	oPSD_180  = oFlux_180 -> Spectrogram( nfft, nshift, $
	                                      /CACHE, $
	                                      NAME   = psd180_vname, $
	                                      WINDOW = 'hanning' )

;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	;BMAG
	oB = MrVar_Get(fgm_b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, mode, level, optdesc], ' ' ) )
	
	;E FIELD
	oE = MrVar_Get(e_vname)
	oE['TITLE'] = 'E!C(mV/m)'
	
	;FLUX RANGE
	oFlux0     = MrVar_Get(f0_vname)
	oFlux180   = MrVar_Get(f180_vname)
	i0         = oFlux0   -> Where(0, /GREATER, COUNT=n0)
	i180       = oFlux180 -> Where(0, /GREATER, COUNT=n180)
	flux_range = [ Min( [oFlux0[i0].min, oFlux180[i180].min] ), $
	               Max( [oFlux0[i0].max, oFlux180[i180].max] ) ]

	;FLUX 0
	iGood = oFlux0 -> Where(0, /GREATER, COUNT=nGood)
	oFlux0['AXIS_RANGE'] = flux_range
	oFlux0['LOG']        = 1B
	oFlux0['PLOT_TITLE'] = 'Electron Flux Parallel to B'
	oFlux0['TITLE']      = 'Flux!Cpa0!C(1/cm^2/s)'
	oFlux0['UNITS']      = '1/cm^2/s'
	
	;FLUX 180
	oFlux180['AXIS_RANGE'] = flux_range
	oFlux180['LOG']        = 1B
	oFlux180['PLOT_TITLE'] = 'Electron Flux Anti-Parallel to B'
	oFlux180['TITLE']      = 'Flux!Cpa180!C(1/cm^2/s)'
	oFlux180['UNITS']      = '1/cm^2/s'
	
	;PSD 0
	oPSD0 = MrVar_Get(psd0_vname)
	oPSD0['TITLE'] = 'PSD Flux!Cpa0!C(1/cm^2/s)^2/Hz'
	
	;PSD 180
	oPSD180 = MrVar_Get(psd180_vname)
	oPSD180['TITLE'] = 'PSD Flux!Cpa180!C(1/cm^2/s)^2/Hz'
	
	;PHI
	oPhi = MrVar_Get(phi0_vname[0])
	oPhi['AXIS_RANGE']   = [-180, 180]
	oPhi['TICKINTERVAL'] = 90.0
	oPhi['TITLE']        = 'Gyrophase!C(deg)'
	
	;THETA
	oTheta = MrVar_Get(theta0_vname[0])
	oTheta['AXIS_RANGE']   = [0, 180]
	oTheta['TICKINTERVAL'] = 90.0
	oTheta['TITLE']        = 'Pitch!C(deg)'
	
	;PAD
	oPAD   = MrVar_Get(pad_vname)
	oPAD['TITLE'] = 'Flux!C(1/cm^2/s)'
	oTheta = oPAD['DEPEND_1']
	oTheta['AXIS_RANGE']   = [0, 180]
	oTheta['TICKINTERVAL'] = 90.0
	oTheta['TITLE']        = 'Pitch!C(deg)'
	
	;GPD
	oGPD   = MrVar_Get(gpd_vname)
	oGPD['TITLE'] = 'Flux!C(1/cm^2/s)'
	oPhi = oGPD['DEPEND_1']
	oPhi['AXIS_RANGE']   = [-180, 180]
	oPhi['TICKINTERVAL'] = 90.0
	oPhi['TITLE']        = 'Gyrophase!C(deg)'
	
;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	trange = MrVar_GetTRange('SSM')
	IF trange[1] - trange[0] LE 60.0 THEN trange -= Floor(trange[0])

	;Plot data
	win = MrVar_PlotTS( [ fgm_b_vname, e_vname, f0_vname, psd0_vname, f180_vname, psd180_vname, $
	                      phi0_vname[0], gpd_vname, theta0_vname[0], pad_vname ], $
	                    /NO_REFRESH, $
	                    XSIZE = 680, $
	                    YSIZE = 700 )

	win = MrVar_OPlotTS( phi0_vname[0],   phi0_vname[1:3] )
	win = MrVar_OPlotTS( phi0_vname[0],   phi180_vname )
	win = MrVar_OPlotTS( theta0_vname[0], theta0_vname[1:3] )
	win = MrVar_OPlotTS( theta0_vname[0], theta180_vname )

	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 14]
	win    -> SetGlobal, XRANGE=trange
	win    -> Refresh

;-------------------------------------------
; Save Figure //////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
		ENDIF ELSE IF ~File_Test(output_dir, /DIRECTORY) THEN BEGIN
			MrPrintF, 'LogText', 'Creating directory: "' + output_dir + '".'
			File_MKDir, output_dir
		ENDIF
		
		;File name
		fname   = StrJoin( [sc, instr, mode, level, optdesc], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------
	RETURN, win
END