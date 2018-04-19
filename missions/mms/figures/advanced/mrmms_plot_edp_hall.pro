; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FGM_4sc
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
;   Generate a plot to provide an overview of reconnection quantities:
;       1. Ex, VixB, VexB, JxB/ne
;       2. Ey, VixB, VexB, JxB/ne
;       3. Ez, VixB, VexB, JxB/ne
;       4. ni, ne
;       5. Bxyz, |B|
;       6. J (moments)
;       7. Rho (Div E)
;
; :Categories:
;   MMS
;
; :Params:
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;
; :Keywords:
;       EPHDESC:    in, optional, type=string, default='ephts04d'
;                   Optional descriptor of the definitive ephemeris datatype to use.
;                       Options are: { 'epht89d' | 'epht89q' | 'ephts04d' | 'defeph' | 'predeph' }
;       LEVEL:      in, optional, type=string, default='l2'
;                   Data quality level. Options are: {'l1a' | 'l1b' | 'ql' | 'l2pre' | 'l2'}
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source CDF files.
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
;       2017/01/05  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_EDP_Hall, sc, mode, $
FGM_INSTR=fgm_instr, $
EPHDESC=ephdesc, $
LEVEL=level, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
NO_LOAD=no_load, $
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
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(trange)    GT 0 THEN MrVar_SetTRange, trange
	
	;EPHDESC
	IF N_Elements(ephdesc) EQ 0 THEN BEGIN
		IF N_Elements(level) EQ 0 $
			THEN ephdesc = 'ephts04d' $
			ELSE ephdesc = level EQ 'ql' ? 'predeph' : 'ephts04d'
	ENDIF
	IF mode EQ 'brst' THEN BEGIN
		IF (mrmms_get_filenames(sc, 'mec', mode, level, OPTDESC=ephdesc))[0] EQ '' $
			THEN mec_mode = 'srvy' $
			ELSE mec_mode = mode
	ENDIF
	
	;Constants
	q  = MrConstants('q')
	e0 = MrConstants('epsilon_0')
	
;-------------------------------------------
; Instrument Parameters ////////////////////
;-------------------------------------------
	;FGM
	CASE fgm_instr OF
		'afg': fgm_level = 'l2pre'
		'dfg': fgm_level = 'l2pre'
		'fgm': fgm_level = 'l2'
		'fsm': fgm_level = 'l3'
		ELSE: Message, 'Invalid FGM instrument: "' + fgm_instr + '".'
	ENDCASE
	fgm_coords  = MrIsMember(['dsl', 'dbcs'], coords) ? 'dmpa' : coords
	fgm_optdesc = fgm_instr EQ 'fsm' ? '8khz' : ''
	
	;EDP
	edp_mode    = mode EQ 'brst' ? mode : 'fast'
	edp_coords  = MrIsMember(['dmpa', 'dbcs'], coords) ? 'dsl' : coords
	edp_optdesc = 'dce'
	
	;FPI
	fpi_mode    = mode EQ 'brst' ? mode : 'fast'
	fpi_coords  = MrIsMember(['dmpa', 'dsl'], coords) ? 'dbcs' : coords
	des_optdesc = 'des-moms'
	dis_optdesc = 'dis-moms'
	
;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	
	;FGM
	b_vname    = StrJoin([sc, fgm_instr, 'b',    fgm_coords, mode, fgm_level], '_')
	bmag_vname = StrJoin([sc, fgm_instr, 'bmag', fgm_coords, mode, fgm_level], '_')
	bvec_vname = StrJoin([sc, fgm_instr, 'bvec', fgm_coords, mode, fgm_level], '_')
	
	;EDP
	all_sc   = 'mms' + ['1', '2', '3', '4']
	e_vnames = all_sc + '_' + StrJoin( ['edp', 'dce', edp_coords, edp_mode, level], '_' )
	
	;MEC
	IF Array_Equal(ephdesc EQ ['defeph', 'predeph'], 0) $
		THEN r_vnames = all_sc + '_' + StrJoin( ['mec', 'r', coords], '_' ) $
		ELSE r_vnames = StrUpCase(sc) + '_' + StrUpCase(ephdesc) + '_' + 'R'
	
	;FPI
	ni_vname = StrJoin([sc, 'dis', 'numberdensity',             fpi_mode], '_')
	vi_vname = StrJoin([sc, 'dis', 'bulkv',         fpi_coords, fpi_mode], '_')
	ne_vname = StrJoin([sc, 'des', 'numberdensity',             fpi_mode], '_')
	ve_vname = StrJoin([sc, 'des', 'bulkv',         fpi_coords, fpi_mode], '_')
	
	;Output names
	CASE StrUpCase(sc) OF
		'MMS1': e_vname = e_vnames[0]
		'MMS2': e_vname = e_vnames[1]
		'MMS3': e_vname = e_vnames[2]
		'MMS4': e_vname = e_vnames[3]
		ELSE: Message, 'Incorrect value for SC: "' + sc + '".'
	ENDCASE
	ex_vname    = e_vname + '_x'
	ey_vname    = e_vname + '_y'
	ez_vname    = e_vname + '_z'
	jxb_vname   = StrJoin( [sc, 'fpi', 'jxb', edp_coords, fpi_mode, level], '_' )
	jxbx_vname  = jxb_vname + '_x'
	jxby_vname  = jxb_vname + '_y'
	jxbz_vname  = jxb_vname + '_z'
	vixb_vname  = StrJoin( [sc, 'dis', 'vxb', edp_coords, edp_mode, level], '_' )
	vixbx_vname = vixb_vname + '_x'
	vixby_vname = vixb_vname + '_y'
	vixbz_vname = vixb_vname + '_z'
	vexb_vname  = StrJoin( [sc, 'des', 'vxb', edp_coords, edp_mode, level], '_' )
	vexbx_vname = vexb_vname + '_x'
	vexby_vname = vexb_vname + '_y'
	vexbz_vname = vexb_vname + '_z'
	j_vname     = StrJoin([sc, 'fpi', 'j', fpi_coords, fpi_mode], '_')
	rho_vname   = StrJoin( ['mms', 'edp', 'rho', mode, level], '_' )
		
;-------------------------------------------
; Get Data /////////////////////////////////
;-------------------------------------------
	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, mode, $
		                     INSTR   = fgm_instr, $
		                     LEVEL   = fgm_level, $
		                     OPTDESC = fgm_optdesc, $
		                     VARFORMAT = '*_b_'+fgm_coords+'_'+mode+'*'
		
		;E-Field
		MrMMS_Load_Data, '', 'edp', edp_mode, level, $
		                 OPTDESC   = edp_optdesc, $
		                 VARFORMAT = ['*dce_'+edp_coords+'*', '*dce_par*']

		;DES
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = des_optdesc, $
		                     VARFORMAT = ['*density_'+fpi_mode, '*bulkv_'+fpi_coords+'_'+fpi_mode]
		                     
		;DIS
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = dis_optdesc, $
		                     VARFORMAT = ['*density_'+fpi_mode, '*bulkv_'+fpi_coords+'_'+fpi_mode]
		
		;Ephemeris
		MrMMS_Load_Data, '', 'mec', mec_mode, level, $
		                 OPTDESC   = ephdesc, $
		                 VARFORMAT = r_vnames
	ENDIF
	
;-------------------------------------------
; Free Charge Density //////////////////////
;-------------------------------------------
	;Grab E & R
	oNe = MrVar_Get(ne_vname)
	oE1 = MrVar_Get(e_vnames[0])
	oE2 = MrVar_Get(e_vnames[1])
	oE3 = MrVar_Get(e_vnames[2])
	oE4 = MrVar_Get(e_vnames[3])
	oR1 = MrVar_Get(r_vnames[0])
	oR2 = MrVar_Get(r_vnames[1])
	oR3 = MrVar_Get(r_vnames[2])
	oR4 = MrVar_Get(r_vnames[3])
	
	;Interpolate to Ne on the selected spacecraft
	oE1 = oE1 -> Interpol(oNe)
	oE2 = oE2 -> Interpol(oNe)
	oE3 = oE3 -> Interpol(oNe)
	oE4 = oE4 -> Interpol(oNe)
	oR1 = oR1 -> Interpol(oNe)
	oR2 = oR2 -> Interpol(oNe)
	oR3 = oR3 -> Interpol(oNe)
	oR4 = oR4 -> Interpol(oNe)
	
	;Create the reciprocal vectors
	oRecipVec = MrVar_RecipVec(oR1, oR2, oR3, oR4)
	
	;Charge Density
	;   - 1e4 converts to C/cm^2
	oRho = (e0 / q * 1e-12) * oRecipVec -> Divergence( oE1, oE2, oE3, oE4 )
	oRho -> SetName, rho_vname
	oRho -> Cache
	oRho['CATDESC'] = 'Free charge density, computed from the divergence of the electric field.'
	oRho['TITLE']   = '$\Delta$#!C(#/cm^3)'
	oRho['UNITS']   = '#/cm^3'
	
;-------------------------------------------
; Current Density //////////////////////////
;-------------------------------------------
	;Get N and V
	oNe = MrVar_Get(ne_vname)
	oNi = MrVar_Get(ni_vname)
	oVi = MrVar_Get(vi_vname)
	oVe = MrVar_Get(ve_vname)
		
	;Interpolate DIS to DES
	oVi_des = oVi -> Interpol(oVe)
	oNi_des = oNi -> Interpol(oVe)
	
	;Compute current density
	;   - 1e15 converts to uA/m^2
	q  = MrConstants('q')
	oJ = q * 1e15 * oNe * (oVi_des - oVe)
	oJ -> SetName, j_vname
	oJ -> Cache
	oJ['CATDESC'] = 'Current density calculated from particle moments: J=q*ne*(Vi-Ve).'
	oJ['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oJ['LABEL']   = ['X', 'Y', 'Z']
	oJ['TITLE']   = 'J!C$\mu$A/m^2'
	oJ['UNITS']   = 'uA/m^2'
	
;-------------------------------------------
; Electric Fields //////////////////////////
;-------------------------------------------
	oB = MrVar_Get(bvec_vname)
	CASE StrUpCase(sc) OF
		'MMS1': oE = MrVar_Get(e_vnames[0])
		'MMS2': oE = MrVar_Get(e_vnames[1])
		'MMS3': oE = MrVar_Get(e_vnames[2])
		'MMS4': oE = MrVar_Get(e_vnames[3])
		ELSE: Message, 'Invalid spacecraft: "' + sc + '".'
	ENDCASE
	
	;Interpolate to DES
	oB_des = oB -> Interpol(oVe)
	oE_des = oE -> Interpol(oVe)
	
	;Hall field
	;   - 1e-18 converts to mV/m
	oJxB = 1.0/(1e18 * q * oNe) * oJ -> Cross(oB_des)
	oJxB -> SetName, jxb_vname
	oJxB -> Cache
	oJxB['CATDESC'] = 'E = (JxB)/(ne)'
	oJxB['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oJxB['LABEL']   = ['X', 'Y', 'Z']
	oJxB['TITLE']   = ['JxB/ne!C(mV/m)']
	oJxB['UNITS']   = ['mV/m']
	
	;Convective field
	;   - 1e-3 converts to mV/m
	oVexB = -1e-3 * oVe -> Cross(oB_des)
	oVexB -> SetName, vexb_vname
	oVexB -> Cache
	oVexB['CATDESC'] = 'E = -VexB'
	oVexB['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oVexB['LABEL']   = ['X', 'Y', 'Z']
	oVexB['TITLE']   = ['-VexB!C(mV/m)']
	oVexB['UNITS']   = ['mV/m']
	
	;Convective field
	;   - 1e-3 converts to mV/m
	oVixB = -1e-3 * oVi_des -> Cross(oB_des)
	oVixB -> SetName, vixb_vname
	oVixB -> Cache
	oVixB['CATDESC'] = 'E = -VixB'
	oVixB['COLOR']   = ['Blue', 'Forest Green', 'Red']
	oVixB['LABEL']   = ['X', 'Y', 'Z']
	oVixB['TITLE']   = ['-VexB!C(mV/m)']
	oVixB['UNITS']   = ['mV/m']
	
;-------------------------------------------
; Split E-Fields into X, Y, Z Components ///
;-------------------------------------------
	oE_des -> Split, oEx, oEy, oEz, /CACHE, NAME=[ex_vname, ey_vname, ez_vname]
	oJxB   -> Split, oJxBx, oJxBy, oJxBz, /CACHE
	oVexB  -> Split, oVexBx, oVexBy, oVexBz, /CACHE
	oVixB  -> Split, oVixBx, oVixBy, oVixBz, /CACHE
	
	;E
	oEx['LABEL'] = 'E'
	oEx['TITLE'] = 'Ex!C(mV/m)'
	
	oEy['LABEL'] = 'E'
	oEy['TITLE'] = 'Ey!C(mV/m)'
	
	oEz['LABEL'] = 'E'
	oEz['TITLE'] = 'Ez!C(mV/m)'
	
	;JXB
	oJxBx['COLOR'] = 'Blue'
	oJxBx['LABEL'] = 'JxB/ne'
	
	oJxBy['COLOR'] = 'Blue'
	oJxBy['LABEL'] = 'JxB/ne'
	
	oJxBz['COLOR'] = 'Blue'
	oJxBz['LABEL'] = 'JxB/ne'
	
	;VIxB
	oVixBx['COLOR'] = 'Forest Green'
	oVixBx['LABEL'] = '-VixB'
	
	oVixBy['COLOR'] = 'Forest Green'
	oVixBy['LABEL'] = '-VixB'
	
	oVixBz['COLOR'] = 'Forest Green'
	oVixBz['LABEL'] = '-VixB'
	
	;VexB
	oVexBx['COLOR'] = 'Red'
	oVexBx['LABEL'] = '-VexB
	
	oVexBy['COLOR'] = 'Red'
	oVexBy['LABEL'] = '-VexB
	
	oVexBz['COLOR'] = 'Red'
	oVexBz['LABEL'] = '-VexB
	
;-------------------------------------------
; Properties ///////////////////////////////
;-------------------------------------------
	oNi = MrVar_Get(ni_vname)
	oNe = MrVar_Get(ne_vname)
	oNi['AXIS_RANGE'] = [ Min([oNi.min, oNe.min]), Max([oNi.max, oNe.max]) ]
	oNi['COLOR']      = 'Blue'
	oNi['LABEL']      = 'Ni'

	oNe['COLOR'] = 'Red'
	oNe['LABEL'] = 'Ne'


	title = StrUpCase(StrJoin([sc, mode, level, 'Hall'], ' '))
	oEx['PLOT_TITLE'] = title

;-------------------------------------------
; Plot Data ////////////////////////////////
;-------------------------------------------
	;Plot data
	win = MrVar_PlotTS( [ex_vname, ey_vname, ez_vname, ni_vname, bvec_vname, j_vname, rho_vname], $
	                    /NO_REFRESH, $
	                    XSIZE = 600, $
	                    YSIZE = 700 )
	win = MrVar_OPlotTS( ex_vname, [vixbx_vname, vexbx_vname, jxbx_vname] )
	win = MrVar_OPlotTS( ey_vname, [vixby_vname, vexby_vname, jxby_vname] )
	win = MrVar_OPlotTS( ez_vname, [vixbz_vname, vexbz_vname, jxbz_vname] )
	win = MrVar_OPlotTS( ni_vname, ne_vname )
	
	;Pretty-up the window
	win[0] -> SetLayout, [1,1]
	win    -> TrimLayout
	win    -> SetProperty, OXMARGIN=[13, 9]
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
		fname   = StrJoin( [sc, 'edp', mode, level, 'hall'], '_' )
		fname   = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done /////////////////////////////////////
;-------------------------------------------

	RETURN, win
END