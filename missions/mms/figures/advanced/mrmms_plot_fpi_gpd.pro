; docformat = 'rst'
;
; NAME:
;       MrMMS_Plot_FPI_GPD
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
;+
;   Plot gyrophase distributions at different energy levels. Complement with plots
;   of Sqrt(Q) calculated within the same energy levels.
;
;       1. Bxyz, |B|
;       2. GPD 20eV
;       3. GPD 60eV
;       4. GPD 250eV
;       5. GPD 500eV
;       6. GPD 1000eV
;       7. Sqrt(Q)
;
; :Params:
;       SC:         in, required, type=string
;                   MMS spacecraft ID. Options are {'mms1' | 'mms2' | 'mms3' | 'mms4'}
;       MODE:       in, required, type=string
;                   Data rate mode. Options are {'srvy' | 'brst'}
;       SPECIES:    in, required, type=string
;                   Particle species. Options are {'e' | 'i'}
;
; :Keywords:
;       NO_LOAD:    in, optional, type=boolean, default=0
;                   If set, data will not be loaded from source files.
;       OUTPUT_DIR: in, optional, type=string, default=pwd
;                   A directory in which to save the figure. If neither `OUTPUT_DIR`
;                       nor `OUTPUT_EXT` are defined, no file is generated.
;       OUTPUT_EXT: in, optional, type=string, default=pwd
;                   File extensions for the output figure. Options include: 'eps', 'gif',
;                       'jpg', 'ps', 'pdf', 'png', 'tiff'. If neither `OUTPUT_DIR` nor
;                       `OUTPUT_EXT` are defined, no file is generated.
;
; :Categories:
;    MMS
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
;       2017/02/07  -   Written by Matthew Argall
;-
FUNCTION MrMMS_Plot_FPI_GPD, sc, mode, species, $
COORDS=coords, $
ENERGIES=energies, $
FGM_INSTR=fgm_instr, $
LEVEL=level, $
NO_LOAD=no_load, $
OUTPUT_DIR=output_dir, $
OUTPUT_EXT=output_ext, $
TAIL=tail
	Compile_Opt idl2
	
	Catch, the_error
	IF the_error NE 0 THEN BEGIN
		Catch, /CANCEL
		MrPrintF, 'LogErr'
		RETURN, Obj_New()
	ENDIF
	
	tf_load = ~Keyword_Set(no_load)
	tf_tail = Keyword_Set(tail)
	IF N_Elements(coords)    EQ 0 THEN coords    = 'gse'
	IF N_Elements(fgm_instr) EQ 0 THEN fgm_instr = 'fgm'
	IF N_Elements(level)     EQ 0 THEN level     = 'l2'
	IF N_Elements(species)   EQ 0 THEN species   = 'e'
	instr   = 'd' + species + 's'
	
	IF N_Elements(energies) EQ 0 THEN BEGIN
		IF tf_tail THEN BEGIN
			energies = [ [ 0, 10], $     ;< 100eV
			             [11, 14], $     ;250 eV
			             [15, 18], $     ;500 eV
			             [19, 25], $     ;2.5 keV
			             [26, 31] ]      ;> 7 keV
		ENDIF ELSE BEGIN
			energies = [ [ 0,  4], $    ;20 eV
			             [ 5,  9], $    ;60 eV
			             [10, 13], $    ;250 eV
			             [14, 17], $    ;500 eV
			             [18, 31] ]     ;> 1 keV
		ENDELSE
	ENDIF
	nEnergies = N_Elements(energies[0,*])
	

;-------------------------------------------
; Variable Names ///////////////////////////
;-------------------------------------------
	fpi_mode  = mode EQ 'brst' ? mode : 'fast'
	fpi_instr = 'd' + species + 's'

	;Source
	b_vname     = StrJoin( [sc, fgm_instr, 'b',     coords, mode, level], '_' )
	bvec_vname  = StrJoin( [sc, fgm_instr, 'bvec',  coords, mode, level], '_' )
	bmag_vname  = StrJoin( [sc, fgm_instr, 'bmag',  coords, mode, level], '_' )
	f_vname     = StrJoin( [sc, instr, 'dist', fpi_mode], '_')
	v_vname     = StrJoin( [sc, fpi_instr, 'bulkv', coords, fpi_mode], '_')
	p_vname     = StrJoin( [sc, fpi_instr, 'prestensor', coords, fpi_mode], '_')
	scpot_vname = StrJoin( [sc, 'edp', 'scpot', 'fast', level], '_' )
	
	;Derived Moments
	pscl_vname    = StrJoin( [sc, fpi_instr, 'p', coords, fpi_mode], '_')
	pscl_calc_vname = StrJoin( [sc, fpi_instr, 'p',    'calc', coords, fpi_mode], '_')
	gpd_vnames    = StrJoin( [sc, fpi_instr, 'gpd',     fpi_mode], '_') + '_' + StrTrim(IndGen(nEnergies+1), 2)
	q_vnames      = StrJoin( [sc, fpi_instr, 'qfactor', fpi_mode], '_') + '_' + StrTrim(IndGen(nEnergies+2), 2)

;-------------------------------------------
; Load Data ////////////////////////////////
;-------------------------------------------

	IF tf_load THEN BEGIN
		;FGM
		MrMMS_FGM_Load_Data, sc, mode, $
		                     INSTR  = fgm_instr, $
		                     LEVEL  = level, $
		                     VARFORMAT = '*b_gse*', $
		                     SUFFIX = suffix
		
		;FPI-DIST
		MrMMS_FPI_Load_Dist3D, sc, fpi_mode, species, $
		                       /APPLY_MODEL, $
		                       COORD_SYS   = coord_sys, $
		                       LEVEL       = level, $
		                       ORIENTATION = orientation
		
		;FPI-MOMS
		MrMMS_FPI_Load_Data, sc, fpi_mode, $
		                     OPTDESC   = fpi_instr + '-moms', $
		                     LEVEL     = level, $
		                     VARFORMAT = ['*pres*'+coords+'*', '*bulkv*'+coords+'*']
		
		;Spacecraft potential
		MrMMS_Load_Data, sc, 'edp', 'fast', 'l2', $
		                 OPTDESC   = 'scpot', $
		                 VARFORMAT = '*scpot*'
	ENDIF
	
;-------------------------------------------
; Energy Bin Labels ////////////////////////
;-------------------------------------------
	
	;Determine axis labels
	oDist = MrVar_Get(f_vname)
	oE    = oDist['DEPEND_3']
	E_labels = StrArr(5)
	E_titles = StrArr(5)
	FOR i = 0, nEnergies-1 DO BEGIN
		E0 = oE['DATA',0,energies[0,i]]
		E1 = oE['DATA',0,energies[1,i]]
		E0 = E0 LT 1000 ? String(E0, FORMAT='(i0)')+'eV' : String(E0/1000.0, FORMAT='(i0)')+'keV'
		E1 = E1 LT 1000 ? String(E1, FORMAT='(i0)')+'eV' : String(E1/1000.0, FORMAT='(i0)')+'keV'
		E_labels[i] = E0 + '-' + E1
		E_titles[i] = E0 + '!C' + E1
	ENDFOR
	
;-------------------------------------------
; GPD & Q //////////////////////////////////
;-------------------------------------------

	;Size of each velocity-space volume element
	theSpecies = species EQ 'i' ? 'H' : species
	oDist4D    = MrDist4D( f_vname, VSC=scpot_vname, SPECIES=theSpecies )
	odV        = oDist4D -> VolumeElement()

	;Rotate the distribution function into field-aligned coordinates
	;   - Because the angles are no longer distributed evenly across the
	;     grid in FAC coordinates, the internal deltas are not correct
	oDist_FAC = MrVar_Dist_Rotate( f_vname, bvec_vname, v_vname, 'VXB' )
	
	;Distribution function
	theSpecies  = species EQ 'i' ? 'H' : species
	oDist4D_FAC = MrDist4D_v1(oDist_FAC, VSC=scpot_vname, SPECIES=theSpecies)
	
	;Gyrophase distribution
	;   - Weight each bin by its original size.
	;   - This compensates for the non-uniform spacing in FAC coordinates
	;   - Do not have to re-bin data into FAC coordinate grid
	oGPD      = ObjArr(nEnergies)
	oPres     = ObjArr(nEnergies+1)
	oQ        = ObjArr(nEnergies+2)
	oDist4D   = MrDist4D(f_vname, VSC=scpot_vname, SPECIES=theSpecies)
	qrange    = [0, -!Values.f_infinity]
	colors = MrDefaultColor(NCOLORS=7)
	
	FOR i = 0, nEnergies - 1 DO BEGIN
		;Gyrophase distributions at each energy
		oGPD[i]  = oDist4D_FAC -> PhiSpec( E_RANGE=energies[*,i], /CACHE, NAME=gpd_vnames[i],  THETA_RANGE=[0.0, 180.0], UNITS='EFLUX', WEIGHT=odV )
		oPhi     = (oGPD[i])['DEPEND_1']
		oPhi['AXIS_RANGE']   = [-180.0, 180.0]
		oPhi['TICKINTERVAL'] = 90.0
		oPhi['TITLE']        = 'Gyrophase!C' + E_titles[0] + '!C(deg)'
		
		;Q-Factor at each energy
		oDist4D -> Moments, PRESSURE = oP, ENERGY_RANGE = energies[*,i]
		oPres[i] = oP
		oQ[i]    = MrVar_Pres_QFactor( bvec_vname, oP, /CACHE, NAME=q_vnames[i] )
		qrange[1]        >= oQ[i].max
		(oQ[i])['COLOR'] = colors[i]
		(oQ[i])['LABEL'] = E_labels[i]
	ENDFOR
	
	;Q-Factor for all energies: DES & Derived
	oDist4D -> Moments, PRESSURE = oP, ENERGY_RANGE = [0, 31]
	oPres[-1] = oP
	
	oQ[-2] = MrVar_Pres_QFactor( bvec_vname, oP, /CACHE, NAME=q_vnames[-2] )
	oQ[-1] = MrVar_Pres_QFactor( bvec_vname, p_vname, /CACHE, NAME=q_vnames[-1] )
	
	Obj_Destroy, [odV, oDist_FAC, oDist4D_FAC, oDist4D]

;-------------------------------------------
; Split Moments ////////////////////////////
;-------------------------------------------
	oP    = MrVar_Get(p_vname)
	oPscl = MrScalarTS( oP['TIMEVAR'], (oP[*,0,0] + oP[*,1,1] + oP[*,2,2]) / 3.0, $
	                    /CACHE, $
	                    NAME = pscl_vname )
	
	oP    = oPres[-1]
	oPscl = MrScalarTS( oP['TIMEVAR'], (oP[*,0,0] + oP[*,1,1] + oP[*,2,2]) / 3.0, $
	                    /CACHE, $
	                    NAME = pscl_calc_vname )

;-------------------------------------------
; Attributes ///////////////////////////////
;-------------------------------------------
	;B
	oB = MrVar_Get(b_vname)
	oB['PLOT_TITLE'] = StrUpCase( StrJoin( [sc, fpi_instr, fpi_mode, level], ' ' ) )
	
	;
	; GPD
	;
	
	
	;
	; Q
	;
	(oQ[0])['AXIS_RANGE'] = qrange
	(oQ[-2])['LABEL']     = 'UNH'
	(oQ[-1])['LABEL']     = 'FPI'
	
	;
	; P
	;
	oP = MrVar_Get(pscl_calc_vname)
	oP['COLOR'] = 'Blue'
	oP['LABEL'] = 'UNH'
	oP['TITLE'] = 'P!C(nPa)'
	
	oP = MrVar_Get(pscl_vname)
	oP['COLOR'] = 'Black'
	oP['LABEL'] = 'FPI'
	oP['TITLE'] = 'P!C(nPa)'

;-------------------------------------------
; Plot /////////////////////////////////////
;-------------------------------------------
	win = MrVar_PlotTS( [b_vname, gpd_vnames, q_vnames[0], pscl_vname], $
	                    /NO_REFRESH, $
	                    YSIZE = 750 )
	
	;Overplot Q-factors
	win = MrVar_OPlotTS( q_vnames[0], q_vnames[1:*] )
	win = MrVar_OPlotTS( pscl_vname, pscl_calc_vname )

	win[0] -> SetLayout, [1,1]
	win -> SetProperty, OXMARGIN=[15,14]
	win -> TrimLayout
	win -> Refresh

;-------------------------------------------
; Save Figure //////////////////////////////
;-------------------------------------------
	IF N_Elements(output_dir) GT 0 || N_Elements(output_ext) GT 0 THEN BEGIN
		;Defaults
		IF N_Elements(output_dir) EQ 0 THEN BEGIN
			CD, CURRENT=output_dir
			MrPrintF, 'LogText', 'Saving file to: "' + output_dir + '".'
		ENDIF
		
		;File name
		fname = StrJoin( [sc, fpi_instr, fpi_mode, level, 'gpd'], '_' )
		fname = FilePath( fname, ROOT_DIR=output_dir )
		
		;Save the figure
		fout = MrVar_PlotTS_Save( win, fname, output_ext )
	ENDIF

;-------------------------------------------
; Done! ////////////////////////////////////
;-------------------------------------------
	RETURN, win
END