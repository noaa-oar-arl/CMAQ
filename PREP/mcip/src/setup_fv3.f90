!------------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in           !
!  continuous development by various groups and is based on information        !
!  from these groups: Federal Government employees, contractors working        !
!  within a United States Government contract, and non-Federal sources         !
!  including research institutions.  These groups give the Government          !
!  permission to use, prepare derivative works of, and distribute copies       !
!  of their work in the CMAQ system to the public and to permit others         !
!  to do so.  The United States Environmental Protection Agency                !
!  therefore grants similar permission to use the CMAQ system software,        !
!  but users are requested to provide copies of derivative works or            !
!  products designed to operate in the CMAQ system to the United States        !
!  Government without restrictions as to use by others.  Software              !
!  that is used with the CMAQ system but distributed under the GNU             !
!  General Public License or the GNU Lesser General Public License is          !
!  subject to their copyright restrictions.                                    !
!------------------------------------------------------------------------------!

SUBROUTINE setup_fv3 (cdfid, cdfid2, ctmlays)

!-------------------------------------------------------------------------------
! Name:     Set Up the FV3 Domain Attributes
! Purpose:  Establishes bounds for FV3 post-processing.
! Revised:  ?? Jun 2004  Modified from MCIP2.2 for FV3. (S.-B. Kim)
!           26 May 2005  Changed vertical dimension to reflect full-layer
!                        dimension in FV3 header.  Added dynamic calculation
!                        of MET_TAPFRQ.  Converted dimensions to X,Y as opposed
!                        to the (former) convention that aligned with MM5.
!                        Included updates from MCIPv2.3.  Added calculation of
!                        cone factor.  Added logic for moist species, 2-m
!                        temperature, and 10-m winds.  Added definitions for
!                        FV3 base state variables.  Added capability to use all
!                        FV3 layers for MCIP without defining a priori.
!                        Cleaned up code.  (T. Otte)
!           15 Jul 2005  Added debugging on variable retrievals.  Changed check
!                        on 3D mixing ratios from rain to ice.  Corrected RADM
!                        seasons for Southern Hemisphere.  Corrected variable
!                        name for retrieval of surface physics option. (T. Otte)
!           18 Aug 2005  Changed internal variable SIGN to FAC to avoid
!                        confusion with F90 intrinsic function.  (T. Otte)
!           10 Apr 2006  Corrected checking of I/O API variables for Mercator
!                        projection.  (T. Otte)
!           12 May 2006  Corrected setting of I/O API variables for polar
!                        stereographic projection.  Revised defining and
!                        setting projection variables for module METINFO.
!                        Added restriction on using Eta/Ferrier microphysics
!                        scheme where QCLOUD represents total condensate.
!                        (T. Otte)
!           20 Jun 2006  Changed setting of IDTSEC from REAL to INTEGER
!                        value.  (T. Otte)
!           27 Jul 2007  Removed settings for RADMdry variable ISESN and for
!                        MET_INHYD.  Updated read of P_TOP to account for new
!                        method of storing "real" scalars in FV3 I/O API with
!                        FV3.  Added checks for fractional land use, leaf
!                        area index, Monin-Obukhov length, aerodynamic and
!                        stomatal resistances, vegetation fraction, canopy
!                        wetness, and soil moisture, temperature, and type in
!                        FV3 file.  Added read for number of land use
!                        categories...new with FV3V2.2.  Added read for number
!                        of soil layers, MET_RELEASE, MET_FDDA_3DAN and
!                        MET_FDDA_OBS.  Set MET_FDDA_SFAN to 0 for now because
!                        that option is not in FV3 ARW as of V2.2.  Changed
!                        MET_RADIATION into MET_LW_RAD and MET_SW_RAD.
!                        (T. Otte)
!           06 May 2008  Changed criteria for setting NUMMETLU when netCDF
!                        dimension "land_cat_stag" does not exist.  Added
!                        checks to determine if 2-m mixing ratio (Q2) and
!                        turbulent kinetic energy (TKE_MYJ) arrays exist, and
!                        set flags appropriately.  Extract nudging coefficients
!                        from header to use in metadata.  Extract whether or
!                        not the urban canopy model was used.  (T. Otte)
!           27 Oct 2009  Cleaned up file opening and logging in FV3 I/O API to
!                        prevent condition with too many files open for long
!                        simulations.  Added MODIFIED IGBP MODIS NOAH and 
!                        NLCD/MODIS as land-use classification options.
!                        Changed MET_UCMCALL to MET_URBAN_PHYS, and allowed
!                        for variable to be set to be greater than 1.  Chnaged
!                        code to allow for surface analysis nudging option
!                        and coefficients to be defined per FV3.  Define
!                        MET_CEN_LAT, MET_CEN_LON, MET_RICTR_DOT, MET_RJCTR_DOT,
!                        and MET_REF_LAT.  Increased MAX_TIMES to 1000.  Compute
!                        MET_XXCTR and MET_YYCTR.  Corrected setting for
!                        DATE_INIT, and fill variable MET_RESTART.  Read number
!                        of land use categories from FV3 global attributes for
!                        FV3V3.1 and beyond.  Allow output from FV3
!                        Preprocessing System (WPS) routine, GEOGRID, to provide
!                        fractional land use output if it is unavailable in FV3
!                        output.  Fill MET_P_ALP_D and MET_P_BET_D here
!                        rather than in setgriddefs.F for Mercator.  Added
!                        new logical variables IFLUFV3OUT and IFZNT.  (T. Otte)
!           12 Feb 2010  Removed unused variables COMM and SYSDEP_INFO.
!                        (T. Otte)
!           18 Mar 2010  Added CDFID as an input argument, and no longer open
!                        and close FV3 history file here.  Added CDFIDG as an
!                        input argument for subroutine CHKWPSHDR.  (T. Otte)
!           15 Dec 2010  Improved support for long MCIP runs from long FV3
!                        runs by increasing MAX_TIMES to 9999.  Added
!                        MET_RAIN_BUCKET.  (T. Otte)
!           23 Feb 2011  Refined error checking for MET_RAIN_BUCKET.  (T. Otte)
!           11 Aug 2011  Added MET_SHAL_CU to input.  Replaced module PARMS3
!                        with I/O API module M3UTILIO.  (T. Otte)
!           24 Aug 2011  Changed name of module FILE to FILES to avoid conflict
!                        with F90 protected intrinsic.  Updated netCDF commands
!                        to F90, and improved error handling.  Replaced calls
!                        to GET_TIMES_CDF with explicit netCDF functions.
!                        (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           21 Nov 2011  Force 2-m water vapor mixing ratio from FV3 with
!                        YSU PBL to be filled with layer 1 QVAPOR to avoid
!                        occasional Q2 < 0 in wintertime.  (T. Otte)
!           07 Dec 2011  Removed requirement to fill nudging coefficient for
!                        moisture when spectral nudging is used in FV3; as of
!                        FV3, spectral nudging toward moisture is not
!                        released in FV3.  Also added provision to collect
!                        nudging coefficient for geopotential when spectral
!                        nudging is used; was added to FV3 header with FV3v3.2.
!                        (T. Otte)
!           21 Aug 2012  Added MET_PCP_INCR for FV3V3.2+.  (T. Otte)
!           10 Sep 2012  Added handling for 40-category 2006 NLCD-MODIS land
!                        use classification as "NLCD40".  Added alternate name
!                        for 50-category 2001 NLCD-MODIS land use classification
!                        as "NLCD50".  (T. Otte)
!           26 Nov 2014  Added reads of ice, lake, and urban land use indices,
!                        and moved those definitions from getluse.f90 to this
!                        routine.  (T. Spero)
!           10 Apr 2015  Determine if 3D cloud fraction is part of FV3 output
!                        and if it represents resolved clouds.  Fill new logical
!                        variable IFCLD3D appropriately so that if resolved
!                        cloud fraction is available, it will be passed through
!                        in output.  (T. Spero)
!           21 Aug 2015  Added flag to capture whether ACM2 was run so that
!                        Monin-Obukhov length can be recalculated following
!                        the "corrector" part of the predictor-corrector in
!                        FV3/ACM2.  (T. Spero)
!           17 Sep 2015  Changed IFMOLACM to IFMOLPX.  (T. Spero)
!           21 Apr 2017  Added MODIS category 21 as "Lake".  (T. Spero)
!           23 Jun 2017  Added a check for FV3's hybrid vertical coordinate
!                        in FV3v3.9 and beyond.  Currently disabled MCIP when
!                        that coordinate is detected.  To be implemented in
!                        a later release of MCIP.  (T. Spero)
!           09 Feb 2018  Added support for hybrid vertical coordinate in FV3
!                        output.  Added capability to read and process data
!                        from the NOAH Mosaic land-surface model.  (T. Spero)
!           26 Jun 2018  Changed name of module with netCDF IO to broaden its
!                        usage.  Now use netCDF tokens for missing data.
!                        (T. Spero)
!           14 Sep 2018  Removed support for MM5v3 input.  (T. Spero)
!           23 Nov 2018  Modify criteria to determine whether incremental
!                        precipitation is available in FV3 output.  FV3v4.0
!                        allows header variable PREC_ACC_DT to appear even if
!                        the accompanying precipitation fields are not in the
!                        output.  (T. Spero)
!           14 Dec 2018  Added flag (IFRCURB) to determine if fraction of urban
!                        area is obtained from urban canopy model.  (T. Spero)
!           10 May 2019  Removed layer collapsing.  (T. Spero)
!           18 Jun 2019  Added a flag (IFPXFV341) to determine of new surface
!                        variables with PX LSM are available to improve dust
!                        simulation in CCTM.  Added a flag (IFCURADFDBK) to
!                        indicate if the convective scheme included radiative
!                        feedbacks.  Added a flag (IFKFRADEXTRAS) for extra
!                        variables available with KF convective scheme with
!                        radiative feedbacks.  (T. Spero)
!           18 Nov 2019  Modified for FV3GFS Capability. (P. C. Campbell)
!-------------------------------------------------------------------------------

  USE metinfo
  USE date_pack
  USE mcipparm
  USE files
  USE netcdf_io
  USE const, ONLY: pi180
  USE netcdf

  IMPLICIT NONE

  INTEGER     ,       INTENT(IN)    :: cdfid, cdfid2
!  INTEGER                           :: cdfid2
  INTEGER                           :: cdfidg
  REAL,               INTENT(OUT)   :: ctmlays     ( maxlays )
  REAL                              :: phalf_lays  ( maxlays )
  REAL                              :: pfull_lays  ( maxlays+1 )
  CHARACTER(LEN=32)                 :: date_init
  INTEGER                           :: dimid
  INTEGER                           :: dimids     ( nf90_max_var_dims )
  REAL,               ALLOCATABLE   :: dum1d      ( : )
  REAL,               ALLOCATABLE   :: dum2d      ( : , : )
  INTEGER,            ALLOCATABLE   :: dum2d_i    ( : , : )
  INTEGER                           :: dx
  INTEGER                           :: dy
  REAL                              :: fac
  CHARACTER(LEN=256)                :: fl
  CHARACTER(LEN=256)                :: fl2
  CHARACTER(LEN=256)                :: flg
  CHARACTER(LEN=256)                :: geofile
  INTEGER                           :: icloud_cu
  INTEGER                           :: id_data
  INTEGER                           :: idtsec
  LOGICAL                           :: ifgeo
  LOGICAL                           :: ifisltyp
  LOGICAL                           :: ifra
  LOGICAL                           :: ifrs
  LOGICAL                           :: ifsmois
  LOGICAL                           :: iftslb
  LOGICAL                           :: ifu10m
  LOGICAL                           :: ifv10m
  INTEGER                           :: it
  INTEGER                           :: ival
  INTEGER                           :: lent
  INTEGER                           :: n_times
  INTEGER                           :: nxm
  INTEGER                           :: nym
  CHARACTER(LEN=16),  PARAMETER     :: pname      = 'SETUP_FV3'
  INTEGER                           :: rcode, rcode2
  REAL                              :: rval
  REAL,  ALLOCATABLE                :: times      ( : )
  INTEGER                           :: varid
  CHARACTER(LEN=80)                 :: fv3version
!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f6000 = "(/, 1x, &
    & '- SUBROUTINE SETUP_FV3 - READING FV3 HEADER')"
  CHARACTER(LEN=256), PARAMETER :: f6100 = "(3x, &
    & 'FV3 GRID DIMENSIONS (X,Y,Z) ', i4, 1x, i4, 1x, i3, //)"

  CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   MISMATCH IN DX AND DY', &
    & /, 1x, '***   DX, DY = ', 2(f7.2), &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9100 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   UNKNOWN LAND USE CLASSIFICATION', &
    & /, 1x, '***   FIRST THREE LETTERS = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9225 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   QCLOUD NOT FOUND IN FV3 OUTPUT...STOPPING', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9250 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ETA/FERRIER SCHEME IS NOT SUPPORTED IN CMAQ', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9275 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   FOUND QCLOUD BUT NOT QRAIN...STOPPING', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9300 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   NQSPECIES SET AT 3',&
    & /, 1x, '***   MCIP NEEDS TO BE MODIFIED FOR THIS CASE', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9410 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING NCF ID FROM FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9420 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR INQUIRING ABOUT VAR IN FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9430 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING DIMS FROM FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9500 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ONLY FOUND ONE FILE WITH ONE TIME PERIOD', &
    & /, 1x, '***   SETTING OUTPUT FREQUENCY TO 1 MINUTE', &
    & /, 1x, 70('*'))" 

  CHARACTER(LEN=256), PARAMETER :: f9550 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   NEED PRECIPITATION ACCUMULATION IN FV3 TO MATCH', &
    & /, 1x, '***   MCIP OUTPUT INTERVAL', &
    & /, 1x, '***   PREC_ACC_DT from FV3: ', i4, &
    & /, 1x, '***   INTVL from MCIP: ', i4, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9600 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR OPENING FV3 NETCDF FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9700 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CLOSING FV3 NETCDF FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9800 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   DID NOT FIND FRACTIONAL LAND USE IN FV3 output', &
    & /, 1x, '***   AND DID NOT FIND GEOGRID FILE' &
    & /, 1x, '***   -- WILL NOT USE FRACTIONAL LAND USE DATA' &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Extract NX, NY, and NZ.
!-------------------------------------------------------------------------------
  WRITE (*,f6000)

  fl = file_mm(1)

  rcode = nf90_get_att (cdfid, nf90_global, 'im', met_nx)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'WEST-EAST_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_get_att (cdfid, nf90_global, 'jm',  &
                        met_ny)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'SOUTH-NORTH_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inq_dimid (cdfid, 'phalf', dimid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'ID for phalf',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
  rcode = nf90_inquire_dimension (cdfid, dimid, len=ival)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'BOTTOM-TOP_GRID_DIMENSION+1',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ELSE

! Set met_nz
    IF ( needlayers ) THEN
     met_nz = MIN(maxlays,ival-1) !If ival > max layers, cap at subset of IOAPI max layers
    ELSE
     met_nz = SIZE(ctmlays)  
    ENDIF

  ENDIF


  WRITE (*,f6100) met_nx, met_ny, met_nz

  met_rictr_dot = FLOAT(met_nx - 1) / 2.0 + 1.0
  met_rjctr_dot = FLOAT(met_ny - 1) / 2.0 + 1.0

!-------------------------------------------------------------------------------
! Read FV3 Pressure layers.
!-------------------------------------------------------------------------------
! Set NLAYS
  nlays = met_nz
  
  rcode = nf90_inq_varid (cdfid, 'phalf', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'phalf',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ELSE
  ALLOCATE ( dum1d ( ival ) ) 
  rcode = nf90_get_var (cdfid, varid, dum1d)

  pfull_lays(1:met_nz+1) = dum1d(ival:ival-met_nz:-1) 

  DEALLOCATE (dum1d)
  ENDIF
  
  rcode = nf90_inq_varid (cdfid, 'pfull', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'pfull',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ELSE
  ALLOCATE ( dum1d ( ival-1 ) )
  rcode = nf90_get_var (cdfid, varid, dum1d)

  phalf_lays(1:met_nz) = dum1d(ival-1:ival-met_nz-1:-1)

  DEALLOCATE (dum1d)
  ENDIF

!  CALL get_var_1d_real_cdf (cdfid, 'pfull', phalf_lays_read, 1, rcode)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'pfull', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF

  
!  CALL get_var_1d_real_cdf (cdfid, 'phalf', pfull_lays(nlays+1:1:-1), 1, rcode)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'phalf', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF




!-------------------------------------------------------------------------------
! Extract domain attributes.
!-------------------------------------------------------------------------------

   rcode = nf90_get_att (cdfid, nf90_global, 'source', fv3version)
   IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'source', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

!  rcode = nf90_get_att (cdfid, nf90_global, 'DX', dx)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'DX', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!
!  rcode = nf90_get_att (cdfid, nf90_global, 'DY', dy)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'DY', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF

!FV3 the user sets the input dx/dy resolutions in namelist
   dx=dx_in
   dy=dy_in
!
  IF (dx == dy) THEN
    met_resoln = dx
  ELSE
    WRITE (*,f9000) TRIM(pname), dx, dy
    CALL graceful_stop (pname)
  ENDIF

  met_nxcoarse = met_nx 
  met_nycoarse = met_ny
  met_gratio   = 1
  met_x_11     = 1
  met_y_11     = 1

!  rcode = nf90_get_att (cdfid, nf90_global, 'grid', met_mapproj)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'MAP_PROJ', TRIM(nf90_strerror(rcode))
!   CALL graceful_stop (pname)
!  ENDIF

!  rcode = nf90_get_att (cdfid, nf90_global, 'STAND_LON', met_proj_clon)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'STAND_LON', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!
!  rcode = nf90_get_att (cdfid, nf90_global, 'MOAD_CEN_LAT', met_proj_clat)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'MOAD_CEN_LAT', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!
!  rcode = nf90_get_att (cdfid, nf90_global, 'CEN_LON', met_cen_lon)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'CEN_LON', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!  met_x_centd = met_cen_lon
!
!  rcode = nf90_get_att (cdfid, nf90_global, 'CEN_LAT', met_cen_lat)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'CEN_LAT', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!  met_y_centd = met_cen_lat
!
!  rcode = nf90_get_att (cdfid, nf90_global, 'TRUELAT1', met_tru1)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'TRUELAT1', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!
!  rcode = nf90_get_att (cdfid, nf90_global, 'TRUELAT2', met_tru2)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'TRUELAT2', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF

!  SELECT CASE ( met_mapproj )

!    CASE (1)    
!      met_p_alp_d  = MIN(met_tru1, met_tru2)  ! true latitude 1  [degrees]
!      met_p_bet_d  = MAX(met_tru1, met_tru2)  ! true latitude 2  [degrees]
!      met_p_gam_d  = met_proj_clon            ! central meridian [degrees]
!      IF ( met_proj_clat < 0.0 ) THEN
!        fac = -1.0  ! Southern Hemisphere
!      ELSE
!        fac =  1.0  ! Northern Hemisphere
!      ENDIF
!      IF ( ABS(met_tru1 - met_tru2) > 1.0e-1 ) THEN
!        met_cone_fac = ALOG10(COS(met_tru1 * pi180)) -  &
!                       ALOG10(COS(met_tru2 * pi180))
!        met_cone_fac = met_cone_fac /                                      &
!                       ( ALOG10(TAN((45.0 - fac*met_tru1/2.0) * pi180)) -  &
!                         ALOG10(TAN((45.0 - fac*met_tru2/2.0) * pi180)) )
!      ELSE
!        met_cone_fac = fac * SIN(met_tru1*pi180)
!      ENDIF
!
!      IF ( wrf_lc_ref_lat > -999.0 ) THEN
!        met_ref_lat = wrf_lc_ref_lat
!      ELSE
!        met_ref_lat = ( met_tru1 + met_tru2 ) * 0.5
!      ENDIF
!
!      CALL ll2xy_lam (met_cen_lat, met_cen_lon, met_tru1, met_tru2,  &
!                      met_proj_clon, met_ref_lat, met_xxctr, met_yyctr)
    
!    CASE (2)  ! polar stereographic
!      met_p_alp_d  = SIGN(1.0, met_y_centd)   ! +/-1.0 for North/South Pole
!      met_p_bet_d  = met_tru1                 ! true latitude    [degrees]
!      met_p_gam_d  = met_proj_clon            ! central meridian [degrees]
!      met_cone_fac = 1.0                      ! cone factor
!      met_ref_lat  = -999.0                   ! not used
!
!      CALL ll2xy_ps (met_cen_lat, met_cen_lon, met_tru1, met_proj_clon,  &
!                     met_xxctr, met_yyctr)
!    
!    CASE (3)  ! Mercator
!      met_p_alp_d  = 0.0                      ! lat of coord origin [deg]
!      met_p_bet_d  = 0.0                      ! (not used)
!      met_p_gam_d  = met_proj_clon            ! lon of coord origin [deg]
!      met_cone_fac = 0.0                      ! cone factor
!      met_ref_lat  = -999.0                   ! not used
!
!      CALL ll2xy_merc (met_cen_lat, met_cen_lon, met_proj_clon,  &
!                       met_xxctr, met_yyctr)
!    
!    CASE DEFAULT
!      met_p_bet_d  = fillreal                 ! missing
!      met_p_alp_d  = fillreal                 ! missing
!      met_p_gam_d  = fillreal                 ! missing
!      met_cone_fac = fillreal                 ! missing
!      met_ref_lat  = fillreal                 ! missing
!  
!  END SELECT

!     FV3 Gaussian Global Grid
      met_mapproj    = 4                      ! FV3 Gaussian Map Projection      
      met_proj_clon  = 0.0
      met_p_alp_d  = 0.0                      ! lat of coord origin [deg]
      met_p_bet_d  = 0.0                      ! (not used)
      met_p_gam_d  = met_proj_clon            ! lon of coord origin [deg]
      met_cone_fac = 0.0                      ! cone factor
      met_ref_lat  = -999.0                   ! not used
      
       met_cen_lat  = met_cen_lat_in          ! set from namelist
       met_cen_lon  = met_cen_lon_in          ! set from namelist                  
!      met_cen_lat  = 89.9103191600434         ! Hardcoded for now GFSv16 Gaussian...
!      met_cen_lon  = 359.882792740194         ! Hardcoded for now GFSv16 Gaussian...
!      met_cen_lat  = 45.0         ! Hardcoded for now GFSv16 Gaussian...
!      met_cen_lon  = -180.0         ! Hardcoded for now GFSv16 Gaussian...

      CALL ll2xy_gau (met_cen_lat, met_cen_lon, met_proj_clon,  &
                       met_xxctr, met_yyctr)

!-------------------------------------------------------------------------------
! Extract model run options.
!-------------------------------------------------------------------------------

!  rcode = nf90_get_att (cdfid, nf90_global, 'MMINLU', met_lu_src)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'MMINLU', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!
!  rcode = nf90_get_att (cdfid, nf90_global, 'ISWATER', met_lu_water)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'ISWATER', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!
!  rcode = nf90_inq_dimid (cdfid, 'soil_layers_stag', dimid)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'ID for soil_layers_stag',  &
!                    TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!  rcode = nf90_inquire_dimension (cdfid, dimid, len=met_ns)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'value for soil_layers_stag',  &
!                    TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF

!   rcode = nf90_get_att (cdfid, nf90_global, 'nsoil', met_ns)
!   IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'value for soil layers', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!   ENDIF

    met_ns=4  ! ytang

  ! Determine if NOAH Mosaic was run and created the correct output fields.
  ! Note that this code is temporarily modified later in this subroutine to
  ! toggle IFMOSAIC to FALSE if the fractional land use arrays are also not
  ! available.  That change is temporary until a subroutine is added to
  ! reconstruct the fractional land use rank if it is missing.

  ! FV3 does not have a Noah Mosaic option, so turn off
  ifmosaic = .FALSE.

!  rcode = nf90_inq_dimid (cdfid, 'mosaic', dimid)
!  IF ( rcode /= nf90_noerr ) THEN
!    ifmosaic = .FALSE.
!  ELSE
!    rcode = nf90_inq_varid (cdfid, 'TSK_MOSAIC', rcode)
!    IF ( rcode /= nf90_noerr ) THEN
!      ifmosaic = .FALSE.
!    ELSE
!      rcode = nf90_inq_varid (cdfid, 'ZNT_MOSAIC', rcode)
!      IF ( rcode /= nf90_noerr ) THEN
!        ifmosaic = .FALSE.
!      ELSE
!        rcode = nf90_inq_varid (cdfid, 'LAI_MOSAIC', rcode)
!        IF ( rcode /= nf90_noerr ) THEN
!          ifmosaic = .FALSE.
!        ELSE
!          rcode = nf90_inq_varid (cdfid, 'RS_MOSAIC', rcode)
!          IF ( rcode /= nf90_noerr ) THEN
!            ifmosaic = .FALSE.
!          ELSE
!            ifmosaic = .TRUE.
!            rcode = nf90_inquire_dimension (cdfid, dimid, len=nummosaic)
!            IF ( rcode /= nf90_noerr ) THEN
!              WRITE (*,f9400) TRIM(pname), 'value for mosaic',  &
!                              TRIM(nf90_strerror(rcode))
!              CALL graceful_stop (pname)
!            ENDIF
!          ENDIF
!        ENDIF
!      ENDIF
!    ENDIF
!  ENDIF


  ! Define number of land use categories.

  rcode = nf90_inq_varid (cdfid2, 'vtype', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'Land Use Type',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
 
  ELSE
  nxm = met_nx - 1
  nym = met_ny - 1
  ALLOCATE ( dum2d_i ( nxm, nym ) )
  rcode = nf90_get_var (cdfid2, varid, dum2d_i)
  nummetlu=MAXVAL(dum2d_i)
  DEALLOCATE (dum2d_i)
  ENDIF

!  IF ( fv3version(18:22) >= "V3.1" ) THEN  ! FV3v3.1 or later
!
!    rcode = nf90_get_att (cdfid, nf90_global, 'NUM_LAND_CAT', nummetlu)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'NUM_LAND_CAT', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!
!    rcode = nf90_get_att (cdfid, nf90_global, 'ISICE', met_lu_ice)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'ISICE', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!
!    rcode = nf90_get_att (cdfid, nf90_global, 'ISLAKE', met_lu_lake)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'ISLAKE', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!
!    rcode = nf90_get_att (cdfid, nf90_global, 'ISURBAN', met_lu_urban)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'ISURBAN', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!
!  ELSE
!    rcode = nf90_inq_dimid (cdfid, 'land_cat_stag', dimid)
!    IF ( rcode /= nf90_noerr ) THEN  ! only exists with fractional land use
!      SELECT CASE ( met_lu_src(1:3) )
!        CASE ( "USG" )  ! USGS -- typically 24, but can be up to 33 in V2.2+
!          IF ( ( fv3version(18:21) == "V2.2" ) .OR.  &
!               ( fv3version(18:19) == "V3"   ) ) THEN
!            nummetlu = 33
!          ELSE
!            nummetlu = 24
!          ENDIF
!          met_lu_water = 16
!          met_lu_ice   = 24
!          met_lu_urban =  1
!          met_lu_lake  = -1
!        CASE ( "OLD" )  ! old MM5 13-category system
!          nummetlu     = 13
!          met_lu_water =  7
!          met_lu_ice   = 11
!          met_lu_urban =  1
!          met_lu_lake  = -1
!        CASE ( "SiB" )  ! SiB 16-category system
!          nummetlu     = 16
!          met_lu_water = 15
!          met_lu_ice   = 16
!          met_lu_urban = -1
!          met_lu_lake  = -1
!        CASE ( "MOD" )  ! Modified IGBP MODIS NOAH 33-category system
!          nummetlu     = 33

! FV3 only  has MODIS 20-category ('MOD'):  Harcoded Ice, water, urban, lake         
          met_lu_src = 'MOD'
          met_lu_water = 17
          met_lu_ice   = 15
          met_lu_urban = 13
          met_lu_lake = -1

!          IF ( fv3version(18:22) >= "V3.8" ) THEN  ! FV3v3.8 or later
!            met_lu_lake = 21
!          ELSE
!            met_lu_lake = -1
!          ENDIF
!        CASE ( "NLC" )  ! NLCD/MODIS combined system
!          IF ( met_lu_src(4:6) == "D40") THEN
!            nummetlu     = 40
!            met_lu_water = 17
!            met_lu_ice   = 15
!            met_lu_urban = 13
!            met_lu_lake  = -1
!          ELSE
!            nummetlu     = 50
!            met_lu_water =  1
!            met_lu_ice   =  2
!            met_lu_urban =  3
!            met_lu_lake  = -1
!          ENDIF
!        CASE DEFAULT
!          WRITE (*,f9100) TRIM(pname), met_lu_src(1:3)
!          CALL graceful_stop (pname)
!      END SELECT
!    ELSE
!      rcode = nf90_inquire_dimension (cdfid, dimid, len=nummetlu)
!      IF ( rcode /= nf90_noerr ) THEN
!        WRITE (*,f9400) TRIM(pname), 'value for land_cat_stag',  &
!                        TRIM(nf90_strerror(rcode))
!        CALL graceful_stop (pname)
!      ENDIF
!    ENDIF
!  ENDIF

!  Initial FV3 version no options/attributes set for different physics options
!  in the met netCDF output file.
!  Initially set all physics = 0, but should match FV3 with WRF settings

!  rcode = nf90_get_att (cdfid, nf90_global, 'RA_LW_PHYSICS', met_lw_rad)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'RA_LW_PHYSICS', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
   met_lw_rad=0
!  rcode = nf90_get_att (cdfid, nf90_global, 'RA_SW_PHYSICS', met_sw_rad)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'RA_SW_PHYSICS', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
   met_sw_rad=0
!  rcode = nf90_get_att (cdfid, nf90_global, 'CU_PHYSICS', met_cumulus)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'CU_PHYSICS', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
   met_cumulus=0   
!  rcode = nf90_get_att (cdfid, nf90_global, 'MP_PHYSICS', met_expl_moist)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'MP_PHYSICS', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
   met_expl_moist=0
!  rcode = nf90_get_att (cdfid, nf90_global, 'BL_PBL_PHYSICS', met_pbl)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'BL_PBL_PHYSICS', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
   met_pbl=0
!  rcode = nf90_get_att (cdfid, nf90_global, 'SF_SFCLAY_PHYSICS', met_sfc_lay)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'SF_SFCLAY_PHYSICS', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
   met_sfc_lay=0
!  rcode = nf90_get_att (cdfid, nf90_global, 'SF_SURFACE_PHYSICS', met_soil_lsm)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'SF_SURFACE_PHYSICS',  &
!                   TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
   met_soil_lsm=0
!
!  ! Determine if an urban model was used.
!  
!  ! FV3 does not have an urban model.
   met_urban_phys = 0
   ifrcurb = .FALSE.
!  IF ( fv3version(18:21) >= "V3.1" ) THEN
!
!    rcode = nf90_get_att (cdfid, nf90_global, 'SF_URBAN_PHYSICS',  &
!                          met_urban_phys)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'SF_URBAN_PHYSICS',  &
!                      TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!
!  ELSE IF ( fv3version(18:21) == "V3.0" ) THEN
!
!    rcode = nf90_get_att (cdfid, nf90_global, 'UCMCALL', met_urban_phys)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'SF_URBAN_PHYSICS',  &
!                      TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!
!  ELSE 
!
!    ! In v2.2, header variable UCMCALL seems to always be 0 for nested runs,
!    ! even when UCM is invoked.  For now, use field TC_URB (canopy temperature)
!    ! as a proxy to determine if the UCM was used.  If the field does not exist,
!    ! then the UCM was not used.  If the field exists, determine if the data are
!    ! "reasonable" (i.e., positive and non-zero); assume that UCM was used if
!    ! the field contains "physical" data.
!
!    nxm = met_nx - 1
!    nym = met_ny - 1
!    it  = 1  ! use first time in file since some files just have one time
!    ALLOCATE ( dum2d ( nxm, nym ) )
!      CALL get_var_2d_real_cdf (cdfid, 'TC_URB', dum2d, it, rcode)
!      IF ( ( rcode == nf90_noerr ) .AND. ( MAXVAL(dum2d) > 100.0 ) ) THEN  ! UCM
!        met_urban_phys = 1
!      ELSE
!        met_urban_phys = 0
!      ENDIF
!    DEALLOCATE ( dum2d )
!
!  ENDIF
!
!  IF ( met_urban_phys >= 1 ) THEN
!    ifrcurb = .TRUE.
!  ELSE
!    ifrcurb = .FALSE.
!  ENDIF
!
!
!  ! Determine if shallow convection was used.
!  ! FV3 does not have shallow convection scheme, or anyway to determine, turn off
   met_shal_cu = 0
!  IF ( fv3version(18:21) >= "V3.3" ) THEN
!
!    rcode = nf90_get_att (cdfid, nf90_global, 'SHCU_PHYSICS', met_shal_cu)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'SHCU_PHYSICS', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!
!    IF ( met_shal_cu == 0 .AND. met_cumulus == 5 ) THEN  ! Grell shallow on?
!      rcode = nf90_get_att (cdfid, nf90_global, 'ISHALLOW', met_shal_cu)
!      IF ( rcode /= nf90_noerr ) THEN
!        WRITE (*,f9400) TRIM(pname), 'ISHALLOW', TRIM(nf90_strerror(rcode))
!        CALL graceful_stop (pname)
!      ENDIF
!    ENDIF
!
!  ELSE  ! no way to easily tell if Grell 3D used shallow convection
!
!    IF ( met_cumulus == 5 ) THEN  ! Grell 3D
!      met_shal_cu = -1
!    ELSE
!      met_shal_cu = 0
!    ENDIF
!
!  ENDIF

  met_snow_opt = 1       ! not used for FV3 yet, but set conistent with original MCIP
  met_rain_bucket = -1.0 ! FV3 includes both bucket and prate, so leave bucket turned off
  met_pcp_incr = 1       ! FV3 has time step average prate, so leave pcp_inc turned on

  !Determine if bucket or time-step precipitation accumulation rate exists

!  rcode = nf90_get_att (cdfid, nf90_global, 'BUCKET_MM', met_rain_bucket)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( fv3version(18:22) >= "V3.2" ) then  ! BUCKET_MM implemented in FV3v3.2
!      WRITE (*,f9400) TRIM(pname), 'BUCKET_MM', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ELSE
!      met_rain_bucket = -1.0
!    ENDIF
!  ENDIF
!
!  rcode = nf90_get_att (cdfid, nf90_global, 'PREC_ACC_DT', rval)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( fv3version(18:22) >= "V3.2" ) then  ! PREC_ACC_DT added in FV3v3.2
!      WRITE (*,f9400) TRIM(pname), 'PREC_ACC_DT', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ELSE
!      met_pcp_incr = 0
!    ENDIF
!  ELSE
!    rcode =         nf90_inq_varid (cdfid, 'PREC_ACC_C',  rcode)
!    rcode = rcode + nf90_inq_varid (cdfid, 'PREC_ACC_NC', rcode)
!    IF ( rcode /= nf90_noerr ) THEN
!      met_pcp_incr = 0
!    ELSE
!      met_pcp_incr = NINT(rval)
!    ENDIF
!  ENDIF
!
!  IF ( met_pcp_incr > 0 ) THEN
!    IF ( met_pcp_incr /= intvl ) THEN  ! can't compute precip for CMAQ
!      WRITE (*,f9550) TRIM(pname), met_pcp_incr, intvl
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF


  ! Determine if radiative feedbacks accompany the convective parameterization
  ! scheme.

  ! Currently there are no radiative feedbacks in FV3 accompaning the convective 
  ! parameterization scheme, turn off.
  ifcuradfdbk = .FALSE.

!  rcode = nf90_get_att (cdfid, nf90_global, 'ICLOUD_CU', ival)
!  IF ( rcode == nf90_noerr ) THEN  ! new enough version `of FV3
!    SELECT CASE ( ival )
!      CASE ( 0 )
!        ifcuradfdbk = .FALSE.
!      CASE ( 1:2 )  ! 1=Grell, 2=KF or MSKF
!        ifcuradfdbk = .TRUE.
!      CASE DEFAULT
!        ifcuradfdbk = .FALSE.
!    END SELECT
!  ELSE
!    ifcuradfdbk = .FALSE.
!  ENDIF

!-------------------------------------------------------------------------------
! Extract FV3 start date and time information.
!-------------------------------------------------------------------------------
  rcode = nf90_inq_varid (cdfid, 'time', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'time',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF 

  rcode = nf90_get_att (cdfid, varid, 'units', date_init)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'SIMULATION_START_DATE',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  met_startdate =  date_init(13:32) // '.0000'

  ! FV3 always turned off for "met restart" in MCIP
  met_restart = 0
 
!  rcode = nf90_get_att (cdfid, nf90_global, 'SIMULATION_START_DATE', date_init)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'SIMULATION_START_DATE',  &
!                    TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!  met_startdate =  date_init(1:19) // '.0000'
!  met_startdate(11:11) = "-"  ! change from "_" to "-" for consistency

!  rcode = nf90_get_att (cdfid, nf90_global, 'START_DATE', date_start)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'START_DATE', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!
!  IF ( date_init == date_start ) THEN
!    met_restart = 0
!  ELSE
!    met_restart = 1
!  ENDIF

!  rcode = nf90_inq_varid (cdfid, 'time', id_data)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9410) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!  rcode = nf90_inquire_variable (cdfid, id_data, dimids=dimids)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9420) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!  rcode = nf90_inquire_dimension (cdfid, dimids(1), len=lent)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9430) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF
!  rcode = nf90_inquire_dimension (cdfid, dimids(2), len=n_times)
!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9430) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF

  rcode = nf90_inquire_variable (cdfid, varid, dimids=dimids)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9420) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inquire_dimension (cdfid, dimids(1), len=n_times)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9430) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  IF ( ALLOCATED ( times ) ) DEALLOCATE ( times )
  ALLOCATE ( times ( n_times ) )
  rcode = nf90_get_var (cdfid, varid, times)
!                        start=(/1/), count=(/n_times/))
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

!  IF ( n_times > 1 ) THEN
!    CALL geth_idts (times(2)(1:19), times(1)(1:19), idtsec)
!  ELSE
!    fl2 = file_mm(2)
!    IF ( fl2(1:10) == '          ' ) THEN
!      WRITE (*,f9500) TRIM(pname)
!      idtsec = 60
!    ELSE
!      rcode = nf90_open (fl2, nf90_nowrite, cdfid2)
!      IF ( rcode == nf90_noerr ) THEN
!        rcode = nf90_inq_varid (cdfid2, 'time', id_data)
!        IF ( rcode /= nf90_noerr ) THEN
!          WRITE (*,f9410) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!          CALL graceful_stop (pname)
!        ENDIF
!        rcode = nf90_inquire_variable (cdfid2, id_data, dimids=dimids)
!        IF ( rcode /= nf90_noerr ) THEN
!          WRITE (*,f9420) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!          CALL graceful_stop (pname)
!        ENDIF
!        rcode = nf90_inquire_dimension (cdfid2, dimids(1), len=lent)
!        IF ( rcode /= nf90_noerr ) THEN
!          WRITE (*,f9430) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!          CALL graceful_stop (pname)
!        ENDIF
!        rcode = nf90_inquire_dimension (cdfid2, dimids(2), len=n_times)
!        IF ( rcode /= nf90_noerr ) THEN
!          WRITE (*,f9430) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!          CALL graceful_stop (pname)
!        ENDIF
!        IF ( ALLOCATED ( times2 ) ) DEALLOCATE ( times2 )
!        ALLOCATE ( times2 ( n_times ) )
!        rcode = nf90_get_var (cdfid2, id_data, times2,   &
!                              start=(/1,1/), count=(/lent,n_times/))
!        IF ( rcode /= nf90_noerr ) THEN
!          WRITE (*,f9400) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
!          CALL graceful_stop (pname)
!        ENDIF
!        CALL geth_idts (times2(1)(1:19), times(1)(1:19), idtsec)
!      ELSE
!        WRITE (*,f9600) TRIM(pname), TRIM(fl2)
!        CALL graceful_stop (pname)
!      ENDIF
!      rcode = nf90_close (cdfid2)
!      IF ( rcode /= nf90_noerr ) THEN
!        WRITE (*,f9700) TRIM(pname), TRIM(fl2)
!        CALL graceful_stop (pname)
!      ENDIF
!    ENDIF
!  ENDIF

   met_tapfrq = (times(1))*60  ! convert hrs --> min
   IF ( met_model == 3 ) THEN !IF FV3 and first 000 forecast hour, assume hourly input = 60 min
    IF (met_tapfrq == 0.0) THEN
     met_tapfrq = 60.0 !minutes
    ENDIF
   ENDIF
   
!-------------------------------------------------------------------------------
! Set variables for non-hydrostatic base state.  There is no option for
! hydrostatic run in FV3.  The base state variables are not currently output
! , so fill in "default" values from FV3 namelist.
!
! Note:  In FV3v2.2 NCAR changed the way "real" scalars (e.g., P_TOP) are
!        stored in the FV3 I/O API.
!-------------------------------------------------------------------------------

!  IF ( (fv3version(18:21) == "V2.2") .OR. (fv3version(18:19) >= "V3") ) THEN
!    CALL get_var_real_cdf (cdfid, 'P_TOP', met_ptop, rcode)
!  ELSE
!    ALLOCATE ( dum1d ( 1 ) )
!    CALL get_var_1d_real_cdf (cdfid, 'P_TOP', dum1d, 1, rcode)
!    met_ptop = dum1d(1)
!  ENDIF

!  IF ( rcode /= nf90_noerr ) THEN
!    WRITE (*,f9400) TRIM(pname), 'P_TOP', TRIM(nf90_strerror(rcode))
!    CALL graceful_stop (pname)
!  ENDIF

!  met_ptop  = phalf_lays(nlays)*100.0 !FV3 set met_ptop to first value of phalf array
!  met_ptop  = pfull_lays(nlays+1)*100.0 ! FV3 top-down hPa --> [Pa]
  met_ptop  = phalf_lays(nlays)*100.0
  met_p00   = 100000.0 ! base state sea-level pressure [Pa]
  met_ts0   =    290.0 ! base state sea-level temperature [K]
  met_tlp   =     50.0 ! base state lapse rate d(T)/d(ln P) from 1000 to 300 mb
  met_tiso  = fillreal ! base state stratospheric isothermal T [K]  ! not used


!-------------------------------------------------------------------------------
! If Eta layer structure (ctmlays) was not defined in user namelist, calculate using FV3 pressure.
!-------------------------------------------------------------------------------

  IF ( needlayers ) THEN
    !FV3 is top down, but CMAQ levels are bottom up, reverse pressure array order:
!    ctmlays = (pfull_lays(nlays+1:1:-1) - pfull_lays(1)) / (MAXVAL(pfull_lays) - pfull_lays(1))
    !Flip again to top down to be consistent with other FV3 vertical grid.
!    ctmlays = ctmlays(nlays+1:1:-1)
!     ctmlays = (pfull_lays - pfull_lays(nlays)) / (MAXVAL(pfull_lays) - pfull_lays(nlays))
     ctmlays = (phalf_lays - phalf_lays(nlays)) / ((phalf_lays(1)) - phalf_lays(nlays))
  ENDIF
!-------------------------------------------------------------------------------
! Determine FV3 release.
!-------------------------------------------------------------------------------

  met_release = ' FV3  '

!  IF ( fv3version(18:18) == "V" ) THEN
!    met_release(1:2) = fv3version(18:19)
!  ENDIF
!
!  IF ( fv3version(20:20) == '.' ) THEN
!    met_release(3:4) = fv3version(20:21)
!  ENDIF
!
!  IF ( fv3version(22:22) == '.' ) THEN
!    met_release(5:6) = fv3version(22:23)
!  ENDIF
!
!  IF ( fv3version(24:24) == '.' ) THEN
!    met_release(7:8) = fv3version(24:25)
!  ENDIF

!-------------------------------------------------------------------------------
! Determine FDDA options.
!-------------------------------------------------------------------------------
  met_fdda_3dan = 0
!  rcode = nf90_get_att (cdfid, nf90_global, 'GRID_FDDA', met_fdda_3dan)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V2.2' ) THEN
!      met_fdda_3dan = 0  ! not implemented until V2.2
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'GRID_FDDA', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF
  met_fdda_gv3d = -1.0 
!  rcode = nf90_get_att (cdfid, nf90_global, 'GUV', met_fdda_gv3d)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V2.2' ) THEN
!      met_fdda_gv3d = -1.0  ! not in header until V2.2
!    ELSE IF ( met_fdda_3dan == 0 ) THEN
!      met_fdda_gv3d = -1.0  ! not in header if analysis nudging is off
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'GUV', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF
  met_fdda_gt3d = -1.0
!  rcode = nf90_get_att (cdfid, nf90_global, 'GT', met_fdda_gt3d)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V2.2' ) THEN
!      met_fdda_gt3d = -1.0  ! not in header until V2.2
!    ELSE IF ( met_fdda_3dan == 0 ) THEN
!      met_fdda_gt3d = -1.0  ! not in header if analysis nudging is off
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'GT', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF
  met_fdda_gq3d = -1.0
!  rcode = nf90_get_att (cdfid, nf90_global, 'GQ', met_fdda_gq3d)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V2.2' ) THEN
!      met_fdda_gq3d = -1.0  ! not in header until V2.2
!    ELSE IF ( met_fdda_3dan /= 1 ) THEN
!      met_fdda_gq3d = -1.0  ! not in header if analysis nudging is off
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'GQ', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF
  met_fdda_gph3d = -1.0
!  rcode = nf90_get_att (cdfid, nf90_global, 'GPH', met_fdda_gph3d)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V3.2' ) THEN
!      met_fdda_gph3d = -1.0  ! not in header until V3.2
!    ELSE IF ( met_fdda_3dan /= 2 ) THEN
!      met_fdda_gph3d = -1.0  ! not in header if spectral nudging is off
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'GPH', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF
!
!  IF ( TRIM(met_release) >= 'V3.1' ) THEN  ! find sfc analysis nudging info
!
!    rcode = nf90_get_att (cdfid, nf90_global, 'GRID_SFDDA', met_fdda_sfan)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'GRID_SFDDA', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!
!    IF ( met_fdda_sfan == 1 ) THEN
!
!      rcode = nf90_get_att (cdfid, nf90_global, 'GUV_SFC', met_fdda_gvsfc)
!      IF ( rcode /= nf90_noerr ) THEN
!        WRITE (*,f9400) TRIM(pname), 'GUV_SFC', TRIM(nf90_strerror(rcode))
!        CALL graceful_stop (pname)
!      ENDIF
!
!      rcode = nf90_get_att (cdfid, nf90_global, 'GT_SFC', met_fdda_gtsfc)
!      IF ( rcode /= nf90_noerr ) THEN
!        WRITE (*,f9400) TRIM(pname), 'GT_SFC', TRIM(nf90_strerror(rcode))
!        CALL graceful_stop (pname)
!      ENDIF
!
!      rcode = nf90_get_att (cdfid, nf90_global, 'GQ_SFC', met_fdda_gqsfc)
!      IF ( rcode /= nf90_noerr ) THEN
!        WRITE (*,f9400) TRIM(pname), 'GQ_SFC', TRIM(nf90_strerror(rcode))
!        CALL graceful_stop (pname)
!      ENDIF
!
!    ELSE
!
!      met_fdda_gvsfc = -1.0
!      met_fdda_gtsfc = -1.0
!      met_fdda_gqsfc = -1.0
!
!    ENDIF
!
!  ELSE
    met_fdda_sfan  =  0  ! sfc analysis nudging not in FV3
    met_fdda_gvsfc = -1.0
    met_fdda_gtsfc = -1.0
    met_fdda_gqsfc = -1.0
!  ENDIF
  met_fdda_obs = 0       ! not implemented in FV3
!  rcode = nf90_get_att (cdfid, nf90_global, 'OBS_NUDGE_OPT', met_fdda_obs)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V2.2' ) THEN
!      met_fdda_obs = 0  ! not implemented until V2.2
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'OBS_NUDGE_OPT', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF
  met_fdda_giv = -1.0    ! not implemented in FV3
!  rcode = nf90_get_att (cdfid, nf90_global, 'OBS_COEF_WIND', met_fdda_giv)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V2.2' ) THEN
!      met_fdda_giv = -1.0  ! not in header until V2.2
!    ELSE IF ( met_fdda_obs == 0 ) THEN
!      met_fdda_giv = -1.0  ! not in header if obs nudging is off
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'OBS_COEF_WIND', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF
  met_fdda_git = -1.0    ! not implemented in FV3
!  rcode = nf90_get_att (cdfid, nf90_global, 'OBS_COEF_TEMP', met_fdda_git)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V2.2' ) THEN
!      met_fdda_git = -1.0  ! not in header until V2.2
!    ELSE IF ( met_fdda_obs == 0 ) THEN
!      met_fdda_git = -1.0  ! not in header if obs nudging is off
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'OBS_COEF_TEMP', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF
  met_fdda_giq = -1.0    ! not implemented in FV3
!  rcode = nf90_get_att (cdfid, nf90_global, 'OBS_COEF_MOIS', met_fdda_giq)
!  IF ( rcode /= nf90_noerr ) THEN
!    IF ( TRIM(met_release) < 'V2.2' ) THEN
!      met_fdda_giq = -1.0  ! not in header until V2.2
!    ELSE IF ( met_fdda_obs == 0 ) THEN
!      met_fdda_giq = -1.0  ! not in header if obs nudging is off
!    ELSE
!      WRITE (*,f9400) TRIM(pname), 'OBS_COEF_MOIS', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ENDIF

!-------------------------------------------------------------------------------
! Determine whether or not fractional land use is available in the output.
! Set the flag appropriately.
!-------------------------------------------------------------------------------
  !Leave check in FV3 MCIP version in case fractional land use becomes available,
  ! and it does not stop model if not available
  rcode2 = nf90_inq_varid (cdfid2, 'LANDUSEF', varid)
  IF ( rcode2 == nf90_noerr ) THEN
    iflufrc    = .TRUE.   ! fractional land use is available
    ifluwrfout = .TRUE.   ! fractional land use is located in FV3 history file
  ELSE
    iflufrc    = .FALSE.   ! fractional land use is not available
    ifluwrfout = .FALSE.  ! fractional land use is not available in FV3 history
    geofile = TRIM( file_geo )
    INQUIRE ( FILE=geofile, EXIST=ifgeo )
    IF ( .NOT. ifgeo ) THEN
      WRITE (*,f9800) TRIM(pname)
      iflufrc = .FALSE.
    ELSE
      flg = file_geo
      rcode = nf90_open (flg, nf90_nowrite, cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9600) TRIM(pname), TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
      CALL chkwpshdr (flg, cdfidg)
      rcode = nf90_inq_varid (cdfidg, 'LANDUSEF', varid)
      IF ( rcode == nf90_noerr ) THEN
        iflufrc = .TRUE.  ! fractional land use is in the file
      ELSE
        iflufrc = .FALSE. ! fractional land use is not in the file
      ENDIF
      rcode = nf90_close (cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9700)  TRIM(pname),TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
    ENDIF
  ENDIF  

  ! For now, require LANDUSEF2 and MOSAIC_CAT_INDEX to process NOAH Mosaic.
  ! IFMOSAIC is toggled to FALSE here if either field is missing.
  ! Can add a subroutine later to reconstruct those fields from LANDUSEF if
  ! LANDUSEF2 and/or MOSAIC_CAT_INDEX is missing.

  IF ( ifmosaic ) THEN
    rcode = nf90_inq_varid (cdfid, 'LANDUSEF2', varid)
    IF ( rcode == nf90_noerr ) THEN
      rcode = nf90_inq_varid (cdfid, 'MOSAIC_CAT_INDEX', varid)
      IF ( rcode == nf90_noerr ) THEN
        iflu2wrfout = .TRUE.   ! lookup for LANDUSEF2 is in FV3 history
      ELSE
        iflu2wrfout = .FALSE.
        ifmosaic    = .FALSE.
      ENDIF
    ELSE
      iflu2wrfout = .FALSE.  ! frac land use 2 is not available in FV3 history
      ifmosaic    = .FALSE.
    ENDIF
  ELSE
    iflu2wrfout = .FALSE.
  ENDIF
!-------------------------------------------------------------------------------
! Determine whether or not the 2-m temperature, the 2-m mixing ratio, the
! 10-m wind components, and the turbulent kinetic energy are in the output,
! and set the flags appropriately.
!-------------------------------------------------------------------------------

  rcode = nf90_inq_varid (cdfid2, 'tmp2m', varid)
  IF ( rcode == nf90_noerr ) THEN
    ift2m = .TRUE.  ! 2-m temperature is in the file
  ELSE
    ift2m = .FALSE. ! 2-m temperature is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'spfh2m', varid)
  IF ( rcode == nf90_noerr ) THEN
    IF ( met_pbl == 1 ) THEN  ! YSU PBL scheme
      ifq2m = .FALSE. ! do not use Q2 from YSU PBL; occasional winter negatives
    ELSE
      ifq2m = .TRUE.  ! 2-m mixing ratio is in the file
    ENDIF
  ELSE
    ifq2m = .FALSE. ! 2-m mixing ratio is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'ugrd10m', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifu10m = .TRUE.  ! 10-m u-component wind is in the file
  ELSE
    ifu10m = .FALSE. ! 10-m u-component wind is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'vgrd10m', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifv10m = .TRUE.  ! 10-m v-component wind is in the file
  ELSE
    ifv10m = .FALSE. ! 10-m v-component wind is not in the file
  ENDIF

  IF ( ( ifu10m ) .AND. ( ifv10m ) ) THEN
    ifw10m = .TRUE.
  ELSE
    ifw10m = .FALSE.
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'TKE_MYJ', varid)  !Left as TKE_MYJ, not in FV3 = False
  IF ( rcode == nf90_noerr ) THEN
    IF ( met_pbl == 2 ) THEN  ! Mellor-Yamada-Janjic (Eta)
      iftke  = .TRUE.  ! turbulent kinetic energy is in the file
      iftkef = .FALSE. ! TKE is not on full-levels; it is on half-layers
    ELSE
      iftke  = .FALSE. ! turbulent kinetic energy is not in the file
      iftkef = .FALSE.
    ENDIF
  ELSE
    iftke  = .FALSE. ! turbulent kinetic energy is not in the file
    iftkef = .FALSE.
  ENDIF
!-------------------------------------------------------------------------------
! Determine whether or not some surface variables are in the output, and set
! the flags appropriately.
!-------------------------------------------------------------------------------

  rcode = nf90_inq_varid (cdfid2, 'LAI', varid)  !not in FV3GFSv16
  IF ( rcode == nf90_noerr ) THEN
    iflai = .TRUE.  ! leaf area index is in the file
  ELSE
    iflai = .FALSE. ! leaf area index is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'RMOL', varid) !not in FV3GFSv16
  IF ( rcode == nf90_noerr ) THEN
    ifmol = .TRUE.  ! (inverse) Monin-Obukhov length is in the file
  ELSE
    ifmol = .FALSE. ! (inverse) Monin-Obukhov length is not in the file
  ENDIF

!  IF ( met_soil_lsm == 7 ) THEN  ! PX LSM is not in FV3. :)
!    ifmolpx = .TRUE.
!    rcode = nf90_inq_varid (cdfid, 'LAI_PX', varid)  ! there are 7 variables
!    rcode = rcode + nf90_inq_varid (cdfid, 'WWLT_PX', varid)
!    rcode = rcode + nf90_inq_varid (cdfid, 'WFC_PX', varid)
!    rcode = rcode + nf90_inq_varid (cdfid, 'WSAT_PX', varid)
!    rcode = rcode + nf90_inq_varid (cdfid, 'CSAND_PX', varid)
!    rcode = rcode + nf90_inq_varid (cdfid, 'FMSAND_PX', varid)
!    rcode = rcode + nf90_inq_varid (cdfid, 'CLAY_PX', varid)
!    IF ( rcode == nf90_noerr ) THEN
!      ifpxwrf41 = .TRUE.  !  all 7 variables are available
!    ELSE
!      ifpxwrf41 = .FALSE.
!    ENDIF
!  ELSE
    ifmolpx   = .FALSE.
    ifpxwrf41 = .FALSE.
!  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'acond', varid) 
  IF ( rcode == nf90_noerr ) THEN
    ifra = .TRUE.  ! aerodynamic resistance is in the file
  ELSE
    ifra = .FALSE. ! aerodynamic resistance is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'RS', varid)  !not in FV3GFSv16??
  IF ( rcode == nf90_noerr ) THEN
    ifrs = .TRUE.  ! stomatal resistance is in the file
  ELSE
    ifrs = .FALSE. ! stomatal resistance is not in the file
  ENDIF

  IF ( ( ifra ) .AND. ( ifrs ) ) THEN
    ifresist = .TRUE.
  ELSE
    ifresist = .FALSE.
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'veg', varid) 
  IF ( rcode == nf90_noerr ) THEN
    ifveg = .TRUE.  ! vegetation fraction is in the file
  ELSE
    ifveg = .FALSE. ! vegetation fraction is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'cnwat', varid) 
  IF ( rcode == nf90_noerr ) THEN
    ifwr = .TRUE.  ! canopy wetness is in the file
  ELSE
    ifwr = .FALSE. ! canopy wetness is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'soill1', varid) 
  IF ( rcode == nf90_noerr ) THEN
    ifsmois = .TRUE.  ! soil moisture is in the file
  ELSE
    ifsmois = .FALSE. ! soil moisture is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'soilt1', varid)
  IF ( rcode == nf90_noerr ) THEN
    iftslb = .TRUE.  ! soil temperature is in the file
  ELSE
    iftslb = .FALSE. ! soil temperature is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'sotyp', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifisltyp = .TRUE.  ! soil type is in the file
  ELSE
    ifisltyp = .FALSE. ! soil type is not in the file
  ENDIF

  If ( ( ifsmois ) .AND. ( iftslb ) .AND. ( ifisltyp ) ) THEN
    ifsoil = .TRUE.
  ELSE
    ifsoil = .FALSE.
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'sfcr', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifznt = .TRUE.  ! roughness length is in the file
  ELSE
    ifznt = .FALSE. ! roughness length is not in the file
  ENDIF

!KF-Radiative Feedback is not in FV3

!  IF ( met_cumulus == 1 .AND. ifcuradfdbk ) THEN  ! KF-Rad was used in FV3
!    rcode = nf90_inq_varid (cdfid, 'QC_CU', varid)  ! there are 4 variables
!    rcode = rcode + nf90_inq_varid (cdfid, 'QI_CU', varid)
!    rcode = rcode + nf90_inq_varid (cdfid, 'CLDFRA_DP', varid)
!    rcode = rcode + nf90_inq_varid (cdfid, 'CLDFRA_SH', varid)
!   IF ( rcode == nf90_noerr ) THEN
!      ifkfradextras = .TRUE.  !  all 4 variables are available
!    ELSE
!      ifkfradextras = .FALSE.
!    ENDIF
!  ELSE
    ifkfradextras = .FALSE.
!  ENDIF

!-------------------------------------------------------------------------------
! Determine the number of 3D cloud moisture species.  Assume that cloud water
! mixing ratio and rain water mixing ratio will occur together.  Also assume
! that cloud ice mixing ratio and cloud snow mixing ratio will occur together,
! but check for availability.  Check for graupel, as well.
! Note:  In FV3v2.1.2 and prior, the Eta/Ferrier microphysics scheme only
! outputs QCLOUD which represents total condensate, not cloud water mixing
! ratio.  CMAQv4.6 and prior cannot handle this field, so MCIP will stop in
! this case.
!-------------------------------------------------------------------------------

  rcode = nf90_inq_varid (cdfid, 'clwmr', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = 1  ! QCLOUD is in the file
  ELSE  ! need hydrometeor fields for CMAQ
    WRITE (*,f9225) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'rwmr', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = nqspecies + 1  ! QRAIN is in the file
  ELSE
    IF ( met_expl_moist == 5 ) THEN  ! Eta/Ferrier scheme
      WRITE (*,f9250) TRIM(pname)
      CALL graceful_stop (pname)
    ELSE
      WRITE (*,f9275) TRIM(pname)
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'icmr', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = nqspecies + 1  ! QICE is in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'snmr', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = nqspecies + 1  ! QSNOW is in the file
  ENDIF

  IF ( nqspecies == 3 ) THEN  ! not set up for QI w/o QS or vice versa
    WRITE (*,f9300) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'grle', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = nqspecies + 1  ! QGRAUP is in the file
  ENDIF

  IF ( nqspecies == 3 ) THEN  ! not set up for QG without QI and QS
    WRITE (*,f9300) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

!-------------------------------------------------------------------------------
! Determine whether 3D resolved cloud fraction is part of FV3 output.  If
! Kain-Fritsch scheme with radiative feedbacks to subgrid clouds is used (new
! in FV3v3.6) or if MSKF is used (new in FV3v3.7) in FV3, then the 3D cloud
! fraction includes both resolved and subgrid clouds.
!-------------------------------------------------------------------------------

  rcode = nf90_inq_varid (cdfid, 'cld_amt', varid)
  IF ( rcode == nf90_noerr ) THEN
!    IF ( TRIM(met_release) >= 'V3.6' ) THEN  !FV3 does not have cloud feedbacks
!      rcode = nf90_get_att (cdfid, nf90_global, 'ICLOUD_CU', icloud_cu)
!      IF ( rcode == nf90_noerr ) THEN
!        IF ( ( ( met_cumulus ==  1 ) .AND. ( icloud_cu == 2 ) ) .OR. &
!               ( met_cumulus == 11 ) ) THEN
!          ifcld3d = .FALSE.  ! 3D resolved cloud fraction is not in the file
!        ELSE
          ifcld3d = .TRUE.  ! 3D resolved cloud fraction is in the file
!        ENDIF
!      ELSE
!        ifcld3d = .TRUE.  ! 3D resolved cloud fraction is in the file
!      ENDIF
!    ELSE
!      ifcld3d = .TRUE.  ! 3D resolved cloud fraction is in the file
!    ENDIF
  ELSE
    ifcld3d = .FALSE. ! 3D cloud fraction is not if the file
  ENDIF
!-------------------------------------------------------------------------------
! Determine if the hybrid vertical coordinate has been used in FV3.  It is
! available as of FV3v3.9.
!-------------------------------------------------------------------------------
   met_hybrid = 2   ! FV3 employs a hybrid vertical coordinate as default = 2 (mimic WRF hybrid flag)
!  IF ( TRIM(met_release) >= "V3.9") THEN
!    rcode = nf90_get_att (cdfid, nf90_global, 'HYBRID_OPT', met_hybrid)
!    IF ( rcode /= nf90_noerr ) THEN
!      WRITE (*,f9400) TRIM(pname), 'HYBRID_OPT', TRIM(nf90_strerror(rcode))
!      CALL graceful_stop (pname)
!    ENDIF
!  ELSE
!    met_hybrid = -1
!  ENDIF
  met_hybrid = 2
END SUBROUTINE setup_fv3
