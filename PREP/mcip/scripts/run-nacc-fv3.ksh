#/bin/ksh -x
APPL=Test_FV3
CoordName=FV3_RPO	      # 16-character maximum
GridName=FV3_CONUS	      # 16-character maximum

InMetDir=/gpfs/hps2/ptmp/Patrick.C.Campbell/fv3gfs_v16_test/12z_hourly
OutDir=/gpfs/hps2/ptmp/Patrick.C.Campbell/fv3gfs_v16_test
ProgDir=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/CMAQ_REPO/PREP/mcip/src
WorkDir=$OutDir

if [ ! -s $InMetDir ]; then
  echo "No such input directory $InMetDir"
  exit 1
fi

if [ ! -d $OutDir ]; then
  echo "No such output directory...will try to create one"
  mkdir -p $OutDir
  if [ $? != 0 ]; then
    echo "Failed to make output directory, $OutDir"
    exit 1
  fi
fi

if [ ! -d $ProgDir ]; then
  echo "No such program directory $ProgDir"
  exit 1
fi

cd $OutDir
cat>namelist.mcip<<!
&FILENAMES
  file_gd    = 'GRIDDESC'
  file_mm    = '$InMetDir/gfs.t12z.atmf','.nc'
  file_sfc   = '$InMetDir/gfs.t12z.sfcf','.nc'
  ioform     =  1
 &END

 &USERDEFS
  inmetmodel =  3
  dx_in      =  12000
  dy_in      =  12000
  met_cen_lat_in =  0.0
  met_cen_lon_in =  0.0
  lpv        =  0
  lwout      =  0
  luvbout    =  1
  mcip_start = "2020-01-12-12:00:00.0000"
  mcip_end   = "2020-01-15-13:00:00.0000"
  intvl      =  60
  coordnam   = "FV3_RPO"
  grdnam     = "FV3_CONUS"
  ctmlays    =  1.000000, 0.995253, 0.990479, 0.985679, 0.980781,
              0.975782, 0.970684, 0.960187, 0.954689, 0.936895, 
              0.930397, 0.908404, 0.888811, 0.862914, 0.829314, 
              0.786714, 0.735314, 0.645814, 0.614214, 0.582114, 
              0.549714, 0.511711, 0.484394, 0.451894, 0.419694, 
              0.388094, 0.356994, 0.326694, 0.297694, 0.270694, 
              0.245894, 0.223694, 0.203594, 0.154394, 0.127094, 0.000000
  btrim      =  -1
  lprt_col   =  0
  lprt_row   =  0
  ntimes     = 72
  wrf_lc_ref_lat = 40.0
  projparm = 2., 33.,45., -97., -97., 40.
  domains = -2508000., -1716000., 12000., 12000., 442, 265
 &END

 &WINDOWDEFS
  x0         =  1
  y0         =  1
  ncolsin    =  442
  nrowsin    =  265
 &END
!

export IOAPI_CHECK_HEADERS=T
export EXECUTION_ID=$PROG

export GRID_BDY_2D=GRIDBDY2D_${APPL}.nc
export GRID_CRO_2D=GRIDCRO2D_${APPL}.nc
export GRID_DOT_2D=GRIDDOT2D_${APPL}.nc
export MET_BDY_3D=METBDY3D_${APPL}.nc
export MET_CRO_2D=METCRO2D_${APPL}.nc
export MET_CRO_3D=METCRO3D_${APPL}.nc
export MET_DOT_3D=METDOT3D_${APPL}.nc
export LUFRAC_CRO=LUFRAC_CRO_${APPL}.nc
export SOI_CRO=SOI_CRO_${APPL}.nc
export MOSAIC_CRO=MOSAIC_CRO_${APPL}.nc

rm -f *.nc 
$ProgDir/mcip.exe
