#------------------------------------------------------------------------------#
#  The Community Multiscale Air Quality (CMAQ) system software is in           #
#  continuous development by various groups and is based on information        #
#  from these groups: Federal Government employees, contractors working        #
#  within a United States Government contract, and non-Federal sources         #
#  including research institutions.  These groups give the Government          #
#  permission to use, prepare derivative works of, and distribute copies       #
#  of their work in the CMAQ system to the public and to permit others         #
#  to do so.  The United States Environmental Protection Agency                #
#  therefore grants similar permission to use the CMAQ system software,        #
#  but users are requested to provide copies of derivative works or            #
#  products designed to operate in the CMAQ system to the United States        #
#  Government without restrictions as to use by others.  Software              #
#  that is used with the CMAQ system but distributed under the GNU             #
#  General Public License or the GNU Lesser General Public License is          #
#  subject to their copyright restrictions.                                    #
#------------------------------------------------------------------------------#

.SUFFIXES:
.SUFFIXES: .o .f90 .F90

MODEL = mcip.exe

#...Portland Group Fortran
#FC     = /usr/local/apps/pgi/linux86-64/17.10/bin/pgf90
#NETCDF = /usr/local/apps/netcdf-4.6.3/pgi-17.10
#IOAPI_ROOT = /usr/local/apps/ioapi-3.2_20181011/pgi-17.10
###FFLAGS = -g -O0 -Ktrap=unf -Ktrap=denorm -Ktrap=inv -Ktrap=ovf -Ktrap=divz  \
###         -Ktrap=fp -pc 32 -Mbounds -Mchkfpstk -Mchkptr -Kieee     \
###         -Minform,inform -Mfree -byteswapio -I$(NETCDF)/include \
###         -I$(IOAPI_ROOT)/Linux2_x86_64pg
#FFLAGS = -O4 -fastsse -pc 32 -Mfree -byteswapio -I$(NETCDF)/include \
#         -I$(IOAPI_ROOT)/Linux2_x86_64pg
#LIBS   = -L$(IOAPI_ROOT)/Linux2_x86_64pg -lioapi \
#	 -L$(NETCDF)/lib -lnetcdff -lnetcdf

#...gfortran
#FC	= gfortran
#NETCDF = /usr/local/apps/netcdf-4.6.3/gcc-6.1.0
#IOAPI_ROOT = /usr/local/apps/ioapi-3.2_20181011/gcc-6.1.0
#FFLAGS	= -O3 -I$(NETCDF)/include -I$(IOAPI_ROOT)/Linux2_x86_64
###FFLAGS	= -g -O0  \
###          -ffpe-trap='invalid','zero','overflow','underflow'  \
###          -I$(NETCDF)/include -I$(IOAPI_ROOT)/Linux2_x86_64
#LIBS    = -L$(IOAPI_ROOT)/Linux2_x86_64 -lioapi  \
#          -L$(NETCDF)/lib -lnetcdff -lnetcdf

#...Intel Fortran
FC      = ifort
#FC      = ftn

#...NetCDF/HDF5 Settings
#NETCDF = /usr/local/apps/netcdf-4.6.3/intel-18.0
#...Serial netCDF I/O

#WCOSS-Cray
#NETCDF =  /opt/cray/netcdf/4.3.2/INTEL/140/
#ARL-AAQEST-Local
#NETCDF =  /opt/netcdf4
#Hera-Princeton
NETCDF = /apps/netcdf/4.7.0/intel/18.0.5.274

#...Parallel netCDF I/O with HDF5
#NETCDF =  /opt/cray/netcdf-hdf5parallel/4.4.1.1/INTEL/150
#NETCDF  = /gpfs/hps/usrx/local/prod/NetCDF/4.2/intel/sandybridge
#HDF5    = /usrx/local/prod/HDF5/1.8.9/parallel/intel/sandybridge

#NETCDF_FORTAN = /gpfs/hps/usrx/local/prod/NetCDF/3.6.3/intel/haswell/src/netcdf-3.6.3-beta1/fortran/
#IOAPI_ROOT = /usr/local/apps/ioapi-3.2/intel-18.0


#...IOAPI Settings
#IOAPI_ROOT = /gpfs/hps3/emc/naqfc/noscrub/Jianping.Huang/nwdev/cmaq.v5.1.8/lib/
#IOAPI_ROOT = /gpfs/hps/nco/ops/nwprod/cmaq.v5.1.4/lib/

#WCOSS-Cray
#IOAPI_ROOT = /gpfs/hps3/emc/naqfc/noscrub/Youhua.Tang/CMAQ/ioapi-3.2/Linux2_x86_64ifort
#ARL-AAQEST-Local
#IOAPI_ROOT = /data/aqf/youhuat/CMAQ/ioapi-3.2/Linux2_x86_64ifort
#Hera-Princeton
IOAPI_ROOT = /scratch2/NCEPDEV/naqfc/Youhua.Tang/CMAQ/ioapi-3.2/Linux2_x86_64ifort

#...MPI Settings
#MPI_ROOT = /opt/cray/mpt/7.2.0/gni/mpich2-intel/140

###FFLAGS  = -g -O0 -check all -C -traceback -FR -I$(NETCDF)/include  \
###          -I$(IOAPI_ROOT)/Linux2_x86_64ifort
#FFLAGS  = -FR -O3 -traceback -I$(NETCDF)/include -I$(IOAPI_ROOT)/Linux2_x86_64ifort

#...Serial
#FFLAGS  = -FR -O3 -traceback -check all -I$(NETCDF)/include -I$(IOAPI_ROOT)/include
#FFLAGS  = -FR -O3 -traceback -I$(NETCDF)/include -I$(IOAPI_ROOT)/include 

FFLAGS  = -FR -O3 -traceback -check all -I$(NETCDF)/include -I$(IOAPI_ROOT)
#FFLAGS  = -FR -O3 -traceback -I$(NETCDF)/include -I$(IOAPI_ROOT)

#...Parallel
#FFLAGS  = -FR -O3 -traceback -I$(NETCDF)/include -I$(IOAPI_ROOT)/include -I$(MPI_ROOT)/include
#FFLAGS  = -FR -O3 -traceback -check all -I$(NETCDF)/include -I$(IOAPI_ROOT)/include -I$(MPI_ROOT)/include
#FFLAGS  = -FR -O3 -traceback -I$(NETCDF)/include -I$(HDF5)/include -I$(IOAPI_ROOT)/include -I$(MPI_ROOT)/include 


#WCOSS-Cray
#LIBS    = -L$(IOAPI_ROOT) -laqm_ioapi \
          -L$(NETCDF)/lib -lnetcdf -lnetcdff -openmp

#ARL-AAQEST-Local
#LIBS    = -L$(IOAPI_ROOT) -lioapi \
          -L$(NETCDF)/lib -lnetcdf -lnetcdff -qopenmp

#Hera-Princeton
LIBS    = -L$(IOAPI_ROOT) -lioapi \
          -L$(NETCDF)/lib -lnetcdf -lnetcdff -qopenmp


#          -L$(HDF5)/lib -lhdf5
#          -L$(NETCDF_FORTRAN) -lnetcdff
#          -L$(NETCDF)/lib -lnetcdf
#          -L$(NETCDF)/lib -lnetcdf -lnetcdff



/EFS   = 


MODULES =\
 const_mod.o \
 const_pbl_mod.o \
 coord_mod.o \
 ctmvars_mod.o \
 date_time_mod.o \
 date_pack_mod.o \
 files_mod.o \
 lucats_mod.o \
 mcipparm_mod.o \
 metinfo_mod.o \
 metvars_mod.o \
 netcdf_io_mod.o \
 vgrd_mod.o \
 xvars_mod.o
 
OBJS =\
 mcip.o \
 alloc_ctm.o \
 alloc_met.o \
 alloc_x.o \
 bcldprc_ak.o \
 blddesc.o \
 chkwpshdr.o \
 chkwrfhdr.o \
 chkfv3hdr.o \
 close_files.o \
 collapx.o \
 comheader.o \
 comheader_lufrac.o \
 comheader_mos.o \
 comheader_soi.o \
 ctmout.o \
 ctmproc.o \
 dealloc_ctm.o \
 dealloc_met.o \
 dealloc_x.o \
 detangle_soil_px.o \
 e_aerk.o \
 dynflds.o \
 getluse.o \
 getmet.o \
 getpblht.o \
 getpsih.o \
 getsdt.o \
 getxyindex.o \
 graceful_stop.o \
 gridout.o \
 gridproc.o \
 init_ctm.o \
 init_io.o \
 init_met.o \
 init_x.o \
 julian.o \
 layht.o \
 ll2xy_lam.o \
 ll2xy_lam_sec.o \
 ll2xy_lam_tan.o \
 ll2xy_merc.o \
 ll2xy_ps.o \
 locate.o \
 mapfac_lam.o \
 mapfac_merc.o \
 mapfac_ps.o \
 metgrid2ctm.o \
 metvars2ctm.o \
 myinterp.o \
 outclog.o \
 outcm3io.o \
 outglog.o \
 outgm3io.o \
 outncf.o \
 outncfbdy.o \
 outncfglobal.o \
 pblsup.o \
 ptemp.o \
 pvs.o \
 qsat.o \
 rdwrfem.o \
 rdfv3.o \
 readnml.o \
 resistcalc.o \
 setgriddefs.o \
 setup.o \
 setup_wrfem.o \
 setup_fv3.o \
 sfclayer.o \
 statflds.o \
 vertarys.o \
 vertnhy_wrf.o \
 vstamp.o \
 vtemp.o \
 wind.o \
 windrotation.o \
 wrfemopts.o \
 wrgdesc.o \
 xy2ll_lam.o \
 xy2ll_merc.o \

all:
	@$(MAKE) $(MODULES)
	@$(MAKE) $(MODEL)

$(MODEL): $(OBJS)
	$(FC) -o $(MODEL) $(FFLAGS) $(OBJS) $(MODULES) $(LIBS)

$(OBJS): $(MODULES)

.f90.o:
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F90.o:
	$(FC) $(FFLAGS) $(INCLUDES) -c $< $(DEFS)

clean:
	rm -f *.o *.mod *.il $(MODEL)
