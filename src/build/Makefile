.SUFFIXES: .o .F90

F90    = gfortran
LIBS   = -L$(NETCDF_DIR) -lnetcdf -lnetcdff -L$(JSONFORTRAN_LIB) -ljsonfortran -L$(NC4FORTRAN_LIB) -lnc4fortran
INCLUDE_MODULES = -I$(NETCDF_DIR) -I$(JSONFORTRAN_INC) -I$(NC4FORTRAN_INC) 

ifeq ($(FC),lf95)
  FFLAGS = --g
else
  FFLAGS = -g
endif

ifeq ($(FC),pgf90)
  FFLAGS += -Mnosave -Ktrap=fp
endif

ifeq ($(FC),gfortran)
  FFLAGS += -ffree-line-length-none -fcheck=bounds,do,pointer -ffpe-trap=zero,overflow,invalid -fbacktrace
endif

ifeq ($(FC),lf95)
  FFLAGS += --nap --chk --pca --nsav --trace --trap --wide
endif

ifeq ($(FC),ifort)
  FFLAGS += -fpe0 -ftrapuv
endif

EXEC = photo_rate_demo

OBJS = constants.o string.o assert.o convert.o environment.o iterator.o config.o \
       wavelength.o netcdf_util.o photo_utils.o abstract.cross_section.type.o \
       base.cross_section.type.o  hno3-oh_no2.cross_section.type.o n2o5-no2_no3.cross_section.type.o cl2-cl_cl.cross_section.type.o \
       abstract.quantum_yield.type.o base.quantum_yield.type.o \
       ch3cho-ch3_hco.cross_section.type.o ch3cho-ch3_hco.quantum_yield.type.o \
       cross_section_factory.o quantum_yield_factory.o photo_kinetics.o photo_driver.o


VPATH = ../core:../project:../project/cross_section_types:../project/quantum_yield_types

##dependencies
#$(OBJECTS) : makefile

${EXEC} :       ${OBJS}
		${F90} -o $@ ${OBJS} ${LIBS} $(CORE_SRC)

.F90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

cleanup:
	 	rm -f ${OBJS} *.mod *.smod

clean:
		rm -f core ${EXEC} ${OBJS} *.mod *.smod
