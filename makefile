#-*- mode: makefile;-*-  
CC=gcc
FC=gfortran
LD=gfortran

ifeq ($(loc),lap)
INC = #-I/home/vishrut.garg/tecplot/include -I/home/chris/usr/bin/openblas/include 
LIB = #-L/home/vishrut.garg/tecplot/bin -L/home/chris/usr/bin/openblas/lib -lopenblas

else

INC = -I/opt/tecplot/current/include #-I/home/atom/a/anthonc/usr/bin/openblas/include
LIB = -L/opt/tecplot/current/bin # -ltecio -lstdc++
#-L/home/atom/a/anthonc/usr/bin/openblas/lib -lopenblas

endif


#pgf90 compiler, options include: debug
ifeq ($(c),pg)
CC=pgcc
FC=pgf90

ifeq ($(opt),debug)
CFLAGS= -Mfree -C -g -Mchkfpstk -Mchkptr -traceback
FFLAGS= $(INC) -module Modules/ -Mfree -C -g -Mchkfpstk -Mchkptr -traceback
LDFLAGS= -module Modules/ -Mfree -C -g -Mchkfpstk -Mchkptr -traceback

else
CFLAGS= -Mfree -O2  -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mvect=sse -Mcache_align -Mflushz -Mpre
FFLAGS= -module Modules/ -Mfree -fast -Mvect=nosse
LDFLAGS= -module Modules/ -Mfree -O2 -Munroll=c:1 -Mnoframe -Mlre -Mautoinline -Mvect=sse -Mcache_align -Mflushz -Mpre
#CFLAGS= -Mfree -fast 
#FFLAGS= -Mfree $(INC) -fast
#LDFLAGS= -Mfree -fast
endif


#ifort compiler, options include: debug
else ifeq ($(c),intel)
CC=icc
FC=ifort

ifeq ($(opt),debug)
CFLAGS=-O3 -w -ipo -xhost 
FFLAGS= $(INC) -O3 -qopenmp -ipo -xhost -FR -u -module Modules -g -traceback -heap-arrays
LDFLAGS= -O3 -qopenmp -ipo -xhost -module Modules -g -traceback -heap-arrays

else
CFLAGS=-O3 -w -ipo -xhost 
FFLAGS= $(INC) -O3 -qopenmp -ipo -xhost -FR -u -module Modules -heap-arrays
LDFLAGS= -O3 -qopenmp -ipo -xhost -module Modules -heap-arrays
endif


#gfortran compiler (no c=), options include: debug (slow through), debug2 (fast debug), prof (for profiling), fast (maybe faster but more dangerous), no_omp (no openmp)
else
CC=gcc
FC=gfortran
LD=gfortran

ifeq ($(opt),debug)
CFLAGS=-pg
FFLAGS = $(INC) -JModules -fcray-pointer -fimplicit-none -ffree-form -fbackslash -pg -g -fbounds-check -Wconversion -ffpe-trap=invalid,zero,overflow #-Wconversion-extra
LDFLAGS = -JModules -pg -fbounds-check -fopenmp -ffpe-trap=invalid,zero,overflow

else ifeq ($(opt),debug2)
CFLAGS=-O2 -fopenmp 
FFLAGS=$(INC) -O3 -g -JModules -fcray-pointer -ffree-form -fopenmp -march=native -funroll-loops -ffpe-trap=invalid,zero,overflow  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-arrays
LDFLAGS= -O3 -g -JModules -fopenmp -ffree-form -march=native -funroll-loops -ffpe-trap=invalid,zero,overflow  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-array

else ifeq ($(opt),no_omp)
CFLAGS=-O2
FFLAGS=$(INC) -O3 -g -JModules -fcray-pointer -ffree-form -march=native -funroll-loops -ffpe-trap=invalid,zero,overflow  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-arrays
LDFLAGS= -O3 -g -JModules -fopenmp -ffree-form -march=native -funroll-loops -ffpe-trap=invalid,zero,overflow  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-array

else ifeq ($(opt),prof)
CFLAGS=-O2 -fopenmp 
FFLAGS=$(INC) -O3 -pg -JModules -fcray-pointer -ffree-form -fopenmp -march=native -funroll-loops  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-arrays
LDFLAGS= -O3 -pg -JModules -fopenmp -ffree-form -march=native -funroll-loops  # -ffpe-trap=invalid,zero,overflow #-ffast-math -fno-protect-parens -flto -fstack-array

else ifeq ($(opt),fast)
CFLAGS=-O2 -fopenmp
FFLAGS=$(INC) -O3 -JModules -fcray-pointer -ffree-form -fopenmp -march=native -funroll-loops -ffast-math -fno-protect-parens -flto
LDFLAGS= -O3 -JModules -fopenmp -ffree-form -march=native -funroll-loops -ffast-math -fno-protect-parens -flto

else
CFLAGS=-O2 -fopenmp
FFLAGS=$(INC) -O3 -JModules -fcray-pointer -ffree-form -fopenmp -march=native -funroll-loops
LDFLAGS= -O3 -JModules -fopenmp -ffree-form -march=native -funroll-loops
endif


endif

SRCDIR=Sources
OBJDIR=Objects
MODDIR=Modules

#Add sources as $(SRCDIR)/file.f or $(SRCDIR)/file.c
CSOURCES= 
FSOURCES=  $(SRCDIR)/AAAkind.f90 $(SRCDIR)/AAdata.f90 $(SRCDIR)/AAdata_local.f90 $(SRCDIR)/A_NOP_mod.f90 $(SRCDIR)/Amultifront_general.f90 $(SRCDIR)/basis_function.f90 $(SRCDIR)/main.f90 $(SRCDIR)/BC.f90 $(SRCDIR)/graph.f90 $(SRCDIR)/flux.f90 $(SRCDIR)/initial_condition.f90 $(SRCDIR)/sj_VI.f90 $(SRCDIR)/sj_SI.f90 $(SRCDIR)/reverse_sj_part.f90 $(SRCDIR)/split_sol.f90 $(SRCDIR)/values_in_an_element.f90 $(SRCDIR)/values_in_sj.f90 $(SRCDIR)/sf.f90 $(SRCDIR)/L2_error.f90 $(SRCDIR)/prediction_and_preparation.f90 $(SRCDIR)/assemble.f90 $(SRCDIR)/jacobian_check.f90 $(SRCDIR)/jacobian_whole.f90 $(SRCDIR)/fsize.f90 $(SRCDIR)/parameter.f90 $(SRCDIR)/special_points.f90 $(SRCDIR)/var_cal.f90  $(SRCDIR)/newton_interation.f90 $(SRCDIR)/initialization.f90 $(SRCDIR)/flag_mesh.f90 $(SRCDIR)/data_folder.f90 $(SRCDIR)/radial_accum.f90 #$(SRCDIR)/drop_volume.f90

FOBJECTS= $(FSOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
COBJECTS= $(CSOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
COBJECTSCMD=$(CSOURCESCMD:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

EXECUTABLE= output

all: $(CSOURCES) $(FSOURCES) $(EXECUTABLE)

cmd: $(FSOURCES) $(CSOURCESCMD) $(EXECUTABLECMD) 

$(EXECUTABLE): $(FOBJECTS) $(COBJECTS)  
	$(FC) $(LDFLAGS) $(FOBJECTS) $(COBJECTS) -o $@ $(LIB)

$(FOBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(COBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *dat #$(MODDIR)/*mod $(EXECUTABLE) $(FOBJECTS) $(COBJECTS)

wipe:
	rm $(MODDIR)/*mod $(EXECUTABLE) $(FOBJECTS) $(COBJECTS) 

wipe_total:
	rm *dat $(MODDIR)/*mod $(EXECUTABLE) $(FOBJECTS) $(COBJECTS)
