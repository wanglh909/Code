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



ifeq ($(opt),debug)

CFLAGS=-pg
FFLAGS = $(INC) -JModules -fimplicit-none -pg -fbounds-check -Wconversion -fopenmp
LDFLAGS = -JModules -pg -fbounds-check -fopenmp

else

CFLAGS=-O2 -fopenmp 
FFLAGS=$(INC) -O3 -JModules -march=native -funroll-loops -fopenmp
LDFLAGS= -O3 -JModules -march=native -funroll-loops -fopenmp

endif

SRCDIR=Sources
OBJDIR=Objects
MODDIR=Modules

#Add sources as $(SRCDIR)/file.f or $(SRCDIR)/file.c
CSOURCES= 
FSOURCES=  $(SRCDIR)/AAAkind.f90 $(SRCDIR)/AAdata.f90 $(SRCDIR)/AAdata_local.f90 $(SRCDIR)/A_NOP_mod.f90 $(SRCDIR)/Amultifront.f90 $(SRCDIR)/basis_function.f90 $(SRCDIR)/main.f90 $(SRCDIR)/BC.f90 $(SRCDIR)/graph.f90 $(SRCDIR)/flux.f90 $(SRCDIR)/initial_condition.f90 $(SRCDIR)/sj_VI.f90 $(SRCDIR)/sj_SI.f90 $(SRCDIR)/reverse_sj_part.f90 $(SRCDIR)/split_sol.f90 $(SRCDIR)/values_in_an_element.f90 $(SRCDIR)/values_in_sj.f90 $(SRCDIR)/sf.f90 $(SRCDIR)/L2_error.f90 $(SRCDIR)/prediction_and_preparation.f90 $(SRCDIR)/assemble.f90 $(SRCDIR)/jacobian_check.f90 $(SRCDIR)/jacobian_whole.f90 $(SRCDIR)/fsize.f90 $(SRCDIR)/parameter.f90 $(SRCDIR)/var_cal.f90  $(SRCDIR)/newton_interation.f90 $(SRCDIR)/initialization.f90 $(SRCDIR)/flag_mesh.f90 $(SRCDIR)/data_folder.f90 #$(SRCDIR)/drop_volume.f90

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
