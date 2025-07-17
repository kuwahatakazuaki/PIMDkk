PROG = pimd.exe
SHELL   = /bin/bash
OBJS    = $(SRCS:%.F90=%.o)
OBJSF77 = $(SRCSF77:%.f=%.o)
DEP     = $(OBJS:.o=.d)
SRCS   += $(SRC_LAMMPS)
OBJS   += $(OBJ_LAMMPS)
FCOPT  += $(DFLAG)
DIR_LAMMPS = LAMMPS
# GFLAG   =  -JModule -IModule
MPI  = True


#ifeq ($(HOSTNAME),genkai)
## ===== Kyushu Super Computer =====
#LMPROOT  = /home/pj24003139/ku40003238/bin
#INCS     = -I$(LMPROOT)/include
#LIBS     = -L$(LMPROOT)/lib64 -llammps_serial
#SRC_LAMMPS = $(DIR_LAMMPS)/LammpsInterface.F90 $(DIR_LAMMPS)/LammpsCalculator.F90 $(DIR_LAMMPS)/force_LAMMPS.F90
## SRC_LAMMPS = $(wildcard $(DIR_LAMMPS)/*.F90)
#OBJ_LAMMPS = $(SRC_LAMMPS:%.F90=%.o)
#DFLAG += -D_LAMMPS_
#FCOPT  = -O2 # -pipe -qmkl
#
#  ifeq ($(MPI),True)
#  FC     = mpiifort
#  DFLAG += -D_mpi_
#  else
#  FC     = ifort
#  endif
## ===== Kyushu Super Computer =====
#
#else

# ===== Other =====
FCOPT  = -O2 -pipe # -JModule
#FCOPT = -g -O0 -pipe
# FCOPT  = -O2 -pipe -JModule -MMD
  ifeq ($(MPI),True)
  FC     = mpif90
  DFLAG += -D_mpi_
  else
  FC     = gfortran
  endif
# ===== Other =====
#endif


SRCS  = \
Parameter.F90                      \
utility.F90                        \
mod_model_force.F90                \
mpi_module.F90                     \
Main_MPI.F90                       \
Broad.F90                          \
read_input.F90                     \
Set_Allocate.F90                   \
Check_Inp.F90                      \
Calc_Constant.F90                  \
Setup_time_mass.F90                \
Normal_Mode.F90                    \
Init_Mass.F90                      \
NM_Position.F90                    \
Init_Velocity.F90                  \
Init_Bath.F90                      \
Init_Bath_Classical.F90            \
Temp_ctr.F90                       \
Kinetic_Energy.F90                 \
NM_Trans.F90                       \
Getforce_Ref.F90                   \
Nhc_Integrate_Cent.F90             \
Ham_Temp.F90                       \
Ham_Temp_Classical.F90             \
Nhc_Integrate.F90                  \
Vupdate_Ref.F90                    \
Update.F90                         \
RandomG.F90                        \
Force_Gaussian_classical.F90       \
Force_MOPAC_MPI.F90                \
Force_Classical.F90                \
GasDev.F90                         \
PI_NEW_MPI.F90                     \
Classical.F90                      \
print_ini.F90                      \
Restart.F90                        \
Remove_TnR_All.F90                 \
Virial_Estimator.F90               \
Set_Deallocate.F90                 \
Force_New_MPI_tk.F90               \
set_pallarel.F90                   \
Set_Gaussian_MPI_tk.F90            \
Force_Gaussian_MPI_tk.F90          \
print_ham.F90                      \
Force_model_Morse.F90              \
Set_siesta.F90                     \
Force_VASP_MPI.F90                 \
print_result_qm.F90                \
calc_umbrella.F90                  \
force_siesta.F90                   \
Set_mopac.F90                      \
Set_VASP.F90                       \
neural_network.F90                 \
force_water.F90                    \
constrain.F90                      \
exit.F90                           \

# LammpsInterface.F90                \
# LammpsCalculator.F90               \
# force_LAMMPS.F90                   \


SRCSF77 =  \
dlarnv.f   \
dlaruv.f   \


all: $(PROG)
	@echo -e '\e[34m Noraml termination!!!\e[m\n'

$(PROG): $(OBJS) $(OBJSF77)
ifeq ($(HOSTNAME),genkai)
	$(FC) $(LIBS) $(OBJS) $(OBJSF77) -o $(PROG)
else
	$(FC) $(FCOPT) $(OBJS) $(OBJSF77) -o $(PROG)
endif


install: $(PROG)
	cp $(PROG) ../bin/$(PROG)


# -- for F90 program --
#.SUFFIXES: .F90 .o
#.F90.o:
%.o:%.F90
	@echo "<< Compiling >>" $<
ifeq ($(HOSTNAME),genkai)
	$(FC) $(FCOPT) $(DFLAG) $(INCS) -c -o $*.o $*.F90
else
	$(FC) $(FCOPT) $(DFLAG) -c -o $*.o $*.F90
endif
	@echo

# .SUFFIXES: .f .o
# .f.o:
%.o:%.f
	@echo "<< Compiling >>" $<
	$(FC) $(FCOPT) -c -o $*.o $*.f
	@echo

-include $(DEP)

clean: 
	rm -rf *.o $(PROG) *.mod *exe *.d
	rm -rf $(DIR_LAMMPS)/*.o
#	del *.o $(program) *.mod 


re: clean all


