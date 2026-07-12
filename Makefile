SHELL   = /bin/bash
OBJS    = $(SRCS:%.F90=%.o)
OBJS   += $(OBJ_LAMMPS)
DIR_LAMMPS = LAMMPS
DIR_MACE = MACE
DIR_PFP = PFP
MACE ?= True
PFP ?= False   # MACE と違いデフォルト無効。有効化は明示的な PFP=True のみ (D1 の安全側フェイル)
# MACE ?= False
# GFLAG   =  -JModule -IModule
# MPI  = True

ifeq ($(MACE),True)
SRC_MACE = $(DIR_MACE)/fortran_mace_isoc.F90 Force_MACE.F90
OBJ_MACE = $(SRC_MACE:%.F90=%.o)
DFLAG += -D_MACE_
else
SRC_MACE =
OBJ_MACE =
endif

ifeq ($(PFP),True)
SRC_PFP = $(DIR_PFP)/fortran_pfp_isoc.F90 Force_PFP.F90
OBJ_PFP = $(SRC_PFP:%.F90=%.o)
DFLAG += -D_PFP_
else
SRC_PFP =
OBJ_PFP =
endif

# MACE か PFP のどちらかが有効なら Python 埋め込みビルドに必要なフラグを設定する
ifeq ($(MACE),True)
NEED_PYTHON = True
endif
ifeq ($(PFP),True)
NEED_PYTHON = True
endif

ifeq ($(NEED_PYTHON),True)
PYINC = $(shell python3-config --includes)
PYLIB = $(shell python3-config --embed --ldflags)
else
PYINC =
PYLIB =
endif


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
  PROG = pimd_mpi.exe
  else
  FC     = gfortran
  PROG = pimd.exe
  endif
# ===== Other =====
#endif


SRCS  = \
Parameter.F90                      \
utility.F90                        \
element_table.F90                  \
mod_model_force.F90                \
mod_lbfgs.F90                      \
mpi_module.F90                     \
main.F90                           \
Broad.F90                          \
read_input.F90                     \
Set_Allocate.F90                   \
$(SRC_MACE)                        \
$(SRC_PFP)                         \
Check_Inp.F90                      \
Calc_Constant.F90                  \
Setup_time_mass.F90                \
Normal_Mode.F90                    \
Init_Mass.F90                      \
NM_Position.F90                    \
Init_Velocity.F90                  \
Init_Bath.F90                      \
Temp_ctr.F90                       \
Kinetic_Energy.F90                 \
NM_Trans.F90                       \
Getforce_Ref.F90                   \
Nhc_Integrate_Cent.F90             \
Ham_Temp.F90                       \
Ham_Temp_Classical.F90             \
Nhc_Integrate.F90                  \
update_pos_vel.F90                 \
Force_MOPAC_MPI.F90                \
Force_Classical.F90                \
PI_NEW_MPI.F90                     \
Classical.F90                      \
print_ini.F90                      \
Restart.F90                        \
remove_trans_rot.F90               \
Virial_Estimator.F90               \
Barostat.F90                       \
Set_Deallocate.F90                 \
Force_New_MPI.F90                  \
set_pallarel.F90                   \
Set_Gaussian_MPI_tk.F90            \
Force_Gaussian.F90                 \
print_ham.F90                      \
Force_model_Morse.F90              \
Force_model_LJ_PBC.F90             \
Force_VASP_MPI.F90                 \
print_result.F90                   \
force_siesta.F90                   \
Set_mopac.F90                      \
Set_VASP.F90                       \
set_Iforce.F90                     \
neural_network.F90                 \
force_water.F90                    \
constrain.F90                      \
PIHMC.F90                          \
exit.F90                           \

# LammpsInterface.F90                \
# LammpsCalculator.F90               \
# force_LAMMPS.F90                   \
# dlarnv.f   \
# dlaruv.f   \



all: $(PROG)
	@echo -e '\e[34m Normal termination!!!\e[m\n'

$(PROG): $(OBJS)
ifeq ($(HOSTNAME),genkai)
	$(FC) $(LIBS) $(OBJS) -o $(PROG) $(PYLIB)
else
	$(FC) $(FCOPT) $(OBJS) -o $(PROG) $(PYLIB)
endif


install: $(PROG)
	cp $(PROG) ../bin/$(PROG)


# -- for F90 program --
#.SUFFIXES: .F90 .o
#.F90.o:
%.o:%.F90
	@echo "<< Compiling >>" $<
ifeq ($(HOSTNAME),genkai)
	$(FC) $(FCOPT) $(DFLAG) $(INCS) $(PYINC) -c -o $*.o $*.F90
else
	$(FC) $(FCOPT) $(DFLAG) $(PYINC) -c -o $*.o $*.F90
endif
	@echo

%.o:%.f
	@echo "<< Compiling >>" $<
	$(FC) $(FCOPT) -c -o $*.o $*.f
	@echo

clean:
	rm -rf *.o $(PROG) *.mod *exe
	rm -rf $(DIR_LAMMPS)/*.o
	rm -rf $(DIR_MACE)/*.o
	rm -rf $(DIR_PFP)/*.o
#	del *.o $(program) *.mod 


re: clean all
