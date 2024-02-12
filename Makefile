PROG = pimd.exe
SHELL   = /bin/bash
OBJS    = $(SRCS:%.f90=%.o)
OBJSF77 = $(SRCSF77:%.f=%.o)

ifeq ($(HOSTNAME),wisteria)
FC = mpifrtpx
else
FC = mpif90
endif

ifeq ($(HOSTNAME),ito)
FC = mpiifort
endif


#fcopt = -g -mcmodel=medium -O3 -no-gcc -traceback
#fcopt = -g -mcmodel=medium -O3 -no-gcc -traceback -cpp

fcopt = -O2 -cpp -pipe

# Debug
#fcopt = -g -check all -cpp

SRCS  = \
Parameter.f90                      \
utility.f90                        \
Main_MPI.f90                       \
Broad.f90                          \
read_input.f90                     \
Set_Allocate.f90                   \
Check_Inp.f90                      \
Calc_Constant.f90                  \
Setup_MPI.f90                      \
Setup.f90                          \
Setup_Classical.f90                \
Normal_Mode.f90                    \
Init_Mass.f90                      \
NM_Position.f90                    \
Init_Velocity.f90                  \
Init_Bath.f90                      \
Init_Bath_Classical.f90            \
Temp_ctr.f90                       \
Kinetic_Energy.f90                 \
NM_Trans.f90                       \
Getforce_Ref.f90                   \
Nhc_Integrate_Cent.f90             \
Ham_Temp.f90                       \
Ham_Temp_Classical.f90             \
Nhc_Integrate.f90                  \
Vupdate_Ref.f90                    \
Update.f90                         \
RandomG.f90                        \
Force_Gaussian_classical.f90       \
Force_MOPAC_MPI.f90                \
Force_Classical.f90                \
GasDev.f90                         \
PI_NEW_MPI.f90                     \
Classical.f90                      \
Print_Ham.f90                      \
Print_Ham_Classical.f90            \
print_ini.f90                      \
Restart.f90                        \
Remove_TnR_All.f90                 \
Virial_Estimator.f90               \
Set_Deallocate.f90                 \
Force_New_MPI_tk.f90               \
Set_etc_MPI_tk.f90                 \
Init_Send_Recv_MPI_tk.f90          \
Init_Recv_Send_MPI_tk.f90          \
Start_Recv_Send_MPI_tk.f90         \
Start_Send_Recv_MPI_tk.f90         \
Set_Gaussian_MPI_tk.f90            \
Force_Gaussian_MPI_tk.f90          \
Print_Ham_tk2.f90                  \
Force_model_Morse.f90              \
Force_model_DoubleWell.f90         \
Set_siesta.f90                     \
Force_VASP_MPI.f90                 \
print_result_qm.f90                \
print_result_cl.f90                \
calc_umbrella.f90                  \
force_siesta.f90                   \
Set_mopac.f90                      \
Set_VASP.f90                       \
neural_network.f90                 \
force_water.f90                    \
exit.f90                           \

# print_ini_cl.f90                   \
# print_ini_qm.f90                   \
# Unset_etc_MPI_tk.f90               \
# Set_Allocate_Classical.f90         \
# Set_Deallocate_Classical.f90       \
#nmtrans_force_r2ur.f90             \
#CMD_NEW_MPI.f90                    \
#RPMD_NEW_MPI.f90                   \
# Read_Inp.f90                       \
# Send_Recv_MPI_tk2.f90              \
# Recv_Send_MPI_tk2.f90              \
# Recv_Send_MPI_tk3.f90              \
# Force_DoubleHarmonic.f90           \


SRCSF77 =  \
dlarnv.f   \
dlaruv.f   \

# -- main program --
#all : $(dis)
#all : $(PROG) $(vir) $(dis) $(bl) $(dist)


all: $(PROG)

$(PROG): $(OBJS) $(OBJSF77)
	$(FC) $(fcopt) $(OBJS) $(OBJSF77) -o $(PROG)


install: $(OBJS)
	$(FC) $(fcopt) $(OBJS)  -o $(PROG)
	cp $(PROG) ../bin/pimd
#	cp $(PROG) ../bin/pimd_muon
#	cp $(PROG) ../bin/$(PROG)


# -- for f90 program --
.SUFFIXES: .f90 .o 
.f90.o: 
	@echo "<< Compiling >>" $<
	$(FC) $(fcopt) -c -o $*.o $*.f90
	@echo

.SUFFIXES: .f .o
.f.o: 
	@echo "<< Compiling >>" $<
	$(FC) $(fcopt) -c -o $*.o $*.f
	@echo

clean: 
	rm -rf *.o $(PROG) *.mod *exe
#	del *.o $(program) *.mod 


re: clean all


