PROG = pimd.exe
SHELL   = /bin/bash
OBJS    = $(SRCS:%.F90=%.o)
OBJSF77 = $(SRCSF77:%.f=%.o)
MPI  = True

#ifeq ($(HOSTNAME),wisteria)
#FC = mpifrtpx
#else
#FC = mpif90
#endif
#ifeq ($(HOSTNAME),ito)
#FC = mpiifort
#endif


#fcopt = -g -mcmodel=medium -O3 -no-gcc -traceback
#fcopt = -g -mcmodel=medium -O3 -no-gcc -traceback -cpp

ifeq ($(MPI),True)
FC = mpif90
fcopt = -O2 -pipe -D_mpi_
else
FC = gfortran
fcopt = -O2 -pipe
#fcopt = -g -O0 -pipe
endif

#fcopt = -O2 -cpp -pipe -Dmpi

# Debug
#fcopt = -g -check all -cpp

SRCS  = \
Parameter.F90                      \
utility.F90                        \
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
Print_Ham_Classical.F90            \
print_ini.F90                      \
Restart.F90                        \
Remove_TnR_All.F90                 \
Virial_Estimator.F90               \
Set_Deallocate.F90                 \
Force_New_MPI_tk.F90               \
set_pallarel.F90                   \
Set_Gaussian_MPI_tk.F90            \
Force_Gaussian_MPI_tk.F90          \
Print_Ham_tk2.F90                  \
Force_model_Morse.F90              \
Force_model_DoubleWell.F90         \
Set_siesta.F90                     \
Force_VASP_MPI.F90                 \
print_result_qm.F90                \
calc_umbrella.F90                  \
force_siesta.F90                   \
Set_mopac.F90                      \
Set_VASP.F90                       \
neural_network.F90                 \
force_water.F90                    \
exit.F90                           \

#Init_Send_Recv_MPI_tk.F90          \
#Init_Recv_Send_MPI_tk.F90          \
#Start_Recv_Send_MPI_tk.F90         \
#Start_Send_Recv_MPI_tk.F90         \


SRCSF77 =  \
dlarnv.f   \
dlaruv.f   \

# -- main program --
#all : $(dis)
#all : $(PROG) $(vir) $(dis) $(bl) $(dist)


all: $(PROG)

$(PROG): $(OBJS) $(OBJSF77)
	$(FC) $(fcopt) $(OBJS) $(OBJSF77) -o $(PROG)


# install: $(OBJS)
install: $(PROG)
#	$(FC) $(fcopt) $(OBJS)  -o $(PROG)
	cp $(PROG) ../bin/pimd.exe
#	cp $(PROG) ../bin/pimd_muon
#	cp $(PROG) ../bin/$(PROG)


# -- for F90 program --
.SUFFIXES: .F90 .o 
.F90.o: 
	@echo "<< Compiling >>" $<
	$(FC) $(fcopt) -c -o $*.o $*.F90
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


