program_main = pimd

#####  Intel Fortran Compiler
#fc          = ifort
##. . . . .   Debug
#fcopt       = -static -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback -c 
##. . . . .   Optimize 
#fcopt       = -i-static  -g  -fast -unroll 
#fcopt       = -i-static  -g  -xW 
#####
#####  Intel Fortran Compiler
fc = mpif90
#fcopt = -g -mcmodel=medium -O3 -no-gcc -traceback -cpp

# By Kuwahata
#fcopt = -g -mcmodel=medium -O3 -cpp
#fcopt = -g -shared-intel -O3 -cpp
# fcopt = -mcmodel=large -O3 -cpp -ffree-line-length-none -pipe
# fcopt = -mcmodel=large -O3 -cpp -ffree-line-length-none
fcopt = -O2 -cpp -ffree-line-length-none

# Debug
#fcopt = -g -check all -cpp

#fc = mpiifort
#fc1 = mpiifort
#fcopt = -g -mcmodel=medium -O2 -mkl -xAVX -DNODIPOLE
#fc = mpif90
#fc1 = mpif77
#fcopt = -g -mcmodel=medium -O2 -xSSE4.2 -DNODIPOLE
#fcopt = -g -mcmodel=medium -openmp -O3
#fcopt = -g -mcmodel=medium -openmp -no-gcc
#fcopt = -g -CB -mpitrace -inline-debug-info -d_lines -debug extended -debug-parameters all -check all -traceback -mcmodel=medium
##. . . . .   Debug
#fcopt       = -fimplicit-none -Wall -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbounds-check  
##. . . . .   Optimize 
#fcopt       = -fimplicit-none -march=core2 -mtune=core2 -O3 
#fcopt       = -O3 -static -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback -c 
#####


obj_main  =          \
global_variable.o    \
utility.o            \
communication.o      \
print_out.o          \
main.o               \
read_input.o         \
set_allocate.o       \
broad_parameters.o   \
check_input.o        \
setup_whole.o        \
random_generator.o   \
setup_indivi.o       \
init_mass.o          \
normal_mode.o        \
restart_read.o       \
print_ini.o          \
init_position.o      \
init_velocity.o      \
init_bath.o          \
get_force.o          \
calc_hamil.o         \
integ_nhc_cent.o     \
getforce_ref.o       \
#estimator.o          \
#setup_mpi.o          \
#normal_mode_matrix.o \
#normal_mode_trans.o  \
#remove_rotation.o   \
#tune_tempe.o        \

module = \
global_variable.mod  \
utility.mod           \
communication.mod  \
print_out.mod  \


%.mod : %.f90 %.o
	@true


#obj_sub = ../gamess_force/getforce_gamess.o  \
#          ../gamess_force/atom.o  \
#          ../gamess_force/utils.o \
#          ../gamess_force/gmscall.o  \
#          ../gamess_force/f_wrapper.o


#vir      =   Virial.go
#obj_vir  =   Virial.o   
#dis      =   Distance.go 
#obj_dis  =   Distance.o   
#bl       =   Block_Average.go 
#obj_bl   =   Block_Average.o
#dist     =   Distribution.go       
#obj_dist =   Distribution.o       
  

# -- main program --
#all : $(dis)
#all : $(program_main) $(vir) $(dis) $(bl) $(dist)
all: $(obj_main) 
#	$(fc) $(obj_main) -o $(program_main) 
	$(fc) $(fcopt) $(obj_main) -o $(program_main) 
#
### Virial Estimator for 2nd and 4th order Trotter expansion
#$(vir): $(obj_vir) 
#	$(fc) $(fcopt)  $(obj_vir)  -o   $(vir)  
### Distance Estimator for 2nd and 4th order Trotter expansion
#$(dis): $(obj_dis) 
#	$(fc) $(fcopt)  $(obj_dis)  -o   $(dis)  
### Block Average 
#$(bl): $(obj_bl) 
#	$(fc) $(fcopt)  $(obj_bl)   -o    $(bl)  
### Block Average 
#$(dist): $(obj_dist) 
#	$(fc) $(fcopt)  $(obj_dist) -o    $(dist)  
  

 
# -- for f90 program --
.SUFFIXES: .f90 .o 
.f90.o: 
	@echo
	@echo "<< Compiling >>" $<
	$(fc) $(fcopt) -c -o $*.o $*.f90


.SUFFIXES: .f .o
.f.o: 
	@echo
	@echo "<< Compiling >>" $<
	$(fc) $(fcopt) -c -o $*.o $*.f
clean: 
	rm -rf *.o $(program_main) *.mod 




#Parameter.o                      \
#Parameter_tk.o                    \
#Main_MPI.o                        \
#Broad.o                           \
#Read_Inp.o                       \
#Set_Allocate.o                   \
#Set_Allocate_Classical.o         \
#Set_Allocate_Shoot.o             \
#Check_Inp.o                      \
#Calc_Constant.o                  \
#Calc_Constant_Shoot.o            \
#Setup_MPI.o                      \
#Setup.o                          \
#Setup_Classical.o                \
#Normal_Mode.o                    \
#Init_Mass.o                      \
#Init_Mass_Shoot.o                \
#NM_Position.o                    \
#NM_Position_Shoot_MPI.o          \
#NM_Position_NewShot.o            \
#Init_Velocity.o                  \
#Init_Velocity_Shoot.o            \
#Init_Velocity_NewShot.o          \
#Init_Bath.o                      \
#Init_Bath_Classical.o            \
#Init_Bath_NewShot.o              \
#Temp_ctr.o                       \
#Temp_ctr_shoot.o                 \
#Temp_ctr_NewShot.o               \
#Kinetic_Energy.o                 \
#Kinetic_Energy_Shoot.o           \
#Kinetic_Energy_NewShot.o         \
#NM_Trans.o                       \
#Print_Out_second.o               \
#Print_Out_fourth.o               \
#Getforce_Ref.o                   \
#Getfnm.o                         \
#Nhc_Integrate_Cent.o             \


#.SUFFIXES: .F90 .o
#.F90.o:
#	@echo
#	@echo "<< Compiling >>" $<
#	$(fc) $(fcopt) -c -o $*.o $*.F90

