##########################################
# CODE DIRECTORIES AND FILES
##########################################
mkfile_path := $(abspath $(firstword $(MAKEFILE_LIST)))
mkfile_dir := $(dir $(mkfile_path))
BIN_DIR := ./bin
SRC_DIR := ./src
LIB_DIR := ./lib
BUILD_DIR = ./obj
LIBS := mod_constants.f90 \
	mod_misc.f90 \
	mod_misc_maths.f90 \
	mod_misc_linalg.f90 \
	mod_rw_geom.f90 \
	mod_rw_vasprun.f90 \
	mod_edit_geom.f90 \
	mod_elements.f90 \
	mod_ml.f90 \
	mod_evolver.f90 \
	mod_evaluator.f90 \
	mod_atom_adder.f90 \
	mod_read_structures.f90
OBJS := $(addprefix $(LIB_DIR)/,$(LIBS))
#$(info VAR is $(OBJS))
SRCS := inputs.f90 \
	generator.f90 \
	main.f90
SRCS := $(OBJS) $(SRCS)
OBJS := $(addprefix $(SRC_DIR)/,$(SRCS))


##########################################
# COMPILER CHOICE SECTION
##########################################
FFLAGS = -O2
#PPFLAGS = -cpp
FC=gfortran
ifeq ($(FC),ifort)
	MPIFLAG = -qopenmp
	MODULEFLAG = -module
	DEVFLAGS = -check all -warn #all
	DEBUGFLAGS = -check all -fpe0 -warn -tracekback -debug extended # -check bounds
else
	MPIFLAG = -fopenmp
	MODULEFLAG = -J
	DEVFLAGS = -g -fbacktrace -fcheck=all -fbounds-check #-g -static -ffpe-trap=invalid
	DEBUGFLAGS = -fbounds-check -Wall -Wno-maybe-uninitialized
	MEMFLAG = #-mcmodel=large
endif


##########################################
# LAPACK SECTION
##########################################
MKLROOT?="/usr/local/intel/parallel_studio_xe_2017/compilers_and_libraries_2017/linux/mkl/lib/intel64_lin"
LLAPACK = $(MKLROOT)/libmkl_lapack95_lp64.a \
	-Wl,--start-group \
	$(MKLROOT)/libmkl_intel_lp64.a \
	$(MKLROOT)/libmkl_sequential.a \
	$(MKLROOT)/libmkl_core.a \
	-Wl,--end-group \
	-lpthread

#$(MKLROOT)/libmkl_scalapack_lp64.a \
#$(MKLROOT)/libmkl_solver_lp64_sequential.a \

ATHENAROOT = $(HOME)/.local/athena
LATHENA = -I$(ATHENAROOT)/include -L$(ATHENAROOT)/lib -lathena

##########################################
# COMPILATION SECTION
##########################################
INSTALL_DIR?=$(HOME)/bin
NAME = raffle
programs = $(BIN_DIR)/$(NAME)

.PHONY: all debug install uninstall dev mpi clean

all: $(programs)

$(BIN_DIR):
	mkdir -p $@

$(BUILD_DIR):
	mkdir -p $@

$(BIN_DIR)/$(NAME): $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(MEMFLAG) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) $(LATHENA) -o $@

install: $(OBJS) | $(INSTALL_DIR) $(BUILD_DIR)
	$(FC) $(MEMFLAG) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) $(LATHENA) -o $(INSTALL_DIR)/$(NAME)

debug: $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(MEMFLAG) $(DEBUGFLAGS) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) $(LATHENA) -o $(programs)

dev: $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(MEMFLAG) $(DEVFLAGS) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) $(LATHENA) -o $(programs)

mpi_dev: $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(MEMFLAG) $(DEVFLAGS) $(MPIFLAG) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) $(LATHENA) -o $(programs)

mpi: $(OBJS) | $(BIN_DIR) $(BUILD_DIR)
	$(FC) $(MEMFLAG) $(MPIFLAG) $(MODULEFLAG) $(BUILD_DIR) $(OBJS) $(LATHENA) -o $(programs)

clean: $(BUILD_DIR) $(BIN_DIR)
	rm -r $(BUILD_DIR)/ $(BIN_DIR)/

uninstall: $(INSTALL_DIR)/$(NAME)
	rm $(INSTALL_DIR)/$(NAME)
