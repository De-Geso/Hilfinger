MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Define color codes
GREEN = '\033[0;32m'
RED = '\033[0;31m'
YELLOW = '\033[0;33m'
RESET = '\033[0m'

# Remove
RM := rm -f

# Directories
SRC := src
MODDIR := ../common
BINDIR := bin

# Compiler and flags
FC := gfortran
FFLAGS := -O3 -fopenmp -I$(MODDIR) -J$(BINDIR)
LD := $(FC)
LDFLAGS := -O3 -fopenmp
LDLIBS :=

# Sources
MODS := kind_parameters.f90 \
	randf.f90 \
	stochastics.f90 \
	utilities.f90

PROGS := inversion_cascade_1.f90 \
	read_binary.f90 \
	single_path_3D.f90 \
	single_path_3D_discrete.f90 \
        cascade_2_benchmark.f90 \
        trace_cascade_2.f90

MODOBJS := $(MODS:%.f90=$(MODDIR)/%.o)
PROGOBJS := $(PROGS:%.f90=$(SRC)/%.o)
TARGETS := $(PROGS:.f90=)

.PHONY: all clean

# Rule for compiling modules
$(MODDIR)/%.o: $(MODDIR)/%.f90
	@echo -e $(YELLOW)Compiling $@$(RESET)
	$(FC) $(FFLAGS) -c $< -o $@

# Rule for compiling program sources
$(SRC)/%.o: $(SRC)/%.f90 $(MODOBJS)
	@echo -e $(YELLOW)Compiling $@$(RESET)
	$(FC) $(FFLAGS) -c $< -o $@

# Linking programs
%: $(MODOBJS) $(SRC)/%.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# All targets
all: $(TARGETS)

inversion_cascade_1: $(MODOBJS) src/inversion_cascade_1.o
	$(FC) $(FFLAGS) -o $@ $^

read_binary: $(MODOBJS) $(SRC)/read_binary.o
	$(FC) $(FFLAGS) -o $@ $^

# Program dependencies (ensure all the necessary .o files are there before trying to link)
$(SRC)/read_binary.o:
$(SRC)/inversion_cascade_1.o: $(MODOBJS)

# Module dependencies
$(MODDIR)/randf.o: $(MODDIR)/kind_parameters.o
$(MODDIR)/stochastics.o: $(MODDIR)/randf.o $(MODDIR)/kind_parameters.o
$(MODDIR)/utilities.o: $(MODDIR)/kind_parameters.o

# Clean up build files
clean:
	$(RM) $(MODOBJS) $(PROGOBJS) $(TARGETS) $(BINDIR)/*.mod
