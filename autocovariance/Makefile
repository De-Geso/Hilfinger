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
STDDIR := ../common
BINDIR := bin

# Compiler and flags
FC := gfortran
FFLAGS := -O3 -fopenmp -I$(MODDIR) -I$(STDDIR) -J$(BINDIR)
LD := $(FC)
LDFLAGS := -O3 -fopenmp
LDLIBS :=

# Sources
MODS := kind_parameters.f90 \
	randf.f90 \
	stochastics.f90 \
	utilities.f90 \
	nr_minmax.f90

PROGS := single_path_3D_discrete.f90 \
	simplex_cascade_3_discrete.f90


MODOBJS := $(MODS:%.f90=$(MODDIR)/%.o)
STDOBJS := $(STDS:%.f90=$(STDDIR)/%.o)
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

single_path_3D_discrete: $(MODOBJS) $(SRC)/single_path_3D_discrete.o
	$(FC) $(FFLAGS) -o $@ $^

simplex_cascade_3_discrete: $(MODOBJS) $(SRC)/simplex_cascade_3_discrete.o
	$(FC) $(FFLAGS) -o $@ $^

$(SRC)/simplex_cascade_3_discrete.o: $(MODOBJS)
$(SRC)/single_path_3D_discrete.o: $(MODOBJS)


# Module dependencies
$(MODDIR)/randf.o: $(MODDIR)/kind_parameters.o
$(MODDIR)/stochastics.o: $(MODDIR)/randf.o $(MODDIR)/kind_parameters.o $(MODDIR)/utilities.o
$(MODDIR)/utilities.o: $(MODDIR)/kind_parameters.o
$(MODDIR)/nr_minmax.o: $(MODDIR)/kind_parameters.o

# Clean up build files
clean:
	$(RM) $(MODOBJS) $(STDOBJS) $(PROGOBJS) $(TARGETS) $(BINDIR)/*.mod
