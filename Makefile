# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Project name

# Configuration settings
FC := gfortran
LD := $(FC)
FFLAGS := -O2
LDFLAGS := $(FCFLAGS)
RM := rm -f

# List of all source files
PATHIN := src/
PATHOUT:= bin/
NAMES := mrna_gene.f90 \
	kind_parameters.f90
SRCS := $(addprefix $(PATHIN), $(NAMES))

# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))

# Declare all public targets
.PHONY: all clean
all: mrna_gene

mrna_gene: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^

$(OBJS): %.o: %
	$(FC) -c -o $@ $<

# define dependencies between object files
src/mrna_gene.f90.o: src/kind_parameters.f90.o

# rebuild all object files in case this Makefile changes
#$(OBJS): $(MAKEFILE_LIST)

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(filter %.exe, $(TEST_EXE)) $(LIB) $(wildcard *.mod)
