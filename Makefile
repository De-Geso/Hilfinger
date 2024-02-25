# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Project name

# Configuration settings
FC := gfortran
LD := $(FC)
FFLAGS := -O3
LDLIBS := -lfftw3
LDFLAGS := -I/usr/include/
RM := rm -f

# List of all source files
PATHIN := src/
BINDIR:= bin/
NAMES := mrna_gene.f90 \
	kind_parameters.f90 \
	init_mrna_gene.f90 \
	randf.f90

# Create lists of the build artefacts in this project
SRCS := $(addprefix $(PATHIN), $(NAMES))
OBJS := $(addsuffix .o, $(SRCS))

# Declare all public targets
.PHONY: all clean

all: mrna_gene

mrna_gene: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(OBJS): %.o: %
	$(FC) $(LDFLAGS) -c -J$(BINDIR) -o $@ $< $(LDLIBS)

# define dependencies between object files
src/mrna_gene.f90.o: src/kind_parameters.f90.o src/init_mrna_gene.f90.o src/randf.f90.o

src/init_mrna_gene.f90.o: src/kind_parameters.f90.o

src/randf.f90.o: src/kind_parameters.f90.o

# rebuild all object files in case this Makefile changes
$(OBJS): $(MAKEFILE_LIST)

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS)) $(wildcard *.mod)
