# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Project name

# Configuration settings
FC := gfortran
LD := $(FC)
FFLAGS := -O3 -fopenmp
LDLIBS :=
LDFLAGS := -I/usr/include/ -O3 -fopenmp
RM := rm -f

# List of all source files
PATHIN := src/
BINDIR:= bin/
NAMES := mrna_gene.f90 \
	mrna_protein_system_parameters.f90 \
	kind_parameters.f90 \
	init_mrna_gene.f90 \
	randf.f90 \
	utilities.f90

# Create lists of the build artefacts in this project
SRCS := $(addprefix $(PATHIN), $(NAMES))
OBJS := $(addsuffix .o, $(SRCS))

# Declare all public targets
.PHONY: all clean

all: mrna_gene

$(OBJS): %.o: %
	$(FC) $(LDFLAGS) -c -J$(BINDIR) -o $@ $< $(LDLIBS)

mrna_gene: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)


# define dependencies between object files
src/mrna_gene.f90.o: src/kind_parameters.f90.o src/init_mrna_gene.f90.o src/randf.f90.o src/utilities.f90.o

src/mrna_protein_system_parameters.f90.o: src/kind_parameters.f90.o

src/init_mrna_gene.f90.o: src/kind_parameters.f90.o

src/randf.f90.o: src/kind_parameters.f90.o

src/utilities.f90.o: src/kind_parameters.f90.o

# rebuild all object files in case this Makefile changes
$(OBJS): $(MAKEFILE_LIST)

run: mrna_gene
	./mrna_gene

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS)) $(wildcard *.mod)
