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

PROGS := mrna_gene.f90 \
	mrna_protein_feedback.f90 \
	mrna_protein_false_rates.f90 \
	single_path_3D_discrete.f90 \
	simplex_cascade_2_discrete.f90 \
	simplex_cascade_3_discrete.f90 \
	acov_discrete_cascade_2D.f90 \
	acov_cascade_2D.f90

MODS := kind_parameters.f90 \
	init_mrna_gene.f90 \
	mrna_protein_system_parameters.f90 \
	randf.f90 \
	utilities.f90 \
	stochastics.f90 \
	nr_minmax.f90

# Create lists of the build artefacts in this project
MODSRCS := $(addprefix $(PATHIN), $(MODS))
MODOBJS := $(addsuffix .o, $(MODSRCS))

PROGSRCS := $(addprefix $(PATHIN), $(PROGS))
PROGOBJS := $(addsuffix .o, $(PROGSRCS))

# Declare all public targets
.PHONY: all clean

all: mrna_protein_false_rates single_path_3D_discrete simplex_cascade_2_discrete simplex_cascade_3_discrete acov_discrete_cascade_2D acov_cascade_2D# mrna_gene # mrna_protein_feedback

$(MODOBJS): %.o: %
	$(FC) $(LDFLAGS) -c -J$(BINDIR) -o $@ $< $(LDLIBS)

$(PROGOBJS): %.o: %
	$(FC) $(LDFLAGS) -c -J$(BINDIR) -o $@ $< $(LDLIBS)


single_path_3D_discrete: $(MODOBJS) src/single_path_3D_discrete.f90.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

mrna_gene: $(MODOBJS) src/mrna_gene.f90.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

mrna_protein_feedback: $(MODOBJS) src/mrna_protein_feedback.f90.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

mrna_protein_false_rates: $(MODOBJS) src/mrna_protein_false_rates.f90.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

simplex_cascade_2_discrete: $(MODOBJS) src/simplex_cascade_2_discrete.f90.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

simplex_cascade_3_discrete: $(MODOBJS) src/simplex_cascade_3_discrete.f90.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

acov_discrete_cascade_2D: $(MODOBJS) src/acov_discrete_cascade_2D.f90.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

acov_cascade_2D: $(MODOBJS) src/acov_cascade_2D.f90.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

	# Define dependencies between object files
# Programs
src/single_path_3D_discrete.f90.o: src/kind_parameters.f90.o src/stochastics.f90.o src/randf.f90.o

src/mrna_gene.f90.o: src/kind_parameters.f90.o src/init_mrna_gene.f90.o src/mrna_protein_system_parameters.f90.o src/randf.f90.o src/utilities.f90.o

src/mrna_protein_feedback.f90.o: src/kind_parameters.f90.o src/mrna_protein_system_parameters.f90.o src/randf.f90.o src/utilities.f90.o

src/mrna_protein_false_rates.f90.o: src/kind_parameters.f90.o src/mrna_protein_system_parameters.f90.o src/randf.f90.o src/utilities.f90.o

src/simplex_cascade_2_discrete.f90.o: src/kind_parameters.f90.o src/mrna_protein_system_parameters.f90.o src/randf.f90.o src/utilities.f90.o src/nr_minmax.f90.o

src/simplex_cascade_3_discrete.f90.o: src/kind_parameters.f90.o src/mrna_protein_system_parameters.f90.o src/randf.f90.o src/utilities.f90.o src/nr_minmax.f90.o

src/acov_discrete_cascade_2D.f90.o: src/kind_parameters.f90.o src/randf.f90.o src/utilities.f90.o src/stochastics.f90.o

src/acov_cascade_2D.f90.o: src/kind_parameters.f90.o src/randf.f90.o src/utilities.f90.o src/stochastics.f90.o

# Modules
src/kind_parameters.f90.o:

src/stochastics.f90.o: src/kind_parameters.f90.o

src/mrna_protein_system_parameters.f90.o: src/kind_parameters.f90.o

src/init_mrna_gene.f90.o: src/kind_parameters.f90.o

src/randf.f90.o: src/kind_parameters.f90.o

src/utilities.f90.o: src/kind_parameters.f90.o

src/nr_minmax.f90.o: src/kind_parameters.f90.o

# Run commands
run_mrna_gene: mrna_gene
	date
	./mrna_gene

run_mrna_protein_feedback: mrna_protein_feedback
	date
	./mrna_protein_feedback

run_mrna_protein_false_rates: mrna_protein_false_rates
	date
	@if [ -z "$(l1)" ] || [ -z "$(l2)" ]; then \
		echo "Usage: make run_mrna_protein_false_rates l1=<value> l2=<value>"; \
	else \
		./mrna_protein_false_rates $(l1) $(l2); \
	fi

# rebuild all object files in case this Makefile changes
$(OBJS): $(MAKEFILE_LIST)

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(MODOBJS), $(PROGOBJS)) $(wildcard *.mod)
