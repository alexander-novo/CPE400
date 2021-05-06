# Upgraded to c++17 to support the Filesystem library
CXXFLAGS     = -std=c++17 -fopenmp -O3 -D_GLIBCXX_PARALLEL
OBJDIR       = obj
DEPDIR       = $(OBJDIR)/.deps
# Flags which, when added to gcc/g++, will auto-generate dependency files
DEPFLAGS     = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d

# Function which takes a list of words and returns a list of unique words in that list
# https://stackoverflow.com/questions/16144115/makefile-remove-duplicate-words-without-sorting
uniq         = $(if $1,$(firstword $1) $(call uniq,$(filter-out $(firstword $1),$1)))

# Source files - add more to auto-compile into .o files
SOURCES      = src/main.cpp
INCLUDES     = -I include/
# Executable targets - add more to auto-make in default 'all' target
EXEC         = proj

SOURCEDIRS   = $(call uniq, $(dir $(SOURCES)))
OBJDIRS      = $(addprefix $(OBJDIR)/, $(SOURCEDIRS))
DEPDIRS      = $(addprefix $(DEPDIR)/, $(SOURCEDIRS))
DEPFILES     = $(SOURCES:%.cpp=$(DEPDIR)/%.d)

.PHONY: all exec clean report required
.SECONDARY:

# By default, make all executable targets and the outputs required for the homework
all: exec Report/report.pdf
exec: $(EXEC)

# Executable Targets
proj: $(OBJDIR)/src/main.o
	$(CXX) $(CXXFLAGS) $^ -o $@

### Experiment Outputs ###

# Figures needed for the report
report:

Report/report.pdf: Report/report.tex report
	latexmk -pdf -cd -use-make -silent -pdflatex='pdflatex -interaction=batchmode -synctex=1' $<

clean:
	rm -rf $(OBJDIR)
	rm -f $(EXEC)
	rm -rf out
	cd Report/; latexmk -c

# Generate .png images from .pgm images. Needed for report, since pdfLaTeX doesn't support .pgm images
%.png: %.pgm
	pnmtopng $< > $@

%.png: %.ppm
	pnmtopng $< > $@

# Auto-Build .cpp files into .o
$(OBJDIR)/%.o: %.cpp
$(OBJDIR)/%.o: %.cpp $(DEPDIR)/%.d | $(DEPDIRS) $(OBJDIRS)
	$(CXX) $(DEPFLAGS) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

# Make generated directories
$(DEPDIRS) $(OBJDIRS) out: ; @mkdir -p $@
$(DEPFILES):
include $(wildcard $(DEPFILES))