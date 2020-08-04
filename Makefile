.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean))) #############

CPPFLAGS := -std=c++20 -Iinclude
CXXFLAGS := -Wall -O3 -flto -fmax-errors=3
# CXXFLAGS := -Wall -O0 -g -fmax-errors=3

# generate .d files during compilation
DEPFLAGS = -MT $@ -MMD -MP -MF .build/$*.d

ROOT_CPPFLAGS := $(shell root-config --cflags \
  | sed 's/-std=[^ ]*//g;s/-I/-isystem /g')
ROOT_LDFLAGS  := $(shell root-config --ldflags) \
  -Wl,-rpath,$(shell root-config --libdir)
ROOT_LDLIBS   := $(shell root-config --libs)

FJ_PREFIX   := $(shell fastjet-config --prefix)
FJ_CPPFLAGS := -I$(FJ_PREFIX)/include
FJ_LDLIBS   := -L$(FJ_PREFIX)/lib -Wl,-rpath=$(FJ_PREFIX)/lib -lfastjet

LHAPDF_PREFIX   := $(shell lhapdf-config --prefix)
LHAPDF_CPPFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LDLIBS   := $(shell lhapdf-config --libs) -Wl,-rpath=$(LHAPDF_PREFIX)/lib

#####################################################################

all: bin/hist

#####################################################################

C_hist := $(ROOT_CPPFLAGS) $(FJ_CPPFLAGS) $(LHAPDF_CPPFLAGS) -I. -DNDEBUG
LF_hist := $(ROOT_LDFLAGS)
L_hist := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS) $(LHAPDF_LDLIBS)
.build/hist.o: .build/punch.hh
bin/hist: .build/reweighter.o .build/Higgs2diphoton.o

C_reweighter := $(ROOT_CPPFLAGS) $(LHAPDF_CPPFLAGS)
C_Higgs2diphoton := $(ROOT_CPPFLAGS)

#####################################################################

.PRECIOUS: .build/%.o

bin/%: .build/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(LF_$*) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

.build/%.o: src/%.cc
	@mkdir -pv $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $(C_$*) -c $(filter %.cc,$^) -o $@

.build/punch.hh: punchcards
	@mkdir -pv $(dir $@)
	@ls $< | sed 's/\.punch$$//;s/^.*/h_(&)/' > $@
	@printf '\nstatic constexpr const char* cards_names[] {\n' >> $@
	@sed -n 's/^h_(\(.*\))$$/  "\1",/p' $@ >> $@
	@printf '};\n' >> $@

-include $(shell find .build -type f -name '*.d' 2>/dev/null)

endif ###############################################################

clean:
	@rm -rfv bin lib .build

