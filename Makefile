CASACORE_INC_DIR := /opt/astro/casacore-2.0.3/include/casacore
CASACORE_LIB_DIR := /opt/astro/casacore-2.0.3/lib
OUT_EXE := dada2ms

dada2ms: *.cc *.h
	g++ -O3 -std=c++11 -I$(CASACORE_INC_DIR) -L$(CASACORE_LIB_DIR) -o $(OUT_EXE) *.cc -lcasa_casa -lcasa_measures -lcasa_ms -lcasa_tables -lcasa_scimath -lcasa_scimath_f -lboost_program_options

clean:
	rm dada2ms
