
dada2ms: *.cc *.h
	g++ -I/scr2/mweastwood/usr/include -L/scr2/mweastwood/usr/lib -O3 -std=c++11 -o dada2ms *.cc -lcasa_casa -lcasa_measures -lcasa_ms -lcasa_tables -lcasa_scimath -lcasa_scimath_f -lboost_program_options

clean:
	rm dada2ms
