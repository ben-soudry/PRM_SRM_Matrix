LDLIBS= -fopenmp -pthread

prmSrmMatrix: prm_srm_matrix.o parse_results_tsv.o io.o ppseq.o theoretical_util.o
	g++ -O3 -std=c++17 prm_srm_matrix.o io.o ppseq.o theoretical_util.o parse_results_tsv.o -o prmSrmMatrix $(LDLIBS)

prm_srm_matrix.o: prm_srm_matrix.cpp prm_srm_matrix.hpp
	g++ -O3 -std=c++17 -c prm_srm_matrix.cpp

parse_results_tsv.o: parse_results_tsv.cpp parse_results_tsv.hpp
	g++ -O3 -std=c++17 -c parse_results_tsv.cpp

io.o: io.cpp io.hpp
	g++ -O3 -std=c++17 -c io.cpp

ppseq.o: ppseq.cpp ppseq.hpp
	g++ -O3 -std=c++17 -c ppseq.cpp

theoretical_util.o: theoretical_util.cpp theoretical_util.hpp
	g++ -O3 -std=c++17 -c theoretical_util.cpp
