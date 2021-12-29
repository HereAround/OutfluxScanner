uninstall:
	( rm -f counter_H1.o && rm -f counter_H2.o && rm -f new_counter.o)
	( rm -f counter_H1 && rm -f counter_H2 && rm -f new_counter)

unzip:
	( cd data_H1 && unzip fluxes_H1.zip )
	( cd data_H2 && unzip fluxes_H2_part1.zip && unzip fluxes_H2_part2.zip )

install: uninstall
	( g++ -std=gnu++11 -c -lboost_thread counter_H1.cpp && g++ -o counter_H1 counter_H1.o -lboost_thread -lpthread )
	( g++ -std=gnu++11 -c -lboost_thread counter_H2.cpp && g++ -o counter_H2 counter_H2.o -lboost_thread -lpthread )
	( g++ -std=gnu++11 -c -lboost_thread new_counter.cpp && g++ -o new_counter new_counter.o -lboost_thread -lpthread )

.PHONY: uninstall install
