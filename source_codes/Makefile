# Copyright 2018, Gurobi Optimization, LLC

TWOUP    = ..
INC      = /home/felipe/Gurobi/gurobi810/linux64/include/
CC       = gcc
CPP      = g++
CARGS    = -m64 -g
CLIB     = -L/home/felipe/Gurobi/gurobi810/linux64/lib -lgurobi81
CPPLIB   = -L/home/felipe/Gurobi/gurobi810/linux64/lib -lgurobi_c++ -lgurobi81

pli-pa: pli-pa.cpp
	$(CPP) $(CARGS) -o $@ $< -I$(INC) $(CPPLIB) -lm

run_pli-pa: pli-pa
	./pli-pa $(args)

clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm kmp