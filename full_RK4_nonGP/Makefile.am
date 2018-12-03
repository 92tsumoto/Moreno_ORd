#INCLUDES = syspara.h

bin_PROGRAMS = ORd_Markov_INa 
#ORd_Markov_INa_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg
ORd_Markov_INa_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -vec-report2 -pg

ORd_Markov_INa_INCLUDE = -I./include
ORd_Markov_INa_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
#ORd_Markov_INa_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ./lib/xhplot.o

ORd_Markov_INa_SOURCES = ./src/ord.c ./src/val_constants.c ./src/mem.c ./src/input_para.c\
	./src/eular.c ./src/runge.c ./src/system.c ./src/currents.c ./src/exptable.c\
	./src/syspara.h
