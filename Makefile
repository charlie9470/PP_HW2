CC = gcc
CXX = g++
LDLIBS = -lpng
CFLAGS = -lm -O3
hw2a: CFLAGS += -pthread
hw2b: CC = mpicc
hw2b: CXX = mpicxx
hw2b: CFLAGS += -fopenmp
CXXFLAGS = $(CFLAGS)
outname = fast01
TARGETS = hw2seq hw2a hw2b

.PHONY: all
all: $(TARGETS)

run_seq:
	srun -n1 -c1 ./hw2seq out/out.png 2602 -3 0.2 -3 0.2 979 2355
run:
	srun -n1 -c5 ./hw2a out/$(outname).png 2602 -3 0.2 -3 0.2 979 2355

.PHONY: clean
clean:
	rm -f $(TARGETS) $(TARGETS:=.o)
