#INC_DIR=weilei_lib
INC_DIR=.
#INC_DIR=~/working/weilei_lib
CXX=g++ -O3 -Wall -std=c++11
### -O2 -O5 -Os
#g++ `pkg-config --cflags itpp` -o hello.out hello.cpp `pkg-config --libs itpp`

START=`pkg-config --cflags itpp`
END=`pkg-config --libs itpp`

ITPP=`pkg-config --cflags itpp` `pkg-config --libs itpp`

#files=$(INC_DIR)/mm_read.cpp $(INC_DIR)/mm_read.h $(INC_DIR)/mmio.c $(INC_DIR)/mmio.h $(INC_DIR)/mm_write.cpp $(INC_DIR)/mm_write.h $(INC_DIR)/lib.cpp $(INC_DIR)/lib.h $(INC_DIR)/dist.cpp $(INC_DIR)/dist.h $(INC_DIR)/concatenation_lib.cpp $(INC_DIR)/concatenation_lib.h  $(INC_DIR)/bp.cpp $(INC_DIR)/bp.h $(INC_DIR)/my_lib.h Makefile
files=$(INC_DIR)/mm_read.cpp $(INC_DIR)/mm_read.h $(INC_DIR)/mmio.c $(INC_DIR)/mmio.h $(INC_DIR)/mm_write.cpp $(INC_DIR)/mm_write.h $(INC_DIR)/lib.cpp $(INC_DIR)/lib.h $(INC_DIR)/dist.cpp $(INC_DIR)/dist.h $(INC_DIR)/product.cpp $(INC_DIR)/product.h  $(INC_DIR)/bp.cpp $(INC_DIR)/bp.h $(INC_DIR)/weilei_lib.h Makefile
command=$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(word 8, $^) $(word 10, $^) $(word 12, $^) $(word 14, $^) $(END)

###include all weilei-written headfiles into weilei_lib.h



all: cpp head
cpp: mmio.o mm_read.o mm_write.o dist.o bp.o lib.o product.o
head: bp_decoder.h.gch weilei_lib.h.gch

#compile object file for cpp 
%.o:%.cpp %.h
	$(CXX) $(ITPP) -c $<
#compile object file for headfile
%.h.gch:%.h
	$(CXX) $(ITPP) -c $<



#mmio:mmio.c mmio.h
#	$(CXX) $(ITPP) -c $<
#mm_read:mm_read.c mm_read.h mmio.c mmio.h
#	$(CXX) $(ITPP) -c $<
#mm_write:mm_write.c mm_write.h mmio.c mmio.h
#	$(CXX) $(ITPP) -c $<
#lib:lib.cpp lib.h
#	$(CXX) $(ITPP) -c $<
#concatenation:concatenation.cpp concatenation.h
#	$(CXX) $(ITPP) -c $<


gnuplot_dist.out:gnuplot_dist.c $(files)
	$(command)
random_concatenation.out:random_concatenation.c $(files)
	$(command)
counter_concatenation.out:counter_concatenation.c $(files)
	$(command)
product.out:product.c $(files)
	$(command) -fopenmp
concatenation.out:concatenation.c $(files)
	$(command)
hypergraph.out:hypergraph.c $(files)
	$(command)

test.out:test.c $(files)
	$(command)

clean:
	rm  *~
	rm \#*

sbatch-dry-run:
	sbatch --test run_prod.sh
sbatch:
	sbatch run_prod.sh
pkill-product:
	pkill .product
