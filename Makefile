TARGET = boruvka

all: $(TARGET)

boruvka: boruvka.h boruvka.cpp graph_tools.cpp seq_generation_utils.h defs.h
	mpic++ -Wall -O3 boruvka.cpp -lm -o boruvka

clean:
	rm -rf *.o $(TARGET)
	rm boruvka
