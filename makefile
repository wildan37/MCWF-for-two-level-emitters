CC = gcc
CFLAGS = -O2
DEPS = header_noise.h

OBJ = main.o noise_process.o init_noisedyn.o mt19937ar.o histogram.o 

all: $(OBJ)
	$(CC) $(CFLAGS) -lm -lgsl -lgslcblas $(OBJ) 

debug: 
	$(MAKE) CFLAGS="-g"

init_noisedyn.o: init_noisedyn.c $(DEPS)
	$(CC) -c $(CFLAGS) init_noisedyn.c

noise_process.o: noise_process.c $(DEPS)
	$(CC) -c $(CFLAGS) noise_process.c

mt19937ar.o: mt19937ar.c mt19937ar.h
	$(CC) -c $(CFLAGS) mt19937ar.c 

histogram.o: histogram.c histogram.h
	$(CC) -c $(CFLAGS) histogram.c 

main.o: main.c $(DEPS) 
	$(CC) -c $(CFLAGS) main.c

clean:
	rm -rf $(OBJ)




