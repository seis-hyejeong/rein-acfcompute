gsldir= your_gsl_directory
CC = gcc
CFLAGS = -O2 -Wall -lm
gsl_lib = -L$(gsldir)/lib -lgsl -lgslcblas
gsl_include = -I$(gsldir)/include

default: all

all: 
	$(CC) -c $(CFLAGS) sacioH.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) readwrite_w.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) complex_pack.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) matrix_pack.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) fft_pack.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) calculation_w.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) chris_only.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) acf_functions.c
	$(CC) -c $(CFLAGS) $(gsl_lib) $(gsl_include) acf_functions_compprint.c
	$(CC) -o acf_water_comp_tauzin $(CFLAGS) $(gsl_lib) $(gsl_include) acf_tele_smoothing_tauzin_comp_main.c readwrite_w.o complex_pack.o matrix_pack.o fft_pack.o calculation_w.o chris_only.o sacioH.o acf_functions_compprint.o
	cp acf_water_comp_tauzin ../bin/
clean:
	rm -rf *.o
