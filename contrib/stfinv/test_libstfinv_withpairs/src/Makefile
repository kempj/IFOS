# Makefile for Test_libstfinv

#--------------------------------------------------------
# edit here:


CC=gcc

#LFLAGS=-lstfinv -lfourierxx -lfftw3 -laff -lstdc++ 
LFLAGS=-lstfinv -lfourierxx -lfftw3 -laff 


CFLAGS= -Wall

all:	Test_libstfinv
	mv Test_libstfinv ..

.c.o:
	$(CC) -c -I$(LOCINCLUDEDIR) $(CFLAGS) $<


FFT_Test= \
	read_parameters.c \
	Test_libstfinv.c  
	
FFT_Test_OBJ = $(FFT_Test:%.c=%.o)

Test_libstfinv:		$(FFT_Test_OBJ) 
	$(CC) $(CFLAGS) -I$(LOCINCLUDEDIR) -L$(LOCLIBDIR) $(FFT_Test_OBJ) -o Test_libstfinv $(LFLAGS)


clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	rm -f ../Test_libstfinv
