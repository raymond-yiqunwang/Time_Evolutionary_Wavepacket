VPATH = ~/Document/program/Fortran/Library
objects = main.o
a.out :  $(objects)
	mpif90 -I ~/Documents/program/Fortran/Library/ ~/Documents/program/Fortran/Library/*.o /usr/local/lib/libfftw3.a -o a.out $(objects)

main.o	:	main.f90
	mpif90 -I ~/Documents/program/Fortran/Library -c main.f90

.PHONY	:	clean
clean	:
	-rm edit $(objects)
