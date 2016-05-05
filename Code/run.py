import os

os.system("python rand.py")
os.system("mpif90 -I ~/Documents/program/Fortran/Library/ ~/Documents/program/Fortran/Library/*.o /usr/local/lib/libfftw3.a  main.f90")
os.system("./a.out")
