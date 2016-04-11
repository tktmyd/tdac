# makefile for tdac

FC=gfortran
FFLAGS=
TARGET:../tdac.x

.SUFFIXES:
.SUFFIXES: .f90 .o

m_params.o: m_params.f90
m_llw2d.o: m_llw2d.f90 m_params.o
m_matrix.o: m_matrix.f90

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@

../tdac.x: m_params.o m_llw2d.o m_matrix.o tdac.f90
	$(FC) $(FFLAGS) -o $@ $^

clean:
	/bin/rm -f *.o *.mod *~ ../tdac.x
