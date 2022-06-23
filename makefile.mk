#make -f makefile.mk

all :: Grafica3D.gif

Grafica3D.gif :: Grafica3D.py SolGeneral.csv
										python Grafica3D.py

SolGeneral.csv :: SolGeneral.py
										python SolGeneral.py
