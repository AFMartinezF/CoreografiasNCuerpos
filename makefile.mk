#make -f Makefile_MartinezAngel_S7C2_EDO2.mk
all :: Resultados.png

Resultados.png :: Archivo.py Resultados.csv
					python Archivo.py

Resultados.csv :: SimulacionMetodoEuler.exe
					./SimulacionMetodoEuler.exe

SimulacionMetodoEuler.exe :: SimulacionMetodoEuler.cpp
					g++ SimulacionMetodoEuler.cpp -o SimulacionMetodoEuler.exe
