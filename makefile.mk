#make -f makefile.mk
all :: df.gif

df.gif :: Grafica.py Resultados.csv
					python Grafica.py

Resultados.csv :: SimulacionMetodoEuler.exe
					./SimulacionMetodoEuler.exe

SimulacionMetodoEuler.exe :: SimulacionMetodoEuler.cpp
					g++ SimulacionMetodoEuler.cpp -o SimulacionMetodoEuler.exe
