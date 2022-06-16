#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;

double distancia (double r_cuerpo1[][2], double r_cuerpo2[][2],int dim);
double diferencia(double r_cuerpo1[][2], double r_cuerpo2[][2]);
double dv_euler_x(double  G, double m, double diferencia, double distancia);
double dy_euler(double v);

int main() {
  double G = 1.2; //Constante de gravitación universal
  int dim = 2;    //Grados de libertad
  int n = 100;
  double t_0 = 0;
  double t_final = 10;
  double dt = t_final/n;
  //Condiciones inciales del problema (posición, velocidad)
  double r1 [n][2] = {{-1,0}};
  double r2 [n][2] = {{1,0}};
  double r3 [n][2] = {{0,0}};
  double v1 [n][2] = {{-0.3,-0.5}};
  double v2 [n][2] = {{-0.3,-0.5}};
  double v3 [n][2] = {{0.7,1.1}};
  double m [3] = {1,1,1}; //Masa de los cuerpos

  //Distancias relativas entre cuerpos
  double d12 = distancia (r1, r2,dim);
  double d13 = distancia (r1, r3,dim);
  double d23 = distancia (r1, r3,dim);
  cout << d12 << endl << d13 << endl << d23 << endl;

  //Ejecución del metodo de Euler (revisar)
  cout<< v1[0][0] << endl;
  v1[1][0] = v1[0][0] + dt*(dv_euler_x(G,m[1],diferencia(r1,r2),distancia(r1, r2,dim))+dv_euler_x(G,m[2],diferencia(r1,r3),distancia(r1, r3,dim)));

  cout<< v1[1][0] << endl;

  return 0;
}

//Función que calcula la distancia entre dos puntos
double distancia (double r_cuerpo1[][2], double r_cuerpo2[][2],int dim){
  double diferencias [dim];
  double distancia;
  for (int i = 0; i< dim; i++){
    diferencias [i] = r_cuerpo1[i][i]-r_cuerpo2[i][i];
    //cout<<dif[i]<<endl;
  }
  distancia = sqrt(pow(diferencias[0],2)+pow(diferencias[1],2));
  return distancia;
}

double diferencia(double r_cuerpo1[][2], double r_cuerpo2[][2]){
  double diferencia;
  int i = 0;
  diferencia = r_cuerpo1[i][i]-r_cuerpo2[i][i];
  return diferencia;
}


//Funciones del metodo de euler para una ecuación de segundo orden
double dv_euler_x(double  G, double m, double diferencia, double distancia){
  double v = (-G*m/pow(distancia,2))*cos(diferencia/distancia); //Unir con las otras funciones? Llamar a las otras funciones?
  return v;//generalizar
}

double dy_euler(double v){
  return  v;
}
