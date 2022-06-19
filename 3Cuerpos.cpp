#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;

double distancia (double r_cuerpo1[][2], double r_cuerpo2[][2],int dim,int iteracion);
//double diferencia(double r_cuerpo1[][2], double r_cuerpo2[][2]);
double dv_euler(double  G, double m, double r_cuerpo1[][2], double r_cuerpo2[][2],int dim,int iteracion);
double dy_euler(double v);

int main() {
  double G = 1.2; //Constante de gravitación universal
  int dim = 2;    //Grados de libertad
  int n = 100;
  double t_0 = 0;
  double t_final = 1;
  double dt = t_final/n;
  int i = 0;
  //Condiciones inciales del problema (posición, velocidad)
  double r1 [n][2] = {{-1,0}};
  double r2 [n][2] = {{1,0}};
  double r3 [n][2] = {{0,0}};
  double v1 [n][2] = {{-0.3,-0.5}};
  double v2 [n][2] = {{-0.3,-0.5}};
  double v3 [n][2] = {{0.7,1.1}};
  double d12;
  double d13;
  double d23;
  double m [3] = {1,1,1}; //Masa de los cuerpos
  double eul12; double eul13; double eul21; double eul23; double eul31; double eul32;
  cout<< v1[0][0] << endl;
  cout<< v1[0][1] << endl;

  for(int t = t_0; t < t_final; t = t + dt){
    //Distancias relativas entre cuerpos
    d12 = distancia (r1, r2,dim,i);
    d13 = distancia (r1, r3,dim,i);
    d23 = distancia (r2, r3,dim,i);
    cout << "Distacia d12: "<< d12 << endl;
    cout << "Distacia d13: "<< d13<< endl;
    cout << "Distacia d23: "<< d23 << endl;
    eul12 = dv_euler(G, m[1],r1, r2, dim, i);
    eul13 = dv_euler(G, m[2],r1, r3, dim, i);
    eul21 = dv_euler(G, m[0],r2, r1, dim, i);
    eul23 = dv_euler(G, m[2],r2, r3, dim, i);
    eul31 = dv_euler(G, m[0],r3, r1, dim, i);
    eul32 = dv_euler(G, m[1],r3, r2, dim, i);

    v1[i+1][0] = v1[i][0] + dt*(eul12*(r1[i][1]-r2[i][1])/(d12)+eul13*(r1[i][1]-r3[i][1])/(d13));
    v1[i+1][1] = v1[i][1] + dt*(eul12*(r1[i][0]-r2[i][0])/(d12)+eul13*(r1[i][0]-r3[i][0])/(d13));
    r1[i+1][0] = r1[i][0] + dt*(dy_euler(v1[i][0])*(r1[i][1]-r2[i][1])/(d12)+dy_euler(v1[i][0])*(r1[i][1]-r3[i][1])/(d13));
    r1[i+1][1] = r1[i][1] + dt*(dy_euler(v1[i][1])*(r1[i][0]-r2[i][0])/(d12)+dy_euler(v1[i][1])*(r1[i][0]-r3[i][0])/(d13));

    v2[i+1][0] = v2[i][0] + dt*(eul21*(r1[i][1]-r2[i][1])/(d12)+eul23*(r2[i][1]-r3[i][1])/(d23));
    v2[i+1][1] = v2[i][1] + dt*(eul21*(r1[i][0]-r2[i][0])/(d12)+eul23*(r2[i][0]-r3[i][0])/(d23));
    r2[i+1][0] = r2[i][0] + dt*(dy_euler(v2[i][0])*(r1[i][1]-r2[i][1])/(d12)+dy_euler(v2[i][0])*(r2[i][1]-r3[i][1])/(d23));
    r2[i+1][1] = r2[i][1] + dt*(dy_euler(v2[i][1])*(r1[i][0]-r2[i][0])/(d12)+dy_euler(v2[i][1])*(r2[i][0]-r3[i][0])/(d23));

    v3[i+1][0] = v3[i][0] + dt*(eul31*(r1[i][1]-r3[i][1])/(d13)+eul32*(r2[i][1]-r3[i][1])/(d23));
    v3[i+1][1] = v3[i][1] + dt*(eul31*(r1[i][0]-r3[i][0])/(d13)+eul32*(r2[i][0]-r3[i][0])/(d23));
    r3[i+1][0] = r3[i][0] + dt*(dy_euler(v3[i][0])*(r1[i][1]-r3[i][1])/(d13)+dy_euler(v2[i][0])*(r2[i][1]-r3[i][1])/(d23));
    r3[i+1][1] = r3[i][1] + dt*(dy_euler(v3[i][1])*(r1[i][0]-r3[i][0])/(d13)+dy_euler(v2[i][1])*(r2[i][0]-r3[i][0])/(d23));

    cout<< i<<endl;
    cout<< v1[i][0] << endl;
    cout<< v1[i][1] << endl;

    //cout << d12 << endl << d13 << endl << d23 << endl;

    i ++;
  }



  //Ejecución del metodo de Euler (revisar)
  //v1[1][0] = v1[0][0] + dt*(dv_euler(G, m[1],r1, r2, dim)*(r1[0]-r2[0])/(d12)+dv_euler(G, m[2],r1, r3, dim)*(r1[0]-r3[0])/(d13));

  //cout<< v1[1][0] << endl;
  //cout<< v1[0][1] << endl;

  return 0;
}

//Función que calcula la distancia entre dos puntos
double distancia (double r_cuerpo1[][2], double r_cuerpo2[][2],int dim, int iteracion){
  double diferencias [dim];
  double distancia;
  for (int i = iteracion; i< iteracion+dim; i++){
    diferencias [i] = r_cuerpo1[i][i]-r_cuerpo2[i][i];
    //cout<<dif[i]<<endl;
  }
  distancia = sqrt(pow(diferencias[0],2)+pow(diferencias[1],2));
  return distancia;
}

//double diferencia(double r_cuerpo1[][2], double r_cuerpo2[][2]){
  //double diferencia;
  //int i = 0;
  //diferencia = r_cuerpo1[i][i]-r_cuerpo2[i][i];
  //return diferencia;
//}

//Funciones del metodo de euler para una ecuación de segundo orden
double dv_euler(double  G, double m, double r_cuerpo1[][2], double r_cuerpo2[][2],int dim,int iteracion){
  double diferencias [dim];
  for (int i = iteracion; i< iteracion+dim; i++){
    diferencias [i] = r_cuerpo1[i][i]-r_cuerpo2[i][i];
    //cout<<dif[i]<<endl;
  }
  double distancia = sqrt(pow(diferencias[0],2)+pow(diferencias[1],2));
  double v = (-G*m/pow(distancia,2)); //Unir con las otras funciones? Llamar a las otras funciones?
  return v;//generalizar
}

double dy_euler(double v){
  return  v;
}
