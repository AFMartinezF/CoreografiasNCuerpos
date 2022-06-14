#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;

double distancia (double r_cuerpo1[][2], double r_cuerpo2[][2],int dim);
double dv_euler(double  G, double cuerpos, double m[], double r[]);
double dy_euler(double v);

int main() {
  double G = 1.2;
  int dim = 2;
  int n = 100;
  double t_0 = 0;
  double t_final = 10;
  double dt = t_final/n;
  double r1 [n][2] = {{-1,0}};
  double r2 [n][2] = {{1,0}};
  double r3 [n][2] = {{0,0}};
  double v1 [n][2] = {{-0.3,-0.5}};
  double v2 [n][2] = {{-0.3,-0.5}};
  double v3 [n][2] = {{0.7,1.1}};
  double m [3] = {1,1,1};

  double d12 = distancia (r1, r2,dim);
  double d13 = distancia (r1, r3,dim);

  v1[1][0] = v1[0][0] + dt*dv_euler(G,3,m,r1d);

  //cout<<v1[1][0]<<endl;
  //cout<<dv_euler(G,3,m,r23)<<endl;

  //for(int i = 1; i<n; i++){
    //for(int j = 1; j<n; j++){
    //Función usando metodo de Euler

    //v1[i,j] = v1[i-1,j-1] + dt*v_euler(G,3,m,r23);

    //}
  //}
  return 0;
}

double distancia (double r_cuerpo1[][2], double r_cuerpo2[][2],int dim){
  double dif [dim];
  double d;
  for (int i = 0; i< dim; i++){
    dif [i] = r_cuerpo1[i][i]-r_cuerpo2[i][i];
    //cout<<r_cuerpo1[i][i]<<endl;
    cout<<dif[i]<<endl;
  }
  d = sqrt(pow(dif[0],2)+pow(dif[1],2));
  return d;
}

double dv_euler(double  G, double cuerpos, double m[], double r[]){
  double v = 0;
  for (int i = 0; i < cuerpos; i++){
    v = v - G*m[i+1]/r[i];
  }
  return v;
}

double dy_euler(double v){
  return  v;
}
  //r[0] = -G*(m[0]/r[0,0]+m[1]/r[1,0])
  //r[1] = -G*(m[0]/r[0,1]+m[1]/r[1,1])
