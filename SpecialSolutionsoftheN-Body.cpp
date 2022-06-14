//Este código fue extraido del artículo Special Solutions of the N-Body Problem: Central Configurations and Choreographies
//http://diposit.ub.edu/dspace/bitstream/2445/181716/2/tfg_victor_sanchez_linan.pdf
//Por tanto pertenece a los autores del mismo
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

double min ( double var1 , double var2 );
int integ_taylor ( double ***z, double *m, double *t, double tmax , double *h,
  double hmin , double hmax , int ordre , int n, int dim , double tol);
/* defineixo :
 - z vector del camp ( z = (q,p); ( ordre +1) *2n* dim )
 - q vector de posicions (( ordre +1) *n* dim ) [ der ][ num ][ coord ]
 - p vector de moments (( ordre +1) *n* dim ) [ der ][ num ][ coord ]
 - m vector de masses (n)
 - t instant de temps
 - tmax temps maxim de l'interval temporal de definicio
 - h pas de Taylor
 - hmin cota inferior pel pas de Taylor
 - hmax cota superior pel pas de Taylor
 - ordre ordre del polinomi de Taylor
 - n nombre de cossos
 - dim dimensio espacial del problema
 - tol tolerancia*/
int integ_taylor ( double ***z, double *m, double *t, double tmax , double *h,
  double hmin , double hmax , int ordre , int n, int dim , double tol ) {
 /* el vector z =(q,p) ve donat per q [0][][] = z[0][0 - >n -1][] i p [0][][] =
z [0][n ->2n -1] i m ve donat per m[0->n -1] */

 /* declaracio de variables */
  int i,j,k,l,s;
  int flag ;
  double alfa = -1.5 , suma , max1 , max2 , h1 , h2;
  double **** difq , *difq2 , *** sumdifq2 , *** norma , ** prod ;

  difq = ( double ****) malloc (( ordre +1)* sizeof ( double ***) );
  difq2 = ( double *) malloc ( dim * sizeof ( double ));
  sumdifq2 = ( double ***) malloc (( ordre +1) * sizeof ( double **) );
  norma = ( double ***) malloc (( ordre +1)* sizeof ( double **) );
  prod = ( double **) malloc (n* sizeof ( double *));
  for (k = 0; k <= ordre ; k ++) {
      difq [k] = ( double ***) malloc (n* sizeof ( double *));
      sumdifq2 [k] = ( double **) malloc (n* sizeof ( double ));
      norma [k] = ( double **) malloc (n* sizeof ( double ));
      for (j = 0; j < n; j ++) {
        difq [k][j] = ( double **) malloc (n* sizeof ( double *));
        sumdifq2 [k][j] = ( double *) malloc (n* sizeof ( double ));
        norma [k][j] = ( double *) malloc (n* sizeof ( double ));
        for (i = 0; i < n; i ++)
        difq [k][j][i] = ( double *) malloc ( dim * sizeof ( double ));
      }
    }
    for (i = 0; i < n; i ++)
      prod [i] = ( double *) malloc ( dim * sizeof ( double ));


 /* Derivada ordre 1 */
 /* Calculem la derivada q [1] */
  for (j = 0; j < n; j ++) {
    for (l = 0; l < dim ; l++)
    z [1][ j][l] = z [0][ j+n][l]/m[j];
  }
 /* Calculem la derivada p [1] */
  for (j = 0; j < n; j ++) {
    for (i = 0; i < j; i ++) {
 /* calculem diferencia entre les q */
      for (l = 0; l < dim ; l++)
        difq [0][ j][i][l] = z [0][ i][l] - z [0][j][l];
 /* calculem el quadrat de les diferencies */
      for (l = 0; l < dim ; l++)
        difq2 [l] = difq [0][ j][i][l]* difq [0][j][i][l];
 /* calculem la suma de les diferencies al quadrat */
      sumdifq2 [0][ j][i] = 0;
      for (l = 0; l < dim ; l++)
        sumdifq2 [0][ j][i] += difq2 [l];
      if ( fabs ( sumdifq2 [0][ j][i]) <= tol )
        return -1;
 /* calculem la norma del denominador */
      norma [0][ j][i] = exp ( alfa *log( sumdifq2 [0][j][i]));
 /* calculem el quocient */
      for (l = 0; l < dim ; l ++)
        prod [i][l] = difq [0][ j][i][l]* norma [0][j][i];
      }
      for (i = j+1; i < n; i ++) {
          for (l = 0; l < dim ; l ++)
            difq [0][ j][i][l] = z [0][ i][l] - z [0][j][l];
          for (l = 0; l < dim ; l ++)
            difq2 [l] = difq [0][ j][i][l]* difq [0][j][i][l];
          sumdifq2 [0][ j][i] = 0;
          for (l = 0; l < dim ; l ++)
            sumdifq2 [0][ j][i] += difq2 [l];
          if ( fabs ( sumdifq2 [0][ j][i]) <= tol )
            return -1;
          norma [0][ j][i] = exp ( alfa *log( sumdifq2 [0][ j][i]));
          for (l = 0; l < dim ; l ++)
            prod [i][l] = difq [0][ j][i][l]* norma [0][ j][i];
          }
 /* per a cada cos j, calculem el seu moment component a component ,
sumant sobre les i */
          for (l = 0; l < dim ; l ++) {
              z [1][ j+n][l] = 0;
              for (i = 0; i < j; i ++)
                z [1][ j+n][l] += m[i]* prod [i][l];
              for (i = j+1; i < n; i ++)
                z [1][ j+n][l] += m[i]* prod [i][l];
              z [1][ j+n][l ]*= m[j];
            }
          }

 /* Derivada ordre k > 1 */
  for (k = 2; k <= ordre ; k ++) {
    for (j = 0; j < n; j ++) {
 /* Calculem les derivades q[k] */
      for (l = 0; l < dim ; l ++)
        z[k][j][l] = z[k -1][ j+n][l ]/( k*m[j]);
 /* Calculem les derivades p[k] */
      for (i = 0; i < j; i ++) {
 /* calculem diferencia entre les q */
      for (l = 0; l < dim ; l ++)
        difq [k -1][ j][i][l] = z[k -1][ i][l] - z[k -1][ j][l];
 /* calculem el quadrat de les diferencies */
      for (l = 0; l < dim ; l ++) {
        difq2 [l] = 0;
        for (s = 0; s <= k -1; s++)
          difq2 [l] += difq [k -1-s][j][i][l]* difq [s][j][i][l];
      }
 /* calculem la suma de les diferencies al quadrat */
        sumdifq2 [k -1][ j][i] = 0;
        for (l = 0; l < dim ; l ++)
          sumdifq2 [k -1][ j][i] += difq2 [l];
 /* calculem la norma , considerem alfa = -1.5 i multipliquem */
        suma = 0;
        for (s = 0; s <= k -2; s++)
          suma += ((k -1) * alfa - s*( alfa +1) )* sumdifq2 [k -1-s][j][i]*
            norma [s][j][i];
          norma [k -1][ j][i] = suma /((k -1)* sumdifq2 [0][j][i]);
/* calculem el quocient */
        for (l = 0; l < dim ; l++) {
          prod [i][l] = 0;
          for (s = 0; s <= k -1; s++)
            prod [i][l] += difq [s][j][i][l]* norma [k -1-s][j][i];
        }
      }
      for (i = j+1; i < n; i++) {
        for (l = 0; l < dim ; l++)
          difq [k -1][ j][i][l] = z[k -1][ i][l] - z[k -1][ j][l];
        for (l = 0; l < dim ; l++) {
          difq2 [l] = 0;
        for (s = 0; s <= k -1; s++)
          difq2 [l] += difq [k -1-s][j][i][l]* difq [s][j][i][l];
        }
        sumdifq2 [k -1][ j][i] = 0;
        for (l = 0; l < dim ; l++)
          sumdifq2 [k -1][ j][i] += difq2 [l];
        suma = 0;
        for (s = 0; s <= k -2; s++)
          suma += ((k -1) * alfa - s*( alfa +1) )* sumdifq2 [k -1-s][j][i]*
            norma [s][j][i];
          norma [k -1][ j][i] = suma /((k -1)* sumdifq2 [0][ j][i]);
        for (l = 0; l < dim ; l++) {
          prod [i][l] = 0;
          for (s = 0; s <= k -1; s++)
            prod [i][l] += difq [s][j][i][l]* norma [k -1-s][j][i];
          }
    }
 /* per a cada cos j, calculem el seu moment component a component ,
sumant sobre les i */
    for (l = 0; l < dim ; l ++) {
      suma = 0;
      for (i = 0; i < j; i ++)
        suma += m[i]* prod [i][l];
        for (i = j+1; i < n; i ++)
        suma += m[i]* prod [i][l];
        z[k][j+n][l] = m[j]* suma /k;
      }
    }
 }

 /* triem pas h */
  max1 = tol; max2 = tol ;
  for (i = 0; i < 2*n; i++) {
    for (l = 0; l < dim ; l++) {
 /* norma subinfinit de z per la derivada d'ordre k=ordre -1 */
      if ( fabs (z[ordre -1][ i][l]) > max1 )
        max1 = fabs (z[ordre -1][ i][l]);
 /* norma subinfinit de z per la derivada d'ordre k= ordre */
      if ( fabs (z[ ordre ][i][l]) > max2 )
        max2 = fabs (z[ ordre ][i][l]);
      }
    }
    h1 = pow ( tol /max1 ,1./( ordre -1) );
    h2 = pow ( tol /max2 ,1./ ordre );
    *h = min (h1 ,h2);
  /* Miro que h estigui en l'interval correcte */
    if (*h > hmax ) *h = hmax ;
    else if (*h < hmin ) {
      *h = hmin ;
      flag = 2;
    }
  /* mirem no passar - nos de temps */
    if (*t + *h > tmax ) {
      *h = tmax - *t;
      flag = 1;
    }
    *t = *t + *h; /* nou temps */

  /* sumem fent Horner pel polinomi de Taylor ( component a component ) */
    for (j = 0; j < 2*n; j ++) {
      for (l = 0; l < dim ; l ++) {
        suma = z[ ordre ][j][l];
        for (k = ordre -1; k >= 0; k --)
          suma = suma *(* h) + z[k][j][l];
        z [0][ j][l] = suma ;
      }
    }

  /* alliberem memoria */
    for (k = 0; k <= ordre ; k ++) {
      for (j = 0; j < n; j ++) {
        for (i = 0; i < n; i ++)
          free ( difq [k][j][i]);
        free ( difq [k][j]); free ( sumdifq2 [k][j]); free ( norma [k][j]);
      }
      free ( difq [k]); free ( sumdifq2 [k]); free ( norma [k]);
    }
    for (i = 0; i < n; i ++)
      free ( prod [i]);
    free ( difq ); free ( difq2 ); free ( sumdifq2 ); free ( norma ); free ( prod );

    return flag ;
  }
