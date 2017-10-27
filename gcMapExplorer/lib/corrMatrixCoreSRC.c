#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>

#include "corrMatrixCoreSRC.h"

int _correlationMatrix(double *in_array, double *out_array, int size, double *maskvalue)
{
  double **in, **out, r = 0.0, random =0;
  int i=0, j=0,  counter =0, nthreads, chunk = 50;

  // Convert 1D array to 2D array by pointers
  in = (double**) malloc(size * sizeof(double*));
  out = (double**) malloc(size * sizeof(double*));
  for(i=0; i<size; i++)
  {
    in[i] = &in_array[i*size];
    out[i] = &out_array[i*size];
  }

  // Determine chunk for multithreading
  nthreads = omp_get_max_threads();
  if(size < 50)
  {
    chunk = size/nthreads;
  }

  // Main loop for calculation
  #pragma omp parallel for schedule(dynamic, chunk) private(i, j, r)
  for(i=0; i<size; i++)
  {
    for(j=0; j<i; j++)
    {
      r = _calculateCorrelation(in[i], in[j], size, maskvalue);
      out[i][j] = r;
      out[j][i] = r;

      if((counter%100)==0)
      {
        printf("\r               ... %d/%d done ...", counter, size);
        fflush(stdout);
      }

    }
    counter += 1;
  }
  printf("\n");


  // Free 2D pointer array, be careful, assign a random address and the free it.
  for(i=0; i<size; i++)
  {
    in[i] = &random;
    out[i] = &random;
  }
  free(in);
  free(out);

  return 0;
}

int _covarianceMatrix(double *in_array, double *out_array, int size, double *maskvalue)
{
  double **in, **out, random=0, cov = 0.0;
  int i=0, j=0, nthreads, chunk = 50;

  // Convert 1D array to 2D array by pointers
  in = (double**) malloc(size * sizeof(double*));
  out = (double**) malloc(size * sizeof(double*));
  for(i=0; i<size; i++)
  {
    in[i] = &in_array[i*size];
    out[i] = &out_array[i*size];
  }

  nthreads = omp_get_max_threads();
  if(size < 50)
  {
    chunk = size/nthreads;
  }

  // Main loop for calculation
  #pragma omp parallel for schedule(dynamic, chunk) private(i, j, cov)
  for(i=0; i<size; i++)
  {
    for(j=0; j<=i; j++)
    {
      cov = _calculateCovariance(in[i], in[j], size, maskvalue);
      out[i][j] = cov;
      out[j][i] = cov;
    }
  }


  // Free 2D pointer array, be careful, assign a random address and the free it.
  for(i=0; i<size; i++)
  {
    in[i] = &random;
    out[i] = &random;
  }
  free(in);
  free(out);

  return 0;
}


double _calculateCorrelation(double *x, double*y, int n, double *maskvalue)
{
  double corr = 0, cov_xy=0, var_x=0, var_y=0;

   /*
  int i=0
  for(i=0; i<n; i++)
  {
    printf(" %f %f \n", x[i], y[i]);
  }
   */

  cov_xy = _calculateCovariance(x, y, n, maskvalue);
  var_x = _calculateCovariance(x, x, n, maskvalue);
  var_y = _calculateCovariance(y, y, n, maskvalue);


  if ( (var_x == 0) || (var_y == 0) )
  {
    corr = 0;
  }
  else
  {
    corr = cov_xy / (sqrt(var_x * var_y));
  }

  //printf("%f %f %f %f \n",cov_xy, var_x, var_y, corr);
  return corr;
}

double _calculateCovariance(double *x, double*y, int n, double *maskvalue)
{
  double xmean=0, ymean=0, counter = 0;
  double mask=0.0, cov=0;
  double xsum=0, ysum=0, dev=0;
  int i = 0;

  if (!(maskvalue == NULL))
  {
    mask = *maskvalue;
  }


  for(i=0; i<n; i++)
  {

    if (maskvalue == NULL)
    {
      xsum +=  x[i];
      ysum +=  y[i];
      counter += 1;
    }
    else
    {
      if ( (x[i] != mask ) && (y[i] != mask) )
      {
        xsum +=  x[i];
        ysum +=  y[i];
        counter += 1;
      }
    }
  }

  if (counter > 3)
  {
    xmean = xsum/counter;
    ymean = ysum/counter;

    for(i=0;i<n;i++)
    {
      if (maskvalue == NULL)
      {
        dev += ( (  x[i]-xmean) * (  y[i]-ymean) );
      }
      else
      {
        if ( (x[i] != mask ) && (y[i] != mask) )
        {
          dev += ( (  x[i]-xmean) * (  y[i]-ymean) );
        }
      }
    }
    cov =  dev / (counter);
  }
  return cov;
}
