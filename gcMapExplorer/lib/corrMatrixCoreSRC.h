#ifndef CORRMATRIXCORE_H
#define CORRMATRIXCORE_H

int _correlationMatrix(double *in_array, double *out_array, int size, double *maskvalue);
int _covarianceMatrix(double *in_array, double *out_array, int size, double *maskvalue);
double _calculateCorrelation(double *x, double*y, int n, double *maskvalue);
double _calculateCovariance(double *x, double*y, int n, double *maskvalue);

#endif
