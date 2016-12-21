#ifndef _RANDOM_H
#define _RANDOM_H

#define Ranuni RandomUniform()
#define Rannor StdNormal()
#include <math.h>

double RandomUniform();
/*inline*/ void Seeds(int ix, int iy, int iz)
{ RandomUniform(); } /* j'ai enleve ix, iy, iz en argument de RandomUniform*/
double StdNormal();
double ChiSquare(int ndf);
double Beta(int v1, int v2);
double Gamma(int num, double rate);
double Weibull(double rate, double shape);

double RandomUniform()
{ 
  static int ix = 10437;
  static int iy = 13568;
  static int iz = 30524;
  double r;
  int rint;
  int newx=0; int newy=0; int newz=0;                    
  if (newx == 0)
  {
    ix = 171*(ix-177*((int) floor(ix/177.0)))
    -  2*((int) floor(ix/177.0));
    iy = 172*(iy-176*((int) floor(iy/176.0)))
      - 35*((int) floor(iy/176.0));
    iz = 170*(iz-178*((int) floor(iz/178.0)))
      - 63*((int) floor(iz/178.0));
    
    if (ix < 0) ix += 30269;
    if (iy < 0) iy += 30307;
    if (iz < 0) iz += 30323;
    
    r = ix/30269.0 + iy/30307.0 + iz/30323.0;
    rint = (int) floor(r);
    return (r - rint);
  }
  else
  {
    ix = newx;
    iy = newy;
    iz = newz;
    return 0;
  }
}

double StdNormal()
{
  static int ir = 0;
  static double an = 0;
  double e, w, v1, v2;
  
  int newr=-1; double newn=0;
  if (newr == -1)
  {
    if (ir == 0) {
      do 
      {
        v1 = 2.0*RandomUniform() - 1.0;
        v2 = 2.0*RandomUniform() - 1.0;
        w  = v1*v1 + v2*v2;
      } while (w > 1);
      e  = sqrt((-2.0*log(w))/w);
      an = v1 * e;
      ir = 1;
      return v2 * e;
    } 
    else 
    {
      ir = 0;
      return an;
    }
  }
  else
  {
    ir = newr;
    an = newn;
    return 0;
  }
}

double ChiSquare(int ndf)
{
  int n;
  double srsq, stdNormal;
  ndf=1;
  srsq = 0.0;
  for (n = 1; n <= ndf; n++)
  {
    stdNormal = StdNormal();
    srsq += (stdNormal * stdNormal);
  }
  return srsq;
}

double Beta(int v1, int v2)
{
  double x, y, totaldf, e, v, logite, randlogite, erle;
  int i, j;
  
  totaldf = v1 + v2;
  if (totaldf < 100) {
    /* Ratio of sigma's method */
    x = 0;
    y = 0;
    for (i = 1; i <= v1; i++)
      x -= log(RandomUniform());
    for (j = 1; j <= v2; j++)
      y -= log(RandomUniform());
    return x / (x + y);
  } else {
    /* Normal Approximation */
    e = v1 / totaldf;
    v = 1.0 / (totaldf * e * (1 - e));
    logite = log(e / (1 - e));
    randlogite = logite + StdNormal() * sqrt(v);
    erle = exp(randlogite);
    return erle / (1 + erle);
  }
}

double Gamma(int num , double rate )
{
  int i;
  double rand;
  
  rand = 0;
  for (i = 1; i <= num; i++)
    rand -= log(RandomUniform());
  return rand / rate;
}

double Weibull(double rate, double shape)
{
  return pow(-log(RandomUniform())/rate, 1.0/shape);
}

double NormalDensity(double dev, double sd)
{
  const double c=0.3989;
  
  dev = dev/sd;
  return exp(-dev*dev/2.0)*c/sd;
}

#endif