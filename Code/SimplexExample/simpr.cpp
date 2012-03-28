void simpr(REAL* y, const int n, REAL* dydx, REAL* dfdx, REAL* dfdy,
	const REAL xs, const REAL htot, const int nstep, REAL* yout,
	void derivs(const REAL, REAL*, REAL*))
{
  int i,j,nn;
  REAL d,h,x;
  
  REAL* a[n*n];
  int* indx[n];
  REAL* del[n],ytemp[n];
  h=htot/nstep;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) a[i][j] = -h*dfdy[i][j];
    ++a[i][i];
  }
  ludcmp(a,indx,d);
  for (i=0;i<n;i++)
    yout[i]=h*(dydx[i]+h*dfdx[i]);
  lubksb(a,indx,yout);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+(del[i]=yout[i]);
  x=xs+h;
  derivs(x,ytemp,yout);
  for (nn=2;nn<=nstep;nn++) {
    for (i=0;i<n;i++)
      yout[i]=h*yout[i]-del[i];
    lubksb(a,indx,yout);
    for (i=0;i<n;i++) ytemp[i] += (del[i] += 2.0*yout[i]);
    x += h;
    derivs(x,ytemp,yout);
  }
  for (i=0;i<n;i++)
    yout[i]=h*yout[i]-del[i];
  lubksb(a,indx,yout);
  for (i=0;i<n;i++)
    yout[i] += ytemp[i];
}
