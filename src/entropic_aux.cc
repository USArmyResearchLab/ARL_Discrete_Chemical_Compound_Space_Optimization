/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 3; tab-width: 3 -*- */
/*
 * This file is part of ARL Discrete Chemical Compound Space Optimization (ARL DCCSO) project.
 *
 * ARL DCSSO constitutes a work of the United States Government and is not
 * subject to domestic copyright protection under 17 USC Sec. 105.
 * Release authorized by the US Army Research Laboratory
 *
 * To the extent possible under law, the author(s) have dedicated all copyright
 * and related and neighboring rights to this software to the public domain
 * worldwide. This software is distributed without any warranty.
 *
 * You should have received a copy of the CC0 Public Domain Dedication along
 * with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
 *
 */

/*! @addtogroup DOChemS Discrete Optimization of Chemical Space */
//! \file entropic_aux.hh Auxiliary functions for *_entropic classes.

#ifndef _ENTROPIC_AUX_CC_
#define _ENTROPIC_AUX_CC_
#include <iostream>
#include <cmath>
#include <BCR_CPP_LA/linear_algebra.h>
#include <entropic_aux.hh>

using namespace std;
using namespace linear_algebra;

//! Conjugate gradients for linear solver.
void linsolve_cg(const mat_sym_full<double>& J,
      const refvector<double>& G,
      refvector<double>& X)
{
   // residual
   refvector<double> r(J*G);
   r-=J*(J*X);
   // conjugate vector
   refvector<double> c;
   c.copy(r);


   double n;
   double n2;
   double a;
   double b;
   double error=1.0;
   n=r*r;
   error=(J*c)*(J*c);
   while (n > 1e-16 && error > 1e-16)
   {
      a=n/error;
      X+=c*a;
      r-=(J*(J*c))*a;
      n2=r*r;
      b=n2/n;
      n=n2;
      c*=b;
      c+=r;
      error=(J*c)*(J*c);
   }

   return;
}

/** Compute the gradient of the distance function \f$ f=\sum_i \ln \sqrt{\sum_j \sin (x_j-s^{(i)}_j\pi/n_j)^2} \f$
 * where \f$n_j\f$ is the number of digits in component j.
 */
refvector<double>& set_gradient(const refvector<long>& b,
      const mat_full<double>& H,
      const refvector<double>& X,
      refvector<double>& G)
      //! We ignore infinities in the derivatives
      {
   string serr="binary_entropic: static refvector<double>& set_gradient(const refvector<long>& b,\
                                                                        const mat_full<double>& H,\
                                                                        const refvector<double>& X,\
                                                                        refvector<double>& G):";
   if(X.dim()!=H.rows())
      throw domain_error(serr+" X and H have incompatible dimensions.");
   if(G.dim()!=H.rows())
      G.resize(H.rows());
   G.zero();
   long i,j;
   bool redundancy=false;
   for(i=0;i<H.cols();i++)
   {
      double n=0;
      for(j=0;j<G.dim();j++)
         n+=sin((X[j]-H[i][j])*M_PI/(double) b[j])*
         sin((X[j]-H[i][j])*M_PI/(double) b[j]);
      if(n>1e-16)
      {
         n=(double) 1/n;
         for(j=0;j<G.dim();j++)
            G[j]-=(
                  cos((X[j]-H[i][j])*M_PI/(double) b[j])*
                  sin((X[j]-H[i][j])*M_PI/(double) b[j])*n*
                  M_PI/(double) b[j]
            );
      }
      else {
         redundancy=true;
         cout << "WARNING: set_gradient has encountered a redundancy: ";
         X.display();
         cout << " at index H[" << i << "]" << endl;
      }
   }

   if(redundancy && G*G<1e-16)
      for(i=0;i<G.dim();i++)
         G[i]-=1.0;
   return G;
      }

/** Compute the hessian of the distance function \f$ f=\sum_i \ln \sqrt{\sum_j \sin (x_j-s^{(i)}_j\pi/n_j)^2} \f$
 * where \f$n_j\f$ is the number of digits in component j.
 */
mat_sym_full<double>& set_hessian(const refvector<long>& b,
      const mat_full<double>& H,
      const refvector<double>& X,
      mat_sym_full<double>& J)
      //! We ignore infinities in the derivatives
      {
   string serr="binary_entropic: static refvector<double>& set_hessian(const refvector<lob>& b,\
                                                                        const mat_full<double>& H,\
                                                                        const refvector<double>& X,\
                                                                        refvector<double>& J):";
   if(X.dim()!=H.rows())
      throw domain_error(serr+" X and H have incompatible dimensions.");
   if(J.cols()!=H.rows())
      throw domain_error(serr+" J and H have incompatible dimensions.");
   J.zero();
   long i,j,k;
   for(i=0;i<H.cols();i++)
   {
      double n=0;
      for(j=0;j<H.rows();j++)
         n+=sin((X[j]-H[i][j])*M_PI/(double) b[j])*
         sin((X[j]-H[i][j])*M_PI/(double) b[j]);
      if(n>1e-16)
      {
         double n_1;
         n_1=(double) 1/n;
         for(j=0;j<J.cols();j++)
         {
            for(k=0;k<=j;k++)
               J[j*(j+1)/2+k]+=(
                     cos((X[j]-H[i][j])*M_PI/(double) b[j])*
                     sin((X[j]-H[i][j])*M_PI/(double) b[j])*n_1*
                     M_PI/(double) b[j]*

                     cos((X[k]-H[i][k])*M_PI/(double) b[k])*
                     sin((X[k]-H[i][k])*M_PI/(double) b[k])*n_1*
                     M_PI/(double) b[k]
               );/*
               J[j*(j+1)/2+k]-=(
                                cos((X[j]-H[i][j])*M_PI/(double) b[j])*
                                sin((X[j]-H[i][j])*M_PI/(double) b[j])*
                                M_PI/(double) b[j]*

                                cos((X[k]-H[i][k])*M_PI/(double) b[k])*
                                sin((X[k]-H[i][k])*M_PI/(double) b[k])*
                                M_PI/(double) b[k]
                                )*(n_1+(double) 1/((double) b.size()-n));*/
            J[j*(j+1)/2+j]-=(
                  cos((X[j]-H[i][j])*M_PI/(double) b[j])*
                  cos((X[j]-H[i][j])*M_PI/(double) b[j])*n_1*
                  M_PI/(double) b[j]*M_PI/(double) b[j]
                                                     -
                                                     sin((X[j]-H[i][j])*M_PI/(double) b[j])*
                                                     sin((X[j]-H[i][j])*M_PI/(double) b[j])*n_1*
                                                     M_PI/(double) b[j]*M_PI/(double) b[j]
            )*0.5;/*
            J[j*(j+1)/2+j]-=(
                             cos((X[j]-H[i][j])*M_PI/(double) b[j])*
                             cos((X[j]-H[i][j])*M_PI/(double) b[j])*
                             M_PI/(double) b[j]*M_PI/(double) b[j]
                             -
                             sin((X[j]-H[i][j])*M_PI/(double) b[j])*
                             sin((X[j]-H[i][j])*M_PI/(double) b[j])*
                             M_PI/(double) b[j]*M_PI/(double) b[j]
                             )*0.5*(log(n)-log((double) b.size()-n));*/
         }
      }
   }
   J*=2.0;
   return J;
      }
//! Create the neighbor if distance 1 around a point X
refvector<refvector<ulong> > create_neighborhood(const refvector<ulong>& Y, const refvector<long>& bases)
{
#ifdef DEBUG
   if(Y.size() != bases.size())
      throw domain_error("create_neighborhood: X and bases have differing dimensions");
#endif
   refvector<ulong> X(Y);
   long count=0;
   for(long dim=0; dim<bases.size();dim++)
      if(bases[dim]>1) count++;
   refvector<refvector<ulong> > N(2*count);
   count=0;
   for(long dim=0;dim<X.size(); dim++)
      if(bases[dim]>1) {
         X[dim]=(X[dim]+1) % bases[dim];
         N[count].copy(X);
         count++;
         if (X[dim]<2) X[dim]=X[dim]-2+bases[dim];
         else X[dim]=X[dim]-2;
         N[count].copy(X);
         X[dim]=(X[dim]+1) % bases[dim];
         count++;
      }
   return N;
}
//! Create the neighbor if distance 1 around a point X
refvector<refvector<ulong> > create_neighborhood(const refvector<refvector<ulong> >& X, const refvector<long> bases)
{
   refvector<refvector<ulong> > N;
   for(long i=0; i<X.size(); i++) {
      N.concat(create_neighborhood(X[i],bases));
   }
   return N;
}

template <class T>
refvector<T> without(const refvector<T> &A,
                     const refvector<T> &B)
{
   refvector<T> C;
   for(long i=0;i<A.size();i++)
      if(B.contains(A[i])<0)
         C.push_back(A[i]);
   return C;
}
template <class T>
refvector<T> Union(const refvector<T> &A,
                   const refvector<T> &B)
{
   refvector<T> C(B);
   for(long i=0;i<A.size();i++)
      if(B.contains(A[i])<0)
         C.push_back(A[i]);
   return C;
}

//! Find a point external to the visited set
refvector<refvector<ulong> > find_external_point(const refvector<refvector<ulong> >& Gamma, const refvector<long>& bases)
{
#ifdef DEBUG
   if(Gamma.size()<1)
      throw domain_error("find_external_point: Gamma is empty");
#endif
   refvector<refvector<ulong> > Gamma2;
   refvector<refvector<ulong> > dGamma;
   refvector<refvector<ulong> > N;
   refvector<refvector<ulong> > R;
   refvector<ulong> X(Gamma[0]);
   dGamma.push_back(X);
   N=create_neighborhood(dGamma,bases);
   R = without(N,Gamma);
   while(R.size()==0)
   {
      Gamma2 = Union(dGamma,Gamma2);
      dGamma = without(N,Gamma2);
      N = create_neighborhood(dGamma,bases);
      R = without(N,Gamma);
   }
   return R;
}

//! This is the sin metric used \f$m(x,Y) = -0.5\sum_{y\in Y} \ln \sum_{d=1}^{N_d} \sin^2 \frac{\pi(x_d-y_d)}{N_d}/\sum N_d\f$
double dmetric_lnsin(const refvector<refvector<ulong> >& Y, const refvector<double>& x, const refvector<long>& bases)
{
   double metric=0.0;
   for(long yi=0;yi<Y.dim();yi++)
   {
      double sumsin=0.0;
      for(long d=0;d<bases.dim();d++)
         sumsin+=pow(sin((x[d]-(double) Y[yi][d])*M_PI/(double) bases[d]),2);
      metric-=0.5*log(sumsin/(double) bases.dim());
   }
   return metric;
}

//! This is the sin metric used \f$m(x,Y) = -0.5\sum_{y\in Y} \ln \sum_{d=1}^{N_d} \sin^2 \frac{\pi(x_d-y_d)}{N_d}/\sum N_d\f$
double metric_lnsin(const refvector<refvector<ulong> >& Y, const refvector<ulong>& x, const refvector<long>& bases)
{
   double metric=0.0;
   for(long yi=0;yi<Y.dim();yi++)
   {
      double sumsin=0.0;
      for(long d=0;d<bases.dim();d++)
         sumsin+=pow(sin(((double) x[d]-(double) Y[yi][d])*M_PI/(double) bases[d]),2);
      metric-=0.5*log(sumsin/(double) bases.dim());
   }
   return metric;
}

//! Return the index of X that has the lowest metric_lnsin with respect to Y
long argmin_lnsin(const refvector<refvector<ulong> >& Y, const refvector<refvector<ulong> >& X, const refvector<long>& bases)
{
   long minindex=0;
   double minvalue=metric_lnsin(Y,X[0],bases);
   cout << "argmin_lnsin: " << minvalue;
   for(long xi=1;xi<X.dim();xi++)
   {
      double value = metric_lnsin(Y,X[xi],bases);
      cout << " " << value;
      if(value < minvalue) {
         minvalue = value;
         minindex = xi;
      }
   }
   cout << endl;
   return minindex;
}

//! Local search around x
refvector<ulong> argmin_lnsin(const refvector<refvector<ulong> >& Y, const refvector<double>& x, const refvector<long>& bases)
{
   refvector<long> xu(bases.dim());
   refvector<long> xl(bases.dim());
   mat_full<long> xb(2,bases.dim());
   for(long i=0; i< bases.dim(); i++) xb[0][i]=(long) ceil(x[i]);
   for(long i=0; i< bases.dim(); i++) xb[1][i]=(long) floor(x[i]);
   long ncandidates=(long) pow(2,bases.dim());
   refvector<refvector<ulong> > candidates(ncandidates);
   for(long i=0; i<ncandidates; i++)
   {
      refvector<ulong> c(bases.dim());
      long d=i;
      for(long j=0;j<bases.dim();j++)
      {
         c[j]=((xb[d % 2][j]  % bases[j]) + bases[j]) % bases[j];
         d=(d-(d%2))/2;
      }
      candidates[i]=c;
   }
   return candidates[argmin_lnsin(Y,candidates,bases)];
}

//! Maximize the entropic distance as declared in set_gradient() using Newton-Raphson.
ulong maximize_entropic_distance(const refvector<ulong>& A, const refvector<long>& b)
/*!
 * \f$ d(x,Y) = \sum_i \ln \sqrt{\sum_j \sin (x_j-Y^{(i)}_j2\pi/n_j)^2} \f$
 */
{
   long i,k;
   mat_full<double> H(A.size(),b.size());
   refvector<refvector<ulong> > ulH(A.size());
   refvector<refvector<ulong> > ulR;
   mat_sym_full<double> J(b.size());
   refvector<double> G(b.size());
   refvector<double> X(b.size());
   refvector<ulong> ulX(b.size());

   for(i=0;i<H.cols();i++)
   {
      long l=1;
      refvector<ulong> dummy(b.size());
      ulH[i]=dummy;
      for(k=0;k<b.size();k++)
      {
         ulH[i][k]=((A[i] -A[i]%l)/l) % b[k];
         H[i][k]=ulH[i][k];
         l*=b[k];
      }
   }
   double error=1.0;
   //X.copy(H[H.cols()-1]);
   ulR = find_external_point(ulH,b);
   ulX = ulR[argmin_lnsin(ulH,ulR,b)];
   cout << "Starting point for entropic search : "; 
   ulX.display();
   for(k=0;k<ulX.dim();k++) X[k]=ulX[k];
   cout << " == ";
   X.display();
   cout << endl;

   double current_metric=dmetric_lnsin(ulH,X,b);
   cout << "Starting metric = " << current_metric << " compared to " << metric_lnsin(ulH,ulX,b) << endl;
   while (error>1e-16)
   {
      G=set_gradient(b,H,X,G);
      J=set_hessian(b,H,X,J);
      //G-=J*X;
      G*=-1.0;

      refvector<double> residual(X.dim());
      linsolve_cg(J,G,residual);

      while(current_metric+2e-16 < dmetric_lnsin(ulH,X+residual,b))
      {
          residual*=0.5;
      }
      error=residual*residual;
      if( error < 2e-16)
      {
          residual = G;
          while(current_metric+2e-16 < dmetric_lnsin(ulH,X+residual,b))
          {
             residual*=0.5;
          }
          error=residual*residual;
      }
      
      X+=residual;
      current_metric = dmetric_lnsin(ulH,X,b);
      cout << "Current metric = " << current_metric << " ";
      X.display();
      cout << " Gradient "; G.display();
      cout << endl;
   }
   ulong conf1=0;
   ulong m=1;
   for(i=0;i<b.size();i++)
   {
      long Xi = (long) lround(X[i]);
      cout << " X[i] rounded = " << Xi << " mod " << b[i] << " = " << (Xi % b[i]) <<  " * " << m << endl;
      Xi = Xi % b[i];
      if(Xi>=0)
         conf1+=Xi*m;
      else
         conf1+=(Xi + b[i])*m;
      m*=b[i];
   }
   cout << endl << conf1 << " ";
   ulX=argmin_lnsin(ulH,X,b);
   m=1;
   conf1=0;
   for(i=0;i<b.size();i++)
   {
      long Xi = (long) ulX[i];
      cout << " X[i] rounded = " << Xi << " mod " << b[i] << " = " << (Xi % b[i]) <<  " * " << m << endl;
      Xi = Xi % b[i];
      if(Xi>=0)
         conf1+=(Xi )*m;
      else
         conf1+=(Xi + b[i])*m;
      m*=b[i];
   }
   cout << endl << conf1 << endl;

   double entropy=0.0;
   for(i=0;i<H.cols();i++)
   {
      double norm=0.0;
      for(m=0;m<H.rows();m++)
         norm+=(
               sin((X[m]-H[i][m])*M_PI/(double) b[m])*
               sin((X[m]-H[i][m])*M_PI/(double) b[m])
         );
      if(norm>=1e-16*(double) b.size())
         entropy-=norm*log(norm/(double) b.size())+((double) b.size() -norm)*log(1.0-norm/(double) b.size());
   }
   // Compute the Manhattan distance to all previous structures
   ulong mindist=0;
   for(i=0;i<b.size();i++) mindist+=b[i];
   for(k=0;k<ulH.dim();k++)
   {
      long dist=0;
      for(i=0;i<b.size();i++)
      {
         long Xi = (long) lround(X[i])-(long) ulH[k][i];
         Xi = Xi % b[i];
         cout << " Xi(" << k << "," << i << ") = " << Xi << " mod " << b[i] << " dist = " << dist;
         if(Xi>=0) {
           if(Xi <= b[i]*0.5)
              dist+=Xi;
           else
            dist+=b[i]-Xi;
         }
         else {
           if(Xi<-b[i]*0.5)
              dist+=b[i]+Xi;
           else
            dist+=-Xi;
         }
      }
      cout << " dist = " << dist << " mindist = " << mindist << endl;
      if(dist<mindist) mindist=dist;
   }
   mindist=0;
   for(i=0;i<b.size();i++) mindist+=b[i];
   for(k=0;k<ulH.dim();k++)
   {
      long dist=0;
      for(i=0;i<b.size();i++)
      {
         long Xi = (long) ulX[i]-(long) ulH[k][i];
         Xi = Xi % b[i];
         cout << " Xi(" << k << "," << i << ") = " << Xi << " mod " << b[i] << " dist = " << dist;
         if(Xi>=0) {
           if(Xi <= b[i]*0.5)
              dist+=Xi;
           else
            dist+=b[i]-Xi;
         }
         else {
           if(Xi<-b[i]*0.5)
              dist+=b[i]+Xi;
           else
            dist+=-Xi;
         }
      }
      cout << " dist = " << dist << " mindist = " << mindist << endl;
      if(dist<mindist) mindist=dist;
   }

   cout << "Entropy: " << entropy/(double) b.size() << " Minimum Manhattan distance to set: " << mindist << endl;
   return conf1;
}
#endif // _ENTROPIC_AUX_CC_
