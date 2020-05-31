#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>
using namespace std;

const double PI = 3.14;

using func = function<double(double)>;                     // Function object

double integral( func f, double a, double b, int n );


//==========================================================


class FourierSeries
{
   double period;                                          // Period
   double omega;                                           // Fundamental frequency
   int nmax;                                               // Highest harmonic
   vector<double> a, b;
public:
   FourierSeries( double (*f)( double ), double T, int N ) : period( T ), nmax( N )
   {
      omega = 2.0 * PI / period;
      a = vector<double>( 1 + nmax, 0 );   b = a;
      int Ndt = max( 10 * N, 100 );                        // Number of intervals for numerical integration
      double factor = 2.0 / period;                        // Normalising factor for basis functions

      for ( int n = 0; n <= nmax; n++ )                    // Compute coefficients
      {
         func c = [=]( double t ){ return f( t ) * cos( n * omega * t ); };
         func s = [=]( double t ){ return f( t ) * sin( n * omega * t ); };
         a[n] = factor * integral( c, 0.0, period, Ndt );
         b[n] = factor * integral( s, 0.0, period, Ndt );
      }
   }

   double evaluate( double t )                                       // Sum Fourier series
   {
      double sum = 0.5 * a[0];
      for ( int n = 1; n <= nmax; n++ ) sum += a[n] * cos( n * omega * t ) + b[n] * sin( n * omega * t );
      return sum;
   }

   void tabulate( ostream &out, double t0, double dt, int nt )       // Output
   {
      #define FMT << fixed << setprecision( 4 ) << setw( 15 ) <<
      for ( int i = 0; i <= nt; i++ )
      {
         double t = t0 + i * dt;
         out FMT t FMT evaluate( t ) << '\n';
      }
   }
};


//==========================================================


double integral( func f, double a, double b, int n )                 // Numerical integration: mid-ordinate rule
{
   double dt = ( b - a ) / n;
   double sum = 0;
   for ( int i = 1; i <= n; i++ ) sum += f( a + ( i - 0.5 ) * dt) * dt;
   return sum;
}


//==========================================================


double squareWave( double t )
{
   t = t - floor( t );
   return t < 0.5 ? 1 : 0;
}


double sawTooth( double t )
{
   t = t - floor( t );
   return t;
}


//==========================================================


int main()
{
   const double period = 1.0;                    // Period of function
   const int N = 20;                             // Highest harmonic

   FourierSeries fs( squareWave, period, N );    // Try other functions besides squareWave
   fs.tabulate( cout, 0.0, period / 100, 400 );
}


//==========================================================
