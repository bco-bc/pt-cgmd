#include "simploce/util/bessel.hpp"
#include <boost/math/special_functions/bessel.hpp>

using namespace simploce::math;
using namespace std;

int main() {
    double x, dx=0.1, l_error, abs_error;
    char space=' ';

    cout.setf(ios::fixed);
    cout.unsetf(ios::scientific);

// Checking quality of series expansion for Jn(z).
// See AS, Table 9.1, p.390, Table 9.2, p. 398.
    cout << "--SERIES EXPANSION FOR Jn(z) WITH REAL NUMBERS---------------------------------\n";
    cout << "   x              J0(x)         J1(x)         J1(2)         J3(x)         J4(x)\n";
    for (int i=0; i<=100; i++) {
         x=i*dx;
         cout << setprecision(1) << setw(4) << x;
         cout.flush();
         double jn;
         for (int n=0; n<5; n++) {
              cout << space;
              jn=B1_Jn_S(n, x);
              if (n==0)
                   cout << setprecision(15) << setw(18) << jn;
              else
                   cout << setprecision(10) << setw(13) << jn;
         }
         cout << "\n";
         cout.flush();
    }

// Checking polynomial approximation for J0(x).
    cout << "\n";
    cout << "--POLYNOMIAL APPROXIMATION FOR J0(x) FROM A&S 9.4.1----------------------\n";
    cout << "   x                 Series            Polyn. App.             Abs. Error\n";
    l_error=0;
    for (int j=-30; j<=30; j++) {
         cout.unsetf(ios::scientific);
         cout.setf(ios::fixed);
         x=j*dx;
         cout << setprecision(1) << setw(4) << x << space;
         cout.flush();
         double jn1, jn2;
         jn1=B1_Jn_S(0, x);    // Exact result.
         cout.setf(ios::scientific);
         cout.unsetf(ios::fixed);
         cout << setprecision(15) << setw(22) << jn1 << space;
         jn2=B1_J0_PA_AS941(x);
         abs_error=fabs(jn1-jn2);
         if (abs_error > l_error) l_error=abs_error;
         cout << setprecision(15) << setw(22) << jn2 << space << setprecision(15) << setw(22) << abs_error << "\n";
         cout.flush();
    }
    cout << "Largest error: " << setprecision(15) << setw(22) << l_error << "\n";

// Checking polynomial approximation for J0(x).
    cout << "\n";
    cout << "--POLYNOMIAL APPROXIMATION FOR J0(x) FROM A&S 9.4.3----------------------\n";
    cout << "   x                 Series            Polyn. App.             Abs. Error\n";
    l_error=0;
    dx=0.1;
    for (int j=30; j<=80; j++) {
         cout.unsetf(ios::scientific);
         cout.setf(ios::fixed);
         x=j*dx;
         cout << setprecision(1) << setw(4) << x << space;
         cout.flush();
         double jn1, jn2;
         jn1=B1_Jn_S(0, x);             // Exact result from series expansion .
         cout.setf(ios::scientific);
         cout.unsetf(ios::fixed);
         cout << setprecision(15) << setw(22) << jn1 << space;
         jn2=B1_J0_PA_AS943(x);
         abs_error=fabs(jn1-jn2);
         if (abs_error > l_error) l_error=abs_error;
         cout << setprecision(15) << setw(22) << jn2 << space << setprecision(15) << setw(22) << abs_error << "\n";
         cout.flush();
    }
    cout << "Largest error: " << setprecision(15) << setw(22) << l_error << "\n";

    // Boost implementation J0(x)
    cout << "\n";
    cout << "--BOOST IMPLEMENTATION FOR J0(x) ----------------------------------------\n";
    cout << "   x                 Series                  Boost             Abs. Error\n";
    l_error=0;
    dx=0.1;
    for (int j=30; j<=80; j++) {
         cout.unsetf(ios::scientific);
         cout.setf(ios::fixed);
         x=j*dx;
         cout << setprecision(1) << setw(4) << x << space;
         cout.flush();
         double jn1, jn2;
         jn1=B1_Jn_S(0, x);             // Exact result from series expansion .
         cout.setf(ios::scientific);
         cout.unsetf(ios::fixed);
         cout << setprecision(15) << setw(22) << jn1 << space;
         jn2=boost::math::cyl_bessel_j(0, x);
         abs_error=fabs(jn1-jn2);
         if (abs_error > l_error) l_error=abs_error;
         cout << setprecision(15) << setw(22) << jn2 << space << setprecision(15) << setw(22) << abs_error << "\n";
         cout.flush();
    }
    cout << "Largest error: " << setprecision(15) << setw(22) << l_error << "\n";

    // MacDonald's function K0(x) and K1(x), see Table on page 417 of AS.
    cout << "\n";
    cout << "--MacDonald's K0(x) and K1(x)-----------------------------------------------------\n";
    cout << "   x               e^z*K0(x)             e^x*K1(x)                                \n";
    l_error=0;
    dx=0.1;
    for (int j=1; j<=100; j++) {
         cout.unsetf(ios::scientific);
         cout.setf(ios::fixed);
         x=j*dx;
         cout << setprecision(1) << setw(4) << x << space;
         cout.flush();
         double jn1, jn2;
         jn1= Bessel_K0(x);
         jn1 *= std::exp(x);
         cout.setf(ios::scientific);
         cout.unsetf(ios::fixed);
         cout << setprecision(15) << setw(22) << jn1 << space;
         jn2=Bessel_K1(x);
         jn2 *= std::exp(x);
         cout.setf(ios::scientific);
         cout.unsetf(ios::fixed);
         cout << setprecision(15) << setw(22) << jn2 << space << std::endl;
    }
}