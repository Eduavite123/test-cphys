#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

#define PI 3.14159265w

using namespace std;

int main(void){
    //------------------------DECLARACION DE VARIABLES-------------------------------------
    double h, s; //discretización: espaciado espacial y temporal, respectivamente
    double phi[][100]; //¡¡¡Dim máxima por determinar!!!!!!!!!!!!!!!!!!!
    //Valores iniciales
    float N, nc; //tamaño pozo y nº ciclos (nc=1,...,N/4 para que haya 4 ciclos)
    float lambda; //proporcional a la altura del pozo (p.e. lambda=0.3)
    //Valores iniciales a generar
    double S; 
    double k_0; 
    double V; //altura de la pared
    //-------------------------------------------------------------------------------------

    //--------------------------CONDICIONES INICIALES--------------------------------------
    h=0.1; //
    s=0.1; //¡¡¡A DETERMINAR!!!!

    N=4; //¡¡¡A DETERMINAR!!!!
    nc=N/4; //lo elegimos así para que haya 4 ciclos

    lambda=0.3;

    k_0=2*PI*nc/N;
    S=1/(4*k_0*k_0);
    //hacer V en función pot
    //-------------------------------------------------------------------------------------
    return 0;
}