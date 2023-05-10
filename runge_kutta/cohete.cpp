#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
#define PI 3.14159265

#define G 6.67e-11 //N m^2 kg^-2
#define MT 5.9736e24 //kg   //masa Tierra
#define ML 0.07349e24 //kg  //masa Luna
#define dTL 3.844e8 //m     //distancia Tierra-Luna
#define w 2.6617e-6 //s^-1  //Velocidad angular Luna??????
#define RT 6.378160e6 //m   //Radio Tierra
#define RL 1.7374e6 //m     //Radio Luna

int main(){
    //------------------------DECLARACION DE VARIABLES-------------------------------------
    double y[4]; //funciones a resolver de las ecuaciones diferenciales
    //   y_0 = r     y_1 = phi     y_2 = p_r     y_3 = p_phi
    float t; //Tiempo
    //-------------------------------------------------------------------------------------
    
    //--------------------------CONDICIONES INICIALES--------------------------------------
    //Aplicamos reescalado a las variables
    y[0]=RT/dTL; //El cohete empieza en el radio de la Tierra
    y[1]=PI/4.0; //Ángulo inicial del cohete
    y[2]=11200.0/dTL; //Componente radial velocidad inicial = velocidad de escape terrestre = 11.2 km/s
    y[3]=0; //Componente angular velocidad inicial nula => Sale normal a la superficie terrestre
    t=0;
    //-------------------------------------------------------------------------------------

    //-------------------------------RUNGE_KUTTA-------------------------------------------
    
    //-------------------------------------------------------------------------------------
}