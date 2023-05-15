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

#define Delta G*MT/pow(dTL,3)
#define mu ML/MT

double f0(double p_r);
double f1(double r, double p_phi);
double f2(double r, double phi, double p_phi, float t);
double f3(double r, double phi, float t);

int main(){
    //------------------------DECLARACION DE VARIABLES-------------------------------------
    double y[4]; //funciones a resolver de las ecuaciones diferenciales
    //   y_0 = r     y_1 = phi     y_2 = p_r     y_3 = p_phi
    float t, h, tmax, iter; //Tiempo y pasos temporales | Tiempo máx y nº iteraciones
    double k1[4], k2[4], k3[4], k4[4];
    double r_prima;
    double aux[4];
    double x_c,y_c; //Posición del cohete en coordenadas cartesianas
    double x_L,y_L; //Posición de la Luna en coordenadas cartesianas
    //-------------------------------------------------------------------------------------
    
    //--------------------------CONDICIONES INICIALES--------------------------------------
    //Paso de tiempo
    h=0.1;
    tmax=50;
    iter=tmax/h;

    //Aplicamos reescalado a las variables
    y[0]=RT/dTL; //El cohete empieza en el radio de la Tierra
    y[1]=PI/4.0; //Ángulo inicial del cohete
    y[2]=11200.0/dTL; //Componente radial velocidad inicial = velocidad de escape terrestre = 11.2 km/s
    y[3]=0; //Componente angular velocidad inicial nula => Sale normal a la superficie terrestre
    t=0;

    //Condiciones iniciales de la Luna
    x_L=1.0*cos(0.0);
    y_L=1.0*sin(0.0);
    //-------------------------------------------------------------------------------------

    //-------------------------------RUNGE_KUTTA-------------------------------------------
    ofstream pos;
    pos.open("cohete_data.dat");
    for (int i = 0; i < iter; i++)
    {
        //Pasamos de coordenadas cilíndricas a cartesianas
        x_c=y[0]*cos(y[1]);
        y_c=y[0]*sin(y[1]); 

        //Escribimos los datos en el archivo
        pos << 0.0 << ", " << 0.0 << endl; //Posición Tierra (fija)
        pos << x_L << ", " << y_L << endl; //Posición Luna
        pos << x_c << ", " << y_c << endl;
        pos << endl;

        //Evaluamos k1
        k1[0]=h*f0(y[2]);
        k1[1]=h*f1(y[0],y[3]);
        k1[2]=h*f2(y[0],y[1],y[3],t);
        k1[3]=h*f3(y[0],y[1],t);

        //Un vector auxiliar facilitará el cálculo
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k1[j]/2.0;
        }

        //Evaluamos k2
        k2[0]=h*f0(aux[2]);
        k2[1]=h*f1(aux[0],aux[3]);
        k2[2]=h*f2(aux[0],aux[1],aux[3],t+h/2.0);
        k2[3]=h*f3(aux[0],aux[1],t+h/2.0);

        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k2[j]/2.0;
        }

        //Evaluamos k3
        k3[0]=h*f0(aux[2]);
        k3[1]=h*f1(aux[0],aux[3]);
        k3[2]=h*f2(aux[0],aux[1],aux[3],t+h/2.0);
        k3[3]=h*f3(aux[0],aux[1],t+h/2.0);        
        
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k3[j];
        }

        //Evaluamos k4
        k4[0]=h*f0(aux[2]);
        k4[1]=h*f1(aux[0],aux[3]);
        k4[2]=h*f2(aux[0],aux[1],aux[3],t+h);
        k4[3]=h*f3(aux[0],aux[1],t+h);   

        //Finalmente evaluamos y_n(t+h)
        for (int k = 0; k < 4; k++)
        {
            y[k]=y[k]+1.0/6.0*(k1[k]+2*k2[k]+2*k3[k]+k4[k]);
        }

        t=t+h; 

        //Movemos la Luna
        x_L=1.0*cos(w*t);
        y_L=1.0*sin(w*t);         

    }
    pos.close();
    //-------------------------------------------------------------------------------------
    return 0;
}

double f0(double p_r){
    double f0;
    f0=p_r;
    return f0;
}

double f1(double r, double p_phi){
    double f1;
    f1=p_phi/pow(r,2);
    return f1;
}

double f2(double r, double phi, double p_phi, float t){
    double f2, r_prima;
    r_prima=sqrt(1+pow(r,2)-2*r*cos(phi-w*t));
    f2=pow(p_phi,2)/pow(r,3)-Delta*(1/pow(r,2)+mu/pow(r_prima,3)*(r-cos(phi-w*t)));
    return f2;
}

double f3(double r, double phi, float t){
    double f3, r_prima;
    r_prima=sqrt(1+pow(r,2)-2*r*cos(phi-w*t));
    f3=-Delta*mu*r*sin(phi-w*t)/pow(r_prima,3);
    return f3;
}