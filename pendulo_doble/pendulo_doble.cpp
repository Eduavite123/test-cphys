#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
#define g 9.81
#define PI 3.14159265

double f0(double O1, double O2, double p_o1, double p_o2);
double f1(double O1, double O2, double p_o1, double p_o2);
double f2(double O1, double O2, double p_o1, double p_o2);
double f3(double O1, double O2, double p_o1, double p_o2);

int main(){
    //------------------------DECLARACION DE VARIABLES-------------------------------------
    double y[4]; //funciones a resolver de las ecuaciones diferenciales
    //   y_0 = \theta_1     y_1 = \theta_2     y_2 = p_{\theta_1}  y_3 = p_{\theta_2}
    float t, h, tmax, iter; //Tiempo y pasos temporales | Tiempo máx y nº iteraciones
    double k1[4], k2[4], k3[4], k4[4];
    double aux[4];
    double dO1, dO2; //Velocidades angulares
    double H; //Hamiltoniano
    //-------------------------------------------------------------------------------------

    //--------------------------CONDICIONES INICIALES--------------------------------------
    //Paso de tiempo y tiempo máximo
    h=0.1;
    tmax=50;
    iter=tmax/h;

    //Posiciones y momentos iniciales de las masas
    y[0]=PI/2.0;
    y[1]=PI/2.0;

    //Energía fija del sistema
    H=20.0;
    if (H-2*g*(1-cos(y[0]))-g*(1-cos(y[2]))<0)
    {
        cout << "La velocidad angular inicial no es un número real. Por favor, elija una energía más alta." << endl;
    }

    //Velocidades angulares 
    dO1=sqrt(H-2*g*(1-cos(y[0]))-g*(1-cos(y[2])));
    dO2=0.0;  

    //Momentos iniciales
    y[2]=(2*(1+pow(sin(y[0]-y[1]),2))*dO1)/(cos(y[0]-y[1]));
    y[3]=dO1*(1+pow(sin(y[0]-y[1]),2));
    //-------------------------------------------------------------------------------------

    //-------------------------------RUNGE_KUTTA-------------------------------------------
    ofstream pos;
    ofstream poin;
    ofstream poin1;
    ofstream poin2;
    ofstream lya;
    pos.open("posiciones.dat");
    poin.open("poincare.dat");
    poin1.open("poincare_O2_dO2.dat");
    poin2.open("poincare_O2_dO1.dat");
    lya.open("lyapunov.txt");
    for (int i = 0; i < iter; i++)
    {
        
        //Escribimos los datos en el archivo
        //Simulación péndulo doble
        pos << sin(y[0]) << ", " << -cos(y[0]) << endl;
        pos << sin(y[0])+sin(y[1]) << ", " << -cos(y[0])-cos(y[1]) << endl;
        pos << endl; 
        //Mapa de Poincaré para ángulos
        poin << y[0] << ", " << y[1] << endl;
        poin << endl;
        //Mapa de Poincaré para O2 y dO2
        poin1 << y[1] << ", " << f1(y[0], y[1], y[2], y[3]) << endl;
        poin1 << endl;
        //Mapa de Poincaré para O2 y dO1
        poin2 << y[1] << ", " << f0(y[0], y[1], y[2], y[3]) << endl;
        poin2 << endl;
        //Exponentes de Lyapunov
        lya << y[1] << "    " << f1(y[0], y[1], y[2], y[3]) << endl;

        //Evaluamos k1
        k1[0]=h*f0(y[0], y[1], y[2], y[3]);
        k1[1]=h*f1(y[0], y[1], y[2], y[3]);
        k1[2]=h*f2(y[0], y[1], y[2], y[3]);
        k1[3]=h*f3(y[0], y[1], y[2], y[3]);

        //Un vector auxiliar facilitará el cálculo
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k1[j]/2.0;
        }

        //Evaluamos k2
        k2[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k2[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k2[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k2[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);

        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k2[j]/2.0;
        }

        //Evaluamos k3
        k3[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k3[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k3[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k3[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);        
        
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k3[j];
        }

        //Evaluamos k4
        k4[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k4[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k4[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k4[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);   

        //Finalmente evaluamos y_n(t+h)
        for (int k = 0; k < 4; k++)
        {
            y[k]=y[k]+1.0/6.0*(k1[k]+2*k2[k]+2*k3[k]+k4[k]);
        }

        t=t+h; 
  
    }
    pos.close();
    poin.close();
    poin1.close();
    poin2.close();
    lya.close();
    //-------------------------------------------------------------------------------------
    return 0;
}

double f0(double O1, double O2, double p_o1, double p_o2){
    double f0;
    f0=(p_o1-p_o2*cos(O1-O2))/(1+pow(sin(O1-O2),2));
    return f0;
}

double f1(double O1, double O2, double p_o1, double p_o2){
    double f1;
    f1=(-p_o1*cos(O1-O2)+2*p_o2)/(1+pow(sin(O1-O2),2));
    return f1;
}

double f2(double O1, double O2, double p_o1, double p_o2){
    double f2, h1, h2;
    h1=(p_o1*p_o2*sin(O1-O2))/(1+pow(sin(O1-O2),2));
    h2=(pow(p_o1,2)+2*pow(p_o2,2)-2*p_o1*p_o2*cos(O1-O2))/(2*pow(1+pow(sin(O1-O2),2),2));
    f2=-2*g*sin(O1)-h1+h2*sin(2*(O1-O2));
    return f2;
}

double f3(double O1, double O2, double p_o1, double p_o2){
    double f3, h1, h2;
    h1=(p_o1*p_o2*sin(O1-O2))/(1+pow(sin(O1-O2),2));
    h2=(pow(p_o1,2)+2*pow(p_o2,2)-2*p_o1*p_o2*cos(O1-O2))/(2*pow(1+pow(sin(O1-O2),2),2));
    f3=-g*sin(O2)+h1-h2*sin(2*(O1-O2));
    return f3;
}