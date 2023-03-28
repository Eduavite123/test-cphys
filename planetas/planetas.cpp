#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

#define Ms 1.99e30 //kg
#define c 1.5e8 //km
#define G 6.67e-11 

void read_file(double vec[9], string nombre_archivo);
double modulo(double rx, double ry);
void acel(double a[9], double r[9], double raux[9], double mass[9]);
void pos(double r[9], double v[9], double a[9], int hop);
void aux(double aux[9], double v[9], double a[9], int hop);
void vel(double v[9], double aux[9], double a[9], int hop);

//programa que representa el sistema solar
int main(){
    //----------------------DECLARACIÓN DE VARIABLES-----------------------------------
    //planetas ordenados del más cercano al Sol al más alejado
    //0: Sol 1:Mercurio 2:Venus 3:Tierra 4:Marte 5:Júpiter 6:Saturno 7:Urano 8:Neptuno
    double m[9]={0}; //masa 
    double x[9]={0}, y[9]={0}; //vector direccion por componentes
    double vx[9]={0}, vy[9]={0}; //velocidad por componentes
    double ax[9]={0}, ay[9]={0}; //aceleración por componentes
    double wx[9]={0}, wy[9]={0}; //auxiliar para calcular la velocidad
    int k=0; //contadores
    float t=0, tmax, h; //tiempo y salto h
    int inter; //numero de iteraciones
    //------------------------------------------------------------------------------------

    //--------------CÁLCULO DEL NÚMERO DE ITERACIONES-----------------
    /* cout << "¿Cuántos años durará la animación?" << endl;
    cin >> tmax;
    tmax=tmax*365*24*60*60; //tiempo en segundos
    cout << "¿Cuántos segundos durarán los pasos temporales?" << endl;
    cin >> h; */
    h=0.01;
    tmax=50;
    inter=tmax/h; 
    //----------------------------------------------------------------

    //-----------------CONDICIONES INICIALES------------------------
    //lectura de ficheros con condiciones iniciales de r, m y v
    /* read_file(m, "m_ini.txt");
    read_file(x, "r_ini.txt");
    read_file(vy, "v_ini.txt"); */

    m[0]=1.99e30;
    m[1]=0.330e24;
    m[2]=4.87e24;
    m[3]=5.97e24;
    m[4]=0.642e24;
    m[5]=1898.0e24;
    m[6]=568.0e24;
    m[7]=86.8e24;
    m[8]=102e24;

    x[0]=0;
    x[1]=57.9e6;
    x[2]=108.2e6;
    x[3]=149.6e6;
    x[4]=228.0e6;
    x[5]=778.5e6;
    x[6]=1432.0e6;
    x[7]=2867.0e6;
    x[8]=4515.0e6;

    vy[0]=0;
    vy[1]=47.4;
    vy[2]=35.0;
    vy[3]=29.8;
    vy[4]=24.1;
    vy[5]=13.1;
    vy[6]=9.7;
    vy[7]=6.8;
    vy[8]=5.4;

    //Reescalado
    for (int i = 0; i < 9; i++)
    {
        m[i]=m[i]/Ms;
        x[i]=x[i]/c; //empezamos desde el eje x
        vy[i]=vy[i]/sqrt(G*Ms/c); //tomamos vx=0 a t=0
    }    

    //Calculamos las aceleraciones iniciales
    acel(ax, x, y, m); 
    //---------------------------------------------------------------

    //--------------------------ALGORITMO DE VERLET-------------------------
    ofstream fich;
    fich.open("planets_data.dat");
    for (k = 0; k<inter; k++)
    {
        //escribir los valores de x e y en un fichero
        for (int i = 0; i < 9; i++)
        {
            fich << x[i] << ", " << y[i] << endl;
        }
        fich << endl;

        t=t+h;

        pos(x,vx,ax,h);
        pos(y,vy,ay,h);

        aux(wx,vx,ax,h);
        aux(wy,vy,ay,h);
        
        acel(ax,x,y,m);
        acel(ay,y,x,m);
        
        vel(vx,wx,ax,h);
        vel(vy,wy,ay,h);

        cout << ax[1] << endl;

    }
    fich.close();
    //----------------------------------------------------------------------

    return 0;
}

//lee ficheros con condiciones iniciales
void read_file(double vec[9], string nombre_archivo){
    ifstream fich;
    fich.open(nombre_archivo);
    for (int i = 0; i<9; i++)
    {
        fich >> vec[i];
    }    
    fich.close();
    return;
}

//calcula el módulo de un vector
double modulo(double rx, double ry)
{
    double mod;
    mod=sqrt(pow(rx,2)+pow(ry,2));
    return mod;
}

/*calcula una componente de la aceleracion a partir de la suma de fuerzas que se 
ejercen sobre cada planeta i*/
void acel(double a[9], double r[9], double raux[9], double mass[9]){
    int i, j=0;
    for (int i = 0; i < 9; i++)
    {
        if (i!=j)
        {
        a[i]=a[i]-mass[j]*((r[i]-r[j]))/(pow(modulo(r[i]-r[j],raux[i]-raux[j]),3));
        }
    }
    return;
}

//evalúa una componente de la dirección r según el algoritmo de Verlet
void pos(double r[9], double v[9], double a[9], int hop){
    int i;
    for (i = 0; i < 9; i++)
    {
        r[i]=r[i]+hop*v[i]+(hop*hop/2.0)*a[i];
    }
    return;
}

//evalúa vector auxiliar que se usará para evaluar la velocidad
void aux(double aux[9], double v[9], double a[9], int hop){
    int i;
    for (i = 0; i < 9; i++)
    {
        aux[i]=v[i]+hop/2.0*a[i];
    }
    return;
}

////evalúa una componente de la velocidad v según el algoritmo de Verlet
void vel(double v[9], double aux[9], double a[9], int hop){
    int i;
    for (i = 0; i < 9; i++)
    {
        v[i]=aux[i]+hop/2.0*a[i];   
    }
    return;
}