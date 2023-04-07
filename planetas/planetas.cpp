#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

#define Ms 1.99e30 //kg
#define c 1.5e8 //km
#define G 6.67e-20 

void read_file(double (&vec)[9], string nombre_archivo);
double modulo(double rx, double ry);
void acel(double (&a)[9], double r[9], double raux[9], double mass[9]);
void pos(double (&r)[9], double v[9], double a[9], float hop);
void aux(double (&aux)[9], double v[9], double a[9], float hop);
void vel(double (&v)[9], double aux[9], double a[9], float hop);
double ener(double mass[9], double vx[9], double vy[9], double x[9], double y[9]);

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
    double e=0; //energia total del sistema
    int k=0; //contadores
    float t=0, tmax, h; //tiempo y salto h
    long inter; //numero de iteraciones
    //------------------------------------------------------------------------------------

    //--------------CÁLCULO DEL NÚMERO DE ITERACIONES-----------------
    /* cout << "¿Cuántos años durará la animación?" << endl;
    cin >> tmax;
    tmax=tmax*365*24*60*60; //tiempo en segundos
    cout << "¿Cuántos segundos durarán los pasos temporales?" << endl;
    cin >> h; */
    h=0.1;
    tmax=50;
    inter=tmax/h; 
    //----------------------------------------------------------------

    //-----------------CONDICIONES INICIALES------------------------
    //lectura de ficheros con condiciones iniciales de r, m y v
    read_file(m, "m_ini.txt");
    read_file(x, "r_ini.txt");
    read_file(vy, "v_ini.txt"); 

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
    ofstream posi;
    ofstream energia;
    posi.open("planets_data.dat");
    energia.open("energia.dat");
    for (k = 0; k<inter; k++)
    {
        //escribir los valores de x e y en un fichero
        for (int i = 0; i < 9; i++)
        {
            posi << x[i] << ", " << y[i] << endl;
        }
        posi << endl;

        t=t+h;

        pos(x,vx,ax,h);
        pos(y,vy,ay,h);

        aux(wx,vx,ax,h);
        aux(wy,vy,ay,h);

        for (int i = 0; i < 9; i++)
        {
            ax[i]=0;
            ay[i]=0;
        }
        
        acel(ax,x,y,m);
        acel(ay,y,x,m);
        
        vel(vx,wx,ax,h);
        vel(vy,wy,ay,h);

        e=ener(m,vx,vy,x,y);
        //escribir los valores de la energía en un fichero
        energia << t << " " << e << endl;
        e=0;

    }
    posi.close();
    energia.close();
    //----------------------------------------------------------------------

    return 0;
}

//lee ficheros con condiciones iniciales
void read_file(double (&vec)[9], string nombre_archivo){
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
void acel(double (&a)[9], double r[9], double raux[9], double mass[9]){
    int i=0, j=0;
    for (i = 0; i < 9; i++)
    {
        for (j = 0; j < 9; j++)
        {
            if (i!=j)
            {
                a[i]=a[i]-mass[j]*((r[i]-r[j]))/(pow(modulo(r[i]-r[j],raux[i]-raux[j]),3));
            }
        }
    }
    return;
}

//evalúa una componente de la dirección r según el algoritmo de Verlet
void pos(double (&r)[9], double v[9], double a[9], float hop){
    int i;
    for (i = 0; i < 9; i++)
    {
        r[i]=r[i]+hop*v[i]+(hop*hop/2.0)*a[i];
    }
    return;
}

//evalúa vector auxiliar que se usará para evaluar la velocidad
void aux(double (&aux)[9], double v[9], double a[9], float hop){
    int i;
    for (i = 0; i < 9; i++)
    {
        aux[i]=v[i]+hop/2.0*a[i];
    }
    return;
}

////evalúa una componente de la velocidad v según el algoritmo de Verlet
void vel(double (&v)[9], double aux[9], double a[9], float hop){
    int i;
    for (i = 0; i < 9; i++)
    {
        v[i]=aux[i]+hop/2.0*a[i];   
    }
    return;
}

//Calcula la energía total del sistema en un instante de tiempo
double ener(double mass[9], double vx[9], double vy[9], double rx[9], double ry[9]){
    double e=0, pot=0;
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if (i!=j)
            {
                pot=pot-mass[j]*mass[i]/modulo(rx[i]-rx[j],ry[i]-ry[j]);
            }
        }
        e=e+0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i])+pot;
    }
    return e;
}