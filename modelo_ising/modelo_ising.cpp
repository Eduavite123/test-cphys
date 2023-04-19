#include <iostream>
#include <cmath>
#include <stdlib.h>  //srand(), radn() (provisional hasta ver cómo funciona gsl)
#include <time.h>
#include <fstream>
using namespace std;

#define MAX 100

void con_ini(int choose, int s[][MAX], int N);
int delta_E(int s[][MAX], int n, int m, int N);
double min(int dE, double temp);

int main(){
    //----------------------DECLARACIÓN DE VARIABLES-----------------------------------
    int N,n,m; //dimensión de la red y coordenadas
    int s[MAX][MAX]; //espin de cada nodo de la red (1 o -1)
    double temp; //temperatura (entre 0 y 5)
    double p, eps;
    int dE; //dE usado en la evalución de p
    int choose, pot_rand; //para elegir espines alineados o desalineados aleatoriamente
    int iter; //nº iteraciones
    //---------------------------------------------------------------------------------

    srand(time(NULL));//semilla para números random

    //----------------------VALORES INICIALES------------------------------------------
    //(0:todos alineados, 1:desalineados aleatoriamente)
    temp=4; //(entre 0 y 5) 
    N=10; //tamaño de la red
    choose=0; //0 o 1. Si temp alta conviene choose=0, si temp baja conviene choose=1
    con_ini(choose,s,N);      
    iter=500; 
    //---------------------------------------------------------------------------------

    //----------------------------MODELO DE ISING--------------------------------------
    ofstream spin;
    spin.open("ising_data.dat");
    for (int k = 0; k < iter; k++)
    {
        //srand(time(NULL));//cambiamos la semilla cada iteración para mejorar la aleatoriedad
        for (int l = 0; l < N*N; l++)
        {
            //elegir nodo al azar
            n=rand() % N;
            m=rand() % N;

            //Calcular dE
            dE=delta_E(s,n,m,N);

            //Calcular p
            p=min(dE,temp);

            //Generación de epsilon
            eps=(double)rand()/RAND_MAX; //entre 0 y 1

            //Cambiar signo del espín
            if (eps<p)
            {
                s[n][m]=-s[n][m];
            }
            
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N-1; j++)
            {
                spin << s[i][j] << ", ";
            }
            spin << s[n][N-1] << endl;   
        }
        spin << endl;
    }
    spin.close();
    //---------------------------------------------------------------------------------

    return 0;
}

void con_ini(int choose, int s[][MAX], int N){
    int pot_rand;
    if (choose==0)
    {
        pot_rand=rand(); //nº aleatorio que será exponente de -1
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                s[i][j]=pow(-1,pot_rand); //todos los espines valen lo mismo
            }
            
        }
    }
    else if (choose==1)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                s[i][j]=pow(-1,rand()); //los espines van arriba o abajo aleatoriamente
            }
        }
    }
}

int delta_E(int s[][MAX], int n, int m, int N){
    int dE;
    if (n!=N-1 && m!=N-1 && n!=0 && m!=0) //si no está en ningún borde
    {
        dE=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
    }
    //CONDICIONES DE CONTORNO PERIÓDICAS
    else if (n==0 && m!=0 && m!=N-1) //borde superior sin esquinas
    {
        dE=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][m-1]);
    }
    else if (n==N-1 && m!=0 && m!=N-1) //borde inferior sin esquinas
    {
        dE=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
    }
    else if (m==0 && n!=0 && n!=N-1) //borde izquierdo sin esquinas
    {
        dE=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
    }
    else if (m==N-1 && n!=0 && n!=N-1) //borde derecho sin esquinas
    {
        dE=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
    }
    else if (n==0 && m==0) //esquina superior izquierda
    {
        dE=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][m+1]+s[n][N-1]);
    }
    else if (n==0 && m==N-1) //esquina superior derecha
    {
        dE=2*s[n][m]*(s[n+1][m]+s[N-1][m]+s[n][0]+s[n][m-1]);
    }
    else if (n==N-1 && m==0) //esquina inferior izquierda
    {
        dE=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][m+1]+s[n][N-1]);
    }
    else if (n==N-1 && m==N-1) //esquina inferior derecha
    {
        dE=2*s[n][m]*(s[0][m]+s[n-1][m]+s[n][0]+s[n][m-1]);
    }
}

double min(int dE, double temp){
    double e;
    e=exp(-dE*1.0/temp);
    if (1<=e) return 1;
    else return e;
}