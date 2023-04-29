#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

#define MAX 200
#define PI 3.14159265

using namespace std;

int main(){
    //------------------------DECLARACION DE VARIABLES-------------------------------------
    complex<double> phi[MAX]; //Solución ec. Schrodinger     ¡¡¡Dim máxima por determinar!!!
    int j,n; //posición y tiempo
    int t; //Tiempo máximo
    //Asumimos los parámetros de discretización espacial y temporal h=1 y s=1, respectivamente

    //Valores iniciales:
    int N; //tamaño pozo y nº ciclos (nc=1,...,N/4 para que haya 4 ciclos)
    float nc;
    float lambda; //proporcional a la altura del pozo (p.e. lambda=0.3)
    double x_0, sigma; //centro y anchura de la gaussiana de la f. onda inicial
    //Valores iniciales a generar:
    double S; //Espaciado temporal reescrito
    double k_0; //Proporcional a la altura del potencial 
    double V[MAX]; //Potencial según posición
    //Parámetros para la resolución del sistema de ecuaciones
    complex<double> alpha[MAX], A_neg, A_0[MAX], A_pos, gamma[MAX], beta[MAX], b[MAX], xi[MAX];
    complex<double> i (0.0,1.0); 
    //-------------------------------------------------------------------------------------

    //--------------------------CONDICIONES INICIALES--------------------------------------
    t=MAX;
    //Parámetros del pozo:
    N=MAX; //Tamaño del pozo    ¡¡¡A DETERMINAR!!!!
    nc=N/4.0; //# ciclos. Lo elegimos así para que haya 4 ciclos

    //Parámetros del potencial:
    lambda=0.3; //Altura del potencial
    k_0=2*PI*nc/(N*1.0);
    //Definimos el potencial
    for (j = 0; j < N; j++)
    {
        if (j<3.0/5.0*(double)N && j>2.0/5.0*(double)N)
        {
            V[j]=lambda*k_0*k_0;
        }
        else V[j]=0;        
    }

    //Condiciones de contorno:
    phi[0]=0.0;
    phi[N]=0.0;
    //Asumimos la función de onda inicial como una gaussiana
    /*Primero decidimos el centro y la anchura de la gaussiana 
    (conviene que la anchura sea mucho menor que el pozo y que esté centrada entre la primera 
    pared y el potencial)*/
    x_0=(double)N/6.0;
    sigma=(double)N/50.0;
    for (j = 1; j < N-1; j++)
    {
        phi[j]=exp(i*k_0*(double)j)*exp((-((double)j-x_0)*((double)j-x_0))/(2*sigma*sigma));
    }

    //Parámetros para la resolución del sistema de ecuaciones:
    S=1.0/(4.0*k_0*k_0); //parámetro de espaciado temporal
    A_neg=A_pos=1;
    for(j=0;j<N;j++)
    {
        A_0[j]=-2.0+2.0*i/S-V[j];
    }
    alpha[N-1]=0.0;
    for(j=0;j<=N-2;j++)
    {
        alpha[(N-2)-j]=-1.0/(A_0[(N-1)-j]+alpha[(N-1)-j]);
    }
    //-------------------------------------------------------------------------------------

    //--------------------------ALGORITMO EC. SCHRODINGER----------------------------------
    ofstream fonda;
    fonda.open("schrodinger_data.dat");
    for(n=0;n<t;n++)
    {
        //escritura de los datos
        for(j=0;j<N;j++)
        {
            fonda << j << ", " << real(phi[j]) /*<< ", " << imag(phi[j][n])*/ << endl;
        }
        fonda << endl;

        //Cálculo de b
        for(j=0;j<N;j++)
        {
            b[j]=4.0*i*phi[j]/S;
        }
        //Cálculo de beta
        beta[N-1]=0; 
        for(j=0;j<=N-2;j++)
        {
            beta[(N-2)-j]=(b[(N-1)-j]-beta[(N-1)-j])/(A_0[(N-1)-j]+alpha[(N-1)-j]);
        }

        //Cálculo de Xi
        xi[0]=xi[N-1]=0;
        for(j=0;j<N-1;j++)
        {
            xi[j+1]=alpha[j]*xi[j]+beta[j];
        }

        //Cálculo de phi en el siguiente instante de tiempo
        for(j=0;j<N;j++)
        {
            phi[j]=xi[j]-phi[j];
        }
    }
    fonda.close();

    return 0;
}