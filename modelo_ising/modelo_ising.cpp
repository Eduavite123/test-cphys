#include <iostream>
#include <cmath>
#include <stdlib.h>  //srand(), radn() (provisional hasta ver cómo funciona gsl)
#include <time.h>
using namespace std;

#define MAX 100

int main(){
    //----------------------DECLARACIÓN DE VARIABLES-----------------------------------
    int N; //dimensión de la red
    int s[MAX][MAX]; //espin de cada nodo de la red (1 o -1)
    //REPASAR DECLARACION DE MATRICES LIBRO BERNIER
    double T; //temperatura (entre 0 y 5)
    double p; //probabilidad 
    double e; //número aleatorio para comparar con p
    int choose; //para elegir espines alineados o desalineados aleatoriamente
    //---------------------------------------------------------------------------------

    //----------------------VALORES INICIALES------------------------------------------
    //El usuario elegirá los valores iniciales:
    cout << "Elija los valores iniciales de los espines (0:todos alineados, 1:desalineados aleatoriamente):" << endl;
    cin >> choose;
    while(choose!=0 || choose!=1)
    {
        cout << "Por favor, elija una opción adecuada (0:todos alineados, 1:desalineados aleatoriamente):" << endl;
        cin >> choose;
    }
    //---------------------------------------------------------------------------------

    return 0;
}