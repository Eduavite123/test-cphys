#include <iostream>

using namespace std;

int main(){
    int num, i;

    cout << "Inserte un numero hasta el que contar: ";
    cin >> num;

    for(i=1;i<=num;i++){
        cout << i << "\n";
    }
    return 0;
}