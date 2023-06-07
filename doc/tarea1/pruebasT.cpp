#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

int main()
{
    // Semilla para generar números aleatorios
    srand(time(0));
    for (int i = 0; i < 10; i++)
    {
        // Generar los números aleatorios
        double tminus = (rand() / static_cast<double>(RAND_MAX)) * 100 - 50; // Número aleatorio entre -50 y 50
        double tplus = (rand() / static_cast<double>(RAND_MAX)) * 100 - 50;  // Número aleatorio entre -50 y 50
        double t = 0;
        //ambos positivos
        if (tminus > 0  && tplus > 0)
        {
            t = min(tminus, tplus);
        }
        // tminus positivo, tplus negativo
        else if (tminus > 0 && tplus < 0)
        {
            t = tminus;
        }
        // tminus negativo, tplus positivo
        else if (tminus < 0 && tplus > 0)
        {
            t = tplus;
        }
        else
        {
            t = 0;
        }

        cout << "tminus: " << tminus << endl;
        cout << "tplus: " << tplus << endl;

        // Imprimir el resultado
        cout << "El menor número positivo es: " << t << endl;
    }
    

    return 0;
}
