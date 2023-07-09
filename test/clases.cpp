#include <iostream>
#include <string>
using namespace std;

class Figuras
{
private:
    /* data */
    string color = "Verde";

public:
    Figuras(string _color)
    {
        color = _color;
    };
    double area();
    double perimetro();
    void setColor(string _color)
    {
        color = _color;
    }
    string getColor()
    {
        return color;
    }
    // ~Figuras();
};

class Cuadrado : public Figuras
{

public:
    double area(double lado)
    {
        return lado * lado;
    }
    double perimetro(double lado)
    {
        return lado * 4;
    }
};

class Triangulo : public Figuras
{
public:
    double area(double lado, double altura)
    {
        return lado * altura / 2;
    }
    double perimetro(double lado)
    {
        return lado * 3;
    }
};

class Circulo : public Figuras
{
private:
    double pi = 3.142592;
    double radio = 0.0;

public:
    Circulo(double _radio)
    {
        radio = _radio;
    }
    double area()
    {
        return pi * radio * radio;
    }
    double perimetro()
    {
        return pi * radio * 2;
    }
};

int main()
{
    Triangulo tri;
    Cuadrado quad;
    Circulo cir(10.0);

    cout << tri.getColor() + "\n";
    quad.setColor("Azul");
    cout << quad.getColor() + "\n";

    cout << "El area del circulo es : \n";
    cout << cir.area();

    return 0;
}