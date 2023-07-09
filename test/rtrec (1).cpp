// rt: un lanzador de rayos minimalista
 // g++ -O3 -fopenmp rt.cpp -o rt
#include <math.h>
#include <stdlib.h>
#include <stdio.h>  
#include <omp.h>
#include <utility>
#include <time.h>       // for clock_t, clock(), CLOCKS_PER_SEC
//#define NUM_THREADS 12
#include <cstdlib>
#include <cmath>
#include<algorithm>
using namespace std;
double pi=3.14159265358979323846; //Creamos el valor de pi para facilitarnos varais cosas.

inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double();
}



class Vector 
{
public:        
	double x, y, z; // coordenadas x,y,z 
  
	// Constructor del vector, parametros por default en cero
	Vector(double x_= 0, double y_= 0, double z_= 0){ x=x_; y=y_; z=z_; }
  
	// operador para suma y resta de vectores
	Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
	Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
	// operator multiplicacion vector y escalar 
	Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }

	//operador para multiplicacion vector vector
	Vector operator*(const Vector &b) const { return Vector(x * b.x , y * b.y , z * b.z); }

	// operator % para producto cruz
	Vector operator%(Vector&b){return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
	
	// producto punto con vector b
	double dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; }

	// producto elemento a elemento (Hadamard product)
	Vector mult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }
	
	// normalizar vector 
	Vector& normalize(){ return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }

	double magnitud() const {
        return sqrt(magnitud2());
        }

    double magnitud2() const {
        return x*x + y*y + z*z;
        }
};


typedef Vector Point;
typedef Vector Color;

inline  Vector cross(const Vector &u, const Vector &v) {
    return Vector(u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x);
}

void coordinateSystem(const Vector &n, Vector &s, Vector &t) { //Esta es la funcion para crear el sistema de coordenadas locales
	if (std::abs(n.x) > std::abs(n.y)) {
		float invLen = 1.0f / std::sqrt(n.x * n.x + n.z * n.z);
		t = Vector(n.z * invLen, 0.0f, -n.x * invLen);
	} else {
		float invLen = 1.0f / std::sqrt(n.y * n.y + n.z * n.z);
		t = Vector(0.0f, n.z * invLen, -n.y * invLen);
	}
	s = cross(t, n);
	}
inline Vector unitVector(const Vector &v){
	return  v * (1.0 / sqrt(v.x * v.x + v.y * v.y + v.z * v.z)); 

};
inline Vector LocalGlobal(const Vector &n,Vector &x){
	Vector s,t;
	coordinateSystem(n,s,t);
	return s*x.x+t*x.y+n*x.z;
	 
}	
//inline Vector sqrtvec(const Vector &v){ return Vector(sqrt(v.x),sqrt(v.y),sqrt(v.z)); }

inline Vector GlobalLocal(const Vector &n,Vector &x){
	Vector s,t;
	coordinateSystem(n,s,t);
	return Vector (s.dot(x),t.dot(x),n.dot(x));
	 
}

inline Vector random_esfera() { //Esta funcion crea direcciones aleatorias esfericas.
    auto r1 = random_double();
    auto r2 = random_double();

	double theta=acos(1.0-2.0*r1);
	double phi=2.0*pi*r2;

    auto x = cos(phi)*sin(theta);
    auto y = sin(phi)*sin(theta);
    auto z = 1.0 - 2.0*r1;

    return Vector(x, y, z);
}
inline void CarteEsfericas(const Vector &w,double &theta,double &phi, double &r){
	r=sqrt(w.x*w.x+w.y*w.y+w.z*w.z);
	theta=acos(w.z/r);
	phi=atan(w.y/w.x);

}

inline Vector random_hemisferio() {//Esta funcion crea direcciones aleatorias en un hemisferio.
    double r1 = random_double();
    double r2 = random_double();

    double theta=acos(r1);
    double phi=2.0*pi*r2;

    double x = sin(theta)*cos(phi);
    double y = sin(theta)*sin(phi);
    double z = r1;

    return Vector(x, y, z);
}

inline Vector random_coseno(double &theta) {//Esta funcion crea direcciones aleatorias con distribucion coseno hemisferico.
    double r1 = random_double();
    double r2 = random_double();

    double phi = 2.0*pi*r2;

    double z = sqrt(1.0-r1);
	theta = acos(z);

    double x = cos(phi)*sin(theta);
    double y = sin(phi)*sin(theta);

    return Vector(x, y, z);
}
inline void random_parametroscoseno(double &theta,double &phi) {
    double r1 = random_double();
    double r2 = random_double();
    phi = 2.0*pi*r2;
	theta = acos(sqrt(1.0-r1));
}

inline Vector esfCarte(double  &theta,double &phi) {
    double x = cos(phi)*sin(theta);
    double y = sin(phi)*sin(theta);
	double z = cos(z);

    return Vector(x, y, z);
}


class Ray 
{ 
public:
	Point o;
	Vector d; // origen y direcccion del rayo
	Ray(){}
	Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};

struct registro{
	Vector n;
	Point x;
	double t;
};

class material {
    public:
        virtual Color Emite(const Point& p) const {
            return Color(0,0,0);
        }

        virtual bool Rebota(
            const Ray &wi, Color &atenuacion,Ray &wo,double &pdf,registro &rec) const = 0;
};

class Luz : public material {
    public:

        Luz(Color c) : emit(c) {}

        virtual bool Rebota(const Ray &wi, Color &atenuacion,Ray &wo,double &pdf,registro &rec) const override {
            return false;
        }

        virtual Color Emite(const Point& p) const override {
            return emit;
        }

    public:
        Color emit;
};

class Abedo : public material {
    public:
        Abedo(const Color& a) : albedo(a) {}

        virtual bool Rebota(const Ray &wi, Color &atenuacion,Ray &wo,double &pdf,registro &rec) const override {
			double ri,thetai,phii;
			CarteEsfericas(wi.d,thetai,phii,ri);
			double ro,thetao,phio;
			wo.o=rec.x;
			wo.d=random_coseno(thetao);

			pdf=(1.0/pi)*cos(thetao);

			CarteEsfericas(wo.d,thetao,phio,ro);	

			double sigma=0.5*0.5;
			double A=1.0-sigma/(2.0*(sigma+0.33));
			double cosio=cos(phii-phio);
			double alpha = ( thetai > thetao) ? thetai : thetao;
			double beta = ( thetai < thetao) ? thetai : thetao;
			
			double B=0;
			if (cosio=!0){
				double B=(0.45*sigma)/(sigma+0.09)*cosio*sin(alpha)*tan(beta);	
			}
			
			wo.d=LocalGlobal(rec.n,wo.d);
			atenuacion = albedo*(1.0/pi)*(A+B)*(1.0/pdf);
			//atenuacion = albedo*(1.0/pi);
            return true;
        }
    public:
        Color albedo;
};

inline double Fresnel(const double &etap,const double &kapap,const double &coste){

	double sente = sin(acos(coste));
	double nksin = etap - kapap - sente * sente;		
	
	double ab = sqrt(nksin*nksin + etap*kapap*4.0);
	double a = sqrt((ab+nksin)/2.0);

	double rpernum=ab+coste*coste-2.0*a*coste;
	double rperdem=ab+coste*coste+2.0*a*coste;
	double rper=rpernum/rperdem;

	double rparnum=ab*coste*coste + sente * sente * sente * sente - 2.0*a*coste *sente * sente;
	double rpardem=ab*coste*coste + sente * sente * sente * sente + 2.0*a*coste *sente * sente;
	double rpar=rper*rparnum/rpardem;

	return (1.0/2.0)*(rper+rpar);
};
inline double Beckmann(double const &costh){
	//printf("%f\n",costh);
	if (costh>0){
	double alpha=0.5*0.5;
	double costh4=costh*costh*costh*costh;
	double tanh = tan(acos(costh));
	return exp(-(tanh*tanh)/alpha)/(pi*alpha*costh4);
	}else return 0.0;
};


inline double Gsmith(const Vector &wi,const Vector &wh,const Vector &wo,const Vector &n){
	double g1=0;
	double g2=0;
	double cosi=wi.dot(wh);
	double coso=wo.dot(wh);
	double alpha=0.5;

	double thetai=acos(cosi);
	double thetao=acos(coso);

	if ( ( cosi / wi.dot(n) ) > 0 ){
		double a1=1.0/(alpha*tan(thetai));
		if (a1>1.6)
			g1=(3.535*a1+2.181*a1*a1)/(1+2.276*a1+2.577*a1*a1);	
		else g1=1.0;
	} 

	if ( (coso / wo.dot(n)) >  0){
		double a2=1.0/(alpha*tan(thetao));
		if (a2>1.6)
			g2=(3.535*a2+2.181*a2*a2)/(1+2.276*a2+2.577*a2*a2);	
		else g2=1.0;
	} 

   return g1*g2;
};
Vector refleccion(const Vector& v, const Vector& n) {
    return  v-n*(v.dot(n)*2.0);
}

class MicroFasetC: public material {
    public:
        MicroFasetC(const Color& a,const Color& b) : eta(a),kappa(b) {}

        virtual bool Rebota(const Ray &wi, Color &atenuacion,Ray &wo,double &pdf,registro &rec) const override {
			wo.o = rec.x;

			Vector v =wi.d;
			v.normalize();
			wo.d = refleccion(v,rec.n);
			Vector wh = wo.d + wi.d;

			double coste = wh.dot(v);

			Color etap = eta*eta*(1.0/(1.00029*1.00029));
			Color kapap = kappa*kappa*(1.0/(1.00029*1.00029));

			//double cosh = rec.n.dot(wh);
			//printf("%f\n",coste);
			double R=Fresnel(etap.x,kapap.x,coste);
			double G=Fresnel(etap.y,kapap.y,coste);
			double B=Fresnel(etap.z,kapap.z,coste);
			//double div = rec.n.dot(v)*rec.n.dot(wo.d)*4.0;
			
			atenuacion = Color(R,G,B)*(1.0/(coste));// *div))*Beckmann(cosh)* Gsmith(v,wh,wo.d,rec.n)

			pdf=1.0;
            return (wo.d.dot(rec.n) > 0);//
        }
    public:
        Color eta;
		Color kappa;
};


class Conductor: public material {
    public:
        Conductor(const Color& a,const Color& b) : eta(a),kappa(b) {}

        virtual bool Rebota(const Ray &wi, Color &atenuacion,Ray &wo,double &pdf,registro &rec) const override {
			wo.o = rec.x;

			Vector v =wi.d;
			wo.d = refleccion(v,rec.n);

			double coste = wo.d.dot(rec.n);
			Color etap = eta*eta*(1.0/(1.00029*1.00029));//*(1.0/(1.00029*1.00029))
			Color kapap = kappa*kappa*(1.0/(1.00029*1.00029));//
			
			double R=Fresnel(etap.x,kapap.x,coste);
			double G=Fresnel(etap.y,kapap.y,coste);
			double B=Fresnel(etap.z,kapap.z,coste);
			
			atenuacion = Color(R,G,B)*(1.0/coste);
			pdf=1.0;
            return (wo.d.dot(rec.n) > 0);//
        }
    public:
        Color eta;
		Color kappa;
};

Vector transmision(Vector &v, const Vector &n,double &eta,double &sint) {
	//v=GlobalLocal(n,v);
	double cost=sqrt(1.0-sint);
	Vector wt=v*(-eta)+n*(v.dot(n)*eta-cost);
	//wt=LocalGlobal(n,wt);
    return wt;
}

Vector Transmision(const Vector& uv, const Vector& n, double etai_over_etat) {
    double cos_theta = fmin(-n.dot(uv), 1.0);
    Vector r_out_perp =  (uv + n*cos_theta)*etai_over_etat  ;
    Vector r_out_parallel =n*( -sqrt(fabs(1.0 - r_out_perp.magnitud2())) );
    return r_out_perp + r_out_parallel;
}
class Dielectrico: public material {
    public:
        Dielectrico(const double &a,const double &b) : eta(a),kappa(b) {}

        virtual bool Rebota(const Ray &wi, Color& atenuacion,Ray &wo,double &pdf,registro &rec) const override {
			wo.o = rec.x;
			

			Vector v = wi.d;
			
			//v.normalize();
			//Vector v = unitVector(wi.d);
			double cosi =-v.dot(rec.n);
			double indicerefrac = cosi > 0.0  ? eta : 1.0/eta;

			//double cost=sqrt(1.0-indicerefrac*indicerefrac*(1.0-cosi*cosi));
			double sint=sqrt(1.0-cosi*cosi);
			
			//bool notrasmite = indicerefrac * sint > 1.0;
			double f=F(cosi,indicerefrac);

			if (f>random_double()){
				wo.d=refleccion(v,rec.n);
				pdf=f;
			}else{ 
				//v =GlobalLocal(rec.n,v);
				//wo.d=transmision(v,rec.n,indicerefrac,sint);
				//wo.d=Transmision(v,rec.n,indicerefrac);
				double cost=sqrt(1-eta*eta*(1-cosi*cosi));
				if (cosi>0.0)
					wo.d=Vector(-indicerefrac*v.x,-indicerefrac*v.y, cost);
			   else wo.d=Vector(-indicerefrac*v.x,-indicerefrac*v.y,-cost);
				//wo.d=LocalGlobal(rec.n,wo.d);
				pdf=1.0-f;}
				//wo.d=refleccion(wi.d,rec.n);
				//pdf=f;
			atenuacion= Color(1.0,1.0,1.0);
			
            return true;//


        }
    public:
        double eta;
		double kappa;
	private:
			static double F(double &cosi,double &eta){
				double r0=(1.0 - eta)/(1.0+eta);
				r0=r0*r0;
				return r0+(1.0-r0)*pow((1.0-cosi),5); 
				//Para F utilize la aproximacion de Sclick que esta en el el libro
				//Raytracing in one weekend la formula esta igual e wikipedia.
			}	
};



class Sphere 
{
public:
	double r;	// radio de la esfera
	Point p;	// posicion 
	material *m;	// Material de la sefera ****Proyecto 2*******

	Sphere(double r_, Point p_,material* m_): r(r_), p(p_), m(m_) {}
  
	double intersect(const Ray &ray) const {
		Vector oc = ray.o-p;

		double a = ray.d.dot(ray.d);
		double b =  oc.dot(ray.d);
		double c = oc.dot(oc)-r*r;
		double discriminant= b*b - a*c;

		if (discriminant==0){
			if (-b>0.0)
				return -b;
			else return 0.0;

		}
		
		if (discriminant<0) {
			return 0.0;
		}
		else{
			double tpositivo = -b + sqrt(discriminant);
			double tnegativo = -b - sqrt(discriminant);
			double t;
			if (tpositivo > 0.0 && tnegativo > 0.0 )
			t = (tpositivo < tnegativo) ? tpositivo : tnegativo;
			else if(tpositivo > 0.0 && tnegativo < 0.0) 
			t = tpositivo;
			else if(tpositivo < 0.0 && tnegativo > 0.0)
			t = tnegativo;
			else if (tpositivo < 0.0 && tnegativo < 0.0)
			t=0.0;
			if (t < 0.001) return 0.0;
			else return t;
		}
	}
};

Luz  Esferaluminoza(Color(10.0, 10.0, 10.0));
Conductor EsAbaDer(Color(1.66058,0.88143,0.531467),Color(9.2282,6.27077,4.83803));//Aluminio
//MicroFasetC EsAbaDer(Color(1.66058,0.88143,0.531467),Color(9.2282,6.27077,4.83803));//Aluminio
//Conductor EsAbaDer(Color(0.143245,0.377423,1.43919),Color(3.98478,2.3847,1.60434));//Oro
//Conductor EsAbaDer(Color(0.208183,0.919438,1.110241),Color(3.92198,2.45627,2.14157));//Cobre
//Luz  Esferaluminoza(Color(1.0, 1.0, 1.0));
Dielectrico Esferacristal(1.5,2.4);
Abedo ParIzq(Color(.75, .25, .25));
Abedo ParDer(Color(.25, .25, .75));
Abedo ParedAt(Color(.25, .75, .25));
Abedo Suelo(Color(.25, .75, .75));
Abedo Techo(Color(.75, .75, .25));
Abedo EsAbaIz(Color(.2, .3, .4));
//Abedo EsAbaDer(Color(.4, .3, .2));


Sphere spheres[] = {
	//Escena: radio, posicion ,material    //Fara facilitarme las cosas quite el Color, y se lo agregue en el material.
        Sphere(1e5,  Point(-1e5 - 49, 0, 0),     &ParIzq), // pared izq
        Sphere(1e5,  Point(1e5 + 49, 0, 0),      &ParDer), // pared der
        Sphere(1e5,  Point(0, 0, -1e5 - 81.6),   &ParedAt), // pared detras
        Sphere(1e5,  Point(0, -1e5 - 40.8, 0),   &Suelo), // suelo
        Sphere(1e5,  Point(0, 1e5 + 40.8, 0),    &Techo), // techo
        Sphere(16.5, Point(-23, -24.3, -34.6),   &EsAbaIz), // esfera abajo-izq
		//Sphere(16.5, Point(-23, -24.3, -34.6),   &Esferaluminoza), // esfera abajo-izq
        Sphere(16.5, Point(23, -24.3, -3.6),     &EsAbaDer), // esfera abajo-der// Para observar las dos fuentes luminosas hay que comentar esta linea
		//Sphere(16.5, Point(23, -24.3, -3.6),     &Esferaluminoza), // esfera abajo-der // Para observar las dos fuentes luminosas hay que descomentar esta linea
        Sphere(10.5, Point(0, 24.3, 0),          &Esferaluminoza), // esfera arriba // esfera iluminada
		Sphere(7.5, Point(-23.0, -33.0, 30.0),          &Esferacristal)
};

// limita el valor de x a [0,1]
inline double clamp(const double x) { 
	if(x < 0.0)
		return 0.0;
	else if(x > 1.0)
		return 1.0;
	return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x) {
	return int( pow( clamp(x), 1.0/2.2 ) * 255 + .5); 
}


inline bool intersect(const Ray &r, double &t, int &id) {
	int NS=8;
	double aux[NS];


	for (int i=0;i<NS;i++){
        aux[i]=spheres[i].intersect(r);
		if ( aux[i]>0){
			t= aux[i];
			id=i;
		};
    };

	for (int i=0;i<NS;i++){
		if ( t>aux[i] && aux[i]>0.001 ){
			t= aux[i];
			id=i;
		};
    };
	if (t>0){
		//printf("%f  %d\t ",t,id);
		return true;
		};

	return false;
}

// Calcula el valor de color para el rayo dado
Color shade(const Ray &r,int prof) { //Agregamos la profundidad para hacer una funcion recursiva, esto nos permite lanzar un segundo rayo desde
	registro rec;
	double t; 						 
	double h;						 
	int id = 0;						 

	if (prof <= 0) // Si ya se ha llegado 
        return Color();

	if (!intersect(r, t, id)){
		return Color();}	// el rayo no intersecto objeto, return Vector() == negro

	const Sphere &obj = spheres[id];

	rec.x=r.d*t+r.o; 
	rec.n=(rec.x-obj.p).normalize();
	rec.t=t;

	Color attenuation;
    Color emite = obj.m->Emite(rec.x);

	Ray rebota;
	double pdf;

    if (!obj.m->Rebota(r, attenuation,rebota,pdf,rec)) 
        return emite;   

	
	double Coseno=rec.n.dot(rebota.d);
	
	return  emite + attenuation*shade(rebota, prof-1)*Coseno;
}


int main(int argc, char *argv[]) {
	double time_spent = 0.0;
	double muestras=32.0;
	int prof=10;
    clock_t begin = clock();
	//sleep(3);
 
	int w = 1024, h = 768; // image resolution
  
	// fija la posicion de la camara y la dirección en que mira
	Ray camera( Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize() );

	// parametros de la camara
	Vector cx = Vector( w * 0.5095 / h, 0., 0.); 
	Vector cy = (cx % camera.d).normalize() * 0.5095;
  
	// auxiliar para valor de pixel y matriz para almacenar la imagen
	Color *pixelColors = new Color[w * h];

	int NUM_THREADS=omp_get_max_threads();
	fprintf(stderr," \r Vamos a trabajar con %d hilos \n",NUM_THREADS);
	omp_set_num_threads(NUM_THREADS); //Lineas de Codigo para paralelizar
	#pragma omp parallel for //schedule(dynamic,1) //Con la paralelizacion se reduce en un 60 % aproxiamadamente el tiempo de ejecucion.

	for(int y = 0; y < h; y++) 
	{ 
		// recorre todos los pixeles de la imagen
		fprintf(stderr,"\r%5.2f%%",100.*y/(h-1));
		for(int x = 0; x < w; x++ ) {

			int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos

			Color pixelValue = Color(); // pixelValue en negro por ahora
			// para el pixel actual, computar la dirección que un rayo debe tener
			for (int i=0; i<muestras;i++)
			{
                
			Vector cameraRayDir = cx * ( double(x)/w - .5) + cy * ( double(y)/h - .5) + camera.d;
			// computar el color del pixel para el punto que intersectó el rayo desde la camara

			pixelValue = pixelValue + shade( Ray(camera.o, cameraRayDir.normalize()) ,prof)*(1.0/muestras);
			// limitar los tres valores de color del pixel a [0,1] 
			}
			//pixelValue = pixelValue;

			pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
		}
		
	
	
	}


	fprintf(stderr,"\n");
	clock_t end = clock();

	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
 
    printf("The elapsed time is %f seconds", time_spent);

	// PROYECTO 1
	// Investigar formato ppm
	FILE *f = fopen("ImplicitoRec.ppm", "w");
	// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int p = 0; p < w * h; p++) 
	{ // escribe todos los valores de los pixeles
    		fprintf(f,"%d %d %d \n", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), 
				toDisplayValue(pixelColors[p].z));
  	}
  	fclose(f);

  	delete[] pixelColors;

	return 0;
}
