#include <iostream>
#include<stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>

/* run this program using the console pauser or add your own getch, system("pause") or input loop */
using namespace std;

//CONSTANTES GLOBALES
const double me(9.10939e-31);
const double kB(1.38066e-23);	// J/°K  *
const double pi(3.141592653589793);

//VARIABLES GLOBALES
//secciones eficaces
double N2[47][2]= {{0, 4.077E-20},{0.1, 4.88e-20},{0.12, 5.13e-20},{0.15, 5.56e-20},{0.17, 5.85e-20},{0.2, 6.25e-20},{0.25, 6.84e-20},{0.3, 7.32e-20},{0.35, 7.72e-20},
{0.4, 8.06e-20},{0.45, 8.33e-20},{0.5, 8.61e-20},{0.6, 8.96e-20},{0.7, 9.25e-20},{0.8, 9.48e-20},{0.9, 9.66e-20},{1.0, 9.85e-20},{1.2, 10.2e-20},{1.5, 11.2e-20},
{1.7, 13.3e-20},{2.0, 25.7e-20},{2.5, 28.5e-20},{3.0, 21.0e-20},{3.5, 14.6e-20},{4.0, 13.2e-20},{4.5, 12.3e-20},{5.0, 11.8e-20},{6.0, 11.4e-20},{7.0, 11.4e-20},
{8.0, 11.5e-20},{9.0, 11.7e-20},{10, 12.0e-20},{12, 12.4e-20},{15, 13.2e-20},{17, 13.5e-20},{20, 13.7e-20},{25, 13.5e-20},{30, 13.0e-20},{35, 12.4e-20},{40, 12.0e-20},
{45, 11.6e-20},{50, 11.3e-20},{60, 10.7e-20},{70, 10.2e-20},{80, 9.72e-20},{90, 9.30e-20},{100, 8.94e-20}};
double N2e[26][2]={{0, 7.874E-20},{0.55, 8.39e-20},{0.7, 9.03e-20},{0.9, 9.62e-20},{1, 9.83e-20},{1.5, 10.53e-20},{2, 17.93e-20},{2.2, 19.5e-20},{2.35, 20.5e-20},
{2.5, 21e-20},{2.7, 17.5e-20},{3, 15e-20},{4, 11.6e-20},{5, 10.75e-20},{6, 10.6e-20},{8, 10.6e-20},{10, 11.4e-20},{15, 11.8e-20},{20, 11.15e-20},{25, 10.25e-20},
{30, 9.65e-20},{40, 8.85e-20},{50, 8.2e-20},{60, 7.4e-20},{80, 6.25e-20},{100, 5.6e-20}};
double O2[32][2]={{0, 0.15},{0.15, 4.96e-20},{0.23, 4.27e-20},{0.26, 4.94e-20},{0.38, 5.07e-20},{0.5, 5.115e-20},{0.6, 5.21e-20},{0.75, 5.77e-20},{0.84, 5.59e-20},
{0.92, 5.85e-20},{0.95, 5.78e-20},{1, 5.8e-20},{2, 6.44667e-20},{3, 6.64667e-20},{4, 6.85e-20},{5, 7.604e-20},{6, 7.9325e-20},{7, 8.43e-20},{8, 8.95e-20},{9, 9.81667e-20},
{10, 1.04875e-19},{12, 1.076e-19},{15, 1.10167e-19},{20, 1.111e-19},{30, 1.1255e-19},{40, 1.091e-19},{50, 1.066e-19},{60, 1.006e-19},{70, 9.7475e-20},{80, 9.245e-20},
{90, 9.025e-20},{100, 8.67e-20}};
double O2e[22][2]={{0, 5.73E-20},{1, 5.95E-20},{2, 6.46667E-20},{3, 7.16667E-20},{4, 6.86E-20},{5, 7.6625E-20},{6, 7.9E-20},{7, 8.275E-20},{8, 8.2E-20},
{9, 8.35E-20},{10, 9.3575E-20},{12, 9.1E-20},{15, 8.9225E-20},{20, 8.605E-20},{30, 8.26667E-20},{40, 7.3E-20},{50, 6.55E-20},{60, 6.1E-20},{70, 5.5E-20},
{80, 5.4E-20},{90, 5E-20},{100, 4.3E-20}};

//DECLARACIÓN DE FUNCIONES
double * vuelo(double aceleracion[], double vel_i[], double densidad);
double interpolacion(double valor, double arreglo[][2], int size);
double colision(double sigma[], double sigma_t);
double * ensayo();

//Main
int main() {
	srand(time(NULL));	/* Toma un número a partir de la fecha y hora actual.*/
	ofstream impactos;
	impactos.open("impactos.csv", ios::out);
	cout <<"Ingrese el número de ensayos: ";
	int tiros;
	cin >> tiros;
	for(int t=0; t< tiros; t++){
		double *salida;	//arreglo que recibe los valores
		salida = ensayo();
		cout << salida[0] << "	" << salida[1] << "	" << salida[2] << "\n";
		impactos << salida[0] << "	" << salida[1] << "	" << salida[2] << "\n";
		impactos.close();
	}
return 0;
}

double * ensayo(){
	double presion = 760;	//preion en mmHg
	presion = presion * (10132.5024/76);	//presion en atm
	double temp = 293.15;	//temperatura en °K
	double densidad;
	densidad = presion * 6.022140857e23 / (8.314472 * temp);	//densidad molécualr - N/m^2*g/mol * molec/mol
	double campo = 1e6;	//	V/m campo electrico
	double Eo; //unidad atómica de energía en Joules
	Eo = me* pow(1.60218e-19,4) / (pow(4*pi*8.85419e-12,2) * pow(1.05457e-34,2));
	//iniciando archivo de salida de la posición
	ofstream posicion;
	posicion.open("posicion.csv", ios::out);
	// Posición inicial
	double posic[3] = {0,0,0};
	posicion << posic[0] << "	" << posic[1] <<"	"<<posic[2]<<"\n";
	// Velocidad inicial
	double vel_i[3] = {0,0,0};
	// inicialización del archivo para la energía cinética
	ofstream energia;
	energia.open("Energia.csv", ios::app);
	double ener = 0;	//inicialización de la energía de impacto del electrón
	energia << ener << "\n";
	//aceleración
	double aceleracion[3] = {0, 0, 1.60218e-19 * campo / me};
	//archivo de salida de los tiempos de vuelo, lapso entre colisiones
	ofstream lapso;
	lapso.open("tiempo de vuelo.csv", ios::out);
	//archivo de salida para los ángulos de dispersión
	ofstream dispersion;
	dispersion.open("dispersion.txt", ios::out);
	// inicialización del tiempos de impacto
	double Tv = 0;
	lapso << Tv << "\n";
	float dist = 0.01; //distancia entre las placas
//	cout<< "Ingrese la distancia entre las placas: ";
//	cin >> dist;
	int k = 0; //contador de lazo de trayectoria
	cout << "comienza el vuelo \n\n";
	while (posic[2] < dist){
		//Función para variables de vuelo libre
		double *V;
		V = vuelo(aceleracion, vel_i, densidad);
		Tv = V[0];	//Tiempo de vuelo
		ener = V[1];	//Energía de impacto
		double sigma[4] = {V[2], V[3], V[4], V[5]};
		double sigma_t = V[6];
//V = {Tv, Energia_k[ind], sigma[ind][0], sigma[ind][1], sigma[ind][2], sigma[ind][3], sigma_t[ind]}
		
		//Función para tipo de colisón
		double elast;
		elast = colision(sigma, sigma_t);
			
		//IMPACTO DEL ELECTRON Y UNA MOLÉCULA
		// actualizacion de la velocidad del electron al impacto
		double vel[3];
		vel[0] = vel_i[0] + aceleracion[0] * Tv;
		vel[1] = vel_i[1] + aceleracion[1] * Tv;
		vel[2] = vel_i[2] + aceleracion[2] * Tv;
		// actualización de la posición
		posic[0] = posic[0] + vel_i[0] *Tv + aceleracion[0] / 2 *pow(Tv,2);
		posic[1] = posic[1] + vel_i[1] *Tv + aceleracion[1] / 2 *pow(Tv,2);
		posic[2] = posic[2] + vel_i[2] *Tv + aceleracion[2] / 2 *pow(Tv,2);
		
		//Solo cuando el salto excede los límites físicos de la simulación
		if (posic[2] > dist){
			double t;
			t = Tv - (vel[2] - sqrt(pow(vel[2],2) + 2*aceleracion[2]*(dist - posic[2]))) / aceleracion[2];
			if (t<0){
			t = Tv - (vel[2] + sqrt(pow(vel[2],2) + 2*aceleracion[2]*(dist - posic[2]))) / aceleracion[2];
			}
			vel[0] = vel[0] - aceleracion[0] * (Tv -t);
			vel[1] = vel[1] - aceleracion[1] * (Tv -t);
			vel[2] = vel[2] - aceleracion[2] * (Tv -t);
			ener = vel[0]*vel[0]+ vel[1]*vel[1]+ vel[2]*vel[2];	//Energía de impacto calculada con vel
			ener = ener * me /(2 * 1.60218e-19);
			posic[0]= posic[0] - vel[0] * (Tv -t) - aceleracion[0] * pow(Tv -t,2) /2;
			posic[1]= posic[1] - vel[1] * (Tv -t) - aceleracion[1] * pow(Tv -t,2) /2;
			posic[2]= posic[2] - vel[2] * (Tv -t) - aceleracion[2] * pow(Tv -t,2) /2;
			Tv = t;
		}
		//ingreso de posición y tiempo en archivos externos.
		posicion << posic[0] << "	" << posic[1] << "	" << posic[2] << "\n";
		lapso << Tv << "\n";
		//guardar la energía en archivo
		energia << ener << "\n";
		
		//CALCULO DE LAS MAGNITUDES LUEGO DE LA COLISIÓN
		double eps = ener * me /(2 * Eo);
		double Ephi; //parametro de corrección para gases diatómicas
		Ephi = (0.065*eps+0.26*sqrt(eps)) / (1+0.05*eps+0.2*sqrt(eps)) - 12*sqrt(eps) / (1 + 40*sqrt(eps));
		//dispersión del electrón al chocar
		double R2 = rand(); //Probabilidad aleatoria de dispersión respecto a la dirección
		R2 = R2 /RAND_MAX;
		double theta; //ángulo axial
		theta = acos(1 - (2*R2*(1-Ephi)) / (1+Ephi*(1-2*R2)));
		dispersion << theta << "\n";
		double fi = rand(); //ángulo en el plano perpendicular entre 0 y 2pi
		fi = fi / RAND_MAX * 2 * pi;
		//el módulo de la velocidad de salida debe multiplicarse por un coeficiente de elasticidad del choque.
		double veloc;
		double vel_sq; //cuadrado de la velocidad
		vel_sq = vel[0]*vel[0]+ vel[1]*vel[1]+ vel[2]*vel[2];
		veloc = sqrt (vel_sq) * elast;	//módulo de la velocidad de salida
		//Velocidad de salida
		if ((vel[0]==0)&&(vel[1]==0)){
			vel_i[0] = veloc * sin(theta) * cos(fi);
			vel_i[1] = veloc * sin(theta) * sin(fi);
			vel_i[2] = veloc * cos(theta);
		}else{
			double cos_psi, sen_psi;
			double cos_alf, sen_alf;
			cos_psi = vel[2] / sqrt(vel_sq);
			sen_psi = sqrt(pow(vel[0],2) +pow(vel[1],2)) / sqrt(vel_sq);
			cos_alf = vel[0] / sqrt(pow(vel[0],2) + pow(vel[1],2));
			sen_alf = vel[1] / sqrt(pow(vel[0],2) + pow(vel[1],2));
		//magnitudes V_z V_r Theta y fi en el sistema rotado (CM) de coordenadas.
			double V_z = veloc * cos(theta);
			double V_r = veloc * sin(theta);
		//Velocidad de salida del choque
			vel_i[0] = (V_r * cos(fi) * cos_psi + V_z * sen_psi) * cos_alf - V_r * sin(fi) * sen_alf;
			vel_i[1] = (V_r * cos(fi) * cos_psi + V_z * sen_psi) * sen_alf + V_r * sin(fi) * cos_alf;
			vel_i[2] =  V_z * cos_psi  -  V_r * cos(fi) * sen_psi;
		}
		k = k+1;
	//	cout << "saltos: "<< k << endl;
		if(k >= 100000){	
		cout<<"\nDemasiado para este pobre ser digital\n" <<"saltos: "<< k << endl;
		break;
		}
	}
	cout<<"el número de choques fue: "<< k <<"\n"
	<< "La posición final: " << posic[0] << " " << posic[1] << " " << posic[2] << "\n";
	static double salida[3] = {posic[0], posic[1], posic[2]};
	
	// Cierre de los archivos externos
	posicion.close(); //cierra archivo de posición
	energia.close();
	lapso.close();
	dispersion.close();
	return salida;
}

//FUNCIONES
//Función de vuelo libre
double * vuelo(double aceleracion[], double vel_i[], double densidad){
	double tiempo[1000]; //intervalo de colisiones
	tiempo[0] = sqrt(2*1e2 * 1.60218e-19 / me) / aceleracion[2];// t = (v-vo)/a zona temporal no relativista
	double dt = tiempo[0] / 1e3;       //tamaño del paso
	double ch =dt;
	double tramo[1000];
	tramo[0] = 0;
	double vel[3] = {vel_i[0], vel_i[1], vel_i[2]};	// vel inicial del electron
	double Energia_K[2000]; // Energía cinética en el tramo
	Energia_K[0] = (pow(vel[0],2) +pow(vel[1],2) +pow(vel[2],2)) * me/(2*1.60218e-19);	//energía en electrón-voltios
	tiempo[0] = 0;
	//vuelo
	int i = 0;
	while (Energia_K[i] <= 100){	//energía en electrón-voltios
		tramo[i+1] = tramo[i] +  sqrt(pow(vel[0]*dt,2) +pow(vel[1]*dt,2) +pow(vel[2]*dt + aceleracion[2] /2 *pow(dt,2),2));
		//actualización de la velocidad sobre la trayectoria
		vel[0] = vel[0] + aceleracion[0] *dt;
		vel[1] = vel[1] + aceleracion[1] *dt;
		vel[2] = vel[2] + aceleracion[2] *dt;
		//La energía en cada punto de la trayectoria
		Energia_K[i+1] = (pow(vel[0],2) +pow(vel[1],2) +pow(vel[2],2)) * me/(2*1.60218e-19);	//energía en electrón-voltios
		tiempo[i+1] = tiempo[i] + dt;
		if(dt != ch){
			cout<<"Error en definicion de intervalo de E_k: " << i
			<< " E_k " << Energia_K[i]<< " Tr " << tramo[i]<<endl
			<< "dt: " << dt << " verificador: " << ch << endl;
			break;
		}
		i = i + 1;
	}
	int limit = i; //longitud de los arreglos de energía y sección eficaz en la región de interés

	//CALCULO DE LA SECCIÓN EFICAZ SOBRE LA TRAYECTORIA, ENTRE COLISIONES
	double sigma[4][limit];	//seccion eficaz en el tramo
	double sigma_t[limit];	//seccion eficaz en el tramo
	int j = 0;
	while (j<limit){
		//Secciones eficaces en función del tiempo
		sigma[0][j] = interpolacion(Energia_K[j], N2,47);
		sigma[1][j] = interpolacion(Energia_K[j], N2e,26);
		sigma[2][j] = interpolacion(Energia_K[j], O2,32);
		sigma[3][j] = interpolacion(Energia_K[j], O2e,22);
		
		if (sigma[1][j] > sigma [0][j]){
			sigma[1][j] = sigma [0][j];
		}
		if (sigma[3][j] > sigma [2][j]){
			sigma[3][j] > sigma [2][j];
		}
		j = j + 1;
	}
	
	//La secion eficaz total	 
	for(int I; I<limit; I++){
		sigma_t[I] = 0.78*sigma[0][I] + 0.22*sigma[2][I];
	}
	
	//CALCULO DEL TIEMPO DE VUELO CON UNA FUNCIÓN ALEATORIA
	double Tv = 0; //tiempo de vuelo
	int indice,ind[100], conteo=0;
	while (Tv <= dt ){
		double arreglo[limit];
		double aleat = rand();	//numero aleatorio
		aleat = aleat/ RAND_MAX; //decimal aleatorio en ]0; 1]
//		cout << "\nAleatorio " << aleat << " log:"<< log(aleat) << endl;
	if(aleat != 0){
		int raices = 0;
		for (int I=0; I<limit; I++){
			arreglo[I] = densidad * sigma_t[I] * tramo[I] + log(aleat);
//if(conteo==200){			trayecto << "I:" << arreglo[I] << "\n";}
		}
		//busco la raiz
		for (int I=1; I<limit;I++){
			if ((arreglo[I]==0)||(arreglo[I-1]<0 && arreglo[I]>=0)){
				ind[raices] = I;
				raices = raices + 1;
//				cout << " raiz"<< I <<": "<<tiempo[I]<<"	";
			}
			if (arreglo[I-1]>0 && arreglo[I]<=0){
				ind[raices] = I;
				raices = raices + 1;
//				cout << " raiz"<< I <<": "<<tiempo[I]<<"	";
			}
		}
		conteo = conteo + 1;
		if(raices >1){
			int iap = rand()%raices;
			indice = ind[iap];
//			cout << "\nRaices "<<raices<<" indice: " << indice << " iap "<< iap + 1 << endl;
		}
		if(raices == 1){
			indice = ind[0];
//			cout << "\nRaices "<<raices<<" indice: " << indice << endl;
		}
		Tv = tiempo[indice];
	}
		if (conteo>200){
//			cout << "Error de busqueda de Tv [ln295]" <<endl;
			break;
		}
	}
//	cout<<"Tv: "<< Tv << " dt: " << dt << " conteo: "<< conteo << endl;
	if(Tv == 0){
		cout << "Error en def de tiempo de vuelo"<< endl;
	}
	static double intervalo[7];	//Arreglo que será exportado
	intervalo[0] = Tv;
	intervalo[1] = Energia_K[indice];
	intervalo[2] = sigma[0][indice];	intervalo[3] = sigma[1][indice];	// N2
	intervalo[4] = sigma[2][indice];	intervalo[5] = sigma[3][indice];	// O2
	intervalo[6] = sigma_t[indice];

	return intervalo;
}

//FUNCIÓN INTERPOLACIÓN
double interpolacion(double valor,double arreglo[][2], int size){
	double X1, X2;
	double Y1, Y2;
	for (int i=0;i < size; i++){
		if (valor >= arreglo[i][0]){
			X1 = arreglo[i][0];		Y1 = arreglo[i][1];
			X2 = arreglo[i+1][0];	Y2 = arreglo[i+1][1];
		}
	}
	double Y;
	Y = (Y2 - Y1)/(X2 - X1)*(valor-X1) + Y1; //interpolación de Taylor de 1er orden
	return Y;
}

//FUNCION TIPO DE COLISIÓN
double colision(double sigma[], double sigma_t){
	//criterios de probabilidad de colisión
	double molec = rand(); //identidad de la molécula
	molec = molec / RAND_MAX;
	double elast;
	// Tipo de molecula nitrogeno
	if (molec<= 0.78*sigma[0] / sigma_t){
		double tipo = rand();
		tipo= tipo / RAND_MAX;
		elast = 1;		//elastica
		if (tipo >= sigma[1] / sigma[0]){
			float aleat = rand();
			elast = (1 - pow(aleat/RAND_MAX,2));// Coeficiente de elasticidad
		}
	}else{ //tipo de molécula oxigeno
		double tipo = rand();
		tipo= tipo/RAND_MAX;
		elast = 1;		//elastica
		if (tipo >= sigma[3] / sigma[2]){
			float aleat = rand();
			elast = (1 - pow(aleat/RAND_MAX,2));// Coeficiente de elasticidad
		}
	}
	return elast;
}