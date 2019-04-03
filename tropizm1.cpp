#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>


//#include "F:\headers\los_isaac.h"

#define LewyPunkt  0
#define PrawyPunkt  1
#define LiczbaPunktow 2

#define domyslny_a  1.0
#define domyslny_l0  1.0
#define domyslny_krokczasowy  80
#define domyslnyN 20
#define domyslnaliczbaargumentow 4

#define domyslny_wzrost_cien 1.004
#define domyslny_wzrost_slonce 1.002


//l-przyrosty
//v- okresla kierunek padania promieni



bool argumenty(int liczbaargumentow, char* args[], double& a, double&l0, int&liczbakrokowczasowych, int& N);
void inicjalizacja(int N, double***& x, int**& iscien, double**& l, double*& v, double*& kunttab, int t);
void zwolnijpamiec(int N, double***& x, int**& iscien, double**& l, double*& v, double*& kunttab);
void wpisaniepoczatkowych(int N, double a, double*** x, int** iscien, double l0, double** l, double* v);

void sprawdzamyczyjestcien(int N, int** iscien, double* v, double*** x);
void robimywzrost(int N, int** iscien, double** l, double wzrost_cien, double wzrost_slonce);
void przesuwanie(int N, double** l, double*** x, double a);
void wykonajkolejnykrok(int N, double***& x, int**& iscien, double**& l, double*& v, double*& kunttab, int t, double a, double wzrost_cien, double wzrost_slonce);

void zapiszdopliku(int N, int** iscien, double*** x);
void zapiszkunt(double* kunttab, int liczbakrokowczasowych);

int main(int argc, char* argv[])
{
    double ***x,*v,**l;
    double a,l0;
    int N,**iscien;
	double *kunttab;
	int liczbakrokowczasowych;


	bool ok = argumenty(argc, argv, a, l0, liczbakrokowczasowych, N);
	if(!ok)
		return 0;

	inicjalizacja(N, x, iscien, l, v,kunttab,liczbakrokowczasowych); // alokuje pamiêæ
	wpisaniepoczatkowych(N, a, x, iscien, l0, l, v);

	for(int t=0;t<liczbakrokowczasowych;t++)  //mega petla
	{
		wykonajkolejnykrok(N, x, iscien, l, v, kunttab, t, a, domyslny_wzrost_cien, domyslny_wzrost_slonce);
	}  //koniec mega petli

	for (int i = 0; i < N; i++)
	{
		printf("d = %f\n", ((l[LewyPunkt][i] - l[PrawyPunkt][i]) - x[LewyPunkt][i][1]));
	}

	zapiszkunt(kunttab, liczbakrokowczasowych);
	zapiszdopliku(N, iscien, x);

	zwolnijpamiec(N, x, iscien, l, v, kunttab); // zwalnia pamiêæ
    return 0;
}

void inicjalizacja(int N, double***& x, int**& iscien, double**& l, double*& v, double*& kunttab, int t)
{
	kunttab = (double *)malloc(t* sizeof(double)); // rozmiar: [liczba kroków czasowych]
	x = (double ***)malloc(LiczbaPunktow * sizeof(double)); // rozmiar: [2][N][2]
	iscien = (int **)malloc(LiczbaPunktow * sizeof(int)); // rozmiar: [2][N]
	l = (double **)malloc(LiczbaPunktow * sizeof(double)); // rozmiar: [2][N]

	for (int k = 0; k<LiczbaPunktow; k++)
	{
		x[k] = (double **)malloc(N*sizeof(double));
		iscien[k] = (int *)malloc(N*sizeof(int));
		l[k] = (double *)malloc(N*sizeof(double));

		for (int i = 0; i<N; i++)
		{
			x[k][i] = (double *)malloc(LiczbaPunktow * sizeof(double));
		}
	}

	v = (double *)malloc(2 * sizeof(double));

}

void zwolnijpamiec(int N, double***& x, int**& iscien, double**& l, double*& v, double*& kunttab)
{
	for (int k = 0; k < LiczbaPunktow; k++)
	{
		for (int i = 0; i<N; i++)
		{
			free(x[k][i]);
		}

		free(x[k]);
		free(iscien[k]);
		free(l);
	}

	free(kunttab);
	free(x);
	free(iscien);
	free(l);
	free(v);

	kunttab = NULL;
	x = NULL;
	iscien = NULL;
	l = NULL;
	v = NULL;
}

void wpisaniepoczatkowych(int N, double a, double*** x, int** iscien, double lpoczatkowe, double** l, double* v)
{
	for (int i = 0; i<N; i++)
	{
		x[0][i][0] = 0.0;
		x[0][i][1] = (double)i*lpoczatkowe;
		x[1][i][0] = a;
		x[1][i][1] = (double)i*lpoczatkowe;

		iscien[0][i] = iscien[1][i] = 0; // zak³adamy, ¿e nie ma cienia

		l[0][i] = l[1][i] = lpoczatkowe; //
	}

	v[0] = -0.4;
	v[1] = -1.5;

}

void sprawdzamyczyjestcien(int N, int** iscien, double* v, double*** x)
{

	for (int i = 1; i<N - 1; i++)
	{
		for (int k = 0; k<2; k++)
		{
			iscien[k][i] = 0;

			for (int j = N - 1; j>i; j--)
			{
				double xrzut;
				xrzut = x[k][j][0] + v[0] / v[1] * (x[k][i][1] - x[k][j][1]);

				if (k == 0 && xrzut<x[k][i][0])
					iscien[k][i] = 1;
				if (k == 1 && xrzut>x[k][i][0])
					iscien[k][i] = 1;
			}
		}
	}
}

void robimywzrost(int N, int** iscien, double** l, double wzrost_cien, double wzrost_slonce)
{
	// czyli: dla ka¿dej komórki zrób to, co tam jest w tym ifie (czyli przesuñ komórkê)

	for (int i = 1; i<N; i++) // dla ka¿dego punktu z tych N punktów
	{
		for (int k = 0; k < LiczbaPunktow; k++) // dla i lewego, i prawego
		{
			if (iscien[k][i] == 0)
				l[k][i] = l[k][i] * wzrost_slonce; // pomnó¿ wysokoœæ, na którym jest ten punkt, przez wzrost dla cienia lub s³oñca
			else
				l[k][i] = l[k][i] * wzrost_cien;
		}
	}
}

void przesuwanie(int N, double** l, double*** x, double a)
{
	for (int i = 1; i<N; i++)
	{
		//std::cout << "punkt nr = " << i << std::endl;
		double* tenlewy = x[LewyPunkt][i];
		double* tenprawy = x[PrawyPunkt][i];
		double* poprzednilewy = x[LewyPunkt][i - 1];
		double* poprzedniprawy = x[PrawyPunkt][i - 1];

		int poprzedni = i - 1;
		int ten = i;


		double fi, jpr, beta;

		fi = (l[LewyPunkt][i] - l[PrawyPunkt][i]) / a; // ró¿nica wysokoœci miêdzy lewym a prawym punktem na i-tym poziomie
		//std::cout << "fi = " << fi << std::endl;
													   // dzielimy przez a -szerokosc stozka stala=1

		jpr = fi / l[PrawyPunkt][i]; // dzielimy przez wysokoœæ prawego punktu (czyli jakby normalizujemy)
		//std::cout << "jpr = " << jpr << std::endl;


		beta = atan(poprzednilewy[1] - poprzedniprawy[1])
					/ (poprzedniprawy[0] - poprzednilewy[0]);

		//std::cout << "beta = " << beta << std::endl;


		if (fabs(jpr) < 0.01) // ró¿nica wysokoœci punktu lewego i prawego jest niewielka
		{
			for (int lewy_prawy = 0; lewy_prawy < LiczbaPunktow; lewy_prawy++) // dla obu punktów // k to na pewno lewy lub prawy punkt
			{
				x[lewy_prawy][ten][1] = x[lewy_prawy][poprzedni][1] + l[lewy_prawy][i] * cos(beta); // wartoœæ z poprzedniego mno¿ymy przez (wysokoœæ * cosinus)
				x[lewy_prawy][ten][0] = x[lewy_prawy][poprzedni][0] + l[lewy_prawy][i] * sin(beta); // wartoœæ z poprzedniego mno¿ymy przez (wysokoœæ * sinus)
			}
		}
		else
		{
			double S0, S1;
			S0 = (poprzedniprawy[0] * (1 / jpr + a) - poprzednilewy[0] / jpr) / a; // do tego z cosinuem
			S1 = (poprzedniprawy[1] * (1 / jpr + a) - poprzednilewy[1] / jpr) / a; // do tego z sinusem
			//

			for (int k = 0; k<LiczbaPunktow; k++) // dla obu punktów
			{
				double rwodz;
				rwodz = sqrt(pow(x[k][i - 1][0] - S0, 2.0) + pow(x[k][i - 1][1] - S1, 2.0)); //
				x[k][i][0] = S0 - rwodz*cos(beta + fi);
				x[k][i][1] = S1 + rwodz*sin(beta + fi);
			}
		}
	}
}

void wykonajkolejnykrok(int N, double***& x, int**& iscien, double**& l, double*& v, double*& kunttab,
	int t, double a, double wzrost_cien, double wzrost_slonce)
{
	//wyszukiwanie cienia
	sprawdzamyczyjestcien(N, iscien, v, x); // wpisuje do iscien[], czy dla danego punktu jest cieñ
	//jest cien

	//sam wzrost, bez przesuwania
	robimywzrost(N, iscien, l, wzrost_cien, wzrost_slonce); // jeœli jest cieñ, to roœnie bardziej, jeœli nie, to roœnie trochê

	//jest wzrost

	//przesuwanie
	przesuwanie(N, l, x, a);

	//poprzesuwane
	double kunt;

	kunt = atan((x[1][N - 1][1] - x[1][N - 5][1]) / (x[1][N - 1][0] - x[1][N - 5][0])); // kunt jest tylko dla prawego
	kunttab[t] = kunt;
}

void zapiszdopliku(int N, int** iscien, double*** x)
{
	FILE* fp;
	fp = fopen( "Sto¿ek_wzrostu.txt", "w");
	for (int k = 0; k<2; k++)
	{
		for (int i = 0; i<N; i++)
		{
			fprintf(fp, "%f\t%f\t%d\n", x[k][i][0], x[k][i][1], iscien[k][i]);
		}
	}
	fclose(fp);

}

void zapiszkunt(double* kunttab, int liczbakrokowczasowych)
{
	FILE* fp;
	fp=fopen( "Kunt.txt", "w");
	for (int t = 0; t < liczbakrokowczasowych; t++)
	{
		fprintf(fp, "%d\t%f\n", t + 1, kunttab[t]);
	}
	fclose(fp);

}

bool argumenty(int liczbaargumentow, char* args[], double& a, double&l0, int&liczbakrokowczasowych, int& N)
{
	bool ok = true;
	if (liczbaargumentow == 1) // tylko nazwa pliku
	{
		a = domyslny_a;
		l0 = domyslny_l0;
		liczbakrokowczasowych = domyslny_krokczasowy;
		N = domyslnyN;
	}
	else if (domyslnaliczbaargumentow + 1 == liczbaargumentow)
	{
		sscanf(args[1], "%d", &liczbakrokowczasowych);
		sscanf(args[2], "%d", &N);
		sscanf(args[3], "%lf", &a);
		sscanf(args[4], "%lf", &l0);
	}
	else
	{
		ok = false;
		printf("Sposob uzycia: liczba krokow czasowych, N, a, l0\n\t(ale bez przecinkow - a double musza miec kropki)");
	}
	return ok;
}




