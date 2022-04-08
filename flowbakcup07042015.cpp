#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <string>
#include <sys/time.h>

# define M_PI 3.14159265358979323846  /* pi */

using namespace std;


class WriteFile{

public:

static void writeMatrix(double **M, int qtdlinhas, int qtdColunas, char nomeArquivo[]){

     ofstream write;
     write.open(nomeArquivo);

     for (int i = 0; i < qtdlinhas; i++ ){

        write << "\n";

        for (int j = 0; j < qtdColunas; j++){

             write.width(10);
             write << M[i][j] << "        ";

        }
    }
             write.close();
    }


static void writeVetAgrup(double *V1, double *V2, int qtdlinhas, char nomeArquivo[]){

    double MatAgrVet[qtdlinhas][2];

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][0]=V1[i];
    }

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][1]=V2[i];
    }

     ofstream write;
     write.open(nomeArquivo);

     for (int i = 0; i < qtdlinhas; i++ ){

        write << "\n";

        for (int j = 0; j < 2; j++){

             write.width(11);
             write << MatAgrVet[i][j] << " ";

        }
    }
             write.close();

}


static void writeVetor(double *V, int qtdColunas, char nomeArquivo[]){

     ofstream write;
     write.open(nomeArquivo);

     for (int i = 0; i < qtdColunas; i++ ){

            write << V[i] <<"\n";

        }
             write.close();
}

};

class Vetor{

    public:

    static void imprimir(double *V, int qtdColunas){

    for (int i = 0; i < qtdColunas; i++ ){

        cout << V[i] << "\t";
    }
        cout << endl;
    }


};

class Matrix{ // Classe para imprimir matrix USAR SÒ nome de classe maisculo

public:

    // Matrix::imprimir()

   static void  imprimir(double **M,int qtdLinhas, int qtdColunas){ // Nome de objeto minusculo

        for ( int i = 0; i < qtdLinhas; i++){

                Vetor::imprimir(M[i],qtdColunas);
        }
    }

};

class esc {

    public:

    int IMAX, JMAX, ILE, ITE;
    double XSF, YSF;

// Criando ponteiro duplos e simples para Alocação dinâmica de memoria e "zeramento" de vetores e matrizes na memoria do computador

    double **LPhi, **Cij, **Phi, **Npj, **u, **v, **uaer, **cp, **X, **Y;
    double  *x_i, *y_j, *dPhi_dy, *dPhi_dy_tresmeios, *LAx, *LBx, *LCx, * LAy, *LBy, *LCy;
    double *DELTAXX, *DELTAYY, *MD, *UP, *LD, *RHS, *A, *B, *C, *D, *cplinha;


void inicializaMatriz_LAx( int IMAX){

    int i;

    LAx = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (i = 0; i < IMAX; i++){

           LAx [i] = 0.0;
              }
}

void inicializaMatriz_LBx( int IMAX){

    int i;

    LBx = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (i = 0; i < IMAX; i++){

           LBx [i] = 0.0;
    }
}

void inicializaMatriz_LCx( int IMAX){

    int i;

    LCx = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (i = 0; i < IMAX; i++){

           LCx [i] = 0.0;
    }
}


void inicializaMatriz_LAy( int JMAX){

    int j;

    LAy = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           LAy [j] = 0.0;
    }
}

void inicializaMatriz_LBy( int JMAX){

    int j;

    LBy = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           LBy [j] = 0.0;
    }
}

void inicializaMatriz_LCy( int JMAX){
   int j;

    LCy = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           LCy [j] = 0.0;
    }
}


void inicializaMatriz_dPhi_dy_tresmeios(int ITE, int ILE, int IMAX, int JMAX, double deltaX, double th, double uinf){ // X do plano fisico é exatamente igual a Xi do plano computacional

// Calculando a malha para a direção x, onde  ILE<=i<=ITE

    int i;

    dPhi_dy_tresmeios = new double [IMAX];

    for (i = 0; i < IMAX; i++){

           dPhi_dy_tresmeios [i] = 0.0;

    }
}

void inicializaMatriz_dPhi_dy(int ITE, int ILE, int IMAX, double deltaX){ // X do plano fisico é exatamente igual a Xi do plano computacional

// Calculando a malha para a direção x, onde  ILE<=i<=ITE

    int i;

    dPhi_dy = new double[IMAX];

    for (i = 0; i < IMAX; i++){
        dPhi_dy [i] = 0.0;
    }
}


void inicializaMatriz_x_i(int ITE, int ILE, int IMAX, double deltaX){ // X do plano fisico é exatamente igual a Xi do plano computacional

// Calculando a malha para a direção x, onde  ILE<=i<=ITE

    int i;

    x_i = new double[IMAX];

    for (i = 0; i < IMAX; i++){
        x_i [i] = 0.0;
    }
}

void inicializaMatriz_y_j(int ITE, int ILE, int JMAX, double deltaX){ // X do plano fisico é exatamente igual a Xi do plano computacional

// Calculando a malha para a direção x, onde  ILE<=i<=ITE

    int j;

    y_j = new double[JMAX];

    for (j = 0; j < JMAX; j++){

            if (j == 0 || j == 1){
             y_j [0] = - deltaX/2.0 ;
             y_j [1] = deltaX/2.0;

                 }else{
                 y_j [j] = 0.0;
                 }
    }
}


void inicializaMatriz_LPhi(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  LPhi = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        LPhi [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 LPhi [i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Phi(int IMAX, int JMAX){ // POTENCIAL Zerando e ja assumindo chute inicial

  int i, j;

  Phi = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Phi [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Phi [i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Cij(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  Cij = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Cij [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                Cij [i][j] = 0.0;
             }
    }
}

void inicializaMatriz_Npj(int IMAX, int JMAX){

  int i, j;

  Npj = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Npj [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                Npj [i][j] = 0.0;
             }
    }
}


void inicializaMatriz_DELTAXX(int IMAX){

// Calculando a malha para a direção x, onde  ILE<=i<=ITE

    int i;

    DELTAXX = new double[IMAX];

    for (i = 0; i < IMAX; i++){
        DELTAXX [i] = 0.0;
    }
}

void inicializaMatriz_DELTAYY(int JMAX){

// Calculando a malha para a direção x, onde  ILE<=i<=ITE

    int j;

    DELTAYY = new double[JMAX];

    for (j = 0; j < JMAX; j++){
        DELTAYY [j] = 0.0;
    }
}


// Alocação de memoria e inicialização de variaiveis para montagem da Matrix Tridiagonal
void inicializaMatriz_MD( int JMAX){
   int j;

    MD = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           MD[j] = 0.0;
    }
}

void inicializaMatriz_UP( int JMAX){
   int j;

    UP = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           UP [j] = 0.0;
    }
}

void inicializaMatriz_LD( int JMAX){
   int j;

    LD = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           LD [j] = 0.0;
    }
}

void inicializaMatriz_RHS( int JMAX){
   int j;

    RHS = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           RHS[j] = 0.0;
    }
}


// Inicializando e alocando espaço para as variaveis calculadas no algoritmo de Thomas
void inicializaMatriz_A( int JMAX){
   int j;

    A = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           A[j] = 0.0;
    }
}

void inicializaMatriz_B( int JMAX){
   int j;

    B = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           B[j] = 0.0;
    }
}

void inicializaMatriz_C( int JMAX){
   int j;

    C = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           C[j] = 0.0;
    }
}

void inicializaMatriz_D( int JMAX){
   int j;

    D = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j <JMAX; j++){

           D[j] = 0.0;
    }
}

void inicializaMatriz_u(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  u = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        u [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                u [i][j] = 0.0;
             }
    }
}

void inicializaMatriz_v(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  v = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        v [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                v [i][j] = 0.0;
             }
    }
}

void inicializaMatriz_uaer(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  uaer = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        uaer [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                uaer [i][j] = 0.0;
             }
    }
}

void inicializaMatriz_cp(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  cp = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        cp [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                cp [i][j] = 0.0;
             }
    }
}

void inicializaMatriz_cplinha( int IMAX){
   int i;

    cplinha = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (i = 0; i < IMAX; i++){

           cplinha[i] = 0.0;
    }
}

void inicializaMatriz_X(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  X = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        X [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                X [i][j] = 0.0;
             }
    }
}

void inicializaMatriz_Y(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  Y = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Y [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                Y [i][j] = 0.0;
             }
    }
}

/*void inicializaMatriz_M( int JMAX){
   int j;

    M = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           M[j] = 0.0;
    }
}*/


void inicializarMatrizes(int ILE, int ITE, int IMAX, int JMAX, double XSF, double YSF, double th, double uinf, double alfa, double r, double eps, double deltaX){

    // Inicialização de matrizes para geração da malha computacional

    inicializaMatriz_LPhi(IMAX, JMAX);
    inicializaMatriz_Phi(IMAX, JMAX);
    inicializaMatriz_Cij(IMAX, JMAX);
    inicializaMatriz_Npj(IMAX, JMAX);

    inicializaMatriz_x_i(ITE, ILE, IMAX, deltaX);
    inicializaMatriz_y_j(ITE, ILE, JMAX, deltaX);

    inicializaMatriz_dPhi_dy(ITE, ILE, JMAX, deltaX);
    inicializaMatriz_dPhi_dy_tresmeios(ITE, ILE,IMAX, JMAX, deltaX, th, uinf);

    inicializaMatriz_LAx(IMAX);
    inicializaMatriz_LCx(IMAX);
    inicializaMatriz_LBx(IMAX);

    inicializaMatriz_LAy(JMAX);
    inicializaMatriz_LCy(JMAX);
    inicializaMatriz_LBy(JMAX);

    inicializaMatriz_DELTAXX(IMAX);
    inicializaMatriz_DELTAYY(JMAX);

    // Inicialização para Montagem da Tridiagonal
    inicializaMatriz_MD(JMAX);
    inicializaMatriz_UP(JMAX);
    inicializaMatriz_LD(JMAX);
    inicializaMatriz_RHS(JMAX);

    // Inicialização para o Algoritmo de THOMAS

    inicializaMatriz_A(JMAX);
    inicializaMatriz_B(JMAX);
    inicializaMatriz_C(JMAX);
    inicializaMatriz_D(JMAX);
    //inicializaMatriz_M(JMAX);


    // Gerar graficos de velocidade e Cp

    inicializaMatriz_u(IMAX,JMAX);
    inicializaMatriz_v(IMAX,JMAX);
    inicializaMatriz_uaer(IMAX,JMAX);
    inicializaMatriz_cp(IMAX,JMAX);
    inicializaMatriz_cplinha(IMAX);


    inicializaMatriz_X(IMAX,JMAX);
    inicializaMatriz_Y(IMAX,JMAX);

}

};

esc  malhacomputacional  (esc esc, int i, int j, int ILE, int ITE, int IMAX, int JMAX, double XSF, double YSF, double deltaX, double uinf, double th);
esc  calculocoeficientes (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX);
esc Jacobimetodo (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX);
esc GAUSSmetodo (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX);
esc SORmetodo (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX, double r);
esc LGS (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX);
esc SLOR (esc esc,int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX, double r);
double MaxResiduo (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX,int kiter, double resmax);
void TDMA(double* a, double* b, double* c, double* d, double *x, int n);
esc Resultados (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX, double uinf);
//void Menu(esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX, double r);


int main(int argc, char*argv[]){


// Contador de Tempo computacional para avaliação dos metodos

clock_t TempoInicial;
TempoInicial = clock();

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         GRID PARAMETERS                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//double res;
double  resmax = 0.0;
double  CMAX = 0.0;
int maxit = 10000;

double *resmaxiter, RMAX = 1.0;
resmaxiter = new double [maxit];

int i,j,kiter = 0;
int ILE = 11; //correspondente ao bordo de ataque do perfil
int ITE = 31; //correspondente ao bordo de fuga do perfil
int IMAX = 41; // número de pontos na direção x
int JMAX = 12; // número de pontos na direção y
double XSF = 1.25; // fator de estiramento ("Stretching") da malha para a direção x
double YSF = 1.25; // fator de estiramento ("Stretching") da malha para a direção y

//CALCULO DO r Médio baseando-se no calculo do r otimo para as coordenadas x e y;

//double rotimo, rotimox, rotimoy, px, py;
//
//px = cos(M_PI/(IMAX-1));
//py = cos(M_PI/(JMAX-1));
//rotimox =  2/(px*px)*(1 - sqrt(1-px*px));
//rotimoy = 2/(py*py)*(1 - sqrt(1-py*py));
//rotimo = (rotimox + rotimoy)/2;
//double r = rotimo;
//cout << rotimo << endl;

 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    DADOS DE ENTRADA    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

// Airfoil thickness ratio
double th = 0.05;
// Free Stream Velocity
double uinf = 1.0;

// Iterations parameters for Numerical Methods
double alfa = 2.0;
double r = 1.9; // Fator de Relaxação do SOR
//double omega = 1.7;

// Convergence Criterion
double eps = 0.100*pow(10,-3); // Limite maquina para o SOR -12

// Calculo do deltaX na direção X apenas no intervalo em cima do perfil biconvexo
double deltaX = 1.0/(ITE-ILE);

esc esc;
//esc.IMAX = IMAX;
//esc.JMAX = JMAX;

esc.inicializarMatrizes(ILE, ITE, IMAX, JMAX, XSF,YSF, th, uinf,alfa,r, eps,deltaX);


esc = malhacomputacional (esc, i, j, ILE, ITE, IMAX, JMAX,  XSF, YSF, deltaX, uinf, th); // Função que retorna a malha computacional
esc = calculocoeficientes (esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF, deltaX);

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Processo Iterativo    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//while (fabs(RMAX) > eps){

for (kiter = 0; kiter < maxit; kiter ++){

     //Atualização da condição de contorno C[i][0] = C[i][1]
    for ( i = 1; i < IMAX-1; i++){
        esc.Cij[i][0] = esc.Cij[i][1]; // Condições para a linha j =0 e j = 1
    }

 //Calculo da correção no potencial
    for ( i = 0; i < IMAX; i++){
        for (j = 0; j < JMAX; j++){
            esc.Phi[i][j] = esc.Phi[i][j] + esc.Cij[i][j];
       }
    }

 //Calculo e atualização do residuo
    for ( i = 1; i < IMAX-1; i++){
        for (j = 1; j < JMAX-1; j++){
            esc.LPhi[i][j] = esc.LAx[i]*esc.Phi[i-1][j] + esc.LBx[i]*esc.Phi[i][j] + esc.LCx[i]*esc.Phi[i+1][j] + esc.LAy[j]*esc.Phi[i][j-1] + esc.LBy[j]*esc.Phi[i][j] + esc.LCy[j]*esc.Phi[i][j+1];
        }
    }

  //Matrix::imprimir(esc.LPhi,IMAX,JMAX);

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Esquemas de Iteração   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//esc = Jacobimetodo(esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF,deltaX);
//esc = GAUSSmetodo(esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF,deltaX);
//esc = SORmetodo(esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF,deltaX,r);
//esc = LGS(esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF,deltaX);
esc = SLOR(esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF,deltaX,r);

RMAX = MaxResiduo(esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF, deltaX,kiter,resmax);
resmaxiter[kiter] =log10(RMAX);

        //Matrix::imprimir(esc.Phi,IMAX,JMAX);
        //Matrix::imprimir(esc.LPhi,IMAX,JMAX);
}

esc = Resultados (esc,i,j,ILE,ITE,IMAX,JMAX, XSF, YSF,deltaX,uinf);

        char nomeArquivo1[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto1/Caso 1/Convergencia Caso 1/USLORcs1.txt";
        WriteFile::writeMatrix(esc.u,IMAX,JMAX,nomeArquivo1);
        //WriteFile::writeVetor(resmaxiter,(maxit-1),nomeArquivo1);

       // Vetor::imprimir(resmaxiter,maxit-1);

double tempo_total_s = (clock()- TempoInicial)/(double)CLOCKS_PER_SEC;
cout << "Tempo Total de Execucao = " << tempo_total_s << endl;

return 0;

};

esc malhacomputacional (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX,double uinf, double th){

    for ( i = (ILE-1); i < ITE; i++){

    esc.x_i[i] = (i-(ILE-1))*deltaX;
    }


    for (i = ITE; i < IMAX; i++){

        esc.x_i[i] = esc.x_i[i-1] + (esc.x_i[i-1] -esc.x_i[i-2])*XSF;
    }

    for (i = ILE-2; i>= 0; i--){

     esc.x_i[i] = esc.x_i[i+1]+(esc.x_i[i+1]-esc.x_i[i+2])*XSF;
    }

    for (j = 2; j < JMAX; j++){

    esc.y_j[j] = esc.y_j[j-1] + (esc.y_j[j-1]-esc.y_j[j-2])*YSF;

    }

    //Vetor::imprimir(esc.y_j,JMAX);

/* Condições de contorno - uinf */

  for(j = 0; j < JMAX; j++){

    esc.Phi[0][j] = uinf*esc.x_i[0]; // INFLOW a valiar a questão do x[i]
    esc.Phi[IMAX-1][j] = uinf*esc.x_i[IMAX-1]; // OUTFLOW

  }

  for(i = 0; i < IMAX; i++){

    esc.Phi[i][JMAX-1] = uinf*esc.x_i[i]; // PARTE SUPERIOR DO DOMINIO

    if ( i >= (ILE-1) && i<=(ITE-1)){ // definição do contorno na parede

       esc.dPhi_dy_tresmeios[i] = uinf*th*2.0*(1.0 - 2.0*esc.x_i[i]);
       esc.Phi[i][0] = - (esc.y_j[1] - esc.y_j[0])*esc.dPhi_dy_tresmeios[i]; //       esc.Phi[i][0] = esc.Phi[i][1] - (esc.y_j[1] - esc.y_j[0])*esc.dPhi_dy_tresmeios[i];
     }
  }

        for ( i = 0; i < JMAX; i++){

            for( j = 0; j < JMAX; j++ ){

            esc.X[i][j] = esc.x_i[i];
            esc.Y[i][j] = esc.y_j[j];

            }
        }

//    for ( i = 1; i < IMAX-1; i++){
//
//            for( j = 1; j < JMAX-1; j++ ){
//
//            esc.Phi[i][j] = 5;
//
//            }
//    }


 return esc;

}

esc calculocoeficientes (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX){

    for (i = 1; i < IMAX-1; i++){ // Calculo para pontos internos da malha

   esc.LAx[i] = 2.0/(esc.x_i[i+1]- esc.x_i[i-1])*1.0/(esc.x_i[i]- esc.x_i[i-1]);
   esc.LCx[i] = 2.0/(esc.x_i[i+1]- esc.x_i[i-1])*1.0/(esc.x_i[i+1] - esc.x_i[i]);
   esc.LBx[i] = -(esc.LAx[i] + esc.LCx[i]);

   esc.DELTAXX[i] = (esc.x_i[i+1] - esc.x_i[i-1])/2.0;

   }

   for (j = 1; j < JMAX-1; j++){

    esc.LAy[j] =2.0/(esc.y_j[j+1] - esc.y_j[j-1])*1.0/(esc.y_j[j]- esc.y_j[j-1]);
    esc.LCy[j] =2.0/(esc.y_j[j+1] - esc.y_j[j-1])*1.0/(esc.y_j[j+1] - esc.y_j[j]);
    esc.LBy[j] = -(esc.LAy[j]+ esc.LCy[j]);

    esc.DELTAYY[j] = (esc.y_j[j+1] - esc.y_j[j-1])/2.0;
   }

    //Vetor::imprimir(esc.LBx,JMAX);

return esc;

}

esc Jacobimetodo (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX){

    for (i =1; i < IMAX-1;i++){
        for (j = 1; j < JMAX-1; j++){
          esc.Npj[i][j] = (-2.0/pow(((esc.x_i[i+1] - esc.x_i[i-1])/2),2) -2.0/pow(((esc.y_j[j+1]- esc.y_j[j-1])/2),2));
        }
   }


    for (i = 1; i <IMAX-1;i++){
       for (j = 1; j < JMAX-1;j++){
         esc.Cij[i][j] = -(esc.LPhi[i][j]/esc.Npj[i][j]);
            }
         }


    return esc;

}

esc GAUSSmetodo (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX){

double NSORX, NSORY, NSORX2, NSORY2;

 for ( i = 1; i < IMAX-1; i++){

          for (j = 1; j < JMAX-1; j++){

                //esc.Cij[0][j] = 0.0;
                //esc.Cij[i][0] = 0.0;

                NSORX = -2.0/(esc.DELTAXX[i]*esc.DELTAXX[i]);
                NSORY = -2.0/(esc.DELTAYY[j]*esc.DELTAYY[j]);

                NSORX2 = (esc.DELTAXX[i]*esc.DELTAXX[i]);
                NSORY2 = (esc.DELTAYY[j]*esc.DELTAYY[j]);

                esc.Cij[i][j] = ((-esc.LPhi[i][j] - esc.Cij[i-1][j]/NSORX2 - esc.Cij[i][j-1]/NSORY2)/(NSORX+NSORY));
         }
 }

return esc;

}

esc SORmetodo (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX, double r){

double NSORX, NSORY, NSORX2, NSORY2;

 for ( i = 1; i < IMAX-1; i++){

          for (j = 1; j < JMAX-1; j++){

                //esc.Cij[0][j] = 0.0;
                //esc.Cij[i][0] = 0.0;

                NSORX = -2.0/(esc.DELTAXX[i]*esc.DELTAXX[i]);
                NSORY = -2.0/(esc.DELTAYY[j]*esc.DELTAYY[j]);

                NSORX2 = (esc.DELTAXX[i]*esc.DELTAXX[i]);
                NSORY2 = (esc.DELTAYY[j]*esc.DELTAYY[j]);

                esc.Cij[i][j] = r*((-esc.LPhi[i][j] - esc.Cij[i-1][j]/NSORX2 - esc.Cij[i][j-1]/NSORY2)/(NSORX+NSORY));
         }
 }

return esc;

}

esc LGS (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX){

     double *x; // Armazenando espaço e zerando o vetor solução do Thomas
     x = new double [JMAX];

     for (i = 0; i < JMAX; i++){

        x[i] =0.0;
     }

    int jj;
    //int ii;

   for ( i = 1; i < IMAX-1; i++){

        //esc.Cij[i][0] = esc.Cij[i][1]; // Condições para a linha j = 1 e j = 2

        for (j = 1; j < JMAX-1; j++){ //

        //esc.Cij[0][j] = 0.0;
        //esc.Cij[i][0] = 0.0;
        //esc.Cij[i][JMAX-1] = 0.0;


        esc.A[j-1] = esc.LAy[j]; // LOW DIAGONAL
        esc.B[j-1] = (esc.LBy[j] -2.0/(esc.DELTAXX[i]*esc.DELTAXX[i])); // MAIN DIAGONAL
        esc.C[j-1] = esc.LCy[j]; // UPER DIAGONAL
        esc.D[j-1] = -esc.LPhi[i][j] -(1.0/(esc.DELTAXX[i]*esc.DELTAXX[i]))*esc.Cij[i-1][j]; //RHS

        if ( j == 1){
        esc.A[0] = 0.0;
        }

        if (j == JMAX-2){
        esc.C[JMAX-2] = 0.0;
        }

        //Vetor::imprimir(esc.A,JMAX);
        //Vetor::imprimir(esc.A,JMAX);
        //Vetor::imprimir(esc.B,JMAX);
        //Vetor::imprimir(esc.C,JMAX);
        //Vetor::imprimir(esc.D,JMAX);
        //Vetor::imprimir(esc.Cij,IMAX,JMAX);
        //Matrix::imprimir(esc.LPhi,IMAX,JMAX);
    }

            TDMA(esc.A, esc.B, esc.C, esc.D, x,(JMAX-2));
        //ThomasAlgorithm(JMAX-2,esc.A,esc.B,esc.C,x, esc.D);
        // Atualização da condição de contorno C[i][0] = C[i][1]

            for ( jj = 1; jj < JMAX-1; jj++){
            esc.Cij[i][jj] = x[jj-1];
            }

}

    return esc;

}

esc SLOR (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX, double r){

     double *x; // Armazenando espaço e zerando o vetor solução do Thomas
     x = new double [JMAX];

     for (i = 0; i < JMAX; i++){

        x[i] =0.0;
     }

    int jj;
    //int ii;

   for ( i = 1; i < IMAX-1; i++){

        //esc.Cij[i][0] = esc.Cij[i][1]; // Condições para a linha j = 1 e j = 2

        for (j = 1; j < JMAX-1; j++){ //

        //esc.Cij[0][j] = 0.0;
        //esc.Cij[i][0] = 0.0;
        //esc.Cij[i][JMAX-1] = 0.0;


        esc.A[j-1] = (1.0/r)*esc.LAy[j]; // LOW DIAGONAL
        esc.B[j-1] = (1.0/r)*(esc.LBy[j] -2.0/(esc.DELTAXX[i]*esc.DELTAXX[i])); // MAIN DIAGONAL
        esc.C[j-1] = (1.0/r)*esc.LCy[j]; // UPER DIAGONAL
        esc.D[j-1] = -esc.LPhi[i][j] -1.0/(esc.DELTAXX[i]*esc.DELTAXX[i])*esc.Cij[i-1][j]; //RHS

        if ( j == 1){
        esc.A[0] = 0.0;
        }

        if (j == JMAX-2){
        esc.C[JMAX-2] = 0.0;
        }

    }

      TDMA(esc.A, esc.B, esc.C, esc.D, x,(JMAX-2));

            for ( jj = 1; jj < JMAX-1; jj++){
            esc.Cij[i][jj] = x[jj-1];
            }
}

    return esc;

}

double MaxResiduo (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX,int kiter, double resmax){

      for (i = 1; i < IMAX-1; i++){

        for (j = 1; j < JMAX-1; j++ ){

            if (fabs(esc.LPhi[i][j])> resmax){
              resmax = fabs(esc.LPhi[i][j]);
            }
        }
      }

return resmax;

}

void TDMA(double* a, double* b, double* c, double* d, double *x, int n){

    int k;
    double m;

    for ( k = 1 ; k <= n-1; k++){

    m = a[k]/b[k-1];
    b[k] = b[k] -m*c[k-1];
    d[k] = d[k] -m*d[k-1];

    }

    x[n-1] = d[n-1]/b[n-1];

    for (k = n-1; k >= 0; k--){

        x[k] = (d[k] - c[k]*x[k+1])/b[k];
    }

return;

}

esc Resultados (esc esc, int i, int j, int ILE, int ITE, int IMAX,int JMAX, double XSF, double YSF, double deltaX, double uinf){

            // Gerando os graficos de Velocidade para Calculo do Cp

    esc.u[0][0] = esc.Phi[0][0]/esc.x_i[0];
    esc.v[0][0] = esc.Phi[0][0]/esc.y_j[0];

    esc.u[IMAX-1][JMAX-1] = esc.Phi[IMAX-1][JMAX-1]/esc.x_i[IMAX-1]; // velocidade u no outflow
    esc.v[IMAX-1][JMAX-1] = esc.Phi[IMAX-1][JMAX-1]/esc.y_j[JMAX-1];

    esc.u[0][JMAX-1] = esc.Phi[0][JMAX-1]/esc.x_i[0];
    esc.v[0][JMAX-1] = esc.Phi[0][JMAX-1]/esc.y_j[JMAX-1];

    esc.u[IMAX-1][0] = esc.Phi[IMAX-1][0]/esc.x_i[IMAX-1];
    esc.v[IMAX-1][0] = esc.Phi[IMAX-1][0]/esc.y_j[0];

    for (i = 1; i< IMAX-1; i++){

        for ( j = 1; j < JMAX-1; j++){
         esc.u[i][0]=(esc.Phi[i+1][0] - esc.Phi[i-1][0])/(esc.x_i[i+1]-esc.x_i[i-1]);
         esc.v[i][0] = (esc.Phi[i][j+1] -esc.Phi[i][j])/(esc.y_j[j]- esc.y_j[j-1]);
        }
    }

    for (i = 1; i< IMAX-1; i++){

        for ( j = 1; j < JMAX-1; j++){
            // Inflow i = 0
            esc.u[0][j] = esc.Phi[0][j]/esc.x_i[0]; // velocidade u no inflow
            esc.v[0][j] = esc.Phi[0][j]/esc.y_j[j];

            // Outflow
            esc.u[IMAX-1][j] = esc.Phi[IMAX-1][j]/esc.x_i[IMAX-1]; // velocidade u no outflow
            esc.v[IMAX-1][j] = esc.Phi[IMAX-1][j]/esc.y_j[j];

            // Topo do dominio computacional
            esc.u[i][JMAX-1] = (esc.Phi[i+1][JMAX-1] - esc.Phi[i-1][JMAX-1])/(esc.x_i[i+1]-esc.x_i[i-1]); //esc.Phi[i][JMAX-1]/esc.x_i[i]; // velocidade u no inflow
            esc.v[i][JMAX-1] = esc.Phi[i][JMAX-1]/esc.y_j[JMAX-1];

            // Velocidade u em j =0

            esc.u[i][j] = (esc.Phi[i+1][j] - esc.Phi[i-1][j])/(esc.x_i[i+1]-esc.x_i[i-1]);
            esc.v[i][j] = (esc.Phi[i][j+1]- esc.Phi[i][j-1])/(esc.y_j[j+1]-esc.y_j[j-1]);
                }

            }

        for (i = 1; i< IMAX-1; i++){

            if ( i >= (ILE-1) && i<=(ITE-1)){ // definição do contorno na parede
                //esc.dPhi_dy_tresmeios[i] = uinf*th*2.0*(1.0 - 2.0*esc.x_i[i]);
                esc.cplinha[i] = - (1.0 - (((esc.u[i][1] + esc.u[i][0])/2.0)*((esc.u[i][1] + esc.u[i][0])/2.0) + (uinf*esc.dPhi_dy_tresmeios[i])*(uinf*esc.dPhi_dy_tresmeios[i]))/(uinf*uinf));
            }

            for ( j = 1; j < JMAX-1; j++){
            esc.cp[i][j] = 1.0 - (esc.u[i][j]*esc.u[i][j] + esc.v[i][j]*esc.v[i][j])/(uinf*uinf);
            esc.cp[i][j] = -esc.cp[i][j];

            }
        }

        Vetor::imprimir(esc.cplinha,IMAX);

        //char nomeArquivo1[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto1/Caso 1/Convergencia Caso 1/cpSLORcs2.txt";
//        WriteFile::writeMatrix(esc.cp,IMAX,JMAX,nomeArquivo1);
        //WriteFile::writeVetor(esc.cplinha,IMAX,nomeArquivo1);

return esc;
}
