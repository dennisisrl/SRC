#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

using namespace std;

void ModMM1();
void ModMM1MM();
void ModMMK();
void ModMMKMM();
void analisisEconomico(double lambda, double miu, double Wq, double W, int k);
void CalcularCosto(double lambda, double miu, double Wq, double W, double HorDia, int k, int tipCostUni, int defTiem);

bool condEsta(double lambda, double miu, int k);

double cambiarValorTiempo(int HorDia, double costo);

double sumatoriaMM1(int a, int b, double c, double d, double e);
double sumatoriaMMK(double lambda, double miu, int k);
double sumatoriaMM1MM(double lambda, double miu, int M);
double sumatoriaPn1(double lambda, double miu, double P0, int min, int max);
double sumatoriaPn2(double lambda, double miu, double P0, int k, int min, int max);
double sumatoriaPn(int min, int max, double M, double lambda, double miu, double P0);
double sumatoriaMMKMM1(double lambda, double miu, int M, int k);
double sumatoriaMMKMM2(double lambda, double miu, int M, int k);
double sumatoriaPn1MMKMM(double lambda,double miu, double P0, int M, int min, int max);
double sumatoriaPn2MMKMM(double lambda, double miu, double P0, int M, int min, int max, int k);
double sumatoriaPn1MMKMM1(double lambda, double miu, double P0, int M, int max, int min);
double sumatoriaPn2MMKMM2(double lambda, double miu, double P0, int M, int max, int min, int k);
double factorial(int a);

//Menu Inicial
int main(){
    int op = 1;
    while (op != 0)
    {
        system("cmd /c cls"); //limpia la pantalla
        cout << "Escoge el modelo" << endl;
        cout << "-------------Menu------------" << endl;
        cout << "1) M/M/1" << endl;
        cout << "2) M/M/K" << endl;
        cout << "3) M/M/1/M/M" << endl;
        cout << "4) M/M/K/M/M" << endl;
        cout << "0) Salir" << endl;
        cout << "Opcion: ";
        cin >> op;
        cout << "-----------------------------" << endl;
        switch (op){
            case 1: ModMM1();
            break;
            case 2: ModMMK();
            break;
            case 3: ModMM1MM();
            break;
            case 4: ModMMKMM();
            break;
        }
    }
    return 0;
}
void ModMM1(){
    double lambda, miu, p, P0, L, Lq, Ln, W, Wq, Wn;
    double costUni, Cte, Ctse, Cs, Ct;
    double Ctte, Ctts, Cttse, Cts;
    int minClien = 0, nClien = 0, res, k = 1;
    system("cmd /c cls"); //limpia la pantalla
    //definir los valores de lambda y miu
    cout << "\tModelo M/M/1" << endl;
    cout << "Ingrese los Valores de " << endl;
    cout << "Lambda: "; 
    cin >> lambda;
    cout << "Miu: ";
    cin >> miu;
    
    if(condEsta(lambda, miu, k)){
        cout << "------------Formulacion Matematica--------------" << endl;
        cout << "Probabilidad de hallar el sistema ocupado" << endl;
        p = lambda/miu;
        cout << "\tp = " << p << endl;
        
        cout << "Probabilidad de hallar el sistema vaicio u ocioso" << endl;
        P0=1-p;
        cout << "\tP0 = " << P0 << endl;
        
        cout << "Probabilidad de hallar exactamente n clientes dentro del sistema" << endl;
        cout << "\tNumero minimo de Clientes" << endl;
        cin >> minClien;
        cout << "\tNumeor Maximo de Clientes" << endl;
        cin >> nClien;
        cout << "Respuesta: ";
        cout << sumatoriaMM1(minClien,nClien,lambda,miu,P0) << endl;
        
        cout << "Numero esperado de clientes en el sistema" << endl;
        L = lambda/(miu-lambda);
        cout << "\tL = " << L << endl;
        
        cout << "Numero esperado de clientes en la cola" << endl;
        Lq = pow(lambda,2)/(miu*(miu-lambda));
        cout << "\tLq = " << Lq <<endl;
        
        cout << "Numero esperado de clientes en la cola vacia" << endl;
        Ln = lambda/(miu - lambda);
        cout << "\tLn = " << Ln << endl;
        
        cout << "Tiempo esperado en el sistema" << endl;
        W = 1/(miu-lambda);
        cout << "\tW = " << W << endl;
        
        cout << "Tiempo esperado en cola" << endl;
        Wq = lambda/(miu*(miu-lambda));
        cout << "\tWq = " << Wq << endl;
        
        cout << "Tiempo esperado en colas no vacias" << endl;
        Wn = lambda/(miu-lambda);
        cout << "\tWn = " << Wn << endl;
        
        system("cmd /c PAUSE");
        while(res != 2){
            cout << "---------------------------------" << endl;
            cout << "\tRealizar Analisis Economico" << endl;
            cout << "1) Si" << endl;
            cout << "2) No" << endl;            
            cout << "Opcion: ";
            cout << "---------------------------------" << endl;
            cin >> res;
            if (res == 1)
                analisisEconomico(lambda, miu, Wq, W, k);
        }
    }
    else
    {
        cout << "Los datos ingresados no cumplen con la Condicion de Estabilidad" << endl;
    }
    

    system("cmd /c PAUSE"); //pausa la pantalla
}

double sumatoriaMM1(int min, int n, double lambda, double miu, double P0){
    double sum = 0;
    for(int i = min;i <= n;i++){
        sum += (P0*pow((lambda/miu),i));
    }
    return sum;
}

//Modelo M/M/K
void ModMMK(){
    int k, min, n, res;
    double lambda, miu, P0, Pk, Pne, Pn, L, Lq, Ln, W, Wq = 0, Wn;

    system("cmd /c cls"); //limpia la pantalla

    //definir los valores de lambda y miu
    cout << "\tModelo M/M/k" << endl;
    cout << "Ingrese los Valores de " << endl;
    cout << "Lambda: "; 
    cin >> lambda;
    cout << "Miu: ";
    cin >> miu;
    cout << "Numero de Servidores: ";
    cin >> k;

    if(condEsta(lambda, miu, k)){
        double mul;
        cout << "------------Formulacion Matematica--------------" << endl;
        mul = (1/factorial(k))*(pow(lambda/miu,k))*((k*miu)/((k*miu)-lambda));
        cout << "Probabilidad de hallar el sistema completamente vacio" << endl;
        P0 = 1/((sumatoriaMMK(lambda,miu,k))+mul);
        cout << "\tP0 = " << P0 << endl;

        cout << "Probabilidad de que un usuario que llega tenga que esperar" << endl;
        Pk = mul*P0;
        cout << "\tPk = " << Pk << endl;

        cout << "Probabilidad de que un usuario que llega no tenga que esperar" << endl;
        Pne = 1 - Pk;
        cout << "\tPne = " << Pne << endl;

        cout << "Probabilidad de Hallar exactamente n clientes" << endl;
        cout << "Cliente: ";
        cin >> n;
        cout << "Iniciando: ";
        cin >> min;

        if(n <= k){
            Pn = sumatoriaPn1(lambda,miu,P0,min,n);
        }
        else{
            Pn = sumatoriaPn2(lambda,miu,P0,k,min,n);
        }
        cout << "\tPn = " << Pn << endl;

        cout << "Numero esperado de clientes en el sistema" << endl;
        L = ((lambda*miu*pow((lambda/miu),k))/((factorial(k-1)*(pow(((k*miu)-lambda),2))))*P0)+(lambda/miu);
        cout << "\tL = " << L << endl;

        cout << "Numero esperado de clientes en la cola" << endl;
        Lq = ((lambda*miu)*(pow((lambda/miu),k)*P0))/((factorial(k-1))*(pow(((k*miu)-lambda),2)));
        cout << "\tLq = " << Lq << endl;

        cout << "Numero esperado de clientes en la cola no vacia" << endl;
        Ln = Lq/Pk;
        cout << "\tLn = " << Ln << endl;

        cout << "Tiempo esperado en el Sistema" << endl;
        W = ((miu*(pow((lambda/miu),k))*P0)/((factorial(k-1))*(pow(((k*miu)-lambda),2))))+(1/miu);
        cout << "\tW = " << W << endl;

        cout << "Tiempo esperado en cola" << endl;
        Wq = (miu*(pow((lambda/miu),k))*P0)/((factorial(k-1)*(pow(((k*miu)-lambda),2))));
        cout << "\tWq = " << Wq << endl;

        cout << "Tiempo esperado en cola para colas no vacias" << endl;
        Wn = Wq/Pk;
        cout << "\tWn = " << Wn << endl;

        system("cmd /c PAUSE");

        while(res != 2){
            cout << "---------------Realizar Analisis Economico------------" << endl;
            cout << "1) Si" << endl;
            cout << "2) No" << endl;
            cout << "---------------------------------" << endl;
            cout << "Opcion: ";
            cin >> res;
            if(res == 1)
                analisisEconomico(lambda, miu, Wq, W, k);
        }
    }
    else{
        cout << "Los valores ingresados no cumplen con la Condicion de Estabilidad" << endl;
    }

    system("cmd /c PAUSE");

}
double sumatoriaPn1(double lambda, double miu, double P0, int min, int max){
    double res = 0;
    for(int n = min; n <= max;n++){
        res = res + ((P0/factorial(n))*(pow((lambda/miu),n)));
    }
    return res;
}

double sumatoriaPn2(double lambda, double miu, double P0, int k, int min, int max){
    double res = 0;
    for(int n = min;n <= max; n++){
        res = res + (P0*(1/((factorial(k)*(pow(k,(n-k))))))*(pow((lambda/miu),n)));
    }
    return res;
}

double sumatoriaMMK(double lambda, double miu, int k){
    double resp = 0;
    for(int n = 0;n <= k-1;n++){
        resp = resp + ((1/(factorial(n)))*(pow((lambda/miu),n)));
    }

    return resp;
}

//Modelo M/M/1/M/M
void ModMM1MM(){
    int n, min, res,k = 1;
    double lambda, miu, M, P0, Pn, Pe, L, Lq, Ln, W, Wq, Wn;

    system("cmd /c cls"); //limpia la pantalla
    cout << "\tModelo M/M/1/M/M" << endl;
    cout << "    Defina los valores Iniciales" << endl;
    cout << "Lambda: "; 
    cin >> lambda;
    cout << "Miu: ";
    cin >> miu;
    cout << "Poblacion: ";
    cin >> M;

    cout << "------------Formulacion Matematica--------------" << endl;
    cout << "Probabilidad de Hallar el Sistema vacio u ocioso" << endl;
    P0 = sumatoriaMM1MM(lambda,miu,M);
    cout << "\tP0 = " << P0 << endl;

    cout << "Probabilidad de hallar el sistema ocupado, utilización del sistema, probabilidad que tienen los usuarios de esperar para ser atendidos";
    Pe = 1 - P0;
    cout << "\tPe = " << Pe << endl;

    cout << "Probabilidad de hallar exactamente n Clientes en el Sistema" << endl;
    cout << "Minimo: ";
    cin >> min;
    cout << "Clientes: ";
    cin >> n;
    Pn = sumatoriaPn(min, n, M, lambda, miu, P0);
    cout << "\tPn = " << Pn << endl;

    cout << "Numero esperado de Clientes en el Sistema" << endl;
    L = M - ((miu/lambda)*(1-P0));
    cout << "\tL = " << L << endl;

    cout << "Numero esperado de clientes en la Cola" << endl;
    Lq = M-(((lambda+miu)/lambda)*(1-P0));
    cout << "\tLq = " << Lq << endl;

    cout << "Numero esperado de clientes en la cola no vacia" << endl;
    Ln = Lq/Pe;
    cout << "\tLn = " << Ln << endl;

    cout << "Tiempo esperado en cola" << endl;
    Wq = Lq/((M-L)*lambda);
    cout << "\tWq = " << Wq << endl;

    cout << "Tiempo esperado en el Sistema" << endl;
    W = Wq+(1/miu);
    cout << "\tW = " << W << endl;

    cout << "Tiempo esperado en cola para colas no vacias" << endl;
    Wn = Wq/Pe;
    cout << "\tWn = " << Wn << endl;

    system("cmd /c PAUSE");

    while(res != 2){
        system("cmd /c cls");
        cout << "---------------------------------" << endl;
        cout << "\tRealizar Analisis Economico" << endl;
        cout << "1) Si" << endl;
        cout << "2) No" << endl;        
        cout << "Opcion: ";
        cout << "---------------------------------" << endl;
        cin >> res;
        if (res == 1)
            analisisEconomico(lambda, miu, Wq, W, k);        
    }
    
    system("cmd /c PAUSE");   
}

double sumatoriaMM1MM(double lambda, double miu, int M){
    double aux = 0;
    int fac1 = factorial(M);
    for(int n = 0;n <= M;n++){
        aux = aux + ((fac1/factorial(M-n))*pow((lambda/miu),n));
    }
    if(aux != 0){
        return 1/aux;
    }
    else{
        return 1;
    }
    
}
double sumatoriaPn(int min, int max, double M, double lambda, double miu, double P0){
    double aux = 0;
    int fac1 = factorial(M);
    for(int n = min;n <= max;n++){
        aux = aux + (((fac1/factorial(M-n))*pow((lambda/miu),n))*P0);
    }
    return aux;
}

//Modelo M/M/K/M/M
void ModMMKMM(){
     int k, min, n, M, res;
    double lambda, miu, P0, Pk, Pne, Pn, Pe, L, Lq, Ln, W, Wq = 0, Wn;

    system("cmd /c cls"); //limpia la pantalla

    //definir los valores de lambda y miu
    cout << "\tModelo M/M/k/M/M" << endl;
    cout << "Ingrese los Valores de " << endl;
    cout << "Lambda: "; 
    cin >> lambda;
    cout << "Miu: ";
    cin >> miu;
    cout << "Poblacion: ";
    cin >> M;
    cout << "Numero de Servidores: ";
    cin >> k;

        cout << "------------Formulacion Matematica--------------" << endl;
        cout << "Probabilidad de hallar el sistema completemente vacio" << endl;
        P0 = 1/(sumatoriaMMKMM1(lambda,miu,M,k)+sumatoriaMMKMM2(lambda,miu,M,k));
        cout << "\tP0 = " << P0 << endl;

        cout << "Probabilidad de hallar exactamente n Clientes dentro del sistema" << endl;
        cout << "Clientes: ";
        cin >> n;
        cout << "Iniciando: ";
        cin >> min;
        if((n >= 0)&&(n < k)){
            Pn = sumatoriaPn1MMKMM(lambda, miu, P0, M, n, min);
        }
        if((n >= k)&&(n <= M)){
            Pn = sumatoriaPn2MMKMM(lambda, miu, P0, M, n, min, k);
        }

        cout << "\tPn = " << Pn << endl;

        cout << "Probabilidad de hallar el sistema completamente ocupado" << endl;
        Pe = 1-sumatoriaPn1MMKMM(lambda,miu, P0,M,(k-1),0);
        //Pe = sumatoriaPn2MMKMM(lambda, miu, P0, M, M, k, k);
        cout << "\tPe = " << Pe << endl;

        cout << "Probabilidad de no esperar" << endl;
        Pne = 1 - Pe;
        cout << "\tPne = " << Pne << endl;

        cout << "Numero esperado de Clientes en el Sistema" << endl;
        L = sumatoriaPn1MMKMM1(lambda,miu,P0,M,(k-1),0)+sumatoriaPn2MMKMM2(lambda,miu,P0,M,M,k,k)+(k*(1-sumatoriaPn1MMKMM(lambda,miu,P0,M,(k-1),0)));
        cout << "\tL = " << L << endl;

        cout << "Numero espera de Clientes en la Cola" << endl;
        Lq = sumatoriaPn2MMKMM2(lambda,miu,P0,M,M,k,k);
        cout << "\tLq = " << Lq << endl;

        cout << "Numero esperado de Clientes en la Cola no Vacia" << endl;
        Ln = Lq/Pe;
        cout << "\tLn = " << Ln << endl;

        cout << "Tiempo esperado en Cola" << endl;
        Wq = Lq/((M-L)*lambda);
        cout << "\tWq = " << Wq << endl;
        
        cout << "Tiempo esperado en el Sistema" << endl;
        W = Wq+(1/miu);
        cout << "\tW = " << W << endl;

        cout << "Tiempo esperado en Cola para Colas no Vacias" << endl;
        Wn = Wq/Pe;
        cout << "\tWn = " << Wn << endl;

        system("cmd /c PAUSE");
        
        while(res != 2){
            system("cmd /c cls");
            cout << "\tRealizar Analisis Economico" << endl;
            cout << "1) Si" << endl;
            cout << "2) No" << endl;
            cout << "---------------------------------" << endl;
            cout << "Opcion: ";
            cin >> res;
            if (res == 1)
            analisisEconomico(lambda, miu, Wq, W, k);        
        }  

    system("cmd /c PAUSE");
}

double sumatoriaMMKMM1(double lambda, double miu, int M, int k){
    double resp = 0;
    int facM = factorial(M);
    for(int n = 0; n <= (k-1);n++){
        resp = resp+(((facM)/((factorial(M-n))*(factorial(n))))*(pow((lambda/miu),n)));
    }
    return resp;
}

double sumatoriaMMKMM2(double lambda, double miu, int M, int k){
    double resp = 0;
    int facM = factorial(M);
    int fack = factorial(k);
    
    for(int n = k;n <= M;n++){
        resp = resp + (((facM)/((factorial(M-n))*(fack)*(pow(k,(n-k)))))*(pow((lambda/miu),n)));
    }
    return resp;
}

double sumatoriaPn1MMKMM(double lambda, double miu, double P0, int M, int max, int min){
    double resp = 0;
    int facM = factorial(M);
    for(int n = min; n <= max;n++){
        resp = resp+((P0)*((facM)/((factorial(M-n))*(factorial(n))))*(pow((lambda/miu),n)));
    }
    return resp;
}

double sumatoriaPn2MMKMM(double lambda, double miu, double P0, int M, int max, int min, int k){
    double resp = 0;
    int facM = factorial(M);
    int fack = factorial(k);

    for(int n = min; n <= max; n++){
        resp = resp + ((P0)*((facM)/((factorial(M-n))*(fack)*(pow(k,(n-k)))))*(pow((lambda/miu),n)));
    }
    return resp;
}

double sumatoriaPn1MMKMM1(double lambda, double miu, double P0, int M, int max, int min){
    double resp = 0;
    int facM = factorial(M);
    for(int n = min; n <= max;n++){
        resp = resp+(n)*((P0)*((facM)/((factorial(M-n))*(factorial(n))))*(pow((lambda/miu),n)));
    }
    return resp;
}

double sumatoriaPn2MMKMM2(double lambda, double miu, double P0, int M, int max, int min, int k){
    double resp = 0;
    int facM = factorial(M);
    int fack = factorial(k);

    for(int n = min; n <= max; n++){
        resp = resp + (n-k)*((P0)*((facM)/((factorial(M-n))*(fack)*(pow(k,(n-k)))))*(pow((lambda/miu),n)));
    }
    return resp;
}
//Analisis Economico
void analisisEconomico(double lambda, double miu, double Wq, double W, int k){
    int defTiem, tipCostUni, HorDia;
    system("cmd /c cls");

    cout << "-------Anlisis Economico---------" << endl;
    cout << "Cantidad de Horas laborales al Dia: ";
    cin >> HorDia;

    system("cmd /c cls");
    while(tipCostUni != 0){
        cout << "-------------Que Costo va a Calcular------------------" << endl;
        cout << "1) Costo Unitario por Tiempo en Cola" << endl;
        cout << "2) Costo Unitario por Tiempo en el Sistema" << endl;
        cout << "3) Costo Unitario por Tiempo de Servicio" << endl;
        cout << "4) Costo Unitario por el Servidor" << endl;
        cout << "5) Costo Total Diario" << endl;
        cout << "0) Ninguno" << endl;
        cout << "Opcion: ";
        cin >> tipCostUni;
        cout << "------------------------------------------------------" << endl;
        system("cmd /c cls");
        if(tipCostUni != 0)
            CalcularCosto(lambda,miu,Wq,W, HorDia,k,tipCostUni,defTiem);
    }
    cout << endl;
}

void CalcularCosto(double lambda, double miu, double Wq, double W, double HorDia, int k, int tipCostUni, int defTiem){
    double costUni, Cte, Ctse, Cs;
    double Ctte, Ctts, Cttse, Cts, Ct;
    system("cmd /c cls");

    cout << "-----------Analisis Economico-----------" << endl;

    if(tipCostUni != 5){
            cout << "\t Costo unitario" << endl;
            cin >> costUni;
            costUni = cambiarValorTiempo(HorDia,costUni);
        }else{
            cout << "Costo Unitario por Tiempo en Cola: ";
            cin >> Cte;
            Cte = cambiarValorTiempo(HorDia,Cte);
            cout << "Costo Unitario por Tiempo en el Sistema: ";
            cin >> Cts;
            Cts = cambiarValorTiempo(HorDia,Cts);
            cout << "Costo Unitario por Tiempo de Servicio: ";
            cin >> Ctse;
            Ctse = cambiarValorTiempo(HorDia, Ctse);
            cout << "Costo Unitario por el Servidor: ";
            cin >> Cs;
            Cs = cambiarValorTiempo(HorDia,Cs);
        }
        switch (tipCostUni){
            case 1:
                cout << "Costo Diario por el Tiempo de Espera en Cola" << endl;
                Ctte = lambda*Wq*costUni;
                cout << "\tCtte = " << Ctte << endl;
            break;
            case 2:
                cout << "Costo Diario por el Tiempo en el Sistema" << endl;
                Ctts = lambda*W*costUni;
                cout << "\tCtts = " << Ctts << endl;
            break;
            case 3:
                cout << "Costo Diario por el Tiempo de Servicio" << endl;
                Cttse = lambda*(1/miu)*costUni;
                cout << "\tCttse = " << Cttse << endl;
            break;
            case 4:
                cout << "Costo Diario por el Servidor" << endl;
                Cts = k*costUni;
                cout << "\tCts = " << Cts << endl;
            break;
            default:
                Ctte = lambda*Wq*Cte;
                Ctts = lambda*W*Cts;
                Cttse = lambda*(1/miu)*Ctse;
                Cts = k*Cs;
                Ct = Ctte+Ctts+Cttse+Cts;
                cout << "Costo Total Diario del Sistema" << endl;
                cout << "\tCt = " << Ct << endl;
                break;
        }
        system("cmd /c PAUSE");
}

double factorial(int a){
    if (a == 0)
        return 1;
    else
        return a*factorial(a-1);
}

double cambiarValorTiempo(int HorDia, double costo){
    int defTiem;
    double resp = costo;
    if(costo > 0){
        cout << "---------------------------" << endl;
        cout << "\tTiempo del Dato" << endl;
        cout << "1) Minutos" << endl;
        cout << "2) Horas" << endl;
        cout << "3) Dias" << endl;
        cout << "Tiempo: ";
        cin >> defTiem;
        cout << "---------------------------" << endl;

        switch(defTiem){
           case 1:
            resp = 60*HorDia*costo;
            break;
            case 2:
            resp = HorDia*costo;
            break;
            default:
            resp = costo;
            break;
        }
    }
    return resp;
}

//Condicion de Estabilidad
bool condEsta(double lambda, double miu, int k){
    if((lambda/(k*miu)) < 1){
        return true;
    }else{
        return false;
    }    
}