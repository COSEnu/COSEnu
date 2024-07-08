#if not defined(__STRUCTURES__)
#define __STRUCTURES__

/*---------------------------------------------------------------------------*/
struct FieldVar {
    double* ee;
    double* xx;
    double* ex_re;
    double* ex_im;
    double* bee;
    double* bxx;
    double* bex_re;
    double* bex_im;

    FieldVar(int size) {
        ee     = new double[size];
        xx     = new double[size];
        ex_re  = new double[size];
        ex_im  = new double[size];

        bee    = new double[size];
        bxx    = new double[size];
        bex_re = new double[size];
        bex_im = new double[size];
    }
    ~FieldVar() {
        delete[] ee;
        delete[] xx;
        delete[] ex_re;
        delete[] ex_im;
        
        delete[] bee;
        delete[] bxx;
        delete[] bex_re;
        delete[] bex_im;
    }
};

/*---------------------------------------------------------------------------*/

struct Flux{
    FieldVar *rflux;
    FieldVar *lflux;
    Flux(int size){
        rflux = new FieldVar(size);
        lflux = new FieldVar(size);
    }
    ~Flux(){
        delete rflux;
        delete lflux;
    }
};
struct Pol{
    double *normP;
    double *P1;
    double *P2;
    double *P3;
    
    double *normbP;
    double *bP1;
    double *bP2;
    double *bP3;

    Pol(unsigned int size){
        normP = new double[size];
        P1    = new double[size];
        P2    = new double[size];
        P3    = new double[size];

        normbP = new double[size];
        bP1    = new double[size];
        bP2    = new double[size];
        bP3    = new double[size];
    }

    ~Pol(){
        delete[] normP;
        delete[] P1;
        delete[] P2;
        delete[] P3;

        delete[] normbP;
        delete[] bP1;
        delete[] bP2;
        delete[] bP3;
    }
};

/*---------------------------------------------------------------------------*/

struct M{

    unsigned int order;
    double e1;
    double e2;
    double e3;
    double norm;

    M(unsigned int order_){

        order = order_;
        e1 = 0.0;
        e2 = 0.0;
        e3 = 0.0;
        norm = 0.0;
    }
    ~M(){}
};

/*---------------------------------------------------------------------------*/

// struct Params{
//     unsigned int nz ;
//     unsigned int nvz;
//     unsigned int gz ;
//     double    z0, z1;
//     double CFL;
//     Params(){}
//     ~Params(){}
// };

/*---------------------------------------------------------------------------*/

struct Profile{
    double *G;
    double *bG;

    Profile(int size){
        G  = new double[size];
        bG = new double[size];
    }

    ~Profile(){
        delete[] G ;
        delete[] bG;
    }
};

/*---------------------------------------------------------------------------*/

#endif // __STRUCTURES__

/*---------------------------------------------------------------------------*/
