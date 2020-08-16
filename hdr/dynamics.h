struct _input load_inputfile();

struct _input{
  double Bp;
  double Prot0;
  double Mej;
  double tmin;
  double tmax;
  double fac_dt;
  double epsB;
};

void calc_dynamics(double Bp, double Prot0, double Mej, double tmin, double tmax, double fac_dt, double epsB);
void sd(double Bp, double Prot0, double Mej, double *tsd, double *Lsd, double *Erot, double *vsd, double *rsd);
double Prot(double t, double tsd, double Prot0);
double Lpsr(double t, double tsd, double lsd);
double Epsr(double t, double tsd, double Erot);
double vej(double t, double tsd, double vsd);
double rej(double t, double tsd, double rsd);
double vnb(double t, double tsd, double vsd);
double rnb(double t, double tsd, double rsd);
double Bnb(double t, double tsd, double Erot, double rsd, double epsB);

