struct _input load_inputfile();

struct _input{
  double Bp;
  double Prot0;
  double Mej;
  double tmin;
  double tmax;
  double fac_dt;
  double epsB;
  double epse;
  double gam_b;
  double gam_max;
  double p1;
  double p2;
  int Nbin_e;
  double gam_ph_min;
  double gam_ph_max;
  int Nbin_ph;
  int output_file_num;
};

void calc_microphysics(double epse, double gam_b, double gam_max, double p1, double p2, int Nbin_e, double gam_ph_min, double gam_ph_max, int Nbin_ph, int output_file_num);
int get_hydro_tstep();
void load_hydro_data(int t_step, double *t, double *dt, double *vej, double *rej, double *vnb, double *rnb, double *Bnb, double *Lpsr, double *dr);
double Ntot(double *dgam, double *dN_dgam, int Nbin_e);
void initialize_e_dis(double *gam, double *dgam, double *dN_dgam, double *dN_dgam_dt, double *dgam_dt, double *tad, double *tsyn, double gam_max, int Nbin_e);
void initialize_ph_dis(double *gam_ph, double *P_nu_syn, double *alpha_nu_syn, double gam_ph_min, double gam_ph_max, int Nbin_ph);
double dNe_dEe_dt_inj(double gam, double Lpsr, double epse, double gam_b, double gam_max, double p1, double p2, int Nbin_e);
void injection(double *gam, double *dgam, double *dN_dgam_dt, double *dgam_dt, double Lpsr, double dt, double *N_inj_tot, double epse, double gam_b, double gam_max, double p1, double p2, int Nbin_e);
double dgam_dt_ad(double gam, double t);
double dgam_dt_syn(double gam, double B);
double tsynb(double B, double gamb);
void cooling(double t, double B, double *dgam_dt, double *gam, double *tad, double *tsyn, int Nbin_e);
void time_evolution_e(double dt, double *gam, double *dgam, double *dN_dgam, double *dN_dgam_dt, double *dgam_dt, int Nbin_e);
double syn_func_fit(double x);
void calc_syn_spec(double B, double r, double dr, double *gam, double *dgam, double *dN_dgam, double gam_max, int Nbin_e, double *gam_ph, double *P_nu_syn, double *alpha_nu_syn, int Nbin_ph);
