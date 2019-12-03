#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../hdr/astro_const.h"
#include "../hdr/ns_pwn_param.h"
#include "../hdr/microphysics.h"

int main()
{
  struct _input in = load_inputfile();
  calc_microphysics(in.epse,in.gam_b,in.gam_max,in.p1,in.p2,in.Nbin_e,in.output_file_num);

  return 0;
}

struct _input load_inputfile()
{
  double dmy[14];
  FILE *ip;
  ip = fopen("../input.dat","r");
  int i=0;
  while (fscanf(ip,"%*s %lf %*[\n]",&dmy[i])!=EOF){
    i++;
  }
  fclose(ip);

  struct _input in;
  in.epse = dmy[7];
  in.gam_b = dmy[8];
  in.gam_max = dmy[9];
  in.p1 = dmy[10];
  in.p2 = dmy[11];
  in.Nbin_e = (int)dmy[12];
  in.output_file_num = (int)dmy[13];

  return in;
}

void calc_microphysics(double epse, double gam_b, double gam_max, double p1, double p2, int Nbin_e, int output_file_num)
{
  int i,j;

  int t_step = get_hydro_tstep();
  //double t[t_step],dt[t_step],vej[t_step],rej[t_step],vnb[t_step],rnb[t_step],Bnb[t_step],Lpsr[t_step];
  double *t,*dt,*vej,*rej,*vnb,*rnb,*Bnb,*Lpsr;
  t = (double *)malloc(t_step*sizeof(double));
  dt = (double *)malloc(t_step*sizeof(double));
  vej = (double *)malloc(t_step*sizeof(double));
  rej = (double *)malloc(t_step*sizeof(double));
  vnb = (double *)malloc(t_step*sizeof(double));
  rnb = (double *)malloc(t_step*sizeof(double));
  Bnb = (double *)malloc(t_step*sizeof(double));
  Lpsr = (double *)malloc(t_step*sizeof(double));
  load_hydro_data(t_step,t,dt,vej,rej,vnb,rnb,Bnb,Lpsr);

  double gam[Nbin_e],dN_dgam[Nbin_e],dgam[Nbin_e],dN_dgam_dt[Nbin_e],dgam_dt[Nbin_e],tad[Nbin_e],tsyn[Nbin_e];
  double N_inj_tot=0.;
  initialize_e_dis(gam,dgam,dN_dgam,dN_dgam_dt,dgam_dt,tad,tsyn,gam_max,Nbin_e);

  FILE *op;
  int output_int = (double)t_step/(double)output_file_num;
  char output_file_name[256]={"\0"},path[256]="../output/",head[256]="ele_dis_",dat[256]=".dat";

  for (j=0;j<t_step;j++){

    if ((j % output_int) == 0){
      sprintf(output_file_name,"%s%s%d%s",path,head,(int)((double)j/(double)output_int),dat);
      op = fopen(output_file_name,"w");
      fprintf(op,"#t = %lf [s], N conservation: %1.3e \n",t[j],Ntot(dgam,dN_dgam,Nbin_e)/N_inj_tot);
      fprintf(op,"#gam, dgam, dN/dgam, Ee*dN, dN/dgam/dt, dgam/dt, tad[s], tsyn[s] \n");

      for (i=0;i<Nbin_e;i++){
	fprintf(op,"%le %le %le %le %le %le %le %le \n",
		gam[i],dgam[i],dN_dgam[i],gam[i]*dgam[i]*dN_dgam[i]*MeC2,dN_dgam_dt[i],dgam_dt[i],tad[i],tsyn[i]);
      }

    }
    fclose(op);

    injection(gam,dgam,dN_dgam_dt,Lpsr[j],dt[j],&N_inj_tot,epse,gam_b,gam_max,p1,p2,Nbin_e);
    cooling(t[j],Bnb[j],dgam_dt,gam,tad,tsyn,Nbin_e);
    time_evolution(dt[j],gam,dgam,dN_dgam,dN_dgam_dt,dgam_dt,Nbin_e);
  }
}

int get_hydro_tstep()
{
  FILE *ip;

  int t_step;
  ip = fopen("../output/dynamics.log","r");
  fscanf(ip,"%*[^\n]%d",&t_step);
  fclose(ip);

  return t_step;
}


void load_hydro_data(int t_step, double *t, double *dt, double *vej, double *rej, double *vnb, double *rnb, double *Bnb, double *Lpsr)
{
  FILE *ip;

  int i;
  ip = fopen("../output/dynamics.dat","r");
  fscanf(ip,"%*[^\n]");
  for (i=0;i<t_step;i++){
    fscanf(ip,"%*d %*le %le %le %le %le %le %le %le %le\n",
	   &t[i],&dt[i],&vej[i],&rej[i],&vnb[i],&rnb[i],&Bnb[i],&Lpsr[i]);
  }
  fclose(ip);
}

double Ntot(double *dgam, double *dN_dgam, int Nbin_e)
{
  int i;
  double tmp=0.;
  for (i=0;i<Nbin_e;i++){
    tmp += dN_dgam[i]*dgam[i];
  }
  return tmp;
}


void initialize_e_dis(double *gam, double *dgam, double *dN_dgam, double *dN_dgam_dt, double *dgam_dt, double *tad, double *tsyn, double gam_max, int Nbin_e)
{
  double dln_gam = log(gam_max)/(double)(Nbin_e-1);
  int i;
  for (i=0;i<Nbin_e;i++){
    gam[i] = exp(dln_gam*i);
    dgam[i] = gam[i]*(exp(dln_gam)-1.);
    dN_dgam[i] = 0.;
    dN_dgam_dt[i] = 0.;
    dgam_dt[i] = 0.;
    tad[i] = 0.;
    tsyn[i] = 0.;
  }
}

double dN_dgam_dt_inj(double gam, double Lpsr, double epse, double gam_b, double gam_max, double p1, double p2)
{
  double fac_bol = 1./(2.-p1)*(1.-pow(gam_b,-2.+p1)) + 1./(2.-p2)*(pow(gam_max,2.-p2)-1.);

  if (gam < gam_b){
    return epse*Lpsr/(MeC2*gam_b*gam_b)/fac_bol*pow(gam/gam_b,-p1);
  } else {
    return epse*Lpsr/(MeC2*gam_b*gam_b)/fac_bol*pow(gam/gam_b,-p2);
  }
}

void injection(double *gam, double *dgam, double *dN_dgam_dt, double Lpsr, double dt, double *N_inj_tot, double epse, double gam_b, double gam_max, double p1, double p2, int Nbin_e)
{
  int i;
  double tmp = 0.;
  for (i=0;i<Nbin_e;i++){
    dN_dgam_dt[i] = dN_dgam_dt_inj(gam[i],Lpsr,epse,gam_b,gam_max,p1,p2);
    tmp += dN_dgam_dt[i]*dt*dgam[i];
  }
  *N_inj_tot += tmp;
}

double dgam_dt_ad(double gam, double t)
{
  return gam/t;
}

double dgam_dt_syn(double gam, double B)
{
  // electron synchrotron energy loss rate (see e.g., Eq. 7.13 of Dermer & Menon)
  // double sin2phi = 2.0/3.0; /* averaging pitch angle */
  // double beta_par = 1.0; /* assuming that particles are relativistic */
    
  return 4.0/3.0*C*SIGMA_T*(B*B/8.0/M_PI)*gam*gam;
}

void cooling(double t, double B, double *dgam_dt, double *gam, double *tad, double *tsyn, int Nbin_e)
{
  int i;
  for (i=0;i<Nbin_e;i++) {
    dgam_dt[i] = dgam_dt_ad(gam[i],t)+dgam_dt_syn(gam[i],B);
    tad[i] = gam[i]/dgam_dt_ad(gam[i],t);
    tsyn[i] = gam[i]/dgam_dt_syn(gam[i],B); 
  }
}

void time_evolution(double dt, double *gam, double *dgam, double *dN_dgam, double *dN_dgam_dt, double *dgam_dt, int Nbin_e)
{
  int i;
    
  dN_dgam[Nbin_e-1] = (dN_dgam[Nbin_e-1]+dN_dgam_dt[Nbin_e-1]*dt)/(1.0+dt/dgam[Nbin_e-1]*dgam_dt[Nbin_e-1]);
  for(i=Nbin_e-2;i>0;i--){
    dN_dgam[i] = (dN_dgam[i]+dN_dgam_dt[i]*dt+dN_dgam[i+1]*dt/dgam[i]*dgam_dt[i+1])/(1.0+dt/dgam[i]*dgam_dt[i]);
  }
  dN_dgam[0] = dN_dgam[0]+dN_dgam_dt[0]*dt+(dN_dgam[1]*dt/dgam[1]*dgam_dt[1])/(1.0+dt/dgam[0]*dgam_dt[0]);    
}
