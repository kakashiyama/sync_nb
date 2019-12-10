#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../hdr/astro_const.h"
#include "../hdr/ns_pwn_param.h"
#include "../hdr/microphysics.h"

int main()
{
  struct _input in = load_inputfile();
  calc_microphysics(in.epse,in.gam_b,in.gam_max,in.p1,in.p2,in.Nbin_e,in.gam_ph_min,in.gam_ph_max,in.Nbin_ph,in.output_file_num);

  return 0;
}

struct _input load_inputfile()
{
  double dmy[17];
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
  in.gam_ph_min = dmy[13];
  in.gam_ph_max = dmy[14];
  in.Nbin_ph = (int)dmy[15];
  in.output_file_num = (int)dmy[16];

  return in;
}

void calc_microphysics(double epse, double gam_b, double gam_max, double p1, double p2, int Nbin_e, double gam_ph_min, double gam_ph_max, int Nbin_ph, int output_file_num)
{
  int i,j;

  int t_step = get_hydro_tstep();
  double *t,*dt,*vej,*rej,*vnb,*rnb,*Bnb,*Lpsr,*dr;
  t = (double *)malloc(t_step*sizeof(double));
  dt = (double *)malloc(t_step*sizeof(double));
  vej = (double *)malloc(t_step*sizeof(double));
  rej = (double *)malloc(t_step*sizeof(double));
  vnb = (double *)malloc(t_step*sizeof(double));
  rnb = (double *)malloc(t_step*sizeof(double));
  Bnb = (double *)malloc(t_step*sizeof(double));
  Lpsr = (double *)malloc(t_step*sizeof(double));
  dr = (double *)malloc(t_step*sizeof(double));
  load_hydro_data(t_step,t,dt,vej,rej,vnb,rnb,Bnb,Lpsr,dr);

  double gam[Nbin_e],dN_dgam[Nbin_e],dgam[Nbin_e],dN_dgam_dt[Nbin_e],dgam_dt[Nbin_e],tad[Nbin_e],tsyn[Nbin_e];
  double N_inj_tot=0.;
  initialize_e_dis(gam,dgam,dN_dgam,dN_dgam_dt,dgam_dt,tad,tsyn,gam_max,Nbin_e);

  double gam_ph[Nbin_ph],P_nu_syn[Nbin_ph],alpha_nu_syn[Nbin_ph];
  initialize_ph_dis(gam_ph,P_nu_syn,alpha_nu_syn,gam_ph_min,gam_ph_max,Nbin_ph);

  FILE *op;
  int output_int = (double)t_step/(double)output_file_num;
  char output_file_name_e[256]={"\0"},path[256]="../output/",head_e[256]="ele_dis_",dat[256]=".dat";
  char output_file_name_ph[256]={"\0"},head_ph[256]="ph_spec_";

  for (j=0;j<t_step;j++){

    if ((j % output_int) == 0){
      calc_syn_spec(Bnb[j],rnb[j],dr[j],gam,dgam,dN_dgam,gam_max,Nbin_e,gam_ph,P_nu_syn,alpha_nu_syn,Nbin_ph);

      sprintf(output_file_name_e,"%s%s%d%s",path,head_e,(int)((double)j/(double)output_int),dat);
      op = fopen(output_file_name_e,"w");
      fprintf(op,"#t = %12.3e [s], r_nb = %12.3e [cm], B_nb = %12.3e [G], L_sd = %12.3e [erg/s]\n",t[j],rnb[j],Bnb[j],Lpsr[j]);
      fprintf(op,"#Number conservation: %1.3e\n",Ntot(dgam,dN_dgam,Nbin_e)/N_inj_tot);
      fprintf(op,"#gam, dgam, dN/dgam, Ee*dN, dN/dgam/dt, dgam/dt, tad[s], tsyn[s] \n");
      for (i=0;i<Nbin_e;i++){
	fprintf(op,"%le %le %le %le %le %le %le %le \n",
		gam[i],dgam[i],dN_dgam[i],gam[i]*dgam[i]*dN_dgam[i]*MeC2,dN_dgam_dt[i],dgam_dt[i],tad[i],tsyn[i]);
      }
      fclose(op);

      sprintf(output_file_name_ph,"%s%s%d%s",path,head_ph,(int)((double)j/(double)output_int),dat);
      op = fopen(output_file_name_ph,"w");
      fprintf(op,"#t = %12.3e [s], r_nb = %12.3e [cm], B_nb = %12.3e [G], L_sd = %12.3e [erg/s]\n",t[j],rnb[j],Bnb[j],Lpsr[j]);
      fprintf(op,"#nu [GHz], Pnu [erg/s/Hz] \n");
      for (i=0;i<Nbin_ph;i++){
	fprintf(op,"%le %le \n",gam_ph[i]*MeC2/H/1.0e9,P_nu_syn[i]);
      }
      fclose(op);

    }

    injection(gam,dgam,dN_dgam_dt,Lpsr[j],dt[j],&N_inj_tot,epse,gam_b,gam_max,p1,p2,Nbin_e);
    cooling(t[j],Bnb[j],dgam_dt,gam,tad,tsyn,Nbin_e);
    time_evolution_e(dt[j],gam,dgam,dN_dgam,dN_dgam_dt,dgam_dt,Nbin_e);

  }
  free(t);
  free(dt);
  free(vej);
  free(rej);
  free(vnb);
  free(rnb);
  free(Bnb);
  free(Lpsr);
  free(dr);
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


void load_hydro_data(int t_step, double *t, double *dt, double *vej, double *rej, double *vnb, double *rnb, double *Bnb, double *Lpsr, double *dr)
{
  FILE *ip;

  int i;
  ip = fopen("../output/dynamics.dat","r");
  fscanf(ip,"%*[^\n]");
  for (i=0;i<t_step;i++){
    fscanf(ip,"%*d %*le %le %le %le %le %le %le %le %le %le\n",
	   &t[i],&dt[i],&vej[i],&rej[i],&vnb[i],&rnb[i],&Bnb[i],&Lpsr[i],&dr[i]);
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

void initialize_ph_dis(double *gam_ph, double *P_nu_syn, double *alpha_nu_syn, double gam_ph_min, double gam_ph_max, int Nbin_ph)
{
  int i;
  double del_ln_gam_ph = (log(gam_ph_max)-log(gam_ph_min))/(double)(Nbin_ph-1);
  for (i=0;i<Nbin_ph;i++){
    gam_ph[i] = gam_ph_min*exp(del_ln_gam_ph*(double)i);
    P_nu_syn[i] = 0.;
    alpha_nu_syn[i] = 0.;
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

void time_evolution_e(double dt, double *gam, double *dgam, double *dN_dgam, double *dN_dgam_dt, double *dgam_dt, int Nbin_e)
{
  int i;
    
  dN_dgam[Nbin_e-1] = (dN_dgam[Nbin_e-1]+dN_dgam_dt[Nbin_e-1]*dt)/(1.0+dt/dgam[Nbin_e-1]*dgam_dt[Nbin_e-1]);
  for(i=Nbin_e-2;i>0;i--){
    dN_dgam[i] = (dN_dgam[i]+dN_dgam_dt[i]*dt+dN_dgam[i+1]*dt/dgam[i]*dgam_dt[i+1])/(1.0+dt/dgam[i]*dgam_dt[i]);
  }
  dN_dgam[0] = dN_dgam[0]+dN_dgam_dt[0]*dt+(dN_dgam[1]*dt/dgam[1]*dgam_dt[1])/(1.0+dt/dgam[0]*dgam_dt[0]);    
}

double syn_func_fit(double x)
{
  /* analytical fitting of synchrotron function F(x) */
  /* see http://arxiv.org/pdf/1301.6908.pdf */
    
  double F1 = M_PI*pow(2.0,5.0/3.0)/sqrt(3.0)/GAMMA13*pow(x,1.0/3.0);
  double F2 = sqrt(M_PI/2.0)*exp(-x)*pow(x,1.0/2.0);
    
  double a1_1 = -0.97947838884478688;
  double a1_2 = -0.83333239129525072;
  double a1_3 = 0.1554179602681624;
  double H_1 = a1_1*pow(x,1.0)+a1_2*pow(x,1.0/2.0)+a1_3*pow(x,1.0/3.0);
  double delta_1 = exp(H_1);
    
  double a2_1 = -0.0469247165562628882;
  double a2_2 = -0.70055018056462881;
  double a2_3 = 0.0103876297841949544;
  double H_2 = a2_1*pow(x,1.0)+a2_2*pow(x,1.0/2.0)+a2_3*pow(x,1.0/3.0);
  double delta_2 = 1.0-exp(H_2);
    
  return F1*delta_1+F2*delta_2;
}

void calc_syn_spec(double B, double r, double dr, double *gam, double *dgam, double *dN_dgam, double gam_max, int Nbin_e, double *gam_ph, double *P_nu_syn, double *alpha_nu_syn, int Nbin_ph)
{
  int i,j,k;
  double nu,x,sin_alpha=2./3.,tau_sa;
  double integ=0.,integ_alpha=0.,del_ln_gam=log(gam_max)/(double)(Nbin_e-1);
  double vol = 4.*M_PI*r*r*dr;

  for (k=0;k<Nbin_ph;k++) {
    integ = 0.0;
    integ_alpha = 0.0;
    nu = gam_ph[k]*MeC2/H;
    for (i=0;i<Nbin_e;i++) {
      x= (2.0*M_PI*nu)/(3.0*ELEC*gam[i]*gam[i]*B/2.0/M_ELE/C*sin_alpha); /* Eq. (6.17c) of Rybicki & Lightman */
      if (i==0 || i==Nbin_e-1) {
	integ += 0.5*dN_dgam[i]*dgam[i]*syn_func_fit(x);
	integ_alpha += -0.5*sin_alpha*pow(gam[i],2.0)*(-dN_dgam[i]/pow(gam[i],2.0))/dgam[i]*syn_func_fit(x)*gam[i]*del_ln_gam/MeC2;
      } else {
	integ += dN_dgam[i]*dgam[i]*syn_func_fit(x);
	integ_alpha += -sin_alpha*pow(gam[i],2.0)*(dN_dgam[i+1]/pow(gam[i+1],2.0)-dN_dgam[i]/pow(gam[i],2.0))/dgam[i]*syn_func_fit(x)*gam[i]*del_ln_gam/MeC2;
      }
    }
        
    P_nu_syn[k] = sqrt(3.0)*pow(ELEC,3.0)*B*sin_alpha/MeC2*integ; /* Eq. (6.33) x (2 pi) of Rybicki & Lightman */
    alpha_nu_syn[k] = C*C/8.0/M_PI/nu/nu*sqrt(3.0)*pow(ELEC,3.0)*B*integ_alpha/vol; /* Eq. (6.52) of Rybicki & Lightman */
    tau_sa = alpha_nu_syn[k]*dr;
        
    if (tau_sa > 1.0e-6){
      P_nu_syn[k] = (1.0-exp(-tau_sa))*P_nu_syn[k]/tau_sa;
    }
        
    integ = 0.0;
    integ_alpha = 0.0;
  }
}
