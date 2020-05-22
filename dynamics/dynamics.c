#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../hdr/astro_const.h"
#include "../hdr/ns_pwn_param.h"
#include "../hdr/dynamics.h"

int main()
{
  struct _input in = load_inputfile();
  calc_dynamics(in.Bp,in.Prot0,in.Mej,in.tmin,in.tmax,in.fac_dt,in.epsB,in.gamb);

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
  in.Bp = dmy[0];
  in.Prot0 = dmy[1];
  in.Mej = dmy[2]*Msun;
  in.tmin = dmy[3];
  in.tmax = dmy[4];
  in.fac_dt = dmy[5];
  in.epsB = dmy[6];
  in.gamb = dmy[8];

  return in;
}

void calc_dynamics(double Bp, double Prot0, double Mej, double tmin, double tmax, double fac_dt, double epsB, double gamb)
{
  FILE *op1,*op2;

  double tsd,Lsd,Erot,vsd,rsd,rnbsd,Bnbsd;
  sd(Bp,Prot0,Mej,&tsd,&Lsd,&Erot,&vsd,&rsd);

  op1 = fopen("../output/dynamics.dat","w");
  fprintf(op1,"#1. No. 2. t/tsd  3. t[s] 4. dt[s] 5. vej[cm/s] 6. rej[cm] 7. vnb[cm/s] 8. rnb[cm] 9. Bnb[G] 10. Lpsr[erg/s] 11. dr[cm] \n");
  int cnt=0;
  double t,dt,dt_dyn,dt_syn,rej_tmp,vej_tmp,rnb_tmp,vnb_tmp,Bnb_tmp,Lpsr_tmp;
  t = tmin;
  dt = tmin;
  while(t < tmax){
    vej_tmp = vej(t,tsd,vsd);
    rej_tmp = rej(t,tsd,rsd);
    vnb_tmp = vnb(t,tsd,vsd);
    rnb_tmp = rnb(t,tsd,rsd);
    Bnb_tmp = Bnb(t,tsd,Erot,rsd,epsB);
    Lpsr_tmp = Lpsr(t,tsd,Lsd);
    fprintf(op1,"%d %le %le %le %le %le %le %le %le %le %le\n",
	    cnt,t/tsd,t,dt,vej_tmp,rej_tmp,vnb_tmp,rnb_tmp,Bnb_tmp,Lpsr_tmp,dt*vnb_tmp);
    cnt++;

    /* Need to fix the time step issue */
    /* Probably I should not split the dynamics and microphysics in order to reduce the memory usage */
    dt_dyn = fac_dt*rej_tmp/vej_tmp;
    dt_syn = tsynb(Bnb_tmp,gamb);
    if (dt_syn <= dt_dyn){
      dt = dt_syn;
    } else {
      dt = dt_dyn;
    }
    t += dt;
    printf("dt= %12.3e [s]\n",dt_dyn);

  }
  fclose(op1);

  op2 = fopen("../output/dynamics.log","w");
  fprintf(op2,"< calculated >\n");
  fprintf(op2,"%d time step \n",cnt);
  fprintf(op2,"Erot = %12.3e [erg] \n",Erot);
  fprintf(op2,"tsd = %12.3e [s] \n",tsd);
  fprintf(op2,"Lsd = %12.3e [erg/s] \n",Lsd);
  fprintf(op2,"vsd = %12.3e [cm/s] \n",vsd);
  fprintf(op2,"rsd = %12.3e [cm]\n",rsd);
  fprintf(op2,"rnbsd = %12.3e [cm] \n",rnb(tsd,tsd,rsd));
  fprintf(op2,"Bnbsd = %12.3e [G]\n",Bnb(tsd,tsd,Erot,rsd,epsB));
  fprintf(op2,"< input >\n");
  fprintf(op2,"Bp = %12.3e [G] \n",Bp);
  fprintf(op2,"Prot0 = %12.3e [s] \n",Prot0);
  fprintf(op2,"Mej = %12.3e [Msun] \n",Mej/Msun);
  fprintf(op2,"tmin = %12.3e [s] \n",tmin);
  fprintf(op2,"tmax = %12.3e [s] \n",tmax);
  fprintf(op2,"dt/tdyn = %12.3e \n",fac_dt);
  fprintf(op2,"epsB = %12.3e \n",epsB);
  fclose(op2);

}

void sd(double Bp, double Prot0, double Mej, double *tsd, double *Lsd, double *Erot, double *vsd, double *rsd)
{
  double mu = .5*Bp*pow(Rns,3.);
  double ome = 2.*M_PI/Prot0;
  double lsd_tmp = pow(mu,2.)*pow(ome,4.)*pow(C,-3.)*(1.+fac_sd*sin(xi_mu)*sin(xi_mu));
  double erot_tmp = .5*Ins*pow(ome,2.);
  double tsd_tmp = .5*erot_tmp/lsd_tmp;
  double vsd_tmp = sqrt(erot_tmp/Mej);
  double rsd_tmp = 3./2.*vsd_tmp*tsd_tmp;

  *tsd = tsd_tmp;
  *Lsd = lsd_tmp;
  *Erot = .5*erot_tmp;
  *vsd = vsd_tmp;
  *rsd = rsd_tmp;
}

double Prot(double t, double tsd, double Prot0)
{
  double x = t/tsd;
  if (x<=1.){
    return Prot0;
  } else {
    return Prot0*pow(x,-0.5);
  }
}

double Lpsr(double t, double tsd, double lsd)
{
  double x = t/tsd;
  if (x<=1.){
      return lsd;
  } else {
      return lsd*pow(x,-2.);
  }
}

double Epsr(double t, double tsd, double Erot)
{
  double x = t/tsd;
  if (x<=1.){
    return 0.5*Erot*x;
  } else {
    return 0.5*Erot*(2.-1./x);
  }
}

double vej(double t, double tsd, double vsd)
{
  double x = t/tsd;
  if (x<=1.){
    return vsd*pow(x,.5);
  } else {
    return vsd*pow(2.-1./x,.5);
  }
}

double rej(double t, double tsd, double rsd)
{
  double x = t/tsd;
  if (x<1.){
    return rsd*pow(x,1.5);
  } else {
    return rsd + 2./3.*rsd*(sqrt(x*(2.*x-1)) -1. /*- log(1-2*(sqrt(4.-2./x)+2)*x/(-3.-pow(2.,1.5)))/pow(2.,1.5)*/);
    /* neglect the log term for similicity */
  }
}

double vnb(double t, double tsd, double vsd)
{
  return fac_nb*vej(t,tsd,vsd);
}

double rnb(double t, double tsd, double rsd)
{
  return fac_nb*rej(t,tsd,rsd);
}

double Bnb(double t, double tsd, double Erot, double rsd, double epsB)
{
  double E = Epsr(t,tsd,Erot);
  double R = rnb(t,tsd,rsd);
  return sqrt(6.*epsB*E/pow(R,3.));
}

double tsynb(double B, double gamb)
{
  double dgam_dt = 4.0/3.0*C*SIGMA_T*(B*B/8.0/M_PI)*gamb*gamb/MeC2;
  return gamb/dgam_dt;
}
