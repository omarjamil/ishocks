/* Original Written by Elme Breedt */

/* Modified to read parameter file read.par  */
/* and not read error for simulation data by Omar Jamil */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(){
  
  char inxray[40],inoptical[40],outfile[40];
  int i,j,k,M[10000],N1,N2,offset;
  double UDCF,DCF[100000],flux_av1,flux_av2,var1,var2,flux1[100000];
  double flux2[100000],time1[100000],time2[100000],error1[100000],error2[100000],tau;
  
  double flux1_2,flux2_2,time,flux,error,kfloat;
  double Deltatau,tau_min;
  FILE *in1,*in2,*out;

  int n;
  char params[40];
  char line[120];
  char par[80];
  

  /*open the parameter file*/
  FILE *fr = fopen ("read.par", "rt");
 
  /* get a line, up to 80 chars from fr.  done if NULL */
  while(fgets(line, 80, fr) != NULL)
    {
      /* convert the string to a long int */
      sscanf(line, "%s %s", par, params);
      if(par[0] == 'F')
        { 
          strcpy (inxray, params);
        }
      else if(par[0] == 'S')
        {
          strcpy (inoptical, params);
        }
      else if(par[0] == 'O')
        {
          strcpy (outfile, params);
        }
      else if(par[0] == 'T')
        {
          tau_min = atof(params);
        }
      else if(par[0] == 'D')
        {
          Deltatau = atof(params);
        }
            
    }
  /*close the parameter file*/
  fclose(fr);
  
  //printf("Enter: xraylc optical_lc output_file tau_min delta_tau")
  //scanf("%s %s %s %lf %lf",&inxray,&inoptical,&outfile,&tau_min,&Deltatau);
  
  in1=fopen(inxray,"r");
  i=0;
  flux_av1=0.0;
  flux1_2=0.0;
  N1=0;
  while(fscanf(in1,"%lf %lf",&time,&flux)!=EOF)
    { 
      time1[i]=time;
      flux1[i]=flux;
      /*error1[i]=error;*/
      //Small arbitrary errors for the simulation data
      //error1[i]=0.001 * flux;
      flux_av1+=flux1[i];
      flux1_2+=flux1[i]*flux1[i];
      N1+=1;
    i++;
    }
  fclose(in1);

  flux_av1=flux_av1/N1;
  var1=flux1_2/(N1-1.0)-flux_av1*flux_av1*N1/(N1-1.0);
  
  in2=fopen(inoptical,"r");
  j=0;
  flux_av2=0.0;
  flux2_2=0.0;
  N2=0;
  
  while(fscanf(in2,"%lf %lf",&time,&flux)!=EOF)
    {
      time2[j] = time;
      flux2[j] = flux;
      /*error2[j] = error;*/
      //Small arbitrary errors for the simulation data
      error2[i]=0.001 * flux;
      flux_av2+=flux2[j];
      flux2_2+=flux2[j]*flux2[j];
      N2+=1;
      j++;
    }
  
  flux_av2=flux_av2/N2;
  var2=flux2_2/(N2-1.0)-flux_av2*flux_av2*N2/(N2-1.0);
  fclose(in2);
  //  printf("flux average1=%lf  var1=%lf\n",flux_av2,var2); 
  
  offset=tau_min/Deltatau;
  
  for(k=0;k<=2*offset+1;k++)
    {
      M[k]=0;
      DCF[k]=0.0;
    }
  
  for(i=0;i<N1;i=i+1)
    {
      for(j=0;j<N2;j=j+1)
        {
          UDCF=(flux1[i]-flux_av1)*(flux2[j]-flux_av2)/sqrt(var1*var2);
          kfloat=(time2[j]-time1[i]+0.5*Deltatau)/Deltatau;
          k=kfloat+offset;
          //           printf("%i %i %i %i %lf %lf\n",k,offset,j,i,time2[j],time1[i]);
          if(k>=0&&k<=2*offset)
            {
              DCF[k]+=UDCF;
              M[k]++;
            }
        }
    }
  
  out=fopen(outfile,"w");
  for(k=0;k<=2*offset+1;k++)
    {
      if(M[k]!=0)
        {
          tau=(k-offset)*Deltatau;
          fprintf(out,"%lf %lf\n",tau,DCF[k]/M[k]);
        }
    }
  fclose(out);
  
}
