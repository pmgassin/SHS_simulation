/*********************************************************************************
******************PROGRAM sphere_SHS***************************************************
**********************************************************************************
   Py_SHS - An open source software about Second Harmonic Scattering 
   developed at Institut Charles Gerhardt Montpellier - ENSCM
   Dr Pierre-Marie GASSIN
   Dr Gassin GASSIN
    (june 2021)  

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

*******************************************************************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<string.h> 


// CONVENTION ZYZ pour matrice Euler!!!

double Txx(double phi,double theta,double psi){
double uu; 
uu=-sin(phi)*sin(psi)+cos(theta)*cos(phi)*cos(psi);
return uu;
} 

double Txy(double phi,double theta,double psi){
double uu; 
uu=-sin(phi)*cos(psi)-cos(theta)*cos(phi)*sin(psi);
return uu;
} 

double Txz(double phi,double theta,double psi){
double uu; 
uu=sin(theta)*cos(phi);
return uu;
} 

double Tyx(double phi,double theta,double psi){
double uu; 
uu=cos(phi)*sin(psi)+cos(theta)*sin(phi)*cos(psi);
return uu;
} 

double Tyy(double phi,double theta,double psi){
double uu; 
uu=cos(phi)*cos(psi)-cos(theta)*sin(psi)*sin(phi);
return uu;
} 

double Tyz(double phi,double theta,double psi){
double uu; 
uu=sin(theta)*sin(phi);
return uu;
} 

double Tzx(double phi,double theta,double psi){
double uu; 
uu=-sin(theta)*cos(psi);
return uu;
} 

double Tzy(double phi,double theta,double psi){
double uu; 
uu=sin(theta)*sin(psi);
return uu;
} 

double Tzz(double phi,double theta,double psi){
double uu; 
uu=cos(theta);
return uu;
} 

double T(double phi,double theta,double psi,int i,int j){
double uu;
if (i==0 && j==0){
uu=Txx(phi,theta,psi);
}
if (i==0 && j==1){
uu=Txy(phi,theta,psi);
}
if (i==0 && j==2){
uu=Txz(phi,theta,psi);
}
if (i==1 && j==0){
uu=Tyx(phi,theta,psi);
}
if (i==1 && j==1){
uu=Tyy(phi,theta,psi);
}
if (i==1 && j==2){
uu=Tyz(phi,theta,psi);
}
if (i==2 && j==0){
uu=Tzx(phi,theta,psi);
}
if (i==2 && j==1){
uu=Tzy(phi,theta,psi);
}
if (i==2 && j==2){
uu=Tzz(phi,theta,psi);
}
return uu;
}

double kimesoijk(double phi,double theta,double psi,int i,int j,int k,double beta[3][3][3]){
double uuu;
uuu=(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[0][0][0]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[0][0][1]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[0][0][2]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[0][1][0]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[0][1][1]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[0][1][2]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[0][2][0]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[0][2][1]);
uuu=uuu+(T(phi,theta,psi,i,0)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[0][2][2]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[1][0][0]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[1][0][1]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[1][0][2]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[1][1][0]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[1][1][1]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[1][1][2]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[1][2][0]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[1][2][1]);
uuu=uuu+(T(phi,theta,psi,i,1)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[1][2][2]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,0)*beta[2][0][0]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,1)*beta[2][0][1]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,0)*T(phi,theta,psi,k,2)*beta[2][0][2]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,0)*beta[2][1][0]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,1)*beta[2][1][1]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,1)*T(phi,theta,psi,k,2)*beta[2][1][2]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,0)*beta[2][2][0]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,1)*beta[2][2][1]);
uuu=uuu+(T(phi,theta,psi,i,2)*T(phi,theta,psi,j,2)*T(phi,theta,psi,k,2)*beta[2][2][2]);
    return uuu;
}


double xyzlabijk(double phi, double theta,double psi, int i, double x, double y,double z){
// c'est i qui dit si c'est x (i=0) y(i=1) ou z (i=2)...
double uu;
uu=T(phi,theta,psi,i,0)*x+T(phi,theta,psi,i,1)*y +T(phi,theta,psi,i,2)*z;
return uu;
}

double complex gxxx(double R,double thetac, double zc,double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}	
uuu=R*kimesoijk(phi,theta,psi,0,0,0,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}

double complex gxyy(double R,double thetac,double zc, double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}		
uuu=R*kimesoijk(phi,theta,psi,0,1,1,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}
	
double complex gxxy(double R,double thetac,double zc, double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}		
uuu=R*kimesoijk(phi,theta,psi,0,0,1,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}



double complex gzxx(double R,double thetac,double zc, double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}	
uuu=R*kimesoijk(phi,theta,psi,2,0,0,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}

double complex gzyy(double R,double thetac,double zc, double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}		
uuu=R*kimesoijk(phi,theta,psi,2,1,1,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}
	
double complex gzxy(double R,double thetac,double zc, double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}		
uuu=R*kimesoijk(phi,theta,psi,2,0,1,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}


double complex gyyy(double R,double thetac,double zc, double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}		
uuu=R*kimesoijk(phi,theta,psi,1,1,1,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}

double complex gyxx(double R,double thetac,double zc, double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}		
uuu=R*kimesoijk(phi,theta,psi,1,0,0,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}

double complex gyxy(double R,double thetac,double zc, double phi, double theta, double psi, double kisurf[3][3][3],double complex petitk,double complex grandk,double grandtheta){
double complex uuu;
double xmeso=R*cos(thetac);
double ymeso=R*sin(thetac);
double kimeso[3][3][3];
int aa,bb,cc;
for (aa=0;aa<3;aa++){
	for (bb=0;bb<3;bb++){
		for (cc=0;cc<3;cc++){
		kimeso[aa][bb][cc]=kimesoijk(thetac,1.570796,0,aa,bb,cc,kisurf);
		}
	}
}		
uuu=R*kimesoijk(phi,theta,psi,1,0,1,kimeso)*sin(theta)*cexp(((-2*petitk+grandk*cos(grandtheta))*(xyzlabijk(phi,theta,psi,2,xmeso,ymeso,zc))+(grandk*sin(grandtheta)*(xyzlabijk(phi,theta,psi,1,xmeso,ymeso,zc))))*I);
    return uuu;
}

int main(int argc, char *argv[])
{
FILE* fichier_beta=NULL;
FILE* fichier_out=NULL;
FILE* fichier_grace1=NULL;
char chaine[121];
char nome[20]; 

printf("********************************************************\n");
printf("****           Py_SHS: An open source software      ****\n");
printf("****         about Second Harmonic Scattering       ****\n");
printf("**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****\n");          
printf("****                  ICGM - ENSCM                  ****\n"); 
printf("****                   june 2021                    ****\n");
printf("********************************************************\n");
printf("                cylinder SHS programm                     \n");
printf("           Don't forget to cite this work:              \n");
printf("    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     \n");
printf("                                                        \n");
printf("********************************************************\n");
printf("********************************************************\n");
printf("********       Enter the input parameters     **********\n");
printf("********************************************************\n");

int i,j,k,kkk,ii,iii,i2,j2,k2;
int nom;
double complex grandk,petitk;
double pi=3.1415926535897;
double pisur2;
double beta[3][3][3];
double av=0,bv=0,cv=0,ah=0,bh=0,ch=0,bv1=0,bv2=0,bh1=0,bh2=0,Iy0,Iz0,Iy90,Iz90,avint,bvint,cvint,ahint,bhint,chint;

short nb_lignes_lues, nb_val_lues;
int pas,nx,ny,nx2,ny2;
double hx,hy,xi,yj,hx2,hy2,xi2,yj2,phii,thetai,psii;
double compteur;
double complex grandgammazxx_moy,grandgammayyy_moy,grandgammazyy_moy,grandgammaxxx_moy,grandgammayxx_moy,grandgammayxy_moy,grandgammaxyy_moy,grandgammaxxy_moy,grandgammazxy_moy;
double I4v,I2v,I2h,I4h;
double Iv,Ih,gammaaa,nomc,noms,deltax,deltay,deltaz;
char str1[20]="polarplot_90";
char str2[20]="polarplot_0";
char str3[20]="angle_scattering";
double grandtheta,grandthetadegree,IH_0,IH_90,IV_0,IV_90;
int pas2,pas3,mmm,mmmf;
double beta_y_moy_0,beta_z_moy_0,beta_y_moy_90,beta_z_moy_90,beta_xxx_moy,beta_xyy_moy;
double complex n;
double realpartn,imagpartn;
double longueuronde=800;
double R,length;
pisur2=pi/2;
pas2=10;
pas=10;




printf("Enter the wavelength in nm? ");
scanf("%lf", &longueuronde);
printf("The wavelength is: %lf \n",longueuronde);
printf("      ******************************     \n");
printf("Enter the real part of the refractive index? ");
scanf("%lf", &realpartn);
printf("      ******************************     \n");
printf("Enter the imaginary part of the refractive index? ");
scanf("%lf", &imagpartn);
printf("      ******************************     \n");
n=realpartn+imagpartn*I;
printf("The refractive index is: %lf + I*%lf \n",creal(n),cimag(n));
printf("      ******************************     \n");
printf("Enter the radius of the cylinder in nm? ");
scanf("%lf", &R);
printf("      ******************************     \n");
printf("Enter the length of the cylinder in nm? ");
scanf("%lf", &length);
printf("      ******************************     \n");

pas3=100;
if( argc == 4 ) {
      printf("The input ki(2) file is %s,  and the output file is %s \n", argv[2], argv[3]);
   }
   else if( argc > 4 ) {
      printf("Too many arguments supplied.\n");
   }
   else {
      printf("Computation failed: \n");
      printf("3 arguments expected : the keyword of the calculation, the Ki(2) inputfile and the outputfile. \n");
      printf("For more details see the documentation \n");
   }

fichier_beta=fopen(argv[2],"rt");
if (fichier_beta==NULL)
    {
        puts("Problem to open the input ki(2) file! \n");
        exit(0) ;
    }

fichier_out = fopen(argv[3], "w");
   if (fichier_out==NULL){
	printf("problem to open the output file to write results!\n");
}
printf("      ******************************     \n");
grandk=4*pi*n/longueuronde;
petitk=2*pi*n/longueuronde;

fgets (chaine, 121, fichier_beta);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][0][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][0][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][0][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][1][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][1][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][1][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][2][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][2][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[0][2][2]);

fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][0][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][0][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][0][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][1][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][1][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][1][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][2][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][2][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[1][2][2]);

fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][0][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][0][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][0][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][1][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][1][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][1][2]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][2][0]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][2][1]);
fgets (chaine, 121, fichier_beta);
nb_val_lues = sscanf (chaine, "%*s%lf",&beta[2][2][2]);

fclose(fichier_beta);


fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"****           Py_SHS: An open source software      ****\n");
fprintf(fichier_out,"****         about Second Harmonic Scattering       ****\n");
fprintf(fichier_out,"**** Dr Pierre-Marie GASSIN and Dr Gaelle GASSIN    ****\n");          
fprintf(fichier_out,"****                  ICGM - ENSCM                  ****\n"); 
fprintf(fichier_out,"****                   june 2021                    ****\n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"                                                        \n");
fprintf(fichier_out,"           Don't forget to cite this work:              \n");
fprintf(fichier_out,"    J. Chem. Inf. Model. 2020, 60, 12, pp 5912–5917     \n");
fprintf(fichier_out,"                                                        \n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"**       SIMULATION RESULTS  for cylinder programm  ****\n");
fprintf(fichier_out,"********************************************************\n");
fprintf(fichier_out,"*********          SURFACE KI(2) TENSOR      ***********\n");
fprintf(fichier_out,"\n");
fprintf(fichier_out,"    xxx=%lf   xxy=%lf   xxz= %lf            \n",beta[0][0][0],beta[0][0][1],beta[0][0][2]);
fprintf(fichier_out,"    xyx=%lf   xyy=%lf   xyz= %lf  \n",beta[0][1][0],beta[0][1][1],beta[0][1][2]);
fprintf(fichier_out,"    xzx=%lf   xzy=%lf   xzz= %lf  \n",beta[0][2][0],beta[0][2][1],beta[0][2][2]);
fprintf(fichier_out,"\n");
fprintf(fichier_out,"    yxx=%lf   yxy=%lf   yxz= %lf  \n",beta[1][0][0],beta[1][0][1],beta[1][0][2]);
fprintf(fichier_out,"    yyx=%lf   yyy=%lf   yyz= %lf  \n",beta[1][1][0],beta[1][1][1],beta[1][1][2]);
fprintf(fichier_out,"    yzx=%lf   yzy=%lf   yzz= %lf  \n",beta[1][2][0],beta[1][2][1],beta[1][2][2]);
fprintf(fichier_out,"\n");
fprintf(fichier_out,"    zxx=%lf   zxy=%lf   zxz= %lf  \n",beta[2][0][0],beta[2][0][1],beta[2][0][2]);
fprintf(fichier_out,"    zyx=%lf   zyy=%lf   zyz= %lf  \n",beta[2][1][0],beta[2][1][1],beta[2][1][2]);
fprintf(fichier_out,"    zzx=%lf   zzy=%lf   zzz= %lf  \n",beta[2][2][0],beta[2][2][1],beta[2][2][2]);
fprintf(fichier_out,"\n");
fprintf(fichier_out,"*******************************************************\n");
if (strcmp(str1,argv[1]) == 0){
printf("********** Begin calculation for cylinder polar plot at 90° ***********\n");
//printf("calculation in progress:");
grandtheta=pisur2;

av=0;
cv=0;
bv=0;
ah=0;
bh=0;
ch=0;

avint=0;
cvint=0;
bvint=0;
ahint=0;
bhint=0;
chint=0;


hx=(2*pi)/(pas);
hy=(length)/(pas3);
nx=pas;
ny=pas3;

hx2=(2*pi)/(pas);
hy2=(pi)/(pas);
nx2=pas;
ny2=pas;

grandgammazxx_moy=0+0*I;
grandgammayyy_moy=0+0*I;
grandgammazyy_moy=0+0*I;
grandgammaxxx_moy=0+0*I;
grandgammayxx_moy=0+0*I;
grandgammayxy_moy=0+0*I;
grandgammaxyy_moy=0+0*I;
grandgammaxxy_moy=0+0*I;
grandgammazxy_moy=0+0*I;

//double boucle sur phi et theta
for (i2=0;i2<nx2;i2++){
	compteur=i2*100/nx2;
	printf("calculation in progress : %lf pourcent \n",compteur);
	for (j2=0;j2<ny2;j2++){
//for (k2=0;k2<nx2;k2++){
phii=(0 + hx2/2 + i2*hx2);
        	thetai=(0 + hy2/2 + j2*hy2);
            	psii=0;//(0 + hx2/2 + k2*hx2);
grandgammazxx_moy=0+0*I;
grandgammayyy_moy=0+0*I;
grandgammazyy_moy=0+0*I;
grandgammaxxx_moy=0+0*I;
grandgammayxx_moy=0+0*I;
grandgammayxy_moy=0+0*I;
grandgammaxyy_moy=0+0*I;
grandgammaxxy_moy=0+0*I;
grandgammazxy_moy=0+0*I;

for (i=0;i<nx;i++){
	for (j=0;j<ny;j++){
        		
            	xi = (0 + hx/2 + i*hx);
            	yj =( (-length/2) + hy/2 + j*hy);
		grandgammazxx_moy=grandgammazxx_moy+hx*hy*gzxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxxx_moy=grandgammaxxx_moy+hx*hy*gxxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammayxx_moy=grandgammayxx_moy+hx*hy*gyxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxyy_moy=grandgammaxyy_moy+hx*hy*gxyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammazyy_moy=grandgammazyy_moy+hx*hy*gzyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammayyy_moy=grandgammayyy_moy+hx*hy*gyyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxxy_moy=grandgammaxxy_moy+hx*hy*gxxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammayxy_moy=grandgammayxy_moy+hx*hy*gyxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammazxy_moy=grandgammazxy_moy+hx*hy*gzxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
			}
           }
avint=avint+hx2*hy2*(grandgammaxxx_moy*conj(grandgammaxxx_moy));     
cvint=cvint+hx2*hy2*(grandgammaxyy_moy*conj(grandgammaxyy_moy));
bvint=bvint+hx2*hy2*(4*grandgammaxxy_moy*conj(grandgammaxxy_moy)+grandgammaxxx_moy*conj(grandgammaxyy_moy)+grandgammaxyy_moy*conj(grandgammaxxx_moy));

Iy0=(grandgammayxx_moy*cos(grandtheta)*cos(grandtheta)-grandgammazxx_moy*cos(grandtheta)*sin(grandtheta))*conj(grandgammayxx_moy*cos(grandtheta)*cos(grandtheta)-grandgammazxx_moy*cos(grandtheta)*sin(grandtheta));
Iz0=(-grandgammayxx_moy*cos(grandtheta)*sin(grandtheta)+grandgammazxx_moy*sin(grandtheta)*sin(grandtheta))*conj(-grandgammayxx_moy*cos(grandtheta)*sin(grandtheta)+grandgammazxx_moy*sin(grandtheta)*sin(grandtheta));
Iy90=(grandgammayyy_moy*cos(grandtheta)*cos(grandtheta)-grandgammazyy_moy*cos(grandtheta)*sin(grandtheta))*conj(grandgammayyy_moy*cos(grandtheta)*cos(grandtheta)-grandgammazyy_moy*cos(grandtheta)*sin(grandtheta));
Iz90=(-grandgammayyy_moy*cos(grandtheta)*sin(grandtheta)+grandgammazyy_moy*sin(grandtheta)*sin(grandtheta))*conj(-grandgammayyy_moy*cos(grandtheta)*sin(grandtheta)+grandgammazyy_moy*sin(grandtheta)*sin(grandtheta));

ahint=ahint+hx2*hy2*(csqrt(Iy0*Iy0+Iz0*Iz0));     
chint=chint+hx2*hy2*(csqrt(Iy90*Iy90+Iz90*Iz90));
bhint=bhint+hx2*hy2*(4*grandgammazxy_moy*conj(grandgammazxy_moy)+grandgammazxx_moy*conj(grandgammazyy_moy)+grandgammazyy_moy*conj(grandgammazxx_moy));
}
           }



//}
/*printf("grandgammaxyy %lf +I%lf \n",creal(grandgammaxyy_moy),cimag(grandgammaxyy_moy));
printf("grandgammayyy %lf +I%lf \n",creal(grandgammayyy_moy),cimag(grandgammayyy_moy));
printf("grandgammazxx %lf +I%lf \n",creal(grandgammazxx_moy),cimag(grandgammazxx_moy));
printf("grandgammaxxx %lf +I%lf \n",creal(grandgammaxxx_moy),cimag(grandgammaxxx_moy));
printf("grandgammayxx %lf +I%lf \n",creal(grandgammayxx_moy),cimag(grandgammayxx_moy));
printf("grandgammazyy %lf +I%lf \n",creal(grandgammazyy_moy),cimag(grandgammazyy_moy));
printf("grandgammayyy %lf +I%lf \n",creal(grandgammayyy_moy),cimag(grandgammayyy_moy));
printf("grandgammaxxy %lf +I%lf \n",creal(grandgammaxxy_moy),cimag(grandgammaxxy_moy));
printf("grandgammazxy %lf +I%lf \n",creal(grandgammazxy_moy),cimag(grandgammazxy_moy));

av=(grandgammaxxx_moy*conj(grandgammaxxx_moy));
cv=(grandgammaxyy_moy*conj(grandgammaxyy_moy));
bv=(4*grandgammaxxy_moy*conj(grandgammaxxy_moy)+grandgammaxxx_moy*conj(grandgammaxyy_moy)+grandgammaxyy_moy*conj(grandgammaxxx_moy));
bh=(4*grandgammazxy_moy*conj(grandgammazxy_moy)+grandgammazxx_moy*conj(grandgammazyy_moy)+grandgammazyy_moy*conj(grandgammazxx_moy));
Iy0=(grandgammayxx_moy*cos(grandtheta)*cos(grandtheta)-grandgammazxx_moy*cos(grandtheta)*sin(grandtheta))*conj(grandgammayxx_moy*cos(grandtheta)*cos(grandtheta)-grandgammazxx_moy*cos(grandtheta)*sin(grandtheta));
Iz0=(-grandgammayxx_moy*cos(grandtheta)*sin(grandtheta)+grandgammazxx_moy*sin(grandtheta)*sin(grandtheta))*conj(-grandgammayxx_moy*cos(grandtheta)*sin(grandtheta)+grandgammazxx_moy*sin(grandtheta)*sin(grandtheta));
ah=(csqrt(Iy0*Iy0+Iz0*Iz0));
Iy90=(grandgammayyy_moy*cos(grandtheta)*cos(grandtheta)-grandgammazyy_moy*cos(grandtheta)*sin(grandtheta))*conj(grandgammayyy_moy*cos(grandtheta)*cos(grandtheta)-grandgammazyy_moy*cos(grandtheta)*sin(grandtheta));
Iz90=(-grandgammayyy_moy*cos(grandtheta)*sin(grandtheta)+grandgammazyy_moy*sin(grandtheta)*sin(grandtheta))*conj(-grandgammayyy_moy*cos(grandtheta)*sin(grandtheta)+grandgammazyy_moy*sin(grandtheta)*sin(grandtheta));
ch=(csqrt(Iy90*Iy90+Iz90*Iz90));
*/

av=avint;
bv=bvint;
cv=cvint;
ah=ahint;
bh=bhint;
ch=chint;

printf("\n"); 
printf("**********calculation complete ************\n");           
printf("av = %lf\n",av);
printf("bv = %lf\n",bv);
printf("cv = %lf\n",cv);
printf("ah = %lf\n",ah);
printf("bh = %lf\n",bh);
printf("ch = %lf\n",ch);


printf("********************************************\n");   
I4v=(av-bv+cv)/((3*av+bv+3*cv));
I2v=(4)*(av-cv)/(3*av+bv+3*cv);
I4h=(ah-bh+ch)/((3*ah+bh+3*ch));
I2h=(4)*(ah-ch)/(3*ah+bh+3*ch);

printf("I2v = %lf \n",I2v);
printf("I4v = %lf \n",I4v);
printf("I2h = %lf \n",I2h);
printf("I4h = %lf \n",I4h);

fichier_grace1 = fopen("out_plot", "w");
   if (fichier_grace1==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire gnuplot!");
}
fprintf(fichier_out,"********************************************\n");
fprintf(fichier_out,"                                            \n");
fprintf(fichier_out,"* Cylinder POLAR PLOT - CONFIGURATION AT 90° *\n");
fprintf(fichier_out,"                                            \n");
fprintf(fichier_out,"********************************************\n");



fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"I2v=%lf  I4v=%lf \n",I2v,I4v);
fprintf(fichier_out,"I2h=%lf  I4h=%lf \n",I2h,I4h);
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"av=%lf   bv=%lf   bv1=%lf  bv2=%lf  cv= %lf \n",av,bv,bv1,bv2,cv);
fprintf(fichier_out,"ah=%lf   bh=%lf   bh1=%lf  bh2=%lf  ch= %lf \n",ah,bh,bh1,bh2,ch);

fclose(fichier_out);

fprintf(fichier_grace1,"gamma IV(gamma) IH(gamma)\n");
for (iii=0;iii<100;iii++){
gammaaa=2*pi*iii/100;
nomc=cos(gammaaa);
noms=sin(gammaaa);
Iv=av*(nomc*nomc*nomc*nomc)+bv*(nomc*nomc*noms*noms)+cv*(noms*noms*noms*noms);
Ih=ah*(nomc*nomc*nomc*nomc)+bh*(nomc*nomc*noms*noms)+ch*(noms*noms*noms*noms);
fprintf(fichier_grace1,"%lf %lf %lf\n",gammaaa,Iv,Ih);
}
printf("-----------------------------------------------------------------------------\n");
printf("calculation finished, to plot the results, open the outplot file with gnuplot\n");
printf("-----------------------------------------------------------------------------\n");
fclose(fichier_grace1);
}


else if (strcmp(str2,argv[1]) == 0){
printf("********** Begin calculation for sphere polar plot at 0° (transmission) ***********\n");
//printf("calculation in progress:");
grandtheta=0;
av=0;
cv=0;
bv=0;
ah=0;
bh=0;
ch=0;

hx=(2*pi)/(pas);
hy=(length)/(pas3);
nx=pas;
ny=pas3;

hx2=(2*pi)/(pas);
hy2=(pi)/(pas);
nx2=pas;
ny2=pas;

grandgammazxx_moy=0+0*I;
grandgammayyy_moy=0+0*I;
grandgammazyy_moy=0+0*I;
grandgammaxxx_moy=0+0*I;
grandgammayxx_moy=0+0*I;
grandgammayxy_moy=0+0*I;
grandgammaxyy_moy=0+0*I;
grandgammaxxy_moy=0+0*I;
grandgammazxy_moy=0+0*I;

//double boucle sur phi et theta
for (i2=0;i2<nx2;i2++){
	compteur=i2*100/nx2;
	printf("calculation in progress : %lf pourcent \n",compteur);
	for (j2=0;j2<ny2;j2++){
for (i=0;i<nx;i++){
	for (j=0;j<ny;j++){
        		phii=(0 + hx2/2 + i2*hx2);
        		thetai=(0 + hy2/2 + j2*hy2);
            	psii=0;
            	xi = (0 + hx/2 + i*hx);
            	yj =( (-length/2) + hy/2 + j*hy);
		    grandgammazxx_moy=grandgammazxx_moy+hx*hy*hy2*hx2*gzxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxxx_moy=grandgammaxxx_moy+hx*hy*hy2*hx2*gxxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammayxx_moy=grandgammayxx_moy+hx*hy*hy2*hx2*gyxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxyy_moy=grandgammaxyy_moy+hx*hy*hy2*hx2*gxyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammazyy_moy=grandgammazyy_moy+hx*hy*hy2*hx2*gzyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammayyy_moy=grandgammayyy_moy+hx*hy*hy2*hx2*gyyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxxy_moy=grandgammaxxy_moy+hx*hy*hy2*hx2*gxxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammayxy_moy=grandgammayxy_moy+hx*hy*hy2*hx2*gyxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammazxy_moy=grandgammazxy_moy+hx*hy*hy2*hx2*gzxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
			}
           }
           }
           }
printf("grandgammazxx %lf +I%lf \n",creal(grandgammazxx_moy),cimag(grandgammazxx_moy));
printf("grandgammaxxx %lf +I%lf \n",creal(grandgammaxxx_moy),cimag(grandgammaxxx_moy));
printf("grandgammayxx %lf +I%lf \n",creal(grandgammayxx_moy),cimag(grandgammayxx_moy));
printf("grandgammazyy %lf +I%lf \n",creal(grandgammazyy_moy),cimag(grandgammazyy_moy));
printf("grandgammayyy %lf +I%lf \n",creal(grandgammayyy_moy),cimag(grandgammayyy_moy));
printf("grandgammaxxy %lf +I%lf \n",creal(grandgammaxxy_moy),cimag(grandgammaxxy_moy));

av=(grandgammaxxx_moy*conj(grandgammaxxx_moy));
cv=(grandgammaxyy_moy*conj(grandgammaxyy_moy));
bv=(4*grandgammaxxy_moy*conj(grandgammaxxy_moy)+grandgammaxxx_moy*conj(grandgammaxyy_moy)+grandgammaxyy_moy*conj(grandgammaxxx_moy));
bh=(4*grandgammazxy_moy*conj(grandgammazxy_moy)+grandgammazxx_moy*conj(grandgammazyy_moy)+grandgammazyy_moy*conj(grandgammazxx_moy));
Iy0=(grandgammayxx_moy*cos(grandtheta)*cos(grandtheta)-grandgammazxx_moy*cos(grandtheta)*sin(grandtheta))*conj(grandgammayxx_moy*cos(grandtheta)*cos(grandtheta)-grandgammazxx_moy*cos(grandtheta)*sin(grandtheta));
Iz0=(-grandgammayxx_moy*cos(grandtheta)*sin(grandtheta)+grandgammazxx_moy*sin(grandtheta)*sin(grandtheta))*conj(-grandgammayxx_moy*cos(grandtheta)*sin(grandtheta)+grandgammazxx_moy*sin(grandtheta)*sin(grandtheta));
ah=(csqrt(Iy0*Iy0+Iz0*Iz0));
Iy90=(grandgammayyy_moy*cos(grandtheta)*cos(grandtheta)-grandgammazyy_moy*cos(grandtheta)*sin(grandtheta))*conj(grandgammayyy_moy*cos(grandtheta)*cos(grandtheta)-grandgammazyy_moy*cos(grandtheta)*sin(grandtheta));
Iz90=(-grandgammayyy_moy*cos(grandtheta)*sin(grandtheta)+grandgammazyy_moy*sin(grandtheta)*sin(grandtheta))*conj(-grandgammayyy_moy*cos(grandtheta)*sin(grandtheta)+grandgammazyy_moy*sin(grandtheta)*sin(grandtheta));
ch=(csqrt(Iy90*Iy90+Iz90*Iz90));

printf("\n"); 
printf("**********calculation complete ************\n");           
printf("av = %lf\n",av);
printf("bv = %lf\n",bv);
printf("cv = %lf\n",cv);
printf("ah = %lf\n",ah);
printf("bh = %lf\n",bh);
printf("ch = %lf\n",ch);


printf("********************************************\n");   
I4v=(av-bv+cv)/((3*av+bv+3*cv));
I2v=(4)*(av-cv)/(3*av+bv+3*cv);
I4h=(ah-bh+ch)/((3*ah+bh+3*ch));
I2h=(4)*(ah-ch)/(3*ah+bh+3*ch);

printf("I2v = %lf \n",I2v);
printf("I4v = %lf \n",I4v);
printf("I2h = %lf \n",I2h);
printf("I4h = %lf \n",I4h);

fichier_grace1 = fopen("out_plot", "w");
   if (fichier_grace1==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire gnuplot!");
}
fprintf(fichier_out,"********************************************\n");
fprintf(fichier_out,"                                            \n");
fprintf(fichier_out,"* Cylinder POLAR PLOT - CONFIGURATION AT 0° *\n");
fprintf(fichier_out,"                                            \n");
fprintf(fichier_out,"********************************************\n");



fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"I2v=%lf  I4v=%lf \n",I2v,I4v);
fprintf(fichier_out,"I2h=%lf  I4h=%lf \n",I2h,I4h);
fprintf(fichier_out,"*************************\n");
fprintf(fichier_out,"av=%lf   bv=%lf   bv1=%lf  bv2=%lf  cv= %lf \n",av,bv,bv1,bv2,cv);
fprintf(fichier_out,"ah=%lf   bh=%lf   bh1=%lf  bh2=%lf  ch= %lf \n",ah,bh,bh1,bh2,ch);

fclose(fichier_out);

fprintf(fichier_grace1,"gamma IV(gamma) IH(gamma)\n");
for (iii=0;iii<100;iii++){
gammaaa=2*pi*iii/100;
nomc=cos(gammaaa);
noms=sin(gammaaa);
Iv=av*(nomc*nomc*nomc*nomc)+bv*(nomc*nomc*noms*noms)+cv*(noms*noms*noms*noms);
Ih=ah*(nomc*nomc*nomc*nomc)+bh*(nomc*nomc*noms*noms)+ch*(noms*noms*noms*noms);
fprintf(fichier_grace1,"%lf %lf %lf\n",gammaaa,Iv,Ih);
}
printf("-----------------------------------------------------------------------------\n");
printf("calculation finished, to plot the results, open the outplot file with gnuplot\n");
printf("-----------------------------------------------------------------------------\n");
fclose(fichier_grace1);
}

else if (strcmp(str3,argv[1]) == 0){

fprintf(fichier_out,"***************************************************************\n");
fprintf(fichier_out,"                                            \n");
fprintf(fichier_out,"****            cylinder SHS angular distribution          ****\n");
fprintf(fichier_out,"                                            \n");
fprintf(fichier_out," ****************OUTPUT****************************************\n");
fprintf(fichier_out,  " grandtheta  IPSS-IH(0°)   IPPP-IH(90°)   ISSS-IV(0°)   ISPP-IV(90°) \n"); 

fichier_grace1 = fopen("out_plot", "w");
   if (fichier_grace1==NULL){
	printf("impossible d'ouvrir le fichier pour ecrire gnuplot!");
}
fprintf(fichier_grace1," grandtheta   IH(0°)  IH(90°)  IV(0°)  IV(90°) \n");
printf("Enter the number of points between 0 and Pi for the angular SHS analysis ?"); 
scanf("%d", &pas2);
mmmf=2*(pas2)+1;
printf("********** Begin Sphere SHS Angular Distribution calculation**************\n");
for (mmm=0;mmm<mmmf;mmm++){
	grandtheta=-pi+mmm*(pi/pas2);
    grandthetadegree=(grandtheta*360)/(2*pi);
    printf("angle value: %lf \n",grandthetadegree);

hx=(2*pi)/(pas);
hy=(length)/(pas3);
nx=pas;
ny=pas3;

hx2=(2*pi)/(pas);
hy2=(pi)/(pas);
nx2=pas;
ny2=pas;

grandgammazxx_moy=0+0*I;
grandgammayyy_moy=0+0*I;
grandgammazyy_moy=0+0*I;
grandgammaxxx_moy=0+0*I;
grandgammayxx_moy=0+0*I;
grandgammayxy_moy=0+0*I;
grandgammaxyy_moy=0+0*I;
grandgammaxxy_moy=0+0*I;
grandgammazxy_moy=0+0*I;

//double boucle sur phi et theta
for (i2=0;i2<nx2;i2++){
	compteur=i2*100/nx2;
	printf("calculation in progress : %lf pourcent \n",compteur);
	for (j2=0;j2<ny2;j2++){
for (i=0;i<nx;i++){
	for (j=0;j<ny;j++){
        		phii=(0 + hx2/2 + i2*hx2);
        		thetai=(0 + hy2/2 + j2*hy2);
            	psii=0;
            	xi = (0 + hx/2 + i*hx);
            	yj =( (-length/2) + hy/2 + j*hy);
		    grandgammazxx_moy=grandgammazxx_moy+hx*hy*hy2*hx2*gzxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxxx_moy=grandgammaxxx_moy+hx*hy*hy2*hx2*gxxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammayxx_moy=grandgammayxx_moy+hx*hy*hy2*hx2*gyxx(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxyy_moy=grandgammaxyy_moy+hx*hy*hy2*hx2*gxyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammazyy_moy=grandgammazyy_moy+hx*hy*hy2*hx2*gzyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammayyy_moy=grandgammayyy_moy+hx*hy*hy2*hx2*gyyy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammaxxy_moy=grandgammaxxy_moy+hx*hy*hy2*hx2*gxxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
       		grandgammayxy_moy=grandgammayxy_moy+hx*hy*hy2*hx2*gyxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
        	grandgammazxy_moy=grandgammazxy_moy+hx*hy*hy2*hx2*gzxy(R,xi,yj,phii,thetai,psii,beta,petitk,grandk,grandtheta);
			}
           }
           }
           }
	
IV_0=(grandgammaxxx_moy*conj(grandgammaxxx_moy));
IV_90=(grandgammaxyy_moy*conj(grandgammaxyy_moy));
//bv=creal(4*grandgammaxxy_moy*conj(grandgammaxxy_moy)+grandgammaxxx_moy*conj(grandgammaxyy_moy)+grandgammaxyy_moy*conj(grandgammaxxx_moy));
//bh=creal(4*grandgammazxy_moy*conj(grandgammazxy_moy)+grandgammazxx_moy*conj(grandgammazyy_moy)+grandgammazyy_moy*conj(grandgammazxx_moy));
Iy0=(grandgammayxx_moy*cos(grandtheta)*cos(grandtheta)-grandgammazxx_moy*cos(grandtheta)*sin(grandtheta))*conj(grandgammayxx_moy*cos(grandtheta)*cos(grandtheta)-grandgammazxx_moy*cos(grandtheta)*sin(grandtheta));
Iz0=(-grandgammayxx_moy*cos(grandtheta)*sin(grandtheta)+grandgammazxx_moy*sin(grandtheta)*sin(grandtheta))*conj(-grandgammayxx_moy*cos(grandtheta)*sin(grandtheta)+grandgammazxx_moy*sin(grandtheta)*sin(grandtheta));
IH_0=creal(csqrt(Iy0*Iy0+Iz0*Iz0));
Iy90=(grandgammayyy_moy*cos(grandtheta)*cos(grandtheta)-grandgammazyy_moy*cos(grandtheta)*sin(grandtheta))*conj(grandgammayyy_moy*cos(grandtheta)*cos(grandtheta)-grandgammazyy_moy*cos(grandtheta)*sin(grandtheta));
Iz90=(-grandgammayyy_moy*cos(grandtheta)*sin(grandtheta)+grandgammazyy_moy*sin(grandtheta)*sin(grandtheta))*conj(-grandgammayyy_moy*cos(grandtheta)*sin(grandtheta)+grandgammazyy_moy*sin(grandtheta)*sin(grandtheta));
IH_90=creal(csqrt(Iy90*Iy90+Iz90*Iz90));
       
printf("IH_0=%lf    IH_90=%lf   IV_0=%lf    IV_90=%lf \n",IH_0,IH_90,IV_0,IV_90);
fprintf(fichier_out," %lf    %lf   %lf   %lf   %lf \n",grandtheta,IH_0,IH_90,IV_0,IV_90);
fprintf(fichier_grace1," %lf    %lf   %lf   %lf   %lf \n",grandtheta,IH_0,IH_90,IV_0,IV_90); 
       
       
       }
printf("-----------------------------------------------------------------------------\n");       
printf("calculation finished, to plot the results, open the outplot file with gnuplot\n");
printf("-----------------------------------------------------------------------------\n");
fclose(fichier_out);
fclose(fichier_grace1);


}

else {
printf("--------------------- --------------------------------------------------\n"); 
printf("             ERROR in the name of the keyword                           \n"); 
printf("USE ONLY polarplot_90 OR polarplot_180 OR angle_scattering  !!  \n"); 
printf(" ! no calculation has been done ! Retry with the correct keyword        \n");
printf("------------------------------------------------------------------------\n"); 
}




}
