// Accepts xyz files
// Open a new file, and creates an input file for CINDO 
// Energy calculation, INDO1/CI semi-empirical
// Singles only, window of 10.
// Correct charge adding
// Asks for charge + Multiplicity
// Atoms = C,N,H,O,F,CL,Si,P, Zn, Ru, Fe, S,Br

#include <stdio.h>
#include <fstream.h>           
#include <stdlib.h> 
#include <iomanip.h>
#include <iostream.h>
#include <math.h>

main(int argc , char *argv[] )
{
  char p[100],namein[50],nameout[50], name[5];
  double X,Y,Z;
  int Atype=0;
  int charge=0;
  int counter=0;
  int multi=0;

  if ( argc != 5 )
      cout<<"Usage: jaguarout_to_CNDO filein fileout charge multiplicity"<<endl;

  sprintf (nameout, "%s" , argv[2] ) ;
  ofstream fpout (nameout);

  sprintf (namein, "%s" , argv[1] ) ;
  ifstream fpin (namein);

  charge =atoi(argv[3]);
  multi=atoi (argv[4]);  

  fpin>>counter;
  fpin.getline(p,80);
  fpin.getline(p,80);

  cout<<"The charge is "<<charge<<" , multiplicity is "<<multi<<endl;
  fpout <<namein<<endl;
  fpout <<"HAMILT=INDO"<<endl;
  fpout <<"STOP=CI"<<endl;
  fpout <<"ROTINV=YES"<<endl;
  fpout <<"BETA=INDO/S"<<endl;
  fpout <<"XCI=100 100"<<endl;
  fpout <<"EX_FROM=300"<<endl;
  fpout <<"MAX_CI=500"<<endl;
  fpout <<"POINTGRP=C1"<<endl;
  fpout <<"CHARGE="<<charge<<endl;
  fpout <<"MULT_CI="<<multi<<endl;
  fpout <<"RESTART=AUTO"<<endl;
  fpout <<"SHIFT=20.0"<<endl;
  fpout <<"MAX_ITS=1000"<<endl; 
  fpout <<endl;
  fpout.close();

  FILE *fp;
  if ((fp=fopen(nameout,"a"))==NULL)
    exit(1);

  for (int i=0;i<counter;i++)
    {
      fpin >>name>>X>>Y>>Z;

      if (name[0]=='H')
	Atype=1;
      if (name[0]=='F')
	  if (name[1]=='e')
            Atype=26;
          else
	    Atype=9;
      if (name[0]=='C')
	  if (name[1]=='l')
	    Atype=17;
	  else
	    Atype=6;
      if (name[0]=='O')
         if (name[1]=='s')
           Atype=76;
         else
	   Atype=8;
      if (name[0]=='B')
	Atype=35;
      if (name[0]=='N')
	Atype=7;
      if (name[0]=='P')
	Atype=15;
      if (name[0]=='Z')
	Atype=30;
      if (name[0]=='R')
	Atype=44;
      if (name[0]=='S')
	  if (name[1]=='i')
	    Atype=14;
	  else
	    Atype=16;
      
      fprintf(fp,"%10.6f%10.6f%10.6f%5d\n",X,Y,Z,Atype);
    }
}
	      
