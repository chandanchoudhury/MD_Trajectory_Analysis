/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id$";

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<float.h>
#include "residence.h"

void distribution(int *, int, int); //function for creating the distributions

struct pdb
{
        char recordtype[7],atnam[7],resnam[5];
        int atno,resno;
        float x,y,z,r,s,t;
        float m,n;
};      
struct pdb a;
int main(int argc,char *argv[])
{ 
//	FILE **fz;
	FILE *fp,*fw,*fs,*fd,*flog;
	char line[121], **info = NULL;
	int llen, counter = 0;
	fp=fopen("53800.pdb","rb");  // SOURCE PDB FILE IS NEEDED. IT CAN BE THE OUTPUT PDB FILE OF THE SIMULATION.
	if(fp==NULL){
                	puts("Cannot Open Source File \n"); exit(1);
                }
	while (fgets(line,120,fp)) {

	// Allocate memory for pointer to line just added
	info = realloc(info,(counter+1) * sizeof(char *));
	// And allocate memory for that line itself!
	llen = strlen(line);
	info[counter] = calloc(sizeof(char),llen+1);
	// Copy the line just read into that memory
	strcpy(info[counter],line);	counter++;
							   }
	free((void*) info); fclose(fp);

	char filename[2000],st[100],f[80];
	int i = 0, t = 0, atno[counter], resno[counter], tot_no_atms = 0;
	int  *water, inc = 0, wat = 0, j = 0, wat_check = 0, count = 0;
	int temp = 0, wat_count = 0, repeat = 0, *wat_count_store;
	int store_size = 0, MAX_FRAME = 0, wat_count_prev = 0;
	float dist_check = 0.0, frames2ps = factor, boundary = ps/frames2ps;
	float *time, tim = 0, dist = 0.0, time_max = 0;
	float x[counter],y[counter],z[counter];
	char atnam[counter][7],resnam[counter][6];
    static char *desc[] = {
    "this is a small test program meant to serve as a template ",
    "when writing your own analysis tools. The advantage of ",
    "using gromacs for this is that you have access to all ",
    "information in the topology, and your program will be ",
    "able to handle all types of coordinates and trajectory ",
    "files supported by gromacs. Go ahead and try it! ",
    "This test version just writes the coordinates of an ",
    "arbitrary atom to standard out for each frame. You can ",
    "select which atom you want to examine with the -n argument."
  						  };
  static int n=1;

  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
  t_pargs pa[] = {
    { "-n", FALSE, etINT, {&n},
      "Plot data for atom number n (starting on 1)"
    }
  };
  
  t_topology top;
  int        ePBC;
  char       title[STRLEN];
  t_trxframe fr;
  rvec       *xtop;
  matrix     box;
  int        status;
  int        flags = TRX_READ_X;

  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  };
  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
   * calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  /* We don't need any topology information to write the coordinates,
   * but to show how it works we start by writing the name and
   * charge of the selected atom. It returns a boolean telling us
   * whether the topology was found and could be read
   */
  
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
  sfree(xtop);

  n=n-1; /* Our enumeration started on 1, but C starts from 0 */
  /* check that this atom exists */
  if(n<0 || n>(top.atoms.nr)) 
  {
    printf("Error: Atom number %d is out of range.\n",n);
    exit(1);
  }
	fp=fopen("53800.pdb","rb");  // SOURCE PDB FILE IS NEEDED. IT CAN BE THE OUTPUT PDB FILE OF THE SIMULATION.
	if(fp==NULL){
                	puts("Cannot Open Source File \n"); exit(1);
                }

	fw=fopen("water.dat","w");   	fs=fopen("water_count.xvg","w");
  	fd=fopen("water_occur.dat","w");
  	flog=fopen("residence.log","w");

//	fz = malloc(sizeof(FILE *) * 900);

    while(fgets(st,100,fp)!=NULL)
    {
    if(strncmp(st,"ATOM",4)==0)
        {
        strncpy(a.recordtype,st,6);
        for(t=0;t<=4;t++)
        f[t]=st[6+t];                        f[5]='\0';	a.atno=atoi(f);
        for(t=0;t<=4;t++)
        f[t]=st[12+t];                       f[6]='\0';	strcpy(a.atnam,f);
        for(t=0;t<=2;t++)
        f[t]=st[17+t];                       f[3]='\0';	strcpy(a.resnam,f);
//      chainnam[i]=st[21];
        for(t=0;t<=3;t++)
        f[t]=st[22+t];                       f[4]='\0';	a.resno=atoi(f);
        for(t=0;t<=7;t++)
        f[t]=st[30+t];                       f[8]='\0';	a.x=atof(f);
        for(t=0;t<=7;t++)
        f[t]=st[38+t];                       f[8]='\0';	a.y=atof(f);
        for(t=0;t<=7;t++)
        f[t]=st[46+t];                       f[8]='\0';	a.z=atof(f);

        atno[i]=a.atno;
        atnam[i][0]=a.atnam[0];atnam[i][1]=a.atnam[1];atnam[i][2]=a.atnam[2];atnam[i][3]=a.atnam[3];
        resnam[i][0]=a.resnam[0];resnam[i][1]=a.resnam[1];resnam[i][2]=a.resnam[2];
        resno[i]=a.resno;        x[i]=a.x;        i++;
        }
    }
    fclose (fp); tot_no_atms = i;

    printf("Tot_no_of_Atoms: %d\n",i);

//  printf("Atom name: %s\n",*(top.atoms.atomname[n]));//  printf("Atom charge:  %f\n",top.atoms.atom[n].q);
  /* The first time we read data is a little special */
    read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
  /* This is the main loop over frames */
  
  water = malloc(sizeof(int) * SIZE); memset(water,0,sizeof(int) * SIZE);
  time  = malloc(sizeof(float) * SIZE); memset(time,0,sizeof(float) * SIZE);
	
	fprintf(flog,"The range of solvation shell from the polymer backbone is %f nm \
	to %f nm\n",solvation_min,solvation_max);

  if (water == NULL || time == NULL)
  {
	  fprintf(stderr,"Could not allocate that much memory\n");
	  return 1;
  }
  counter = 0, temp = 0, MAX_FRAME = 0;

    do {
		inc = 0;
		for(i = 0; i < tot_no_atms; i++){
		if((strncmp(atnam[i]," N  ",4) == 0)|| (strncmp(atnam[i]," CA ",4) == 0)
		|| (strncmp(atnam[i]," CB ",4) == 0)){ // considering the backbone atoms
//		of the polymer
			for(t = 162; t < 39327; t++){
				if((strncmp(atnam[t]," OW ",4) == 0)){
				dist = sqrt(pow((fr.x[i][XX]-fr.x[t][XX]),2) +
				pow((fr.x[i][YY]-fr.x[t][YY]),2) +
				pow((fr.x[i][ZZ]-fr.x[t][ZZ]),2));

				if(dist > solvation_min && dist <= solvation_max){ // Checking the distance criteria, dist
//				value taken from the rdf plot of backbone atoms with OW

					wat_check = t; inc = 0; j = temp;

					for(j = temp; j < counter; j++){ 
						if(wat_check == water[j]){ // If a particular water
//                  (OW) is being repated in a frame, then it's discared, i.e;
//                  counting only once.
						inc++;
						}
					}
					if(inc == 0){
					water[counter] = t; time[counter] = fr.time;  // Unique
//					Waters of each frame are stored in the water array.
					time_max = fr.time;	counter++; 
					}
				}
			}
			}
		}
		else
			continue;
		}
		temp = counter;
		MAX_FRAME++;
   }while(read_next_frame(status,&fr));


	printf("\nMAX_FRAME  = %6d \n\n",MAX_FRAME);
	fprintf(flog,"\ntotal no. of frames = %6d \n",MAX_FRAME);
	fprintf(flog,"\nsoft boundary of %f ps which implies  %f frames \n",ps,boundary);

/*  counter represents the water array size */

	printf("SIZE %6d\n",counter); fprintf(flog,"\ntotal no. of loops %6d\n",counter);
	
	printf("\nWait ... ... Execution is in Progress ... ...\n\n");

	fprintf(flog,"\nProgress of loops ... ... ...\n\n");
 
 	wat_count_store = malloc(sizeof(int) * SIZE_2);
	memset(wat_count_store,0,sizeof(int) * SIZE_2);
	j = 0;
	
	for(i = 0; i < counter; i++)
	{
		fprintf(flog,"%8d\n",i);
		if(time[i] == time[i+1])
		{
			count++;	
		}

		else
		{
			fprintf(fs,"%8.3f\t%6d\n",time[i], count+1);			
			count = 0;
		}

		fprintf(fw,"%8.3f\t%6d\n",time[i],water[i]);

		tim = time[i]; wat = water[i]; wat_count = 0; repeat = 0;

		for(t = 0; t < i; t++){ // Checking for repeatation in the array
//          (water array).
        	if(wat == water[t])
            	repeat++;
        }
		
		if(repeat == 0)
		{

			wat_count_prev = 0; wat_count = 0; inc = 0;
		
			for(t = i; t < counter; t++)
			{
				if(time[t] >= tim )
				{
					if(time[t] == time[t+1])
					{
						if(water[t] == wat)
						{
//							wat_count_prev = wat_count; 
							wat_count++; 
//							inc = 0; 
						}
					}
	
					else
					{
						if(water[t] == wat)
						{
//							wat_count_prev = wat_count; 
							wat_count++;
//							inc = 0;
						}
	
						else if((wat_count != 0) && (wat_count_prev == wat_count))
						{
//							printf("%d %d \n",wat_count,inc);
							inc++; 
						}

						else if(inc < boundary)
						{
//							printf("%d %d \n",wat_count,inc);
//							printf("%d %d %d\n",water[i], wat_count,inc);
							wat_count += inc;	inc = 0;
						}
						
						if(inc == (int)boundary)
						{
							if(wat_count != 0){
							fprintf(fd,"%6d%8d\n",wat, wat_count);

							wat_count_store[j] = wat_count; j++;}
							wat_count_prev = 0; wat_count = 0; inc = 0;
						}
						
						if(time[t] == time_max)
						{
							if(wat_count != 0){
							fprintf(fd,"%6d%8d\n",wat, wat_count);

							wat_count_store[j] = wat_count; j++;
							}
							wat_count_prev = 0; wat_count = 0; inc = 0;
						}

						wat_count_prev = wat_count;
					}
				}
			}
		}
	}

	fprintf(flog,"\n\nLoops scanning completed.\nDistribution function is called for\n");

	distribution(&wat_count_store[0],j,MAX_FRAME); //calling the distribution function

	printf("\nExecution Completed.\n\n"); fprintf(flog,"\n\nExecution Completed.\n");

	free(water); free(time); free(wat_count_store);

	return 0;
}

void distribution(int *array, int size, int MAX_FRAME)
{
    FILE *fz, *flog;
  	flog=fopen("distri.log","w+");

	if(flog==NULL)
	{
		puts("Cannot Open Source File \n"); exit(1);
	}

    fz=fopen("distribution.xvg","w");

    int i = 0, count = 0, sum_count = 0, t = 0;
    float x = 0, y = 0, min = 0, max = 0, occur[500000], time[500000];
   //Calulate Min and max
    min = *array; max = *array;

    for(i = 0; i < size; i++)
    {
        if(min > array[i])
            min = array[i];
        if(max < array[i])
            max = array[i];
    }
	
	printf("%f %f \n\n",min,max); fprintf(flog,"\nminimum %f, maximum %f, factor\
	%f	\n",min,max, factor);

fprintf(fz,"@    xaxis  label \"stay of water molecules, soft boundary %5.2f ps, ps\"\n \
@ yaxis  label \"Probability \"\n@    target G0.S0 \n@   type \
xy\n@ altxaxis off\n@    altyaxis  off \n@    xaxis  tick place normal\n@   yaxis \
tick place normal\n@    frame type 1\n@    s0 hidden false\n@    s0 type xy\n@   s0 \
symbol 1\n@    s0 symbol size 0.400000\n@    s0    symbol color 1\n@    s0    symbol \
pattern 1\n@    s0 symbol fill color 1\n@    s0 symbol fill pattern 0\n@ s0   symbol \
linewidth 1.0\n@    s0 symbol  linestyle 1\n@    s0 symbol char 65\n@  s0 symbol   \
char font 0\n@    s0    symbol skip 0\n@    s0 line type 0\n@    s0 line  linestyle \
1\n@    s0 line  linewidth 1.0\n@    s0 line color 1\n",ps);

    for(x = min; x <= max; x += bin)
    {	
//		printf("\n%f %d",x, bin);
        for(i = 0; i < size; i++)
        {
            if((array[i] >= x) && (array[i] < (x + bin))){
            count ++;
			}
        }
//		printf("\ncount %d\n",count);
		sum_count += count;
		
        if(count != 0 && x !=0)
		{
			time[t] = x; occur[t] = (float)count; t++;
		}

        count = 0;
    }
	printf("Total Occurance %d\n ",sum_count);
	
	fprintf(flog,"\n\nTotal occurance %d \n",sum_count);

	for(i = 0; i < t; i++)
	{
//		fprintf(flog,"\n\n %d %f %f \n\n",i,time[i],occur[i]);
		fprintf(fz,"%f\t%E\n", time[i]*factor, occur[i]/sum_count);
	}

}
