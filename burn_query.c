/*
	Burn Query Constructor
	Because writing a new program for terribly similar tasks is a waste of time
	
	Takes a list of HDF5 files from the command prompt and then performs abundance queries
	based on the input you specify on run time. Still has a long way to go at this point. 
	~ Jack M. Sexton - 2014
	
	
	
	
	
	Planned Features
	
	* All of a species sorting - Just specify a proton count and catch all isotopes 
		above the fractional mass thresh-hold
		
	*Unified Query Entry - Enter nz nn frac-mass threshold in one line of input for 
		faster query construction
	
	*Improved Performance - Dog-slow when it works, lots of tweaks to help that though. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <SE.h>


struct linked_list
{
    int n;
	int z;
	double cutoff;
    struct linked_list *next ;
} ;

typedef struct linked_list NODE ;

NODE * getnode();
NODE * insert(NODE *head , int nn_sort, int nz_sort, double fmass_sort);
NODE * delete(NODE *head);
void display(NODE *head);
int get_nn();
int get_nz();
double get_fmass();


int main(int argc, char *argv[])
{

	int i, j, k;
	int sefp, nspecies, nread, printflag, nobj, *nn, *nz, *ids;
	char buff[50];
	char *pos;
	char tmp[50];
	FILE *filevar;
	NODE *head = NULL ;
	NODE *curr = NULL ;
	int input, nn_sort, nz_sort;
	double fmass_sort, *total_mass, mass, *frac_mass, mtot = 0.0; 
	
	if (argc < 2) 
    {
        fprintf(stderr, "Usage: %s HDF5_file [HDF5_file ...]\n", argv[0]);
        exit(-1);
    }
	
	while(1)
	{
		printf("\n 1. Insert isotope \n");
        printf("\n 2. Delete last isotope added \n");
		printf("\n 3. View current Query \n");
		printf("\n 4. Perform current Query on loaded files\n");
		printf("\n 5. Exit \n");
		
		printf("\n Please type the appropriate option: \n\n");
		
		scanf("%d" , &input);
		switch(input)
        {
			case 1:
				if (head == NULL)
				{
					printf("\n First Isotope is the primary sort parameter");
					printf("\n I.E. 26Al with a mass fraction above 10E_6 would be");
					printf("\n 13 13 6\n\n");
					
					nn_sort = get_nn();
					nz_sort = get_nz();
					fmass_sort = get_fmass();
					head = insert(head , nn_sort, nz_sort, fmass_sort);
					break;
            
				}
				else
				{
					printf("\n Secondary Isotope to return upon flagging of primary sort parameter \n");
					printf("\n I.E. In particles containing 26Al with 10E-6 fractional abundance");
					printf("\n (to stick with the example from the first particle)");
					printf("\n output 28Si above 10E-6 would be\n");
					printf("\n 14 14 6 \n");
					
					nn_sort = get_nn();
					nz_sort = get_nz();
					fmass_sort = get_fmass();
					head = insert(head, nn_sort, nz_sort, fmass_sort);
					break;
				}
				
			case 2:
				head = delete(head);
				break;
				
			case 3:
				display(head);
				break;
				
			case 4:
				fgets(tmp, 50, stdin);
				printf("\nAre you sure this is the query you wish to perform?\n");
				display(head);
				printf("\n(Y/y)es to continue (N/n)o to go back\n");
				fgets(buff, sizeof(buff), stdin);
				if ((pos=strchr(buff, '\n')) != NULL)
							*pos = '\0';
				if (buff == NULL) 
				{
					printf ("\nNo input\n");
					break;
				}
				else if (strcmp(buff,"No") == 0||strcmp(buff,"no") == 0||strcmp(buff,"N") == 0||strcmp(buff,"n") == 0) 
				{
					break;
				}
				else if (strcmp(buff,"Yes") == 0||strcmp(buff,"yes") == 0||strcmp(buff,"Y") == 0||strcmp(buff,"y") == 0)
				{
					printf("\nEnter output filename\n");
					fgets(buff, sizeof(buff), stdin);
					if ((pos=strchr(buff, '\n')) != NULL)
							*pos = '\0';
					if (buff == NULL) 
					{
						// Extra NL since my system doesn't output that on EOF.
						printf ("\nNo input\n");
						break;
					}
					else
					{
						sefp = SEopen(argv[1]);
						SEreadIArrayAttr(sefp, -1, "nn", &nn, &nspecies);
						SEreadIArrayAttr(sefp, -1, "nz", &nz, &nspecies);
						SEclose(sefp);
						printf("%s contains %d species\n", argv[1], nspecies);

						total_mass = (double *)calloc(nspecies, sizeof(double));

						filevar = fopen(buff,"a"); /* open the output file */

						for (i = 1; i < argc; ++i) /* for each i from 1 to argc */
						{
							sefp = SEopen(argv[i]); /* open the current file */
							printf("%s opened\n",argv[i]);

							/* For each file, determine the list of particles in the file. */
							nobj = SEncycles(sefp);
							ids = (int *)malloc(nobj * sizeof(int));
							printf("%u bytes Allocated for %i particles\n", nobj * sizeof(int), nobj);
							
							SEcycles(sefp, ids, nobj);
							
							for (j = 0; j < nobj; ++j)
							{
								mass = SEreadDAttr(sefp, ids[j], "mass");
								SEreadDArrayAttr(sefp, ids[j], "fmass", &frac_mass, &nread);
								
								if (nread != nspecies)
								{
									fprintf(stderr, "%s particle %d: nread (%d) != nspecies (%d)\n", 
										argv[i], ids[j], nread, nspecies);
									exit(-1);
								}
								
								printflag = 0;
								
								
								for (k = 0; k < nspecies; ++k)
								{						
									if ( (nz[k] == head->z) && (nn[k]== head->n) && (frac_mass[k] >= head->cutoff) )
									{
										printflag = 1;
										printf("particle %d flagged for %d:%d at %g \n", ids[j], nn[k], nz[k], frac_mass[k]);
										break;
									}
								}
								for (k = 0; k < nspecies; ++k)
								{ 
									total_mass[k] += mass * frac_mass[k];
									if(printflag == 1)
									{
										curr = head;
										while (curr != NULL)
										{
											if((nz[k] == curr->z) && (nn[k] == curr->n) && (frac_mass[k] >= curr->cutoff))
											{
												fprintf(filevar, "%d %d %e\n", nz[k], nn[k], frac_mass[k]);
												printf("%d:%d at %g saved to file\n", nn[k], nz[k], frac_mass[k]);
											}
											curr=curr->next;
											if(curr == NULL)
											break;
										}
									}
								}
							free(frac_mass);
							}
						free(ids);
						SEclose(sefp);
						}
						fclose(filevar);
						for (k = 0; k < nspecies; ++k) 
						{
							mtot += total_mass[k];
						}

						for (k = 0; k < nspecies; ++k)
						{
							printf("nn = %d\tnz = %d\tmass = %e (%.2f%%)\n", nn[k], nz[k],
								   total_mass[k], total_mass[k]/mtot * 100.0);
						}
						free(nn);
						free(nz);
						free(total_mass);
						while (head != NULL)
						{
							head=delete(head);
						}
						exit(1);
					}
				}
				else
				{
				printf("last input: |%s|", buff);
				exit(-1);
				}
				
			case 5:
				while (head != NULL)
				{
					head=delete(head);
				}
				exit(1);
		}
	}
}

NODE * getnode()
{
    NODE *create ;
    create = (NODE *)(malloc(sizeof(NODE)));
    create->next = NULL ;
    return create;
}

NODE *insert(NODE *head ,  int nn_sort, int nz_sort, double fmass_sort)
{
    NODE *makenode;
    NODE *prev = NULL, *curr = NULL ;
    curr = head ;
    if(head==NULL)
    {
        makenode = getnode();
        makenode->n = nn_sort;
		makenode->z = nz_sort;
		makenode->cutoff = fmass_sort;
		makenode->next  = NULL ;
        head = makenode ;
		printf("\nNN:NZ %d:%d added \n" , head->n, head->z ) ;
        return head ;
    }
    while(curr != NULL)
    {
        prev = curr ;
		curr = curr->next ;
    }
    makenode = getnode();
    makenode->n = nn_sort;
	makenode->z = nz_sort;
	makenode->cutoff = fmass_sort;
	makenode->next  = NULL ;
	prev->next = makenode ;
	printf("\nNN:NZ %d:%d added \n" , makenode->n, makenode->z ) ;
	return head;
}

NODE * delete(NODE *head)
{
    if (head == NULL)
    {
        printf("\n Deleting Not Possible, List Empty \n");
    }
    else if (head->next == NULL )
    {
        printf("\nNN:NZ %d:%d removed \n" , head->n, head->z ) ;
		free(head);
        return NULL;
    }
    else
    {
        NODE *prev = NULL, *curr = head ;
        while(curr->next != NULL)
        {
            prev = curr ;
			curr = curr->next ;
        }
		
        prev->next = NULL;
		printf("\nNN:NZ %d:%d removed \n" , curr->n, curr->z ) ;
        free(curr);
    }
    return head;
}

void display(NODE *head)
{
    NODE *q;
    q = head;
    if(q == NULL)
    {
        printf("\n  List Empty \n");
        return;
    }
    while(q != NULL)
    {
        if(q->next == NULL)
        {
            printf("\nNN:NZ %d:%d above %g \n" , q->n, q->z, q->cutoff );
            break;
        }
        else
        {
            printf("\nNN:NZ %d:%d above %g \n" , q->n, q->z, q->cutoff );
            q = q->next ;
        }
    }
}

int get_nn()
{
	int nn_get;
	printf("\nNeutron Count?\n\n");
	scanf("%d" , &nn_get);
	return(nn_get);
}

int get_nz()
{
	int nz_get;
	printf("\nProton Count?\n\n");
	scanf("%d" , &nz_get);
	return(nz_get);
}

double get_fmass()
{
	int fget ;
	double ftemp = 0 ;
	printf("\nFractional Mass Threshold? \n");
	printf("\nin the form of: cut-off values below 10 to the -n");
	printf("\nwhere n is what you enter now\n\n");
	scanf("%d" , &fget);
	switch(fget)
	{
		case 1: 
			ftemp = 10E-1;
			break;
		case 2:
			ftemp = 10E-2;
			break;
		case 3:
			ftemp = 10E-3;
			break;
		case 4:
			ftemp = 10E-4;
			break;
		case 5:
			ftemp = 10E-5;
			break;
		case 6:
			ftemp = 10E-6;
			break;
		case 7:
			ftemp = 10E-7;
			break;
		case 8:
			ftemp = 10E-8;
			break;
		case 9:
			ftemp = 10E-9;
			break;
		case 10:
			ftemp = 10E-10;
			break;
	}
	return(ftemp);
}
