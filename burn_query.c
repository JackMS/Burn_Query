/*
	Burn Query Constructor
	Because writing a new program for terribly similar tasks is a waste of time
	
	Takes a list of HDF5 files from the command prompt and then performs abundance queries
	based on the input you specify on run time. Still has a long way to go at this point. 
	~ Jack M. Sexton - 2014
	
	1/29/2015  
	COMMENTS!
	~ JMS
	
	3/18/2015
	Now polls a SDF file as well and puts the abundance data into it's spatial context.
	Has the same restrictions as SDF_reader though, as this functionality is inherited.
	More machine-friendly output styling.
	~ JMS
	
	Planned Features
	
	* All of a species sorting - Just specify a proton count and catch all isotopes 
		above the fractional mass thresh-hold
		
	*Unified Query Entry - Enter nz nn frac-mass threshold in one line of input for 
		faster query construction
	
	*Re-write SDF functionality using Tree16's library
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <SE.h>

//The particle structure from SNSPH, defaulting to jet3b et al's styling
//WIL NOT WORK WITH CCO_2, without major hand holding for each file anyway
typedef struct
{
    double x, y, z;		/* position of body */
    float mass;			/* mass of body */
    float vx, vy, vz;		/* velocity of body */
    float u;			/* internal energy */
    float h;			/* smoothing length */
    float rho;			/* density */
    float drho_dt;              /* time derivative of rho */
    float udot;			/* time derivative of u */
    float ax, ay, az;		/* acceleration */
    float lax, lay, laz;	/* acceleration at tpos-dt */
    float phi;			/* potential */
    float idt;			/* timestep */
    unsigned int nbrs;          /* number of neighbors */
    unsigned int ident;		/* unique identifier */
    unsigned int windid;        /* wind id */
    float temp;                 /* temperature */
    float Y_el;                  /* for alignment */
    float f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22; 
    int p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22; 
    int m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22; 
} particle; 

// Linked List Node Object
struct linked_list
{
    int n;
	int z;
	double cutoff;
    struct linked_list *next ;
} ;
// C doesn't understand linked lists so you have to define it in terms it does
typedef struct linked_list NODE ;

// Declaring the linked list functions and container variables
NODE * getnode();
NODE * insert(NODE *head , int nn_sort, int nz_sort, double fmass_sort);
NODE * delete(NODE *head);
void display(NODE *head);
int get_nn();
int get_nz();
double get_fmass();

// The Main Loop - The Actual Program
int main(int argc, char *argv[])
{

// Declarations
	int i, j, k, h, sz, nobj_sdf;
	int offset=1600;
	int sefp, nspecies, nread, printflag, nobj, *nn, *nz, *ids;
	char buff[50];
	char *pos;
	char tmp[50];
	particle part;
	FILE *filevar;
	NODE *head = NULL ;
	NODE *curr = NULL ;
	int input, nn_sort, nz_sort;
	double fmass_sort, *total_mass, mass, *frac_mass, mtot = 0.0; 
	
	// Naked command test
	if (argc < 2) 
    {
        fprintf(stderr, "Usage: %s SDF_file HDF5_file [HDF5_file ...]\n", argv[0]);
        exit(-1);
    }
	
	// User Interface 
	while(1)
	{
		printf("\n 1. Insert isotope \n");
        printf("\n 2. Delete last isotope added \n");
		printf("\n 3. View current Query \n");
		printf("\n 4. Perform current Query on loaded files\n");
		printf("\n 5. Exit \n");
		
		printf("\n Please type the appropriate option: \n\n");
		
		// Interface with the User 
		scanf("%d" , &input);
		switch(input)
        {
			case 1:
				if (head == NULL)
				{
					printf("\n First Isotope is the primary sort parameter");
					printf("\n I.E. 26Al with a mass fraction above 10E_6 would be");
					printf("\n 13 13 6\n\n");
					
					// Fix this Ugly Shit later 
					nn_sort = get_nn();
					nz_sort = get_nz();
					fmass_sort = get_fmass();
					// Insert the node at the top of the list 
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
					// Insert the node at the next position of the list 
					head = insert(head, nn_sort, nz_sort, fmass_sort);
					break;
				}
				
			case 2:
				// Deletes the last node, not necessarily the head but including the head 
				head = delete(head);
				break;
				
			case 3:
				// Read the List Node by Node and display*/
				display(head);
				break;
				
			case 4:
				// Purge the buffer to remove trailing newline from using scanf 
				fgets(tmp, 50, stdin);
				// Warn the User before embarking 
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
					// Get the filename
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
						//Starting the Query
						
						//Read the SDF file
						FILE *fp = fopen(argv[1], "rb");
						
						//Read the first HDF5 file
						sefp = SEopen(argv[2]);
						
						//Get the species
						SEreadIArrayAttr(sefp, -1, "nn", &nn, &nspecies);
						SEreadIArrayAttr(sefp, -1, "nz", &nz, &nspecies);
						
						//Close it for now to keep things tidy
						SEclose(sefp);
						printf("%s contains %d species\n", argv[2], nspecies);

						//Build a buffer for summing the masses
						total_mass = (double *)calloc(nspecies, sizeof(double));

						//Open the output file
						filevar = fopen(buff,"a"); 
						
						//Inspect the SDF for details
						//Will be replaced later with proper LibSDF calls
						fseek(fp, 0L, SEEK_END);
						sz = ftell(fp);
						nobj_sdf = (sz-offset)/sizeof(particle);

						//Write the header of the output file
						fprintf(filevar, "id, x, y, z, rho");
						curr = head;
						while (curr != NULL)
						{
							fprint(filevar, ", %d:%d", curr->z, curr->n);
							curr=curr->next;
							if(curr == NULL) 
							{	
								fprint(filevar, "\n");
								break;
							}
						}
						
						//Reset the file-pointer for the SDF file
						fseek(fp,offset,SEEK_SET);
						
						//Iterate through the HDF5 files
						for (i = 2; i < argc; ++i) 
						{
							//Open the current file
							sefp = SEopen(argv[i]); 
							printf("%s opened\n",argv[i]);

							// For each HDF5 file, determine the list of particles in the file
							nobj = SEncycles(sefp);
							ids = (int *)malloc(nobj * sizeof(int));
							printf("%u bytes Allocated for %i particles\n", nobj * sizeof(int), nobj);
							
							SEcycles(sefp, ids, nobj);
							
							for (j = 0; j < nobj; ++j)
							{
								//Read the mass attributes for each
								mass = SEreadDAttr(sefp, ids[j], "mass");
								SEreadDArrayAttr(sefp, ids[j], "fmass", &frac_mass, &nread);
								
								//Look for mismatch
								if (nread != nspecies)
								{
									fprintf(stderr, "%s particle %d: nread (%d) != nspecies (%d)\n", 
										argv[i], ids[j], nread, nspecies);
									exit(-1);
								}
								
								//Null out the flag for each particle
								printflag = 0;
								
								for (k = 0; k < nspecies; ++k)
								{						
									//Test against the head of the linked list for the isotope we want
									if ( (nz[k] == head->z) && (nn[k]== head->n) && (frac_mass[k] >= head->cutoff) )
									{
										//Flag it for later
										printflag = 1;
										printf("particle %d flagged for %d:%d at %g \n", ids[j], nn[k], nz[k], frac_mass[k]);
										break;
									}
								}
								for (k = 0; k < nspecies; ++k)
								{ 
									//Count this isotope towards the total mass
									total_mass[k] += mass * frac_mass[k];
									//Check that print flag from earlier
									if(printflag == 1)
									{
										
										
										fprintf(filevar, "%d, %g, %g, %g, %g", part.ident, part.x, part.y, part.z, part.rho);
										
										
										
										//Go to the top of the list
										curr = head;
										//For the whole list
										while (curr != NULL)
										{
											//Check sort parameters
											if((nz[k] == curr->z) && (nn[k] == curr->n) && (frac_mass[k] >= curr->cutoff))
											{
												//Output to the file
												fprintf(filevar, "%d %d %e\n", nz[k], nn[k], frac_mass[k]);
												printf("%d:%d at %g saved to file\n", nn[k], nz[k], frac_mass[k]);
											}
											//Next item on list
											curr=curr->next;
											//Unless it's not
											if(curr == NULL)
											break;
										}
									}
								}
							//Tidy Up
							free(frac_mass);
							}
						//Tidy Up
						free(ids);
						SEclose(sefp);
						}
						fclose(filevar);
						//Update final totals
						for (k = 0; k < nspecies; ++k) 
						{
							mtot += total_mass[k];
						}
						//Console Output of Final results
						for (k = 0; k < nspecies; ++k)
						{
							printf("nn = %d\tnz = %d\tmass = %e (%.2f%%)\n", nn[k], nz[k],
								   total_mass[k], total_mass[k]/mtot * 100.0);
						}
						//Tidy Up
						free(nn);
						free(nz);
						free(total_mass);
						//Don't forget the linked list
						while (head != NULL)
						{
							head=delete(head);
						}
						exit(1);
					}
				}
				else
				{
				//It's dead, spit out the last input so you can diagnose why
				printf("last input: |%s|", buff);
				exit(-1);
				}
				
			case 5:
			`//For the whole list
				while (head != NULL)
				{
					//Remove
					head=delete(head);
				}
				exit(1);
		}
	}
}
//Make a Node
NODE * getnode()
{
    NODE *create ;
    create = (NODE *)(malloc(sizeof(NODE)));
    create->next = NULL ;
    return create;
}

//Add a node to the list
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

//Remove a Node
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

//Display the List
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

//Get User Input
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
