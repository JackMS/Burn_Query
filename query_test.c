#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

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
	NODE *head = NULL ;
	int input, nn_sort, nz_sort ;
	double fmass_sort ;
	
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
					printf("\nFirst Isotope is the primary sort parameter");
					printf("\nI.E. 26Al with a mass fraction above 10E_6 would be");
					printf("\n13 13 6\n");
					
					nn_sort = get_nn();
					nz_sort = get_nz();
					fmass_sort = get_fmass();
					head = insert(head , nn_sort, nz_sort, fmass_sort);
					break;
            
				}
				else
				{
					printf("\nSecondary Isotope to return upon flagging of primary sort parameter \n");
					printf("\nI.E. In particles containing 26Al with 10E-6 fractional abundance");
					printf("\n(to stick with the example from the first particle)");
					printf("\noutput 28Si above 10E-6 would be");
					printf("\n14 14 6 ");
					
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
				exit(1);
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
