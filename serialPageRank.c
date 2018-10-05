/* Parallel and Distributed Systems 4th assignment
 * 
 * Serial PageRank Algorith with Gauss Seidel
 * 
 * Author: Argiris Papoudakis, papoudaa@auth.gr
 * 
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>


// error computation for the convergence of the algorithm
double errorComputation(double *array1, double *array2, int N)
{
    
    double max = -500;
    int i;

    for (i = 0; i < N; ++i)
    {
    double absDiff = fabs(array1[i] - array2[i]);

    if (max < absDiff)
      max = absDiff;
    }

    //return the maximum difference
    return max;
}


void pageRankGaussSeidel(int N,int **Graph,int *outBoundLinks,int *inBoundLinks,double *pageRankVector,double *E,int *iterations)
{
    int i,j;
    double e=1,sum1=0,sum2=0;
    
    //damping factor
    double d = 0.85;
    
    //save the pagerank vector of previous iteration in z
    double *previous;
    
    previous = malloc(sizeof(double) * N);
    if (previous == NULL)
    {
        printf("Could not allocate memory for the z array");
        exit(-1);
    }
    
    // initaliaze the pagerank vector
    for (i=0;i<N;i++)
    {
        pageRankVector[i] = 1.0 / N;
    }

    (*iterations) = 0;
    
    while (e > 1e-6)
    {
        
        sum1=0;
        
        //dangling nodes
        for(j = 0; j < N; j++)
        {
            
            if(outBoundLinks[j] == 0)
            {
                sum1 = sum1 + pageRankVector[j] / N;
            }
        
        }        
        
        for(i = 0; i < N; i++)
        {
            
            // Initialize the necessary vectors for the current iteration.
            previous[i] = pageRankVector[i];
            sum2 = 0;
            
            
            for (j=0;j<inBoundLinks[i];j++)
            {    
                
                sum2 = sum2  + pageRankVector[Graph[i][j]] / outBoundLinks[Graph[i][j]];
                
            }
            
            //update pageRank vector
            pageRankVector[i] = d * (sum1+sum2) + (1-d)* E[i];
            
        }
        
        //increase the counter
        (*iterations)++;
        
        //compute error
        e = errorComputation(pageRankVector, previous, N);
    }

    
}


//test the results
void test(int N,int *pageRankVector,char *filename)
{
    
    int i,n,counter=0;
    FILE *file = fopen(filename, "r");
    
    double *test;
    
    test = malloc(sizeof(double)*N);
    
    for (i=0;i<N;i++)
    {
        n = fscanf(file,"%lf",&test[i]);
    }
    
    for (i=0;i<N;i++)
    {
        if (fabs(pageRankVector[i]-(double)test[i])/fabs((double)test[i]) > 0.1)
        {
            counter++;
        }
    }
    
    if((double)(N-counter) / N * 100 < 95.0)
		printf("Test Failed!\n");
	else
		printf("Test Passed!\n");
        
}

int main(int argc , char **argv)
{
    
    int i=0,j=0,k=0,iter,**Graph; 
    int NumberOfNodes,numEdges,N;
    
    struct timeval startwtime, endwtime;
    double time;
    
    char *filename,line[1000],s[2] = " ", *token;
    
    
    if (argc != 2) 
    {
        printf("Usage: %s <filename> \n", argv[0]);
        exit(1);
    }
    
    filename = argv[1];

    //open the txt file containing the graph
    FILE *file = fopen(filename, "r");
    
    if(file == NULL)
    {
        printf("Error opening file \n");   
        exit(1);             
    }
    
   
    //find the nodes and the edges of the graph
    while (fgets(line,1000,file) != NULL)
    {
       
        //check if the line start with #
        if (line[0] == '#')
        {
            //find the line that have nodes and edges
            token = strtok(line,s);
            token = strtok(NULL,s);
             
            if ( !strcmp(token,"Nodes:") )
            {
                token = strtok(NULL,s);
                
                //convert string to integer
                NumberOfNodes = atoi(token);
                
                token = strtok(NULL,s);
                token = strtok(NULL,s);
                
                //convert string to integer
                numEdges = atoi(token); 
            }
          
        } 
       
    }
    
    printf("Number of nodes: %d \n",NumberOfNodes);
    printf("Number of edges: %d \n",numEdges);
    
    N = NumberOfNodes;
    
    //resetting pointer to the start of file
    fseek(file, 0, SEEK_SET);
   
    //check the case that the file have nodes with ids bigger than the number of nodes (web-Google.txt)
    while (fgets(line,1000,file) != NULL)
    {
        token = strtok(line,s);
       
        //check if the line dont start with #
        if (strcmp(token,"#"))
        {
            sscanf(line, "%d %d \n", &i, &j);
           
            if (i>N)
            {
                N = i;
            }
            
            if (j>N)
            {
                N = j;
            }
           
           
        } 
       
    }

    // inbound links for the nodes
    Graph = malloc(sizeof(int*) * N);
    if (Graph==NULL)
    {
        printf("Could not allocate memory for the Graph array\n");
        exit(-1);
    }
    
    // number of the outbound links for every nodes 
    int *outBoundLinks;
    
    outBoundLinks = malloc(sizeof(int)*N);
    if (outBoundLinks==NULL)
    {
        printf("Could not allocate memory for the outBoundLinks array \n");
        exit(-1);
    
    }
    

    //resetting pointer to the start of file
    fseek(file, 0, SEEK_SET);

    while (fgets(line,1000,file) != NULL)
    {
        
        token = strtok(line,s);
    
              
        // check that the line dont start with #
        if (strcmp(token,"#"))
        {
            
            sscanf(line, "%d %d \n", &i, &j);
            
            outBoundLinks[i]++;
            
        }
        
       
    }
   
    //resetting pointer to the start of file
    fseek(file, 0, SEEK_SET);
    
    int *inBoundLinks;
    
    inBoundLinks = malloc(sizeof(int) * N);
    if (inBoundLinks == NULL)
    {
        printf("Could not allocate memory for the inBoundLinks array \n");
        exit(-1);
    }
    
    //initialize the inBoundLinks array
    for (i=0;i<N;i++)
    {
        inBoundLinks[i] = 0;
    }
    
    
    //iterate the file line by line
    while (fgets(line,1000,file) != NULL)
    {
        
        token = strtok(line,s);
    
        
        // check that the line dont start with #
        if (strcmp(token,"#"))
        {
            
            sscanf(line, "%d %d \n", &i, &j);
            
            inBoundLinks[j]++;
        
        }

    }
    


    for (i=0;i<=N;i++)
    {
        Graph[i] = malloc(sizeof(int) * inBoundLinks[i] );
    }
    
    //resetting pointer to the start of file
    fseek(file, 0, SEEK_SET);
    
    
    int *counter;
    
    counter = malloc(sizeof(int) * N);
    
    if (counter == NULL)
    {
        printf("Could not allocate memory for the counter array \n");
        exit(-1);
    }
    
    //initialize counter array
    for (i=0;i<N;i++)
    {
        counter[i] = 0;
    }
    
    while (fgets(line,1000,file) != NULL)
    {
        
        token = strtok(line,s);
    
        
        // check that the line dont start with #
        if (strcmp(token,"#"))
        {
            
            sscanf(line, "%d %d \n", &i, &j);
            
            Graph[j][counter[j]] = i;
        
            counter[j]++;
        }

    }
    
    
    double *pageRankVector;
    
    pageRankVector = malloc(sizeof(double) * N);
    if (pageRankVector==NULL)
    {
        printf("Could not allocate memory for the pageRankVector array \n");
        exit(-1);
    }
    
    // The matrix with the random jump probabilities.
    double *E;
    
    E = malloc(sizeof(double) * N);
    if (E == NULL)
    {
        printf("Could not allocate memory for the E array \n");
        exit(-1);
    }
    
    
    for (i=0;i<N;i++)
    {
        E[i] = 1.0 / N ;
    }
    
    gettimeofday (&startwtime, NULL);
    
    pageRankGaussSeidel(N,Graph,outBoundLinks,inBoundLinks,pageRankVector,E,&iter);
    
    gettimeofday (&endwtime, NULL);  
    
    
    time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

    /*
    for (i=0;i<30;i++)
    {
        printf("%.20f \n",pageRankVector[i]);
    }
    */
    printf("Serial PageRank time: %f \n", time);
    
    printf("Number of iterations: %d \n",iter);
    
    //test the results
    //test(N,pageRankVector,testFile);	
	
    // deallocate memory
    free(inBoundLinks);
    free(outBoundLinks);
    free(pageRankVector);
    free(counter);
    free(E);
    
    for (i=0;i<N;i++)
    {
        free(Graph[i]);
    }
    
    free(Graph);
    
   
    return 0;
    
}
