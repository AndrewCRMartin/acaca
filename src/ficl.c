/*************************************************************************

   Program:    ficl
   File:       ficl.c
   
   Version:    V3.6
   Date:       09.01.96
   Function:   Find the cluster into which a PDB loop fits
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995-6
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 387 7050 X 3284
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  26.07.95 Original
   V3.4  10.09.95 Skipped
   V3.5  06.11.95 Skipped
   V3.6  09.01.96 Skipped
   V3.6a 30.01.09 Compile cleanups

*************************************************************************/
/* Includes
*/
#define MAIN
#include "acaca.h"
#include "bioplib/matrix.h"

/************************************************************************/
/* Defines and macros
*/
#define HUGEBUFF 640


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *datafile, char *pdbfile, 
                  char *startres, char *lastres, BOOL *CATorsions,
                  BOOL *Verbose);
REAL **ReadClusterFile(char *datafile, BOOL CATorsions, int *pMethod,
                       int *pNData, int *pVecLength, 
                       CLUSTER **ppClusters, int *NClusters,
                       CLUSTER **ppMedians,  int *NMedians);
BOOL ReadData(FILE *fp, REAL **data, int NLoops, int VecLength);
BOOL ReadHeader(FILE *fp, int *pMethod, int *pNLoops, int *pMaxLen);
REAL **AllocateDataArrays(int NLoops, int VecLength, CLUSTER **ppClusters);
void Usage(void);
void CleanUp(REAL **data1, int NData1, int VecLen1,
             REAL **data2, int NData2, int VecLen2);
int ReadClusters(FILE *fp, CLUSTER *clusters);
int MatchCluster(REAL **data, int NData, int VecLength, 
                 CLUSTER *clusters, int NClusters, REAL *LoopData, 
                 BOOL CATorsions, int method, BOOL *pError);
int ConfirmCluster(REAL **data, int NVec, int VecLen, 
                   CLUSTER *clusters, int TheCluster, REAL *vector,
                   BOOL *pError);
int InClusterBounds(REAL **data, int NVec, int VecDim, CLUSTER *clusters,
                    int ClusNum, REAL *vector);
REAL MinDistInCluster(REAL **data, int NVec, int VecLen, 
                      CLUSTER *clusters, REAL *vector, int ClusNum);
int FindNearestMedian(REAL **data, int NVec, int VecLen, 
                      CLUSTER *clusters, int NClusters, REAL *vector);
REAL *FindMedian(REAL **data, int NVec, int VecLen, CLUSTER *clusters, 
                 int ClusNum);
int ReadMedians(FILE *fp, CLUSTER **ppMedians);
void PrintClusterInfo(int TheCluster, CLUSTER *MedianData, int NMedians,
                      REAL dist, BOOL Verbose);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for scanning a loop against a set of cluster definitions

   26.07.95 Original    By: ACRM
   31.07.95 Moved all result printing into PrintClusterInfo()
*/
int main(int argc, char **argv)
{
   char    datafile[MAXBUFF],
           pdbfile[MAXBUFF],
           startres[16],
           lastres[16];
   int     retval = 0,
           NData,
           NLoopData,
           NMedians,
           method,
           NClusters,
           VecLength,
           TheCluster;
   BOOL    CATorsions = TRUE,
           Error      = FALSE,
           Verbose    = FALSE;
   REAL    **data     = NULL,
           **LoopData = NULL,
           dist;
   CLUSTER *clusters  = NULL,
           *medians   = NULL;
   
   gOutfp = stdout;

   if(ParseCmdLine(argc, argv, datafile, pdbfile, startres, lastres, 
                   &CATorsions, &Verbose))
   {
      if((data=ReadClusterFile(datafile, CATorsions, &method, 
                               &NData, &VecLength, 
                               &clusters, &NClusters,
                               &medians,  &NMedians))!=NULL)
      {
         if(HandleLoopSpec(pdbfile, startres, lastres, CATorsions, FALSE))
         {
            if((LoopData=ConvertData(gDataList,&NLoopData,CATorsions))
               !=NULL)
            {
               TheCluster = MatchCluster(data, NData, VecLength, 
                                         clusters, NClusters, 
                                         LoopData[0], CATorsions, 
                                         method, &Error);
               if(TheCluster == 0 && Error)
               {
                  fprintf(stderr,"Cluster matching failed (memory)\n");
                  retval = 1;
               }
               else
               {
                  dist = MinDistInCluster(data, NData, VecLength, 
                                          clusters, LoopData[0], 
                                          ABS(TheCluster));
                  PrintClusterInfo(TheCluster, medians, NMedians, dist,
                                   Verbose);
               }
            }
            else
            {
               fprintf(stderr,"Unable to get torsion data from loop\n");
               retval = 1;
            }
         }
         else
         {
            fprintf(stderr,"Failure in reading loop\n");
            retval = 1;
         }
      }
      else
      {
         fprintf(stderr,"Error reading cluster file: %s\n",datafile);
         retval = 1;
      }
   }
   else
   {
      Usage();
   }

   CleanUp(data, NData, VecLength, LoopData, NLoopData, VecLength);

   return(retval);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *datafile, 
                     char *pdbfile, char *startres, char *lastres, 
                     BOOL *CATorsions, BOOL *Verbose)
   ---------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *datafile    Cluster data file
            char   *startres    Start residue spec
            char   *lastres     Last residue spec
            BOOL   *CATorsions  Do CA torsions?
            BOOL   *Verbose     Print verbose information
   Returns: BOOL                Success?

   Parse the command line
   
   26.07.95 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *datafile, char *pdbfile, 
                  char *startres, char *lastres, BOOL *CATorsions,
                  BOOL *Verbose)
{
   argc--;
   argv++;
   
   datafile[0] = '\0';
   pdbfile[0]  = '\0';
   startres[0] = '\0';
   lastres[0]  = '\0';
   *CATorsions = TRUE;

   /* Handle the switches                                               */
   while(argc && argv[0][0] == '-')
   {
      switch(argv[0][1])
      {
      case 't':
         *CATorsions = FALSE;
         break;
      case 'v':
         *Verbose = TRUE;
         break;
      default:
         return(FALSE);
         break;
      }
      argc--;
      argv++;
   }

   /* Check there are four additional arguments                         */
   if(argc != 4)
      return(FALSE);
         
   /* Copy the arguments to variables                                   */
   strcpy(datafile, argv[0]);
   strcpy(pdbfile,  argv[1]);
   strcpy(startres, argv[2]);
   strcpy(lastres,  argv[3]);

   return(TRUE);
}


/************************************************************************/
/*>REAL **ReadClusterFile(char *datafile, BOOL CATorsions, int *pMethod,
                          int *pNData, int *pVecLength, 
                          CLUSTER **ppClusters, int *NClusters,
                          CLUSTER **ppMedians,  int *NMedians)
   ---------------------------------------------------------------------
   Read the file produced by CLAN which defines the known clusters

   Returns: REAL   **     Array of cluster data
                          NULL on error

   26.07.95 Original    By: ACRM
*/
REAL **ReadClusterFile(char *datafile, BOOL CATorsions, int *pMethod,
                       int *pNData, int *pVecLength, 
                       CLUSTER **ppClusters, int *NClusters,
                       CLUSTER **ppMedians,  int *NMedians)
{
   FILE *fp;
   BOOL ok = TRUE;
   int  NLoops,
        MaxLen;
   REAL **data = NULL;
   
   /* Open the datafile for reading                                     */
   if((fp=fopen(datafile,"r"))==NULL)
      return(NULL);
   
   /* Read in the required header info                                  */
   if(ReadHeader(fp, pMethod, &NLoops, &MaxLen))
   {
      *pNData     = NLoops;
      *pVecLength = MaxLen * 2 * (CATorsions ? 1 : 3);

      if((data=AllocateDataArrays(NLoops, *pVecLength, ppClusters))!=NULL)
      {
         /* Read in the clustering data                                 */
         if(ReadData(fp, data, NLoops, *pVecLength))
         {
            /* Read the cluster table                                   */
            if((*NClusters = ReadClusters(fp, *ppClusters))==0)
            {
               fprintf(stderr,"Unable to read CLUSTABLE section in \
CLAN output\n");
               ok = FALSE;
            }
            else
            {
               if((*NMedians = ReadMedians(fp, ppMedians))==0)
               {
                  fprintf(stderr,"Unable to read MEDIANS section in \
CLAN output\n");
                  ok = FALSE;
               }
            }
         }
         else
         {
            fprintf(stderr,"Unable to read DATA section in CLAN \
output\n");
            ok = FALSE;
         }
      }
      else
      {
         fprintf(stderr,"Unable to allocate data arrays\n");
         ok = FALSE;
      }
   }
   else
   {
      fprintf(stderr,"Unable to read HEADER section in CLAN output\n");
      ok = FALSE;
   }

   /* Close the file and return                                         */
   fclose(fp);
   if(!ok)
   {
      if(data!=NULL)
         FreeArray2D((char **)data, NLoops, 
                     (MaxLen * 2 * (CATorsions ? 1 : 3)));
      data = NULL;
   }
      
   return(data);
}


/************************************************************************/
/*>BOOL ReadData(FILE *fp, REAL **data, int NLoops, int VecLength)
   ---------------------------------------------------------------
   Read the DATA section from the CLAN output file

   Returns: BOOL               Success?

   26.07.95 Original    By: ACRM
*/
BOOL ReadData(FILE *fp, REAL **data, int NLoops, int VecLength)
{
   char buffer[HUGEBUFF],
        word[MAXBUFF],
        *p;
   BOOL InSection = FALSE;
   int  i,
        LoopCount = 0;

   rewind(fp);
   
   while(fgets(buffer,HUGEBUFF,fp))
   {
      TERMINATE(buffer);

      /* Check for end of data section                                  */
      if(!strncmp(buffer,"END DATA",8))
         return(TRUE);

      /* If in section, check for required data                         */
      if(InSection)
      {
         for(i=0, p=buffer; p!=NULL && i<VecLength; i++)
         {
            p = GetWord(p,word,MAXBUFF);
            if(!sscanf(word,"%lf",&(data[LoopCount][i])))
               return(FALSE);
         }
         
         if(++LoopCount > NLoops)
            return(FALSE);
      }

      /* Check for beginning of data section                            */
      if(!strncmp(buffer,"BEGIN DATA",10))
         InSection = TRUE;
   }
   
   return(FALSE);
}


/************************************************************************/
/*>BOOL ReadHeader(FILE *fp, int *pMethod, int *pNLoops, int *pMaxLen)
   -------------------------------------------------------------------
   Read the HEADER section from the CLAN output file

   Returns: BOOL               Success?

   26.07.95 Original    By: ACRM
*/
BOOL ReadHeader(FILE *fp, int *pMethod, int *pNLoops, int *pMaxLen)
{
   char buffer[MAXBUFF],
        word[MAXBUFF],
        *p;
   BOOL InSection = FALSE,
        GotMethod = FALSE,
        GotLength = FALSE,
        GotNLoops = FALSE,
        GotScheme = FALSE;
   int  i;

   rewind(fp);
   
   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);

      /* Check for beginning and end of data section                    */
      if(!strncmp(buffer,"BEGIN HEADER",12))
         InSection = TRUE;
      if(!strncmp(buffer,"END HEADER",10))
         return(GotMethod && GotLength && GotNLoops && GotScheme);

      /* If in section, check for required data                         */
      if(InSection)
      {
         p = GetWord(buffer,word,MAXBUFF);

         if(!strncmp(word,"METHOD",6))
         {
            if(sscanf(p,"%d",pMethod))
               GotMethod = TRUE;
         }
         else if(!strncmp(word,"NLOOPS",6))
         {
            if(sscanf(p,"%d",pNLoops))
               GotNLoops = TRUE;
         }
         else if(!strncmp(word,"MAXLENGTH",9))
         {
            if(sscanf(p,"%d",pMaxLen))
               GotLength = TRUE;
         }
         else if(!strncmp(word,"SCHEME",6))
         {
            /* Read the scheme out of the following values              */
            for(i=0;p!=NULL && i<MAXLOOPLEN;i++)
            {
               p = GetWord(p,word,MAXBUFF);
               sscanf(word,"%d",&(gScheme[i]));
            }
            gMaxLoopLen = i;
            GotScheme = TRUE;
         }
      }
   }
   
   return(GotMethod && GotLength && GotNLoops && GotScheme);
}


/************************************************************************/
/*>REAL **AllocateDataArrays(int NLoops, int VecLength, 
                             CLUSTER **ppClusters)
   ----------------------------------------------------
   Allocate array to store the data section from the CLAN output

   Returns: REAL    **   Data matrix storage
                         NULL on error

   26.07.95 Original    By: ACRM
*/
REAL **AllocateDataArrays(int NLoops, int VecLength, CLUSTER **ppClusters)
{
   REAL **data;

   /* Allocate and check 2D array                                       */
   if((data = (REAL **)Array2D(sizeof(REAL),NLoops,VecLength))==NULL)
   {
      fprintf(stderr,"No memory for data array.\n");
      return(NULL);
   }

   /* Allocate memory for clusters                                      */
   if((*ppClusters = (CLUSTER *)malloc(NLoops * sizeof(CLUSTER)))==NULL)
   {
      FreeArray2D((char **)data, NLoops, VecLength);
      return(NULL);
   }
   
   return(data);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Print a usage message

   26.07.95 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nficl (c) 1995 Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: ficl [-t] clusterfile pdbfile startres \
lastres\n");
   fprintf(stderr,"       -t    Use true torsions rather than CA \
pseudo-torsions\n");

   fprintf(stderr,"\nTakes the output from CLAN and compares a loop in \
a PDB file with the\n");
   fprintf(stderr,"clusters defined in the CLAN file. Outputs the \
cluster into which this\n");
   fprintf(stderr,"loop falls or an indication this loop does not match \
any of the existing\n");
   fprintf(stderr,"clusters.\n");
}


/************************************************************************/
/*>void CleanUp(REAL **data1, int NData1, int VecLen1, 
                REAL **data2, int NData2, int VecLen2)
   ---------------------------------------------------
   Free allocated data

   26.07.95 Original    By: ACRM
*/
void CleanUp(REAL **data1, int NData1, int VecLen1, 
             REAL **data2, int NData2, int VecLen2)
{
   FreeArray2D((char **)data1, NData1, VecLen1);
   FreeArray2D((char **)data2, NData2, VecLen2);
}


/************************************************************************/
/*>int ReadClusters(FILE *fp, CLUSTER *clusters)
   ---------------------------------------------
   Reads in the cluster table

   Returns: int                  Number of clusters read
                                 0 on error

   27.07.95 Original    By: ACRM
   30.01.09 Initialize some variables
*/
int ReadClusters(FILE *fp, CLUSTER *clusters)
{
   char    buffer[HUGEBUFF],
           word[MAXBUFF],
           *p;
   int     NClusters = 0,
           SecPos    = 0,
           i,
           loopnum = 0;

   rewind(fp);
   
   /* First find the number of clusters                                 */
   while(fgets(buffer,HUGEBUFF,fp))
   {
      TERMINATE(buffer);

      if(!strncmp(buffer,"BEGIN MEDIANS",13))
      {
         sscanf(buffer+13,"%d",&NClusters);
         break;
      }
   }
   
   if(NClusters == 0)
   {
      fprintf(stderr,"Unable to find BEGIN MEDIANS statement with number \
of clusters\n");
      return(0);
   }

   /* Rewind the file and find the appropriate cluster numbers from the
      cluster table
   */
   rewind(fp);
   
   while(fgets(buffer,HUGEBUFF,fp))
   {
      TERMINATE(buffer);
      p = buffer;

      /* Check for end of data section                                  */
      if(!strncmp(buffer,"END CLUSTABLE",13))
         return(NClusters);

      /* If in section, check for required data                         */
      if(SecPos)
      {
         if(SecPos == 3)
            loopnum=0;
         
         if(SecPos > 2)
         {
            p = GetWord(p,(clusters[loopnum].loopid),MAXLOOPID);
            for(i=0; i<NClusters; i++)
               p = GetWord(p,word,MAXBUFF);

            sscanf(word,"%d",&(clusters[loopnum].clusnum));
         }
         SecPos++;
         loopnum++;
      }  
   
      /* Check for beginning of data section                            */
      if(!strncmp(buffer,"BEGIN CLUSTABLE",15))
         SecPos = 1;
   }
   
   return(NClusters);
}


/************************************************************************/
/*>int MatchCluster(REAL **data, int NData, int VecLength, 
                    CLUSTER *clusters, int NClusters, REAL *LoopData, 
                    BOOL CATorsions, int method, BOOL *pError)
   ---------------------------------------------------------------
   Returns: int             Appropriate cluster 
                            -ve if cluster is a singleton
                            0 if not a member of any cluster 
                            or if error (check the flag).

   Match a vector in LoopData against the clusters using method

   Algorithm is:

      C = FindNearestMedian();
      TheCluster = C;
      Dmin = DistToNearestVectorInCluster(C);
      for(each cluster, D, except C)
      {  if(InBoundsOfCluster(D))
         {  d = DistToNearestVectorInCluster(D);
            if(d < DMin)
            {  DMin = d;
               TheCluster = D;
            }
         }
      }

   28.07.95 Original (Note `method' is not used)   By: ACRM
   31.07.95 Added error check from InClusterBounds()
            Added error check from FindNearestMedian()
*/
int MatchCluster(REAL **data, int NData, int VecLength, 
                 CLUSTER *clusters, int NClusters, REAL *LoopData, 
                 BOOL CATorsions, int method, BOOL *pError)
{
   int  i,
        TheCluster = 0,
        C, D,
        ok;
   REAL DMin,
        dist;

   *pError = FALSE;

   /* Find the closest cluster median and the distance to the nearest
      item in that cluster
   */
   TheCluster = C = FindNearestMedian(data, NData, VecLength, clusters, 
                                      NClusters, LoopData);
   if(C==0)
   {
      *pError = TRUE;
      return(0);
   }
   DMin = MinDistInCluster(data, NData, VecLength, clusters, LoopData, C);

   /* Test each other cluster to see if we are in the bounds of that
      cluster
   */
   for(i=0; i<NClusters; i++)
   {
      D = i+1;
      
      if(D != C)
      {
         if((ok=InClusterBounds(data, NData, VecLength, clusters, D, 
                                LoopData)) > 0)
         {
            /* If we are in this cluster's bounds, then see if there
               is a point in this cluster which is closer than the
               nearest point in the previous cluster. If so, then make
               this our cluster
            */
            if((dist = MinDistInCluster(data, NData, VecLength, clusters, 
                                        LoopData, D)) < DMin)
            {
               DMin = dist;
               TheCluster = D;
            }
         }
         else if(ok<0)
         {
            fprintf(stderr,"InClusterBounds() out of memory\n");
            *pError = TRUE;
            return(0);
         }
      }
   }
   
   /* Having established the most likely cluster, test whether our
      vector is really a member of this cluster or whether it simply
      happens to be the closest
   */
   TheCluster = ConfirmCluster(data, NData, VecLength, clusters, 
                               TheCluster, LoopData, pError);
   return(TheCluster);
}


/************************************************************************/
/*>int ConfirmCluster(REAL **data, int NVec, int VecLen, 
                      CLUSTER *clusters, int TheCluster, REAL *vector,
                      BOOL *pError)
   -------------------------------------------------------------------
   Returns: int                   The cluster number
                                  0 if cluster mis-match
                                  0 if error (check error flag)

   Confirms that a vector really is a member of a cluster by checking
   that there are points closer to this vector than the median is
   and seeming that the bounding box does not increase by more than
   50% in any dimension.

   28.07.95 Original    By: ACRM
   30.01.09 Initialize some variables
*/
int ConfirmCluster(REAL **data, int NVec, int VecLen, 
                   CLUSTER *clusters, int TheCluster, REAL *vector,
                   BOOL *pError)
{
   REAL *median,
        DistMedian,
        DistNearest,
        dist,
        MaxVal = 0.0,
        MinVal = 0.0;
   int  i, j, ok,
        NMembers;

   *pError = FALSE;
   
   /* If we are out of the bounds of the cluster                        */
   if((ok=InClusterBounds(data,NVec,VecLen,clusters,TheCluster,vector))
      ==0)
   {
      /* First ensure that we are closer to one of the points in the 
         cluster than we are to the median
      */
      if((median=FindMedian(data,NVec,VecLen,clusters,TheCluster))==NULL)
      {
         *pError = TRUE;
         return(0);
      }
      
      DistMedian  = VecDist(vector, median, VecLen);
      free(median);
      median = NULL;
      
      DistNearest = MinDistInCluster(data, NVec, VecLen, clusters, vector,
                                     TheCluster);
      if(DistNearest > DistMedian)
         return(0);

      /* Now ensure that we are not expanding the cluster's bounding
         box by more than 50%
      */
      for(j=0; j<VecLen; j++)     /* For each dimension                 */
      {
         NMembers = 0;
         
         /* Find bounds of cluster                                      */
         for(i=0; i<NVec; i++)
         {
            if(clusters[i].clusnum==TheCluster)
            {
               if(!NMembers)
               {
                  MinVal = MaxVal = data[i][j];
                  NMembers++;
               }
               else
               {
                  if(data[i][j] < MinVal)
                     MinVal = data[i][j];
                  else if(data[i][j] > MaxVal)
                     MaxVal = data[i][j];
                  NMembers++;
               }
            }
         }

         /* If the cluster has only one member return the nigative 
            version of the cluster number
         */
         if(NMembers == 1)
            return(-TheCluster);

         /* Find distance between bounds                                */
         dist = MaxVal-MinVal;
         
         /* See if the new bounds exceed the old by more than 50%       */
         if(vector[j] > MaxVal)
         {
            if((vector[j] - MinVal) > (REAL)1.5*dist)
               return(0);
         }
         else if(vector[j] < MinVal)
         {
            if((MaxVal - vector[j]) > (REAL)1.5*dist)
               return(0);
         }
      }
      
   }
   else if(ok<0)
   {
      *pError = TRUE;
      return(0);
   }

   return(TheCluster);
}


/************************************************************************/
/*>int InClusterBounds(REAL **data, int NVec, int VecDim, 
                       CLUSTER *clusters, int ClusNum, REAL *vector)
   -----------------------------------------------------------------
   Sees if a vector is within the bounds of a cluster.
   The bounds are extended by 10% to account for rounding error resulting
   from reading the cluster data from a file rather than calculating
   true values.

   31.07.95 Original    By: ACRM
*/
int InClusterBounds(REAL **data, int NVec, int VecDim, CLUSTER *clusters,
                    int ClusNum, REAL *vector)
{
   int      i, j,
            retval = 1;
   REAL     *minval, 
            *maxval,
            dist;
   BOOL     Done   = FALSE;

   /* Allocate arrays to store min and max values in each dimension     */
   minval=(REAL *)malloc(VecDim*sizeof(REAL));
   maxval=(REAL *)malloc(VecDim*sizeof(REAL));

   if(minval==NULL || maxval==NULL)
   {
      if(minval!=NULL) free(minval);
      if(maxval!=NULL) free(maxval);

      return(-1);
   }

   /* Find the min and max values in each dimension                     */
   for(i=0, Done=FALSE; i<NVec; i++)
   {
      if(clusters[i].clusnum == ClusNum)
      {
         /* On first item, just copy in the data                        */
         if(!Done)
         {
            for(j=0; j<VecDim; j++)
               minval[j] = maxval[j] = data[i][j];
            Done = TRUE;
         }
         else
         {
            for(j=0; j<VecDim; j++)
            {
               if(data[i][j] < minval[j])
                  minval[j] = data[i][j];
               if(data[i][j] > maxval[j])
                  maxval[j] = data[i][j];
            }
         }
      }
   }

   /* Expand the bounds by 10%                                          */
   for(j=0; j<VecDim; j++)
   {
      dist = maxval[j] - minval[j];
      if(dist==(REAL)0.0)
         dist = ABS(minval[j]);
         
      dist /= (REAL)10.0;
      minval[j] -= dist;
      maxval[j] += dist;
   }

   /* Check whether our vector falls within these bounds                */
   for(j=0; j<VecDim; j++)
   {
      if((vector[j] < minval[j]) || (vector[j] > maxval[j]))
      {
         retval = 0;
         break;
      }
   }
      
   /* Free up the arrays                                                */
   free(minval);
   free(maxval);

   return(retval);
}


/************************************************************************/
/*>REAL MinDistInCluster(REAL **data, int NVec, int VecLen, 
                         CLUSTER *clusters, REAL *vector, int ClusNum)
   -------------------------------------------------------------------
   Returns the minimum distance from the vector to a member of the
   cluster.

   31.07.95 Original    By: ACRM
*/
REAL MinDistInCluster(REAL **data, int NVec, int VecLen, 
                      CLUSTER *clusters, REAL *vector, int ClusNum)
{
   int  i;
   REAL DMin = INF,
        dist;
   
   for(i=0; i<NVec; i++)
   {
      if(clusters[i].clusnum == ClusNum)
      {
         dist = VecDist(vector,data[i],VecLen);
         if(dist < DMin)
            DMin = dist;
      }
   }
   
   return(DMin);
}


/************************************************************************/
/*>int FindNearestMedian(REAL **data, int NVec, int VecLen, 
                         CLUSTER *clusters, int NClusters, REAL *vector)
   ---------------------------------------------------------------------
   Finds the cluster with the median closest to the vector.
   Returns 0 on error.

   31.07.95 Original    By: ACRM
*/
int FindNearestMedian(REAL **data, int NVec, int VecLen, 
                      CLUSTER *clusters, int NClusters, REAL *vector)
{
   int  i, j,
        ClusNum = 0;
   REAL **medians,
        dist,
        DMin = INF;
   
   /* Allocate memory to store median arrays                            */
   if((medians = (REAL **)malloc(NClusters * sizeof(REAL *)))==NULL)
      return(0);

   /* For each cluster, find the median                                 */
   for(i=0; i<NClusters; i++)
   {
      if((medians[i] = FindMedian(data, NVec, VecLen, clusters, i+1))
         ==NULL)
      {
         /* Free previous medians if allocation failed                  */
         for(j=0; j<i; j++)
            free(medians[j]);
         free(medians);
         return(0);
      }
   }
   
   /* Find the nearest median                                           */
   for(i=0; i<NClusters; i++)
   {
      if((dist = VecDist(vector, medians[i], VecLen)) < DMin)
      {
         DMin    = dist;
         ClusNum = i+1;
      }
   }

   /* Free allocated memory                                             */
   for(i=0; i<NClusters; i++)
      free(medians[i]);
   free(medians);

   /* Return the nearest cluster                                        */
   return(ClusNum);
}


/************************************************************************/
/*>REAL *FindMedian(REAL **data, int Nvec, int VecLen, 
                    CLUSTER *clusters, int ClusNum)
   ----------------------------------------------------
   Find the median of cluster ClusNum.
   Returns a REAL pointer (array) containing the median vector. This
   must be freed after use!
   Returns NULL on error.

   28.07.95 Original based on code from cluster.c    By: ACRM
*/
REAL *FindMedian(REAL **data, int NVec, int VecLen, CLUSTER *clusters, 
                 int ClusNum)
{
   int      i, j;
   REAL     *minval, 
            *maxval,
            *medval;
   BOOL     Done = FALSE;

   /* Allocate arrays to store min and max values in each dimension and 
      the median which will be returned.
   */
   minval=(REAL *)malloc(VecLen*sizeof(REAL));
   maxval=(REAL *)malloc(VecLen*sizeof(REAL));
   medval=(REAL *)malloc(VecLen*sizeof(REAL));

   if(minval==NULL || maxval==NULL || medval==NULL)
   {
      if(minval!=NULL) free(minval);
      if(maxval!=NULL) free(maxval);
      if(medval!=NULL) free(medval);

      return(NULL);
   }

   /* Find the min and max values in each dimension                     */
   for(i=0, Done=FALSE; i<NVec; i++)
   {
      if(clusters[i].clusnum == ClusNum)
      {
         /* On first item, just copy in the data                        */
         if(!Done)
         {
            for(j=0; j<VecLen; j++)
               minval[j] = maxval[j] = data[i][j];
            Done = TRUE;
         }
         else
         {
            for(j=0; j<VecLen; j++)
            {
               if(data[i][j] < minval[j])
                  minval[j] = data[i][j];
               if(data[i][j] > maxval[j])
                  maxval[j] = data[i][j];
            }
         }
      }
   }

   /* Now store the median values                                       */
   for(j=0; j<VecLen; j++)
      medval[j] = (minval[j] + maxval[j]) / (REAL)2.0;
      
   /* Free up the arrays                                                */
   free(minval);
   free(maxval);

   /* Return the medians array                                          */
   return(medval);
}


/************************************************************************/
/*>int ReadMedians(FILE *fp, CLUSTER **ppMedians)
   ----------------------------------------------
   Read the MEDIANS section from the CLAN output file.
   Returns the number of clusters or 0 on error

   31.07.95 Original    By: ACRM
*/
int ReadMedians(FILE *fp, CLUSTER **ppMedians)
{
   char    buffer[HUGEBUFF],
           *p;
   int     NClusters = 0,
           count     = 0;
   BOOL    InSection = FALSE;

   *ppMedians = NULL;

   rewind(fp);
   
   /* First find the number of clusters                                 */
   while(fgets(buffer,HUGEBUFF,fp))
   {
      TERMINATE(buffer);
      p = buffer;

      /* Check for end of data section                                  */
      if(!strncmp(buffer,"END MEDIANS",11))
         return(NClusters);

      /* If in section, read the median data                            */
      if(InSection)
      {
         if(*ppMedians == NULL)
            return(0);
         
         if((sscanf(buffer,"%d %s", &(((*ppMedians)[count]).clusnum),
                    ((*ppMedians)[count]).loopid))!=2)
         {
            fprintf(stderr,"This is an old CLAN file without cluster \
numbers in MEDIANS\n");
            free(*ppMedians);
            return(0);
         }
         
         count++;
      }  
   
      /* Check for beginning of data section                            */
      if(!strncmp(buffer,"BEGIN MEDIANS",13))
      {
         sscanf(buffer+13,"%d",&NClusters);

         InSection = TRUE;

         if(NClusters == 0)
         {
            fprintf(stderr,"Unable to find BEGIN MEDIANS statement with \
number of clusters\n");
            return(0);
         }

         /* Allocate memory for median storage                          */
         if((*ppMedians = (CLUSTER *)malloc(NClusters * sizeof(CLUSTER)))
            ==NULL)
            return(0);
      }
   }
   
   return(NClusters);
}


/************************************************************************/
/*>void PrintClusterInfo(int TheCluster, CLUSTER *MedianData, 
                         int NMedians, REAL dist, BOOL Verbose)
   ------------------------------------------------------------
   Prints details about the best cluster match

   28.07.95 Original    By: ACRM
   31.07.95 Moved printing of distance in here and added verbose version
*/
void PrintClusterInfo(int TheCluster, CLUSTER *MedianData, int NMedians,
                      REAL dist, BOOL Verbose)
{
   int  i;

   if(TheCluster == 0)
   {
      if(Verbose)
      {
         printf("No cluster found\n");
      }
      else
      {
         printf("Best: 0 Representitive: (none) NOMATCH Distance: \
9999.000\n");
      }
   }
   else
   {
      if(Verbose)
      {
         printf("Cluster %d\n", ABS(TheCluster));
         
         for(i=0; i<NMedians; i++)
         {
            if(MedianData[i].clusnum == ABS(TheCluster))
            {
               printf("Representitive for this cluster is: %s\n",
                      MedianData[i].loopid);
               break;
            }
         }
         
         if(TheCluster < 0)
         {
            printf("Note, however, that there is only one structure in \
this cluster, so\n");
            printf("it is not possible to see how well the conformation \
fits into the\n");
            printf("cluster.\n");
         }
         printf("The distance of this conformation (in cluster space) \
from the nearest\n");
         printf("member of the cluster is %f\n",dist);
      }
      else
      {
         printf("Cluster: %d ", ABS(TheCluster));
         
         for(i=0; i<NMedians; i++)
         {
            if(MedianData[i].clusnum == ABS(TheCluster))
            {
               printf("Representitive: %s ",
                      MedianData[i].loopid);
               break;
            }
         }
         
         printf("%s ",((TheCluster < 0)?"SINGLETON":"CLUSTER"));
         printf(" Distance: %f\n",dist);
      }
   }
}

