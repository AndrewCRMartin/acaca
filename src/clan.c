/*************************************************************************

   Program:    clan
   File:       clan.c
   
   Version:    V3.8
   Date:       16.01.23
   Function:   Perform cluster analysis on loop conformations
   
   Copyright:  (c) UCL, Prof. Andrew C. R. Martin 1995-2023
   Author:     Prof. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      +44 (0)207 679 7034
   EMail:      andrew@bioinf.org.uk
               
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
   V1.0  04.07.95 Original
   V1.1  05.07.95 Added true torsion option
   V1.2  06.07.95 Added post-clustering RMS check
   V1.3  07.07.95 Added Keywords to control post-clustering
   V1.4  24.07.95 Prints loop id's in cluster table
                  Prints raw clustering data
   V2.0  25.07.95 Rearranged source to share code between parts of the
                  acaca package
   V2.1  16.08.95 Some code reorganisation to maintain list of cluster
                  assignments being used in a 1D rather than 2D array.
   V2.2  12.09.95 During postclustering, if there are only 2 members in
                  a cluster, both must be less than the cutoffs in order
                  to merge the clusters.
   V3.0  13.09.95 Clustering now considers a set of distances from the
                  N-ter C-alpha as well as torsions
   V3.1  21.09.95 Added DISTANCE/ANGLE/NOANGLE keywords
   V3.2  02.10.95 Modifications to decr.c
   V3.3  04.10.95 Modified critical residue definition to show residues
                  which are conserved in at least one cluster.
   V3.4  10.10.95 Various changes to make deleted residues work with 
                  residues conserved in at least one cluster
   V3.5  06.11.95 Added EXCLUDE command for excluding from template 
                  analysis
   V3.6  09.01.96 Added code to free up unneeded storage when not doing
                  critical residues
   V3.7  14.03.96 Cluster merging now considers CB as well
   V3.7a 30.01.09 Compile cleanups
   V3.8  16.01.23 Fixed some bugs running under Linux

*************************************************************************/
/* Includes
*/
#define MAIN
#include "acaca.h"
#include "decr.h"


/************************************************************************/
/* Defines and macros
*/
#define KEY_METHOD           0
#define KEY_LOOP             1
#define KEY_OUTPUT           2
#define KEY_MAXLENGTH        3
#define KEY_SCHEME           4
#define KEY_DENDOGRAM        5
#define KEY_TABLE            6
#define KEY_POSTCLUSTER      7
#define KEY_DATA             8
#define KEY_CRITICAL         9
#define KEY_INFO             10
#define KEY_NODISTANCE       11
#define KEY_DISTANCE         12
#define KEY_NOANGLE          13
#define KEY_ANGLE            14
#define KEY_TRUETORSIONS     15
#define KEY_PSEUDOTORSIONS   16
#define KEY_EXCLUDE          17
#define PARSER_NCOMM         18
#define PARSER_MAXSTRPARAM   3
#define PARSER_MAXSTRLEN     80
#define PARSER_MAXREALPARAM  MAXLOOPLEN

#define UP '|'
#define ACROSS '-'
#define BLANK ' '

/*  Map row I and column J of upper half diagonal symmetric matrix 
    onto vector.
*/
#define IOFFSET(n,i,j) (j+(i-1)*n-(i*(i+1))/2)


/************************************************************************/
/* Globals
*/
static MKeyWd sKeyWords[PARSER_NCOMM];         /* Parser keywords       */
static char   *sStrParam[PARSER_MAXSTRPARAM];  /* Parser string params  */
static REAL   sRealParam[PARSER_MAXREALPARAM]; /* Parser real params    */
static int    sInfoLevel = 0;                  /* Info level            */


/************************************************************************/
/* Prototypes
*/
#include "clan.p"
#include "acaca.p"
#include "decr.p"
#include "decr2.p"


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for cluster analysis on PDB loops

   27.06.95 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF];
   FILE *fp=NULL;
   int  retval = 0;

   gOutfp = stdout;

   InitProperties();

   gPClusCut[0] = RMSCUT;
   gPClusCut[1] = MAXDEV;
   gPClusCut[2] = MAXCBDEV;

   if(ParseCmdLine(argc, argv, infile, &gCATorsions))
   {
      if((fp = fopen(infile, "r"))!=NULL)
      {
         if(ReadInputFile(fp,gCATorsions))
         {
            if(!DoClustering(gCATorsions))
            {
               fprintf(stderr,"Clustering failed\n");
               retval = 1;
            }
         }
         else
         {
            fprintf(stderr,"Error while reading input file\n");
            retval = 1;
         }
      }
      else
      {
         fprintf(stderr,"Unable to open input file: %s\n",infile);
         retval = 1;
      }
   }
   else
   {
      Usage();
   }

   CleanUp();

   return(retval);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, 
                     BOOL *CATorsions)
   ------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            BOOL   *CATorsions  Do pseudo-CA torsions rather than true
                                torsions.
   Returns: BOOL                Success?

   Parse the command line
   
   26.06.95 Original    By: ACRM
   05.07.95 Added -t
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, BOOL *CATorsions)
{
   argc--;
   argv++;
   
   infile[0] = '\0';

   /* Handle the switches                                               */
   while(argc && argv[0][0] == '-')
   {
      switch(argv[0][1])
      {
      case 't':
         *CATorsions = FALSE;
         break;
      default:
         return(FALSE);
         break;
      }
      argc--;
      argv++;
   }

   /* Check there is one additional argument                            */
   if(argc != 1)
      return(FALSE);
         
   /* Copy the first to infile                                          */
   strcpy(infile, argv[0]);

   return(TRUE);
}


/************************************************************************/
/*>BOOL ReadInputFile(FILE *fp, BOOL CATorsions)
   ---------------------------------------------
   Input:   FILE  *fp         Input file pointer
            BOOL  CATorsions  Do CA-pseudo-torsions rather than true
                              torsions
   Returns: BOOL              Success?

   Calls routines to set up the command parser and to read the control 
   file.

   27.06.95 Original   By: ACRM
*/
BOOL ReadInputFile(FILE *fp, BOOL CATorsions)
{
   if(SetupParser())
   {
      if(DoCmdLoop(fp, CATorsions))
         return(TRUE);
   }
   
   return(FALSE);
}


/************************************************************************/
/*>BOOL SetupParser(void)
   ----------------------
   Returns:   BOOL         Success of memory allocations

   Set up the command parser

   26.06.95 Original    By: ACRM
   07.07.95 Added postcluster
   08.08.95 Added Criticalresidues
   16.08.95 Postcluster takes 1 or 2 parameters
   18.08.95 Added infolevel
   13.09.95 Added nodistance
   21.09.95 Added angle/noangle
   06.11.95 Added exclude
*/
BOOL SetupParser(void)
{
   int i;
   
   /* Allocate memory for the string parameters                         */
   for(i=0; i<PARSER_MAXSTRPARAM; i++)
   {
      if((sStrParam[i] = 
         (char *)malloc(PARSER_MAXSTRLEN * sizeof(char)))==NULL)
      {
         int j;
         
         for(j=0;j<i;j++) 
            free(sStrParam[j]);
         fprintf(stderr,"No memory for parser string array\n");
         
         return(FALSE);
      }
   }
   
   /* Set up the keywords                                               */
   MAKEMKEY(sKeyWords[KEY_METHOD],        "METHOD",          STRING,1,1);
   MAKEMKEY(sKeyWords[KEY_LOOP],          "LOOP",            STRING,3,3);
   MAKEMKEY(sKeyWords[KEY_OUTPUT],        "OUTPUT",          STRING,1,1);
   MAKEMKEY(sKeyWords[KEY_MAXLENGTH],     "MAXLENGTH",       NUMBER,1,1);
   MAKEMKEY(sKeyWords[KEY_SCHEME],        "SCHEME",          NUMBER,1,
                                                             MAXLOOPLEN);
   MAKEMKEY(sKeyWords[KEY_DENDOGRAM],     "DENDOGRAM",       STRING,0,0);
   MAKEMKEY(sKeyWords[KEY_TABLE],         "TABLE",           STRING,0,0);
   MAKEMKEY(sKeyWords[KEY_POSTCLUSTER],   "POSTCLUSTER",     NUMBER,1,3);
   MAKEMKEY(sKeyWords[KEY_DATA],          "DATA",            STRING,0,0);
   MAKEMKEY(sKeyWords[KEY_CRITICAL],      "CRITICALRESIDUES",
                                                             STRING,0,0);
   MAKEMKEY(sKeyWords[KEY_INFO],          "INFOLEVEL",       NUMBER,1,1);
   MAKEMKEY(sKeyWords[KEY_NODISTANCE],    "NODISTANCE",      STRING,0,0);
   MAKEMKEY(sKeyWords[KEY_DISTANCE],      "DISTANCE",        STRING,0,0);
   MAKEMKEY(sKeyWords[KEY_NOANGLE],       "NOANGLE",         STRING,0,0);
   MAKEMKEY(sKeyWords[KEY_ANGLE],         "ANGLE",           STRING,0,0);

   MAKEMKEY(sKeyWords[KEY_TRUETORSIONS],  "TRUETORSIONS",    STRING,0,0);
   MAKEMKEY(sKeyWords[KEY_PSEUDOTORSIONS],"PSEUDOTORSIONS",  STRING,0,0);

   MAKEMKEY(sKeyWords[KEY_EXCLUDE],       "EXCLUDE",         STRING,1,1);
   
   /* Check all allocations OK                                          */
   for(i=0; i<PARSER_NCOMM; i++)
   {
      if(sKeyWords[i].name == NULL)
      {
         fprintf(stderr,"No memory for keywords, or keyword undefined\n");
         return(FALSE);
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL DoCmdLoop(FILE *fp, BOOL CATorsions)
   -----------------------------------------
   Input:   FILE  *fp          Input file pointer
            BOOL  CATorsions   Do CA rather than true torsions
   Returns: BOOL               Success? Fails if illegal input encountered

   Main loop to handle the command parser for the control file

   27.06.95 Original   By: ACRM
   07.07.95 Added postcluster
   08.08.95 Added Criticalresidues
   18.08.95 Added infolevel
   13.09.95 Added nodistance
   21.09.95 Added distance, angle, noangle
   26.09.95 Added truetorsions/pseudotorsions and GotLoop checking
   06.11.95 Added exclude
*/
BOOL DoCmdLoop(FILE *fp, BOOL CATorsions)
{
   char buffer[MAXBUFF],
        loopid[MAXBUFF];
   int  NParams,
        i,
        key;
   BOOL GotLoop = FALSE;
   
   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);
      
      key = mparse(buffer, PARSER_NCOMM, sKeyWords, sRealParam,
                   sStrParam, &NParams);

      switch(key)
      {
      case PARSE_COMMENT:
         break;
      case PARSE_ERRC:
         fprintf(stderr,"Error in command: %s\n",buffer);
         break;
      case PARSE_ERRP:
         fprintf(stderr,"Error in parameters: %s\n",buffer);
         break;
      case KEY_METHOD:
         if(!SetClusterMethod(sStrParam[0]))
            return(FALSE);
         break;
      case KEY_LOOP:
         GotLoop = TRUE;
         if(!HandleLoopSpec(sStrParam[0], sStrParam[1], sStrParam[2],
                            CATorsions, TRUE))
            fprintf(stderr,"Loop skipped!\n");
         break;
      case KEY_OUTPUT:
         if(!SetOutputFile(sStrParam[0]))
            return(FALSE);
         break;
      case KEY_MAXLENGTH:
         if(gMaxLoopLen)   /* Scheme has already been defined           */
         {
            if((int)sRealParam[0] != gMaxLoopLen)
            {
               fprintf(stderr,"The number of items in your scheme \
definition does not match the number\n");
               fprintf(stderr,"specified by MAXLENGTH\n");
               return(FALSE);
            }
         }
         else
         {
            CreateDefaultScheme((int)sRealParam[0]);
         }
         gMaxLoopLen = (int)sRealParam[0];
         break;
      case KEY_SCHEME:
         if(gMaxLoopLen)
         {
            if(gMaxLoopLen != NParams)
            {
               fprintf(stderr,"The number of items in your scheme \
definition does not match the number\n");
               fprintf(stderr,"specified by MAXLENGTH\n");
               return(FALSE);
            }
         }
         gMaxLoopLen = NParams;
         
         for(i=0; i<NParams; i++)
            gScheme[i] = (int)sRealParam[i];
         break;
      case KEY_DENDOGRAM:
         gDoDendogram = TRUE;
         break;
      case KEY_TABLE:
         gDoTable = TRUE;
         break;
      case KEY_POSTCLUSTER:
         gPClusCut[0] = sRealParam[0];
         if(NParams>1)
            gPClusCut[1] = sRealParam[1];
         if(NParams>2)
            gPClusCut[2] = sRealParam[2];
         break;
      case KEY_DATA:
         gDoData = TRUE;
         break;
      case KEY_CRITICAL:
         gDoCritRes = TRUE;
         break;
      case KEY_INFO:
         sInfoLevel = (int)sRealParam[0];
         break;
      case KEY_NODISTANCE:
         if(GotLoop)
         {
            fprintf(stderr,"Error: %s command must appear \
before all LOOP commands\n",sKeyWords[key].name);
            return(FALSE);
         }
         gDoDistance = FALSE;
         break;
      case KEY_DISTANCE:
         if(GotLoop)
         {
            fprintf(stderr,"Error: %s command must appear \
before all LOOP commands\n",sKeyWords[key].name);
            return(FALSE);
         }
         gDoDistance = TRUE;
         break;
      case KEY_NOANGLE:
         if(GotLoop)
         {
            fprintf(stderr,"Error: %s command must appear \
before all LOOP commands\n",sKeyWords[key].name);
            return(FALSE);
         }
         gDoAngles = FALSE;
         break;
      case KEY_ANGLE:
         if(GotLoop)
         {
            fprintf(stderr,"Error: %s command must appear \
before all LOOP commands\n",sKeyWords[key].name);
            return(FALSE);
         }
         gDoAngles = TRUE;
         break;
      case KEY_TRUETORSIONS:
         if(GotLoop)
         {
            fprintf(stderr,"Error: %s command must appear \
before all LOOP commands\n",sKeyWords[key].name);
            return(FALSE);
         }
         gCATorsions = FALSE;
         break;
      case KEY_PSEUDOTORSIONS:
         if(GotLoop)
         {
            fprintf(stderr,"Error: %s command must appear \
before all LOOP commands\n",sKeyWords[key].name);
            return(FALSE);
         }
         gCATorsions = TRUE;
         break;
      case KEY_EXCLUDE:
         sprintf(loopid,"%s-%s-%s",
                 sStrParam[0],sStrParam[1],sStrParam[2]);
         if((gStringList=StoreString(gStringList, loopid))==NULL)
         {
            fprintf(stderr,"Error: No memory for string: %s\n", loopid);
            return(FALSE);
         }
         break;
      default:
         break;
      }
   }  

   return(TRUE);
}


/************************************************************************/
/*>BOOL ShowClusters(FILE *fp, REAL **data, int NVec, int VecDim, 
                     int Method, BOOL ShowTable, BOOL ShowDendogram)
   -----------------------------------------------------------------
   Input:   FILE *fp              Output file pointer
            REAL **data           2D array of data to cluster
            int  NVec             Number of vectors to cluster
            int  VecDim           Dimension of the vectors
            int  Method           Clustering method (1--7)
            BOOL ShowTable        Display clustering table?
            BOOL ShowDendogram    Display clustering dendogram?
   Returns: BOOL                  Success of memory allocations (local,
                                  in clustering or in critical residue
                                  definition)

   Allocate temporary arrays, call clustering, post-clustering, 
   data-writing and critical residue code and free the temporary arrays.

   21.06.95 Original    By: ACRM
   29.06.95 Initialise all pointers to NULL
   06.07.95 Added PostCluster() call
   24.07.95 Added WriteClusData() call
   08.08.95 Added DefineCriticalResidues() call
   15.08.95 Modified to use a single array to store the clusters at the
            clustering level of ineterest rather than always looking in
            the 2D clusters array
   16.08.95 PostCluster() now returns the new number of clusters
   25.09.95 Passes VecDim to ClusterDendogram if Method==1 (else passes
            1.0)
*/
BOOL ShowClusters(FILE *fp, REAL **data, int NVec, int VecDim, 
                  int Method, BOOL ShowTable, BOOL ShowDendogram)
{
   int  *ia          = NULL,
        *ib          = NULL,
        lev          = NVec,
        *iorder      = NULL,
        *height      = NULL,
        **clusters   = NULL,
        *TheClusters = NULL,
        NClus,
        OldNClus;
   REAL *crit        = NULL,
        *critval     = NULL;
   BOOL ok           = TRUE;
   char **out        = NULL;

   ia          = (int *)malloc(NVec * sizeof(int));
   ib          = (int *)malloc(NVec * sizeof(int));
   crit        = (REAL *)malloc(NVec * sizeof(REAL));
   TheClusters = (int *)malloc(NVec * sizeof(int));

   iorder      = (int *)malloc(lev * sizeof(int));
   height      = (int *)malloc(lev * sizeof(int));
   critval     = (REAL *)malloc(lev * sizeof(REAL));

   if(ia==NULL      ||
      ib==NULL      ||
      crit==NULL    ||
      critval==NULL ||
      iorder==NULL  ||
      height==NULL  ||
      TheClusters==NULL)
   {
      if(ia          != NULL) free(ia);
      if(ib          != NULL) free(ib);
      if(iorder      != NULL) free(iorder);
      if(height      != NULL) free(height);
      if(crit        != NULL) free(crit);
      if(critval     != NULL) free(critval);
      if(TheClusters != NULL) free(TheClusters);
      return(FALSE);
   }

   /* Write a header for the clustering                                 */
   WriteHeader(fp,Method,NVec,VecDim,gScheme);

   /* Write the raw cluster data                                        */
   if(gDoData)
      WriteClusData(fp,NVec,VecDim,data);
   
   /* Do the clustering                                                 */
   if(HierClus(NVec,VecDim,Method,data,ia,ib,crit))
   {
      /* Assign data to clusters                                        */
      if((clusters = ClusterAssign((ShowTable?fp:NULL),NVec,ia,ib,crit,
                                   lev,iorder,critval,height))!=NULL)
      {
         /* Print the dendogram                                         */
         if(ShowDendogram)
         {
            if((out = ClusterDendogram(fp,lev,iorder,height,critval,
                                       ((Method==1)?VecDim:1.0)))
               ==NULL)
               ok = FALSE;
         }


         /* Find number of distinct clusters                            */
         NClus = FindNumTrueClusters(crit, lev, VecDim);

         /* Setup the TheClusters array with these values               */
         FillClusterArray(clusters, NVec, NClus, TheClusters);

         WriteResults(fp, TheClusters, NClus, data, NVec, VecDim, crit, 
                      FALSE);
         OldNClus = NClus;
         if((NClus = PostCluster(fp, TheClusters, data, NVec, VecDim, 
                                 crit, NClus))==0)
         {
            ok = FALSE;
         }
         else
         {
            WriteResults(fp, TheClusters, NClus, data, NVec, VecDim, 
                         crit, TRUE);
            if(gDoCritRes)
            {
               if(!DefineCriticalResidues(fp,TheClusters,data,NVec,VecDim,
                                          crit, NClus))
                  ok = FALSE;
            }
         }
      }
      else
      {
         ok = FALSE;
      }
   }
   else
   {
      ok = FALSE;
   }

   /* Free allocated memory                                             */
   if(clusters != NULL) FreeArray2D((char **)clusters,NVec,VecDim);
   if(out      != NULL) FreeArray2D((char **)out,lev*3,lev*3);
   if(ia       != NULL) free(ia);
   if(ib       != NULL) free(ib);
   if(iorder   != NULL) free(iorder);
   if(height   != NULL) free(height);
   if(crit     != NULL) free(crit);
   if(critval  != NULL) free(critval);
   
   return(ok);
}


/************************************************************************/
/*>BOOL HierClus(int NVec, int VecDim, int ClusterMethod, REAL **data, 
                 int *ia, int *ib, REAL *crit)
   -------------------------------------------------------------------
   Input:   int  NVec                 Number of vectors to cluster
            int  VecDim               Dimension of each vector
            REAL data[NVec][VecDim]   Input data matrix
            int  ClusterMethod        Clustering criterion to be used
   Output:  int  ia[NVec]             \
                 ib[NVec]             | History of allomerations
            REAL crit[NVec]           /
   Returns: BOOL                      Success of memory allocations

   Hierarchical clustering using user-specified criterion. 
                                                             
   20.06.95 Original By: ACRM
            Based on FORTRAN code by F. Murtagh, ESA/ESO/STECF, Garching,
            February 1986 available in STATLIB.
   26.06.95 Fixed frees on error
   29.06.95 Fix the pointers in the data array when finished
   30.01.09 Initialize some variables
*/
BOOL HierClus(int NVec, int VecDim, int ClusterMethod, REAL **data, 
              int *ia, int *ib, REAL *crit)
{
   int  ind, 
        ind1, 
        ind2, 
        ind3, 
        NClusters, 
        i, 
        j, 
        k, 
        i2, 
        j2, 
        jj = 0, 
        im = 0, 
        jm = 0,
        *NearNeighb = NULL;
   REAL DMin, 
        x, 
        xx,
        *DissimNearNeighb = NULL,
        *LDDissim = NULL,
        *membr = NULL;
   BOOL *Flag = NULL;

   /* Indicate agglomerable object/clusters                             */
   Flag  = (BOOL *)malloc(NVec * sizeof(BOOL));
   /* Current nearest neighbour storage                                 */
   NearNeighb    = (int  *)malloc(NVec * sizeof(int));
   /* Cluster cardinalities                                             */
   membr = (REAL *)malloc(NVec * sizeof(REAL));
   /* Dissimilarity of nearest neighbour                                */
   DissimNearNeighb = (REAL *)malloc(NVec * sizeof(REAL));
   /* Stores dissimilarities in lower half diagonal                     */
   LDDissim  = (REAL *)malloc(NVec*(NVec-1)/2 * sizeof(REAL));

   /* Check allocations                                                 */
   if(Flag             == NULL || 
      NearNeighb       == NULL || 
      membr            == NULL || 
      DissimNearNeighb == NULL || 
      LDDissim         == NULL)
   {
      if(Flag!=NULL)             free(Flag);
      if(NearNeighb!=NULL)       free(NearNeighb);
      if(membr!=NULL)            free(membr);
      if(DissimNearNeighb!=NULL) free(DissimNearNeighb);
      if(LDDissim!=NULL)         free(LDDissim);

      return(FALSE);
   }

   /* For all arrays, move pointer back one so we count FORTRAN-style
      from 1 rather than from 0
   */
   Flag--;                  /* Local arrays                             */
   NearNeighb--;
   membr--;
   DissimNearNeighb--;
   LDDissim--;

   crit--;                  /* Passed parameter arrays                  */
   ib--;
   ia--;
   for(i=0; i<NVec; i++)
      (data[i])--;
   data--;
   
   /* Initializations                                                   */
   for(i=1; i<=NVec; i++) 
   {
      membr[i] = (REAL)1.0;
      Flag[i]  = TRUE;
   }
   NClusters = NVec;
   
   /* Construct dissimilarity matrix                                    */
   for(i=1; i<=NVec-1; i++) 
   {
      for(j=i+1; j<=NVec; j++) 
      {
         ind = IOFFSET(NVec, i, j);
         LDDissim[ind] = (REAL)0.0;
         for(k=1; k<=VecDim; k++) 
         {
            LDDissim[ind] += (data[i][k] - data[j][k]) * 
                             (data[i][k] - data[j][k]);
         }
         
         /* For the case of the min. var. method where merging criteria 
            are defined in terms of variances rather than distances. 
         */
         if (ClusterMethod == 1) 
         {
            LDDissim[ind] /= (REAL)2.0;
         }
      }
   }
   
   /* Carry out an agglomeration - first create list of near neighbours */
   for(i=1; i<=NVec-1; i++) 
   {
      DMin = INF;
      for(j=i+1; j<=NVec; j++) 
      {
         ind = IOFFSET(NVec, i, j);
         if (LDDissim[ind] < DMin) 
         {
            DMin = LDDissim[ind];
            jm = j;
         }
      }
      NearNeighb[i] = jm;
      DissimNearNeighb[i] = DMin;
   }
   
   /* Next, determine least dissimilar using list of near neighbours    */
   do
   {
      DMin = INF;
      for(i=1; i<=NVec-1; i++) 
      {
         if(Flag[i] && (DissimNearNeighb[i] < DMin))
         {
            DMin = DissimNearNeighb[i];
            im   = i;
            jm   = NearNeighb[i];
         }
      }
      NClusters--;
      
      /* This allows an agglomeration to be carried out                 */
      i2 = MIN(im,jm);
      j2 = MAX(im,jm);
      ia[NVec   - NClusters] = i2;
      ib[NVec   - NClusters] = j2;
      crit[NVec - NClusters] = DMin;
      
      /* Update dissimilarities from new cluster                        */
      Flag[j2] = FALSE;
      DMin     = INF;
      for(k=1; k<=NVec-1; k++) 
      {
         if(Flag[k] && (k != i2))
         {
            x = membr[i2] + membr[j2] + membr[k];
            
            if (i2 < k) 
               ind1 = IOFFSET(NVec, i2, k);
            else 
               ind1 = IOFFSET(NVec, k, i2);
            
            if (j2 < k)
               ind2 = IOFFSET(NVec, j2, k);
            else
               ind2 = IOFFSET(NVec, k, j2);
            
            ind3 = IOFFSET(NVec, i2, j2);
            xx   = LDDissim[ind3];
            
            switch(ClusterMethod)
            {
            case 1:
               /*  Ward's minimum variance method                       */
               LDDissim[ind1] = (membr[i2] + membr[k]) * LDDissim[ind1] + 
                                (membr[j2] + membr[k]) * LDDissim[ind2] -
                                membr[k] * xx;
               LDDissim[ind1] /= x;
               break;
            case 2:
               /*  Single link method                                   */
               LDDissim[ind1] = MIN(LDDissim[ind1], LDDissim[ind2]);
               break;
            case 3:
               /*  Complete link method                                 */
               LDDissim[ind1] = MAX(LDDissim[ind1], LDDissim[ind2]);
               break;
            case 4:
               /*  Average link (or group average) method               */
               LDDissim[ind1] = (membr[i2] * LDDissim[ind1] + 
                                 membr[j2] * LDDissim[ind2]) / 
                                (membr[i2] + membr[j2]);
               break;
            case 5:
               /*  McQuitty's method                                    */
               LDDissim[ind1] = LDDissim[ind1] * (REAL).5 + 
                                LDDissim[ind2] * (REAL).5;
               break;
            case 6:
               /*  Median (Gower's) method                              */
               LDDissim[ind1] = LDDissim[ind1] * (REAL).5 + 
                                LDDissim[ind2] * (REAL).5 - 
                                xx         * (REAL).25;
               break;
            case 7:
               /*  Centroid method                                      */
               LDDissim[ind1] = (membr[i2] * LDDissim[ind1] + 
                                 membr[j2] * LDDissim[ind2] - 
                                 membr[i2] * membr[j2]  * xx / 
                                 (membr[i2] + membr[j2])) / 
                                (membr[i2] + membr[j2]);
               break;
            }
            
            if((i2 <= k) && (LDDissim[ind1] < DMin))
            {
               DMin = LDDissim[ind1];
               jj = k;
            }
         }
      }
      
      membr[i2] += membr[j2];
      DissimNearNeighb[i2] = DMin;
      NearNeighb[i2] = jj;
      
      /* Update list of nearest neighbours as required.                 */
      for(i=1; i<=NVec-1; i++) 
      {
         if(Flag[i]) 
         {
            if(NearNeighb[i]==i2 || NearNeighb[i]==j2) 
            {
               /* Redetermine nearest neighbour of I                    */
               DMin = INF;
               for(j=i+1; j<=NVec; j++) 
               {
                  ind = IOFFSET(NVec, i, j);
                  if(Flag[j] && (i!=j) && (LDDissim[ind] < DMin))
                  {
                     DMin = LDDissim[ind];
                     jj = j;
                  }
               }
               NearNeighb[i] = jj;
               DissimNearNeighb[i] = DMin;
            }
         }
      }
   }  while(NClusters>1);
   /* Repeat previous steps until N-1 agglomerations carried out.       */

   free(++DissimNearNeighb);
   free(++LDDissim);
   free(++membr);
   free(++NearNeighb);
   free(++Flag);

   data++;
   for(i=0; i<NVec; i++)
      (data[i])++;

   return(TRUE);
}


/************************************************************************/
/*>int **ClusterAssign(FILE *fp, int NVec, int *ia, int *ib, REAL *crit, 
                       int lev, int *iorder, REAL *critval, int *height)
   ---------------------------------------------------------------------
   Input:   FILE *fp          File for output (or NULL)
            int  NVec         Number of vectors
            int  ia[NVec]     \
                 ib[NVec]     | History of allomerations
            REAL crit[NVec]   /
            int  lev          Number of clusters in largest partition
   Output:  int  iorder[lev]  \
            REAL critval[lev] | vectors describing the dendrogram
            int  height[lev]  /
   Returns: int  **           NVec X lev array describing clusters

   Given a HIERARCHIC CLUSTERING, described as a sequence of    
   agglomerations, derive the assignments into clusters for the 
   top LEV-1 levels of the hierarchy.                           
   Prepare also the required data for representing the          
   dendrogram of this top part of the hierarchy.                

   Pick out the clusters which the N objects belong to, 
   at levels N-2, N-3, ... N-LEV+1 of the hierarchy. 
   The clusters are identified by the lowest seq. no. of 
   their members. 
   There are 2, 3, ... LEV clusters, respectively, for the 
   above levels of the hierarchy. 

   20.06.95 Original By: ACRM
            Based on FORTRAN code by F. Murtagh, ESA/ESO/STECF, Garching,
            February 1986 available in STATLIB.     
   22.06.95 Size of hvals array should be [lev+2], not [lev]
   26.06.95 Fixed potential off-bottom of array accesses
   24.07.95 Added code to print loop identifiers
*/
int **ClusterAssign(FILE *fp, int NVec, int *ia, int *ib, REAL *crit, 
                    int lev, int *iorder, REAL *critval, int *height)
{
   int      ilev, 
            i, 
            j, 
            k, 
            level, 
            icl, 
            NClusters, 
            loc,
            *hvals = NULL,
            **clusters = NULL;
   BOOL     BreakOut;
   DATALIST *p;
   char     *loopid;

   /* Allocate memory                                                   */
   clusters = (int **)Array2D(sizeof(int),NVec,lev);
   hvals  = (int *)malloc((2+lev) * sizeof(int));

   /* Check allocations; free and return if failed                      */
   if(clusters==NULL || hvals==NULL)
   {
      if(clusters!=NULL)
         FreeArray2D((char **)clusters, NVec, lev);
      if(hvals!=NULL)
         free(hvals);
      return(NULL);
   }

   /* For all arrays, move pointer back one so we count FORTRAN-style
      from 1 rather than from 0
   */
   crit--;
   ia--;
   ib--;
   height--;
   critval--;
   iorder--;
   hvals--;
   for(i=0; i<NVec; i++)
      (clusters[i])--;
   clusters--;


   hvals[1] = 1;
   hvals[2] = ib[NVec - 1];
   loc = 3;
   for(i=NVec-2; i>=NVec-lev && i>0; i--) 
   {
      for(j=1,BreakOut=FALSE; j<=loc-1; j++) 
      {
         if(ia[i] == hvals[j]) 
         {
            BreakOut=TRUE;
            break;
         }
      }
      if(!BreakOut)
      {
         hvals[loc] = ia[i];
         loc++;
      }
      
      for(j=1,BreakOut=FALSE; j<=loc-1; j++) 
      {
         if(ib[i]==hvals[j]) 
         {
            BreakOut=TRUE;
            break;
         }
      }
      if(!BreakOut)
      {
         hvals[loc] = ib[i];
         loc++;
      }
   }

   for(level=NVec-lev; level<=NVec-2; level++) 
   {
      for(i=1; i<=NVec; i++) 
      {
         icl = i;
         for(ilev=1; ilev<=level; ilev++) 
         {
            if(ib[ilev] == icl) 
            {
               icl = ia[ilev];
            }
         }
         NClusters = NVec-level;
         clusters[i][NClusters-1] = icl;
      }
   }

   for(i=1; i<=NVec; i++) 
   {
      for(j=1; j<=lev-1; j++) 
      {
         for(k=2; k<=lev; k++) 
         {
            if(clusters[i][j] == hvals[k]) 
            {
               clusters[i][j] = k;
               break;
            }
         }
      }
   }

   if(fp!=NULL)
   {
      fprintf(fp, "\nBEGIN CLUSTABLE\n");
      
      fprintf(fp,"     SEQ NOS 2CL 3CL 4CL 5CL 6CL 7CL 8CL 9CL\n");
      fprintf(fp,"     ------- --- --- --- --- --- --- --- --- ----\n");
      for(i=1,p=gDataList; i<=NVec && p!=NULL; i++, NEXT(p))
      {
         loopid = FNam2PDB(p->loopid);
         
         if(loopid==NULL)
            fprintf(fp,"%11d",i);
         else
            fprintf(fp,"   %4s%4d",loopid,i);

         for(j=1; j<=lev-1; j++)
         {
            fprintf(fp,"%4d",clusters[i][j]);
         }
         fprintf(fp,"\n");
      }
      fprintf(fp, "END CLUSTABLE\n");
   }
   
   /* Determine an ordering of the LEV clusters (at level LEV-1) 
      for later representation of the dendrogram. 
      These are stored in IORDER. 
      Determine the associated ordering of the criterion values 
      for the vertical lines in the dendrogram. 
      The ordinal values of these criterion values may be used in 
      preference, and these are stored in HEIGHT. 
      Finally, note that the LEV clusters are renamed so that they 
      have seq. nos. 1 to LEV.
   */ 
   iorder[1] = ia[NVec - 1];
   iorder[2] = ib[NVec - 1];
   critval[1] = (REAL)0.;
   critval[2] = crit[NVec - 1];
   height[1] = lev;
   height[2] = lev - 1;
   loc = 2;
   for(i=NVec-2; i>=NVec-lev+1; i--) 
   {
      for(j=1; j<=loc; j++) 
      {
         if(ia[i] == iorder[j]) 
         {
            /* Shift rightwards and insert IB(I) beside IORDER(J)       */
            for(k=loc+1; k>=j+1; k--) 
            {
               iorder[k]  = iorder[k-1];
               critval[k] = critval[k-1];
               height[k]  = height[k-1];
            }
            iorder[j + 1]  = ib[i];
            critval[j + 1] = crit[i];
            height[j + 1]  = i - (NVec - lev);
            loc++;
         }
      }
   }

   for(i=1; i<=lev; i++) 
   {
      for(j=1; j<=lev; j++) 
      {
         if(hvals[i]==iorder[j]) 
         {
            iorder[j] = i;
            break;
         }
      }
   }


   /* Fix iorder[] array to give the correct numbers along the bottom   */
   iorder[1] = 1;
   iorder[2] = 2;
   
   for(j=2; j<=lev-1; j++)
   {
      int parent;
      
      for(i=1; i<=NVec; i++)
      {
         if(clusters[i][j] == j+1)
         {
            /* This is the new cluster; see what it's parent was        */
            parent = clusters[i][j-1];

            /* Insert this cluster number to the right of the parent in 
               iorder
            */
            InsertIorder(iorder,lev,j+1,parent);
            
            break;
         }
      }
      
   }

   /* Fix pointers into arrays                                          */
   clusters++;
   for(i=0; i<NVec; i++)
      (clusters[i])++;

   /* Free temporary storage                                            */
   free(++hvals);

   /* Return cluster array                                              */
   return(clusters);
}


/************************************************************************/
/*>char **ClusterDendogram(FILE *fp, int lev, int *iorder, int *height, 
                           REAL *critval, REAL DivFactor) 
   ----------------------------------------------------------------------
   Input:   FILE *fp          Output file pointer (or NULL)
            int  lev          Number of clustering levels to display
            int  iorder[lev]  Ordering of objects along the bottom of the
                              dendrogram     
            int  height[lev]  Height of the vertical above each object, 
                              in ordinal values   
            REAL critval[lev] Height in real values
            REAL DivFactor    Value to divide critical values by before
                              displaying them
   Returns: char **           The dendogram [lev*3][lev*3]
                                                  
   Construct a dendrogram of the top lev levels of a hierarchic 
   clustering.                       

   20.06.95 Original By: ACRM
            Based on FORTRAN code by F. Murtagh, ESA/ESO/STECF, Garching,
            February 1986 available in STATLIB.     
   25.09.95 Divides critical values by DivFactor before display
*/
char **ClusterDendogram(FILE *fp, int lev, int *iorder, int *height, 
                        REAL *critval, REAL DivFactor) 
{
   int  i,  j,  k,  l,
        i2, j2, i3, 
        ic, idum;
   char **out;
   
   if((out = (char **)Array2D(sizeof(char),lev*3,lev*3))==NULL)
   {
      return(NULL);
   }
   
   /* Blank the dendogram                                               */
   for(i=0; i<lev*3; i++)
   {
      for(j=0; j<lev*3; j++)
      {
         out[i][j]=BLANK;
      }
   }
   
   /* Build the dendogram                                               */
   for(i=3; i<=lev*3; i+=3)
   {
      i2=i/3;
      
      j2=(lev*3+1)-3*height[i2-1];
      for(j=lev*3; j>=j2; j--)
      {
         out[j-1][i-1]=UP;
      }
      
      for(k=i; k>=3; k--)
      {
         i3=(int)((k+2)/3);
         if(((lev*3+1)-height[i3-1]*3)<j2)
            break;
         out[j2-1][k-1]=ACROSS;
      }
   }
   
   /* Print the dendogram                                               */
   if(fp != NULL)
   {
      fprintf(fp, "\nBEGIN DENDOGRAM\n");
      
      ic=3;
      for(i=1; i<=lev*3; i++)
      {
         if(i==ic+1)
         {
            idum=ic/3;
            idum=lev-idum;
            for(l=1; l<=lev; l++)
            {
               if(height[l-1]==idum)
                  break;
            }
            idum=l;
            
            fprintf(fp,"         %12.2g    ",critval[idum-1]/DivFactor);
            for(j=0; j<lev*3; j++)
            {
               fprintf(fp,"%c",out[i-1][j]);
            }
            fprintf(fp,"\n");
            
            ic+=3;
         }
         else
         {
            fprintf(fp,"                         ");
            for(j=0; j<lev*3; j++)
            {
               fprintf(fp,"%c",out[i-1][j]);
            }
            fprintf(fp,"\n");
         }
      }
      
      fprintf(fp,"\n                         ");
      for(i=0; i<lev; i++)
         fprintf(fp,"%3d",iorder[i]);
      fprintf(fp,"\n\n");
      
      fprintf(fp,"              CRITERION        CLUSTERS 1 TO LEV\n");
      fprintf(fp,"              VALUES.      (TOP LEV-1 LEVELS OF \
HIERARCHY).\n");
      fprintf(fp, "END DENDOGRAM\n");
   }

   return(out);
}


/************************************************************************/
/*>BOOL DoClustering(BOOL CATorsions)
   ----------------------------------
   Input:   BOOL      CATorsions    Do CA rather than true torsions?
   Returns: BOOL                    Success of memory allocations.
   Globals: DATALIST *gDataList     Linked list of data
            FILE     *gOutfp        Output file pointer
            int      gMaxLoopLen    Maximum loop length
            int      gClusterMethod Clustering method
            BOOL     gDoTable       Print the cluster table?
            BOOL     gDoDendogram   Print the dendogram?

   Calls ConvertData() which converts data stored as a linked list into 
   a 2D array and calls the ShowClusters() clustering code. Finally frees
   the memory allocated for the 2D array.

   27.06.95 Original   By: ACRM
   13.09.95 Corrected dimensionality of vectors (was 2*gMaxLoopLen: 
            wrong when doing true torsions) and account for distances.
   21.09.95 Modified calculation of VecDim
*/
BOOL DoClustering(BOOL CATorsions)
{
   REAL **data;
   int  NData,
        VecDim;
   BOOL retval;

   VecDim = 2;
   if(!CATorsions) VecDim += 4;
   if(gDoAngles)   VecDim += 1;
   if(gDoDistance) VecDim += 1;
   VecDim *= gMaxLoopLen;

   if((data = ConvertData(gDataList, &NData, CATorsions))==NULL)
      return(FALSE);

   retval = ShowClusters(gOutfp, data, NData, VecDim, 
                         gClusterMethod, gDoTable, gDoDendogram);
   FreeArray2D((char **)data, NData, VecDim);

   return(retval);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   27.06.95 Original   By: ACRM
   05.07.95 V1.1
   06.07.95 V1.2
   07.07.95 V1.3
   24.07.95 V1.4
   25.07.95 V2.0
   10.10.95 V3.4
   06.11.95 V3.5
   09.01.96 V3.6
*/
void Usage(void)
{
   fprintf(stderr,"\nCLAN V3.6 (c) 1995, Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: clan [-t] <datafile>\n");
   fprintf(stderr,"       -t Do true torsions\n");

   fprintf(stderr,"\nCLAN (CLuster ANalysis of Loops) performs cluster \
analysis to examine\n");
   fprintf(stderr,"loops in proteins. See the documentation for details \
of the data file\n");
   fprintf(stderr,"format.\n\n");
}


/************************************************************************/
/*>void CreateDefaultScheme(int maxres)
   ------------------------------------
   Input:   int   maxres         Number of residues
   Globals: int   gScheme[]      Scheme array

   Generates a default scheme for where insertions should be placed
   (i.e. in the middle of the loop)

   27.06.95 Original   By: ACRM
*/
void CreateDefaultScheme(int maxres)
{
   int i, j;
   
   /* Number forwards from 1, stepping by 2                             */
   for(i=1, j=0; 
       i<=maxres; 
       i+=2, j++)
   {
      gScheme[j] = i;
   }

   /* Number backwards from 2, stepping by 2                            */
   for(i=2, j=maxres-1;
       i<=maxres;
       i+=2, j--)
   {
      gScheme[j] = i;
   }
}


/************************************************************************/
/*>void InsertIorder(int *iorder, int lev, int cluster, int parent)
   ----------------------------------------------------------------
   I/O:     int    *iorder     Array of numbers representing the cluster
                               ordering on the final dendogram. Array
                               positions ARE NUMBERED FROM 1, not 0
   Input:   int    lev         Length of iorder array
            int    cluster     Value to be inserted
            int    parent      Value to search for in iorder to find
                               insertion point.

   Inserts the number cluster into iorder to the right of where parent
   is found. There are lev positions in iorder which is NUMBERED FROM 1

   03.07.95 Original   By: ACRM
*/
BOOL InsertIorder(int *iorder, int lev, int cluster, int parent)
{
   int i,j;
   
   iorder++;
   
   for(i=0; i<lev; i++)
   {
      if(iorder[i]==parent)
      {
         if(i==(lev-1))
            return(FALSE);
         
         for(j=lev-1; j>=i+2; j--)
            iorder[j] = iorder[j-1];
         
         iorder[i+1] = cluster;
         
         return(TRUE);
      }
   }
   return(FALSE);
}


/************************************************************************/
/*>void WriteHeader(FILE *fp, int Method, int NVec, int VecDim, 
                    int *Scheme)
   ------------------------------------------------------------
   Input:   FILE    *fp           Output file pointer
            int     Method        Clustering method (1--7)
            int     NVec          Number of vectors to cluster
            int     VecDim        Vector dimension
            int    *Scheme        Loop insertion scheme

   Writes a header containing details of the clustering method and 
   vectors.

   04.07.95 Original    By: ACRM
   13.09.95 Prints NODISTANCE if appropriate.
   21.09.95 Also print DISTANCE and now does (NO)ANGLE as well
   26.09.95 Added TRUETORSIONS/PSEUDOTORSION
   14.03.96 Added gPClusCut[2] to POSTCLUSTER
*/
void WriteHeader(FILE *fp, int Method, int NVec, int VecDim, int *Scheme)
{
   int i;
   
   fprintf(fp,"BEGIN HEADER\n");
   fprintf(fp,"   METHOD %d\n",Method);
   fprintf(fp,"   NLOOPS %d\n",NVec);
   fprintf(fp,"   POSTCLUSTER %f %f %f\n",gPClusCut[0], gPClusCut[1],
           gPClusCut[2]);
   fprintf(fp,"   MAXLENGTH %d\n",gMaxLoopLen);
   fprintf(fp,"   SCHEME ");
   for(i=0; i<gMaxLoopLen; i++)
      fprintf(fp,"%d ",Scheme[i]);
   fprintf(fp,"\n");

   fprintf(fp,"   %s\n", (gDoDistance) ? "DISTANCE" : "NODISTANCE");
   fprintf(fp,"   %s\n", (gDoAngles)   ? "ANGLES"   : "NOANGLES");
   fprintf(fp,"   %s\n", (gCATorsions) ? 
                          "PSEUDOTORSIONS" : "TRUETORSIONS");
   
   fprintf(fp,"END HEADER\n");
}


/************************************************************************/
/*>BOOL WriteResults(FILE *fp, int *clusters, int NClus, REAL **data, 
                     int NVec, int VecDim, REAL *crit, BOOL PostClus)
   ------------------------------------------------------------------
   Input:   FILE     *fp          output file pointer
            int      *clusters    Cluster array
            int      NClus        Number of clusters
            REAL     **data       Data which is clustered
            int      NVec         Number of vectors
            int      VecDim       Vector dimension
            REAL     *crit        Critical values for clustering
            BOOL     PostClus     These are post-cluster results
   Globals: DATALIST *gDataList   Linked list of raw data

   Write detailed clustering data results

   04.07.95 Original   By: ACRM
   26.07.95 Prints number of clusters in BEGIN statement
   02.08.95 Writes ASSIGNMENTS block containing the assignment for
            each loop
   15.08.95 Now takes a vector of the cluster assignments rather than the
            whole clustering matrix
   17.08.95 Added postcluster handling
   12.09.95 Added nmemb parameter to FindMedian
*/
BOOL WriteResults(FILE *fp, int *clusters, int NClus, REAL **data, 
                  int NVec, int VecDim, REAL *crit, BOOL PostClus)
{
   int      i,
            nmemb;
   DATALIST *repres;
   
   /* Write the cluster assignment for each loop                        */
   fprintf(fp,"\nBEGIN %sASSIGNMENTS\n",(PostClus?"":"RAW"));
   for(repres=gDataList,i=0; repres!=NULL; NEXT(repres),i++)
   {
      fprintf(fp,"%3d %s\n", clusters[i], repres->loopid);
   }
   fprintf(fp,"END %sASSIGNMENTS\n",(PostClus?"":"RAW"));
   

   /* For each of the clusters, find a representitive                   */
   fprintf(fp,"\nBEGIN %sMEDIANS %d\n",(PostClus?"":"RAW"), NClus);
   for(i=1; i<=NClus; i++)
   {
      if((repres = FindMedian(clusters,data,NVec,VecDim,i,&nmemb))==NULL)
      {
         fprintf(fp,"END %sMEDIANS (failed!)\n",(PostClus?"":"RAW"));
         fprintf(stderr,"FindMedian() failed\n");
         return(FALSE);
      }
      if(repres != (DATALIST *)(-1))
         fprintf(fp,"%3d %s\n",i, repres->loopid);
   }
   fprintf(fp,"END %sMEDIANS\n",(PostClus?"":"RAW"));

   return(TRUE);
}


/************************************************************************/
/*>int FindNumTrueClusters(REAL *crit, int lev, int VecDim)
   --------------------------------------------------------
   Input:   REAL   *crit    Array of critical values in clustering
            int    lev      Number of clustering levels (length of crit[]
                            array)
            int    VecDim   Number of dimensions in vector
   Returns: int             Number of truely different clusters

   Finds the number of really different clusters. 

   This is currently done by dividing the critical value by the Vector
   dimensionality and selecting values greater than 0.06 which seems
   to work well for Ward's minimum variance method.

   03.07.95 Original    By: ACRM
   25.09.95 Added VecDim parameter and modified code to use it
*/
int FindNumTrueClusters(REAL *crit, int lev, int VecDim)
{
   int i;
   
   for(i=0; i<lev-1; i++)
   {
      if(crit[i]/(REAL)VecDim > (REAL)0.06)
         return(lev-i);
   }
   
   return(1);
}


/************************************************************************/
/*>void CleanUp(void)
   ------------------
   Frees up memory allocated within the global datalist linked list, then
   frees the data list

   06.07.95 Original   By: ACRM
   08.08.95 Added freeing of truestart linked list
*/
void CleanUp(void)
{
   DATALIST *p;
   for(p=gDataList; p!=NULL; NEXT(p))
   {
      if(p->torsionpdb != NULL)
         FREELIST(p->torsionpdb, PDB);
      if(p->allatompdb != NULL)
         FREELIST(p->allatompdb, PDB);
   }
   FREELIST(gDataList, DATALIST);
}


/************************************************************************/
/*>void WriteClusData(FILE *fp, int NVec, int VecDim, REAL **data)
   ---------------------------------------------------------------
   Input:   FILE  *fp           Output file pointer
            int   NVec          Number of vectors
            int   VecDim        Dimension of vectors
            REAL  **data        Data for clustering

   Writes the raw clustering data to the output file

   24.07.95 Original    By: ACRM
   26.07.95 Corrected inner loop to count to VecDim not NVec!
*/
void WriteClusData(FILE *fp, int NVec, int VecDim, REAL **data)
{
   int i, j;
   
   fprintf(fp,"\nBEGIN DATA\n");
   
   for(i=0; i<NVec; i++)
   {
      for(j=0; j<VecDim; j++)
      {
         fprintf(fp,"%10.4f",data[i][j]);
      }
      fprintf(fp,"\n");
   }
   fprintf(fp,"END DATA\n");
}


/************************************************************************/
/*>DATALIST *FindMedian(int *clusters, REAL **data, int NVec, 
                        int VecDim, int ClusNum, int *NMemb)
   ----------------------------------------------------------
   Input:   int   *clusters     Clustring vector
            REAL  **data        Data which is clustered
            int   NVec          Number of vectors
            int   VecDim        Dimension of vectors
            int   ClusNum       The cluster number we are interested in
   Output:  int   *NMemb        Number of members of this cluster
   Returns: DATALIST *          Pointer to median for this cluster.
                                NULL on memory error
                                (-1) if no members of this cluster

   Calculates the median of cluster ClusNum out of NClus clusters. Then
   finds the vector closest to the median and returns a DATALIST
   pointer to the associated loop identifier.

   03.07.95 Original    By: ACRM
   06.07.95 Changed to return DATALIST pointer rather than the actual
            loopid string
   02.08.95 Added code to handle NClus==1. Fixed array access code
            to [NClus-2] rather than [NClus]
   15.08.95 Uses 1D cluster vector rather than the 2D matrix
   12.09.95 Added NMemb output parameter
   30.01.09 Initialize some variables
*/
DATALIST *FindMedian(int *clusters, REAL **data, int NVec, int VecDim, 
                     int ClusNum, int *NMemb)
{
   DATALIST *p;
   int      i, j,
            best = 0;
   REAL     *minval, 
            *maxval,
            *medval,
            mindist = 10000.0,
            dist;
   BOOL     Done = FALSE;

   *NMemb = 0;

   /* Allocate arrays to store min and max values in each dimension     */
   if((minval=(REAL *)malloc(VecDim*sizeof(REAL)))==NULL)
      return(NULL);
   if((maxval=(REAL *)malloc(VecDim*sizeof(REAL)))==NULL)
   {
      free(minval);
      return(NULL);
   }

   /* We just use the same storage space for the median values          */
   medval = minval;

   /* Find the min and max values in each dimension                     */
   for(i=0; i<NVec; i++)
   {
      if(clusters[i]==ClusNum)
      {
         /* On first item, just copy in the data                        */
         if((*NMemb) == 0)
         {
            for(j=0; j<VecDim; j++)
               minval[j] = maxval[j] = data[i][j];
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

         /* Increment member count                                      */
         (*NMemb)++;
      }
   }

   /* If nothing found in this cluster, return NULL                     */
   if(*NMemb == 0)
      return((DATALIST *)(-1));

   /* Now store the median values                                       */
   for(j=0; j<VecDim; j++)
      medval[j] = (minval[j] + maxval[j]) / (REAL)2.0;
      
   /* Now run through again and find which is closest to the medval     */
   for(i=0, Done=FALSE; i<NVec; i++)
   {
      if(clusters[i]==ClusNum)
      {
         if(!Done)
         {
            best    = i;
            mindist = (REAL)0.0;
            for(j=0; j<VecDim; j++)
               mindist += (data[i][j] - medval[j]) *
                          (data[i][j] - medval[j]);
            
            Done = TRUE;
         }
         else
         {
            dist = (REAL)0.0;
            for(j=0; j<VecDim; j++)
               dist += (data[i][j] - medval[j]) *
                       (data[i][j] - medval[j]);
            if(dist < mindist)
            {
               mindist = dist;
               best    = i;
            }
         }
      }
   }

   /* Now run through the global data linked list to find the `best'
      example's loop pointer
   */
   for(p=gDataList, i=0; p!=NULL && i<best; NEXT(p), i++) ;

   /* Free up the arrays                                                */
   free(minval);
   free(maxval);

   /* Return the loop pointer                                           */
   return(p);
}


/************************************************************************/
/*>void FillClusterArray(int **clusters, int NVec, int NClus,
                         int *TheClusters)
   ----------------------------------------------------------
   Fills in the `TheClusters' array with the column corresponding to NClus
   from the clusters matrix.

   15.08.95 Original    By: ACRM
*/
void FillClusterArray(int **clusters, int NVec, int NClus,
                      int *TheClusters)
{
   int i;
   
   if(NClus > 1)
   {
      for(i=0; i<NVec; i++)
      {
         TheClusters[i] = clusters[i][NClus-2];
      }
   }
   else
   {
      for(i=0; i<NVec; i++)
      {
         TheClusters[i] = 1;
      }
   }
}


/************************************************************************/
/*>REAL RmsPDB(PDB *pdb1, PDB *pdb2, int length)
   ---------------------------------------------
   Input:   PDB   *pdb1        Reference PDB linked list
            int   length       Number of residues to fit
   I/O:     PDB   *pdb2        Mobile PDB linked list
                               Note that this will be moved in space
   Returns: REAL               RMS deviation

   Returns the RMS over length residues of the two PDB linked lists.
   If length is zero, all residues will be used.

   Note that pdb2 will be moved in space at the end of this.

   06.07.95 Original    By: ACRM
   26.09.95 If FitPDB() returned an error, the lists weren't being
            reassembled
   30.01.09 Initialize some variables
*/
REAL RmsPDB(PDB *pdb1, PDB *pdb2, int length)
{
   PDB  *end1 = NULL, 
        *end2 = NULL, 
        *p;
   REAL rms = 0.0;
   BOOL ok = TRUE;
  
   if(length)
   {
      /* Terminate each PDB linked list after length residues           */
      end1 = TermPDB(pdb1, length);
      end2 = TermPDB(pdb2, length);
   }
   
   if(FitPDB(pdb1, pdb2, NULL))
   {
      rms = CalcRMSPDB(pdb1, pdb2);
   }
   else
   {
      ok = FALSE;
   }

   if(length)
   {
      /* Rejoin PDB linked lists                                        */
      p=pdb1;
      LAST(p);
      p->next = end1;
      p=pdb2;
      LAST(p);
      p->next = end2;
   }

   return((ok)?rms:(REAL)9999.0);
}


/************************************************************************/
/*>REAL RmsCAPDB(PDB *pdb1, PDB *pdb2, int length)
   -----------------------------------------------
   Input:   PDB   *pdb1        Reference PDB linked list
            int   length       Number of residues to fit
   I/O:     PDB   *pdb2        Mobile PDB linked list
                               Note that this will be moved in space
   Returns: REAL               RMS deviation

   Returns the CA-RMS over length residues of the two PDB linked lists.
   If length is zero, all residues will be used.

   Note that pdb2 will be moved in space at the end of this.

   26.09.95 Original based on RmsPDB   By: ACRM
   30.01.09 Initialize some variables
*/
REAL RmsCAPDB(PDB *pdb1, PDB *pdb2, int length)
{
   PDB  *end1 = NULL, 
        *end2 = NULL, 
        *p,
        *pdbca1 = NULL,
        *pdbca2 = NULL;
   REAL rms = (REAL)9999.0;
   int  natoms;
   BOOL ok = TRUE;
   char *sel[2];

   SELECT(sel[0], "CA  ");
   if(sel[0]==NULL)
      return((REAL)9999.0);
  
   if(length)
   {
      /* Terminate each PDB linked list after length residues           */
      end1 = TermPDB(pdb1, length);
      end2 = TermPDB(pdb2, length);
   }

   if((pdbca1 = SelectAtomsPDB(pdb1, 1, sel, &natoms))==NULL)
      ok = FALSE;
   if((pdbca2 = SelectAtomsPDB(pdb2, 1, sel, &natoms))==NULL)
      ok = FALSE;
   free(sel[0]);

   if(ok)
   {
      if(!FitPDB(pdbca1, pdbca2, NULL))
         ok = FALSE;
   }

   if(ok)
      rms = CalcRMSPDB(pdbca1, pdbca2);

   if(length)
   {
      /* Rejoin PDB linked lists                                        */
      p=pdb1;
      LAST(p);
      p->next = end1;
      p=pdb2;
      LAST(p);
      p->next = end2;
   }

   if(pdbca1 != NULL)
      FREELIST(pdbca1,PDB);
   if(pdbca2 != NULL)
      FREELIST(pdbca2,PDB);

   return((ok)?rms:(REAL)9999.000);
}


/************************************************************************/
/*>REAL MaxCADeviationPDB(PDB *pdb1, PDB *pdb2, int length)
   --------------------------------------------------------
   Input:   PDB   *pdb1        Reference PDB linked list
            int   length       Number of residues to fit
   I/O:     PDB   *pdb2        Mobile PDB linked list
                               Note that this will be moved in space
   Returns: REAL               RMS deviation

   Returns the max CA-CA deviation over length residues of the two PDB 
   linked lists.
   If length is zero, all residues will be used.

   Note that pdb2 will be moved in space at the end of this.

   06.07.95 Original    By: ACRM
   26.09.95 Modified so fitting is done on CAs only
            If FitPDB() returned an error, the lists weren't being
            reassembled
   14.03.96 Changed to use FitCaPDB() which fits the complete PDB linked
            lists on CAs.
   30.01.09 Removed redundant variables
   30.01.09 Initialize some variables
*/
REAL MaxCADeviationPDB(PDB *pdb1, PDB *pdb2, int length)
{
   PDB  *end1 = NULL, 
        *end2 = NULL, 
        *p,   *q;
   REAL dev,
        maxdev = (REAL)0.0;
   BOOL ok = TRUE;

   if(length)
   {
      /* Terminate each PDB linked list after length residues           */
      end1 = TermPDB(pdb1, length);
      end2 = TermPDB(pdb2, length);
   }
   
   if(!FitCaPDB(pdb1, pdb2, NULL))
      ok = FALSE;

   /* Check max CA deviation                                            */
   if(ok)
   {
      for(p=pdb1, q=pdb2; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->atnam,"CA  ",4))
         {
            /* Step q, until we hit a CA                                */
            while((q!=NULL) && (strncmp(q->atnam,"CA  ",4))) NEXT(q);
            
            if(q!=NULL)
            {
               dev = DISTSQ(p,q);
               if(dev > maxdev)
                  maxdev = dev;

               /* Step q on by one                                      */
               NEXT(q);
            }
            else
            {
               fprintf(stderr,"MaxCADeviationPDB(): second list \
expired!\n");
               ok=FALSE;
               break;
            }
         }
      }
   }
   
   if(length)
   {
      /* Rejoin PDB linked lists                                        */
      p=pdb1;
      LAST(p);
      p->next = end1;
      p=pdb2;
      LAST(p);
      p->next = end2;
   }

   return((ok)?(REAL)sqrt((double)maxdev):(REAL)9999.0);
}


/************************************************************************/
/*>REAL MaxCBDeviationPDB(PDB *pdb1, PDB *pdb2, int length)
   --------------------------------------------------------
   Input:   PDB   *pdb1        Reference PDB linked list
            int   length       Number of residues to fit
   I/O:     PDB   *pdb2        Mobile PDB linked list
                               Note that this will be moved in space
   Returns: REAL               RMS deviation

   Returns the max CB-CB deviation over length residues of the two PDB 
   linked lists.
   If length is zero, all residues will be used.

   Note that pdb2 will be moved in space at the end of this.

   14.03.96 Original based on MaxCADeviationPDB()   By: ACRM
   30.01.09 Removed redundant variables
   30.01.09 Initialize some variables
*/
REAL MaxCBDeviationPDB(PDB *pdb1, PDB *pdb2, int length)
{
   PDB  *end1 = NULL, 
        *end2 = NULL, 
        *p,   *q,
        *pcb, *qcb;
   REAL dev,
        maxdev = (REAL)0.0;
   BOOL ok = TRUE;

   if(length)
   {
      /* Terminate each PDB linked list after length residues           */
      end1 = TermPDB(pdb1, length);
      end2 = TermPDB(pdb2, length);
   }
   
   if(!FitCaPDB(pdb1, pdb2, NULL))
      ok = FALSE;

   /* Check max CB deviation                                            */
   if(ok)
   {
      for(p=pdb1, q=pdb2; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->atnam,"N   ",4))
         {
            /* Step q, until we hit a N                                 */
            while((q!=NULL) && (strncmp(q->atnam,"N   ",4))) NEXT(q);
            
            if(q!=NULL)
            {
               /* Now check CB                                          */
               if(strncmp(p->resnam,"GLY ",4) && 
                  strncmp(q->resnam,"GLY ",4))
               {
                  pcb=FindAtomInRes(p,"CB  ");
                  qcb=FindAtomInRes(q,"CB  ");

                  if(pcb!=NULL && qcb!=NULL)
                  {
                     dev = DISTSQ(pcb,qcb);
                     if(dev > maxdev)
                        maxdev = dev;
                  }
               }
               
               /* Step q on by one                                      */
               NEXT(q);
            }
            else
            {
               fprintf(stderr,"MaxCBDeviationPDB(): second list \
expired!\n");
               ok=FALSE;
               break;
            }
         }
      }
   }
   
   if(length)
   {
      /* Rejoin PDB linked lists                                        */
      p=pdb1;
      LAST(p);
      p->next = end1;
      p=pdb2;
      LAST(p);
      p->next = end2;
   }

   return((ok)?(REAL)sqrt((double)maxdev):(REAL)9999.0);
}


/************************************************************************/
/*>int RenumClusters(int *clusters, int NVec)
   ------------------------------------------
   Renumber the clusters from 1.
   Returns the number of clusters or 0 if allocations failed.

   16.08.95 Original    By: ACRM
*/
int RenumClusters(int *clusters, int NVec)
{
   int  i, j,
        ClusNum;
   BOOL *Flags;

   /* Allocate and blank a flags array                                  */
   if((Flags=(BOOL *)malloc(NVec * sizeof(BOOL)))==NULL)
      return(0);
   for(i=0; i<NVec; i++)
      Flags[i] = FALSE;

   /* Set all used cluster number flags to TRUE                         */
   for(i=0; i<NVec; i++)
   {
      Flags[clusters[i]-1] = TRUE;
   }

   /* Run through the Flags array looking for used cluster numbers and
      renumber. This only works because renumbering will only ever
      reduce a cluster number not increase it.
   */
   for(i=0, ClusNum=0; i<NVec; i++)
   {
      if(Flags[i])
      {
         ClusNum++;
         for(j=0; j<NVec; j++)
         {
            if(clusters[j] == i+1)
            {
               clusters[j] = ClusNum;
            }
         }
      }
   }

   free(Flags);
   return(ClusNum);
}


/************************************************************************/
/*>int PostCluster(FILE *fp, int *clusters, REAL **data, int NVec, 
                   int VecDim, REAL *crit, int NClus)
   ---------------------------------------------------------------
   Input:   FILE   *fp          Output file pointer
            int    *clusters    Cluster vector
            REAL   **data       Data used for clustering
            int    NVec         Number of vectors
            int    VecDim       Dimension of vectors
            REAL   *crit        Critical values
            int    NClus        Number of clusters
   Returns: int                 Revised number of clusters   
                                (0 if memory allocations failed)
   Globals: REAL   gPClusCut    Clustering RMS cutoffs

   Show results of PostClustering on RMS deviation. This simply does LSQ
   fits of representatives from each cluster to look for those with low
   RMS deviation and merges these clusters.

   06.07.95 Original   By: ACRM
   15.08.95 Changed to use cluster vector rather than full matrix. No
            longer requires the member array or the FindMultiMedian() code
   16.08.95 Now returns revised number of clusters and renumbers the 
            clusters
   18.08.95 RenumClusters() now returns the new number of clusters
            (NClus-NMerge) did not give the correct number because of
            multiple merges into a single cluster
   12.09.95 If there are only 2 members in a cluster, both must be within
            the cutoffs in order to do a merge
   15.09.95 Removed unused variables
*/
int PostCluster(FILE *fp, int *clusters, REAL **data, int NVec, 
                int VecDim, REAL *crit, int NClus)
{
   int      NMerge = 0,
            NewNClus = NClus,
            *NMembers,
            *NewNumbers,
            i, j;
   DATALIST **repres,
            *rep_i,
            *rep_j,
            *rep_k,
            *rep_l;
   REAL     rms, rms1, rms2, rms3, rms4,
            CADev, CADev1, CADev2, CADev3, CADev4,
            CBDev, CBDev1, CBDev2, CBDev3, CBDev4;
   
   /* Allocate memory for array of representitives                      */
   if((repres = (DATALIST **)malloc(NClus * sizeof(DATALIST *)))==NULL)
      return(0);

   /* Allocate memory for new cluster numbers                           */
   if((NewNumbers = (int *)malloc(NClus * sizeof(int)))==NULL)
   {
      free(repres);
      return(0);
   }

   /* Allocate memory to store number of members of each cluster        */
   if((NMembers = (int *)malloc(NClus * sizeof(int)))==NULL)
   {
      free(repres);
      free(NewNumbers);
      return(0);
   }
   
   /* For each of these clusters, find a representitive and initialise
      the NewNumbers array
   */
   for(i=1; i<=NClus; i++)
   {
      if((repres[i-1]=FindMedian(clusters,data,NVec,VecDim,i,
                                 &(NMembers[i-1])))==NULL)
      {
         free(repres);
         free(NewNumbers);
         return(0);
      }
      NewNumbers[i-1] = i;
   }

   fprintf(fp,"\nBEGIN POSTCLUSTER\n");

   /* Compare each representitive against each other and if any matches 
      with RMS < gPClusCut[0] and max CA deviation < gPClusCut[1]
      and max CB deviation < gPClusCut[2] merge the clusters, printing 
      a message.
   */
   for(i=0; i<NClus-1; i++)
   {
      for(j=i+1; j<NClus; j++)
      {
         if((NMembers[i] != 2) && (NMembers[j] != 2))
         {
            if(TestMerge(repres[i], repres[j], &rms, &CADev, &CBDev))
            {
               NMerge++;
               DoMerge(fp,i,repres[i],j,repres[j],rms,CADev,CBDev,
                       NewNumbers,NClus);
            }
         }
         else if((NMembers[i] == 2) && (NMembers[j] != 2))
         {
            rep_i = FindLoop(clusters,NVec,i+1,0);
            rep_j = FindLoop(clusters,NVec,i+1,1);
            
            if((rep_i != NULL) && (rep_j != NULL))
            {
               if(TestMerge(rep_i, repres[j], &rms1, &CADev1, &CBDev1) &&
                  TestMerge(rep_j, repres[j], &rms2, &CADev2, &CBDev2))
               {
                  NMerge++;
                  rms = (rms1 + rms2) / (REAL)2.0;
                  CADev = (CADev1 + CADev2) / (REAL)2.0;               
                  CBDev = (CBDev1 + CBDev2) / (REAL)2.0;               
                  DoMerge(fp,i,repres[i],j,repres[j],rms,CADev,CBDev,
                          NewNumbers,NClus);
               }
            }
            else
            {
               if(rep_i == NULL)
                  fprintf(stderr,"INTERR: Loop 0 not found in cluster \
%d\n",i);
               if(rep_j == NULL)
                  fprintf(stderr,"INTERR: Loop 1 not found in cluster \
%d\n",i);
            }
         }
         else if((NMembers[i] != 2) && (NMembers[j] == 2))
         {
            rep_i = FindLoop(clusters,NVec,j+1,0);
            rep_j = FindLoop(clusters,NVec,j+1,1);
            
            if((rep_i != NULL) && (rep_j != NULL))
            {
               if(TestMerge(repres[i], rep_i, &rms1, &CADev1, &CBDev1) &&
                  TestMerge(repres[i], rep_j, &rms2, &CADev2, &CBDev2))
               {               
                  NMerge++;
                  rms = (rms1 + rms2) / (REAL)2.0;
                  CADev = (CADev1 + CADev2) / (REAL)2.0;
                  CBDev = (CBDev1 + CBDev2) / (REAL)2.0;               
                  DoMerge(fp,i,repres[i],j,repres[j],rms,CADev,CBDev,
                          NewNumbers,NClus);
               }
            }
            else
            {
               if(rep_i == NULL)
                  fprintf(stderr,"INTERR: Loop 0 not found in cluster \
%d\n",j);
               if(rep_j == NULL)
                  fprintf(stderr,"INTERR: Loop 1 not found in cluster \
%d\n",j);
            }
         }
         else /* Both clusters have 2 members                           */
         {
            rep_i = FindLoop(clusters,NVec,i+1,0);
            rep_j = FindLoop(clusters,NVec,i+1,1);
            rep_k = FindLoop(clusters,NVec,j+1,0);
            rep_l = FindLoop(clusters,NVec,j+1,1);

            if((rep_i != NULL) && (rep_j != NULL) && 
               (rep_k != NULL) && (rep_l != NULL))
            {
               if(TestMerge(rep_i, rep_k, &rms1, &CADev1, &CBDev1) &&
                  TestMerge(rep_i, rep_l, &rms2, &CADev2, &CBDev2) &&
                  TestMerge(rep_j, rep_l, &rms3, &CADev3, &CBDev3) &&
                  TestMerge(rep_j, rep_k, &rms4, &CADev4, &CBDev4))
               {
                  NMerge++;
                  rms = (rms1 + rms2 + rms3 + rms4) / (REAL)4.0;
                  CADev = (CADev1 + CADev2 + CADev3 + CADev4) / (REAL)4.0;
                  CBDev = (CBDev1 + CBDev2 + CBDev3 + CBDev4) / (REAL)4.0;
                  DoMerge(fp,i,repres[i],j,repres[j],rms,CADev,CBDev,
                          NewNumbers,NClus);
               }
            }
            else
            {
               if(rep_i == NULL)
                  fprintf(stderr,"INTERR: Loop 0 not found in cluster \
%d\n",i);
               if(rep_j == NULL)
                  fprintf(stderr,"INTERR: Loop 1 not found in cluster \
%d\n",i);
               if(rep_k == NULL)
                  fprintf(stderr,"INTERR: Loop 0 not found in cluster \
%d\n",j);
               if(rep_l == NULL)
                  fprintf(stderr,"INTERR: Loop 1 not found in cluster \
%d\n",j);
            }
         }
      }
   }
   
   /* If some merging has occured, update the actual clusters array and
      print message
   */
   if(NMerge)
   {
      for(i=0; i<NClus; i++)
      {
         for(j=0; j<NVec; j++)
         {
            if(clusters[j] == (i+1))
               clusters[j] = NewNumbers[i];
         }
      }

      if((NewNClus=RenumClusters(clusters, NVec))==0)
      {
         fprintf(stderr,"Warning: Cluster renumbering out of memory\n");
         fprintf(stderr,"         Number of clusters not corrected, so \
expect strange results!\n");
      }
   }
   else
   {
      fprintf(fp,"No merges were performed\n");
   }

   fprintf(fp,"END POSTCLUSTER\n");

   free(repres);
   free(NewNumbers);
   
   return(NewNClus);
}


/************************************************************************/
/*>BOOL TestMerge(DATALIST *loop1, DATALIST *loop2, REAL *rms, 
                  REAL *CADev, REAL *CBDev)
   -----------------------------------------------------------
   Input:   DATALIST *loop1         First loop to test
            DATALIST *loop2         Second loop to test
   Output:  REAL     *rms           RMS between loops
            REAL     *CADev         Max CA deviation between loops
            REAL     *CBDev         Max CB deviation between loops
   Returns: BOOL                    Should they be merged?

   Tests whether two loop examples should be merged into one cluster

   12.09.95 Original    By: ACRM
   26.09.95 Changed to use RmsCAPDB()
            Changed comparison from < to <= 
   14.03.96 Checks CB deviation as well as CA. Required code to duplicate
            the main PDB linked list.
            Added CBDev parameter
   15.04.96 This causes a problem if critical residues have not been
            requested, as these data are not available. Prompts with
            a message in this case. Only gives one warning message.
*/
BOOL TestMerge(DATALIST *loop1, DATALIST *loop2, REAL *rms, REAL *CADev,
               REAL *CBDev)
{
   static BOOL Warned = FALSE;
   PDB         *dupe1 = NULL,
               *dupe2 = NULL,
               *p, *end;
   BOOL        DupeDone = TRUE;
   
   if((loop1 != (DATALIST *)(-1)) && (loop2 != (DATALIST *)(-1)))
   {
      /* Only bother trying to merge clusters if loops are of the same 
         length.
      */
      if(loop1->length == loop2->length)
      {

         /* If less than cutoff, merge clusters. Note that this must
            be done on a copy of the all-atom linked list, otherwise
            the search for critical residues will fail since these
            routines move the loops in space 
         */
         if((p=FindResidueSpec(loop1->allatompdb,loop1->start))!=NULL)
         {
            /* Duplicate from the start res on                          */
            if((dupe1 = DupePDB(p))!=NULL)
            {
               if((end = TermPDB(dupe1, loop1->length + 1))!=NULL)
                  FREELIST(end,PDB);
            }
         }
         if((p=FindResidueSpec(loop2->allatompdb,loop2->start))!=NULL)
         {
            /* Duplicate from the start res on                          */
            if((dupe2 = DupePDB(p))!=NULL)
            {
               if((end = TermPDB(dupe2, loop2->length + 1))!=NULL)
                  FREELIST(end,PDB);
            }
         }
         
         if(dupe1==NULL || dupe2==NULL)
         {
            if(dupe1!=NULL) FREELIST(dupe1, PDB);
            if(dupe2!=NULL) FREELIST(dupe2, PDB);
            if(!Warned)
            {
               fprintf(stderr,"Warning: Unable to duplicate PDB linked \
lists.\n");
               fprintf(stderr,"         Max deviations in merging will \
only be done on CA, not CB\n");

               if(loop1->allatompdb == NULL || loop2->allatompdb == NULL)
               {
                  fprintf(stderr,"         You can solve this by using \
the CRITICAL keyword *before* the LOOP specifications.\n");
               }
               Warned = TRUE;
            }
            
            dupe1    = loop1->pdbloop;
            dupe2    = loop2->pdbloop;
            DupeDone = FALSE;
         }
            
         *rms = RmsCAPDB(dupe1, 
                         dupe2, 
                         loop1->length);
         *CADev = MaxCADeviationPDB(dupe1, 
                                    dupe2,
                                    loop1->length);
         *CBDev = MaxCBDeviationPDB(dupe1, 
                                    dupe2,
                                    loop1->length);

         /* If we managed to duplicate our PDB linked lists, free them  */
         if(DupeDone)
         {
            FREELIST(dupe1, PDB);
            FREELIST(dupe2, PDB);
         }

         if(sInfoLevel)
         {
            fprintf(stderr,"Test %s with %s. RMS=%.3f MAXCA=%.3f \
MAXCB=%.3f\n", loop1->loopid, loop2->loopid, *rms, *CADev, *CBDev);
         }

         if(((gPClusCut[0] == 0.0) || (*rms   <= gPClusCut[0])) &&
            ((gPClusCut[1] == 0.0) || (*CADev <= gPClusCut[1])) &&
            ((gPClusCut[2] == 0.0) || (*CBDev <= gPClusCut[2])))
         {
            return(TRUE);
         }
      }
   }
   return(FALSE);
}


/************************************************************************/
/*>void DoMerge(FILE *fp, int i, DATALIST *loop1, int j, DATALIST *loop2,
                REAL rms, REAL CADev, REAL CBDev, int *NewNumbers, 
                int NClus)
   ----------------------------------------------------------------------
   Input:   FILE     *fp          File to write information to
            int      i            First cluster number
            DATALIST *loop1       Representative of first cluster
            int      j            Second cluster number
            DATALIST *loop2       Representative of second cluster
            REAL     rms          RMS between clusters
            REAL     CADev        Max CA deviation between clusters
            REAL     CBDev        Max CB deviation between clusters
            int      NClus        Initial number of clusters
   I/O:     int      *NewNumbers  Array of new cluster numbers  
           
   Actually merges two clusters. The high numbered cluster will always
   be given the number of the low numbered cluster.

   The input i,j are the raw cluster numbers -1 and are also offsets
   into the cluster number array

   12.09.95 Original    By: ACRM
   14.03.96 Added CBDev printing & parameter
*/
void DoMerge(FILE *fp, int i, DATALIST *loop1, int j, DATALIST *loop2, 
             REAL rms, REAL CADev, REAL CBDev, int *NewNumbers, int NClus)
{
   int OldClusNum,
       NewClusNum,
       k;

   fprintf(fp,"MERGED cluster %d (%s) with %d (%s), rmsd = %f, \
max CA deviation = %f, max CB deviation = %f\n",
           i+1,
           loop1->loopid,
           j+1,
           loop2->loopid,
           rms, CADev, CBDev);
   
   
   /* Merge these two clusters                                          */
   OldClusNum = MAX(NewNumbers[i], NewNumbers[j]);
   NewClusNum = MIN(NewNumbers[i], NewNumbers[j]);

   for(k=0; k<NClus; k++)
   {
      if(NewNumbers[k] == OldClusNum)
         NewNumbers[k] = NewClusNum;
   }


#ifdef DEBUG
   for(k=0; k<NClus; k++)
      fprintf(stderr,"%2d ",NewNumbers[k]);
   fprintf(stderr,"\n");
#endif
}


/************************************************************************/
/*>DATALIST *FindLoop(int *clusters, int NVec, int ClusNum, int loopnum)
   ---------------------------------------------------------------------
   Input:   int   *clusters     Clustring vector
            int   NVec          Number of vectors
            int   ClusNum       Cluster of interest
            int   loopnum       Example number within this cluster
   Returns: DATALIST *          Pointer to loop example in this cluster.
                                NULL if not found

   Finds the loopnum'th example within this cluster (counts from 0)

   12.09.95 Original    By: ACRM
*/
DATALIST *FindLoop(int *clusters, int NVec, int ClusNum, int loopnum)
{
   DATALIST *p;
   int      i, j,
            example = 0;

   for(i=0; i<NVec; i++)
   {
      if(clusters[i]==ClusNum)
      {
         if(example == loopnum)
         {
            /* Run through the global data linked list to find the
               example's loop pointer
            */
            for(p=gDataList, j=0; p!=NULL && j<i; NEXT(p), j++) ;
            return(p);
         }
         example++;
      }
   }
   
   return(NULL);
}


/************************************************************************/
/*>BOOL DefineCriticalResidues(FILE *fp, int *clusters, REAL **data, 
                               int NVec, int VecDim, REAL *crit, 
                               int NClus)
   -----------------------------------------------------------------
   Input:   FILE   *fp         Output file pointer
            int    *clusters   Clustering table
            REAL   **data      Data for clustering
            int    NVec        Number of vectors
            int    VecDim      Dimesion of vector
            REAL   *crit       Critical values
            int    NClus       Number of clusters

   Does all set up from definition of clusters to call routines which
   analyse critical residues.

   08.08.95 Original    By: ACRM
   15.08.95 Modified to use cluster vector rather than matrix
   04.10.95 Modified additionally to show details of residues conserved
            in at least one cluster. This required changing the old
            code to handle cinfo as an array.
   05.10.95 Removed First checking for MergeAllProperties() as this
            is now done on a per-residue basis within that routine.
   06.11.95 Added check on exclude list before processing
*/
BOOL DefineCriticalResidues(FILE *fp, int *clusters, REAL **data, 
                            int NVec, int VecDim, REAL *crit, int NClus)
{
   int         clusnum,
               i,  j,
               *NMembers,
               NCons,
               InfoPos   = 0,
               InfoStart = 0;
   LOOPINFO    *loopinfo;
   DATALIST    *p;
   CLUSTERINFO *cinfo;
   PDB         *pdb,
               *pdb_start,
               *pdb_end;
   RESSPEC     *ConsList = NULL;

   /* Allocate memory for maximum possible amount of loop data          */
   if((loopinfo=(LOOPINFO *)malloc(NVec * sizeof(LOOPINFO)))==NULL)
      return(FALSE);

   /* Blank the loopinfo structures                                     */
   for(i=0; i<NVec; i++)
      BlankLoopInfo(&(loopinfo[i]));

   /* Allocate memory for the cluster info structures                   */
   if((cinfo=(CLUSTERINFO *)malloc(NClus * sizeof(CLUSTERINFO)))==NULL)
   {
      free(loopinfo);
      return(FALSE);
   }

   /* Allocate array to store number of members of each cluster         */
   if((NMembers=(int *)malloc((NClus+1) * sizeof(int)))==NULL)
   {
      free(loopinfo);
      free(cinfo);
      return(FALSE);
   }

   for(i=0; i<=NClus; i++)
      NMembers[i] = 0;
   
   /* Initialise the properties lookup tables                           */
   InitProperties();
   
   /* Blank the cluster info structures                                 */
   for(i=0; i<NClus; i++)
      BlankClusterInfo(&(cinfo[i]));
   
   /* Print a header                                                    */
   fprintf(fp, "\nBEGIN CRITICALRESIDUES %d\n", NClus);

   /* For each cluster in turn                                          */
   for(clusnum=1; clusnum<=NClus; clusnum++)
   {
      /* Currently no members in this cluster                           */
      NMembers[clusnum] = 0;
      
      /* Run through the cluster table to find members of this cluster  */
      for(i=0; i<NVec; i++)
      {
         /* If this is the current cluster                              */
         if(clusters[i] == clusnum)
         {
            /* Find the PDB linked list for this example                */
            for(j=0, p=gDataList; j<i && p!=NULL; j++, NEXT(p));

            /* Find PDB pointers for the structure and the start and end 
               of the loop itself.
            */
            pdb       = p->allatompdb;
            pdb_start = FindResidueSpec(pdb, p->start);
            pdb_end   = FindResidueSpec(pdb, p->end);
            if(pdb_end != NULL)
               pdb_end   = FindNextResidue(pdb_end);
            
            if(pdb_start != NULL)
            {
               /* If this loop is not in the list of loops to be ignored
                  in sequence template analysis
               */
               if(!InStringList(gStringList, p->loopid))
               {
                  /* Store the loop properties in the array             */
                  if(!FindNeighbourProps(pdb, pdb_start, pdb_end, clusnum,
                                         &(loopinfo[InfoPos])))
                  {
                     free(loopinfo);
                     return(FALSE);
                  }
                  
                  InfoPos++;
                  (NMembers[clusnum])++;
               }
            }
            else
            {
               fprintf(stderr,"Unable to find start residue (%s) in PDB \
file (%s)\n",p->start,p->loopid);
               return(FALSE);
            }
         }  /* Correct cluster                                          */
      }  /* End of for(each vector)                                     */

      InfoStart += NMembers[clusnum-1];

      /* If there were some members in this cluster                     */
      if(NMembers[clusnum])
      {
         /* We've accumulated loop info for each member of this cluster, 
            so merge the properties, filling in a cluster info structure.
            Note that we use the number of members of the previous
            clusters as an offset for the loopinfo array. clusnum
            starts at 1 and NMembers[0] is set to 0
         */
         if(!MergeProperties(NMembers[clusnum], 
                             loopinfo+InfoStart, clusnum, 
                             &(cinfo[clusnum-1])))
         {
            fprintf(fp, "END CRITICALRESIDUES (failed!)\n");
            fprintf(stderr,"MergeProperties() failed\n");
            return(FALSE);
         }
         
         /* Print merged properties for this cluster                    */
         PrintMergedProperties(fp, clusnum, cinfo[clusnum-1], 
                               NMembers[clusnum]);
         fprintf(fp,"\n");
      }
   }  /* End of for(each cluster)                                       */
   fprintf(fp, "END CRITICALRESIDUES\n");

   /* Build a list of all the residues which are conserved in any one
      cluster
   */
   InfoStart = 0;
   if((ConsList = BuildConservedList(cinfo, NClus, &NCons))!=NULL)
   {
      /* Print a header                                                 */
      fprintf(fp, "\nBEGIN ALLCRITICALRESIDUES %d\n", NClus);

      /* Clean out our current cluster information                      */
      for(i=0; i<NClus; i++)
         CleanClusInfo(&(cinfo[i]));

      /* For each cluster in turn                                       */
      for(clusnum=1; clusnum<=NClus; clusnum++)
      {
         InfoStart += NMembers[clusnum-1];

         /* If there were some members in this cluster                  */
         if(NMembers[clusnum])
         {
            for(i=0; i<NVec; i++)
            {
               /* If this is our cluster                                */
               if(clusters[i] == clusnum)
               {
                  /* Find the PDB linked list for this example          */
                  for(j=0, p=gDataList; j<i && p!=NULL; j++, NEXT(p));

                  if(!InStringList(gStringList, p->loopid))
                  {
                     /* Find PDB pointer for the structure.             */
                     pdb = p->allatompdb;
                     
                     if(!MergeAllProperties(pdb, ConsList, NCons,
                                            &(cinfo[clusnum-1])))
                     {
                        fprintf(fp,"END ALLCRITICALRESIDUES (failed!)\n");
                        fprintf(stderr,"MergeAllProperties() failed\n");
                        return(FALSE);
                     }
                  }
                  
               }  /* Correct cluster                                    */
            }  /* End of for(each vector)                               */
         }  /* There were members of this cluster                       */
         
         /* Print merged properties for this cluster                    */
         PrintMergedProperties(fp, clusnum, cinfo[clusnum-1], 
                               NMembers[clusnum]);
         PrintDeletedResidues(fp, cinfo[clusnum-1], ConsList, NCons);
         fprintf(fp,"\n");
      }
      fprintf(fp, "END ALLCRITICALRESIDUES\n");
   }
   
   

   /* Clean up allocated memory in the loopinfo and clusinfo structures */
   CleanLoopInfo(loopinfo, InfoPos);  /* 07.11.95 Changed to InfoPos
                                         rather than NVec
                                      */
   for(i=0; i<NClus; i++)
      CleanClusInfo(&(cinfo[i]));
      
   free(loopinfo);
   free(cinfo);
   free(NMembers);

   return(TRUE);
}


