/* Use absolute residue conservation                                    */
#define USE_ABSCONS

/* Use Loop--f/w HBonds                                                 */
#define USE_HBONDS

/* Use buried HPhobs in loop                                            */
#define USE_HPHOB

/* Use f/w h/phob partners                                              */
#define USE_HPHOB_PARTNERS 

/* Use conserved Gly/Pro                                                */
#define USE_GLYPRO             

/* Use cis Pro even if there's only one                                 */
#define USE_CISPRO

/* Use Loop--loop s/c-m/c HBonds                                        */
#define USE_LOOP_SM_HBONDS     

/* When unifying SDR lists, choose positions from all clusters of the
   same loop length
*/
#define UNIFY_ON_LENGTH

/* When unifying SDR lists, choose positions from all large clusters    */
#define UNIFY_ON_LARGE_CLUSTER

/* When unifying SDR lists on length, exclude added residues if they do
   not give added descriminatory power
*/
#define EXCLUDE_NONINFORM                                            




/* Report reasons for residues                                          */
/* #define REPORT_REASONS */

/* Various debugging                                                    */
/* #define DEBUG */

/*************************************************************************

   Program:    FindSDRs
   File:       FindSDRs.c
   
   Version:    V1.0
   Date:       02.02.96
   Function:   Find SDRs in a set of loops
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1996
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
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
   Reads the output of CLAN and attempts to define SDRs using the
   PDB files and sequence templates for each cluster.

   The algorithm is as follows:

   For each cluster:

   1. If a residue is absolutely conserved and the cluster has at least
      MINABSCONS (5) members it is defined as key
   2. If a Gly/Pro is absolutely conserved and the cluster has at least
      MINGLYPRO (2) members it is defined as key
   3. Any residues which make sidechain HBonds between loop and framework
      in every member of the cluster are defined as key
   4. Any residues which make sidechain/backbone HBonds within the loop
      in every member of the cluster are defined as key
   5. Any residues in the loop which are buried (mean SA < SACUT (=3.0))
      hydrophobics in every member of the cluster are defined as key
   6. Framework hydrophobic residues which make sidechain interactions
      (atom distance < sqrt(HPHOBCONTDISTSQ) (=5.0)) with loop key
      hydrophobics in every member of the cluster are defined as key


   To report unified SDRs:

   7. A list of key positions defined above in any cluster (of any loop 
      length) with at least MINCLUSSIZE (5) members is assembled.
   9. For each cluster, the key residues defined in step 7 are appended
      to the list generated in steps 1--6.
   8. For each cluster, key positions from small clusters (< MINCLUSSIZE)
      are appended to the list if the loop length matches
      [OPTIONALLY: There must also be some ``added value'' (i.e. the amino 
      acid at this position descriminates between the conformations)].


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1  02.02.96 Original development version
   V1.0  29.04.96 Finds rogues when parent isn't the largest cluster
                  of a given length
   V1.0a 30.01.09 Compile cleanups

*************************************************************************/
/* Includes
*/
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "bioplib/general.h"
#include "bioplib/fsscanf.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/hbond.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"

#include "resprops.h"
#include "decr2.h"

/************************************************************************/
/* Defines and macros
*/

/* Buffer size                                                          */
#define MAXBUFF         160

/* Max length of a word to pull out of the buffer                       */
#define MAXWORD         32

/* malloc() step                                                        */
#define ALLOCQUANTUM    16  

/* Max amino acid types                                                 */
#define MAXRES          24  

/* Max SA for a buried res                                              */
#define SACUT           3.0 

/* The copy command                                                     */
#define CPCOMMAND       "cp"

/* Command to create an access file from PDB                            */
#define SOLVACC         "pdbsolv %s | pdbsumbval -a -q > %s" 

/* Dir for temp files; may need to be a blank string                    */
#define TEMPDIR         ""

/* Square distance considered to be a hphob contact                     */
#define HPHOBCONTDISTSQ ((REAL)25.0)          

/* Min number of members of a cluster when reporting unified SDR lists  */
#define MINCLUSSIZE     5

/* Min number of members if a cluster when considering absolute 
   conservation in defining SDRs
*/
#define MINABSCONS      5

/* Min number of members of a cluster when considering absolute 
   conservation of Gly/Pro in defining SDRs       
*/
#define MINGLYPRO       2

/* Position codes used by SDRLIST->position                             */
#define POS_NOCONTACT   0
#define POS_CONTACT     1
#define POS_LOOP        2

/* Onlength codes used by SDRLIST->onlength                             */
#define OL_FALSE        0
#define OL_ONLENGTH     1
#define OL_DELETABLE    2

/* Structure (linked list) used to store unified SDR list               */
typedef struct _sdrlist
{
   struct _sdrlist *next;
   int             resnum,
                   nobsres,
                   position;
   char            chain,
                   insert,
                   obsres[MAXRES+1];
   int             onlength;
}  SDRLIST;

/* Structure which defines each loop name and the cluster to which it
   belongs
*/
typedef struct
{
   int  cluster;
   char filename[MAXBUFF],
        firstres[16],
        lastres[16];
}  LOOPCLUS;

/* Structure which defines the characteristics of a cluster             */
typedef struct
{
   SDRLIST *sdrlist;
   int     *resnum;
   char    *chain,
           *insert;
   PROP_T  *props;
   int     NRes,                 /* Number of common residue ids        */
           length,               /* Length of loop itself               */
           NMembers,
           ArraySize,
           *count,
           *PartnerCount,
           rogue;
   BOOL    *absolute,
           *deleted,
           *key,
           *flagged;
   char    *ConsRes;
}  CLUSINFO;

/* Structure which defines mean and sd Ooi values and hydrophobicity flag
   for each residue name
*/
typedef struct
{
   char   resnam[8];
   REAL   mean, sd;
   BOOL   hphob;
}  OOIDATA;


/************************************************************************/
/* Globals
*/
LOOPCLUS *gLoopClus = NULL;  /* Store loop name and cluster number      */
CLUSINFO *gClusInfo = NULL;  /* Store sequence templates for clusters   */
OOIDATA  gOoiData[MAXRES];   /* Store mean/sd Ooi data for general loops*/
int      gMinLoopLength,
         gMaxLoopLength;


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
int ReadClanFile(FILE *in, int *NLoops);
BOOL ReadAssignments(FILE *fp);
void FreeGlobalStorage(int nclus);
void StorePDBNameCluster(char *inbuff, int LoopNum);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *KeepSA);
void Usage(void);
void BlankTemplates(int nclus);
BOOL ExpandTemplateArrays(CLUSINFO *ClusInfo);
BOOL ReadTemplates(FILE *in);
BOOL FindSDRs(int nclus, int nloops, BOOL KeepSA);
void ReportSDRs(FILE *out, int nclus);
void Report(CLUSINFO *ClusInfo, int residx, char *reason);
void FillOoiData(void);
BOOL IsInRange(char *resspec, char *firstres, char *lastres);
PDB *ReadPDBAsSA(char *filename, BOOL KeepSAFile);
void MarkPartners(CLUSINFO *ClusInfo, PDB *pdb, PDB *res, 
                  char *firstres, char *lastres);
BOOL MakeSCContact(PDB *res1, PDB *res2);
BOOL MarkHPhob(CLUSINFO *ClusInfo, int clusnum, int nloops, 
               BOOL KeepSA);
BOOL MarkHBonders(CLUSINFO *ClusInfo, int clusnum, int nloops);
SDRLIST *InSDRList(SDRLIST *sdrlist, char chain, int resnum, char insert);
BOOL FillSDRsForCluster(SDRLIST *sdrlist, int clusnum, int nloops);
BOOL ReportUnifiedSDRs(FILE *out, int nclus, int nloops);
void PrintSDRList(FILE *out, SDRLIST *sdrlist);
void indexint(int n, int *arrin, int *indx);
void FlagNonInformativeSDRs(int nclus);
BOOL ValueIsAdded(SDRLIST *s1, SDRLIST *s2);
void FlagRogueClusters(int nclus, int nloops);
BOOL IsRogue(int clus, int LargestClus);
BOOL IsCisProline(CLUSINFO *ClusInfo, int clusnum, int resoffset, 
                  int nloops);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for defining SDRs from output of CLAN

   02.02.96 Original   By: ACRM
   30.01.09 Initialize some variables
*/
int main(int argc, char **argv)
{
   char InFile[MAXBUFF],
        OutFile[MAXBUFF];
   FILE *in  = stdin,
        *out = stdout;
   int  nclus = 0, nloops;
   BOOL KeepSA;

   if(ParseCmdLine(argc, argv, InFile, OutFile, &KeepSA))
   {
      if(OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if((nclus=ReadClanFile(in, &nloops))!=0)
         {
            /* Initialise various data structures                       */
            FillOoiData();
            InitProperties();
            
            if(FindSDRs(nclus, nloops, KeepSA))
            {
               ReportSDRs(out, nclus);
               ReportUnifiedSDRs(out, nclus, nloops);
            }
            else
            {
               fprintf(stderr,"Unable to run FindSDR code\n");
            }
         }
         else
         {
            fprintf(stderr,"Unable to read data from CLAN file\n");
            return(0);
         }
      }
   }
   else
   {
      Usage();
   }

   FreeGlobalStorage(nclus);
   
   return(0);
}


/************************************************************************/
/*>int ReadClanFile(FILE *in, int *NLoops)
   ---------------------------------------
   Read in the CLAN output file
   Stores the cluster assignments (i.e. each filename and cluster number)
   in the global gLoopClus[] array and stores the sequence templates
   for each cluster in the gClusInfo[] array.

   02.02.96 Original   By: ACRM
   06.02.96 Outputs NLoops
*/
int ReadClanFile(FILE *in, int *NLoops)
{
   int  nloops   = 0,
        nclus    = 0;
   char word[MAXWORD],
        *buffer,
        *buffp;
   BOOL GotAssignments = FALSE,
        GotTemplates   = FALSE;

   *NLoops = 0;
   
   while((buffer=fgetsany(in))!=NULL)
   {
      if(nloops == 0)
      {
         if(strstr(buffer,"NLOOPS"))
         {
            buffp = buffer;
            
            /* Pull out NLOOPS                                          */
            buffp = GetWord(buffp,word,MAXWORD);
            /* Pull out the actual number                               */
            buffp = GetWord(buffp,word,MAXWORD);
            if(!sscanf(word,"%d",&nloops) || nloops==0)
            {
               fprintf(stderr,"Unable to read NLOOPS from clan file\n");
               return(0);
            }
            /* Allocate loop storage                                    */
            if((gLoopClus = (LOOPCLUS *)
                malloc(nloops * sizeof(LOOPCLUS)))==NULL)
            {
               fprintf(stderr,"No memory to store loop/cluster list\n");
               return(0);
            }
         }
      }
      else /* We've now got the number of loops, look for data          */
      {
         if(strstr(buffer,"BEGIN ASSIGNMENTS"))
         {
            if(!ReadAssignments(in))
            {
               fprintf(stderr,"Failed to read ASSIGNMENTS\n");
               return(0);
            }
            GotAssignments = TRUE;
         }

         if(strstr(buffer,"BEGIN CRITICALRESIDUES"))
         {
            buffp = buffer;
            
            /* Pull out BEGIN                                           */
            buffp = GetWord(buffp,word,MAXWORD);
            /* Pull out CRITICALRESIDUES                                */
            buffp = GetWord(buffp,word,MAXWORD);
            /* Pull out the actual number                               */
            buffp = GetWord(buffp,word,MAXWORD);
            if(!sscanf(word,"%d",&nclus) || nclus==0)
            {
               fprintf(stderr,"Unable to read number of clusters from \
clan file\n");
               return(0);
            }

            /* Allocate storage for sequence templates for each cluster */
            if((gClusInfo = (CLUSINFO *)
                malloc(nclus * sizeof(CLUSINFO)))==NULL)
            {
               fprintf(stderr,"No memory to store cluster sequence \
templates\n");
               return(0);
            }

            BlankTemplates(nclus);

            if(!ReadTemplates(in))
            {
               fprintf(stderr,"Failed to read CRITICALRESIDUES\n");
               return(0);
            }
            GotTemplates = TRUE;
         }
      }
      
      free(buffer);
   }
   if(errno==ENOMEM)
   {
      fprintf(stderr,"No memory for file buffer\n");
      return(0);
   }

   if(!GotAssignments)
   {
      fprintf(stderr,"Failed to find BEGIN ASSIGNMENTS record\n");
      return(0);
   }
   if(!GotTemplates)
   {
      fprintf(stderr,"Failed to find BEGIN CRITICALRESIDUES record\n");
      return(0);
   }
   
   *NLoops = nloops;
   return(nclus);
}


/************************************************************************/
/*>BOOL ReadAssignments(FILE *fp)
   ------------------------------
   Read the ASSIGNMENTS section of the CLAN file

   02.02.96 Original   By: ACRM
   07.02.96 Added parentheses around while()
*/
BOOL ReadAssignments(FILE *fp)
{
   int  i = 0;
   char *buffer;
   
   while((buffer=fgetsany(fp))!=NULL)
   {
      if(strstr(buffer,"END ASSIGNMENTS"))
         break;

      StorePDBNameCluster(buffer,i++);
      free(buffer);
   }
   if(buffer==NULL)
      return(FALSE);
   
   free(buffer);
   
   return(TRUE);
}

         
/************************************************************************/
/*>void FreeGlobalStorage(int nclus)
   ---------------------------------
   Free any globally allocated storage

   02.02.96 Original   By: ACRM
   08.02.96 Added count and PartnerCount
   09.02.96 Added flagged
*/
void FreeGlobalStorage(int nclus)
{
   int i;
   
   if(gLoopClus != NULL)
      free(gLoopClus);

   if(gClusInfo != NULL)
   {
      for(i=0; i<nclus; i++)
      {
         free(gClusInfo[i].resnum);
         free(gClusInfo[i].chain);
         free(gClusInfo[i].insert);
         free(gClusInfo[i].props);
         free(gClusInfo[i].absolute);
         free(gClusInfo[i].count);
         free(gClusInfo[i].PartnerCount);
         free(gClusInfo[i].deleted);
         free(gClusInfo[i].flagged);
         free(gClusInfo[i].ConsRes);
      }
      free(gClusInfo);
   }
}


/************************************************************************/
/*>void StorePDBNameCluster(char *inbuff, int LoopNum)
   ---------------------------------------------------
   Store the PDB name, first and last residue and cluster number from
   withing the ASSIGNMENTS section of the file

   02.02.96 Original   By: ACRM
*/
void StorePDBNameCluster(char *inbuff, int LoopNum)
{
   char buff[MAXBUFF],
        *buffp,
        *chp;

   gLoopClus[LoopNum].firstres[0] = '\0';
   gLoopClus[LoopNum].lastres[0]  = '\0';
   
   sscanf(inbuff,"%d %s", &(gLoopClus[LoopNum].cluster),
                          buff);
   buffp = buff;
   
   if((chp = strchr(buffp,'-'))!=NULL)
   {
      *chp = '\0';
   }
   strcpy(gLoopClus[LoopNum].filename,buffp);

   if(chp!=NULL)
   {
      buffp=chp+1;
      if((chp = strchr(buffp,'-'))!=NULL)
      {
         *chp = '\0';
      }
      strcpy(gLoopClus[LoopNum].firstres,buffp);

      if(chp!=NULL)
      {
         buffp=chp+1;
         strcpy(gLoopClus[LoopNum].lastres,buffp);
      }
   }
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *KeepSA)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            BOOL   *KeepSA      Should generated SA files be kept?
   Returns: BOOL                Success?

   Parse the command line
   
   02.02.96 Original    By: ACRM
   08.02.96 Added KeepSA flag
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *KeepSA)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *KeepSA = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'k':
            *KeepSA = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   02.02.96 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nFindSDRs V1.0 (c) 1996, Dr. Andrew C.R. Martin, \
UCL\n");

   fprintf(stderr,"\nUsage: findsdrs [-k] [clanfile [outfile]]\n");
   fprintf(stderr,"       -k  Keep any generated SA files\n");

   fprintf(stderr,"\nTakes the output from the Clan loop clustering \
program and reads the\n");
   fprintf(stderr,"PDB files specified within this file to attempt to \
define key residues.\n\n");
}

/************************************************************************/
/*>void BlankTemplates(int nclus)
   ------------------------------
   Sets all pointers in the template structures to NULL and zeros the
   counts.

   02.02.96 Original   By: ACRM
   06.02.96 Added key and count
   08.02.96 Added PartnerCount
   09.02.96 Added flagged
   15.03.96 Added sdrlist
*/
void BlankTemplates(int nclus)
{
   int i;

   for(i=0; i<nclus; i++)
   {
      gClusInfo[i].sdrlist        = NULL;
      gClusInfo[i].resnum         = NULL;
      gClusInfo[i].chain          = NULL,
      gClusInfo[i].insert         = NULL;
      gClusInfo[i].props          = NULL;
      gClusInfo[i].absolute       = NULL,
      gClusInfo[i].deleted        = NULL;
      gClusInfo[i].flagged        = NULL;
      gClusInfo[i].ConsRes        = NULL;
      gClusInfo[i].key            = NULL;
      gClusInfo[i].count          = NULL;
      gClusInfo[i].PartnerCount   = NULL;
      gClusInfo[i].ArraySize      = 0;
      gClusInfo[i].NRes           = 0;
   }
}


/************************************************************************/
/*>BOOL ExpandTemplateArrays(CLUSINFO *ClusInfo)
   ---------------------------------------------
   Allocate more space in the template arrays for a given cluster

   02.02.96 Original   By: ACRM
   06.02.96 Added key and count
   08.02.96 Added PartnerCount
   09.02.96 Added flagged
*/
BOOL ExpandTemplateArrays(CLUSINFO *ClusInfo)
{
   int size;
   
   size = ClusInfo->ArraySize + ALLOCQUANTUM;

   if((ClusInfo->resnum = 
       (int *)realloc(ClusInfo->resnum, size * sizeof(int)))==NULL)
      return(FALSE);
   if((ClusInfo->chain  = 
       (char *)realloc(ClusInfo->chain,  size * sizeof(char)))==NULL)
      return(FALSE);
   if((ClusInfo->insert = 
       (char *)realloc(ClusInfo->insert, size * sizeof(char)))==NULL)
      return(FALSE);
   if((ClusInfo->props = 
       (PROP_T *)realloc(ClusInfo->props, size * sizeof(PROP_T)))==NULL)
      return(FALSE);
   if((ClusInfo->absolute = 
       (BOOL *)realloc(ClusInfo->absolute, size * sizeof(BOOL)))==NULL)
      return(FALSE);
   if((ClusInfo->deleted = 
       (BOOL *)realloc(ClusInfo->deleted, size * sizeof(BOOL)))==NULL)
      return(FALSE);
   if((ClusInfo->flagged = 
       (BOOL *)realloc(ClusInfo->flagged, size * sizeof(BOOL)))==NULL)
      return(FALSE);
   if((ClusInfo->ConsRes = 
       (char *)realloc(ClusInfo->ConsRes, size * sizeof(char)))==NULL)
      return(FALSE);
   if((ClusInfo->key = 
       (BOOL *)realloc(ClusInfo->key, size * sizeof(BOOL)))==NULL)
      return(FALSE);
   if((ClusInfo->count = 
       (int *)realloc(ClusInfo->count, size * sizeof(int)))==NULL)
      return(FALSE);
   if((ClusInfo->PartnerCount = 
       (int *)realloc(ClusInfo->PartnerCount, size * sizeof(int)))==NULL)
      return(FALSE);

   ClusInfo->ArraySize = size;
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL ReadTemplates(FILE *in)
   ----------------------------
   Read the CRITICALRESIDUES (template) section of the CLAN file

   02.02.96 Original   By: ACRM
*/
BOOL ReadTemplates(FILE *in)
{
   char *buffer,
        *buffp,
        *chp,
        word[32];
   int  clusnum,
        length,
        offset = (-1),
        count,
        nmembers;
   
   while((buffer=fgetsany(in))!=NULL)
   {
      /* If we hit a CLUSTER record, extract clusnum, length, members   */
      if(strstr(buffer, "CLUSTER"))
      {
         /* Remove punctuation from the record                          */
         if((buffp=strchr(buffer,','))!=NULL)
            *buffp = ' ';
         if((buffp=strchr(buffer,'('))!=NULL)
            *buffp = ' ';
         if((buffp=strchr(buffer,')'))!=NULL)
            *buffp = ' ';
         
         if(sscanf(buffer,"%s %d %s %s %d %s %s %d",
                   word,
                   &clusnum,
                   word, word,
                   &length,
                   word, word,
                   &nmembers)!=8)
         {
            fprintf(stderr,"Unable to read CLUSTER record:\n%s\n",
                    buffer);
            return(FALSE);
         }

         if(clusnum > 0)
         {
            offset = clusnum-1;
            
            /* Store these in the CLUSINFO array and allocate storage 
               space
            */
            gClusInfo[offset].length   = length;
            gClusInfo[offset].NMembers = nmembers;
            gClusInfo[offset].NRes     = 0;
            if(!ExpandTemplateArrays(&(gClusInfo[offset])))
            {
               fprintf(stderr,"No memory to store cluster sequence \
templates\n");
               return(FALSE);
            }
         }
      }
      else if(strstr(buffer,"WARNING") || (strlen(buffer) < 13))
      {
         /* Do nothing                                                  */
         ;
      }
      else if(strstr(buffer,"END CRITICALRESIDUES"))
      {
         free(buffer);
         break;
      }
      else /* It's an actual residue specification record               */
      {
         count = gClusInfo[offset].NRes;

         /* Make more storage space if needed                           */
         if(count >= gClusInfo[offset].ArraySize)
         {
            if(!ExpandTemplateArrays(&(gClusInfo[offset])))
            {
               fprintf(stderr,"No memory to store cluster sequence \
templates\n");
               return(FALSE);
            }
         }
            
         /* Read the required values out of the buffer                  */
         fsscanf(buffer,"%c%3d%c%1x%6s",
                 &(gClusInfo[offset].chain[count]),
                 &(gClusInfo[offset].resnum[count]),
                 &(gClusInfo[offset].insert[count]),
                 word);
         sscanf(word,"%hx",&(gClusInfo[offset].props[count]));

         /* Other flags                                                 */
         gClusInfo[offset].absolute[count] = FALSE;
         gClusInfo[offset].deleted[count]  = FALSE;
         gClusInfo[offset].ConsRes[count]  = ' ';
         
         if((buffp=strstr(buffer,"deleted"))!=NULL)
         {
            gClusInfo[offset].deleted[count] = TRUE;
            gClusInfo[offset].ConsRes[count] = '-';
         }
         else if((buffp=strstr(buffer,"CONSERVED"))!=NULL)
         {
            gClusInfo[offset].absolute[count] = TRUE;
            if((chp=strchr(buffp,')'))!=NULL)
            {
               *chp = '\0';
               if((chp=strchr(buffp,'('))!=NULL)
               {
                  buffp = chp+1;
                  gClusInfo[offset].ConsRes[count] = *buffp;
               }
            }
         }
         (gClusInfo[offset].NRes)++;
      }
         
      free(buffer);
   }

   if(buffer==NULL && errno==ENOMEM)
      return(FALSE);

   return(TRUE);
}


/************************************************************************/
/*>BOOL FindSDRs(int nclus, int nloops, BOOL KeepSA)
   -------------------------------------------------
   Main routine for doing the work of finding the SDRs

   06.02.96 Original   By: ACRM
   08.02.96 Added KeepSA flag
*/
BOOL FindSDRs(int nclus, int nloops, BOOL KeepSA)
{
   int clus,
       i;
   
   /* For each cluster                                                  */
   for(clus=0; clus<nclus; clus++)
   {
      /* Start off assuming nothing is a key residue                    */
      for(i=0; i<gClusInfo[clus].NRes; i++)
         gClusInfo[clus].key[i] = FALSE;


#ifdef USE_ABSCONS
      /* If there are more than MINABSCONS members in the cluster, residues
         which are absolutely conserved are marked
      */
      if(gClusInfo[clus].NMembers >= MINABSCONS)
      {
         for(i=0; i<gClusInfo[clus].NRes; i++)
         {
            if(gClusInfo[clus].absolute[i])
            {
               gClusInfo[clus].key[i] = TRUE;
#ifdef REPORT_REASONS
               Report(&(gClusInfo[clus]), i, "Absolute Conservation");
#endif
            }
         }
      }
#endif
      
#ifdef USE_GLYPRO
      /* If there are more than MINGLYPRO members in the cluster, residues
         which are absolutely conserved are marked
      */
      if(gClusInfo[clus].NMembers >= MINGLYPRO)
      {
         for(i=0; i<gClusInfo[clus].NRes; i++)
         {
            if(gClusInfo[clus].absolute[i] &&
               ((gClusInfo[clus].ConsRes[i] == 'G') ||
                (gClusInfo[clus].ConsRes[i] == 'P')))
            {
               gClusInfo[clus].key[i] = TRUE;
#ifdef REPORT_REASONS
               Report(&(gClusInfo[clus]), i, "Conserved G/P");
#endif
            }
         }
      }
#ifdef USE_CISPRO
      else
      {
         for(i=0; i<gClusInfo[clus].NRes; i++)
         {
            if(gClusInfo[clus].absolute[i] &&
               gClusInfo[clus].ConsRes[i]  == 'P')
            {
               if(IsCisProline(&(gClusInfo[clus]), clus+1, i, nloops))
               {
                  gClusInfo[clus].key[i] = TRUE;
#ifdef REPORT_REASONS
                  Report(&(gClusInfo[clus]), i, "Cis-Pro");
#endif
               }
            }
         }
      }
#endif
#endif
      
#ifdef USE_HBONDS
      /* Any residues which make s/c HBonds are marked as key        */
      if(!MarkHBonders(&(gClusInfo[clus]), clus+1, nloops))
         return(FALSE);
#endif
      
#ifdef USE_HPHOB
      /* Any buried hydrophobics are marked as key                   */
      if(!MarkHPhob(&(gClusInfo[clus]), clus+1, nloops, KeepSA))
         return(FALSE);
#endif
   }

   return(TRUE);
}



/************************************************************************/
/*>void ReportSDRs(FILE *out, int nclus)
   -------------------------------------
   Reports the SDRs. No unification of which residues are key between
   the clusters is performed

   06.02.96 Original   By: ACRM
*/
void ReportSDRs(FILE *out, int nclus)
{
   int i, j;
   
   for(i=0; i<nclus; i++)
   {
      fprintf(out,"\nCLUSTER %d (Length = %d, Members = %d)\n",
              i+1,gClusInfo[i].length,gClusInfo[i].NMembers);
      
      for(j=0; j<gClusInfo[i].NRes; j++)
      {
         if(gClusInfo[i].key[j])
         {
            fprintf(out,"%c %4d %c  0x%04x ",
                    gClusInfo[i].chain[j],
                    gClusInfo[i].resnum[j],
                    gClusInfo[i].insert[j],
                    gClusInfo[i].props[j]);

            PrintProps(out,gClusInfo[i].props[j], FALSE);
            
            if(gClusInfo[i].absolute[j])
            {
               fprintf(out," [CONSERVED] (%c)",gClusInfo[i].ConsRes[j]);
            }
            else
            {
               PrintSampleResidues(out,gClusInfo[i].props[j], FALSE);
            }
            
            fprintf(out,"\n");
         }
      }
   }
}


/************************************************************************/
/*>void Report(CLUSINFO *ClusInfo, int residx, char *reason)
   ---------------------------------------------------------
   Reports the reson why a residue has been defined as key

   09.02.96 Original   By: ACRM
*/
void Report(CLUSINFO *ClusInfo, int residx, char *reason)
{
   fprintf(stderr,"Residue %c%d%c %s\n",
           ClusInfo->chain[residx],
           ClusInfo->resnum[residx],
           ClusInfo->insert[residx],
           reason);
}


/************************************************************************/
/*>void FillOoiData(void)
   ----------------------
   Sets up Ooi(6.5,resmean) data from analysis of all protein loops

   06.02.96 Original   By: ACRM
*/
void FillOoiData(void)
{
   strcpy(gOoiData[0].resnam, "ALA ");
   gOoiData[0].mean   = 47.769764;
   gOoiData[0].sd     = 12.190481;
   gOoiData[0].hphob  = TRUE;
   
   strcpy(gOoiData[1].resnam, "CYS ");
   gOoiData[1].mean   = 54.118146;
   gOoiData[1].sd     = 9.775772;
   gOoiData[1].hphob  = TRUE;
   
   strcpy(gOoiData[2].resnam, "ASP ");
   gOoiData[2].mean   = 44.637253;
   gOoiData[2].sd     = 11.047142;
   gOoiData[2].hphob  = FALSE;
   
   strcpy(gOoiData[3].resnam, "GLU ");
   gOoiData[3].mean   = 41.081560;
   gOoiData[3].sd     = 10.852106;
   gOoiData[3].hphob  = FALSE;
   
   strcpy(gOoiData[4].resnam, "PHE ");
   gOoiData[4].mean   = 52.152497;
   gOoiData[4].sd     = 9.640339;
   gOoiData[4].hphob  = TRUE;
   
   strcpy(gOoiData[5].resnam, "GLY ");
   gOoiData[5].mean   = 47.194653;
   gOoiData[5].sd     = 12.340285;
   gOoiData[5].hphob  = FALSE;

   strcpy(gOoiData[6].resnam, "HIS ");
   gOoiData[6].mean   = 49.390870;
   gOoiData[6].sd     = 12.144934;
   gOoiData[6].hphob  = FALSE;
   
   strcpy(gOoiData[7].resnam, "ILE ");
   gOoiData[7].mean   = 50.533408;
   gOoiData[7].sd     = 9.875437;
   gOoiData[7].hphob  = TRUE;
   
   strcpy(gOoiData[8].resnam, "LYS ");
   gOoiData[8].mean   = 39.634909;
   gOoiData[8].sd     = 9.690623;
   gOoiData[8].hphob  = FALSE;
   
   strcpy(gOoiData[9].resnam, "LEU ");
   gOoiData[9].mean   = 50.446448;
   gOoiData[9].sd     = 9.488016;
   gOoiData[9].hphob  = TRUE;
   
   strcpy(gOoiData[10].resnam, "MET ");
   gOoiData[10].mean  = 48.193995;
   gOoiData[10].sd    = 12.442053;
   gOoiData[10].hphob = TRUE;
   
   strcpy(gOoiData[11].resnam, "ASN ");
   gOoiData[11].mean  = 45.770615;
   gOoiData[11].sd    = 11.542167;
   gOoiData[11].hphob = FALSE;
   
   strcpy(gOoiData[12].resnam, "PRO ");
   gOoiData[12].mean  = 46.169103;
   gOoiData[12].sd    = 10.762206;
   gOoiData[12].hphob = FALSE;
   
   strcpy(gOoiData[13].resnam, "GLN ");
   gOoiData[13].mean  = 43.548330;
   gOoiData[13].sd    = 11.127800;
   gOoiData[13].hphob = FALSE;
   
   strcpy(gOoiData[14].resnam, "ARG ");
   gOoiData[14].mean  = 44.059901;
   gOoiData[14].sd    = 11.499307;
   gOoiData[14].hphob = FALSE;
   
   strcpy(gOoiData[15].resnam, "SER ");
   gOoiData[15].mean  = 46.334175;
   gOoiData[15].sd    = 12.230367;
   gOoiData[15].hphob = FALSE;
   
   strcpy(gOoiData[16].resnam, "THR ");
   gOoiData[16].mean  = 47.319994;
   gOoiData[16].sd    = 11.508988;
   gOoiData[16].hphob = FALSE;
   
   strcpy(gOoiData[17].resnam, "VAL ");
   gOoiData[17].mean  = 50.048591;
   gOoiData[17].sd    = 10.195003;
   gOoiData[17].hphob = TRUE;
   
   strcpy(gOoiData[18].resnam, "TRP ");
   gOoiData[18].mean  = 54.772689;
   gOoiData[18].sd    = 9.012706;
   gOoiData[18].hphob = TRUE;
   
   strcpy(gOoiData[19].resnam, "TYR ");
   gOoiData[19].mean  = 51.645258;
   gOoiData[19].sd    = 9.760537;
   gOoiData[19].hphob = TRUE;
   
   strcpy(gOoiData[20].resnam, "UNK ");
   gOoiData[20].mean  = 42.753234;
   gOoiData[20].sd    = 12.584526;
   gOoiData[20].hphob = FALSE;
   
   strcpy(gOoiData[21].resnam, "GLX ");
   gOoiData[21].mean  = 43.943333;
   gOoiData[21].sd    = 7.774092;
   gOoiData[21].hphob = FALSE;
   
   strcpy(gOoiData[22].resnam, "ASX ");
   gOoiData[22].mean  = 41.491000;
   gOoiData[22].sd    = 7.759922;
   gOoiData[22].hphob = FALSE;
   
   strcpy(gOoiData[23].resnam, "PCA ");
   gOoiData[23].mean  = 27.750000;
   gOoiData[23].sd    = 11.199888;
   gOoiData[23].hphob = FALSE;
}


/************************************************************************/
/*>BOOL IsInRange(char *resspec, char *firstres, char *lastres)
   ------------------------------------------------------------
   Takes 3 residue specs as <chain><resnum><insert> and returns TRUE
   if the first is within the range defined by the other two.

   07.02.96 Original   By: ACRM
*/
BOOL IsInRange(char *resspec, char *firstres, char *lastres)
{
   char chain,  firstchain,  lastchain,
        insert, firstinsert, lastinsert;
   int  resnum, firstresnum, lastresnum;
   
   if(ParseResSpec(resspec, &chain, &resnum, &insert))
   {
      if(ParseResSpec(firstres, &firstchain, &firstresnum, &firstinsert))
      {
         if(ParseResSpec(lastres, &lastchain, &lastresnum, &lastinsert))
         {
            /* If chains match                                          */
            if((chain     == firstchain) &&
               (lastchain == firstchain))
            {
               /* If residue number is *within* the range, return TRUE  */
               if((resnum > firstresnum) && (resnum < lastresnum))
                  return(TRUE);
               
               /* If the range has a single residue number, check both
                  inserts
               */
               if((resnum == firstresnum) && (resnum == lastresnum))
               {
                  if(((int)insert >= (int)firstinsert) &&
                     ((int)insert <= (int)lastinsert))
                     return(TRUE);
               }

               /* If residue number matches ends of range check insert  */
               if(((resnum == firstresnum) && 
                   ((int)insert >= (int)firstinsert)) ||
                  ((resnum == lastresnum) && 
                   ((int)insert <= (int)lastinsert)))
                  return(TRUE);
            }
         }
      }
   }
   
   return(FALSE);
}


/************************************************************************/
/*>PDB *ReadPDBAsSA(char *filename, BOOL KeepSAFile)
   -------------------------------------------------
   Given a filename, looks to see if that file with a .sa extension
   exists in TEMPDIR. If so reads it as PDB. If not, calls appropriate
   programs to create it from the specified PDB file and deletes it
   again unless the KeepSAFile flag is set

   07.02.96 Original   By: ACRM
   08.02.96 Added KeepSAFile option
            Added code to check whether file already exists in TEMPDIR
*/
PDB *ReadPDBAsSA(char *filename, BOOL KeepSAFile)
{
   char filestem[MAXBUFF],
        tempfile[MAXBUFF],
        safile[MAXBUFF],
        buffer[MAXBUFF],
        *ext;
   int  natom;
   BOOL SAExists  = FALSE,
        PDBExists = FALSE;
   PDB  *pdb;
   FILE *fp;
   
   /* Get the filename extension                                        */
   ext = filename + strlen(filename);
   while((ext > filename) && (*(ext-1) != '.'))
      ext--;

   /* Get the filestem out of the filename                              */
   GetFilestem(filename, filestem);
   /* Create the name for the temp file                                 */
   sprintf(tempfile,"%s%s",TEMPDIR,filestem);
   strcpy(safile, tempfile);
   
   /* Append the .sa extension                                          */
   SetExtn(safile,"sa");
   /* Set the extension to that of the input file                       */
   SetExtn(tempfile,ext);

   /* See if the filename specified is the same as what it would be
      if it's in TEMPDIR
   */
   if(!strcmp(filename,tempfile))
      PDBExists = TRUE;

   /* See if the .sa file already exists                                */
   if(access(safile,F_OK)==0)
      SAExists = TRUE;
   
   /* If the SA file doesn't exist...                                   */
   if(!SAExists)
   {
      /* If the filename for the original PDB file and the temp file
         differ, copy the original to the temp
      */
      if(!PDBExists)
      {
         sprintf(buffer,"%s %s %s",CPCOMMAND,filename,tempfile);
         system(buffer);
      }
      
      /* Create the command for building the SA file and run it         */
      sprintf(buffer,SOLVACC,tempfile,safile);
      system(buffer);

      /* Remove our temp file if it didn't exists before we started     */
      if(!PDBExists)
         unlink(tempfile);
   }
   
   /* Open the SA file                                                  */
   if((fp=fopen(safile,"r"))==NULL)
   {
      fprintf(stderr,"Warning: Unable to open solvent \
accessibility file %s for reading\n", safile);
      return(NULL);
   }
   
   /* Read the SA file as PDB                                           */
   if((pdb=ReadPDB(fp, &natom))==NULL)
   {
      fclose(fp);
      fprintf(stderr,"No atoms read from accessibility file: %s\n",
              safile);
      return(NULL);
   }

   /* Close the file and remove if we created it                        */
   fclose(fp);
   if(!SAExists && !KeepSAFile)
      unlink(safile);
   
   return(pdb);
}


/************************************************************************/
/*>void MarkPartners(CLUSINFO *ClusInfo, PDB *pdb, PDB *res, 
                     char *firstres, char *lastres)
   ---------------------------------------------------------------------
   Given a hydrophobic buried in the loop (res), looks for other
   hydrophobics not in the loop which might make s/c--s/c contacts 
   with the loop one

   Note that this has to be done in 2 stages. First we flag those
   residues which are partners, then, at the end, we increment the
   partner count for each flagged residue. This is necessary because
   a residue could be a partner for more than one hydrophobic in the
   loop so simply incrementing the partner count the first time through
   could cause a residue to get incremented more than once.

   08.02.96 Original   By: ACRM
   09.02.96 Added flagged code
*/
void MarkPartners(CLUSINFO *ClusInfo, PDB *pdb, PDB *res, 
                  char *firstres, char *lastres)
{
   int  i,
        j;
   char resspec[16];
   PDB  *partner;

   for(i=0; i<ClusInfo->NRes; i++)
   {
      sprintf(resspec,"%c%d%c", 
              ClusInfo->chain[i],
              ClusInfo->resnum[i],
              ClusInfo->insert[i]);

      /* If this residue is not in the loop                             */
      if(!IsInRange(resspec, firstres, lastres))
      {
         /* Find this residue in the PDB linked list                    */
         if((partner = FindResidue(pdb, 
                                   ClusInfo->chain[i],
                                   ClusInfo->resnum[i],
                                   ClusInfo->insert[i]))!=NULL)
         {
            /* Search for this res type in the Ooi data                 */
            for(j=0; j<MAXRES; j++)
            {
               /* If found...                                           */
               if(!strncmp(partner->resnam, gOoiData[j].resnam, 4))
               {
                  /* If it's hydrophobic                                */
                  if(gOoiData[j].hphob)
                  {
                     /* If they make contact, increment the PartnerCount 
                        variable
                     */
                     if(MakeSCContact(res, partner))
                     {
                        ClusInfo->flagged[i] = TRUE;
                     }
                  }
                  break;
               }
            }
         }
      }
   }   
}


/************************************************************************/
/*>BOOL MakeSCContact(PDB *res1, PDB *res2)
   ----------------------------------------
   Determines whether 2 residues make a sidechain contact

   09.02.96 Original   By: ACRM
*/
BOOL MakeSCContact(PDB *res1, PDB *res2)
{
   PDB *end1,
       *end2,
       *p, *q;

   end1 = FindNextResidue(res1);
   end2 = FindNextResidue(res2);

   for(p=res1; p!=end1; NEXT(p))
   {
      if(strncmp(p->atnam,"N   ",4) &&
         strncmp(p->atnam,"CA  ",4) &&
         strncmp(p->atnam,"C   ",4) &&
         strncmp(p->atnam,"O   ",4))
      {
         for(q=res2; q!=end2; NEXT(q))
         {      
            if(strncmp(q->atnam,"N   ",4) &&
               strncmp(q->atnam,"CA  ",4) &&
               strncmp(q->atnam,"C   ",4) &&
               strncmp(q->atnam,"O   ",4))
            {
               if(DISTSQ(p,q) <= HPHOBCONTDISTSQ)
                  return(TRUE);
            }
         }
      }
   }
   return(FALSE);
}


/************************************************************************/
/*>BOOL MarkHPhob(CLUSINFO *ClusInfo, int clusnum, int nloops,
                  BOOL KeepSA)
   ---------------------------------------------------------------
   Marks residues as key if they are in the loop and are buried 
   hydrophobics or if they are in the framework and interact with the
   key residues in the loop.

   07.02.96 Original   By: ACRM
   08.02.96 Calls library FindResidue()
            Added KeepSA flag
            Fixed bug: Had 2 variables called i (one is now LoopNum)
            Added call to MarkPartners() and support code
   09.02.96 Added code to transfer partners from the flagged array to
            the PartnerCount array
*/
BOOL MarkHPhob(CLUSINFO *ClusInfo, int clusnum, int nloops, 
               BOOL KeepSA)
{
   PDB  *pdb,
        *res1;
   char resspec[16];
   int  i, j,
        LoopNum,
        NRequired;
#ifdef DEBUG
   BOOL L48;
#endif

   NRequired = ClusInfo->NMembers;
   
   /* Zero the counts for each residue                                  */
   for(i=0; i<ClusInfo->NRes; i++)
   {
      ClusInfo->count[i]        = 0;
      ClusInfo->PartnerCount[i] = 0;
   }

   /* Run through all the loops                                         */
   for(LoopNum=0; LoopNum<nloops; LoopNum++)
   {
      /* If we've found a loop in this cluster                          */
      if(gLoopClus[LoopNum].cluster == clusnum)
      {
         if((pdb=ReadPDBAsSA(gLoopClus[LoopNum].filename, KeepSA))==NULL)
         {
            fprintf(stderr,"Warning: Unable to create or read solvent \
accessibility file from %s\n",
                    gLoopClus[LoopNum].filename);
            if(!(--NRequired))
               break;
            continue;
         }

#ifdef DEBUG
         fprintf(stderr,"Marking HPhobs for %s\n",
                 gLoopClus[LoopNum].filename);
#endif

         /* Clear the flags which indicate partner residues             */
         for(i=0; i<ClusInfo->NRes; i++)
            ClusInfo->flagged[i] = FALSE;

         /* Run through the template residues in the CLUSINFO structure,
            seeing if they are in the loop and buried
         */
         for(i=0; i<ClusInfo->NRes; i++)
         {
            sprintf(resspec,"%c%d%c", 
                    ClusInfo->chain[i],
                    ClusInfo->resnum[i],
                    ClusInfo->insert[i]);

            if(IsInRange(resspec,
                         gLoopClus[LoopNum].firstres,
                         gLoopClus[LoopNum].lastres))
            {
               /* Find this residue in the PDB linked list              */
               if((res1 = FindResidue(pdb, 
                                      ClusInfo->chain[i],
                                      ClusInfo->resnum[i],
                                      ClusInfo->insert[i]))!=NULL)
               {
                  /* Search for this res type in the Ooi data           */
                  for(j=0; j<MAXRES; j++)
                  {
                     /* If found...                                     */
                     if(!strncmp(res1->resnam, gOoiData[j].resnam, 4))
                     {
#ifdef DEBUG
                        if(res1->chain[0] == 'L' && res1->resnum == 48)
                        {
                           fprintf(stderr,"\nIn L48:\t");
                           L48 = TRUE;
                        }
                        else
                        {
                           L48 = FALSE;
                        }
#endif
                           
                        /* If it's hydrophobic                          */
                        if(gOoiData[j].hphob)
                        {
#ifdef DEBUG
                           if(L48) fprintf(stderr,"Hydrophobic\t");
#endif
                           /* If the SA is < SACUT, it is buried        */
                           if(res1->bval < SACUT)
                           {
#ifdef DEBUG
                              if(L48) fprintf(stderr,"Buried\t");
#endif
                              (ClusInfo->count[i])++;
#ifdef USE_HPHOB_PARTNERS
                              /* Mark any hydrophobic partner residues  */
                              MarkPartners(ClusInfo, pdb, res1,
                                           gLoopClus[LoopNum].firstres,
                                           gLoopClus[LoopNum].lastres);
#endif
                           }
                        }
                        break;
                     }
                  }
               }
            }
         }

         /* Transfer residues flagged as partners to the partner count  
            (which counts over structures)
         */
         for(i=0; i<ClusInfo->NRes; i++)
         {
            if(ClusInfo->flagged[i])
               (ClusInfo->PartnerCount[i])++;
         }

         /* Free the PDB linked list                                    */
         FREELIST(pdb, PDB);
      }  /* In the correct cluster                                      */
   }  /* For each loop                                                  */

   /* Run through the list of counts marking as key all those which are 
      buried in every loop
   */
   if(NRequired)
   {
      for(i=0; i<ClusInfo->NRes; i++)
      {
         if((ClusInfo->count[i]        == NRequired) ||
            (ClusInfo->PartnerCount[i] == NRequired))
         {
            ClusInfo->key[i] = TRUE;
#ifdef REPORT_REASONS
            Report(ClusInfo, i, (ClusInfo->count[i] == NRequired) ?
                                "Buried Hydrophobic"              :
                                "Partner Hydrophobic");
#endif
            
         }
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL MarkHBonders(CLUSINFO *ClusInfo, int clusnum, int nloops)
   --------------------------------------------------------------
   Identifies residues which make sidechain H-bonds in every loop
   in a cluster

   Allowed HBonds are:
      S/C--ANY if one is in the loop and the other in the f/w
      S/C--B/B if both are in the loop (if USE_LOOP_SM_HBONDS is defined)

   06.02.96 Original   By: ACRM
   08.02.96 Calls library FindResidue()
            Fixed bug: Had 2 variables called i (one is now LoopNum)
            One and only one of the two residues must be in the loop
   09.02.96 Added code to handle S/C--B/B HBonds if both residues
            are in the loop. 
   30.01.09 Initialize some variables
*/
BOOL MarkHBonders(CLUSINFO *ClusInfo, int clusnum, int nloops)
{
   PDB  *pdb,
        *res1 = NULL,
        *res2 = NULL;
   int  i, j,
        natom,
        LoopNum,
        LoopResCount,
        NRequired;
   FILE *fp;
   char resspec[16];

   NRequired = ClusInfo->NMembers;
   
   /* Zero the counts for each residue                                  */
   for(i=0; i<ClusInfo->NRes; i++)
      ClusInfo->count[i]   = 0;

   /* Run through all the loops                                         */
   for(LoopNum=0; LoopNum<nloops; LoopNum++)
   {
      /* If we've found a loop in this cluster                          */
      if(gLoopClus[LoopNum].cluster == clusnum)
      {
#ifdef DEBUG
         fprintf(stderr,"Marking HBonds for %s\n",
                 gLoopClus[LoopNum].filename);
#endif

         /* Open the PDB file                                           */
         if((fp=fopen(gLoopClus[LoopNum].filename, "r"))==NULL)
         {
            fprintf(stderr,"Warning: Unable to open %s for reading\n",
                    gLoopClus[LoopNum].filename);
            if(!(--NRequired))
               break;
            continue;
         }

         /* Read the PDB file and close it                              */
         if((pdb=ReadPDB(fp, &natom))==NULL)
         {
            fprintf(stderr,"No atoms read from PDB file: %s\n",
                    gLoopClus[LoopNum].filename);
            fclose(fp);
            return(FALSE);
         }
         fclose(fp);

         /* Run through the template residues in the CLUSINFO structure,
            seeing if they make sidechain H-bonds to any other template
            residue
         */
         for(i=0; i<ClusInfo->NRes; i++)
         {
            for(j=0; j<ClusInfo->NRes; j++)
            {
               if(i!=j)
               {
                  /* Count how many of these 2 residues are in the loop */
                  LoopResCount = 0;
                  sprintf(resspec,"%c%d%c", 
                          ClusInfo->chain[i],
                          ClusInfo->resnum[i],
                          ClusInfo->insert[i]);
                  if(IsInRange(resspec,
                               gLoopClus[LoopNum].firstres,
                               gLoopClus[LoopNum].lastres))
                     LoopResCount++;
                  sprintf(resspec,"%c%d%c", 
                          ClusInfo->chain[j],
                          ClusInfo->resnum[j],
                          ClusInfo->insert[j]);
                  if(IsInRange(resspec,
                               gLoopClus[LoopNum].firstres,
                               gLoopClus[LoopNum].lastres))
                     LoopResCount++;

                  /* If at least one is in the loop, identify the 2 
                     residues
                  */
                  if(LoopResCount > 0)
                  {
                     res1 = FindResidue(pdb, 
                                        ClusInfo->chain[i],
                                        ClusInfo->resnum[i],
                                        ClusInfo->insert[i]);
                     res2 = FindResidue(pdb, 
                                        ClusInfo->chain[j],
                                        ClusInfo->resnum[j],
                                        ClusInfo->insert[j]);
                  }

                  if(res1!=NULL && res2!=NULL)
                  {
                     if(LoopResCount == 1)
                     {
                        /* If it's just one of them, then test for HBond 
                           between sidechain and anything

                           If the first makes a s/c HBond to the second, 
                           then increment its count
                        */
                        if(IsHBonded(res1, res2, HBOND_SIDECHAIN))
                        {
                           (ClusInfo->count[i])++;
                           /* We now break out of the inner loop to stop
                              us counting more than one HBond involving i
                           */
                           break;
                        }                           
                     }
#ifdef USE_LOOP_SM_HBONDS
                     else if(LoopResCount == 2)
                     {
                        /* If both residues are in the loop, check for 
                           HBond between sidechain and backbone

                           If the first makes a s/c HBond to the second, 
                           then increment its count
                        */
                        if(IsHBonded(res1, res2, HBOND_SB))
                        {
                           (ClusInfo->count[i])++;
                           /* We now break out of the inner loop to stop
                              us counting more than one HBond involving i
                           */
                           break;
                        }
                     }
#endif
                  }
               }
            }
         }
         
         /* Free the PDB linked list                                    */
         FREELIST(pdb, PDB);
      }  /* In the correct cluster                                      */
   }  /* For each loop                                                  */

   /* Run through the list of counts marking all those which do
      make a s/c HBond in every loop as key
   */
   if(NRequired)
   {
      for(i=0; i<ClusInfo->NRes; i++)
      {
         if(ClusInfo->count[i] == ClusInfo->NMembers)
         {
            ClusInfo->key[i] = TRUE;
#ifdef REPORT_REASONS
            Report(ClusInfo, i, "Conserved Hydrogen Bond");
#endif
         }         
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>SDRLIST *InSDRList(SDRLIST *sdrlist, char chain, int resnum,
                      char insert)
   ------------------------------------------------------------
   Sees whether a residue is in the linked list of SDRs
   13.02.96 Original   By: ACRM
*/
SDRLIST *InSDRList(SDRLIST *sdrlist, char chain, int resnum, char insert)
{
   SDRLIST *s;
   
   for(s=sdrlist; s!=NULL; NEXT(s))
   {
      if((s->resnum == resnum) &&
         (s->chain  == chain)  &&
         (s->insert == insert))
         return(s);
   }
   
   return(NULL);
}


/************************************************************************/
/*>BOOL FillSDRsForCluster(SDRLIST *sdrlist, int clusnum, int nloops)
   ------------------------------------------------------------------
   For a given cluster, runs through the PDB files and stores the
   amino acids seen at each key residue position.

   13.02.96 Original   By: ACRM
*/
BOOL FillSDRsForCluster(SDRLIST *sdrlist, int clusnum, int nloops)
{
   int     LoopNum,
           natom;
   PDB     *pdb,
           *p;
   FILE    *fp;
   SDRLIST *s;
   char    res;
   BOOL    Found;

   /* Clear the number of observed residues for each SDR position       */
   for(s=sdrlist; s!=NULL; NEXT(s))
      s->nobsres = 0;

   /* Run through all the loops                                         */
   for(LoopNum=0; LoopNum<nloops; LoopNum++)
   {
      /* If we've found a loop in this cluster                          */
      if(gLoopClus[LoopNum].cluster == clusnum)
      {
         /* Open the PDB file                                           */
         if((fp=fopen(gLoopClus[LoopNum].filename, "r"))==NULL)
         {
            fprintf(stderr,"Warning: Unable to open %s for reading\n",
                    gLoopClus[LoopNum].filename);
            continue;
         }

         /* Read the PDB file and close it                              */
         if((pdb=ReadPDB(fp, &natom))==NULL)
         {
            fprintf(stderr,"No atoms read from PDB file: %s\n",
                    gLoopClus[LoopNum].filename);
            fclose(fp);
            return(FALSE);
         }
         fclose(fp);

         /* Run through the PDB file, checking if each residue is in the
            key residue list
         */
         for(p=pdb; p!=NULL; NEXT(p))
         {
            if(!strncmp(p->atnam,"CA  ",4))
            {
               /* If it is in the SDR list                              */
               if((s=InSDRList(sdrlist, p->chain[0], p->resnum, 
                               p->insert[0]))!=NULL)
               {
                  /* See if it's in the list of observed residues and,
                     if not, add it.
                  */
                  res = throne(p->resnam);
                  TESTINARRAY(s->obsres, s->nobsres, res, Found);
                  if(!Found)
                     s->obsres[(s->nobsres)++] = res;
               }
            }
         }

         FREELIST(pdb, PDB);
      }  /* Correct cluster                                             */
   }  /* For each loop                                                  */

   /* Terminate all SDR lists as strings                                */
   for(s=sdrlist; s!=NULL; NEXT(s))
      s->obsres[s->nobsres] = '\0';

   return(TRUE);
}


/************************************************************************/
/*>void ReportUnifiedSDRs(FILE *out, int nclus, int nloops)
   ------------------------------------------------
   For each cluster store the observed residues at each key position
   within this cluster's loop or framework and at positions in the 
   framework defined as key for all other clusters (of any loop length) 
   with at least MINCLUSSIZE members. Also add positions from smaller 
   clusters of the same SDR length. These are flagged so that 
   non-informative extras positions can be removed.

   06.02.96 Original   By: ACRM
   15.02.96 Added code to unify on loops of same length
   16.02.96 Added elimination of non-informative SDRs added on length
            and notification of rogue clusters
   30.01.09 Initialize some variables
*/
BOOL ReportUnifiedSDRs(FILE *out, int nclus, int nloops)
{
   int     clus,
           clus2,
           res,
           i;
   SDRLIST *sdrlist = NULL,
           *s = NULL;
   char    resspec[16],
           firstres[16],
           lastres[16];
   BOOL    UseCluster,
           InFramework;

   
   fprintf(out,"%cObserved residues for each cluster at unified SDR \
positions:\n",(char)12);
   
   /* For each cluster, build a linked list of the SDR positions in this
      cluster and all f/w SDR positions in all the other clusters
   */
   for(clus=0; clus<nclus; clus++)
   {
      /* Find an example from this cluster to define the firstres and
         lastres
      */
      for(i=0; i<nloops; i++)
      {
         if(gLoopClus[i].cluster == clus+1)
         {
            strcpy(firstres, gLoopClus[i].firstres);
            strcpy(lastres,  gLoopClus[i].lastres);
            break;
         }
      }

      /* Store SDR positions for this cluster from key residues for
         this cluster's sequence template                              
      */
      for(res=0; res<gClusInfo[clus].NRes; res++)
      {
         /* If it's marked as a key residue                             */
         if(gClusInfo[clus].key[res])
         {
            /* If it's not already in the SDR list                      */
            if(InSDRList(sdrlist,
                         gClusInfo[clus].chain[res],
                         gClusInfo[clus].resnum[res],
                         gClusInfo[clus].insert[res])==NULL)
            {
               /* Allocate space in the linked list                     */
               if(sdrlist==NULL)
               {
                  INIT(sdrlist,SDRLIST);
                  s = sdrlist;
                  /* Store this pointer in the ClusInfo structure       */
                  gClusInfo[clus].sdrlist = sdrlist;
               }
               else
               {
                  ALLOCNEXT(s,SDRLIST);
               }
               if(s==NULL)
               {
                  fprintf(stderr,"No memory for SDR list\n");
                  return(FALSE);
               }
               
               /* Insert the data                                       */
               s->chain    = gClusInfo[clus].chain[res];
               s->resnum   = gClusInfo[clus].resnum[res];
               s->insert   = gClusInfo[clus].insert[res];
               s->nobsres  = 0;
               s->onlength = OL_FALSE;

               sprintf(resspec,"%c%d%c",
                       gClusInfo[clus].chain[res],
                       gClusInfo[clus].resnum[res],
                       gClusInfo[clus].insert[res]);
                  
               /* If it's not in the loop (i.e. is f/w)                 */
               if(IsInRange(resspec,firstres,lastres))
                  s->position = POS_LOOP;
               else
                  s->position = POS_CONTACT;
            }  /* If it's not already in the SDR list                   */
         }  /* If it's a key residue                                    */
      }  /* For each residue in sequence template                       */

      /* Unify the SDR list by adding key residues for other clusters.
         This is done for all large clusters and/or for all clusters
         with the same loop length depending on the defines
      */
      for(clus2=0; clus2<nclus; clus2++)
      {
         UseCluster = FALSE;
         if(clus2 != clus)
         {
#ifdef UNIFY_ON_LARGE_CLUSTER
            if(gClusInfo[clus2].NMembers >= MINCLUSSIZE)
               UseCluster = TRUE;
#endif
#ifdef UNIFY_ON_LENGTH
            if(gClusInfo[clus].length == gClusInfo[clus2].length)
               UseCluster = TRUE;
#endif
         }

         /* If it's not the current cluster and it is big enough        */
         if(UseCluster)
         {
            /* Run through the sequence template                        */
            for(res=0; res<gClusInfo[clus2].NRes; res++)
            {
               /* If it's marked as a key residue                       */
               if(gClusInfo[clus2].key[res])
               {
                  sprintf(resspec,"%c%d%c",
                          gClusInfo[clus2].chain[res],
                          gClusInfo[clus2].resnum[res],
                          gClusInfo[clus2].insert[res]);
                  
                  /* If it's not in the loop (i.e. is f/w) or this loop
                     is of the same length, add it to the SDR list
                  */
                  InFramework = (!IsInRange(resspec,firstres,lastres));
                  
                  if(InFramework ||
                     (gClusInfo[clus].length == gClusInfo[clus2].length))
                  {
                     /* If it's not already in the SDR list             */
                     if(InSDRList(sdrlist,
                                  gClusInfo[clus2].chain[res],
                                  gClusInfo[clus2].resnum[res],
                                  gClusInfo[clus2].insert[res])==NULL)
                     {
                        /* Allocate space in the linked list            */
                        if(sdrlist==NULL)
                        {
                           INIT(sdrlist,SDRLIST);
                           s = sdrlist;
                        }
                        else
                        {
                           ALLOCNEXT(s,SDRLIST);
                        }
                        if(s==NULL)
                        {
                           fprintf(stderr,"No memory for SDR list\n");
                           return(FALSE);
                        }
               
                        /* Insert the data                              */
                        s->chain    = gClusInfo[clus2].chain[res];
                        s->resnum   = gClusInfo[clus2].resnum[res];
                        s->insert   = gClusInfo[clus2].insert[res];
                        s->nobsres  = 0;
                        s->position = InFramework?POS_NOCONTACT:POS_LOOP;
                        s->onlength = (gClusInfo[clus2].NMembers <
                                       MINCLUSSIZE) ? 
                                       OL_ONLENGTH  : 
                                       OL_FALSE;
                     }  /* If it's not already in the SDR list          */
                  }  /* If it's not in the loop                         */
               }  /* If it's a key residue                              */
            }  /* For each residue in sequence template                 */
         }  /* Not the current cluster and has enough members           */
      }  /* For each cluster (looking at f/w residues)                  */

      /* Fill in the observed residues from the PDB files               */
      FillSDRsForCluster(sdrlist, clus+1, nloops);

      /* We don't free the linked list of SDRs as we keep the pointer in
         the ClusInfo structure. However, we do set sdrlist to NULL so
         we will initialise the linked list for the next cluster
      */
      sdrlist = NULL;
   }  /* For each cluster (main loop)                                   */

#ifdef EXCLUDE_NONINFORM
   FlagNonInformativeSDRs(nclus);
#endif

   FlagRogueClusters(nclus, nloops);

   /* Print out the allowed residue lists                               */
   for(clus=0; clus<nclus; clus++)
   {
      char buffer[MAXBUFF];
      
      if(gClusInfo[clus].rogue)
         sprintf(buffer,", ROGUE - Matches cluster %d",
                 gClusInfo[clus].rogue);
      else
         buffer[0] = '\0';
      
      /* Print a header for this cluster                                */
      fprintf(out,"\nCLUSTER %d (Length = %d, Members = %d%s)\n",
              clus+1,
              gClusInfo[clus].length,
              gClusInfo[clus].NMembers,
              buffer);

      PrintSDRList(out, gClusInfo[clus].sdrlist);
      FREELIST(gClusInfo[clus].sdrlist, SDRLIST);
      gClusInfo[clus].sdrlist = NULL;
   }

   return(TRUE);
}

      
/************************************************************************/
/*>void PrintSDRList(FILE *out, SDRLIST *sdrlist)
   ----------------------------------------------
   Prints the SDR list in order sorted by residue number (note that we
   don't bother to sort by chain label or by insert...)

   13.02.96 Original   By: ACRM
   16.02.96 Doesn't print if the residue has been flagged as 
            non-informative
*/
void PrintSDRList(FILE *out, SDRLIST *sdrlist)
{
   SDRLIST *s,
           **sdrlistidx;
   int     NumSDR,
           *resnum,
           *idx,
           i;
   BOOL    ok = FALSE;
   
   /* First count the number of SDRs                                    */
   for(s=sdrlist, NumSDR=0; s!=NULL; NEXT(s))
      NumSDR++;

   /* Allocate array to index the linked list                           */
   if((sdrlistidx=(SDRLIST **)malloc(NumSDR * sizeof(SDRLIST *)))!=NULL)
   {
      /* Create an array to store the residue numbers                   */
      if((resnum=(int *)malloc(NumSDR * sizeof(int)))!=NULL)
      {
         /* Create a sort index for the res num array                   */
         if((idx=(int *)malloc(NumSDR * sizeof(int)))!=NULL)
         {
            /* All allocations succeeded                                */
            ok = TRUE;
         
            /* Create index into the linked list                        */
            for(s=sdrlist, i=0; s!=NULL; NEXT(s), i++)
               sdrlistidx[i] = s;
            
            /* Copy the residue numbers into the resnum array           */
            for(i=0; i<NumSDR; i++)
               resnum[i] = sdrlistidx[i]->resnum;

            /* Index sort the resnum array                              */
            indexint(NumSDR, resnum, idx);

            /* Print out the SDRs in index sorted order                 */
            for(i=0; i<NumSDR; i++)
            {
               if(sdrlistidx[idx[i]]->onlength != OL_DELETABLE)
               {
                  fprintf(out,"%c %4d %c : %-20s (%s%s)\n",
                          sdrlistidx[idx[i]]->chain, 
                          sdrlistidx[idx[i]]->resnum, 
                          sdrlistidx[idx[i]]->insert, 
                          (sdrlistidx[idx[i]]->obsres[0] ?
                           sdrlistidx[idx[i]]->obsres    :
                           "-"),
                          (sdrlistidx[idx[i]]->position==POS_NOCONTACT  ?
                           "No contact"                                 :
                           (sdrlistidx[idx[i]]->position==POS_CONTACT   ?
                            "Makes contact"                             :
                            "In loop")),
                          (sdrlistidx[idx[i]]->onlength                 ?
                           ", Added on length"                          :
                           ""));
               }
            }
            free(idx);
         }
         free(resnum);
      }
      free(sdrlistidx);
   }

   /* If allocations failed, just print in the current order            */
   if(!ok)
   {
      for(s=sdrlist; s!=NULL; NEXT(s))
      {
         fprintf(out,"%c %4d %c : %-20s (%s%s)\n",
                 s->chain, s->resnum, s->insert, 
                 (s->obsres[0]?s->obsres:"-"),
                 (s->position==POS_NOCONTACT    ?
                  "No contact"                  :
                  (s->position==POS_CONTACT     ?
                   "Makes contact"              :
                   "In loop")),
                 (s->onlength?", Added on length" : ""));
      }
   }
}


/************************************************************************/
/*>void indexint(int n, int *arrin, int *indx)
   ---------------------------------------------
   Input:  n      int      Number of elements in array
           arrin  *int    Array to be indexed
   Output: indx   *int     Index array
   
   Index a int array by Heapsort.
   
   03.06.90 Original
   01.06.92 ANSIed and autodoc'd
   13.02.96 Tideid and changed from double to int
*/
void indexint(int n, int *arrin, int *indx)
{
   int i, j, l, ir, 
       indxt;
   int q;

   /* Special cases where we don't need to do any sorting!              */
   if(n==0)
      return;
   if(n==1)
   {
      indx[0] = 0;
      return;
   }
   
   /* Modify the array starts so we can count from 1                    */
   indx--;
   arrin--;
   
   for(j=1; j<=n; j++)
      indx[j] = j;

   l  = n/2+1;
   ir = n;
   
   for(;;)
   {
      if(l>1)
      {
         indxt = indx[--l];
         q     = arrin[indxt];
      }
      else
      {
         indxt    = indx[ir];
         q        = arrin[indxt];
         indx[ir] = indx[1];

         if((--ir) == 1)
         {
            indx[1] = indxt;

            /* Correct the offsets, so they run from 0                  */
            for(j=1; j<=n; j++)
               indx[j]--;

            return;
         }
      }

      i = l;
      j = l+l;
      while(j <= ir)
      {
         if(j < ir)
         {
            if(arrin[indx[j]] < arrin[indx[j+1]])
               j++;
         }
         if(q < arrin[indx[j]])
         {
            indx[i] = indx[j];
            i       = j;
            j      += j;
         }
         else
         {
            j = ir+1;
         }
      }
      indx[i] = indxt;
   }
}


/************************************************************************/
/*>void FlagNonInformativeSDRs(int nclus)
   --------------------------------------
*/
void FlagNonInformativeSDRs(int nclus)
{
   int     clus1,
           clus2,
           MaxAllowed,
           ClusWithMaxAllowed;
   SDRLIST *sdrlist1,
           *sdrlist2,
           *s1,
           *s2,
           *MaxAllowedSDR;
   BOOL    AddedValue;
   
   /* For each cluster                                                  */
   for(clus1=0; clus1<nclus; clus1++)
   {
      sdrlist1 = gClusInfo[clus1].sdrlist;
      
      /* For each SDR                                                   */
      for(s1=sdrlist1; s1!=NULL; NEXT(s1))
      {
         MaxAllowed         = s1->nobsres;
         ClusWithMaxAllowed = clus1;
         MaxAllowedSDR      = s1;
         
         /* If this SDR was added on loop length match alone            */
         if(s1->onlength)
         {
            /* See which cluster of this length has the most allowed
               amino acids for this SDR
            */

            /* For each other cluster of the same loop length           */
            for(clus2=0; clus2<nclus; clus2++)
            {
               sdrlist2 = gClusInfo[clus2].sdrlist;

               if((clus1 != clus2) &&
                  (gClusInfo[clus1].length == gClusInfo[clus2].length))
               {
                  /* Find the equivalent SDR in this second cluster     */
                  for(s2=sdrlist2; s2!=NULL; NEXT(s2))
                  {
                     if((s2->resnum == s1->resnum) &&
                        (s2->chain  == s1->chain)  &&
                        (s2->insert == s1->insert))
                     {
                        if(s2->nobsres > MaxAllowed)
                        {
                           MaxAllowed = s2->nobsres;
                           ClusWithMaxAllowed = clus2;
                           MaxAllowedSDR = s2;
                        }
                        break;
                     }
                  }
               }
            }

            /* Run through the clusters with the same length loop a
               second time to see if any of them gives "added value"
               i.e. all the allowed residues are excluded from the
               allowed residues in the biggest group.
            */
            AddedValue = FALSE;
            for(clus2=0; clus2<nclus; clus2++)
            {
               sdrlist2 = gClusInfo[clus2].sdrlist;

               /* Note we check it's not the cluster with the max allowed
                  this time rather than clus1
               */
               if((clus2 != ClusWithMaxAllowed) &&
                  (gClusInfo[clus2].length == 
                   gClusInfo[ClusWithMaxAllowed].length))
               {
                  /* Find the equivalent SDR in this second cluster     */
                  for(s2=sdrlist2; s2!=NULL; NEXT(s2))
                  {
                     if((s2->resnum == s1->resnum) &&
                        (s2->chain  == s1->chain)  &&
                        (s2->insert == s1->insert))
                     {
                        if(ValueIsAdded(MaxAllowedSDR, s2))
                        {
                           AddedValue = TRUE;
                           break;
                        }
                     }
                  }
               }
            }

            /* If no value was added, mark this SDR from the main cluster
               loop as deletable
            */
            if(!AddedValue)
               s1->onlength = OL_DELETABLE;
         }  /* Was added on length                                      */
      }  /* For each SDR                                                */
   }  /* For each cluster1                                              */
}


/************************************************************************/
/*>BOOL ValueIsAdded(SDRLIST *s1, SDRLIST *s2)
   -------------------------------------------
   Sees if any one of the amino acid types listed in s2 occur in s1.
   If any one does, then returns FALSE as no value is added; if none
   does, returns TRUE as value is added

   16.02.96 Original   By: ACRM
*/
BOOL ValueIsAdded(SDRLIST *s1, SDRLIST *s2)
{
   int i;
   
   for(i=0; i<s2->nobsres; i++)
   {
      if(strchr(s1->obsres, (s2->obsres)[i]))
         return(FALSE);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void FlagRogueClusters(int nclus, int nloops)
   ---------------------------------------------
   Flag `rogue' clusters i.e. small clusters which have no distinguishing
   characteristics at the SDR positions

   16.02.96 Original   By: ACRM
   29.04.96 Now also searches for clusters which are rogues wrt clusters
            other than the largest one
*/
void FlagRogueClusters(int nclus, int nloops)
{
   int     clus1, clus2, clus3,
           clusa, clusb,
           LoopLength,
           LargestClus,
           LargestClusSize,
           MinLoopLength = 0,
           MaxLoopLength = 0;

   /* Assume all clusters are not rogues and find min and max loop
      lengths
   */
   MinLoopLength = MaxLoopLength = gClusInfo[0].length;
   for(clus1=0; clus1<nclus; clus1++)
   {
      gClusInfo[clus1].rogue   = 0;
      if(gClusInfo[clus1].length > MaxLoopLength)
         MaxLoopLength = gClusInfo[clus1].length;
      if(gClusInfo[clus1].length < MinLoopLength)
         MinLoopLength = gClusInfo[clus1].length;
   }

   /* For each possible loop length                                     */
   for(LoopLength=MinLoopLength; LoopLength<=MaxLoopLength; LoopLength++)
   {
      LargestClus     = (-1);
      LargestClusSize = (-1);
      
      /* For each cluster with this loop length, find the largest cluster
      */
      for(clus1=0; clus1<nclus; clus1++)
      {
         if(gClusInfo[clus1].length == LoopLength)
         {
            if(gClusInfo[clus1].NMembers > LargestClusSize)
            {
               LargestClusSize = gClusInfo[clus1].NMembers;
               LargestClus     = clus1;
            }
         }
      }

      /* If there were some clusters with this loop length              */
      if(LargestClus != (-1))
      {
         /* Run through all the clusters of this length and see if any
            are rogues c.f. the largest cluster
         */
         for(clus1=0; clus1<nclus; clus1++)
         {
            if((clus1 != LargestClus) &&
               (gClusInfo[clus1].length == LoopLength))
            {
               if(IsRogue(clus1, LargestClus))
               {
                  gClusInfo[clus1].rogue = LargestClus+1;
               }
            }
         }

         /* Now see if there are any rogues vs. non-largest cluster     

            For each cluster of this length which is not already flagged
            as a rogue
         */
         for(clus1=0; clus1<nclus; clus1++)
         {
            if((gClusInfo[clus1].length == LoopLength) &&
               (!gClusInfo[clus1].rogue))
            {
               /* For each other cluster of the same length             */
               for(clus2=clus1+1; clus2<nclus; clus2++)
               {
                  if(gClusInfo[clus2].length == LoopLength)
                  {
                     /* Set clusa to the larger and clusb to the smaller*/
                     clusa = clus1;
                     clusb = clus2;
                     if(gClusInfo[clusa].NMembers <
                        gClusInfo[clusb].NMembers)
                     {
                        int temp = clusa;
                        clusa    = clusb;
                        clusb    = temp;
                     }
                        
                     if(IsRogue(clusb, clusa))
                     {
                        /* Update any rogue clusters for which clusb was
                           the "parent" and switch the parent to clusa

                           14.06.07 Did say
                              gClusInfo[clus3].rogue -= clusa+1;
                        */
                        for(clus3=0; clus3<nclus; clus3++)
                        {
                           if(gClusInfo[clus3].rogue == clusb+1)
                              gClusInfo[clus3].rogue = clusa+1;
                        }
                        
                        /* Flag clusb as rogue, setting clusa as its
                           parent
                        */
                        gClusInfo[clusb].rogue = clusa+1;
                     }
                  }
               }
            }
         }  /* Finding rogues vs. non-largest cluster                   */
      }  /* Have some clusters for this loop length                     */
   }  /* For each loop length                                           */
}


/************************************************************************/
/*>BOOL IsRogue(int clus, int LargestClus)
   ---------------------------------------
   Checks whether a cluster (clus) is a rogue compared with the
   LargestClus cluster i.e. it has no distinguishing characteristics at 
   the SDR positions

   16.02.96 Original   By: ACRM
*/
BOOL IsRogue(int clus, int LargestClus)
{
   SDRLIST *sdrlist1,
           *sdrlist2,
           *s1, *s2;

   sdrlist1 = gClusInfo[LargestClus].sdrlist;
   sdrlist2 = gClusInfo[clus].sdrlist;

   for(s1=sdrlist1; s1!=NULL; NEXT(s1))
   {
      for(s2=sdrlist2; s2!=NULL; NEXT(s2))
      {
         if((s2->resnum == s1->resnum) &&
            (s2->chain  == s1->chain)  &&
            (s2->insert == s1->insert))
         {
            if(ValueIsAdded(s1, s2))
            {
               return(FALSE);
            }
         }
      }
   }

   /* There were no `added-value' SDRs, so this is a rogue              */
   return(TRUE);
}






/************************************************************************/
/*>BOOL IsCisProline(CLUSINFO *ClusInfo, int clusnum, int resoffset, 
                     int nloops)
   -----------------------------------------------------------------
   Determines whether a residue from the first PDB file in a specified
   cluster is a cis proline.

   22.03.96 Original   By: ACRM
*/
BOOL IsCisProline(CLUSINFO *ClusInfo, int clusnum, int resoffset, 
                  int nloops)
{
   PDB  *pdb,
        *prev,
        *ResPro,
        *ResPrev,
        *ResNext,
        *p,
        *CA1 = NULL,
        *C1  = NULL,
        *N2  = NULL,
        *CA2 = NULL;
   int  natom,
        LoopNum;
   FILE *fp;
   REAL angle;
   

   /* Run through all the loops                                         */
   for(LoopNum=0; LoopNum<nloops; LoopNum++)
   {
      /* If we've found a loop in this cluster                          */
      if(gLoopClus[LoopNum].cluster == clusnum)
      {
#ifdef DEBUG
         fprintf(stderr,"Marking cis-prolines for %s\n",
                 gLoopClus[LoopNum].filename);
#endif

         /* Open the PDB file                                           */
         if((fp=fopen(gLoopClus[LoopNum].filename, "r"))==NULL)
         {
            fprintf(stderr,"Warning: Unable to open %s for reading\n",
                    gLoopClus[LoopNum].filename);
            continue;
         }

         /* Read the PDB file and close it                              */
         if((pdb=ReadPDB(fp, &natom))==NULL)
         {
            fprintf(stderr,"No atoms read from PDB file: %s\n",
                    gLoopClus[LoopNum].filename);
            fclose(fp);
            return(FALSE);
         }
         fclose(fp);

         /* Find the PDB pointer for the resoffset residue from the 
            ClusInfo structure
         */
         ResPro = FindResidue(pdb,
                              ClusInfo->chain[resoffset],
                              ClusInfo->resnum[resoffset],
                              ClusInfo->insert[resoffset]);

         /* Find the previous residue. This is a little messy since the
            PDB linked list is not doubly linked
         */
         for(prev = pdb; prev->next != ResPro; NEXT(prev));
         ResPrev = FindResidue(pdb,
                               prev->chain[0],
                               prev->resnum,
                               prev->insert[0]);

         /* Find the next residue                                       */
         ResNext = FindNextResidue(ResPro);

         /* Find the atoms describing the omega torsion angle           */
         for(p=ResPrev; p!=ResPro; NEXT(p))
         {
            if(!strncmp(p->atnam,"CA  ",4))
               CA1 = p;
            if(!strncmp(p->atnam,"C   ",4))
               C1  = p;
         }
         for(p=ResPro; p!=ResNext; NEXT(p))
         {
            if(!strncmp(p->atnam,"N   ",4))
               N2  = p;
            if(!strncmp(p->atnam,"CA  ",4))
               CA2 = p;
         }
         
         if(CA1==NULL || C1==NULL || N2==NULL || CA2==NULL)
         {
            fprintf(stderr,"Warning: Missing atom around possible \
cis-proline in %s %c%d%c\n",
                    gLoopClus[LoopNum].filename,
                    ClusInfo->chain[resoffset],
                    ClusInfo->resnum[resoffset],
                    ClusInfo->insert[resoffset]);
         }
         else
         {
            angle = phi(CA1->x, CA1->y, CA1->z,
                        C1->x,  C1->y,  C1->z,
                        N2->x,  N2->y,  N2->z,
                        CA2->x, CA2->y, CA2->z);
            FREELIST(pdb, PDB);
            if((angle > (-PI/2.0)) && (angle < (PI/2.0)))
               return(TRUE);
            else
               return(FALSE);
         }

         /* Free the PDB linked list and see if there's another struct. */
         FREELIST(pdb, PDB);
      }  /* In the correct cluster                                      */
   }  /* For each loop                                                  */

   return(FALSE);
}


