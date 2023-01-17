/*************************************************************************

   Program:    clan/ficl
   File:       acaca.c
   
   Version:    V3.6
   Date:       09.01.96
   Function:   Perform cluster analysis on loop conformations
   
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
   EMail:      martin@biochem.ucl.ac.uk
               
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
   V3.3  04.10.95 Skipped
   V3.4  10.09.95 Skipped
   V3.5  06.11.95 Skipped
   V3.6  09.01.96 Frees memory not needed when not doing critical residues
                  Truncated structures were causing the whole data list to
                  be freed.
   V3.6a 30.01.09 Compile cleanups

*************************************************************************/
/* Includes
*/
#include "acaca.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>BOOL SetClusterMethod(char *method)
   -----------------------------------
   Input:   char *method         Clustering method string
   Returns: BOOL                 Success? (FALSE if unknown method)
   Globals: int  gClusterMethod  Set to specified method

   Set the clustering method global variable based on the supplied
   text.

   27.06.95 Original   By: ACRM
*/
BOOL SetClusterMethod(char *method)
{
   if(!blUpstrncmp(method,"WAR",3) || method[0] == '1')
      gClusterMethod = 1;
   else if(!blUpstrncmp(method,"SIN",3) || method[0] == '2')
      gClusterMethod = 2;
   else if(!blUpstrncmp(method,"COM",3) || method[0] == '3')
      gClusterMethod = 3;
   else if(!blUpstrncmp(method,"AVE",3) || 
           !blUpstrncmp(method,"GRO",3) || method[0] == '4')
      gClusterMethod = 4;
   else if(!blUpstrncmp(method,"MCQ",3) || method[0] == '5')
      gClusterMethod = 5;
   else if(!blUpstrncmp(method,"MED",3) || 
           !blUpstrncmp(method,"GOW",3) || method[0] == '6')
      gClusterMethod = 6;
   else if(!blUpstrncmp(method,"CEN",3) || method[0] == '7')
      gClusterMethod = 7;
   else
   {
      fprintf(stderr,"Unknown clustering method: %s\n",method);
      return(FALSE);
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL SetOutputFile(char *filename)
   ----------------------------------
   Input:   char  *filename         Output filename
   Returns: BOOL                    Success?
   Globals: FILE  *gOutfp           Output file pointer

   Open an output file other then stdout

   27.06.95 Original   By: ACRM
*/
BOOL SetOutputFile(char *filename)
{
   if((gOutfp = fopen(filename,"w"))==NULL)
   {
      fprintf(stderr,"Unable to open output file: %s\n",filename);
      return(FALSE);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL HandleLoopSpec(char *filename, char *start, char *end,
                       BOOL CATorsions, BOOL Verbose)
   -----------------------------------------------------------
   Input:   char   *filename      PDB structure filename
            char   *start         Loop start residue spec
            char   *end           Loop end residue spec
            BOOL   CATorsions     Do CA rather than true torsions
            BOOL   Verbose        Print the filename being processed?
   Returns: BOOL                  Success of opening file and allocating
                                  memory

   Takes a PDB/loop specification and reads the PDB file, selects the
   CAs and calls routine to store the pseudo-torsions and related data
   in global linked list.

   27.06.95 Original   By: ACRM
   06.07.95 No longer frees pdbca list as this is stored in the
            global data list
   26.07.95 Frees sel[] array items correctly
   31.07.95 Added Verbose flag
   08.08.95 StoreTorsions() now stores the whole PDB linked list as
            well as the selected atoms linked list. Consequently, the 
            main linked list is no longer freed.
   15.08.95 Initialise the sel[] array
   09.01.95 Added !gDoCritRes handling; we can free up the PDB linked
            lists
*/
BOOL HandleLoopSpec(char *filename, char *start, char *end, 
                    BOOL CATorsions, BOOL Verbose)
{
   FILE *fp    = NULL;
   PDB  *pdb   = NULL, 
        *pdbca = NULL,
        *p_start,
        *p_end;
   int  natom,
        resnum1, resnum2;
   char chain1,  chain2, 
        insert1, insert2;
   char *sel[3];
   BOOL retval = TRUE;

   /* Initialse so checks on free are OK                                */
   sel[0] = sel[1] = sel[2] = NULL;
   
   /* Open the specified PDB file                                       */
   if((fp=fopen(filename,"r"))==NULL)
   {
      fprintf(stderr,"Unable to open file: %s\n",filename);
      retval = FALSE;
   }
   else
   {
      /* Read in the file                                               */
      if((pdb=blReadPDBAtoms(fp,&natom))==NULL)
      {
         fprintf(stderr,"Unable to read atoms from file: %s\n",filename);
         retval = FALSE;
      }
      else
      {
         if(Verbose)
            fprintf(stderr,"Processing file: %s\n",filename);
         
         /* Select just the CAs                                         */
         SELECT(sel[0],"CA  ");
         SELECT(sel[1],"N   ");
         SELECT(sel[2],"C   ");
         if((sel[0] != NULL) && 
            (sel[1] != NULL) && 
            (sel[2] != NULL))
         {
            /* Parse the resspecs for start and end                     */
            blParseResSpec(start, &chain1, &resnum1, &insert1);
            blParseResSpec(end,   &chain2, &resnum2, &insert2);
               
            if(CATorsions)
            {
               /* Do CA-pseudo torsions                                 */
               if((pdbca = blSelectAtomsPDBAsCopy(pdb, 1, sel, &natom))==NULL)
               {
                  fprintf(stderr,"Unable to select CA atoms\n");
                  retval = FALSE;
               }
               else
               {
                  /* Find one back from the startres and 2 on from last 
                     res
                  */
                  if(!FindCAResidues(pdbca,
                                     chain1,resnum1,insert1,
                                     chain2,resnum2,insert2,
                                     &p_start,&p_end))
                  {
                     retval = FALSE;
                  }
                  else
                  {
                     /* Calculate and store the CA pseudo-torsions      */
                     if(!StoreTorsions(pdb,pdbca,p_start,p_end,filename,
                                       start,end))
                     {
                        retval = FALSE;
                     }
                  }
               }
            }
            else    /* Do true torsions                                 */
            {
               if((pdbca = blSelectAtomsPDBAsCopy(pdb, 3, sel, &natom))==NULL)
               {
                  fprintf(stderr,"Unable to select backbone atoms\n");
                  retval = FALSE;
               }
               else
               {
                  if(!FindBBResidues(pdbca,
                                     chain1,resnum1,insert1,
                                     chain2,resnum2,insert2,
                                     &p_start,&p_end))
                  {
                     retval = FALSE;
                  }
                  else
                  {
                     /* Calculate and store the true torsions           */
                     if(!StoreTorsions(pdb,pdbca,p_start,p_end,filename,
                                       start,end))
                     {
                        retval = FALSE;
                     }
                  }
               }
            }
         }
         else
         {
            fprintf(stderr,"No memory for selection list\n");
            retval = FALSE;
         }
      }
   }

   if(fp     != NULL) fclose(fp);
   if(sel[0] != NULL) free(sel[0]);
   if(sel[1] != NULL) free(sel[1]);
   if(sel[2] != NULL) free(sel[2]);
   if(!gDoCritRes)
   {
      if(pdb   != NULL) FREELIST(pdb,   PDB);
      /* Note that we can't free the CA linked list if we're doing
         postclustering

         if(pdbca != NULL) FREELIST(pdbca, PDB);
      */
   }
   
   return(retval);
}


/************************************************************************/
/*>BOOL FindCAResidues(PDB *pdbca, char chain1, int resnum1, char insert1,
                       char chain2, int resnum2, char insert2,
                       PDB **pp_start, PDB **pp_end)
   ---------------------------------------------------------------------
   Input:   PDB  *pdbca      CA PDB linked list
            char chain1      Start of loop chain spec
            int  resnum1     Start of loop residue number
            char insert1     Start of loop insert code
            char chain2      End of loop chain spec
            int  resnum2     End of loop residue number
            char insert2     End of loop insert code
   Output:  PDB  **pp_start  Pointer to atom before start of loop in 
                             linked list
            PDB  **pp_end    Pointer to atom before end of loop in 
                             linked list
   Returns: BOOL             Success of finding specfied residues

   Finds pointers to the atoms before the specified start and end 
   residues.

   N.B. This routine assumes that only CA atoms are in the linked list

   27.06.95 Original   By: ACRM
*/
BOOL FindCAResidues(PDB *pdbca, char chain1, int resnum1, char insert1,
                    char chain2, int resnum2, char insert2,
                    PDB **pp_start, PDB **pp_end)
{
   PDB *p;

   *pp_start = NULL;
   *pp_end   = NULL;
   
   /* Search for the residue before the one specified by ID 1           */
   for(p=pdbca; p!=NULL; NEXT(p))
   {
      if((p->next != NULL)                &&
         (p->next->resnum    == resnum1)  &&
         (p->next->chain[0]  == chain1)   &&
         (p->next->insert[0] == insert1))
      {
         *pp_start = p;
         break;
      }
   }
   
   if(*pp_start==NULL)
   {
      fprintf(stderr,"Unable to find residue (before) %c%d%c\n",
              chain1, resnum1, insert1);
      return(FALSE);
   }
   
   /* Search for the residue before the one specified by ID 2           */
   for(p=pdbca; p!=NULL; NEXT(p))
   {
      if((p->next != NULL)                &&
         (p->next->resnum    == resnum2)  &&
         (p->next->chain[0]  == chain2)   &&
         (p->next->insert[0] == insert2))
      {
         *pp_end = p;
         break;
      }
   }
   
   if(*pp_end==NULL)
   {
      fprintf(stderr,"Unable to find residue (before) %c%d%c\n",
              chain2, resnum2, insert2);
      return(FALSE);
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL StoreTorsions(PDB *allatompdb, PDB *pdb, PDB *p_start, PDB *p_end,
                      char *filename, char *start, char *end)
   -----------------------------------------------------------------------
   Input:   PDB   *allatompdb   The start of the raw all-atom PDB linked 
                                list
            PDB   *pdb          The start of the CA/backbone PDB linked 
                                list
            PDB   *p_start      Atom before start of the loop in the 
                                CA/backbone PDB linked list
            PDB   *p_end        Atom vefore end of the loop in the 
                                CA/backbone PDB linked list
            char  *filename     PDB filename
            char  *start        Residue spceification for start of loop
            char  *end          Residue spceification for end of loop
   Returns: BOOL                Success of memory allocation
   Globals: DATALIST *gDataList Created linked list of loop 
                                specifications and associated torsion data
            
   Calculate and store torsions in a global linked list

   27.06.95 Original   By: ACRM
   06.07.95 Stores start of PDB linked list and first interesting
            residue in DATALIST
   26.07.95 Moved setting of loopid so it isn't done on every loop!
   08.08.95 New parameter, allatompdb, now stored in the global linked 
            list. Start and end strings stored explicitly in structure
   13.09.95 Also stores distances between CAs if gDoDistances is TRUE
   21.09.95 Also stores the angles at each CA if gDoAngles is TRUE
   09.01.96 Added !gDoCritRes handling
            Truncated structures were causing the whole data list to
            be freed rather than just this entry.
   30.01.09 Initialize some variables
*/
BOOL StoreTorsions(PDB *allatompdb, PDB *pdb, PDB *p_start, PDB *p_end, 
                   char *filename, char *start, char *end)
{
   PDB             *p1, *p2 = NULL, 
                   *p3, *p4;
   static DATALIST *p = NULL,
                   *prev = NULL;
   int             i;

   /* Allocate space in data linked list                                */
   if(gDataList == NULL)
   {
      INIT(gDataList, DATALIST);
      p=gDataList;
   }
   else
   {
      prev = p;
      ALLOCNEXT(p, DATALIST);
   }
   if(p==NULL)
   {
      FREELIST(gDataList, DATALIST);
      fprintf(stderr,"No memory for storing torsions\n");
      return(FALSE);
   }

   if(!gDoCritRes)
   {
      p->allatompdb  = NULL;
      p->torsionpdb  = NULL;
   }
   else
   {
      p->allatompdb  = allatompdb;
      p->torsionpdb  = pdb;
   }
   
   p->length      = 0;
   p->pdbloop     = p_start->next; /* Start of loop itself since p_start
                                      is the atom before the loop
                                   */
   sprintf(p->loopid,"%s-%s-%s",filename,start,end);
   strcpy(p->start, start);
   strcpy(p->end,   end);

   for(p1=p_start, i=0; 
       p1!=p_end->next && i<MAXLOOPLEN*TORPERRES; 
       NEXT(p1), i++)
   {
      p2 = p1->next;
      p3 = p2->next;
      p4 = p3->next;
   
      /* Check all are valid atoms                                      */
      if(p1==NULL || p2==NULL || p3==NULL || p4==NULL)
      {
         /* 09.01.95 Fixed destruction of whole list...                 
---      FREELIST(gDataList, DATALIST);
         */
         free(p);
         p = prev;
         
         fprintf(stderr,"Structure is truncated, unable to calculate all \
torsions.\n");
         return(FALSE);
      }
      
      /* Calculate torsion                                              */
      p->torsions[i] = blPhi(p1->x, p1->y, p1->z,
                             p2->x, p2->y, p2->z,
                             p3->x, p3->y, p3->z,
                             p4->x, p4->y, p4->z);

      /* Calculate angle                                                */
      p->angles[i] = blAngle(p1->x, p1->y, p1->z,
                             p2->x, p2->y, p2->z,
                             p3->x, p3->y, p3->z);

      /* Update loop length when p3 is a CA (i.e. for the phi angle)    */
      if(!strncmp(p3->atnam,"CA  ",4))
         (p->length)++;
   }

   /* See if we ran out of storage space                                */
   if(p1 != p_end->next)
   {
      fprintf(stderr,"Loop length exceeded maximum of %d\n",MAXLOOPLEN);
      return(FALSE);
   }
   
   /* Store distances if required                                       */
   if(gDoDistance)
   {
      BOOL FirstCA = TRUE;
      
      for(p1=p->pdbloop, i=0;
          p1!=p_end->next && i<MAXLOOPLEN;
          NEXT(p1))
      {
         if(!strncmp(p1->atnam,"CA  ",4))
         {
            if(FirstCA)
            {
               FirstCA = FALSE;
               p2 = p1;
            }
            p->dist[i] = DIST(p1,p2);
            i++;
         }
      }
   }

   /* See if we ran out of storage space                                */
   if(p1 != p_end->next)
   {
      fprintf(stderr,"Loop length exceeded maximum of %d\n",MAXLOOPLEN);
      return(FALSE);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL FindBBResidues(PDB *pdbbb, char chain1, int resnum1, char insert1,
                       char chain2, int resnum2, char insert2,
                       PDB **pp_start, PDB **pp_end)
   ---------------------------------------------------------------------
   Input:   PDB  *pdbbb      Backbone PDB linked list
            char chain1      Start of loop chain spec
            int  resnum1     Start of loop residue number
            char insert1     Start of loop insert code
            char chain2      End of loop chain spec
            int  resnum2     End of loop residue number
            char insert2     End of loop insert code
   Output:  PDB  **pp_start  Pointer to atom before start of loop in 
                             linked list
            PDB  **pp_end    Pointer to atom before end of loop in 
                             linked list
   Returns: BOOL             Success of finding specfied residues

   Finds pointers to the C before the specified start and the CA in
   the specified end residues.

   N.B. This assumes N,CA,C ordering within the PDB file.

   27.06.95 Original   By: ACRM
*/
BOOL FindBBResidues(PDB *pdbbb, char chain1, int resnum1, char insert1,
                    char chain2, int resnum2, char insert2,
                    PDB **pp_start, PDB **pp_end)
{
   PDB *p;

   *pp_start = NULL;
   *pp_end   = NULL;
   
   /* Search for the C atom before the residue specified by ID 1        */
   for(p=pdbbb; p!=NULL; NEXT(p))
   {
      if((!strncmp(p->atnam,"C   ",4))    &&
         (p->next != NULL)                &&
         (p->next->resnum    == resnum1)  &&
         (p->next->chain[0]  == chain1)   &&
         (p->next->insert[0] == insert1))
      {
         *pp_start = p;
         break;
      }
   }
   
   if(*pp_start==NULL)
   {
      fprintf(stderr,"Unable to find C in residue before %c%d%c\n",
              chain1, resnum1, insert1);
      return(FALSE);
   }
   
   /* Search for the CA in the residue specified by ID 2                */
   for(p=pdbbb; p!=NULL; NEXT(p))
   {
      if((!strncmp(p->atnam,"CA  ",4)) &&
         (p->resnum    == resnum2)     &&
         (p->chain[0]  == chain2)      &&
         (p->insert[0] == insert2))
      {
         *pp_end = p;
         break;
      }
   }
   
   if(*pp_end==NULL)
   {
      fprintf(stderr,"Unable to find CA in residue after %c%d%c\n",
              chain2, resnum2, insert2);
      return(FALSE);
   }

   return(TRUE);
}


/************************************************************************/
/*>REAL **ConvertData(DATALIST *indata, int *NData, BOOL CATorsions)
   -----------------------------------------------------------------
   Input:   DATALIST *indata      Linked list of loop data structures
            BOOL     CATorsions   Do CA rather than true torsions
   Output:  int      *NData       Number of items in linked list
            
   Converts the linked list of torsions into a 2D array of sin and cos
   values of the torsions.

   27.06.95 Original   By: ACRM
   04.07.94 No longer frees up the input data linked list
   10.08.95 Changed explicit 9999.0 to DUMMY
   13.09.95 Added storage of distances
   21.09.95 Modified method to calc storage requirements & handling of
            angle data
*/
REAL **ConvertData(DATALIST *indata, int *NData, BOOL CATorsions)
{
   REAL     **data;
   DATALIST *p;
   int      i,
            n,
            count,
            pos,
            maxval,
            ArrayDim,
            AngleOffset = 2,
            DistOffset  = 2;

   /* The number of array positions used per residue and offsets to
      different parts of the data
   */
   maxval = 2;
   if(!CATorsions) 
   {
      maxval      += 4;
      AngleOffset += 4;
      DistOffset  += 4;
   }
   if(gDoAngles)
   {
      maxval      += 1;
      DistOffset  += 1;
   }
   if(gDoDistance) 
   {
      maxval      += 1;
   }
   
   
   /* Count items in data linked list                                   */
   for(p=indata, n=0; p!=NULL; NEXT(p), n++);
   *NData = n;
   
   /* Calculate array dimension                                         */
   ArrayDim = gMaxLoopLen * maxval;

   /* Allocate and check 2D array                                       */
   if((data = (REAL **)blArray2D(sizeof(REAL),*NData,ArrayDim))==NULL)
   {
      fprintf(stderr,"No memory for data array.\n");
      return(NULL);
   }

   /* Walk the linked list, copying in the data                         */
   for(p=indata,n=0; p!=NULL; NEXT(p),n++)
   {
      /* Set everything to DUMMY                                        */
      for(count=0; count<gMaxLoopLen; count++)
      {
         for(i=0; i<maxval; i++)
            data[n][count*maxval + i] = DUMMY;
         if(gDoDistance)
            data[n][count*maxval + DistOffset] = DUMMY2;
      }

      /* Insert from the start of the scheme till we hit something bigger
         than the number of residues we have.
      */
      for(count=0; 
          (gScheme[count] <= p->length) && (count<gMaxLoopLen); 
          count++)
      {
         if(CATorsions)
         {
            data[n][count*maxval]     = sin(p->torsions[count]);
            data[n][count*maxval + 1] = cos(p->torsions[count]);
         }
         else
         {
            data[n][count*maxval]     = sin(p->torsions[count*3]);
            data[n][count*maxval + 1] = cos(p->torsions[count*3]);
            data[n][count*maxval + 2] = sin(p->torsions[count*3 + 1]);
            data[n][count*maxval + 3] = cos(p->torsions[count*3 + 1]);
            data[n][count*maxval + 4] = sin(p->torsions[count*3 + 2]);
            data[n][count*maxval + 5] = cos(p->torsions[count*3 + 2]);
         }
         
         if(gDoDistance)
            data[n][count*maxval + DistOffset] = p->dist[count];

         if(gDoAngles)
         {
            data[n][count*maxval + AngleOffset]   = 2.0*(p->angles[count])/PI - 1.0;
/*            data[n][count*maxval + AngleOffset+1] = cos(p->angles[count]); */
         }
      }
         
      /* Insert from the end of the scheme till we hit something bigger
         than the number of residues we have.
      */
      for(count=gMaxLoopLen-1, pos=p->length-1; 
          (gScheme[count] <= p->length) && (count>=0); 
          count--, pos--)
      {
         if(CATorsions)
         {
            data[n][count*maxval]     = sin(p->torsions[pos]);
            data[n][count*maxval + 1] = cos(p->torsions[pos]);
         }
         else
         {
            data[n][count*maxval]     = sin(p->torsions[pos*3]);
            data[n][count*maxval + 1] = cos(p->torsions[pos*3]);
            data[n][count*maxval + 2] = sin(p->torsions[pos*3 + 1]);
            data[n][count*maxval + 3] = cos(p->torsions[pos*3 + 1]);
            data[n][count*maxval + 4] = sin(p->torsions[pos*3 + 2]);
            data[n][count*maxval + 5] = cos(p->torsions[pos*3 + 2]);
         }
         
         if(gDoDistance)
            data[n][count*maxval + DistOffset] = p->dist[pos];

         if(gDoAngles)
         {
            data[n][count*maxval + AngleOffset]   = 2.0*(p->angles[pos])/PI - 1.0;
/*            data[n][count*maxval + AngleOffset+1] = cos(p->angles[pos]); */
         }
      }
   }

#ifdef DEBUG
   PrintArray(data, *NData, ArrayDim);
#endif

   return(data);
}


/************************************************************************/
/*>void PrintArray(REAL **data, int NData, int width)
   --------------------------------------------------
   Input:   REAL  **data         Data array for clustering
            int   NData          Number of vectors
            int   width          Dimension of vectors

   Prints the 2D data array for debugging

   27.06.95 Original   By: ACRM
*/
void PrintArray(REAL **data, int NData, int width)
{
   int i,j;
   
   for(i=0; i<NData; i++)
   {
      for(j=0; j<width; j++)
      {
         fprintf(stderr,"%5.2f ",data[i][j]);
      }
      fprintf(stderr,"\n");
   }
}


