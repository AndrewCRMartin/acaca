/*************************************************************************

   Program:    
   File:       decr.c
   
   Version:    V3.7
   Date:       06.02.96
   Function:   DEfine Critical Residues
   
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
   V1.0  01.08.95 Original
   V3.2  02.10.95 Modified for completely conserved residues
   V3.3  04.10.95 Modified to show residues which are conserved in at
                  least one conformation
   V3.4  10.10.95 Various changes to make deleted residues work with 
                  residues conserved in at least one cluster
   V3.5  06.11.95 Skipped
   V3.6  09.01.96 Skipped
   V3.7  06.02.96 Separated out bits for findsdrs
   V3.7a 30.01.09 Fixed initial check on same residue

*************************************************************************/
/* Includes
*/
#include "decr.h"
#include "decr2.h"


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
/*>BOOL FindNeighbourProps(PDB *pdb, PDB *start, PDB *stop, int clusnum,
                           LOOPINFO *loopinfo)
   ---------------------------------------------------------------------
   Input:   PDB      *pdb        PDB linked list for whole structure
            PDB      *start      Start of loop in pdb
            PDB      *stop       Record after end of loop in pdb
            int      clusnum     Cluster number
   Output:  LOOPINFO *loopinfo   Details of loop and contacting residues
                                 (Space is allocated within this 
                                 structure)
   Returns: BOOL                 Success of memory allocations

   Peforms all allocations in a LOOPINFO structure and fills it in with
   details of the loop and contacting residues.

   01.08.95 Original    By: ACRM
   02.08.95 Added AALoop and AAContact, clusnum parameter
   10.10.95 Changed USHORT to PROP_T
*/
BOOL FindNeighbourProps(PDB *pdb, PDB *start, PDB *stop, int clusnum,
                        LOOPINFO *loopinfo)
{
   PDB  *p, *p_next,
        *q, *q_next,
        **contacts;
   int  ncontacts   = 0,
        maxcontacts = ALLOCQUANTUM,
        looplen     = 0;
   BOOL InArray;

   /* Allocate an array to store the contact PDB pointers               */
   if((contacts=(PDB **)malloc(maxcontacts*sizeof(PDB)))==NULL)
      return(FALSE);
   
   /* Step through the loop                                             */
   for(p=start,looplen=0; p!=stop; p=p_next)
   {
      /* Find the following residue                                     */
      p_next = blFindNextResidue(p);
      looplen++;

      /* For each residue N-ter to the loop                             */
      for(q=pdb; q!=start && q!=NULL; q=q_next)
      {
         q_next = blFindNextResidue(q);

         /* If loop residue makes contact with this non-loop residue    */
         if(ResidueContact(p,p_next,q,q_next,CONTACTDIST))
         {
            /* If this non-loop residue is not already stored           */
            TESTINARRAY(contacts,ncontacts,q,InArray);
            if(!InArray)
            {
               /* Store the contact residue pointer                     */
               contacts[ncontacts++] = q;

               /* Increase the contacts array size                      */
               if(ncontacts==maxcontacts)
               {
                  maxcontacts += ALLOCQUANTUM;
                  if((contacts=(PDB **)
                      realloc(contacts, maxcontacts*sizeof(PDB)))==NULL)
                     return(FALSE);
               }
            }
         }
      }
                  
      /* For each residue C-ter to the loop                             */
      for(q=stop; q!=NULL; q=q_next)
      {
         q_next = blFindNextResidue(q);

         /* If loop residue makes contact with this non-loop residue    */
         if(ResidueContact(p,p_next,q,q_next,CONTACTDIST))
         {
            /* If this non-loop residue is not already stored           */
            TESTINARRAY(contacts,ncontacts,q,InArray);
            if(!InArray)
            {
               /* Store the contact residue pointer                     */
               contacts[ncontacts++] = q;

               /* Increase the contacts array size                      */
               if(ncontacts==maxcontacts)
               {
                  maxcontacts += ALLOCQUANTUM;
                  if((contacts=(PDB **)
                      realloc(contacts, maxcontacts*sizeof(PDB)))==NULL)
                     return(FALSE);
               }
            }
         }
      }
   }

   /* Set values in the loopinfo structure                              */
   loopinfo->ncontacts = ncontacts;
   loopinfo->contacts  = contacts;
   loopinfo->length    = looplen;

   /* Allocate an array to store the loop residue pointers              */
   if((loopinfo->residues=(PDB **)malloc(looplen*sizeof(PDB)))==NULL)
   {
      free(contacts);
      loopinfo->contacts = NULL;
      return(FALSE);
   }
   
   /* Copy residue pointers into this array                             */
   for(p=start,looplen=0; p!=stop; p=p_next)
   {
      p_next = blFindNextResidue(p);
      loopinfo->residues[looplen++] = p;
   }

   /* Allocate memory for property arrays                              */
   loopinfo->ResProps = 
      (PROP_T *)malloc(looplen*sizeof(PROP_T));
   loopinfo->ContactProps = 
      (PROP_T *)malloc(ncontacts*sizeof(PROP_T));
   loopinfo->AALoop = 
      (char *)malloc(looplen*sizeof(char));
   loopinfo->AAContact = 
      (char *)malloc(ncontacts*sizeof(char));
   loopinfo->ResFlag = 
      (BOOL *)malloc(looplen*sizeof(BOOL));
   loopinfo->ContactFlag = 
      (BOOL *)malloc(ncontacts*sizeof(BOOL));
   

   if(loopinfo->ResProps     == NULL || 
      loopinfo->ContactProps == NULL ||
      loopinfo->AALoop       == NULL ||
      loopinfo->AAContact    == NULL ||
      loopinfo->ResFlag      == NULL ||
      loopinfo->ContactFlag  == NULL)
   {
      if(loopinfo->ResProps     != NULL) free(loopinfo->ResProps);
      if(loopinfo->ContactProps != NULL) free(loopinfo->ContactProps);
      if(loopinfo->AALoop       != NULL) free(loopinfo->AALoop);
      if(loopinfo->AAContact    != NULL) free(loopinfo->AAContact);
      if(loopinfo->contacts     != NULL) free(loopinfo->contacts);
      if(loopinfo->residues     != NULL) free(loopinfo->residues);
      if(loopinfo->ResFlag      != NULL) free(loopinfo->ResFlag);
      if(loopinfo->ContactFlag  != NULL) free(loopinfo->ContactFlag);

      loopinfo->ResProps     = NULL;
      loopinfo->ContactProps = NULL;
      loopinfo->AALoop       = NULL;
      loopinfo->AAContact    = NULL;
      loopinfo->contacts     = NULL;
      loopinfo->residues     = NULL;
      loopinfo->ResFlag      = NULL;
      loopinfo->ContactFlag  = NULL;
      
      return(FALSE);
   }

   FillLoopInfo(loopinfo);
   loopinfo->clusnum = clusnum;
    
   return(TRUE);
}


/************************************************************************/
/*>BOOL ResidueContact(PDB *p_start, PDB *p_stop, PDB *q_start, 
                       PDB *q_stop, REAL dist)
   ------------------------------------------------------------
   Input:   PDB  *p_start    Start of first residue
            PDB  *p_stop     Record after end of first residue
            PDB  *q_start    Start of second residue
            PDB  *q_stop     Record after end of second residue
            REAL dist        Maximum distace to define contact
   Returns: BOOL             In contact?

   See if a contact of <= dist Angstroms is made between atoms in the 
   residue bounded by pointers p_start/p_stop and sidechain atoms
   bounded by q_start/q_stop

   01.08.95 Original    By: ACRM
   30.01.09 Corrected initial self-contact check
*/
BOOL ResidueContact(PDB *p_start, PDB *p_stop, PDB *q_start, PDB *q_stop,
                    REAL dist)
{
   PDB *p, *q;

   /* Ignore contact with itself                                        */
   if(p_start==q_start)  /* 30.01.09 Corrected from p and q             */
      return(FALSE);

   /* Square the distance to save on doing square roots                 */
   dist *= dist;
   
   for(p=p_start; p!=p_stop; NEXT(p))
   {
      for(q=q_start; q!=q_stop; NEXT(q))
      {
         if(strncmp(q->atnam,"N   ",4) &&
            strncmp(q->atnam,"CA  ",4) &&
            strncmp(q->atnam,"C   ",4) &&
            strncmp(q->atnam,"O   ",4))
         {
            if(DISTSQ(p,q) <= dist)
               return(TRUE);
         }
      }
   }

   return(FALSE);
}


/************************************************************************/
/*>void FillLoopInfo(LOOPINFO *loopinfo)
   -------------------------------------
   I/O:     LOOPINFO  *loopinfo     Input with PDB pointer arrays filled
                                    in listing the loop and contact 
                                    residues.
                                    Output with the residue property and
                                    contact property arrays completed.

   Fill in residue property flags in a loopinfo structure for both the
   loop and contacting residues.

   01.08.95 Original    By: ACRM
   02.08.95 Added AALoop and AAContact
*/
void FillLoopInfo(LOOPINFO *loopinfo)
{
   int  i;
   char res;
   
   
   for(i=0; i<loopinfo->length; i++)
   {
      res = blThrone((loopinfo->residues[i])->resnam);
      loopinfo->AALoop[i] = res;
      loopinfo->ResProps[i] = SetProperties(res);
   }
   for(i=0; i<loopinfo->ncontacts; i++)
   {
      res = blThrone((loopinfo->contacts[i])->resnam);
      loopinfo->AAContact[i] = res;
      loopinfo->ContactProps[i] = SetProperties(res);
   }
}


/************************************************************************/
/*>BOOL MergeProperties(int NLoops, LOOPINFO *loopinfo, int clusnum, 
                        CLUSTERINFO *clusterinfo)
   -----------------------------------------------------------------
   Input:   int         NLoops       Number of loops in a cluster
            LOOPINFO    *loopinfo    Array of completed structures for
                                     loops in this cluster
            int         clusnum      The number of the cluster in which
                                     we are interested
   Output:  CLUSTERINFO *clusterinfo Compiled data about this cluster
                                     (Memory allocated within this
                                     structure)
   Returns: BOOL                     Success of memory allocations

   Allocate memory in and complete a clusterinfo structure with merged
   property data for the residues ids common to all loops.

   03.08.95 Original    By: ACRM
   08.08.95 Added setting of clusterinfo->length from first loop's
            length
   02.09.95 Added ->ConsRes[] and ->absolute[] handling
            Corrected check on NULLs for each free() to != rather 
            than ==
   05.10.95 Added clusterinfo->First = NULL;
   10.10.95 Added clusterinfo->deletable = NULL;
            Changed USHORT to PROP_T
*/
BOOL MergeProperties(int NLoops, LOOPINFO *loopinfo, int clusnum,
                     CLUSTERINFO *clusterinfo)
{
   int  i,j,k,
        NRes,
        first;
   BOOL GotResid;

   /* Find the residue IDs common to all loops in this cluster          */
   if((NRes = FlagCommonResidues(NLoops, loopinfo, clusnum)) < 0)
      return(FALSE);
   clusterinfo->NRes = NRes;
   clusterinfo->length = loopinfo[0].length;

   /* If there weren't any common residues (or FlagCommonResidues() found
      no loops in this cluster, just return
   */
   if(NRes == 0)
      return(TRUE);

   /* Allocate memory in the clusterinfo structure to store properties
      for these residues
   */
   clusterinfo->resnum         = (int *)malloc(NRes * sizeof(int));
   clusterinfo->chain          = (char *)malloc(NRes * sizeof(char));
   clusterinfo->insert         = (char *)malloc(NRes * sizeof(char));
   clusterinfo->ConservedProps = (PROP_T *)malloc(NRes * sizeof(PROP_T));
   clusterinfo->RangeOfProps   = (PROP_T *)malloc(NRes * sizeof(PROP_T));
   clusterinfo->absolute       = (BOOL *)malloc(NRes * sizeof(BOOL));
   clusterinfo->ConsRes        = (char *)malloc(NRes * sizeof(char));
   clusterinfo->First          = NULL;
   clusterinfo->deletable      = NULL;

   if((clusterinfo->resnum         == NULL) ||
      (clusterinfo->chain          == NULL) ||
      (clusterinfo->insert         == NULL) ||
      (clusterinfo->ConservedProps == NULL) ||
      (clusterinfo->RangeOfProps   == NULL) ||
      (clusterinfo->absolute       == NULL) ||
      (clusterinfo->ConsRes        == NULL))
   {
      if(clusterinfo->resnum != NULL)
         free(clusterinfo->resnum);
      if(clusterinfo->chain  != NULL)
         free(clusterinfo->chain);
      if(clusterinfo->insert != NULL)
         free(clusterinfo->insert);
      if(clusterinfo->ConservedProps != NULL)
         free(clusterinfo->ConservedProps);
      if(clusterinfo->RangeOfProps != NULL)
         free(clusterinfo->RangeOfProps);
      if(clusterinfo->ConsRes != NULL)
         free(clusterinfo->ConsRes);
      if(clusterinfo->absolute != NULL)
         free(clusterinfo->absolute);
      return(FALSE);
   }
   
   for(i=0; i<NRes; i++)
   {
      clusterinfo->ConsRes[i]  = '-';
      clusterinfo->absolute[i] = TRUE;
   }
   
   /* Find the first loop in the specified cluster                      */
   first = (-1);
   for(i=0; i<NLoops; i++)
   {
      if(loopinfo[i].clusnum == clusnum)
      {
         first = i;
         break;
      }
   }
   
   /* This shouldn't happen as it's checked for in FlagCommonResidues() */
   if(first == (-1))
   {
      fprintf(stderr,"INTERR: FlagCommonResidues() found cluster %d, but \
MergeProperties() can't\n",clusnum);
      return(FALSE);
   }

   /* Copy in the flagged residue ids from the first loop               */
   for(j=0,k=0; j<loopinfo[first].length; j++)   /* Loop residues       */
   {
      if(loopinfo[first].ResFlag[j])
      {
         clusterinfo->chain[k]  = loopinfo[first].residues[j]->chain[0];
         clusterinfo->resnum[k] = loopinfo[first].residues[j]->resnum;
         clusterinfo->insert[k] = loopinfo[first].residues[j]->insert[0];
         clusterinfo->ConservedProps[k] = loopinfo[first].ResProps[j];
         clusterinfo->RangeOfProps[k]   = loopinfo[first].ResProps[j];
         clusterinfo->ConsRes[k]        = loopinfo[first].AALoop[j];
         
         k++;
      }
   }
   for(j=0; j<loopinfo[first].ncontacts; j++)    /* Contacting residues */
   {
      if(loopinfo[first].ContactFlag[j])
      {
         clusterinfo->chain[k]  = loopinfo[first].contacts[j]->chain[0];
         clusterinfo->resnum[k] = loopinfo[first].contacts[j]->resnum;
         clusterinfo->insert[k] = loopinfo[first].contacts[j]->insert[0];
         clusterinfo->ConservedProps[k] = loopinfo[first].ContactProps[j];
         clusterinfo->RangeOfProps[k]   = loopinfo[first].ContactProps[j];
         clusterinfo->ConsRes[k]        = loopinfo[first].AAContact[j];

         k++;
      }
   }

   /* Examine the loopinfo array for each loop                          */
   for(i=0; i<NLoops; i++)
   {
      /* Is it the cluster of interest?                                 */
      if(loopinfo[i].clusnum == clusnum)
      {
         /* Work our way through the loop properties array              */
         for(j=0; j<loopinfo[i].length; j++)
         {
            /* If it's a flagged (common) residue                       */
            if(loopinfo[i].ResFlag[j])
            {
               /* Find it in the residue id arrays                      */
               GotResid = FALSE;
               for(k=0; k<clusterinfo->NRes; k++)
               {
                  if((clusterinfo->chain[k] == 
                      loopinfo[i].residues[j]->chain[0]) &&
                     (clusterinfo->resnum[k] ==
                      loopinfo[i].residues[j]->resnum)   &&
                     (clusterinfo->insert[k] ==
                      loopinfo[i].residues[j]->insert[0]))
                  {
                     GotResid = TRUE;
                     break;
                  }
               }               
               
               if(GotResid)
               {
                  /* If this residue id is already stored, combine 
                     properties from this loop into the array  
                  */
                  clusterinfo->ConservedProps[k] &= 
                     loopinfo[i].ResProps[j];
                  clusterinfo->RangeOfProps[k]   |= 
                     loopinfo[i].ResProps[j];
                  if(clusterinfo->ConsRes[k] != loopinfo[i].AALoop[j])
                     clusterinfo->absolute[k] = FALSE;
               }
               else
               {
                  /* This shouldn't happen!                             */
                  fprintf(stderr,"Cockup! Residue in loop %d flagged as \
common but not in loop 0\n",i);
                  blWritePDBRecord(stderr,loopinfo[i].residues[j]);
               }
            }
         }

         /* Repeat for the contacting residue array                     */
         for(j=0; j<loopinfo[i].ncontacts; j++)
         {
            /* If it's a flagged (common) residue                       */
            if(loopinfo[i].ContactFlag[j])
            {
               /* Find it in the residue id arrays                      */
               GotResid = FALSE;
               for(k=0; k<clusterinfo->NRes; k++)
               {
                  if((clusterinfo->chain[k] == 
                      loopinfo[i].contacts[j]->chain[0]) &&
                     (clusterinfo->resnum[k] ==
                      loopinfo[i].contacts[j]->resnum)   &&
                     (clusterinfo->insert[k] ==
                      loopinfo[i].contacts[j]->insert[0]))
                  {
                     GotResid = TRUE;
                     break;
                  }
               }               

               if(GotResid)
               {
                  /* If this residue id is already stored, combine 
                     properties from this loop into the array  
                  */
                  clusterinfo->ConservedProps[k] &= 
                     loopinfo[i].ContactProps[j];
                  clusterinfo->RangeOfProps[k]   |= 
                     loopinfo[i].ContactProps[j];
                  if(clusterinfo->ConsRes[k] != loopinfo[i].AAContact[j])
                     clusterinfo->absolute[k] = FALSE;
               }
               else
               {
                  /* This shouldn't happen!                             */
                  fprintf(stderr,"Cockup! Residue in loop %d flagged as \
common but not in loop 0\n",i);
                  blWritePDBRecord(stderr,loopinfo[i].contacts[j]);
               }
            }
         }
      }
   }
   return(TRUE);
}


/************************************************************************/
/*>void BlankClusterInfo(CLUSTERINFO *clusterinfo)
   -----------------------------------------------
   Output:  CLUSTERINFO  *clusterinfo    Cleared structure.

   Set all pointers in a clusterinfo structure to NULL

   02.08.95 Original    By: ACRM
   08.08.95 Added length
   02.10.95 Added absolute and ConsRes
   04.10.95 Doesn't zero ->length, as this isn't set by the second
            critical residue phase.
*/
void BlankClusterInfo(CLUSTERINFO *clusterinfo)
{
   clusterinfo->resnum         = NULL;
   clusterinfo->chain          = NULL;
   clusterinfo->insert         = NULL;
   clusterinfo->ConservedProps = NULL;
   clusterinfo->RangeOfProps   = NULL;
   clusterinfo->absolute       = NULL;
   clusterinfo->ConsRes        = NULL;
   clusterinfo->First          = NULL;
   clusterinfo->NRes           = 0;
}


/************************************************************************/
/*>void BlankLoopInfo(LOOPINFO *loopinfo)
   --------------------------------------
   I/O:     LOOPINFO  *loopinfo     A loopinfo structure pointer

   Clear all info in a loop info structure (assumes memory has been
   freed or not yet allocated).

   08.08.95 Original    By: ACRM
*/
void BlankLoopInfo(LOOPINFO *loopinfo)
{
   loopinfo->contacts     = NULL;
   loopinfo->residues     = NULL;
   loopinfo->ContactProps = NULL;
   loopinfo->ResProps     = NULL;
   loopinfo->AALoop       = NULL;
   loopinfo->AAContact    = NULL;
   loopinfo->ResFlag      = NULL;
   loopinfo->ContactFlag  = NULL;
   loopinfo->length       = 0;
   loopinfo->clusnum      = 0;
   loopinfo->ncontacts    = 0;
}


/************************************************************************/
/*>int FlagCommonResidues(int NLoops, LOOPINFO *loopinfo, int clusnum)
   -------------------------------------------------------------------
   Input:   int      NLoops     Number of loops
            int      clusnum    Cluster number of interest
   I/O:     LOOPINFO *loopinfo  Array of structures containing loop info
                                On output, the ResFlag and ContactFlag
                                arrays are filled in with the residues
                                common to all loops in the array
   Returns: int                 Number of residues in common
                                (-1 if a memory allocation failed)

   Runs through the NLoops loops in the loopinfo structure array, looking
   at cluster clusnum. Finds residue ids common to all loops in this
   cluster and flags them in the loopinfo structures.

   Returns the total number of common residues.

   03.08.95 Original   By: ACRM
*/
int FlagCommonResidues(int NLoops, LOOPINFO *loopinfo, int clusnum)
{
   int  i, j, k,
        looplen,
        *resnum,
        *count,
        nres,
        retval = (-1),
        first;
   char *chain,
        *insert;

   /* Find how many residues there are in the first example and allocate
      memory of this size
   */
   nres = loopinfo[0].length + loopinfo[0].ncontacts;

   chain  = (char *)malloc(nres * sizeof(char));
   insert = (char *)malloc(nres * sizeof(char));
   resnum = (int *)malloc(nres * sizeof(int));
   count  = (int *)malloc(nres * sizeof(int));

   /* If allocations all OK, continue                                   */
   if((chain  != NULL) && (insert != NULL) && 
      (resnum != NULL) && (count  != NULL))
   {
      /* Initialise retval to 0 to say that allocations succeeded       */
      retval = 0;
      
      /* Find the first loop which is in the required cluster           */
      first = (-1);
      for(i=0; i<NLoops; i++)
      {
         if(loopinfo[i].clusnum == clusnum)
         {
            first = i;
            break;
         }
      }
      
      /* If a loop is found in the required cluster                     */
      if(first != (-1))
      {
         /* Copy the residue ids into the arrays for first loop         */
         looplen = loopinfo[first].length;
         for(j=0; j<looplen; j++)                /* Loop residues       */
         {
            chain[j]  = loopinfo[first].residues[j]->chain[0];
            resnum[j] = loopinfo[first].residues[j]->resnum;
            insert[j] = loopinfo[first].residues[j]->insert[0];
            count[j]  = 1;
         }
         for(j=0; j<loopinfo[0].ncontacts; j++)  /* Contacting residues */
         {
            chain[looplen+j]  = loopinfo[first].contacts[j]->chain[0];
            resnum[looplen+j] = loopinfo[first].contacts[j]->resnum;
            insert[looplen+j] = loopinfo[first].contacts[j]->insert[0];
            count[looplen+j]  = 1;
         }
         
         
         /* Run through the rest of the loops incrementing the count if
            the residue label is found in this loop
         */
         for(i=first+1; i<NLoops; i++)
         {
            if(loopinfo[i].clusnum == clusnum)
            {
               for(j=0; j<loopinfo[i].length; j++)    /* Loop residues  */
               {
                  for(k=0; k<nres; k++)
                  {
                     if((chain[k] ==loopinfo[i].residues[j]->chain[0]) &&
                        (resnum[k]==loopinfo[i].residues[j]->resnum)   &&
                        (insert[k]==loopinfo[i].residues[j]->insert[0]))
                        (count[k])++;
                  }
               }
               
               for(j=0; j<loopinfo[i].ncontacts; j++) /* Contacting res */
               {
                  for(k=0; k<nres; k++)
                  {
                     if((chain[k] ==loopinfo[i].contacts[j]->chain[0]) &&
                        (resnum[k]==loopinfo[i].contacts[j]->resnum)   &&
                        (insert[k]==loopinfo[i].contacts[j]->insert[0]))
                        (count[k])++;
                  }
               }
            }
         }
         
         /* Run through the loops again flagging all residues which are 
            common to all loops.
         */
         for(i=first; i<NLoops; i++)
         {
            if(loopinfo[i].clusnum == clusnum)
            {
               for(j=0; j<loopinfo[i].length; j++)    /* Loop residues  */
               {
                  loopinfo[i].ResFlag[j] = FALSE;
                  for(k=0; k<nres; k++)
                  {
                     if((chain[k] ==loopinfo[i].residues[j]->chain[0])  &&
                        (resnum[k]==loopinfo[i].residues[j]->resnum)    &&
                        (insert[k]==loopinfo[i].residues[j]->insert[0]) &&
                        (count[k] ==NLoops))
                     {
                        loopinfo[i].ResFlag[j] = TRUE;
                        break;
                     }
                  }
               }
               
               for(j=0; j<loopinfo[i].ncontacts; j++) /* Contacting res */
               {
                  loopinfo[i].ContactFlag[j] = FALSE;
                  for(k=0; k<nres; k++)
                  {
                     if((chain[k] ==loopinfo[i].contacts[j]->chain[0])  &&
                        (resnum[k]==loopinfo[i].contacts[j]->resnum)    &&
                        (insert[k]==loopinfo[i].contacts[j]->insert[0]) &&
                        (count[k] ==NLoops))
                     {
                        loopinfo[i].ContactFlag[j] = TRUE;
                        break;
                     }
                  }
               }
            }
         }
         
         /* Count the common residues as the return value               */
         retval = 0;
         for(k=0; k<nres; k++)
         {
            if(count[k] == NLoops)
               retval++;
         }
      }  /* There was a loop in the specified cluster                   */
   }  /* Memory allocations OK                                          */
   
   /* Free up allocated memory                                          */
   if(chain  != NULL) free(chain);
   if(insert != NULL) free(insert);
   if(resnum != NULL) free(resnum);
   if(count  != NULL) free(count);
   
   return(retval);
}


/************************************************************************/
/*>void CleanLoopInfo(LOOPINFO *loopinfo, int NMembers)
   ----------------------------------------------------
   I/O:     LOOPINFO  *loopinfo     Array of loopinfo structures
   Input:   int       NMembers      Number of items in array

   Free up memory allocated within the loopinfo[] array

   08.08.95 Original    By: ACRM
*/
void CleanLoopInfo(LOOPINFO *loopinfo, int NMembers)
{
   int i;
   
   for(i=0; i<NMembers; i++)
   {
      if(loopinfo[i].contacts     != NULL) free(loopinfo[i].contacts);
      if(loopinfo[i].residues     != NULL) free(loopinfo[i].residues);
      if(loopinfo[i].ContactProps != NULL) free(loopinfo[i].ContactProps);
      if(loopinfo[i].ResProps     != NULL) free(loopinfo[i].ResProps);
      if(loopinfo[i].AALoop       != NULL) free(loopinfo[i].AALoop);
      if(loopinfo[i].AAContact    != NULL) free(loopinfo[i].AAContact);
      if(loopinfo[i].ResFlag      != NULL) free(loopinfo[i].ResFlag);
      if(loopinfo[i].ContactFlag  != NULL) free(loopinfo[i].ContactFlag);
      
      loopinfo[i].contacts     = NULL;
      loopinfo[i].residues     = NULL;
      loopinfo[i].ContactProps = NULL;
      loopinfo[i].ResProps     = NULL;
      loopinfo[i].AALoop       = NULL;
      loopinfo[i].AAContact    = NULL;
      loopinfo[i].ResFlag      = NULL;
      loopinfo[i].ContactFlag  = NULL;
   }
}


/************************************************************************/
/*>void CleanClusInfo(CLUSTERINFO *cinfo)
   --------------------------------------
   I/O:     CLUSTERINFO  *cinfo     Pointer to cluster infor structure
                                    to be cleaned

   Frees up memory alliocated in a cluster information structure and
   calls BlankClusterInfo() to clear everything.

   08.08.95 Original    By: ACRM
   02.10.95 Added absolute and ConsRes
   05.10.95 Added First
*/
void CleanClusInfo(CLUSTERINFO *cinfo)
{
   if(cinfo->resnum         != NULL) free(cinfo->resnum);
   if(cinfo->chain          != NULL) free(cinfo->chain);
   if(cinfo->insert         != NULL) free(cinfo->insert);
   if(cinfo->ConservedProps != NULL) free(cinfo->ConservedProps);
   if(cinfo->RangeOfProps   != NULL) free(cinfo->RangeOfProps);
   if(cinfo->absolute       != NULL) free(cinfo->absolute);
   if(cinfo->ConsRes        != NULL) free(cinfo->ConsRes);
   if(cinfo->First          != NULL) free(cinfo->First);

   BlankClusterInfo(cinfo);
}


/************************************************************************/
/*>void PrintMergedProperties(FILE *fp, int clusnum, CLUSTERINFO cinfo,
                              int NMembers)
   --------------------------------------------------------------------
   Input:   FILE        *fp      Output file pointer
            int         clusnum  Cluster number
            CLUSTERINFO cinfo    Cluster information structure
            int         NMembers Number of members of the cluster

   Print property information from merged properties for a cluster

   08.08.95 Original    By: ACRM
   14.08.95 Always prints number of members
   02.10.95 Handles printing of absolutely conserved residues
   10.10.95 Modified printing for deletable as separate flag
            Changed USHORT to PROP_T
*/
void PrintMergedProperties(FILE *fp, int clusnum, CLUSTERINFO cinfo,
                           int NMembers)
{
   int    i;
   PROP_T props;
   
   fprintf(fp,"CLUSTER %d (Length = %d, Members = %d)\n",
           clusnum,cinfo.length,NMembers);

   if(NMembers < 2)
      fprintf(fp, "WARNING: This cluster has%s%s %s!\n",
              (NMembers?" only ":" "),
              (NMembers?"one":"no"),
              (NMembers?"member":"members"));

   if(cinfo.NRes == 0)
   {
      fprintf(fp, "WARNING: No common residues identified for this \
cluster!\n");
   }
   else
   {
      for(i=0;i<cinfo.NRes;i++)
      {
         props = cinfo.ConservedProps[i];
         if(cinfo.deletable!=NULL && cinfo.deletable[i])
            props |= DELETED_FLAG;
         
         if(!(cinfo.absolute[i] && cinfo.ConsRes[i] == '-'))
         {
            fprintf(fp,"%c%3d%c 0x%04x ",
                    cinfo.chain[i],
                    cinfo.resnum[i],
                    cinfo.insert[i],
                    props);
            PrintProps(fp,cinfo.ConservedProps[i],
                       ((cinfo.deletable==NULL)?FALSE:
                        cinfo.deletable[i]));
            if(cinfo.absolute[i])
            {
               if(cinfo.deletable!=NULL && cinfo.deletable[i])
                  fprintf(fp," [CONSERVED/deletable] (%c-)",
                          cinfo.ConsRes[i]);
               else
                  fprintf(fp," [CONSERVED] (%c)",cinfo.ConsRes[i]);
            }
            else
            {
               PrintSampleResidues(fp,cinfo.ConservedProps[i],
                                   ((cinfo.deletable==NULL)?FALSE:
                                    cinfo.deletable[i]));
            }
            fprintf(fp,"\n");
         }
      }
   }
}


/************************************************************************/
/*>RESSPEC *BuildConservedList(CLUSTERINFO *cinfo, int NClus, int *NCons)
   ----------------------------------------------------------------------
   Input:   CLUSTERINFO  *cinfo      Array of cluster info structures
            int          NClus       Number of clusters
   Output:  int          *NCons      Number of conserved residues
   Returns: RESSPEC      *           Allocated array of residue
                                     sepcifications

   Build an array of unique residues specified in the CLUSTERINFO
   structure array.

   04.10.95 Original    By: ACRM
*/
RESSPEC *BuildConservedList(CLUSTERINFO *cinfo, int NClus, int *NCons)
{
   int     i, j,
           MaxCons   = ALLOCQUANTUM;
   RESSPEC *ConsList = NULL;

   *NCons = 0;

   /* Initialise the conserved list                                     */
   if((ConsList=(RESSPEC *)malloc(MaxCons * sizeof(RESSPEC)))==NULL)
      return(NULL);

   /* Step through each cluster                                         */
   for(i=0; i<NClus; i++)
   {
      /* For each residue in this cluster                               */
      for(j=0; j<cinfo[i].NRes; j++)
      {
         /* If this residue is not already stored, then store it        */
         if(InConsList(ConsList, *NCons, 
                       cinfo[i].chain[j], 
                       cinfo[i].resnum[j], 
                       cinfo[i].insert[j]) == (-1))
         {
            /* Allocate more space if necessary                         */
            if(*NCons == MaxCons)
            {
               MaxCons += ALLOCQUANTUM;
               if((ConsList=(RESSPEC *)
                   realloc(ConsList, MaxCons*sizeof(RESSPEC)))==NULL)
               {
                  *NCons = (-1);
                  return(NULL);
               }
            }
            
            /* Store this residue in the conserved list                 */
            ConsList[*NCons].chain  = cinfo[i].chain[j];
            ConsList[*NCons].resnum = cinfo[i].resnum[j];
            ConsList[*NCons].insert = cinfo[i].insert[j];

            /* Increment the number of conserved residues counter       */
            (*NCons)++;
         }
      }
   }

   return(ConsList);
}


/************************************************************************/
/*>int InConsList(RESSPEC *ConsList, int NCons, char *chain, int resnum, 
                  char *insert)
   ----------------------------------------------------------------------
   Input:   RESSPEC    *ConsList    Array of residue specifications
            int        NCons        Number of specs in ConsList
            char       *chain       Chain specification to look for
            int        resnum       Res num specification to look for
            char       *insert      Insert specification to look for
   Returns: int                     Position in ConsList (-1 if not found)

   Test whether the specified residue is in the ConsList

   04.10.95 Original    By: ACRM
   05.10.95 Returns int rather than BOOL
*/
int InConsList(RESSPEC *ConsList, int NCons, char chain, int resnum, 
                char insert)
{
   int i;
   
   for(i=0; i<NCons; i++)
   {
      if((ConsList[i].chain  == chain)  &&
         (ConsList[i].resnum == resnum) &&
         (ConsList[i].insert == insert))
      {
         ConsList[i].flag = TRUE;
         return(i);
      }
   }
   
   return(-1);
}


/************************************************************************/
/*>BOOL MergeAllProperties(int NLoops, LOOPINFO *loopinfo, int clusnum, 
                           RESSPEC *ConsList, int NRes,
                           CLUSTERINFO *clusterinfo)
   --------------------------------------------------------------------
   Input:   int         NLoops       Number of loops in a cluster
            LOOPINFO    *loopinfo    Array of completed structures for
                                     loops in this cluster
            int         clusnum      The number of the cluster in which
                                     we are interested
            RESSPEC     *ConsList    Conserved residue list
            int         NRes         Number of items in ConsList
   Output:  CLUSTERINFO *clusterinfo Compiled data about this cluster
                                     (Memory allocated within this
                                     structure)
   Returns: BOOL                     Success of memory allocations

   Allocate memory in and complete a clusterinfo structure with merged
   property data for the residues ids conserved in any one cluster.

   04.10.95 Original   By: ACRM
   05.10.95 Corrected checking of absolute conservation when not first
            loop (Was coming out with everythings as conserved).
            Removed First parameter and changed CLUSTERINFO structure
            such that each residue has a `First' parameter. Added handling
            of this.
   06.10.95 Added deleted property handling
   09.10.95 Wasn't clearing `first' flag if - occurred first time so -
            information was written over on the next go rather than ORed
   10.10.95 Handles deletable as a separate array rather than as a residue
            property flag. Simplifies various logic!
            Changed USHORT to PROP_T
*/
BOOL MergeAllProperties(PDB *pdb,
                        RESSPEC *ConsList, int NRes,
                        CLUSTERINFO *clusterinfo)
{
   int  i,k;
   char res;
   PDB *p;

   if(NRes==(-1))
      return(FALSE);
   
   /* If there weren't any common residues, return                      */
   if(NRes == 0)
      return(TRUE);

   /* See if we've done the allocations for this structure              */
   if(clusterinfo->resnum == NULL)
   {
      /* Not yet allocated, so allocate memory in the clusterinfo 
         structure to store properties for these residues
      */
      clusterinfo->resnum         = (int    *)malloc(NRes*sizeof(int));
      clusterinfo->chain          = (char   *)malloc(NRes*sizeof(char));
      clusterinfo->insert         = (char   *)malloc(NRes*sizeof(char));
      clusterinfo->ConservedProps = (PROP_T *)malloc(NRes*sizeof(PROP_T));
      clusterinfo->RangeOfProps   = (PROP_T *)malloc(NRes*sizeof(PROP_T));
      clusterinfo->absolute       = (BOOL   *)malloc(NRes*sizeof(BOOL));
      clusterinfo->ConsRes        = (char   *)malloc(NRes*sizeof(char));
      clusterinfo->First          = (BOOL   *)malloc(NRes*sizeof(BOOL));
      clusterinfo->deletable      = (BOOL   *)malloc(NRes*sizeof(BOOL));
      
      
      if((clusterinfo->resnum         == NULL) ||
         (clusterinfo->chain          == NULL) ||
         (clusterinfo->insert         == NULL) ||
         (clusterinfo->ConservedProps == NULL) ||
         (clusterinfo->RangeOfProps   == NULL) ||
         (clusterinfo->absolute       == NULL) ||
         (clusterinfo->First          == NULL) ||
         (clusterinfo->deletable      == NULL) ||
         (clusterinfo->ConsRes        == NULL))
      {
         if(clusterinfo->resnum != NULL)
            free(clusterinfo->resnum);
         if(clusterinfo->chain  != NULL)
            free(clusterinfo->chain);
         if(clusterinfo->insert != NULL)
            free(clusterinfo->insert);
         if(clusterinfo->ConservedProps != NULL)
            free(clusterinfo->ConservedProps);
         if(clusterinfo->RangeOfProps != NULL)
            free(clusterinfo->RangeOfProps);
         if(clusterinfo->ConsRes != NULL)
            free(clusterinfo->ConsRes);
         if(clusterinfo->absolute != NULL)
            free(clusterinfo->absolute);
         if(clusterinfo->First != NULL)
            free(clusterinfo->First);
         if(clusterinfo->deletable != NULL)
            free(clusterinfo->deletable);
         return(FALSE);
      }

      for(i=0; i<NRes; i++)
      {
         clusterinfo->ConsRes[i]   = '-';
         clusterinfo->absolute[i]  = TRUE;
         clusterinfo->First[i]     = TRUE;
         clusterinfo->deletable[i] = FALSE;
      }

      clusterinfo->NRes = NRes;
   }
   
   /* Clear all found flags in the conserved list                       */
   for(i=0; i<NRes; i++)
      ConsList[i].flag = FALSE;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"CA  ",4))
      {
         if((k=InConsList(ConsList, NRes, 
                          p->chain[0], 
                          p->resnum,
                          p->insert[0]))!=(-1))
         {
            res = blThrone(p->resnam);
            
            clusterinfo->chain[k]          = p->chain[0];
            clusterinfo->resnum[k]         = p->resnum;
            clusterinfo->insert[k]         = p->insert[0];

            if(clusterinfo->First[k])
            {
               clusterinfo->ConservedProps[k] = SetProperties(res);
               clusterinfo->RangeOfProps[k]   = SetProperties(res);
               clusterinfo->ConsRes[k]        = res;

               clusterinfo->First[k]          = FALSE;
            }
            else
            {
               clusterinfo->ConservedProps[k] &= SetProperties(res);
               clusterinfo->RangeOfProps[k]   |= SetProperties(res);
               if(res != clusterinfo->ConsRes[k])
                  clusterinfo->absolute[k] = FALSE;
            }
         }
      }
   }

   /* For each residue in our conserved list                            */
   for(i=0; i<NRes; i++)
   {
      /* If we haven't looked at this residue, then it must be deleted
         in this protein. Strictly, the deleted property is not a
         conserved property in the same was as others. In effect it
         represents an additional amino acid type allowed at a given
         position, so it is ORed with both Conserved and RangeOf
         properties.
      */
      if(!ConsList[i].flag)  /* i.e. This residue is deleted            */
      {
         clusterinfo->deletable[i] = TRUE;
      }
   }
   
   return(TRUE);
}




/************************************************************************/
/*>void PrintDeletedResidues(FILE *fp, CLUSTERINFO cinfo, 
                             RESSPEC *ConsList, int NCons)
   -------------------------------------------------------
   Input:   FILE         *fp           Output file pointer
            CLUSTERINFO  cinfo         Information on this cluster
            RESSPEC      *ConsList     Array of residue specifications
                                       forming the conserved residues list
            int          NCons         Number of members of ConsList

   Print any residues which appear in the conserved residues list but
   which haven't been flagged.

   04.10.95 Original    By: ACRM
*/
void PrintDeletedResidues(FILE *fp, CLUSTERINFO cinfo, 
                          RESSPEC *ConsList, int NCons)
{
   int  i;
   BOOL found;
   
   /* Clear all found flags in the conserved list                       */
   for(i=0; i<NCons; i++)
      ConsList[i].flag = FALSE;
   
   /* For each residues in the cluster info structure, call InConsList()
      which (internally) sets the flags in the ConsList array
   */
   for(i=0; i<cinfo.NRes; i++)
   {
      found = (InConsList(ConsList, NCons, 
                          cinfo.chain[i], 
                          cinfo.resnum[i], 
                          cinfo.insert[i]) != (-1));
   }

   /* Run through the conserved list again, printing missing residues   */
   for(i=0; i<NCons; i++)
   {
      if(!ConsList[i].flag)
      {
         fprintf(fp,"%c%3d%c 0xFFFF /deleted/ (-)\n",
                 ConsList[i].chain,
                 ConsList[i].resnum,
                 ConsList[i].insert);
      }
   }      
}

