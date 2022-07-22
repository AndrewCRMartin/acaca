/*************************************************************************

   Program:    
   File:       decr.h
   
   Version:    V3.6
   Date:       09.01.96
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
   V0.1  01.08.95 Original
   V3.4  10.10.95 Various changes to make deleted residues work with 
                  residues conserved in at least one cluster
   V3.5  06.11.95 Skipped
   V3.6  09.01.96 Skipped

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/macros.h"

#include "resprops.h"

/************************************************************************/
/* Defines and macros
*/

#define ALLOCQUANTUM 16
#define CONTACTDIST  ((REAL)4.0)

typedef struct
{
   PDB    **contacts,
          **residues;
   PROP_T *ContactProps,
          *ResProps;
   char   *AALoop,
          *AAContact;
   int    length,
          clusnum,
          ncontacts;
   BOOL   *ResFlag,
          *ContactFlag;
}  LOOPINFO;

typedef struct
{
   int    *resnum;
   char   *chain,
          *insert;
   PROP_T *ConservedProps,
          *RangeOfProps;
   int    NRes,                  /* Number of common residue ids        */
          length;                /* Length of loop itself               */
   BOOL   *absolute,
          *First,
          *deletable;
   char   *ConsRes;
}  CLUSTERINFO;

typedef struct
{
   int  resnum;
   char chain,
        insert;
   BOOL flag;
}  RESSPEC;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

