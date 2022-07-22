/*************************************************************************

   Program:    
   File:       resprops.h
   
   Version:    V1.0
   Date:       02.02.96
   Function:   Residue property flags
   
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
   V1.0  02.02.96 Original extracted from decr.h

*************************************************************************/
/* Includes
*/
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines and macros
*/

/* This is the variable type used to store property flags. It must be
   sufficiently wide for all the flags we use. (Currently 2 bytes.)
*/
typedef USHORT PROP_T;

/* Bit-wise flags for residue properties                                */
#define HPHOB_FLAG     0x0001
#define HPHIL_FLAG     0x0002

#define NEGATIVE_FLAG  0x0004
#define POSITIVE_FLAG  0x0008
#define UNCHARGED_FLAG 0x0010

#define AROMATIC_FLAG  0x0020
#define ALIPHATIC_FLAG 0x0040

#define SMALL_FLAG     0x0080
#define MEDIUM_FLAG    0x0100
#define LARGE_FLAG     0x0200

#define GLY_FLAG       0x0400
#define PRO_FLAG       0x0800
#define OTHER_FLAG     0x1000

#define HBOND_FLAG     0x2000
#define NOHBOND_FLAG   0x4000

#define DELETED_FLAG   0x8000

