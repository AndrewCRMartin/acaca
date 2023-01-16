/*************************************************************************

   Program:    acaca suite
   File:       acaca.h
   
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
   V3.0  13.09.95 Added gDoDistance
   V3.4  10.09.95 Skipped
   V3.5  06.11.95 Added gStringList
   V3.6  09.01.96 Skipped
   V3.7  14.03.96 gPClusCut[] now 3 long rather than 2
   V3.7a 30.01.09 Increased MAXLOOPLEN and added comment

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <values.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/parse.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/angle.h"
#include "bioplib/array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF              160
#define MAXLOOPLEN           40  /* Note this must be actual maximum +2 */
#define TORPERRES            3
#define MAXLOOPID            32

#define RMSCUT               1.0
#define MAXDEV               1.5
#define MAXCBDEV             1.9

#define DUMMY                ((REAL)10.0)  /* For non-existent torsions */
#define DUMMY2               ((REAL)100.0) /* For non-existent dists.   */
#define INF                  ((REAL)MAXDOUBLE)

typedef struct _datalist
{
   struct _datalist *next;
   PDB    *allatompdb,      /* The native PDB linked list from the file */
          *pdbloop,         /* Start of loop in pdbstart linked list    */
          *torsionpdb;      /* C-alpha or backbone atoms of whole struc */
   REAL   torsions[MAXLOOPLEN * TORPERRES],
          angles[MAXLOOPLEN],
          dist[MAXLOOPLEN];
   int    length;
   char   loopid[MAXBUFF],
          start[16],
          end[16];
}  DATALIST;

typedef struct _cluster
{
   int    clusnum;
   char   loopid[MAXLOOPID];
}  CLUSTER;

/************************************************************************/
/* Globals
*/
#ifdef MAIN /* ------------------- Define externals ------------------- */
int        gMaxLoopLen    = 0,
           gScheme[MAXLOOPLEN],
           gClusterMethod = 1;
BOOL       gDoDendogram   = FALSE,
           gDoTable       = FALSE,
           gDoData        = FALSE,
           gDoCritRes     = FALSE,
           gDoDistance    = FALSE,    /* Handle dists. in clustering    */
           gDoAngles      = FALSE,    /* Handle angles in clustering    */
           gCATorsions    = FALSE;    /* Do CA pseudo torsions          */
FILE       *gOutfp        = NULL;
DATALIST   *gDataList     = NULL;
STRINGLIST *gStringList   = NULL;
REAL       gPClusCut[3];
#else       /* ----------------- Reference externals ------------------ */
extern int        gMaxLoopLen,
                  gScheme[MAXLOOPLEN],
                  gClusterMethod;
extern BOOL       gDoDendogram,
                  gDoTable,
                  gDoData,
                  gDoCritRes,
                  gDoDistance,
                  gDoAngles,
                  gCATorsions;
extern FILE       *gOutfp;
extern DATALIST   *gDataList;
extern STRINGLIST *gStringList;
extern REAL       gPClusCut[2];
#endif

/************************************************************************/
/* Prototypes
*/

BOOL SetClusterMethod(char *method)
;
BOOL SetOutputFile(char *filename)
;
BOOL HandleLoopSpec(char *filename, char *start, char *end, 
                    BOOL CATorsions, BOOL Verbose)
;
BOOL FindCAResidues(PDB *pdbca, char chain1, int resnum1, char insert1,
                    char chain2, int resnum2, char insert2,
                    PDB **pp_start, PDB **pp_end)
;
BOOL StoreTorsions(PDB *allatompdb, PDB *pdb, PDB *p_start, PDB *p_end, 
                   char *filename, char *start, char *end)
;
BOOL FindBBResidues(PDB *pdbbb, char chain1, int resnum1, char insert1,
                    char chain2, int resnum2, char insert2,
                    PDB **pp_start, PDB **pp_end)
;
REAL **ConvertData(DATALIST *indata, int *NData, BOOL CATorsions)
;
void PrintArray(REAL **data, int NData, int width)
;
