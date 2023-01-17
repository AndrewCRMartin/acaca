/*************************************************************************

   Program:    GetLoops
   File:       getloops.c
   
   Version:    V3.6
   Date:       09.01.96
   Function:   Get loops specified in a clan input file
   
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
   V1.0  03.07.95 Original
   V3.4  10.09.95 Skipped
   V3.5  06.11.95 Skipped
   V3.6  09.01.96 Filenames have start and end residues


*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define MAXRES  16

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL GetLoop(char *filename, char *firstres, char *lastres);
char *MakeOutFilename(char *filename);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   03.07.95 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   char filename[MAXBUFF],
        buffer[MAXBUFF],
        word[MAXBUFF],
        firstres[MAXRES],
        lastres[MAXRES],
        *p;
   FILE *fp;

   if(argc != 2 || argv[1][0] == '-')
   {
      Usage();
      exit(0);
   }
   
   if((fp=fopen(argv[1],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open CLAN input file: %s\n",argv[1]);
      exit(1);
   }
   
   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);

      if((p = blGetWord(buffer,word,MAXBUFF))!=NULL)
      {
         if(!blUpstrncmp(word,"LOOP",4))
         {
            p = blGetWord(p,filename,MAXBUFF);
            p = blGetWord(p,firstres,MAXRES);
            p = blGetWord(p,lastres,MAXRES);
            
            if(!GetLoop(filename,firstres,lastres))
            {
               fprintf(stderr,"Unable to open ore read file: %s\n",
                       filename);
            }
         }
      }
   }

   fclose(fp);

   return(0);
}


/************************************************************************/
/*>BOOL GetLoop(char *filename, char *firstres, char *lastres)
   -----------------------------------------------------------
   03.07.95 Original    By: ACRM
   09.01.96 Filename now contains start and end residues
*/
BOOL GetLoop(char *filename, char *firstres, char *lastres)
{
   FILE *fp,     *fpout;
   PDB  *pdb,    *p;
   int  res1,    res2,
        natoms;
   char chain1,  chain2,
        insert1, insert2,
        *outfile,
        namebuffer[MAXBUFF];
   BOOL InLoop = FALSE,
        InLast = FALSE;
   
   if((fp=fopen(filename,"r")) == NULL)
      return(FALSE);

   if((pdb=blReadPDBAtoms(fp, &natoms))==NULL)
   {
      fclose(fp);
      return(FALSE);
   }
   fclose(fp);

   outfile = MakeOutFilename(filename);
   sprintf(namebuffer,"%s-%s-%s",outfile,firstres,lastres);
   if((fpout=fopen(namebuffer,"w")) == NULL)
   {
      FREELIST(pdb, PDB);
      fclose(fp);
      return(FALSE);
   }

   blParseResSpec(firstres, &chain1, &res1, &insert1);
   blParseResSpec(lastres, &chain2, &res2, &insert2);
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    == res1) &&
         (p->chain[0]  == chain1) &&
         (p->insert[0] == insert1))
         InLoop = TRUE;

      if((p->resnum    == res2) &&
         (p->chain[0]  == chain2) &&
         (p->insert[0] == insert2))
         InLast = TRUE;

      if(InLast)
      {
         if((p->resnum    != res2) ||
            (p->chain[0]  != chain2) ||
            (p->insert[0] != insert2))
            InLoop = FALSE;
      }
      
      if(InLoop)
      {
         blWritePDBRecord(fpout, p);
      }
   }

   fclose(fpout);
   FREELIST(pdb, PDB);

   if(!InLast)
   {
      fprintf(stderr,"%s skipped! Last residue (%s) not found\n",filename,lastres);
      unlink(namebuffer);
   }
   
   return(TRUE);
}

/************************************************************************/
/*>char *MakeOutFilename(char *filename)
   -------------------------------------
   03.07.95 Original    By: ACRM
*/
char *MakeOutFilename(char *filename)
{
   char *p, *q;
   
   q = filename;
   
   while((p=strchr(q,'/')) != NULL)
      q = p+1;

   return(q);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   03.07.95 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nGetLoops V1.0 (c) Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: getloops <filename>\n");

   fprintf(stderr,"\nExtracts loops from a set of PDB files specified in \
a file of the form\n");
   fprintf(stderr,"used as input to CLAN.\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"N.B. If the input files are in teh current directory, \
they will be\n");
   fprintf(stderr,"OVER-WRITTEN by the loop files.\n");
}
