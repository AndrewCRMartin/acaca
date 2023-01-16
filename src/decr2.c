/*************************************************************************

   Program:    
   File:       decr2.c
   
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

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "resprops.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "decr2.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXPROPAA 21             /* 20 aa's plus - for deletion         */


/************************************************************************/
/* Globals
*/
PROP_T sPropsArray[MAXPROPAA];   /* Properties for the 20 aa's and -    */
char   sResArray[MAXPROPAA];     /* 1-letter codes in same order        */


/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>void PrintSampleResidues(FILE *fp, PROP_T props, BOOL deletable)
   ----------------------------------------------------------------
   Input:   FILE   *fp        Output file pointer
            PROP_T props      Properties flags
            BOOL   deletable  Is the residue deletable?

   Prints sample amino acids which possess a set of properties.

   03.08.95 Original    By: ACRM
   10.08.95 Prints all 20 aas rather than a - if no conserved properties
   06.10.95 Added handling of the deleted property which is effectively
            just an additional amino acid rather than a residue
            property as such
   10.10.95 Deletable now handled with separate parameter
            Changed USHORT to PROP_T
*/
void PrintSampleResidues(FILE *fp, PROP_T props, BOOL deletable)
{
   int    i;

   fprintf(fp,"  (");

   if(props == (PROP_T)0)
   {
      fprintf(fp,"ACDEFGHIKLMNPQRSTVWY");
   }
   else
   {
      for(i=0; i<20; i++)
      {
         if((sPropsArray[i] & props) == props)
         {
            fprintf(fp,"%c",sResArray[i]);
         }
      }
   }

   /* If the deleted flag was set in our copy of the properties, then
      print a -
   */
   if(deletable)
      fprintf(fp,"-");

   fprintf(fp,")");
}


/************************************************************************/
/*>void PrintProps(FILE *fp, PROP_T props, BOOL deletable)
   ------------------------------------------------------
   Input:   FILE   *fp         Output file pointer
            PROP_T props       Properties value
            BOOL   deletable   Is the residue deletable?

   Print the properties associated with the props value as moderately
   verbose, but parsable, text.

   03.08.95 Original    By: ACRM
   06.10.95 Added deleted handling
   10.10.95 Deleted now handled as separate array
            Changed USHORT to PROP_T
*/
void PrintProps(FILE *fp, PROP_T props, BOOL deletable)
{
   /* Note that we and this with everything but the DELETED_FLAG first.
      This effectively makes sure we switch off the DELETED_FLAG before
      making this comparison
   */
   if(props==(PROP_T)0)
   {
      fprintf(fp,"No conserved properties");

      if(deletable)
         fprintf(fp,"/deletable/");

      return;
   }      
   
   if(ISSET(props, GLY_FLAG))
   {
      fprintf(fp,"glycine");
      return;
   }
   
   if(ISSET(props, PRO_FLAG))
   {
      fprintf(fp,"proline");
      return;
   }

   fprintf(fp,"/");

   if(ISSET(props, HPHOB_FLAG))
   {
      fprintf(fp,"hydrophobic/");

      if(ISSET(props, AROMATIC_FLAG))
      {
         fprintf(fp,"aromatic/");

         if(ISSET(props, HBOND_FLAG))
            fprintf(fp,"H-bonding/");
         if(ISSET(props, NOHBOND_FLAG))
            fprintf(fp,"non-H-bonding/");
      }
   }
   else
   {
      if(ISSET(props, UNCHARGED_FLAG))
      {
         fprintf(fp,"uncharged/");
      }

      if(ISSET(props, NEGATIVE_FLAG))
         fprintf(fp,"negative/");
      if(ISSET(props, POSITIVE_FLAG))
         fprintf(fp,"positive/");

      if(!ISSET(props, NEGATIVE_FLAG) &&
         !ISSET(props, POSITIVE_FLAG))
      {
         if(ISSET(props, HPHIL_FLAG))
            fprintf(fp,"hydrophilic/");

         if(ISSET(props, HBOND_FLAG))
            fprintf(fp,"H-bonding/");
         if(ISSET(props, NOHBOND_FLAG))
            fprintf(fp,"non-H-bonding/");
      }
   }
   
   if(!ISSET(props, AROMATIC_FLAG))
   {   
      if(ISSET(props, SMALL_FLAG))
         fprintf(fp,"small/");
      if(ISSET(props, MEDIUM_FLAG))
         fprintf(fp,"medium/");
      if(ISSET(props, LARGE_FLAG))
         fprintf(fp,"large/");
   }
   
   if(ISSET(props, ALIPHATIC_FLAG))
      fprintf(fp,"aliphatic/");

   if(ISSET(props, OTHER_FLAG))
      fprintf(fp,"not glycine or proline/");

   if(deletable)
      fprintf(fp,"deletable/");
}


/************************************************************************/
/*>void InitProperties(void)
   -------------------------
   Initialise static global property flags tables.
   
   01.08.95 Original    By: ACRM
   06.10.95 Added -
*/
void InitProperties(void)
{
   sResArray[0] = 'A';
   SET(sPropsArray[0], HPHOB_FLAG);
   SET(sPropsArray[0], UNCHARGED_FLAG);
   SET(sPropsArray[0], ALIPHATIC_FLAG);
   SET(sPropsArray[0], SMALL_FLAG);
   SET(sPropsArray[0], OTHER_FLAG);
   SET(sPropsArray[0], NOHBOND_FLAG);

   sResArray[1] = 'C';
   SET(sPropsArray[1], HPHOB_FLAG);
   SET(sPropsArray[1], UNCHARGED_FLAG);
   SET(sPropsArray[1], ALIPHATIC_FLAG);
   SET(sPropsArray[1], SMALL_FLAG);
   SET(sPropsArray[1], OTHER_FLAG);
   SET(sPropsArray[1], NOHBOND_FLAG);

   sResArray[2] = 'D';
   SET(sPropsArray[2], HPHIL_FLAG);
   SET(sPropsArray[2], NEGATIVE_FLAG);
   SET(sPropsArray[2], ALIPHATIC_FLAG);
   SET(sPropsArray[2], SMALL_FLAG);
   SET(sPropsArray[2], OTHER_FLAG);
   SET(sPropsArray[2], NOHBOND_FLAG);

   sResArray[3] = 'E';
   SET(sPropsArray[3], HPHIL_FLAG);
   SET(sPropsArray[3], NEGATIVE_FLAG);
   SET(sPropsArray[3], ALIPHATIC_FLAG);
   SET(sPropsArray[3], MEDIUM_FLAG);
   SET(sPropsArray[3], OTHER_FLAG);
   SET(sPropsArray[3], NOHBOND_FLAG);

   sResArray[4] = 'F';
   SET(sPropsArray[4], HPHOB_FLAG);
   SET(sPropsArray[4], UNCHARGED_FLAG);
   SET(sPropsArray[4], AROMATIC_FLAG);
   SET(sPropsArray[4], LARGE_FLAG);
   SET(sPropsArray[4], OTHER_FLAG);
   SET(sPropsArray[4], NOHBOND_FLAG);

   sResArray[5] = 'G';
   SET(sPropsArray[5], HPHOB_FLAG);
   SET(sPropsArray[5], UNCHARGED_FLAG);
   SET(sPropsArray[5], ALIPHATIC_FLAG);
   SET(sPropsArray[5], SMALL_FLAG);
   SET(sPropsArray[5], GLY_FLAG);
   SET(sPropsArray[5], NOHBOND_FLAG);

   sResArray[6] = 'H';
   SET(sPropsArray[6], HPHIL_FLAG);
   SET(sPropsArray[6], POSITIVE_FLAG);
   SET(sPropsArray[6], ALIPHATIC_FLAG);
   SET(sPropsArray[6], LARGE_FLAG);
   SET(sPropsArray[6], OTHER_FLAG);
   SET(sPropsArray[6], HBOND_FLAG);

   sResArray[7] = 'I';
   SET(sPropsArray[7], HPHOB_FLAG);
   SET(sPropsArray[7], UNCHARGED_FLAG);
   SET(sPropsArray[7], ALIPHATIC_FLAG);
   SET(sPropsArray[7], MEDIUM_FLAG);
   SET(sPropsArray[7], OTHER_FLAG);
   SET(sPropsArray[7], NOHBOND_FLAG);

   sResArray[8] = 'K';
   SET(sPropsArray[8], HPHIL_FLAG);
   SET(sPropsArray[8], POSITIVE_FLAG);
   SET(sPropsArray[8], ALIPHATIC_FLAG);
   SET(sPropsArray[8], LARGE_FLAG);
   SET(sPropsArray[8], OTHER_FLAG);
   SET(sPropsArray[8], NOHBOND_FLAG);

   sResArray[9] = 'L';
   SET(sPropsArray[9], HPHOB_FLAG);
   SET(sPropsArray[9], UNCHARGED_FLAG);
   SET(sPropsArray[9], ALIPHATIC_FLAG);
   SET(sPropsArray[9], MEDIUM_FLAG);
   SET(sPropsArray[9], OTHER_FLAG);
   SET(sPropsArray[9], NOHBOND_FLAG);

   sResArray[10] = 'M';
   SET(sPropsArray[10], HPHOB_FLAG);
   SET(sPropsArray[10], UNCHARGED_FLAG);
   SET(sPropsArray[10], ALIPHATIC_FLAG);
   SET(sPropsArray[10], LARGE_FLAG);
   SET(sPropsArray[10], OTHER_FLAG);
   SET(sPropsArray[10], NOHBOND_FLAG);

   sResArray[11] = 'N';
   SET(sPropsArray[11], HPHIL_FLAG);
   SET(sPropsArray[11], UNCHARGED_FLAG);
   SET(sPropsArray[11], ALIPHATIC_FLAG);
   SET(sPropsArray[11], SMALL_FLAG);
   SET(sPropsArray[11], OTHER_FLAG);
   SET(sPropsArray[11], HBOND_FLAG);

   sResArray[12] = 'P';
   SET(sPropsArray[12], HPHIL_FLAG);
   SET(sPropsArray[12], UNCHARGED_FLAG);
   SET(sPropsArray[12], ALIPHATIC_FLAG);
   SET(sPropsArray[12], MEDIUM_FLAG);
   SET(sPropsArray[12], PRO_FLAG);
   SET(sPropsArray[12], NOHBOND_FLAG);

   sResArray[13] = 'Q';
   SET(sPropsArray[13], HPHIL_FLAG);
   SET(sPropsArray[13], UNCHARGED_FLAG);
   SET(sPropsArray[13], ALIPHATIC_FLAG);
   SET(sPropsArray[13], MEDIUM_FLAG);
   SET(sPropsArray[13], OTHER_FLAG);
   SET(sPropsArray[13], HBOND_FLAG);

   sResArray[14] = 'R';
   SET(sPropsArray[14], HPHIL_FLAG);
   SET(sPropsArray[14], POSITIVE_FLAG);
   SET(sPropsArray[14], ALIPHATIC_FLAG);
   SET(sPropsArray[14], LARGE_FLAG);
   SET(sPropsArray[14], OTHER_FLAG);
   SET(sPropsArray[14], NOHBOND_FLAG);

   sResArray[15] = 'S';
   SET(sPropsArray[15], HPHIL_FLAG);
   SET(sPropsArray[15], UNCHARGED_FLAG);
   SET(sPropsArray[15], ALIPHATIC_FLAG);
   SET(sPropsArray[15], SMALL_FLAG);
   SET(sPropsArray[15], OTHER_FLAG);
   SET(sPropsArray[15], HBOND_FLAG);

   sResArray[16] = 'T';
   SET(sPropsArray[16], HPHIL_FLAG);
   SET(sPropsArray[16], UNCHARGED_FLAG);
   SET(sPropsArray[16], ALIPHATIC_FLAG);
   SET(sPropsArray[16], MEDIUM_FLAG);
   SET(sPropsArray[16], OTHER_FLAG);
   SET(sPropsArray[16], HBOND_FLAG);

   sResArray[17] = 'V';
   SET(sPropsArray[17], HPHOB_FLAG);
   SET(sPropsArray[17], UNCHARGED_FLAG);
   SET(sPropsArray[17], ALIPHATIC_FLAG);
   SET(sPropsArray[17], MEDIUM_FLAG);
   SET(sPropsArray[17], OTHER_FLAG);
   SET(sPropsArray[17], NOHBOND_FLAG);

   sResArray[18] = 'W';
   SET(sPropsArray[18], HPHOB_FLAG);
   SET(sPropsArray[18], UNCHARGED_FLAG);
   SET(sPropsArray[18], AROMATIC_FLAG);
   SET(sPropsArray[18], LARGE_FLAG);
   SET(sPropsArray[18], OTHER_FLAG);
   SET(sPropsArray[18], NOHBOND_FLAG);

   sResArray[19] = 'Y';
   SET(sPropsArray[19], HPHOB_FLAG);
   SET(sPropsArray[19], UNCHARGED_FLAG);
   SET(sPropsArray[19], AROMATIC_FLAG);
   SET(sPropsArray[19], LARGE_FLAG);
   SET(sPropsArray[19], OTHER_FLAG);
   SET(sPropsArray[19], HBOND_FLAG);

   sResArray[20] = '-';
   SET(sPropsArray[20], DELETED_FLAG);
}


/************************************************************************/
/*>PROP_T SetProperties(char res)
   ------------------------------
   Input:   char    res        Residue 1-letter code
   Returns: PROP_T             Property flags for residue
                               (0 if not found)

   Set the props variable from a 1-letter code residue by looking it up 
   in the static residue properties table.

   03.08.95 Original    By: ACRM
   08.08.95 Changed so it returns the properties rather than outputting
            them.
   06.10.95 Changed to use MAXPROPAA constant rather than 20
   10.10.95 Changed USHORT to PROP_T
*/
PROP_T SetProperties(char res)
{
   int   i;
   
   for(i=0; i<MAXPROPAA; i++)
   {
      if(sResArray[i] == res)
      {
         return(sPropsArray[i]);
      }
   }
   return((PROP_T)0);
}


