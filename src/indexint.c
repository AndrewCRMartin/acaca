/*************************************************************************

   Program:    
   File:       indexint.c
   
   Version:    V1.3
   Date:       13.02.96
   Function:   Index heapsort an int array
   
   Copyright:  (c) SciTech Software 1991-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   This routine uses a heapsort to index a floating point array
   such that arrin[indx[j]] is in ascending order with j.
   It is modified from the FORTRAN version in 'Numerical Recipes'
   Page 233. This version correctly sorts from array element 0
   as opposed to 1 in the FORTRAN version.

**************************************************************************

   Usage:
   ======
   indexf(n,arrin,indx)
   Input:  n      int      Number of elements in array
           arrin  *double   Array to be indexed
   Output: indx   *int     Index array

**************************************************************************

   Revision History:
   =================
   V1.0  23.06.90 Original
   V1.1  01.06.92 ANSIed and autodoc'd
   V1.2  19.07.93 Corrected bug (said j=l+1 rather than j=l+l). Oops!

*************************************************************************/
/* Includes
*/

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
/*>void indexint(int n, int *arrin, int *indx)
   ---------------------------------------------
   Input:  n      int      Number of elements in array
           arrin  *int    Array to be indexed
   Output: indx   *int     Index array
   
   Index a int array by Heapsort.
   
   03.06.90 Original
   01.06.92 ANSIed and autodoc'd
   13.02.96 Rewritten
*/
void indexint(int n, int *arrin, int *indx)
{
   int i, j, l, ir, 
       indxt;
   int q;
   
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
#ifdef DEMO_CODE
#include <stdio.h>
int main(void)
{
   int values[10];
   int  indx[10],
        j;

   values[0] = 1;
   values[1] = 15;
   values[2] = 25;
   values[3] = 12;
   values[4] = 2;
   values[5] = 3;
   values[6] = 26;
   values[7] = 290;
   values[8] = 5;
   values[9] = 7;

   indexint(10,values,indx);

   for(j=0; j<10; j++)
      printf("%d %d\n",indx[j],values[indx[j]]);

   return(0);
}
#endif
