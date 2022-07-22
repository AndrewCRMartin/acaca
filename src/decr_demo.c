#define PRINT_LOOP_PROPS








/*************************************************************************
***                                                                    ***
***                             DEMO CODE                              ***
***                                                                    ***
*************************************************************************/

#include <stdio.h>
/************************************************************************/
BOOL ProcessPDBFile(char *pdbfile, LOOPINFO *loopinfo)
{
   PDB  *pdb, *start, *stop;
   FILE *fp;
   int  natom, i,
        clusnum = 4;
   
   
   fp = fopen(pdbfile,"r");
   pdb = ReadPDBAtoms(fp,&natom);
   fclose(fp);

   for(start=pdb,natom=0; start!=NULL; NEXT(start))
   {
      if(!strncmp(start->atnam,"N   ",4))
         natom++;
      if(natom == 24)
         break;
   }
   for(stop=pdb,natom=0; stop!=NULL; NEXT(stop))
   {
      if(!strncmp(stop->atnam,"N   ",4))
         natom++;
      if(natom == 34)
         break;
   }

   return(FindNeighbourProps(pdb, start, stop, clusnum, loopinfo));
}



/************************************************************************/
int main(void)
{
   LOOPINFO    loopinfo[4];
   int         i, j;
   CLUSTERINFO cinfo;

   InitProperties();


   if(!ProcessPDBFile("/pdb/p2hfl.pdb", &(loopinfo[0])))
   {
      fprintf(stderr,"No memory for p2hfl");
      return(1);
   }
   if(!ProcessPDBFile("/pdb/p2fbj.pdb", &(loopinfo[1])))
   {
      fprintf(stderr,"No memory for p2fbj");
      return(1);
   }
   if(!ProcessPDBFile("/pdb/p1for.pdb", &(loopinfo[2])))
   {
      fprintf(stderr,"No memory for p1for");
      return(1);
   }
   if(!ProcessPDBFile("/pdb/p1baf.pdb", &(loopinfo[3])))
   {
      fprintf(stderr,"No memory for p1baf");
      return(1);
   }

   /* Blank the cluster data                                            */
   BlankClusterInfo(&cinfo);
   
   if(!MergeProperties(4, loopinfo, 4, &cinfo))
   {
      fprintf(stderr,"MergeProperties() failed\n");
      return(1);
   }

#ifdef PRINT_LOOP_PROPS
   for(j=0; j<4; j++)
   {
      printf("Loop properties:\n");
      for(i=0;i<loopinfo[j].length;i++)
      {
         printf("%c%4d%c %c 0x%04x %d\n",
                loopinfo[j].residues[i]->chain[0],
                loopinfo[j].residues[i]->resnum,
                loopinfo[j].residues[i]->insert[0],
                loopinfo[j].AALoop[i],
                loopinfo[j].ResProps[i],
                (int)loopinfo[j].ResFlag[i]
                );
      }
      printf("Contact properties:\n");
      for(i=0;i<loopinfo[j].ncontacts;i++)
      {
         printf("%c%4d%c %c 0x%04x %d\n",
                loopinfo[j].contacts[i]->chain[0],
                loopinfo[j].contacts[i]->resnum,
                loopinfo[j].contacts[i]->insert[0],
                loopinfo[j].AAContact[i],
                loopinfo[j].ContactProps[i],
                (int)loopinfo[j].ContactFlag[i]
                );
      }
   }
#endif

   printf("MERGED PROPERTIES:\n");
   for(i=0;i<cinfo.NRes;i++)
   {
/*
      printf("%c%3d%c   0x%04x   0x%04x\n",
             cinfo.chain[i],
             cinfo.resnum[i],
             cinfo.insert[i],
             cinfo.ConservedProps[i],
             cinfo.RangeOfProps[i]);
*/
      printf("%c%3d%c 0x%04x ",
             cinfo.chain[i],
             cinfo.resnum[i],
             cinfo.insert[i],
             cinfo.ConservedProps[i]);
      PrintProps(stdout,cinfo.ConservedProps[i]);
      PrintSampleResidues(stdout,cinfo.ConservedProps[i]);
      printf("\n");
   }

   return(0);
}
#endif


