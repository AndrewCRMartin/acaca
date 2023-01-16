BOOL FindNeighbourProps(PDB *pdb, PDB *start, PDB *stop, int clusnum,
                        LOOPINFO *loopinfo)
;
BOOL ResidueContact(PDB *p_start, PDB *p_stop, PDB *q_start, PDB *q_stop,
                    REAL dist)
;
void FillLoopInfo(LOOPINFO *loopinfo)
;
BOOL MergeProperties(int NLoops, LOOPINFO *loopinfo, int clusnum,
                     CLUSTERINFO *clusterinfo)
;
void BlankClusterInfo(CLUSTERINFO *clusterinfo)
;
void BlankLoopInfo(LOOPINFO *loopinfo)
;
int FlagCommonResidues(int NLoops, LOOPINFO *loopinfo, int clusnum)
;
void CleanLoopInfo(LOOPINFO *loopinfo, int NMembers)
;
void CleanClusInfo(CLUSTERINFO *cinfo)
;
void PrintMergedProperties(FILE *fp, int clusnum, CLUSTERINFO cinfo,
                           int NMembers)
;
RESSPEC *BuildConservedList(CLUSTERINFO *cinfo, int NClus, int *NCons)
;
int InConsList(RESSPEC *ConsList, int NCons, char chain, int resnum, 
                char insert)
;
BOOL MergeAllProperties(PDB *pdb,
                        RESSPEC *ConsList, int NRes,
                        CLUSTERINFO *clusterinfo)
;
void PrintDeletedResidues(FILE *fp, CLUSTERINFO cinfo, 
                          RESSPEC *ConsList, int NCons)
;
