int main(int argc, char **argv)
;
BOOL ParseCmdLine(int argc, char **argv, char *infile, BOOL *CATorsions)
;
BOOL ReadInputFile(FILE *fp, BOOL CATorsions)
;
BOOL SetupParser(void)
;
BOOL DoCmdLoop(FILE *fp, BOOL CATorsions)
;
BOOL ShowClusters(FILE *fp, REAL **data, int NVec, int VecDim, 
                  int Method, BOOL ShowTable, BOOL ShowDendogram)
;
BOOL HierClus(int NVec, int VecDim, int ClusterMethod, REAL **data, 
              int *ia, int *ib, REAL *crit)
;
int **ClusterAssign(FILE *fp, int NVec, int *ia, int *ib, REAL *crit, 
                    int lev, int *iorder, REAL *critval, int *height)
;
char **ClusterDendogram(FILE *fp, int lev, int *iorder, int *height, 
                        REAL *critval, REAL DivFactor) 
;
BOOL DoClustering(BOOL CATorsions)
;
void Usage(void)
;
void CreateDefaultScheme(int maxres)
;
BOOL InsertIorder(int *iorder, int lev, int cluster, int parent)
;
void WriteHeader(FILE *fp, int Method, int NVec, int VecDim, int *Scheme)
;
BOOL WriteResults(FILE *fp, int *clusters, int NClus, REAL **data, 
                  int NVec, int VecDim, REAL *crit, BOOL PostClus)
;
int FindNumTrueClusters(REAL *crit, int lev, int VecDim)
;
void CleanUp(void)
;
void WriteClusData(FILE *fp, int NVec, int VecDim, REAL **data)
;
DATALIST *FindMedian(int *clusters, REAL **data, int NVec, int VecDim, 
                     int ClusNum, int *NMemb)
;
void FillClusterArray(int **clusters, int NVec, int NClus,
                      int *TheClusters)
;
REAL RmsPDB(PDB *pdb1, PDB *pdb2, int length)
;
REAL RmsCAPDB(PDB *pdb1, PDB *pdb2, int length)
;
REAL MaxCADeviationPDB(PDB *pdb1, PDB *pdb2, int length)
;
REAL MaxCBDeviationPDB(PDB *pdb1, PDB *pdb2, int length)
;
int RenumClusters(int *clusters, int NVec)
;
int PostCluster(FILE *fp, int *clusters, REAL **data, int NVec, 
                int VecDim, REAL *crit, int NClus)
;
BOOL TestMerge(DATALIST *loop1, DATALIST *loop2, REAL *rms, REAL *CADev,
               REAL *CBDev)
;
void DoMerge(FILE *fp, int i, DATALIST *loop1, int j, DATALIST *loop2, 
             REAL rms, REAL CADev, REAL CBDev, int *NewNumbers, int NClus)
;
DATALIST *FindLoop(int *clusters, int NVec, int ClusNum, int loopnum)
;
BOOL DefineCriticalResidues(FILE *fp, int *clusters, REAL **data, 
                            int NVec, int VecDim, REAL *crit, int NClus)
;
