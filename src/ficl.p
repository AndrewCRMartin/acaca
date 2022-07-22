int main(int argc, char **argv)
;
BOOL ParseCmdLine(int argc, char **argv, char *datafile, char *pdbfile, 
                  char *startres, char *lastres, BOOL *CATorsions,
                  BOOL *Verbose)
;
REAL **ReadClusterFile(char *datafile, BOOL CATorsions, int *pMethod,
                       int *pNData, int *pVecLength, 
                       CLUSTER **ppClusters, int *NClusters,
                       CLUSTER **ppMedians,  int *NMedians)
;
BOOL ReadData(FILE *fp, REAL **data, int NLoops, int VecLength)
;
BOOL ReadHeader(FILE *fp, int *pMethod, int *pNLoops, int *pMaxLen)
;
REAL **AllocateDataArrays(int NLoops, int VecLength, CLUSTER **ppClusters)
;
void Usage(void)
;
void CleanUp(REAL **data1, int NData1, int VecLen1, 
             REAL **data2, int NData2, int VecLen2)
;
int ReadClusters(FILE *fp, CLUSTER *clusters)
;
int MatchCluster(REAL **data, int NData, int VecLength, 
                 CLUSTER *clusters, int NClusters, REAL *LoopData, 
                 BOOL CATorsions, int method, BOOL *pError)
;
int ConfirmCluster(REAL **data, int NVec, int VecLen, 
                   CLUSTER *clusters, int TheCluster, REAL *vector,
                   BOOL *pError)
;
int InClusterBounds(REAL **data, int NVec, int VecDim, CLUSTER *clusters,
                    int ClusNum, REAL *vector)
;
REAL MinDistInCluster(REAL **data, int NVec, int VecLen, 
                      CLUSTER *clusters, REAL *vector, int ClusNum)
;
int FindNearestMedian(REAL **data, int NVec, int VecLen, 
                      CLUSTER *clusters, int NClusters, REAL *vector)
;
REAL *FindMedian(REAL **data, int NVec, int VecLen, CLUSTER *clusters, 
                 int ClusNum)
;
int ReadMedians(FILE *fp, CLUSTER **ppMedians)
;
void PrintClusterInfo(int TheCluster, CLUSTER *MedianData, int NMedians,
                      REAL dist, BOOL Verbose)
;
