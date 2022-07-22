int main(int argc, char **argv)
;
int ReadClanFile(FILE *in, int *NLoops)
;
BOOL ReadAssignments(FILE *fp)
;
void FreeGlobalStorage(int nclus)
;
void StorePDBNameCluster(char *inbuff, int LoopNum)
;
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *KeepSA)
;
void Usage(void)
;
void BlankTemplates(int nclus)
;
BOOL ExpandTemplateArrays(CLUSINFO *ClusInfo)
;
BOOL ReadTemplates(FILE *in)
;
BOOL FindSDRs(int nclus, int nloops, BOOL KeepSA)
;
void ReportSDRs(FILE *out, int nclus)
;
void Report(CLUSINFO *ClusInfo, int residx, char *reason)
;
void FillOoiData(void)
;
BOOL IsInRange(char *resspec, char *firstres, char *lastres)
;
PDB *ReadPDBAsSA(char *filename, BOOL KeepSAFile)
;
void MarkPartners(CLUSINFO *ClusInfo, PDB *pdb, PDB *res, 
                  char *firstres, char *lastres)
;
BOOL MakeSCContact(PDB *res1, PDB *res2)
;
BOOL MarkHPhob(CLUSINFO *ClusInfo, int clusnum, int nloops, 
               BOOL KeepSA)
;
BOOL MarkHBonders(CLUSINFO *ClusInfo, int clusnum, int nloops)
;
SDRLIST *InSDRList(SDRLIST *sdrlist, char chain, int resnum, char insert)
;
BOOL FillSDRsForCluster(SDRLIST *sdrlist, int clusnum, int nloops)
;
BOOL ReportUnifiedSDRs(FILE *out, int nclus, int nloops)
;
void PrintSDRList(FILE *out, SDRLIST *sdrlist)
;
void indexint(int n, int *arrin, int *indx)
;
void FlagNonInformativeSDRs(int nclus)
;
BOOL ValueIsAdded(SDRLIST *s1, SDRLIST *s2)
;
void FlagRogueClusters(int nclus, int nloops)
;
BOOL IsRogue(int clus, int LargestClus)
;
BOOL IsCisProline(CLUSINFO *ClusInfo, int clusnum, int resoffset, 
                  int nloops)
;
