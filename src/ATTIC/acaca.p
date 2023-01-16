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
