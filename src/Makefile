CC       = cc
#COPT     = -g -Wall -ansi -I./bioplib -Wno-unused-but-set-variable -Wno-unused-function -Wno-stringop-truncation
COPT     = -O3 -Wall -ansi -I./bioplib -Wno-unused-but-set-variable -Wno-unused-function -Wno-stringop-truncation

EXE       = clan ficl getloops findsdrs
CLANFILES = clan.o acaca.o decr.o decr2.o
FICLFILES = ficl.o acaca.o
GETLFILES = getloops.o
FSDRFILES = FindSDRs.o decr2.o
LIBS      = -lm 
HFILES    = acaca.h decr.h
CBFILES   = bioplib/CalcRMSPDB.o       \
            bioplib/FindResidueSpec.o  \
            bioplib/TermPDB.o          \
            bioplib/fit.o              \
            bioplib/DupePDB.o          \
            bioplib/InStringList.o     \
            bioplib/FreeStringList.o   \
            bioplib/WritePDB.o         \
            bioplib/BuildConect.o      \
            bioplib/PDB2Seq.o          \
            bioplib/hash.o             \
            bioplib/prime.o            \
            bioplib/GetPDBChainLabels.o \
            bioplib/parse.o            \
            bioplib/stringutil.o       \
            bioplib/IndexPDB.o         \
            bioplib/FindAtomInRes.o    \
            bioplib/SelAtPDB.o         \
            bioplib/angle.o            \
            bioplib/throne.o           \
            bioplib/FindNextResidue.o  \
            bioplib/StoreString.o      \
            bioplib/array2.o           \
            bioplib/upstrncmp.o        \
            bioplib/FNam2PDB.o         \
            bioplib/FitPDB.o           \
            bioplib/FitCaPDB.o         \
            bioplib/ParseRes.o         \
            bioplib/phi.o              \
            bioplib/CopyPDB.o          \
            bioplib/KillLeadSpaces.o   \
            bioplib/StringToUpper.o    \
            bioplib/padterm.o          \
            bioplib/GetCGPDB.o         \
            bioplib/OriginPDB.o        \
            bioplib/GetPDBCoor.o       \
            bioplib/ApMatPDB.o         \
            bioplib/MatMult3_33.o      \
            bioplib/ReadPDB.o          \
            bioplib/FindResidue.o      \
            bioplib/fsscanf.o          \
            bioplib/chindex.o          \
            bioplib/PDBHeaderInfo.o    \
            bioplib/strcatalloc.o      \
            bioplib/GetWord.o          \
            bioplib/stringcat.o        \
            bioplib/TranslatePDB.o 

FBFILES =   bioplib/array2.o           \
            bioplib/FindResidue.o      \
            bioplib/GetWord.o          \
	    bioplib/VecDist.o          \
	    bioplib/upstrncmp.o        \
	    bioplib/ReadPDB.o          \
            bioplib/BuildConect.o      \
	    bioplib/ParseRes.o         \
	    bioplib/SelAtPDB.o         \
            bioplib/IndexPDB.o         \
	    bioplib/phi.o              \
	    bioplib/fsscanf.o          \
	    bioplib/chindex.o          \
	    bioplib/WritePDB.o         \
	    bioplib/padterm.o          \
	    bioplib/CopyPDB.o          \
	    bioplib/FindNextResidue.o  \
            bioplib/StoreString.o      \
            bioplib/FreeStringList.o   \
            bioplib/hash.o             \
            bioplib/prime.o            \
            bioplib/stringutil.o       \
            bioplib/PDB2Seq.o          \
            bioplib/throne.o           \
            bioplib/GetPDBChainLabels.o \
            bioplib/PDBHeaderInfo.o    \
            bioplib/strcatalloc.o      \
            bioplib/stringcat.o        \
	    bioplib/angle.o

GBFILES =   bioplib/upstrncmp.o        \
            bioplib/GetWord.o          \
	    bioplib/ReadPDB.o          \
	    bioplib/ParseRes.o         \
	    bioplib/padterm.o          \
	    bioplib/FindNextResidue.o  \
	    bioplib/chindex.o          \
	    bioplib/hash.o             \
	    bioplib/prime.o            \
	    bioplib/fsscanf.o          \
            bioplib/throne.o           \
            bioplib/stringutil.o       \
            bioplib/FindResidue.o      \
            bioplib/GetPDBChainLabels.o \
            bioplib/PDB2Seq.o          \
            bioplib/IndexPDB.o         \
            bioplib/StoreString.o      \
            bioplib/FreeStringList.o   \
            bioplib/BuildConect.o      \
            bioplib/array2.o           \
            bioplib/PDBHeaderInfo.o    \
            bioplib/strcatalloc.o      \
            bioplib/stringcat.o        \
	    bioplib/WritePDB.o

SBFILES =   bioplib/hbond.o            \
	    bioplib/padterm.o          \
	    bioplib/FindNextResidue.o  \
	    bioplib/angle.o            \
	    bioplib/GetFilestem.o      \
	    bioplib/ParseRes.o         \
	    bioplib/OpenStdFiles.o     \
	    bioplib/GetWord.o          \
	    bioplib/fgetsany.o         \
	    bioplib/throne.o           \
	    bioplib/SetExtn.o          \
	    bioplib/chindex.o          \
	    bioplib/FindResidue.o      \
	    bioplib/phi.o              \
	    bioplib/fsscanf.o          \
	    bioplib/WritePDB.o         \
            bioplib/array2.o           \
            bioplib/GetPDBChainLabels.o \
            bioplib/BuildConect.o      \
            bioplib/PDB2Seq.o          \
            bioplib/IndexPDB.o         \
            bioplib/StoreString.o      \
            bioplib/FreeStringList.o   \
            bioplib/stringutil.o       \
	    bioplib/hash.o             \
	    bioplib/prime.o            \
            bioplib/PDBHeaderInfo.o    \
            bioplib/strcatalloc.o      \
            bioplib/stringcat.o        \
	    bioplib/ReadPDB.o          

all : $(EXE)

.c.o : $(HFILES)
	$(CC) $(COPT) -o $@ -c $<

clan : $(CLANFILES) $(CBFILES)
	$(CC) $(COPT) -o clan $(CLANFILES) $(CBFILES) $(LIBS)
ficl : $(FICLFILES) $(FBFILES)
	$(CC) $(COPT) -o ficl $(FICLFILES) $(FBFILES) $(LIBS)
getloops : $(GETLFILES) $(GBFILES)
	$(CC) $(COPT) -o getloops $(GETLFILES) $(GBFILES) $(LIBS)
findsdrs : $(FSDRFILES) $(SBFILES)
	$(CC) $(COPT) -o findsdrs $(FSDRFILES) $(SBFILES) $(LIBS)




clean :
	\rm -f $(CLANFILES) $(FICLFILES) $(GETLFILES) $(FSDRFILES)
	\rm -f $(CBFILES) $(FBFILES) $(GBFILES) $(SBFILES)

install :
	mkdir -p ../bin
	cp $(EXE) ../bin


