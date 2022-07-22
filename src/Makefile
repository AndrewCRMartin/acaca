CC       = cc
COPT     = -g -Wall -ansi -I./bioplib
PROTOIZE = /home/bsm/martin/bin/ansi -p

CLANFILES = clan.o acaca.o decr.o decr2.o
FICLFILES = ficl.o acaca.o
GETLFILES = getloops.o
FSDRFILES = FindSDRs.o decr2.o
LIBS      = -lm 
PFILES    = clan.p acaca.p ficl.p decr.p FindSDRs.p decr2.p
HFILES    = acaca.h decr.h
CBFILES   = bioplib/CalcRMSPDB.c       \
            bioplib/FindResidueSpec.c  \
            bioplib/TermPDB.c          \
            bioplib/fit.c              \
            bioplib/DupePDB.c          \
            bioplib/InStringList.c     \
            bioplib/WritePDB.c         \
            bioplib/parse.c            \
            bioplib/FindAtomInRes.c    \
            bioplib/SelAtPDB.c         \
            bioplib/angle.c            \
            bioplib/throne.c           \
            bioplib/FindNextResidue.c  \
            bioplib/StoreString.c      \
            bioplib/array2.c           \
            bioplib/upstrncmp.c        \
            bioplib/FNam2PDB.c         \
            bioplib/FitPDB.c           \
            bioplib/FitCaPDB.c         \
            bioplib/ParseRes.c         \
            bioplib/phi.c              \
            bioplib/CopyPDB.c          \
            bioplib/KillLeadSpaces.c   \
            bioplib/StringToUpper.c    \
            bioplib/padterm.c          \
            bioplib/GetCGPDB.c         \
            bioplib/OriginPDB.c        \
            bioplib/GetPDBCoor.c       \
            bioplib/ApMatPDB.c         \
            bioplib/MatMult3_33.c      \
            bioplib/ReadPDB.c          \
            bioplib/FindResidue.c      \
            bioplib/fsscanf.c          \
            bioplib/chindex.c          \
            bioplib/TranslatePDB.c 

FBFILES =   bioplib/array2.c           \
            bioplib/GetWord.c          \
	    bioplib/VecDist.c          \
	    bioplib/upstrncmp.c        \
	    bioplib/ReadPDB.c          \
	    bioplib/ParseRes.c         \
	    bioplib/SelAtPDB.c         \
	    bioplib/phi.c              \
	    bioplib/fsscanf.c          \
	    bioplib/chindex.c          \
	    bioplib/WritePDB.c         \
	    bioplib/padterm.c          \
	    bioplib/CopyPDB.c          \
	    bioplib/FindNextResidue.c  \
	    bioplib/angle.c

GBFILES =   bioplib/upstrncmp.c        \
            bioplib/GetWord.c          \
	    bioplib/ReadPDB.c          \
	    bioplib/ParseRes.c         \
	    bioplib/padterm.c          \
	    bioplib/FindNextResidue.c  \
	    bioplib/chindex.c          \
	    bioplib/fsscanf.c          \
	    bioplib/WritePDB.c

SBFILES =   bioplib/hbond.c            \
	    bioplib/padterm.c          \
	    bioplib/FindNextResidue.c  \
	    bioplib/angle.c            \
	    bioplib/GetFilestem.c      \
	    bioplib/ParseRes.c         \
	    bioplib/OpenStdFiles.c     \
	    bioplib/GetWord.c          \
	    bioplib/fgetsany.c         \
	    bioplib/throne.c           \
	    bioplib/SetExtn.c          \
	    bioplib/chindex.c          \
	    bioplib/FindResidue.c      \
	    bioplib/phi.c              \
	    bioplib/fsscanf.c          \
	    bioplib/ReadPDB.c          


all : clan ficl getloops findsdrs

.c.o : $(HFILES)
	$(CC) $(COPT) -o $@ -c $<




clan : $(CLANFILES)
	$(CC) $(COPT) -o clan $(CLANFILES) $(CBFILES) $(LIBS)
ficl : $(FICLFILES)
	$(CC) $(COPT) -o ficl $(FICLFILES) $(FBFILES) $(LIBS)
getloops : $(GETLFILES)
	$(CC) $(COPT) -o getloops $(GETLFILES) $(GBFILES) $(LIBS)
findsdrs : $(FSDRFILES)
	$(CC) $(COPT) -o findsdrs $(FSDRFILES) $(SBFILES) $(LIBS)




clean :
	rm -f $(CLANFILES) $(FICLFILES) $(GETLFILES) $(FSDRFILES)

protos : $(PFILES)

.c.p :
	$(PROTOIZE) $< $@
