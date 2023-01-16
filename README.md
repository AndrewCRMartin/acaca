ACACA - Automatic Canonical Assignment by Cluster Analysis
==========================================================

Reads a set of PDB files and a control file which specifies the loops to be
analysed. For each loop writes a file containing the CA pseudo-torsions
for the loops. Shorter loops are padded with torsions of 9999.0

For each loop in turn:
- Performs cluster analysis.
- Process results of cluster analysis to identify a sensible number of clusters.
- Assigns each conformation to a cluster.
- Looks for features in each cluster. e.g. Conserved buried residues

Cluster analysis must provide the centre of the cluster, the dimensions of 
the cluster and the distance to the nearest other cluster such that new
structures can be scanned against the clusters.

Three programs:

## CLAN - CLuster ANalysis

Performs cluster analysis on a loop in a set of PDB files. Generates
information on the clusters.

Takes an input file of the following syntax:

| Keyword     | Parameters                 | Details                                    |
| ----------- | -------------------------- | ------------------------------------------ |
| METHOD      | `clustering-method`        | Ward, single, multiple, etc.               |
| OUTPUT      | `outfile`                  | or stdout if not specified                 |
| MAXLENGTH   | `length`                   | Max length of a loop in analysis           |
| SCHEME      | `insert-scheme`            | Order in which positions in the loop should be assigned from the actual loop i.e. where insertions should be considered |
| DENDOGRAM   |                            | Show the clustering tree dendogram         |
| TABLE       |                            | Show the cluster table                     |
| DATA        |                            | Show the data which is used for clustering |
| POSTCLUSTER | `cutoff`                   | Specify RMS cutoff for post-cluster        |
| LOOP        | `pdb` `startres` `lastres` | Multiple records (must come last)          |

The output file contains the METHOD, MAXLENGTH and SCHEME information as
well as the clustering data which includes the centre and size of each
cluster and the distance to the nearest neighbouring cluster.


## FICL - FInd CLuster

Takes a loop in a PDB file and matches it against a set of clusters to
find which cluster (if any) it matches.

Is run with the following syntax:

```
ficl datafile pdb startres lastres
```

where:

- `pdb startres lastres`    is the loop to be tested
- `datafile`                is the output file from CLAN

FICL will pick up the cluster method from the output of CLAN

CLAN *must* be run with TABLE and DATA switched on!


## FindSDRs - Find structure determining (key) residues

Takes the output of CLAN and identifies key residues for the clusters

Is run with the following syntax:

```
findsdrs [-k] [clanfile [outfile]]
```

