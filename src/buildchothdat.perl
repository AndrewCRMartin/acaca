#!/bin/perl
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 1995
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   Phone:      (Home) +44 (0)1372 275775
#               (Work) +44 (0)171 387 7050 X 3284
#   EMail:      INTERNET: martin@biochem.ucl.ac.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************

if($#ARGV < 0)
{
    print "\nBuildChothDat V1.0 (c) 1995, Dr. Andrew C.R. Martin, UCL\n";
    print "Usage: buildchothdat <loopid> <clanfile>\n";

    print "\nBuilds a Chothia program loop template file from a Clan output file\n\n";

    exit;
}


$indata = 0;
$loop = $ARGV[0];
shift;

while(<>)
{
    if(/END/)
    {
        $indata = 0;
    }
    elsif($indata)
    {
        if(/CLUSTER/)
        {
            # Parse out the cluster number & length
            ($j1,$clusnum,$j1,$j1,$length,$j1) = split;
            chop $length;
            # Write a header line
            printf "LOOP $loop $clusnum $length\n";
            printf "SOURCE Cluster analysis\n";
        }
        else
        {
            if(/WARNING/)
            {
                ;
            }
            else
            {
                chop;
                chop;

                # Parse out the chain name, resid
                $chain = substr($_,0,1);
                $resid = substr($_,1,4);

                # Remove leading spaces
                while(substr($resid,0,1) eq " ")
                {
                    $resid = substr($resid,1);
                }

                # Remove trailing spaces
                while(length($resid) && 
                      substr($resid,length($resid)-1,1) eq " ")
                {
                    chop $resid;
                }


                # Parse out the allowed residues
                ($j1,$allowed) = split(/\(/);

                if($allowed ne "-")
                {
                    # Write the residue record
                    printf "%s%s %s\n", $chain, $resid, $allowed;
                }
            }
        }
    }
    elsif(/BEGIN CRITICALRESIDUES/)
    {
        $indata = 1;
    }
}
