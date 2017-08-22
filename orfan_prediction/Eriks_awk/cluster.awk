#!/usr/local/bin/gawk -f

# Make single-linkage clusters from a list of pairwise links (matches)
# Script was originally made for blast output (see Gerald/THE_NEW_VH_UNIVERSE)

# Requires file 'names' with all sequence names in current dir.
# Most recent version adapted to work with output from selfcomp.crossmatch


BEGIN {

  # seq[seqname] = cluster number of 'seqname'

  while (getline < "names") {
    seq[$1] = ++n		# pointer seq -> cluster
    names[$1] = $1		# all seqs
  }
}

{
  s1 = $2
  s2 = $5

  # Join clusters of s1 and s2 by setting all s2 clusters to the s1 cluster

  clus1 = seq[s1]
  clus2 = seq[s2]

  for (i in names) {
    if (seq[i] == clus2) seq[i] = clus1
  }
}

END {

  # Count cluster sizes
  clus = 0
  for (j=1; j<=n; j++) {
    for (i in names) {
      if (seq[i] == j) {
	clussize[j]++
      }
    }
  }


  clus = 0
  for (j=1; j<=n; j++) {
    first = 1
    for (i in names) {
      if (seq[i] == j && clussize[j] > 1) {
	if (first) {
	  print "Cluster "++clus":"
	  first = 0
	}
	print i
      }
    }
    if (!first) print ""
  }
}
