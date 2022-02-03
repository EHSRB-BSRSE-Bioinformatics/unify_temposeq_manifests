# unify_temposeq_manifests
Code to standarize the format of the TempO-Seq probe manifests and fill in any missing ID values where possible

Output files consist of Probe ID, Probe Name, Gene Symbol, Ensembl ID, Entrez (NCBI) ID, Probe Sequence, and Targeted Transcripts.
Note that the Probe Name consists of the Probe ID + Gene Symbol, *but* it is not modified from BioSpyder's original IDs, whereas the gene symbols have occasionally been updated and those updates are captured in this manifest.

Cases where the Ensembl/NCBI ID are missing include genes that have been withdrawn from repositories, and cases where it was unclear precisely which gene is being targeted.


## Stats

### Human S1500 1.2 (2981 total)

Before: no Ensembl IDs or Entrez IDs

After: 34 blank Ensembl IDs, 54 blank Entrez IDs

### Human S1500 2.0 (3386 total)

Before: 69 blank Ensembl IDs, 132 blank Entrez IDs

After: 3 blank Ensembl IDs, 2 blank Entrez IDs

### Human WT 2.0 (22538 total)

Before: 122 blank Ensembl IDs, 77 blank Entrez IDs

After: 34 blank Ensembl IDs, 29 blank Entrez IDs

