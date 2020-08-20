# CCHFV-PrimerChecker
Bioinformatic tool to evaluate detection efficacy of existing molecular assays and to propose specific sets of assays to detect as much as possible CCHFV strains from all geographic regions.

**Reference:** Gruber CEM, Bartolini B, Castilletti C, et al. Geographical Variability Affects CCHFV Detection by RT-PCR: A Tool for In-Silico Evaluation of Molecular Assays. Viruses. 2019;11(10):953. Published 2019 Oct 16. doi:10.3390/v11100953

Instructions:

The software is ready-to-use and no need any compilation.

- If you use CCHFV-PrimerChecker in your PC for the first time:

1) clone the folders "bins" and "inputfiles" in your local directory;
2) open a shell (Linux) or a prompt (Windows) in that directory;
3) digit: "python bins/CCHFVPrimerChecker.py";
4) search for the file "PrimersMatches_RESULTS.xls" in your local directory: here you can find the result of all the anayzed in   the reference work.

- If you want to compare a new assay with the previously analyzed:

    a) open the Primers_table.csv file in folder "inputfiles";
    
    b) compile the file adding all the informations of the primers/probes of the new assay;

    c) follow the instructions from 2) to 4).

- If you want to add a new sequence to the S segment database:

    a) open the CCHFVsequences.fasta file in folder "inputfiles";   
    
    b) add the new sequence to the file;

    c) open the CCHFVsequencesAligned.fasta file in folder "inputfiles";

    d) align the CCHFVsequencesAligned.fasta using your preferred aligner tool;

    e) save the new file as "CCHFVsequencesAligned.fasta";

    f) chek the refence sequence "NC_005302_S_Aligned.fasta" is well aligned with the new "CCHFVsequencesAligned.fasta";
    
    g) follow the instructions from 2) to 4).
