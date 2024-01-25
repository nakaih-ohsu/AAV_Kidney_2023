Files for Barcode-Seq analysis to demultiplex and detect barcodes in a fastq file. 
Script_BarcodeExtraction.cpp is the script file for analysis. 
Sample.fq is a sample FASTQ file to run the script.
List_of_filenames.txt, Primers.txt, Sample_name_SBC_ID_SBC.txt, VBC.txt are Config files.
Sample_name_SBC_ID_SBC.txt is for the analysis of Sample.fq. Use Sample_name_SBC_ID_SBC_RV_RP_DNA.txt or Sample_name_SBC_ID_SBC_RV_RNA.txt for the actual data analysis.
Output files that show how many times VBC sequences are detected in each sample will be created separately for BC1 and BC2 and will be accumulated in "Pool" folder. Output file names are determined by Sample_name_SBC_ID_SBC.txt.
VBC_capsid.txt is a lookup table between VBCs and capsids.
