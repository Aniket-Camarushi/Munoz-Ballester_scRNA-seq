# Munoz-Ballester_scRNA-seq

This branch contain my experience with scRNA-seq. The main file(s) that run the scRNA-seq pipelines are: "Mouse Brain 10k - Auto Annotations.R" and "Mouse Brain 10k - CellChat Integration.R". Before the main scripts are run, please set your working directories according to the below file structure.

I have organized my directories as such:

./Data
    
    |-- 10k_Adult_Mouse_Brain
        
        |-- sample_filtered_featured_bc_matrix
        
            |-- barcodes.tsv.gz
         
            |-- features.tsv.gz
         
            |-- matrix.mtx.gz

./Scripts
    
    |-- Mouse Brain 10k - Auto Annotations.R
    
    |-- Mouse Brain 10k - CellChat Integration.R
    
    |-- Mouse Brain 10k Functions.R
    
    |-- Setup.R
