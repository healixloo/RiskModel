Establishment and evaluation of module-based immune-associated gene signature to predict overall survival in patients of colon adenocarcinoma

ReadMe
Jing Lu, 10.10.2022

Usage of R scripts in the folder original data:
Put the R script in the folder that contained the input files, and run as "Rscript script.R", then outputs will be generated.

====================================================================================

Packages and their versions in R used in this project are listed as below:

R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.15.7

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils    
[7] datasets  methods   base     

other attached packages:
 [1] dplyr_0.8.0.1           GSEABase_1.44.0        
 [3] graph_1.60.0            annotate_1.58.0        
 [5] XML_3.98-1.19           AnnotationDbi_1.44.0   
 [7] IRanges_2.16.0          S4Vectors_0.20.1       
 [9] Biobase_2.40.0          BiocGenerics_0.28.0    
[11] GSVA_1.30.0             survminer_0.4.9        
[13] nomogramFormula_1.2.0.0 rms_5.1-4              
[15] SparseM_1.77            Hmisc_4.3-1            
[17] Formula_1.2-3           lattice_0.20-38        
[19] survival_3.2-10         beeswarm_0.2.3         
[21] pheatmap_1.0.12         survivalROC_1.0.3      
[23] estimate_1.0.13         limma_3.36.5           
[25] circlize_0.4.10         ggpubr_0.2.3           
[27] magrittr_1.5            ggplot2_3.3.3          

loaded via a namespace (and not attached):
 [1] nlme_3.1-139        bitops_1.0-6        bit64_0.9-7        
 [4] RColorBrewer_1.1-2  tools_3.5.1         backports_1.1.4    
 [7] R6_2.4.0            rpart_4.1-15        DBI_1.0.0          
[10] colorspace_1.4-1    nnet_7.3-12         withr_2.1.2        
[13] tidyselect_0.2.5    gridExtra_2.3       bit_1.1-14         
[16] compiler_3.5.1      quantreg_5.51       htmlTable_1.13.1   
[19] sandwich_3.0-0      scales_1.0.0        checkmate_1.9.1    
[22] survMisc_0.5.5      polspline_1.1.17    mvtnorm_1.1-0      
[25] stringr_1.4.0       digest_0.6.18       foreign_0.8-71     
[28] base64enc_0.1-3     pkgconfig_2.0.2     htmltools_0.4.0    
[31] fastmap_1.0.1       htmlwidgets_1.5.1   rlang_0.4.5        
[34] GlobalOptions_0.1.2 RSQLite_2.1.1       rstudioapi_0.10    
[37] shiny_1.4.0         shape_1.4.5         generics_0.0.2     
[40] zoo_1.8-6           acepack_1.4.1       RCurl_1.95-4.12    
[43] Matrix_1.2-17       Rcpp_1.0.7.2        munsell_0.5.0      
[46] stringi_1.4.3       multcomp_1.4-16     yaml_2.2.0         
[49] MASS_7.3-51.4       blob_1.1.1          grid_3.5.1         
[52] promises_1.1.0      crayon_1.3.4        splines_3.5.1      
[55] knitr_1.22          pillar_1.3.1        ggsignif_0.6.0     
[58] geneplotter_1.58.0  codetools_0.2-16    glue_1.3.1         
[61] latticeExtra_0.6-28 data.table_1.12.2   httpuv_1.5.2       
[64] MatrixModels_0.4-1  gtable_0.3.0        purrr_0.3.2        
[67] tidyr_0.8.3         km.ci_0.5-2         assertthat_0.2.1   
[70] xfun_0.6            mime_0.6            xtable_1.8-3       
[73] broom_0.5.2         later_1.0.0         tibble_2.1.1       
[76] shinythemes_1.1.2   memoise_1.1.0       KMsurv_0.1-5       
[79] cluster_2.0.8       TH.data_1.0-10     
