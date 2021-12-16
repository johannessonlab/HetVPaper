# The fertility effects of Pa_5_12720 point mutations

In here you will find a very simply R script and the data required for Fig. S11b. Unlike other figures in the pipelines, it uses two additional R packages: `rstatix` v. 0.7.0 and `ggpubr` 0.4.0.

Here is a session of the last environment where I ran the script and it worked well:

    > sessionInfo()
    R version 4.1.1 (2021-08-10)
    Platform: aarch64-apple-darwin20 (64-bit)
    Running under: macOS Big Sur 11.6

    Matrix products: default
    LAPACK: /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRlapack.dylib

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] ggpubr_0.4.0  rstatix_0.7.0 cowplot_1.1.1 ggplot2_3.3.5 dplyr_1.0.7  

    loaded via a namespace (and not attached):
     [1] tidyselect_1.1.1   purrr_0.3.4        haven_2.4.3        carData_3.0-4      colorspace_2.0-2   vctrs_0.3.8        generics_0.1.0    
     [8] utf8_1.2.2         rlang_0.4.12       pillar_1.6.4       foreign_0.8-81     glue_1.5.1         withr_2.4.3        DBI_1.1.1         
    [15] RColorBrewer_1.1-2 readxl_1.3.1       lifecycle_1.0.1    munsell_0.5.0      ggsignif_0.6.3     gtable_0.3.0       cellranger_1.1.0  
    [22] zip_2.2.0          labeling_0.4.2     rio_0.5.27         forcats_0.5.1      curl_4.3.2         fansi_0.5.0        broom_0.7.10      
    [29] Rcpp_1.0.7         backports_1.2.1    scales_1.1.1       abind_1.4-5        farver_2.1.0       digest_0.6.29      hms_1.1.1         
    [36] stringi_1.7.4      openxlsx_4.2.4     grid_4.1.1         cli_3.1.0          tools_4.1.1        magrittr_2.0.1     tibble_3.1.6      
    [43] crayon_1.4.2       car_3.0-11         tidyr_1.1.3        pkgconfig_2.0.3    ellipsis_0.3.2     data.table_1.14.0  assertthat_0.2.1  
    [50] rstudioapi_0.13    R6_2.5.1           compiler_4.1.1  
