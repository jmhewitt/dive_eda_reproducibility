## Quick Start Guide

To get up and rolling with the code just to see what happens, you need
to do these 5 things:

1.  clone the [dive_eda repo from
    GitHub](https://github.com/jmhewitt/dive_eda)
2.  make sure you have all the packages installed (see below for list)
3.  In RStudio, type `targets::tar_make()` into the R console and run
    the entire analysis pipeline (it is structured to run using the
    [targets package in R](https://books.ropensci.org/targets/))
4.  Load the output with `targets::tar_load(eda_results)`. This will
    load a 12-element list that contains the summary output for each tag
    (note that you can view the summary manifest with
    `targets::tar_load(eda_results_manifest)`.
5.  View the results as a summary table of the significant results,
    e.g. `eda_results[[7]]$df` where the 7 is the index in the
    `eda_results_manifest` list

These packages are required to run the analysis pipeline:

-   `targets`
-   `future`
-   `future.batchtools`
-   `dplyr`
-   `lubridate`
-   `ggplot2`
-   `ggthemes`
-   `stringr`
-   `tidyr`

A note on run time - as of this writing, the pipeline takes about 76
minutes to run over the 8 example tags, each of which contains
approximately 14 days’ worth of diving data at 5-minute intervals.

Processing was performed and timed using a 2019 MacBook Pro with a
2.3GHz 8-Core Intel Core i9 with 64GB of RAM. This can be sped up with
use of multiple cores; e.g. `targets::tar_make_future(workers = 4)` see
`targets::tar_make_future()` for guidance. See the end of the document
for full session info, including the R version, and R packages version.

Assuming you are in the R project, and have the `targets` library
loaded, in Rcode, you would source this code to run the analysis:

``` r
targets::tar_make()
targets::tar_load(eda_results) 
eda_results
```

And that summary table of the conditional tests for an example tag
(e.g., `ZcTag095_DUML`) looks as follows:

| stat                 | 0–6h | 0–3h | 3–6h     | 0–1h     | 1–2h | 2–3h | 3–4h | 4–5h | 5–6h     |
|:---------------------|:-----|:-----|:---------|:---------|:-----|:-----|:-----|:-----|:---------|
| prop_surface         |      |      |          |          |      |      |      |      |          |
| iddi_obs             |      |      |          | 0.09 (-) |      |      |      |      |          |
| prop_upward          |      |      |          | 0.09 (-) |      |      |      |      | 0.08 (-) |
| prop_downward        |      |      |          |          |      |      |      |      |          |
| prop_no_change       |      |      |          | 0.05 (+) |      |      |      |      | 0.08 (+) |
| total_absolute       |      |      |          |          |      |      |      |      |          |
| total_upward         |      |      |          |          |      |      |      |      |          |
| total_downward       |      |      |          |          |      |      |      |      |          |
| depth_seq            |      |      | 0.05 (+) |          |      |      |      |      | 0.02 (+) |
| depth_diff_seq       |      |      | 0.06 (+) |          |      |      |      |      |          |
| type_seq             |      |      |          |          |      |      |      |      | 0.1 (+)  |
| direction_change_seq |      |      |          | 0.07 (-) |      |      |      |      |          |

Table 1. Results from one example whale (ZcTag095_DUML, conditional)
indicates sigificance for each test-statistic across the 9 temporal
windows. Results are significant at the 0.1 level if printed in table.
The conditional tests represent pre-post pairs, where the ‘pre’ window
must match the context of the diving state of the animal at the time of
the exposure. See Section 2.2 of the manuscript for details.

The `eda_results_manifest` allows you to glimpse the different animals
under the different conditions:

-   To get the results for unconditional tests for `ZcTag069`, type
    `eda_results[[2]]$df`
-   To get conditional results for `ZcTag097_DUML`, type
    `eda_results[[11]]$df`
-   The full manifest for the example tags is:

``` r
targets::tar_load(eda_results_manifest)
eda_results_manifest
```

    ##         tag conditional index
    ## 1  ZcTag069        TRUE     1
    ## 2  ZcTag069       FALSE     2
    ## 3  ZcTag085        TRUE     3
    ## 4  ZcTag085       FALSE     4
    ## 5  ZcTag093        TRUE     5
    ## 6  ZcTag093       FALSE     6
    ## 7  ZcTag095        TRUE     7
    ## 8  ZcTag095       FALSE     8
    ## 9  ZcTag096        TRUE     9
    ## 10 ZcTag096       FALSE    10
    ## 11 ZcTag097        TRUE    11
    ## 12 ZcTag097       FALSE    12

## Data

The code base is set up to run using the
[targets](https://books.ropensci.org/targets/) package. We describe the
structure of the targets-based workflow, as well as the analytical R
code. We assume you have cloned the [dive_eda
repo](https://github.com/jmhewitt/dive_eda) and opened up the project in
R. We wrote this code to work with series data that arise from a
[SPLASH-10 tag from WildLife
Computers](https://wildlifecomputers.com/our-tags/splash-archiving-tags/splash10/).
In our case, the tags were programmed to record the depth every 5
minutes for approximately 2 weeks:

``` r
dive_dat <- read_csv(here::here('data/sattag', 'ZcTag084_DUML_series_20200108.csv'))

dive_dat %>% 
  select(DeployID, Day, Time, Depth, DRange)
```

    ## # A tibble: 3,878 × 5
    ##    DeployID Day         Time   Depth DRange
    ##    <chr>    <chr>       <time> <dbl>  <dbl>
    ##  1 ZcTag084 23-May-2019 12:50    267   27.5
    ##  2 ZcTag084 23-May-2019 12:55    267   27.5
    ##  3 ZcTag084 23-May-2019 13:00    170   26.5
    ##  4 ZcTag084 23-May-2019 13:05     12   12.5
    ##  5 ZcTag084 23-May-2019 13:10    534   54.5
    ##  6 ZcTag084 23-May-2019 13:15   1165  110. 
    ##  7 ZcTag084 23-May-2019 13:20   1360  112. 
    ##  8 ZcTag084 23-May-2019 13:25   1360  112. 
    ##  9 ZcTag084 23-May-2019 13:30   1360  112. 
    ## 10 ZcTag084 23-May-2019 13:35   1360  112. 
    ## # … with 3,868 more rows

Thus, we are concerned with the time series of values in the `Depth`
column, which as described in the User Manual for the tags, is the
center of a depth bin (in meters); `DRange` is the +/- of the error
(also in meters). In the manuscript, we run the code over a series of
animals, e.g. 8 tags; the code is structured to run over all of these
tags.

## References

Quick, N. J., Cioffi, W. R., Shearer, J., & Read, A. J. (2019). Mind the
gap—Optimizing satellite tag settings for time series analysis of
foraging dives in Cuvier’s beaked whales (*Ziphius cavirostris*). Animal
Biotelemetry, 7(1), 5. <https://doi.org/10.1186/s40317-019-0167-5>

Shearer, J. M., Quick, N. J., Cioffi, W. R., Baird, R. W., Webster, D.
L., Foley, H. J., Swaim, Z. T., Waples, D. M., Bell, J. T., & Read, A.
J. (2019). Diving behaviour of Cuvier’s beaked whales (*Ziphius
cavirostris*) off Cape Hatteras, North Carolina. R Soc Open Sci, 6(2),
181728. <https://doi.org/10.1098/rsos.181728>

## R Environment

The current R environment on which this was run is here:

``` r
sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] targets_0.6.0   forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7    
    ##  [5] purrr_0.3.4     readr_2.0.2     tidyr_1.1.3     tibble_3.1.5   
    ##  [9] ggplot2_3.3.5   tidyverse_1.3.1 here_1.0.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7        lubridate_1.7.10  ps_1.6.0          assertthat_0.2.1 
    ##  [5] rprojroot_2.0.2   digest_0.6.28     utf8_1.2.2        R6_2.5.1         
    ##  [9] cellranger_1.1.0  backports_1.2.1   reprex_2.0.1      evaluate_0.14    
    ## [13] httr_1.4.2        highr_0.9         pillar_1.6.3      rlang_0.4.11     
    ## [17] readxl_1.3.1      rstudioapi_0.13   data.table_1.14.0 callr_3.7.0      
    ## [21] rmarkdown_2.9     bit_4.0.4         igraph_1.2.6      munsell_0.5.0    
    ## [25] broom_0.7.9       compiler_4.1.0    modelr_0.1.8      xfun_0.26        
    ## [29] pkgconfig_2.0.3   htmltools_0.5.1.1 tidyselect_1.1.1  codetools_0.2-18 
    ## [33] fansi_0.5.0       crayon_1.4.1      tzdb_0.1.2        dbplyr_2.1.1     
    ## [37] withr_2.4.2       grid_4.1.0        jsonlite_1.7.2    gtable_0.3.0     
    ## [41] lifecycle_1.0.1   DBI_1.1.1         magrittr_2.0.1    scales_1.1.1     
    ## [45] vroom_1.5.5       cli_3.0.1         stringi_1.7.5     fs_1.5.0         
    ## [49] xml2_1.3.2        ellipsis_0.3.2    generics_0.1.0    vctrs_0.3.8      
    ## [53] tools_4.1.0       bit64_4.0.5       glue_1.4.2        hms_1.1.1        
    ## [57] parallel_4.1.0    processx_3.5.2    yaml_2.2.1        colorspace_2.0-2 
    ## [61] rvest_1.0.1       knitr_1.36        haven_2.4.3
