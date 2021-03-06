library(dplyr)
library(xtable)
library(stringr)

targets::tar_load(eda_results)
targets::tar_load(eda_results_manifest)

# 88, 89, 93, 95, 96, 97

fdr_table = TRUE

tag_id = 'ZcTag093'

ind = eda_results_manifest %>% 
  filter(tag == tag_id, conditional == TRUE) %>% 
  select(index) %>% 
  as.numeric()

# select table to format
if(fdr_table) {
  df = eda_results[[ind]]$df.bivariate_fdr
} else {
  df = eda_results[[ind]]$df.bivariate
}

# print (+) and (-) using latex math formatting
df = apply(df, 2, function(col) {
  txt = gsub(pattern = '\\(', replacement = '$(', col)
  txt = gsub(pattern = '\\)', replacement = ')$', txt)
  txt = gsub(pattern = '_', replacement = ' ', txt)
})

# update principal column name
colnames(df)[1] = '$t_k(w)$'

# pair the sequential and non-sequential test results
df = df[c(rbind(1:(nrow(df)/2),(nrow(df)/2+1):nrow(df))), ]

# name substitution
sub_names = c(
  "time on surface" = 'Time on surface', 
  "avg depth" = 'Average depth', 
  "time deep diving" = 'Time in deep dives', 
  "time ascending" = 'Time spent ascending', 
  "time descending" = 'Time spent descending', 
  "time level" = 'Time spent without vertical movement', 
  "total vertical" = 'Total vertical distance traveled', 
  "total ascending" = 'Total ascent distance', 
  "total descending" = 'Total descent distance', 
  "direction changes" = 'Total vertical direction changes'
)

# update row names wrt. pairing
df[seq(from = 2, to = nrow(df), by = 2), 1] = ''
df[seq(from = 1, to = nrow(df), by = 2), 1] = paste(
  '\\multirow{2}{*}{',
  paste(1:(nrow(df)/2), '.', sep =''),
  sub_names[df[seq(from = 1, to = nrow(df), by = 2), 1]],
  '}'
)

# add cell color coding and spacing
df[seq(from = 2, to = nrow(df), by = 2),-1] = apply(
  df[seq(from = 2, to = nrow(df), by = 2),-1], 2, function(col) {
  paste('\\emph{\\color{gray}', col, '}')
})

cutoff_fmt = eda_results[[ind]]$cutoffs$df.bivariate %>% 
  filter(tag == tag_id) %>% 
  select(cutoff) %>% 
  as.numeric() %>% 
  format(scientific = TRUE, digits = 2) %>% 
  toupper()

if(fdr_table) {
  caption = paste(
    eda_results[[ind]]$tag_id, ', tail probabilities conditional on dive ',
    'state at exposure without order sensitivity $t_k(w)$ (black), and ',
    'with order sensitivity $t_k^D(w)$ (grey, emphasized). Probabilities ', 
    'are included if smaller than the tag\'s FDR threshold of ', cutoff_fmt, 
    ', and are annotated with $(-)$ or ',
    '$(+)$ to indicate they are for the left or right tail, respectively.',
    '\\linebreak',
    sep = ''
  )
} else {
  caption = paste(
    eda_results[[ind]]$tag_id, ', tail probabilities conditional on dive ',
    'state at exposure without order sensitivity $t_k(w)$ (black), and ',
    'with order sensitivity $t_k^D(w)$ (grey, emphasized). Probabilities ', 
    'are excluded if greater than 0.1, and are annotated with $(-)$ or',
    '$(+)$ to indicate they are for the left or right tail, respectively.',
    '\\linebreak',
    sep = ''
  )
}



x = print.xtable(
  xtable(
    x = df, 
    align = 'rl|l|ll|llllll', 
    caption = caption
  ), 
  booktabs = TRUE, 
  sanitize.text.function = identity, 
  include.rownames = FALSE,
  hline.after = unique(c(-1,0,seq(from = 2, to = nrow(df), by = 2), nrow(df))), 
  caption.placement = 'top', 
  comment = FALSE
)

x = stringr::str_replace(
  string = x, 
  pattern = fixed('\\bottomrule'),
  replacement = '\\arrayrulecolor{black}\\bottomrule'
)
x = stringr::str_replace(
  string = x, 
  pattern = fixed('\\end{tabular}'),
  replacement = '\\end{tabular}\n}'
)
x = stringr::str_replace(
  string = x, 
  pattern = fixed('\\begin{tabular}'),
  replacement = '\\resizebox{\\textwidth}{!}{\n\\begin{tabular}'
)
x = stringr::str_replace(
  string = x, 
  pattern = fixed('\\multirow'),
  replacement = '\\arrayrulecolor{black!12}\\multirow'
)

clipr::write_clip(x)
