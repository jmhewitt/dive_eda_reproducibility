library(dplyr)
library(xtable)
library(stringr)
library(tidyr)

targets::tar_load(eda_results)
targets::tar_load(eda_results_manifest)

# work toward total number of small tail probs. per metric
aggregate_results = do.call(rbind, lapply(eda_results, function(x) {
  if(x$conditional) {
    x$df.bivariate.raw
  } else {
    NULL
  }
}))

# restrict results to subset of 2019 animals
aggregate_results = aggregate_results %>% 
  filter(tag %in% c('ZcTag083', 'ZcTag085', 'ZcTag087', 'ZcTag088',
                    'ZcTag089', 'ZcTag093', 'ZcTag095', 'ZcTag096',
                    'ZcTag097'))

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

aggregate_results = aggregate_results %>% 
  mutate(
    description = factor(
      x = gsub('--', '-', description),
      levels = c('0-6h', '0-3h', '3-6h', '0-1h', '1-2h', '2-3h', 
                 '3-4h', '4-5h', '5-6h')
    ),
    stat = factor(
      x = stat,
      levels = rev(c('time_on_surface', 'avg_depth', 'time_deep_diving',
                     'time_ascending', 'time_descending', 'time_level',
                     'total_vertical', 'total_ascending', 'total_descending',
                     'direction_changes', 'timing_time_on_surface', 
                     'timing_avg_depth', 'timing_time_deep_diving', 
                     'timing_time_ascending', 
                     'timing_time_descending', 'timing_time_level',
                     'timing_total_vertical', 'timing_total_ascending',
                     'timing_total_descending', 'timing_direction_changes'))
    ),
    sequential = grepl('timing', stat)
  )

aggregate_results = aggregate_results %>% 
  group_by(stat, description, sequential) %>% 
  summarise(num_signif_animals = sum(p < .1))


df = aggregate_results %>% 
  select(stat, description, num_signif_animals) %>% 
  pivot_wider(
    id_cols = stat, 
    names_from = description, 
    values_from = num_signif_animals
  ) %>% 
  mutate(
    stat = as.character(stat),
    across(everything(), as.character)
  )
  
# name order and substitution
sub_names = c(
  "time_on_surface" = 'Time on surface', 
  "avg_depth" = 'Average depth', 
  "time_deep_diving" = 'Time in deep dives', 
  "time_ascending" = 'Time spent ascending', 
  "time_descending" = 'Time spent descending', 
  "time_level" = 'Time spent without vertical movement', 
  "total_vertical" = 'Total vertical distance traveled', 
  "total_ascending" = 'Total ascent distance', 
  "total_descending" = 'Total descent distance', 
  "direction_changes" = 'Total vertical direction changes'
)
timing_names = sub_names
names(timing_names) = paste('timing_', names(timing_names), sep='') 
sub_names_tgts = as.vector(rbind(names(sub_names), names(timing_names)))
sub_names = as.vector(rbind(sub_names, timing_names))
names(sub_names) = sub_names_tgts

# re-order and pair the sequential and non-sequential test results
df = df[
  sapply(names(sub_names), function(nom) {
    which(unname(nom) == df$stat)
  }),
]

# update principal column name
colnames(df)[1] = '$t_k(w)$'

# remove 0's from display
df[df==0] = NA

# update row names wrt. pairing
df[seq(from = 2, to = nrow(df), by = 2), 1] = ''
df[seq(from = 1, to = nrow(df), by = 2), 1] = paste(
  '\\multirow{2}{*}{',
  paste(1:(nrow(df)/2), '.', sep =''),
  sub_names[df[seq(from = 1, to = nrow(df), by = 2), 1] %>% unlist()],
  '}'
)

# add cell color coding and spacing
df[seq(from = 2, to = nrow(df), by = 2),-1] = apply(
  df[seq(from = 2, to = nrow(df), by = 2),-1], 2, function(col) {
      ifelse(
        is.na(col), 
        rep(NA,length(col)), 
        paste('\\emph{\\color{gray}', col, '}')
      )
  })

x = print.xtable(
  xtable(
    x = df, 
    align = 'rl|l|ll|llllll', 
    caption = paste(
      'Number of animals with either left or right tail probabilities less ',
      'than 0.1, conditional on dive state at exposure without order ',
      'sensitivity $t_k(w)$ (black), and with order sensitivity $t_k^D(w)$ ', 
      '(grey, emphasized).\\linebreak',
      sep = ''
    )
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
