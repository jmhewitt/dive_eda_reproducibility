summary_targets = list(

  # collapse raw results for printing
  tar_target(
    name = eda_results, 
    command = {
      
      # store FDR pvalue thresholds
      cutoffs = list()
      
      # standard significance value cutoff
      cutoff_std = .1
      
      # FDR control level
      qstar = .05
      
      # extract all difference-based results into a single data frame
      df = do.call(rbind, lapply(sattag_eda, function(x) {
        cbind(
          data.frame(tag = x$tag$deployid),
          x$window,
          data.frame(stat = colnames(x$null_samples),
                     p = pmin(x$p.left, x$p.right)) %>% 
            mutate(right_tail = p == x$p.right,
                   pval = punif(q = p, min = 0, max = .5)),
          conditional_results = x$conditional
        )
      })) %>% 
        filter(conditional_results == conditional) %>% 
        select(-conditional_results)
      
      # FDR significant value cutoffs
      cutoffs$df = df %>% 
        group_by(tag) %>% 
        summarise(
          cutoff = {
            o = order(pval)
            signif = pval[o] <= seq_along(pval) / length(pval) * qstar
            ifelse(sum(signif) > 0, max(pval[o][signif]), 0)
          }
        )
      
      # simplify display of p-values
      df$p_fmt = apply(
        df %>% left_join(cutoffs$df, by = 'tag') %>% 
          select(p, right_tail, cutoff), 1, 
        function(r) {
          if(r['p'] <= cutoff_std) {
            paste(
              sprintf(fmt = '%.2f', round(r['p'],2)), 
              ifelse(r['right_tail'], '(+)', '(-)'))
          } else {
            ''
          }
        })
      
      # simplify display of FDR-corrected p-values
      df$p_fmt_fdr = apply(
        df %>% left_join(cutoffs$df, by = 'tag') %>% 
          select(pval, right_tail, cutoff), 1, 
        function(r) {
        if(r['pval'] <= r['cutoff']) {
          ifelse(r['right_tail'], '(+)', '(-)')
        } else {
          ''
        }
      })
      
      # filter results for a single tag and wide-format the table
      df_print = df %>% 
        filter(tag == tag_id) %>%
        select(description, stat, p_fmt) %>% 
        pivot_wider(names_from = description, values_from = p_fmt)
      
      # filter results for a single tag and wide-format the table
      df_print_fdr = df %>% 
        filter(tag == tag_id) %>%
        select(description, stat, p_fmt_fdr) %>% 
        pivot_wider(names_from = description, values_from = p_fmt_fdr)
      
      if(nrow(df_print) == 0) {
        return(NULL)
      }
      
      #
      # format results for bivariate tests
      #
      
      # extract all difference-based results into a single data frame
      df.bivariate = do.call(rbind, lapply(sattag_eda, function(x) {
        cbind(
          data.frame(tag = x$tag$deployid),
          x$window,
          x$p.bivariate %>% 
            filter(type == 'bivariate_kde') %>% 
            mutate(p = pmin(p.left, p.right),
                   pval = punif(q = p, min = 0, max = 0.5),
                   right_tail = p == p.right),
          conditional_results = x$conditional
        ) %>% select(tag, description, response_lag, window_length, stat, p, 
                     pval, right_tail, conditional_results)
      })) %>% 
        filter(conditional_results == conditional) %>% 
        select(-conditional_results)
      
      # FDR significant value cutoffs
      cutoffs$df.bivariate = df.bivariate %>% 
        group_by(tag) %>% 
        summarise(
          cutoff = {
            o = order(pval)
            signif = pval[o] <= seq_along(pval) / length(pval) * qstar
            ifelse(sum(signif) > 0, max(pval[o][signif]), 0)
          }
        )
      
      # simplify display of FDR-corrected p-values
      df.bivariate$p_fmt_fdr = apply(
        df.bivariate %>% left_join(cutoffs$df.bivariate, by = 'tag') %>% 
          select(pval, right_tail, cutoff), 1, 
        function(r) {
          if(r['pval'] <= r['cutoff']) {
            ifelse(r['right_tail'], '(+)', '(-)')
          } else {
            ''
          }
        })
      
      # simplify display of p-values
      df.bivariate$p_fmt = apply(
        df.bivariate %>% left_join(cutoffs$df.bivariate, by = 'tag') %>% 
          select(p, right_tail, cutoff), 1, 
        function(r) {
          if(r['p'] <= cutoff_std) {
            paste(
              sprintf(fmt = '%.2f', round(r['p'],2)), 
              ifelse(r['right_tail'], '(+)', '(-)'))
          } else {
            ''
          }
        })
      
      # filter results for a single tag and wide-format the table
      df.bivariate_print = df.bivariate %>% 
        filter(tag == tag_id) %>%
        select(description, stat, p_fmt) %>% 
        pivot_wider(names_from = description, values_from = p_fmt)
      
      # filter results for a single tag and wide-format the table
      df.bivariate_print_fdr = df.bivariate %>% 
        filter(tag == tag_id) %>%
        select(description, stat, p_fmt_fdr) %>% 
        pivot_wider(names_from = description, values_from = p_fmt_fdr)
      
      #
      # format results for mahalanobis tests
      #
      
      # extract all difference-based results into a single data frame
      df.mahalanobis = do.call(rbind, lapply(sattag_eda, function(x) {
        cbind(
          data.frame(tag = x$tag$deployid),
          x$window,
          x$p.bivariate %>% 
            filter(type == 'bivariate_mdist') %>% 
            mutate(p = p.right,
                   pval = p,
                   right_tail = TRUE),
          conditional_results = x$conditional
        ) %>% select(tag, description, response_lag, window_length, stat, p, 
                     pval, right_tail, conditional_results)
      })) %>% 
        filter(conditional_results == conditional) %>% 
        select(-conditional_results)
      
      # FDR significant value cutoffs
      cutoffs$df.mahalanobis = df.mahalanobis %>% 
        group_by(tag) %>% 
        summarise(
          cutoff = {
            o = order(pval)
            signif = pval[o] <= seq_along(pval) / length(pval) * qstar
            ifelse(sum(signif) > 0, max(pval[o][signif]), 0)
          }
        )
      
      # simplify display of FDR-corrected p-values
      df.mahalanobis$p_fmt_fdr = apply(
        df.mahalanobis %>% left_join(cutoffs$df.mahalanobis, by = 'tag') %>% 
          select(pval, right_tail, cutoff), 1, 
        function(r) {
          if(r['pval'] <= r['cutoff']) {
            ifelse(r['right_tail'], '(+)', '(-)')
          } else {
            ''
          }
        })
      
      # simplify display of p-values
      df.mahalanobis$p_fmt = apply(
        df.mahalanobis %>% left_join(cutoffs$df.mahalanobis, by = 'tag') %>% 
          select(p, right_tail, cutoff), 1, 
        function(r) {
          if(r['p'] <= cutoff_std) {
            paste(
              sprintf(fmt = '%.2f', round(r['p'],2)), 
              ifelse(r['right_tail'], '(+)', '(-)'))
          } else {
            ''
          }
        })
      
      # filter results for a single tag and wide-format the table
      df.mahalanobis_print = df.mahalanobis %>% 
        filter(tag == tag_id) %>%
        select(description, stat, p_fmt) %>% 
        pivot_wider(names_from = description, values_from = p_fmt)
      
      # filter results for a single tag and wide-format the table
      df.mahalanobis_print_fdr = df.mahalanobis %>% 
        filter(tag == tag_id) %>%
        select(description, stat, p_fmt_fdr) %>% 
        pivot_wider(names_from = description, values_from = p_fmt_fdr)
      
      
      # 
      # display the formatted table with the desired row order
      #
      
      # row_order = c('prop_surface', 'iddi_obs', 'prop_upward', 'prop_downward', 
      #               'prop_no_change', 'total_absolute', 'total_upward', 
      #               'total_downward', 'depth_seq', 'depth_diff_seq', 'type_seq',
      #               'direction_change_seq')
      
      row_order = c('time_on_surface', 'avg_depth', 'time_deep_diving',
                    'time_ascending', 'time_descending', 'time_level',
                    'total_vertical', 'total_ascending', 'total_descending',
                    'direction_changes', 'timing_time_on_surface', 
                    'timing_avg_depth', 'timing_time_deep_diving', 
                    'timing_time_ascending', 
                    'timing_time_descending', 'timing_time_level',
                    'timing_total_vertical', 'timing_total_ascending',
                    'timing_total_descending', 'timing_direction_changes')
      
      list(list(
        df = df_print[
          sapply(row_order, function(x) { which(x == df_print$stat) }),
        ],
        df_fdr = df_print_fdr[
          sapply(row_order, function(x) { which(x == df_print_fdr$stat) }),
        ], 
        df.bivariate = df.bivariate_print[
          sapply(row_order, function(x) { which(x == df.bivariate_print$stat) }),
        ],
        df.bivariate_fdr = df.bivariate_print_fdr[
          sapply(row_order, function(x) { which(x == df.bivariate_print_fdr$stat) }),
        ],
        df.mahalanobis = df.mahalanobis_print[
          sapply(row_order, function(x) { which(x == df.mahalanobis_print$stat) }),
        ],
        df.mahalanobis_fdr = df.mahalanobis_print_fdr[
          sapply(row_order, function(x) { which(x == df.mahalanobis_print_fdr$stat) }),
        ],
        df.bivariate.raw = df.bivariate %>% filter(tag == tag_id),
        conditional = conditional, 
        tag_id = tag_id,
        cutoffs = cutoffs
      ))
      
      #
      # aggregate number of baseline samples used to train results
      #
      
      # df_n_baseline = do.call(rbind, lapply(sattag_eda, function(x) {
      #   cbind(
      #     data.frame(tag = x$tag$deployid,
      #                n = nrow(x$null_samples)),
      #     x$window
      #   )
      # }))
      # 
      # # filter results for a single tag and wide-format the table
      # df_n_print = df_n_baseline %>% 
      #   filter(tag == 'ZcTag069') %>%  
      #   select(-response_lag, -window_length, -tag) %>% 
      #   pivot_wider(names_from = description, values_from = n)
      # 
      # df_n_baseline %>% 
      #   filter(window_length == 6)   
      
    },
    pattern = cross(tag_id, map(conditional)), 
    deployment = 'worker',
    memory = 'transient',
    storage = 'worker'
  ),
  
  tar_target(
    name = eda_results_manifest,
    command = {
      do.call(rbind, lapply(1:length(eda_results), function(ind) { 
        x = eda_results[[ind]]
        data.frame(tag = x$tag_id, conditional = x$conditional, index = ind)
      }))
    }
  ),
  
  tar_target(
    name = population_summary_table,
    command = {
      
      aggregate_results = do.call(rbind, lapply(eda_results, function(x) {
        if(x$conditional) {
          x$df.bivariate.raw
        } else {
          NULL
        }
      }))
      
      aggregate_results %>% 
        group_by(tag) %>% 
        summarise(signif_level = .1,
                  total_signif = sum(p < signif_level),
                  prop_signif = total_signif / length(p),
                  total_signif_chance = n() * signif_level,
                  excess_signif = total_signif / total_signif_chance)
    }
  )
)
