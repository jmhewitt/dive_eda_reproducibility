sim_plot_targets = list(
  
  # significance levels to make simulation results plots for
  tar_target(
    name = alpha,
    command = c(.01, .05, .1)
  ),
  
  # simulation results plots
  tar_target(
    name = sim_summary,
    command = {
      
      # load and merge simulation output
      f = dir(path = file.path('output', 'sim', 'samples'), full.names = TRUE)
      samples = do.call(rbind, lapply(f, readRDS))
      
      # enrich with final pvalues and rejections
      samples = samples %>% 
        mutate(
          p = pmin(p.left, p.right),
          pval = punif(q = p, min = 0, max = .5),
          reject = pval < alpha
        )
      
      # compute statistical power function
      res = samples %>% 
        group_by(sim, window_len, lag_len, stat) %>% 
        summarise(
          preject = mean(reject)
        )
      
      # refine plot labels: name order and substitution
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
      timing_names = paste('(Timing of)', sub_names)
      names(timing_names) = paste('timing_', names(sub_names), sep='') 
      levels(res$stat) = c(sub_names, timing_names)[levels(res$stat)]
      res$stat = forcats::fct_rev(res$stat)
      
      res$window_len = paste(res$window_len, 'h window', sep = '')
      res$lag_len = paste(res$lag_len, 'h lag', sep = '')
      
      # # view a good-looking density estimator
      # x = samples %>% 
      #   filter(
      #     window_len == 6,
      #     lag_len == 0,
      #     stat == 'time_deep_diving',
      #     sim == 'baseline'
      #   ) %>% 
      #   select(pval) %>% 
      #   unlist()
      # plot(density(x))
      
      # # view a poor-looking density estimator
      # x = samples %>% 
      #   filter(
      #     window_len == 1,
      #     lag_len == 0,
      #     stat == 'time_level',
      #     sim == 'baseline'
      #   ) %>% 
      #   select(pval) %>% 
      #   unlist()
      # plot(density(x))
      
      # plot rejection rates for simulation configurations
      pl = res %>% 
        mutate(
          sim = ifelse(sim == 'baseline', 'Yes', 'No')
        ) %>% 
        ggplot(aes(y = stat, x = preject, col = sim)) + 
        geom_vline(xintercept = alpha, col = 'grey60') +
        geom_point() + 
        scale_color_brewer(expression(H[0]~' true (No response)'), 
                           type = 'qual', palette = 'Dark2') + 
        facet_grid(window_len ~ lag_len) + 
        xlab(expression(P(Reject~H[0]))) +
        ylab('Summary statistic') +
        theme_few() +
        theme(
          axis.title.y = element_blank(),
          panel.grid.major.y = element_line(colour = 'grey90'),
          panel.grid.major.x = element_line(colour = 'grey95'),
          strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
      
      #
      # save figures
      #
      
      f = file.path('output', 'figures', 'simulation')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      sc = 1.5
      
      ggsave(pl, filename = file.path(f, paste(tar_name(), '_', alpha, '.pdf',
                                               sep = '')), 
             width = sc*11, height = sc*8)
    }, 
    pattern = map(alpha)
  )
  
)
