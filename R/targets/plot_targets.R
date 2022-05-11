plot_targets = list(

  tar_target(
    name = statistic_demo,
    command = {

      tgt_tag = 'ZcTag093'
      
      
      #
      # identify window to analyze
      #
      
      manifest = do.call(rbind, lapply(pre_post_pair_defs, function(w) {
        data.frame(tag = w$tag$deployid,
                   response_lag = w$window$response_lag,
                   window_length = w$window$window_length,
                   conditional = w$baseline_pairs_are_conditional,
                   response_start = w$response_window[1],
                   response_end = w$response_window[2],
                   pre_exposure_start = w$pre_exposure_window[1],
                   pre_exposure_end = w$pre_exposure_window[2])
      }))
      manifest$ind = 1:nrow(manifest)
      
      window_selected = manifest %>% 
        filter(tag == tgt_tag,
               response_lag == 1,
               window_length == 1,
               conditional == TRUE)
      
      #
      # load tag data
      #
      
      # identify depth record file for tag
      tar_load(tag_info)
      tag_info = tag_info %>% filter(deployid == window_selected$tag)
      f = grep(pattern = tag_info$deployid, x = depth_files, value = TRUE)
      
      # load data
      d = read.csv(f)
      d$Date = as.POSIXct(d$Date, origin = '1970-01-01 00:00.00 UTC', 
                          tz = 'UTC')
      
      # map all depths to standardized bins
      d$depth.bin = sapply(d$Depth, function(depth) {
        which.min(abs(depth - template_bins$center))
      })
      d$depth.standardized = template_bins$center[d$depth.bin]
      
      # add dive segmentation to depths time series
      d$diveId = dive_labels[[
        which(sapply(dive_labels, function(x) x$tag == tag_info$deployid))
      ]]$labels
      
      # classify dives by max observed depth
      diveTypes = d %>%
        dplyr::group_by(diveId) %>%
        dplyr::summarise(
          diveType = ifelse(any(diveId == 0), 'Unknown', ifelse(
            max(depth.standardized) > deep_dive_depth, 'Deep', 'Shallow'
          ))
        )
      
      # merge dive type information with sattag record
      d = d %>% left_join(diveTypes, by = 'diveId')
      
      # reformat column names and add derived data to sattag record
      d = d %>% mutate(
        depths = depth.standardized,
        ddepths = c(0, diff(depths)),
        ddepths.sign = sign(ddepths),
        descending = ddepths.sign > 0,
        ascending = ddepths.sign < 0,
        no_change = ddepths.sign == 0
      )
      
      #
      # plot data, demonstrating the raw statistics
      #
      
      pl = ggplot(d, aes(x = Date, y = depths, group = 1)) + 
        # pre/post windows
        geom_rect(aes(xmin = window_selected$pre_exposure_start,
                  xmax = window_selected$pre_exposure_end, 
                  ymin = 0,
                  ymax = Inf),
                  data = data.frame(x = 1),
                  inherit.aes = FALSE,
                  alpha = .15) +
        geom_rect(aes(xmin = window_selected$response_start,
                      xmax = window_selected$response_end, 
                      ymin = 0,
                      ymax = Inf),
                  data = data.frame(x = 1),
                  inherit.aes = FALSE,
                  alpha = .15) + 
        # horizontal surface line
        geom_hline(yintercept = template_bins$center[1], lty = 3) +
        # depth differences
        geom_segment(mapping = aes(x = Date, xend = Date, 
                                   y = depths, yend = depths_next, 
                                   group = Date), 
                     data = d %>% 
                       mutate(depths_next = lead(depths)) %>% 
                       filter(depths != depths_next), 
                     arrow = arrow(angle = 30, length = unit(.016,'npc'),
                                   type = 'closed'),
                     col = '#de2d26', alpha = .6) + 
        # time series plot of depth bins
        geom_line() +
        geom_point(mapping = aes(alpha = diveType)) + 
        # cee start time
        geom_vline(xintercept = tag_info$cee_start, lty = 2, alpha = .5) + 
        # axis and plot formatting
        scale_x_datetime('Time (UTC)', 
                         date_breaks = 'hour',
                         date_labels = '%H%M',
                         limits = c(window_selected$pre_exposure_start, 
                                    window_selected$response_end) + 
                                    c(-1,1) * 3600) + 
        scale_alpha_manual(values = c(Deep = 1, Shallow = 0, Unknown = 0)) + 
        scale_y_reverse() +
        guides(col = 'none', alpha = 'none') + 
        ylab('Depth (m)') + 
        theme_few()
      
      f = file.path('output','figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      ggsave(pl, filename = file.path(f, paste(tar_name(), '_', tgt_tag, '.pdf', 
                                               sep ='')),
             width = 8, height = 4, dpi = 'print')
    }
  ),
  
  tar_target(
    name = tag_exposure_plot,
    command = {
      
      tag_id = tags_2019
      
      #
      # identify window to analyze
      #
      
      manifest = do.call(rbind, lapply(pre_post_pair_defs, function(w) {
        data.frame(tag = w$tag$deployid,
                   response_lag = w$window$response_lag,
                   window_length = w$window$window_length,
                   conditional = w$baseline_pairs_are_conditional,
                   response_start = w$response_window[1],
                   response_end = w$response_window[2],
                   pre_exposure_start = w$pre_exposure_window[1],
                   pre_exposure_end = w$pre_exposure_window[2])
      }))
      manifest$ind = 1:nrow(manifest)
      
      window_selected = manifest %>% 
        filter(tag == tag_id,
               response_lag == 0,
               window_length == 6,
               conditional == TRUE)
      
      window_selected2 = manifest %>%
        filter(tag == tag_id,
               response_lag %in% c(1,3,5),
               window_length == 1,
               conditional == TRUE)
      
      #
      # load tag data
      #
      
      # identify depth record file for tag
      tar_load(tag_info)
      tag_info = tag_info %>% filter(deployid == window_selected$tag)
      f = grep(pattern = tag_info$deployid, x = depth_files, value = TRUE)
      
      # load data
      d = read.csv(f)
      d$Date = as.POSIXct(d$Date, origin = '1970-01-01 00:00.00 UTC', 
                          tz = 'UTC')
      
      # map all depths to standardized bins
      d$depth.bin = sapply(d$Depth, function(depth) {
        which.min(abs(depth - template_bins$center))
      })
      d$depth.standardized = template_bins$center[d$depth.bin]
      
      # add dive segmentation to depths time series
      d$diveId = dive_labels[[
        which(sapply(dive_labels, function(x) x$tag == tag_info$deployid))
      ]]$labels
      
      # classify dives by max observed depth
      diveTypes = d %>%
        dplyr::group_by(diveId) %>%
        dplyr::summarise(
          diveType = ifelse(any(diveId == 0), 'Unknown', ifelse(
            max(depth.standardized) > deep_dive_depth, 'Deep', 'Shallow'
          ))
        )
      
      # merge dive type information with sattag record
      d = d %>% left_join(diveTypes, by = 'diveId')
      
      # reformat column names and add derived data to sattag record
      d = d %>% mutate(
        depths = depth.standardized,
        ddepths = c(0, diff(depths)),
        ddepths.sign = sign(ddepths),
        descending = ddepths.sign > 0,
        ascending = ddepths.sign < 0,
        no_change = ddepths.sign == 0
      )
      
      #
      # plot data, demonstrating the raw statistics
      #
      
      pl = ggplot(d, aes(x = Date, y = depths, group = 1)) + 
        # horizontal surface line
        geom_hline(yintercept = template_bins$center[1], lty = 3) +
        # time series plot of depth bins
        geom_line() +
        geom_point() + 
        # cee start time
        geom_vline(xintercept = tag_info$cee_start, lty = 2, alpha = .5) + 
        # axis and plot formatting
        scale_x_datetime('Time (UTC)', 
                         date_breaks = 'hour',
                         date_labels = '%H%M',
                         limits = c(window_selected$pre_exposure_start, 
                                    window_selected$response_end) + 
                           c(-1,1) * 3600) + 
        scale_alpha_manual(values = c(Deep = 1, Shallow = 0, Unknown = 0)) + 
        scale_y_reverse() +
        guides(col = 'none', alpha = 'none') + 
        ylab('Depth (m)') + 
        theme_few()
      
      #
      # plot data, demonstrating the time windows
      #
      
      pl2 = ggplot(d, aes(x = Date, y = depths, group = 1)) + 
        # time windows
        geom_rect(aes(xmin = window_selected$pre_exposure_start,
                      xmax = window_selected$pre_exposure_end, 
                      ymin = 0,
                      ymax = Inf),
                  data = data.frame(x = 1),
                  inherit.aes = FALSE,
                  alpha = .15) +
        geom_rect(aes(xmin = window_selected$response_start,
                      xmax = window_selected$response_end, 
                      ymin = 0,
                      ymax = Inf),
                  data = data.frame(x = 1),
                  inherit.aes = FALSE,
                  alpha = .15) +
        geom_rect(aes(xmin = response_start,
                      xmax = response_end, 
                      ymin = 0,
                      ymax = Inf),
                  data = window_selected2,
                  inherit.aes = FALSE,
                  alpha = .15) +
        geom_rect(aes(xmin = pre_exposure_start,
                      xmax = pre_exposure_end, 
                      ymin = 0,
                      ymax = Inf),
                  data = window_selected2[1,],
                  inherit.aes = FALSE,
                  alpha = .15) +
        # horizontal surface line
        geom_hline(yintercept = template_bins$center[1], lty = 3) +
        # time series plot of depth bins
        geom_line() +
        geom_point() + 
        # cee start time
        geom_vline(xintercept = tag_info$cee_start, lty = 2, alpha = .5) + 
        # axis and plot formatting
        scale_x_datetime('Time (UTC)', 
                         date_breaks = 'hour',
                         date_labels = '%H%M',
                         limits = c(window_selected$pre_exposure_start, 
                                    window_selected$response_end) + 
                           c(-1,1) * 3600) + 
        scale_alpha_manual(values = c(Deep = 1, Shallow = 0, Unknown = 0)) + 
        scale_y_reverse() +
        guides(col = 'none', alpha = 'none') + 
        ylab('Depth (m)') + 
        theme_few()
      
      f = file.path('output','figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      ggsave(pl, filename = file.path(f, paste(tar_name(), '_', tag_id, '.pdf', 
                                               sep ='')),
             width = 8, height = 4, dpi = 'print')
      ggsave(pl2, filename = file.path(f, paste(tar_name(), 
                                                '_windows_', tag_id, '.pdf', 
                                               sep ='')),
             width = 8, height = 4, dpi = 'print')
    },
    pattern = map(tags_2019),
    memory = 'transient'
  ),
  
  tar_target(
    name = null_distribution_plots,
    command = {
      
      # plot settings
      tag = 'ZcTag093'
      stat = 'time_deep_diving'
      
      #
      # load and identify data for plotting
      #
      
      tar_load(sattag_eda)
      
      manifest = do.call(rbind, lapply(sattag_eda, function(x) {
        data.frame(x$tag, x$window, conditional = x$conditional)
      })) %>% mutate(ind = 1:n())
      
      ind = manifest %>% 
        filter(
          deployid == tag, 
          response_lag == 0, 
          window_length == 6,
          conditional == TRUE
        ) %>% 
        select(ind) %>% 
        unlist()
    
      df = data.frame(
        pre = sattag_eda[[ind]]$null_samples_pre[[stat]],
        post = sattag_eda[[ind]]$null_samples_post[[stat]]
      ) * 5 / 60
      
      observed_pre = sattag_eda[[ind]]$observed_pre[[stat]] * 5 / 60
      observed_post = sattag_eda[[ind]]$observed_post[[stat]] * 5 / 60
      
      #
      # build marginal distribution and Mahalanobis distance components
      #
      
      # bandwidths for kde, modified s.t. results are always non-zero if data
      # are not degenerate
      h = c(MASS::bandwidth.nrd(df$pre),
            MASS::bandwidth.nrd(df$post))
      if(h[1] == 0) {
        x = df$pre[,ind]
        h[1] = 4 * 1.06 * sqrt(var(x)) * length(x)^(-1/5)
      }
      if(h[2] == 0) {
        x = df$post[,ind]
        h[2] = 4 * 1.06 * sqrt(var(x)) * length(x)^(-1/5)
      }
      
      # bivariate kernel density estimate (KDE)
      dens = MASS::kde2d(
        x = df$pre, 
        y = df$post, 
        n = 500, 
        lims = c(range(c(df$pre, observed_pre)), 
                 range(c(df$post, observed_post))),
        h = h
      )
      
      # indexes to retrieve conditional distributions for KDE 
      kde_ind_pre = which.min(abs(dens$x - observed_pre))
      kde_ind_post = which.min(abs(dens$y - observed_post))
      
      # conditional distribution P(Post | Pre) and complement
      C = sum(dens$z[kde_ind_pre,])
      cdf.conditional = cumsum(dens$z[kde_ind_pre,]) / C
      ccdf.conditional = rev(cumsum(rev(dens$z[kde_ind_pre,]))) / C
      
      # mahalanobis distance calculation components
      mdist_mean = colMeans(df)
      mdist_var = var(df)
      md = mahalanobis(
        x = rbind(df, cbind(pre = observed_pre, post = observed_post)), 
        center = mdist_mean, 
        cov = mdist_var
      )
      
      # 
      # visualize tests and components!
      #
      
      # bivariate data plot
      pl = ggplot(df, aes(x = pre, y = post)) + 
        # overplotted raw data
        geom_point(col = 'grey', alpha = .3) +
        # bivariate density
        stat_density_2d(col = 'black') +
        # conditional baseline distribution used for test
        geom_vline(xintercept = observed_pre, col = 'darkgreen') +
        # pre/post pair observed during CEE
        geom_point(
          data = data.frame(
            pre = observed_pre, 
            post = observed_post
          ),
          col = 'darkgreen'
        ) + 
        # axis labels and formatting
        # xlab(expression(t(w[pre]))) + 
        # ylab(expression(t(w[post]))) +
        xlab(expression(s[31]^(i)))+ 
        ylab(expression(s[32]^(i)))+ 
        theme_few() + 
        theme(axis.title.y = element_text(angle = 0, vjust = .5)) + 
        coord_equal() # + 
        # ggtitle('Time (h) deep diving within 6h window')
      
      # univariate data plot
      df_uni = data.frame(post = dens$y, d = dens$z[kde_ind_pre,] / C)
      pl_uni = ggplot(df_uni, aes(x = post, y = d)) + 
        # baseline observations associated with pre-exposure conditions
        geom_rug(
          sides = 'b', 
          data = df %>% filter(pre == observed_pre), 
          inherit.aes = FALSE, 
          mapping = aes(x = post),
          col = 'grey',
          alpha = .4
        ) + 
        # estimated density 
        geom_line() +
        # observed response value
        geom_point(
          data = data.frame(
            post = dens$y[kde_ind_post],
            d = dens$z[kde_ind_pre, kde_ind_post] / C
          ),
          col = 'darkgreen'
        ) +
        # shade the observed p-value
        geom_area(
          data = df_uni %>% filter(post <= dens$y[kde_ind_post]),
          col = 'darkgreen',
          fill = 'darkgreen'
        ) + 
        # formatting
        xlab(expression(t(w[post]))) + 
        # ylab(expression('['~t(w[post])~'|'~t(w[pre])~']')) +
        ylab(expression(f(t(w[post])~'|'~t(w[pre])))) +
        theme_few() + 
        theme(axis.title.y = element_text(angle = 0, vjust = .5)) + 
        ggtitle('Time (h) deep diving within 6h window')
      
      
      #
      # mahalanobis distance contour comparison
      #
      
      ellipsePoints = function(mu, Sigma, segments = 100, scale = 1) {
        angles = (0:segments) * 2 * pi / segments
        unit.circle = cbind(cos(angles), sin(angles)) * scale
        t(mu + t(unit.circle %*% Sigma))
      }
      
      pl_mdist = pl +
        # center
        geom_point(
          data = data.frame(pre = mdist_mean['pre'],
                            post = mdist_mean['post']), 
          col = 'darkred'
        ) +
        # contours 
        geom_polygon(
          data = do.call(rbind, lapply(c(1,3,5), function(sc) {
            data.frame(
              ellipsePoints(mu = mdist_mean, Sigma = mdist_var, scale = sc),
              scale = sc
            )
          })),
          mapping = aes(group = scale),
          fill = NA,
          col = 'darkred'
        )
      
      #
      # univariate Mahalanobis distribution
      #
      
      dens_mdist = density(md[1:nrow(df)])
      
      md_ind = which.min(abs(dens_mdist$x - tail(md, 1)))
      
      pl_mdist_uni = ggplot(
        data.frame(mdist = dens_mdist$x, d = dens_mdist$y), 
        aes(x = mdist, y = d)
      ) + 
        # baseline distribution
        geom_line() + 
        # data points
        geom_rug(
          data = data.frame(mdist = md[1:nrow(df)]),
          inherit.aes = FALSE,
          mapping = aes(x = mdist),
          sides = 'b',
          col = 'grey',
          alpha = .4
        ) + 
        # p-value and shading
        geom_area(
          data = data.frame(mdist = dens_mdist$x, d = dens_mdist$y) %>% 
            slice(md_ind : n()),
          col = 'darkgreen',
          fill = 'darkgreen'
        ) + 
        geom_point(
          data = data.frame(mdist = dens_mdist$x, d = dens_mdist$y)[md_ind,],
          col = 'darkgreen'
        ) + 
        # formatting
        xlab(expression(D[M](t(w[pre]),t(w[post])))) + 
        ylab(expression('['~D[M](symbol("\327"))~']')) +
        theme_few() + 
        ggtitle('Mahalanobis distance of time deep diving within 6h window')
        
      #
      # save plots
      #
      
      f = file.path('output','figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      ggsave(pl + theme(text = element_text(size = 20), 
                        axis.title = element_text(size = 30)), 
             filename = file.path(f, paste(tar_name(), '_bivariate_', tag, 
                                               '.pdf', sep ='')),
             width = 8, height = 8, dpi = 'print')
      
      ggsave(pl_uni + theme(text = element_text(size = 13), 
                             axis.title = element_text(size = 30)), 
             filename = file.path(f, paste(tar_name(), '_univariate_', 
                                                   tag, '.pdf', sep ='')),
             width = 8, height = 6, dpi = 'print')
      
      ggsave(pl_mdist + theme(text = element_text(size = 20), 
                              axis.title = element_text(size = 30)), 
             filename = file.path(f, paste(tar_name(), '_mdist_', tag, 
                                               '.pdf', sep ='')),
             width = 8, height = 8, dpi = 'print')
      
      ggsave(pl_mdist_uni  + theme(text = element_text(size = 13), 
                                   axis.title = element_text(size = 30)), 
             filename = file.path(f, paste(tar_name(), '_mdist_univariate_', 
                                                   tag, '.pdf', sep ='')),
             width = 8, height = 6, dpi = 'print')
     }
  ),
  
  tar_target(
    name = population_signif_trends, 
    command = {
      
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
        summarise(prop_signif_animals = sum(p < .1) / length(tag))
      
      levels(aggregate_results$stat) = sub_names[
        gsub('_', ' ', gsub('timing_', '', levels(aggregate_results$stat)))
      ]
      
      pl = ggplot(aggregate_results %>% filter(sequential == FALSE) %>% 
                    mutate(prop_signif_animals = prop_signif_animals + 
                             prop_signif_animals / prop_signif_animals - 1),
             aes(x = description, y = stat, fill = prop_signif_animals)) + 
        geom_raster() + 
        scale_fill_distiller(
          'Prop. animals with small tail probs. (p<.1)',
          palette = 'YlOrBr',
          direction = 1, 
          na.value = '#ffffff'
        ) +
        scale_x_discrete(position = 'top') + 
        theme_few() + 
        theme(axis.title = element_blank()) + 
        ggtitle(expression('Aggregate statistics: '~t[i](w)==sum(Y[j]^(i),j))) +
        coord_equal()
      
      pl_seq = ggplot(aggregate_results %>% filter(sequential == TRUE) %>% 
                    mutate(prop_signif_animals = prop_signif_animals + 
                             prop_signif_animals / prop_signif_animals - 1),
                  aes(x = description, y = stat, fill = prop_signif_animals)) + 
        geom_raster() + 
        scale_fill_distiller(
          'Prop. responding animals (p<.1)',
          palette = 'YlOrBr',
          direction = 1, 
          na.value = '#ffffff'
        ) +
        scale_x_discrete(position = 'top') + 
        theme_few() + 
        theme(axis.title = element_blank()) + 
        ggtitle(
          expression('Sequential statistics: '~t[i](w)==sum(j*Y[j]^(i),j))
        ) +
        coord_equal()
      
      #
      # save plots
      #
      
      f = file.path('output','figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      ggsave(pl + theme(text = element_text(size = 20)), 
             filename = file.path(f, paste(tar_name(), '_aggregate.png',
                                           sep ='')),
             width = 16, height = 8, dpi = 'print')
      
      ggsave(pl_seq + theme(text = element_text(size = 20)), 
             filename = file.path(f, paste(tar_name(), '_sequential.png',
                                           sep ='')),
             width = 16, height = 8, dpi = 'print')
      
      
      #
      # reformat results in tabular format
      #
      
      
      

      NULL
    }
  ),
  
  tar_target(
    name = result_similarity_trends,
    command = {

      sim_subset = result_similarity[tags_2019, tags_2019]
      
      sim_subset = result_similarity
      
      df = data.frame(sim_subset, from = rownames(sim_subset)) %>%
        pivot_longer(cols = starts_with('ZcTag'), names_to = 'to', 
                     values_to = 'sim') %>% 
        mutate(from = gsub(pattern = 'ZcTag', replacement = '', x = from),
               to = gsub(pattern = 'ZcTag', replacement = '', x = to))
      
      #
      # confusion matrix plot
      #
      
      pl = ggplot(df, aes(x = from, y = to, fill = sim)) + 
        geom_tile() + 
        scale_fill_viridis(expression(italic(D)(i,j))) + 
        coord_equal() + 
        xlab('ZcTag (i)') + 
        ylab('ZcTag (j)') +
        theme_few() + 
        theme(
          axis.title.y = element_text(angle = 0, vjust = .5),
          panel.border = element_blank()
        )
      
      f = file.path('output','figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      ggsave(pl, filename = file.path(f, paste(tar_name(), '.pdf', sep = '')))
      
      #
      # nearest neighbor distance plot
      #
      
      df2 = df %>% 
        group_by(from) %>% 
        arrange(from, sim) %>% 
        mutate(nn = cumsum(from != to)) %>% 
        ungroup()
      
      pl2 = ggplot(df2 %>% filter(nn > 0), aes(x = nn, y = sim, col = from)) + 
        # geom_point() + 
        geom_line(alpha = .6) + 
        scale_color_discrete('ZcTag') + 
        scale_x_continuous('Nth most similar neighbor', 
                           breaks = 1:max(df2$nn)) + 
        ylab('Dissimilarity') + 
        theme_few()
      
      ggsave(pl2, 
             filename = file.path(f, paste(tar_name(), '_line.pdf', sep = '')),
             width = 8, height = 6)
        
       NULL
    }
  )
)
