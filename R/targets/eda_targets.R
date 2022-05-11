eda_targets = list(
  
  tar_target(
    name = ts_plots,
    command = {

      # skip processing if no information
      if(is.null(pre_post_pair_defs)) {
        return(NULL)
      }
      
      # identify depth record file for tag
      tag_info = pre_post_pair_defs[[1]]$tag
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
      
      # time series of test statistics
      ts_samples = do.call(rbind, lapply(
        pre_post_pair_defs[[1]]$all_pairs, function(x) {
          data.frame(
            time = x$context_window[2],
            sattag_summary(d[x$study_inds,]) - 
            sattag_summary(d[x$context_inds,])
          )
        }
      ))
      
      # output location
      f = file.path('output', 'eda', 'series', tag_info$deployid)
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      # plot the results; split plots by statistic
      manifest = do.call(rbind, lapply(colnames(ts_samples)[-1], function(stat) {
        # build plot
        pl = ggplot(ts_samples, aes_string(x = 'time', y = stat)) + 
          geom_line() + 
          geom_point() + 
          geom_vline(xintercept = tag_info$cee_start, lty = 3) + 
          scale_x_datetime('Time (UTC)', 
                           date_breaks = '6 hours', 
                           date_labels = '%m-%d %Hh') + 
          ylab(expression(s(bold(X)))) + 
          ggtitle(paste(pre_post_pair_defs[[1]]$tag$deployid, ', ',
                        pre_post_pair_defs[[1]]$window$response_lag, 'h lag, ', 
                        pre_post_pair_defs[[1]]$window$window_length, 'h window ', 
                        sep = ''), 
                  subtitle = paste(stat, ', exposure at dotted line', 
                                   ifelse(pre_post_pair_defs[[1]]$baseline_pairs_are_conditional,
                                          ' (conditional)', ' (unconditional)'),
                                   sep = '')) + 
          theme_few() + 
          theme(axis.text.x = element_text(angle = 90),
                axis.title.y = element_text(angle = 0, vjust = .5))
        # save plot
        fname = file.path(
          f, 
          paste(pre_post_pair_defs[[1]]$tag$deployid, '_',
                pre_post_pair_defs[[1]]$window$description, '_',
                stat, 
                ifelse(pre_post_pair_defs[[1]]$baseline_pairs_are_conditional,
                       '_conditional_', '_unconditional_'),
                '.png', 
                sep = '')
        )
        ggsave(pl, filename = fname, dpi = 'print', width = 14, height = 7)
        
        data.frame(deployid = pre_post_pair_defs[[1]]$tag$deployid, 
                   response_lag = pre_post_pair_defs[[1]]$window$response_lag, 
                   window_length = pre_post_pair_defs[[1]]$window$window_length,
                   stat = stat,
                   conditional = pre_post_pair_defs[[1]]$baseline_pairs_are_conditional,
                   path = fname,
                   type = 'tsplots')
      }))
      
      manifest
    },
    pattern = map(pre_post_pair_defs),
    deployment = 'worker'
  )
)
  