table_targets = list(
  
  tar_target(
    name = ms_tags,
    command = {
      
      x = print.xtable(
        xtable(
          tag_meta %>% 
            mutate(
              Tag = gsub(pattern = '_DUML', replacement = '', x = Tag),
              Npre = paste('$', format(Npre, big.mark = ','), '$', sep = ''),
              Ntot = paste('$', format(Ntot, big.mark = ','), '$', sep = ''),
              Date = strftime(x = Time, format = '%d-%b'),
              Time = strftime(x = Time, format = '%H%M'),
              Start = gsub(pattern = '-2019', replacement = '', x = Start),
              End = gsub(pattern = '-2019', replacement = '', x = End)
            ) %>% 
            filter(sapply(Tag, function(x) x %in% tags_2019)) %>% 
            select(
              ID=Tag, Start, End, "$U$"=Npre, '$T$'=Ntot, "ID "=Number, Date, 
              "Time (UTC)" = Time
            ), 
          caption = paste(
            'Summary of tags and CEEs analyzed in study, including start and ',
            'end dates ',
            'for data collection, and start date and time for CEE.  All data ',
            'were collected in 2019.  The variables $U$ and $T$ denote the ',
            'number of pre-exposure observations used to estimate baseline ',
            'distributions, and the total number of observations for the tag, ',
            'respectively.',
            sep = ''
          ), 
          align = 'cccccc|ccc', 
          label = 'table:tag_info'
        ), 
        booktabs = TRUE, 
        include.rownames = FALSE,
        sanitize.text.function = identity
      )
      
      x = stringr::str_replace(
        string = x, 
        pattern = fixed('\\toprule'),
        replacement = paste(
          '\\toprule\n\\multicolumn{5}{c|}{Tag info.} & ',
          '\\multicolumn{3}{c}{CEE info.} \\\\ \n',
          sep = ''
        )
      )
      
      write_clip(x)
      
      x
    }
  ),
  
  tar_target(
    name = result_similarity,
    command = {

      meta = do.call(rbind, lapply(eda_results, function(x) 
        data.frame(tag = x$tag_id, conditional = x$conditional)
      ))
      
      tags = unique(meta$tag)
      
      sim = matrix(nrow = length(tags), ncol = length(tags))
      rownames(sim) = tags
      colnames(sim) = rownames(sim)
      
      for(i in 1:nrow(sim)) {
        for(j in 1:ncol(sim)) {
          t1ind = which(meta$tag == tags[i] & meta$conditional == TRUE)
          t2ind = which(meta$tag == tags[j] & meta$conditional == TRUE)
          sim[i, j] = eda_results[[t1ind]]$df.bivariate.raw %>% 
            mutate(p_std = ifelse(right_tail, 1 - p, p)) %>% 
            left_join(
              eda_results[[t2ind]]$df.bivariate.raw %>% 
                mutate(p_std = ifelse(right_tail, 1 - p, p)),
              by = c('response_lag', 'window_length', 'stat')
            ) %>% 
            summarize(mean((p.x-p.y)^2)) %>% 
            as.numeric()
        }
      }
      
      sim
    }
  )

)
