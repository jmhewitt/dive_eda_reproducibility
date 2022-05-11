data_targets = list(

  # seconds between depth observations
  tar_target(sattag_timestep, 300),
  
  # deep dive threshold (m)
  tar_target(deep_dive_depth, 800),
  
  # location of sattag files
  tar_target(sattag_dir, file.path('data', 'sattag')),
  
  # sattag series data
  tar_target(
    name = depth_files, 
    command = dir(path = sattag_dir, pattern = 'Series\\.csv', 
                  full.names = TRUE),
    format = 'file'
  ),
  
  # sattag messages
  tar_target(
    name = message_files, 
    command = dir(path = sattag_dir, pattern = 'SeriesRange\\.csv', 
                  full.names = TRUE),
    format = 'file'
  ),
  
  # formats of which tag metadata are stored in
  tar_target(
    name = date_formats, 
    command = c('mdy IMS p', 'ymd HMS')
  ),
  
  
  tar_target(
    name = tag_info_file,
    command = file.path(sattag_dir, 'tag_info.csv'),
    format = 'file'
  ),
  
  # tag metadata: sex and CEE info
  tar_target(
    name = tag_info, 
    command = read.csv(tag_info_file, colClasses = 'factor') %>%
        dplyr::mutate(
          baseline_end = lubridate::parse_date_time(
            x = baseline_end, orders = date_formats, tz = 'UTC'
          ),
          cee_start = lubridate::parse_date_time(
            x = cee_start, orders = date_formats, tz = 'UTC'
          )
        )
  ),
  
  # tag deployid
  tar_target(
    name = tag_id,
    command = tag_info$deployid
  ),
  
  # target 2019 animals to analyze
  tar_target(
    name = tags_2019,
    command = c('ZcTag083', 'ZcTag085', 'ZcTag087', 'ZcTag088', 'ZcTag089',
                'ZcTag093', 'ZcTag095', 'ZcTag096', 'ZcTag097')
  ),
  
  tar_target(
    name = tag_meta,
    command = {
        do.call(rbind, lapply(depth_files, function(f) {
          
          d = read.csv(f)
          d$Date = as.POSIXct(d$Date, origin = '1970-01-01 00:00.00 UTC', 
                              tz = 'UTC')
          
          id = gsub(
            pattern = '_DUML', 
            replacement = '', 
            x = d$DeployID[1]
          )
          
          cee_start = tag_info %>% 
            filter(deployid == id) %>% 
            select(cee_start) %>% 
            unname()
          
          if(nrow(cee_start) == 0) {
            cee_start = NA
            Npre = nrow(d)
            cee_num = NA
          } else {
            Npre = d %>% filter(Date < as.numeric(cee_start)) %>% nrow()
            cee_num = tag_info %>% 
              filter(deployid == id) %>% 
              select(cee_id) %>% 
              unname()
          }
          
          data.frame(
            Tag = id,
            # Data collection columns
            Start = range(d$Day)[1], 
            End = range(d$Day)[2],
            # CEE columns
            Number = cee_num,
            Time = cee_start,
            # data information
            Npre = Npre,
            Ntot = nrow(d)
          )
        }))
    }
  )
  
)
