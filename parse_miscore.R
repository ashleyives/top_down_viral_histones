
library(tidyverse)

#function that will parse MIScores into a more useable series of list-based columns
#needs columns Proteform, firstAA, lastAA, protLength (from mapping fasta)
#note that position in MIScore is based on proteoform sequence not gene sequence from FASTA, so need to correct that 
#also classifies as proteoform level, if site is known = 1, if localized range = 2
#assumes there is only one mod in MIScore 
#assumes miscores look like "IronAdduct[E41:19.9%;E56:19.9%;E61:19.9%;E64:19.9%;D72:19.9%]" OR "Phosph[T16:92.0%]"

parse_miscore <- function(x) {
  x <- x %>%
    mutate(mi_mod = sub("\\[.*", "", MIScore),  # Extract everything before "["
           mi_location = sub(".*\\[(.*)\\].*", "\\1", MIScore), # Extract everything inside "[...]"
           mi_location_list = strsplit(mi_location, ";")        # Split mi_location into a list using ";"
    )
  
  #remove percentages from MISCore
  x$mi_location_list <- lapply(x$mi_location_list, function(loc) {
    sapply(loc, function(entry) sub(":.*", "", entry))
  })
  
  # Create mi_adj_location_list by adjusting the numerical value using firstAA - 1
  x$mi_adj_location_list <- lapply(seq_along(x$mi_location_list), function(idx) {
    loc_list <- x$mi_location_list[[idx]]
    firstAA <- x$firstAA[idx]
    sapply(loc_list, function(entry) {
      letter <- substr(entry, 1, 1)  # Extract the letter (first character)
      number <- as.numeric(substr(entry, 2, nchar(entry))) # Extract the number (rest of the string)
      adjusted_number <- number + firstAA - 1
      paste0(letter, adjusted_number)  # Combine the letter and adjusted number
    })
  })
  
  #adds a range for a given mod as start and stop site 
  x <- x %>%
    mutate(
      mi_adj_location_start = sapply(mi_adj_location_list, function(loc) {
        # Extract numeric values, ignoring letters
        numeric_entries <- as.numeric(gsub("[^0-9]", "", loc))
        if (all(!is.na(numeric_entries))) {
          min(numeric_entries)  # Return the smallest value
        } else {
          NA  # Handle cases with no valid numbers
        }
      }),
      mi_adj_location_stop = sapply(mi_adj_location_list, function(loc) {
        # Extract numeric values, ignoring letters
        numeric_entries <- as.numeric(gsub("[^0-9]", "", loc))
        if (all(!is.na(numeric_entries))) {
          max(numeric_entries)  # Return the largest value
        } else {
          NA  # Handle cases with no valid numbers
        }
      })
    )
  
  #add proteoform level 
  x <- x %>%
    mutate(mi_level = sapply(mi_location_list, function(loc) {
      if (length(loc) == 1) {
        return(1)  # Level 1 if there's only one entry
      } else {
        return(2)  # Level 2 if there are multiple entries
      }
    }))
  
  
  return(x)
}

