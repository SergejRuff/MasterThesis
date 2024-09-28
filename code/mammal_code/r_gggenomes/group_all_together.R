# # Function to check overlap and update positions
# update_overlaps <- function(group) {
#   # Check for overlaps and update positions
#   new_start <- min(group$start)
#   new_end <- max(group$end)
#   
#   # Determine the representative name
#   rep_name <- group %>%
#     count(interpro_description) %>%
#     arrange(desc(n)) %>%
#     slice(1) %>%
#     pull(interpro_description)
#   
#   # Return a single row with updated values
#   return(group %>%
#            slice(1) %>%
#            mutate(start = new_start, end = new_end, interpro_description = rep_name))
# }
# 
# # Apply the function to each group of seq_id and orf_start
# selected_rows <- selected_rows %>%
#   group_by(seq_id, orf_start) %>%
#   group_split() %>%
#   lapply(update_overlaps) %>%
#   bind_rows()
# 
# # View the updated dataframe
# print(selected_rows_updated)


# Update specific rows
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(seq_id == 'SRR13364363_cap3_Contig-11', 27, start),
    end = if_else(seq_id == 'SRR13364363_cap3_Contig-11', 420, end),
    interpro_description = if_else(seq_id == 'SRR13364363_cap3_Contig-11', 'DNA/RNA polymerase superfamily', interpro_description)
  )%>%
  group_by(seq_id, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


selected_rows <- selected_rows %>%
  mutate(
    start = if_else(seq_id == 'SRR13364364_cap3_Contig-2', 31, start),
    end = if_else(seq_id == 'SRR13364364_cap3_Contig-2', 781, end),
    interpro_description = if_else(seq_id == 'SRR13364364_cap3_Contig-2', 'DNA/RNA polymerase superfamily', interpro_description)
  )%>%
  group_by(seq_id, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()



# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(seq_id == 'SRR17345012_cap3_Contig-13' & orf_start == 2356, 2530, start),
    end = if_else(seq_id == 'SRR17345012_cap3_Contig-13' & orf_start == 2356, 2770, end),
    interpro_description = if_else(seq_id == 'SRR17345012_cap3_Contig-13' & orf_start == 2356, "DNA/RNA polymerase superfamily", interpro_description)
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()



# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(seq_id == 'SRR17345013_cap3_Contig-21' & orf_start == 2521, 2648, start),
    end = if_else(seq_id == 'SRR17345013_cap3_Contig-21' & orf_start == 2521, 2941, end),
    interpro_description = if_else(seq_id == 'SRR17345013_cap3_Contig-21' & orf_start == 2521, "DNA/RNA polymerase superfamily", interpro_description)
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(seq_id == 'SRR17345012_cap3_Contig-11' & orf_start == 2716, 2971, start),
    end = if_else(seq_id == 'SRR17345012_cap3_Contig-11' & orf_start == 2716, 3256, end),
    interpro_description = if_else(seq_id == 'SRR17345012_cap3_Contig-11' & orf_start == 2716, "DNA/RNA polymerase superfamily", interpro_description)
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()

##########################################


# Define the new values, specific orf_start, and interval
new_start <- 4
new_end <- 338
new_description <- 'Alphavirus-like methyltransferase (MT) domain'
specific_orf_start <- 3
start_interval <- c(3, 16)  # Example interval for start
end_interval <- c(198, 340)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


###################################################################


# Define the new values, specific orf_start, and interval
new_start <- 1338
new_end <- 1637
new_description <- 'RNA-directed RNA polymerase,  catalytic domain'
specific_orf_start <- 3
start_interval <- c(1330, 1424)  # Example interval for start
end_interval <- c(1533, 1638)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()




#################################################################


# Define the new values, specific orf_start, and interval
new_start <- 903
new_end <- 1185
new_description <- '(+) RNA virus helicase core domain'
specific_orf_start <- 3
start_interval <- c(902, 944)  # Example interval for start
end_interval <- c(1152, 1186)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()

#############################################################


# Define the new values, specific orf_start, and interval
new_start <- 936
new_end <- 1171
new_description <- 'P-loop containing nucleoside triphosphate hydrolase'
specific_orf_start <- 3
start_interval <- c(936-1, 1052+1)  # Example interval for start
end_interval <- c(1038-1, 1171+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()

#########################################################################


# Define the new values, specific orf_start, and interval
new_start <- 744
new_end <- 890
new_description <- 'Macro domain'
specific_orf_start <- 3
start_interval <- c(744-1, 756+1)  # Example interval for start
end_interval <- c(871-1, 890+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR17345012_cap3_Contig-6' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()

##############################################################################


# Define the new values, specific orf_start, and interval
new_start <- 4212
new_end <- 4510
new_description <- 'P-loop containing nucleoside triphosphate hydrolase'
specific_orf_start <- 4118
start_interval <- c(4212-1, 4363+1)  # Example interval for start
end_interval <- c(4352-1, 4510+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR18779479_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR18779479_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR18779479_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


#############################################################################

# Define the new values, specific orf_start, and interval
new_start <- 4639
new_end <- 5066
new_description <- 'DNA/RNA polymerase superfamily'
specific_orf_start <- 4118
start_interval <- c(4639-1, 4837+1)  # Example interval for start
end_interval <- c(4949-1, 5066+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR18779479_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR18779479_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR18779479_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


##################################################################################

# Define the new values, specific orf_start, and interval
new_start <- 2766
new_end <- 2793
new_description <- 'Zinc finger, CCCH-type'
specific_orf_start <- 2747
start_interval <- c(2766-1, 2769+1)  # Example interval for start
end_interval <- c(2791-1, 2793+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR16151768_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR16151768_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR16151768_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


################################################################################


# Define the new values, specific orf_start, and interval
new_start <- 3693
new_end <- 5837
new_description <- 'DNA/RNA polymerase superfamily'
specific_orf_start <- 3669
start_interval <- c(3693-1, 5497+1)  # Example interval for start
end_interval <- c(4554-1, 5837+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR16151768_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR16151768_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR16151768_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


#######################################################################


# Define the new values, specific orf_start, and interval
new_start <- 3619
new_end <- 4201
new_description <- 'DNA/RNA polymerase superfamily'
specific_orf_start <- 338
start_interval <- c(3619-1, 3863+1)  # Example interval for start
end_interval <- c(3909-1, 4201+1)    # Example interval for end

selected_rows$interpro_description <- gsub("Reverse transcriptase/Diguanylate cyclase domain", "DNA/RNA polymerase superfamily", 
                                           selected_rows$interpro_description)

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()

##############################################################################


# Define the new values, specific orf_start, and interval
new_start <- 2121
new_end <- 2487
new_description <- 'P-loop containing nucleoside triphosphate hydrolase'
specific_orf_start <- 338
start_interval <- c(2121-1, 2345+1)  # Example interval for start
end_interval <- c(2291-1, 2487+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


#############################################################################


# Define the new values, specific orf_start, and interval
new_start <- 1038
new_end <- 1408
new_description <- 'Pestivirus envelope glycoprotein E2'
specific_orf_start <- 338
start_interval <- c(1038-1, 1277+1)  # Example interval for start
end_interval <- c(1124-1, 1408+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()

######################################################################


# Define the new values, specific orf_start, and interval
new_start <- 1785
new_end <- 2144
new_description <- 'Pestivirus NS3, peptidase S31'
specific_orf_start <- 338
start_interval <- c(1785-1, 2082+1)  # Example interval for start
end_interval <- c(1933-1, 2144+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR15647801_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()


#############################################################################


# Define the new values, specific orf_start, and interval
new_start <- 58
new_end <- 1635
new_description <- 'Nonstructural protein 2, N-terminal domain, coronavirus'
specific_orf_start <- 3
start_interval <- c(58-1, 787+1)  # Example interval for start
end_interval <- c(420-1, 1635+1)    # Example interval for end

# Update the rows and ensure only one row remains
selected_rows <- selected_rows %>%
  mutate(
    start = if_else(
      seq_id == 'SRR14579882_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_start, start
    ),
    end = if_else(
      seq_id == 'SRR14579882_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_end, end
    ),
    interpro_description = if_else(
      seq_id == 'SRR14579882_cap3_Contig-1' & 
        orf_start == specific_orf_start &
        start >= start_interval[1] & start <= start_interval[2] &
        end >= end_interval[1] & end <= end_interval[2],
      new_description, interpro_description
    )
  ) %>%
  group_by(seq_id, orf_start, start, end, interpro_description) %>%
  filter(row_number() == 1) %>%  # Keep only the first row in each group
  ungroup()