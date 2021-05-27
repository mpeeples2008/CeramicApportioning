require(msm)

apportion <-
  function(Site,      # vector of site names with one entry for every site/type/count
           Type,      # vector of type names with one entry for every site/type/count
           ct,        # vector of counts with one entry for every site/type/count
           beg.date,  # beginning date of the interval to be considered
           end.date,  # end date of the interval to be considered
           interval,  # interval to define period length for apportioning
           site.beg,  # vector of site beginning dates with one entry for every site/type/count
           site.end,  # vector of site ending dates with one entry for every site/type/count
           type.beg,  # vector of type beginning dates with one entry for every site/type/count
           type.end)  # vector of type ending dates with one entry for every site/type/count
                    {
    
    
    # calculate the number of periods of length "interval" and round up to nearest whole number
    nperiods <- ceiling((end.date - beg.date) / interval)
    
    # Create column labels for apportioned  matrix for each interval
    c.lab <- NULL
    for (i in 1:nperiods) {
      c.lab[i] <- beg.date + ((i - 1) * interval)
    }
    
    # create short variable names to satisfy my need for tidiness (and to match variable names in Roberts et al. 2012 paper)
    t0 <- type.beg
    t1 <- type.end
    s0 <- site.beg
    s1 <- site.end
    
    # set truncation points for normal distribution at 2 standard deviations
    zsd <- 2
    
    # Create empty matrix of 0s with a row for each Site/Type/Count and columns for each period defined above and apply column labels
    mat <- matrix(0, length(Site), nperiods)
    colnames(mat) <- c.lab
    
    # Create an additional vector for ceramics not apportioned to any period considered here
    # This vector will contain counts for ceramics outside of the site occupation span or outside of the user defined period.
    unapportioned <- rep(0, length(Site))
    
    # round site occupation start dates to multiples of interval
    for (i in 1:length(s0)) {
      if ((s0[i] / interval) != (round(s0[i] / interval, 0))) {
        s0[i] <- interval * floor(s0[i] / interval)
      }
    }
    
    # round site occupation end dates to multiples of interval
    for (i in 1:length(s1)) {
      if ((s1[i] / interval) != (round(s1[i] / interval, 0))) {
        s1[i] <- interval * ceiling(s1[i] / interval)
      }
    }
    
    # Apply apportioning procedure in a loop across every Site/Type/Count
    for (j in 1:length(Site)) {
      
      # create a sequence by 1 for every year from the type begin date to end date
      type_seq <- seq(t0[j], t1[j], by = 1)
      # create a sequence from -2 to 2 at intervals of 0.0001 to define positions in the truncated normal distribution
      zscore <- seq(-zsd, zsd, by = 0.0001)
      # create a vector of probabilities for the truncated normal distribution at positions defined by "zscore"
      pdist <- ptnorm(zscore, 0, 1, -zsd, zsd)
      
      # define length of period for type
      g <- t1[j] - t0[j]
      # define midpoint type date
      m <- mean(type_seq)
      # define a sequence by 1 from site occupation start to end date
      site_seq <- seq(s0[j], s1[j])
      
      
      # determine the length of intersection between the type date and the site date
      ws_overlap <- intersect(type_seq, site_seq)
      
      # If type dates and sites dates intersect, do the following:
      if (length(ws_overlap) > 0) {
        vj0 <- min(ws_overlap) # minimum date of overlap
        vj1 <- max(ws_overlap) # maximum date of overlap
        
        z20 <- (vj0 - m) / (g / (2 * zsd)) # transform start date of overlap into z values
        z21 <- (vj1 - m) / (g / (2 * zsd)) # transform end date of overlap into z values
        i0 <- which(abs(zscore - z20) == min(abs(zscore - z20))) # find the closest z score value to z20 defined above
        i1 <- which(abs(zscore - z21) == min(abs(zscore - z21))) # find the closest z score value to z21 defined above
        z_phi20 <- pdist[i0] # assign the probability associated with the selected z score value for start of overlap
        z_phi21 <- pdist[i1] # assign the probability associated with the selected z score value for end of overlap
        
        
        # Loop through each period and assign apprtioned counts as ct * phi value 
        for (i in 1:nperiods) {
          # if statement determines if a particular interval should be considered
          if (length(intersect(seq(c.lab[i], (c.lab[i] + interval - 1)), ws_overlap[1:length(ws_overlap) - 1])) > 0) {
            zj0 <- (c.lab[i] - m) / (g / (2 * zsd))
            zj1 <- ((c.lab[i] + interval) - m) / (g / (2 * zsd))
            i0 <- which(abs(zscore - zj0) == min(abs(zscore - zj0)))
            i1 <- which(abs(zscore - zj1) == min(abs(zscore - zj1)))
            z_phi0 <- pdist[i0]
            z_phi1 <- pdist[i1]
            pjt_1a <- (z_phi1 - z_phi0) / (z_phi21 - z_phi20)
            mat[j, i] <- as.numeric(pjt_1a * ct[j])
          } # end if statement
        } # end for loop
      } # end if statement
      
      # add any portion of the ceramics that extends beyond the period considered to the unapportioned vector
      if (s0[j] < beg.date ||
          s1[j] > end.date) {
        unapportioned[j] <- ct[j] - sum(mat[j, ])
      }
      
      # assign any type count that doesn't overlap with site occupation to unapportioned vector
      if (length(ws_overlap) == 0) {
        unapportioned[j] <- unapportioned[j] + ct[j]
      }
    } # end for loop
    
    # round numbers to 3 digits combine vectors and output
    mat <- round(mat, 3)
    unapportioned <- round(unapportioned, 3)
    mat <- cbind(Site, Type, t0, t1, unapportioned, mat)
    mat <- as.data.frame(mat)
    colnames(mat)[1:4] <- c('Site','Type','Start_Date','End_Date')
    mat[,3:ncol(mat)] <- do.call(cbind,lapply(mat[,3:ncol(mat)],as.numeric))
    return(mat)
} # end function

