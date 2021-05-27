## This function is modified based on the iterative proportional fitting procedure posted by
## Eddie Hunsinger: https://edyhsgr.github.io/datafitting.html

ipf_apportion <-
  function(apportion_output,
           pop,
           maxiter = 50,
           closure = 0.0000001,
           debugger = FALSE) {
    
    temp.attr <- apportion_output[,1:5]
    seed <- apportion_output[,6:ncol(apportion_output)]
    rowcontrol <- as.vector(rowSums(seed))
    colcontrol <- pop/sum(pop) * sum(seed)
    
    # input data checks:
    if (debugger)
      print("checking inputs")
    #sum of marginal totals equal and no zeros in marginal totals
    if (debugger) {
      print("checking rowsum=colsum")
    }
    if (round(sum(rowcontrol), 4) != round(sum(colcontrol), 4))
      stop("sum of rowcontrol must equal sum of colcontrol")
    if (debugger) {
      print("checking rowsums for zeros")
    }
    if (any(rowcontrol == 0)) {
      numzero <- sum(rowcontrol == 0)
      rowcontrol[rowcontrol == 0] <- 0.0000001
      warning(paste(
        numzero,
        "zeros in rowcontrol argument replaced with 0.000001",
        sep = " "
      ))
    }
    if (debugger) {
      print("Checking colsums for zeros")
    }
    if (any(colcontrol == 0)) {
      numzero <- sum(colcontrol == 0)
      colcontrol[colcontrol == 0] <- 0.0000001
      warning(paste(
        numzero,
        "zeros in colcontrol argument replaced with 0.0000001",
        sep = " "
      ))
    }
    if (debugger) {
      print("Checking seed for zeros")
    }
    if (any(seed == 0)) {
      numzero <- sum(seed == 0)
      seed[seed == 0] <- 0.0000001
      warning(paste(numzero, "zeros in seed argument replaced with 0.0000001", sep =
                      " "))
    }
    
    # set initial values
    result <- seed
    rowcheck <- 1
    colcheck <- 1
    checksum <- 1
    iter <- 0
    
    # successively proportion rows and columns until closure or iteration criteria are met
    ##########
    if (debugger) {
      print(checksum > closure)
      print(iter < maxiter)
    }
    while ((checksum > closure) && (iter < maxiter))
    {
      #########
      if (debugger) {
        print(paste("(re)starting the while loop, iteration=", iter))
      }
      
      coltotal <- colSums(result)
      colfactor <- colcontrol / coltotal
      result <- sweep(as.matrix(result), 2, as.matrix(colfactor), "*")
      if (debugger) {
        print(paste("column factor = ", colfactor))
        print(result)
      }
      
      rowtotal <- rowSums(result)
      rowfactor <- rowcontrol / rowtotal
      result <- sweep(as.matrix(result), 1, as.matrix(rowfactor), "*")
      if (debugger) {
        print(paste("row factor = ", rowfactor))
        print(result)
      }
      
      rowcheck <- sum(abs(1 - rowfactor))
      colcheck <- sum(abs(1 - colfactor))
      checksum <- max(rowcheck, colcheck)
      iter <- iter + 1
      #print(paste("Ending while loop, checksum > closure",checksum > closure,"iter < maxiter",iter < maxiter))
    }#End while loop
    
    out <-
      list(
        fitted.table = result,
        number.iterations = iter,
        tolerance = checksum
      )
    result <- round(result,3)
    
    res.bin <- result
    res.bin[res.bin>0] <- 1
    
    orig.bin <- apportion_output[,6:ncol(apportion_output)]
    orig.bin[orig.bin>0] <- 1
    
    if (min(res.bin-orig.bin) != 0) {
      result <- apportion_output[,6:ncol(apportion_output)]
      print("IPF produced impossible results and original apportioned data are returned")
    }
    return(cbind(temp.attr,result))
  }
