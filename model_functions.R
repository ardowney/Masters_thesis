##Run this entire script to create the necessary functions into your working environment

eem_parafac <- function(eem_list, comps, maxit = 7000, normalise = TRUE, const = c("nonneg","nonneg","nonneg"), nstart = 30, ctol = 10^-8, strictly_converging = FALSE, cores = parallel::detectCores(logical=FALSE), verbose = FALSE, output = "best",...){
  eem_array <- eem2array(eem_list)
  if(normalise){
    if(verbose) cat("EEM matrices are normalised!",fill=TRUE)
    eem_array <- eem_array %>% norm_array()
  }
  if(verbose) cat(paste0(cores," cores are used for the calculation."),fill=TRUE)
  res <- lapply(comps,function(comp){
    #comp <- 6
    if(verbose) cat(paste0("calculating ",comp," components model..."),fill=TRUE)
    cl <- NULL
    if(cores > 1){
      cl <- makeCluster(min(cores,nstart), type="PSOCK")
      clusterExport(cl, c("eem_array","comp","maxit","const","ctol","cores"), envir=environment())
      clusterEvalQ(cl, library(multiway))
    }
    if(strictly_converging){
      cpresult <- parafac_conv(eem_array, nfac = comp, const = const, maxit = maxit, parallel = (cores > 1), cl = cl, ctol = ctol, nstart = nstart, output = "all", verbose = verbose, ...)
    } else {
      cpresult <- parafac(eem_array, nfac = comp, const = const, maxit = maxit, parallel = (cores > 1), cl = cl, ctol = ctol, nstart = nstart, output = "all",...)#, ...
    }
    Rsqs <- lapply(cpresult,`[[`,"Rsq") %>% unlist()
    cpresult1 <- cpresult[[which.max(Rsqs)]]
    if(cores > 1){
      stopCluster(cl)
    }
    cflags <- lapply(cpresult,`[[`,"cflag") %>% unlist()
    converged <- sum(cflags == 0)/nstart
    if(converged <= 0.5){
      warning("Calculating the ",comp," component",ifelse(comp > 1, "s","")," model, only ",sum(cflags == 0)," out of ",nstart," models converged! You might want to increase the number of initialisations (nstart) or iterations (maxit).")
    }
    if(cpresult1$cflag != 0){
      warning("The PARAFAC model with ",comp," component",ifelse(comp > 1, "s", "")," did not converge! Increasing the number of initialisations (nstart) or iterations (maxit) might solve the problem.")
    }
    if(output == "all"){
      cpresult1$models <- lapply(cpresult,.trans_parafac, em = eem_list[[1]]$em, ex = eem_list[[1]]$ex, samples = eem_list %>% eem_names(), comp = comp, const = const, norm_factors = attr(eem_array,"norm_factors"))
    }
    cpresult1 <- .trans_parafac(cpresult1, em = eem_list[[1]]$em, ex = eem_list[[1]]$ex, samples = eem_list %>% eem_names(), comp = comp, const = const, norm_factors = attr(eem_array,"norm_factors"))
    cpresult1$converged <- converged
    cpresult1
  })
  mostattributes(res) <- attributes(eem_array)
  return(res)
}

.trans_parafac <- function(parafac, em, ex, samples, comp, const, norm_factors){
  attr(parafac,"norm_factors") <- norm_factors
  rownames(parafac$B) <- eem_list[[1]]$em
  rownames(parafac$C) <- eem_list[[1]]$ex
  rownames(parafac$A) <- eem_list %>% eem_names()
  labComp <- paste("Comp.",1:comp,sep="")
  colnames(parafac$A) <- labComp
  colnames(parafac$B) <- labComp
  colnames(parafac$C) <- labComp
  # small issue with multiway: slightly negative values are possible despite using nonnegative constraints
  non <- grepl("no|unsmpn",const)
  if(non[1]) parafac$A[parafac$A < 0] <- 0
  if(non[2]) parafac$B[parafac$B < 0] <- 0
  if(non[3]) parafac$C[parafac$C < 0] <- 0
  parafac
}

parafac_conv <- function(X, nstart, verbose = FALSE, output = c("best", "all"), cl = NULL, ...){
  nmod <- 0
  ntot <- 0
  cpresult_all <- list()
  while(nmod < nstart & ntot <= 10 * nstart){
    pred_factor <- ifelse(ntot == 0, 1, ifelse(nmod == 0, 3, ntot/nmod/2))
    if(verbose) cat("Due to previous model calculations, pred_factor was set", pred_factor, fill = TRUE)
    nmiss <- ceiling((nstart - nmod) * pred_factor / 8) * 8
    if(verbose) cat("start run with",nmiss,"models...",fill = TRUE)
    if(!is.null(cl)){
      clusterExport(cl, c("nmiss"), envir = environment())
    }
    cpresult <- parafac(X, nstart = nmiss, output = "all", cl = cl, ...) #, ...
    cpresult_conv <- cpresult[lapply(cpresult,`[[`,"cflag") %>% unlist() == 0]
    cpresult_all <- c(cpresult_all,cpresult_conv)
    nmod <- length(cpresult_all)
    ntot <- ntot + nmiss
    if(verbose) cat(length(cpresult_conv),"models converged successfully in this run!", fill = TRUE)
    if(verbose) cat(length(cpresult_all),"models calculated!", fill = TRUE)
  }
  if(verbose) cat(nmod,"out of",ntot,"models converged!",fill=TRUE)
  if(output[1] == "best"){
    sses <- cpresult_all %>% lapply(`[[`,"SSE") %>% unlist()
    cpresult_all <- cpresult_all[[which.min(sses)]]
  } else if (output[1] == "all"){
    sses <- cpresult_all %>% lapply(`[[`,"SSE") %>% unlist()
    cpresult_all <- cpresult_all[[sses[sort(order(-sses)[1:nstart])]]]
  }
  if(ntot >= 10 * nstart) warning("Maximum number of starts reached without generating the desired number of valid models.")
  cpresult_all
}

eempf_rescaleBC <- function(pfmodel,newscale = "Fmax"){
  nf <- attr(pfmodel,"norm_factors")
  comp <- ncol(pfmodel$A)
  if(newscale == "Fmax"){
    Bmax <- pfmodel$B %>% abs() %>% matrixStats::colMaxs()
    Cmax <- pfmodel$C %>% abs() %>% matrixStats::colMaxs()
    pfmodel$B <- pfmodel$B %*% diag(1/Bmax, nrow=comp) %>%
      `colnames<-`(colnames(pfmodel$B))
    pfmodel$C <- pfmodel$C %*% diag(1/Cmax, nrow=comp) %>%
      `colnames<-`(colnames(pfmodel$B))
    pfmodel$A <- pfmodel$A %*% diag(Bmax, nrow=comp) %*% diag(Cmax, nrow=comp) %>%
      `colnames<-`(colnames(pfmodel$B))
  } else {
    pfmodel <- rescale(pfmodel,mode = "C", newscale = newscale, absorb = "A")
    pfmodel <- rescale(pfmodel,mode = "B", newscale = newscale, absorb = "A")
  }
  attr(pfmodel,"norm_factors") <- nf
  labComp <- paste("Comp.",1:comp,sep="")
  colnames(pfmodel$A) <- labComp
  colnames(pfmodel$B) <- labComp
  colnames(pfmodel$C) <- labComp
  pfmodel
}

eempf_comp_names <- function(pfmodel){
  if(class(pfmodel) == "parafac") {
    colnames(pfmodel$A)
  }else if(class(pfmodel) == "list" & class(pfmodel[[1]]) == "parafac"){
    lapply(pfmodel, function(pfm) colnames(pfm$A))
  } else{
    stop("pfmodel is not a parafac model or a list of parafac models!")
  }
}

`eempf_comp_names<-` <- function(pfmodel, value){
  if(class(pfmodel) == "parafac") {
    colnames(pfmodel$A) <- value
    colnames(pfmodel$B) <- value
    colnames(pfmodel$C) <- value
    pfmodel %>% `class<-`("parafac")
  }else if(class(pfmodel) == "list" & class(pfmodel[[1]]) == "parafac"){
    if(!is.list(value) | (length(value) == 1 & length(pfmodel) > 1)) value <- lapply(1:length(pfmodel), function(x) value)
    lapply(1:length(pfmodel), function(pfn){
      colnames(pfmodel[[pfn]]$A) <- value[[pfn]][1:ncol(pfmodel[[pfn]]$A)]
      colnames(pfmodel[[pfn]]$B) <- value[[pfn]][1:ncol(pfmodel[[pfn]]$A)]
      colnames(pfmodel[[pfn]]$C) <- value[[pfn]][1:ncol(pfmodel[[pfn]]$A)]
      pfmodel[[pfn]] %>% `class<-`("parafac")
    })
  } else{
    stop("pfmodel is not a parafac model or a list of parafac models!")
  }
}

eem2array <- function(eem_list){
  eem_matrices <- eem_list %>%
    sapply("[", "x")
  dv <- lapply(eem_matrices,dim) %>% bind_cols()
  if(all(dv[1,1] %>% unlist() == dv[1,]) & all(dv[2,1] %>% unlist() == dv[2,])) dim_eem <- c(dv[,1] %>% unlist(), eem_matrices %>% length()) else dim_eem <- NA
  if(is.na(dim_eem[1])) stop("dimensions mismatch!")
  
  eem_array <- array(eem_matrices %>% unlist, dim=dim_eem)
  eem_array <- eem_array %>% aperm(perm = c(3,1,2), resize = TRUE)
  attr(eem_array,"em") <- eem_list[[1]]$em
  attr(eem_array,"ex") <- eem_list[[1]]$ex
  attr(eem_array,"samples") <- eem_list %>% sapply("[","sample") %>% unlist()
  attr(eem_array,"mdim") <- dim(eem_array)
  eem_array
}

norm_array <- function(eem_array){
  norm_factors <- lapply(1:(dim(eem_array)[1]),function(s) {
    sd(eem_array[s,,], na.rm=TRUE)
  }) %>% unlist()
  eem_array <- eem_array / norm_factors
  attr(eem_array,"norm_factors") <- norm_factors
  eem_array
}

eempf_comp_mat <- function(pfmodel,gather=TRUE){
  mat <- lapply(seq(1:ncol(pfmodel$A)), function(comp){
    m <- matrix(pfmodel$B[,comp]) %*% t(matrix(pfmodel$C[,comp])) %>%
      data.frame()
    colnames(m) <- pfmodel$C %>% rownames()
    rownames(m) <- pfmodel$B %>% rownames()
    if(gather==TRUE) m <- m %>% tibble::rownames_to_column("em") %>% gather(ex,value,-em)
    m
  })
  names(mat) <- colnames(pfmodel$B)
  mat
}

eempf_leverage <- function(pfmodel){
  cpl <- lapply(pfmodel[c("A","B","C")], function(M) diag(M %*% pinv(t(M) %*% M) %*% t(M)))
  names <- list(rownames(pfmodel$A),rownames(pfmodel$B),rownames(pfmodel$C))
  cpl <- lapply(1:(cpl %>% length()),function(i){ cpl[[i]] %>% setNames(names[[i]])}) %>%
    setNames(c("A","B","C"))
}

eempf_mleverage <- function(pfres_comps,ecdf = FALSE, stats = FALSE){
  cpls <- lapply(pfres_comps,eempf_leverage) %>%
    lapply(unlist) %>%
    lapply(function(ll) data.frame(parameter=names(ll),value=ll)) %>%
    list_join(by = "parameter") %>%
    `colnames<-`(c("parameter",paste0("comps",lapply(pfres_comps,function(cpout){ncol(cpout$A)}) %>% unlist())))
  if(ecdf){
    cpls <- cpls %>%
      mutate_if(is.numeric,function(col) ecdf(col)(col))
  }
  if(stats){
    cpls <- cpls %>%
      mutate(mean = rowMeans(select(., -parameter)), stdev=rowSds(select(., -parameter) %>% as.matrix()))
  }
  cpls
}

eempf_leverage_data <- function(cpl,qlabel=0.1){
  cpl <- cpl %>%
    lapply(. %>% data.frame() %>% rownames_to_column("x") %>% setNames(c("x","leverage")))
  mode_name <- c("sample","em","ex")
  cpl <- lapply(1:3,function(i){
    M <- cpl[[i]] %>% mutate(mode = mode_name[i])
  }) %>%
    bind_rows()
  pl <- cpl %>%
    group_by(mode) %>%
    mutate(q = quantile(leverage,1 - qlabel)) %>%
    mutate(label = ifelse(leverage > q,x,NA)) %>%
    ungroup()
}

norm2A <- function(pfmodel){
  if(class(pfmodel) != "parafac") stop("pfmodel must be an object of class parafac!")
  if(!is.null(attr(pfmodel,"norm"))){
    pfmodel$A <- pfmodel$A * attr(pfmodel,"norm_factors")
    attr(pfmodel,"norm_factors") <- NULL
  }
  pfmodel
}


eempf_cortable <- function(pfmodel,normalisation = FALSE, method="pearson",...){
  if(normalisation) pfmodel <- norm2A(pfmodel)
  pfmodel %>%
    .$A %>%
    cor(method=method,...)
}


maxlines <- function(pfmodel){
  maxl <- lapply(colnames(pfmodel$C),function(comp){
    em = (pfmodel$B[,comp]*max(pfmodel$C[,comp])) %>%
      data.frame(e="em", wavelength = as.numeric(names(.)),value=.)
    ex = (pfmodel$C[,comp]*max(pfmodel$B[,comp])) %>%
      data.frame(e="ex", wavelength = as.numeric(names(.)),value=.)
    res <- bind_rows(em,ex) %>%
      setNames(c("e","wavelength",comp))
  }) %>%
    list_join(by=c("e","wavelength"))
}


eempf_residuals <- function(pfmodel,eem_list,select=NULL, cores = parallel::detectCores(logical = FALSE)/2){
  pfmodel <- norm2A(pfmodel)
  if(!is.null(select)){
    eem_list <- eem_extract(eem_list,sample = select ,keep=TRUE,verbose = FALSE)
  }
  if(!all(eem_names(eem_list) %in% rownames(pfmodel$A)) | length(eem_list) == 0){
    pfmodel <- A_missing(eem_list,pfmodel,cores=cores)
  }
  what <- which(rownames(pfmodel$A) %in% (eem_list %>% eem_names()))
  pfmodel$A <- as.data.frame(pfmodel$A)[what,]
  res_data <- lapply(pfmodel$A %>% rownames(),function(sample){
    comps <- lapply(pfmodel$A %>% colnames(),function(component){
      pfmodel$B[,component] %*% t(pfmodel$C[,component]) * pfmodel$A[sample,component]
    })
    names(comps) <- pfmodel$C %>% colnames()
    fit <- comps %>%
      Reduce('+', .)
    eem <- eem_list[[which(eem_list %>% eem_names == sample)]]
    samp <- eem$x[eem$em %in% rownames(pfmodel$B),eem$ex %in% rownames(pfmodel$C)]
    res <- samp - fit
    
    comps <- lapply(pfmodel$A %>% colnames(),function(component){
      comps[[component]] %>% data.frame() %>% mutate(type = component, em = rownames(pfmodel$B)) %>%
        gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    }) %>%
      bind_rows()
    colnames(samp) <- rownames(pfmodel$C)
    samp <- samp %>% data.frame() %>% mutate(type = "sample", em = rownames(pfmodel$B)) %>%
      gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    rownames(res) <- rownames(pfmodel$B)
    colnames(res) <- rownames(pfmodel$C)
    res <- res %>% data.frame() %>% mutate(type = "residual", em = rownames(pfmodel$B)) %>%
      gather(ex,value,-em,-type) %>% mutate(ex = substr(ex,2,4))
    bind_rows(list(comps,samp,res)) %>%
      mutate(Sample = sample)
  }) %>%
    bind_rows()
}


A_missing <- function(eem_list,pfmodel = NULL,cores = parallel::detectCores(logical = FALSE),components = NULL, const = NULL, control = NULL, ...){
  eem_list <- eem_red2smallest(eem_list)
  if(is.null(pfmodel) & is.null(components)) stop("You must either specify a model or components as a base for the newly generated model!")
  
  exclude <- list("ex" = eem_list[[1]]$ex[!(eem_list[[1]]$ex %in% rownames(pfmodel$C))],
                  "em" = eem_list[[1]]$em[!(eem_list[[1]]$em %in% rownames(pfmodel$B))],
                  "sample" = c()
  )
  x <- eem_list %>%
    eem_exclude(exclude)
  if(!is.null(components)){
    if(!is.null(pfmodel)) warning("The base model is ignored since you provided components manually!")
    if(class(components[[1]]) == "parafac_components") components <- eempf_bindxc(components)
    if(class(components) == "parafac_components"){
      Bfixed <- components$B
      Cfixed <- components$C
      comps <- ncol(components$B)
      if(is.null(const)) const <- c("nonneg", "nonneg", "nonneg")
    } else {
      stop("The list of components you supplied is invalid!")
    }
  } else {
    Bfixed <- pfmodel$B
    Cfixed <- pfmodel$C
    comps <- pfmodel$A %>% ncol()
    normalise = (!is.null(attr(pfmodel,"norm_factors")))
    control = pfmodel$control
    const = pfmodel$const
  }
  missingAs <- eem_parafac(x,comps = comps,normalise = (!is.null(attr(pfmodel,"norm_factors"))),Bfixed = Bfixed, Cfixed = Cfixed,cores = cores,const = const, control = control, ...)
  missingAs[[1]]
}


splithalf <- function(eem_list, comps, splits = NA, rand = FALSE, normalise = TRUE, nstart = 20, cores = parallel::detectCores(logical = FALSE), maxit = 2500, ctol = 10^(-7), rescale = TRUE, verbose = FALSE, ...){
  a <- seq(1,eem_list %>% length())
  if(rand){
    a <- a %>% sample()
  }
  if(is.na(splits[1])) splits <- lapply(seq(1:4),function(sp) a[seq(sp,length(a),by=4)] %>% sort())
  names(splits) <- LETTERS[1:length(splits)]
  
  spl_eems <- lapply(combn(seq(1:length(splits)),2) %>% split(rep(1:ncol(.), each = nrow(.))), function(co){
    eem_list %>% eem_extract(splits[co] %>% unlist() %>% sort(),keep=TRUE, verbose=FALSE)
  })
  
  split_designations <- c("AB","AC","AD","BC","BD","CD")
  
  names(spl_eems) <- split_designations
  
  if(verbose){
    cat(paste0(cores," cores are used for the calculation."),fill=TRUE)
    if(normalise){
      if(verbose) cat("EEM matrices are normalised!",fill=TRUE)
    }
    cat(paste0("Calculating PARAFAC models with split-half data..."),fill=TRUE)
    pb <- txtProgressBar(max = length(spl_eems), style=3)
  }
  fits <- lapply(1:length(spl_eems),function(i){
    # i <- 1
    mod <- eem_parafac(spl_eems[[i]], comps = comps, normalise = normalise, maxit = maxit, nstart = nstart, cores = cores,ctol = ctol, verbose = FALSE)#,...
    if(rescale) mod <- lapply(mod,eempf_rescaleBC,newscale="Fmax")
    if(verbose) setTxtProgressBar(pb, i)
    mod
  }) #
  if(verbose) close(pb)
  
  sscs <- eempf_ssc(fits, tcc = TRUE)
  
  sscs2 <- sscs %>%
    lapply(lapply,ssc_max)
  
  C_sort <- sscs2 %>%
    .[grepl("1vs",names(.))] %>%
    lapply(`[[`,2) %>%
    lapply(attr, "order")
  
  fits <- lapply(1:length(spl_eems),function(sel){
    lapply(fits[[sel]],function(f){
      f$A <- f$A[,C_sort[[sel]]]
      f$B <- f$B[,C_sort[[sel]]]
      f$C <- f$C[,C_sort[[sel]]]
      f
    })
  })
  
  sel_comb <- lapply(1:(length(fits)/2), function(i){
    paste0(split_designations[i],"vs",split_designations[length(fits) + 1 - i]) %>%
      setNames(paste0(i,"vs",length(fits) + 1 - i))
  }) %>%
    unlist()
  
  attr(fits,"tcc_table") <- sscs2 %>%
    .[names(sel_comb)] %>%
    lapply(lapply,data.frame) %>%
    lapply(bind_cols) %>%
    lapply(setNames, c("tcc_ex", "tcc_em")) %>%
    lapply(mutate, component = paste0("Comp.",1:n())) %>%
    bind_rows(.id = "comb") %>%
    mutate(comb = sel_comb[comb]) %>%
    select(component,comb,tcc_em,tcc_ex) %>%
    arrange(component,comb)
  attr(fits,"splits") <- lapply(spl_eems,eem_names) %>% setNames(c("AB","AC","AD","BC","BD","CD"))
  
  fits
}




tcc_find_pairs <- function(fits){
  warning("This function is deprecated! Please use eempf_ssc and ssc_max.")
  sel <- 0
  problem <- FALSE
  table <- lapply(fits,function(fit){
    sel <<- sel + 1
    c <- fit %>% lapply(eempf_comp_mat)
    tab <- lapply(c,function(c1){
      nc1 <- length(c1)
      nc2 <- 0
      lapply(c1,function(c2){
        nc2 <<- nc2 + 1
        c2 <- c2 %>%
          mutate(comps = nc1, comp = paste0("Comp.",nc2))
      }) %>%
        bind_rows()
    }) %>%
      bind_rows() %>%
      mutate(selection = sel, ex = as.numeric(ex), em = as.numeric(em)) %>%
      group_by(comp) %>%
      filter(em == em[which.max(value)] | ex == ex[which.max(value)]) %>%
      mutate(em2 = em, em = ifelse(ex != ex[which.max(value)],NA,em),ex = ifelse(em2 != em2[which.max(value)],NA,ex)) %>%
      select(-em2) %>%
      ungroup()
  }) %>%
    bind_rows()
  
  comps <- max(table$comps)
  
  ct <- lapply(c("ex","em"),function(e){
    t <- table %>%
      mutate(o = paste0(selection,comp)) %>%
      filter_(.dots=paste0("!is.na(",e,")")) %>%
      select(one_of(!!e),value,o) %>%
      arrange_(.dots=e) %>%
      spread(o,value)
    comb_nam <- colnames(t) %>% .[. != e]
    t <- t %>%
      select(one_of(comb_nam)) %>%
      multiway::congru()
    rownames(t) <- comb_nam
    colnames(t) <- comb_nam
    attr(t,"e") <- e
    t
  }) %>%
    setNames(c("ex","em"))
  
  arr <- ct %>%
    .[["em"]] %>%
    .[1:comps,(comps+1):(length(fits)*comps)] %>%
    data.frame()
  
  a <- data.frame("comp1" = arr %>% rownames(),stringsAsFactors = FALSE)
  for(i in 2:6){
    b <- c()
    for(j in 1:(arr %>% nrow())){
      sel <- (arr %>% colnames() %>% substr(2,2) == i) & (!arr %>% colnames() %in% b)
      if(sum(sel) == 1){
        d <- arr %>% colnames() %>% .[sel]
      }else{
        d <- arr[j,which(arr %>% colnames() %>% substr(2,2) == i & !arr %>% colnames() %in% b)] %>% which.max() %>% names()
      }
      b <- c(b,d)
    }
    a <- a %>% bind_cols(data.frame(i=b,stringsAsFactors = FALSE))
  }
  
  arrange_emex <- a %>%
    mutate_all(as.character()) %>%
    gather(i,comp2,-comp1) %>%
    mutate(selection = substr(comp2,2,2), comp2 = substr(comp2,3,8))%>%
    select(comp1,selection,comp2) %>% #View()
    rowwise() %>%
    mutate(tcc_ex = ct[[1]][comp1,paste0(selection,comp2)],tcc_em = ct[[2]][comp1,paste0(selection,comp2)]) %>%
    ungroup()
  
  if(length(fits) == 6){
    tt <- arrange_emex %>%
      select(comp1) %>%
      distinct(comp1) %>%
      mutate(selection = "1", comp2 = substr(comp1,2,7)) %>%
      bind_rows(arrange_emex) %>%
      select(-tcc_em,-tcc_ex) %>%
      group_by(comp1) %>%
      spread(selection,comp2) %>%
      ungroup()
    
    split_designations <- c("AB","AC","AD","BC","BD","CD")
    
    ttt <- lapply(1:3,function(i){
      tt %>%
        select(comp1,one_of(paste0(i),paste0(7-i))) %>%
        mutate(i=i, comb = paste0(split_designations[i],"vs",split_designations[7-i])) %>%
        setNames(c("set","x","y","i","comb"))
    }) %>%
      bind_rows() %>%
      mutate(component = substr(set,2,7)) %>%
      rowwise() %>%
      mutate(tcc_ex = ct[[1]][paste0(i,x),paste0(7-i,y)],tcc_em = ct[[2]][paste0(i,x),paste0(7-i,y)]) %>%
      arrange(component) %>%
      select(component,comb,tcc_ex,tcc_em)
    attr(arrange_emex,"tcc_table") <- ttt
  }
  arrange_emex
}


splithalf_tcc <- function(fits){
  attr(fits,"tcc_table")
}


splithalf_splits <- function(fits){
  attr(fits,"splits")
}


tcc <- function(maxl_table,na.action="na.omit"){
  c <- lapply(c("em","ex"), function(E) {
    c2 <- maxl_table %>% filter(e == E) %>% select(-e) %>% arrange(wavelength)
    if(na.action == "na.omit") c2 <- c2 %>% na.omit()
    c2 <- c2 %>%
      select(-wavelength) %>%
      congru()
    na <- maxl_table %>% select(-e,-wavelength) %>% names()
    rownames(c2) <- na
    colnames(c2) <- na
    c2
  }) %>%
    setNames(c("em","ex"))
  c
}


eempf_openfluor <- function(pfmodel, file, Fmax = TRUE){
  if(!dir.exists(dirname(file.path(file)))){
    stop("The path to your file does not contain an existing directory. Please enter a correct path!")
  }
  factors <- rbind(pfmodel$C %>%
                     data.frame(mode = "Ex", wl = rownames(.),.),
                   pfmodel$B %>%
                     data.frame(mode = "Em", wl = rownames(.),.))
  template <- system.file("openfluor_template.txt",package="staRdom")
  template <- readLines(template)
  template <- stringr::str_replace_all(template,"toolbox","toolbox\tstaRdom")
  template <- stringr::str_replace_all(template,"nSample",paste0("nSample\t",pfmodel$A %>% nrow()))
  write(template,file)
  write.table(factors,file=file,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
  message("An openfluor file has been successfully written. Please fill in missing header fields manually!")
}


eempf4analysis <- function(pfmodel,eem_list = NULL, absorbance = NULL, cuvl = NULL, n = 4, export = NULL,...){
  loadings <- pfmodel %>%
    norm2A() %>%
    .$A %>%
    data.frame() %>%
    tibble::rownames_to_column("sample")
  if(!is.null(eem_list)){
    eem4peaks <- eem_list %>% eem_smooth(n=n)
    indices_peaks <- eem4peaks %>% eem_biological_index() %>%
      full_join(eem4peaks %>% eem_coble_peaks(), by="sample")  %>%
      full_join(eem4peaks %>% eem_fluorescence_index(), by="sample") %>%
      full_join(eem4peaks %>% eem_humification_index(scale=TRUE), by="sample")
    loadings <- full_join(loadings,indices_peaks,by="sample")
  }
  if(!is.null(absorbance)){
    if(is.null(cuvl)){
      warning("Because of missing cuvette length, absorbance slope parameters were not calculated!")
    } else {
      abs_parameters <- abs_parms(absorbance %>%
                                    select(one_of(names(absorbance))), cuvl)
      loadings <- full_join(loadings,abs_parameters,by="sample")
    }
  }
  if(!is.null(export)){
    write.table(loadings, file=export, row.names=FALSE, ...)
  }
  loadings
}

eempf_export<- function(pfmodel,export = NULL, Fmax = TRUE,...){
  pfmodel <- norm2A(pfmodel)
  if(Fmax) pfmodel <- eempf_rescaleBC(pfmodel)
  tabs <- list(pfmodel$A %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("sample") #%>%
               #`colnames<-`(colnames(.) %>% stringr::str_replace_all("Comp.","Fmax"))
               ,
               pfmodel$B %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("Em"),
               pfmodel$C %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("Ex"))
  rows <- lapply(tabs,nrow) %>% unlist() %>% max()
  tabs <- lapply(tabs,function(tab){
    if(nrow(tab) < rows){
      tab <- rbind(tab,as.data.frame(matrix(nrow = rows - nrow(tab), ncol = ncol(tab))) %>% `colnames<-`(colnames(tab)))
    }
    tab <- cbind(tab,as.data.frame(matrix(nrow = rows, ncol = 1)) %>% `colnames<-`(" "))
  })
  tabs <- do.call(cbind,tabs)
  if(!is.null(export)) write.table(tabs,file=export,row.names = FALSE,na="",...)
  tabs %>% invisible()
}


eempf_corcondia <- function(pfmodel,eem_list,divisor="core"){
  arr <- eem_list %>% eem2array()
  if(!is.null(attr(eem_list,"norm_factors"))) arr <- arr %>% norm_array()
  corcondia(arr,pfmodel,divisor=divisor)
}

eempf_eemqual <- function(pfmodel,eem_list,splithalf = NULL, ...){
  corec <- eempf_corcondia(pfmodel,eem_list)
  fit <- pfmodel$Rsq
  if(is.null(splithalf)) sh <- splithalf(eem_list,comps = pfmodel$A %>% ncol(),...) else sh <- splithalf
  splth <- splithalf_tcc(sh) %>%
    group_by(component) %>%
    summarise_at(vars(contains("tcc")),max,na.rm=TRUE) %>%
    .[c(2,3)] %>%
    prod()
  tab <- data.frame(components = ncol(pfmodel$A) ,fit = fit, corec = corec, splithalf = splth, eemqual = fit*corec*splth)
  attr(tab,"shmodel") <- sh
  tab
}


eempf_varimp <- function(pfmodel, eem_list, cores = parallel::detectCores(logical=FALSE),...){
  
  exclude <- list("ex" = eem_list[[1]]$ex[!(eem_list[[1]]$ex %in% rownames(pfmodel$C))],
                  "em" = eem_list[[1]]$em[!(eem_list[[1]]$em %in% rownames(pfmodel$B))],
                  "sample" = c()
  )
  x <- eem_list %>%
    eem_red2smallest() %>%
    eem_exclude(exclude)
  mods <- lapply(1:(pfmodel$A %>% ncol()), function(c){
    eem_parafac(x, comps = ncol(pfmodel$A)-1,normalise = (!is.null(attr(pfmodel,"norm_factors"))),Bfixed = pfmodel$B[,-c], Cfixed = pfmodel$C[,-c],cores = cores,control = pfmodel$control, const = pfmodel$const, ...)
  })
  Rsq_al <- pfmodel$Rsq - (mods %>% lapply(function(mod) {mod[[1]]$Rsq}) %>% unlist())
}




eempf_reorder <- function(pfmodel,order,decreasing = FALSE){
  if(!(order[1] == "em" | order[1] == "ex" | is.vector(order))) stop("no valid data suppli ed for order!")
  if(order[1] == "em") order <- apply(pfmodel$B,2,which.max) %>% sort.list(decreasing = decreasing)
  if(order[1] == "ex") order <- apply(pfmodel$C,2,which.max) %>% sort.list(decreasing = decreasing)
  if(ncol(pfmodel$A) != length(order)) stop("the length of the order vector does not fit the number of components")
  
  mod <- try(reorder.parafac(pfmodel,neworder = order), silent=TRUE)
  if(class(mod) =="try-error") stop(mod) else {
    mod
  }
}


eempf_excomp <- function(pfmodel,comps){
  list(pfmodel$B[,comps], pfmodel$C[,comps]) %>%
    `names<-`(c("B","C")) %>%
    `class<-`("parafac_components")
}


empf_ssc <- function(pfmodels, tcc = FALSE, m = FALSE, cores = parallel::detectCores(logical = FALSE)){
  classes <- unlist(lapply(unlist(pfmodels, recursive = FALSE),class))
  if(any(classes == "parafac") & !is.null(classes)){ ## Results from splithalf
    pfmodels %>%
      unlist(recursive = FALSE) %>%
      lapply(function(mod){
        list(B=mod$B,C=mod$C)
      }) %>%
      eempf_ssc(tcc = tcc, m = m, cores = cores)
  } else if(all(unlist(lapply(pfmodels,class)) == "parafac") & !is.null(unlist(lapply(pfmodels,class)))){ ## PARAFAC models
    pfmodels %>% lapply(function(mod){
      list(B=mod$B,C=mod$C)
    }) %>%
      eempf_ssc(pfmodels = ., tcc = tcc, m = m, cores = cores)
  } else if(all(classes == "matrix")){ ## matrices
    
    cl <- makePSOCKcluster(min(cores,length(pfmodels)))
    clusterExport(cl, c("pfmodels","tcc"), envir = environment())
    clusterEvalQ(cl,require(staRdom))
    
    SSCs <- parLapply(cl,1:length(pfmodels),function(k){
      lapply(k:length(pfmodels), function(l){
        B = ssc(pfmodels[[k]][[1]], pfmodels[[l]][[1]], tcc = tcc)
        C = ssc(pfmodels[[k]][[2]], pfmodels[[l]][[2]], tcc = tcc)
        list(B=B, C=C)
      }) %>%
        setNames(paste0(k,"vs",k:length(pfmodels)))
    }) %>% unlist(recursive = FALSE)
    stopCluster(cl)
    
    if(m){
      SSCs <- lapply(SSCs, function(mats){
        sqrt(mats[[1]]*mats[[2]])
      })
    }
    SSCs
  } else {
    stop("No suitable data supplied! Please refer to the eempf_ssc help.")
  }
}


ssc <- function(mat1, mat2, tcc = FALSE){
  if(any(is.null(mat1),is.na(mat1),is.null(mat2), is.na(mat2))){
    a <- NA
  } else {
    a <- lapply(1:ncol(mat1),function(nc){
      col1 <- mat1[,nc]
      apply(mat2,2,function(col2){
        tcc_cal <- sum(col1*col2)/sqrt(sum(col1^2)*sum(col2^2))
        if(!tcc){
          wl <- as.numeric(names(col1))
          if(any(is.na(wl)) | pracma::isempty(wl)){
            stop("SSCs cannot be calculated. Please add wavelengths as rownames of the matrices!")
          }
          alpha <- abs((wl[which.max(col1)]-wl[which.max(col2)]) / diff(range(wl)))
          beta <- abs((sum(col1/max(col1)) - sum(col2/max(col2))) / diff(range(wl)))
          ssc <- tcc_cal -alpha - beta
        } else {
          tcc_cal
        }
      })
    }) %>% setNames(colnames(mat1)) %>%
      do.call(rbind,.)
  }
  attr(a,"method") <- ifelse(tcc, "TCC", "SSC")
  a
}


eempf_ssccheck <- function(pfmodels, best = length(pfmodels), tcc = FALSE, cores = parallel::detectCores(logical = FALSE)){
  #pfmodels <- pf3
  Rsqs <- lapply(pfmodels,`[[`,"Rsq") %>% unlist()
  checkmods <- pfmodels[order(Rsqs, decreasing = TRUE)[1:best]]
  not_conv <- checkmods %>%
    lapply(`[[`,"cflag") %>%
    unlist() %>%
    sapply(identical, 0) %>%
    sapply(`!`) %>%
    sum()
  if(not_conv){
    warning(paste0(not_conv," of the best ", best, " chosen models ",ifelse(not_conv == 1, "is","are")," not converging!"))
  }
  #Rsqs <- lapply(checkmods,`[[`,"Rsq") %>% unlist()
  sscs <- eempf_ssc(pfmodels = checkmods, tcc = tcc, cores = cores)
  
  cl <- makePSOCKcluster(min(cores, length(sscs)))
  clusterExport(cl, c("sscs"), envir=environment())
  clusterEvalQ(cl,require(staRdom))
  
  maxs <- parLapply(cl, sscs, lapply, ssc_max)
  
  stopCluster(cl)
  
  a <- names(maxs) %>%
    lapply(function(na){
      str_extract(na,"^[0-9]{1,2}") != str_extract(na,"[0-9]{1,2}$")
    }) %>%
    unlist() %>%
    maxs[.] %>%
    lapply(bind_rows)
  ssccheck <- a %>%
    lapply(mutate, comp = 1:n()) %>%
    bind_rows() %>%
    mutate(comparison = names(a)[cumsum(comp == 1)])
  attr(ssccheck,"method") <- ifelse(tcc, "TCC", "SSC")
  ssccheck
}

ssc_max <- function(mat){
  n <- min(dim(mat))
  p <- permutations(n , n)
  combinations <- lapply(1:nrow(p),function(row){
    per <- p[row,]
    matrix(c(1:n,per), ncol = n,byrow = TRUE)
  })
  
  best_comb <- lapply(combinations, function(c){
    pair <- c[,2] %>% unlist()
    apply(c,2,function(pair){
      mat[pair[1],pair[2]] %>% `^`(2)
    }) %>%
      sum() %>%
      sqrt()
  }) %>%
    which.max() %>%
    combinations[[.]]
  
  res <- apply(best_comb,2,function(pair){
    mat[pair[1],pair[2]]
  })
  attr(res,"order") <- best_comb[2,]
  res
}


eempf_convergence <- function(pfmodel, print = TRUE){
  if(!is.list(pfmodel$models)){
    stop("The supplied PARAFAC model does not contain the whole model set used in the calculation! Please rerun eem_parafac setting output = 'all'")
  } else {
    n <- length(pfmodel$models)
    sses <- lapply(pfmodel$models,`[[`,"SSE") %>% unlist()
    conv <- lapply(pfmodel$models,`[[`,"cflag") %>% unlist()
    res <- list(models = n, converging = sum(conv == 0), nc_itlim = sum(conv == 1), nc_other = sum(conv == 2), conv = conv, sses = sses)
    if(print){
      cat("Calculated models: ", n, fill = TRUE)
      cat("Converging models: ", sum(conv == 0), fill = TRUE)
      cat("Not converging Models, iteration limit reached: ", sum(conv == 1), fill = TRUE)
      cat("Not converging models, other reasons: ", sum(conv == 2), fill = TRUE)
      cat("Best SSE: ", min(sses), fill = TRUE)
      cat("Summary of SSEs of converging models:",fill = TRUE)
      summary(sses[conv == 0]) %>%
        print()
    }
    invisible(res)
  }
}


#####plot fubctions
eempf_compare <- function(pfres,...){
  #pfres <- pf4
  p1 <- eempf_fits(pfres,...)
  p2 <- eempf_plot_comps(pfres,type=1,...)
  p3 <- eempf_plot_comps(pfres,type=2,...)
  p1 %>% print()
  p2 %>% print()
  p3 %>% print()
  return(list(p1,p2,p3)) %>% invisible()
}



eempf_fits <- function(pfres,...){
  #pfres <- pf1
  if(is.null(names(pfres))) names <- paste0("model",1:length(pfres)) else names <- names(pfres) #paste0("model",1:length(pfres))
  pl <- data.frame(comps=lapply(pfres,"[[","A") %>% lapply(ncol) %>% unlist(),
                   fit=lapply(pfres,"[[","Rsq") %>% unlist(), mod_name = names) %>%
    rowwise() %>%
    mutate(comps = ifelse(is.na(mod_name),paste0(comps, " comps"),paste0(mod_name," (",comps, " comps)"))) %>%
    ggplot(aes(x=comps,y=fit),...)+
    labs(x="model", y="model fit (Rsq)")+
    geom_point(size=5,shape=4)
  pl
}


eempf_plot_comps <- function(pfres,type=1,names=TRUE,contour = FALSE,...){
  #pfres <- pf3
  c <- pfres %>% lapply(eempf_comp_mat)
  if(is.null(names(c))) names(c) <- paste0("model",seq(1:length(c)))
  tab <- lapply(1:length(c),function(n){
    c1 <- c[[n]]
    mod_name <- names(c)[n]
    nc1 <- length(c1)
    nc2 <- 0
    lapply(c1,function(c2){
      nc2 <<- nc2 + 1
      c2 <- c2 %>%
        mutate(comps = nc1, comp = paste0("Comp.",nc2), modname = mod_name)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    mutate(modname = factor(modname,levels = names(c))) %>%
    mutate_at(vars(ex,em,value),as.numeric)
  fill_max <- tab$value %>% max(na.rm=TRUE)
  vals <- seq(from=0,to=fill_max,length.out = 51)
  vals <- (vals - min(vals))/diff(range(vals))
  if(type==2){
    plot <- tab %>%
      group_by(modname,comp) %>%
      mutate(max_pos = which.max(value), max_em = em[max_pos],max_ex = ex[max_pos]) %>%
      mutate(exn = ifelse(em == max_em,ex,NA),emn = ifelse(ex == max_ex,em,NA)) %>%
      filter(!is.na(emn) | !is.na(exn)) %>%
      ungroup() %>%
      ggplot()+
      geom_line(aes(x=exn,y=value),colour="lightblue",group="excitation", na.rm=TRUE)+
      geom_line(aes(x=emn,y=value),colour="darkblue",group="emission", na.rm=TRUE)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      facet_grid(comp ~ modname)+
      labs(x="wavelength (nm)")
  } else {
    diffs <- tab %>%
      select(-value,-comps) %>%
      gather("spec","wl", -comp, -modname) %>%
      group_by(comp,modname,spec) %>%
      unique() %>%
      summarise(slits = diff(wl) %>% n_distinct()) %>% #View()
      .$slits != 1
    
    plot <- tab %>%
      ggplot(aes(x = ex, y = em, z = value))
    
    if(any(diffs)){
      plot <- plot +
        #geom_raster(aes(fill = value), interpolate = interpolate)
        layer(mapping = aes(colour = value, fill = value),
              geom = "tile", stat = "identity", position = "identity")
    } else {
      plot <- plot +
        layer(mapping = aes(fill = value),
              geom = "raster", stat = "identity", position = "identity")
    }
    
    plot <- plot +
      scale_fill_gradientn(colours=rainbow(75)[51:1],values=vals,limits = c(tab$value %>% min(na.rm=TRUE),fill_max), aesthetics = c("fill", "colour"))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      labs(x = "Excitation (nm)", y = "Emission (nm)") +
      facet_grid(comp ~ modname)
    if(contour){
      plot <- plot +
        geom_contour(colour = "black", size = 0.3)
    }
  }
  #print(pl)
  plot
}


eempf_leverage_plot <- function(cpl, qlabel = 0.1){
  cpl <- eempf_leverage_data(cpl,qlabel=qlabel) %>%
    ungroup() %>%
    mutate(mode = str_replace(mode,"em","Emission (nm)") %>% str_replace("ex","Excitation (nm)") %>% str_replace("sample","Sample"))
  breaks <- seq(0,1000,50)
  breaks2 <- cpl$x[cpl$mode == "Emission (nm)"] %>%
    as.numeric() %>%
    range(na.rm = TRUE)
  breaks2 <- breaks[breaks >= breaks2[1] & breaks <= breaks2[2]]
  vals2 <- cpl$x[cpl$mode == "Emission (nm)"][findInterval(breaks2, cpl$x[cpl$mode =="Emission (nm)"])]
  breaks1 <- cpl$x[cpl$mode == "Excitation (nm)"] %>%
    as.numeric() %>%
    range(na.rm = TRUE)
  breaks1 <- breaks[breaks >= breaks1[1] & breaks <= breaks1[2]]
  vals <- cpl$x[cpl$mode == "Excitation (nm)"][findInterval(breaks1, cpl$x[cpl$mode =="Excitation (nm)"])]
  pl <- cpl %>%
    ggplot(aes(x=x,y=leverage))+
    geom_point(alpha = 0.4)+
    geom_text(aes(label=label),vjust="inward",hjust="inward", na.rm=TRUE, check_overlap = TRUE)+
    scale_x_discrete(labels = c(breaks1,breaks2), breaks = c(vals, vals2)) +
    labs(x="Variables (wavelengths or samples)", y = "Leverage") +
    facet_wrap( ~ mode, scales = "free")
  pl
}

eempf_leverage_ident <- function(cpl,qlabel=0.1){
  pl <- eempf_leverage_data(cpl,qlabel=qlabel) %>%
    mutate(label = ifelse(is.na(label),"",label))
  exclude <- lapply(pl$mode %>% unique(),function(mod){
    data <- pl %>% filter(mode == mod)
    plot(data$x %>% factor(), data$leverage, xlab = mod, ylab = "leverage")
    text(data$x %>% factor(), data$leverage, label = data$label)
    ide <- identify(data$x %>% factor(), data$leverage) #, labels = data$label
    ide <- data$x %>% factor() %>% .[ide] %>% as.character()
  }) %>%
    setNames(pl$mode %>% unique())
}



eempf_comp_load_plot <- function(pfmodel,...){
  pl1 <- ggeem(pfmodel,...)
  pl2 <- eempf_load_plot(pfmodel)
  #pl1 %>% print()
  #pl2 %>% print()
  list(pl1,pl2)
}


eempf_load_plot <- function(pfmodel){
  pfmodel <- norm2A(pfmodel)
  (pfmodel$A) %>%
    data.frame() %>%
    rownames_to_column("sample") %>%
    #mutate(sample = names[[3]]) %>%
    gather(comp,amount,-sample) %>%
    ggplot()+
    geom_bar(aes(x=sample,y=amount,fill=comp),stat="identity",width=0.8)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}



eempf_comps3D <- function(pfmodel,which=NULL){
  data <- pfmodel %>% eempf_comp_mat()
  z <- lapply(data,function(mat){
    mat %>%
      data.frame() %>%
      spread(em,value) %>%
      remove_rownames() %>%
      column_to_rownames("ex") %>%
      as.matrix()
  })
  ex <- lapply(data,function(mat){
    mat$ex %>% unique() %>% as.numeric() %>% as.vector()
  })
  em <- lapply(data,function(mat){
    mat$em %>% unique() %>% as.numeric() %>% as.vector()
  })
  scene <- list(xaxis=list(title="em"),
                yaxis=list(title="ex"))
  lapply(1:length(ex),function(comp){
    if(is.null(which) | comp %in% which){
      plotly::plot_ly(x=em[[comp]],y=ex[[comp]], z=z[[comp]],colors = rainbow(12)[9:1]) %>%
        plotly::layout(scene=scene) %>%
        plotly::add_surface()
    }
  })
}


eempf_corplot <- function(pfmodel,normalisation=FALSE,lower=list(continuous="smooth"),mapping=aes(alpha=0.2),...){
  if(normalisation) pfmodel <- norm2A(pfmodel)
  pfmodel %>%
    .$A %>%
    data.frame() %>%
    ggpairs(lower=lower,mapping=mapping,...)
}


eempf_residuals_plot <- function(pfmodel,eem_list,res_data = NULL, spp = 5, select=NULL, residuals_only = FALSE , cores = parallel::detectCores(logical = FALSE), contour = FALSE){
  #pfmodel,eem_list,select=eem_names(eem_list)[10:19]
  #pfmodel <- pf4[[1]]
  #eem_list <- eem_ex
  if(is.null(res_data)){
    res_data <- eempf_residuals(pfmodel,eem_list,select=select,cores = cores)
  }
  if (!is.null(select)){
    res_data <- res_data %>% filter(Sample %in% select)
  }
  res_data <- res_data %>%
    mutate_at(vars(ex,em,value),as.numeric)
  #if(!is.numeric(fill_max)){
  if(residuals_only){
    res_data <- res_data %>%
      filter(type == "residual")
  }
  fill_max <- res_data$value %>% max(na.rm=TRUE)
  #}
  vals <- c(res_data$value %>% min(na.rm=TRUE),seq(from=0,to=fill_max,length.out = 50))
  vals <- (vals - min(vals))/diff(range(vals))
  ppp <- res_data$Sample %>% unique() %>% length() /spp
  ov_plot <- lapply(1:ceiling(ppp),function(pos){
    #pos <- 1
    pl <- res_data %>%
      filter(Sample %in% (res_data$Sample %>% unique() %>% .[(spp*(pos-1)+1):(spp*pos)])) %>%
      ggplot(aes(x=ex,y=em,z=value))
    
    diffs <- res_data %>%
      select(-value, -type) %>%
      gather("spec","wl", -Sample) %>%
      group_by(Sample,spec) %>%
      unique() %>%
      #arrange(sample,spec,wl) %>%
      #mutate(diffs = wl - lag(wl))
      summarise(slits = diff(wl) %>% n_distinct()) %>%
      .$slits != 1
    
    if(any(diffs)){
      pl <- pl +
        #geom_raster(aes(fill = value), interpolate = interpolate)
        layer(mapping = aes(colour = value, fill = value),
              geom = "tile", stat = "identity", position = "identity") +
        scale_colour_gradientn(colours=c(rainbow(70)[62],rainbow(70)[50:1]),values=vals,limits = c(res_data$value %>% min(na.rm=TRUE),fill_max))
    } else {
      pl <- pl +
        layer(mapping = aes(fill = value),
              geom = "raster", stat = "identity", position = "identity")
    }
    pl <- pl +
      #geom_raster(aes(fill=value))+ #,interpolate=TRUE
      scale_fill_gradientn(colours=c(rainbow(70)[62],rainbow(70)[50:1]),values=vals,limits = c(res_data$value %>% min(na.rm=TRUE),fill_max))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      labs(x = "Excitation (nm)", y = "Emission (nm)", fill = "fluorescence")
    if(contour){
      pl <- pl +
        geom_contour(colour = "black", size = 0.3)
    }
    if(residuals_only){
      pl <- pl +
        facet_wrap(~ Sample)
    } else {
      pl <- pl +
        facet_grid(type ~ Sample)
    }
    pl
  })
  ov_plot
}


splithalf_plot <- function(fits){
  sel <- 0
  table <- lapply(fits,function(fit){
    sel <<- sel + 1
    c <- fit %>% lapply(eempf_comp_mat)
    tab <- lapply(c,function(c1){
      nc1 <- length(c1)
      nc2 <- 0
      lapply(c1,function(c2){
        nc2 <<- nc2 + 1
        c2 <- c2 %>%
          mutate(comps = nc1, comp = paste0("Comp.",nc2))
      }) %>%
        bind_rows()
    }) %>%
      bind_rows() %>%
      mutate(selection = sel) %>%
      group_by(comps,comp) %>%
      mutate(max_pos = which.max(value), max_em = em[max_pos],max_ex = ex[max_pos]) %>%
      mutate(exn = ifelse(em == max_em,ex,NA),emn = ifelse(ex == max_ex,em,NA)) %>%
      filter(!is.na(emn) | !is.na(exn)) %>%
      mutate(ex=exn,em=emn) %>%
      select(-exn,-emn,-max_pos,-max_em,-max_ex) %>%
      ungroup() %>%
      mutate_at(vars(em,ex,value),as.numeric)
  }) %>%
    bind_rows()
  pl1 <- table %>%
    mutate(selection = factor(selection,ordered=FALSE)) %>%
    ggplot()+
    geom_line(data = . %>% filter(!is.na(ex)),aes(x=ex,y=value,colour=selection,group=selection),linetype=2)+
    geom_line(data = . %>% filter(!is.na(em)),aes(x=em,y=value,colour=selection,group=selection),linetype=1)+
    labs(x="Wavelength (nm)",y="Loading") +
    theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_grid(. ~ comp)
  pl1 %>% print()
}



eempf_report <- function(pfmodel, export, eem_list = NULL, absorbance = NULL, meta = NULL, metacolumns = NULL, splithalf = FALSE, shmodel = NULL, performance = FALSE, residuals = FALSE, spp = 5, ...){
  #shmodel <- sh
  # rm(shmodel)
  #splithalf = TRUE
  rmdfile <- system.file("PARAFAC_report.Rmd",package="staRdom")
  dir <- dirname(export)
  file <- basename(export)
  if(splithalf | performance){
    if(!is.null(eem_list) & is.null(shmodel)){
      splithalf <- splithalf(eem_list, comps = ncol(pfmodel$A), normalise = !is.null(attr(pfmodel,"norm_factors")),...)
      tcc <- splithalf_tcc(shmodel)
    } else if (!is.null(shmodel)){
      splithalf <- shmodel
      tcc <- splithalf_tcc(shmodel)
    } else {
      tcc <- NULL
      warning("Split-half analysis and/or model performance could not be incorporated due to missing EEM data or an already calculated split-half analysis.",fill=TRUE)
    }
  }
  if(performance){
    if(!is.null(eem_list)){
      performance <- eempf_eemqual(pfmodel,eem_list,splithalf)
    } else{
      warning("For a performance calculation, EEM data is needed!")
    }
  }
  imgwidth <- nrow(pfmodel$A)/8
  rmarkdown::render(rmdfile, output_file = file, output_dir = dir, params = list(pfmodel = pfmodel, eem_list = eem_list, absorbance = absorbance, meta = meta, tcc = tcc, metacolumns = metacolumns, splithalf = splithalf, performance = performance, residuals = residuals, spp = spp, imgwidth = imgwidth))
  TRUE
}


eempf_plot_ssccheck <- function(ssccheck){
  ssccheck %>%
    mutate(excitation = B,  emission = C) %>%
    select(-B,-C) %>%
    gather("spectrum", "TCC", excitation, emission) %>%
    group_by(comp, spectrum) %>%
    mutate(mean = mean(TCC), min = min(TCC), max = max(TCC)) %>%
    ggplot()+
    geom_point(aes(x=comp + 0.1 * ifelse(spectrum == "excitation",1,-1),y=mean, colour = comp, group = comp), shape = 21, size = 3)+
    geom_point(aes(x=comp + 0.1 * ifelse(spectrum == "excitation",1,-1),y=TCC, colour = comp, group = comp), alpha = 0.4)+
    geom_errorbar(aes(x=comp + 0.1 * ifelse(spectrum == "excitation",1,-1),ymin = min, ymax = max, colour = comp, group = comp, linetype = spectrum))+
    scale_color_viridis_c(guide = FALSE)+
    labs(x = "Component", y = attr(ssccheck,"method"), linetype = "", shape = "")+
    scale_x_continuous(breaks = c(1:max(ssccheck$comp)))
}




