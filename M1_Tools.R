# Toolbox for M1 Beta version
# (c) Zijian HUANG 2023

ccm.nullCreate = function(source, v, 
                          # pass to SurrogateData function
                          method = c('seasonal', 'ebisuzaki'),
                          # length of seasonal trend to be maintain
                          Tp = 12,
                          # how many models to create
                          rep = 1000,
                          # What degree of noise to add (0 = none)
                          noise = 3)
  # This function will call `SurrogateData` from CCM package
  # and create data of null model based on given method
  # Output: a list where data of each null model is an element
{
  nullSurr = SurrogateData(source[, v],
                           method = method,
                           T_period = Tp,
                           num_surr = rep, 
                           alpha = noise)
  nullList = as.list(as.data.frame(nullSurr))
  names(nullList) = c(paste('Mod', as.character(1:rep), sep = ''))
  return(nullList)
}

ccm.perform = function(tr, lt, 
                       lib = '36 236 3', 
                       variable = NULL, 
                       nullMod = c('both', 'seas', 'ebi'), 
                       nullModRep = 500)
  # Function to perform CCM (version 3.0)
  # Required: baseLookUp, edm_dmTbl, climScale, edmE
  # @param tr: trap number
  # @param lt: litter type
  # @param lib: step of library size to CCM
  # @param variable: variables to perform CCM. If NULL select all variables in table
  # @param nullMod: null models to perform, seas - seasonal, ebi - Ebisuzaki
  # @param nullModRep: number of surrogate data to create in null models
{
  # Identify species
  sp = baseLookUp$species[which(baseLookUp$trapID == tr)]
  # Table input to CCM
  baseTbl = edm_dmTbl[[as.character(tr)]][[lt]] %>%
    left_join(climScale, by = c('time' = 'Date'))
  # Variables to perform CCM
  if (is.null(variable)) {
    varList = names(climScale)[-1]
  } else {
    varList = variable
  }
  # Optimal E
  optE = edmE[[sp]][[lt]][['E']]
  # ---
  # Null models
  ## Surrogate data
  seasonalNullModels = lapply(varList, function(Var){
    ccm.nullCreate(source = baseTbl,
                   v = Var,
                   method = 'seasonal',
                   rep = nullModRep)
  }); names(seasonalNullModels) = varList
  ebiNullModels = lapply(varList, function(Var){
    ccm.nullCreate(source = baseTbl,
                   v = Var,
                   method = 'ebisuzaki',
                   rep = nullModRep)
  }); names(ebiNullModels) = varList
  nullModelList = list('seas' = seasonalNullModels,
                       'ebi' = ebiNullModels)
  ## Perform CCM on surrogate data
  message('>>> Performing CCM on Null Models')
  if (nullMod == 'both') {
    nullMode = c('seas', 'ebi')
  } else {nullMode = nullMod}
  nullCCM = lapply(nullMode, function(nmd){
    # For each type of null model
    byVar = pbapply::pblapply(varList, function(Var){
      surList = nullModelList[[nmd]][[Var]]
      # Perform CCM on each surrogate data set (i.e., each null model)
      nullPerf = lapply(surList, function(nm){
        # Create table input to CCM
        nullTbl = data.frame(time = baseTbl$time,
                             DM = baseTbl$DM,
                             env = nm)
        nullPerform = CCM(dataFrame = nullTbl,
                          columns = 'env',
                          target = 'DM',
                          libSizes = lib,
                          E = optE,
                          sample = 100)
        result = nullPerform[, 'DM:env']
        return(result)
      }) %>% do.call(cbind.data.frame, args = .)
      return(nullPerf)
    }); names(byVar) = varList
    return(byVar)
  }); names(nullCCM) = nullMode
  message('>>> Null model CCM performance success.')
  ## Rearrange and convert Null Model CCM performance
  nullCCMRes = lapply(nullMode, function(nm){
    # Retrieve time steps in CCM
    timeStepsInfo = strsplit(lib, split = ' ', fixed = T)[[1]]
    timeSteps = eval(parse(text = paste('seq(', 
                                        paste(timeStepsInfo, collapse = ','),
                                        ')', sep = '')))
    # Retrieve performance of null models
    rhoList = lapply(nullCCM[[nm]], function(x){
      resTbl = do.call(rbind.data.frame, x)
      names(resTbl) = timeSteps
      return(resTbl)
    })
    # Convert null model performance to long table
    rhoTbl = rhoList %>% do.call(rbind.data.frame, args = .) %>%
      mutate(Var = rep(varList, each = nullModRep)) %>%
      select(Var, everything()) %>%
      tidyr::gather(key = 'Size', value = 'Rho',
                    2:(length(timeSteps)+1)) %>%
      mutate(Size = as.numeric(Size))
    # Calculate mean for null model performance
    rhoMean = rhoTbl %>% group_by(Var, Size) %>% 
      summarise(nullRho = mean(Rho)) %>% as.data.frame()
    # Calculate 95% CI for null models
    rhoCI = lapply(rhoList, function(x){
      apply(x, 2, quantile, probs = c(.025, .975)) %>%
        as.data.frame()
    }) %>% do.call(rbind.data.frame, args = .) %>%
      # Add marks for variables
      mutate(Var = rep(varList, each = 2),
             Boundary = rep(c('lower', 'upper'), length(varList))) %>%
      select(Var, Boundary, everything()) %>%
      # Convert to long table
      tidyr::gather(key = 'Size', value = 'Rho',
                    3:(length(timeSteps)+2)) %>%
      # Spread CI boundaries to two columns
      tidyr::spread(Boundary, Rho) %>%
      mutate(Size = as.numeric(Size)) %>% 
      # Join Mean
      left_join(rhoMean, by = c('Var', 'Size'))
    ## Save to output
    output = list('data' = rhoTbl,
                  'CI' = rhoCI)
    return(output)
  }); names(nullCCMRes) = nullMode
  # ---
  # CCM of abiotic factors and litter biomass
  message('>>> Performing CCM...')
  factorPerf = lapply(varList, function(Var){
    ccmPerformance = CCM(dataFrame = baseTbl %>% select(time, DM, all_of(Var)),
                         columns = Var,
                         target = 'DM',
                         libSizes = lib,
                         E = optE,
                         sample = 100)
    ccmResTbl = ccmPerformance %>%
      select(Size = LibSize,
             Env2DM = paste(Var, ':DM', sep = ''),
             DM2Env = paste('DM:', Var, sep = '')) %>%
      tidyr::gather(key = 'Direction',
                    value = 'Rho',
                    Env2DM:DM2Env) %>%
      mutate(Var = Var)
    return(ccmResTbl)
  }) %>% do.call(rbind.data.frame, args = .) %>%
    mutate(species = sp, litter = lt, trap = tr)
  message('>>> CCM process finished.')
  # ---
  # Performance arrangement
  ## Join results
  if (length(nullMode) == 2) {
    fullResTbl = factorPerf %>%
      subset(Direction == 'DM2Env') %>%
      select(Var, Size, Rho) %>%
      left_join(nullCCMRes[['seas']][['CI']], by = c('Size', 'Var')) %>%
      select(Var, Size, Rho,
             seasRho = nullRho, seasLower = lower, seasUpper = upper) %>%
      left_join(nullCCMRes[['ebi']][['CI']], by = c('Size', 'Var')) %>%
      select(Var, Size, Rho, seasRho, seasLower, seasUpper, 
             ebiRho = nullRho, ebiLower = lower, ebiUpper = upper)
  } else {
    fullResTbl = factorPerf %>%
      subset(Direction == 'DM2Env') %>%
      select(Var, Size, Rho) %>%
      left_join(nullCCMRes[[nullMode]][['CI']], by = c('Size', 'Var')) %>%
      select(Var, Size, Rho, nullRho, lower, upper)
    names(fullResTbl) = c('Var', 'Size', 'Rho', 
                          paste(nullMode, 'Rho', sep = ''),
                          paste(nullMode, 'Lower', sep = ''),
                          paste(nullMode, 'Upper', sep = ''))
  }
  ## Significance test
  sigTest = lapply(nullMode, function(nm){
    byRow = lapply(1:dim(fullResTbl)[1], function(rowID){
      targetRow = fullResTbl[rowID,]
      nullDataBase = nullCCMRes[[nm]][['data']]
      nullData = subset(nullDataBase, Var == targetRow[['Var']] &
                          Size == targetRow[['Size']])$Rho
      p = 1 - ecdf(nullData)(targetRow[['Rho']])
      return(p)
    }) %>% as.numeric()
    return(byRow)
  }); names(sigTest) = nullMode
  if (length(nullMode) == 2) {
    fullResTbl = fullResTbl %>% 
      mutate(pSEAS = sigTest[['seas']],
             pEBI = sigTest[['ebi']],
             sigSEAS = ifelse(pSEAS >= 0.05, 'NS',
                              ifelse(pSEAS >= 0.01, '*',
                                     ifelse(pSEAS >= 0.001, '**', '***'))),
             sigEBI = ifelse(pEBI >= 0.05, 'NS',
                             ifelse(pEBI >= 0.01, '*',
                                    ifelse(pEBI >= 0.001, '**', '***'))),
             species = sp,
             litter = lt, 
             trap = tr) %>%
      select(species, trap, litter, Var, Size, Rho, 
             seasRho, seasLower, seasUpper, pSEAS, sigSEAS,
             ebiRho, ebiLower, ebiUpper, pEBI, sigEBI)
  } else {
    fullResTbl = fullResTbl %>% 
      mutate(p = sigTest[[nullMode]],
             sig = ifelse(p >= 0.05, 'NS',
                          ifelse(p >= 0.01, '*',
                                 ifelse(p >= 0.001, '**', '***'))),
             species = sp,
             litter = lt, 
             trap = tr) %>%
      select(species, trap, litter, everything())
  }
  ## Monotonic test
  kendallTbl = fullResTbl %>% nest_by(Var) %>%
    mutate(kentau = list(cor.test(data$Size, data$Rho, method = 'kendall'))) %>%
    summarise(tau = kentau$estimate,
              p = kentau$p.value) %>%
    mutate(sig = ifelse(p >= 0.05, 'NS',
                        ifelse(p >= 0.01, '*',
                               ifelse(p >= 0.001, '**', '***'))),
           species = sp,
           litter = lt,
           trap = tr) %>%
    select(species, trap, litter, Var, tau, p, sig)
  ## Fishers Z test
  ## Test whether skill of maximal library length significantly higher than the minimal
  fishersTbl = lapply(varList, function(v){
    timeStepsInfo = strsplit(lib, split = ' ', fixed = T)[[1]]
    infoTbl = fullResTbl %>% subset(Var == v) %>% arrange(Size)
    rho1 = infoTbl[['Rho']][1] # Rho in minimal length
    rho2 = infoTbl[['Rho']][dim(infoTbl)[1]] # Rho in maximal length
    # Fisher's Z transformation
    fi1 = 0.5*log((1+rho1)/(1-rho1)); fi2 = 0.5*log((1+rho2)/(1-rho2))
    sdExpected = sqrt(1/(dim(infoTbl)[1] - 3))
    pValue = 2*(1-pnorm(abs(fi1 - fi2), sd = sdExpected))
    # Output
    outputTbl = data.frame(species = sp, 
                           trap = tr,
                           litter = lt,
                           Var = v,
                           fisherMin = fi1,
                           fisherMax = fi2,
                           ExpSD = sdExpected,
                           p = pValue) %>%
      mutate(sig = ifelse(p >= 0.05, 'NS',
                          ifelse(p >= 0.01, '*',
                                 ifelse(p >= 0.001, '**', '***'))))
    return(outputTbl)
  }) %>% do.call(rbind.data.frame, args = .)
  # Store to output
  output = list('curve' = factorPerf, # CCM Curve table of factors on litter fall
                'tbl' = fullResTbl, # Table of rho, null model CI & significant test
                'mono' = kendallTbl, # Kendall table of monotonic test by variable
                'fishers' = fishersTbl, # Fishers Z test for rho improvement
                'nullTbl' = nullCCMRes, # A list of null model rho & CI
                'nullData' = nullModelList) # A list of surrogate data
  return(output)
}
corrTest <- function(mat, ...)
{
  mat = as.matrix(mat)
  n = ncol(mat)
  pMat = matrix(NA, n, n)
  diag(pMat) = 0
  for (i in 1:(n-1)){
    for (j in 1:(n-1)){
      tmp = cor.test(mat[,i], mat[,j], ...)
      pMat[i,j] = pMat[j, i] = tmp$p.value
    }
  }
  colnames(pMat) = rownames(pMat) = colnames(mat)
  return(pMat)
}

envCleaning <- function(data, lb = .001, ub = .999)
{
  q1 = quantile(data, lb, na.rm = T)
  q99 = quantile(data, ub, na.rm = T)
  if (length(data[data < q1]) != 0) data[data < q1] = q1
  if (length(data[data > q99]) != 0) data[data > q99] = q99
  return (data)
}

ltCleaning <- function(datafile, lb = .001, ub = .999)
  # A function to remove outsiders and invalid values from the original data set
  # @param datafile: data frame named 'dfsp' loaded from 'litter_allsp9919.Rdata'
  # @output ds: data frame with cleaned data and standard species names
{
  # Filter outliners
  ds = lapply(unique(datafile$species), function(sp){
    lapply(c('FDM', 'LDM', 'TDM'), function(lt){
      d = datafile %>% subset(species == sp & litter == lt)
      res = lapply(unique(d$trapID), function(ID){
        df = d %>% subset(trapID == ID) %>% select(date, DM)
        # Setting range border
        q1 = quantile(df$DM, lb, na.rm = T)
        q99 = quantile(df$DM, ub, na.rm = T)
        # Find if there are points outside of range
        if (dim(df[(df$DM < q1),])[1] != 0) {
          df[(df$DM < q1),]$DM = q1
        }
        if (dim(df[(df$DM > q99),])[1] != 0) {
          df[(df$DM > q99),]$DM = q99
        }
        res = data.frame(df,
                         species = sp,
                         litter = lt,
                         trapID = ID) %>%
          select(date, species, litter, trapID, DM)
        return(res)
      }) %>% do.call(rbind.data.frame, args = .)
      return(res)
    }) %>% do.call(rbind.data.frame, args = .)
  }) %>% do.call(rbind.data.frame, args = .)
  ds = ds %>% mutate(species = factor(species, 
                                      levels = c('83a Am', '83a Ko', 
                                                 '28a Ko', '28a Sa', '26a Sc')))
  return(ds)
}

# A function to assign significance level based on p values
markSignificance = function(n) {
  if (n >= 0.05) {
    mark = 'ns'
  } else if (n >= 0.01) {
    mark = '*'
  } else if (n >= 0.001) {
    mark = '**'
  } else {
    mark = '***'
  }
  return(mark)
}

# A function to perform calculations with moving-window slices
# Note: Slices for litter should be named by trap ID & litter type in sequence
# Note: Slices for environment should be named by variables
movCal = function(slices,
                  type = c('litter', 'env'),
                  fun = c('mean'),
                  ltList = c('LDM', 'TDM')) {
  res = lapply(fun, function(param){
    if (param == 'mean'){
      calFun = 'mean'
    } else {
      stop("Function not defined. Please modify the code of the function.")
    }
    if (type == 'litter') {
      output = lapply(names(slices), function(tr){
        byLT = lapply(ltList, function(lt){
          targetList = slices[[tr]][[lt]]
          result = lapply(1:length(targetList), function(dateID){
            eval(parse(text = paste('calRes = ',
                                    calFun,
                                    '(targetList[[dateID]])',
                                    sep = '')))
            outRow = data.frame(trapID = as.numeric(tr),
                                litter = lt,
                                date = dateID,
                                value = calRes,
                                param = param)
            return(outRow)
          }) %>% do.call(rbind.data.frame, args = .)
        }) %>% do.call(rbind.data.frame, args = .)
      }) %>% do.call(rbind.data.frame, args = .)
    } else if (type == 'env') {
      output = lapply(names(slices), function(v){
        targetList = slices[[v]]
        lapply(1:length(targetList), function(dateID){
          eval(parse(text = paste('calRes = ',
                                  calFun,
                                  '(slices[[v]][[dateID]])',
                                  sep = '')))
          outTbl = data.frame(var = v,
                              date = dateID,
                              value = calRes,
                              param = param)
          return(outTbl)
        }) %>% do.call(rbind.data.frame, args = .)
      }) %>% do.call(rbind.data.frame, args = .)
    } else {
      stop("Type not defined. Please modify the code of the function.")
    }
  }) %>% do.call(rbind.data.frame, args = .)
  return(res)
}

pwLook = function(meanTbl, compareResult,
                  lt = c('FDM', 'LDM')){
  # Order by mean
  baseTbl = meanTbl %>% subset(litter == lt) %>% arrange(desc(biomass))
  spOrder = baseTbl$species
  # Comparison table
  cpTbl = compareResult %>% subset(litter == lt)
  cpTblMirror = cpTbl
  names(cpTblMirror) = c('litter', 'sp2', 'sp1', 'pvalue', 'sig')
  difMtx = bind_rows(cpTbl, cpTblMirror) %>%
    dplyr::select(sp1, sp2, sig) %>%
    reshape2::dcast(sp1~sp2, drop = F, fill = NA, value.var = 'sig')
  rownames(difMtx) = colnames(difMtx)[-1]
  difMtx = difMtx[match(spOrder, difMtx$sp1),] # reorder
  # initial mark matrix
  markMtx = matrix(data = '', nrow = length(spOrder), ncol = length(spOrder))
  rownames(markMtx) = colnames(markMtx) = spOrder
  # assign significance mark
  for (spID in 1:length(spOrder)) {
    targetSp = spOrder[spID]
    cpCol = difMtx[targetSp]
    if (length(which(cpCol == 'ns')) == 0) {
      markMtx[spID, spID] = 1
    } else if (length(which(cpCol == 'ns')) < length(spOrder)) {
      markMtx[c(spID, which(cpCol == 'ns')), spID] = 1
    }
  }
  return(markMtx)
}

scaleArray <- function(x, na.rm = FALSE)
{
  (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
}
