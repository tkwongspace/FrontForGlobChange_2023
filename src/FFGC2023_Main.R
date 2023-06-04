# Time-series Analysis of Litter Production in Futian Mangrove
# Gamma version
#
# (c) Zijian HUANG 2023

# ********************
# 0. Preparations ----
# ********************

library(dplyr); library(lubridate); library(ggpubr); library(zoo)

# # Path for input data
# dataPath = "/Volumes/TKssd/Graduate/dataBackup/Manuscript1Backup/data"

# Tools
source('FFGC2023_Tools.R')

# Colors
base.Colors = list('White' = '#F2F2F2',
                   'LightGrey' = '#9BAFB6',
                   'Blue' = '#2E6171',
                   'SeaBlue' = '#27187E',
                   'Cyan' = '#00798C',
                   'Orange' = '#F79824',
                   'Red' = '#772014',
                   'Purple' = '#42253B',
                   'Violet' = '#A05C7B',
                   'Eggplant' = '#673C4F',
                   'Pink' = '#E26D5A',
                   'Aqua' = '#4CE0B3',
                   'LightGreen' = '#B6C649',
                   'Green' = '#488B49',
                   'DarkLava' = '#524632',
                   'Black' = '#07020D')

## 0.1 Load data ----
# ## Harvested litter biomass by traps
# {
#   # Mind that the first column could cause error when reading in Windows OS
#   tp0 = read.csv(paste(dataPath, "litterByTrap.csv", sep = "/"))
#   # Fill missing date by linear interpolation
#   tp = tp0 %>%
#     mutate(species = paste(Age+1, "a ", Species, sep = ""),
#            litter = paste(Type, "DM", sep = ""),
#            date = ym(paste(Year, Month, sep = "-")) %m+%
#              months(1) - 1) %>%
#     select(date, species, litter, trapID = TrapID, DM, 
#            Year, Month) %>%
#     subset(species != "29a Sc")
#   tp = lapply(unique(tp$species), function(sp){
#     lapply(unique(tp$litter), function(lt){
#       df = tp %>% subset(species == sp & litter == lt)
#       res = lapply(unique(df$trapID), function(ID){
#         gapCheck = data.frame(date = seq(as.Date('1999-01-01'),
#                                          as.Date('2019-10-01'),
#                                          by = 'month')) %>%
#           mutate(date = date %m+% months(1) - 1) %>%
#           left_join(df %>% subset(trapID == ID),
#                     by = 'date')
#         gapFilled = data.frame(
#           date = gapCheck$date,
#           species = sp,
#           litter = lt,
#           trapID = ID,
#           DM = na.approx(gapCheck$DM)
#         )
#         return(gapFilled)
#       }) %>% do.call(rbind.data.frame, args = .) 
#       return(res)
#     }) %>% do.call(rbind.data.frame, args = .)
#   }) %>% do.call(rbind.data.frame, args = .) %>%
#     ltCleaning() %>%
#     select(date, species, litter, trapID, DM)
#   rm(tp0)
# }
# 
# ## Environmental data
# {
#   # Temperature & Rainfall
#   tr = read.csv(paste(dataPath, "LFS_weather.csv", sep = "/")) %>%
#     mutate(date = ymd(paste(year, month, '01', sep = '-')) %m+% months(1) -1) %>%
#     select(Date = date,
#            AVGT = temp, 
#            PRCP = prcp)
#   tr = data.frame(Date = tr$Date, 
#                   lapply(tr[,-1], envCleaning, lb = .05, ub = .95))
#   # Sea levels
#   sl = read.csv(paste(dataPath, "TBT_tides.csv", sep = "/")) %>%
#     mutate(date = ymd(paste(year, month, '01', sep = '-')) %m+% months(1) -1) %>%
#     select(Date = date, SL = msl)
#   sl = data.frame(Date = sl$Date, 
#                   SL = envCleaning(na.approx(sl[,-1]), lb = .05, ub = .95))
#   # Salinity & Nutrients
#   wq = read.csv(paste(dataPath, "HKEPD_quality.csv", sep = "/")) %>%
#     mutate(SampleDate = as_date(Date),
#            Date = ymd(paste(Year, Month, '01', sep = '-')) %m+% months(1) -1)
#   wqDateRange = data.frame(date = seq(as.Date('1998-01-01'), 
#                                       as.Date('2020-10-01'), by = 'month')) %>%
#     mutate(Date = date %m+% months(1) - 1) %>% select(Date)
#   wqGapFilled = wqDateRange %>%
#     left_join(wq %>% subset(Station == 'DM2'), by = 'Date') %>%
#     select(-Station, -Date, -Year, -Month, -SampleDate) %>%
#     na.approx()
#   wqGapCheck = data.frame(Date = wqDateRange, wqGapFilled)
#   wq = data.frame(Date = wqGapCheck$Date, 
#                   lapply(wqGapCheck %>% select(-Date),
#                          envCleaning, lb = .05, ub = .95))
#   # Join all environmental table
#   wt = left_join(tr, sl, by = 'Date') %>% left_join(wq, by = 'Date')
#   rm(tr, sl, wq, wqDateRange, wqGapCheck, wqGapFilled)
# }

# *************************
# 1. Data Descriptions ----
# *************************
# Some pre-arrangements
# List of environmental drivers
varList = c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP')
# Annual Summation
# >> litter production
annLit = tp %>% mutate(Year = year(date)) %>%
  group_by(species, trapID, litter, Year) %>%
  summarise(annSum = sum(DM) * 10) %>%  # x10 to convert unit from g/m2 -> kg/ha
  as.data.frame()
annLitSum = annLit %>%
  group_by(species, litter) %>%
  summarise(biomass = mean(annSum), SE = sd(annSum)/sqrt(n())) %>%
  as.data.frame()
# >> environment
annEnv = wt %>% subset(Date >= ym('1999-01') & Date <= ym('2019-12')) %>%
  mutate(Year = year(Date)) %>%
  group_by(Year) %>%
  summarise(Tmax = max(AVGT), Tmin = min(AVGT), Tmean = mean(AVGT), Tvar = var(AVGT),
            Rmax = max(PRCP), Rmin = min(PRCP), Rsum = sum(PRCP), Rvar = var(PRCP),
            Lmax = max(SL), Lmin = min(SL), Lmean = mean(SL), Lvar = var(SL),
            Smax = max(SAL), Smin = min(SAL), Smean = mean(SAL), Svar = var(SAL),
            Nmax = max(TN), Nmin = min(TN), Nmean = mean(TN), Nvar = var(TN),
            Pmax = max(TP), Pmin = min(TP), Pmean = mean(TP), Pvar = var(TP)) %>%
  dplyr::select(-Year) %>%
  scale() %>%
  data.frame()

## 1.1 Between-assemblage annual comparisons ----
# one-way ANOVA
spComb = combn(unique(annLit$species), 2, simplify = F)
#
onewayTbl = lapply(spComb, function(sp){
  lapply(c("FDM", "LDM", "TDM"), function(lt){
    compr = oneway.test(DM ~ species, 
                        data = annLit %>%
                          subset(litter == lt & species %in% sp) %>%
                          mutate(DM = log(annSum)))
    resTbl = data.frame(litter = lt,
                        sp1 = sp[1],
                        sp2 = sp[2],
                        pvalue = compr$p.value)
    return(resTbl)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  rowwise() %>% mutate(sig = markSignificance(pvalue))

## 1.2 Temporal trends ----
# sliding windows
# >> litterfall
tpMovWin = lapply(baseLookUp$trapID, function(tr){
  res = lapply(c('LDM', 'FDM'), function(lt){
    lb = tp %>% subset(trapID == tr & litter == lt)
    mv = runner::runner(lb$DM, f = c, k = 36)[-c(1:35)]
    return(mv)
  })
  names(res) = c('LDM', 'FDM')
  return(res)
}); names(tpMovWin) = baseLookUp$trapID
# Calculate mean for each slice
tpMovMean = movCal(slices = tpMovWin,
                   type = 'litter',
                   fun = 'mean',
                   ltList = c('LDM', 'FDM')) %>%
  # date corresponds to the middle of the moving window
  mutate(date = ym('1999-01') %m+% months(date + 17)) %>%
  left_join(baseLookUp, by = 'trapID')
# Analyse temporal trends
library(nlme)
tpMovTrend = tpMovMean %>% 
  mutate(time = scale(date)) %>%
  nest_by(species, litter) %>%
  mutate(mod = list(lme(value~time, data = data,
                        random = ~1 | trapID, method = 'ML',
                        correlation = corAR1()))) %>%
  summarise(slope = summary(mod)$coefficients$fixed[[2]],
            p = round(anova(mod)$`p-value`[2], digits = 4)) %>%
  mutate(sig = ifelse(p < 0.05, 'solid', 'dashed'))
detach("package:nlme")
#
# >> environment
# Slice by moving windows
wtMovWin = lapply(varList, function(v){
  colID = which(names(wt) == v)
  lb = wt %>% subset(Date >= ym('1999-01') & Date <= ym('2019-12')) %>%
    dplyr::select(1, all_of(colID))
  names(lb) = c('Date', 'Value')
  mv = runner::runner(lb$Value, f = c, k = 36)[-c(1:35)]
  return(mv)
}); names(wtMovWin) = varList
# Calculate mean for each slice
wtMovMean = movCal(slices = wtMovWin,
                   type = 'env',
                   fun = 'mean') %>%
  # date corresponds to the middle of the moving window
  mutate(date = ym('1999-01') %m+% months(date + 17))
# Analyse temporal trends
library(nlme)
wtMovTrend = wtMovMean %>% 
  mutate(time = scale(date)) %>%
  nest_by(var) %>%
  mutate(mod = list(gls(value~time, data = data,
                        method = 'ML',
                        correlation = corAR1()))) %>%
  summarise(slope = summary(mod)$tTable[2],
            p = round(summary(mod)$tTable[8], digits = 4)) %>%
  mutate(sig = ifelse(p < 0.05, 'solid', 'dashed'))
detach("package:nlme")

# ***********************************
# 2. Empirical Dynamic Modelling ----
# ***********************************
#
library(rEDM)
library(hash)
#
## 2.1 Create table for EDM ----
# Differentiated version of litter biomass for stationary
# To test stationary in 'Stationarity_Test.R')
edm_dmTbl = lapply(unique(tp$trapID), function(tr){
  byLT = lapply(c("FDM", "LDM"), function(lt){
    dfBase = data.frame(
      time = subset(tp, trapID == tr & litter == lt)$date,
      DM = c(NA, diff(subset(tp, trapID == tr & litter == lt)$DM))
    )
    dfBase = dfBase[-1,]
    return(dfBase)
  }); names(byLT) = c("FDM", "LDM")
  return(byLT)
}); names(edm_dmTbl) = unique(tp$trapID)
#
## 2.2 Optimization of embedding size ----
#
# Identify optimal E for each assemblage and litter type
# Optimal E is selected by uni-variate SSR model
edm_findE = lapply(unique(baseLookUp$trapID), function(tr){
  lapply(c("FDM", "LDM"), function(lt){
    findE = EmbedDimension(dataFrame = edm_dmTbl[[as.character(tr)]][[lt]],
                           lib = '1 179',  # 1999-02 to 2013-12
                           pred = '180 249',  # 2014-01 to 2019-10
                           maxE = 13,
                           columns = 'DM',
                           target = 'DM',
                           embedded = F,
                           showPlot = F)
    # add species and litter marks
    findE = mutate(findE,
                   trap = tr,
                   species = baseLookUp$species[which(baseLookUp$trapID == tr)],
                   litter = lt)
    return(findE)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .)
#
# Find rows with max rho for each assemblage
edm_Ecandidates = lapply(levels(edm_findE$species), function(sp){
  byLT = lapply(c("FDM", "LDM"), function(lt){
    eSummary = edm_findE %>%
      group_by(trap, litter) %>%
      mutate(maxRho = max(rho),
             optimal = ifelse(rho == maxRho, T, F)) %>%
      subset(optimal)
    output = unique(subset(eSummary, species == sp & litter == lt)$E)
    return(output)
  }); names(byLT) = c("FDM", "LDM"); return(byLT)
}); names(edm_Ecandidates) = levels(edm_findE$species)
#
# Select optimal E with maximum prediction skill
edm_Eopt = pbapply::pblapply(levels(edm_findE$species), function(sp){
  byLT = lapply(c("FDM", "LDM"), function(lt){
    # locate E candidates
    eList = edm_Ecandidates[[sp]][[lt]]
    byE = lapply(eList, function(testE){
      # locate traps
      trapList = baseLookUp[which(baseLookUp$species == sp), 2]
      byTrap = lapply(trapList, function(tr){
        # input to Simplex
        inputdf = edm_dmTbl[[as.character(tr)]][[lt]]
        models = Simplex(dataFrame = inputdf,
                         lib = '1 179',
                         pred = '180 249',
                         E = testE,
                         columns = 'DM',
                         target = 'DM',
                         embedded = F,
                         showPlot = F)
        rho = ComputeError(models$Observations,
                           models$Predictions)$rho
        return(rho)
      }) # byTrap
      # create data frame
      resTbl = data.frame(species = sp,
                          litter = lt,
                          E = testE,
                          rho = mean(unlist(byTrap)),
                          se = sd(unlist(byTrap))/sqrt(length(byTrap)))
      return(resTbl)
    }) %>% do.call(rbind.data.frame, args = .)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  group_by(species, litter) %>%
  mutate(maxRho = max(rho),
         optimal = ifelse(rho == maxRho, T, F)) %>%
  ungroup() %>%
  dplyr::select(species, litter, E, rho, se, optimal)
#
# Create list for optimal E (easy to use)
edmE = pbapply::pblapply(levels(edm_findE$species), function(sp){
  byLT = lapply(c("FDM", "LDM"), function(lt){
    optimalRow = edm_Eopt %>%
      subset(optimal & species == sp & litter == lt)
    output = list('E' = optimalRow$E,
                  'rho' = optimalRow$rho)
    return(output)
  }); names(byLT) = c("FDM", "LDM"); return(byLT)
}); names(edmE) = levels(edm_findE$species)

rm(edm_findE, edm_Eopt, edm_Ecandidates)

## 2.3 CCM ----
# Warning: This section takes a long time
edm_CCM = lapply(levels(baseLookUp$species), function(sp){
  trapList = baseLookUp$trapID[which(baseLookUp$species == sp)]
  byLT = lapply(c("FDM", "LDM"), function(lt){
    byTrap = lapply(trapList, function(tr){
      message(paste('> Process starts for ', sp, ' ', lt,
                    ' (trap #', tr, ') at ', Sys.time(), sep = ''))
      ccmPerform = ccm.perform(tr = tr,
                               lt = lt,
                               lib = '36 236 3',
                               nullMod = 'ebi',
                               nullModRep = 50)
      # save to disk
      fileName = paste('CCM', tr, lt, sep = '_')
      save(ccmPerform, file = paste(dataPath, '/ccm/',
                                    fileName, '.Rdata', sep = ''))
      message(paste('> File output of ', sp, ' ', lt,
                    ' (trap #', tr, ') finished at ', Sys.time(), '.', sep = ''))
      return(ccmPerform)
    }); names(byTrap) = trapList; return(byTrap)
  }); names(byLT) = c("FDM", "LDM"); return(byLT)
}); names(edm_CCM) = levels(baseLookUp$species)

# => Summary of CCM results please refer to M1_Visualization.R Sec.4.1

## 2.4 S-Map models ----
#
# Table for SMap input
edmSmapBase = lapply(names(edm_dmTbl), function(tr){
  lapply(c('FDM', 'LDM'), function(lt){
    edm_dmTbl[[tr]][[lt]] %>% 
      mutate(trap = as.numeric(tr), 
             litter = lt,
             species = baseLookUp$species[which(baseLookUp$trapID == as.numeric(tr))])
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  left_join(climScale, by = c('time' = 'Date'))
#
## Test predictive power by using different E with Simplex
edmSmapFindE = pbapply::pblapply(levels(edmSmapBase$species), function(sp){
  lapply(c('FDM', 'LDM'), function(lt){
    ## trap list
    trapList = baseLookUp$trapID[which(baseLookUp$species == sp)]
    ## process by trap & variables
    byTrap = lapply(trapList, function(tr){
      byV = lapply(varList, function(targetV){
        ### basic table for E finding (remove useless columns)
        inputBase = edmSmapBase %>% subset(trap == tr & litter == lt) %>%
          select(time, DM, Env = all_of(targetV))
        eToTest = 1:13
        ### perform Simplex on each E candidate
        findE = lapply(eToTest, function(optE){
          #### embed tables
          if (optE == 1) {
            inputTbl = inputBase
            colInput = 'Env'
          } else {
            embBase = inputBase %>% select(time, Env)
            inputEmb = Embed(dataFrame = embBase,
                             E = optE,
                             columns = 'Env')
            colInput = c('Env', paste('Env', 1:(optE-1), sep = ''))
            names(inputEmb) = colInput
            inputTbl = inputBase %>% select(time, DM) %>%
              cbind(inputEmb) %>% slice(-(1:(optE-1)))
          }
          #### Simplex
          trainLib = '1 179'
          testLib = paste(180, dim(inputTbl)[1], collapse = ' ')
          mod = Simplex(dataFrame = inputTbl,
                        lib = trainLib,
                        pred = testLib,
                        E = optE,
                        columns = paste(colInput, collapse = ' '),
                        target = 'DM',
                        embedded = T,
                        showPlot = F)
          rho = ComputeError(mod$Observations, mod$Predictions)$rho
          resTbl = data.frame(species = sp,
                              litter = lt,
                              trapID = tr,
                              Var = targetV, 
                              E = optE, 
                              rho = rho)
          return(resTbl)
        }) %>% do.call(rbind.data.frame, args = .)
        return(findE)
      }) %>% do.call(rbind.data.frame, args = .) #byV
    }) %>% do.call(rbind.data.frame, args = .) #byTrap
  }) %>% do.call(rbind.data.frame, args = .) #lt
}) %>% do.call(rbind.data.frame, args = .) #sp
#
# Candidate E's for each assemblage under given contexts
edmSmapEcand = edmSmapFindE %>%
  group_by(litter, trapID, Var) %>%
  mutate(maxRho = max(rho), optimal = ifelse(rho == maxRho, T, F)) %>%
  subset(optimal) %>%
  as.data.frame() %>%
  group_by(species, litter, Var) %>%
  summarise(candE = paste(unique(E), collapse = ', '),
            rho = mean(rho))
#
# Identify optimal E 
edmSmapEopt = lapply(unique(edmSmapEcand$species), function(sp){
  lapply(c('FDM', 'LDM'), function(lt){
    baseInfo = edmSmapEcand %>% subset(species == sp & litter == lt)
    trapList = baseLookUp$trapID[which(baseLookUp$species == sp)]
    varToLook = baseInfo$Var
    ## search each variable
    byV = lapply(varToLook, function(targetV){
      paramRow = baseInfo %>% subset(Var == targetV)
      eToLook = as.numeric(strsplit(paramRow$candE, split = ', ', fixed = T)[[1]])
      if (length(eToLook) >1){
        byE = lapply(eToLook, function(testE){
          #### perform Simplex under same E on each trap
          meanRho = lapply(trapList, function(tr) {
            inputBase = edmSmapBase %>% subset(trap == tr & litter == lt) %>%
              select(time, DM, Env = all_of(targetV))
            if (testE == 1){
              colInput = 'Env'
              inputTbl = inputBase
            } else {
              embBase = inputBase %>% select(time, Env)
              inputEmb = Embed(dataFrame = embBase,
                               E = testE,
                               columns = 'Env')
              colInput = c('Env', paste('Env', 1:(testE-1), sep = ''))
              names(inputEmb) = colInput
              inputTbl = inputBase %>% select(time, DM) %>%
                cbind(inputEmb) %>% slice(-(1:(testE-1)))
            }
            ##### Simplex
            trainLib = '1 179'
            testLib = paste(180, dim(inputTbl)[1], collapse = ' ')
            mod = Simplex(dataFrame = inputTbl,
                          lib = trainLib,
                          pred = testLib,
                          E = testE,
                          columns = paste(colInput, collapse = ' '),
                          target = 'DM',
                          embedded = T,
                          showPlot = F)
            rho = ComputeError(mod$Observations, mod$Predictions)$rho
            return(rho)
          }) %>% unlist() %>% mean()
          resTbl = data.frame(species = sp, litter = lt, Var = targetV,
                              optE = testE, rho = meanRho)
          return(resTbl)
        }) %>% do.call(rbind.data.frame, args = .)
      } else {
        byE = data.frame(species = sp, litter = lt, Var = targetV,
                         optE = paramRow$candE, rho = paramRow$rho)
      }
      return(byE)
    }) %>% do.call(rbind.data.frame, args = .) #byV
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  group_by(species, litter, Var) %>%
  mutate(maxRho = max(rho), optimal = ifelse(rho == maxRho, T, F)) %>%
  ungroup() %>% subset(optimal) %>% 
  select(species, litter, Var, optE, rho) %>% as.data.frame()
#
# Select optimal theta based on optimal E for each trap
edmSmapTheta = pbapply::pblapply(levels(baseLookUp$species), function(sp){
  trapList = baseLookUp$trapID[which(baseLookUp$species == sp)]
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    byV = lapply(varList, function(targetV){
      paramRow = edmSmapEopt %>% 
        subset(species == sp & litter == lt & Var == targetV)
      optE = as.numeric(paramRow$optE)
      #### process by trap
      byTrap = lapply(trapList, function(tr){
        baseTbl = edmSmapBase %>% 
          subset(trap == tr & litter == lt) %>%
          dplyr::select(time, DM, Env = all_of(targetV))
        if (optE == 1){
          colInput = 'Env'
          inputTbl = baseTbl
        } else {
          embBase = baseTbl %>% select(time, Env)
          inputEmb = Embed(dataFrame = embBase,
                           E = optE,
                           columns = 'Env')
          colInput = c('Env', paste('Env', 1:(optE-1), sep = ''))
          names(inputEmb) = colInput
          inputTbl = baseTbl %>% select(time, DM) %>%
            cbind(inputEmb) %>% slice(-(1:(optE-1)))
        }
        ##### identify optimal theta
        thetaTbl = PredictNonlinear(dataFrame = inputTbl,
                                    lib = '1 179',
                                    pred = paste(180, dim(inputTbl)[1], sep = ' '),
                                    E = optE,
                                    columns = paste(colInput, collapse = ' '),
                                    target = 'DM',
                                    embedded = T,
                                    showPlot = F)
        optTheta = thetaTbl$Theta[which.max(thetaTbl$rho)]
        optRho = thetaTbl$rho[which.max(thetaTbl$rho)]
        ##### Arrange to output
        outputList = list('data' = inputTbl,
                          'col' = paste(colInput, collapse = ' '),
                          'E' = optE,
                          'theta' = optTheta,
                          'rho' = optRho)
        return(outputList)
      }); names(byTrap) = trapList; return(byTrap)
    }); names(byV) = varList; return(byV)
  }); names(byLT) = c('FDM', 'LDM'); return(byLT)
}); names(edmSmapTheta) = levels(baseLookUp$species)
# Extract optimal model parameters (by traps) to table
# note: the following data has been saved to file
edmSmapParam = lapply(levels(baseLookUp$species), function(sp){
  lapply(c('FDM', 'LDM'), function(lt){
    ## Retrieve optimal parameters for all traps of the species
    listToCheck = edmSmapTheta[[sp]][[lt]]
    varToCheck = names(listToCheck)
    tableConverted = lapply(varToCheck, function(V){
      lapply(1:length(listToCheck[[V]]), function(tID){
        baseInfo = listToCheck[[V]][[tID]]
        outTbl = data.frame(species = sp,
                            litter = lt,
                            var = V,
                            trap = as.numeric(names(listToCheck[[V]])[tID]),
                            colInput = baseInfo[['col']],
                            E = baseInfo[['E']],
                            theta = baseInfo[['theta']],
                            rho = baseInfo[['rho']])
        return(outTbl)
      }) %>% do.call(rbind.data.frame, args = .)
    }) %>% do.call(rbind.data.frame, args = .)
    return(tableConverted)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  group_by(species, litter, var, colInput, E) %>%
  summarise(avg = mean(rho),
            theta = paste(unique(theta), collapse = ', ')) %>%
  ungroup() %>% as.data.frame()
#
save(edmSmapParam, file = paste(dataPath, "/SMap/param.Rdata", sep = ''))
# Optimize parameter and construct S-Map model for each assemblage
# note: the following data has been saved to file
edmSmapMods = pbapply::pblapply(levels(baseLookUp$species), function(sp){
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    trapList = baseLookUp$trapID[which(baseLookUp$species == sp)]
    byVar = lapply(varList, function(V){
      paramRow = edmSmapParam %>% 
        subset(species == sp & litter == lt & var == V)
      optE = paramRow$E
      colInput = paramRow$colInput
      thetaToLook = as.numeric(strsplit(paramRow$theta, split = ', ', fixed = T)[[1]])
      #### perform by theta
      byTheta = lapply(thetaToLook, function(optT){
        #### perform S-Map by trap
        byTrap = lapply(trapList, function(tr){
          inputTbl = edmSmapTheta[[sp]][[lt]][[V]][[as.character(tr)]][['data']]
          mod = SMap(dataFrame = inputTbl,
                     lib = '1 179',
                     pred = paste(180, dim(inputTbl)[1], collapse = ' '),
                     E = optE,
                     theta = optT,
                     columns = colInput,
                     target = 'DM',
                     embedded = T)
          smapTbl = data.frame(trap = tr, mod$predictions)
          return(smapTbl)
        }); names(byTrap) = trapList
        #### summary rho under given theta
        thetaSummary = lapply(byTrap, function(eachTrap){
          outTbl = data.frame(theta = optT,
                              rho = ComputeError(eachTrap$Observations,
                                                 eachTrap$Predictions)$rho)
          return(outTbl)
        }) %>% do.call(rbind.data.frame, args = .) %>%
          group_by(theta) %>% summarise(rho = mean(rho))
        #### save to output
        outList = list('base' = byTrap, 
                       'rho' = thetaSummary)
      }); names(byTheta) = thetaToLook
      ### optimize theta
      evalTheta = lapply(1:length(byTheta), function(ID){
        byTheta[[ID]][['rho']]
      }) %>% do.call(rbind.data.frame, args = .) %>%
        mutate(maxRho = max(rho), optimal = ifelse(rho == maxRho, T, F)) %>%
        subset(optimal)
      ### optimal theta & SMap prediction table
      outList = list('theta' = evalTheta$theta,
                     'rho' = evalTheta$rho,
                     'base' = byTheta[[as.character(evalTheta$theta)]][['base']])
      return(outList)
    }); names(byVar) = varList; return(byVar)
  }); names(byLT) = c('FDM', 'LDM'); return(byLT)
}); names(edmSmapMods) = levels(baseLookUp$species)
#
save(edmSmapMods, file = paste(dataPath, "/SMap/optModels.Rdata", sep = ''))
# => Now can move to visualization of uni-variate S-Map

## 2.5 Scenario exploration ----
# note: the following data has been saved to file
edmScenExp = pbapply::pblapply(levels(baseLookUp$species), function(sp){
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    trapList = baseLookUp$trapID[which(baseLookUp$species == sp)]
    ## parameters of the system
    baseParam = edmSmapParam %>% subset(species == sp & litter == lt)
    varToLook = baseParam$var
    ## perform by driver
    byVar = lapply(varToLook, function(V){
      optE = subset(baseParam, var == V)$E
      optTheta = edmSmapMods[[sp]][[lt]][[V]][['theta']]
      colInput = subset(baseParam, var == V)$colInput
      ### by trap
      byTrap = lapply(trapList, function(tr){
        baseTbl = edmSmapBase %>% 
          subset(trap == tr & litter == lt) %>%
          dplyr::select(time, DM, Env = all_of(V))
        #### modify environmental driver
        valueToChange = baseTbl$Env
        dX = 0.1 * sd(valueToChange)
        xPlus = xMinus = valueToChange
        maxLength = length(valueToChange)
        xPlus[(180+(optE-1)):maxLength] = xPlus[(180+(optE-1)):maxLength] + dX
        xMinus[(180+(optE-1)):maxLength] = xMinus[(180+(optE-1)):maxLength] - dX
        plusTbl = data.frame(time = baseTbl$time, Env = xPlus)
        minusTbl = data.frame(time = baseTbl$time, Env = xMinus)
        #### embed table
        plusEmb = Embed(dataFrame = plusTbl, E = optE, columns = 'Env')
        minusEmb = Embed(dataFrame = minusTbl, E = optE, columns = 'Env')
        names(plusEmb) = names(minusEmb) = strsplit(colInput, split = ' ', fixed = T)[[1]]
        plusTbl = bind_cols(baseTbl %>% select(time, DM), plusEmb) %>%
          slice(-(1:(optE-1)))
        minusTbl = bind_cols(baseTbl %>% select(time, DM), minusEmb) %>%
          slice(-(1:(optE-1)))
        rm(plusEmb, minusEmb)
        #### S-Map
        dxPlus = SMap(dataFrame = plusTbl,
                      lib = '1 179',
                      pred = paste(180, dim(plusTbl)[1], collapse = ' '),
                      E = optE,
                      theta = optTheta,
                      columns = colInput,
                      target = 'DM',
                      embedded = T)
        dxMinus = SMap(dataFrame = minusTbl,
                       lib = '1 179',
                       pred = paste(180, dim(minusTbl)[1], collapse = ' '),
                       E = optE,
                       theta = optTheta,
                       columns = colInput,
                       target = 'DM',
                       embedded = T)
        scneOut = bind_rows(data.frame(mode = 'plus',
                                       dxPlus$predictions),
                            data.frame(mode = 'minus',
                                       dxMinus$predictions)) %>%
          mutate(dX = dX) %>%
          select(time, dX, Y1 = Predictions, mode)
        #### join S-Map predictions with unchanged environment
        smapBase = edmSmapMods[[sp]][[lt]][[V]][['base']][[as.character(tr)]]
        scneRes = smapBase %>% select(time, Y0 = Predictions) %>% 
          left_join(scneOut, by = 'time') %>%
          subset(!is.na(Y0)) %>%
          mutate(species = sp, litter = lt, trap = tr, var = V,
                 dY = Y1 - Y0,
                 sens = ifelse(mode == 'plus', dY/dX, -dY/dX)) %>%
          select(species, litter, trap, time, var, sens, mode, Y1, Y0, dY, dX)
        return(scneRes)
      }) %>% do.call(rbind.data.frame, args = .)
    }) %>% do.call(rbind.data.frame, args = .) 
  }) %>% do.call(rbind.data.frame, args = .) 
}) %>% do.call(rbind.data.frame, args = .)
#
save(edmScenExp, file = paste(dataPath, "/SMap/scneExpFull.Rdata", sep = ''))
# => Now can move to visualization for scenario exploration

# 3. Structural Equation Modelling ----
library(lavaan)
## 3.1 Reproduction - Exotic ----
inputTbl = annLit %>% 
  subset(species %in% c('28a Sa', '26a Sc') & litter == 'FDM') %>%
  left_join(annEnv, by = 'Year') %>%
  mutate(DM = log(annSum))
#
optimalMod <- '
  # direct effects
  DM ~ 1 + tb*Tmean + rb*Rsum + sb*Smean + nb*Nmean
  # environmental interactions
  Nmean ~ 1 + rn*Rsum + ln*Lmean + sn*Smean
  Smean ~ 1 + ts*Tmean + rs*Rsum
  # covariance
  Lmean ~~ Rsum
  Lmean ~~ Tmean
  # indirect effects
  Lind := ln*nb
  Sind := sn*nb
  Tind := ts*(Sind+sb)
  Rind := rs*(Sind+sb) + rn*nb
  # total effects
  Ttot := tb + Tind
  Rtot := rb + Rind
  Stot := sb + Sind
'
#
fitSEM = sem(optimalMod, data = inputTbl)
summary(fitSEM, standardized = T, ci = T, fit.measures = T)
#
## 3.2 Reproduction - Native ----
inputTbl = annLit %>% 
  subset(litter == 'FDM' & !species %in% c('28a Sa', '26a Sc')) %>%
  left_join(annEnv, by = 'Year') %>%
  mutate(DM = log(annSum))
#
optimalMod <- '
  # direct effects
  DM ~ 1 + sb*Smean
  # environmental interactions
  Smean ~ 1 + ts*Tmean + rs*Rsum
  # indirect effects
  Tind := ts*sb
  Rind := rs*sb
'
#
fitSEM = sem(optimalMod, data = inputTbl)
summary(fitSEM, standardized = T, ci = T, fit.measures = T)
#
##  3.3 Leaves - Exotic ----
inputTbl = annLit %>% 
  subset(species %in% c('28a Sa', '26a Sc') & litter == 'LDM') %>%
  left_join(annEnv, by = 'Year') %>%
  mutate(DM = log(annSum))
#
optimalMod <- '
  # direct effects
  DM ~ 1 + tb*Tmean + rb*Rsum + sb*Smean + lb*Lmean
  # environmental interactions
  Smean ~ 1 + ts*Tmean + rs*Rsum
  # covariance
  Lmean ~~ Rsum
  Lmean ~~ Tmean
  # indirect effects
  Tind := ts*sb
  Rind := rs*sb
  # total effects
  Ttot := tb + Tind
  Rtot := rb + Rind
'
#
fitSEM = sem(optimalMod, data = inputTbl)
summary(fitSEM, standardized = T, ci = T, fit.measures = T)

## 3.4 Leaves - Native ----
inputTbl = annLit %>% 
  subset(litter == 'LDM' & !species %in% c('28a Sa', '26a Sc')) %>%
  left_join(annEnv, by = 'Year') %>%
  mutate(DM = log(annSum))
#
optimalMod <- '
  # direct impacts
  DM ~ 1 + tb*Tmean + sb*Smean + nb*Nmean
  # environmental interactions
  Smean ~ 1 + ts*Tmean + rs*Rsum
  Nmean ~ 1 + ln*Lmean + sn*Smean
  # covariance
  Lmean ~~ Tmean
  Lmean ~~ Rsum
  Rsum ~~ Nmean
  # indirect effects
  Sind := sn*nb
  Tind := ts*(Sind + sb)
  Rind := rs*(Sind + sb)
  Lind := ln*nb
  # total effects
  Ttot := tb + Tind
  Stot := sb + Sind
'
#
fitSEM = sem(optimalMod, data = inputTbl)
summary(fitSEM, standardized = T, ci = T, fit.measures = T)

#
# 4. Redundancy Analysis ----
#
library(vegan)
#
## 4.1 Annual reproductive litterfall ----
# Arrange table
rdaR = annLit %>% subset(litter == 'FDM') %>%
  group_by(species, trapID) %>%
  mutate(value = scale(annSum)[,1],
         trapMark = paste(species, trapID)) %>%
  ungroup() %>%
  dplyr::select(Year, trapMark, value) %>%
  tidyr::spread(key = trapMark, value = value) %>%
  dplyr::select(-Year)
rownames(annEnv) = 1999:2019; rownames(rdaR) = 1999:2019
# Model the effect of all environmental drivers on FL ratio across traps
rdaMod = rda(rdaR ~ ., data = annEnv %>% select(all_of(selectedVar)))
# Forward selection of variables
rdaSel = ordiR2step(rda(rdaR ~ 1, 
                        data = annEnv %>% select(all_of(selectedVar))),
                    scope = formula(rdaMod),
                    direction = 'forward',
                    R2scope = F,
                    pstep = 1000,
                    trace = T)
#
# ... Omit the section of manual selections ...
#
# Build optimal model
rdaSMod = rda(rdaR ~ Tmin + Rmin + Smin + Nmean + Lmax + Smean, 
              data = annEnv %>% select(all_of(selectedVar)))
# Calculate adjusted R-squared
RsquareAdj(rdaSMod)
# Check weights 
summary(rdaSMod)
# Significance of RDA & variables
anova.cca(rdaSMod, step = 1000)
anova.cca(rdaSMod, step = 1000, by = "term")

## 4.2 Annual leaf litterfall ----
# Arrange table
rdaL = annLit %>% subset(litter == 'LDM') %>%
  group_by(species, trapID) %>%
  mutate(value = scale(annSum)[,1],
         trapMark = paste(species, trapID)) %>%
  ungroup() %>%
  dplyr::select(Year, trapMark, value) %>%
  tidyr::spread(key = trapMark, value = value) %>%
  dplyr::select(-Year)
rownames(annEnv) = 1999:2019; rownames(rdaL) = 1999:2019
# Model the effect of all environmental drivers on FL ratio across traps
rdaMod = rda(rdaL ~ ., data = annEnv %>% select(all_of(selectedVar)))
# Forward selection of variables
rdaSel = ordiR2step(rda(rdaL ~ 1, 
                        data = annEnv %>% select(all_of(selectedVar))),
                    scope = formula(rdaMod),
                    direction = 'forward',
                    R2scope = F,
                    pstep = 1000,
                    trace = T)
#
# ... Omit the section of manual selections ...
#
# Build optimal model
rdaSMod = rda(rdaL ~ Rmax + Lmax + Smax + Pvar + Nmean + Lmin + Smean, 
              data = annEnv %>% select(all_of(selectedVar)))
# Calculate adjusted R-squared
RsquareAdj(rdaSMod)
# Check weights 
summary(rdaSMod)$cont
# Significance of RDA & variables
anova.cca(rdaSMod, step = 1000)
anova.cca(rdaSMod, step = 1000, by = "term")

# ARCHIVE - CCA ----
## CCA on annual reproduction ----
rdaOri = annLit %>% subset(litter == 'FDM') %>%
  group_by(species, trapID) %>%
  mutate(value = annSum,
         trapMark = paste(species, trapID)) %>%
  ungroup() %>%
  dplyr::select(Year, trapMark, value) %>%
  tidyr::spread(key = trapMark, value = value) %>%
  dplyr::select(-Year)
rdaEnv = annEnv %>% select_all(all_of(selectedVar))

cc1 = CCA::cc(rdaOri, rdaEnv)
cc2 = cca(rdaOri ~., data = rdaEnv)
cc3 = ordiR2step(cca(rdaOri~1, data = rdaEnv), scope = formula(cc2), 
                 R2scope = F, pstep = 1000, trace=T)

## CCA on annual leaf ----
rdaOri = annLit %>% subset(litter == 'LDM') %>%
  group_by(species, trapID) %>%
  mutate(value = annSum,
         trapMark = paste(species, trapID)) %>%
  ungroup() %>%
  dplyr::select(Year, trapMark, value) %>%
  tidyr::spread(key = trapMark, value = value) %>%
  dplyr::select(-Year)
rdaEnv = annEnv %>% select_all(all_of(selectedVar))

cc1 = CCA::cc(rdaOri, rdaEnv)
cc2 = cca(rdaOri ~., data = rdaEnv)
cc3 = ordiR2step(cca(rdaOri~1, data = rdaEnv), scope = formula(cc2), 
                 R2scope = F, pstep = 1000, trace=T)

