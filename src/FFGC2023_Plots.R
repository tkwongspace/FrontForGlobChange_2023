# Visualization for Manuscript 1 (Gamma & Published version)
#
# (c) Zijian HUANG 2023

plotPath = "/Volumes/TKssd/Graduate/dataBackup/Manuscript1Backup/plots"

plotSpeciesLabels = c(expression(paste('28a ', italic('S. apetala'))),
                      expression(paste('26a ', italic('S. caseolaris'))),
                      expression(paste('28a ', italic('K. obovata'))),
                      expression(paste('83a ', italic('K. obovata'))),
                      expression(paste('83a ', italic('A. marina'))))

# *************************************************
# Figure 1 - Temporal variations in litterfall ----
# *************************************************
# Require: tp, tpMovTrend from Section 1.2
{
  plotDM = tp %>%
    group_by(species, litter, date) %>%
    summarise(avg = mean(DM), se = sd(DM)/sqrt(n())) %>%
    ungroup() %>% rowwise() %>%
    mutate(lower = avg - se, upper = avg + se) %>%
    left_join(tpMovTrend %>% dplyr::select(species, litter, sig),
              by = c('species', 'litter')) %>%
    mutate(species = factor(species,
                            levels = c('28a Sa', '26a Sc', '28a Ko',
                                       '83a Ko', '83a Am'))) %>%
    subset(!is.na(sig)) # remove TDM rows
  plotLitterSubs = lapply(1:10, function(fID){
    if (fID > 5) plotLT = 'FDM' else plotLT = 'LDM'
    if (fID > 5) markLT = 'Reproductive Parts' else markLT = 'Leaves'
    if (fID%%5 == 0) plotSp = '83a Am' else plotSp = levels(plotDM$species)[fID%%5]
    # subplots
    if (fID < 5) {
      p = ggplot(plotDM %>% subset(litter == plotLT & species == plotSp), 
                 aes(x = date)) +
        geom_line(aes(y = avg), color = base.Colors[['LightGrey']]) +
        # geom_ribbon(aes(ymin = lower, ymax = upper),
        #             fill = base.Colors[['LightGrey']], alpha = .5) +
        geom_smooth(aes(y = avg, linetype = sig), method = 'lm', se = F,
                    color = base.Colors[['Blue']]) +
        scale_linetype_identity() +
        scale_x_continuous(breaks = ym(c('1999-01','2019-01')),
                           labels = c('1999', '2019')) +
        facet_grid(cols = vars(species)) +
        theme_classic() +
        theme(text = element_text(size = 22),
              strip.background = element_blank(),
              strip.text = element_text(size = 22),
              axis.title = element_blank())
    } else if (fID == 5) {
      p = ggplot(plotDM %>% subset(litter == plotLT & species == plotSp) %>%
                   mutate(litter = markLT), 
                 aes(x = date)) +
        geom_line(aes(y = avg), color = base.Colors[['LightGrey']]) +
        # geom_ribbon(aes(ymin = lower, ymax = upper),
        #             fill = base.Colors[['LightGrey']], alpha = .5) +
        geom_smooth(aes(y = avg, linetype = sig), method = 'lm', se = F,
                    color = base.Colors[['Blue']]) +
        scale_linetype_identity() +
        scale_x_continuous(breaks = ym(c('1999-01','2019-01')),
                           labels = c('1999', '2019')) +
        facet_grid(cols = vars(species),
                   rows = vars(litter)) +
        theme_classic() +
        theme(text = element_text(size = 22),
              strip.background = element_blank(),
              strip.text = element_text(size = 22),
              axis.title = element_blank())
    } else if (fID < 10) {
      p = ggplot(plotDM %>% subset(litter == plotLT & species == plotSp), 
                 aes(x = date)) +
        geom_line(aes(y = avg), color = base.Colors[['LightGrey']]) +
        # geom_ribbon(aes(ymin = lower, ymax = upper),
        #             fill = base.Colors[['LightGrey']], alpha = .5) +
        geom_smooth(aes(y = avg, linetype = sig), method = 'lm', se = F,
                    color = base.Colors[['Blue']]) +
        scale_linetype_identity() +
        scale_x_continuous(breaks = ym(c('1999-01','2019-01')),
                           labels = c('1999', '2019')) +
        theme_classic() +
        theme(text = element_text(size = 22),
              strip.background = element_blank(),
              strip.text = element_text(size = 22),
              axis.title = element_blank())
    } else {
      p = ggplot(plotDM %>% subset(litter == plotLT & species == plotSp) %>%
                   mutate(litter = markLT), 
                 aes(x = date)) +
        geom_line(aes(y = avg), color = base.Colors[['LightGrey']]) +
        # geom_ribbon(aes(ymin = lower, ymax = upper),
        #             fill = base.Colors[['LightGrey']], alpha = .5) +
        geom_smooth(aes(y = avg, linetype = sig), method = 'lm', se = F,
                    color = base.Colors[['Blue']]) +
        scale_linetype_identity() +
        scale_x_continuous(breaks = ym(c('1999-01','2019-01')),
                           labels = c('1999', '2019')) +
        facet_grid(rows = vars(litter)) +
        theme_classic() +
        theme(text = element_text(size = 22),
              strip.background = element_blank(),
              strip.text = element_text(size = 22),
              axis.title = element_blank())
    }
    return(p)
  })
}
ggarrange(plotLitterSubs[[1]], plotLitterSubs[[2]],
          plotLitterSubs[[3]], plotLitterSubs[[4]],
          plotLitterSubs[[5]], plotLitterSubs[[6]],
          plotLitterSubs[[7]], plotLitterSubs[[8]],
          plotLitterSubs[[9]], plotLitterSubs[[10]],
          nrow = 2, ncol = 5,
          labels = c(LETTERS[1:5], '', '', '', '', ''), 
          heights = c(0.53, 0.47),
          widths = c(0.195, 0.195, 0.195, 0.195, 0.22)) %>%
  annotate_figure(bottom = text_grob("Time", size = 22),
                  left = text_grob(expression("Litter Biomass ("~g/m^2~")"),
                                   size = 22, rot = 90)) %>%
  ggexport(filename = paste(plotPath, "Figure1.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 90 * 25.3999,  # mm to inch
           res = 300)

rm(plotDM, plotLitterSubs)

# **********************************************
# Figure 2 - Seasonal pattern of litterfall ----
# **********************************************
# Require: tp
plotBaseTP = tp %>% mutate(Month = month(date)) %>% subset(litter != 'TDM') %>%
  group_by(litter, Month, species) %>%
  summarise(value = mean(DM), se = sd(DM)/sqrt(n())) %>%
  mutate(lower = value - se, upper = value + se,
         xMonth = ifelse(species == '28a Sa', Month - 0.1,
                         ifelse(species == '26a Sc', Month - 0.05,
                                ifelse(species == '83a Ko', Month + 0.05,
                                       ifelse(species == '83a Am', 
                                              Month + 0.1, Month)))),
         species = factor(species,
                          levels = c('28a Sa', '26a Sc', '28a Ko',
                                     '83a Ko', '83a Am')),
         litter = factor(litter,
                         levels = c('LDM', 'FDM'),
                         labels = c('Leaves', 'Reproductive Parts')))
plotTpSubs = lapply(1:2, function(id){
  if (id == 1) {
    param = 'Leaves'
  } else if (id == 2) {
    param = 'Reproductive Parts'
  }
  p = ggplot(plotBaseTP %>% subset(litter == param),
             aes(x = xMonth, y = value, color = species)) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    scale_x_continuous(breaks = 1:12) +
    ggalt::geom_xspline() +
    ggsci::scale_color_npg(labels = plotSpeciesLabels) +
    labs(color = 'Assemblages') +
    facet_grid(cols = vars(litter)) +
    theme_classic() +
    theme(text = element_text(size = 22),
          strip.background = element_blank(),
          strip.text = element_text(size = 22),
          axis.title = element_blank(),
          legend.position = 'top') +
    guides(color = guide_legend(nrow = 2))
  return(p)
})
ggarrange(plotTpSubs[[1]], plotTpSubs[[2]],
          labels = LETTERS[1:2], nrow = 1, ncol = 2,
          common.legend = T) %>%
  annotate_figure(left = text_grob(expression('Litter Biomass ('~g/m^2~')'),
                                   size = 22, rot = 90),
                  bottom = text_grob("Months", size = 22)) %>%
  ggexport(filename = paste(plotPath, "Figure2.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 90 * 25.3999,  # mm to inch
           res = 300)

rm(plotBaseTP, plotTpSubs)

# *********************************************
# Figure 3 - Annual litterfall comparisons ----
# *********************************************
# Require: annLit, annLitSum, onewayTbl
# Assign significance mark
annSig = lapply(c("FDM", "LDM"), function(lt){
  baseTbl = annLitSum %>% subset(litter == lt) %>% arrange(desc(biomass))
  if (lt == 'LDM') {
    baseTbl$mark = c('a', 'b', 'c', 'c', 'd')
  } else {
    baseTbl$mark = c('a', 'b', 'bc', 'c', 'd')
  }
  # calculate coordinate of label
  maxDM = max(annLit$annSum[which(annLit$litter == lt)])
  minDM = min(annLit$annSum[which(annLit$litter == lt)])
  baseTbl$ylabel = maxDM + 0.2 * (maxDM - minDM)
  return(baseTbl)
}) %>% do.call(rbind.data.frame, args = .) %>%
  dplyr::select(species, litter, ylabel, mark) %>%
  as.data.frame()
# join table for plotting
plotTbl = annLit %>%
  mutate(DM = log(annSum),
         species = factor(species,
                          levels = c('28a Sa', '26a Sc', 
                                     '28a Ko', '83a Ko', '83a Am')),
         litter = factor(litter,
                         levels = c('FDM', 'LDM', 'TDM'),
                         labels = c('Reproductive Parts',
                                    'Leaves', 'Twigs')))
annSig = annSig %>%
  mutate(ylabel = log(ylabel),
         species = factor(species,
                          levels = c('28a Sa', '26a Sc', 
                                     '28a Ko', '83a Ko', '83a Am')),
         litter = factor(litter,
                         levels = c('FDM', 'LDM', 'TDM'),
                         labels = c('Reproductive Parts',
                                    'Leaves', 'Twigs')))
annPlotSpeciesLabels = c(expression(atop(NA,
                                         atop(textstyle(italic('S. apetala')), 
                                              '(Exotic, 28a)'))),
                         expression(atop(NA,
                                         atop(textstyle(italic('S. caseolaris')), 
                                              '(Exotic, 26a)'))),
                         expression(atop(NA,
                                         atop(textstyle(italic('K. obovata')), 
                                              '(Native, 28a)'))),
                         expression(atop(NA,
                                         atop(textstyle(italic('K. obovata')), 
                                              '(Native, 83a)'))),
                         expression(atop(NA,
                                         atop(textstyle(italic('A. marina')), 
                                              '(Native, 83a)'))))

annPlotList = lapply(1:2, function(id){
  if (id == 1) plotLT = 'Reproductive Parts' else plotLT = 'Leaves'
  p = ggplot(plotTbl %>% subset(litter == plotLT)) +
    geom_jitter(aes(x = species, y = DM), 
                color = base.Colors[['LightGrey']], alpha = 0.3) +
    geom_boxplot(aes(x = species, y = DM),
                 color = base.Colors[['Blue']], fill = NA) +
    geom_text(data = annSig %>% subset(litter == plotLT),
              aes(x = species, y = ylabel, label = mark), inherit.aes = F,
              size = 5) +
    scale_x_discrete(labels = annPlotSpeciesLabels) +
    labs(x = '', y = 'Log-Transformed Annual Litter Biomass') +
    facet_grid(cols = vars(litter)) +
    theme_classic() +
    theme(text = element_text(size = 22),
          strip.background = element_blank(),
          strip.text = element_text(size = 22),
          axis.text.x = element_text(hjust = 1, vjust = 1, angle = 30,
                                     margin = margin(b = -5)),
          axis.title = element_blank())
  return(p)
})
ggarrange(annPlotList[[2]], annPlotList[[1]],
          labels = LETTERS[1:2], nrow = 1, ncol = 2) %>%
  annotate_figure(left = text_grob('Log-Transformed Annual Litterfall',
                                   size = 22, rot = 90)) %>%
  ggexport(filename = paste(plotPath, "Figure3.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 90 * 25.3999,  # mm to inch
           res = 300)

rm(annSig, plotTbl, annPlotSpeciesLabels, annPlotList)

# ****************************************************
# Figure 4 S7-8 - Seasonal patterns of effect sizes ----
# ****************************************************
#
load(file = paste(dataPath, "/SMap/scneExpFull.Rdata", sep = ''))
# Driver significance from CCM models
# significance are illustrated by point shapes: 16 = significant, 1 = ns
edmDriverSign = bind_rows(
  # 83a Am
  data.frame(species = '83a Am',
             litter = rep(c('LDM', 'FDM'), each = 6),
             var = rep(varList,2),
             pointShape = c(16, 16, 1, 16, 1, 1, 
                            16, 16, 1, 16, 1, 1)),
  # 83a Ko
  data.frame(species = '83a Ko',
             litter = rep(c('LDM', 'FDM'), each = 6),
             var = rep(varList,2),
             pointShape = c(16, 16, 16, 16, 1, 1,
                            16, 16, 16, 16, 16, 16)),
  # 28a Ko
  data.frame(species = '28a Ko',
             litter = rep(c('LDM', 'FDM'), each = 6),
             var = rep(varList,2),
             pointShape = c(16, 1, 1, 1, 16, 1,
                            16, 16, 1, 16, 16, 1)),
  # 28a Sa
  data.frame(species = '28a Sa',
             litter = rep(c('LDM', 'FDM'), each = 6),
             var = rep(varList,2),
             pointShape = c(17, 17, 17, 17, 2, 2,
                            2, 2, 2, 2, 17, 17)),
  # 26a Sc
  data.frame(species = '26a Sc',
             litter = rep(c('LDM', 'FDM'), each = 6),
             var = rep(varList,2),
             pointShape = c(2, 2, 17, 2, 2, 2,
                            17, 2, 2, 17, 2, 2))
)
# Summary by months
sensSum = edmScenExp %>% mutate(Month = month(time)) %>% 
  subset(!is.na(sens)) %>%
  group_by(species, litter, var, Month) %>% 
  summarise(avg = mean(sens), sd = sd(sens), count = n(),
            margin = qt(0.975, df = count-1)*sd/sqrt(count)) %>% 
  as.data.frame() %>% 
  mutate(upper = avg + margin, lower = avg - margin)  %>%
  left_join(edmDriverSign, by = c('species', 'litter', 'var')) %>%
  mutate(var = factor(var,
                      levels = c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP'),
                      labels = c('Temperature', 'Rainfall', 'Sea Level',
                                 'Salinity', 'Nitrogen', 'Phosphorus')),
         species = factor(species, 
                          levels = c('28a Sa', '26a Sc',
                                     '28a Ko', '83a Ko', '83a Am')),
         litter = factor(litter,
                         levels = c('LDM', 'FDM'),
                         labels = c('Leaves', 'Reproductive Parts')))
edmDriverSign = edmDriverSign %>% 
  mutate(lineType = ifelse(pointShape %in% c(16, 17),
                           'solid', 'dashed'))
# Summary by months
sensPlain = edmScenExp %>% mutate(Month = month(time)) %>% 
  subset(!is.na(sens)) %>%
  group_by(species, litter, var, mode, Month) %>% 
  summarise(pred = mean(Y1), predMargin = qt(0.975, df = n()-1)*sd(Y1)/sqrt(n()),
            ori = mean(Y0), oriMargin = qt(0.975, df = n()-1)*sd(Y0)/sqrt(n())) %>% 
  mutate(predUp = pred + predMargin, predLow = pred - predMargin,
         oriUp = ori + oriMargin, oriLow = ori - oriMargin)  %>%
  as.data.frame() %>% 
  dplyr::select(species, litter, var, mode, Month, 
                pred, predUp, predLow, ori, oriUp, oriLow) %>%
  left_join(edmDriverSign, by = c('species', 'litter', 'var')) %>%
  mutate(var = factor(var,
                      levels = c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP'),
                      labels = c('Temperature', 'Rainfall', 'Sea Level',
                                 'Salinity', 'Nitrogen', 'Phosphorus')),
         species = factor(species, 
                          levels = c('28a Sa', '26a Sc',
                                     '28a Ko', '83a Ko', '83a Am')),
         litter = factor(litter,
                         levels = c('LDM', 'FDM'),
                         labels = c('Leaves', 'Reproductive Parts')),
         mode = factor(mode, 
                       levels = c('plus', 'minus'),
                       labels = c('Increase', 'Decrease')))

# Plot by months
ltFig4 = 'Leaves'; ltFigS8 = 'Reproductive Parts'
spList = levels(sensPlain$species)

## Figure 4 ----
plotPlainSubsF = lapply(1:30, function(id){
  if (id%%5==0) {
    plotSp = spList[5]
    plotV = levels(sensPlain$var)[(id - id%%5)/5]
    p = ggplot(sensPlain %>% 
                 subset(species == plotSp & litter == ltFig4 & var == plotV)) + 
      geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
                 linetype = 'dashed') +
      geom_ribbon(aes(x = Month, ymin = predLow, ymax = predUp, fill = mode),
                  alpha = .5) +
      # geom_line(aes(x = Month, y = pred, color = mode)) + 
      geom_line(aes(x = Month, y = ori, linetype = lineType), 
                color = 'black') +
      scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
      # ggsci::scale_color_npg() +
      ggsci::scale_fill_npg() +
      scale_linetype_identity() +
      labs(fill = 'Change in Target Driver') +
      facet_grid(cols = vars(species),
                 rows = vars(var)) +
      theme_classic() +
      theme(text = element_text(size = 22),
            axis.title = element_blank(),
            strip.text = element_text(size = 22),
            strip.background = element_blank(),
            legend.position = 'bottom')
  } else {
    plotSp = spList[id%%5]
    plotV = levels(sensPlain$var)[(id - id%%5)/5+1]
    p = ggplot(sensPlain %>% 
                 subset(species == plotSp & litter == ltFig4 & var == plotV)) + 
      geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
                 linetype = 'dashed') +
      geom_ribbon(aes(x = Month, ymin = predLow, ymax = predUp, fill = mode),
                  alpha = .5) +
      # geom_line(aes(x = Month, y = pred, color = mode)) + 
      geom_line(aes(x = Month, y = ori, linetype = lineType), 
                color = 'black') +
      scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
      # ggsci::scale_color_npg() +
      ggsci::scale_fill_npg() +
      scale_linetype_identity() +
      labs(fill = 'Change in Target Driver') +
      facet_grid(cols = vars(species)) +
      theme_classic() +
      theme(text = element_text(size = 22),
            axis.title = element_blank(),
            strip.text = element_text(size = 22),
            strip.background = element_blank(),
            legend.position = 'bottom')
  }
  return(p)
})

# Temperature on leaf litter
ggarrange(plotPlainSubsF[[1]] + 
            scale_x_continuous(breaks = 1:12) +
            labs(fill = 'Temperature'), 
          plotPlainSubsF[[3]] + 
            scale_x_continuous(breaks = 1:12) +
            labs(fill = 'Temperature'), 
          plotPlainSubsF[[4]] + 
            scale_x_continuous(breaks = 1:12) +
            labs(fill = 'Temperature'), 
          plotPlainSubsF[[5]] + 
            scale_x_continuous(breaks = 1:12) +
            labs(fill = 'Temperature') + 
            theme(strip.text.y = element_blank()),
          labels = LETTERS[1:4], nrow = 2, ncol = 2,
          common.legend = T) %>%
  annotate_figure(bottom = text_grob("Month", size = 22),
                  left = text_grob(expression("Predicted"~Delta*"Biomass (g/m"^2*")"), 
                                   size = 22, rot = 90)) %>%
  ggexport(filename = paste(plotPath, "Figure4.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 120 * 25.3999,  # mm to inch
           res = 300)

rm(plotPlainSubsF, ltFig4)
rm(edmDriverSign, sensSum, sensPlain, edmScenExp)

## Figure S7 ----
plotScneExp = edmScenExp %>%
  subset(!is.na(sens)) %>%
  group_by(species, litter, var) %>%
  summarise(avg = mean(sens), sd = sd(sens), count = n(),
            margin = qt(0.975, df = count - 1) * sd / sqrt(count)) %>%
  as.data.frame() %>%
  mutate(upper = avg + margin, lower = avg - margin) %>%
  left_join(edmDriverSign, by = c('species', 'litter', 'var')) %>%
  mutate(var = factor(var,
                      levels = c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP'),
                      labels = c('Temperature', 'Rainfall', 'Sea Level',
                                 'Salinity', 'Nitrogen', 'Phosphorus')),
         species = factor(species, 
                          levels = c('28a Sa', '26a Sc',
                                     '28a Ko', '83a Ko', '83a Am')),
         litter = factor(litter,
                         levels = c('LDM', 'FDM'),
                         labels = c('Leaves', 'Reproductive Parts')))
plotSESubs = lapply(levels(plotScneExp$litter), function(lt){
  plotTbl = plotScneExp %>% subset(litter == lt)
  p = ggplot(plotTbl, aes(x = species, y = avg, color = species)) +
    geom_hline(yintercept = 0, 
               color = base.Colors[['LightGrey']],
               linetype = 'dashed') +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    geom_point(aes(shape = pointShape),
               size = 2) +
    scale_x_discrete(limits = rev) +
    scale_shape_identity() +
    ggsci::scale_color_npg() +
    labs(x = 'Mangrove Assemblages', 
         y = 'Effect Size') +
    facet_grid(cols = vars(var),
               rows = vars(litter))  +
    theme_classic() +
    theme(text = element_text(size = 22),
          strip.background = element_blank(),
          strip.text = element_text(size = 22),
          axis.title = element_blank(),
          legend.position = 'none')
  if (lt == 'Leaves') {
    p = p + 
      scale_y_continuous(breaks = seq(-50, 50, 50)) + 
      coord_flip()
  } else {
    p = p + coord_flip()
  }
  return(p)
})

ggarrange(plotSESubs[[1]], plotSESubs[[2]],
          labels = LETTERS[1:2], nrow = 2, ncol = 1) %>%
  annotate_figure(bottom = text_grob("Effect Size", size = 18)) %>%
  ggexport(filename = paste(plotPath, "FigureS7.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 90 * 25.3999,  # mm to inch
           res = 300)

rm(plotScneExp, plotSESubs)
rm(edmDriverSign, sensSum, sensPlain, edmScenExp)

## Figure S8 ----
plotPlainSubsL = lapply(1:30, function(id){
  if (id%%5==0) {
    plotSp = spList[5]
    plotV = levels(sensPlain$var)[(id - id%%5)/5]
    p = ggplot(sensPlain %>% 
                 subset(species == plotSp & litter == ltFigS8 & var == plotV)) + 
      geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
                 linetype = 'dashed') +
      geom_ribbon(aes(x = Month, ymin = predLow, ymax = predUp, fill = mode),
                  alpha = .5) +
      # geom_line(aes(x = Month, y = pred, color = mode)) + 
      geom_line(aes(x = Month, y = ori, linetype = lineType), 
                color = 'black') +
      scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
      # ggsci::scale_color_npg() +
      ggsci::scale_fill_npg() +
      scale_linetype_identity() +
      labs(fill = 'Change in Target Driver') +
      facet_grid(cols = vars(species),
                 rows = vars(var)) +
      theme_classic() +
      theme(text = element_text(size = 22),
            axis.title = element_blank(),
            strip.text = element_text(size = 22),
            strip.background = element_blank(),
            legend.position = 'bottom')
  } else {
    plotSp = spList[id%%5]
    plotV = levels(sensPlain$var)[(id - id%%5)/5+1]
    p = ggplot(sensPlain %>% 
                 subset(species == plotSp & litter == ltFigS8 & var == plotV)) + 
      geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
                 linetype = 'dashed') +
      geom_ribbon(aes(x = Month, ymin = predLow, ymax = predUp, fill = mode),
                  alpha = .5) +
      # geom_line(aes(x = Month, y = pred, color = mode)) + 
      geom_line(aes(x = Month, y = ori, linetype = lineType), 
                color = 'black') +
      scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
      # ggsci::scale_color_npg() +
      ggsci::scale_fill_npg() +
      scale_linetype_identity() +
      labs(fill = 'Change in Target Driver') +
      facet_grid(cols = vars(species)) +
      theme_classic() +
      theme(text = element_text(size = 22),
            axis.title = element_blank(),
            strip.text = element_text(size = 22),
            strip.background = element_blank(),
            legend.position = 'bottom')
  }
  return(p)
})
# Nitrogen on reproduction
ggarrange(plotPlainSubsL[[21]] + 
            scale_x_continuous(breaks = 1:12) +
            labs(fill = 'Coastal Nitrogen'), 
          plotPlainSubsL[[23]] + 
            scale_x_continuous(breaks = 1:12) +
            labs(fill = 'Coastal Nitrogen'), 
          plotPlainSubsL[[24]] + 
            scale_x_continuous(breaks = 1:12) +
            labs(fill = 'Coastal Nitrogen'),
          labels = LETTERS[1:3], nrow = 1, ncol = 3,
          common.legend = T) %>%
  annotate_figure(bottom = text_grob("Month", size = 22),
                  left = text_grob(expression("Predicted"~Delta*"Biomass (g/m"^2*")"), 
                                   size = 22, rot = 90)) %>%
  ggexport(filename = paste(plotPath, "FigureS8.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 63 * 25.3999,  # mm to inch
           res = 300)

rm(plotPlainSubsL, ltFigS8)
rm(edmDriverSign, sensSum, sensPlain, edmScenExp)

# *************************************************
# Figure S2 - Temporal patterns of environment ----
# *************************************************
# require: wt
plotBaseTW = wt %>% mutate(Month = month(Date)) %>%
  subset(Date >= ym('1999-01') & Date <= ym('2019-12')) %>%
  dplyr::select(Month, AVGT, PRCP, SL, SAL, TN, TP) %>%
  tidyr::gather(key = 'Var', value = 'Value', AVGT:TP) %>%
  group_by(Var) %>%
  mutate(Value = GGally::range01(Value)) %>%
  group_by(Month, Var) %>%
  summarise(value = mean(Value), 
            se = sd(Value)/sqrt(n())) %>%
  mutate(lower = value - se, upper = value + se,
         Var = factor(Var,
                      levels = c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP'),
                      labels = c('Temperature', 'Rainfall', 'Sea Level',
                                 'Salinity', 'Nitrogen', 'Phosphorus')))
plotTW = ggplot() +
  geom_col(data = plotBaseTW %>% subset(Var == 'Rainfall'), 
           aes(x = as.factor(Month), y = value, fill = Var),
           color = base.Colors[['Blue']]) +
  ggalt::geom_xspline(data = plotBaseTW %>% subset(Var != 'Rainfall'), 
                      aes(x = Month, y = value, color = Var),
                      size = 1) +
  scale_fill_manual(values = base.Colors[['White']]) +
  ggsci::scale_color_d3() +
  labs(x = 'Month', y = 'Transformed Mean') +
  theme_classic() +
  theme(text = element_text(size = 22),
        axis.title.x = element_text(margin = margin(t = 10)),
        # axis.title.y = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())
ggsave(filename = "FigureS2.tiff",
       plot = plotTW,
       device = "tiff",
       path = plotPath,
       scale = 1.7,
       width = 180,
       height = 90,
       units = "mm",
       dpi = 300,
       limitsize = F)
dev.off()

rm(plotBaseTW, plotTW)

# ***************************************************
# Figure S4 - Temporal variations of environment ----
# ***************************************************
# require: wt
{
  plotTbl = wt %>%
    dplyr::select(Date, all_of(varList)) %>%
    tidyr::gather(key = 'var', value = 'value', AVGT:TP) %>%
    subset(Date >= ym('1999-01') & Date <= ym('2019-12')) %>%
    left_join(wtMovTrend %>% dplyr::select(var, slope, sig),
              by = c('var')) %>%
    mutate(var = factor(var, 
                        levels = varList,
                        labels = c('Temperature', 'Rainfall', 'Sea Level',
                                   'Salinity', 'Nitrogen', 'Phosphorus')))
  plotWtSubs = lapply(levels(plotTbl$var), function(plotV){
    p = ggplot(plotTbl %>% subset(var == plotV), aes(x = Date)) +
      geom_line(aes(y = value), color = base.Colors[['LightGrey']]) +
      geom_smooth(aes(y = value, linetype = sig), method = 'lm', se = F,
                  color = base.Colors[['Blue']]) +
      scale_linetype_identity() +
      scale_x_continuous(breaks = ym(c('1999-01','2019-01')),
                         labels = c('1999', '2019')) +
      facet_grid(cols = vars(var)) +
      theme_classic() +
      theme(text = element_text(size = 22),
            strip.background = element_blank(),
            axis.title = element_blank())
    return(p)
  })
}
ggarrange(plotWtSubs[[1]], plotWtSubs[[2]], plotWtSubs[[3]],
          plotWtSubs[[4]], plotWtSubs[[5]], plotWtSubs[[6]],
          labels = LETTERS[1:6], nrow = 2, ncol = 3) %>%
  annotate_figure(bottom = text_grob('Time', size = 22),
                  left = text_grob('Values', size = 22, rot = 90)) %>%
  ggexport(filename = paste(plotPath, "FigureS4.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 90 * 25.3999,  # mm to inch
           res = 300)

rm(plotTbl, plotWtSubs)

#
# Figure S5 - Seasonal pattern of diff litterfall ----
#
# require: difLit
plotDif = difLit %>%
  mutate(Month = month(time)) %>% 
  group_by(litter, Month, species) %>%
  summarise(value = mean(DM), se = sd(DM)/sqrt(n())) %>%
  mutate(lower = value - se, upper = value + se,
         xMonth = ifelse(species == '28a Sa', Month - 0.1,
                         ifelse(species == '26a Sc', Month - 0.05,
                                ifelse(species == '83a Ko', Month + 0.05,
                                       ifelse(species == '83a Am', 
                                              Month + 0.1, Month)))),
         species = factor(species,
                          levels = c('28a Sa', '26a Sc', '28a Ko',
                                     '83a Ko', '83a Am')),
         litter = factor(litter,
                         levels = c('LDM', 'FDM'),
                         labels = c('Leaves', 'Reproductive Parts')))
plotDifSubs = lapply(1:2, function(id){
  if (id == 1) {
    param = 'Leaves'
  } else if (id == 2) {
    param = 'Reproductive Parts'
  }
  p = ggplot(plotDif %>% subset(litter == param),
             aes(x = xMonth, y = value, color = species)) +
    geom_hline(yintercept = 0, linetype = 'dashed', 
               color = base.Colors[['LightGrey']]) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    scale_x_continuous(breaks = 1:12) +
    geom_line() +
    ggsci::scale_color_npg(labels = plotSpeciesLabels) +
    labs(color = 'Assemblages') +
    facet_grid(cols = vars(litter)) +
    theme_classic() +
    theme(text = element_text(size = 22),
          strip.background = element_blank(),
          axis.title = element_blank(),
          legend.position = 'top') +
    guides(color = guide_legend(nrow = 2))
  return(p)
})
# left: (a) Monthly LDM, right: (b) Monthly FDM
ggarrange(plotDifSubs[[1]], plotDifSubs[[2]],
          labels = LETTERS[1:2], nrow = 1, ncol = 2,
          common.legend = T) %>%
  annotate_figure(left = text_grob(expression('Monthly Increment in Litter Biomass ('~g/m^2~')'),
                                   size = 22, rot = 90),
                  bottom = text_grob("Months", size = 22)) %>%
  ggexport(filename = paste(plotPath, "FigureS5.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 90 * 25.3999,  # mm to inch
           res = 300)

rm(plotDif, plotDifSubs)

#
# Figure S6 - Phenophases ----
#
# require: tp
{
  # Set up phase table and assign seasons to each month
  phaseSubTbl = tp %>% 
    group_by(trapID, litter) %>%
    mutate(DM = GGally::range01(DM),
           Month = month(date)) %>%
    ungroup() %>%
    group_by(Month, species, trapID, litter) %>%
    summarise(DM = mean(DM)) %>%
    ungroup() %>%
    group_by(Month, species, litter) %>%
    summarise(Mean = mean(DM), SE = sd(DM)/sqrt(n())) %>%
    mutate(upper = Mean + SE, lower = Mean - SE) %>%
    as.data.frame()
  
  phaseTbl = full_join(phaseSubTbl %>% subset(litter == 'FDM') %>%
                         select(species, Month, Mean, upper, lower),
                       phaseSubTbl %>% subset(litter == 'LDM') %>%
                         select(species, Month, Mean, upper, lower),
                       by = c('species', 'Month'),
                       suffix = c('F', 'L')) %>%
    select(species, Month, MeanF, upperF, lowerF, MeanL, upperL, lowerL) %>%
    mutate(Season = ifelse(Month %in% 3:5, 'Spring',
                           ifelse(Month %in% 6:8, 'Summer',
                                  ifelse(Month %in% 9:11, 'Autumn', 'Winter'))),
           Season = factor(Season, levels = c('Spring', 'Summer', 'Autumn', 'Winter'))) %>%
    arrange(Month)
  
  ## Subplots
  plotList = lapply(1:6, function(SPID){
    if (SPID == 1) {  # 28a Sa
      pointTbl = phaseTbl %>% subset(species == '28a Sa') %>% 
        mutate(species = as.character(species),
               species = 'Exotic: 28a Sa')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.0396 + 0.025, 0), 
               labelY = c(0.00926 + 0.015, 0)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    } else if (SPID == 2) {  # 26a Sc
      pointTbl = phaseTbl %>% subset(species == '26a Sc') %>% 
        mutate(species = as.character(species),
               species = 'Exotic: 26a Sc')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.299 + 0.015, 0.206 - 0.015), 
               labelY = c(0.139 + 0.015, 0.163 - 0.015)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    } else if (SPID == 3) {  # 28a Ko
      pointTbl = phaseTbl %>% subset(species == '28a Ko') %>% 
        mutate(species = as.character(species),
               species = 'Native: 28a Ko')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.117 + 0.015, 0.127 - 0.015), 
               labelY = c(0.046 + 0.01, 0.0234 - 0.015)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    } else if (SPID == 4) {  # 83a Ko
      pointTbl = phaseTbl %>% subset(species == '83a Ko') %>% 
        mutate(species = as.character(species),
               species = 'Native: 83a Ko')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.165 + 0.015, 0.154 - 0.015), 
               labelY = c(0.0324 + 0.01, 0.0186 - 0.02)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    } else {  # 83a Am
      pointTbl = phaseTbl %>% subset(species == '83a Am') %>% 
        mutate(species = as.character(species),
               species = 'Native: 83a Am')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.0926 + 0.015, 0.1364 + 0.02), 
               labelY = c(0.00755 + 0.005, 0.00337 + 0.01)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    }
    ## PLOTTING
    ## SPID == 6 is specially for legend
    if (SPID %in% 1:5) {
      p = ggplot(pointTbl, aes(x = MeanL, y = MeanF)) + 
        geom_path(data = segmentTbl, aes(color = Season), 
                  linewidth = 1, alpha = .7) + 
        ## Line ranges of F and L
        geom_linerange(aes(xmin = lowerL, xmax = upperL)) +
        geom_linerange(aes(ymin = lowerF, ymax = upperF)) +
        ## Labels
        geom_text(data = labelTbl, aes(label = Month, x = labelX, y = labelY), 
                  check_overlap = T, size = 6.5) +
        ## Assign special point marks and colors
        geom_point(data = pointTbl %>% subset(Month == 1),
                   color = base.Colors[['Purple']], size = 5, shape = 12) +
        geom_point(data = pointTbl %>% subset(Month == 12),
                   color = base.Colors[['Purple']], size = 5, shape = 5) +
        scale_color_manual(values = c(base.Colors[['Pink']], 
                                      base.Colors[['Green']], 
                                      base.Colors[['Orange']], 
                                      base.Colors[['LightGrey']])) +
        scale_y_continuous(breaks = seq(0, 0.5, 0.1)) +
        scale_x_continuous(breaks = seq(0, 0.5, 0.1)) +
        facet_wrap(~species) +
        theme_classic() +
        theme(text = element_text(size = 22),
              strip.text = element_text(size = 22),
              strip.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              panel.background = element_rect(color = 'black'),
              legend.position = 'none')
    } else { ## Get legend
      p = ggplot(pointTbl, aes(x = MeanL, y = MeanF)) + 
        geom_path(data = segmentTbl, aes(color = Season), 
                  linewidth = 1, alpha = .7) +
        scale_color_manual(values = c(base.Colors[['Pink']], 
                                      base.Colors[['Green']], 
                                      base.Colors[['Orange']], 
                                      base.Colors[['LightGrey']])) +
        labs(color = 'Seasons') +
        theme_classic() +
        theme(text = element_text(size = 22),
              legend.text = element_text(size = 22),
              legend.spacing.y = unit(0.5, 'cm')) +
        guides(color = guide_legend(byrow = T))
      p = get_legend(p)
    }
    return(p)
  })
}
# Plot
ggarrange(plotList[[1]], plotList[[2]], plotList[[6]],
          plotList[[3]], plotList[[4]], plotList[[5]], 
          labels = c(LETTERS[1:2], '', LETTERS[3:5])) %>%
  annotate_figure(left = text_grob('Reproductive Litterfall', 
                                   size = 22, rot = 90),
                  bottom = text_grob('Leaf Litterfall', size = 22)) %>%
  ggexport(filename = paste(plotPath, "FigureS6.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 18*7 * 25.3999,  # mm to inch
           res = 300)

rm(phaseTbl, phaseSubTbl, plotList)

#
# Figure S9 - Temporal pattern of annual litterfall ----
#
# require: annLit
plotTbl = annLit %>% subset(litter != 'TDM') %>%
  group_by(species, trapID, litter) %>%
  mutate(value = scale(annSum)[,1]) %>%
  group_by(species, litter, Year) %>%
  summarise(avg = mean(value), se = sd(value)/sqrt(n())) %>%
  ungroup() %>% 
  mutate(upper = avg + se, lower = avg - se,
         xYear = ifelse(species == '28a Sa', Year - 0.1,
                        ifelse(species == '26a Sc', Year - 0.05,
                               ifelse(species == '83a Ko', Year + 0.05,
                                      ifelse(species == '83a Am', 
                                             Year + 0.1, Year)))),
         species = factor(species,
                          levels = c('28a Sa', '26a Sc', '28a Ko',
                                     '83a Ko', '83a Am')))

scaledSubPlots = lapply(1:2, function(id){
  if (id == 1) {plotLT = 'LDM'; mark = 'Leaves'} else {plotLT = 'FDM'; mark = 'Reproductive Parts'}
  p = ggplot(plotTbl %>% subset(litter == plotLT) %>% mutate(mark = mark),
             aes(x = xYear, y = avg, color = species)) +
    geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line() +
    geom_point(aes(shape = species), fill = 'white') +
    scale_shape_manual(values = c(24, 24, 21, 21, 21),
                       labels = plotSpeciesLabels) +
    ggsci::scale_color_npg(labels = plotSpeciesLabels) +
    labs(color = 'Assemblages', shape = 'Assemblages') +
    facet_grid(rows = vars(mark)) +
    theme_classic() +
    theme(text = element_text(size = 22),
          strip.background = element_blank(),
          axis.title = element_blank(),
          legend.position = 'top')
  return(p)
})

ggarrange(scaledSubPlots[[1]], scaledSubPlots[[2]],
          nrow = 2, ncol = 1, labels = LETTERS[1:2],
          common.legend = T) %>%
  annotate_figure(bottom = text_grob("Year", size = 22),
                  left = text_grob("Scaled Litter Biomass", 
                                   size = 22, rot = 90)) %>%
  ggexport(filename = paste(plotPath, "FigureS9.tiff", sep = "/"),
           width = 180 * 25.3999,  # mm to inch
           height = 90 * 25.3999,  # mm to inch
           res = 300)

rm(scaledSubPlots, plotLT)

#
# Figure S10 - Bivariate correlation of annual proxies ----
#
# require: wt
# summarization of Environment
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
  data.frame() %>% mutate(Year = 1999:2019)
corMatrix = annEnv %>% 
  dplyr::select(Tmean, Rsum, Lmean, Smean, Nmean, Pmean) %>%
  cor()
corPlot = ggcorrplot::ggcorrplot(corMatrix,
                                 type = 'lower', 
                                 method = 'square',
                                 hc.order = T,
                                 outline.color = "white",
                                 ggtheme = theme_bw(),
                                 lab = T,
                                 lab_size = 6,
                                 tl.cex = 18) +
  theme(legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.key.width = unit(5, 'mm'),
        legend.key.height = unit(20, 'mm'))
ggsave(filename = "FigureS10.tiff",
       plot = corPlot,
       device = "tiff",
       path = plotPath,
       scale = 2,
       width = 85,
       height = 80,
       units = "mm",
       dpi = 300,
       limitsize = F)
dev.off()

