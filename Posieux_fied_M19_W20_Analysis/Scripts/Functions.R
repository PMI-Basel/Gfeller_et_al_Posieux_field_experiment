
# plot soil benzoxazinoids at the end of maize growth
ggplot_soil_box <- function(data, x, y, fill, facet) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    geom_boxplot(outlier.colour = NA, alpha = 0.3, 
                 color = "grey68", show.legend = FALSE) + 
    geom_quasirandom(size = 1, shape = 1, width = 0.15, groupOnX = TRUE,
                     color = "black", alpha = 0.6, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", 
                 position = position_dodge(width = 0.75), show.legend = FALSE) +
    stat_summary(fun.data = mean_se, geom = "errorbar", size = 1, width = 0.25,
                 position = position_dodge(width = 0.75), show.legend = FALSE)  +
    scale_fill_manual(values = c(W = "gold2", b = "darkgreen"))  +
    scale_x_discrete(name = "Maize genotype",
                     labels = c(W = "WT", b = expression(italic(bx1)))) +
    facet_wrap(vars({{facet}}), nrow = 3, scales = "free")
}


# plot soil benzoxazinoids at wheat emergence
ggplot_soil_box_res <- function(data, x, y, fill, facet) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    geom_boxplot(outlier.colour = NA, alpha = 0.3, 
                 color = "grey68", show.legend = FALSE) + 
    geom_quasirandom(size = 1, shape = 1, width = 0.15, groupOnX = TRUE,
                     color = "black", alpha = 0.6, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", 
                 position = position_dodge(width = 0.75), show.legend = FALSE) +
    stat_summary(fun.data = mean_se, geom = "errorbar", size = 1, width = 0.25,
                 position = position_dodge(width = 0.75), show.legend = FALSE)  +
    scale_fill_manual(values = c(W = "gold2", b = "darkgreen"))  +
    scale_x_discrete(name = "Soil conditioning",
                     labels = c(W = "WT", b = expression(italic(bx1)))) +
    facet_wrap(vars({{facet}}), nrow = 3, scales = "free")
  
}


# plot model assumptions for linear models
plot.mod.vg <- function(mod){
  par(mfrow = c(1,2))
  plot1 <- plot(fitted(mod), resid(mod), xlab = "Fitted values", ylab = "Residuals", main = "Tukey-Anscombe plot")
  plot2 <- car::qqPlot(resid(mod), dist = "norm", mean = mean(resid(mod)), sd = sd(resid(mod)),xlab = "Theoretical quantiles", ylab = "Empirical quantiles", main = "Q-Q plot of residuals")
}


# tidy ANOVA table and create label to print ANOVA in ggplot as text for lm and gls models
#get table: t.anova_seed <- create_anova_table(model)[[1]]
#get label: lab_anova_seed <- create_anova_table(model)[[2]]
create_anova_table <- function(mod){
  t.anova <- Anova(mod) %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()
  
  if(colnames(t.anova)[4] == "Pr(>Chisq)") {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:wheat_var", "C x V:"), 
             Variable = str_replace_all(Variable, "trt", "Cond:"), 
             Variable = str_replace_all(Variable, "wheat_var", "Var:"),
             Variable = str_replace_all(Variable, "pos_y", "Pos:"),
             p_new = if_else(`Pr(>Chisq)` < 0.001, paste0("p < 0.001 "), 
                             paste0("p = ", sprintf("%.3f", `Pr(>Chisq)`))),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., "")))
    
    t.anova_red <- t.anova
    
    an_title <- "**ANOVA**<Br>"
  } else {
    t.anova <- t.anova %>% 
      mutate(Variable = str_replace_all(Variable, "trt:wheat_var", "C x V:"),
             Variable = str_replace_all(Variable, "trt", "Cond:"),
             Variable = str_replace_all(Variable, "wheat_var", "Var:"),
             Variable = str_replace_all(Variable, "pos_y", "Pos:"),
             p_new = if_else(`Pr(>F)` < 0.001, paste0("p < 0.001 "), 
                             paste0("p = ", sprintf("%.3f", `Pr(>F)`))),
             across(where(is.numeric), ~ as.character(signif(., 2))),
             across(everything(), ~ replace_na(., ""))) 
    
    t.anova_red <- t.anova %>% 
      filter(Variable != "Residuals")
    
    an_title <- "**ANOVA**<Br>"
  }
  
  if(dim(t.anova_red)[1] == 3){
    lab_anova <- paste(an_title, 
                       t.anova_red$Variable[1], t.anova_red$p_new[1], "<Br>", 
                       t.anova_red$Variable[2], t.anova_red$p_new[2], "<Br>",
                       t.anova_red$Variable[3], t.anova_red$p_new[3])
  } else {
    lab_anova <- paste(an_title, 
                       t.anova_red$Variable[1], t.anova_red$p_new[1], "<Br>", 
                       t.anova_red$Variable[2], t.anova_red$p_new[2], "<Br>",
                       t.anova_red$Variable[3], t.anova_red$p_new[3], "<Br>",
                       t.anova_red$Variable[4], t.anova_red$p_new[4])
  }
  
  an.list <-  list(t.anova, lab_anova)
  return(an.list)
}



# plot soil benzoxazinoids in the response phase (wheat)
ggplot_resp_box_bx <- function(data, x, y, fill, t.emmeans = TRUE, ...) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    geom_boxplot(outlier.colour = NA, alpha = 0.3, 
                 color = "grey68") + 
    geom_quasirandom(size = 1, shape = 1, width = 0.08, dodge.width = 0.75,
                     color = "black", alpha = 0.6, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point",shape = 18, size = 3, color = "black", 
                 position = position_dodge(width = 0.75), show.legend = FALSE) +
    stat_summary(fun.data = mean_se, geom = "errorbar", size = 1, width = 0.25,
                 position = position_dodge(width = 0.75), show.legend = FALSE)  +
    scale_fill_manual(name = "Soil conditioning", 
                      labels = c(W = "WT", b = expression(italic(bx1))), 
                      values = c(W = "gold2", b = "darkgreen"))  +
    scale_x_discrete(labels = c(claro = "Claro", fiorina = "Fiorina", 
                                sailor = "Sailor")) +
    xlab("Wheat variety") +
    {if(t.emmeans)geom_text(data = tab.emm, 
                            aes(y = max.y, group = 0.5,
                                label = if_else(p.value < 0.001, paste0("p < 0.001 "), 
                                                paste0("p = ", sprintf("%.3f", p.value)))),
                            size = 3)} +
    theme(legend.position = "right") 
  
}

# plot soil wheat phenotypes in the response phase
ggplot_resp_box <- function(data, x, y, fill, t.emmeans = TRUE, ...) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    geom_boxplot(outlier.colour = NA, alpha = 0.3, 
                 color = "grey68") + 
    geom_quasirandom(size = 1, shape = 1, width = 0.08, dodge.width = 0.75,
                     color = "black", alpha = 0.6, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point",shape = 18, size = 3, color = "black", 
                 position = position_dodge(width = 0.75), show.legend = FALSE) +
    stat_summary(fun.data = mean_se, geom = "errorbar", size = 1, width = 0.25,
                 position = position_dodge(width = 0.75), show.legend = FALSE)  +
    scale_fill_manual(name = "Soil conditioning", 
                      labels = c(W = "WT", b = expression(italic(bx1))), 
                      values = c(W = "gold2", b = "darkgreen"))  +
    scale_x_discrete(labels = c(claro = "Claro", fiorina = "Fiorina", 
                                sailor = "Sailor")) +
    xlab("Wheat variety") +
    {if(t.emmeans)geom_text(data = tab.emm, 
                            aes(y = max.y, group = 0.5,
                                label = if_else(p.value < 0.001, paste0("p < 0.001 "), 
                                                paste0("p = ", sprintf("%.3f", p.value)))))} +
    theme(legend.position = "top") 
  
}


# Tidy ANOVA tables for element plot; label in one line
output_anova_ele <- function(data) {
  t.anova <- data %>% 
    rownames_to_column(var = "Variable") %>%
    tibble()%>%
    mutate(Variable = str_replace_all(Variable, "trt:wheat_var", "C x V:"), 
           Variable = str_replace_all(Variable, "trt", "Cond:"),
           Variable = str_replace_all(Variable, "wheat_var", "Var:"),
           p_new = if_else(`Pr(>Chisq)` < 0.001, paste0("p < 0.001 "), 
                           paste0("p = ", sprintf("%.2f", `Pr(>Chisq)`))),
           across(where(is.numeric), ~ as.character(signif(., 2))),
           across(everything(), ~ replace_na(., ""))) 
  
  # tidy anova table for plot
  lab_anova <- paste("**ANOVA:**", t.anova$Variable[1], " ", t.anova$p_new[1], ",", 
                     t.anova$Variable[2], " ", t.anova$p_new[2], ",",
                     t.anova$Variable[3], " ", t.anova$p_new[3])
  return(lab_anova)
}

