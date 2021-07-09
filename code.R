library(ape)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(nortest)
library(car)

library(cvms)
library(broom)    # tidy()
library(tibble)   # tibble()

mi_IM <- function(xy, resp, cov, trt, blk){
  d = as.matrix(dist(xy))
  di = 1/d
  diag(di) = 0
  W = di
  
  mod1 = aov(resp ~ trt)
  mod2 = aov(resp ~ blk + trt)
  mod3 = aov(resp ~ cov + trt)
  mod4 = aov(resp ~ cov + blk + trt)
  
  smod1 = summary(mod1)
  smod2 = summary(mod2)
  smod3 = summary(mod3)
  smod4 = summary(mod4)
  
  # MORAN INDEX
  IMy = Moran.I(resp, W)
  IM1 = Moran.I(mod1$residuals, W)
  IM2 = Moran.I(mod2$residuals, W)
  IM3 = Moran.I(mod3$residuals, W)
  IM4 = Moran.I(mod4$residuals, W)
  
  # NORMALITY
  # Shapiro-Wilk, Shapiro Francia, Anderson-Darlin
  nt_sh_1 = shapiro.test(mod1$residuals)
  nt_sh_2 = shapiro.test(mod2$residuals)
  nt_sf_1 = sf.test(mod1$residuals)
  nt_sf_2 = sf.test(mod2$residuals)
  nt_ad_1 = ad.test(mod1$residuals)
  nt_ad_2 = ad.test(mod2$residuals)
  
  # HOMOSCEDASTICITY
  # Barlett, Levine y Fligner-Killeen
  ht_bt_1 = bartlett.test(mod1$residuals, trt)
  ht_bt_2 = bartlett.test(mod2$residuals, trt)
  ht_lv_1 = leveneTest(mod1$residuals, trt)
  ht_lv_2 = leveneTest(mod2$residuals, trt)
  ht_fl_1 = fligner.test(mod1$residuals, trt)
  ht_fl_2 = fligner.test(mod2$residuals, trt)
  
  # RETURN
  return(list(
    pvIMy = IMy$p.value,
    pvIM1 = IM1$p.value,
    pvIM2 = IM2$p.value,
    pvIM3 = IM3$p.value,
    pvIM4 = IM4$p.value,
    IMy = IMy$observed,
    IM1 = IM1$observed,
    IM2 = IM2$observed,
    IM3 = IM3$observed,
    IM4 = IM4$observed,
    H = smod2[[1]][1,4],
    pvmod2_blk = smod2[[1]][1,5],
    
    pvsh1 = ifelse(!nt_sh_1$p.value<0.05,'Normal','Non-normal'),
    pvsh2 = ifelse(!nt_sh_2$p.value<0.05,'Normal','Non-normal'),
    pvsf1 = ifelse(!nt_sf_1$p.value<0.05,'Normal','Non-normal'),
    pvsf2 = ifelse(!nt_sf_2$p.value<0.05,'Normal','Non-normal'),
    pvad1 = ifelse(!nt_ad_1$p.value<0.05,'Normal','Non-normal'),
    pvad2 = ifelse(!nt_ad_2$p.value<0.05,'Normal','Non-normal'),
    
    pvbt1 = ifelse(!ht_bt_1$p.value<0.05,'Homoscedasticity','Heteroscedasticity'),
    pvbt2 = ifelse(!ht_bt_2$p.value<0.05,'Homoscedasticity','Heteroscedasticity'),
    pvlv1 = ifelse(!ht_lv_1$`Pr(>F)`[1]<0.05,'Homoscedasticity','Heteroscedasticity'),
    pvlv2 = ifelse(!ht_lv_2$`Pr(>F)`[1]<0.05,'Homoscedasticity','Heteroscedasticity'),
    pvfl1 = ifelse(!ht_fl_1$p.value<0.05,'Homoscedasticity','Heteroscedasticity'),
    pvfl2 = ifelse(!ht_fl_2$p.value<0.05,'Homoscedasticity','Heteroscedasticity')))
}


# GENERACION DE DATOS -----------------------------------------------------

set.seed(123)
rep_clorofila = replicate(1000, sort.int(rnorm(30, 300, 30), 1))

densidad = factor(as.vector(replicate(10, sample(c('s1','s2','s3'), 3, F))))
bloq = factor(rep(paste0('b',1:10), each = 3))
table(densidad, bloq)

xy = expand.grid(x=1:3, y = 1:10)

# COVARIABLE --------------------------------------------------------------

#################################
set.seed(2020)
covar=matrix(NA,nrow = 10, ncol = 3)
for(k in 1:3){
  for(i in 1:10){
    covar[i,k]=(runif(1,i/5,i/5+rnorm(1,0.75,0.01)))
  }
}
covar

covarb = as.vector(t(covar))
df = data.frame(xy, densidad, bloq, covarb, rep_clorofila)

summary(aov(df$X1~df$bloq+df$densidad))
#########################

df %>% 
  ggplot()+
  aes(x,y, fill = covarb)+
  geom_tile()

# DESIGN
df %>% 
  ggplot()+
  aes(x,y, fill = X2)+
  geom_tile()+
  geom_text(aes(label = interaction(densidad, bloq)),
            color = 'white')+
  labs(fill = 'Response', 
       x='X-coordinate', 
       y='Y-coordinate')+
  theme_bw()+
  geom_segment(x = 3.6, xend = 3.6,
               y = 1, yend = 10,
               arrow = arrow(length = unit(5, "pt")))+
  annotate(geom = 'text', label=c('+','-'),
           x = c(3.6,3.6), y = c(10.5,0.8), size=c(6,10))

# ANALISIS ----------------------------------------------------------------
pre_mod = aov(rep_clorofila[,1] ~ bloq + densidad)
sum_pre_mod = summary(pre_mod)
sum_pre_mod

df[,1:6]
rep_IM_pre = apply(rep_clorofila, 2, function(y){
  mi_IM(xy = xy, resp = y, cov = covarb, trt=densidad, blk=bloq)
})

rep_IM = do.call(rbind, lapply(rep_IM_pre, as.data.frame))
rep_IM[1:10,]
class(rep_IM)

rep_IM %>% 
  ggplot()+
  aes(pvIM3, pvIM2)+
  geom_point(color=gray(0.3), size = 0.6)+
  geom_hline(yintercept = 0.05, color = 'black', linetype='dotted')+
  geom_vline(xintercept = 0.05, color = 'black', linetype='dotted')+
  theme_bw()

table('CRD'=rep_IM$pvIM1<0.05, 'RCB'=rep_IM$pvIM2<0.05)
table('CRD'=rep_IM$pvIM1<0.05, 'COV'=rep_IM$pvIM3<0.05)
table('CRD'=rep_IM$pvIM1<0.05, 'COV-RCB'=rep_IM$pvIM4<0.05)
table('RCB'=rep_IM$pvIM2<0.05, 'COV'=rep_IM$pvIM3<0.05)


rep_IM %>% 
  ggplot()+
  aes(pvIM2, pvIM1)+
  geom_point(shape = 21, color='white', fill = '#08306B', size = 2)+
  geom_hline(yintercept = 0.05, color = '#6BAED6', size = 0.6)+
  geom_vline(xintercept = 0.05, color = '#6BAED6', size = 0.6)+
  labs(x='MI (p-value) - RCB', y='MI (p-value) - CRD')+
  annotate(geom = 'text', label=c('A','B','C','D'),
           x=c(0,0.5,0.5,0), y=c(-0.04,-0.04, 0.5,0.5))+
  theme_bw()

table(pvIM1 = rep_IM$pvIM1<0.05, 
      pvIM2 = rep_IM$pvIM2<0.05)
pv = (rep_IM[,c('pvIM1','pvIM2')]<0.05)
pv[pv==T] = 'Spatial D'
pv[pv==F] = 'Spatial I'
tbl_pv = table(CRD = pv[,1], RCB = pv[,2])
cfm <- tidy(tbl_pv)
plot_confusion_matrix(cfm, 
                      target_col = "RCB", 
                      prediction_col = "CRD",
                      counts_col = "n", add_sums = T,
                      palette = 'Blues',
                      rotate_y_text = F)+
  labs(x='RCB',y='CRD')+
  ggsave(filename = 'table_1.jpeg', device = 'jpeg',
         path = fig_path,
         scale = 1.2, width = 7, height = 7, units = 'cm')

f3a = rep_IM %>%
  ggplot()+
  aes(pvIM1, pvIMy)+
  geom_point(shape = 21, color='white', fill = '#08306B', size = 2)+
  geom_hline(yintercept = 0.05, color = '#6BAED6', size = 0.6)+
  geom_vline(xintercept = 0.05, color = '#6BAED6', size = 0.6)+
  labs(x='MI (p-value) - CRD', y='MI (p-value) - Response')+
  annotate(geom = 'text', label=c('A','B','C','D'),
           x=c(0,0.2,0.2,0), y=c(-0.04,-0.04, 0.5,0.5))+
  theme_bw()
table(pvIMy = rep_IM$pvIMy<0.05,
      pvIM1 = rep_IM$pvIM1<0.05)

f3b = rep_IM %>%
  ggplot()+
  aes(pvIM2,pvIMy)+
  geom_point(shape = 21, color='white', fill = '#08306B', size = 2)+
  geom_hline(yintercept = 0.05, color = '#6BAED6', size = 0.6)+
  geom_vline(xintercept = 0.05, color = '#6BAED6', size = 0.6)+
  labs(x='MI (p-value) - CRD', y='MI (p-value) - Response')+
  annotate(geom = 'text', label=c('A','B','C','D'),
           x=c(0,0.5,0.5,0), y=c(-0.04,-0.04, 0.5,0.5))+
  theme_bw()
table(pvIMy = rep_IM$pvIMy<0.05,
      pvIM2 = rep_IM$pvIM2<0.05)

f3 = arrangeGrob(f3a+annotate(geom='label',x = 1, y = 1, label='a'),
                 f3b+annotate(geom='label',x = 1, y = 1, label='b'), 
                 nrow=1)
grid.arrange(f3)

# TABLAS DE NORMALIDAD ----------------------------------------------------

tbl2a = table(CRD = rep_IM$pvsh1, RCB = rep_IM$pvsh2) %>% 
  tidy() %>% 
  plot_confusion_matrix(target_col = "RCB", 
                        prediction_col = "CRD",
                        counts_col = "n", add_sums = T,
                        palette = 'Blues',
                        rotate_y_text = F)+
  labs(x='RCB',y='CRD')

tbl2b = table(CRD = rep_IM$pvsf1, RCB = rep_IM$pvsf2) %>% 
  tidy() %>% 
  plot_confusion_matrix(target_col = "RCB", 
                        prediction_col = "CRD",
                        counts_col = "n", add_sums = T,
                        palette = 'Blues',
                        rotate_y_text = F)+
  labs(x='RCB',y='CRD')

tbl2c = table(CRD = rep_IM$pvad1, RCB = rep_IM$pvad2) %>% 
  tidy() %>% 
  plot_confusion_matrix(target_col = "RCB", 
                        prediction_col = "CRD",
                        counts_col = "n", add_sums = T,
                        palette = 'Blues',
                        rotate_y_text = F)+
  labs(x='RCB',y='CRD')

tbl2 = arrangeGrob(tbl2a+
                     annotate(geom='label',x = Inf, y = Inf, label='S-W', 
                              hjust = 1, vjust = -1)+
                     theme(plot.margin = unit(rep(0.1,4), "lines"))+
                     coord_cartesian(clip = "off"),
                   tbl2b+
                     annotate(geom='label',x = Inf, y = Inf, label='S-F', 
                              hjust = 1, vjust = -1)+
                     theme(plot.margin = unit(rep(0.1,4), "lines"))+
                     coord_cartesian(clip = "off"),
                   nrow=1)
grid.arrange(tbl2)

# HOMOSCEDASTICITY TABLES ----------------------------------------------------

tbl3a = table(CRD = rep_IM$pvbt1, RCB = rep_IM$pvbt2) %>% 
  tidy() %>% 
  plot_confusion_matrix(target_col = "RCB", 
                        prediction_col = "CRD",
                        counts_col = "n", add_sums = T,
                        palette = 'Blues',
                        rotate_y_text = F)+
  labs(x='RCB',y='CRD')

tbl3b = table(CRD = rep_IM$pvlv1, RCB = rep_IM$pvlv2) %>% 
  tidy() %>% 
  plot_confusion_matrix(target_col = "RCB", 
                        prediction_col = "CRD",
                        counts_col = "n", add_sums = T,
                        palette = 'Blues',
                        rotate_y_text = F)+
  labs(x='RCB',y='CRD')

tbl3c = table(CRD = rep_IM$pvfl1, RCB = rep_IM$pvfl2) %>% 
  tidy() %>% 
  plot_confusion_matrix(target_col = "RCB", 
                        prediction_col = "CRD",
                        counts_col = "n", add_sums = T,
                        palette = 'Blues',
                        rotate_y_text = F)+
  labs(x='RCB',y='CRD')

tbl3 = arrangeGrob(tbl3a+
                     annotate(geom='label',x = Inf, y = Inf, label='B', 
                              hjust = 1, vjust = -1)+
                     theme(plot.margin = unit(rep(0.1,4), "lines"))+
                     coord_cartesian(clip = "off"),
                   tbl3c+
                     annotate(geom='label',x = Inf, y = Inf, label='F-K', 
                              hjust = 1, vjust = -1)+
                     theme(plot.margin = unit(rep(0.1,4), "lines"))+
                     coord_cartesian(clip = "off"), 
                   nrow=1)
grid.arrange(tbl3)


# BLOCK EFFICIENCY ---------------------------------------------------
tbl_h_blk = as.vector(table('H' = rep_IM$H>1, 'pvmod2_blk' = rep_IM$pvmod2_blk<0.05))

table(rep_IM$H>1, rep_IM$pvmod2_blk<0.05,
      dnn = list('H', 'pvmod2_blk'))

rep_IM %>%
  ggplot()+
  aes(H, pvmod2_blk)+
  geom_point(shape = 21, color='white', fill = '#08306B', size = 2)+
  geom_hline(yintercept = 0.05, color = '#6BAED6', size = 0.6)+
  geom_vline(xintercept = 0.05, color = '#6BAED6', size = 0.6)+
  labs(x='MI (p-value) - CRD', y='MI (p-value) - Response')+
  labs(x='H', y='Block p-value - RCB')+
  annotate(geom = 'text', label=tbl_h_blk,
           x=c(-3,20,-3, 20), y=c(0.25, 0.25, -0.05, -0.05))+
  theme_bw()

rep_IM %>%
  ggplot()+
  aes(H, pvIM2)+
  geom_point(shape = 21, color='white', fill = '#08306B', size = 2)+
  geom_hline(yintercept = 0.05, color = '#6BAED6', size = 0.6)+
  geom_vline(xintercept = 0.05, color = '#6BAED6', size = 0.6)+
  labs(x='H', y='MI (p-value) - RCB')+
  theme_bw()
