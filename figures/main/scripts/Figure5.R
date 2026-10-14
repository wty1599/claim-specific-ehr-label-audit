script_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1]), winslash='/'))
source(file.path(script_dir, 'common.R'))

# Presentation-only rendering from the archived aggregate inputs.
read_fig5 <- function(name) read_source(file.path('fig5', name))
rate <- read_fig5('Figure6_panelB_treatment_prevalence_v2.csv')
ess <- read_fig5('arm_ess_concentration_by_target.csv')
balance <- read_fig5('covariate_balance_fixed_scale.csv')
flow <- read_fig5('Figure6_panelA_landmark_flow_v2.csv')
rate <- rate[order(rate$phenotype), ]
ess <- ess[ess$target == 'ATE', ]
ess <- ess[order(ess$phenotype, -ess$treatment), ]
c2 <- balance[balance$phenotype == 'C2 lower-risk' &
              balance$target == 'ATO' &
              balance$variable != 'urine_output_24h_ml', ]
c2 <- c2[order(abs(c2$smd_fixed_unweighted_scale), decreasing=TRUE), ]
stopifnot(identical(as.integer(rate$n_treated), c(265L,49L)),
          identical(as.integer(rate$n), c(3151L,14374L)),
          identical(round(100*rate$treatment_prevalence,2), c(8.41,.34)),
          identical(round(ess$arm_ess,1), c(54.9,2279.3,23.9,14324.5)),
          sum(rate$n)==tail(flow$n_remaining,1), nrow(c2)==15L,
          sum(abs(c2$smd_fixed_unweighted_scale)>.1)==8L,
          c2$variable[1]=='creatinine_max',
          round(abs(c2$smd_fixed_unweighted_scale[1]),3)==.927)

labels <- c(age='Age', sofa_score='SOFA score', aki_stage_0_24h='AKI stage',
            sapsii='SAPS II', creatinine_max='Creatinine, maximum',
            bun_max='BUN, maximum', potassium_max='Potassium, maximum',
            bicarbonate_min='Bicarbonate, minimum', ph_min='pH, minimum',
            lactate_max='Lactate, maximum', mbp_min='Mean BP, minimum',
            heart_rate_max='Heart rate, maximum',
            resp_rate_max='Respiratory rate, maximum', gcs_min='GCS, minimum',
            gender_male='Male sex')
stopifnot(setequal(names(labels), c2$variable))

rate$group <- factor(sub(' .*','',rate$phenotype),
                     levels=rev(sub(' .*','',rate$phenotype)))
rate$label <- sprintf('%s%% (%s/%s)',fmt_f(100*rate$treatment_prevalence,2),
                      fmt_n(rate$n_treated),fmt_n(rate$n))
p_a <- ggplot(rate,aes(100*treatment_prevalence,group)) +
  geom_col(width=.23,fill=BLUE) +
  geom_text(aes(x=0,label=label),hjust=0,nudge_y=.25,
            size=FS/ggplot2::.pt) +
  scale_x_continuous(limits=c(0,10),breaks=c(0,5,10),
                     labels=function(x) paste0(x,'%'),
                     expand=expansion(mult=c(0,.04))) +
  scale_y_discrete(expand=expansion(add=c(.5,.65))) +
  labs(x='RRT record proportion',y=NULL) +
  theme(plot.margin=margin(4,8,3,4))

ess$row <- rev(seq_len(nrow(ess)))
ess$row_label <- paste(sub(' .*','',ess$phenotype),
                       ifelse(ess$treatment==1,'RRT+','RRT-'))
counts <- rbind(data.frame(row=ess$row,value=ess$n,Series='Raw count'),
                data.frame(row=ess$row,value=ess$arm_ess,Series='ATE ESS'))
counts$Series <- factor(counts$Series,levels=c('Raw count','ATE ESS'))
p_b <- ggplot(ess,aes(y=row)) +
  geom_segment(aes(x=arm_ess,xend=n,yend=row),linewidth=LW,colour=GREY) +
  geom_point(data=counts,aes(x=value,shape=Series,colour=Series),
             size=2,stroke=.45) +
  scale_shape_manual(values=c(1,16)) +
  scale_colour_manual(values=c(INK,BLUE)) +
  scale_x_log10(limits=c(10,20000),breaks=c(10,100,1000,10000),labels=fmt_n) +
  scale_y_continuous(limits=c(.35,5.25),breaks=ess$row,
                     labels=ess$row_label,expand=c(0,0)) +
  labs(x='Raw count or ESS',y=NULL) +
  theme(legend.position='none',axis.ticks.y=element_blank(),
        plot.margin=margin(4,4,3,2))
p_btext <- ggplot(ess,aes(y=row)) +
  geom_text(aes(x=.02,label=fmt_n(n)),hjust=0,size=FS/ggplot2::.pt) +
  geom_text(aes(x=.57,label=fmt_f(arm_ess)),hjust=0,size=FS/ggplot2::.pt) +
  annotate('text',x=c(.02,.57),y=4.85,
           label=c('Raw n','ATE ESS'),hjust=0,family=FONT,size=FS/ggplot2::.pt) +
  scale_x_continuous(limits=c(0,1.15),expand=c(0,0)) +
  scale_y_continuous(limits=c(.35,5.25),expand=c(0,0)) +
  labs(x=' ',y=NULL) +
  theme_classic(base_size=FS,base_family=FONT) +
  theme(axis.text=element_blank(),axis.ticks=element_blank(),
        axis.line=element_blank(),plot.margin=margin(4,3,3,0),
        axis.title.x=element_text(colour='white'))
p_both <- p_b + p_btext + patchwork::plot_layout(widths=c(1.8,1))

balance_long <- rbind(
  data.frame(variable=c2$variable,value=abs(c2$unweighted_smd_fixed_scale),
             series='Unweighted'),
  data.frame(variable=c2$variable,value=abs(c2$smd_fixed_unweighted_scale),
             series='Overlap weighted'))
balance_long$label <- factor(unname(labels[balance_long$variable]),
                             levels=rev(unname(labels[c2$variable])))
balance_long$series <- factor(balance_long$series,
                              levels=c('Unweighted','Overlap weighted'))
p_c <- ggplot(balance_long,aes(value,label,shape=series,colour=series)) +
  geom_vline(xintercept=.1,linetype='dashed',colour=GREY,linewidth=LW) +
  geom_point(position=position_dodge(width=.42),size=2.2,stroke=.45) +
  annotate('text',x=.13,y=16,label='SMD = 0.1',hjust=0,
           family=FONT,size=FS/ggplot2::.pt,colour=INK) +
  scale_shape_manual(values=c(1,16)) +
  scale_colour_manual(values=c(GREY,BLUE)) +
  scale_x_continuous(limits=c(0,1.4),breaks=c(0,.1,.5,1.0),
                     expand=expansion(mult=c(0,.02))) +
  coord_cartesian(ylim=c(.5,16.4),clip='off') +
  labs(x='Absolute standardized mean difference',y=NULL) +
  theme(legend.position='top',legend.title=element_blank(),
        plot.margin=margin(4,7,3,4),axis.ticks.y=element_blank())

draw <- function() {
  grid.newpage()
  heading <- function(letter,text,x,y) {
    grid.text(letter,x=unit(x,'mm'),y=unit(y,'mm'),just=c('left','top'),
              gp=gpar(fontfamily=FONT,fontsize=TAG,fontface='bold',col=INK))
    grid.text(text,x=unit(x+6,'mm'),y=unit(y-.5,'mm'),just=c('left','top'),
              gp=gpar(fontfamily=FONT,fontsize=FS,col=INK))
  }
  heading('A','RRT records',2,158)
  heading('B','Effective sample size',62,158)
  heading('C','C2 covariate balance',2,103)
  print(p_a,newpage=FALSE,vp=viewport(x=unit(1,'mm'),y=unit(106,'mm'),
        width=unit(59,'mm'),height=unit(46,'mm'),just=c('left','bottom')))
  print(p_both,newpage=FALSE,vp=viewport(x=unit(61,'mm'),y=unit(106,'mm'),
        width=unit(112,'mm'),height=unit(46,'mm'),just=c('left','bottom')))
  for (i in seq_len(2)) {
    x <- c(88,116)[i]
    grid.points(x=unit(x,'mm'),y=unit(149.5,'mm'),pch=c(1,16)[i],
                size=unit(1.8,'mm'),gp=gpar(col=c(INK,BLUE)[i],lwd=.8))
    grid.text(c('Raw count','ATE ESS')[i],x=unit(x+3,'mm'),y=unit(149.5,'mm'),
              just='left',gp=gpar(fontfamily=FONT,fontsize=FS,col=INK))
  }
  print(p_c,newpage=FALSE,vp=viewport(x=unit(1,'mm'),y=unit(1,'mm'),
        width=unit(172,'mm'),height=unit(98,'mm'),just=c('left','bottom')))
}
export_figure(draw,'Figure5',160)
write.csv(balance_long,file.path(QA,'Figure5C_displayed_balance.csv'),row.names=FALSE)
write.csv(ess,file.path(QA,'Figure5B_displayed_ESS.csv'),row.names=FALSE)
pass('Figure5',c('Source aggregates are unchanged.',
 'Both RRT record proportions and all four arm-specific ESS values match the source.',
 'Fifteen C2 covariates displayed; eight overlap-weighted |SMD| values exceed 0.1.',
 'The dashed 0.1 line is labeled within panel C.'))
