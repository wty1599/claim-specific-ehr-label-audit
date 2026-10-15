script_dir <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1]), winslash='/'))
source(file.path(script_dir, 'common.R'))
# Keep the raster device local to this renderer.
export_figure <- function(plot,name,height_mm,width_mm=SPEC$width_mm) {
  draw <- function() { if(is.function(plot)) plot() else if(inherits(plot,'grob')) {grid.newpage();grid.draw(plot)} else print(plot) }
  grDevices::cairo_pdf(file.path(OUT,paste0(name,'.pdf')),width=width_mm/25.4,height=height_mm/25.4,family=FONT,onefile=TRUE)
  draw();dev.off()
  ragg::agg_tiff(file.path(OUT,paste0(name,'.tiff')),width=width_mm,height=height_mm,units='mm',res=600,compression='lzw',background='white')
  draw();dev.off()
  ragg::agg_png(file.path(QA,paste0(name,'.png')),width=width_mm,height=height_mm,units='mm',res=180,background='white')
  draw();dev.off()
  writeLines(c(paste('Width_mm',width_mm),paste('Height_mm',height_mm),paste('Body_pt',FS),paste('Tag_pt',TAG)),file.path(QA,paste0(name,'_dimensions.txt')))
}
# R's graphics lwd unit is 1/96 inch; allow for conversion in the PDF device.
LW <- max(LW, .28)
theme_update(axis.line=element_line(linewidth=LW), axis.ticks=element_line(linewidth=LW))

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
threshold <- .1
stopifnot(identical(as.integer(rate$n_treated), c(265L,49L)),
          identical(as.integer(rate$n), c(3151L,14374L)),
          identical(round(100*rate$treatment_prevalence,2), c(8.41,.34)),
          identical(round(ess$arm_ess,1), c(54.9,2279.3,23.9,14324.5)),
          all(rate$n_treated+rate$n_control==rate$n),
          sum(rate$n)==tail(flow$n_remaining,1),
          all(abs(rate$treatment_prevalence-rate$n_treated/rate$n)<1e-14),
          nrow(c2)==15L, !anyDuplicated(c2$variable),
          all(is.finite(c2$smd_fixed_unweighted_scale)),
          sum(abs(c2$smd_fixed_unweighted_scale)>threshold)==8L,
          c2$variable[1]=='creatinine_max',
          round(abs(c2$smd_fixed_unweighted_scale[1]),3)==.927,
          SPEC$width_mm==174, FS==8, TAG==10, LW>=.2, FONT=='Arial')
for (i in seq_len(nrow(rate))) {
  z <- ess[ess$phenotype==rate$phenotype[i], ]
  stopifnot(z$n[z$treatment==1]==rate$n_treated[i],
            z$n[z$treatment==0]==rate$n_control[i])
}

labels <- c(age='Age', sofa_score='SOFA score', aki_stage_0_24h='AKI stage',
            sapsii='SAPS II', creatinine_max='Creatinine, maximum',
            bun_max='BUN, maximum', potassium_max='Potassium, maximum',
            bicarbonate_min='Bicarbonate, minimum', ph_min='pH, minimum',
            lactate_max='Lactate, maximum', mbp_min='Mean BP, minimum',
            heart_rate_max='Heart rate, maximum',
            resp_rate_max='Respiratory rate, maximum', gcs_min='GCS, minimum',
            gender_male='Male sex')
stopifnot(setequal(names(labels), c2$variable))

# Supply existing aggregate SMDs directly; do not call a balance estimator.
B <- data.frame(Type=ifelse(c2$variable=='gender_male','Binary','Contin.'),
                Diff.Un=c2$unweighted_smd_fixed_scale,
                Diff.Adj=c2$smd_fixed_unweighted_scale,
                row.names=unname(labels[c2$variable]))
bt <- structure(list(Balance=B), class=c('bal.tab.bin','bal.tab'))
attr(bt,'print.options') <- list(type='bin', stats='mean.diffs', abs=FALSE,
                                 continuous='std', binary='std', weight.names='Adj',
                                 nweights=1L, un=TRUE, quick=FALSE)
p_c <- cobalt::love.plot(bt, stats='mean.diffs', abs=TRUE, var.order=NULL,
                         thresholds=c(m=threshold), colors=c(GREY,BLUE), shapes=c(1,16),
                         size=2.2, line=FALSE, sample.names=c('Unweighted','Overlap weighted'),
                         position='top', title=NULL, wrap=60) +
  labs(x='Absolute standardized mean difference',y=NULL) +
  theme_classic(base_size=FS,base_family=FONT) +
  theme(text=element_text(colour=INK),axis.text=element_text(size=FS,colour=INK),
        axis.line=element_line(linewidth=LW),axis.ticks=element_line(linewidth=LW),
        legend.position='top',legend.title=element_blank(),legend.text=element_text(size=FS),
        legend.margin=margin(0,0,0,0),plot.margin=margin(0,7,3,4)) +
  annotate('text',x=.12,y=15.5,label='0.1 threshold',hjust=0,vjust=1,
           colour=INK,family=FONT,size=FS/ggplot2::.pt)
for (i in seq_along(p_c$layers)) {
  if (inherits(p_c$layers[[i]]$geom,'GeomVline')) p_c$layers[[i]]$aes_params$linewidth <- LW
  if (inherits(p_c$layers[[i]]$geom,'GeomPoint')) {
    p_c$layers[[i]]$mapping$stroke <- NULL
    p_c$layers[[i]]$aes_params$stroke <- .45
  }
}
stopifnot(identical(tail(levels(p_c$data$var),1),unname(labels[c2$variable[1]])))
for (i in seq_len(nrow(c2))) {
  z <- p_c$data[as.character(p_c$data$var)==labels[c2$variable[i]], ]
  stopifnot(z$stat[z$Sample=='Unweighted']==abs(c2$unweighted_smd_fixed_scale[i]),
            z$stat[z$Sample=='Overlap weighted']==abs(c2$smd_fixed_unweighted_scale[i]))
}
point_layer <- which(vapply(p_c$layers,function(x) inherits(x$geom,'GeomPoint'),logical(1)))
rendered_points <- ggplot_build(p_c)$data[[point_layer]]
stopifnot(nrow(rendered_points)==2L*nrow(c2),all(rendered_points$stroke>0),
          sum(rendered_points$shape==1)==nrow(c2),sum(rendered_points$shape==16)==nrow(c2),
          identical(sort(rendered_points$x[rendered_points$shape==1]),sort(abs(c2$unweighted_smd_fixed_scale))),
          identical(sort(rendered_points$x[rendered_points$shape==16]),sort(abs(c2$smd_fixed_unweighted_scale))),
          round(max(rendered_points$x[rendered_points$shape==1]),4)==1.2098)

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
stopifnot(file.copy(file.path(QA,'Figure5.png'),file.path(OUT,'Figure5.png'),overwrite=TRUE))
write.csv(c2,file.path(QA,'Figure5C_displayed_balance.csv'),row.names=FALSE)
write.csv(p_c$data,file.path(QA,'figure5_cobalt_plotted_values.csv'),row.names=FALSE)
write.csv(rendered_points,file.path(QA,'figure5_cobalt_rendered_points.csv'),row.names=FALSE)
write.csv(ess,file.path(QA,'Figure5B_displayed_ESS.csv'),row.names=FALSE)
