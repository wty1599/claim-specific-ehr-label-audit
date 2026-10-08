script_file <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1]), winslash='/')
source(file.path(dirname(script_file), 'common.R'))
if (.Platform$OS.type=='windows') invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
# R device line units render at roughly 0.75 of the nominal ggplot width.
LW <- max(LW,.28)
theme_update(axis.line=element_line(linewidth=LW),axis.ticks=element_line(linewidth=LW))
src <- function(f) figure_reader(file.path(DATA,'fig23',f),check.names=FALSE,
                            stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM')
sep <- src('Figure4_panelA_separation_v2.csv')
archived <- src('Figure4_panelA_shape_tiles_v2.csv')
tile <- src('Figure3_panelB_shape_summary.csv')
d2 <- src('Figure4_panelB_delta_auc_v3_lineage_resolved.csv')
oracle <- src('Figure4_panelC_oracle_v2.csv')
bridge <- src('Table_D1_D2_bridge_method_corrected_v2.csv')
fixed <- src('SF4_locked_variants.csv')
stopifnot(SPEC$width_mm == 174, FS == 8, TAG == 10, LW >= .2,
 round(sep$estimate[sep$evidence_type=='Reference'],6)==2.716950,
 round(sep$estimate[sep$evidence_type=='Separation alert'],6)==3.521073,
 nrow(tile)==9L, !anyNA(tile$median_p),
 all(tile$median_p>=0 & tile$median_p<=1),
 all(tile$shape_alert_threshold==.05),
 !anyDuplicated(tile[c('condition','n_eval')]),
 all(tile$n==100L),nrow(d2)==8L,nrow(oracle)==12L)

# This is a presentation summary of the held-out shape diagnostics.
archived$stage[archived$stage=='G7 discrete control'] <- 'G2 discrete control'
idx <- match(paste(tile$condition,tile$n_eval),paste(archived$stage,archived$n_eval))
stopifnot(nrow(tile)==9L,!anyNA(idx),!anyDuplicated(idx),all(tile$n==100L),
 identical(as.integer(tile$count),as.integer(archived$shape_alert_count[idx])),
 identical(as.integer(tile$n),as.integer(archived$shape_alert_n[idx])),
 identical(sort(unique(tile$n_eval)),c(1000L,5000L,10025L)))
tile$count_label <- paste0(tile$count,'/',tile$n)
stopifnot(identical(tile$count_label,archived$count_label[idx]))
write.csv(tile,file.path(QA,'Figure3_B_tile_values.csv'),row.names=FALSE)
condition_order <- c('N1 continuous control','Empirical Lloyd refits','G2 discrete control')
tile$row <- factor(tile$condition,levels=rev(condition_order))
tile$column <- factor(tile$n_eval,levels=sort(unique(tile$n_eval)))
tile$positive <- factor(ifelse(tile$count==0,'None','At least one'),levels=c('None','At least one'))
threshold_label <- paste0('P < ',unique(tile$shape_alert_threshold))
tile$median_label <- sprintf('P = %.3f',tile$median_p)
tile$median_label[tile$median_p==0] <- 'P = 0'
shape_plot <- ggplot(tile,aes(column,row)) +
 geom_tile(aes(fill=positive),width=.94,height=.87,colour=GREY,linewidth=LW) +
 geom_text(aes(label=count_label),nudge_y=.12,fontface='bold') +
 geom_text(aes(label=median_label),nudge_y=-.17) +
 scale_fill_manual(values=c(None=SPEC$neutral_fill,`At least one`=SPEC$signal_fill),name=threshold_label) +
 scale_x_discrete(labels=fmt_n(sort(unique(tile$n_eval))),expand=expansion(add=.5),position='top') +
 scale_y_discrete(labels=c('Discrete control','Empirical K2','Continuous Gaussian-\ncopula reference'),expand=expansion(add=.6)) +
 labs(x='Evaluation sample size',y=NULL) +
 theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='bottom',
       legend.justification='right',legend.margin=margin(0,0,0,0),
       plot.margin=margin(2,2,0,0),axis.title.x=element_text(margin=margin(b=5)),
       axis.text.y=element_text(hjust=1),legend.key.width=unit(3,'mm'))

orders <- c('Base + K2','Locked K2','Imputation-specific K2','Fold-wise K2')
stopifnot(!anyNA(match(d2$analysis,orders)),
 all(d2$label_object[d2$analysis=='Locked K2']=='fixed_full_cohort_MICE_imputation1_plain_kmeans_K2'),
 all(d2$label_object[d2$analysis %in% c('Base + K2','Imputation-specific K2')]=='MICE_imputation_specific_K2_with_Rubin_pooling'),
 all(d2$label_object[d2$analysis=='Fold-wise K2']=='outcome_specific_foldwise_cross_fitted_K2'),
 all(is.na(d2$low)==is.na(d2$high)),
 all(is.na(d2$low)==(d2$analysis %in% c('Base + K2','Locked K2'))))
for(outcome in c('Mortality','Renal')) {
  suffix <- if(outcome=='Mortality') 'mortality_30d' else 'make30'
  mi <- bridge[bridge$evidence_id==paste0('D2_mice_',suffix),]
  cf <- bridge[bridge$evidence_id==paste0('D2_crossfit_',suffix),]
  getrow <- function(a) d2[d2$outcome==outcome & d2$analysis==a,]
  stopifnot(getrow('Base + K2')$estimate==mi$baseline_plus_label_increment,
    getrow('Locked K2')$estimate==if(outcome=='Mortality') fixed$mort_dAUC else fixed$make_dAUC,
    all(unlist(getrow('Imputation-specific K2')[c('estimate','low','high')])==unlist(mi[c('estimate','ci_low','ci_high')])),
    all(unlist(getrow('Fold-wise K2')[c('estimate','low','high')])==unlist(cf[c('estimate','ci_low','ci_high')])))
}
ok <- !is.na(d2$low)
stopifnot(all(d2$low[ok]<=d2$estimate[ok] & d2$estimate[ok]<=d2$high[ok]),
 all(oracle$mean_ci_low<=oracle$mean_delta_auc & oracle$mean_delta_auc<=oracle$mean_ci_high),
 all(abs(oracle$mean_ci_low-(oracle$mean_delta_auc-qt(.975,oracle$n_repeats-1)*oracle$mcse_mean))<1e-12),
 all(abs(oracle$mean_ci_high-(oracle$mean_delta_auc+qt(.975,oracle$n_repeats-1)*oracle$mcse_mean))<1e-12))
oracle$outcome_display <- ifelse(oracle$outcome=='mortality_30d','30-day mortality','Renal composite')
oracle$outcome_display <- factor(oracle$outcome_display,levels=c('30-day mortality','Renal composite'))
colours <- c(`30-day mortality`=BLUE,`Renal composite`=ORANGE)
symbols <- c(`30-day mortality`=21,`Renal composite`=22)
oracle_plot <- ggplot(oracle,aes(oracle_accuracy,mean_delta_auc,colour=outcome_display,shape=outcome_display,
                                linetype=outcome_display)) +
 geom_hline(yintercept=0,colour=GREY,linewidth=LW,linetype=2) +
 geom_line(linewidth=LW) +
 geom_errorbar(aes(ymin=mean_ci_low,ymax=mean_ci_high),width=.012,linewidth=LW) +
 geom_point(fill='white',size=1.8,stroke=LW*1.5) +
 scale_colour_manual(values=colours) + scale_shape_manual(values=symbols) +
 scale_linetype_manual(values=c(1,2)) +
 scale_x_continuous(breaks=sort(unique(oracle$oracle_accuracy)),limits=c(.48,1.02)) +
 scale_y_continuous(breaks=c(0,.02,.04,.06),limits=c(-.0045,.067)) +
 labs(x='Fixed noisy-oracle accuracy',y='Mean delta AUC (95% MC CI)') +
 theme(legend.position='none',plot.margin=margin(5,7,7,3))

text3 <- function(label,x,y,just='left',bold=FALSE,size=FS,col=INK) {
 grid.text(label,x=x,y=y,just=just,gp=gpar(fontfamily=FONT,fontsize=size,col=col,
                                        fontface=if(bold) 'bold' else 'plain'))
}
line3 <- function(x0,x1,y0,y1,col=INK,lty=1) grid.segments(x0=x0,x1=x1,y0=y0,y1=y1,
 gp=gpar(col=col,lwd=LW*72/25.4,lty=lty))
point3 <- function(x,y,col=BLUE,pch=21,fill='white') grid.points(x,y,pch=pch,size=unit(1.8,'mm'),
 gp=gpar(col=col,fill=fill,lwd=LW*72/25.4))
axis3 <- function(lim,ticks,labels,left,right,bottom,top,title,reference=numeric()) {
 pos <- function(v) left+(v-lim[1])/diff(lim)*(right-left)
 for(v in reference) line3(pos(v),pos(v),bottom,top,col=GREY,lty=2)
 line3(left,right,bottom,bottom)
 for(i in seq_along(ticks)) {
   line3(pos(ticks[i]),pos(ticks[i]),bottom,bottom-.015)
   text3(labels[i],pos(ticks[i]),bottom-.05,just='centre')
 }
 text3(title,(left+right)/2,bottom-.12,just='centre')
 pos
}
draw_a <- function() {
 ref <- sep[sep$evidence_type=='Reference',]
 empirical <- sep[sep$evidence_type=='Separation alert',]
 text3('Fixed K2 (full cohort)',.05,.85)
 text3('95th-percentile copula reference',.05,.54)
 pos <- axis3(c(2.62,3.82),c(2.7,3.0,3.3,3.6),c('2.7','3.0','3.3','3.6'),.08,.89,.23,.78,
               'Pooled-within Mahalanobis separation')
 line3(pos(ref$estimate),pos(ref$estimate),.23,.47,col=GREY,lty=2)
 line3(pos(ref$estimate),pos(ref$estimate),.66,.78,col=GREY,lty=2)
 point3(pos(empirical$estimate),.73,BLUE)
 text3(sprintf('%.3f',empirical$estimate),pos(empirical$estimate)+.04,.73)
 line3(pos(ref$ci_low),pos(ref$ci_high),.42,.42,INK)
 line3(pos(ref$ci_low),pos(ref$ci_low),.40,.44,INK)
 line3(pos(ref$ci_high),pos(ref$ci_high),.40,.44,INK)
 point3(pos(ref$estimate),.42,INK,23)
 text3(sprintf('%.3f',ref$estimate),pos(ref$estimate)+.04,.42)
}
draw_c <- function() {
 ys <- c(.80,.60,.43,.26)
  labels <- c('Base + K2 vs Base*\n(imputation-specific)',
             'Fixed K2 vs Raw-EN*\n(full cohort)',
             'Imputation-specific K2\nvs Raw-EN',
             'Fold-wise K2\nvs Raw-EN')
 text3('Limited clinical comparator',.05,.91,bold=TRUE)
 line3(.05,.97,.72,.72,col=GREY)
 text3('All source variables (Raw-EN)',.05,.68,bold=TRUE)
  text3(labels,.50,ys,just='right')
 pos <- axis3(c(-.006,.057),c(0,.02,.04),c('0','0.02','0.04'),.56,.97,.16,.9,'Delta AUC',0)
 for(i in seq_len(nrow(d2))) {
   y <- ys[match(d2$analysis[i],orders)]+if(d2$outcome[i]=='Mortality') .027 else -.027
   col <- if(d2$outcome[i]=='Mortality') BLUE else ORANGE
   if(!is.na(d2$low[i])) {
     line3(pos(d2$low[i]),pos(d2$high[i]),y,y,col)
     line3(pos(d2$low[i]),pos(d2$low[i]),y-.01,y+.01,col)
     line3(pos(d2$high[i]),pos(d2$high[i]),y-.01,y+.01,col)
     text3(sprintf('[%.4f, %.4f]',d2$low[i],d2$high[i]),.97,y,
           just='right',size=FS,col=col)
   }
   point3(pos(d2$estimate[i]),y,col,if(d2$outcome[i]=='Mortality') 21 else 22,
          if(is.na(d2$low[i])) 'white' else col)
 }
}
draw_figure3 <- function() {
 grid.newpage()
 pushViewport(viewport(width=unit(170,'mm'),height=unit(157,'mm'),
                       layout=grid.layout(2,1,heights=unit(c(65,92),'mm'))))
 pushViewport(viewport(layout.pos.row=1,layout=grid.layout(1,2,widths=unit(c(69,101),'mm'))))
 pushViewport(viewport(layout.pos.col=1))
 text3('A',.005,.975,bold=TRUE,size=TAG)
 pushViewport(viewport(y=.46,height=.87));draw_a();popViewport(2)
 pushViewport(viewport(layout.pos.col=2))
 text3('B',.005,.975,bold=TRUE,size=TAG)
 print(shape_plot,newpage=FALSE,vp=viewport(y=.44,height=.86,width=.97))
 popViewport(2)
 pushViewport(viewport(layout.pos.row=2,layout=grid.layout(1,2,widths=unit(c(97,73),'mm'))))
 pushViewport(viewport(layout.pos.col=1))
  text3('C',.005,.97,bold=TRUE,size=TAG)
  point3(.26,.91,BLUE,21);text3('30-day mortality',.285,.91)
  point3(.65,.91,ORANGE,22);text3('Renal composite',.675,.91)
  text3('Open symbols: point estimates without available intervals',.04,.875,size=FS)
  pushViewport(viewport(y=.425,height=.84));draw_c();popViewport(2)
 pushViewport(viewport(layout.pos.col=2))
 text3('D',.005,.97,bold=TRUE,size=TAG)
 print(oracle_plot,newpage=FALSE,vp=viewport(y=.445,height=.88,width=1))
 popViewport(3)
}
export_figure(draw_figure3,'Figure3',164)
svglite::svglite(file.path(OUT,'Figure3.svg'),width=SPEC$width_mm/25.4,height=164/25.4)
draw_figure3();dev.off()
if (requireNamespace('pdftools',quietly=TRUE) &&
    requireNamespace('xml2',quietly=TRUE) &&
    requireNamespace('magick',quietly=TRUE)) {
pdf_text <- pdftools::pdf_text(file.path(OUT,'Figure3.pdf'))
forbidden <- 'audit|module|alert|inadequate|inconclusive|surrogate|locked|MAKE[- ]?30|kidney[[:space:]]*composite'
stopifnot(!any(grepl(forbidden,pdf_text,ignore.case=TRUE)),any(grepl('Fixed K2 \\(full cohort\\)',pdf_text)))
writeLines(pdf_text,file.path(QA,'Figure3_text.txt'))
pass('Figure3',c('Separation full-cohort 3.521073; continuous reference 2.716950.',
 'Nine condition/size summary cells; each represents 100 repeats and matches the archived alert counts.',
 'All nine BH P < 0.05 numerator/denominator pairs exactly match the archived tile CSV.',
 'Panel C: all eight estimates verified against bridge/fixed-label CSVs; all four available intervals match.',
 'Panel C: four absent intervals remain absent; asterisks identify the two affected comparisons.',
 'Three label objects retained: fixed full cohort, imputation-specific, outcome-specific fold-wise.',
 'Panel D: all 12 points and interval endpoints read unchanged; MC CI matches mean +/- t(0.975, n-1)*MCSE.',
 'No prohibited figure text in extracted PDF. Width 174 mm; Arial 8 pt; tags 10 pt; lines >=0.20 mm.'))

svg <- xml2::read_xml(file.path(OUT,'Figure3.svg'))
styles <- xml2::xml_attr(xml2::xml_find_all(svg,'//*[@style]'),'style')
font_sizes <- as.numeric(sub('.*font-size: ([0-9.]+)px.*','\\1',styles[grepl('font-size:',styles)]))
line_styles <- styles[grepl('stroke-width:',styles) & !grepl('stroke: none',styles)]
line_pt <- as.numeric(sub('.*stroke-width: ([0-9.]+).*','\\1',line_styles))
stopifnot(all(font_sizes %in% c(8,10)),min(line_pt)*25.4/72>=.2)
pdf_boxes <- pdftools::pdf_data(file.path(OUT,'Figure3.pdf'))[[1]]
stopifnot(all(pdf_boxes$x>=0 & pdf_boxes$y>=0),
 all(pdf_boxes$x+pdf_boxes$width<=SPEC$width_mm/25.4*72+1),
 all(pdf_boxes$y+pdf_boxes$height<=164/25.4*72+1))
writeLines(c(paste('SVG font sizes (pt):',paste(sort(unique(font_sizes)),collapse=', ')),
 sprintf('Minimum visible SVG stroke: %.3f mm',min(line_pt)*25.4/72),
 'PDF text boxes are within page bounds.'),file.path(QA,'Figure3_export_checks.txt'))
magick::image_write(magick::image_convert(magick::image_read(file.path(QA,'Figure3.png')),
                                         colorspace='gray'),file.path(QA,'Figure3_grayscale.png'))
}
