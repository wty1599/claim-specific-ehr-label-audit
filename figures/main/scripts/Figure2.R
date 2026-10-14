script_file <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1]), winslash='/')
source(file.path(dirname(script_file), 'common.R'))
if (.Platform$OS.type=='windows') invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
# R device line units render at roughly 0.75 of the nominal ggplot width.
LW <- max(LW,.28)
src <- function(f) figure_reader(file.path(DATA,'fig23',f),check.names=FALSE,
                            stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM')
loadings <- src('Figure3_panelA_top_loadings_v2.csv')
recon <- src('Figure3_panelB_identity_auc_v2.csv')
ablation <- src('Figure3_panelC_ablation_v2.csv')
kdigo <- src('KDIGO_criterion4_locked_K2.csv')
loadings <- loadings[order(loadings$loading_rank), ]
recon <- recon[match(c('sofa_nonrenal','renal_core','acid_base','renal_acid_base'), recon$model_id), ]
ablation <- ablation[match(c('A3','A4','A5'), ablation$variant_id), ]
kdigo <- kdigo[kdigo$cluster == 'C1 higher-risk', ]
kdigo <- kdigo[match(c('Neither creatinine nor urine output','Urine output only',
                       'Creatinine only','Creatinine + urine output'), kdigo$criterion_pattern), ]
stopifnot(SPEC$width_mm == 174, FS == 8, TAG == 10, LW >= .2,
          nrow(loadings) == 10L, identical(loadings$loading_rank, 1:10),
          identical(loadings$standardized_variable_name, c('bun_max','creatinine_max','aniongap_max',
            'bicarbonate_min','lactate_max','glucose_max','ph_min','potassium_max','sodium_min','inr_max')),
          isTRUE(all.equal(loadings$abs_w, c(.420942121600222,.420512198689902,.327791357704267,
            .312351019582630,.306883558401731,.260284706186111,.197950002037362,
            .191384285609545,.174758045317930,.160246373725941), tolerance=1e-12)),
          all(abs(loadings$signed_w) == loadings$abs_w),
          identical(round(recon$auc,3), c(.589,.886,.944,.972)),
          identical(round(ablation$mean,3), c(.572,.284,.927)),
          identical(round(100*kdigo$proportion_cluster_given_pattern,1), c(11,13,41.8,59.3)),
          sum(kdigo$n) == 3992L, sum(kdigo$pattern_total) == 20049L,
          all(recon$auc_ci_low <= recon$auc & recon$auc <= recon$auc_ci_high),
          all(ablation$low <= ablation$mean & ablation$mean <= ablation$high),
          all(kdigo$cluster_given_pattern_wilson_low <= kdigo$proportion_cluster_given_pattern),
          all(kdigo$proportion_cluster_given_pattern <= kdigo$cluster_given_pattern_wilson_high),
          all(abs(kdigo$n/kdigo$pattern_total-kdigo$proportion_cluster_given_pattern) < 1e-12))

# Fixed column geometry keeps labels and numerical columns readable at print size.
text2 <- function(label, x, y, just='left', bold=FALSE, size=FS, col=INK) {
  grid.text(label, x=x, y=y, just=just,
            gp=gpar(fontfamily=FONT, fontsize=size, fontface=if(bold) 'bold' else 'plain', col=col))
}
line2 <- function(x0,x1,y0,y1,col=INK,lty=1) {
  grid.segments(x0=x0,x1=x1,y0=y0,y1=y1,gp=gpar(col=col,lwd=LW*72/25.4,lty=lty))
}
point2 <- function(x,y,col,pch=21,fill=col) {
  grid.points(x,y,pch=pch,size=unit(1.8,'mm'),gp=gpar(col=col,fill=fill,lwd=LW*72/25.4))
}
axis2 <- function(lim,ticks,labels,left,right,bottom,top,title,reference=numeric()) {
  pos <- function(v) left+(v-lim[1])/diff(lim)*(right-left)
  for(v in reference) line2(pos(v),pos(v),bottom,top,col=GREY,lty=2)
  line2(left,right,bottom,bottom)
  for(i in seq_along(ticks)) {
    line2(pos(ticks[i]),pos(ticks[i]),bottom,bottom-.012)
    text2(labels[i],pos(ticks[i]),bottom-.035,just='centre')
  }
  text2(title,(left+right)/2,bottom-.09,just='centre')
  pos
}
interval2 <- function(est,lo,hi,y,pos,col,pch=21,fill=col) {
  line2(pos(lo),pos(hi),y,y,col)
  line2(pos(lo),pos(lo),y-.012,y+.012,col)
  line2(pos(hi),pos(hi),y-.012,y+.012,col)
  point2(pos(est),y,col,pch,fill)
}
domain_col <- c(Renal=BLUE,`Acid-base`=ORANGE,Other=GREY)
domain_pch <- c(Renal=22,`Acid-base`=24,Other=21)
draw_a <- function() {
  ys <- seq(.83,.24,length.out=nrow(loadings))
  text2('|w|',.97,.90,just='right',bold=TRUE)
  pos <- axis2(c(0,.45),c(0,.2,.4),c('0','0.2','0.4'),.58,.88,.17,.85,'Absolute loading')
  text2(loadings$display_label,.52,ys,just='right')
  for(i in seq_len(nrow(loadings))) {
    col <- domain_col[[loadings$domain_group[i]]]
    grid.rect(x=pos(0),y=ys[i],width=pos(loadings$abs_w[i])-pos(0),height=.034,
              just='left',gp=gpar(col=NA,fill=col))
  }
  text2(sprintf('%.3f',loadings$abs_w),.97,ys,just='right')
  for(i in seq_along(domain_col)) {
    x <- c(.17,.43,.76)[i]
    grid.rect(x=x,y=.01,width=.018,height=.018,
              gp=gpar(col=NA,fill=domain_col[i]))
    text2(names(domain_col)[i],x+.025,.01)
  }
}
draw_b <- function() {
  ys <- seq(.76,.28,length.out=nrow(recon))
  text2('AUC [95% CI]',.98,.90,just='right',bold=TRUE)
  text2(c('Non-renal\nSOFA','Renal core','Acid-base','Renal +\nacid-base'),.30,ys,just='right')
  pos <- axis2(c(.5,1),c(.5,.75,1),c('0.50','0.75','1.00'),.35,.60,.17,.85,'Out-of-fold AUC',.5)
  cols <- c(GREY,BLUE,ORANGE,INK); shapes <- c(21,22,24,23)
  for(i in seq_len(nrow(recon))) interval2(recon$auc[i],recon$auc_ci_low[i],recon$auc_ci_high[i],
    ys[i],pos,cols[i],shapes[i],if(i==4) 'white' else cols[i])
  text2(sprintf('%.3f [%.3f-%.3f]',recon$auc,recon$auc_ci_low,recon$auc_ci_high),.98,ys,just='right')
}
draw_c <- function() {
  ys <- seq(.75,.31,length.out=nrow(ablation))
  text2('ARI',.96,.90,just='right',bold=TRUE)
  text2(c('Remove acid-base','Remove renal\n+ acid-base','Remove non-renal\nSOFA'),.45,ys,just='right')
  pos <- axis2(c(0,1),c(0,.5,1),c('0','0.5','1'),.50,.81,.17,.85,
               'ARI with fixed labels',1)
  cols <- c(ORANGE,INK,GREY); shapes <- c(24,23,21)
  for(i in seq_len(nrow(ablation))) interval2(ablation$mean[i],ablation$low[i],ablation$high[i],
    ys[i],pos,cols[i],shapes[i],if(i==2) 'white' else cols[i])
  text2(sprintf('%.3f',ablation$mean),.96,ys,just='right')
}
draw_d <- function() {
  ys <- seq(.76,.28,length.out=nrow(kdigo))
  text2('C1 % (n/N)',.99,.90,just='right',bold=TRUE)
  text2(c('Neither criterion','Urine output\nonly','Creatinine only','Both criteria'),.30,ys,just='right')
  pos <- axis2(c(0,.70),c(0,.3,.6),c('0','30','60'),.35,.62,.17,.85,'C1 prevalence (%)')
  for(i in seq_len(nrow(kdigo))) interval2(kdigo$proportion_cluster_given_pattern[i],
    kdigo$cluster_given_pattern_wilson_low[i],kdigo$cluster_given_pattern_wilson_high[i],
    ys[i],pos,BLUE)
  text2(sprintf('%.1f (%s/%s)',100*kdigo$proportion_cluster_given_pattern,fmt_n(kdigo$n),
                fmt_n(kdigo$pattern_total)),.99,ys,just='right')
}
draw_figure2 <- function() {
  grid.newpage()
  pushViewport(viewport(width=unit(170,'mm'),height=unit(151,'mm'),
                        layout=grid.layout(2,2,widths=unit(c(1,1),'null'),heights=unit(c(1.08,.92),'null'))))
  for(i in 1:4) {
    pushViewport(viewport(layout.pos.row=(i-1)%/%2+1,layout.pos.col=(i-1)%%2+1))
    text2(LETTERS[i],.005,.975,bold=TRUE,size=TAG)
    pushViewport(viewport(x=.5,y=.455,width=.98,height=.88))
    list(draw_a,draw_b,draw_c,draw_d)[[i]]()
    popViewport(2)
  }
  popViewport()
}
export_figure(draw_figure2,'Figure2',158)
svglite::svglite(file.path(OUT,'Figure2.svg'),width=SPEC$width_mm/25.4,height=158/25.4)
draw_figure2(); dev.off()
if (requireNamespace('pdftools',quietly=TRUE) &&
    requireNamespace('xml2',quietly=TRUE) &&
    requireNamespace('magick',quietly=TRUE)) {
pdf_text <- pdftools::pdf_text(file.path(OUT,'Figure2.pdf'))
forbidden <- 'audit|module|alert|inadequate|inconclusive|surrogate|locked|MAKE[- ]?30|kidney[[:space:]]*composite'
stopifnot(!any(grepl(forbidden,pdf_text,ignore.case=TRUE)))
writeLines(pdf_text,file.path(QA,'Figure2_text.txt'))
pass('Figure2',c('Top ten loading identities, ranks, and full-precision values match.',
 'AUC: 0.589, 0.886, 0.944, 0.972; ARI: 0.572, 0.284, 0.927.',
 'KDIGO C1 prevalence: 11.0%, 13.0%, 41.8%, 59.3%; numerators and denominators match.',
 'All intervals are imported unchanged and contain their point estimates.',
 'No prohibited figure text in extracted PDF. Width 174 mm; Arial 8 pt; tags 10 pt; lines >=0.20 mm.'))

svg <- xml2::read_xml(file.path(OUT,'Figure2.svg'))
styles <- xml2::xml_attr(xml2::xml_find_all(svg,'//*[@style]'),'style')
font_sizes <- as.numeric(sub('.*font-size: ([0-9.]+)px.*','\\1',styles[grepl('font-size:',styles)]))
line_styles <- styles[grepl('stroke-width:',styles) & !grepl('stroke: none',styles)]
line_pt <- as.numeric(sub('.*stroke-width: ([0-9.]+).*','\\1',line_styles))
stopifnot(all(font_sizes %in% c(8,10)),min(line_pt)*25.4/72>=.2)
pdf_boxes <- pdftools::pdf_data(file.path(OUT,'Figure2.pdf'))[[1]]
stopifnot(all(pdf_boxes$x>=0 & pdf_boxes$y>=0),
 all(pdf_boxes$x+pdf_boxes$width<=SPEC$width_mm/25.4*72+1),
 all(pdf_boxes$y+pdf_boxes$height<=158/25.4*72+1))
writeLines(c(paste('SVG font sizes (pt):',paste(sort(unique(font_sizes)),collapse=', ')),
 sprintf('Minimum visible SVG stroke: %.3f mm',min(line_pt)*25.4/72),
 'PDF text boxes are within page bounds.'),file.path(QA,'Figure2_export_checks.txt'))
magick::image_write(magick::image_convert(magick::image_read(file.path(QA,'Figure2.png')),
                                         colorspace='gray'),file.path(QA,'Figure2_grayscale.png'))
}
