script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(script_path), "../../.."), winslash="/"))
source(file.path(repo,"R/common/02_figure_paths.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript Figure4.R [authorized_figure4_aggregate_dir]")
figure_paths <- ehr_audit_figure_paths(repo,
  figure4_source = if (length(args)) args[[1]] else NULL,
  require_sources = "figure4", create_outputs = TRUE)
src <- figure_paths$figure4_source
out <- figure_paths$figure4_output
qa <- figure_paths$figure4_qa
figure_reader <- ehr_audit_figure_reader(src, file.path(qa,"Figure4_input_fingerprints.csv"), repo)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(qa,recursive=TRUE,showWarnings=FALSE)
read_src <- function(name) figure_reader(file.path(src,name),stringsAsFactors=FALSE,check.names=FALSE)
fidelity <- read_src("L0_L1_L2_fidelity.csv")
hospital <- read_src("Figure4_panelB_public_funnel_points.csv")
summary <- read_src("Figure4_panelB_funnel_summary.csv")
prev <- read_src("Figure4_pooled_L1_prevalence.csv")
curve <- read_src("heldout_calibration_curves.csv")
risk <- read_src("heldout_risk_distribution_screened.csv")
point <- read_src("heldout_performance_point.csv")
split <- read_src("hospital_holdout_split_summary.csv")
performance <- read_src("Figure4_panelC_L1_transport.csv")
calibration <- read_src("Figure4_panelC_L1_calibration.csv")
ranges <- read_src("Figure4_reference_ranges.csv")

stopifnot(identical(fidelity$layer,c("L0","L1","L2")),
          all(fidelity$n==20049L),
          identical(round(fidelity$concordance,3),c(1,.989,.986)),
          identical(round(fidelity$ari,3),c(1,.949,.936)),
          sum(hospital$hospitals_aggregated)==199L,
          sum(hospital$n_stays)==17465L,
          sum(hospital$point_type=="hospital")==96L,
          sum(hospital$point_type=="pooled")==5L,
          all(hospital$hospitals_aggregated[hospital$point_type=="pooled"]>=2L),
          all(round(hospital$n_stays*hospital$c1_prevalence)[hospital$point_type=="hospital"]>=10L),
          all((hospital$n_stays-round(hospital$n_stays*hospital$c1_prevalence))[hospital$point_type=="hospital"]>=10L),
          summary$outside_99_8_mimic==9L,
          nrow(curve)==1206L,
          nrow(risk)==56L, all(risk$stays>=10L),
          nrow(point)==6L,
          all(point$evaluation_stays==6478L),
          all(point$evaluation_events==1068L),
          split$stay_n[split$side=="evaluation"]==6478L,
          split$hospital_n[split$side=="evaluation"]==60L,
          nrow(performance)==8L, nrow(calibration)==4L,
          all(performance$n_external_stays[performance$outcome=="mortality"]==17284L),
          all(performance$n_external_stays[performance$outcome=="make"]==17333L),
          all(calibration$bootstrap_valid==1000L))
for (key in unique(paste(risk$model,risk$method))) {
  rows <- paste(risk$model,risk$method)==key
  stopifnot(sum(risk$stays[rows])==6478L)
}
mimic <- prev$c1_prevalence[prev$cohort=="MIMIC-IV"]
eicu <- prev$c1_prevalence[prev$cohort=="eICU"]
stopifnot(length(mimic)==1L,length(eicu)==1L,
          abs(mimic-3902/20049)<1e-12,
          abs(eicu-3319/17465)<1e-12)

width_mm <- 174
height_mm <- 276
font <- "Arial"
fs <- 8
tag <- 10
ink <- "#424242"
blue <- "#0072B2"
orange <- "#C15B29"
grey <- "#9B9B9B"
line_pt <- .22*96/25.4
mm <- function(x) grid::unit(x,"mm")
shown <- character()
txt <- function(label,x,y,bold=FALSE,just="left",size=fs,col=ink,rot=0) {
  shown <<- c(shown,as.character(label))
  grid::grid.text(label,x=mm(x),y=mm(height_mm-y),just=just,rot=rot,
                  gp=grid::gpar(fontfamily=font,fontsize=size,
                                fontface=if(bold) "bold" else "plain",col=col))
}
seg <- function(x1,y1,x2,y2,col=ink,lty=1) {
  grid::grid.segments(mm(x1),mm(height_mm-y1),mm(x2),mm(height_mm-y2),
                      gp=grid::gpar(col=col,lwd=line_pt,lty=lty,lineend="butt"))
}
path <- function(x,y,col=ink,lty=1) {
  grid::grid.lines(mm(x),mm(height_mm-y),
                   gp=grid::gpar(col=col,lwd=line_pt*1.45,lty=lty,lineend="round"))
}
mark <- function(x,y,pch=16,col=ink,size=1.5) {
  grid::grid.points(mm(x),mm(height_mm-y),pch=pch,size=mm(size),
                    gp=grid::gpar(col=col,fill=col,lwd=line_pt))
}
rect <- function(x1,y1,x2,y2,fill=grey,alpha=1) {
  grid::grid.rect(mm((x1+x2)/2),mm(height_mm-(y1+y2)/2),
                  width=mm(x2-x1),height=mm(y2-y1),
                  gp=grid::gpar(col=NA,fill=grDevices::adjustcolor(fill,alpha.f=alpha)))
}
axis_x <- function(ticks,xmap,y,labels=format(ticks,trim=TRUE)) {
  seg(xmap(min(ticks)),y,xmap(max(ticks)),y)
  for(i in seq_along(ticks)) {
    seg(xmap(ticks[i]),y,xmap(ticks[i]),y+1)
    txt(labels[i],xmap(ticks[i]),y+4,just="centre")
  }
}
draw_a <- function() {
  txt("A",2,4,bold=TRUE,size=tag)
  txt("Estimate",78,10,bold=TRUE,just="centre")
  ax <- function(v) 45+(v-.90)/.10*21
  rows <- c(20,27,34,47,54,61)
  vals <- c(fidelity$concordance,fidelity$ari)
  labels <- rep(paste(c("Original matrix","Deployment rule","Earlier reconstruction"),
                      paste0("(",fidelity$layer,")")),2)
  txt("Label agreement",3,13,bold=TRUE)
  txt("Adjusted Rand index",3,40,bold=TRUE)
  seg(ax(1),17,ax(1),63,col=grey,lty=2)
  for(i in seq_along(vals)) {
    txt(labels[i],3,rows[i])
    mark(ax(vals[i]),rows[i],pch=c(16,17,15)[(i-1)%%3+1],
         col=c(ink,blue,grey)[(i-1)%%3+1],size=1.7)
    txt(sprintf("%.3f",vals[i]),78,rows[i],just="centre")
  }
  axis_x(c(.90,.95,1),ax,66,labels=c("0.90","0.95","1.00"))
  txt("Agreement with fixed labels",44,77,just="centre")
}
draw_b <- function() {
  txt("B",89,4,bold=TRUE,size=tag)
  bx <- function(v) 103+log10(v)/log10(800)*66
  by <- function(v) 66-v*51
  n <- seq_len(746L)
  for(coverage in c(.95,.998)) {
    alpha <- 1-coverage
    for(lim in list(stats::qbinom(alpha/2,n,mimic)/n,
                    stats::qbinom(1-alpha/2,n,mimic)/n)) {
      path(bx(n),by(lim),col=if(coverage==.95) grey else ink,
           lty=if(coverage==.95) 2 else 3)
    }
  }
  seg(bx(1),by(mimic),bx(746),by(mimic),col=blue)
  seg(bx(1),by(eicu),bx(746),by(eicu),col=orange,lty=4)
  ind <- hospital[hospital$point_type=="hospital",]
  pooled <- hospital[hospital$point_type=="pooled",]
  for(i in seq_len(nrow(ind))) mark(bx(ind$plot_x[i]),by(ind$c1_prevalence[i]),size=.75)
  for(i in seq_len(nrow(pooled))) mark(bx(pooled$plot_x[i]),by(pooled$c1_prevalence[i]),
                                        pch=5,col=blue,size=2)
  seg(103,66,169,66);seg(103,15,103,66)
  for(v in c(0,.25,.5,.75,1)) {
    seg(102,by(v),103,by(v))
    txt(sprintf("%.0f%%",v*100),100,by(v),just="right")
  }
  for(v in c(1,10,100,700)) {
    seg(bx(v),66,bx(v),67)
    txt(format(v,scientific=FALSE),bx(v),71,just="centre")
  }
  seg(108,8,114,8,col=blue);txt(sprintf("MIMIC-IV %.2f%%",mimic*100),116,8,col=blue)
  seg(146,8,152,8,col=orange,lty=4);txt(sprintf("eICU %.2f%%",eicu*100),154,8,col=orange)
  txt("C1 prevalence",89,43,just="centre",rot=90)
  txt("Assigned stays per hospital (log scale)",136,76,just="centre")
  seg(101,80,107,80,col=grey,lty=2);txt("95%",109,80)
  seg(123,80,129,80,col=ink,lty=3);txt("99.8%",131,80)
  mark(147,80,pch=5,col=blue,size=1.7);txt("Pooled bands",151,80)
}
models <- c("base","base_gen","feat_pen","feat_rf")
model_labels <- c("Base","Base + K2","Raw-EN","Raw-RF")
metric_keys <- c("auc","slope","intercept")
metric_labels <- c("External AUC","Calibration slope","Calibration intercept")
transport <- do.call(rbind,lapply(c("mortality","make"),function(outcome) {
  do.call(rbind,lapply(metric_keys,function(metric) {
    z <- performance[performance$outcome==outcome,]
    z <- z[match(models,z$model),]
    value <- z[[if(metric=="auc") "external_auc" else paste0("calibration_",metric)]]
    low <- high <- rep(NA_real_,nrow(z))
    if(metric=="auc") {
      low <- z$external_auc_ci_low; high <- z$external_auc_ci_high
    } else {
      ci <- calibration[calibration$outcome==outcome & calibration$parameter==metric,]
      idx <- match(z$model,ci$model)
      low <- ci$ci_low[idx]; high <- ci$ci_high[idx]
      stopifnot(max(abs(value[!is.na(idx)]-ci$estimate[idx[!is.na(idx)]]))<1e-12)
    }
    gate <- ranges[ranges$quantity==paste0("calibration_",metric),]
    outside <- if(nrow(gate)) value<gate$lower | value>gate$upper else rep(FALSE,nrow(z))
    data.frame(outcome=outcome,model=z$model,metric=metric,value=value,
               low=low,high=high,outside=outside)
  }))
}))
stopifnot(nrow(transport)==24L,sum(!is.na(transport$low))==12L,
          abs(transport$value[transport$outcome=="mortality" & transport$model=="base_gen" & transport$metric=="slope"]-.779840685646913)<1e-12)
draw_c <- function() {
  seg(2,83,172,83,col=grey)
  txt("C",2,88,bold=TRUE,size=tag)
  mark(15,88,pch=17,col=orange,size=1.4)
  txt("Outside prespecified calibration range",19,88,size=7)
  starts <- c(37,83,129); widths <- c(37,37,40)
  limits <- list(c(.72,.83),c(.64,2.08),c(-.62,.22))
  ticks <- list(c(.74,.78,.82),c(.8,1,1.5,2),c(-.6,-.2,0,.2))
  ticklabels <- list(c("0.74","0.78","0.82"),c("0.8","1.0","1.5","2.0"),
                     c("-0.6","-0.2","0","0.2"))
  for(j in seq_along(metric_keys))
    txt(metric_labels[j],starts[j]+widths[j]/2,94,bold=TRUE,just="centre")
  for(i in seq_along(c("mortality","make"))) {
    outcome <- c("mortality","make")[i]
    top <- c(102,136)[i]
    yy <- top+c(0,5,10,15)
    txt(if(outcome=="mortality") "In-hospital death: 17,284 stays" else
          "eICU composite: 17,333 stays",3,top-5,bold=TRUE,size=7)
    for(k in seq_along(models)) txt(model_labels[k],3,yy[k],size=7)
    for(j in seq_along(metric_keys)) {
      z <- transport[transport$outcome==outcome & transport$metric==metric_keys[j],]
      z <- z[match(models,z$model),]
      lim <- limits[[j]]
      xx <- function(v) starts[j]+(v-lim[1])/diff(lim)*widths[j]
      stopifnot(all(z$value>=lim[1] & z$value<=lim[2]),
                all(z$low[!is.na(z$low)]>=lim[1]),
                all(z$high[!is.na(z$high)]<=lim[2]))
      gate <- ranges[ranges$quantity==paste0("calibration_",metric_keys[j]),]
      if(nrow(gate)) {
        rect(xx(gate$lower),top-2,xx(gate$upper),top+17,fill=grey,alpha=.13)
        seg(xx(gate$reference),top-2,xx(gate$reference),top+17,col=grey,lty=2)
      }
      for(k in seq_len(nrow(z))) {
        col <- if(z$outside[k]) orange else ink
        if(!is.na(z$low[k])) {
          seg(xx(z$low[k]),yy[k],xx(z$high[k]),yy[k],col=col)
          seg(xx(z$low[k]),yy[k]-.7,xx(z$low[k]),yy[k]+.7,col=col)
          seg(xx(z$high[k]),yy[k]-.7,xx(z$high[k]),yy[k]+.7,col=col)
        }
        mark(xx(z$value[k]),yy[k],pch=if(z$outside[k]) 17 else 16,col=col,size=1.35)
      }
      axis_x(ticks[[j]],xx,top+18,labels=ticklabels[[j]])
    }
  }
}
draw_hist <- function(model,method,xmap,top,col) {
  h <- risk[risk$model==model & risk$method==method,]
  stopifnot(nrow(h)==14L,sum(h$stays)==6478L,all(h$stays>=10L))
  density <- h$stays/(h$bin_high-h$bin_low)
  scale <- max(density)
  for(i in seq_len(nrow(h))) {
    left <- xmap(h$bin_low[i])+.08
    right <- xmap(h$bin_high[i])-.08
    bar_h <- 4.3*density[i]/scale
    rect(left,top-bar_h,right,top,fill=col,alpha=.9)
  }
  seg(xmap(0),top,xmap(1),top,col=grey)
}
draw_d_one <- function(model,x0,title) {
  xmap <- function(v) x0+45*v
  ymap <- function(v) 247-45*v
  txt(title,x0+22.5,192,bold=TRUE,just="centre")
  txt("n = 6,478 stays; 1,068 deaths",x0+22.5,197,just="centre")
  seg(xmap(0),ymap(0),xmap(1),ymap(0));seg(xmap(0),ymap(0),xmap(0),ymap(1))
  path(xmap(c(0,1)),ymap(c(0,1)),col=grey,lty=3)
  for(v in seq(0,1,.25)) {
    seg(xmap(v),247,xmap(v),248)
    seg(xmap(0)-1,ymap(v),xmap(0),ymap(v))
    txt(sprintf("%.2f",v),xmap(v),251,just="centre")
    txt(sprintf("%.2f",v),xmap(0)-2,ymap(v),just="right")
  }
  for(method in c("no_update","intercept_slope")) {
    z <- curve[curve$model==model & curve$method==method,]
    z <- z[order(z$grid_lp),]
    stopifnot(nrow(z)==201L,
              all(z$predicted_risk>=0 & z$predicted_risk<=1),
              all(z$observed_smooth>=0 & z$observed_smooth<=1))
    path(xmap(z$predicted_risk),ymap(z$observed_smooth),
         col=if(method=="no_update") ink else blue,
         lty=if(method=="no_update") 1 else 2)
  }
  draw_hist(model,"no_update",xmap,256,grey)
  draw_hist(model,"intercept_slope",xmap,264,blue)
}
draw_d <- function() {
  seg(2,166,172,166,col=grey)
  txt("D",2,173,bold=TRUE,size=tag)
  txt("Held-out hospital mortality calibration",13,173,bold=TRUE)
  seg(39,182,48,182,col=ink);txt("No update",51,182)
  seg(89,182,98,182,col=blue,lty=2);txt("Intercept + slope update",101,182)
  seg(151,182,158,182,col=grey,lty=3);txt("Ideal",160,182)
  draw_d_one("base",21,"Base")
  draw_d_one("base_l1_k2",109,"Base + K2")
  txt("Observed death proportion",4,226,just="centre",rot=90)
  txt("Predicted in-hospital death risk",93,273,just="centre")
  txt("Risk distributions: gray = original; blue = updated; 0.65-1.00 pooled",93,269,
      just="centre",size=7)
}
draw <- function() { grid::grid.newpage();draw_a();draw_b();draw_c();draw_d() }

grDevices::cairo_pdf(file.path(out,"Figure4.pdf"),width=width_mm/25.4,
                     height=height_mm/25.4,family=font,onefile=TRUE)
draw();grDevices::dev.off()
grDevices::tiff(file.path(out,"Figure4.tiff"),width=width_mm,height=height_mm,
                units="mm",res=600,compression="lzw",type="cairo",bg="white")
draw();grDevices::dev.off()
grDevices::png(file.path(qa,"Figure4.png"),width=width_mm,height=height_mm,
               units="mm",res=180,type="cairo",bg="white")
draw();grDevices::dev.off()
stopifnot(all(!grepl("patient|hospitalid|pass|unlocked",shown,ignore.case=TRUE)),
          abs(point$brier[point$model=="base_l1_k2" & point$method=="no_update"]-.121935797684142)<1e-10,
          abs(point$brier[point$model=="base_l1_k2" & point$method=="intercept_slope"]-.118241987606364)<1e-10)
writeLines(c("PASS", "Panel A: source fidelity; panel B: disclosure-screened hospital funnel.",
             "Panel C: full-eICU mortality and composite performance, four models.",
             "Panel D: same 6,478 stays before and after mortality-risk updating.",
             "Risk histogram bins each represent at least 10 stays; top tail pooled.",
             "No patient-level data or hospital identifiers were read."),
           file.path(qa,"Figure4_assertions.txt"))
writeLines(capture.output(sessionInfo()),file.path(qa,"Figure4_sessionInfo.txt"))
cat("Figure4 complete: source and disclosure assertions passed.\n")
