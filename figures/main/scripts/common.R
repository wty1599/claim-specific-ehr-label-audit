argv <- commandArgs(FALSE)
sf <- sub('^--file=', '', grep('^--file=', argv, value=TRUE)[1])
ROOT <- normalizePath(file.path(dirname(sf), '..'), winslash='/', mustWork=TRUE)
repo <- Sys.getenv('EHR_AUDIT_REPO_ROOT', unset =
  normalizePath(file.path(ROOT, '../..'), winslash='/', mustWork=TRUE))
source(file.path(repo, 'R/common/02_figure_paths.R'))
figure_paths <- ehr_audit_figure_paths(repo, require_sources='main', create_outputs=TRUE)
DATA <- figure_paths$main_source
OUT <- figure_paths$main_output
QA <- figure_paths$main_qa
figure_reader <- ehr_audit_figure_reader(DATA,
  file.path(QA, paste0(tools::file_path_sans_ext(basename(sf)), '_input_fingerprints.csv')), repo)
suppressPackageStartupMessages({library(ggplot2); library(grid); library(patchwork)})
SPEC <- jsonlite::fromJSON(file.path(ROOT,'figure_spec.json'))
INK <- SPEC$ink; BLUE <- SPEC$blue; ORANGE <- SPEC$orange; GREY <- SPEC$grey
FONT <- SPEC$font; FS <- SPEC$body_pt; TAG <- SPEC$tag_pt; LW <- SPEC$line_mm
theme_set(theme_classic(base_size=FS, base_family=FONT) + theme(
  text=element_text(colour=INK), axis.text=element_text(size=FS,colour=INK),
  axis.title=element_text(size=FS), plot.title=element_text(size=FS,face='bold'),
  plot.tag=element_text(size=TAG,face='bold'),
  legend.text=element_text(size=FS),legend.title=element_text(size=FS),
  strip.text=element_text(size=FS,face='bold'),strip.background=element_blank(),
  axis.line=element_line(linewidth=LW), axis.ticks=element_line(linewidth=LW),
  panel.grid=element_blank(),plot.margin=margin(5,5,5,5),
  legend.key.height=unit(3,'mm'),legend.key.width=unit(5,'mm')))
update_geom_defaults('text',list(family=FONT,size=FS/ggplot2::.pt,colour=INK))
read_source <- function(name) figure_reader(file.path(DATA,name),check.names=FALSE,stringsAsFactors=FALSE)
fmt_n <- function(x) format(x,big.mark=',',scientific=FALSE,trim=TRUE)
fmt_f <- function(x,d=1) formatC(x,format='f',digits=d,big.mark=',')
export_figure <- function(plot,name,height_mm,width_mm=SPEC$width_mm) {
  draw <- function() { if(is.function(plot)) plot() else if(inherits(plot,'grob')) {grid.newpage();grid.draw(plot)} else print(plot) }
  grDevices::cairo_pdf(file.path(OUT,paste0(name,'.pdf')),width=width_mm/25.4,height=height_mm/25.4,family=FONT,onefile=TRUE)
  draw();dev.off()
  grDevices::tiff(file.path(OUT,paste0(name,'.tiff')),width=width_mm,height=height_mm,units='mm',res=600,compression='lzw',bg='white',type='cairo')
  draw();dev.off()
  grDevices::png(file.path(QA,paste0(name,'.png')),width=width_mm,height=height_mm,units='mm',res=300,bg='white',type='cairo')
  draw();dev.off()
  writeLines(c(paste('Width_mm',width_mm),paste('Height_mm',height_mm),paste('Body_pt',FS),paste('Tag_pt',TAG)),file.path(QA,paste0(name,'_dimensions.txt')))
}
pass <- function(name,details) writeLines(c('PASS',details),file.path(QA,paste0(name,'_assertions.txt')))
