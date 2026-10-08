this <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),winslash='/'))
source(file.path(this,'common.R'))
d <- read_source('fig1/nodes.csv')
for(nm in names(d)) d[[nm]] <- gsub('\\n','\n',d[[nm]],fixed=TRUE)
stopifnot(identical(d$id[d$section=='row'],paste0('D',1:4)))
draw <- function() {
  grid.newpage()
  pushViewport(viewport(xscale=c(0,174),yscale=c(126,0)))
  txt <- function(t,x,y,bold=FALSE,just='left') grid.text(t,x=unit(x,'native'),y=unit(y,'native'),just=just,
    gp=gpar(fontfamily=FONT,fontsize=FS,fontface=if(bold)'bold' else 'plain',col=INK,lineheight=1.1))
  line <- function(x,y,arrow=FALSE) grid.lines(x=unit(x,'native'),y=unit(y,'native'),
    arrow=if(arrow)grid::arrow(length=unit(1.6,'mm'),type='closed')else NULL,
    gp=gpar(col=INK,lwd=0.7))
  txt('A',4,5,TRUE)
  xx <- c(4,64,124)
  for(i in 1:3) {
    grid.rect(x=unit(xx[i],'native'),y=unit(12,'native'),width=unit(46,'mm'),height=unit(19,'mm'),just=c('left','top'),gp=gpar(fill=if(i==2)SPEC$neutral_fill else 'white',col=INK,lwd=0.7))
    txt(d$question[i],xx[i]+23,21.5,just='centre')
  }
  line(c(50,64),c(21.5,21.5),TRUE);line(c(110,124),c(21.5,21.5),TRUE)
  txt('B',4,42,TRUE)
  txt('Question',10,49,TRUE);txt('Label used',64,49,TRUE);txt('Comparison or reference',119,49,TRUE)
  line(c(10,170),c(54,54))
  rows <- d[d$section=='row',]
  for(i in 1:4) {
    y <- 62+(i-1)*17.5
    txt(rows$question[i],10,y)
    txt(rows$label_used[i],64,y)
    txt(rows$comparison[i],119,y)
    if(i<4) grid.lines(x=unit(c(10,170),'native'),y=unit(rep(y+8.75,2),'native'),gp=gpar(col=GREY,lwd=0.7))
  }
  line(c(147,147,6,6),c(31,37,37,114.5))
  for(y in 62+(0:3)*17.5) line(c(6,9),c(y,y),TRUE)
  popViewport()
}
export_figure(draw,'Figure1',126)
pass('Figure1','Design text adapted from the approved Figure 1 script; no numerical analysis. Four questions; D4 compares RRT records; shared visual specification.')
