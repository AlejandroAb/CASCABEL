library('plotly')
library('webshot')

#Read arguments from command line.
args <- commandArgs(trailingOnly = T)
#args[1] = path for seting WD to use current dir use $PWD
#args[2] = file with sequence distribution
#args[3] = path to storee the files
#arg[4] = chart type pie or bar
#play at https://plot.ly/create/?fid=RPlotBot:3157
setwd(args[1])
file <- args[2]
#histTable <- read.table(file)
tabla <- read.table(file, header = FALSE, sep = "\t")

if (args[4] == "pie"){
  m <- list(
  l = 10,
  r = 0,
  b = 150,
  t = 150,
  pad = 10

  )
  le <- list(orientation = "v",   # show entries horizontally
              yanchor = "middle",  # use center of legend as anchor
              xanchor = "right",
              y= 0,
              x=1

              )

  p <- plot_ly(tabla, labels = tabla[,2], values = tabla[,1], type ='pie') %>%
  layout(title = 'Sample distribution                   ',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend = T,
         margin = m,
         legend = le)
  export(p, file=args[3])
#print(p)
}else if (args[4] == "bar"){
    p <- plot_ly(tabla, y = tabla[,2], x = tabla[,1], type ='bar',orientation = 'h')
    export(p, file=args[3])
}else { #do both!
  m <- list(l = 10,  r = 0,  b = 150,  t = 150,  pad = 10  )
  le <- list(orientation = "v",   # show entries horizontally
              yanchor = "middle",  # use center of legend as anchor
              xanchor = "right",
              y= 0, x=1 )

  p <- plot_ly(tabla, labels = tabla[,2], values = tabla[,1], type ='pie') %>%
  layout(title = 'Sample distribution                   ',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend = T,
         margin = m,
         legend = le)
  export(p, file=args[3])

  p2 <- plot_ly(tabla, y = tabla[,2], x = tabla[,1], type ='bar',orientation = 'h') %>%
  layout(title = 'Sample distribution',
         bargap=0.1)
  bar_file=gsub("seqs_fw_rev_filtered.dist.png","seqs_fw_rev_filtered.dist.bar.png",args[3])
  export(p2, file=bar_file)

}
#export(p, file=paste(args[3],"seqs_fw_rev_filtered.dist.png", sep=''))
export(p, file=args[3])
