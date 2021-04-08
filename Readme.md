# Toolbox - FeaturePlot

Please, to run this shiny app, execute following lines in a R session:

```R
if (!require('shiny') ) {
  install.packages('shiny')
}

if (!require('Seurat') ) {
  install.packages('Seurat')
}

shiny::runGitHub("martin-jeremy/toolbox-featureplot", ref = "main")
```

Enjoy !
