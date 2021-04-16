# Toolbox - FeaturePlot

Please, to run this shiny app, execute following lines in a R session:

```R
if (!require('shiny') ) {
  install.packages('shiny')
}

if (!require('shinythemes') ) {
  install.packages('shinythemes')
}

if (!require('Seurat') ) {
  install.packages('Seurat')
}

shiny::runGitHub("martin-jeremy/toolbox-featureplot", ref = "main")
```

Enjoy !

Exemple:

![gif](https://github.com/martin-jeremy/toolbox-featureplot/blob/main/FeaturePlot_presentation.gif)

PS: If you want to try it with a more consistent dataset, you can save dataset from [SeuratData packages (satijalab/seurat-data)](https://github.com/satijalab/seurat-data) as Rds file, then load it in the app :

```R
if (!require('shiny') ) { install.packages('shiny') }
if (!require('Seurat') ) { install.packages('Seurat') }
if (!require('devtools') ) { install.packages('devtools') }
if (!require('SeuratData') ) { devtools::install_github('satijalab/seurat-data') }

InstallData("pbmc3k")
saveRDS(pbmc3k.final, "pbmc3k.Rds")

shiny::runGitHub("martin-jeremy/toolbox-featureplot", ref = "main")
```
