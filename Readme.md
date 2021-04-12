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

Exemple:
![gif]("FeaturePlot_presentation.gif")

PS: If you want to try it with a more consistent dataset, you can save dataset from [SeuratData packages]() as Rds file, then load it in the shiny app

```R
if (!require('shiny') ) { install.packages('shiny') }
if (!require('Seurat') ) { install.packages('Seurat') }
if (!require('devtools') ) { install.packages('devtools') }
if (!require('SeuratData') ) { devtools::install_github('satijalab/seurat-data') }

InstallData("pbmc3k")
saveRDS(pbmc3k.final, "pbmc3k.Rds")

shiny::runGitHub("martin-jeremy/toolbox-featureplot", ref = "main")
```
