# This is a little exercise, meant for my own educational purpose and that of others
# The goals are:
#   - download a dataset from the excellent bodenmiller et al. 2012
#   - load the data into a flowSet object
#   - extract a sample of the data
#   - run a SPADE analysis
#   - display the SPADE tree as a d3.js graph object
# script is intended to run from the root of the repository

# libraries
require(flowCore)
require(flowViz)
require(spade)
library(stringr)
library(igraph)
library(rjson)

# functions

# parameters
base.dir <- 'input'
out.dir <- 'output'

# the following are taken from the publication
stimulations <- c('Vanadate','IL-3','IL-2','IL-12','Reference','G-CSF','GM-CSF','BCR/FcR-XL','IFN-g','IFN-a','LPS','PMA/Ionomycin')
start.conc <- 2.5 # Staurosporine
dil <- 4

pheno.markers <- c('CD20','IgM','CD4','HLA-DR','CD14','CD7','CD3','CD123') # c('CD20','IgM','CD4','CD33','HLA-DR','CD14','CD7','CD3','CD123')
func.markers <- c('pStat1','pSlp76','pBtk','pPlcg2','pErk','pLat','pS6','pNFkB','pp38','pStat5','pAkt','pSHP2','pZap70','pStat3')

# check that some data is there
if(!'Staurosporine' %in% dir(base.dir)) {
  # if not download it from cytobank
  download.file('https://s3.amazonaws.com/reports.public.cytobank.org/105/Staurosporine.zip',
                base.dir)
  unzip(file.path(base.dir,'Staurosporine.zip'),exdir=base.dir)
} else {
  # check that we have the right number of files
  cat(length(dir(file.path(base.dir,'Staurosporine'))),'files found!')
}

# load all data into a flowSet
staurosporine <- read.flowSet(files=list.files(path=file.path(base.dir,'Staurosporine'),pattern = 'fcs'),
                              path = file.path(base.dir,'Staurosporine'))

# annotate the flowSet
annots <- lapply(sampleNames(staurosporine),function(x) {
    conds <- str_split(x,'\\.')[[1]][1]
    annot <- str_split(conds,'_')[[1]]
    annot <- data.frame(name=x,Compound=annot[1],Cells=annot[2],WellID=annot[3])
    return(annot)
  }
)
annots <- do.call('rbind',annots)


# reorder the cell lines
annots$Cells <- factor(annots$Cells,
                     levels=c('cd4+', 'cd8+', 'igm-', 'igm+', 
                              'cd14-hladr-', 'cd14-hladrmid', 'cd14-hladrhigh',
                              'cd14+hladr-', 'cd14+hladrmid', 'cd14+hladrhigh', 
                              'cd14-surf-', 'cd14+surf-', 'dendritic', 'nk')
)

well.annot <- data.frame(Row=rep(LETTERS[1:8],each=12),Column=rep(seq(1,12),8))
well.annot$WellID <- apply(well.annot,1,function(x) paste(x['Row'],sprintf('%02i',as.integer(x['Column'])),sep=''))
well.annot$Concentration <- start.conc/dil^(match(well.annot$Row,LETTERS)-1)
well.annot[well.annot$Row=='H','Concentration'] <- 0
well.annot$pConcentration <- -log10(well.annot$Concentration*1E6)
well.annot[well.annot$Row=='H','pConcentration'] <- NA
well.annot$Stimulation <- stimulations[well.annot$Column]

annots <- merge(annots,well.annot,by='WellID')

row.names(annots) <- annots$name
annots <- annots[sampleNames(staurosporine),]

pData(staurosporine) <- annots
varMetadata(phenoData(staurosporine))$labelDescription <- colnames(annots)

# check the parameters
pData(parameters(staurosporine[[1]]))

# generate a SPADE tree for untreated samples only
# remove samples where no cells would be left after downsampling
unstimulated.samples <- with(pData(staurosporine),Row=='H') & ( fsApply(staurosporine,nrow)*0.01>1 )

unstimulated.files <- sapply(subset(pData(staurosporine),unstimulated.samples,name),
                             function(x) file.path(base.dir,'Staurosporine',x),
                             simplify=FALSE)
unstimulated.files <- unlist(unstimulated.files)

# SPADE.driver(files = unstimulated.files,
#              out_dir = file.path(out.dir,'Staurosporine','SPADE','unstimulated'),
#              cluster_cols = subset(pData(parameters(staurosporine[[1]])),
#                                    desc %in% c(pheno.markers,func.markers),
#                                    name,drop=TRUE),
#              transforms = arcsinhTransform(transformationId="defaultArcsinhTransform", a=0, b=1/5, c=0)
# )

# load SPADE tree results
graph <- read.graph(file.path(out.dir,'Staurosporine','SPADE','unstimulated',"mst.gml"),format="gml")
layout <- read.delim(file.path(out.dir,'Staurosporine','SPADE','unstimulated',"layout.table"),header=F,sep=' ')

# check your results
plot(graph,layout=as.matrix(layout))

# now load the median values for all samples read.flowset( )
clusterUnstimSet <- read.flowSet(path = file.path(out.dir,'Staurosporine','SPADE','unstimulated'),
                                  files=paste(sampleNames(staurosporine)[unstimulated.samples],'density.fcs.cluster.fcs',sep='.') )

unstim.annots <- annots[unstimulated.samples,]
unstim.annots$name <- paste(unstim.annots$name,'density.fcs.cluster.fcs',sep='.')
row.names(unstim.annots) <- unstim.annots$name

pData(clusterUnstimSet) <- unstim.annots[sampleNames(clusterUnstimSet),]

# extract raw data
cluster.mat <- fsApply(clusterUnstimSet,exprs)
cluster.mat <- cluster.mat[,with(pData(parameters(clusterUnstimSet[[1]])),name[desc %in% c(pheno.markers,func.markers)])]
dimnames(cluster.mat)[[2]] <- with(pData(parameters(clusterUnstimSet[[1]])),desc[desc %in% c(pheno.markers,func.markers)])
cluster.mat <- asinh(cluster.mat/5)

cells <- lapply(sampleNames(clusterUnstimSet),function(x) rep(pData(clusterUnstimSet)[x,'Cells'],
                                                              nrow(clusterUnstimSet[[x]]) )
)
cells <- unlist(cells)

stims <- lapply(sampleNames(clusterUnstimSet),function(x) rep(pData(clusterUnstimSet)[x,'Stimulation'],
                                                                     nrow(clusterUnstimSet[[x]]) )
)
stims <- unlist(stims)
stims <- factor(stims, levels=stimulations)
stims <- relevel(stims,'Reference')

clusters <- lapply(sampleNames(clusterUnstimSet),function(x) exprs(clusterUnstimSet[[x]])[,'cluster'])
clusters <- unlist(clusters)

# check
heatmap(table(clusters,cells),scale='none')
heatmap(table(clusters,stimulations),scale='none')

names(layout) <- c('x_fixed','y_fixed')

nodes <- lapply(sort(unique(clusters)),function(cur.clus) {
    cat('Cluster:',cur.clus,'\n')
    res <- lapply(levels(stims),function(cur.stim) {
        #cat(cur.clus," | ",cur.stim,'\n')
        cur.mat <- cluster.mat[stims==cur.stim & clusters==cur.clus,,drop=FALSE]
        if(nrow(cur.mat)>1) {
          vals <- as.list(apply(cur.mat,2,median))
        } else if(nrow(cur.mat)==1) {
          vals <- as.list(cur.mat)
          names(vals) <- dimnames(cluster.mat)[[2]]
        } else {
          vals <- as.list(setNames(rep(NA,ncol(cluster.mat)),dimnames(cluster.mat)[[2]]))
        }
        vals$percenttotal <- sum(stims==cur.stim & clusters==cur.clus )/sum(stims==cur.stim)
        return(vals)
      }
    )
    names(res) <- levels(stims)
    res <- c(res,as.list(layout[cur.clus,]))
    res$name <- cur.clus
    return(res)
  }
)

edges <- as.data.frame(get.edgelist(graph))
names(edges) <- c('source','target')
edges$source <- as.integer(as.character(edges$source))-1
edges$target <- as.integer(as.character(edges$target))-1
edges$weight <- get.edge.attribute(graph,'weight')

spade <- list('nodes'=nodes,
              'links'=apply(edges,1,as.list), #edges,
              'reagents'=c(pheno.markers,func.markers),
              'stimulations'=levels(stims)
)

cat(toJSON(spade),file=file.path('web','spade.json'))
