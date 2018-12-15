
# install required libraries
installLibraries <- function(){
  library(limma)
  library(edgeR)
  library(gplots)
  library(RColorBrewer)
  source('https://bioconductor.org/biocLite.R')
  biocLite('org.Hs.eg.db')
  library('org.Hs.eg.db')
  library(Rsubread)
  library(Glimma)

}

#make counts file
makeCounts <- function(files){ # files <- list of bam files
  fileNames <- c()
  for(i in 1:length(files)){
    
    name <- paste( sub('\\..*', '', sub('.*bam.','',files[i])),sep = "", ".txt")
    fileNames <- name
    counts <- featureCounts(files= files[i], annot.inbuilt = "hg19",GTF.featureType="exon", GTF.attrType="gene_id",useMetaFeatures=TRUE, nthreads = 12)
    sink(name) # create file
    print(counts)
    sink() # save
    }
}

# sub function for preproscessing to assign geneID and geneName
geneNameCheck <- function(geneNames, rowName){
  x = 1
  for(i in geneNames){
    if(is.na(i)){
      geneNames[ as.character(rowName[x]) ] = as.character(rowName[x])
    }
    x = x+1
  }
  geneNames <<- geneNames
  return(geneNames)
}

# create DataFrame from bam files
createCountData <- function(numberOfGenes, countsFiles) { # how many genes(rows), 
  #find a way to get only one column
  run = TRUE
  rowName <-c()
  
  for(i in 1:length(countsFiles) ){
    # get column Names
    colName <- sub('*.txt','', countsFiles[i] )
    seqdata <- read.delim(file = countsFiles[i], stringsAsFactors = FALSE) # read count file
    
    #for first time create rows from count file
    # assign geneID to geneName
    if(run) { # create row names
      processedData <- data.frame(matrix(seqdata, nrow=numberOfGenes, ncol=0))
      rowName <- as.numeric( sub( ' .*', '', seqdata[1:numberOfGenes+1,]) )
      geneNames <- mapIds(org.Hs.eg.db,as.character(rowName), 'SYMBOL', 'ENTREZID') # assign geneID to geneName
      rownames(processedData) <- geneNameCheck(geneNames, rowName) # check if all geneID succesfully assigned and use it for rownames
      run = FALSE
    }
    col = seqdata[2:(numberOfGenes+1), ]
    col =sub('.* ','', col) 
    processedData[colName] = as.numeric(col)
    
  }
  Mylabels <<- paste(sampleinfo$FileName, sampleinfo$CellType, sampleinfo$Status)
  group <<- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
  group <<- factor(group)
  print("Created globa variabled Mylabels and group")
  return(processedData)
}


# preprocess data(takes counts and number of genes)
preprocessData <- function(numberOfGenes, seqdata, numberOfSamples) {
  #find a way to get only one column
  processedData <- data.frame(matrix(seqdata, nrow=numberOfGenes, ncol=0))
  run = TRUE
  rowName <-c()
  
  for(i in 0:numberOfSamples-1) {
    print("entered")
    # get column Names
    colName <- sub('.*bam.','', seqdata[ (1+(numberOfGenes+1)*i),] )
    print(colName)
    colName <- sub('\\..*', '', colName)
    
    
    # assign geneID to geneName
    if(run) { # create row names
      rowName <- as.numeric(sub(' .*', '', seqdata[(2+i*(numberOfGenes+1)):((numberOfGenes+1)+i*(numberOfGenes+1)),]))
      geneNames <- mapIds(org.Hs.eg.db,as.character(rowName), 'SYMBOL', 'ENTREZID') # assign geneID to geneName
      rownames(processedData) <- geneNameCheck(geneNames, rowName) # check if all geneID succesfully assigned and use it for rownames
      run = FALSE
    }
    
    col = sub('.* ','', seqdata[ ( 2+i*(numberOfGenes+1) ):( (numberOfGenes+1)+i*(numberOfGenes+1) ),] )
    col = seqdata[(2+i*(numberOfGenes+1)):((numberOfGenes+1)+i*(numberOfGenes+1)),]
    col =sub('.* ','', col) 
    processedData[colName] = as.numeric(col)
    
  }
  labels <<- paste(sampleinfo$FileName, sampleinfo$CellType, sampleinfo$Status)
  group <<- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
  group <<- factor(group)
  glMDSPlot(y, labels=labels, groups=group, folder="mds")
  print("Created globa variabled labels and group")
  return(processedData)
}

# countdata <- data frame of genes, rate <- rate of removing
removeLowlyExpressedGenes <- function(countdata, rate = 0.5){
  View(countdata)
  #check for numeric values in DF
  myCPM <- cpm(countdata)
  thresh <- myCPM > rate
  keep <- rowSums(thresh) >= 2 # used to be 2
  counts.keep <- countdata[keep,]
  y <<- DGEList(counts.keep)
  logcounts <- cpm(y,log=TRUE)
 # y <<- DGEList(counts.keep, norm.factors = rep(1,ncol(counts.keep)))
 # logcounts <- cpm(y,log=TRUE,prior.count=2)
  print("Global variable 'y' was created")
  View(logcounts)
  return(logcounts)
  
}

#build hierarchical Clustering heatmap (dataFrame, colors, numberOfSamples, method of clustering, decreasing, reorder Columns, build all possible graphs, give name for graph file)
hierarchicalClusteringHeatmap <- function(logCounts, colors, numberOfSamples, method = 'Default', decreasing = TRUE, Colv =TRUE, buildAll = FALSE, name = ""){
  # methods :
  #            :default
  #           1: 1-Pearson
  #           12: 1-Pearson other
  #           E: Euclidean
  #           E2: Euclidean other
  #           T: Test option
  #           
  par(mfrow=c(1,2), mar=c(5,4,10,2))
  mypalette <- brewer.pal(11,"RdYlBu")
  morecols <- colorRampPalette(mypalette)
  
#  for(i in 2:numberOfSamples){
   # stri = sprintf("Top %d  most variable \n genes across samples", i)
  stri = sprintf("Top %d  most variable \n genes across samples", numberOfSamples)
    var_genes <- apply(logCounts , 1, var)
   # select_var <- names(sort(var_genes, decreasing = decreasing))[1:i]
    select_var <- names(sort(var_genes, decreasing = decreasing))[1:numberOfSamples]
    highly_variable_lcpm <<- logCounts[select_var,]
    
    #build all possible graphs
    if(buildAll == TRUE){
      png(file = paste(name,sep = "_", paste(numberOfSamples, sep = "", "_1.png") ) )
      heatmap.2(highly_variable_lcpm, col=rev(morecols(2)), Colv = Colv, offsetRow=0, offsetCol = -0.2, 
                cexCol = 0.6, trace="none", main=stri,ColSideColors=colors,scale="row",
                reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
                distfun=function(x) as.dist(1-cor(t(x))),
                hclustfun=function(x) hclust(x, method="complete"))
      dev.off()
      
      png(file = paste(name,sep = "_", paste(numberOfSamples, sep = "", "_E.png") ) )
      heatmap.2(highly_variable_lcpm, col=rev(morecols(50)) , Colv = Colv, offsetRow=0, offsetCol = -0.2, 
                cexCol = 0.6, trace="none", main=stri,ColSideColors=colors,scale="row",
                reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
                distfun=function(y) dist(y, method="euclidean"),
                hclustfun=function(y) hclust(y, method="ward.D2"))
      dev.off()
      
      png(file = paste(name,sep = "_", paste(numberOfSamples, sep = "", "_T.png") ) )
      heatmap.2(highly_variable_lcpm, col=rev(morecols(50)),trace="none", Colv = Colv, offsetRow=0, 
                cexCol = 0.6, main=stri, ColSideColors=colors)
      dev.off()
      
      png(file = paste(name,sep = "_", paste(numberOfSamples, sep = "", "_E2.png") ) )
      heatmap.2(highly_variable_lcpm, col=rev(morecols(50)) , Colv = Colv, offsetRow=0, offsetCol = -0.2, 
                cexCol = 0.6, trace="none", main=stri,ColSideColors=colors,scale="row",
                distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="ward.D2"))
      dev.off()
      
      png(file = paste(name,sep = "_", paste(numberOfSamples, sep = "", "_12.png") ) )
      heatmap.2(highly_variable_lcpm, col=rev(morecols(50)) , Colv = Colv, offsetRow=0, offsetCol = -0.2, 
                cexCol = 0.6, trace="none", main=stri,ColSideColors=colors,scale="row",
                distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="average"))
      dev.off()
      
      png(file = paste(name,sep = "_", paste(numberOfSamples, sep = "", "_Default.png") ) )
      heatmap.2(highly_variable_lcpm, col=rev(morecols(50)), Colv = Colv, 
                cexCol = 0.6, trace="none", main=stri, ColSideColors=colors,scale="row")
      dev.off()
   
    }
    # build only spesific graph
    else{
      if(method == '1'){ # 1-Pearson
        heatmap.2(highly_variable_lcpm, col=rev(morecols(2)), Colv = Colv, offsetRow=0, offsetCol = -0.2, 
                cexCol = 0.6, trace="none", main=stri,ColSideColors=colors,scale="row",
                reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
                distfun=function(x) as.dist(1-cor(t(x))),
                hclustfun=function(x) hclust(x, method="complete"))
      }
      else if (method == 'E'){ # Euclidean
        heatmap.2(highly_variable_lcpm, col=rev(morecols(50)) , Colv = Colv, offsetRow=0, offsetCol = -0.2, 
                  cexCol = 0.6, trace="none", main=stri,ColSideColors=colors,scale="row",
                  reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
                  distfun=function(y) dist(y, method="euclidean"),
                  hclustfun=function(y) hclust(y, method="ward.D2"))
      }
      else if (method == "T"){ # Test
        heatmap.2(highly_variable_lcpm, col=rev(morecols(50)),trace="none", Colv = Colv, offsetRow=0, 
                  cexCol = 0.6, main=stri, ColSideColors=colors)
      }
      else if(method == "E2"){
        heatmap.2(highly_variable_lcpm, col=rev(morecols(50)) , Colv = Colv, offsetRow=0, offsetCol = -0.2, 
                  cexCol = 0.6, trace="none", main=stri,ColSideColors=colors,scale="row",
                  distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="ward.D2"))
      }
      else if(method == "12"){
        heatmap.2(highly_variable_lcpm, col=rev(morecols(50)) , Colv = Colv, offsetRow=0, offsetCol = -0.2, 
                  cexCol = 0.6, trace="none", main=stri,ColSideColors=colors,scale="row",
                  distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="average"))
      }
      else{
        heatmap.2(highly_variable_lcpm, col=rev(morecols(50)), Colv = Colv, 
                  cexCol = 0.6, trace="none", main=stri, ColSideColors=colors,scale="row")
      }
    }
    
  View(highly_variable_lcpm)
}

# build voom transoform data and create global variable v (page36)
voomTransoftData <- function(group, indexY, normalize.method = "none"){
  indexY <- calcNormFactors(indexY)
  indexY$samples
  par(mfrow=c(1,1))
  print(levels(group))
  design <- model.matrix(~ 0 + group)
  #design <- design[,indexForLevel] 
  colnames(design) <- levels(group)
  
  par(mar=c(5,5,5,5))
  v <- voom(indexY,design,plot = TRUE)
 # v <- voom(indexY,design,plot = TRUE, normalize.method = normalize.method)
  
  print("Global var Y was changed")
  print("Global var Design was created")
  indexY <<- indexY
  Design <<- design
  return(v)
}

#build to graphs unnormalise and normalise(page 39)
normalize <- function(logcounts,v, ylim = c(5,15), cex.axis = 0.7 ) { # (dataframe, v, ylim(axes))
  par(mfrow=c(1,1), mar=c(5,5,4,2), cex.axis = cex.axis)
  boxplot(logcounts, xlab="", ylim=ylim,  ylab="Log2 counts per million", las=2,main="Unnormalise d logCPM",offsetRow=0, offsetCol = -0.2)
  abline(h=median(logcounts),col="blue")
  boxplot(v$E, xlab="",ylim=ylim, ylab="Log2 counts per million",las=2,main="Voom transformed logCPM",offsetRow=0, offsetCol = -0.2)
  abline(h=median(v$E),col="blue")
  
}


#testing for differential expression (page 39)
testingForDifferentialExpression <- function(v, design, name, displayPlots =TRUE){ # voom for sample, design for sample, name for plots
  # TODO find a way to spesify my own name and parameters
  #                       in cont.matrix(B.Scrambled = HeLa.Scrambled, L.shoct4 = HeLa.shoct4, levels=design)
  
  # returns summa.fit for sample
  fit <- lmFit(v)
  names(fit)
  cont.matrix <- makeContrasts(Scrambled = HeLa.shoct4 - HeLa.Scrambled, levels=design)
  fit.cont <- contrasts.fit(fit, cont.matrix)
  fit.cont <- eBayes(fit.cont)
  print("created fit.cont")
  summa.fit <- decideTests(fit.cont)
  
  # VennDiagram Plot
  plotName = paste(name, sep = "_vennDiagram", ".png")
  vennDiagram(summa.fit,include=c("up", "down"),
              counts.col=c("red", "blue"),
              circle.col = c("red", "blue", "green3"),
              cex = 1)
  title(name)
  dev.off()
  
  if(displayPlots){
    vennDiagram(summa.fit,include=c("up", "down"),
                counts.col=c("red", "blue"),
                circle.col = c("red", "blue", "green3"),
                cex = 1)
    title(name)
  }
  
  topTable(fit.cont,coef=1,sort.by="p")
  
  limma.res <- topTable(fit.cont,coef="Scrambled",sort.by="p",n="Inf")
  write.csv(limma.res,file=paste(name, sep = "", "Shoct-scrambled.csv") )
  
 # limma.res <- topTable(fit.cont,coef="shoct4",sort.by="p",n="Inf")
 # write.csv(limma.res,file="shoct4.csv")
  
  # PlotMD 1
  png(paste(name,sep = "_Scrambled", ".png"))
  plotMD(fit.cont,coef=1,status=summa.fit[,"Scrambled"], values = c(-1, 1), title("Scrambled"))
  dev.off()
  
  if(displayPlots){
    plotMD(fit.cont,coef=1,status=summa.fit[,"Scrambled"], values = c(-1, 1), title("ScrambledVsShoct"))
  }
  
  # Volcano Plot
  
  png(paste(name,sep = "_volcano_shOct4-Scrambled", ".png"))
  volcanoplot(fit.cont,coef=1,highlight=200,names=rownames(fit.cont),ylim=c(0,3) )
  title(name)
  dev.off()
  
  if(displayPlots){
   # cont.matrix <- makeContrasts(Scrambled = C33A.Scrambled - C33A.shoct4, shoct4 = C33A.shoct4-C33A.Scrambled, levels=design)
    volcanoplot(fit.cont,coef=1,highlight=20,names=rownames(fit.cont), ylim=c(0,4) )
    title(name)
  }
  
  # PlotMD 2
  png(paste(name,sep = "_shoct", ".png"))
 # plotMD(fit.cont,coef=1,status=summa.fit[,"shoct4"], values = c(-1, 1), title("shoct4"))
  dev.off()
  
  if(displayPlots){
#    plotMD(fit.cont,coef=2,status=summa.fit[,"shoct4"], values = c(-1, 1), title("shoct4"))
  }
  
  return(summa.fit)
}




#steps from tutorial
tutorialExample <- function(){
  seqdata <- read.delim("mouseCounts.txt")
  sampleinfo <- read.delim("SampleInfo.txt")
  countdata <- seqdata[,-(1:2)]
  View(countdata)
  rownames(countdata) <- seqdata[,1]
  colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
  table(colnames(countdata)==sampleinfo$SampleName)
  
  View(countdata)
  
  # remove Lowly exper
  myCPM <- cpm(countdata)
  thresh <- myCPM > 0.5
  keep <- rowSums(thresh) >= 2
  counts.keep <- countdata[keep,]
  y <- DGEList(counts.keep)
  logcounts <- cpm(y,log=TRUE)
  
  View(logcounts)
  
  col.cell <- c("purple","orange")[sampleinfo$CellType]
  col.status <- c("blue","red","dark green")[sampleinfo$Status]
  sampleinfo <- read.delim("SampleInfo_Corrected.txt")
  
  col.cell <- c("purple","orange")[sampleinfo$CellType]
  col.status <- c("blue","red","dark green")[sampleinfo$Status]
  labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Status)
  group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
  group <- factor(group)
  
  # hierarchical Clustering Heatmap
  var_genes <- apply(logcounts, 1, var)
  select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
  highly_variable_lcpm <- logcounts[select_var,]
  mypalette <- brewer.pal(11,"RdYlBu")
  morecols <- colorRampPalette(mypalette)
  col.cell <<- c("purple","orange")[sampleinfo$CellType]
  heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
  
  print("Global variables created : ")
  print("col.cell - for colors in heatmap ")
        
  
  return(logcounts)
}


