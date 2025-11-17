library(scales)

# ENTER YOUR DIRECTORY PATH---------
PathToData=""
# define function-----------
iClassification.class.data.to.subset.cols <- function (data, col.name, values){
  data.subset = data$col.annot[,col.name] %in% values
  if(sum(data.subset)==1){ #cast data values as matrices in case number of cols = 1 (otherwise the matrices become vectors)
    if(nrow(data$row.annot)==1){
      stop(paste("[iClassification.class.data.to.subset.cols]: can't subset to a single column if data has only 1 row!", sep=""), call. = TRUE)	
    }
    # subset the data		
    for(i in 1:length(data$data)){
      data$data[[i]] = matrix(data$data[[i]][,data.subset], ncol=1)				
    }		
    # subset the annotation		
    data$col.annot = data$col.annot[data.subset,]				
  }else{	
    # subset the data		
    for(i in 1:length(data$data)){
      data$data[[i]] = data$data[[i]][,data.subset]				
    }		
    # subset the annotation		
    data$col.annot = data$col.annot[data.subset,]				
  }
  
  return(data)
}

# CO PCA----------
data = get(load(paste0(PathToData,"/updated_organoids.RData")))
data = iClassification.class.data.to.subset.cols(data, "origin", "CO") 
color.points=(data$col.annot$color)
names(color.points)=data$col.annot$genotype
exprs=data$data$quantile


exprs=data$data$quantile
pc=prcomp(t(exprs))
standardError=rep(1.1,length(unique(data$col.annot$genotype)))
names(standardError)=unique(data$col.annot$genotype)
for (i in unique(data$col.annot$genotype)) {
  if(sum(data$col.annot$genotype==i)>1){
  x=pc$x[data$col.annot$genotype==i,1:2]
  standardError[i]=(sd(x[,1])^2/nrow(x)+sd(x[,2])^2/nrow(x))^(0.5)
  }
}
exprs=exprs[,1:length(unique(data$col.annot$genotype))]
colnames(exprs)=unique(data$col.annot$genotype)

for (i in unique(data$col.annot$genotype)) {
  if(sum(data$col.annot$genotype==i)>1){
  x=data$data$quantile[,data$col.annot$genotype==i]
  exprs[,i]=apply(x, 1, mean)
  }else{
    x=data$data$quantile[,data$col.annot$genotype==i]
    exprs[,i]=x
  }
}
pc=prcomp(t(exprs))
plot(pc$x[,1:2],col=alpha(color.points[rownames(pc$x)], 0.4), pch=19,
     cex=(standardError)/3
     )
text(pc$x[,1:2],col=alpha("gray10", 1), pch=19, rownames(pc$x))

