#' @title package
#' @export

d<- function(vector){
  vector<- unlist(strsplit(vector, ","))
}

#' @title seurat featureplot loop
#' @export
d.featureplot<- function(data, feature, split.by, max.cutoff, pt.size, scale, order){

  if(missing(max.cutoff)){max.cutoff =  "q90"}
  if(missing(pt.size)){pt.size = 1}
  if(missing(split.by)){split.by = "sample_type"}
  if(missing(scale)){scale = "feature"}
  if(missing(order)){order = TRUE}

  feature<- fun.l(feature)
  pal<- colorRampPalette(c("#3983A0","lightgreen", "yellow", "red"))

  x<-for (i in seq_along(feature)){

    t = try(FeaturePlot(data, features = feature[i]),silent = TRUE)

    if(isTRUE(class(t) ==  "try-error")) {
      print(feature[i])
      next }else{ print(i)

        print(FeaturePlot(data, features = feature[i], cols = pal(10), max.cutoff = max.cutoff, pt.size = pt.size, order = order, split.by = split.by, keep.scale = scale)&
                theme(title = element_text(size = 25),
                      axis.title = element_text(size = 20),
                      legend.text = element_text(size = 15)))
      }
  }
}

#' @title Seurat dotplot loop
#' @export
d.dotplot<- function(data, feature, ident, cols, split.by){

  feature<- fun.l(feature)

  #need to have a complete list of idents to work, so adds the missing ones (checks the whole ident vector vs the ident vector) back but orders the specified ones
  ident<- fun.l(ident)
  whole.ident<- levels(data)
  x<- ident
  if(!all(x == whole.ident)){x<- append(x, whole.ident)}
  x<- x[!duplicated(x)]
  data@active.ident<- factor(data@active.ident, levels = x)

  #for generting the plot titles
  w<- strsplit(ident, "[.]") # sqaure brackets to allow me to split by "." as otherwise "." is treated as any character
  w<- unlist(w)

  #checks more than just the first ident location: splits the string vector, takes only the even number (back/dorsal...), deletes duplicated values and reoders the string with "&"
  x = 0
  w.II<- NULL

  n<- 1:length(w)
  n<- subset(n, n %% 2 == 0)

  for(i in seq_along(n)){w.II<- append(w.II,w[n[i]])}
  w.II<- w.II[!duplicated(w.II)]

  n<- seq_along(w.II)
  if(sum(n)>1){
    #bad fix, to be more robust implement inside a for loop
    t<- paste(w.II, w.II[1], sep = " & ")
    if(sum(n) > 3){ t<- paste(t[2], w.II[3], sep = " & ")}else{t<- t[2]}

  }else{t = w.II[1]}

  if(missing(cols)){cols = fun.l("#3983A0,white,red")}else{cols<- fun.l(cols)}

  DotPlot(object = data, features = feature, idents = ident, dot.scale = 7) + scale_colour_gradient2(low = cols[1], mid = cols[2], high = cols[3])+
    scale_shape_manual(values = "grey")+#outline code
    theme_bw()+
    theme(panel.grid= element_blank(),
          panel.background = element_rect(colour = "black", size= NULL),
          axis.text.x = element_text(angle = 90))+
    labs(title = t)+
    geom_point(aes(size=pct.exp), shape = 21, colour="#525659", stroke=0.5)

}

#' @title Feature plot in ggplot2
#' @export
d.ggfeatureplot<- function(data, feature, split.by, reduction, max.cutoff){
  feature<- fun.l(feature)
  if(missing(reduction)){reduction = "UMAP"}
  if(missing(max.cutoff)){max.cutoff = 1}

  for(i in seq_along(feature)){

    t = try(FeaturePlot(data, features = feature[i]),
            silent = TRUE)

    if(isTRUE(class(t) ==  "try-error")) {
      print(feature[i])
      next
    }else{

      print(i)

      if(!missing(split.by)){
        x<- FetchData(data, vars = paste(reduction,"_1", sep = ""))
        y<- FetchData(data, vars = paste(reduction,"_2", sep = ""))
        z<- FetchData(data, vars = feature[i])
        split<- FetchData(data, vars = split.by)

        df<- data.frame(x,y,z,split)

        names(df)[1]<- "Reduction_1"
        names(df)[2]<- "Reduction_2"
        names(df)[3]<- "Expression"
        names(df)[4]<- "split"

        cut<- quantile(df$Expression, max.cutoff) #calculates quartile
        df$Expression<- ifelse(df$Expression > cut, max(df$Expression), df$Expression) #replaces values above the max.cutoff line with 0

        df<- df[order(df$Expression, decreasing = FALSE),]

        p<- ggplot(df, aes(Reduction_1,Reduction_2, colour = Expression))+
          geom_point()+
          theme_minimal()+
          scale_colour_gradientn(colors = c("#3983A0","lightgreen", "yellow", "red"))+
          ggtitle(feature[i])+
          xlab(paste(reduction,"_1", sep = ""))+
          ylab(paste(reduction,"_2", sep = ""))+
          facet_wrap(vars(split))

        print(p)
      }else{
        x<- FetchData(data, vars = paste(reduction,"_1", sep = ""))
        y<- FetchData(data, vars = paste(reduction,"_2", sep = ""))
        z<- FetchData(data, vars = feature[i])

        df<- data.frame(x,y,z)

        names(df)[1]<- "Reduction_1"
        names(df)[2]<- "Reduction_2"
        names(df)[3]<- "Expression"

        cut<- quantile(df$Expression, max.cutoff) #calculates quartile
        df$Expression<- ifelse(df$Expression > cut, max(df$Expression), df$Expression)

        #df<- df[df$Expression < cut,]

        df<- df[order(df$Expression, decreasing = FALSE),]#so most expressed nuceli are on top

        p<- ggplot(df, aes(Reduction_1,Reduction_2, colour = Expression))+
          geom_point()+
          theme_minimal()+
          scale_colour_gradientn(colors = c("#3983A0","lightgreen", "yellow", "red"))+
          ggtitle(feature[i])+
          xlab(paste(reduction,"_1", sep = ""))+
          ylab(paste(reduction,"_2", sep = ""))

        print(p)
      }
    }
  }
}

#' @title 3d plot baaaaby
#' @export
d.3d_featureplot<- function(data, feature, reduction, type, ident, size, max.cutoff, point){
  if(missing(reduction)){reduction = "UMAP"}
  if(missing(type)){type = 1}
  if(missing(ident)){ident = "umap.idents"}
  if(missing(size)){size = 3}
  if(missing(max.cutoff)){max.cutoff = 1}
  if(missing(point)){point = "p"}

  x<- FetchData(data, vars = paste(reduction,"_1", sep = ""))
  y<- FetchData(data, vars = paste(reduction,"_2", sep = ""))
  z<- FetchData(data, vars = feature)
  ident<- FetchData(data, vars = ident)

  df<- data.frame(x,y,z,ident)
  names(df)[1]<- "reduction_1"
  names(df)[2]<- "reduction_2"
  names(df)[3]<- "feature"
  names(df)[4]<- "ident"

  cut<- quantile(df$feature, max.cutoff)#max.cutoff argument
  df$feature<- ifelse(df$feature > cut, max(df$feature), df$feature)#NEEDS FIXING, turns max values white

  if(type == 1){ #each ident = z axis
    df$colour<- as.numeric(as.character(df$feature))#copies ptch1 data to the colours column
    range<- range(df$colour)
    range<- seq(from = range[1]-0.01, to = range[2]+0.5, by = max(range)/10) #-0.01 to the min range value to ensure that 0 values have a col, wont always meet max val so add +0.5

    cols<- colorRampPalette(c("#3983A0","lightgreen", "yellow", "red"))
    cols<- cols(length(range)-1)

    df$colour<- cut(df$colour, breaks = range, labels = cols)#allow me to assign a discrete value to a continious spectrum

    plot3d(x = df$reduction_1, y=df$reduction_2, z = df$ident, col = df$colour, size = size, type = point,
           xlab = paste(reduction,"_1", sep = ""), ylab = paste(reduction,"_2", sep = ""), zlab = "Cluster")

    #Title stuff
    bgplot3d({
      plot.new()
      title(main = feature, line = 2)
    })
  }else{ #gene expression = z axis

    df$colour<- as.numeric(as.character(df$ident))

    identities<- df$ident
    identities<- identities[!duplicated(identities)]
    identities<- identities[order(identities, decreasing = FALSE)]


    cols<- colorRampPalette(fun.l("#4F2577,#1F4995,#0090A8,#229548,#8EB71B,#F08080,#FFD700,#00FF00,#00BFFF,#FF69B4,#708090,#800000,#9400D3,#00CED1,#3CB371,#FF8C00,#FFDAB9,#98FB98,#F3E500"))
    cols<- cols(length(identities))
    df$colour<- factor(df$colour, levels = identities)
    levels(df$colour)<- cols  #now each ident is now its assigned colour

    plot3d(x = df$reduction_1, y=df$reduction_2, z = df$feature, col = df$colour, size = size, type = point,
           xlab = paste(reduction,"_1", sep = ""), ylab = paste(reduction,"_2", sep = ""), zlab = "Expression")

    #Title stuff
    bgplot3d({
      plot.new()
      title(main = feature, line = 2)
    })
  }
}

#' @title Pseudotime plot
#' @export
d.pseudotimeplot<- function(data1, data2, feature, split.by, points, lineage){
  gene<- FetchData(data1, vars = feature)
  if(missing(data2)){data2<- sce}
  if(missing(points)){points = 1}
  if(missing(lineage)){lineage = 1}

  lineage<- paste0(lineage)#to add quotation marks to it
  lineage<- paste("slingPseudotime_",lineage,sep = "")
  df<- data.frame(sce[[lineage]],gene)
  names(df)[1]<- "time"
  names(df)[2]<- "gene"

  if(missing(split.by)){

    ggplot(df, aes(x = time, y = gene))+
      geom_point(alpha = points)+
      geom_smooth()+
      xlab("Pseudotime ordering")+
      ylab(paste(feature,"expression", sep = " "))

  }else{

    ident<- FetchData(data1, vars = split.by)
    df<- data.frame(df,ident)
    names(df)[3]<- "Identity"

    ggplot(df, aes(x = time, y = gene))+
      geom_point(alpha = points)+
      geom_smooth()+
      ggtitle("Aggregated")+
      xlab("Pseudotime ordering")+
      ylab(paste(feature,"expression", sep = " "))+

      ggplot(df, aes(x = time, y = gene, colour = Identity))+
      geom_point(alpha = points)+
      geom_smooth()+
      xlab("Pseudotime ordering")+
      ylab(paste(feature,"expression", sep = " "))+
      facet_wrap(vars(Identity))

  }
}

#' @title qausi compoent plot
#' @export
d.quasi_componentplot<- function(data, comp1, comp2, split.by, col){
  if(missing(col)){if(!missing(split.by)){col = split.by}else{col = "sample_type"}} #needs fixing

  #data<- NormalizeData(data)
  data<- FindVariableFeatures(data)
  #data<- ScaleData(data, features = rownames(data))

  genes<- VariableFeatures(data)
  df<- FetchData(data, vars = genes) #to keep only interesting genes
  df<- t(df)
  lef1.exp<- FetchData(data, vars = comp1)
  ptch1.exp<- FetchData(data, vars = comp2)

  #want to find genes, across the sample, that are most correlated -> rows in the dataframe

  ####correlation test####
  x.lef<- NULL
  y.ptch<- NULL
  for(i in seq_along(df[,1])){
    print(i)
    x<- cor(lef1.exp, df[i,], method = "pearson")
    y<- cor(ptch1.exp, df[i,], method = "pearson")

    x.lef<- append(x.lef, x)
    y.ptch<- append(y.ptch, y)
  }
  names<- names(df[,1])
  ####sort the output and take only the top 30####
  df<- data.frame(names,x.lef,y.ptch)
  df<- df[order(df$x.lef, decreasing = TRUE),]
  lef.genes<- head(df$names,30)
  df<- df[order(df$y.ptch, decreasing = TRUE),]
  ptch.genes<- head(df$names,30)

  ####take the top 30 most correlated genes, use PCA on these genes only, plot####
  data<- RunPCA(data, features = lef.genes)
  pc.lef<- FetchData(data, vars = "PC_1")

  data<- RunPCA(data, features = ptch.genes)
  pc.ptch<- FetchData(data, vars = "PC_1")

  if(!missing(split.by)){
    df<- data.frame(pc.lef, pc.ptch, data[[split.by]],FetchData(data, vars = col))
    names(df)[1]<- "wnt_component"
    names(df)[2]<- "sonic_component"
    names(df)[3]<- "split"
    names(df)[4]<- "col"

    ggplot(df, aes(x = wnt_component, y = sonic_component, colour = col))+
      geom_point()+
      facet_wrap(vars(split))+
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)+
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
      xlab(label = paste(comp1,"Component", sep = " "))+
      ylab(label = paste(comp2,"Component", sep = " "))+
      labs(colour = col)
  }else{
    df<- data.frame(pc.lef, pc.ptch, FetchData(data, vars = col))
    names(df)[1]<- "wnt_component"
    names(df)[2]<- "sonic_component"
    names(df)[3]<- "col"

    ggplot(df, aes(x = wnt_component, y = sonic_component, colour = col))+
      geom_point()+
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)+
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)+
      xlab(label = paste(comp1,"Component", sep = " "))+
      ylab(label = paste(comp2,"Component", sep = " "))+
      labs(colour = col)
  }
}
