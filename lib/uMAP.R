
require("umap")
require("data.table") # for the shift function

umap_main <- function(data_path, filename, lam, maximal_cluster_number, on_off_smoothing){

  set.seed(123)
  
  print(data_path)
  print(filename)
  print(lam)
  print(maximal_cluster_number)
  print(on_off_smoothing)
  
  ############
  ##  Input ##
  ############
  
  if (dir.exists(paste0(data_path, "/clustering_", filename)) == FALSE ){
    dir.create(paste0(data_path, "/clustering_", filename))
  }
  output_path <- paste0(data_path, "/clustering_", filename, "/")
  
  # get data
  data <- read.table(paste0(data_path, "/data_classification_", filename, ".tsv"), header=FALSE)
  data2 <- matrix(as.numeric(unlist(data)), nrow=nrow(data), ncol=ncol(data))
  rownames(data2) <- seq(1, nrow(data2))
  
  outputmatrix <- read.table(paste0(data_path, "/CV_tab_", filename , ".bed"), header=FALSE)
  
  ################
  ##  Normalize ##
  ################
  
  min_max_noramlize <- function(x){
    x_new <- (x-min(x)) / (max(x)-min(x)) 
    return(x_new)
  }
  
  data_normalized <- t(apply(data2, 1, min_max_noramlize))
  
  # Set NaNs to zero
  data_normalized[which(is.na(data_normalized))] <- 0
  
  # Remove duplicates for uMAP
  duplicates <- which(duplicated(data_normalized))
  data_removed_duplicates <- data_normalized[which(duplicated(data_normalized) == FALSE),]
  
  ############################
  ##  Check for Translocate ##
  ############################
  
  end <- ncol(data_removed_duplicates)
  
  left <- which(apply(data_removed_duplicates[,1:3], 1, mean) < 0.8)
  right <- which(apply(data_removed_duplicates[,(end-2):end], 1, mean) < 0.8)
  
  peaks_to_translocate <- intersect(left, right)
  
  #################
  ##  Tranlocate ##
  #################
  
  translocate_forward <- function(v, s){
    end <- length(v)
    v_new <- shift(v, s)
    v_new[1:s] <- v[(end-s+1):end]
    return(v_new)
  }
  
  translocate_backward <- function(v, s){
    end <- length(v)
    v_new <- shift(v, s)
    v_new[(end-abs(s)+1):end] <- v[1:abs(s)]
    return(v_new)
  }
  
  center_peak <- function(y_peak, x_dist){
    s <- which.max(convolve(x_dist , y_peak, type="open")) - length(y_peak)
    if ( s < 0 ){
      return(translocate_backward(y_peak, s))
    }
    if ( s > 0 ){
      return(translocate_forward(y_peak, s))
    } 
    return(y_peak)
  }
  
  # Center Peak
  # Use a Gaussian distribution to center the peaks around the middle of the diagram.
  dist <- dnorm(seq(-4, 4, 8/ncol(data_removed_duplicates)), mean=0, sd=1)
  dist <- dist/max(dist)
  d <- length(dist) - ncol(data_removed_duplicates)
  dist <- dist[1:ncol(data_removed_duplicates)]
  
  if ( d > 0 ){
    dist <- translocate_backward(dist, -d)
  }
  
  if ( d < 0 ){
    dist <- translocate_forward(dist, d)
  } 
  
  pre_data_ready <- t(apply(data_removed_duplicates, 1, center_peak, x_dist=dist))
  data_ready <- data_removed_duplicates
  data_ready[peaks_to_translocate, ] <- pre_data_ready[peaks_to_translocate, ]
  
  #################
  ##  Smooth Out ##
  #################
  
  smoothing <- function(y, lambda, dim){
    data_points <- length(y)
    
    # Test if we have just a vector of constant values.
    which(y != max(y))
    if ( length(which(y != max(y))) != 0 ){
      x <- c(1:data_points)
      y.smooth <- smooth.spline(x, y, spar=lambda)
      x_new <- seq(1, data_points, data_points/dim)
      y_new <- predict(y.smooth, x_new)$y
      y_new <- y_new/max(y_new)
      return(y_new)
    } 
    
    x_new <- seq(1, data_points, data_points/dim)
    y_new <- rep( max(y), length(x_new) )
    return(y_new)
  }
  
  smoothing_dim <- 150
  
  if ( ncol(data_ready) > smoothing_dim ) {
    smoothing_dim <- ncol(data_ready)
  }
    
  data_smoothed <- data_ready
  if( on_off_smoothing ) {
    print("[NOTE] Smooth data")
    data_smoothed <- t(apply(data_ready, 1, smoothing, lambda=lam, dim=smoothing_dim))
  }
  
  #######################
  ##  kmeans Functions ##
  #######################
  
  kmeansAIC <- function(fit){
    m <- ncol(fit$centers)
    k <- nrow(fit$centers)
    D <- fit$tot.withinss
    return(D + 2*m*k)
  }
  
  # function to find optimal number of clusters
  kopt <- function(data, k.max){
    
    num_centroids <- c(2:k.max) 
    
    kmeans <- sapply(num_centroids, function(k){kmeans(data, centers = k, nstart = 100, iter.max = 1000)})
    
    # get total within-cluster sum of squares
    twss <- numeric(ncol(kmeans))
    for ( i in 1:ncol(kmeans) ){
      twss[i] <- kmeans[,i]$tot.withinss
    }
    
    # get percentage of variance explained
    perc <- numeric(ncol(kmeans))
    for ( i in 1:ncol(kmeans) ){
      perc[i] <- (kmeans[,i]$betweenss / kmeans[,i]$totss) * 100
    }
    
    aic <- apply(kmeans, 2, kmeansAIC)
    
    opt <- which.min(aic) - 1
    
    differences <- numeric()
    for( i in 1:(length(twss)-1) ) {
      differences <- c(differences, abs(twss[i] - twss[i+1]))
    }
    percantage_diff <- differences/sum(differences)
    
    if ( which( percantage_diff < 0.1 )[1] < opt ) 
      opt <- which( percantage_diff < 0.1 )[1] + 2
    
    pdf(paste0(output_path, "/kmeans_Optimization.pdf"), 15, 15)
    par(cex = 2.0, mfrow = c(2,2), cex = 1.5, family = 'serif')
    
    plot(num_centroids, twss,
         type="b", pch = 19, frame = FALSE,
         xlab="Number of Clusters",
         ylab="Total within Cluster Sum of Squares",
         xlim=c(1,k.max))
    abline(v = opt, lty =2)
    
    plot(num_centroids, perc,
         type="b", pch = 19, frame = FALSE,
         xlab="Number of Clusters",
         ylab=" ",
         ylim = c(0,100),
         xlim=c(1,k.max))
    text( (opt - 0.5) , perc[opt], round(perc[opt], digits = 1))
    abline(v = opt, lty =2)
    
    plot(num_centroids, log(aic),
         type="b", pch = 19, frame = FALSE,
         xlab="Number of Clusters",
         ylab="AIC (log)",
         xlim=c(1,k.max))
    abline(v = opt, lty =2)
    dev.off()
    
    return(c(opt, perc[opt]))
  }
  
  #############
  ### uMAP  ###
  #############
  
  d <- 2
  custom.config = umap.defaults
  custom.config$n_epochs = 5000
  custom.config$n_components = d
  custom.config$min_dist = 0.01
  custom.config$n_neighbors = 5
  
  data_umap = umap(data_smoothed, config = custom.config)
  new_data <- data_umap$layout
  
  ###############
  ### kmeans  ###
  ###############
  
  returned_optimisation_values <- kopt(new_data, maximal_cluster_number)
  optimal_num_centroids <- returned_optimisation_values[1]
  optimal_perc <- returned_optimisation_values[2]
  
  ## use kmeans clustering with optimal number of centroids
  clusters <- kmeans(new_data, centers=optimal_num_centroids, nstart = 100, iter.max = 1000)$cluster
  
  #########################
  ##  Visualize Clusters ##
  #########################
  colors <- rainbow(length(unique(clusters)))
  cluster_colors <- unlist(lapply(clusters, function(x) { return(colors[x]) }))
  seq_clusters <- c(1:length(unique(clusters)))
  
  plot_uMAP <- function(c1, c2, data){
    par(mar=c(4.1, 4.1, 4.1, 4.1), xpd=TRUE)
    plot(data[,c1], data[,c2], xlab=paste0("Dim ",c1), ylab=paste0("Dim ",c2), col="white")
    legend("topright", inset=c(-0.2,0), col=colors[seq_clusters], legend=seq_clusters, pch=seq_clusters, title="Group")
    
    for(k in 1:length(unique(clusters))){
      datapoints <- which(clusters == k)
      points(data[datapoints,c1], data[datapoints,c2], col=colors[k], pch=k)
    }
  }
  
  pdf(paste0(output_path,"/uMAP.pdf"))
  par(family = 'serif', cex = 1.5)
  
  plot_uMAP(1, 2, new_data)
  
  if (d >= 3){
    plot_uMAP(1, 3, new_data)
    plot_uMAP(2, 3, new_data)
  }
  
  dev.off()
  
  ##########################
  ## Create Output Matrix ##
  ##########################
  
  # match duplicates with non duplicates 
  uMAP_peaks_corresponding_to_duplicate <- match(data.frame(t(data_normalized[duplicates,])), data.frame(t(data_removed_duplicates)))
  clusters_of_duplicates <- clusters[uMAP_peaks_corresponding_to_duplicate]
  
  # Generate cluster column
  cluster_col <- matrix(ncol = 1, nrow=nrow(data))
  cluster_col[which(duplicated(data_normalized) == FALSE), 1] <- clusters
  cluster_col[duplicates, 1] <- clusters_of_duplicates
  
  # Fill output matrix 
  outputmatrix <- cbind(outputmatrix, rep(-1, nrow(outputmatrix)))
  outputmatrix[match(duplicates, outputmatrix[,13]), ncol(outputmatrix)] <- clusters_of_duplicates
  outputmatrix[match(rownames(data_removed_duplicates), outputmatrix[,13]), ncol(outputmatrix)] <- clusters
  
  write.table(outputmatrix, file=paste0(data_path, "/final_tab_", filename, ".bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
  
  ########################################
  ## Generate Peak Profiles of Clusters ##
  ########################################
  
  end <- ncol(data_normalized)
  left <- which(apply(data_normalized[,1:3], 1, mean) < 0.8)
  right <- which(apply(data_normalized[,(end-2):end], 1, mean) < 0.8)
  peaks_to_translocate <- intersect(left, right)
  pre_testing1 <- t(apply(data_normalized, 1, center_peak, x_dist=dist))
  testing1 <- data_normalized
  testing1[peaks_to_translocate, ] <- pre_testing1[peaks_to_translocate,]
  
  testing2 <- t(apply(testing1, 1, smoothing, lambda=lam, dim=150))
  
  unique_clusters <- unique(clusters)
  for( i in 1:length(unique_clusters) ){
    peaks <- which(cluster_col[,1] == unique_clusters[i])
    num_peaks <- length(peaks)
    
    if ( num_peaks >= 20 ){
      num_peaks <- 20
    }
    
    if( on_off_smoothing ) {
      pdf(paste0(output_path, "/cluster_smoothed", unique_clusters[i],".pdf"))
      par(family = 'serif', cex = 1.5)
      for( j in 1:num_peaks ){
          plot(testing2[peaks[j],], ylab = "Normalized Read Count", xlab = "Nucleotide Position")
      }
      dev.off()
    }
    
    pdf(paste0(output_path, "/cluster_", unique_clusters[i],".pdf"))
    par(family = 'serif', cex = 1.5)
    for( j in 1:num_peaks ){
      plot(data_normalized[peaks[j],], ylab = "Normalized Read Count", xlab = "Nucleotide Position")
    }
    dev.off()
    
  }
}

#################
##  Parameters ##
#################

args = commandArgs(trailingOnly=TRUE)

if ( length(args) < 5 ) {
  print("[TEST] Not enough arguments.")
} else {
  umap_main(args[1], args[2], as.numeric(args[3]), as.numeric(args[4]), args[5])
}
