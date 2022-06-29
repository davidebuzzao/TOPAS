TOPAS = function(network,
                 seeds,
                 expansion_steps = 2,
                 cores = 1){
  
  if (is.null(network) | ncol(network) < 2){
    stop('Network file is missing.', call. = FALSE)
  } 
  if (is.null(seeds)) {
    stop('Seeds are missing.', call. = FALSE)
  } 
  if (is.null(expansion_steps)) {
    expansion_steps = 2
  }
  if (is.null(cores)) {
    cores = parallel::detectCores()
  }
  if (class(cores) != 'numeric' | parallel::detectCores() < cores){
    stop('Please, introduce a numeric for cores.')
  }
 
  
  #########################################
  ########## STEP 1 - Full Seed Network
  message(paste0('STEP 1: Full Seed Network'))
  seeds_graph = igraph::simplify(igraph::graph_from_data_frame(network, directed = FALSE),
                                 remove.multiple = TRUE, remove.loops = TRUE)
  seeds = as.character(seeds)
  
  #########################################
  ########## STEP 2 - Largest Connected Module
  message(paste0('STEP 2: Largest Connected Module (using ',cores,ifelse(cores>1,' cores)', ' core)')))
  
  .sp_compute = function(source_v,seeds,graph){
    ## This function computes all shortest paths between 
    ## a seed and all other seeds and returns only those
    ## which are shorter than input "expansion_steps"
    dest_v = setdiff(seeds,source_v)
    sp = igraph::all_shortest_paths(
      graph,
      from = which(igraph::V(graph)$name == source_v),
      to = which(igraph::V(graph)$name %in% dest_v),
      weights = NULL)[['res']]
    
    return(unlist(lapply(sp,function(s){
      # print(s)
      if(length(names(s))>2 &
         length(names(s))<=(expansion_steps+2)){
        names(s[[2:(length(s)-1)]])
        }
      })
    )) 
  }
  
  while(igraph::count_components(seeds_graph)>1){
    ## Loop until 1 only connected component is retrieved
    best_cc = 1
    max_seeds = 0
    graph_lcc = igraph::components(seeds_graph)
    
    for (i in 1:igraph::count_components(seeds_graph)){
      seeds.cc = sum(seeds %in% names(graph_lcc$membership[graph_lcc$membership==i]))
      if (seeds.cc > max_seeds){
        best_cc = i
        max_seeds = seeds.cc
      }
    }
    
    graph_lcc = igraph::induced_subgraph(seeds_graph, which(igraph::V(seeds_graph)$name %in%
                                                      names(graph_lcc$membership[graph_lcc$membership==best_cc])))
    
    seeds = seeds[seeds %in% igraph::V(graph_lcc)$name]
    
    if (cores>1){
      ## Set up parallelization
      pb <- utils::txtProgressBar(min=0, max=length(seeds), style = 3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      cl <- parallel::makeCluster(cores)
      doSNOW::registerDoSNOW(cl)
      boot <- foreach::foreach(i = seeds, .options.snow = opts)
      system.time(connectors <- unique(
        unlist(
          foreach::`%dopar%`(boot, .sp_compute(i, seeds, graph_lcc)))
        ))
      parallel::stopCluster(cl)
      
    } else {
      connectors = unique(
        unlist(
          lapply(
            seeds,
            .sp_compute, seeds, graph_lcc)
        )
      )
    }
    
    seeds_graph = 
      igraph::induced_subgraph(
        graph_lcc, 
        unique(c(which(igraph::V(graph_lcc)$name %in% c(seeds,connectors)))))
    
  }
  
  if(igraph::ecount(seeds_graph)<1){
    ## If no edge is retrieved, the output is NULL
    message('\nNo module found!')
    return(NULL)
  }
  
  #########################################
  ########## STEP 3 - Pruned Module
  message(paste0('\nSTEP 3: Pruned Module'))
  
  ## Execute Random Walk With restart from seeds
  df = data.frame(v=igraph::V(seeds_graph)$name,val=0)
  df$val = ifelse(df$v %in% seeds, 1, 0)
  restart_list = as.matrix(df$val)
  rownames(restart_list) = df$v
  
  set.seed(1996)
  stationary_probability = 
    dnet::dRWR(seeds_graph,
               normalise = 'row',
               setSeeds = restart_list,
               restart = 0.75,
               verbose=F)
  
  df$stationary_probability = as.numeric(stationary_probability)
  
  set.seed(1996)
  df = df[sample(nrow(df)),] ## Randomize the alphabetical order first, it adds nothing if all probabilities vary
  df = df[order(df$stationary_probability),]; rownames(df) = NULL ## order by probability afterwards
  
  for (v in df[df$val==0,'v']){
    if (!v %in% igraph::V(seeds_graph)$name) {
      next
      }
    
    ## Temporary delete a candidate connector from the module.
    tmp_graph = igraph::delete_vertices(
      seeds_graph,
      which(igraph::V(seeds_graph)$name==v))
    
    ## If removing a vertex generated more than 1 connected components,
    ## then check if seed coverage decreases. 
    ## Else remove vertex permanently.
    if (igraph::count_components(tmp_graph)>1){
      best_cc = 1
      max_seeds = 0
      tmp_graph_lcc = igraph::components(tmp_graph)
      
      for (i in 1:igraph::count_components(tmp_graph)){
        seeds.cc = sum(seeds %in% names(tmp_graph_lcc$membership[tmp_graph_lcc$membership==i]))
        if (seeds.cc > max_seeds){
          best_cc = i
          max_seeds = seeds.cc
        }
      }
      
      if ((length(seeds) - max_seeds)==0){
        tmp_graph_lcc = 
          igraph::induced_subgraph(
            tmp_graph, 
            which(igraph::V(tmp_graph)$name %in% names(tmp_graph_lcc$membership[tmp_graph_lcc$membership==best_cc])))
        
        seeds_graph = rlang::duplicate(tmp_graph_lcc)
      } 
    } else {
      seeds_graph = rlang::duplicate(tmp_graph)
    }
  }
  return(as.data.frame(igraph::get.edgelist(seeds_graph)))
}
