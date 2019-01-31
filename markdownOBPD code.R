
library("Rmpfr");

##### Simulations #####
SIS <- function(
  obs, #observed table
  n.sim = 1000, # Tables to generate per simulation
  dist = c("unif") #default proposal dist. is uniform but can also take hypergeometric
){
  r.sums <- rowSums(obs); #= c(10,62,13,11,39), #vector of row sums
  c.sums <- colSums(obs); #= c(65,25,45), #vector of column sums
  r.dim <- length(r.sums);#dimensions row
  c.dim <- length(c.sums);# dimensions column
  
  w <- as(1:n.sim,"mpfr"); # weights for simulated table I
  l.x <- as(1:n.sim,"mpfr"); #proportional dist. to pi(x). 
  
  for(z in 1:n.sim){
    X <- matrix(nrow = r.dim,ncol = c.dim) #Empty Matrix representing a table
    g <- matrix(nrow = (r.dim-1),ncol = (c.dim-1)) #simulated proposal dist.
    
    #saving the lower and upper values to compute q(T)
    lower <- matrix(nrow = r.dim-1,ncol = c.dim - 1);
    upper <- matrix(nrow = r.dim-1,ncol = c.dim - 1);
    for(j in 1:(c.dim-1)){
      for(i in 1:(r.dim-1)){
        #finding the lower and upper bound for each cell
        lower[i,j]<-max(0,c.sums[j]-sum(X[(1:i),j],na.rm = TRUE)-#total column sum
                          sum(
                            r.sums[(i+1):r.dim], #row sums of non simulated values
                            -sum(X[(i+1):r.dim,1:j],na.rm = TRUE),#row sum of already sim. values
                            na.rm = TRUE))
        upper[i,j]<- min(
          c.sums[j]- sum(X[1:i,j],na.rm = TRUE),#column sum - values already simulated in the col.
          sum(r.sums[i],-sum(X[i,1:j],na.rm = TRUE),na.rm = TRUE))#row sum - values sim. in row
        
        if(dist == c("unif")){#simulate from uniform distribution
          X[i,j] <- round(runif(1,min = lower[i,j],max =upper[i,j]))
          g[i,j] <- 1/(upper[i,j]-lower[i,j] + 1);
        }
        if(dist == c("hyper")){#simulate from hypergeometric distribution
          X[i,j] <- rhyper(1,m = upper[i,j], n = upper[i,j], k = lower[i,j]+upper[i,j]); 
          g[i,j] <- dhyper(X[i,j],m = upper[i,j], n = upper[i,j], k = lower[i,j]+upper[i,j]);
        }
      }
      X[r.dim,j] <- c.sums[j] - sum(X[1:(r.dim-1),j]);
    }
    
    X[,c.dim]<- (r.sums - apply(X,1,function(x)sum(x,na.rm = TRUE)));
    
    
    l.x[z]<-1/prod(gamma(as(X+1,"mpfr")))
    w[z]<- prod(g)/prod(gamma(as(X+1,"mpfr"))); #weight l.x/g.x 
    #g.x is the product of conditionals and l.x is proportional to multinomial
    
    
  }
  p.val <- sum((1/prod(factorial(obs))>=l.x)*w)/sum(w); #p-value
  #let h(x) = indicator function. w = weight.
  #use normal importance sampling to get an estimate of p-value
  
  list(weight = w, #vector with length (n.sim) of weights
       last.weight = g, #example of proposal dist of last table
       upper = upper, #example of the upper bounds of the last table
       lower = lower, #example of the lower bounds of the last table
       X = X, #example of last table
       #n.sim = n.sim, # number of simulated tables
       p.value = p.val #simulated p-value
  )
}


#list of files in dir.
contigency.table.dir <- "/Users/omachowda/Google Drive/StatCom 3/OBPD/Contingency tables/Income folder"
cont.tables.list <- list.files(contigency.table.dir)
#take only .csv files
cont.tables.list <- cont.tables.list[grep(cont.tables.list,pattern = ".csv")]
#get table names
table.name <- gsub(cont.tables.list,pattern = ".csv",replacement = "");
#assign table to name
for(i in 1:length(table.name)){
  tab <- read.csv(cont.tables.list[i]);
  tab <- tab[-c(nrow(tab)),-c(1,ncol(tab))]; #get rid of row and column sums and first 2 rows
  assign(table.name[i],tab)
}

SIS <- dget("/Users/omachowda/Google Drive/StatCom 3/OBPD/Robin/Sequential Importance Sampling Function.R") #import SIS function
#facilities to age output
SIS(get(table.name[1]),dist = "hyper")$upper
SIS(get(table.name[1]),dist = "hyper")$lower
SIS(get(table.name[1]),dist = "hyper")$X
SIS(get(table.name[1]),dist = "hyper")$p.value

#reference character items as objects
SIS.results <- lapply(sapply(table.name,get),function(x) SIS(x,dist = "hyper")$p.value)
SIS.results


##### M^2 Test ######

corr_calc <- function(df){
  columns = length(df)
  print(length(df))
  new_df = data.frame()
  v1= c()
  v2= c()
  z = 0
  for (i in 1:4){
    for (j in 1:columns){
      new_df <- rbind(new_df, c(i,j, df[i,j]))
      z = z + df[i,j]
    }
  }
  
  for (i in 1:(4*columns)){
    v1 = append(v1, c(rep(new_df[i,1],new_df[i,3])))
  }
  
  for (i in 1:(4*columns)){
    v2 = append(v2, c(rep(new_df[i,2],new_df[i,3])))
  }
  
  fit = cor(x = v1,y = v2)
  M = (z-1)*(fit)
  p = pchisq(M,12)

  return (p)
}

resulted = lapply(sapply(table.name,get),FUN= corr_calc)

#### Mosaic Plots #####

counter=1;
par(mfrow=c(1,2))
for( i in 1:12){
  
  mosaicplot(x = get(table.name[i]), shade = TRUE,color = TRUE, main= table.name[i])
  
}

