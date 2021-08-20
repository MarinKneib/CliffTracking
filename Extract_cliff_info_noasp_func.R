ExtractCliffInfo_noasp<-function(img,cliff_out,pond_out){
  # Function to extract necessary information for tracking of cliffs & ponds from shapefiles + raster
  #
  # Input data
  #   img: raster image used to map cliffs & ponds
  #   cliff_out: outline of cliffs from shapefile - in the shapefile, the cliffs are identified by a random number, their 'rank'
  #   pond_out: outline of ponds from shapefile - in the shapefile, the ponds are identified by a random number, their 'rank'
  #
  # Output data
  #   Cliff_df: data frame containing coordinates & aspect of each pixel of each cliff
  #   Pond_df: data frame containing coordinates of each pixel of each pond
  #   cliff_df_info: data frame containing area, center coordinates, mean aspect, aspect standard deviation, rank of each cliff 
  #                   along with radius and center coordinates of the circle fitted to it
  #   pond_df_info: data frame containing area, center coordinates, and rank of each pond in the shapefile
  #
  #Author: Marin Kneib
  #Work address: Swiss Federal Research Institute WSL
  #Email: marin.kneib@wsl.ch

  ###### Extract xy coordinates of each cliff & pond pixel
  cell_ids<-init(img,fun='cell') # Give a number to each of the image cells
  cliff_cell<-extract(cell_ids,cliff_out) # Extract cell number for each cliff
  cliff_xy<-list()
  for (i in 1:length(cliff_out)){
    cliff_xy[[i]]<-xyFromCell(img,cliff_cell[[i]]) # Retrieve coordinates of each cell for each cliff
  }
  
  pond_cell<-extract(cell_ids,pond_out) # Extract cell number for each pond
  pond_xy<-list()
  for (i in 1:length(pond_out)){
    pond_xy[[i]]<-xyFromCell(img,pond_cell[[i]]) # Retrieve coordinates of each cell for each cliff
  }
  
  ###### Create list with all the data for each cliff & pond
  Cliff_df<-list()
  ###### General info on cliff & ponds
  cliff_df_info<-list()
  
  ###### Calculate curvature & 'aspect' by fitting circle to the cliff's points and calculate angle by looking for each point at angle 
  # between vector from point to circle center and North 
  
  ### Function to optimize circle fit
  fitCircle <- function(xy,
                    a0=mean(xy[,1]),
                    b0=mean(xy[,2]),
                    r0=mean(sqrt((xy[,1]-a0)^2 + (xy[,2]-b0)^2)),
                    ...){
    SS <- function(abr){
      sum((abr[3] - sqrt((xy[,1]-abr[1])^2 + (xy[,2]-abr[2])^2))^2)
    }
    optim(c(a0,b0,r0), SS, lower=c(0,0,10), upper=c(Inf,Inf,5000), method="L-BFGS-B")
  }
  
  #### For each cliff calculate mean & std aspect
  for (i in 1:length(cliff_out)){
    # Calculate circle charateristics
    f = fitCircle(cliff_xy[[i]])
    a<-f$par[1]
    b<-f$par[2]
    R<-f$par[3]
    
    ### Calculate 'aspect' of each pixel
    aspect<-cliff_xy[[i]][,1]
    for (p in 1:length(cliff_xy[[i]][,1])){
      x<-cliff_xy[[i]][p,1]
      y<-cliff_xy[[i]][p,2]
      if ((y!=b) & (x!=a)){
        alpha<-acos((y-b)/(sqrt((y-b)^2+(x-a)^2)))*180/pi
        if(a>x){
          alpha<- -alpha
        }
        aspect[p]<-180+alpha
        if (aspect[p]>180){
          aspect[p]<-aspect[p]-360
        }
      }
      # If the pixel is at the same position than the center of the circle, aspect should be 0
      if ((y==b) & (x==a)){
        aspect[p]<-0
      }
    }
    m_aspect<-circ.mean(aspect*pi/180)*180/pi
    std_aspect<-sd.circular(aspect*pi/180)*180/pi
    
    #append aspect information to cliff pixel data frame
    Cliff_df[[i]]<-as.data.frame(cbind(cliff_xy[[i]],aspect))
    colnames(Cliff_df[[i]])<-c('x','y','Aspect')
    
    #Fill in general information on cliff
    area<-as.numeric(cliff_out@data$Area[i]) # Cliff area
    order<-cliff_out@data$rank[i]  # Cliff rank
    CenterX<-mean(cliff_xy[[i]][,1])  # centroid coordinates
    CenterY<-mean(cliff_xy[[i]][,2])
    mAspect<-m_aspect # aspect mean & std (in radians)
    stdAspect<-std_aspect
    cliff_df_info[[i]]<-as.data.frame(cbind(as.numeric(area),as.numeric(CenterX),as.numeric(CenterY),as.numeric(mAspect),as.numeric(stdAspect),as.numeric(order),as.numeric(R),as.numeric(a),as.numeric(b)))
    colnames(cliff_df_info[[i]])<-c('Area','CenterX','CenterY','meanAspect','stdAspect','Rank','Radius','CircleX','CircleY')
  }
  
  # For ponds keep x & y of each pixel
  Pond_df<-list()
  pond_df_info<-list()
  
  for (i in 1:length(pond_out)){
    Pond_df[[i]]<-as.data.frame(cbind(pond_xy[[i]]))
    colnames(Pond_df[[i]])<-c('x','y')
    area<-as.numeric(pond_out@data$Area[i]) # Pond area
    order<-pond_out@data$rank[i]  # Pond rank
    CenterX<-mean(pond_xy[[i]][,1])  # centroid coordinates
    CenterY<-mean(pond_xy[[i]][,2])
    pond_df_info[[i]]<-as.data.frame(cbind(area,CenterX,CenterY,order))
    colnames(pond_df_info[[i]])<-c('Area','CenterX','CenterY','Rank')
  }
  
  return(list(Cliff_df,cliff_df_info,Pond_df,pond_df_info))
}