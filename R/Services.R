#' service class
#'
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' @useDynLib MV4
#' @export
#' @import stringr XML misc3d
services<-function() {
  # ------------------------------------------------
  # rotateMatrix
  # ------------------------------------------------
  rotateMatrix<-function( m , rotations = 1) {
    if(rotations == 1 ) m<-t(m[nrow(m):1,])
    if(rotations == 2 ) m<-m[nrow(m):1,ncol(m):1]
    if(rotations == 3 ) m<-t(m)[ncol(m):1,]
    return( m )
  }
  # ------------------------------------------------
  # getPointPlaneDistance
  # ------------------------------------------------
  getPointPlaneDistance<-function(Punto,Piano) {
    return(abs(Piano[1]*Punto[1]+Piano[2]*Punto[2]+Piano[3]*Punto[3]+Piano[4])/sqrt(Piano[1]^2+Piano[2]^2+Piano[3]^2))
  }
  # ------------------------------------------------
  # get3DPosFromNxNy
  # ------------------------------------------------
  get3DPosFromNxNy<-function(Nx,Ny,oM) {
    return(xy<-t(oM%*%c(Nx,Ny,0,1)))
  }
  # ------------------------------------------------
  # getPlaneEquationBetween3Points
  # ------------------------------------------------
  getPlaneEquationBetween3Points<-function(Pa,Pb,Pc) {
    ac<-(Pb[2]-Pa[2])*(Pc[3]-Pa[3])-(Pc[2]-Pa[2])*(Pb[3]-Pa[3])
    bc<-(Pb[3]-Pa[3])*(Pc[1]-Pa[1])-(Pc[3]-Pa[3])*(Pb[1]-Pa[1])
    cc<-(Pb[1]-Pa[1])*(Pc[2]-Pa[2])-(Pc[1]-Pa[1])*(Pb[2]-Pa[2])
    dc<--(ac*Pa[1]+bc*Pa[2]+cc*Pa[3])
    return(c(ac,bc,cc,dc))
  }

  # ===============================================================
  # getXMLStructureFromDICOMFile
  # ===============================================================
  getXMLStructureFromDICOMFile<-function(fileName, folderCleanUp = TRUE) {
    # browser()
    # aggiungi estensione xml
    fileNameXML<-paste(fileName,".xml")
    fileNameXML<-str_replace_all(string = fileNameXML , pattern = " .xml",replacement = ".xml")
    # salva solo path senza nome oggetto DICOM
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    # browser() # rg-im
    # se file con estensione xml gia' esiste nella cartella non fare nulla altrimenti lo aggiunge
    if(!file.exists( fileNameXML ) | folderCleanUp==TRUE) {
      # browser()
      stringa1<-"dcm2xml";
      if ( Sys.info()["sysname"] == "Windows") {
        options(warn=-1)
        # -im
        fileNameXML <- str_replace_all(string = fileNameXML , pattern = ".dcm.xml",replacement =  ".xml")
        # -fm
        p.fileNameXML <- paste( c('"',fileNameXML,'"') ,collapse = '')
        p.fileName <- paste( c('"',fileName,'"') ,collapse = '')
        stringa2<-paste(" +M  ",p.fileName,p.fileNameXML,collapse='')
        system2(stringa1,stringa2,stdout=NULL)
        options(warn=0)        
      } else {
        # -im
        fileNameXML <- str_replace_all(string = fileNameXML , pattern = ".dcm.xml",replacement =  ".xml")
        # -fm        
        p.fileNameXML <- paste( c("'",fileNameXML,"'") ,collapse = '')
        p.fileName <- paste( c("'",fileName,"'") ,collapse = '')
        stringa2<-paste(" +M  ",p.fileName,p.fileNameXML,collapse='')
        options(warn=-1)
        system2(stringa1,stringa2,stdout=NULL)
        options(warn=0)        
      }
    }
    # Load the XML file: restituisce il file xml nella variabile doc
    doc = xmlInternalTreeParse(fileNameXML)
    return(doc);
  }

  # ===============================================================
  # getDICOMTag
  # ===============================================================
  getDICOMTag<-function(tag=tag, fileName="", folderCleanUp = TRUE) {
    obj.S<-services();
# browser()
    # exemption: you want an Image!
    if(tag == "7fe0,0010") stop("Not available Yet! (for tag 7fe0,0010)")

    # build the XML file and get the XML structure
    doc<-getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = folderCleanUp)

    # build the QUERY
    stringaQuery<-''
    if(tag=="0018,0031" | tag=="0018,1072" | tag=="0018,1074" | tag=="0018,1075") {
      stringaQuery<-paste(c('/file-format/data-set/sequence[@tag="0054,0016"]/item/element[@tag="',tag,'"]'),collapse='');
    }
    
    if(stringaQuery=='') stringaQuery<-paste(c('/file-format/data-set/element[@tag="',tag,'"]'),collapse='');

    # execute the QUERY
    valore<-xpathApply(doc,stringaQuery,xmlValue);

    if(length(valore)==2) logObj$handle( "error" , "a tag in DICOM file seems to be duplicated"  );
    if(length(valore)==0) return(NA);

    valore<-valore[[1]]
    return(valore);
  }
  # ===============================================================
  # setDICOMTag
  # ===============================================================
  setDICOMTag<-function(tag, value, fileName ) {
    obj.S<-services();
    if(!file.exists(fileName)) stop(" the file does not exist")
    stringa1<-"dcmodify";
    stringa2<-paste(c(" -nb -m  '",tag,"'='",value,"' ",fileName),collapse='')
    options(warn=-1)
    system2(stringa1,stringa2,stdout=NULL)
    options(warn=0)
  }
  # ===============================================================
  # splittaTAG
  # ===============================================================
  splittaTAG<-function(stringa) {
    return( as.numeric(strsplit(stringa,split = "\\\\")[[1]])   )
  }

  # ===============================================================
  # new.trilinearInterpolator
  # ===============================================================
  new.trilinearInterpolator<-function( voxelCube , pixelSpacing.new  ,pixelSpacing.old  ) {
    Nx.old<-dim(voxelCube)[1];	Ny.old<-dim(voxelCube)[2];	Nz.old<-dim(voxelCube)[3]
    xDim.old<-pixelSpacing.old[1];	yDim.old<-pixelSpacing.old[2];	zDim.old<-pixelSpacing.old[3]
    xDim.new<-pixelSpacing.new[1];	yDim.new<-pixelSpacing.new[2];	zDim.new<-pixelSpacing.new[3]

    fattoreDiScalaX<-pixelSpacing.old[1]/pixelSpacing.new[1];
    fattoreDiScalaY<-pixelSpacing.old[2]/pixelSpacing.new[2];
    fattoreDiScalaZ<-pixelSpacing.old[3]/pixelSpacing.new[3];

    Nx.new<-ceiling(Nx.old * fattoreDiScalaX)
    Ny.new<-ceiling(Ny.old * fattoreDiScalaY)
    Nz.new<-ceiling(Nz.old * fattoreDiScalaZ)

    result<-array(rep( 0 , Nx.new * Ny.new * Nz.new ))

    res<-.C("newnewtrilinearInterpolator",
            as.integer(Nx.old),as.integer(Ny.old),as.integer(Nz.old),
            as.integer(Nx.new),as.integer(Ny.new),as.integer(Nz.new),
            as.double(pixelSpacing.old[1]),as.double(pixelSpacing.old[2]),as.double(pixelSpacing.old[3]),
            as.double(pixelSpacing.new[1]),as.double(pixelSpacing.new[2]),as.double(pixelSpacing.new[3]),
            as.double(voxelCube),as.double(result) );
    result<-array( res[[14]] , dim=c(Nx.new,Ny.new,Nz.new) )
    return( result )
  }
  # ========================================================================================
  # cropCube: crop a voxel cube in order to limit its dimension to the needs
  # ========================================================================================
  cropCube<-function( bigCube ) {
    matPos<-which(bigCube!=0,arr.ind = T)

    min.x<-min(matPos[,1]);     max.x<-max(matPos[,1])
    min.y<-min(matPos[,2]);     max.y<-max(matPos[,2])
    min.z<-min(matPos[,3]);     max.z<-max(matPos[,3])
    newCube<-bigCube[ min.x:max.x, min.y:max.y , min.z:max.z]
    location<-list( "min.x"=min.x, "max.x"=max.x, "min.y"=min.y, "max.y"=max.y, "min.z"=min.z, "max.z"=max.z  )
    return( list ( "voxelCube"=newCube, "location"=location) )
  }
  # ========================================================================================
  # expandCroppedCube: expand a cropped cube
  # ========================================================================================
  expandCroppedCube<-function( croppedCube, toDim, def.val.for.expanded.space = NA) {
    min.x <- croppedCube$info$location$min.x
    min.y <- croppedCube$info$location$min.y
    min.z <- croppedCube$info$location$min.z
    max.x <- croppedCube$info$location$max.x
    max.y <- croppedCube$info$location$max.y
    max.z <- croppedCube$info$location$max.z

    littleCube <- croppedCube$voxelCube
    bigCube <- array( def.val.for.expanded.space, dim = toDim )
    bigCube[ min.x:max.x,  min.y:max.y, min.z:max.z ] <- littleCube

    return( bigCube )
  }
  # ========================================================================================
  # applyErosion.2D: crop a voxel cube in order to limit its dimension to the needs
  # ========================================================================================
  erosion.2D<-function(  imageSlice, margin.x=1, margin.y=1 ) {
    erodedVoxelCube<- imageSlice
    nX<-dim(erodedVoxelCube)[1];    nY<-dim(erodedVoxelCube)[2];
    mx<-margin.x; my<-margin.y;
    iterator<-0; # this is just to avoid infinite loops...
    
    # erode it!
    minValue<-min(erodedVoxelCube[which(!is.na(erodedVoxelCube),arr.ind = T)])-100;
    erodedVoxelCube[which(is.na(erodedVoxelCube),arr.ind = T)]<-minValue
    
    aa<-.C("erosion",as.double(erodedVoxelCube), as.integer(nX), as.integer(nY),
           as.integer(1),as.integer(margin.x),as.integer(margin.y),
           as.integer(0), as.integer(iterator), as.integer(minValue))
    
    erodedVoxelCube<-array(aa[[1]], dim=c(nX,nY))
    erodedVoxelCube[which(erodedVoxelCube==minValue,arr.ind = T)]<-NA
    return(erodedVoxelCube)
  }
  # ========================================================================================
  # regionGrowing
  # ========================================================================================
  regionGrowing<-function(  imageVC, c.pos, c.threshold ) {
    # browser()
    nX<-dim(imageVC)[1];    nY<-dim(imageVC)[2];  nZ<-dim(imageVC)[2];
    iterator<-0; # this is just to avoid infinite loops...
    
    # erode it!
    minValue<-min(imageVC[which(!is.na(imageVC),arr.ind = T)])-100;
    imageVC[which(is.na(imageVC),arr.ind = T)] <- minValue
    maskedCube <- imageVC * 0;
    maskedCube[c.pos] <- 1
    xPos <- c.pos[1]; yPos <- c.pos[2]; zPos <- c.pos[3];
    
    aa<-.C("regionGrowing",as.double(imageVC) )
    
    # aa<-.C("regionGrowing",as.double(imageVC), as.double(maskedCube), 
    #        as.integer(nX), as.integer(nY),as.integer(nZ),
    #        as.integer(xPos), as.integer(yPos),as.integer(zPos),
    #        as.integer(iterator), as.double(c.threshold))
    # browser()
    # erodedVoxelCube<-array(aa[[1]], dim=c(nX,nY))
    # erodedVoxelCube[which(erodedVoxelCube==minValue,arr.ind = T)]<-NA
    # return(erodedVoxelCube)
  }  
  # # ========================================================================================
  # # applyErosion.2D: crop a voxel cube in order to limit its dimension to the needs
  # # ========================================================================================
  # PointInPolygon.Axial.2D<-function( axialImg, xVertex, yVertex   ) {
  #   aa<-.C("bidimentionalFlatPointInPolygon",
  #          as.integer(dim(axialImg)[1]), as.integer(dim(axialImg)),
  #          as.integer(1),as.integer(margin.x),as.integer(margin.y),
  #          as.integer(0), as.integer(iterator), as.integer(minValue))
  #   
  #   erodedVoxelCube<-array(aa[[1]], dim=c(nX,nY))
  #   erodedVoxelCube[which(erodedVoxelCube==minValue,arr.ind = T)]<-NA
  #   return(erodedVoxelCube)
  # } 
  # int *nX, int *nY,
  # int *numVertex,
  # double *xVertex, 
  # double *yVertex,
  # int *PIPvector 
  # ========================================================================================
  # triangle2mesh
  # ========================================================================================  
  triangle2mesh <- function(x) {
    v <- list()
    n <- nrow(x$v1)
    nit <- 1:n
    v$vb <- t(cbind(rbind(x$v1,x$v2,x$v3),1))
    v$it <- rbind(nit,nit+n,nit+2*n)
    class(v) <- "mesh3d"
    return(v)
  }
  
  # ========================================================================================
  # StructureVolume
  # ========================================================================================
  StructureVolume<-function(mesh, measure.unit=c("cm3", "mm3")) {
    if (class(x = mesh)!="mesh3d") stop("mesh isn't a mesh3d object")
    measure.unit<-match.arg(arg = measure.unit)
    if (measure.unit=="mm3") mu<-1
    if (measure.unit=="cm3") mu<-1000
    Volume=0
    return(.C("MeshVolume", as.double(mesh$vb[1,]), as.double(mesh$vb[2,]), as.double(mesh$vb[3,]),
              as.integer(ncol(mesh$it)), as.integer(mesh$it[1,]-1), as.integer(mesh$it[2,]-1),
              as.integer(mesh$it[3,]-1), as.double(Volume))[[8]]/mu)
  }
  
  # ========================================================================================
  # StructureSurface
  # ========================================================================================
  StructureSurface<-function(mesh, measure.unit=c("cm2", "mm2")) {
    if (class(x = mesh)!="mesh3d") stop("mesh isn't a mesh3d object")
    measure.unit<-match.arg(arg = measure.unit)
    if (measure.unit=="mm2") mu<-1
    if (measure.unit=="cm2") mu<-100
    Surface=0
    return(.C("MeshSurface", as.double(mesh$vb[1,]), as.double(mesh$vb[2,]), as.double(mesh$vb[3,]),
              as.integer(ncol(mesh$it)), as.integer(mesh$it[1,]-1), as.integer(mesh$it[2,]-1),
              as.integer(mesh$it[3,]-1), as.double(Surface))[[8]]/mu)
  }  
  # ========================================================================================
  # get.HOT_IRON
  # ========================================================================================  
  get.HOT_IRON<-function( alpha ) {
    hotIron <- c(0,0,0,2,0,0,4,0,0,6,0,0,8,0,0,10,0,0,12,0,0,14,0,0,16,0,0,18,0,0,20,0,0,22,0,0,24,0,0,26,0,0,28,0,0,30,0,0,32,0,0,34,0,0,36,0,0,38,0,0,40,0,0,42,0,0,44,0,0,46,0,0,48,0,0,50,0,0,52,0,0,54,0,0,56,0,0,58,0,0,60,0,0,62,0,0,64,0,0,66,0,0,68,0,0,70,0,0,72,0,0,74,0,0,76,0,0,78,0,0,80,0,0,82,0,0,84,0,0,86,0,0,88,0,0,90,0,0,92,0,0,94,0,0,96,0,0,98,0,0,100,0,0,102,0,0,104,0,0,106,0,0,108,0,0,110,0,0,112,0,0,114,0,0,116,0,0,118,0,0,120,0,0,122,0,0,124,0,0,126,0,0,128,0,0,130,0,0,132,0,0,134,0,0,136,0,0,138,0,0,140,0,0,142,0,0,144,0,0,146,0,0,148,0,0,150,0,0,152,0,0,154,0,0,156,0,0,158,0,0,160,0,0,162,0,0,164,0,0,166,0,0,168,0,0,170,0,0,172,0,0,174,0,0,176,0,0,178,0,0,180,0,0,182,0,0,184,0,0,186,0,0,188,0,0,190,0,0,192,0,0,194,0,0,196,0,0,198,0,0,200,0,0,202,0,0,204,0,0,206,0,0,208,0,0,210,0,0,212,0,0,214,0,0,216,0,0,218,0,0,220,0,0,222,0,0,224,0,0,226,0,0,228,0,0,230,0,0,232,0,0,234,0,0,236,0,0,238,0,0,240,0,0,242,0,0,244,0,0,246,0,0,248,0,0,250,0,0,252,0,0,254,0,0,255,0,0,255,2,0,255,4,0,255,6,0,255,8,0,255,10,0,255,12,0,255,14,0,255,16,0,255,18,0,255,20,0,255,22,0,255,24,0,255,26,0,255,28,0,255,30,0,255,32,0,255,34,0,255,36,0,255,38,0,255,40,0,255,42,0,255,44,0,255,46,0,255,48,0,255,50,0,255,52,0,255,54,0,255,56,0,255,58,0,255,60,0,255,62,0,255,64,0,255,66,0,255,68,0,255,70,0,255,72,0,255,74,0,255,76,0,255,78,0,255,80,0,255,82,0,255,84,0,255,86,0,255,88,0,255,90,0,255,92,0,255,94,0,255,96,0,255,98,0,255,100,0,255,102,0,255,104,0,255,106,0,255,108,0,255,110,0,255,112,0,255,114,0,255,116,0,255,118,0,255,120,0,255,122,0,255,124,0,255,126,0,255,128,4,255,130,8,255,132,12,255,134,16,255,136,20,255,138,24,255,140,28,255,142,32,255,144,36,255,146,40,255,148,44,255,150,48,255,152,52,255,154,56,255,156,60,255,158,64,255,160,68,255,162,72,255,164,76,255,166,80,255,168,84,255,170,88,255,172,92,255,174,96,255,176,100,255,178,104,255,180,108,255,182,112,255,184,116,255,186,120,255,188,124,255,190,128,255,192,132,255,194,136,255,196,140,255,198,144,255,200,148,255,202,152,255,204,156,255,206,160,255,208,164,255,210,168,255,212,172,255,214,176,255,216,180,255,218,184,255,220,188,255,222,192,255,224,196,255,226,200,255,228,204,255,230,208,255,232,212,255,234,216,255,236,220,255,238,224,255,240,228,255,242,232,255,244,236,255,246,240,255,248,244,255,250,248,255,252,252,255,255,255)
    hotIron <- hotIron /max(hotIron)
    newHI <- c()
    for( i in seq(1,length(hotIron),by = 3) ) {
      newHI <- c(newHI,rgb(  hotIron[i], hotIron[i+1], hotIron[i+2],alpha = alpha))
    }
    return(newHI)
  }
  # ========================================================================================
  # getSOPClassUIDsTable
  # ========================================================================================  
  getSOPClassUIDsTable <- function() {
    SOPClassUIDs <- c(
      "0008,0016"="ComprehensiveSRStorage",
      "0008,0016"="CTImageStorage"
    )
    return( SOPClassUIDs )    
  } 
  # ========================================================================================
  # getSOPClassUIDsTable
  # ========================================================================================  
  anonymizeFolder <- function( folderPath ) {
    ooo <- geoLet()
    ooo$openDICOMFolder(pathToOpen = folderPath)
    
    arr.DICOMFiles <- ooo$getFilesInfo()[,"fileName"]
    
    cat("\n now anonymizing")
    for( fileName in arr.DICOMFiles) {
      cat("\n\t anonymizing: ",fileName)
      setDICOMTag(tag = "0010,0010",value = "xxx",fileName = fileName)
      setDICOMTag(tag = "0010,0020",value = "xxx",fileName = fileName)  
      setDICOMTag(tag = "0010,0030",value = "19211212",fileName = fileName)  
      setDICOMTag(tag = "0010,1010",value = "075Y",fileName = fileName)  
    }
    rm(ooo)
    gc()
  }
  binarize <- function(VC , from , to, reverse = FALSE ) {
    BVC <- VC
    BVC <- BVC * 0
    if(reverse == FALSE ) {
      BVC[ which( VC>=from & VC<=to) ] <- 1
      BVC[ which( VC<from | VC>to) ] <- 0
      
    } else {
      BVC[ which( VC>=from & VC<=to) ] <- 0
      BVC[ which( VC<from | VC>to) ] <- 1
    }
    return(BVC)
  }
  Skeletonize <- function( VC ,  binarization.window=c(300,1000) )  {
    bau <- VC
    bau <- binarize(bau, binarization.window[1], binarization.window[2]) * bau
    matrice.1 <- apply(  bau  , c(1,3), sd );  matrice.2 <- apply(  bau  , c(1,3), mean );   matrice <- (matrice.2 * 2 - matrice.1)
    cippa <-  binarize(matrice, -0.05, -0.03 ,reverse = TRUE) * matrice
    return(cippa)
  }

  # -------------------------------------------------------------------------------------------------
  return( list(
    "get3DPosFromNxNy"=get3DPosFromNxNy,
    "getPlaneEquationBetween3Points"=getPlaneEquationBetween3Points,
    "getPointPlaneDistance"=getPointPlaneDistance,
    "getXMLStructureFromDICOMFile" = getXMLStructureFromDICOMFile,
    "getDICOMTag"=getDICOMTag,
    "setDICOMTag"=setDICOMTag,
    "splittaTAG"=splittaTAG,
    "new.trilinearInterpolator"=new.trilinearInterpolator,
    "rotateMatrix"=rotateMatrix,
    "cropCube"=cropCube,
    "expandCroppedCube"=expandCroppedCube,
    "triangle2mesh"=triangle2mesh,
    "StructureVolume"=StructureVolume,
    "StructureSurface"=StructureSurface,
    "get.HOT_IRON"=get.HOT_IRON,
    "getSOPClassUIDsTable"=getSOPClassUIDsTable,
    "erosion.2D"=erosion.2D,
    "anonymizeFolder"=anonymizeFolder,
    "binarize"=binarize,
    "Skeletonize"=Skeletonize,
    "regionGrowing"=regionGrowing
  ))
}


