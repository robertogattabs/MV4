#' class for loading and presenting DICOM data
#'
#' @description  Instantiate an object of the class \code{geoLet}.This represents just the classname,
#'               methods are exposed with the technique of 'closure'.
#'               In order to see manuals for the single mathods, consider the vignette or use the
#'               available for the following wrapping functions:
#'               \itemize{
#'               \item \code{openDICOMFolder( );} : to load a DICOM series into an geoLet object
#'               \item \code{getImageVoxelCube( );} : to get the ImageVoxelCube stored into a geoLet object
#'               \item \code{getPixelSpacing( );} : to get the pixelSpacing (x,y,z) of the main ImageVoxelCube stored into a geoLet object
#'               \item \code{getROIList( );} : to get the list of the ROI defined in a geoLet object
#'               \item \code{getTag( );} : to get a single DICOM-tag of a DICOM file loaded into a geoLet object
#'               \item \code{getROIVoxels( );} : to get the IMAGE Voxels geometrically located into a ROI, for a given geoLet object
#'               \item \code{extractDoseVoxels( );} : to get the DOSE Voxels geometrically located into a ROI, for a given geoLet object
#'               \item \code{calculateDVH( );} : to get the DVH calculated from a geoLet object
#'               }
#'               The original methods for the class geoLet can also be invocked using the same name without the previx 'GTL.', i.e.:
#'               misc3d rgl Rvcg oce rmarkdown moments
#' @export
#' @useDynLib MV4
#' @import stringr XML progress
#' @importFrom mgcv in.out
geoLet<-function( use.ROICache = FALSE ) {
  
  # global variables
  internalAttributes<-list()                             # Attributes
  cacheArea <- list()                                    # Cache
  use.cacheArea <- FALSE                                  # do you have to use the cache?
  logObj<-logHandler()                                   # log/error handler Object
  dataStorage<-list()                                    # memory data structure
  SOPClassUIDList<-c()
  global_tableROIPointList<-c()
  global_openedPath <- c()
  
  #=================================================================================
  # openDICOMFolder
  # Loads a Folder containing one or more DICOM Studies
  #=================================================================================
  # Open a folder and load the content
  openDICOMFolder<-function( pathToOpen ) {

    if(!dir.exists(pathToOpen) & !file.exists(pathToOpen)) logObj$handle( "error" , "The indicate Path does not exist"  );
    # ----------------------------------------------
    # get the dcm file type
    # ----------------------------------------------
    if( internalAttributes$verbose == TRUE ) cat("\n Dir scouting:")
    SOPClassUIDList<<-getFolderContent( pathToOpen );
    
    # ----------------------------------------------
    # Load CT/RMN Scans
    # ----------------------------------------------
    if( internalAttributes$verbose == TRUE ) cat("\n Image Loading:\n ")
    loadCTRMNRDScans( );
    
    # ----------------------------------------------
    # Load US
    # ----------------------------------------------
    if( internalAttributes$verbose == TRUE ) cat("\n Image Loading ( US and SC ):\n ")
    loadUS( );    
    loadSecondaryCapture();

        # ----------------------------------------------
    # Carica l'RTStruct, se presente
    # ----------------------------------------------
    if( internalAttributes$verbose == TRUE ) cat("\n RTStruct Loading: ")
    a <- loadRTStructFiles()
    if( internalAttributes$verbose == TRUE ) cat( a$quantity," structures loaded" )
    
    # ----------------------------------------------
    # Carica l'RTDose, se presente
    # ----------------------------------------------
    # if( internalAttributes$verbose == TRUE ) cat("\n RTDose Loading: ")
    # loadRTDoseFiles()
    # if( internalAttributes$verbose == TRUE ) cat( a$quantity," dose loaded" )

    # ----------------------------------------------
    # Carica i nifti, se presenti
    # ----------------------------------------------
    # if( internalAttributes$verbose == TRUE ) cat("\n nifti files Loading: ")
    # a <- loadNIFTIFileDescription()
    # if( internalAttributes$verbose == TRUE ) cat( a$quantity," structures loaded" )
    
    global_openedPath <<- pathToOpen
    
  }
  
  loadRTDoseFiles <- function() {
    
    qualiRTDose <- which( SOPClassUIDList[  , "kind"] == "RTDoseStorage" )
    if( length(qualiRTDose) > 0 ) {
      for( riga in qualiRTDose ) {
        fileName <- SOPClassUIDList[qualiRTDose,"fileName" ]
        righe <- as.numeric(getTag(tag = "0028,0010" , fileName = fileName ))
        colonne <- as.numeric(getTag(tag = "0028,0011" , fileName = fileName ))
        PixelSpacing <- getTag(tag = "0028,0030" , fileName = fileName )
        BitsAllocated <- getTag(tag = "0028,0100" , fileName = fileName )
        BitsStored <- getTag(tag = "0028,0101" , fileName = fileName )
        HighBit <- getTag(tag = "0028,0102" , fileName = fileName )
        PixelRepresentation <- getTag(tag = "0028,0103" , fileName = fileName )
        DoseUnits <- getTag(tag = "3004,0002" , fileName = fileName )
        DoseType <- getTag(tag = "3004,0004" , fileName = fileName )
        GridFrameOffsetVector <- getTag(tag = "3004,000c" , fileName = fileName )
        FrameOfReferenceUID <- getTag(tag = "0020,0052" , fileName = fileName )
        SamplesPerPixel <- getTag(tag = "0028,0002" , fileName = fileName )
        DoseGridScaling <- as.numeric(getTag(tag = "3004,000e" , fileName = fileName ))
        
        numSlices <- length(unlist(strsplit( GridFrameOffsetVector , "\\\\")))        
        
        if( SamplesPerPixel != "1") { stop(" SamplesPerPixel is expected to be '1' in the current version of moddicom") }
        if( DoseType != "PHYSICAL") { stop(" SamplesPerPixel is expected to be 'PHYSICAL' in the current version of moddicom") }
        if( DoseUnits != "GY") { stop(" SamplesPerPixel is expected to be 'GY' in the current version of moddicom") }
        if( PixelRepresentation != "0") { stop(" PixelRepresentation is expected to be '0' in the current version of moddicom") }
        if( HighBit != "31") { stop(" HighBit is expected to be '31' in the current version of moddicom") }
        if( BitsStored != "32") { stop(" BitsStored is expected to be '32' in the current version of moddicom") }
        if( BitsAllocated != "32") { stop(" HighBit is expected to be '32' in the current version of moddicom") }
        
        
        
        # rn <- readBin(con = fileName, what = "integer", size = 4, endian = "little",n = file.size(fileName))
        rn <- readBin(con = fileName, what = "integer", size = 1, endian = "little",n = file.size(fileName))
        rn <- rn[ length(rn):1 ]
        
        aaa <- unlist(lapply( seq(1, ( righe * colonne * numSlices * 4) , by = 4) , function( pos ) {
          rn[ (pos+3) ] + rn[ (pos+2) ] * 2^8 + rn[ (pos+1) ] * 2^12 + rn[ (pos+0) ] * 2^16  
        }))
        rn <- aaa
        
        # rn <- rn[ length(rn):1 ]
        # browser()
        # oppa <- rn[ (length(rn)-( righe * colonne *numSlices)):length(rn)  ]
        oppa <- rn[ 1:( righe * colonne * numSlices)   ]
        oppa <- oppa[ length(oppa):1 ]
        matRN <- array(0,c(righe,colonne,numSlices))
        
        ct<-1
        for( z in seq(1,numSlices)) {
          for(x in seq(1,righe)) {
            for(y in seq(1,colonne)) {
              matRN[x,colonne-y,z]<-oppa[ct]
              ct<-ct+1 
            }
          }
        } 
        browser()
        # matRN <- matRN * DoseGridScaling
        
        SOPInstanceUID <- as.character( SOPClassUIDList[riga,"SOPInstanceUID"] )
        ImageOrientationPatient <-  SOPClassUIDList[riga,"ImageOrientationPatient"] 
        
        newMatRN <- array( 0 , c(colonne,righe,dim(matRN)[3]))
        for(z in 1:(dim(matRN)[3]) ) {
          immagine <- t(matRN[(dim(matRN)[1]):1,,z])
          newMatRN[,,z] <- immagine
        }        
        
        obj.S <- services()
        doc <- obj.S$getXMLStructureFromDICOMFile( fileName = fileName, folderCleanUp = FALSE )
        TransferSyntaxUID <- xpathApply(doc,'//element[@tag="0002,0010" and @name="TransferSyntaxUID"]',xmlValue)[[1]]
        if( !(TransferSyntaxUID %in% c("1.2.840.10008.1.2.1","1.2.840.10008.1.2") )) {
          stop("TransferSyntaxUID not compatible with the current version")
        }
        browser()
        if( TransferSyntaxUID == "1.2.840.10008.1.2"){
          # stop("TransferSyntaxUID not compatible with the current version")
          cat("\n\t TransferSyntaxUID is not little Endian Explicit: conversion could take some minutes..")
          toBits <- function (x, nBits = 32){ tail(rev(as.numeric(intToBits(x))),nBits) }
          soglia.cache <- 5
          voxel.2.consider <- which(newMatRN!=0,arr.ind = T)
          # tabella <- table(voxel.2.consider)
          # tabella <- tabella[order(tabella,decreasing = T)]
          # arr.valori.indici <- names(which(tabella > soglia.cache))
          tmp <- lapply(1:nrow(voxel.2.consider),function(riga) {
            voxval <- newMatRN[ voxel.2.consider[riga,1], voxel.2.consider[riga,2], voxel.2.consider[riga,3] ] 
            bitty <- toBits( voxval , nBits = 32)
            newbitty <- c(bitty[17:32],bitty[1:16])
            newMatRN[ voxel.2.consider[riga,1], voxel.2.consider[riga,2], voxel.2.consider[riga,3] ] <<- sum(unlist(lapply(1:length(newbitty) ,  function(i){ (2*newbitty[i])^(length(newbitty)-i)}  )))-1
          })
          browser()
        }
        newMatRN <- newMatRN * DoseGridScaling
        
        
        dataStorage$doses[[ SOPInstanceUID ]] <<- list()
        dataStorage$doses[[ SOPInstanceUID ]]$dose <<- newMatRN
        dataStorage$doses[[ SOPInstanceUID ]]$rows <<- righe
        dataStorage$doses[[ SOPInstanceUID ]]$cols <<- colonne
        dataStorage$doses[[ SOPInstanceUID ]]$PixelSpacing <<- PixelSpacing
        dataStorage$doses[[ SOPInstanceUID ]]$DoseUnits <<- DoseUnits
        dataStorage$doses[[ SOPInstanceUID ]]$DoseType <<- DoseType
        dataStorage$doses[[ SOPInstanceUID ]]$GridFrameOffsetVector <<- GridFrameOffsetVector
        dataStorage$doses[[ SOPInstanceUID ]]$FrameOfReferenceUID <<- FrameOfReferenceUID
        dataStorage$doses[[ SOPInstanceUID ]]$FrameOfReferenceUID <<- FrameOfReferenceUID

      }
    }
  }
  
  #=================================================================================
  # loadNIFTIFiles
  # Loads the nifti files in the folder
  #=================================================================================
  loadNIFTIFileDescription<-function( ) {
    total.number <- 0
    for( riga in 1:nrow(SOPClassUIDList) ) {
      if( SOPClassUIDList[ riga , "kind"] == "nifti" ) {
        fileNameWithPath <- SOPClassUIDList[ riga , "fileName"]
        lastBS <- rev(unlist(str_locate_all(fileNameWithPath,"/")))[1]
        ROIName <- str_trim(str_replace_all(str_trim(str_sub(fileNameWithPath, lastBS+1)),".nii.gz",""))
        ROIName <- paste(c(ROIName,".nii"), collapse = '')
        
        # fileNameWithPath<-SOPClassUIDList[ riga , "fileName"]
        # aaa <- readNIfTI(fname = fileNameWithPath)
        # tmpVC <- slot(aaa,".Data")
        # if( slot(aaa,"reoriented") == FALSE ) logObj$sendLog(  "In the NIFTI file, 'reoriented' is set to FALSE" ,"ERR" );
        # if( slot(aaa,"scl_slope") != 1 ) logObj$sendLog(  "In the NIFTI file, 'slope' is no 1" ,"ERR" );
        # if( slot(aaa,"scl_inter") != 0 ) logObj$sendLog(  "In the NIFTI file, 'intercept' is no 1" ,"ERR" );
        # dim.x <- slot(aaa,"dim_")[2]
        # dim.y <- slot(aaa,"dim_")[3]
        # dim.z <- slot(aaa,"dim_")[4]
        dataStorage$structures[[ROIName]] <<- list()
        if( !("structures" %in% names(dataStorage$info)) ) dataStorage$info$structures <<- list()
        dataStorage$info$structures[[ROIName]]$type <<- "NIFTI"
        dataStorage$info$structures[[ROIName]]$fileName <<- fileNameWithPath
        dataStorage$info$structures[[ROIName]]$loaded <<- FALSE
        total.number <- total.number + 1
      }
    }
    return( list("quantity"=total.number))
  }
  #=================================================================================
  # loadRTStructFiles
  # Loads a DICOM RT Struct (one x folder)
  #=================================================================================
  loadRTStructFiles<-function( ) {
    
    imageSerie<-list();   listaPuntiROI<-list()
    explicitRTStructFileName <- internalAttributes$explicitRTStructFileName
    
    TMP<-list()
    if( is.na(explicitRTStructFileName)) {
      righe.RTStruct <- which(SOPClassUIDList[,"kind"] == "RTStructureSetStorage")
      for( riga in righe.RTStruct ) {
        fileName <- SOPClassUIDList[riga ,"fileName"]
        SOPInstanceUID <- SOPClassUIDList[riga ,"SOPInstanceUID"]
        TMP[[ SOPInstanceUID ]]<-getStructuresFromXML( fileName );
      }
    } else {
      TMP[[ explicitRTStructFileName ]]<-getStructuresFromXML( explicitRTStructFileName );
    }
    # browser()
    # now let me use some more easy to handle variable names
    matrice2<-c(); matrice3<-c(); FORUID.m<-NA;
    for(i in names(TMP)) {
      matrice2<-cbind(matrice2,TMP[[i]]$IDROINameAssociation)
      matrice3<-rbind(matrice3,TMP[[i]]$tableROIPointList)
      # Aggiungi le informazioni relative al FrameOfReferenceUID delle ROI caricate
      # ed il ReferencedROINumber!!!!! (xè non è il numero della ROI di moddicom, possono essere diversi)
      if( !("structures" %in% names(dataStorage$info)) ) dataStorage$info$structures<<-list();
      for( nomeROI in TMP[[i]]$IDROINameAssociation[2,] ) {
        dataStorage$info$structures[[nomeROI]]<<-list();
        dataStorage$info$structures[[nomeROI]]$FrameOfReferenceUID<<-TMP[[i]]$FORUID.m
        dataStorage$info$structures[[nomeROI]]$SeriesInstanceUID<<-TMP[[i]]$RTStructSeriesInstanceUID
        dataStorage$info$structures[[nomeROI]]$type<<-"DICOMRTStruct"
        dataStorage$info$structures[[nomeROI]]$associatedSlices<<-c()
      }
    }
    listaROI<-list()
    
    # for each ROI
    # browser()
    for(i in matrice2[2,]) {
      
      # get the points
      subMatrix<-matrice3[which(matrice3[,2]==i,arr.ind = TRUE),]
      # if some points exist
      quantiElementiTrovati <- -1
      
      if(is.list(subMatrix) & !is.array(subMatrix)) quantiElementiTrovati<-1
      if(length(subMatrix)==4 & !is.array(subMatrix) & is.matrix(subMatrix)==FALSE) subMatrix <- t(subMatrix)
      if(is.matrix(subMatrix) & is.array(subMatrix)) quantiElementiTrovati<-dim(subMatrix)[1]
      
      if(quantiElementiTrovati==-1) {
        logObj$sendLog( "Unexpected error in loading slices. No slices found.", "ERR"  );
      }
      
      if( quantiElementiTrovati >0 ) {
        listaROI[[i]]<-list()
        # add properly the points to the 'listaROI' structure
        for(contatore in seq(1,quantiElementiTrovati) ) {
          if( quantiElementiTrovati == 1) {
            ROIPointStringList<-subMatrix[[3]]
            SOPInstance<-subMatrix[[4]]
          }
          else {
            ROIPointStringList<-subMatrix[contatore,3][[1]]
            SOPInstance<-subMatrix[contatore,4][[1]]
          }
          listaCoords<-strsplit(ROIPointStringList,"\\\\");
          listaCoords<-as.numeric(listaCoords[[1]])
          # if a ROI already exists for the slice, well append it to the list
          if( !( SOPInstance  %in% names(listaROI[[i]])  ) )  listaROI[[i]][[   SOPInstance  ]]<-list()
          listaROI[[i]][[   SOPInstance  ]][[ length(listaROI[[i]][[   SOPInstance  ]])+1  ]]<-matrix(listaCoords,ncol=3,byrow=T)
          # Add the first one as last (close the loop)
          listaROI[[i]][[   SOPInstance  ]][[ length(listaROI[[i]][[   SOPInstance  ]]) ]]<-rbind(listaROI[[i]][[   SOPInstance  ]][[length(listaROI[[i]][[   SOPInstance  ]])]],listaROI[[i]][[   SOPInstance  ]][[length(listaROI[[i]][[   SOPInstance  ]])]][1,])
        }
      } else {
        listaROI[[i]]<-NA
      }
    }
    
    for( tmpSOPIUID in names(TMP)) {
      tableROIPointList <- TMP[[tmpSOPIUID]]$tableROIPointList
      ROINamesDaSistemare <- unique(tableROIPointList[,"ROIName"])
      for(ROIName in ROINamesDaSistemare ) {
        SOTTOMATRICE <- tableROIPointList[ which(tableROIPointList[,"ROIName"]==ROIName), ]
        dataStorage$info$structures[[ROIName]]$associatedSlices <<- SOTTOMATRICE
      }
    }
    
    dataStorage[["structures"]]<<-listaROI
    return( list("quantity"=length(listaROI)))
  }
  #####################################################################################
  # getStructuresFromXML: carica il file xml del RT struct
  #
  # INPUT:
  #   - fileName: il nome del file del RT struct
  # OUTPUT:
  #   -
  #################################################################################
  getStructuresFromXML<-function( fileName ) {
    
    obj.S<-services();
    massimo<-0
    folderCleanUp <- internalAttributes$attr_folderCleanUp
    arr.SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
    # browser()
    # Load the XML file if not in cache
    doc <- obj.S$getXMLStructureFromDICOMFile( fileName = fileName, folderCleanUp = folderCleanUp )
    # browser()
    # prima di tutto controlla che la FrameOfReferenceUID sia la stessa OVUNQUE e che punti
    # ad una serie di immagini ESISTENTE!
    # E' un chiodo ma .... ragionevole, almeno per ora
    # Estrae la seriesIstanceUID, la frame of reference UID e il Referenced Frame of Reference UID
    RTStructSeriesInstanceUID <- xpathApply(doc,'//element[@tag="0020,000e" and @name="SeriesInstanceUID"]',xmlValue)[[1]]
    FORUID.m <- xpathApply(doc,'//element[@tag="0020,0052" and @name="FrameOfReferenceUID"]',xmlValue)[[1]]
    FORUID.d <- xpathApply(doc,'//element[@tag="3006,0024" and @name="ReferencedFrameOfReferenceUID"]',xmlValue)
    # Analizza tutti i valori del Referenced frame of reference UID e controlla se tutti i valori coincidono con il
    # frame of reference UID
    for( FORUID.d_index in seq( 1, length(FORUID.d) ) ) {
      if( FORUID.d[[ FORUID.d_index ]] !=  FORUID.m ) {
        logObj$sendLog(  "FrameOfReferenceUID not aligned in RTStruct file" , "ERR" );
      }
    }
    
    if( length(unique(unlist(FORUID.d))) > 1 ) logObj$sendLog(  "more than 1 FrameOfReferenceUID in RTStruct file" , "ERR" );
    if( (unique(unlist(FORUID.d)) == FORUID.m ) == FALSE) ogObj$sendLog(  "FrameOfReferenceUID not aligned (?) in RTStruct file" , "ERR" );
    referencedFORUID <-  unique(unlist(FORUID.d))[1]
    
    # Guarda se ci sono delle immagini con quel FrameOfReferenceUID
    immagini.associabili <- which(SOPClassUIDList[ ,"FrameOfReferenceUID"] == FORUID.m & SOPClassUIDList[ ,"type"] =="IMG")
    if( length(immagini.associabili) == 0 ) {
      logObj$sendLog(  "the FrameOfReferenceUID of the RTStruct is not associated to an image" , "ERR"  );
    }
    
    # SEQUENCES: the one with the attribute  tag="3006,0020"  and name="StructureSetROISequence"
    # is the one with association NAME<->ID
    # Estrazione della parte di xml contenente il nome,numero delle varie ROI
    n2XML<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3006,0020" and @name="StructureSetROISequence"]/item')
    
    # SEQUENCES: now get the true coords
    # Estrazione delle coordinate di ciascuna ROI
    n3XML<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3006,0039" and @name="ROIContourSequence"]/item')
    
    # ROI Names
    matrice2<-c()
    # Per ciascuna ROI viene estratto il numero, il nome e organizzati in una matrice
    for(i in n2XML) {
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0022"]',xmlValue)[[1]]
      ROIName<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0026"]',xmlValue)[[1]]
      if( str_trim(ROIName) == "" ) ROIName <- ROINumber
      matrice2<-rbind(matrice2,c(ROINumber,ROIName))
    }
    
    matrice2<-t(matrice2)
    # ROI Point list
    massimo<-0
    matrice3<-c()
    
    # esegue una interazione per ogni ROI
    # Non fa altro che estrarre da 'i' la parte che contiene i vertici della ROI.
    for(i in n3XML) {
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0084"]',xmlValue)[[1]]
      ROIName<-matrice2[2,which(matrice2[1,]==ROINumber)]
      
      # la funzione getNodeSet() restituisce l'insieme delle coordinate dei vertici di una sola ROI
      # (sto intendendo per ROI l'insieme di tutti i contorni con lo stesso ROIName)
      listaPuntiDaRavanare<-getNodeSet(xmlDoc(i),'/item/sequence/item')
      numero.Punti.semiperimetro.massimo<-0
      
      # logObj$sendLog( paste( c( "\n Loading the ROI '".ROIName."' : ",length(listaPuntiDaRavanare)," polylines ") , collapse = '') )
      
      # Estrae solo un punto (fPoint.x, fPoint.y, fPoint.z) da listaPuntiDaRavanare in quanto viene assunto
      # che la ROI rispetto alla immagine sono tra loro paralleli (non sempre necessariamente vero)
      for(i2 in listaPuntiDaRavanare)   {
        
        # ReferencedSOPInstanceUID<-xpathApply(xmlDoc(i2),'//element[@tag="0008,1155"]',xmlValue)[[1]]
        # salva le coordinate delle ROI in una lista
        ROIPointList<-xpathApply(xmlDoc(i2),'/item/element[@tag="3006,0050"]',xmlValue)
        # organizza le coordinate in un vettore
        splittedROIPointList<-as.numeric(strsplit(ROIPointList[[1]],split = "\\\\")[[1]])
        # Separa in vettori diversi le coordinate dei punti x, y, z
        fPoint.x<-splittedROIPointList[1]
        fPoint.y<-splittedROIPointList[2]
        fPoint.z<-splittedROIPointList[3]
        
        # Calcola di quanti punti consta il "semiperimetro" della ROI
        # (questo serve per avere un'idea di dove prendere 3 punti abbstanza distanti per calcolare
        # il piano su cui giace)
        numero.Punti.semiperimetro <- as.integer((length(splittedROIPointList)/3)/2)
        
        # Se sei in presenza del numero di punti massimo, visto finora, prendi 3 punti
        if(numero.Punti.semiperimetro>numero.Punti.semiperimetro.massimo &
           numero.Punti.semiperimetro>10) {
          
          # prendi i tre punti sperabilmente più "distanti"
          # E' una stima, lo sa il cielo quali siano in realtà: dovrei
          # calcolare tutte le distanze reciproche! (ma anche no...)
          # Considero il primo punto della sequenza, quello ad un quarto e quello a metà
          p1 <- 1;    p2 <- numero.Punti.semiperimetro/2;      p3 <- numero.Punti.semiperimetro;
          # Costruisci la matrice di tutti i punti (così è più facile estrarre i 3 punti)
          matrice.punti <- matrix(splittedROIPointList,ncol=3,byrow = TRUE)
          p1 <- matrice.punti[p1,];
          p2 <- matrice.punti[p2,];
          p3 <- matrice.punti[p3,];
          # memorizza l'equazione del piano ed il numero di punti massimo della ROI
          numero.Punti.semiperimetro.massimo<-numero.Punti.semiperimetro
        }
        
        arr.ReferencedSOPInstanceUID<-c()
        if( length(arr.SeriesInstanceUID) == 0 )  {
          logObj$sendLog( "giveBackImageSeriesInstanceUID(); gave back nothing" , "ERR" );
        }
        
        # Cerca di assegnarlo ad una slice di immagine
        # Considerando un punto nella matrice della ROI calcola per ogni slice la distanza punto piano
        # in modo che se la distanza risulta minore di 0.2 assegna la ROI a quella determinata slice
        # cicla su tutte le immagini con quel FrameOfReferenceUID e cerca di capire a quale potresti associarla
        righe <- which(SOPClassUIDList[,"FrameOfReferenceUID"] == referencedFORUID & SOPClassUIDList[,"type"] == "IMG")
        for( riga in righe ) {
          tmp.SerInstUID <- SOPClassUIDList[riga,"SeriesInstanceUID"]
          slice.index <- as.character(SOPClassUIDList[riga,"SOPInstanceUID"])
          distanza <- obj.S$getPointPlaneDistance( c( fPoint.x, fPoint.y, fPoint.z ), dataStorage$info[[tmp.SerInstUID]][[slice.index]]$planeEquation)
          if( abs( distanza ) < internalAttributes$maxDistanceForImageROICoupling ) {
            arr.ReferencedSOPInstanceUID <- c( arr.ReferencedSOPInstanceUID, slice.index )
          }
        }
        
        # Assumento le slice di immagine fra loro parallele, vedi se è su un piano parallelo ad esse
        # se 'numero.Punti.semiperimetro.massimo'>0 significa che i tre punti della ROI sono stati estratti
        RSOPIUID.da.rimuovere <- c()
        for( RSOPIUID in arr.ReferencedSOPInstanceUID) {
          tmp.SerInstUID <- SOPClassUIDList[which(SOPClassUIDList[,"SOPInstanceUID"] == RSOPIUID),"SeriesInstanceUID"]
          
          if(numero.Punti.semiperimetro.massimo>0) {
            distanza1<-obj.S$getPointPlaneDistance(p1,dataStorage$info[[tmp.SerInstUID]][[RSOPIUID]]$planeEquation)
            distanza2<-obj.S$getPointPlaneDistance(p2,dataStorage$info[[tmp.SerInstUID]][[RSOPIUID]]$planeEquation)
            distanza3<-obj.S$getPointPlaneDistance(p3,dataStorage$info[[tmp.SerInstUID]][[RSOPIUID]]$planeEquation)
            # calcola il massimo gap fra i 3 punti.
            tot <- c(distanza1,distanza2,distanza3)
            tot <- max(tot)-min(tot)
            # Se è maggiore di un errore indicato ( .1 mm ) dichiara non paralleli i piani
            # (in realtà punti troppo vicini potrebbero fregarmi. Speriamo di no!)
            if(tot > .1 ) {
              # annovera fra le SOPInstanceUID da togliere dalla lista delle associabili
              logObj$sendLog ( paste( c("\n ROI ",ROIName," is not coplanar with images)" ),collapse = '')  );
              RSOPIUID.da.rimuovere <- c(RSOPIUID.da.rimuovere,RSOPIUID)
            }
          }
        }
        # Trattieni solo le SOPClassUID delle immagini che sono risultate paralele
        # (togli quelle non paralle dalla lista costruita precedentemente)
        arr.ReferencedSOPInstanceUID <- arr.ReferencedSOPInstanceUID[  which(!(arr.ReferencedSOPInstanceUID %in% RSOPIUID.da.rimuovere)) ]
        
        # Aggiorna la tabella delle associazioni ROI, immagini
        if( length(arr.ReferencedSOPInstanceUID) > 0 ) {
          for( tmp.RSOPUID in arr.ReferencedSOPInstanceUID ) {
            matrice3 <- rbind( matrice3, unlist(c( ROINumber, ROIName, ROIPointList, tmp.RSOPUID )) )
          }
        } else {
          logObj$sendLog (   paste(c("ROI ",ROIName,": the point (",fPoint.x,",",fPoint.y,",",fPoint.z,") has no image slice! "),collapse='')  );
        }
        # matrice3<-rbind(matrice3,c(ROINumber,ROIName,ROIPointList,ReferencedSOPInstanceUID))
      }
    }
    
    if( length(matrice3)>0) {
      if(is.matrix(matrice3)) {
        colnames(matrice3)<-c( "ROINumber", "ROIName", "ROIPointList", "ReferencedSOPInstanceUID"  )
      } else {
        names(matrice3)<-c( "ROINumber", "ROIName", "ROIPointList", "ReferencedSOPInstanceUID"  )
      }
    }
    
    colnames(matrice3) <- c("ROINumber", "ROIName", "ROIPointList", "ReferencedSOPInstanceUID")
    global_tableROIPointList <<- matrice3
    return(list("IDROINameAssociation"=matrice2,"tableROIPointList"=matrice3,"FORUID.m"=FORUID.m,"RTStructSeriesInstanceUID"=RTStructSeriesInstanceUID))
  }
  #=================================================================================
  # getFolderContent
  # Entra in una cartella ed estrai la lista dei file DICOM al suo interno, ricavandone
  # il tipo. Inoltre filtra gli oggetti DICOM solo alle SOPClassUID note.
  #=================================================================================
  getFolderContent <- function( pathToOpen ) {
    objS <- services()
    # browser()
    # if no path is given, use the set one
    if(!dir.exists(pathToOpen) & !file.exists(pathToOpen))  logObj$sendLog(  "The indicate Path does not exist" ,"ERR" );
    
    # salva in un array tutti i DICOM presenti nella cartella
    if( dir.exists(pathToOpen) == TRUE ) {
      DCMFilenameArray<-list.files(pathToOpen,internalAttributes$defaultExtension.dicom)
      NIFTIFilenameArray<-list.files(pathToOpen,internalAttributes$defaultExtension.nifti)
    } else {
      if( file.exists(pathToOpen) == TRUE ) {
        DCMFilenameArray <- substr( pathToOpen, max(unlist(str_locate_all( pathToOpen, "/")))+1, str_length(pathToOpen))   
        NIFTIFilenameArray <- c()
        pathToOpen <- substr( pathToOpen, 1,max(unlist(str_locate_all( pathToOpen, "/"))))
      }
    }
    
    # lista con la SOP Class UID di ciascun DICOM
    SOPClassUIDList<-list()
    nomiColonne <- c("fileName","tag","kind","type","IPP.x","IPP.y","IPP.z","FrameOfReferenceUID","ImageOrder","field2Order","p.x","p.y","p.z","SOPInstanceUID","ImageOrientationPatient","SeriesInstanceUID")
    MMatrix <- matrix(ncol=length(nomiColonne),nrow=0)
    colnames(MMatrix)<-nomiColonne
    # browser()
    ImagingPositionArray <- c()
    Iteration <- 0
    ImageOrder <- 1
    
    if( internalAttributes$verbose == TRUE )    pb <- progress_bar$new(total = length(DCMFilenameArray))
    
    for(i in 1:length(DCMFilenameArray) ) {
      
      if( internalAttributes$verbose == TRUE ) pb$tick()
      
      fileNameWithPath<-paste(pathToOpen,"/",DCMFilenameArray[i] , sep="")
      # browser()
      # Nel caso in cui sia DICOM
      # if( substr(fileNameWithPath,nchar(fileNameWithPath)-3,nchar(fileNameWithPath))=='.dcm' ) {
      if( substr( fileNameWithPath, nchar( fileNameWithPath ) - 3,nchar(fileNameWithPath))!='.xml' &
          substr( fileNameWithPath, nchar( fileNameWithPath ) - 3,nchar(fileNameWithPath))!='.raw' &
          substr( fileNameWithPath, nchar( fileNameWithPath ) - 3,nchar(fileNameWithPath))!='.txt' &
          substr( fileNameWithPath, nchar( fileNameWithPath ) - (str_length(internalAttributes$defaultExtension.nifti)-1),nchar(fileNameWithPath))!=internalAttributes$defaultExtension.nifti
      ) {
        
        # if( internalAttributes$verbose == TRUE ) cat(".")
        
        riga <- nrow(MMatrix) + 1
        MMatrix <- rbind(MMatrix, rep("",ncol(MMatrix))  )
        valore<-getDICOMTag( fileName = fileNameWithPath, tag = "0008,0016")
        MMatrix[riga, "fileName"] <- fileNameWithPath
        MMatrix[riga, "tag"] <- valore
        FrameOfReferenceUID<-getDICOMTag( fileName = fileNameWithPath, tag = "0020,0052")
        
        MMatrix[riga, "FrameOfReferenceUID"]<-FrameOfReferenceUID
        MMatrix[riga, "kind"]<-"Unknown"
        MMatrix[riga, "type"]<-"Unknown"
        
        if( valore == "1.2.840.10008.5.1.4.1.1.2" ) MMatrix[riga, "kind"]<-"CTImageStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.481.2" ) MMatrix[riga, "kind"]<-"RTDoseStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.481.3" ) MMatrix[riga, "kind"]<-"RTStructureSetStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.481.3" ) MMatrix[riga, "type"]<-"RTStruct"
        if( valore == "1.2.840.10008.5.1.4.1.1.481.5" ) MMatrix[riga, "kind"]<-"RTPlanStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.4" ) MMatrix[riga, "kind"]<-"MRImageStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.128" ) MMatrix[riga, "kind"]<-"PositronEmissionTomographyImageStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.2.1" ) MMatrix[riga, "kind"]<-"CTImageStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.6.1" ) MMatrix[riga, "kind"]<-"USImageStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.88.22" ) MMatrix[riga, "kind"]<-"EnhancedSR"
        if( valore == "1.2.840.10008.5.1.4.1.1.7" ) MMatrix[riga, "kind"]<-"SecondaryCaptureImageStorage"        
        
        
        if( MMatrix[riga, "kind"] == "RTDoseStorage") {
          MMatrix[riga, "ImageOrientationPatient"]<-getDICOMTag(fileName = fileNameWithPath , tag = "0020,0037")
          GridFrameOffsetVector <- getTag(tag = "3004,000c" , fileName = fileNameWithPath )
          PixelSpacing <- getTag(tag = "0028,0030" , fileName = fileNameWithPath )
          psp <- unlist(str_split(PixelSpacing,"\\\\"))
          GFOV <- as.numeric(unlist(str_split( GridFrameOffsetVector ,"\\\\")))
          GFOV <- unique(as.numeric(unlist(lapply( 2:length(GFOV) , function(i) {  GFOV[i]-GFOV[i-1] })))    )
          if( length(GFOV) > 1 ) {
            stop("GridFrameOffsetVector has more delta grid")
          }
          MMatrix[riga,"p.x"] <- psp[1]; MMatrix[riga,"p.y"] <- psp[2]; MMatrix[riga,"p.z"] <- GFOV
          # MMatrix[riga, "ImagePositionPatient"]<-getDICOMTag(fileName = fileNameWithPath , tag = "0020,0032")
        }
        if( MMatrix[riga, "kind"] == "USImageStorage") {
          MMatrix[riga, "type"]<-"IMG.US"  
        }
        if( MMatrix[riga, "kind"] == "SecondaryCaptureImageStorage") {
          MMatrix[riga, "type"]<-"IMG.SecondCapt"  
          pixelSpacing<-objS$splittaTAG(getDICOMTag(fileName = fileNameWithPath,tag = "0028,0030"))    
          MMatrix[riga, "p.x"]<-pixelSpacing[1]
          MMatrix[riga, "p.y"]<-pixelSpacing[2]          
        }        
        if( MMatrix[riga, "kind"] == "CTImageStorage" |
            MMatrix[riga, "kind"] == "MRImageStorage" |
            MMatrix[riga, "kind"] == "PositronEmissionTomographyImageStorage" )  {
          MMatrix[riga, "type"]<-"IMG"
          MMatrix[riga, "ImageOrder"]<-ImageOrder
          ImageOrder <- ImageOrder + 1
          
          # Prendi il Pixel Spacing e lo slice thickness
          pixelSpacing<-objS$splittaTAG(getDICOMTag(fileName = fileNameWithPath,tag = "0028,0030"))
          sliceThickness<-objS$splittaTAG(getDICOMTag(fileName = fileNameWithPath,tag = "0018,0050"))
          # MMatrix[riga, "ImagePositionPatient"]<-getDICOMTag(fileName = fileNameWithPath , tag = "0020,0032")
          MMatrix[riga, "p.x"]<-pixelSpacing[1]
          MMatrix[riga, "p.y"]<-pixelSpacing[2]
          MMatrix[riga, "p.z"]<-sliceThickness
          
          # ImageOrientation
          ImageOrientationPatient <- getDICOMTag(fileName = fileNameWithPath , tag = "0020,0037")
          MMatrix[riga, "ImageOrientationPatient"] <- ImageOrientationPatient
        }
        
        SOPInstanceUID <- getDICOMTag( tag = "0008,0018", fileName = fileNameWithPath)
        MMatrix[riga, "SOPInstanceUID"] <- SOPInstanceUID
        
        ImagePositionPatient <- getDICOMTag( fileName = fileNameWithPath, tag = "0020,0032")
        if( !is.na(ImagePositionPatient) ) {
          ImagingPosition <- objS$splittaTAG(ImagePositionPatient)
          MMatrix[riga, "IPP.x"] <- ImagingPosition[1]
          MMatrix[riga, "IPP.y"] <- ImagingPosition[2]
          MMatrix[riga, "IPP.z"] <- ImagingPosition[3]
        }
        SeriesInstanceUID<-getDICOMTag(tag = "0020,000e", fileName = fileNameWithPath )
        MMatrix[riga, "SeriesInstanceUID"] <- SeriesInstanceUID
      }
    }
    # Ora gestisti eventuali NIFTI
    for(i in 1:length(NIFTIFilenameArray) ) {
      fileNameWithPath<-paste(pathToOpen,"/",NIFTIFilenameArray[i] , sep="")
      # if( substr(fileNameWithPath,nchar(fileNameWithPath)-str_length("nii.gz"),nchar(fileNameWithPath))=='.nii.gz' ) {
      if(substr( fileNameWithPath, nchar( fileNameWithPath ) - (str_length(internalAttributes$defaultExtension.nifti)-1),nchar(fileNameWithPath))==internalAttributes$defaultExtension.nifti) {
        riga <- nrow(MMatrix)+1
        MMatrix <- rbind(MMatrix, rep("",ncol(MMatrix))  )
        MMatrix[riga, "fileName"] <- fileNameWithPath
        MMatrix[riga, "tag"] <- ""
        MMatrix[riga, "kind"] <- "nifti"
        MMatrix[riga, "type"] <- "nifti"
      }
    }
    
    # Ora devi ordinare le immagini in funzione delle loro coordinate!
    arr.SeriesInstanceUID <- unique(MMatrix[  which( MMatrix[,"type"]=="IMG" ), "SeriesInstanceUID" ])
    for( SIUID in arr.SeriesInstanceUID ) {
      righe <- which( MMatrix[,"type"]=="IMG" & MMatrix[,"SeriesInstanceUID"]==SIUID )
      sottomatrice <- MMatrix[  righe,  c("SOPInstanceUID","IPP.x","IPP.y","IPP.z") ]
      FOV.z <- diff(range(as.numeric(sottomatrice[,"IPP.z"])))
      FOV.y <- diff(range(as.numeric(sottomatrice[,"IPP.y"])))
      FOV.x <- diff(range(as.numeric(sottomatrice[,"IPP.x"])))
      
      # Ordinale per il gap maggiore (para-assiale/coronale/saggittale)
      winner <- order(c(FOV.x,FOV.y,FOV.z),decreasing = T)[1]
      if( winner == 1 ) field2Order <- "IPP.x"
      if( winner == 2 ) field2Order <- "IPP.y"
      if( winner == 3 ) field2Order <- "IPP.z"
      nuovo.ordine <- order(as.numeric(sottomatrice[,field2Order]),decreasing = F)
      tmp.SOPIUID <- sottomatrice[nuovo.ordine,"SOPInstanceUID"]
      for( i in 1:length(tmp.SOPIUID )) {
        MMatrix[ which( MMatrix[,"SOPInstanceUID"]==tmp.SOPIUID[i]  ) , "ImageOrder"] <- i
        MMatrix[ which( MMatrix[,"SOPInstanceUID"]==tmp.SOPIUID[i]  ) , "field2Order"] <- field2Order
      }
    }
    
    SOPClassUIDList <- MMatrix
    return(SOPClassUIDList);
  }
  
  loadUS <- function() {
    objS <- services()
    # Cicla sulle sole righe corrispondenti a delle immagini US
    righe.con.immagini <- which(SOPClassUIDList[,"type"]=="IMG.US") 
    for( riga in righe.con.immagini ) {
      fileName <- SOPClassUIDList[riga,"fileName"]
      SeriesInstanceUID<-SOPClassUIDList[riga, "SeriesInstanceUID"]
      
      if( internalAttributes$verbose == TRUE ) cat("\n\t ",fileName)
      SOPInstanceUID <- SOPClassUIDList[riga,"SOPInstanceUID"]
      
      if( !(SeriesInstanceUID %in% names(dataStorage[["info"]])) )  dataStorage[["info"]][[SeriesInstanceUID]] <<- list() 
      
      SingleSliceLoader <- loadCTRMNRDScans.SingleSlice( fileName = fileName )
      immagine <- SingleSliceLoader$immagine
      
      # now update the structure in memory
      if( length( dataStorage$img ) == 0 ) dataStorage$img <<- list()
      if( length( dataStorage$img[[SeriesInstanceUID]] ) == 0 ) dataStorage$img[[ SeriesInstanceUID ]] <<- list()
      dataStorage$img[[ SeriesInstanceUID ]][[ SOPInstanceUID ]] <<- immagine
    }
  }
  
  loadSecondaryCapture <- function() {
    objS <- services()
    # Cicla sulle sole righe corrispondenti a delle immagini US
    righe.con.immagini <- which(SOPClassUIDList[,"type"]=="IMG.SecondCapt") 
    for( riga in righe.con.immagini ) {
      fileName <- SOPClassUIDList[riga,"fileName"]
      SeriesInstanceUID<-SOPClassUIDList[riga, "SeriesInstanceUID"]
      
      if( internalAttributes$verbose == TRUE ) cat("\n\t ",fileName)
      SOPInstanceUID <- SOPClassUIDList[riga,"SOPInstanceUID"]
      
      if( !(SeriesInstanceUID %in% names(dataStorage[["info"]])) )  dataStorage[["info"]][[SeriesInstanceUID]] <<- list() 
      
      SingleSliceLoader <- loadCTRMNRDScans.SingleSlice( fileName = fileName )
      immagine <- SingleSliceLoader$immagine
      
      # now update the structure in memory
      if( length( dataStorage$img ) == 0 ) dataStorage$img <<- list()
      if( length( dataStorage$img[[SeriesInstanceUID]] ) == 0 ) dataStorage$img[[ SeriesInstanceUID ]] <<- list()
      dataStorage$img[[ SeriesInstanceUID ]][[ SOPInstanceUID ]] <<- immagine
    }
  }  
  #=================================================================================
  # loadCTRMRDNScans
  # Loads a DICOM CT/MR Scans
  #=================================================================================
  loadCTRMNRDScans<-function( ) {
    objS <- services()
    # Cicla sulle sole righe corrispondenti a delle immagini
    righe.con.immagini <- which(SOPClassUIDList[,"type"]=="IMG")
    
    for( riga in righe.con.immagini ) {
      
      fileName <- SOPClassUIDList[riga,"fileName"]
      SeriesInstanceUID<-SOPClassUIDList[riga, "SeriesInstanceUID"]
      
      if( internalAttributes$verbose == TRUE ) cat("\n\t ",fileName)
      
      FrameOfReferenceUID<-SOPClassUIDList[riga,"FrameOfReferenceUID"]
      ImageOrder <- SOPClassUIDList[riga,"ImageOrder"]
      SOPInstanceUID <- SOPClassUIDList[riga,"SOPInstanceUID"]
      kind.of.SOPClassUID <- SOPClassUIDList[riga, "kind"]
      IPP.x <- SOPClassUIDList[riga,"IPP.x"]
      IPP.y <- SOPClassUIDList[riga,"IPP.y"]
      IPP.z <- SOPClassUIDList[riga,"IPP.z"]
      ImageOrientationPatient <- SOPClassUIDList[riga,"ImageOrientationPatient"]
      ImagePositionPatient <- c(SOPClassUIDList[riga,"IPP.x"],SOPClassUIDList[riga,"IPP.y"],SOPClassUIDList[riga,"IPP.z"])
      pixelSpacing <- as.numeric(c(SOPClassUIDList[riga,"p.x"],SOPClassUIDList[riga,"p.y"]))
      
      # Se non era ancora stato settato, setta il FrameOfReferenceUID di riferimento
      if(is.na(internalAttributes$attr_mainFrameOfReferenceUID)) internalAttributes$attr_mainFrameOfReferenceUID <<- FrameOfReferenceUID
      # cat("\n\t",FrameOfReferenceUID)
      
      if(is.na(internalAttributes$attr_mainFrameOfReferenceUID)) stop("Error #39847")
      
      # Verifica che il FORUID sia compatibile con quello caricato (se no, dai errore)
      if( FrameOfReferenceUID != internalAttributes$attr_mainFrameOfReferenceUID ) stop("FORUID differente!")
      
      # if( SeriesInstanceUID %in% names(dataStorage[["info"]]) )  dataStorage[["info"]][[SeriesInstanceUID]] <- list()
      # if( SeriesInstanceUID %in% names(dataStorage[["info"]]) )
      if( !(SeriesInstanceUID %in% names(dataStorage[["info"]])) )  dataStorage[["info"]][[SeriesInstanceUID]] <<- list()
      
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]]<<-list()
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["PatientPosition"]]<<-c(SOPClassUIDList[riga,"IPP.x"],SOPClassUIDList[riga,"IPP.y"],SOPClassUIDList[riga,"IPP.z"])
      
      # Carica l'immagine
      SingleSliceLoader <- loadCTRMNRDScans.SingleSlice( fileName = fileName )
      immagine <- SingleSliceLoader$immagine
      
      # now update the structure in memory
      if( length( dataStorage$img ) == 0 ) dataStorage$img <<- list()
      if( length( dataStorage$img[[SeriesInstanceUID]] ) == 0 ) dataStorage$img[[ SeriesInstanceUID ]] <<- list()
      dataStorage$img[[ SeriesInstanceUID ]][[ SOPInstanceUID ]] <<- immagine
      
      # Costruisci la matrice per le trasformazioni affini
      iPP<-as.numeric(c(IPP.x,IPP.y,IPP.z))
      iOP<-objS$splittaTAG( ImageOrientationPatient )
      oM<-matrix(c(iOP[1],iOP[2],iOP[3],0,iOP[4],iOP[5],iOP[6],0,0,0,0,0,iPP[1],iPP[2],iPP[3],1),ncol=4);
      
      dataStorage[[ "info" ]][[ SeriesInstanceUID ]][[ SOPInstanceUID ]][[ "ImagePositionPatient" ]] <<- ImagePositionPatient
      dataStorage[[ "info" ]][[ SeriesInstanceUID ]][[ SOPInstanceUID ]][[ "ImageOrientationPatient" ]] <<- ImageOrientationPatient
      oM[1,1]<-oM[1,1]*pixelSpacing[1]
      oM[2,1]<-oM[2,1]*pixelSpacing[1]
      oM[3,1]<-oM[3,1]*pixelSpacing[1]
      oM[1,2]<-oM[1,2]*pixelSpacing[2]
      oM[2,2]<-oM[2,2]*pixelSpacing[2]
      oM[3,2]<-oM[3,2]*pixelSpacing[2]
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["orientationMatrix"]]<<-oM
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["pixelSpacing"]]<<-pixelSpacing
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["Rows"]]<<-getDICOMTag(tag = "0028,0010", fileName = fileName)
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["Columns"]]<<-getDICOMTag(tag = "0028,0011", fileName = fileName)
      if( kind.of.SOPClassUID == "MRImageStorage" ) {
        dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["RepetitionTime"]]<<-getDICOMTag( tag = "0018,0080", fileName = fileName)
        dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["EchoTime"]]<<-getDICOMTag(tag = "0018,0081", fileName = fileName)
      }
      # Aggiungi eventuali campi specifici per quel tipo di imaging (vero sopratutto per le PET/CT)
      for( campo in names(SingleSliceLoader$fields)) {
        dataStorage[[ "info" ]][[ SeriesInstanceUID ]][[ SOPInstanceUID ]][[ campo ]] <<- SingleSliceLoader$fields[[ campo ]]
      }
      
      if( kind.of.SOPClassUID == "PositronEmissionTomographyImageStorage" ) {
        # Fai qualche controllo di qualita', sulla reale possibilita' di importare le immagini PET
        rescale.type <-  SingleSliceLoader$fields$rescale.type
        UM <-  SingleSliceLoader$fields$UM
        CountsSource <- SingleSliceLoader$fields$CountsSource
        DecayCorrection <- SingleSliceLoader$fields$DecayCorrection
        deltaT <- SingleSliceLoader$fields$deltaT
        
        if(is.na(rescale.type)) rescale.type <- UM
        # browser()
        if( UM != rescale.type) {  logObj$sendLog(  "in PET image the rescale slope/intercept have different UM than the one used in the image (0054,1001) vs (0028,1054)" , "ERR"  )  }
        # if(CountsSource!='EMISSION' | DecayCorrection!='START') { logObj$sendLog(  c("\n ERROR: CountsSource!='EMISSION' or DecayCorrection!='START' ! This modality is not yet supported") , "ERR")  }
        if( DecayCorrection!='START') { logObj$sendLog(  c("\n ERROR: CountsSource!='EMISSION' or DecayCorrection!='START' ! This modality is not yet supported") , "ERR")  }
        if(is.na(deltaT) | is.null(deltaT) | deltaT==0) {  logObj$sendLog( "\n Error: deltaT between RadiopharmaceuticalStartTime and AcquisitionTime seems to be invalid" , "ERR" )  }
      }
      
      # three points to find out plane equation
      Pa<-c(oM[1,4],oM[2,4],oM[3,4])
      Pb<-objS$get3DPosFromNxNy(1000,0,oM)
      Pc<-objS$get3DPosFromNxNy(0,1000,oM)
      
      abcd<-objS$getPlaneEquationBetween3Points(Pa,Pb,Pc)
      piano<-matrix(abcd,nrow=1)
      colnames(piano)<-c("a","b","c","d")
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["planeEquation"]]<<-piano
      
      if(  kind.of.SOPClassUID == "RTDoseStorage" ) {
        stop("not yet implemented")
      }
    }
  }
  
  #=================================================================================
  # loadCTRMNRDScans.SingleSlice
  # Si occupa del caricamento di una singola slice. E' stato scorporato in quanto deve
  # essere potenzialmente evocato da più punti del programma ( gestione cache )
  #=================================================================================
  loadCTRMNRDScans.SingleSlice<-function( fileName ) {
    
    objServ<-services()
    fields <- list()
    
    # get the image data
    immagine<-getDICOMTag(tag = "7fe0,0010", fileName = fileName);
    SOPClassUID <- SOPClassUIDList[ SOPClassUIDList[,"fileName"]==fileName, "kind"]
    # browser()
    # apply rescaleSlope and rescaleIntercept, if needed
    rescale.intercept<-as.numeric(getDICOMTag(fileName = fileName,tag ="0028,1052" )); # rescale Intercept
    rescale.slope<-as.numeric(getDICOMTag(fileName = fileName,tag ="0028,1053" )); # rescale Slope
    rescale.type<-getDICOMTag(fileName = fileName,tag ="0028,1054" ); # rescale Type
    
    if(is.na(rescale.intercept)) rescale.intercept = 0;
    if(is.na(rescale.slope)) rescale.slope = 1;

    # se e' un US, non moltiplicare, esci subito'
    if( SOPClassUID == "USImageStorage" ) {
      if( rescale.intercept != 0 ) stop("\n Error: rescale.intercept not 0 for an US")
      if( rescale.slope != 1 ) stop("\n Error: rescale.slope not 1 for an US")      
      
      fields$rescale.intercept <- rescale.intercept
      fields$rescale.slope <- rescale.slope
      fields$rescale.type <- rescale.type
      
      return( list( "immagine" = immagine,
                    "fields" = fields ) )      
    }
    # se e' un US, non moltiplicare, esci subito'
    if( SOPClassUID == "SecondaryCaptureImageStorage" ) {
      # browser()
      if( rescale.intercept != 0 ) stop("\n Error: rescale.intercept not 0 for an US")
      if( rescale.slope != 1 ) stop("\n Error: rescale.slope not 1 for an US")      
      
      fields$rescale.intercept <- rescale.intercept
      fields$rescale.slope <- rescale.slope
      fields$rescale.type <- rescale.type
      
      return( list( "immagine" = immagine,
                    "fields" = fields ) )      
    }    
    
    immagine <- immagine * rescale.slope + rescale.intercept
    
    immagine <- objServ$rotateMatrix( immagine, rotations = 1 )
    
    if( SOPClassUID == "PositronEmissionTomographyImageStorage" ) {
      # browser()
      res <- calculate.SUVCoefficient.BW( fileName = fileName)
      SUVCoefficient.BW <- res$SUVCoefficient.BW
      immagine <- SUVCoefficient.BW * immagine
      
      # copia i campi 'fields' acquisisti dalla funzione per il calcolo del SUV e copiali nei campi
      # fields della lista che verra' restituita
      for( campo in names(res$fields)) {
        fields[[ campo ]] <- res$fields[[ campo ]]
      }
    }
    # browser()
    fields$rescale.intercept <- rescale.intercept
    fields$rescale.slope <- rescale.slope
    fields$rescale.type <- rescale.type
    
    return( list( "immagine" = immagine,
                  "fields" = fields ) )
  }
  #=================================================================================
  # calculate.SUVCoefficient.BW
  # Calcola il coefficiente moltiplicativo per il SUV
  #=================================================================================
  calculate.SUVCoefficient.BW<-function( fileName ) {
    # browser()
    fields <- list()
    AcquisitionTime<-getDICOMTag(fileName = fileName,tag ="0008,0032" );
    RadiopharmaceuticalStartTime<-getDICOMTag(fileName = fileName,tag ="0018,1072" );
    PatientWeight<-as.numeric(getDICOMTag(fileName = fileName,tag ="0010,1030" ));
    RadionuclideTotalDose<-as.numeric(getDICOMTag(fileName = fileName,tag ="0018,1074" )); # RadionuclideTotalDose
    RadionuclideHalfLife<-as.numeric(getDICOMTag(fileName = fileName,tag ="0018,1075" )); # RadionuclideHalfLife
    UM<-getDICOMTag(fileName = fileName, tag ="0054,1001" ); # UM of voxel Cube
    CountsSource<-getDICOMTag(fileName = fileName, tag ="0054,1002" ); # CountsSource
    DecayCorrection<-getDICOMTag(fileName = fileName, tag ="0054,1102" ); # DecayCorrection
    rescale.type<-getDICOMTag(fileName = fileName,tag ="0028,1054" ); # rescale Type)
    
    # deltaT<-as.numeric(difftime(as.POSIXct(AcquisitionTime, format = "%H:%M:%S"),as.POSIXct(RadiopharmaceuticalStartTime, format = "%H:%M:%S"),units = 'secs'))
    deltaT<-as.numeric(difftime(as.POSIXct(AcquisitionTime, format = "%H%M%S"),as.POSIXct(RadiopharmaceuticalStartTime, format = "%H%M%S"),units = 'secs'))
    
    # Una bella considerazione sul rescale.type ci potrebbe anche stare. Nel frattempo pongo il rescale a 1
    rescaleDueToUM<-1
    
    #SUVCoefficient.BW<-PatientWeight/( RadionuclideTotalDose * exp( -deltaT *log(2)/(RadionuclideHalfLife) ) )
    SUVCoefficient.BW<-PatientWeight/( RadionuclideTotalDose * 2^( -deltaT / RadionuclideHalfLife ) ) * 1000
    SUVCoefficient.BW<-SUVCoefficient.BW*rescaleDueToUM
    
    fields$AcquisitionTime <- AcquisitionTime
    fields$RadiopharmaceuticalStartTime <- RadiopharmaceuticalStartTime
    fields$PatientWeight <- PatientWeight
    fields$RadionuclideTotalDose <- RadionuclideTotalDose
    fields$RadionuclideHalfLife <- RadionuclideHalfLife
    fields$UM <- UM
    fields$CountsSource <- CountsSource
    fields$DecayCorrection <- DecayCorrection
    fields$rescale.type <- rescale.type
    fields$deltaT <- deltaT
    fields$SUVCoefficient.BW <- SUVCoefficient.BW
    
    return( list( "SUVCoefficient.BW" = SUVCoefficient.BW,
                  "fields" = fields ) );
  }
  #=================================================================================
  # getImageFromRAW
  # build a row data from a DICOM file stored on filesystem and load it
  # into memory (using DCMTK)
  #=================================================================================
  getImageFromRAW<-function(fileName) {
    
    # browser()
    
    objSV<-services()
    fileNameRAW<-paste(fileName,".0.raw")
    fileNameRAW<-str_replace_all(string = fileNameRAW , pattern = " .0.raw",replacement = ".0.raw")
    
    if(!file.exists(fileName)) logObj$sendLog( " the fileName is missing in geoLet::getImageFromRAW()", "ERR"  );
    # browser()
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    if(!file.exists( fileNameRAW )  | internalAttributes$attr_folderCleanUp==TRUE) {
      stringa1<-"dcmdump";
      if ( Sys.info()["sysname"] == "Windows") {
        fileNameFS<-chartr("\\","/",fileName);
        stringa2<-chartr("/","\\\\",stringa1)
      }
      else fileNameFS<-fileName;
      fileNameFS <- paste(c("'",fileNameFS,"'"),collapse = '')
      stringa2<-paste(" +W  ",pathToStore,fileNameFS,collapse='')
      options(warn=-1)
      stringone<-as.character(paste( c(stringa1," ",stringa2),collapse=''))
      # browser()
      # gestisci le system call in maniera diversa in funzione che sia WINDOWS o LINUX
      if ( Sys.info()["sysname"] == "Windows") {
        stringa2 <- chartr("/","\\",stringa2)
        stringa2 <- chartr("'",'"',stringa2)
        system2(stringa1,stringa2,stdout=NULL)
        # res<-.C("executeCMDLine",  as.character(stringone), as.integer(str_length(stringone))  )
      }
      else {
        system2(stringa1,stringa2,stdout=NULL)
      }
      options(warn=0)
    }
    rowsDICOM<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0010'))
    columnsDICOM<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0011'))
    bitsAllocated<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0100'))
    pixelRepresentation<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0100'))
    SOPClassUID <- SOPClassUIDList[ SOPClassUIDList[,"fileName"]==fileName, "kind"]
    # browser()
    
    if( SOPClassUID == "SecondaryCaptureImageStorage" ) {
      PhotometricInterpretation <- getDICOMTag(fileName = fileName,tag = '0028,0004')
      SamplesPerPixel <- getDICOMTag(fileName = fileName,tag = '0028,0002')        
      BitsAllocated <- as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0100'))
      BitsStored <- getDICOMTag(fileName = fileName,tag = '0028,0101')
      HighBit <- getDICOMTag(fileName = fileName,tag = '0028,0102')
      PixelRepresentation <- getDICOMTag(fileName = fileName,tag = '0028,0103')    
      
      # browser()
      if(! PhotometricInterpretation %in% c("MONOCHROME1","MONOCHROME2","RGB") ) {
        stop("\n PhotometricInterpretation not valid: #MaiUnaGioiaError")
      }
      if(! (BitsAllocated == 8  & BitsStored==8 & HighBit == 7 ) ) {
        stop("\n BitsAllocated, BitsStored and HighBit are expected to be 8,8,7: #MaiUnaGioiaError")
      }
      if(! (SamplesPerPixel == 1) ) {
        stop("\n SamplesPerPixel not valid, should be 1: #MaiUnaGioiaError")
      }
      if ( Sys.info()["sysname"] == "Windows") {
        # browser()
        fileNameRAWFS<-chartr("\\","/",fileNameRAW);
        fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
      }
      else fileNameRAWFS<-fileNameRAW;
      
      if(!file.exists(fileNameRAWFS)) logObj$sendLog( "problem in creating image binary file in geoLet::getImageFromRAW()", "ERR"  );
      
      rn <- readBin(con = fileNameRAWFS, what="integer", size=1, endian="little",n=rowsDICOM*columnsDICOM, signed = FALSE)
      us.R <- matrix(rn,ncol=columnsDICOM, byrow = TRUE)      
      us.R <- t(us.R)
      us.R <- us.R[,ncol(us.R):1]

      rn <- us.R
    }
    if( SOPClassUID == "USImageStorage" ) {
      UltrasoundColorDataPresent <- as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0014'))
      PhotometricInterpretation <- getDICOMTag(fileName = fileName,tag = '0028,0004')
      SamplesPerPixel <- getDICOMTag(fileName = fileName,tag = '0028,0002')  
      # browser()
      if(! PhotometricInterpretation %in% c("MONOCHROME1","MONOCHROME2","RGB") ) {
        stop("\n PhotometricInterpretation not valid: #MaiUnaGioiaError")
      }
      if(! (UltrasoundColorDataPresent == 1)  ) {
        stop("\n UltrasoundColorDataPresent not valid: #MaiUnaGioiaError")
      }
      if(! (SamplesPerPixel == 3) ) {
        stop("\n SamplesPerPixel not valid: #MaiUnaGioiaError")
      }
      if ( Sys.info()["sysname"] == "Windows") {
        # browser()
        fileNameRAWFS<-chartr("\\","/",fileNameRAW);
        fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
      }
      else fileNameRAWFS<-fileNameRAW;
      
      if(!file.exists(fileNameRAWFS)) logObj$sendLog( "problem in creating image binary file in geoLet::getImageFromRAW()", "ERR"  );
      
      # -im per Maria
      # if(file.info(fileNameRAWFS)$size == 0 & PhotometricInterpretation == "YBR_FULL_422") {
      #   fileNameRAWFS <- gsub("\\.0\\.raw","\\.1\\.raw",fileNameRAWFS)
      # }
      # -fm
# browser()
      rn <- readBin(con = fileNameRAWFS, what="integer", size=1, endian="little",n=rowsDICOM*columnsDICOM * 3, signed = FALSE)
      
      us.R <- rn[which(1:length(rn) %% 3 == 1)]
      us.G <- rn[which(1:length(rn) %% 3 == 2)]
      us.B <- rn[which(1:length(rn) %% 3 == 0)]
      
      us.R <- matrix(us.R,ncol=columnsDICOM, byrow = TRUE)
      us.G <- matrix(us.G,ncol=columnsDICOM, byrow = TRUE)
      us.B <- matrix(us.B,ncol=columnsDICOM, byrow = TRUE)
      
      us.R <- t(us.R)[,nrow(us.R):1]
      us.G <- t(us.G)[,nrow(us.G):1]
      us.B <- t(us.B)[,nrow(us.B):1]      
      
      rn <- list( "us.R"=us.R , "us.G"=us.G, "us.B"=us.B )
    }
    if( !( SOPClassUID %in% c("RTDoseStorage","USImageStorage","SecondaryCaptureImageStorage") ) ) {
      if( bitsAllocated!=16 & !(SOPClassUID %in% c("USImageStorage","SecondaryCaptureImageStorage")) )
        logObj$sendLog( "16bit pixel are allowed only for non-RTDoseStorage", "ERR"  );
      
      if( bitsAllocated!=8 & SOPClassUID == "USImageStorage" )
        logObj$sendLog( "US has been tested for 8 bit only", "ERR"  );      
      
      if ( Sys.info()["sysname"] == "Windows") {
        # browser()
        fileNameRAWFS<-chartr("\\","/",fileNameRAW);
        fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
      }
      else fileNameRAWFS<-fileNameRAW;
      
      if(!file.exists(fileNameRAWFS)) logObj$sendLog( "problem in creating image binary file in geoLet::getImageFromRAW()", "ERR"  );
      if( pixelRepresentation == 1 ) {
        stop("first time here")
        rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM)  
      } else {
        rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM, signed = FALSE)  
      }
      rn<-matrix(rn,ncol=columnsDICOM, byrow = TRUE)
    }
    if( SOPClassUID == "RTDoseStorage" ) {
      if(bitsAllocated==32) {
        if( SOPClassUID != "RTDoseStorage" )
          logObj$sendLog(  "32bit pixel are allowed only for RTDoseStorage" ,"ERR" );
        numberOfFrames<-as.numeric(getDICOMTag(fileName = fileName, tag = '0028,0008'))
        if ( Sys.info()["sysname"] == "Windows") {
          fileNameRAWFS<-chartr("\\","/",fileNameRAW);
          fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
        }
        else fileNameRAWFS<-fileNameRAW;
        
        if(!file.exists(fileNameRAWFS)) logObj$sendLog( "problem in creating image binary file in geoLet::getImageFromRAW()" ,"ERR" );
        
        rn<-readBin(con = fileNameRAWFS, what="integer", size=4, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames)
        # per ora va via come ciclo FOR, poi ci ragioniamo.... (per le performances)
        matRN<-array(0,c(rowsDICOM,columnsDICOM,numberOfFrames))
        ct<-1
        for( z in seq(1,numberOfFrames)) {
          for(x in seq(1,rowsDICOM)) {
            for(y in seq(1,columnsDICOM)) {
              matRN[x,columnsDICOM-y,z]<-rn[ct]
              ct<-ct+1
            }
          }
        }
        new_atRN<-array(0,c(columnsDICOM,rowsDICOM,numberOfFrames))
        for(ct in seq(1:dim(matRN)[3]  )) {
          new_atRN[,,ct]<-t(objSV$rotateMatrix( matRN[,,ct], rotations=2 ))
        }
        rn<-new_atRN
      } else  {
        numberOfFrames<-as.numeric(getDICOMTag(fileName = fileName, tag = '0028,0008'))
        
        if ( Sys.info()["sysname"] == "Windows") {
          fileNameRAWFS<-chartr("\\","/",fileNameRAW);
          fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS)
        }
        else fileNameRAWFS<-fileNameRAW;
        
        if(!file.exists(fileNameRAWFS)) logObj$sendLog( "problem in creating image binary file in geoLet::getImageFromRAW()", "ERR"  );
        
        rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames)
        matRN<-array(0,c(rowsDICOM,columnsDICOM,numberOfFrames))
        ct<-1
        for( z in seq(1,numberOfFrames)) {
          for(x in seq(1,rowsDICOM)) {
            for(y in seq(1,columnsDICOM)) {
              matRN[x,columnsDICOM-y,z]<-rn[ct]
              ct<-ct+1
            }
          }
        }
        new_atRN<-array(0,c(columnsDICOM,rowsDICOM,numberOfFrames))
        for(ct in seq(1:dim(matRN)[3]  )) {
          new_atRN[,,ct]<-t(objSV$rotateMatrix( matRN[,,ct], rotations=2 ))
        }
        rn<-new_atRN
      }
    }
    
    return(rn)
  }
  #=================================================================================
  # getROIList
  # restituisce la lista delle ROI
  #=================================================================================
  getROIList<-function() {
    if( length(dataStorage$structures) == 0 ) return( c() )
    mat2Ret<-matrix( c(seq(1,length(names(dataStorage$structures))),names(dataStorage$structures)),nrow=2 ,byrow=T )
    return(mat2Ret[2,])
  }
  #=================================================================================
  # getImageVoxelCube
  # give back the greyLevel voxel cube. If no ps.x/y/z are specified it gives back
  # the voxelCube of the original dimensions, otherwise it gives back the interpolated
  # voxelCube according to the wished pixelSpacing along x,y or z
  #=================================================================================
  getImageVoxelCube<-function( ps.x=NA, ps.y=NA, ps.z=NA , SeriesInstanceUID = NA , 
                               RGB.US.format = FALSE) {
    objS<-services();
    
    if( length(giveBackImageSeriesInstanceUID()) > 1 &
        is.na(SeriesInstanceUID) ) {
      logObj$sendLog(  "There are more than one Series, please specify which SeriesInstanceUID you want" ,"ERR" );
    }
    # browser()
    
    arr.img.type <- SOPClassUIDList[ which(SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID), "type"]
    # Se si tratta di una Secondary Capture
    if( "IMG.SecondCapt" %in% arr.img.type ) {
      if(all.equal( arr.img.type , rep("IMG.SecondCapt",length(arr.img.type))) == FALSE) {
        stop("Error.. some img from the same SeriesInstanceUID is not  IMG.SecondCapt")
      }    
      # browser()
      SOPInstanceUID <- names(dataStorage$img[[SeriesInstanceUID]])
      if( length(SOPInstanceUID) > 1 ) {
        stop("Error.. too many SOPInstanceUID for a given serie of IMG.SecondCapt")
      }  
      # browser()
      # -im
      mat.all <- dataStorage$img[[SeriesInstanceUID]][[SOPInstanceUID]]
      new.mat.all <- list()
      new.mat.all[[ SOPInstanceUID ]] <- mat.all
      # -fm
      return( new.mat.all )
      
    }
    
    # Se si tratta di US
    if( "IMG.US" %in% arr.img.type ) {
      if(all.equal( arr.img.type , rep("IMG.US",length(arr.img.type))) == FALSE) {
        stop("Error.. some img from the same SeriesInstanceUID is not an IMG.US")
      } 
      
      # mat.all <- dataStorage$img[[SeriesInstanceUID]]
      new.mat.all <- list()
      # browser()
      for( SOPInstanceUID in names(dataStorage$img[[SeriesInstanceUID]])) {
        mat.all <- dataStorage$img[[SeriesInstanceUID]][[SOPInstanceUID]]
        if( RGB.US.format == TRUE ) {
          mat.all$us.R <- t(mat.all$us.R[,ncol(mat.all$us.R):1])
          mat.all$us.G <- t(mat.all$us.G[,ncol(mat.all$us.G):1])
          mat.all$us.B <- t(mat.all$us.B[,ncol(mat.all$us.B):1])
          
          rgb.mat <- array( dim=c(dim(mat.all[[1]]),3) )
          rgb.mat[ , ,1 ] <- mat.all$us.R
          rgb.mat[ , ,2 ] <- mat.all$us.G
          rgb.mat[ , ,3 ] <- mat.all$us.B
          # rgb.mat <- rgb.mat / max(rgb.mat)        
        } else {
          rgb.mat <- array( dim=c(dim(mat.all[[1]]),3) )
          rgb.mat[ , ,1 ] <- mat.all$us.R
          rgb.mat[ , ,2 ] <- mat.all$us.G
          rgb.mat[ , ,3 ] <- mat.all$us.B        
        }
        new.mat.all[[ SOPInstanceUID ]] <- rgb.mat
      }
      # 
      # plot.new()
      # rasterImage( aaa[[1]]/max( aaa[[1]]) , 0,0,1,1, interpolate=FALSE )
      return( new.mat.all )
    }
        
    # Se no, prendi il cubone
    voxelCube <- createImageVoxelCube( SeriesInstanceUID = SeriesInstanceUID)
    
    # se non  server interpolare
    if(is.na(ps.x) && is.na(ps.y) && is.na(ps.z) ) return(voxelCube)
    
    # se invece serve interpolare: prendi i pixelSpacing lungo la X, la Y e la Z (slice thickness)
    oldPixelSpacing<-getPixelSpacing( SeriesInstanceUID = SeriesInstanceUID);
    
    if(is.na(ps.x))  ps.x <- oldPixelSpacing[1];
    if(is.na(ps.y))  ps.y <- oldPixelSpacing[2];
    if(is.na(ps.z))  ps.z <- oldPixelSpacing[3];
    
    voxelCube<-objS$new.trilinearInterpolator(
      voxelCube = voxelCube,
      pixelSpacing.new = c(ps.x,ps.y,ps.z),
      pixelSpacing.old = oldPixelSpacing )
    
    invisible(gc())
    
    return( voxelCube )
  }
  #=================================================================================
  # NAME: getROIVoxels
  # restituisce i voxel interni ad una data ROI
  #=================================================================================
  getROIVoxels<-function( Structure  , new.pixelSpacing=c(), SeriesInstanceUID = NA, croppedCube  = TRUE, 
                          onlyVoxelCube = FALSE, voxel.inclusion.threshold = 0.5,
                          giveBackOriginalImageToo = FALSE) {
    # browser()
    if( is.na(SeriesInstanceUID) ) {
      SeriesInstanceUID <- giveBackImageSeriesInstanceUID()  
    }
    
    if( !(Structure %in% names(dataStorage$info$structures)) ) {
      logObj$sendLog(  "Error: the mentioned ROI does not exist", "ERR"  );
      return( NA )
    }
    
    if( dataStorage$info$structures[[Structure]]$type == "DICOMRTStruct" ) {
      cat("\n=================================================================================")
      res <- getROIVoxels.DICOM(Structure = Structure , new.pixelSpacing=new.pixelSpacing,
                                SeriesInstanceUID = SeriesInstanceUID, croppedCube  = croppedCube,
                                onlyVoxelCube = onlyVoxelCube, 
                                giveBackOriginalImageToo = giveBackOriginalImageToo)
      cat("\n=================================================================================")
    }
    if( dataStorage$info$structures[[Structure]]$type == "NIFTI" ) {
      stop("\n Sorry, the NIFTI format has been pissed off :) \n")
      # res <- getROIVoxels.NIFTI(Structure = Structure , new.pixelSpacing=new.pixelSpacing,
      #                           SeriesInstanceUID = SeriesInstanceUID, croppedCube  = croppedCube,
      #                           onlyVoxelCube = onlyVoxelCube, voxel.inclusion.threshold = voxel.inclusion.threshold,
      #                           giveBackOriginalImageToo = giveBackOriginalImageToo)
    }
    
    return( res )
  }
  #=================================================================================
  # NAME: getROIFeatures
  # restituisce le features di una ROI
  #=================================================================================
  getROIFeatures<-function( Structure , statistical = TRUE, morphological = TRUE) {
    ROI <- getROIVoxels( Structure = Structure)
    arr.features <- c()
    if( statistical == TRUE ) {
      a <- statisticalFeatures(imgObj = ROI$voxelCube) 
      arr.features <- c( arr.features , a)
    }
    if( morphological == TRUE ) {
      a <- morphologicalFeatures(imgObj = ROI$voxelCube,px = getPixelSpacing()[1], py = getPixelSpacing()[2], pz = getPixelSpacing()[3])   
      arr.features <- c( arr.features , a)
    }    
    arr.features <- unlist(arr.features)
    return( arr.features )
  }  
  #=================================================================================
  # NAME: getROIVoxels.dicom
  # restituisce i voxel interni ad una data ROI - DICOM
  #=================================================================================
  getROIVoxels.DICOM<-function( Structure  , new.pixelSpacing=c(), SeriesInstanceUID = NA, croppedCube  = TRUE, 
                                onlyVoxelCube = FALSE, giveBackOriginalImageToo = FALSE) {
    objS<-services();
    # browser()
    if(!(Structure %in% getROIList())) logObj$sendLog(  paste(c( Structure," not present."  ),collapse = ''), "ERR"  );
    if( length(SeriesInstanceUID) > 1 ) logObj$sendLog(  paste( "Error, too many SeriesInstanceUIDs. No more than one is admitted."  ), "ERR"  );
    
    # try to find out which Series is the CT/MR serie
    if(is.na(SeriesInstanceUID)) {
      arr.SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
      if( length(arr.SeriesInstanceUID) == 0 ) {
        logObj$sendLog( "No image series seem to be available.. do you have any CT/MRI/PET in the indicated folder?" , "ERR" );
      }
      if( length(arr.SeriesInstanceUID) >1 ) {
        logObj$sendLog( "too many series are available: please pick up one" , "ERR" );
      }
      SeriesInstanceUID <- arr.SeriesInstanceUID[1]
    }
    # Questo va fatto solo se e' chiaro che non si tratta di un nifti'
    risultato.ROI <- getROIVoxelsFromCTRMN( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                            new.pixelSpacing = new.pixelSpacing, 
                                            giveBackOriginalImageToo = giveBackOriginalImageToo)  
    voxelCube <- risultato.ROI$ROIVoxelCube
    original.ROIVoxelCube <- risultato.ROI$original.ROIVoxelCube

    if(sum(!is.na(voxelCube)) == 0 ) return(list())
    
    info <- list()
    if( croppedCube == TRUE ) {
      croppedVC <- objS$cropCube( bigCube = voxelCube )
      voxelCube <- croppedVC$voxelCube
      info$cropped <- TRUE
      info$location <- croppedVC$location
    } else {
      info$cropped <- FALSE
      info$location$fe<-dim(voxelCube)[1]; info$location$se<-dim(voxelCube)[2]; info$location$te<-dim(voxelCube)[3]
      info$location$min.x <- 1; info$location$min.y <- 1; info$location$min.z <- 1
      info$location$max.x <- dim(voxelCube)[1]; info$location$max.y <- dim(voxelCube)[2]; info$location$max.z <- dim(voxelCube)[3];
    }
    
    # Se si vuole indietro SOLO il VOXELCube, prevediamo un ritorno semplificato
    if( onlyVoxelCube == TRUE ) return( voxelCube )
    
    invisible(gc())
    toReturn <- list( "voxelCube" = voxelCube, "info"= info  , "original.voxelCube"=original.ROIVoxelCube)
    class(toReturn)<-"ROIVoxel.struct"
    
    return( toReturn )
  }
  #=================================================================================
  # NAME: getROIVoxelsFromCTRMN
  # Estrae i voxel da scansioni CT,MR
  #=================================================================================
  getROIVoxelsFromCTRMN<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                   new.pixelSpacing=c() , giveBackOriginalImageToo = FALSE) {
    objService<-services()
    
    # se esiste gia' in cache, usa quanto gia' c'e'
    # browser()
    # if( "ROI" %in% names(cacheArea)) {
    #   if( Structure %in% names(cacheArea[["ROI"]]) & use.cacheArea == TRUE ) {
    #     if ( cacheArea[["ROI"]][[Structure]]$SeriesInstanceUID == SeriesInstanceUID & length(new.pixelSpacing)==0 ) {
    #       # return( cacheArea[["ROI"]][[Structure]]$voxelCube  )
    #       return( cacheArea[["ROI"]][[Structure]]  )
    #     }
    #   }
    # }
    
    # define some variables to make more clear the code
    numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
    numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);
    numberOfSlices<-length(dataStorage$img[[SeriesInstanceUID]]);
    
    # in case of interpolation
    old.ps <- getPixelSpacing( SeriesInstanceUID =  SeriesInstanceUID)
    if( length(new.pixelSpacing)>0 ) {
      if(new.pixelSpacing[1] == old.ps[1] & new.pixelSpacing[2]==old.ps[2]) new.pixelSpacing<-c()
    }
    if(  length(new.pixelSpacing)>0 ){
      new.pixelSpacing <- c(new.pixelSpacing,old.ps[3])
      ratio.ps.x <- old.ps[1] / new.pixelSpacing[1];   ratio.ps.y <- old.ps[2] / new.pixelSpacing[2]
    } else {
      new.numberOfRows <- numberOfRows;   new.numberOfColumns <- numberOfColumns
    }
    
    # initialize the image array with the right dimension
    # - im
    image.arr<-array( data = NA, dim = c(numberOfColumns, numberOfRows, numberOfSlices ) )
    # image.arr<-array( data = 0, dim = c(numberOfColumns, numberOfRows, numberOfSlices ) )
    # - fm
    
    # prendi le immagini cui e' associata la ROI (referenziata)'
    tabellaAssociazioni <- dataStorage$info$structures[[Structure]]$associatedSlices
    if( length(tabellaAssociazioni) == 4 ) {tmpAN <- names(tabellaAssociazioni) } else { tmpAn <- colnames(tabellaAssociazioni) }
    
    pb <- progress_bar$new(total = numberOfSlices)
    for( n in 1:numberOfSlices ) {
      tmpSOPIUID <- SOPClassUIDList[  which(SOPClassUIDList[,"ImageOrder"]==n  & SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID ),  "SOPInstanceUID" ]
      field2Order <- as.character(SOPClassUIDList[which( SOPClassUIDList[,"SOPInstanceUID"] == tmpSOPIUID),"field2Order"])
      
      # per il caso in cui la tabella abbia solo una riga... (castala a marix)
      if( length(tabellaAssociazioni) == 4 ) {  tabellaAssociazioni <- matrix(tabellaAssociazioni,ncol=4); colnames(tabellaAssociazioni)<-tmpAN    }
      
      # se questa slice risulta associata ad una ROI, allora calcola l'eventuale punto nel poligono
      if ( tmpSOPIUID %in% tabellaAssociazioni[,"ReferencedSOPInstanceUID"] ) {
        # browser()
        # Calcola la DOM di quella SLICE
        DOM <- dataStorage$info[[SeriesInstanceUID]][[tmpSOPIUID]]$orientationMatrix[c(1:3,5:7,13:15)]
        if(all(new.pixelSpacing==old.ps)==FALSE) {
          DOM[1:3] <- DOM[1:3]* (new.pixelSpacing[1]/old.ps[1])
          DOM[4:6] <- DOM[4:6]* (new.pixelSpacing[2]/old.ps[2])
        }
        DOM <- matrix(c(DOM[1:3],0,DOM[4:6],0,0,0,0,0,DOM[7:9],1),ncol=4)
        # expand grid dei punti della slice
        mtr.punti <- expand.grid( 1:numberOfColumns , 1:numberOfRows )
        # calcolo delle coordinate di ogni punto
        risultato <- t(apply( mtr.punti,MARGIN = 1, function(x) {  DOM %*% c( unlist(x),0,1)  }))
        risultato <- cbind( risultato, c( "Nx" = rep( NA, nrow(risultato) ) ) )
        risultato <- cbind( risultato, c( "Ny" = rep( NA, nrow(risultato) ) ) )
        risultato <- cbind( risultato, c( "Out" = rep( NA, nrow(risultato) ) ) )
        
        # cicla, perche' piu' poliline di una stessa ROI potrebbe essere associata alla slice
        polylineDaAnalizzare <- tabellaAssociazioni[tabellaAssociazioni[,"ReferencedSOPInstanceUID"]==tmpSOPIUID,"ROIPointList"]
        
        # if( length(polylineDaAnalizzare) > 1 ) browser()
        
        for( pol.num in 1:length(polylineDaAnalizzare)) {
          # browser()
          # -im
          tmp.risultato <- risultato
          # -fm
          # if( tmpSOPIUID == "1.3.12.2.1107.5.1.4.98879.30000019111423541678700003550" ) browser()
          
          ROI <- matrix(as.numeric(unlist(str_split(polylineDaAnalizzare[pol.num],"\\\\")[[1]])),ncol=3, byrow=T)
          xlim <- range(ROI[,1]); ylim <- range(ROI[,2]);  zlim <- range(ROI[,3])
          
          # estrai le righe valide (interne al bounding box di 'risultato')
          if(field2Order == "IPP.z" ) {
            righe.valide <- which((risultato[,1] >= xlim[1] & risultato[,1] <= xlim[2]) & (risultato[,2] >= ylim[1] & risultato[,2] <= ylim[2])  )
          }
          if(field2Order == "IPP.y" ) {
            righe.valide <- which((risultato[,1] >= xlim[1] & risultato[,1] <= xlim[2]) & (risultato[,3] >= zlim[1] & risultato[,3] <= zlim[2]) )            
          }          
          if(field2Order == "IPP.x" ) {
            righe.valide <- which((risultato[,2] >= ylim[1] & risultato[,2] <= ylim[2]) & (risultato[,3] >= zlim[1] & risultato[,3] <= zlim[2]) )
          }        
          
          # copia le righe valide all'interno della tabella
          risultato[righe.valide,5] <- mtr.punti[righe.valide,1]
          risultato[righe.valide,6] <- mtr.punti[righe.valide,2]
          
          if(field2Order == "IPP.z" ) {
            points2Test <- risultato[righe.valide,c(1,2)] # not yet validated
            ROI <- ROI[,c(1,2)]
          }
          if(field2Order == "IPP.y" ) {
            points2Test <- risultato[righe.valide,c(1,3)] # not yet validated
            ROI <- ROI[,c(1,3)]
          }
          if(field2Order == "IPP.x" ) {
            points2Test <- risultato[righe.valide,c(2,3)] # not yet validated
            ROI <- ROI[,c(2,3)]
          }
          
          # Ora fai il point-in-polygon
          # browser()
          punti.interni <- which(in.out(ROI,points2Test))
          
          # filtra risultato sui soli punti che ha senso scorrere per costruire la maschera di '1'
          # if( n >= 193) browser()
          # cat("\n N=",n)
          if( length(righe.valide) > 1 & length(punti.interni) > 1) {
            # -im
            # risultato <- risultato[righe.valide,][punti.interni,]
            tmp.risultato <- risultato[righe.valide,][punti.interni,]
            # -fm
            # posiziona gli '1' della maschera (sempre che il count dei punti sopra sia > 0)
            if( length(risultato) > 0 ) {
              # -im
              # tmp <- apply( risultato[,c(5,6)], MARGIN = 1, function(x){ image.arr[ x[1], x[2], n ] <<- 1} )
              tmp <- apply( tmp.risultato[,c(5,6)], MARGIN = 1, function(x){ image.arr[ x[1], x[2], n ] <<- 1} ) 
              # -fm
            }
          }
        }
      } 
      pb$tick()
    }
    # browser()
    # ROIVoxelCube <- VC[,,] * image.arr[,dim(image.arr)[2]:1,]
    # browser()
    ROIVoxelCube <- getImageVoxelCube( SeriesInstanceUID = SeriesInstanceUID)
    if(giveBackOriginalImageToo==TRUE) {
      original.ROIVoxelCube <- ROIVoxelCube
    } else {
      original.ROIVoxelCube <- NA
    }
    
    # -im 
    # ROIVoxelCube[which( image.arr[,dim(image.arr)[2]:1,] == 0, arr.ind = T)] <- NA
    ROIVoxelCube <- ROIVoxelCube * image.arr[,dim(image.arr)[2]:1,]
    # -fm 
    
    # aggiorna la cache
    da.restituire <- list( "ROIVoxelCube" = ROIVoxelCube,
                           "original.ROIVoxelCube" = original.ROIVoxelCube 
    )
    
    # browser()
    # if(  length(new.pixelSpacing)==0 & use.cacheArea == TRUE) {
    #   cacheArea[["ROI"]][[Structure]] <<- list()
    #   cacheArea[["ROI"]][[Structure]]$SeriesInstanceUID <<- SeriesInstanceUID
    #   if(length(new.pixelSpacing)>0)  cacheArea[["ROI"]][[Structure]]$pixelSpacing <<- new.pixelSpacing
    #   if(length(old.ps)>0)  cacheArea[["ROI"]][[Structure]]$pixelSpacing <<- old.ps
    #   cacheArea[["ROI"]][[Structure]]$voxelCube <<- ROIVoxelCube
    #   cacheArea[["ROI"]][[Structure]]$original.ROIVoxelCube <<- original.ROIVoxelCube
    # }
    
    return( da.restituire   )
  }
  #=================================================================================
  # getPixelSpacing
  #=================================================================================
  getPixelSpacing <- function( SeriesInstanceUID  = NA) {
    if(is.na(SeriesInstanceUID)) {
      SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
      if( length(SeriesInstanceUID) > 1) stop("Too many SeriesIntanceUID have been found: which one?")
    }
    # -im 
    # browser()
    # riga <- which(SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID & SOPClassUIDList[ ,"type"]=="IMG")[1]
    riga <- which(SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID & (SOPClassUIDList[ ,"type"] %in% c("IMG","IMG.US","IMG.SecondCapt")) )[1]
    # -fm
    p.x <- SOPClassUIDList[ riga, "p.x"]
    p.y <- SOPClassUIDList[ riga, "p.y"]
    p.z <- SOPClassUIDList[ riga, "p.z"]
    return( as.numeric(c( p.x , p.y, p.z ) ))
  }
  #=================================================================================
  # giveBackImageSeriesInstanceUID
  #=================================================================================
  giveBackImageSeriesInstanceUID<-function() {
    # arr.SeriesInstanceUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[ ,"type"]=="IMG")    ,"SeriesInstanceUID"])
    arr.SeriesInstanceUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[ ,"type"] %in% c("IMG","IMG.US","IMG.SecondCapt"))    ,"SeriesInstanceUID"])
    # browser()
    # if( length(arr.SeriesInstanceUID) > 1 ) stop("Error, two SeriesInstanceUID seems to be present")
    if( length(arr.SeriesInstanceUID) ==0 ) stop("Error, no images seems to be loaded")
    # return(arr.SeriesInstanceUID[1])
    return( arr.SeriesInstanceUID )
  }
  #=================================================================================
  # get3DPosFromNxNy
  #=================================================================================
  get3DPosFromNxNy<-function(  Nx, Ny, Nz, SeriesInstanceUID = NA , get.OM.back = FALSE) {
    
    # if not passed, get the series Instance UID of the images
    if(is.na(SeriesInstanceUID)) {
      SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
      if( length(SeriesInstanceUID) > 1) stop("Too many SeriesIntanceUID have been found: which one?")
    }    
    # browser()
    SOPInstanceUID <- as.character(SOPClassUIDList[ SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID & SOPClassUIDList[ ,"type"]=="IMG" & SOPClassUIDList[ ,"ImageOrder"]== Nz, "SOPInstanceUID"])
    oppa <- services()
    OM <- dataStorage$info[[ SeriesInstanceUID ]][[SOPInstanceUID]]$orientationMatrix
    res <- oppa$get3DPosFromNxNy(Nx,Ny,OM)[1:3]
    
    if( get.OM.back == TRUE) res <- list( "res"=res,"OM"=OM )
    return(res)
  }
  #=================================================================================
  # createImageVoxelCube
  # create the imageVoxelCube for the current obj and for the image stored
  #=================================================================================
  createImageVoxelCube<-function( SeriesInstanceUID = NA ) {
    # if not passed, get the series Instance UID of the images
    if(is.na(SeriesInstanceUID)) {
      SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
      if( length(SeriesInstanceUID) > 1) stop("Too many SeriesIntanceUID have been found: which one?")
    }
    
    Rows <- dim(dataStorage$img[[SeriesInstanceUID]][[1]])[1]
    Columns <- dim(dataStorage$img[[SeriesInstanceUID]][[1]])[2]
    Rows <- dim(dataStorage$img[[SeriesInstanceUID]][[1]])[2]
    Columns <- dim(dataStorage$img[[SeriesInstanceUID]][[1]])[1]
    
    Slices <- length(which(SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID))
    
    cubone<-array(data = 0,dim = c(Columns,Rows,Slices))
    
    for( i in seq( 1 , Slices ) )  {
      SOPInstanceUID <- SOPClassUIDList[ which(SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID & SOPClassUIDList[ ,"ImageOrder"]==i), "SOPInstanceUID"]
      cubone[,,i] <- dataStorage$img[[SeriesInstanceUID]][[SOPInstanceUID]]
    }
    return(cubone)
  }
  #=================================================================================
  # getDICOMTag
  # tag = which tag
  # fileName = nome del file (se presente)
  #=================================================================================
  getDICOMTag<-function( tag = tag, fileName ) {
    obj.S<-services();
    if(tag == "7fe0,0010") return( getImageFromRAW(fileName)  );
    
    if( internalAttributes$getTagXMLCacheList$fileName == fileName  ) {
      doc <- internalAttributes$getTagXMLCacheList$doc  
    } else {
      doc <- obj.S$getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = internalAttributes$attr_folderCleanUp)
      internalAttributes$getTagXMLCacheList$fileName  <<- fileName  
      internalAttributes$getTagXMLCacheList$doc <<- doc
    }
    
    a <- obj.S$getDICOMTag(tag = tag,fileName = fileName, 
                           folderCleanUp = internalAttributes$attr_folderCleanUp,
                           cached.Doc = doc)
    return(a)
  }
  #=================================================================================
  # getROIImageImageAssociations
  #=================================================================================
  getROIImageImageAssociations<-function( ROIName, details = FALSE  ) {
    ROISlicePositions <- c()
    if( length( ROIName ) != 1 ) stop("Please, specify the interested ROI")
    kkk <- global_tableROIPointList[  which( global_tableROIPointList[,"ROIName"] %in% ROIName ) ,"ReferencedSOPInstanceUID"]
    aaa <- table(SOPClassUIDList[  which( SOPClassUIDList[,"SOPInstanceUID"] %in% kkk ) , "kind"])
    
    if( details == TRUE ) {
      OIVC <- getROIVoxels(Structure = ROIName,croppedCube = F,onlyVoxelCube = T )  
      ROISlicePositions <- as.numeric((1:dim(OIVC)[3]) %in% unique(which(!is.na(OIVC),arr.ind = T)[,3] ))
    }
    
    return( list("recap" = aaa , "ROISlicePositions"= ROISlicePositions) )
  }
  #=================================================================================
  # getFilesInfo
  #=================================================================================
  getFilesInfo<-function( ) {
    return( SOPClassUIDList )
  }
  #=================================================================================
  # get.CT.SeriesInstanceUID
  #=================================================================================
  get.CT.SeriesInstanceUID<-function( ) {
    SIUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[,"kind"] == "CTImageStorage"), "SeriesInstanceUID" ])
    return( SIUID )
  }
  #=================================================================================
  # get.MR.SeriesInstanceUID
  #=================================================================================
  get.MRI.SeriesInstanceUID<-function( ) {
    SIUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[,"kind"] == "MRImageStorage"), "SeriesInstanceUID" ])
    return( SIUID )
  }
  #=================================================================================
  # get.MR.SeriesInstanceUID
  #=================================================================================
  get.PET.SeriesInstanceUID<-function( ) {
    SIUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[,"kind"] == "PositronEmissionTomographyImageStorage"), "SeriesInstanceUID" ])
    return( SIUID )
  }
  #=================================================================================
  # get.US.SeriesInstanceUID
  #=================================================================================
  get.US.SeriesInstanceUID<-function( ) {
    SIUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[,"kind"] == "USImageStorage"), "SeriesInstanceUID" ])
    return( SIUID )
  }  
  #=================================================================================
  # get.US.SeriesInstanceUID
  #=================================================================================
  get.SC.SeriesInstanceUID<-function( ) {
    SIUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[,"kind"] == "SecondaryCaptureImageStorage"), "SeriesInstanceUID" ])
    return( SIUID )
  }    
  
  #=================================================================================
  # class
  #=================================================================================
  class<-function( what="class") {
    if(what=="class") return("geoLet");
    if(what=="version") return("3");
    return("")
  }
  #=================================================================================
  # setAttribute
  #=================================================================================
  setAttribute<-function( folderCleanUp = NA, verbose = NA, RTStructFileName = NA, maxROIPlaneDistance = NA ) {
    if(!is.na(folderCleanUp)) internalAttributes$attr_folderCleanUp <<- folderCleanUp
    if(!is.na(verbose)) internalAttributes$verbose <<- verbose
    if(!is.na(RTStructFileName)) internalAttributes$explicitRTStructFileName <<- RTStructFileName
    if(!is.na(maxROIPlaneDistance))  internalAttributes$maxDistanceForImageROICoupling <<- maxROIPlaneDistance
  }
  # ===============================================================
  # getInterpolatedSlice
  # ===============================================================
  getInterpolatedSlice<-function( SIUID.pty , SIUID.ref , n.slice = NA , z.slice = NA , debug.mode = FALSE  ) {
    
    if( is.na(n.slice) & is.na(z.slice) ) stop(" n.slice OR z.slice should have a value")
    
    if( debug.mode == TRUE ) browser()
    CT.SIUID <- SIUID.pty
    PT.SIUID <- SIUID.ref
    
    CT.VC <- getImageVoxelCube(SeriesInstanceUID = CT.SIUID )
    PT.VC <- getImageVoxelCube(SeriesInstanceUID = PT.SIUID )    
    
    maxCT <- max( dim(CT.VC) ); maxPT <- max( dim(PT.VC) )
    mis.x.CT <- length( maxCT - dim(CT.VC)[1] )
    mis.y.CT <- length( maxCT - dim(CT.VC)[2] )
    mis.z.CT <- length( maxCT - dim(CT.VC)[3] )
    
    CT.x <- matrix(unlist(lapply( 1:dim(CT.VC)[1], function(x) {get3DPosFromNxNy(Nx = x,Ny = 1,Nz = 1,SeriesInstanceUID = CT.SIUID ) } )),ncol=3,byrow = T)[,1]
    PT.x <- matrix(unlist(lapply( 1:dim(PT.VC)[1], function(x) {get3DPosFromNxNy(Nx = x,Ny = 1,Nz = 1,SeriesInstanceUID = PT.SIUID ) } )),ncol=3,byrow = T)[,1]
    CT.y <- matrix(unlist(lapply( 1:dim(CT.VC)[2], function(x) {get3DPosFromNxNy(Nx = 1,Ny = x,Nz = 1,SeriesInstanceUID = CT.SIUID ) } )),ncol=3,byrow = T)[,2]
    PT.y <- matrix(unlist(lapply( 1:dim(PT.VC)[2], function(x) {get3DPosFromNxNy(Nx = 1,Ny = x,Nz = 1,SeriesInstanceUID = PT.SIUID ) } )),ncol=3,byrow = T)[,2]
    CT.z <- matrix(unlist(lapply( 1:dim(CT.VC)[3], function(x) {get3DPosFromNxNy(Nx = 1,Ny = 1,Nz = x,SeriesInstanceUID = CT.SIUID ) } )),ncol=3,byrow = T)[,3]
    PT.z <- matrix(unlist(lapply( 1:dim(PT.VC)[3], function(x) {get3DPosFromNxNy(Nx = 1,Ny = 1,Nz = x,SeriesInstanceUID = PT.SIUID ) } )),ncol=3,byrow = T)[,3]
    
    if(!is.na(n.slice))   posizione.PT.from <- between( CT.z[n.slice], PT.z , debug.mode = debug.mode)
    if(!is.na(z.slice))   posizione.PT.from <- between( z.slice, PT.z , debug.mode = debug.mode)
    posizione.PT.to <- posizione.PT.from + 1
    delta.z <- CT.z[2]-CT.z[1]
    if(!is.na(n.slice)) zRelativePosition <- CT.z[n.slice]-PT.z[posizione.PT.from]
    if(!is.na(z.slice)) zRelativePosition <- z.slice-PT.z[posizione.PT.from]
    # CT.z <- c(CT.z,seq(maxCT - length(CT.z)) * delta.z + CT.z[ length(CT.z) ]) 
    
    minVal <- min(PT.VC)-1000
    res <- CT.VC[,,1]
    res[1:dim(res)[1],1:dim(res)[2]] <- minVal
    
    if( !is.na(n.slice) & (CT.z[n.slice] %in% PT.z) ) {
      browser()
    } else {
      aaa <- .C("c_getInterpolatedSlice2D",
                as.integer(length(PT.x)), as.integer(length(PT.y)), as.double(delta.z), as.double(zRelativePosition),
                as.integer(length(CT.x)), as.integer(length(CT.y)),
                as.double(PT.x), as.double(PT.y),
                as.double(CT.x), as.double(CT.y),
                as.double(as.array(PT.VC[,,posizione.PT.from])),
                as.double(as.array(PT.VC[,,posizione.PT.to])),
                as.double(as.array(res)), as.double(minVal) );      
    }
    
    
    
    
    res <- matrix(aaa[13][[1]],nrow=dim(CT.VC)[1])
    return( res )
  }    
  #=================================================================================
  # getTag
  #=================================================================================    
  getTag<-function(tag=tag, fileName=NA, whichFile='first', SeriesInstanceUID = NA) {   
    
    if(is.na(SeriesInstanceUID) & is.na(fileName)) {
      SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
      if( length(SeriesInstanceUID) > 1) stop("Too many SeriesIntanceUID have been found: which one?")
    }
    
    
    if(!is.na(fileName)) return(getDICOMTag(tag = tag , fileName = fileName) )
    
    fileName <- SOPClassUIDList[ SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID,"fileName"][1]
    return(getDICOMTag(tag = tag , fileName = fileName) )
  }  
  # ===============================================================
  # between
  # ===============================================================
  between <- function(value, arr , debug.mode = FALSE) {
    if( debug.mode == TRUE ) browser()
    position <- unlist(lapply( 1:(length(arr - value)-1), function(x){  if( (arr - value)[x]*(arr - value)[x+1]<0  ) return(x) } ))
    return(position)
  }  
  # ===============================================================
  # getInterpolatedSlice
  # ===============================================================
  c_getInterpolatedSlice<-function( SIUID.pty , SIUID.ref , slice  ) {
    
    CT.SIUID <- SIUID.pty
    PT.SIUID <- SIUID.ref
    
    CT.VC <- getImageVoxelCube(SeriesInstanceUID = CT.SIUID )
    PT.VC <- getImageVoxelCube(SeriesInstanceUID = PT.SIUID )    
    
    CT.z <- matrix(unlist(lapply( 1:dim(CT.VC)[3], function(x) {get3DPosFromNxNy(Nx = 1,Ny = 1,Nz = x,SeriesInstanceUID = CT.SIUID ) } )),ncol=3,byrow = T)[,3]
    PT.z <- matrix(unlist(lapply( 1:dim(PT.VC)[3], function(x) {get3DPosFromNxNy(Nx = 1,Ny = 1,Nz = x,SeriesInstanceUID = PT.SIUID ) } )),ncol=3,byrow = T)[,3]
    
    zRiferimento <- CT.z[slice]
    slice.b <- which(unlist(lapply( 1:(length(PT.z)-1), function(x) {  sign((PT.z[x]-zRiferimento)*(PT.z[x+1]-zRiferimento)) } ))==-1)
    slice.t <- slice.b + 1
    
    SOPInstanceUID.b.PT <- SOPClassUIDList[which(SOPClassUIDList[,"SeriesInstanceUID"] == PT.SIUID & SOPClassUIDList[,"ImageOrder"] == slice.b),"SOPInstanceUID"]
    SOPInstanceUID.t.PT <- SOPClassUIDList[which(SOPClassUIDList[,"SeriesInstanceUID"] == PT.SIUID & SOPClassUIDList[,"ImageOrder"] == slice.t),"SOPInstanceUID"]
    SOPInstanceUID.CT <- SOPClassUIDList[which(SOPClassUIDList[,"SeriesInstanceUID"] == CT.SIUID & SOPClassUIDList[,"ImageOrder"] == slice.t),"SOPInstanceUID"]
    
    IOM.b.PT <- array(dataStorage$info[[PT.SIUID]][[SOPInstanceUID.b.PT]]$orientationMatrix)
    IOM.t.PT <- array(dataStorage$info[[PT.SIUID]][[SOPInstanceUID.t.PT]]$orientationMatrix)
    IOM.CT <- array(dataStorage$info[[CT.SIUID]][[SOPInstanceUID.CT]]$orientationMatrix)    
    nx.PT <- dim( PT.VC )[1]; ny.PT <- dim( PT.VC )[2];
    nx.CT <- dim( CT.VC )[1]; ny.CT <- dim( CT.VC )[2];
    slice.b.PT <- slice.b; slice.t.PT <- slice.t
    slice.CT <- slice
    image.b.PT <- PT.VC[,,slice.b]; image.t.PT <- PT.VC[,,slice.t];
    result <- CT.VC[,,slice] * 0 
    
    browser()
    
    res<-.C("c_getInterpolatedSlice",
            as.integer(nx.PT), as.integer(ny.PT),
            as.integer(slice.b.PT), as.integer(slice.t.PT),
            as.double(IOM.b.PT), as.double(IOM.t.PT),
            
            as.integer(nx.CT), as.integer(ny.CT),
            as.integer(slice.CT),
            as.double(IOM.CT),
            
            as.double(image.b.PT),as.double(image.t.PT),
            as.double(result) );
    return( res )
  }   
  getDataStorage  <- function() {
    return(dataStorage)
  }
  #=================================================================================
  # Constructor
  #=================================================================================
  constructor<-function( use.ROICache , verbose ) {
    
    # Attributes - set by user
    internalAttributes$maxDistanceForImageROICoupling<<-0.2
    internalAttributes$explicitRTStructFileName<<-NA
    internalAttributes$verbose<<-TRUE
    internalAttributes$attr_folderCleanUp<<-FALSE                      # force to re-dump DICOM files
    internalAttributes$attr_ROIVoxelMemoryCache<<-TRUE                 # force to cache ROI Voxel
    internalAttributes$attr_ROIVoxelMemoryCacheArray<<-list();
    internalAttributes$attr_arrayXMLCache<<-list();                    # array containint XML files
    internalAttributes$attr_dataChache<<-list();
    internalAttributes$attr_attributeList<<-''
    internalAttributes$attr_loadXMLInCache<<-FALSE
    internalAttributes$attr_loadRAWInCache<<-TRUE
    internalAttributes$attr_arrayXMLCache<<-list();
    internalAttributes$attr_arrayRAWCache<<-list();
    internalAttributes$attr_ROI.non.compl<<-c();
    internalAttributes$attr_neglectedROIs<<-c()
    internalAttributes$attr_withRTStruct<<-TRUE
    internalAttributes$attr_mainFrameOfReferenceUID<<-NA              # frameOfReference (geometry)
    internalAttributes$defaultExtension.dicom <<- ""
    internalAttributes$defaultExtension.nifti<<- ".nii.gz"
    internalAttributes$threshold.4.niftiROI<<- 0.4
    internalAttributes$getTagXMLCacheList$fileName  <<- ""
    internalAttributes$getTagXMLCacheList$doc <<- ""
    
    # Internal Structures and objs
    logObj <<- logHandler()                                   # log/error handler Object
    dataStorage <<- list()                                    # memory data structure
    dataStorage$info <<- list()
    SOPClassUIDList <<- list()
    global_tableROIPointList<<-c()
    cacheArea <<- list( "ROI" = list() )
    use.cacheArea <<- use.ROICache
  }
  constructor( use.ROICache = use.ROICache )
  return( list(
    "openDICOMFolder"=openDICOMFolder,
    "setAttribute"=setAttribute,
    "getImageVoxelCube"=getImageVoxelCube,
    "getPixelSpacing"=getPixelSpacing,
    "getROIList"=getROIList,
    "getFilesInfo"=getFilesInfo,
    "getROIVoxels"=getROIVoxels,
    "getROIFeatures"=getROIFeatures,
    "getROIImageImageAssociations"=getROIImageImageAssociations,
    "get.CT.SeriesInstanceUID"=get.CT.SeriesInstanceUID,
    "get.MRI.SeriesInstanceUID"=get.MRI.SeriesInstanceUID,
    "get.PET.SeriesInstanceUID"=get.PET.SeriesInstanceUID,
    "get.US.SeriesInstanceUID"=get.US.SeriesInstanceUID,  
    "get.SC.SeriesInstanceUID"=get.SC.SeriesInstanceUID,
    "get3DPosFromNxNy"=get3DPosFromNxNy,
    "getInterpolatedSlice"=getInterpolatedSlice,
    "getTag"=getTag,
    "getDataStorage"=getDataStorage
  ))
}
