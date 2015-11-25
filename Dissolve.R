#
# Version of gUnaryUnion from rgeos that returns a spatial
# polygons data frame with the polygon IDs as attributes in 
# the @data slot.
#
# Arguments:
#   SDF:      spatial polygons data frame for dissolve
#   Id:       either a scalar or a vector of new polygon IDs
#              - if a scalar (or omitted) all internal polygons removed
#              - vector must be same length as SDF
#
# Returns:
#             Spatial Polygons Data Frame - the @data slot has
#             two columns:
#               - SeqNum: sequence number 1...N
#               - rowID:  new polygon IDs (from input Id vector)
#
# Version 1.0
#
Dissolve <- function(SDF, Id=1) {

   if (length(SDF) == length(Id) | length(Id) == 1) {
      require(rgeos)

      if(length(Id) == 1) {
         temp <- gUnaryUnion(SDF,id=rep(Id,length(SDF)))
      } else {
         temp <- gUnaryUnion(SDF,id=Id)
      }

      ntemp <- length(temp)
      rowID  <- rep(NA,ntemp)
      for (i in 1:ntemp) rowID[i] <- temp@polygons[[i]]@ID
      temp.df <- data.frame(1:ntemp,rowID)
      colnames(temp.df) <- c("SeqNum","rowID")
      rownames(temp.df) <- rowID

      return(SpatialPolygonsDataFrame(temp,temp.df))
      } else {
      cat("Dissolve: SDF and Id length mismatch: aborting\n")
      NA
      }
}
