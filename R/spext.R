setMethod("coordinates",signature(obj="SpatialPolygons"),function(obj){
    coord_vals <- t(sapply(1:length(obj),function(i) obj@polygons[[i]]@Polygons[[1]]@labpt))
    colnames(coord_vals) <- colnames(obj@polygons[[1]]@Polygons[[1]]@coords)
    coord_vals
})
