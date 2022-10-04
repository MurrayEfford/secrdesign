maskEdgePolygon <- function (mask) {
    mat <- plotMaskEdge(mask, plt = FALSE)
    oneline <- function(x) sf::st_linestring(matrix(mat[,x], nrow=2,byrow=T))
    regions <- lapply(1:ncol(mat), oneline) 
    regions <- sf::st_multilinestring (regions)
    regions <- sf::st_polygonize(regions)
    sf::st_sfc(regions[[1]])
}
# region <- maskEdgePolygon(as.mask(alltraps))