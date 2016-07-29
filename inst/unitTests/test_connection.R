test.flowCL.connection <- function()
{
    prefix.info                 <- flowCL:::flowCL_query_data_prefix.info()
    que.hasProperSynonym        <- flowCL:::flowCL_query_data_hasProperSynonym()
    test.res <- flowCL:::queryMarker ( marker = "CD8+", query.file = que.hasProperSynonym, prefix.info = prefix.info,
                                       endpoint="http://cell.ctde.net:8080/openrdf-sesame/repositories/CL")
    return( checkTrue(nrow(test.res) >= 1, "Connection Check") )
}
