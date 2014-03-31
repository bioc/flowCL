test.flowCL.connection <- function()
{
    prefix.info                 <- flowCL:::flowCL_query_data_prefix.info()
    que.hasProperSynonym        <- flowCL:::flowCL_query_data_hasProperSynonym()
    test.res <- flowCL:::queryMarker ( marker = "CD8+", query.file = que.hasProperSynonym, prefix.info = prefix.info )
    return( checkTrue(nrow(test.res) >= 1, "Connection Check") )
}
