#################################################################################
# query functions (flowCL: Semantic labelling of flow cytometric cell populations)
# Author: Justin Meskas (jmeskas@bccrc.ca)
# Date last modified: March 14, 2014
# Author: Radina Droumeva (radina.droumeva@gmail.com)
# Date last modified: July 19, 2013
#################################################################################
# Note:
#   The key regular expression is of the form: "CD19$|CD19[^0-9a-zA-Z][^ abes]".
#   This selects entries where the match is either of the exact marker passed 
#   followed by the end of line, or followed by anything other than a number
#   or letter or followed by a space which is followed by an 'a', 'b', 'e' or 's'. 
#   This will exclude alpha, beta, epsilon and single, however, it will exclude 
#   any other words starting with a, b, e or s (A potential bug).
#   Research into proper use of Regular Expression" yielded a \b command 
#   (i.e. \bbeta\b). Unfortunately this did not work. 
#   The following scenarios will no longer occur:
#   Case 1: Looking for CD19 returns CD193, etc.
#   Case 2: Looking for Ly6 returns Ly6g, etc.
#   Case 3: Looking for CD8 returns CD8 alpha chain
#################################################################################

#################################################################################
# function used by ontologyLabel
ontolExceptionsNeg1 <- function ( marker.list, q1 ) { 
    
    updatedException <- FALSE
    if ( marker.list[q1] == "CD8" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "IgD" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "CD3" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "HLA-DR" ) { updatedException <- TRUE }  
    return ( updatedException )
}
#################################################################################
# function used by ontologyLabel
ontolExceptionsPos1 <- function ( marker.list, q1 ) { 
  
    updatedException <- FALSE
    if ( marker.list[q1] == "CD8" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "IgD" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "CD3" )    { updatedException <- TRUE }
    if ( marker.list[q1] == "HLA-DR" ) { updatedException <- TRUE }  
    return ( updatedException )
}
#################################################################################
# function used by ontologyLabel
ontolExceptionsNeg2 <- function ( marker.list, q1, CompInfo ) { 
  
    if ( marker.list[q1] == "CD8" ) {
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"T cell receptor co-receptor CD8\"\n" ) }
        marker.list[q1] <-"T cell receptor co-receptor CD8"
    }
    if(marker.list[q1]=="IgD"){
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"IgD immunoglobulin complex\"\n" ) }
        marker.list[q1] <-"IgD immunoglobulin complex"
    }
    if(marker.list[q1]=="CD3"){
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"CD3 epsilon\"\n" ) }
        marker.list[q1] <-"CD3 epsilon"
    }
    if(marker.list[q1]=="HLA-DR"){
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"MHC class II histocompatibility antigen alpha chain HLA-DRA\"\n" ) }
        marker.list[q1] <-"MHC class II histocompatibility antigen alpha chain HLA-DRA"  
    }
    return(marker.list)
}
#################################################################################
# function used by ontologyLabel
ontolExceptionsPos2 <- function(marker.list, q1, CompInfo){ 
  
    if(marker.list[q1]=="CD8"){
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"T cell receptor co-receptor CD8\"\n" ) }
        marker.list[q1] <-"T cell receptor co-receptor CD8"
    }
    if(marker.list[q1]=="IgD"){
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"IgD immunoglobulin complex\"\n" ) }
        marker.list[q1] <-"IgD immunoglobulin complex"
    }
    if(marker.list[q1]=="CD3"){
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"alpha-beta T cell receptor complex\"\n" ) }
        marker.list[q1] <-"alpha-beta T cell receptor complex" 
    }
    if(marker.list[q1]=="HLA-DR"){
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[q1], "has been FORCED to update to \"MHC class II histocompatibility antigen alpha chain HLA-DRA\"\n" ) }
        marker.list[q1] <-"MHC class II histocompatibility antigen alpha chain HLA-DRA"
    }
    return ( marker.list )
}
#################################################################################
# Searches the ontology to find the correct label for each marker
ontologyLabel <- function ( marker.list, marker.list.short, CompInfo="", save.dir="", que.hasProperLabel="", que.hasProperSynonym="", prefix.info="" ) {
  
    if ( length ( marker.list[["Positive"]] ) != 0 ) {            
    for ( q1 in 1:length ( marker.list[["Positive"]] ) ) {
    
    # skips query if there is a file named "markers_ShortName_OntologyName" with all the marker ontology labels
    fname <-  paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
    if(file.exists(fname)){
        markers_ShortName_OntologyName <- read.csv ( fname , header = FALSE )    
        markers_ShortName_OntologyName <- as.matrix ( markers_ShortName_OntologyName )

        temp.location <- which ( markers_ShortName_OntologyName[,1] == marker.list.short[["Positive"]][q1] )
        if ( length ( temp.location ) == 1 ) {
            if ( CompInfo == TRUE ) {  cat ( "Marker", marker.list[["Positive"]][q1], "is called", markers_ShortName_OntologyName[temp.location,2], "\n" ) }
            marker.list[["Positive"]][q1] <- markers_ShortName_OntologyName[temp.location,2]
            next
        }
    }
      
    # Change the short name markers in marker.list to the marker labels in the ontology.
    
    # This is only needed for the ones that the code has trouble finding.
    # Hopefully with an updated ontology these next three line can be deleted.
    updatedException <- ontolExceptionsPos1 ( marker.list[["Positive"]], q1 )
    marker.list[["Positive"]] <- ontolExceptionsPos2 ( marker.list[["Positive"]], q1, CompInfo )
    if ( updatedException == TRUE ) { next }
    
    temp.marker <- marker.list[["Positive"]][q1]
    # Define the cell.ctde.net SPARQL endpoint
    endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"

    # Concatenate the query preceded by all prefix information as a single string to be passed to SPARQL
    query <- paste ( c ( prefix.info, que.hasProperLabel ), collapse="\n" )

    # Prepare marker for query by ensuring the marker is either followed by the  
    # end of the line or it has a symbol other than a letter, a number or a space with 
    # 'a', 'b', 'e' or 's' after it, as that may indicate a different marker.
    temp.marker <- paste ( temp.marker, "$|", temp.marker, "[^0-9a-zA-Z-+/][^ abes]", collapse="", sep = "" )
    query <- gsub ( "\\$marker", temp.marker, query )
    
    # Execute query
    res <- SPARQL ( url = endpoint, query )$results

    if ( nrow ( res ) == 1 ) {
        temp.marker <- res[1,'label']          
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Positive"]][q1], "has been updated to", temp.marker, "\n" ) }
        marker.list[["Positive"]][q1] <- temp.marker
      
        # Update the markers_ShortName_OntologyName.csv file
        fname <-  paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
        temp.table <- read.csv ( fname, header = FALSE ) ;   temp.table <- as.matrix ( temp.table )
        temp.matrix <- matrix ( 0, length ( temp.table[,1] ) + 1, 2 )
        temp.matrix[1:length ( temp.table[,1] ) , ] <- temp.table
        temp.table <- temp.matrix
        temp.table[length ( temp.table[,1] ), 1] <- as.character ( marker.list.short[["Positive"]][q1] )    
        temp.table[length ( temp.table[,1] ), 2] <- paste ( temp.marker )
        write.table(temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE)  
    }
    if ( nrow ( res ) >= 2 ) {
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Positive"]][q1],"has not been changed, multiple possible markers in label\n" ) }
    }
    if ( nrow ( res ) == 0 | nrow ( res ) >= 2 ) {
        temp.marker <- marker.list[["Positive"]][q1]
        # Define the cell.ctde.net SPARQL endpoint
        endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"

        # Concatenate the query preceded by all prefix information as a single string to be passed to SPARQL
        query <- paste ( c ( prefix.info, que.hasProperSynonym ), collapse="\n" )

        # Prepare marker for query by ensuring the marker is either followed by the  
        # end of the line or it has a symbol other than a letter, a number or a space with 
        # 'a', 'b', 'e' or 's' after it, as that may indicate a different marker.
        temp.marker <- paste ( temp.marker, "$|", temp.marker, "[^0-9a-zA-Z-+/][^ abes]", collapse="", sep = "" )
        query <- gsub ( "\\$marker", temp.marker, query )

        # Execute query
        res <- SPARQL ( url=endpoint, query )$results
 
        # Small loop to check if the query is giving multiple label names. In this case the marker will not be changed
        if ( nrow ( res ) >= 1 ) {
            temp = res[1,'label']
            for ( q2 in 1:nrow ( res ) ) {
                if ( temp == res[q2,'label'] ) {
                    temp = res[q2,2]
                    sameLabels <- TRUE;
                } else {
                    sameLabels <- FALSE;
                    break
                }
            }
          
            if ( sameLabels == TRUE ) { # label exists and there is only one label
                temp.marker <- res[1,'label']
                if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Positive"]][q1], "has been updated to", temp.marker, "\n" ) }
                marker.list[["Positive"]][q1] <-temp.marker  
            
                # Update the markers_ShortName_OntologyName.csv file
                fname <-  paste ( save.dir,"markers_ShortName_OntologyName.csv",sep = "" )
                temp.table <- read.csv ( fname, header = FALSE);       temp.table <- as.matrix ( temp.table)
                temp.matrix <- matrix ( 0,length ( temp.table[,1] ) + 1, 2 )
                temp.matrix[1:length ( temp.table[,1] ), ] <- temp.table
                temp.table <- temp.matrix
                temp.table[length ( temp.table[,1] ), 1] <- as.character ( marker.list.short[["Positive"]][q1] )    
                temp.table[length ( temp.table[,1] ), 2] <-paste ( temp.marker )
                write.table ( temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE )  
            }
            if ( sameLabels == FALSE ) { # label exists however there is two or more labels
                if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Positive"]][q1], "has not been changed, multiple possible markers in synonyms\n" ) }
            }
        }
        if ( nrow ( res ) == 0 ) { # label does not exist
            if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Positive"]][q1], "could not be found\n" ) }
        }
    }
    } # end of for-loop
    } # end of if statement
  
    # The negative case
    if ( length ( marker.list[["Negative"]] ) != 0 ) {
    for ( q1 in 1:length ( marker.list[["Negative"]] ) ) {

    # skips query if there is a file named "markers_ShortName_OntologyName" with all the marker ontology labels
    fname <-  paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
    if(file.exists(fname)){
        markers_ShortName_OntologyName <- read.csv ( fname , header = FALSE )    
        markers_ShortName_OntologyName <- as.matrix ( markers_ShortName_OntologyName )
        temp.location <- which ( markers_ShortName_OntologyName[,1] == marker.list.short[["Negative"]][q1] )
        if ( length ( temp.location ) == 1 ) {
            if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Negative"]][q1], "is called", markers_ShortName_OntologyName[temp.location,2], "\n" ) }
            marker.list[["Negative"]][q1] <- markers_ShortName_OntologyName[temp.location,2]
            next
        }
    }
    
    # Change the short name markers in marker.list to the marker labels in the ontology.
    
    # This is only needed for the ones that the code has trouble finding.
    # Hopefully with an updated Ontology these next three line can be deleted.
    updatedException <- ontolExceptionsNeg1 ( marker.list[["Negative"]], q1 )
    marker.list[["Negative"]] <- ontolExceptionsNeg2 ( marker.list[["Negative"]], q1, CompInfo )
    if ( updatedException == TRUE ) { next }
    
    temp.marker <- marker.list[["Negative"]][q1]
    # Define the cell.ctde.net SPARQL endpoint
    endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"
    
    # Concatenate the query preceded by all prefix information as a single string
    # to be passed to SPARQL
    query <- paste ( c ( prefix.info, que.hasProperLabel), collapse="\n" )
    
    # Prepare marker for query by ensuring the marker is either followed by the  
    # end of the line or it has a symbol other than a letter, a number or a space with 
    # 'a', 'b', 'e' or 's' after it, as that may indicate a different marker.
    temp.marker <- paste ( temp.marker, "$|", temp.marker, "[^0-9a-zA-Z-+/][^ abes]", collapse="", sep = "" )
    query <- gsub ( "\\$marker", temp.marker, query )

    # Execute query    
    res <- SPARQL ( url=endpoint, query )$results
    
    if ( nrow ( res ) == 1 ) {
        temp.marker <- res[1,'label']
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Negative"]][q1], "has been updated to", temp.marker, "\n" ) }
        marker.list[["Negative"]][q1] <-temp.marker
       
        # Update the markers_ShortName_OntologyName.csv file
        fname <-  paste ( save.dir,"markers_ShortName_OntologyName.csv",sep = "" )
        temp.table <- read.csv ( fname, header = FALSE ) ;    temp.table <- as.matrix ( temp.table )
        temp.matrix <- matrix ( 0,length ( temp.table[,1] ) + 1, 2 )
        temp.matrix[1:length ( temp.table[,1] ), ] <- temp.table
        temp.table <- temp.matrix
        temp.table[length ( temp.table[,1] ), 1] <- as.character ( marker.list.short[["Negative"]][q1] );    
        temp.table[length ( temp.table[,1] ), 2] <- paste ( temp.marker )
        write.table ( temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE )
    }
    if ( nrow ( res ) >= 2 ) {
        if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Negative"]][q1],"has not been changed, multiple possible markers in label\n" ) }
    }
    if ( nrow ( res ) == 0 | nrow ( res ) >= 2 ) {
        temp.marker <- marker.list[["Negative"]][q1]
        # Define the cell.ctde.net SPARQL endpoint.
        endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"
      
        # Concatenate the query preceded by all prefix information as a single string
        # to be passed to SPARQL
        query <- paste ( c ( prefix.info, que.hasProperSynonym ), collapse="\n" )
      
        # Prepare marker for query by ensuring the marker is either followed by the  
        # end of the line or it has a symbol other than a letter, a number or a space with 
        # 'a', 'b', 'e' or 's' after it, as that may indicate a different marker.
        temp.marker <- paste ( temp.marker, "$|", temp.marker, "[^0-9a-zA-Z-+/][^ abes]", collapse="", sep = "" )
        query <- gsub ( "\\$marker", temp.marker, query )
      
        # Execute query
        res <- SPARQL ( url=endpoint, query )$results
      
        # small loop to check if the multiple hits are giving the same result
        if ( nrow ( res ) >= 1 ) {
            temp = res[1,'label']
            for ( q2 in 1:nrow ( res ) ) {
                if ( temp == res[q2,'label'] ) {
                    temp = res[q2,'label']
                    sameLabels <- TRUE;
                } else {
                    sameLabels <- FALSE;
                    break
                }
            }
            if ( sameLabels == TRUE ) { # label exists and there is only one label
                temp.marker <- res[1,'label']
                if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Negative"]][q1], "has been updated to", temp.marker, "\n" ) }
                marker.list[["Negative"]][q1] <-temp.marker
          
                # Update the markers_ShortName_OntologyName.csv file
                fname <-  paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
                temp.table <- read.csv ( fname, header = FALSE );   temp.table <- as.matrix ( temp.table )
                temp.matrix <- matrix ( 0,length ( temp.table[,1] ) + 1, 2 )
                temp.matrix[1:length ( temp.table[,1] ), ] <- temp.table
                temp.table <- temp.matrix
                temp.table[length ( temp.table[,1] ), 1] <- as.character ( marker.list.short[["Negative"]][q1] )
                temp.table[length ( temp.table[,1] ), 2] <- paste ( temp.marker )
                write.table ( temp.table, fname, sep = ",", col.names = FALSE, row.names = FALSE )
            }
            if ( sameLabels == FALSE ) { # label exists however there is two or more labels
                if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Negative"]][q1], "has not been changed, multiple possible markers in synonyms\n" ) }
            }
        }
        if ( nrow ( res ) == 0 ) { # label does not exist
            if ( CompInfo == TRUE ) { cat ( "Marker", marker.list[["Negative"]][q1], "could not be found\n" ) }
        }
    }
    } # end of for-loop
    } # end of if statement

    return ( marker.list )  
}

#################################################################################
# Similar to queryMarker, but for finding the parents for a specific cell population
# by cell label.
# NOTE: I tried matching by ID and for some reason I could do it on the HE group
# but not here for some reason, that's why I'm matching by label...
parentQuery <- function ( child.label = "common myeloid progenitor", 
                        query.file = "getParentClasses.txt", prefix.info="" ) {
  
    # Define the cell.ctde.net SPARQL endpoint
    endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"

    # Concatenate the query preceded by all prefix information as a single string
    # to be passed to SPARQL.
    query <- paste ( c ( prefix.info, query.file ), collapse="\n" )
    # Add "^" to beginning and "$" to end to find an exact match for the label.
    child.label <- paste ( "^", child.label, "$", sep = "", collapse = "" )
    query <- gsub ( "\\$label", child.label, query )
    # Execute query
    res <- SPARQL ( url=endpoint, query )$results 
    return ( res )
}

#################################################################################
# Makes a query with SPARQL to the CL Ontology, and returns the results.
queryMarker <- function ( marker = NULL, query.file = "getMatchingSynonyms.txt",
                        celllabel = NULL, prefix.info="" ) {

    # TO DO: improve input check
    # For now, a NULL, NA or length 0 phenotype will have empty results
    #  matrix returned.
    if ( !is.character ( marker ) || length ( marker ) == 0 ) {
        warning ( "No marker found." )
        res <- matrix ( ncol = 2, nrow = 0 )
        colnames ( res ) <- c ( "ID", "Synonym Match" )
        return ( res )
    }

    # Define the cell.ctde.net SPARQL endpoint.
    endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"

    # Concatenate the query preceded by all prefix information as a single string
    # to be passed to SPARQL.
    query <- paste ( c ( prefix.info, query.file ), collapse="\n" )

    # Prepare marker for query by ensuring the marker is either followed by the  
    # end of the line or it has a symbol other than a letter, a number or a space with 
    # 'a', 'b', 'e' or 's' after it, as that may indicate a different marker.
    marker <- paste ( marker, "$|", marker, "[^0-9a-zA-Z-+/][^ abes]", collapse="", sep = "" )
    query <- gsub ( "\\$marker", marker, query )
    if ( !is.null ( celllabel ) ) {
        celllabel <- paste ( "^", celllabel, "$", sep = "", collapse = "" )
        query <- gsub ( "\\$celllabel", celllabel, query )
    }
    # Execute query
    res <- SPARQL ( url = endpoint, query )$results
    return ( res )
}


