flowCL <- function ( MarkerList = "HIPC", Indices = NULL, CompInfo = FALSE, KeepArch = TRUE,   
                     MaxHitsPht = 5, OntolNamesTD = FALSE, ResetArch = FALSE, VisualSkip = FALSE ) {

# flowCL (Semantic labelling of flow cytometric cell populations)
# Author: Justin Meskas (jmeskas@bccrc.ca)
# Date last modified: April 30, 2014
# Author: Radina Droumeva (radina.droumeva@gmail.com)
# Date last modified: July 19, 2013
    
initialTime <- Sys.time ( )

VisualizationSkip <- VisualSkip       # Skips the making of the visualization (tree diagram)
penalty.skip      <- TRUE             # Skips the penalty calculation (slowest part of the code, and not that important)
OntolNamesTD      <- OntolNamesTD     # Shows either the common marker name or the ontology marker name in the tree diagram

finalResults <- list ( )

# Load some required data from supportingFunctions.
listPhenotypes_flowCL <- listPhenotypes_flowCL()
listColours_flowCL    <- listColours_flowCL()

# Load more required data from supportingFunctions. 
prefix.info                     <- flowCL_query_data_prefix.info()
que.getParentClasses            <- flowCL_query_data_getParentClasses()
que.hasPlasmaMembranePart       <- flowCL_query_data_hasPlasmaMembranePart()
que.hasProperLabel              <- flowCL_query_data_hasProperLabel()
que.hasPMPsingle                <- flowCL_query_data_hasPMPsingle()
que.hasProperSynonym            <- flowCL_query_data_hasProperSynonym()
que.lacksPMPsingle              <- flowCL_query_data_lacksPMPsingle()
que.lacksPlasmaMembranePart     <- flowCL_query_data_lacksPlasmaMembranePart()
que.hasLowPlasmaMembraneAmount  <- flowCL_query_data_hasLowPlasmaMembraneAmount()
que.hasHighPlasmaMembraneAmount <- flowCL_query_data_hasHighPlasmaMembraneAmount()

# Function to test if the user has input cd in non-capitals. If so, the code will terminate.
if ( cdTest ( MarkerList ) == TRUE)
    return()

# Check date of Ontology update
if ( MarkerList == "Date"||MarkerList == "date"||MarkerList == "DATE" ) {
    que.getDate <- flowCL_query_date()    
    endpoint <- "http://cell.ctde.net:8080/openrdf-sesame/repositories/CL"
    res <- SPARQL ( url = endpoint, paste(que.getDate, collapse="\n") )$results
    return ( c(paste(strsplit(res[1,2],split="/")[[1]][7],": year-month-day")))
}

# Load phenotypes.
suppressWarnings ( if ( MarkerList == "HIPC" ) {
    listPhenotypes <- listPhenotypes_flowCL
} else {
    listPhenotypes <- MarkerList
} )
listPhenotypes <- as.matrix ( listPhenotypes )

# Start and end of the iterations from MarkerList.
if ( is.null ( Indices ) ) {
    IterStart = 1
    IterEnd   = length ( listPhenotypes )
    markersToQuery <- IterStart:IterEnd
} else {
    listPhenotypes <- listPhenotypes[Indices]
    markersToQuery <- 1:length(Indices)
}

# Preallocate lists.
listPhenotypeUpdate  <- listPhenotypeOriginal <- listPhenotypeID <- listCellID <- listPhenotypes
listPhenotypesuccess <- listMarkerLabels      <- listCellLabels                <- listPhenotypes

# Default listPhenotypesuccess as a "No"
for ( q in markersToQuery ) {    
    listPhenotypesuccess[[q]] <-   "No"
}

# Removes the flowCL_results directory if ResetArch is TRUE
if ( ResetArch == TRUE ) {
    unlink ( "flowCL_results", recursive = TRUE ) 
}

# Define directories for storing results.
save.dir             <- 'flowCL_results/'
save.dirResults      <- 'flowCL_results/results/'
save.dirParents      <- 'flowCL_results/parents/'
save.dirParentsQuery <- 'flowCL_results/parents_query/'
# Create directories.
dir.create ( save.dir,             showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirResults,      showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirParents,      showWarnings=FALSE, recursive=TRUE )
dir.create ( save.dirParentsQuery, showWarnings=FALSE, recursive=TRUE )

# Create a file for future quick matching of the short name with the ontology name of each marker.
fname <- paste ( save.dir, "markers_ShortName_OntologyName.csv", sep = "" )
if( !file.exists ( fname ) ) {
    temp.table <- matrix ( )
    temp.table[1] <- paste ( 'CD4+' );    temp.table[2] <- paste ( 'CD4 molecule' )
    write.table ( t ( temp.table ), fname, sep = ",", col.names = FALSE, row.names = FALSE )  
}

# Create null variables for memory allocation and easy additions of indices
rEG <- nAttrs <- eAttrs <- attrs <- list ( )

#################################################################################
# Main body starts here. This for-loop goes through each phenotype one at a time.
for ( q in markersToQuery ) {
  
# Get starting time to keep track of each phenotype's running time.
start <- Sys.time ( )

# Define a phenotype to test here.
phenotype <- listPhenotypes[[q]]
# Change the phenotype to have only hi and lo and not other variations
phenotype <- sub("bright", "hi", phenotype)
phenotype <- sub("bri",    "hi", phenotype)
phenotype <- sub("high",   "hi", phenotype)
phenotype <- sub("dim",    "lo", phenotype)
phenotype <- sub("low",    "lo", phenotype)
    
# Output for the user
if ( CompInfo == TRUE ) { cat("\nThe phenotype of interest is", phenotype, "\n" ) }

# The following breaks up the phenotype into single markers and sorts them
# by their expression (only positive or negative expression implemented for now)
# see supportingFunctions.R for details on how 'phenoParse' works.
marker.list <- phenoParse ( phenotype )

# Change the . in HLA.DR to a - since the + and - signs are reserved for splitting the string.
marker.list <- changeHLADR ( marker.list )
# Creates another copy which will have the + and - signs. Used by treeDiagram and when searching files.
marker.list.short <- marker.list

# Put the +, -, hi and lo signs back on marker.list.short. Used by treeDiagram and when searching files.
AddOn <- c("+", "-", "hi", "lo")
for ( q3 in 1:length(marker.list)){
    if ( length ( marker.list.short[[q3]] ) != 0 ) {
        for ( q2 in 1:( length ( marker.list.short[[q3]] ) ) ) {
            marker.list.short[[q3]][q2] <- paste ( marker.list.short[[q3]][q2], AddOn[q3], sep = "" )
        }
    }
}

# Update the marker list with the full label names from the ontology.
marker.list <- ontologyLabel ( marker.list, marker.list.short, CompInfo = CompInfo, save.dir = save.dir, 
                                que.hasProperLabel = que.hasProperLabel, que.hasProperSynonym = que.hasProperSynonym, prefix.info = prefix.info )
# Make a list of the ontology names for each phenotype searched for, which will be exported to a table in .csv form
listPhenotypeUpdate[[q]]   <- phenoUnparse ( phenotype, marker.list )

# Default skip Query to TRUE. If it changes to FALSE then querying will have to be done.
skipQuery <- TRUE
# Create res as a NULL variable for easy additions of indices
res <- NULL
# Cycle through positive and negative markers.
# For positive, use a SPARQL query using the property 'has plasma membrane part'
# of each marker. Then for the negative markers, query 'lacks plasma membrane part'.
for ( i in unlist ( marker.list.short ) ) {
    fname <- paste ( save.dirResults, 'results_', i, '.csv', sep = "" )
    if ( file.exists ( fname ) ) {
        tempCSV <- read.csv ( fname, as.is=TRUE )
        if ( length ( which ( tempCSV[,"Number.Of.Hits"] > 1 )) >= 1 ) { # Number Of Hits column
            warning ( "You are receiving more marker labels than markers that were input. Most likely the input markers are not defined correctly." )
        }
        # delete the last two columns to combine and redue the results of the last two columns with other markers
        tempCSV <- tempCSV[ - c ( ncol ( tempCSV ) , ncol ( tempCSV ) - 1 ) ] 

        res <- rbind ( res, tempCSV )
    } else {
        if ( CompInfo == TRUE ) { cat ( "At least one marker was not previously queried. Querying all.\n" ) }
        skipQuery <- FALSE
        break;
    }               
}

# If all results files exist then no querying needs to be done
if ( skipQuery == TRUE ) {
    
    clean.res <- tabulateResults ( res )
  
} else {
    # Initialize result collector.
    res <- NULL
    for ( q3 in 1:length(marker.list)){
        for ( m in marker.list[[q3]] ) {
            if ( CompInfo == TRUE ) { cat( "Locating marker", m, "\n" ) }
            # Get relevant information about the marker - the population names which
            # have "plasma membrane part" of this marker.
            if ( q3 == 1 ) { cur.res <- queryMarker ( marker = m, query.file = que.hasPlasmaMembranePart,       prefix.info = prefix.info ) }
            if ( q3 == 2 ) { cur.res <- queryMarker ( marker = m, query.file = que.lacksPlasmaMembranePart,     prefix.info = prefix.info ) }
            if ( q3 == 3 ) { cur.res <- queryMarker ( marker = m, query.file = que.hasHighPlasmaMembraneAmount, prefix.info = prefix.info ) }
            if ( q3 == 4 ) { cur.res <- queryMarker ( marker = m, query.file = que.hasLowPlasmaMembraneAmount,  prefix.info = prefix.info ) }
            if ( nrow ( cur.res ) == 0 ) {
                if ( CompInfo == TRUE ) { cat( "No cell populations found which have plasma membrane part", m, "!\n" ) }
            }
    
            # VERY INEFFICIENTLY, cycle through the other markers and assign a penalty 
            # according to a contradiction count.
            penalties <- rep ( 0, nrow ( cur.res ) )
            # Skips the penalty calculator if TRUE (very slow part of the code).
    #         if ( penalty.skip == FALSE ) {
    #         for ( other.marker in marker.list[["Negative"]] ) {
    #             for ( i in 1:nrow ( cur.res ) ) {
    #                 temp <- queryMarker ( marker = other.marker, query.file = que.hasPMPsingle,
    #                                         celllabel = cur.res[i, 'celllabel'], prefix.info = prefix.info )
    #                 if ( nrow ( temp ) > 0 ) {
    #                   # If one of the negatively expressed markers was found for the current
    #                   # cell type as "has plasma membrane part", that is a contradiction which
    #                   # should be penalized, proportionally to the number of markers in total.
    #                   penalties[i] <- penalties[i] - 1 / length ( unlist ( marker.list ) )
    #                 }
    #             }
    #         }
    #         # Similarly for the positive markers showing as negative.
    #         for ( other.marker in setdiff ( marker.list[["Positive"]], m ) ) {
    #             for ( i in 1:nrow ( cur.res ) ) {
    #                 temp <- queryMarker ( marker = other.marker, query.file = que.lacksPMPsingle,
    #                                         celllabel = cur.res[i, 'celllabel'], prefix.info = prefix.info )
    #                 if ( nrow ( temp ) > 0 ) {
    #                     # If one of the other positively expressed markers was found for the current
    #                     # cell type as "lacks plasma membrane part", that is a contradiction which
    #                     # should be penalized, proportionally to the number of markers in total.
    #                     penalties[i] <- penalties[i] - 1 / length ( unlist ( marker.list ) )
    #                 }
    #             }
    #         }
    #         } # end of penalty skip if statement.
            
            
            cur.res <- cbind ( cur.res, penalties )
            temp.clean.res <- tabulateResults ( cur.res )
            temp.name <- marker.list.short[[q3]][which ( m == marker.list[[q3]] )]
            fname <- paste ( save.dirResults, 'results_', temp.name, '.csv', sep = "" )
            # save the marker's data to skip the query next time
            write.csv ( temp.clean.res, fname, row.names = FALSE )
            res <- rbind ( res, cur.res )
            if ( length ( which ( temp.clean.res[,"Number Of Hits"]  > 1 )) >= 1 ) { # Number Of Hits column
                warning ( "You are receiving more marker labels than markers that were input. Most likely the input markers are not defined correctly." )
            }
        }
    }

    # For each returned owl object ID, tabulate how many times the result was 
    # returned. This is essentially the "hits" part of the score -- telling us
    # how many of the markers we have were matched to each population. From this,
    # the penalty tally will be subtracted to get a final overall score for each
    # population label.
    clean.res <- tabulateResults ( res )

}
    
# Save the results. 
fname <- paste ( save.dirResults, 'results_', phenotype, '.csv', sep = "" )
write.csv ( clean.res, fname, row.names = FALSE )
if ( CompInfo == TRUE ) { cat( 'Initial query results saved in', fname, "\n" ) }

# If there were no hits the code will skip to the next iteration instead of 
# proceeding and getting an error
if ( nrow ( clean.res ) == 0 ) {
    if ( CompInfo == TRUE ) { 
        cat ( "Time elapsed:", timeOutput(start), "\n" ) 
        cat ( "Iterations at", which ( markersToQuery == q ) , "out of", length ( markersToQuery ), "\n" )
    }  
    next
}

# Added in to make the code faster when the user only wants to know if the markers 
# are in the ontology or not and does not want the tree diagrams
if ( VisualizationSkip == FALSE ) {

# Refine list of parents by focusing on highest scores:
scores <- as.numeric ( as.character ( clean.res[ , "Score"] ) )
       
# If at least 3 perfect scores are found, only report them:
if ( length ( which ( scores == length ( unlist ( marker.list ) ) ) ) >= 3 ) {
    cutoff.score <- length ( unlist ( marker.list ) )
} else {
    # Typically there will be less than perfect scores (due to marker missing/
    # marker not declared in the population/non-typical phenotype)
  
    # Only show the matches with the most hits. If there are only a few that are off by 1 then show them as well.
    if ( length ( which ( as.integer ( clean.res[,"Number Of Hits"] ) >= ( as.integer ( clean.res[1,"Number Of Hits"] ) - 1 ) ) ) <= 4 ) {
        cutoff.score <- as.integer ( clean.res[1,"Number Of Hits"] ) - 1
    } else {
        cutoff.score <- as.integer ( clean.res[1,"Number Of Hits"] )
    }
    #cutoff.score <- getScoreCutoff(scores) #quantile(scores, 0.85) # (Old version of cutoff.score)
}
if (length ( which ( scores == length ( unlist ( marker.list ) ) ) ) >= 1 ) {
    if ( CompInfo == TRUE ) { cat ( "Perfect match(es) found.\n" ) } 
} else {
    if ( CompInfo == TRUE ) { cat ( "Using a cutoff score of", cutoff.score, "\n" ) }
}

# Extract highly scored parents only:
clean.res <- clean.res[which ( scores >= cutoff.score ), ]

# Fixes small bug when there is only one result
if ( length ( which ( scores >= cutoff.score ) ) == 1 ) {
    clean.res <- as.matrix ( clean.res )
    clean.res <- t ( clean.res )
}

# Identify all parents of the matches.
parent.res <- matrix ( nrow = 0, ncol = 5 )
colnames ( parent.res ) <- c ( "x", "celllabel", "parent", "parentlabel", "score" )
for ( i in 1:nrow ( clean.res ) ) {
    # Fixes a bug when trying to save cell type of "F4/80-negative adipose macrophage".
    # the "sub" replaces the "/" with a "_" to avoid R thinking that F4 is a folder.
    fname <- paste ( save.dirParentsQuery, sub ( "/","_", ( clean.res[i, "celllabel"] ) ), '.csv', sep = "" )
    if ( file.exists ( fname ) ) {  # read from file if there is one, instead of querying the same thing again
        res <- read.csv ( fname, as.is = TRUE )
    } else { 
        res <- parentQuery ( child.label =  clean.res[i, "celllabel"], query.file = que.getParentClasses, prefix.info )
        write.csv ( res, fname, row.names = FALSE )      
    }
    res <- cbind ( res, rep ( clean.res[i, "Score"], nrow ( res ) ) )
    colnames ( res ) [ ncol ( res ) ] <- "Score"
    parent.res <- rbind ( parent.res, res )  
}
# Save parent information.
fname <- paste ( save.dirParents, 'parent.res', phenotype, '.csv', sep = "" )
write.csv ( parent.res, fname, row.names = FALSE )
if ( CompInfo == TRUE ) { cat ( 'Parent information saved in', fname, "\n" ) }

summary <- table ( parent.res [ , 'parentlabel'] )

scores2 <- sapply ( names ( summary ), function ( p ) mean ( as.numeric ( as.character ( parent.res [ which ( parent.res [ , 'parentlabel'] == p ), 'Score'] ) ) ) )

# The higher up in the tree -- the more times a population is called a parent to 
# others -- the more certain we are that the label applies (e.g. 'cell' is usually
# a parent to most population labels under investigation, and we are pretty sure
# whatever phenotype we are working with is at least a cell!)
# The higher the Score is, the more markers in our phenotype matched. Combining 
# these two measures gives an overall estimate of how specific and how reliable
# the hits are.
scored.summary <- summary / max ( summary ) * scores2 / max ( scores2 )
s <- order ( scores2, scored.summary, decreasing = TRUE )
  
# Create a list of the populations and their parents for visualization purposes.
parent.analysis <- NULL
for ( q1 in 1:length ( s ) ) {
    # Fixes a bug when trying to save cell type of "F4/80-negative adipose macrophage".
    # the "sub" replaces the "/" with a "_" to avoid R thinking that F4 is a folder.    
    fname <- paste ( save.dirParentsQuery, sub ( "/", "_", names ( summary ) [ s[q1] ] ), '.csv', sep = "" )
    if ( file.exists ( fname ) ) { # read from file if there is one, instead of querying the same thing again
        res <- read.csv ( fname, as.is = TRUE )
    } else { 
        res <- parentQuery ( child.label = names ( summary ) [ s[q1] ], query.file = que.getParentClasses, prefix.info )
        write.csv ( res, fname, row.names = FALSE )      
    }
    parent.analysis [ length ( parent.analysis ) + 1 ] <- list ( res [ which ( is.element ( res [ , 'parentlabel'], names ( summary ) ) ), 'parentlabel'] )
}
names ( parent.analysis ) <- names ( summary ) [s]
 
# Also create a list of the populations and their children for visualization purposes.
child.analysis <- lapply ( names ( summary ) [s], function(x) { 
    children <- c()
    for ( y in setdiff ( names ( summary ), x ) ) {
        if ( is.element ( x, parent.analysis[[y]] ) ) {
            children <- c ( children, y )
        }
    }
    return ( children )
})
names ( child.analysis ) <- names ( summary ) [s]
           
# Create and save a pdf file of a tree diagram
finalResults[[listPhenotypes[q]]] <- treeDiagram ( child.analysis, clean.res, phenotype, OntolNamesTD, marker.list.short, marker.list, save.dir, listColours_flowCL = listColours_flowCL )

}# end of visualization if statement

# Check to see if all markers were hits (A perfect match).
if ( max ( clean.res[,"Number Of Hits"] ) == ( length ( unlist ( marker.list ) ) ) ) {
    if ( length ( which ( clean.res[,"Number Of Hits"] == length ( unlist ( marker.list ) ) ) ) == 1 ) {
        listPhenotypesuccess[[q]] <- ( "Yes" )
    } else {
        listPhenotypesuccess[[q]] <- paste( "Yes - ", length ( which ( clean.res[,"Number Of Hits"] == length ( unlist ( marker.list ) ) ) ), " hits", sep = "" )
    }
}

# Compile the new row of the .csv file by updating all the list functions. 
r <- updateLists ( clean.res, MaxHitsPht )
listMarkerLabels[[q]] <- r[1] ; listCellLabels[[q]] <- r[2]
listPhenotypeID[[q]]  <- r[3] ; listCellID[[q]]     <- r[4]

# The time for one marker/phenotype to be queried
if ( CompInfo == TRUE ) { cat ( "Time elapsed:", timeOutput(start), "\n" ) }
# Shown so the user knows how far the code has run
if ( CompInfo == TRUE ) { cat ( "Iterations at", which(markersToQuery==q) , "out of",length(markersToQuery), "\n" ) }

} # end of main body for-loop

#################################################################################

# Creating a .csv file with the short name phenotypes and the ontology name phenotypes with a "Yes" or "No" indicating if it is in the ontology or not.
# Plus a list of the markers and the cell types with their ontology IDs for the cases with the maximum number of hits.
if ( length ( listPhenotypes ) == 1 ) {
    transpose.listP <- TRUE # if there is only one row
} else {
    transpose.listP <- FALSE # if there is two or more rows
}

listPhenotypes <- cbind ( listPhenotypes, listPhenotypeUpdate, listPhenotypesuccess, listPhenotypeID, listMarkerLabels, listCellID, listCellLabels )
colnames(listPhenotypes) <- c("Short marker names","Ontology marker names","Successful Match?","Marker ID","Marker Label","Cell ID","Cell Label")
# If there is only one marker queried the results will only have one row and this cases a small error. "transpose.listP" fixes this
if ( transpose.listP == FALSE ) { # if there is two or more rows
    listPhenotypes <- listPhenotypes[markersToQuery,]
}

# Save the results into a .csv file for outside of R viewing
fname <- paste ( save.dir, 'listPhenotypes.csv', sep = "" )
write.csv ( listPhenotypes, fname, row.names = FALSE )

if ( transpose.listP == FALSE ) { # if there is two or more rows
    # Format the .csv file into a table format to show in R
    listPhenotypes[,2] <- gsub ( "\n", ", ", listPhenotypes[,2] )
    for ( q2 in 4:7 ) {
        for ( q1 in  1:MaxHitsPht ) {
            listPhenotypes[,q2] <- gsub ( paste ( "\n", q1, sep = "" ) , paste ( " ", q1, sep = "" ), listPhenotypes[,q2] )
        }
        listPhenotypes[,q2] <- gsub ( "\n"  , ", ", listPhenotypes[,q2] )
        listPhenotypes[,q2] <- gsub ( ",  +", " ", listPhenotypes[,q2] )
    }
}
if ( transpose.listP == TRUE ) { # if there is only one row
    # Format the .csv file into a table format to show in R
    listPhenotypes[2] <- gsub ( "\n", ", ", listPhenotypes[2] )
    for ( q2 in 4:7 ) {
        for ( q1 in  1:MaxHitsPht ) {
            listPhenotypes[q2] <- gsub ( paste ( "\n", q1, sep = "" ) , paste ( " ", q1, sep = "" ), listPhenotypes[q2] )
        }
        listPhenotypes[q2] <- gsub ( "\n"  , ", ", listPhenotypes[q2] )
        listPhenotypes[q2] <- gsub ( ",  +", " ", listPhenotypes[q2] )
    }    
}

SegmentTotal <- list()
# Create the output if []$Table is called for multiple phenotypes
if ( transpose.listP == FALSE ) { # if there is two or more rows
    for ( l1 in 1:length ( listPhenotypes[,1] ) ) {
        listTem <- t ( listPhenotypes[l1,] )

        columnof_listPhenotypes <- colnames ( listPhenotypes ) 
        segment <- strwrap ( listTem, 70 )
        
        segIndices <- c( 1, 2, union ( which ( regexpr ( "No", segment ) == 1 ), which ( regexpr ( "Yes", segment ) == 1 ) ), which ( regexpr ( "1)", segment ) == 1 ) )
        temp <- rep( " ", length ( segment ) )

        temp[segIndices] <- columnof_listPhenotypes[1:length(segIndices)]

        segment <- t ( segment )
        colnames(segment) <- temp

        if ( l1 == 1 ) {
            SegmentTotal[segment[1]] <- list( t ( segment ) )
        } else {
            SegmentTotal[segment[1]] <- list( t ( segment ) )
        }
    }
    
    # Store results in a table
    finalResults["Table"] <- list(( SegmentTotal ))
}
# Create the output if []$Table is called for one phenotype
if ( transpose.listP == TRUE ) { # if there is only one row
    
    listTem <- t ( listPhenotypes )
    
    columnof_listPhenotypes <- colnames ( listPhenotypes ) 
    segment <- strwrap ( listTem, 70 )
    
    segIndices <- c( 1, 2, union ( which ( regexpr ( "No", segment ) == 1 ), which ( regexpr ( "Yes", segment ) == 1 ) ), which ( regexpr("1)", segment ) == 1 ) )

    temp <- rep( " ", length ( segment ) )
    
    temp[segIndices] <- columnof_listPhenotypes[1:length(segIndices)]
    
    segment <- t ( segment )
    colnames(segment) <- temp
    SegmentTotal <-  segment
    # Store results in a table
    finalResults["Table"] <- list ( t ( SegmentTotal ) )
}

# Removes the flowCL_results directory if KeepArch is FALSE
if ( KeepArch == FALSE ) {
    unlink ( "flowCL_results", recursive = TRUE ) 
}

if ( transpose.listP == FALSE ) { # if there is two or more rows
    temp.celllable <- listPhenotypes[,7]
}
if ( transpose.listP == TRUE ) { # if there is only one row
    temp.celllable <- listPhenotypes[7]
}
for ( q4 in markersToQuery ) {
    temp.celllable[q4] <- gsub ( " [+] more", "", temp.celllable[q4] )
    temp.string <- list()
    # Look in the Label of the marker ID for "GO" then pulls out the GO and the numbers that follow
    temp.index <- gregexpr ( ")",  temp.celllable[q4] )
    if ( ( temp.index[[1]][1] != - 1 ) ) {
        for ( q7 in 1:length ( temp.index[[1]] ) ) {
            if ( q7 < length ( temp.index[[1]] ) ) { 
                temp.string[q7] <- paste (substr ( temp.celllable[q4], temp.index[[1]][q7] + 2, temp.index[[1]][q7+1] - 3 ), sep = "" )
            }
            if ( q7 == length ( temp.index[[1]] ) ) { 
                temp.string[q7] <- paste (substr ( temp.celllable[q4], temp.index[[1]][q7] + 2, nchar(temp.celllable[q4])), sep = "" )
            }      
        }
    }
    finalResults[["Cell_Label"]][[listPhenotypes[q4]]] <- temp.string
}

# Output for the user
if ( CompInfo == TRUE ) { 
    cat ( "\nTotal time was: ", timeOutput ( initialTime ), "\n" ) 
    cat ( "Archive saved in \"[current directory]/", save.dir, "\"\n", sep = "" )
}

return ( finalResults )
}

