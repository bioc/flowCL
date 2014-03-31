#################################################################################
# supporting functions (flowCL: Semantic labelling of flow cytometric cell populations)
# Author: Justin Meskas (jmeskas@bccrc.ca)
# Date last modified: March 14, 2014
# Author: Radina Droumeva (radina.droumeva@gmail.com )
# Date last modified: July 19, 2013
#################################################################################

#################################################################################
# Change the "." in HLA.DR to a - since the + and - signs are reserved for splitting the string
# Should change this to allow for any marker with a period to change to a dash
changeHLADR <- function ( marker.list ) {
if ( length ( marker.list[["Positive"]] ) >= 1 ) {
    for ( q1 in 1:length ( marker.list[["Positive"]] ) ) {
        if ( marker.list[["Positive"]][q1] == "HLA.DR" ) {
            marker.list[["Positive"]][q1] <- "HLA-DR"
        }
    }
}
if ( length ( marker.list[["Negative"]] ) >= 1 ) {
    for ( q1 in 1:length ( marker.list[["Negative"]] ) ) {
        if ( marker.list[["Negative"]][q1] == "HLA.DR" ) {
            marker.list[["Negative"]][q1] <- "HLA-DR"
        }
    }
}
return ( marker.list )
}

#################################################################################
# Removes all the non first generation children from the child.analysis and then
# produces a flow chart.

treeDiagram <- function ( child.analysis, clean.res, phenotype, OntolNamesTD, marker.list.short, marker.list, save.dir, listColours_flowCL = "" ) {
    
    # Sort the child.analysis by starting with the cell population which has
    # the most children to the one with the least-- i.e. is the most likely to be the 
    # root parent, such as 'cell' or 'native cell':
    child.lengths <- unlist ( lapply(child.analysis, length ) )
    sort.child <- child.analysis[order ( child.lengths, decreasing = TRUE )]
    labels <- names ( sort.child )
    arrows.colour <- arrows.label <- arrows.dashed.solid <- arrows.head <- NULL
  
    if ( OntolNamesTD == TRUE ) { marker.list.short <- marker.list }
  
    # A quadruple for loop (probably not the most elegant, but it is still fast ). Start with q8,
    # which selects the label with the most children first. q9 stores an element from the q8 label. q10 
    # and q11 look through all other elements in the q8 label and detects if the stored 
    # element (q9) is a parent of other elements in the q8 label. If this is true, then the element
    # q11 is removed from the q8 label.
    # In short, something that is not a first generation child will be removed from the list.
    for ( q8 in 1:length ( labels ) ) {
        deletePoints <- NULL
        for ( q9 in 1:length ( child.analysis[labels[q8]][[1]] ) ) {
      
            if ( is.null ( child.analysis[labels[q8]][[1]][q9] ) ) { next }
      
            temp.label <- child.analysis[labels[q8]][[1]][q9]
      
            for ( q10 in 1:length ( child.analysis[temp.label][[1]] ) ) {
                for ( q11 in 1:length ( child.analysis[labels[q8]][[1]] ) ) {          
          
                    if ( !is.null ( child.analysis[temp.label][[1]][q10] ) ) {
                        if ( ( child.analysis[temp.label][[1]][q10] == child.analysis[labels[q8]][[1]][q11] ) & ( q9!=q11 ) ) {
                            # stores points for removal
                            deletePoints <- c ( deletePoints, q11 )
                        }
                    }
                }
            }
        }
        if ( !is.null ( deletePoints ) ) {
            # removes all non first generation children
            child.analysis[labels[q8]][[1]] <- child.analysis[labels[q8]][[1]][-deletePoints] 
        }
    }  
    # Sets up Rgraphviz with all the node titles in "labels"
    rEG <- new ( "graphNEL", nodes = c ( labels, marker.list.short[["Positive"]], marker.list.short[["Negative"]] ), edgemode="directed" )
      
    for ( q13 in 1:length ( labels ) ) {
        for ( q14 in 1:length ( child.analysis[labels[q13]][[1]] ) ) {
            if ( !is.null ( child.analysis[labels[q13]][[1]][q14] ) ) {
                # Adds an arrow from parent to child
                rEG <- addEdge ( child.analysis[labels[q13]][[1]][q14], labels[q13], rEG, 1 )
                arrows.colour <- c ( arrows.colour, "black" )
                arrows.label  <- c ( arrows.label , paste ( child.analysis[labels[q13]][[1]][q14], "~", labels[q13], sep = "" ) )
                arrows.dashed.solid <- c ( arrows.dashed.solid, "solid" )    
                arrows.head <- c ( arrows.head, "open" )                   
            }
        }
    }
        
    labels.colour <- labels
  
    # Load all predetermined colours
    listColours_flowCL <- as.matrix ( listColours_flowCL )
  
    # Make lists of arrows start and end locations, colour and label
    count.markers <- 0
    for ( q15 in 1:length ( labels ) ) {
        for ( q16 in 1:length ( clean.res[,'celllabel'] ) ) {
            if ( clean.res[q16,'celllabel'] == ( labels[q15] ) ) {
                if ( length ( marker.list[["Positive"]] ) != 0 ) {
                    for ( q17 in 1:length ( marker.list[["Positive"]] ) ) {
                        if ( grepl ( marker.list[["Positive"]][q17],clean.res[q16,'markerlabel'] ) == TRUE ) {
                            count.markers <- count.markers + 1
                            rEG <- addEdge ( marker.list.short[["Positive"]][q17], labels[q15], rEG, 1 )
                            arrows.colour <- c ( arrows.colour, listColours_flowCL[q17] )
                            arrows.label  <- c ( arrows.label, paste ( marker.list.short[["Positive"]][q17],"~",labels[q15], sep = "" ) )
                            arrows.dashed.solid <- c ( arrows.dashed.solid, "dashed" )       
                            arrows.head <- c ( arrows.head, "none" )                
                        }
                    }
                }
                if ( length ( marker.list[["Negative"]] ) != 0 ) {
                    for ( q18 in 1:length ( marker.list[["Negative"]] ) ) {
                        if ( grepl ( marker.list[["Negative"]][q18],clean.res[q16,'markerlabel'] ) == TRUE ) {
                            count.markers <- count.markers + 1
                            rEG <- addEdge(marker.list.short[["Negative"]][q18],labels[q15], rEG, 1 )
                            arrows.colour <- c ( arrows.colour, listColours_flowCL[length ( marker.list.short[["Positive"]] )+q18] )
                            arrows.label  <- c ( arrows.label , paste ( marker.list.short[["Negative"]][q18],"~",labels[q15], sep = "" ) )
                            arrows.dashed.solid <- c ( arrows.dashed.solid, "dashed" )
                            arrows.head <- c ( arrows.head, "none" )         
                        }
                    }
                }
                # make perfect matches green  
                if ( count.markers == length ( unlist ( marker.list ) ) ) {        
                    labels.colour[q15] <- "lightgreen"
                    count.markers <- 0
                } else  { # colour the partial matches a beige colour
                    labels.colour[q15] <- "bisque"
                    count.markers <- 0     
                }
                break
            } else  { # colour all the non-important parents white
                labels.colour[q15] <- "white"
            }
        }
    }
    if ( length ( marker.list[["Positive"]] ) != 0 ) { # colour positive markers sky blue
        for ( q19 in 1:( length ( marker.list[["Positive"]] ) ) ) {
            labels.colour[length ( labels ) + q19] <- "skyblue"
        }
    }
    if ( length ( marker.list[["Negative"]] ) != 0 ) {
        for (q19 in 1:(length ( marker.list[["Negative"]] ) ) ) { # colour negative markers a light red colour
            labels.colour[length ( labels ) + length ( marker.list[["Positive"]] ) + q19] <- "lightcoral"
        }
    }
    nAttrs <- list ( )
    eAttrs <- list ( )

    nAttrs$fillcolor <- structure(c (labels.colour ), .Names = c (labels, marker.list.short[["Positive"]], marker.list.short[["Negative"]] ) )
    eAttrs$color <- structure(c (arrows.colour ), .Names = c (arrows.label ) )
  
#   # Even though it looks like arrows.dashed.solid and arrows.head are being used, they are not. 
#   # This is how the vignette for Rgraphiz implements these, however it does not work currently.
#   eAttrs$style <- structure(c ( arrows.dashed.solid ), .Names = c ( arrows.label ) )
#   eAttrs$arrowhead <- structure(c ( arrows.head ), .Names = c ( arrows.label ) )
    
    attrs <- list( node = list ( shape="ellipse", fontsize = 14, fixedsize=FALSE ), graph = list ( rankdir = "BT" ) )

    # create a pdf flow chart
    child.file.name <- paste ( save.dir, "tree_", phenotype, ".pdf", sep = "" )     
    pdf ( file = child.file.name, width = 10.5, height = 8 )    
    suppressWarnings ( plot ( rEG, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs ) )
    dev.off ( )

    return ( list ( rEG, nAttrs, eAttrs, attrs ) )
}

#################################################################################
### Break up a phenotype into individual markers and 'signs'
#   TO DO: After discussion this should be changed to also support expression
#   levels other than 'positive' and 'negative' -- such as 'dim', 'low', 'lo', 'high' and ??

phenoParse <- function ( phenotype )  {
    
    if ( !is.character ( phenotype ) ) {
        warning ( "Phenotype is not a valid string!" )
        try ( phenotype <- as.character ( phenotype ) )
    }
  
    # First split up the string based on + or - to get the markers
    markers <- unlist ( strsplit ( x = phenotype, split = "\\+|\\-" ) )
    # Next, split the original string based on the markers found above, leaving
    # only the signs (remove first result as there is no leading sign )
  
    signs <- unlist ( strsplit ( x = phenotype, split = paste ( markers, sep = "", collapse="|" ) ) )[-1]
  
    # Return a list of positive and negative markers
    res <- list ( `Positive` = markers[signs == "+"], `Negative` = markers[signs == "-"] )
  
    return ( res )
}

#################################################################################
# Creates a string with the updated phenotypes to have better formatting for listPhenotypes.csv
phenoUnparse <- function ( phenotype, marker.list ) {

posit.pos <- gregexpr ( pattern="[+]", phenotype )
posit.neg <- gregexpr ( pattern="[-]", phenotype )
posit.pos <- posit.pos[[1]]
posit.neg <- posit.neg[[1]]
temp <- c ( )
q1 <- 1; q2 <- 1
# create the vector temp containing 1's and 0's to record the order of negatives and positives.
for ( k in 1:( length ( posit.pos ) + length ( posit.neg ) - as.numeric ( posit.pos[q1] == - 1 ) - as.numeric ( posit.neg[q2] == - 1 ) ) ) {
    if ( posit.pos[q1] != - 1 && posit.neg[q2] != - 1 && !is.na ( posit.pos[q1] ) && !is.na ( posit.neg[q2] ) ) {
        if ( posit.pos[q1] < posit.neg[q2] ) {
            temp[k] <- 1; q1 <- q1 + 1 # pos occurred earlier therefore pos
        } else {
            temp[k] <- 0; q2 <- q2 + 1 # neg occurred earlier therefore neg 
        }
    } else { # this else statement only comes into effect if there is no more +'s or no -'s left in the lists.
        if ( posit.pos[q1] == - 1 || is.na ( posit.pos[q1] ) ) {
            temp[k] <- 0 # no more pos therefore neg
        }  
        if ( posit.neg[q2] == - 1 || is.na ( posit.neg[q2] ) ) {
            temp[k] <- 1 # no more neg therefore pos
        }  
    }
}

# add +'s and -'s onto marker.list according to the order of the 1's and 0's from temp
temp.string <- c ( )
q1 <- 1; q2 <- 1
for ( k in 1:length ( temp ) ) {
    if ( temp[k] == 1  )  {
        temp.string <-  paste ( temp.string, marker.list[["Positive"]][q1], "\n", sep = "" ) ; q1 <- q1 +1
    }
    if ( temp[k] == 0  )  {
        temp.string <-  paste ( temp.string, marker.list[["Negative"]][q2], "\n", sep = "" ) ; q2 <- q2 +1
    }
}
# removes the \n from the last line for formatting reasons in the .csv file
temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) 
return ( temp.string )
}

#################################################################################
### Condense a results table to the unique entries only and tabulate repeated hits
#   TO DO: For now, this relies on the first column having the unique IDs of
#   the owl objects returned.
tabulateResults <- function ( res ) {
    number.of.hits <- table ( res[ , 1] )
    res <- cbind ( res, number.of.hits[res[, 1]] )
    colnames ( res )[ncol ( res )] <- "Number Of Hits"
    unique.ids <- unique ( res[, 1] )
    clean.res <- matrix ( "", nrow = length ( unique.ids ), ncol = ncol ( res ) + 1 )
    colnames ( clean.res ) <- c ( colnames ( res ), "Score" )
    rownames ( clean.res ) <- unique.ids
    for ( id in unique.ids ) {
        locate.entries <- which ( res[, 1] == id )
        clean.res[id, ] <- c ( id, res[locate.entries[1], 2], 
                                apply(res[locate.entries, 3:( ncol ( res ) - 2 )], 2, paste, 
                                collapse = "\n" ), sum ( res[locate.entries, 'penalties'] ),
                                res[locate.entries[1], ncol(res )], 0 )
        clean.res[id, "Score"] <- as.numeric ( clean.res[id, "Number Of Hits"] ) + 
        as.numeric ( clean.res[id, "penalties"] )
    }
    # Required because of a bug cause from only having 1 row (f[1,] notation was a problem )
    if ( nrow ( clean.res ) >= 2 ) {
        sort.scores <- sort ( as.numeric ( clean.res[, "Score"] ), 
                                decreasing = TRUE, index.return = TRUE )$ix
        return ( clean.res[sort.scores, ] )
    } else  {
        return ( clean.res )
    }
}

#################################################################################
# Prints out the time since start_time. Used for optimizing code and for informing the user how long certain processes take.
timeOutput <- function ( start_time )  {
    start_time <- as.POSIXct ( start_time )
    dt <- difftime ( Sys.time ( ), start_time, units = "secs" )
    # Since you only want the H:M:S, we can ignore the date...
    # but you have to be careful about time-zone issues
    format ( .POSIXct ( dt, tz = "GMT" ), "%H:%M:%S" )
}
timeOutput ( Sys.Date ( ) )

#################################################################################
# Finds and stores information for display in a .csv file
updateLists <- function ( clean.res, MaxHitsPht ) {
    # Creates the lists MarkerLabels, CellLabels, PhenotypeID and CellID
    BreakTrue <- FALSE
    if ( max ( clean.res[ , 'Number Of Hits'] ) >= 1 ) {
        for ( q2 in 1:min ( length ( clean.res[ ,'Number Of Hits'] ), MaxHitsPht ) ) {
            if ( q2 == 1 & ( clean.res[q2, 'Number Of Hits'] ) == max ( clean.res[ , 'Number Of Hits'] ) ) {
                listMarkerLabels.temp <- paste ( q2, ") ", ( clean.res[q2,'markerlabel'] ), "\n", sep = "" ) 
                listCellLabels.temp   <- paste ( q2, ") ", ( clean.res[q2,'celllabel'] ), "\n", sep = "" )
                # Look in the Label of the marker ID for "PR" then pulls out the PR and the numbers that follow 
                temp.index <- gregexpr ( "PR", ( clean.res[q2,'marker'] ) )
                temp.string <- ""
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    for (q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr(clean.res[q2,'marker'],temp.index[[1]][q7],temp.index[[1]][q7]+11 ),"\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) {temp.string <- substr(temp.string,1,nchar ( temp.string )-1 )} # Removes the last \n 
                    }
                }
                # Look in the Label of the marker ID for "GO" then pulls out the GO and the numbers that follow
                temp.index <- gregexpr ( "GO", ( clean.res[q2,'marker'] ) )
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    if ( temp.string != "" ) {
                        temp.string <- paste ( temp.string, "\n", sep = "" )
                    }
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'marker'], temp.index[[1]][q7], temp.index[[1]][q7] + 9 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n      
                    }
                }
                # Stores this string into the lists
                listPhenotypeID.temp  <- paste ( q2, ") ", temp.string, "\n", sep = "" )
        
                # Look in the Label of the Cell ID for "CL" then pulls out the CL and the numbers that follow
                temp.index <- gregexpr ( "CL", ( clean.res[q2,'x'] ) )
                temp.string <- ""
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'x'], temp.index[[1]][q7], temp.index[[1]][q7] + 9 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n 
                    }
                }
                # Stores this string into the lists
                listCellID.temp  <- paste ( q2,") ", temp.string, "\n", sep = "" )
            }
            if ( q2 > 1 & ( clean.res[q2,'Number Of Hits'] ) == max ( clean.res[,'Number Of Hits'] ) ) {
                listMarkerLabels.temp <- paste ( listMarkerLabels.temp , q2, ") ", ( clean.res[q2,'markerlabel'] ), "\n", sep = "" )
                listCellLabels.temp   <- paste ( listCellLabels.temp   , q2, ") ", ( clean.res[q2,'celllabel'] ), "\n", sep = "" )
        
                # Look in the Label of the marker ID for "PR" then pulls out the PR and the numbers that follow 
                temp.index <- gregexpr ( "PR", ( clean.res[q2,'marker'] ) )
                temp.string <- ""
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'marker'], temp.index[[1]][q7], temp.index[[1]][q7] + 11 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n 
                    }
                }
                # Look in the Label of the marker ID for "GO" then pulls out the GO and the numbers that follow
                temp.index <- gregexpr ( "GO", ( clean.res[q2,'marker'] ) )
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    if ( temp.string != "" ) {
                        temp.string <- paste ( temp.string, "\n", sep = "" )
                    }
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'marker'], temp.index[[1]][q7], temp.index[[1]][q7] + 9 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n        
                    }
                }
                # Stores this string into the lists
                listPhenotypeID.temp  <- paste ( listPhenotypeID.temp  , q2, ") ", temp.string,"\n", sep = "" )
        
                # Look in the Label of the Cell ID for "CL" then pulls out the CL and the numbers that follow
                temp.index <- gregexpr ( "CL", ( clean.res[q2,'x'] ) )
                temp.string <- ""
                if ( ( temp.index[[1]][1] != - 1 ) ) {
                    for ( q7 in 1:length ( temp.index[[1]] ) ) {
                        temp.string <- paste ( temp.string, substr ( clean.res[q2,'x'], temp.index[[1]][q7], temp.index[[1]][q7] + 9 ), "\n", sep = "" )
                        if ( q7 == length ( temp.index[[1]] ) ) { temp.string <- substr ( temp.string, 1, nchar ( temp.string ) - 1 ) } # Removes the last \n 
                    }
                }
                # Stores this string into the lists
                listCellID.temp  <- paste ( listCellID.temp  , q2, ") ", temp.string, "\n", sep = "" )
            }
            if ( ( clean.res[q2,'Number Of Hits'] ) != max ( clean.res[ , 'Number Of Hits'] ) ) {
                BreakTrue <- TRUE
                break
            }  
        }# end of for loop
    }# end of if statement
  
    # If there is more than 5 elements to store, then the rest is cut off and a "+ more" is displayed
    if ( length ( clean.res[ , 'Number Of Hits'] ) > MaxHitsPht & BreakTrue == FALSE ) {
        listMarkerLabels.temp <- paste ( listMarkerLabels.temp, "+ more" )
        listCellLabels.temp   <- paste ( listCellLabels.temp,   "+ more" )
        listPhenotypeID.temp  <- paste ( listPhenotypeID.temp,  "+ more" )
        listCellID.temp       <- paste ( listCellID.temp,       "+ more" )
    } else  {  # Removes the last \n from the string for better formatting int the .csv file
        listMarkerLabels.temp <- substr ( listMarkerLabels.temp, 1, nchar ( listMarkerLabels.temp ) - 1 )
        listCellLabels.temp   <- substr ( listCellLabels.temp,   1, nchar ( listCellLabels.temp )   - 1 )
        listPhenotypeID.temp  <- substr ( listPhenotypeID.temp,  1, nchar ( listPhenotypeID.temp )  - 1 )
        listCellID.temp       <- substr ( listCellID.temp,       1, nchar ( listCellID.temp )       - 1 )
    }
  
    return ( c ( listMarkerLabels.temp, listCellLabels.temp, listPhenotypeID.temp, listCellID.temp ) )
}
###############################################################3
# List of phenotypes commonly used by HIPC
listPhenotypes_flowCL <- function ( ) {

return <- c(
"CD3+","CD4+","CD8+","CCR7+","CD45RA+","CD38+","HLA.DR+","CD127+","CD25+","CCR4+","CD45RO+","CXCR3+",
"CCR6+","CD19+","CD20+","CD27+","IgD+","CD24+","CD14+","CD11c+","CD123+","CD16+","CD56+",
"CD3-","CD4-","CD8-","CCR7-","CD45RA-","CD38-","HLA.DR-","CD127-","CD25-","CCR4-","CD45RO-","CXCR3-",
"CCR6-","CD19-","CD20-","CD27-","IgD-","CD24-","CD14-","CD11c-","CD123-","CD16-","CD56-",
"CD3-CD19+CD20+IgD+CD27-",
"CD3-CD19+CD20+IgD-CD27+",
"CD3-CD19+CD20+IgD+CD27+",
"CD3-CD19+CD20+CD38+CD24+",
"CD3-CD19+CD20-CD38+CD27+",
"CD3+CD4+CD8-CCR7+CD45RA+",
"CD3+CD4+CD8-CCR7+CD45RA-",
"CD3+CD4+CD8-CCR7-CD45RA-",
"CD3+CD4+CD8-CCR7-CD45RA+",
"CD3+CD4+CD8-CD38+HLA.DR+",
"CD3+CD4-CD8+CCR7+CD45RA+",
"CD3+CD4-CD8+CCR7+CD45RA-",
"CD3+CD4-CD8+CCR7-CD45RA-",
"CD3+CD4-CD8+CCR7-CD45RA+",
"CD3+CD4-CD8+CD38+HLA.DR+",
"CD3+CD4+CD25+CD127-CCR4+CD45RO+",
"CD3+CD4+CD25+CD127-CCR4+CD45RO-",
"CD3+CD4+CD25+CD127-CCR4+HLA.DR+",
"CD3+CD4+CD8-CXCR3-CCR6-",
"CD3+CD4+CD8-CXCR3+CCR6-",
"CD3+CD4+CD8-CXCR3-CCR6+",
"CD3+CD4-CD8+CXCR3-CCR6-",
"CD3+CD4-CD8+CXCR3+CCR6-",
"CD3+CD4-CD8+CXCR3-CCR6+",
"CD14-CD3-CD19-CD20-CD16-CD56-HLA.DR+CD11c-CD123+",
"CD14-CD3-CD19-CD20-CD16-CD56-HLA.DR+CD11c+CD123-",
"CD3+CD4+CD127-CD25+",
"CD3+CD4+CD127-CD25+CCR4+",
"CD3-CD19-CD20-CD14+",
"CD3-CD19-CD20-CD14+CD16+",
"CD3-CD19-CD20-CD14+CD16-",
"CD3-CD19-CD20-CD14-HLA.DR-",
"CD3-CD19-CD20-CD14-HLA.DR-CD16+CD56+",
"CD3-CD19-CD20-CD14-HLA.DR-CD16+CD56-",
"CD3-CD19-CD20-CD14-CD16-CD56-HLA.DR+")

}
####################################################
# List of colours used for the tree diagram
listColours_flowCL <- function ( ) {
    return <- c(
"red", "blue", "forestgreen", "darkviolet", "gold2", "darkorange3", "aquamarine4", "aquamarine2", "khaki3", "deeppink", "tan4", 
"darkmagenta", "coral2", "burlywood4", "hotpink", "lightblue2", "lightpink4", "mediumseagreen", "navyblue", "olivedrab1", 
"orangered", "purple", "peru", "plum2", "tan3", "bisque1", "darkslategray3", "darkslategray1", 
"gray23", "gray24", "gray25", "gray26", "gray27", "gray28", "gray29", "gray30", "gray31", "gray32", "gray33", "gray34", 
"gray35", "gray36", "gray37", "gray38", "gray39", "gray40", "gray41", "gray42", "gray43", "gray44", "gray45", "gray46", 
"gray47", "gray48", "gray49", "gray50", "gray51", "gray52", "gray53", "gray54", "gray55", "gray56", "gray57", "gray58", 
"gray59", "gray60", "gray61", "gray62", "gray63", "gray64", "gray65", "gray66", "gray67", "gray68", "gray69", "gray70", 
"gray71", "gray72", "gray73", "gray74", "gray75", "gray76", "gray77", "gray78")
}
##################################
# function for loading hasProperLabel data
flowCL_query_data_hasProperLabel <- function(){
return <- c("select distinct ?x ?label",
"{",
"?x a owl:Class.",                         
"?x rdfs:label ?label. ",                  
"FILTER regex(?label, \"$marker\", \"i\")",
"}")                                       
}
##################################
# function for loading prefix.info data
flowCL_query_data_prefix.info <- function(){

return <- c("# Common prefix and abbreviation",                                                       
"prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>",                               
"prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>",                                    
"prefix owl: <http://www.w3.org/2002/07/owl#>",                                            
"prefix obo: <http://purl.obolibrary.org/obo/>",                                           
"prefix oboinowl: <http://www.geneontology.org/formats/oboInOwl#>",                        
"prefix probe: <http://purl.obolibrary.org/obo/CL_0000903>",                               
"# Has plasma membrane part (has_pmp)",                                                    
"prefix has_pmp: <http://purl.obolibrary.org/obo/RO_0002104>",                             
"prefix lacks_pmp: <http://purl.obolibrary.org/obo/cl#lacks_plasma_membrane_part>",        
"prefix has_high_pma: <http://purl.obolibrary.org/obo/cl#has_high_plasma_membrane_amount>",
"prefix has_low_pma: <http://purl.obolibrary.org/obo/cl#has_low_plasma_membrane_amount>",  
"prefix definition: <http://purl.obolibrary.org/obo/IAO_0000115>")                         
}

##################################
# function for loading hasProperSynonym data
flowCL_query_data_hasProperSynonym <- function(){
    
return <- c("select distinct ?x ?label ?synonym ",
"where",                                     
"{",                                         
"?x a owl:Class.",                           
"?x rdfs:label ?label.",                     
"?x oboinowl:hasExactSynonym ?synonym. ",    
"FILTER regex(?synonym, \"$marker\", \"i\")",
"}")                                         
}

##################################
# function for loading getParentClasses data
flowCL_query_data_getParentClasses <- function(){
    
return <- c("# Find all parent classes of the cell type of interest. Note that for some reason,",
"# matching on ?x does not work, but matching on ?celllabel (x's label) does.",      
"# Matching directly on ?x works on http://sparql.hegroup.org/sparql !",             
"select distinct ?x ?celllabel ?parent ?parentlabel",                                
"where",                                                                             
"{",                                                                                 
"  ?parent a owl:Class.",                                                            
"  ?x a owl:Class.",                                                                 
"  ?x rdfs:label ?celllabel.",                                                       
"  ?x rdfs:subClassOf ?parent.",                                                     
"  ?parent rdfs:label ?parentlabel.",                                                
"  FILTER regex(?celllabel, \"$label\", \"i\")",                                     
"}")                                                                                 
}

##################################
# function for loading hasPlasmaMembranePart data
flowCL_query_data_hasPlasmaMembranePart <- function(){
    
return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",                                                     
"{",                                                         
"  ?x a owl:Class.",                                         
"  ?x rdfs:label ?celllabel.",                               
"  ?x rdfs:subClassOf ?sub.",                                
"  ?sub rdf:type owl:Restriction.",                          
"  ?sub owl:onProperty has_pmp:.",                           
"  ?sub owl:someValuesFrom ?marker.",                        
"  ?marker rdfs:label ?markerlabel.  ",                      
"  has_pmp: rdfs:label ?plabel.",                            
"  FILTER regex(?markerlabel, \"$marker\", \"i\")",          
"}")                                                         
}

##################################
# function for loading lacksPlasmaMembranePart data
flowCL_query_data_lacksPlasmaMembranePart <- function(){
    
return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",                                                     
"{",                                                         
"  ?x a owl:Class.",                                         
"  ?x rdfs:label ?celllabel.",                               
"  ?x rdfs:subClassOf ?sub.",                                
"  ?sub rdf:type owl:Restriction.",                          
"  ?sub owl:onProperty lacks_pmp:.",                         
"  ?sub owl:someValuesFrom ?marker.",                        
"  ?marker rdfs:label ?markerlabel. ",                       
"  lacks_pmp: rdfs:label ?plabel.",                          
"  FILTER regex(?markerlabel, \"$marker\", \"i\")",          
"}")                                                         
}

##################################
# function for loading hasPMPsingle data
flowCL_query_data_hasPMPsingle <- function(){
    
return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",                                                     
"{",                                                         
"  ?x a owl:Class.",                                         
"  ?x rdfs:label ?celllabel.",                               
"  ?x rdfs:subClassOf ?sub.",                                
"  ?sub rdf:type owl:Restriction.",                          
"  ?sub owl:onProperty has_pmp:.",                           
"  ?sub owl:someValuesFrom ?marker.",                        
"  ?marker rdfs:label ?markerlabel. ",                       
"  has_pmp: rdfs:label ?plabel.",                            
"  FILTER regex(?markerlabel, \"$marker\", \"i\")",          
"  FILTER regex(?celllabel, \"$celllabel\", \"i\")",         
"}")                                                         
}

##################################
# function for loading lacksPMPsingle data
flowCL_query_data_lacksPMPsingle <- function(){
    
return <- c("select distinct ?x ?celllabel ?plabel ?marker ?markerlabel",
"where",                                                     
"{",                                                         
"  ?x a owl:Class.",                                         
"  ?x rdfs:label ?celllabel.",                               
"  ?x rdfs:subClassOf ?sub.",                                
"  ?sub rdf:type owl:Restriction.",                         
"  ?sub owl:onProperty lacks_pmp:.",                         
"  ?sub owl:someValuesFrom ?marker.",                        
"  ?marker rdfs:label ?markerlabel. ",                       
"  lacks_pmp: rdfs:label ?plabel.",                          
"  FILTER regex(?markerlabel, \"$marker\", \"i\")",          
"  FILTER regex(?celllabel, \"$celllabel\", \"i\")",         
"}")  
}

##################################
# function for loading lacksPMPsingle data
flowCL_query_date <- function(){
    
return <- c("PREFIX :<http://purl.obolibrary.org/obo/cl.owl#>",
"PREFIX geo-pos:<http://www.w3.org/2003/01/geo/wgs84_pos#>",
"PREFIX uberon:<http://purl.obolibrary.org/obo/uberon#>",
"PREFIX umbel-ac:<http://umbel.org/umbel/ac/>",
"PREFIX sw-vocab:<http://www.w3.org/2003/06/sw-vocab-status/ns#>",
"PREFIX ff:<http://factforge.net/>",
"PREFIX music-ont:<http://purl.org/ontology/mo/>",
"PREFIX dc-term:<http://purl.org/dc/terms/>",
"PREFIX om:<http://www.ontotext.com/owlim/>",
"PREFIX opencyc-en:<http://sw.opencyc.org/2008/06/10/concept/en/>",
"PREFIX factbook:<http://www.daml.org/2001/12/factbook/factbook-ont#>",
"PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>",
"PREFIX pext:<http://proton.semanticweb.org/protonext#>",
"PREFIX ot:<http://www.ontotext.com/>",
"PREFIX dc:<http://purl.org/dc/elements/1.1/>",
"PREFIX onto:<http://www.ontotext.com/>",
"PREFIX foaf:<http://xmlns.com/foaf/0.1/>",
"PREFIX yago:<http://mpii.de/yago/resource/>",
"PREFIX umbel:<http://umbel.org/umbel#>",
"PREFIX pkm:<http://proton.semanticweb.org/protonkm#>",
"PREFIX wordnet16:<http://xmlns.com/wordnet/1.6/>",
"PREFIX owl:<http://www.w3.org/2002/07/owl#>",
"PREFIX gr:<http://purl.org/goodrelations/v1#>",
"PREFIX wordnet:<http://www.w3.org/2006/03/wn/wn20/instances/>",
"PREFIX opencyc:<http://sw.opencyc.org/concept/>",
"PREFIX wordn-sc:<http://www.w3.org/2006/03/wn/wn20/schema/>",
"PREFIX nytimes:<http://data.nytimes.com/>",
"PREFIX dbp-prop:<http://dbpedia.org/property/>",
"PREFIX geonames:<http://sws.geonames.org/>",
"PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>",
"PREFIX dbpedia:<http://dbpedia.org/resource/>",
"PREFIX oasis:<http://psi.oasis-open.org/iso/639/#>",
"PREFIX geo-ont:<http://www.geonames.org/ontology#>",
"PREFIX umbel-en:<http://umbel.org/umbel/ne/wikipedia/>",
"PREFIX ubprop:<http://purl.obolibrary.org/obo/ubprop#>",
"PREFIX bbc-pont:<http://purl.org/ontology/po/>",
"PREFIX ptop:<http://proton.semanticweb.org/protontop#>",
"PREFIX lingvoj:<http://www.lingvoj.org/ontology#>",
"PREFIX fb:<http://rdf.freebase.com/ns/>",
"PREFIX dbtune:<http://dbtune.org/bbc/peel/work/>",
"PREFIX obo:<http://purl.obolibrary.org/obo/>",
"PREFIX psys:<http://proton.semanticweb.org/protonsys#>",
"PREFIX umbel-sc:<http://umbel.org/umbel/sc/>",
"PREFIX dbp-ont:<http://dbpedia.org/ontology/>",
"PREFIX xsd:<http://www.w3.org/2001/XMLSchema#>",
"PREFIX ub:<http://www.lehigh.edu/~zhp2/2004/0401/univ-bench.owl#>",
"PREFIX oboInOwl:<http://www.geneontology.org/formats/oboInOwl#>",
"PREFIX skos:<http://www.w3.org/2004/02/skos/core#>",

"select ?s  ?p",
"where",
"{",
"?s owl:versionIRI ?p.",
"}"
)  
}

test.flowCL.connection <- function()
{   
    require("RUnit")
    testsuite <- defineTestSuite("flowCL.check",
                    dirs = system.file("unitTests", package="flowCL"),
                    testFileRegexp = "^test_.*\\.R$", testFuncRegexp = "^test.+")
    testResult <- runTestSuite(testsuite)
    printTextProtocol(testResult)
}

