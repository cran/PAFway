#' Compare the frequency of observed edges of a specific type and the expected frequency of these edges, given the presence of a different edge type
#'
#' @param GO A vector of Strings, equal to the length of all the nodes.  The names of the vector should be the names of the nodes.  The values should either be the functional annotations, or a concatenation of the functional annotations, separated by a "_" symbol.
#' @param edges A matrix of Strings, with at least two columns.  Each row will represent an edge, linking the node in the first column to the node in the second column.  Please make sure the node names are the same as those in "GO"
#' @param GOtypes This is a vector that contains the functional annotations or GO terms that are of interest
#' @return A matrix that has the same number of rows and columns as length(GOtypes)*length(GOtypes).  The value at position [i,j] will contain ratio between the observed and expected frequency of edges of type i, based on the frequency of edges of type j.
#' @export
#' @examples
#' GO=c(rep('A,C', 5), rep('A', 5), rep('C', 5), rep('B,D', 5), rep('B', 5), rep('D', 5))
#' names(GO)=paste("node", 1:length(GO))
#' edges=cbind(names(GO)[1:(length(GO)/2)], names(GO)[(length(GO)/2+1):length(GO)])
#' GOtypes=c('A', "B", "C", "D")
#' pafway_meta(GO, edges, GOtypes)
pafway_meta <- function(GO, edges, GOtypes) {
    
    GOinNetwork = GO[unique(c(edges[, 1], edges[, 2]))]
    
    freq_together=sapply(GOtypes, function(i) { #how often a gene has both terms
        sapply(GOtypes, function(j) {
            length(which(grepl(i, GOinNetwork) & grepl(j, GOinNetwork)))
        })
    })
    
    freq_from=sapply(GOtypes, function(i){
        length(which(grepl(i, GOinNetwork[edges[, 1]])))
    })
    
    freq_to=sapply(GOtypes, function(i){
        length(which(grepl(i, GOinNetwork[edges[, 2]])))
    })
    
    freq_mat=t(sapply(GOtypes, function(i) { #how often an edge is found
        sapply(GOtypes, function(j) {
            length(which(grepl(i, GOinNetwork[edges[, 1]]) & grepl(j, GOinNetwork[edges[, 2]])))
        })
    }))
    
    possible_edges=strsplit(sapply(GOtypes, function(i) { 
        sapply(GOtypes, function(j) { 
            paste(i,j,sep=",")
            })}), ",")
    

    
    names_possible_edges=sapply(possible_edges, function(i){paste(i[1], i[2], sep="->")})
    
    expected_freq_mat=sapply(possible_edges, function(i){
        sapply(possible_edges, function(j){
            
            a=i[1]
            b=i[2]
            c=j[1]
            d=j[2]
            
            a_to_b=freq_mat[a,b]
            a_to_nb=(freq_from[a]-freq_mat[a,b])
            na_to_b=(freq_to[b]-freq_mat[a,b])
            na_to_nb=(length(edges[,1])-freq_from[a]-freq_to[b]+freq_mat[a,b])
            
            p_a_and_c=freq_together[a,c]/freq_together[a,a]
            p_na_and_c=(freq_together[c,c]-freq_together[a,c])/(length(GOinNetwork)-freq_together[a,a])
            p_b_and_d=freq_together[b,d]/freq_together[b,b]
            p_nb_and_d=(freq_together[d,d]-freq_together[b,d])/(length(GOinNetwork)-freq_together[b,b])
            
            a_to_b*p_a_and_c*p_b_and_d+
                a_to_nb*p_a_and_c*p_nb_and_d+
                   na_to_b*p_na_and_c*p_b_and_d+
                       na_to_nb*p_na_and_c*p_nb_and_d
        })
    })

    rownames(expected_freq_mat)= names_possible_edges
    colnames(expected_freq_mat)= names_possible_edges
    true_freq_mat=sapply(possible_edges, function(i){
        sapply(possible_edges, function(j){
            c=j[1]
            d=j[2]
            freq_mat[c,d]
        })
    })
    true_freq_mat/expected_freq_mat
            
            }

        