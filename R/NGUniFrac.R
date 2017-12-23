NGUniFrac_cum <- function (otu.tab, tree) {
	# Calculate Generalized UniFrac distances. Unweighted and
	# Variance-adjusted UniFrac distances will also be returned.
	#
	# Args:
	#		otu.tab: OTU count table, row - n sample, column - q OTU
	#		tree: rooted phylogenetic tree of R class "phylo"
	#		alpha: parameter controlling weight on abundant lineages
	#
	# Returns:
	# 	unifracs: three dimensional array containing the generalized
	#							UniFrac distances, unweighted UniFrac distance and
	#							variance adjusted UniFrac distances.
	#
	if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")

	# Convert into proportions
	otu.tab <- as.matrix(otu.tab)
	#row.sum <- rowSums(otu.tab)
	#otu.tab <- otu.tab / row.sum

	####ADD@2017.12.23
	library(metagenomeSeq)
	otu.tab <- cumNormMat(t(otu.tab), p = cumNormStatFast(t(otu.tab)), sl = 1000)
	otu.tab <- t(otu.tab)
	####

	n <- nrow(otu.tab)

	# Construct the returning array
	if (is.null(rownames(otu.tab))) {
		rownames(otu.tab) <- paste("comm", 1:n, sep="_")
	}

	# Check OTU name consistency
	if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
		stop("The OTU table contains unknown OTUs! OTU names
					in the OTU table and the tree should match!" )
	}

	# Get the subtree if tree contains more OTUs
	absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
	if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
		warning("The tree has more OTU than the OTU table!")
	}

	# Reorder the otu.tab matrix if the OTU orders are different
	tip.label <- tree$tip.label
	otu.tab <- otu.tab[, tip.label]


	ntip <- length(tip.label)
	nbr <- nrow(tree$edge)
	edge <- tree$edge
	edge2 <- edge[, 2]
  br.len <- tree$edge.length  # branch length


    #  Accumulate OTU proportions up the tree
	cum <- matrix(0, nbr, n)							# Branch abundance matrix
	for (i in 1:ntip) {
		tip.loc <- which(edge2 == i)
		cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
		node <- edge[tip.loc, 1]						# Assume the direction of edge
		node.loc <- which(edge2 == node)
		while (length(node.loc)) {
			cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
			node <- edge[node.loc, 1]
			node.loc <- which(edge2 == node)
		}
	}

	out = list(cum = cum, br.len= br.len)
	return(out)
}






