#=
	Constant for defining location of vertices
	Each array specifies the positioning of the
	in terms of relative I-II-coordinates.
=#
const vertices = Dict("SW"=>[0,0],
					  "NW"=>[0,1],
					  "SE"=>[1,0],
					  "NE"=>[1,1])

#=
	Constant for defining location of edges
	Each array specifies the positioning of the
	in terms of relative I-II-coordinates.
	The first number indicates the axis in which
	the edge is oriented, the remaining numbers
	define the coordinates of the edge normal to
	the orientation axis.
=#
const edges = Dict("S"=>[1,0],
				   "N"=>[1,1],
				   "W"=>[2,0],
				   "E"=>[2,1])
#=
	Define which coordinates are defined by
	the second number of the arrays
	in edges.
=#
const edgeCS = [2;1]

"""

	isEqual(num1::Number,num2::Number,tolerance::Number)

"""
function isEqual(num1::Number,num2::Number,tolerance::Number)
	val = num1 > num2 - tolerance &&
		  num1 < num2 + tolerance
	return val
end

"""

	checkNum(d::Dict,k1::AbstractString,k2::AbstractString)

"""
function checkNum(d::Dict,k1::AbstractString, k2::AbstractString)
	if length(d[k1]) != length(d[k2])
		println("Node numbers of $k1 and $k2 are not equal!")
	end
	return
end

"""

	sortNodes(nodes::Array{GlobNode,1})

Function for sorting of an array containing nodes, according to their coordinates.
Sorting considers x-coordinates first, y-coordinates second, and z-coordinates last.
"""
function sortNodes(nodes::Array{GlobNode,1})
	# Array with weighting factors for sorting of nodes
	weights = [1e5 1]
	# Initiate array, soon-to-hold the values decisive for sorting
	sortVals = Array{Float64,1}()
	# Fill array sortVals with values, calculated by coordinates
	for n in nodes
		sV = weights * n.node.coords
		append!(sortVals,sV)
	end
	# Obtain array containing the order of the sorted elements
	sortOrder = sortperm(sortVals)
	# Return sorted array of nodes
	return nodes[sortOrder]
end


"""

	findVertex(abq::AbqModel, name::AbstractString)

"""
function findVertex(abq::AbqModel, name::AbstractString)
	nodeBool = map(abq.nodes) do n
		!(n.instance in abq.ecc) && # Is instance containing node in ecceptions?
		isEqual(n.node.coords[abq.csys[1]], abq.minC[abq.csys[1]] + vertices[name][1] * abq.dim[abq.csys[1]], abq.tol) &&
		isEqual(n.node.coords[abq.csys[2]], abq.minC[abq.csys[2]] + vertices[name][2] * abq.dim[abq.csys[2]], abq.tol)
	end
	if length(abq.nodes[nodeBool]) != 1
		throw(ToleranceError)
	end
	return abq.nodes[nodeBool][1]
end

"""

	findVertices!(abq::AbqModel)

"""
function findVertices!(abq::AbqModel)
	for v in keys(vertices)
		abq.vertices[v] = findVertex(abq, v)
	end
	println("Vertices written to AbqModel.")
	return
end

"""

	findEdge(abq::AbqModel, name::AbstractString)

"""
function findEdge(abq::AbqModel, name::AbstractString)
	axis1 = abq.csys[edges[name][1]]
	axis2 = abq.csys[edgeCS[edges[name][1]]]
	nodeBool = map(abq.nodes) do n
		!(n.instance in abq.ecc) && # Is instance containing node in ecceptions?
		n.node.coords[axis1] > abq.minC[axis1] + abq.tol &&
		n.node.coords[axis1] < abq.minC[axis1] + abq.dim[axis1] - abq.tol &&
		isEqual(n.node.coords[axis2], abq.minC[axis2] + edges[name][2]*abq.dim[axis2], abq.tol)
	end
	return abq.nodes[nodeBool]
end

"""

	findEdges!(abq::AbqModel)

"""
function findEdges!(abq::AbqModel)
	for e in keys(edges)
		abq.edges[e] = sortNodes(findEdge(abq, e))
	end
	# Check if corresponding edges have equal node numbers
	checkNum(abq.edges,"S","N")
	checkNum(abq.edges,"W","E")
	println("Edges written to AbqModel.")
	return
end

"""

	nset(name::AbstractString, node::Int, inst::AbstractString)

Generate an array containing the definition of a node set in a Abaqus input file.
"""
function nset(name::AbstractString, node::Int, inst::AbstractString)
	return ["*Nset, nset=$(name), instance=$(inst)",
			" $(node),"]
end


"""

	nodeDesignation!(abq::AbqModel)

"""
function nodeDesignation!(abq::AbqModel)
	abq.defRA ? println("Reference axis not explicitly set. Default $(abq.refAxis) is used.") : print("")
	findVertices!(abq)
	findEdges!(abq)
	println("Node designation written to AbqModel.")
	return
end


