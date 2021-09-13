"""

	pbc1d!(abq::AbqModel)

Generate periodic boundary conditions for onedimensional periodicity and
append them to the given AbqModel.
"""
function pbc1D!(abq::AbqModel)
	# abbreviation for coordinate system
	c = abq.csys
	# abbreviation for unit cell dimensions
	l = abq.dim
	# abbreviation for minimum coordinates of unit cell
	min_c = abq.minC
	# initiate counter for equation numbering
	i=1
	################
	# VERTEX NODES #
	################
	# no equations for vertices needed!
	##############
	# EDGE NODES #
	##############
	for n = 1:length(abq.edges["E"])
		this_node = abq.edges["E"][n]
		x = this_node.node.coords
		for a = 1:2
			i+=1
			push!(abq.eqns,Equation(i,["E-$(n)","W-$(n)","SE","SW","NE","NW"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [+1.0,
									   -1.0,
									   -1.0 + (x[c[2]] - min_c[c[2]])/l[c[2]],
									   +1.0 - (x[c[2]] - min_c[c[2]])/l[c[2]],
									   -(x[c[2]] - min_c[c[2]])/l[c[2]],
									   +(x[c[2]] - min_c[c[2]])/l[c[2]]
									  ]))
		end
	end
	println("$(i) equations written to AbqModel.")
	return
end
