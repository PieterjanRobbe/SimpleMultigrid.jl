#
# UTILS.JL
#
# Utility functions for Multigrid

# "discrete" L2 norm
discr_norm(u,h,d) = sqrt(h^d*sum(u.^2))

# interpolation
function interpolate(x, f, xp)
	x_ = Tuple([0:x[i][1]:1 for i = 1:length(x)])
	sizes = length.(x)
	f_ = zeros(sizes.+2)
	ranges = Tuple([2:sizes[i]+1 for i = 1:length(sizes)])
	f_[ranges...] = reshape(f,sizes)	
	ip = interpolate(f_, BSpline(Cubic(Line())), OnGrid())
	sip = scale(ip, x_...)

	xp_ = zeros(length.(xp))
	for k in Base.product(Tuple([1:length(xp[i]) for i=1:length(xp)])...)
		xp_[k...] = sip[[xp[l][k[l]] for l = 1:length(xp)]...]
	end

	return xp_[:]	
end