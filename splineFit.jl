# SplineFit.jl v0.1 
# 2022 Andreas Carlberger

using Plots
gr()



"""

	splineCoeffs = splineFit(x, y, nSplines ; x1=x[1], xn=x[end], curved_ends=true, doplot=false)

	splineCoeffs, yp = splineFit(x, y, nSplines, xp... ; x1=x[1], xn=x[end], curved_ends=true, doplot=false)

Function approximation with descrete data using cubic splines and least squares fit. Option to use the splines for interpolaton. Use `evalSpline`(`splineCoeffs`,`xp`) for futher use of the spline functions. 

`x`: 			Vector with data x-coordinates.  

`y`: 			Vector with data y-coordinates.

`nSplines`: 	Number of equally spaced cubic splines. 

`xp`: 			Optional: x-coordinates to be interpolated and returned.  

`x1`: 			First spline start point, default = `x`[1].

`xn`: 			Last spline end point, default = `x[end].

`curved_ends`: 	If the splines can have non-zero curvature at `x1` and `xn`, default = `true`.

`doplot`: 		plot results.

`splineCoeffs`: coefficients of spline functions.

`yp`: 			interpolated results.

# Example
```julia-repl
julia> x = collect(1:100);
julia> y = 3*log.(x)+0.2*randn(length(x))+2*sin.(x/10);
julia> nSplines = 6;
julia> xp = collect(-0:3:100);
julia> splineCoeffs,yp = splineFit(x, y, nSplines, xp, curved_ends=true, doplot=true);
julia> xp2 = [5;15];
julia> yp2 = evalSpline(splineCoeffs,xp2)
2-element Vector{Float64}:
  5.476211419705948
 10.17353202236643
julia> scatter!(xp2,yp2, mc=:orange,ms=7)
```
"""
function splineFit(x, y, nSplines, xp... ; x1=x[1], xn=x[end], curved_ends=true, doplot=false)
	allow_curvature_at_ends = curved_ends
	Δx = (xn-x1)/nSplines

	SS = createSS(nSplines, Δx, x1, allow_curvature_at_ends)
	splineID, xi = splineID2(x, x1, xn, Δx, nSplines)
	eqs = createEQS(x, y, splineID, xi, SS, allow_curvature_at_ends)
	splineCoeffs = solveEQS(eqs, y, nSplines, Δx, xi, allow_curvature_at_ends)
	
	if doplot
		plotSplines(x, y, splineCoeffs, xi, Δx)
	end
	if ~(xp == ())
		yp = evalSpline(splineCoeffs,xp[1])
		if doplot
			scatter!(xp[1],yp, markercolor=:magenta,ms=4)
			gui()
		end
		return splineCoeffs, yp
	end
	return splineCoeffs
end

mutable struct SplineCoeffs
	a
	b
	c
	d
	x
end
function solveEQS(eqs,y,nSplines,Δx,xi,allow_curvature_at_ends)
	coeffReduced=eqs\y # a1 b1 d1 d2 d3 dn-1
	tmp = zeros(nSplines)
	splineCoeffs = SplineCoeffs(copy(tmp),copy(tmp),copy(tmp),copy(tmp),xi)
	splineCoeffs.a[1] = coeffReduced[1]
	splineCoeffs.b[1] = coeffReduced[2]
	if allow_curvature_at_ends
		splineCoeffs.c[1] = coeffReduced[3]	
		splineCoeffs.d[1] = coeffReduced[4]
		nNormalSubstitutions = nSplines
		addon_for_d = 3
	else
		splineCoeffs.c[1] = 0	
		splineCoeffs.d[1] = coeffReduced[3]
		nNormalSubstitutions = nSplines-1
		addon_for_d = 2
	end
	for i in 2:nNormalSubstitutions
		splineCoeffs.a[i] = splineCoeffs.a[i-1]+splineCoeffs.b[i-1]*Δx+splineCoeffs.c[i-1]*Δx^2+splineCoeffs.d[i-1]*Δx^3
		splineCoeffs.b[i] = splineCoeffs.b[i-1]+2*splineCoeffs.c[i-1]*Δx+3*splineCoeffs.d[i-1]*Δx^2
		splineCoeffs.c[i] = splineCoeffs.c[i-1]+3*splineCoeffs.d[i-1]*Δx
		splineCoeffs.d[i] = coeffReduced[i+addon_for_d]
	end
	if ~allow_curvature_at_ends
		i = nSplines
		splineCoeffs.a[i] = splineCoeffs.a[i-1]+splineCoeffs.b[i-1]*Δx+splineCoeffs.c[i-1]*Δx^2+splineCoeffs.d[i-1]*Δx^3
		splineCoeffs.b[i] = splineCoeffs.b[i-1]+2*splineCoeffs.c[i-1]*Δx+3*splineCoeffs.d[i-1]*Δx^2
		splineCoeffs.c[i] = splineCoeffs.c[i-1]+3*splineCoeffs.d[i-1]*Δx
		splineCoeffs.d[i] = -splineCoeffs.c[i]/(3*Δx)
	end
	return splineCoeffs
end
"""

	yp = evalSpline(splineCoeffs,xp)


Return spline function values of splines created with splineCoeffs = splineFit(x, y, nSplines).

`splineCoeffs`: coefficients of spline functions.

`xp`: 			x-coordinates.

`yp`: 			corresponding y-coordinates.

# Example
```julia-repl
julia> x = collect(1:100);
julia> y = 3*log.(x)+0.2*randn(length(x))+2*sin.(x/10);
julia> nSplines = 6;
julia> xp = collect(-0:3:100);
julia> splineCoeffs,yp = splineFit(x, y, nSplines, xp, curved_ends=true, doplot=true);
julia> xp2 = [5;15];
julia> yp2 = evalSpline(splineCoeffs,xp2)
2-element Vector{Float64}:
  5.476211419705948
 10.17353202236643
julia> scatter!(xp2,yp2, mc=:orange,ms=7)
```
"""
function evalSpline(splineCoeffs,xp)
	nSplines = length(splineCoeffs.a)
	x1 = splineCoeffs.x[1]
	xn = splineCoeffs.x[end]
	Δx = (xn-x1)/(nSplines-1)
	splineID::Array{Int} = fill(-1,size(xp))
	splineID=Int.(ceil.((xp.-x1)/Δx))
	splineID[splineID.>nSplines] .= nSplines
	splineID[splineID.<1] .= 1
	yp = evalSpline(splineCoeffs,splineID,xp)
	return yp
end

function evalSpline(splineCoeffs,sid,xp)
	xloc = xp .- splineCoeffs.x[sid]
	yp = splineCoeffs.a[sid] .+ splineCoeffs.b[sid].*xloc .+ splineCoeffs.c[sid].*xloc.^2 .+ splineCoeffs.d[sid].*xloc.^3
	return yp
end

function plotSplines(x,y,splineCoeffs,xi,Δx)
	scatter(x,y,legend=false, mc=:grey)
	endmarker = 0.05*(maximum(y)-minimum(y))
	for i in 1:nSplines
		xp = collect(range(xi[i],xi[i]+Δx,length=30))
		yp = evalSpline(splineCoeffs,i,xp)
		plot!(xp,yp, lw = 2, linecolor = RGB(0,0,0),showlegend=false)
		delim_color = RGB(0,0.1,0.9)
		plot!([xp[1];xp[1]],[yp[1]+endmarker;yp[1]-endmarker], lw = 1, lc=delim_color,showlegend=false)
		plot!([xp[end];xp[end]],[yp[end]+endmarker;yp[end]-endmarker], lw = 1, lc=delim_color,showlegend=false)
	end
	gui()
end

function createEQS(x,y,splineID,xi,SS,allow_curvature_at_ends)
	n = length(x)
	if allow_curvature_at_ends
		nVariables = nSplines+3
	else
		nVariables = nSplines+1
	end
	eqs = zeros(n,nVariables)
	for i in 1:n
		sid = splineID[i]
		xloc = x[i]-xi[sid]
		eqs[i,:] = SS.a[sid,:] + SS.b[sid,:]*xloc + SS.c[sid,:]*xloc^2 + SS.d[sid,:]*xloc^3
	end
	return eqs
end

function splineID2(x,x1,xn,Δx,nSplines)
	splineID::Array{Int} = fill(-1,size(x))
	splineID=Int.(ceil.((x.-x1)/Δx))
	splineID[splineID.>nSplines] .= nSplines
	splineID[splineID.<1] .= 1
	xi = collect(range(x1,stop=xn,length=nSplines+1))
	pop!(xi)
	return splineID, xi
end

mutable struct SSmat
	a
	b
	c
	d
end

function createSS(nSplines,Δx,x1,allow_curvature_at_ends)
	xn = x1 + nSplines*Δx
	if allow_curvature_at_ends
		nVariables = nSplines+3	# a1 b1 c1 d1 d2 d3 ... dn
	else
		nVariables = nSplines+1 # a1 b1 d1 d2 d3 ... dn-1
	end
	SSa = zeros(nSplines,nVariables)
	SSb = copy(SSa)
	SSc = copy(SSa)
	SSd = copy(SSa)

	SSa[1,1] = 1
	SSb[1,2] = 1
	if allow_curvature_at_ends
		SSc[1,3] = 1
		SSd[1,4] = 1
	else
		if nSplines > 1
			SSd[1,3] = 1
		end
	end
	

	for i in 2:nSplines
		for j in 1:nVariables
			SSa[i,j] += SSa[i-1,j]
			SSa[i,j] += SSb[i-1,j]*Δx
			SSa[i,j] += SSc[i-1,j]*Δx^2
			SSa[i,j] += SSd[i-1,j]*Δx^3

			SSb[i,j] += SSb[i-1,j]
			SSb[i,j] += SSc[i-1,j]*2*Δx
			SSb[i,j] += SSd[i-1,j]*3*Δx^2

			SSc[i,j] += SSc[i-1,j]
			SSc[i,j] += SSd[i-1,j]*3*Δx

			#SSd[i,j] += SSd[i-1,j]
		end
		if i < nSplines && ~allow_curvature_at_ends
			SSd[i,i+2] = 1
		elseif allow_curvature_at_ends
			SSd[i,i+3] = 1
		end
	end
	if nSplines > 1 && ~allow_curvature_at_ends
		i = nSplines
		for j in 1:nSplines+1
			SSd[i,j] -= SSc[i-1,j]/(3*Δx)
		end
	end
	
	SS = SSmat(SSa,SSb,SSc,SSd)
	return SS
end

















