

using Revise

include("/Users/Andreas/Desktop/Julia/projects/curveFit/splineFit.jl")



test_ID = 1

x = collect(1:200)
y = 3*log.(x)+0.2*randn(length(x))+2*sin.(x/10)
nSplines = 9
#y[100:200] .= 10
xp = collect(-0:12:200)
xp2 = [5,15]		
if test_ID == 1
	splineCoeffs,yp = splineFit(x, y, nSplines, xp, curved_ends=true, doplot=true)
	yp2 = evalSpline(splineCoeffs,xp2)
	scatter!(xp2,yp2, mc=:orange,ms=7)
	gui()
elseif test_ID == 2
	splineCoeffs = splineFit(x, y, nSplines, curved_ends=false, doplot=true)
	yp2 = evalSpline(splineCoeffs,xp2)
	scatter!(xp2,yp2, mc=:orange,ms=7)
	gui()
end

nothing