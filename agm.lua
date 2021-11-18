-- https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean

local math = require 'ext.math'

--[[
a = a0
g = g0
epsilon = stop epsilon
maxiter = max iterations
callback = callback every iteration: callback(a, g, i), starting with (a0, g0)
	return true to stop early
--]]
local function arithmeticGeometricMean(a, g, epsilon, maxiter, callback)
	epsilon = epsilon or 1e-7
	maxiter = maxiter or 100
	for i=0,maxiter-1 do
		if callback then 
			if callback(a, g, i) then
				return a
			end
		end
		local an = .5 * (a + g)
		local delta = math.abs(a - g)
		if delta < epsilon then return an end
		a, g = an, math.sqrt(a * g)
	end
	return math.nan
end

return arithmeticGeometricMean
