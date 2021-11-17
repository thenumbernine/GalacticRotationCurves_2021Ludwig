-- https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean

local math = require 'ext.math'

local epsilon = 1e-7
local function arithmeticGeometricMean(a,g)
	for iter=1,100 do
		local an = .5 * (a + g)
		local delta = math.abs(a - g)
		if delta < epsilon then return an end
		a, g = an, math.sqrt(a * g)
	end
	return math.nan
end

return arithmeticGeometricMean
