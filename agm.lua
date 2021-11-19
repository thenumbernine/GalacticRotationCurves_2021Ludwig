-- https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean

local math = require 'ext.math'

local agm = {}
setmetatable(agm, agm)

-- overrideable defaults
agm.epsilon = 1e-7
agm.maxiter = 100

--[=[ debug -- watching the values
local ffi = require 'ffi'
ffi.cdef[[
typedef union {
	double d;
	long l;
} dl;
]]

local x = ffi.new'dl'
require 'ffi.c.stdio'
--]=]

--[[
a = a0
g = g0
epsilon = stop epsilon
maxiter = max iterations
callback = callback every iteration: callback(a, g, i), starting with (a0, g0)
	return true to stop early
--]]
local function arithmeticGeometricMean(a, g, callback, epsilon, maxiter)
	epsilon = epsilon or agm.epsilon
	maxiter = maxiter or agm.maxiter
	local lastdelta = math.huge
	for i=0,maxiter-1 do
		if callback then 
			if callback(a, g, i) then
				return a
			end
		end
		local an = .5 * (a + g)
		local delta = math.abs(a - g)

--[[
x.d = delta
ffi.C.printf('δ=%016lx,%.30e\t', x.l, x.d)
x.d = a
ffi.C.printf('a=%016lx,%.30e\t', x.l, x.d)
x.d = g
ffi.C.printf('g=%016lx,%.30e\n', x.l, x.d)
--]]

		-- delta == 0? we're done
		if delta == 0 
		-- delta has converged to some small value, it isn't changing, we're done
		or (delta <= epsilon and delta == lastdelta) 
		then 
--print'done'
			return an 
		end
		lastdelta = delta
		a, g = an, math.sqrt(a * g)
	end
--print'maxiter'
	return math.nan
end
agm.agm = arithmeticGeometricMean

-- how much will this slow things down?
-- ah well, just use agm.agm to avoid the tail call
function agm:__call(...)
	return arithmeticGeometricMean(...)
end

return agm
