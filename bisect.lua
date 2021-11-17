local math = require 'ext.math'

--[[
if this is just a bisection method then it should initialize with f(xL) and f(xR) on either side of the root
--]]
local function bisectRootFind(f, xL, xR, maxiter, ...)
	maxiter = maxiter or 100
	local fL = f(xL, ...)
	local fR = f(xR, ...)
	if (fL < 0) == (fR < 0) then
		return false,
			"can't do bisection without initializing with samples on opposite evaluations of the root\n"
			.." f(xL="..xL..") = " .. fL .. "\n"
			.." f(xR="..xR..") = " .. fR .. "\n"
			.." for args = " .. require 'ext.tolua'{...}
	end
	for i=1,maxiter do
		local dx = xR - xL
		local xmid = xL + .5 * dx 
		if dx < 1e-15 then return xmid end
		local fmid = f(xmid, ...)
		
		if (fmid < 0) == (fL < 0) then
			xL, fL = xmid, fmid
		else
			xR, fR = xmid, fmid
		end
	end
	return math.nan
end

return bisectRootFind
