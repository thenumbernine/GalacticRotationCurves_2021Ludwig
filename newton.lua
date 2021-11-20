local math = require 'ext.math'

local epsilon = 1e-7
local maxiter = 100

local function newtonRootFind(f, df_dx, x, ...)
	local lastabsdx = math.huge
	for iter=1,maxiter do
		local f_val = f(x, ...)
		local df_dx_val = df_dx(x, ...)
		local dx = -f_val / df_dx_val
		local absdx = math.abs(dx)
		if absdx < epsilon and absdx == lastabsdx then return x end 
		lastabsdx = absdx
		x = x + dx
--[[
		print(
			'f = '..f_val
			..' df_dx = '..df_dx_val
			..' dx = '..dx
			..' x = '..x
		)
--]]
	end
	return math.nan
end

return newtonRootFind
