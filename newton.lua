local function newtonRootFind(f, df_dx, x0, maxiter, ...)
	local x = x0
	maxiter = maxiter or 100
	for iter=1,maxiter do
		local f_val = f(x, ...)
		local df_dx_val = df_dx(x, ...)
		local dx = -f_val / df_dx_val
		if math.abs(dx) < 1e-7 then break end 
		x = x + dx
--[[
		print(
			'f = '..f_val
			..' df_dz = '..df_dz_val
			..' dz = '..dz
			..' z = '..z
		)
--]]
		if iter == maxiter then return math.nan end
	end
	return x
end

return newtonRootFind
