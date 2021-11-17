local function integrateTrapezoid(f, xL, xR, n)
	n = n or 100
	local dx = (xR - xL) / n
	local sum = .5 * (f(xL) + f(xR))
	for i=1,n-1 do
		local x = xL + dx * i
		sum = sum + f(x)
	end
	return sum * dx
end

return integrateTrapezoid
