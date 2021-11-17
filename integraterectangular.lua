local function integrateRectangular(f, xL, xR, n)
	n = n or 100
	local dx = (xR - xL) / n
	local sum = 0
	for i=1,n do
		local x = xL + dx * (i - .5)
		sum = sum + f(x)
	end
	return sum * dx
end

return integrateRectangular
