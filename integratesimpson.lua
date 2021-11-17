local function integrateSimpson(f, xL, xR, n)
	if n then
		n = math.ceil(n / 2) * 2
	else
		n = 100
	end
	local dx = (xR - xL) / n
	local sum = f(xL) + f(xR)
	local sum2 = 0
	for i=1,n/2-1 do
		local x = xL + dx * 2 * i
		sum2 = sum2 + f(x)
	end
	sum2 = sum2 * 2
	local sum4 = 0
	for i=1,n/2 do
		local x = xL + dx * (2 * i - 1)
		sum4 = sum4 + f(x)
	end
	sum4 = sum4 * 4
	return (sum + sum2 + sum4) * dx / 3
end

return integrateSimpson
