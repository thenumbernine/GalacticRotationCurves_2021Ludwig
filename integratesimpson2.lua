local function integrateSimpson3_8(f, xL, xR, n)
	if n then
		n = math.ceil(n / 3) * 3
	else
		n = 99
	end
	local dx = (xR - xL) / n
	local sum = f(xL) + f(xR)
	local sum2 = 0
	for i=1,n/3-1 do
		local x = xL + dx * 3 * i
		sum2 = sum2 + f(x)
	end
	sum2 = sum2 * 2
	local sum3 = 0
	for i=1,n-1 do
		if i % 3 ~= 0 then
			local x = xL + dx * i
			sum3 = sum3 + f(x)
		end
	end
	sum3 = sum3 * 3
	return (sum + sum2 + sum3) * dx * (3 / 8)
end

return integrateSimpson3_8
