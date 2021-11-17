local function Sign(X, Y)
	if Y < 0 then return -math.abs(X) end
	return math.abs(X)
end

-- Reference: From Numath Library By Tuan Dang Trong in Fortran 77.
-- C++ Release 1.0 By J-P Moreau, Paris.
-- (www.jpmoreau.fr)
local function BesselFunctionFirstKindOrder1(X)
--[[
	***********************************************************************
	This subroutine calculates the First Kind Bessel Function of
	order 1, for any real number X. The polynomial approximation by
	series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	REFERENCES:
	M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	VOL.5, 1962.
	***********************************************************************
--]]
	local P1 = 1.0
	local P2 = 0.183105E-2
	local P3 = -0.3516396496E-4
	local P4 = 0.2457520174E-5
	local P5 = -0.240337019E-6
	local P6 = 0.636619772
	local Q1 =  0.04687499995
	local Q2 = -0.2002690873E-3
	local Q3 = 0.8449199096E-5
	local Q4 = -0.88228987E-6
	local Q5 =  0.105787412E-6
	local R1 =  72362614232.0
	local R2 = -7895059235.0
	local R3 = 242396853.1
	local R4 = -2972611.439
	local R5 = 15704.48260
	local R6 = -30.16036606
	local S1 = 144725228442.0
	local S2 = 2300535178.0
	local S3 = 18583304.74
	local S4 = 99447.43394
	local S5 = 376.9991397
	local S6 = 1.0

	local AX = math.abs(X)
	if AX < 8 then
		local Y = X*X
		local FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
		local FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
		return X*(FR/FS)
	else
		local Z = 8.0/AX
		local Y = Z*Z
		local XX = AX-2.35619491
		local FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
		local FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
		return math.sqrt(P6/AX)*(math.cos(XX)*FP-Z*math.sin(XX)*FQ)*Sign(S6,X)
	end
end

return BesselFunctionFirstKindOrder1
