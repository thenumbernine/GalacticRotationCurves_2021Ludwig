-- https://en.wikipedia.org/wiki/Elliptic_integral#Complete_elliptic_integral_of_the_first_kind
-- "See Carlson (2010, 19.8) for details"

local math = require 'ext.math'
local agm = require 'agm'

local function completeEllipticIntegralFirstKind(k)
	local ksq = k * k
	if ksq > 1 then return math.nan end
	if ksq == 1 then return math.huge end
	return math.pi / (2 * agm(1, math.sqrt(1 - ksq)))
end

return completeEllipticIntegralFirstKind
