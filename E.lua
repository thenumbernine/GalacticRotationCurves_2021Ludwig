-- complete elliptic integral of the second kind

local math = require 'ext.math'
local agm = require 'agm'

local function completeEllipticIntegralSecondKind(k)
	local ksq = k * k
	if ksq > 1 then return math.nan end
	if ksq == 1 then return 1 end
	local sum = 0
	local _2n_minus_1 = .5
	local a0 = 1
	local g0 = math.sqrt(1 - ksq)
	local cn = math.abs(k)
	local ainf = agm(a0, g0, nil, nil, function(an,gn,i)
		local cnsq = cn * cn
		sum = sum + _2n_minus_1 * cnsq 
		_2n_minus_1 = _2n_minus_1 * 2
		local an_plus_1 = .5 * (an + gn)	-- maybe pass this?  though idk , you'd have to restructure too much
		-- TODO rearrange this by storing a previous cn_minus_1
		-- and then just redo the sum for cn = (cn_minus_1^2) / (4 * an)
		-- but that would take skipping this for i==0
		-- what's slower? an add-mul or a branch?  i'm thinking the branch.
		cn = cnsq / (4 * an_plus_1)
	end)
	return math.pi / (2 * ainf) * (1 - sum)
end

return completeEllipticIntegralSecondKind
