#!/usr/bin/env lua128

local range = require 'ext.range'
local table = require 'ext.table'
local timer = require 'ext.timer'
local math = require 'ext.math'
local gnuplot = require 'gnuplot'
local matrix = require 'matrix'
local symmath = require 'symmath'
local J1 = require 'J1'
		
symmath.fixVariableNames = true
symmath.tostring = symmath.export.SingleLine

local Gstorage = {}
setmetatable(_G, {
	__newindex = function(...)
		local t, k, v = ...
		local oldv = Gstorage[k]
--		print(debug.traceback())
		print(tostring(k)..' = '..tostring(v)..(oldv ~= nil and oldv ~= v and (' (from '..oldv..')') or ''))
		Gstorage[k] = v
	end,
	__index = Gstorage,
})


-- Section 5, after eqn 5.3 (sloppiest research paper ever)
local Msun = 1.9891e+30	-- kg 
-- wikipedia:
--Msun = 1.98847e+30	-- kg
print('Msun = '..Msun)

-- Appendix D, after eqn D11
local Lsun = 3.828e+26 	-- W
-- wikipedia:
--Lsun = 3.828e+26	-- W
print('Lsun = '..Lsun)

-- mass-to-light ratio
local UpsilonSun = Msun / Lsun	-- kg / W ~ wikipedia: 5133 kg / W

-- me
local m_in_pc = 648000 / math.pi * 149597870700
print('m_in_pc = '..m_in_pc)

-- gravitational "constant":
local G_in_m3_per_kg_s2 = 6.67384e-11
print('G_in_m3_per_kg_s2 = '..G_in_m3_per_kg_s2)

-- speed of light:
local c_in_m_per_s = 299792458
print('c_in_m_per_s = '..c_in_m_per_s)

-- natural unit proportionality between kilograms and meters:
local G_over_c2_in_m_per_kg = G_in_m3_per_kg_s2 / (c_in_m_per_s * c_in_m_per_s)
print('G_over_c2_in_m_per_kg = '..G_over_c2_in_m_per_kg)

local n = 1000
local xvec = matrix{n}:lambda(function(i) return (i-1)/(n-1) end)	-- [0,1] edge coords

-- assumes xvec has values from 0 to 1
local function makePow10Range(rmin, rmax)
	local log10rmin = math.log(rmin) / math.log(10)
	local log10rmax = math.log(rmax) / math.log(10)
	return xvec:map(function(x) return 10^((log10rmax - log10rmin) * x + log10rmin) end)
end


--[[
where is this found ...
Appendix D:
in text after eqn D11
in text after eqn D12, repeated twice for (r0,alpha0) and (r1,alpha1)
in text after eqn D21, last sentence of the paper.
too bad it never has a formal equation label ... that'd make this paper a bit easier to follow.

'd' is in ... what units? kpc I think
'r' is in ... what units? kpc I think
'alpha' is in ... what units? unitless I think
--]]
local function r_for_d_alpha(d, alpha)
	local r = d * (math.pi / (180 * 3600)) * alpha
	return r
end
local function d_for_r_alpha(r, alpha)
	local d = (180 * 3600 / math.pi) * r / alpha
	return d
end
local function alpha_for_r_d(r, d)
	local alpha = (180 * 3600 / math.pi) * r / d
	return alpha
end

-- eqn D14 def of sb ... requires b0 which doesn't exist
-- sb(alpha) = b0 + alpha * (b1 + alpha * (b2 + alpha * b3))	
-- but just after D14 in text it says b0 == s0, so we get the eqn D17 def of sb(alpha):
local function sb_for_alpha_eqn_D_17(alpha)
	--[[
	return s0 + alpha * (b1 + alpha * (b2 + alpha * (b3 + alpha * b4)))
	--]]
	--[[
	return s0 + b1 * alpha + b2 * alpha^2 + b3 * alpha^3 + b4 * alpha^4
	--]]
	-- [[
	local sum = 0
	for i=bPolyOrder,1,-1 do
		sum = sum + _G['b'..i]
		sum = sum * alpha
	end
	return sum + s0
	--]]
end
local function dSb_dalpha(alpha) 
	--[[
	return b1 + 2 * b2 * alpha + 3 * b3 * alpha^2 + 4 * b4 * alpha^3
	--]]
	-- [[
	local sum = 0
	for i=bPolyOrder,2,-1 do
		sum = sum + i * _G['b'..i]
		sum = sum * alpha
	end
	return sum + b1
	--]]
end

-- eqn D14 def of sd ... requires d0, d1 which doesn't exist
-- sd(alpha) = d0 + alpha * (d1 + alpha * (d2 + alpha * d3))	
-- how is this derived?
local function sd_for_alpha_eqn_D_17(alpha) 
	-- eqn D17 def of sd
	local dalpha = alphae - alpha
	--[[
	return se + dalpha * dalpha * (d2 + dalpha * (d3 + dalpha * (d4)))
	--]]
	--[[
	return se + d2 * dalpha^2 + d3 * dalpha^3 + d4 * dalpha^4
	--]]
	-- [[
	local sum = 0
	for i=dPolyOrder,2,-1 do
		sum = sum + _G['d'..i]
		sum = sum * dalpha
	end
	return sum * dalpha + se
	--]]
end
local function dSd_dalpha(alpha) 
	local dalpha = alphae - alpha
	--[[
	return se + 2 * d2 * dalpha + 3 * d3 * dalpha^2 + 4 * d4 * dalpha^3
	--]]
	-- [[
	local sum = 0
	for i=dPolyOrder,2,-1 do
		sum = sum + i * _G['d'..i]
		sum = sum * dalpha
	end
	return sum + se
	--]]
end

-- eqn D13 
local function s_for_alpha_eqn_D_13(alpha) 
	if alpha <= alpha0 then
		return sb_for_alpha_eqn_D_17(alpha) 
	end
	if alpha <= alphae then
		return sd_for_alpha_eqn_D_17(alpha)
	end
	return se
end


-- eqn D.18.a
-- evaluate b1 according to se, s0, alphae, alpha0, b2 ... bn, d3 ... dn
local function b1_eqn_D_18_a()
	return 1 / (alphae + alpha0) * (
		2 * (se - s0)
		- range(2,bPolyOrder):mapi(function(i)
			return (i * (alphae - alpha0) + 2 * alpha0) * alpha0 ^ (i - 1) * _G['b'..i]
		end):sum()
		- range(3,dPolyOrder):mapi(function(j)
			return (j - 2) * (alphae - alpha0) ^ j * _G['d'..j]
		end):sum()
	)
end

-- eqn D.18.b
-- evaluate d2 according to the same
local function d2_eqn_D_18_b()
	return -1 / (alphae * alphae - alpha0 * alpha0) * ( 
		se - s0
		+ range(2,bPolyOrder):mapi(function(i)
			return (i - 1) * alpha0 ^ i * _G['b'..i]
		end):sum()
		+ range(3,dPolyOrder):mapi(function(j)
			return (alphae + (j - 1) * alpha0) * (alphae - alpha0) ^ (j - 1) * _G['d'..j]
		end):sum()
	)
end

-- eqn D.20
local function se_eqn_D_20(mu_for_alpha_fn)
	return math.log(alphae / alphaeff) / math.log(2 * math.log(10) / 5 * (mu_for_alpha_fn(alphae) - mu0))	
end

-- equation C.2
-- where A, B, R, Z are normalized somehow
-- I'm not using this one at the moment
local function phi_for_R_Z_eqn_C_2(R,Z)
	return -G_in_m3_per_kg_s2 * M / math.sqrt(R*R + (A + math.sqrt(B * B + Z * Z))^2)
end

-- eqn 8.2
local function r_for_i(i)
	-- 'i' is the index, right? the subscript of 'i' right? and not the imaginary unit *and* the denotation of r ...
	return rspiral * math.exp(2 * math.pi * i * kspiral)
end

-- eqn 8.3
local function y_for_i(i)
	-- to the v?  or is that a reference? what is v?
	-- why does the eqn say exp(...)^v and not exp(v * ...) ?
	return y0 * math.exp(2 * math.pi * i * kspiral * v)
end

-- eqn 8.1
local function Y_for_r_eqn_8_1(r)
	local sum = 1
	for i=0,YOrder do
-- I guess this is the constant they provided?
-- then again, I am only ever seeing gammai = gamma0 in the text description of section 8 on NGC 3198, sooo ...
-- maybe 'i' is a subscript and it just so happens that all gamma[i] = gamma0 ?
--		local gammai = gamma_for_i(i)
		sum = sum + (y_for_i(i) * gammai / math.pi) / ((r - r_for_i(i))^2 + gammai^2)
	end
	return sum
end

-- eqn 8.5
-- basically same as S(alpha) from eqn D13, but with (r, r0, re) swapped with (alpha, alpha0, alphae)
local function s_for_r_eqn_8_5(r)
	--[[ 
	-- since r and alpha are linear, can we just do this? 
	-- no.  
	-- the coefficient of difference will scale with the polynomials ...
	if r <= r0 then
		return sb_for_alpha_eqn_D_17(r)
	end
	if r <= re then
		return sd_for_alpha_eqn_D_17(r)
	end
	return se
	--]]
	-- [[ convert to alpha and evaluate accordingly
	local alpha = alpha_for_r_d(r, d)
	return s_for_alpha_eqn_D_13(alpha)
	--]]
end

-- eqn 8.4
-- ... has both these equations ... and they aren't even equal ... smh
local function normrho_z_eq_0_eqn_8_4(r)
	--[[ can't do this without s1 and s1 of the galaxy being defined ... which it is not for NGC 3198
	local alpha = alpha_for_r_d(r, d)
	return Y_for_r_eqn_8_1(r) * 10 ^ (-2/5 * (mu_for_alpha_eqn_D_11(alpha) - mu0))
	--]]
	-- [[ looks good
	return 
		Y_for_r_eqn_8_1(r) 
		* math.exp(-(r / reff) ^ (1 / s_for_r_eqn_8_5(r)))
	--]]
	-- "where s(r) is the Sérsic profile adjusted to the luminosity:"
end

-- eqn D8, doesn't require alpha1 (derived from r1) or alpha0 (derived from r0) or s2, like eqn D11 requires
-- but what is "s" ?
-- after eqn D1: "where s is the Sersic index"
-- but this 's' varies with alpha ... correct?
local function mu_for_alpha_eqn_D_8(alpha)
	local s = s_for_alpha_eqn_D_13(alpha)
	return mu0 + 5 / (2 * math.log(10)) * (alpha / alphaeff) ^ (1 / s)
end

local function normrho_for_r_z_eq_0_eqn_9_1_a(r)
	local alpha = alpha_for_r_d(r, d)
	return Y_for_r_eqn_8_1(r) * 10 ^ (-2/5 * (mu_for_alpha_eqn_D_8(alpha) - mu0))
end

local function normrho_for_r_z_eq_0_eqn_9_1_b(r)
	local alpha = alpha_for_r_d(r, d)
	if alpha < alpha0 then
		return Y_for_r_eqn_8_1(r) * math.exp(-(r / r0) ^ (1 / sb_for_alpha_eqn_D_17(alpha)))
	else
		return Y_for_r_eqn_8_1(r) * math.exp(-(r / r0) ^ (1 / sd_for_alpha_eqn_D_17(alpha)))
	end
end



-- expressions ... putting all in one place so the symbolic vars will be too
local normrho_for_r_z_eqn_5_2_b	-- also a candidate for eqn C 9 lhs over lambda
local normphi_for_r_z_eqn_5_2_a 
local eqn_C_9_root_expr_based_on_5_2_b, dz_eqn_C_9_root_expr_based_on_5_2_b	-- internal to the newton root-finder
local galacticWidth_for_r_eqn_C_9_based_on_eqn_5_2_b
local eqn_C_9_rhs_over_lambda		-- purely for debugging
local dr_beta_for_r_eqn_5_3
local g_for_r_eqn_5_3
timer("deriving root-finding", function()

	-- let the newton func and deriv exist in this scope
	do	-- use symmath vars in this scope
		local Constant = symmath.Constant
		local var = symmath.var
		local frac = symmath.frac
		local sqrt = symmath.sqrt
		local exp = symmath.exp
		-- these are in the same order the arguments are provided in the compiled function
		local argvars = table{symmath.vars('r', 'z')}
		local r, z = argvars:unpack()
		-- globals:
		local a, b, l, lambda, mu0 = symmath.vars('a', 'b', 'l', 'lambda', 'mu0')
	

		local beta = var('beta', {r})
		local f = var('f', {r})
		local g = var('g', {r})
		local rs = var'rs'
		local normrho = var('ϱ', {r})
		local normphi = var('φ', {r})

		--eqn 5.1 is derived from the Abel equation 4.13:
		local eqn_5_1 = ((beta + r * beta:diff(r)) * g):eq(beta * ((beta - r * beta:diff(r)) * beta + f * (1 - beta^2)))
		print('eqn 5.1:', eqn_5_1)
	
		-- tmp for eqn 5.2
		local tmp = sqrt(b*b + z*z)

		-- eqn 5.2.a
		-- normalized potential
		local normphi_eqn_5_2_a_expr = normphi:eq(-rs / (2 * sqrt(r^2 + (a + tmp)^2)))
		print('eqn 5.2.a:', normphi_eqn_5_2_a_expr)

		-- this is assuming rs is a point mass, right?  and not distributed?
		local dr_normphi_eqn_5_2_a_expr = normphi_eqn_5_2_a_expr:diff(r):prune() 
		print('d/dr of eqn 5.2.a:', dr_normphi_eqn_5_2_a_expr)

		-- eqn 5.2.b
		-- normalized density
		-- also eqn C.4 but with lambda's definition of C.5 substituted
		local normrho_eqn_5_2_b_expr = normrho:eq(
			(b*b / (3 * lambda)) * (
				a*r*r 
				+ (a + 3 * tmp) * (a + tmp)^2
			) / (
				(r*r + (a + tmp)^2)^frac(5,2)
				* tmp^3
			)
		)
		print('eqn 5.2.b:', normrho_eqn_5_2_b_expr)

		-- eqn 5.3 requires substitutions (text after eqn 5.2):
		local f_eqn_5_3_expr = f:eq( frac(3,2) * lambda * rs * normrho * r^2 )
		print('f for eqn 5.3:', f_eqn_5_3_expr)

		f_eqn_5_3_expr = f_eqn_5_3_expr:subst(normrho_eqn_5_2_b_expr, z:eq(0))()
		print('f for eqn 5.3:', f_eqn_5_3_expr)

		local g_eqn_5_3_expr = g:eq(r * normphi:diff(r))
		print('g for eqn 5.3:', g_eqn_5_3_expr)
		
		g_eqn_5_3_expr = g_eqn_5_3_expr:subst(normphi_eqn_5_2_a_expr, z:eq(0))()
		print('g for eqn 5.3:', g_eqn_5_3_expr)

		-- deriving eqn 5.3 ...
		local eqn_5_3 = eqn_5_1:clone()
		print('eqn 5.1 with substitutions:', eqn_5_3) 
		
		eqn_5_3 = eqn_5_3()
		print('eqn 5.1 with substitutions:', eqn_5_3) 
		
		eqn_5_3 = (eqn_5_3 - beta * g + r * beta^2 * beta:diff(r))()
		print('eqn 5.1 with substitutions:', eqn_5_3) 
		
		eqn_5_3 = (eqn_5_3 / (r * (g + beta^2)))()
		print('eqn 5.1 with substitutions:', eqn_5_3) 
		
		eqn_5_3 = eqn_5_3:subst(f_eqn_5_3_expr, g_eqn_5_3_expr)
		print('eqn 5.1 with substitutions:', eqn_5_3) 

		-- now it looks like they use explicit integration to solve this, starting at some "lbeta" point ......?
		-- it says the point of evaluation being separate from r=0
		-- it never says implicit/explicit/whatever type of integration.

		-- ok, eqn_5_3 is now a def of dbeta/dr, ready for integration 

--[[ TODO
		local mu_for_alpha_eqn_D_11_expr = 

		-- eqn D.12
		local normrho_eqn_D_12_a_expr = 10 ^ (frac(-2, 5) * (mu_for_alpha_eqn_D_11_expr - mu0))
--]]

		-- eqn C.9 
		-- uses a normalized density function based on eqn C.4
		local eqn_C_9_rhs_over_lambda_expr = exp(-l^2 / 2)

		-- TODO eqn C.9 modified with a normalized density function based on eqn D.12


		local eqn_C_9_root_expr_based_onn_5_2_b_expr = normrho_eqn_5_2_b_expr:rhs() - eqn_C_9_rhs_over_lambda_expr
		print(eqn_C_9_root_expr_based_onn_5_2_b_expr)
	
		-- TODO - Variable evaluateDerivative is commented out in Variable because it is also duplciated in Derivative's prune() rule ...
		
		local dz_eqn_C_9_root_expr = eqn_C_9_root_expr_based_onn_5_2_b_expr:diff(z):prune()
		print(dz_eqn_C_9_root_expr)
		
		-- assert d/dx rhs == 0, 
		-- such that if width = width_lhs - width_rhs then d/dz (width) == d/dz (width_lhs)
		assert(Constant.isValue(eqn_C_9_rhs_over_lambda_expr:diff(z):prune(), 0))
		
		-- eqn C.4
		-- using the substitution of eqn C5 to replace R0 with lambda
		-- ok in eqn C 9 the z is replaced with "delta(z)" ... but delta is really just z. .. and they are solving for z and callign that variable "delta"
		normrho_for_r_z_eqn_5_2_b = symmath.export.Lua:toFunc{output={normrho_eqn_5_2_b_expr:rhs()}, input=argvars}

		normphi_for_r_z_eqn_5_2_a = symmath.export.Lua:toFunc{output={normphi_eqn_5_2_a_expr:rhs()}, input=argvars}

		eqn_C_9_rhs_over_lambda = symmath.export.Lua:toFunc{output={eqn_C_9_rhs_over_lambda_expr}, input=argvars}
		eqn_C_9_root_expr_based_on_5_2_b = symmath.export.Lua:toFunc{output={eqn_C_9_root_expr_based_onn_5_2_b_expr}, input=argvars}
		dz_eqn_C_9_root_expr_based_on_5_2_b = symmath.export.Lua:toFunc{output={dz_eqn_C_9_root_expr}, input=argvars}
	
		dr_beta_for_r_eqn_5_3 = symmath.export.Lua:toFunc{output={eqn_5_3:rhs()}, input={r, beta}}
	
		g_for_r_eqn_5_3 = symmath.export.Lua:toFunc{output={g_eqn_5_3_expr:rhs()}, input={r}}
	end	

	
	-- equation C9
	-- this is solving the nonlinear equation normrho(r,z) = exp(-1/2 l^2) for z, for fixed z
	galacticWidth_for_r_eqn_C_9_based_on_eqn_5_2_b = function(r)
--print('for r = '..r)
		local z = 1
		-- newton method
		local maxiter = 100
		for iter=1,maxiter do
			local f_val = eqn_C_9_root_expr_based_on_5_2_b(r, z)
			local df_dz_val = dz_eqn_C_9_root_expr_based_on_5_2_b(r, z)
			local dz = -f_val / df_dz_val
			if math.abs(dz) < 1e-7 then break end 
			z = z + dz
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
		return z
	end
end)

local function normrho_for_r_z_eq_0_eqn_5_2_b(r)
	return normrho_for_r_z_eqn_5_2_b(r,0)
end

--[[
also equation C.7
φ = ϕ/c^2 ... WHY DO YOU HAVE TO USE TWO DIFFERENT PHIS!?!?!?!?
r, rs, a, b, z are in kpc ... or equivalently normalized by R0 = 1 kpc ...
--]]
local function normphi_for_r_z_eq_0_eqn_5_2_a(r)
	return normphi_for_r_z_eqn_5_2_a(r,0)
end



print[[


NGC 1560

]]


-- section 6, on spheroid fitting to NGC 1560 and then on solving the Abel equation for beta (= v/c?)
a = .373					-- kpc
b = .300					-- kpc
rs = 0.00000701				-- kpc
lbeta = 8.29				-- kpc ... this is the starting point of integration according to the text after eqn 5.3
v_at_lbeta = 80				-- km/s
beta_at_lbeta = 0.000267	-- unitless <-> v / c
M = 7.3e+10 * Msun			-- kg

-- section 7 but they are normalized anyways so it shouldn't matter right? 
lambda = 0.134
-- used for section 6 graphs but defined in Appendix D:
d = 3.0 * 1000 -- kpc

--[[
section 5 after eqn 5.3:
reference point: beta(lbeta) = beta_l ~ 0.000267 <-> v ~ 80 km/s
sure enough, on figure 4, looks like the sample point and the graph overlap at r = lbeta.
however the graph extents go further in each direction.
i'm guessing they explicit integrated it forwards and backwards.

so this makes the smooth graphs of a spheroid rotation curve.  how to get the bumpy ones? maybe replace the rho and phi with the bumpy ones, and re-derive the Abel function?
--]]
local function makeNormalizedRotationCurve()
	local r = lbeta
	local beta = beta_at_lbeta
	local dr = lbeta / 1000
	local rAndBetaVec = table()

	rAndBetaVec:insert{r, beta}

	while r >= 0 do
		local dbeta_dr = dr_beta_for_r_eqn_5_3(r, beta)
		beta = beta + dbeta_dr * -dr
		r = r + -dr
		rAndBetaVec:insert{r, beta}
	end
	
	local r = lbeta
	local beta = beta_at_lbeta

--[[ only figure 4 rhs has graph data from r=0 through r=lbeta=8.29 to r=rmax=12
-- but this graph appears to be closest to fig 1 lhs, which goes from r=0 to r=lbeta=8.29
	while r < 10 do
		local dbeta_dr = dr_beta_for_r_eqn_5_3(r, beta)
		beta = beta + dbeta_dr * dr
		r = r + dr
		rAndBetaVec:insert{r, beta}
	end
--]]	
	rAndBetaVec:sort(function(a,b) return a[1] < b[1] end) 

	local rvec, betavec = matrix(rAndBetaVec):T():unpack()
	return rvec, betavec
end

-- col 1 = dist from galaxy center, in arcseconds
-- col 2 = mu_B = luminosity
local _1992_Broeils_Table_3 = table{
	{0,		22.27},
	{2,		22.30},
	{4,		22.31},
	{6,		22.29},
	{8,		22.31},
	{10,	22.27},
	{12,	22.27},
	{14,	22.27},
	{16,	22.27},
	{18,	22.25},
	{20,	22.26},
	{22,	22.30},
	{24,	22.31},
	{26,	22.34},
	{28,	22.35},
	{30,	22.42},
	{32,	22.36},
	{34,	22.36},
	{36,	22.38},
	{38,	22.36},
	{40,	22.36},
	{42,	22.41},
	{44,	22.47},
	{46,	22.47},
	{48,	22.50},
	{50,	22.51},
	{52,	22.50},
	{54,	22.54},
	{56,	22.54},
	{58,	22.57},
	{60,	22.60},
	{62,	22.63},
	{64,	22.65},
	{66,	22.68},
	{68,	22.73},
	{70,	22.75},
	{72,	22.81},
	{74,	22.81},
	{76,	22.83},
	{78,	22.89},
	{80,	22.95},
	{82,	22.92},
	{84,	22.98},
	{86,	23.01},
	{88,	23.02},
	{90,	23.03},
	{92,	23.09},
	{94,	23.11},
	{96,	23.14},
	{98,	23.11},
	{100,	23.13},
	{102,	23.20},
	{104,	23.22},
	{106,	23.23},
	{108,	23.24},
	{110,	23.32},
	{112,	23.32},
	{114,	23.33},
	{116,	23.28},
	{118,	23.34},
	{120,	23.39},
	{122,	23.43},
	{124,	23.39},
	{126,	23.46},
	{128,	23.49},
	{130,	23.52},
	{132,	23.53},
	{134,	23.57},
	{136,	23.62},
	{138,	23.60},
	{140,	23.62},
	{142,	23.63},
	{144,	23.61},
	{146,	23.67},
	{148,	23.70},
	{150,	23.73},
	{152,	23.75},
	{154,	23.79},
	{156,	23.84},
	{158,	23.80},
	{160,	23.84},
	{162,	23.90},
	{164,	23.90},
	{166,	23.95},
	{168,	24.01},
	{170,	24.02},
	{172,	24.02},
	{174,	24.04},
	{176,	24.12},
	{178,	24.13},
	{180,	24.11},
	{182,	24.18},
	{184,	24.18},
	{186,	24.25},
	{188,	24.24},
	{190,	24.23},
	{192,	24.29},
	{194,	24.30},
	{196,	24.30},
	{198,	24.34},
	{200,   24.32},
	{202,	24.36},
	{204,   24.41},
	{206,	24.44},
	{208,	24.43},
	{210,	24.41},
	{212,	24.54},
	{214,	24.55},
	{216,	24.47},
	{218,	24.51},
	{220,   24.57},
	{222,	24.56},
	{224,	24.57},
	{226,	24.59},
	{228,	24.67},
	{230,	24.73},
	{232,	24.75},
	{234,	24.71},
	{236,	24.82},
	{238,	24.86},
	{240,	24.86},
	{242,	24.91},
	{244,	24.89},
	{246,	24.99},
	{248,	24.93},
	{250,	24.97},
	{253,	24.99},
	{257,	25.05},
	{261,	25.11},
	{265,	25.08},
	{269,	25.16},
	{273,	25.18},
	{277,	25.37},
	{281,	25.33},
	{285,	25.43},
	{289,	25.37},
	{293,	25.39},
	{297,	25.44},
	{301,	25.41},
	{305,	25.61},
	{309,	25.63},
	{313,	25.72},
	{317,	25.73},
	{321,	25.75},
	{325,	25.60},
	{329,	25.80},
	{333,	25.88},
	{337,	25.81},
	{341,	25.94},
	{345,	25.99},
	{349,	25.91},
	{353,	26.13},
}
local alpha_in_arcsec_1992_Broeils_table_3, luminosity_1992_Broeils = matrix(_1992_Broeils_Table_3):T():unpack()
local r_in_kpc_1992_Broeils_table_3 = alpha_in_arcsec_1992_Broeils_table_3:map(function(alpha)
	return r_for_d_alpha(d, alpha)
end)


-- 1992 Broeis, table 4
-- col 1 is distance from galaxy center, in arcseconds
-- col 2 is velocity, in km/s
-- col 3 is velocity corrected for asymmetric drifts
local _1992_Broeils_Table_4 = table{
	{15, 4.5, 5.0},
	{30, 8.4, 8.9},
	{45, 14.0, 14.5},
	{60, 26.1, 26.4},
	{75, 28.7, 28.9},
	{90, 27.9, 27.8},
	{105, 32.1, 31.8},
	{120, 43.0, 42.8},
	{135, 47.8, 48.2},
	{150, 47.6, 48.4},
	{165, 49.8, 50.6},
	{180, 52.6, 53.5},
	{195, 56.4, 57.2},
	{210, 58.3, 59.1},
	{225, 58.8, 59.8},
	{240, 59.4, 60.3},
	{255, 59.1, 60.7},
	{270, 59.9, 62.1},
	{285, 61.9, 63.6},
	{300, 60.9, 62.0},
	{315, 60.6, 60.5},
	{330, 62.1, 60.3},
	{345, 64.5, 63.8},
	{360, 65.6, 66.1},
	{375, 67.1, 67.7},
	{390, 68.7, 70.4},
	{405, 70.4, 73.0},
	{420, 71.3, 74.2},
	{435, 72.0, 75.1},
	{450, 72.3, 75.2},
	{465, 73.2, 76.3},
	{480, 74.1, 77.2},
	{510, 74.4, 76.9},
	{540, 75.2, 77.5},
	{570, 76.6, 78.7},
}
for _,row in ipairs(_1992_Broeils_Table_4) do
	-- convert distance from arcseconds to kpc
	row[1] = r_for_d_alpha(d, row[1])
	-- convert from km/s to beta = v/c
	row[2] = row[2] * 1000 / c_in_m_per_s
	row[3] = row[3] * 1000 / c_in_m_per_s
end
local r_in_kpc_1992_Broeils_table_4, beta_1992_Broeils, beta_corr_1992_Broeils = matrix(_1992_Broeils_Table_4):T():unpack()

-- now try plotting it ...
-- should be figure 1 lhs, figure 2 lhs, or figure 4 rhs ... yeah they all have the same label
-- CHECK
-- this is figure 1 lhs
-- looks like the paper's samples for Broeils 1992 use the corrected velocities
local rvec, betavec = makeNormalizedRotationCurve()
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_1a_NGC_1560_normalized_rotation_curve_eqn_5.3.svg",
	style = 'data lines',
	xlabel = "r (kpc)",
	ylabel = 'v/c',
	title = "Normalized rotation curve of NGC 1560",
	xrange = {rvec[1], rvec[#rvec]},
	yrange = {(table.inf(betavec)), (table.sup(betavec))},
	data = {
		rvec,
		betavec,
		r_in_kpc_1992_Broeils_table_4, 
		beta_1992_Broeils, 
		beta_corr_1992_Broeils
	},
	key = 'bottom right',
	{using='1:2', title='beta fitted'},
	{using='3:4', title='beta 1992 Broeils', with='points'},			-- *NOT* in 2021 Ludwig paper
	{using='3:5', title='beta 1992 Broeils, corrected', with='points'},	-- is in the 2021 Ludwig paper
}

-- eqn 3.17:
-- v^2 ~ R dphi/dR
-- eqn 4.14:
-- beta^2 = r d/dr normphi(r,z=0) = g(r)
local function betacirc_for_r_eqn_4_14(r)
	return math.sqrt(g_for_r_eqn_5_3(r))
end
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_2a_NGC_1560_rotation_velocity_of_spheriod_eqn_4.14.svg",
	style = 'data lines',
	xlabel = 'r (kpc)',
	xrange = {rvec[1], rvec[#rvec]},
	ylabel = 'v/c',
	yrange = {0, .0015},
	title = 'Rotation velocity of a spheroid',
	data = {
		rvec,
		rvec:map(betacirc_for_r_eqn_4_14),
		betavec,
		r_in_kpc_1992_Broeils_table_4, 
		beta_1992_Broeils, 
		beta_corr_1992_Broeils
	},
	{using='1:2', title='β circ'},
	{using='1:3', title='β'},
	{using='4:5', title='beta 1992 Broeils', with='points'},			-- *NOT* in 2021 Ludwig paper
	{using='4:6', title='beta 1992 Broeils, corrected', with='points'},	-- is in the 2021 Ludwig paper
}
-- CHECK


local rvec = xvec * lbeta
local phivec = rvec:map(normphi_for_r_z_eq_0_eqn_5_2_a)
phivec = phivec / math.abs(phivec[1])
local rhovec = rvec:map(normrho_for_r_z_eq_0_eqn_5_2_b)
rhovec = rhovec / math.abs(rhovec[1])
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_1b_NGC_1560_mass_density_and_potential_eqn_5.2.svg",
	style = 'data lines',
	xlabel = "r (kpc)",
	xrange = {rvec[1], rvec[#rvec]},
	yrange = {-1, 1},
	title = 'Mass density and potential',
	data = {rvec, phivec, rhovec},
	{using='1:2', title=''},
	{using='1:3', title=''},
	{'0', title='', lc='rgb "grey"'},	-- zero line
}
-- CHECK.


--[[ 
text after eqn 4.5: 
beta = v(R,0) / c
R0 = 1 kpc
r = R / R0
z = Z / R0
normrho = rho / rho0
rho0 = rho(r=0,z=0)
normphi = phi / c^2

text from section 7:
"Using the variable Sérsic index profile, 
the calculated value of the absolute magnitude is Ms = −15.3,
the total luminosity is Ls = 1.02 × 108 Lsun
and the apparent magnitude is md = 12.1.
The total luminosity is calculated up to the maximum galactic radius rmax = 12.2 kpc 
that is estimated by the rotation velocity model described in the next paragraph."

from Appendix D, after eqn D11
"As an example, the surface brightness μ B of the dwarf galaxy NGC 1560 analyzed in Sect. 7 can be adjusted to the observed values listed by Broeils [32] by taking..."
so I guess that means these values *aren't* supposed to match figures 3 and 4 in Section 7 ...
but then why do they match? and why do Section 7's numbers not match?

--]]
mu0 = 22.28
alpha0 = 61.46	-- arcsec
s1 = 0.435
alpha1 = 99.05	-- arcsec
s2 = 1.144
-- r = d * (math.pi / (180. * 3600.)) * alpha
Ms = -16.5
Ls = 3.06e+8 * Lsun	-- W
d = 3.0 * 1000 -- kpc
md = 10.9
-- after eqn D12:
r0 = d * (math.pi / (180 * 3600)) * alpha0	-- r0 = 0.00089389946522976 ... 
r1 = d * (math.pi / (180 * 3600)) * alpha1	-- r1 = 0.001440623853417 ... 
-- r0 seems different from section 7 r0
-- but idk where these are used


-- eqn D11:
-- depends on s1, s2, which we only have for NGC 1560 (right?)
local function mu_for_alpha_eqn_D_11(alpha) 
	return mu0 + (
		alpha < alpha0
		and 5 / (2 * math.log(10)) * (alpha / alpha1) ^ (1 / s1)
		or 5 / (2 * math.log(10)) * (alpha0 / alpha1) ^ (1 / s1) * (1 - s2 / s1 + s2 / s1 * (alpha / alpha0) ^ (1 / s2))
	)
end


-- Figure 3 subtext:
-- "The profiles extend tot he last measurement taken at lrho = 5.13 kpc"
lrho = 5.13	-- kpc
alphae_check = alpha_for_r_d(lrho, d)
-- alphae_check = 352.71281868253 ... where the fig 3 xrange comes from
alphae = alphae_check

local alphavec = xvec * alphae

-- CHECK
-- looks like figure 3 lhs
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_3a_NGC_1560_luminosity_eqn_D11_using_appendix_D_numbers.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ (mag arcsec^-2)",
	style = 'data lines',
	title = "Luminosity profile of NGC 1560 (using appendix D adjusted values)",
	xrange = {0, alphae},
	yrange = {[3] = 'reverse'},
	data = {
		alphavec,
		alphavec:map(mu_for_alpha_eqn_D_11),
		alpha_in_arcsec_1992_Broeils_table_3,
		luminosity_1992_Broeils,
	},
	{using='1:2', title='2021 Ludwig'}, 
	{using='3:4', title='1992 Broeils', with='points'}, 
}


-- -- section 7 1st paragraph:
mu0 = 22.27
-- "alphaeff = 116.5 (reff = 1.69 kpc)" -- does this mean that the alphaeff should be derivable from the reff?
-- after D12: r0 = d (pi / (180 * 3600)) * alpha0, and then repeated for r1 and alpha1 ...

alphaeff = 116.5
reff = 1.69 -- kpc
d_check = d_for_r_alpha(reff, alphaeff)
-- d_check = 2992.1675756017 --  ... what units?
-- if the units of 'd' were 'kpc' then this would set 'd' to ~ 3.0 Mpc, which is what our Appendix D values for NGC 1560 claim
-- soo I think maybe that's a win?
-- though what about our changed alpha0 and r0 values?

-- from Appendix D, text after eqn D8 ...
-- ... is Reff the same as reff?
-- man this paper is a mess
-- "and d is the distance to the galaxy" -- at least I can recover that for NGC 1560 since the paper doesn't list its 'd' anywhere.
--alphaeff_check = (180 * 3600 / math.pi) * bs^-s * Reff / d
-- this equation, matching it with the eqn after D12 def of r for d and alpha ... 
-- ... makes it look like bs^-s Reff == r 
-- ... but what 'r' ?

alpha0 = 140.3	-- but didn't we just define an alpha0? and doesn't it work better?
r0 = 2.04 		-- kpc ... different from the Appendix D r0
-- should alpha0 be calculated from r0, or should r0 be calculated from alpha0?
alpha0_check = alpha_for_r_d(r0, d_check)	-- alpha0_check = 140.62721893491, which is close
r0_check = r_for_d_alpha(d, alpha0)		-- r0_check = 2.04058078379 ... which actually matches
-- so since the paper says alpha0 = ... (r0 = ...), 
-- maybe it is saying the parenthesis values should be derivable from the non-parenthesis variables?

-- non-parenthesis:
-- these are 4th order poly coeffs that are lin reg fitted to Sersec index profile
s0 = 0.360		-- s0 = b0
b2 = 0.0000344
b3 = -1.41e-7
b4 = -4.05e-10
d3 = -1.11e-9
d4 = 5.73e-11
bPolyOrder = 4
dPolyOrder = 4

-- parenthesis: (derivable?)
alphae = 353.0		-- derived from alpha of d and lrho (lrho for NGC 1560)
se = 0.874			-- derived from alphae using eqn D20
b1 = 0.00245		-- derived from eqn D18
d2 = -3.22e-6		-- derived from eqn D18
--[[
why can't this paper at least provide the d and b coeffs sequentially?
probably because b1 & d2 are derived using eqn D18

so how are they derived?  seems I'm having trouble verifying them.
also, where does alphae come from?  
it seems the text after eqn D14 says "sd(alphae) = se" ... 
... but both se and alphae area listed as derived values ... 
... so how do you derive se then?
--]]

lrho_check = r_for_d_alpha(d, alphae)		-- lrho_check = 5.13417688295 vs 5.13 ... close
alphae_check = alpha_for_r_d(lrho, d)		-- alphae_check = 352.71281868253 - close enough

alphae = alphae_check

local alphavec = xvec * alphae

--[[
ok going over coefficients / polynomials / constraints in Appendix D, eqns D14-D18 ...
from eqn D14 to eqn D17, they go from using 
	sd(alpha) = sum j=0..m d[j] * alpha^j	<- eqn D14
to using
	sd(alpha) = sum j=0..m d[j] * (alphae - alpha)^j	<- eqn D17
... but ... 
the d's of D14 will not match the d's of D17.
unless we were remapping only d0 or d1, as the paper does using the constraints:
	sb(alpha0) = sd(alpha0),
	sb'(alpha0) = sd'(alpha0),
	sd'(alphae) = 0
but this will only allow us to equate a fixed number of coefficients ... 3 right?
and the equation in D17 for sd(alpha) assumes that an arbitrary number of d's will match ... (n-3) I think, right?
soooo ...
it looks like the paper is reassigning the values of its d_i's without specifying it ...
yeah sure enough, recentering a polynomial changes its coefficients, though the paper acts like the coefficients will stay the same ...


compared to Fig 3a this is wrong -- inflection too far right, the value raises too high (low on the flipped graph ... why is the graph flipped in the paper?!?!?!)
but compared to Fig 2b ...
--]]

-- fig 2b subtext: "maximum radial distance 5.13 kpc" which happens to match lrho for NGC 1560 ...
-- but the graphs extend beyond r=8 ...
-- FAIL - what's the rmax of this?  shape doesn't look right either ...
local rvec = xvec * lbeta
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_2b_NGC_1560_luminosity_eqn_D11_using_section_7_numbers.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ_B (mag arcsec^-2)",	-- why is it μ_B? what does the "B" mean? "Broeils"?  Displaying the sampled data.
	style = 'data lines',
	title = "Luminosity profile of NGC 1560",
	--xrange = {rvec[1], rvec[#rvec]},	-- graph goes to 6.4, xrange goes to 8 ... hmm ...
	xrange = {0, lbeta},
	yrange = {[3] = 'reverse'},
	data = {
		rvec,
		rvec:map(function(r)
			local alpha = alpha_for_r_d(r, d)
			return mu_for_alpha_eqn_D_11(alpha)
		end),
		r_in_kpc_1992_Broeils_table_3,
		luminosity_1992_Broeils,
	},
	{using='1:2', title='2021 Ludwig'}, 
	{using='3:4', title='1992 Broeils', with='points'}, 
}




-- does eqn D20 match up with se's definition in section 7?
-- eqn D20
se_check = se_eqn_D_20(mu_for_alpha_eqn_D_11)
-- using the section 7 alpha0: se_check = 0.49270654394305
-- using the Appendix D alpha0: se_check = 0.88439918653437 ... which matches section 7
-- ... WHY DOESN'T THE SECTION 7 ALPHAE MATCH THE SECTION 7 SE?!?!?!?

-- ... annnnd no.  eqn D20's se is half of what section 7's variable listing's se is.  but substutitnig the eqn D20 se makes things worse.
-- however, if we use the alpha0 from Appendix D, then we get se = 0.88439918653437 which is close to correct. 
se = se_check

-- how about eqn D18 vs b1 and d2 in section 7?
b1_check = b1_eqn_D_18_a()
-- b1_check = 0.00091103173224744 ... wrong.  half of the given value.
-- and if we change the denom from 1 / (alphae + alpha0) to 1 / (alphae - alpha0) ... it gets a bit closer ... 
-- b1_check = 0.0021128911777981 ... close
--b1 = b1_check

d2_check = d2_eqn_D_18_b()
-- d2_check = 4.2717212743056e-07
-- ... and this is at least close 
-- one order off ...
d2 = d2_check


-- -- for the sake of the luminosity profile, section 7 still has no mention of s1, s2, or alpha1 ...
-- -- now the mu0 matches, but if I use this alpha0 then the graphs no longer match.
-- -- and using D8:
-- mu(alpha) = mu0 + (5 / (2 * math.log(10)) * (alpha / alphaeff) ** (1 / s)
-- -- ... but what is 's' ?  or would we use eqn D3 bs ~ 2s - 1/3 <=> s = (bs + 1/3) / 2
-- -- and even with eqn D3 we still aren't given a value of bs for NGC 1560 ... but we do have a way to calculate it ... what a rabbit hole
-- set output "luminosity_NGC_1560_eqn_D8.svg"
-- plot [0:350] mu(x) notitle
-- bleh .. I'm just going to assume Appendix D is used for NGC 1560.  I can only hope Appendix D also has hidden in the text somewhere the equivalent information for NGC 3198 and NGC 3115 ... because I doubt their section info does.



-- after eqn D14, verify sb(0) = s0
print('sb(0) = '..sb_for_alpha_eqn_D_17(0)..' should equal s0 = '..s0)
-- sb(0) = 0.36 should equal s0 = 0.36
-- CHECK

-- after eqn D14, verify sd(alphae) = se
print('sd(alphae) = '..sd_for_alpha_eqn_D_17(alphae)..' should equal se = '..se)
-- sd(alphae) = 0.49270654394305 should equal se = 0.49270654394305
-- CHECK

-- verify that sb(alpha0) = sd(alpha0) (eqn D15)
print('sb(alpha0) = '..sb_for_alpha_eqn_D_17(alpha0)..' should equal sd(alpha0) = '..sd_for_alpha_eqn_D_17(alpha0))
-- sb(alpha0) = 0.7872520779912 should equal sd(alpha0) = 0.61863119778044
-- close ...
-- but using the Appendix D alpha0 instead of Section 7 alpha0:
-- sb(alpha0) = 0.32668194699093 should equal sd(alpha0) = 0.36367889131819
-- it isn't so far off

-- verify that sb'(alpha0) = sd'(alpha0) (eqn D16)
print('dSb_dalpha(alpha0) = '..dSb_dalpha(alpha0)..' should equal dSd_dalpha(alpha0) = '..dSd_dalpha(alpha0))
-- dSb_dalpha(alpha0) = -0.0010347569719419 should equal dSd_dalpha(alpha0) = 0.49494316036054
-- FAIL

-- and verify that sd'(alphae) = 0
print('dSd_dalpha(alphae) = '..dSd_dalpha(alphae)..' should equal 0')
-- dSd_dalpha(alphae) = 0.49270654394305 should equal 0
-- FAIL


gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_Sersic_index_eqn_D13.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = 'Sersic index of NGC 1560',
	xrange = {0, alphae},
	data = {alphavec, alphavec:map(s_for_alpha_eqn_D_13)},
	{using='1:2', title=''},
}
-- FAIL. this looks incorrect.


-- how about using D19
local function S_eqn_D_19(alpha) 
	-- eqn D19 ... looks bad
	return math.log(alpha / alphaeff) / math.log(2 * math.log(10) / 5 * (mu_for_alpha_eqn_D_11(alpha) - mu0))	
end
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_Sersic_index_eqn_D19.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = 'Sersic index of NGC 1560',
	xrange = {0, alphae},
	data = {alphavec, alphavec:map(S_eqn_D_19)},
	{using='1:2', title=''},
}


-- R0's def:
-- in the text after eqn 4.5
-- and in the text after eqn C.3
local R0_in_kpc = 1 -- kpc
local R0_in_m = R0_in_kpc * 1000 * m_in_pc

-- eqn 4.9.a
-- also eqn C.8
-- rs = 2 * G * M / (c^2 * R0)
-- notice, typically the Schwarzschild radius rs is defined as rs = 2 M * G / c^2 ... 
-- ... so where does this R0 come into play?
local function rs_eqn_4_9_a()
	return 2 * G_over_c2_in_m_per_kg * M / R0_in_m
end

-- assuming rho0 is in kg/m^3 ... means we need to use R0_in_m^3 ... and M is in kg ...
local function lambda_eqn_4_9_b()
	return (4/3) * math.pi * R0_in_m^3 * rho0 / M
end


-- section 7 text:
lambda = 0.134
M = 1.52e+10 * Msun -- kg
rho0 = 3.31e-20	-- kg/m^3
l = 3
rs = 1.46e-6	-- kpc ... or, by eqn C.8, unitless because it is divided by R0 = 1 kpc
a = 7.19		-- kpc
b = 0.567		-- kpc
rmax = 12.2		-- kpc
alphamax = alpha_for_r_d(rmax, d)
-- alphamax = 838.81021207153
-- hmm, not the NGC 1560 graph alphamax for sure ... which is around 350
-- instead it turns out the alpha (arcsec) associated with the distance of lrho = 5.13 kpc , 
-- which is the furthest distance specified in the figure 3 subtext.


local rvec = xvec * rmax
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_5b_NGC_1560_Gravitational_potential.svg",
	style = 'data lines',
	xlabel = "r (kpc)",
	xrange = {rvec[1], rvec[#rvec]},
	format = {y = '%.2e'},
	title = 'Gravitatioanl potential of NGC 1560',
	data = {
		rvec,
		rvec:map(normphi_for_r_z_eq_0_eqn_5_2_a)
	},
	{using='1:2', title=''},
	{'0', title='', lc='rgb "grey"'},	-- zero line
}
-- fail
-- seems I need a different normphi(r,z) ...
-- if I use the old a, b, rs then this has too big of an amplitude
-- if I use the new a, b, rs then this has too small of an amplitude


-- eqn C5: lambda = 4 pi R0^3 rho0 / (3 M)
R0_in_m_check = (3 * M * lambda / (4 * math.pi * rho0)) ^ (1/3) 	-- m
-- R0_in_m_check = 3.08008569838469e+19
R0_in_kpc_check = R0_in_m_check / m_in_pc / 1000	-- kpc
-- R0_in_kpc_check = 9.98187794103888e-01 kpc
-- very close to the paper's stated R0 = 1 kpc

rs_check = rs_eqn_4_9_a()
-- rs_check = 1.4551710314407e-06 ... close to 1.46e-6

lambda_check = lambda_eqn_4_9_b()
-- lambda_check = 0.13473115517544 ... close to 0.134

-- Appendix D, after eqn D12
-- for NGC 1560, there is a d = 3 Mpc ... defined only in the Appendix D ...
-- for NGC 3198, there is a d = 9.2 Mpc
-- for NGC 3115, there is a d = 10 Mpc
-- ... so what do we use for NGC 1560?
r0_check = r_for_d_alpha(d, alpha0)

r1 = r_for_d_alpha(d, alpha1)

-- eqn D12
-- depends upon s1 and s2 to be defined
-- and those are only defined for NGC 1560 in Appendix D
-- why are the labeled equations combining two equations into one?  how do you denote them?
local function normrho_z_eq_0_eqn_D_12_a(r) 
	local alpha = alpha_for_r_d(r, d)
	return 10 ^ ((-2/5) * (mu_for_alpha_eqn_D_11(alpha) - mu0))
end
local function normrho_z_eq_0_eqn_D_12_b(r) 
	if r <= r0 then
		return math.exp(-(r / r1) ^ (1 / s1))
	end
	return math.exp(-(r0 / r1) ^ (1 / s1) * (1 - s2 / s1 + s2 / s1 * (r / r0) ^ (1 / s2)))
end


local rvec = makePow10Range(.01, rmax)

-- [[ CHECK
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_4a_NGC_1560_normalized_mass_density_eqn_D12_a.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 1560',
	log = 'x',
	data = {rvec, rvec:map(normrho_z_eq_0_eqn_D_12_a)},
	{using='1:2', title=''},
}
--]]

-- [[ CHECK
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_4a_NGC_1560_normalized_mass_density_eqn_D12_b.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 1560',
	log = 'x',
	data = {rvec, rvec:map(normrho_z_eq_0_eqn_D_12_b)},
	{using='1:2', title=''},
}
--]]

--[[ FAIL - inflection is too far to the right ... until I changed something, and now it's 100% wrong.
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_normalized_mass_density_eqn_5.2b.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 1560',
	log = 'x',
	data = {rvec, rvec:map(normrho_for_r_z_eq_0_eqn_5_2_b)},
	{using='1:2', title=''},
}
--]]


-- TODO rotation curve
-- TODO while you're here, also use the sampled rotation curve points: r_in_kpc_1992_Broeils_table_4, beta_1992_Broeils, beta_corr_1992_Broeils = matrix(_1992_Broeils_Table_4):T():unpack()
-- where tf is the formula for the rotation curve?


-- trying to get this ... can't get it fro 1560 or 3198
local function makeGalaxyWidthGraphs(name)
	do
		local r = 1
		local deltavec = xvec
		gnuplot{
			terminal = 'svg size 1024,768 background rgb "white"',
			output = "NGC_"..name.."_eqn_C9_f_lhs_rhs.svg",
			style = 'data lines',
			xlabel = "δ (kpc)",
			title = 'eqn C.9 f / lambda lhs & rhs for root finding delta based on eqn 5.2b',
			data = {
				deltavec,
				-- lhs
				deltavec:map(function(delta)
					return normrho_for_r_z_eqn_5_2_b(r, delta)
				end),
				-- rhs
				deltavec:map(function(delta)
					return eqn_C_9_rhs_over_lambda(r, delta)
				end),
			},
			{using='1:2', title='lhs'},
			{using='1:3', title='rhs'},
		}
		gnuplot{
			terminal = 'svg size 1024,768 background rgb "white"',
			output = "NGC_"..name.."_eqn_C9_f.svg",
			style = 'data lines',
			xlabel = "δ (kpc)",
			title = 'f for root finding delta based on norm rho of eqn 5.2b',
			data = {
				deltavec,
				deltavec:map(function(delta)
					return eqn_C_9_root_expr_based_on_5_2_b(r, delta)
				end),
			},
			{using='1:2', title=''},
		}
		gnuplot{
			terminal = 'svg size 1024,768 background rgb "white"',
			output = "NGC_"..name.."_eqn_C9_df.svg",
			style = 'data lines',
			xlabel = "δ (kpc)",
			title = 'df for root finding delta based on norm rho of eqn 5.2b',
			data = {
				deltavec,
				deltavec:map(function(delta)
					return dz_eqn_C_9_root_expr_based_on_5_2_b(r, delta)
				end),
			},
			{using='1:2', title=''},
		}
	end

	local rvec = xvec * rmax
	timer('plotting galacticWidth_for_r_eqn_C_9_based_on_eqn_5_2_b', function()
		gnuplot{
			terminal = 'svg size 1024,768 background rgb "white"',
			output = "NGC_"..name.."_galactic_width_eqn_C9.svg",
			style = 'data lines',
			xlabel = "r (kpc)",
			ylabel = "δ (kpc)",
			title = "Galactic width of NGC "..name,
			data = {rvec, rvec:map(galacticWidth_for_r_eqn_C_9_based_on_eqn_5_2_b)},
			{using='1:2', title=''},
		}
	end)
end

makeGalaxyWidthGraphs'1560'

-- TODO gravitational potential

for _,k in ipairs(table.keys(Gstorage)) do Gstorage[k] = nil end
print[[


NGC 3198

]]


-- Alright, after enough nonsense and contradictions, and being split between two sections (Section 7 and Appendix D) that have contradicting information, 
-- lets try starting fresh with Section 8, NGC 3198 ...
--
-- "The analysis of the spiral galaxy NGC 3198 [11,34,35] follows the same procedure established for NGC 1560"
-- Already we're off to a bad start.

mu0 = 19.49
alphaeff = 22.4
reff = 1.00	-- kpc. parenthesis.
alpha0 = 154.0
r0 = 6.87	-- kpc. parenthesis.
s0 = 0.586
b2 = -0.000433
b3 = 5.86e-8
b4 = 3.29e-9
b5 = -1.06e-11
b6 = 1.52e-13
b7 = 2.90e-15
b8 = -1.75e-17
d3 = -4.96e-7
d4 = 1.85e-9
d5 = 1.07e-11
d6 = 2.04e-14
d7 = -1.75e-16
d8 = -1.20e-18
bPolyOrder = 8
dPolyOrder = 8

-- parenthesis (derived)

--[[
how do you derive this? 
wait ... is this alpha of lrho?
why don't they call it l_e? or alpha_lp? or anything else to connect the two? gah
--]]
alphae = 316.8		
se = 1.49			-- TODO how do you derive this?
b1 = 0.0499
d2 = 0.0000181

b1_check = b1_eqn_D_18_a()	-- b1_check = 0.050231295947568 -- close
d2_check = d2_eqn_D_18_b()	-- d2_check = 1.9517551593474e-05 -- close
-- we can do this to fix the tiny gap in the piecewise graph:
b1 = b1_check
d2 = d2_check

rmax = 30.7	-- kpc
Ms = -21
Ls = 1.93e+10 * Lsun
d = 9.2 * 1000 -- kpc
md = 8.85 -- apparent magnitude
rspiral = 4.0	-- kpc
kspiral = 0.1
gamma0 = 0.95 -- kpc
gammai = gamma0	-- is this a ... function or a constant? is ri and yi functions or constants?
y0 = 8.0
v = 1.4

-- where did I get this from? 
-- somewhere near eqn 8.1 and "This population is tentatively represented by a function Y defined by a sum of n + 1 Lorentzian (Cauchy) distributions"
-- but where does it say "n=3" for NGC 3198?
-- either way, this corresponds to the number of peaks in the curves, and for NGC 3198 it has 4 peaks
YOrder = 4

alphaeff_check = alpha_for_r_d(reff, d)			-- alphaeff_check = 22.420087635554 vs 22.4 ... check
reff_check = r_for_d_alpha(d, alphaeff)			-- reff_check = 0.99910403403053 vs 1.00 ... check
alpha0_check = alpha_for_r_d(r0, d)				-- alpha0_check = 154.02600205626 vs 154.0 ... check
r0_check = r_for_d_alpha(d, alpha0)				-- r0_check = 6.8688402339599 vs 6.87 ... check
se_check = se_eqn_D_20(mu_for_alpha_eqn_D_8)	-- se_check = 1.49 vs 1.49 ... check

arctan_kspiral_in_deg = math.deg(math.atan(kspiral))
-- arctan_kspiral_in_deg = 5.7105931374996
-- theta = math.atan(ksprial) = math.rad(5.71)




rs = 1.20e-5		-- kpc ... or, by eqn C.8, unitless because it is divided by R0 = 1 kpc

-- used for normrho_for_r_z_eqn_5_2_b
a = 9.10			-- kpc - major radius bulge
b = 2.64			-- kpc - minor radius bulge
l = 3.0				-- range parameter
lambda = 0.00323	-- coefficient

rmax = 30.7			-- kpc - maximum radius
M = 1.25e+11 * Msun	-- kg - total mass of galaxy
rho0 = 6.54e-21		-- kg/m^3 - central density
Upsilon = 6.50 * UpsilonSun	-- kg/W ? ... total mass-to-light ratio

alphamax = alpha_for_r_d(rmax, d)	
-- alphamax = 688.29669041151
-- but these graphs end at about alpha=315 ...

-- Figure 6 subtext:
-- "The profiles extend to the last measurement taken at lrho = 14.1 kpc"
lrho = 14.1	-- kpc
alphae_check = alpha_for_r_d(lrho, d)
-- alphae_check = 316.12323566131	-- tada - there's where the alpha rhs range of the graph comes from.  which is also alphae.  fuckin paper.
local alphavec = xvec * alphae

-- this is the last dot on the normalized velocity graph. its xmax is rmax, which is 30.7 though.
lbeta = 29.4	-- kpc


-- not used, unless we do evaluate s as a poly of 'r', but that takes remapping the coeffs, right?
-- wait, if alphae is derived from alpha of lrho then ... re = lrho ... so ... why are they even separate variables?
-- I'm pretty sure the answer is because the other two galaxies have fully different variable labels, which smh why....
lrho_check = r_for_d_alpha(d, alphae)	-- lrho_check = 14.130185624146 vs 14.1 ... close

rs_check = rs_eqn_4_9_a()			-- rs_check = 1.1966867034874e-05 vs 1.20e-5 ... check 
lambda_check = lambda_eqn_4_9_b()	-- lambda_check = 0.0032370645736991 vs 0.00323 ... check

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_6a_NGC_3198_luminosity_eqn_D11.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ (mag arcsec^-2)",
	style = 'data lines',
	title = "Luminosity profile of NGC 3198",
	xrange = {0, alphae},
	yrange = {[3] = 'reverse'},
	data = {alphavec, alphavec:map(mu_for_alpha_eqn_D_8)},
	{using='1:2', title=''}, 
}
-- CHECK.  good.
-- though there is a hiccup in the sd->sb transition at alpha0
-- solution? re-derive the parenthesis values like b1

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_6b_NGC_3198_Sersic_index_eqn_D13.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = "Sersic index of NGC 3198",
	xrange = {0, alphae},
	data = {alphavec, alphavec:map(s_for_alpha_eqn_D_13)},
	{using='1:2', title=''}, 
}
-- CHECK.  good.
-- same as prev graph: fix the hiccup.

-- "The mass density profile extends to the last measurement taken at lrho = 1.41 kpc"
-- so looks like this mass density graph uses lrho as the rmax instead of rmax (like I think figure 3 used)
local rvec = makePow10Range(.1, lrho)
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_7a_NGC_3198_normalized_mass_density_eqn_5.2b.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3198',
	log = 'xy',
	xrange = {.1, lrho},
	yrange = {.002, 2},
	data = {rvec, rvec:map(normrho_z_eq_0_eqn_8_4)},
	{using='1:2', title=''},
}

--[[
Trying to reproduce the rotation curve data points for NGC 3198. 
The paper says they are from 1989 Begeman
I see his figure 2, whose radial distance is in arcminutes
The last entry is alpha = 11.0 arcminutes.
Also the distance he uses is 9.4 Mpc, as opposed to Ludwig's 9.2 Mpc

Using d=9.2 Mpc, r=29.4 kpc, we get alpha=659.15057648529 arcseconds
very close to 660 = 11.0 arcminutes.

Using the 1989 Begeman paper's d=9.4 Mpc r=29.4 kpc we get alpha=645.12609613454 arcseconds
aka alpha=10.752101602242 arcminutes
So d=9.4 Mpc is bad.

Actually the paper doesn't say it uses d=9.4 Mpc, it says that its table is from a previous paper that uses d=9.4 Mpc.

row #1 is arcminutes
row #2 is velocity in km/s
--]]
local _1989_Begeman_Table_2 = table{	-- r is in arcmin, v is in km/s
	{0.25,	55},
	{0.5,	92},
	{0.75,	110},
	{1.0,	123},
	{1.25,	134},
	{1.5,	142},
	{1.75,	145},
	{2.0,	147},
	{2.25,	148},
	{2.5,	152},
	{2.75,	155},
	{3.0,	156},
	{3.5,	157},
	{4.0,	153},
	{4.5,	153},
	{5.0,	154},
	{5.5,	153},
	{6.0,	150},
	{6.5,	149},
	{7.0,	148},
	{7.5,	146},
	{8.0,	147},
	{8.5,	148},
	{9.0,	148},
	{9.5,	149},
	{10.0,	150},
	{10.5,	150},
	{11.0,	149},
}
for _,row in ipairs(_1989_Begeman_Table_2) do
	-- convert from arcminutes into kpc
	row[1] = r_for_d_alpha(d, row[1] * 60)
	-- convert from km/s to m/s
	row[2] = row[2] * 1e+3 / c_in_m_per_s
end

local r_in_kpc_1989_Begeman, beta_1989_Begeman = matrix(_1989_Begeman_Table_2):T():unpack()

-- verify that the last sample point in kpc is about what lbeta in the paper is
lbeta_check = table.last(r_in_kpc_1989_Begeman)		-- lbeta_check = 29.437886716971 vs 29.4 ... close

-- now use beta_at_lbeta from the paper ... SINCE THE 2021 LUDWIG PAPER DIDN'T PROVIDE ITTTTTT>@I$#@JK$@#J:LK$@
beta_at_lbeta = table.last(beta_1989_Begeman)		-- beta_at_lbeta = 0.00049701050184525

-- Looks like beta(lbeta) does match with the sample data, fig 7b says came from Begeman source
-- but ... the value is not provided in this paper.
-- well at least the 1989 Begeman paper looks correct
local rvec, betavec = makeNormalizedRotationCurve()
local rvec = xvec * rmax
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_7b_NGC_3198_normalized_rotation_curve.svg",
	xlabel = "r (kpc)",
	ylabel = "v / c",
	style = 'data lines',
	title = 'Normalized rotation curve of NGC 3198',
	xrange = {0, rmax},
	data = {
		rvec,
		betavec,
		r_in_kpc_1989_Begeman,
		beta_1989_Begeman,
	},
	{using='1:2', title='2021 Ludwig'},
	{using='3:4', title='1989 Begeman', with='points'},
}



--[[
going through sources
2006 Cooperstock et al

Figure 1
eqn 17:
Phi = sum_n Cn exp(-kn |z|) J_0(kn r)
V = -c sum_n Cn kn exp(-kn |z|) J1(kn r)

eqn 14:
V = c N / r = c dPhi/dr
so if this 2021 Ludwig paper changes its Phi def, maybe we can use V = c dPhi/dr ?

−Cn kn			kn
0.00093352334660 0.07515079869
0.00020761839560 0.17250244090
0.00022878035710 0.27042899730
0.00009325578799 0.3684854512
0.00007945062639 0.4665911784
0.00006081834319 0.5647207491
0.00003242780880 0.6628636447
0.00003006457058 0.7610147353
0.00001687931928 0.8591712228
0.00003651365250 0.9573314522
Table 3: Curve-fitted coefficients for NGC 3198

finally an organized paper
--]]
local function v_for_r(r)
	local sum = 0
	for _,coeff in ipairs{
		{0.00093352334660, 0.07515079869},
		{0.00020761839560, 0.17250244090},
		{0.00022878035710, 0.27042899730},
		{0.00009325578799, 0.3684854512},
		{0.00007945062639, 0.4665911784},
		{0.00006081834319, 0.5647207491},
		{0.00003242780880, 0.6628636447},
		{0.00003006457058, 0.7610147353},
		{0.00001687931928, 0.8591712228},
		{0.00003651365250, 0.9573314522},
	} do
		local neg_Cn_kn, kn = table.unpack(coeff)
		sum = sum + neg_Cn_kn * J1(kn * r)
	end
	-- the eqn says times speed of light, but 2021 Ludwig graphs v/c = Lorentz beta
	-- but for reproducing 2006 Cooperstock et al, return beta*c = v
	return c_in_m_per_s * sum
end
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_3198_normalized_rotation_curve_2006_Cooperstock_Fig_3.svg",
	xlabel = "r (kpc)",
	ylabel = "v",
	style = 'data lines',
	title = 'Normalized rotation curve of NGC 3198',
	xrange = {0, rmax},
	data = {rvec, rvec:map(v_for_r) / c_in_m_per_s},
	{using='1:2', title=''},
}
--[[
CHECK
it looks just like the 2006 Cooperstock et al paper
but this graph doesnt' look like the 2021 Ludwig paper
--]]

--[[
hmm, where's the equation for delta(r)?
-- all I seem to find are equalities using delta(r) ... so then you solve the nonlinear equation?
-- 6.9, C14 ...
--]]


makeGalaxyWidthGraphs'3198'

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_3198_gravitational_potential_eqn_C7.svg",
	style = 'data lines',
	xlabel = "r (kpc)",
	format = {y = '%.2e'},
	title = 'Gravitational potential of NGC 3198',
	data = {rvec, rvec:map(normphi_for_r_z_eq_0_eqn_5_2_a)},
	{using='1:2', title=''},
}
--[[
FAIL
this does not look like the graph at all
the graph has waves like the density profile of the NGC 3198 ....
so i'm thinking it's a different galaxy width
Y(r) = 1 + sum ... this is the source of the waves in the profile ...
so where else could Y(r) influence the gravitational potential?
it is used in normrho(r,z) in eqn 8.4
sure enough, from eqn C.4 and C.9 it looks like the rhs of the root-finding is when the z in rho(r,z) is replaced with delta(r) ... and they cross-multiply lambda ...
--]]

-- if I could match phi(r), then I bet I could match d/dr phi(r)
-- and then I bet I could match the rotation curve v/c = beta(r) = sqrt(r d/dr phi) ... ??? or is it the Abel function solution?


for _,k in ipairs(table.keys(Gstorage)) do Gstorage[k] = nil end
print[[


NGC 3115

]]


mu0 = 15.21
alphaeff = 1.88
reff = 0.0911	-- kpc
alpha0 = 18.24
r0 = 0.884		-- kpc
s0 = 0.845
b2 = -0.00291
b3 = -5.37e-6
b4 = -1.13e-8
b5 = 2.84e-10
b6 = 5.38e-12
b7 = 7.62e-15
b8 = -1.07e-16
d3 = 8.94e-9
d4 = -1.42e-11
d5 = -4.82e-16
d6 = 2.04e-17
d7 = -1.31e-20
bPolyOrder = 8
dPolyOrder = 7

alphae = 983.45
se = 2.43
b1 = 0.116
d2 = -2.21e-6

rmax = 54.5	-- kpc
Ms = -21.9
Ls = 4.63 * Lsun	-- W
d = 10 * 1000	-- kpc
md = 8.08

kspiral = 0.16	-- kpc?
rspiral = 1.0	-- kpc
gamma0 = 0.4
-- section 8 says "gamma0 = gammai = ...", section 9 says "gamma0 = 0.4" ... did section 9 mean to imply that "gammai = gamm0" still / always?
gammai = gamma0
y0 = 0.4
v = 0.3	-- v? nu? upsilon? 

-- Section 9, just before eqn 9.1: "These current rings can be represented taking two terms (i = 0, 1) in Eq. (8.1)." 
YOrder = 1

rs = 6.79e-6	-- kpc
a = 5.83		-- kpc
b = 0.282		-- kpc
l = 4.3
lambda = 0.432
rmax = 54.5		-- kpc

M = 7.28e+10 * Msun	-- kg
rho0 = 5.09e-19		-- kg / m^3
Upsilon = 1.57 * UpsilonSun	-- kg / W

lrho = 47.7


alphaeff_check = alpha_for_r_d(reff, d)	-- alphaeff_check = 1.879072384911 vs 1.88 ... close
reff_check = r_for_d_alpha(d, alphaeff)	-- reff_check = 0.091144972048593 vs 0.0911 ... close
alpha0_check = alpha_for_r_d(r0, d)		-- alpha0_check = 18.233808872243 vs 18.24 ... close
r0_check = r_for_d_alpha(d, alpha0)		-- r0_check = 0.88430015434379 vs 0.884 ... close
alphae_check = alpha_for_r_d(lrho, d)		-- alphae_check = 983.88312579865 vs 983.45 ... close
se_check = se_eqn_D_20(mu_for_alpha_eqn_D_8)	-- se_check = 2.43 vs 2.43 ... check
b1_check = b1_eqn_D_18_a()				-- b1_check = 0.11562347077137 vs 0.116 ... close
d2_check = d2_eqn_D_18_b()				-- d2_check = -2.2101537775354e-06 vs -2.21e-6 ... close
d_check = d_for_r_alpha(reff, alphaeff)	-- d_check = 9995.0658771864 vs 10000 ... close

-- we can do this to fix the tiny gap in the piecewise graph, which for NGC 3115 only appears in the Sersic index graph
b1 = b1_check
d2 = d2_check

-- why is lrho the name of the dist of alphae, when r0 is the dist of alpha0 and reff is the dist of alphaeff?
-- why not call it 're' instead of 'lrho' ?
lrho_check = r_for_d_alpha(d, alphae)		-- lrho_check = 47.679001468717 vs 47.7 ... close

local alphavec = xvec * alphae
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_9a_NGC_3115_luminosity_eqn_D11.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ (mag arcsec^-2)",
	style = 'data lines',
	title = "Luminosity profile of NGC 3115",
	xrange = {0, alphae},
	yrange = {[3] = 'reverse'},
	data = {alphavec, alphavec:map(mu_for_alpha_eqn_D_8)},
	{using='1:2', title=''}, 
}
-- CHECK

local alphavec = makePow10Range(0.1, alphae)	-- where does the rmax range come from?
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_9b_NGC_3115_Sersic_index_eqn_D13.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = "Sersic index of NGC 3115",
	xrange = {0, alphae},
	yrange = {0, 2.5},
	log = 'x',
	data = {alphavec, alphavec:map(s_for_alpha_eqn_D_13)},
	{using='1:2', title=''}, 
}
-- CHECK

local rvec = makePow10Range(.01, lrho)
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_3115_normalized_mass_density_eqn_9.1a.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3115 (eqn 9.1a)',
	log = 'xy',
	xrange = {rvec[1], rvec[#rvec]},
	yrange = {1e-6, 1},
	format = {y = '%.2e'},
	-- eqn 9.1 or 8.4 depends on Y(r), which depends on gammai ... which isn't defined for NGC 3115 ... so I'm guessing
	data = {rvec, rvec:map(normrho_for_r_z_eq_0_eqn_9_1_a)},
	{using='1:2', title=''},
}
-- FAIL ...
-- xrange and yrange are good, peaks are in proper place, but peaks are wrong amplitude.
-- maybe cuz nowhere is 'gammai' specified?

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_3115_normalized_mass_density_eqn_9.1b.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3115 (eqn 9.1b)',
	log = 'xy',
	xrange = {rvec[1], rvec[#rvec]},
	yrange = {1e-6, 1},
	format = {y = '%.2e'},
	-- eqn 9.1 or 8.4 depends on Y(r), which depends on gammai ... which isn't defined for NGC 3115 ... so I'm guessing
	data = {rvec, rvec:map(normrho_for_r_z_eq_0_eqn_9_1_b)},
	{using='1:2', title=''},
}
-- FAIL, wrong yrange, but right peaks 


-- TODO normalized rotation curve

-- TODO galactic width

-- TODO gravitational potential

-- TODO normalized mass density with small increase in rs

-- TODO normalized rotation curve with small increase in rs
