#!/usr/bin/env lua128

local range = require 'ext.range'
local table = require 'ext.table'
local timer = require 'ext.timer'
local math = require 'ext.math'
local gnuplot = require 'gnuplot'

local matrix = require 'matrix'

local J1 = require 'J1'
local K = require 'K'
local E = require 'E'

local newtonRootFind = require 'newton'
local bisectRootFind = require 'bisect'

-- integration function api: integrate(f, xL, xR, n)
local integrate = require 'integrate'
--integrate.method = 'adaptlob'		-- goes slow with lua128
integrate.method = 'simpson'	

-- ivp function api: ivp(t, x, dt, f, ...)
local ivp = require 'integrate.ivp.euler'

local symmath = require 'symmath'
symmath.fixVariableNames = true
symmath.tostring = symmath.export.SingleLine


-- goes slow, don't do this unless you need to
local calcGravPotGraphs = false

-- also goes slow:
local calcRotCurveGraphs = false


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

-- for when 'd' is defined as a global
-- (useful for function concatenation)
local function r_for_alpha(alpha)
	return r_for_d_alpha(d, alpha)
end
local function alpha_for_r(r)
	return alpha_for_r_d(r, d)
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

-- eqn D8, doesn't require alpha1 (derived from r1) or alpha0 (derived from r0) or s2, like eqn D11 requires
-- but what is "s" ?
-- after eqn D1: "where s is the Sersic index"
-- but this 's' varies with alpha ... correct?
local function mu_for_alpha_eqn_D_8(alpha)
	local s = s_for_alpha_eqn_D_13(alpha)
	return mu0 + 5 / (2 * math.log(10)) * (alpha / alphaeff) ^ (1 / s)
end

-- eqn D.20
-- for NGC 1560 using eqn D.11 for mu(alpha) this messes up, but using eqn D.8 it works.
-- for NGC 3198 and 3115 it works for mu(alpha) with eqn D.8
local function se_eqn_D_20()
	return math.log(alphae / alphaeff) / math.log(2 * math.log(10) / 5 * (mu_for_alpha_eqn_D_8(alphae) - mu0))	
end

-- eqn 8.5
-- basically same as S(alpha) from eqn D13, but with (r, r0, re) swapped with (alpha, alpha0, alphae)
local function s_for_r_eqn_8_5(r)
	local alpha = alpha_for_r(r)
	--[[ 
	if r <= r0 then
		return sb_for_alpha_eqn_D_17(alpha)
	end
	if r <= re then
		return sd_for_alpha_eqn_D_17(alpha)
	end
	return se
	--]]
	-- [[ convert to alpha and evaluate accordingly
	return s_for_alpha_eqn_D_13(alpha)
	--]]
end

-- unlabeled equation?
-- Same as eqn 8.4a, 8.4b, 9.1b except with Y(r) set to 1
-- (You can set YOrder = -1 and use these others and it should also work)
-- This is the same normalized mass density as NGC 3198 and NGC 3115 except with Y(r) set to 1.
-- (Since there was no Y(r) defined for NGC 1560).
-- I can't find it anywhere in the paper that this equation is given.
local function normrho_for_r_NGC_1560(r)
	local alpha = alpha_for_r(r)
	local s = s_for_alpha_eqn_D_13(alpha)
	return math.exp(-(r / reff) ^ (1 / s))
end

--[[
eqn 8.4a
also eqn 9.1a (which I got to match eqn 9.1b)
does this mean that eqn 8.4b == eqn 9.1b ?
if mu(alpha) is of the form in eqn D.8 then everything cancels down to look like eqn. 8.4b
--]]
local function normrho_for_r_z_eq_0_eqn_8_4_a(r)
	local alpha = alpha_for_r(r)
	return Y_for_r_eqn_8_1(r) * 10 ^ (-2/5 * (mu_for_alpha_eqn_D_8(alpha) - mu0))
end

--[[
Ok so what happens when we rewrite 8.4b in terms of mu(r)?
I'm suspicious Fig 2.b uses something like this.
Especially since Fig 4.a uses a similar trick: use the next section's equation and set Y(r) to 1 ...
	normrho(r) = Y(r) 10^(-2/5 (mu(r) - mu0))
	let Y(r) = 1 ...
	mu(r) = mu0 - 5/(2 log(10)) log(normrho(r)) 
	looks similar to D.8, if log(normrho(r)) = (alpha / alphaeff) ^ (1 / s(alpha))
--]]
local function mu_for_r_eqn_8_4_a(r, normrho)
	return mu0 - 5 / (2 * math.log(10)) * math.log(normrho)
end

--[[
eqn 8.4b
... has both these equations ... and they aren't even equal ... smh
--]]
local function normrho_for_r_z_eq_0_eqn_8_4_b(r)
	--[[ can't do this without s1 and s1 of the galaxy being defined ... which it is not for NGC 3198
	local alpha = alpha_for_r(r)
	return Y_for_r_eqn_8_1(r) * 10 ^ (-2/5 * (mu_for_alpha_eqn_D_11(alpha) - mu0))
	--]]
	-- [[ looks good
	local alpha = alpha_for_r(r)
	local s = s_for_alpha_eqn_D_13(alpha)
	return 
		Y_for_r_eqn_8_1(r)
		* math.exp(-(r / reff) ^ (1 / s))
	--]]
end


--[[
seems 9.1b is not equal to 9.1a
eqn 9.1b: 
	normrho(r) = Y(r) * ...
		... exp( -(r/r0) ^ (1/sb(r)) )	for 0 <= r <= r0
		... exp( -(r/r0) ^ (1/sd(r)) )	for r0 < r
looks wrong
eqn 9.1a: 
	normrho(r) = Y(r) 10^(-2/5 (mu(r) - mu0))
looks correct ... for mu(r) defined in eqn D.8:
	mu(r) = mu0 + 5 / (2 * log(10)) * (r / reff) ^ (1 / s(r))

combining eqn 9.1a (good) with eqn D.8:
	normrho(r) = Y(r) exp(-(r/reff) ^ (1/s(r)))

now my current eqn D.8 is using the s(r(alpha)) defintiion from eqn D.13:
	s(alpha) = 
		... sb(alpha) 	for 0 <= alpha <= alpha0
		... sd(alpha) 	for alpha0 < alpha <= alphae
		... se 			for alphae < alpha

can you spot the difference?
yup, 
the good equation uses reff in the denominator 
the bad equation uses r0 in the denominator 
smh, another typo.
--]]
local function normrho_for_r_z_eq_0_eqn_9_1_b(r)
	local alpha = alpha_for_r(r)
	--[[ what the paper says ... THIS IS WRONG AND THE GRAPH DOESN'T MATCH 9.1a
	if alpha <= alpha0 then
		return Y_for_r_eqn_8_1(r) * math.exp(-(r / r0) ^ (1 / sb_for_alpha_eqn_D_17(alpha)))
	else
		return Y_for_r_eqn_8_1(r) * math.exp(-(r / r0) ^ (1 / sd_for_alpha_eqn_D_17(alpha)))
	end
	--]]
	-- [[ what the math says.  THIS IS RIGHT.  AND THE GRAPH MATCHES 9.1a
	if alpha <= alpha0 then
		return Y_for_r_eqn_8_1(r) * math.exp(-(r / reff) ^ (1 / sb_for_alpha_eqn_D_17(alpha)))
	elseif alpha <= alphae then
		return Y_for_r_eqn_8_1(r) * math.exp(-(r / reff) ^ (1 / sd_for_alpha_eqn_D_17(alpha)))
	else
		return Y_for_r_eqn_8_1(r) * math.exp(-(r / reff) ^ (1 / se))
	end
	--]]
end

-- expressions ... putting all in one place so the symbolic vars will be too
local normrho_for_r_z_eqn_5_2_b	-- also a candidate for eqn C 9 lhs over lambda
local normphi_for_r_z_eqn_5_2_a 
local dz_eqn_6_9_root_expr_based_on_5_2_b
local eqn_6_9_rhs_over_lambda		-- purely for debugging
local dr_beta_for_r_beta_f_g_eqn_5_1 		-- dbeta_dr(r, beta, f, g)
local f_for_r_eqn_5_3
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
		local a = var'a'
		local b = var'b'
		local l = var'l'
		local lambda = var'lambda'
		local mu0 = var'mu0'
	

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
		print('∂/∂r of eqn 5.2.a:', dr_normphi_eqn_5_2_a_expr)

		-- eqn 5.2.b
		-- normalized density
		-- also eqn C.4 but with lambda's definition of C.5 substituted
		-- for root finding:
		-- also equation C.9
		-- also equation C.14
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

		-- eqn 5.3 requires substitutions (eqn 4.12, text after eqn 5.2):
		local f_eqn_4_12_a_expr = f:eq( frac(3,2) * lambda * rs * normrho * r^2 )
		print('f for eqn 4.12a:', f_eqn_4_12_a_expr)
		
		local g_eqn_4_12_b_expr = g:eq(r * normphi:diff(r))
		print('g for eqn 4.12b:', g_eqn_4_12_b_expr)

		local f_z_eq_0_eqn_5_3_expr = f_eqn_4_12_a_expr:subst(normrho_eqn_5_2_b_expr, z:eq(0))()
		print('f(r,z=0) for eqn 5.3 with normrho eqn 5.2b:', f_z_eq_0_eqn_5_3_expr)
		
		local g_z_eq_0_eqn_5_3_expr = g_eqn_4_12_b_expr:subst(normphi_eqn_5_2_a_expr, z:eq(0))()
		print('g(r,z=0) for eqn 5.3 with normphi eqn 5.2a:', g_z_eq_0_eqn_5_3_expr)

		
		-- solving eqn 5.1 for dbeta/dr, but preserving f(r) and g(r) instead of analytically simplifying their replacement in, like eqn 5.3 does ...
		local dr_beta_eqn_5_1_expr = eqn_5_1()
		print('eqn 5.1 for ∂β/∂r:', dr_beta_eqn_5_1_expr)

		dr_beta_eqn_5_1_expr = (dr_beta_eqn_5_1_expr - beta * g + r * beta^2 * beta:diff(r))()
		print('eqn 5.1 for ∂β/∂r:', dr_beta_eqn_5_1_expr)

		dr_beta_eqn_5_1_expr = (dr_beta_eqn_5_1_expr / (r * (g + beta^2)))()
		print('eqn 5.1 for ∂β/∂r:', dr_beta_eqn_5_1_expr)

--[[ TODO
		local mu_for_alpha_eqn_D_11_expr = 

		-- eqn D.12
		local normrho_eqn_D_12_a_expr = 10 ^ (frac(-2, 5) * (mu_for_alpha_eqn_D_11_expr - mu0))
--]]

		-- eqn C.9 
		-- also eqn 6.9
		-- also eqn 6.11 is the same when z=0
		-- uses a normalized density function based on eqn C.4
		-- this is just setting eqn 5.2b * lambda = eqn 6.9 * lambda
		
		local eqn_6_9_rhs_over_lambda_expr = exp(-l^2 / 2)

		-- TODO eqn C.9 modified with a normalized density function based on eqn D.12


		local eqn_6_9_root_expr_based_on_5_2_b_expr = normrho_eqn_5_2_b_expr:rhs() - eqn_6_9_rhs_over_lambda_expr
		print('eqn 6.9 root using eqn 5.2.b:', eqn_6_9_root_expr_based_on_5_2_b_expr)
	
		-- TODO - Variable evaluateDerivative is commented out in Variable because it is also duplciated in Derivative's prune() rule ...
		
		local dz_eqn_6_9_root_expr = eqn_6_9_root_expr_based_on_5_2_b_expr:diff(z):prune()
		print('∂/∂z eqn 6.9 root using eqn 5.2.b:', dz_eqn_6_9_root_expr)
		
		-- assert ∂/∂x rhs == 0, 
		-- such that if width = width_lhs - width_rhs then d/dz (width) == d/dz (width_lhs)
		assert(Constant.isValue(eqn_6_9_rhs_over_lambda_expr:diff(z):prune(), 0))
		
		-- eqn C.4
		-- using the substitution of eqn C5 to replace R0 with lambda
		-- ok in eqn C 9 the z is replaced with "delta(z)" ... but delta is really just z. .. and they are solving for z and callign that variable "delta"
		normrho_for_r_z_eqn_5_2_b = symmath.export.Lua:toFunc{output={normrho_eqn_5_2_b_expr:rhs()}, input=argvars}

		normphi_for_r_z_eqn_5_2_a = symmath.export.Lua:toFunc{output={normphi_eqn_5_2_a_expr:rhs()}, input=argvars}

		eqn_6_9_rhs_over_lambda = symmath.export.Lua:toFunc{output={eqn_6_9_rhs_over_lambda_expr}, input=argvars}
		
		dz_eqn_6_9_root_expr_based_on_5_2_b = symmath.export.Lua:toFunc{output={dz_eqn_6_9_root_expr}, input=argvars}
	
		dr_beta_for_r_beta_f_g_eqn_5_1 = symmath.export.Lua:toFunc{output={dr_beta_eqn_5_1_expr:rhs()}, input={r, beta, f, g}}

		f_for_r_eqn_5_3 = symmath.export.Lua:toFunc{output={f_z_eq_0_eqn_5_3_expr:rhs()}, input={r}}
		
		g_for_r_eqn_5_3 = symmath.export.Lua:toFunc{output={g_z_eq_0_eqn_5_3_expr:rhs()}, input={r}}
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

local function galacticWidth_newton_f_lhs_based_on_eqn_5_2_b(r, z)
	return normrho_for_r_z_eqn_5_2_b(r, z)
end

local function galacticWidth_newton_f_rhs(r, z)
	return eqn_6_9_rhs_over_lambda(r, delta)
end

-- f = eqn 5.2b - exp(-l^2/2)
local function galacticWidth_newton_f_root_based_on_eqn_5_2_b(z, r) 
	return galacticWidth_newton_f_lhs_based_on_eqn_5_2_b(r, z) - galacticWidth_newton_f_rhs(r, z)
end
	
-- df/dz = d/dz (eqn 5.2b)
local function galacticWidth_newton_df_dz_based_on_eqn_5_2_b(z, r)
	return dz_eqn_6_9_root_expr_based_on_5_2_b(r, z)
end

--[[
equation C9
this is solving the nonlinear equation:
	normrho(r,z) = exp(-1/2 l^2) for z, for fixed r
where 5.2b = C.9 / lambda = normrho(r,z)

so that means we need to solve this but replace normrho(r,z) with whatever eqn we want, usu 5.2b?
nah, it looks good enough.
--]]
local function galacticWidth_for_r_eqn_6_9_based_on_eqn_5_2_b(r)
	--[[ I think the function generation is screwing up, esp for df/dz ...
	return newtonRootFind(
		galacticWidth_newton_f_root_based_on_eqn_5_2_b,
		galacticWidth_newton_df_dz_based_on_eqn_5_2_b, 
		1, r)
	--]]
	-- [[ looks good, faster than newton, but how do you determine deltamax?
	return bisectRootFind(galacticWidth_newton_f_root_based_on_eqn_5_2_b, 0, 20, nil, r) or math.nan
	--]]
end

-- all the galaxy width functions seem to be dependent on normrho(r,z) eqn 5.2.b ... but that is specific to spheroids ... 
-- I'm surprised these graphs don't only work for the normrho(r,z) defined elsewhere, like 8.4.a 
local function makeGalaxyWidthGraphs(figname, name)
	local rvec = xvec * rmax
	timer('plotting galacticWidth_for_r_eqn_6_9_based_on_eqn_5_2_b', function()
		gnuplot{
			terminal = 'svg size 1024,768 background rgb "white"',
			output = figname.."_NGC_"..name.."_galactic_width_eqn_C9.svg",
			style = 'data lines',
			xlabel = "r (kpc)",
			ylabel = "δ (kpc)",
			title = "Galactic width of NGC "..name,
			data = {rvec, rvec:map(galacticWidth_for_r_eqn_6_9_based_on_eqn_5_2_b)},
			{using='1:2', title=''},
		}
	end)
end

local function mmin_for_r_eqn_6_7(r)
	local tmp = r + rmax
	return 4 * r * rmax / (tmp * tmp)
end

-- this is the boundary epsilon, not the integral error tolerance epsilon
local boundaryEpsilon = 1e-7

local sqrt_pi = math.sqrt(math.pi)
local sqrt_2 = math.sqrt(2)

--TODO both of these integral functions can go faster if I use the galacticWidth vector instead of re-evaluating it every time

-- eqn C.17
-- normphi(r,0) is based on normrho(r,0)
local function normphi_for_r_z_eq_0_eqn_C_17(r)
	local delta = galacticWidth_for_r_eqn_6_9_based_on_eqn_5_2_b

	-- normrho -- normalized density
	-- NGC 3198 uses eqn 8.4a
	-- NGC 3115 shows 8.4a == 8.4b == 9.1b
	local normrho = normrho_for_r_z_eq_0_eqn_8_4_a

	local mmin = mmin_for_r_eqn_6_7(r)

	-- eqn C.17:
	local normphi =
		-(sqrt_2 / sqrt_pi) * (3/2) * lambda * rs * r * 
			integrate(function(m)
				local tmp = (2 - m - 2 * math.sqrt(1 - m)) / m
				return
					K(m) / (2 * math.sqrt(m * (1 - m)))
					* tmp^(3/2)
					* delta(tmp * r)
					* normrho(tmp * r)
			end, boundaryEpsilon, 1 - boundaryEpsilon)
		-(sqrt_2 / sqrt_pi) * (3/2) * lambda * rs * r
			* integrate(function(m)
				local tmp = (2 - m + 2 * math.sqrt(1 - m)) / m
				return 
					K(m) / (2 * math.sqrt(m * (1 - m)))
					* tmp^(3/2)
					* delta(tmp * r)
					* normrho(tmp * r)
			end, mmin, 1 - boundaryEpsilon)
print(r, mmin, normphi)
	return normphi
end

-- eqn C.18
local function dr_normphi_for_r_z_eq_0_eqn_C_18(r)
	local delta = galacticWidth_for_r_eqn_6_9_based_on_eqn_5_2_b
	local normrho = normrho_for_r_z_eq_0_eqn_8_4_a
	local mmin = mmin_for_r_eqn_6_7(r)
	local dr_normphi = 
		1 / (sqrt_2 * sqrt_pi) * (3/2 * lambda * rs) * 
			integrate(function(m)
				local tmp = (2 - m - 2 * math.sqrt(1 - m)) / m
				return
					tmp * tmp
					* delta(tmp * r)
					* normrho(tmp * r)
					* (
						K(m) / (2 * math.sqrt(1 - m) * (1 - math.sqrt(1 - m)))
						- E(m) / (2 * math.sqrt(1 - m) * (1 - m - math.sqrt(1 - m)))
					)
			end, boundaryEpsilon, 1 - boundaryEpsilon)
		+ 1 / (sqrt_2 * sqrt_pi) * (3/2 * lambda * rs) *
			integrate(function(m)
				local tmp = (2 - m + 2 * math.sqrt(1 - m)) / m
				return
					tmp * tmp
					* delta(tmp * r)
					* normrho(tmp * r)
					* (
						K(m) / (2 * math.sqrt(1 - m) * (1 + math.sqrt(1 - m)))
						- E(m) / (2 * math.sqrt(1 - m) * (1 - m + math.sqrt(1 - m)))
					)
			end, mmin, 1 - boundaryEpsilon)
print(r, mmin, dr_normphi)
	return dr_normphi
end


--[[
section 5 after eqn 5.3:
it says the point of evaluation being separate from r=0
reference point: beta(lbeta) = beta_l ~ 0.000267 <-> v ~ 80 km/s
sure enough, on figure 4, looks like the sample point and the graph overlap at r = lbeta.
however the graph extents go further in each direction.

it never says implicit/explicit/whatever type of integration.
it looks like they use explicit integration to solve this
i'm guessing they explicit integrated it forwards and backwards.

so this makes the smooth graphs of a spheroid rotation curve.  

how to get the bumpy ones? 
maybe replace the rho and phi with the bumpy ones, and re-derive the Abel function?
yes, sure enough, (down at the Fig__4b graphs) they replace the f(r) and g(r)

hmm, how about a better integrator? would that make things go faster?

∂β_∂r = (β * (β^2 * (1 - f) + f - g)) / (r * (g + β^2))

forward-euler:
β(t+dt) = β(t) + dt * dβ_dr(t)

backward-euler:
β(t+dt) = β(t) + dt * dβ_dr(t+dt)
β(t+dt) - dt * dβ_dr(t+dt) - β(t) = 0

-- solve the cubic for β(t+dt):
	β(t+dt)^3 * (dt * (f - 1) + r)
	+ β(t+dt)^2 * -β(t) * r
	+ β(t+dt) * ((r + dt) * g - dt * f)
	+ -β(t) * r * g
	= 0

β(t+dt) - dt * dβ_dr(t+dt) = β(t)
... assume dβ_dr(t+dt) is linear and equal to A * β(t+dt)
β(t+dt) - dt * A * β(t+dt) = β(t)
β(t+dt) * (1 - dt * A) = β(t)
β(t+dt) = β(t) / (1 - dt * A)

			1000	500		200		100
euler		good	bad				bad
midpoint	good	bad				bad
heun		bad		bad				bad
rk2alpha	good	bad				bad
rk4			bad		bad				bad
rkf45		good	bad		bad		bad
--]]
local function makeRotationCurve(dr_beta_fun, rL, rR)
	local r = lbeta
	local beta = beta_at_lbeta
	local n = 1000			-- make 1000 points
	local dr = (rR - rL) / n
	local rAndBetaVec = table()
--print('dr', dr, 'r', r, 'beta', beta)

	rAndBetaVec:insert{r, beta}

	while r > rL do
local oldbeta = beta
		r, beta = r - dr, ivp(r, beta, -dr, dr_beta_fun)
print('r='..r..' β='..beta..' ∂β/∂r='..((beta - oldbeta) / -dr))
		rAndBetaVec:insert{r, beta}
	end
	
	local r = lbeta
	local beta = beta_at_lbeta

-- [=[ only figure 4 rhs has graph data from r=0 through r=lbeta=8.29 to r=rmax=12
-- but this graph appears to be closest to fig 1 lhs, which goes from r=0 to r=lbeta=8.29
	while r < rR do
local oldbeta = beta
		r, beta = r + dr, ivp(r, beta, dr, dr_beta_fun)
print('r='..r..' β='..beta..' ∂β/∂r='..((beta - oldbeta) / dr))
		rAndBetaVec:insert{r, beta}
	end
--]=]
	rAndBetaVec:sort(function(a,b) return a[1] < b[1] end) 

	local rvec, betavec = matrix(rAndBetaVec):T():unpack()
	return rvec, betavec
end

-- eqn D11:
-- depends on s1, s2, which we only have for NGC 1560 (right?)
local function mu_for_alpha_eqn_D_11(alpha) 
	return mu0 + (
		alpha < alpha0
		and 5 / (2 * math.log(10)) * (alpha / alpha1) ^ (1 / s1)
		or 5 / (2 * math.log(10)) * (alpha0 / alpha1) ^ (1 / s1) * (1 - s2 / s1 + s2 / s1 * (alpha / alpha0) ^ (1 / s2))
	)
end

-- eqn D.12 but as a function of mu:
local function normrho_for_mu_eqn_D_12_a(mu)
	return 10 ^ ((-2/5) * (mu - mu0))
end

-- eqn D12
-- depends upon s1 and s2 to be defined
-- and those are only defined for NGC 1560 in Appendix D
-- why are the labeled equations combining two equations into one?  how do you denote them?
local function normrho_for_r_z_eq_0_eqn_D_12_a(r) 
	local alpha = alpha_for_r(r)
	local mu = mu_for_alpha_eqn_D_11(alpha)
	return normrho_for_mu_eqn_D_12_a(mu)
end
local function normrho_for_r_z_eq_0_eqn_D_12_b(r) 
	if r <= r0 then
		return math.exp(-(r / r1) ^ (1 / s1))
	end
	return math.exp(-(r0 / r1) ^ (1 / s1) * (1 - s2 / s1 + s2 / s1 * (r / r0) ^ (1 / s2)))
end

-- eqn C.15
-- how to calculate rmax
-- only used for validation right now
local function rmax_eqn_C_15()
	-- calculate rmax based on a, b, l, lambda (which itself is based a, b, l)
	return assert(bisectRootFind(
		function(r)	-- f
			return galacticWidth_newton_f_root_based_on_eqn_5_2_b(0, r)
		end,
		1,			-- xL
		100,		-- xR
		nil			-- maxiter
	))
end

-- eqn 4.12a
local function f_eqn_4_12_a(r, normrho)
	return 3/2 * lambda * rs * r * r * normrho
end

-- eqn 4.12b
local function g_eqn_4_12_b(r, dr_normphi)
	return r * dr_normphi
end

-- f from eqn 4.12a with normrho from ... 
local function f_for_r_eqn_4_12_a(r)
	--local normrho = normrho_for_r_z_eq_0_eqn_5_2_b(r)
	local normrho = normrho_for_r_z_eq_0_eqn_8_4_b(r)
	--local normrho = normrho_for_r_z_eq_0_eqn_D_12_a(r)
	--local normrho = normrho_for_r_z_eq_0_eqn_D_12_b(r)
	return f_eqn_4_12_a(r, normrho)
end

-- g from eqn 4.12b with dr_normphi from eqn C.18
local function g_for_r_eqn_4_12_b(r)
	local dr_normphi = dr_normphi_for_r_z_eq_0_eqn_C_18(r)
	return g_eqn_4_12_b(r, dr_normphi)
end

-- eqn 4.14 .a=.c, also eqn 3.17
local function beta_for_g_eqn_4_14(g)
	return math.sqrt(g)
end

-- eqn 3.17:
-- v^2 ~ R dphi/dR
-- eqn 4.14 .b=.c:
-- beta^2 = r d/dr normphi(r,z=0) = g(r)
-- (using the 5.2a definition of normphi(r,z))
local function beta_circ_for_r_eqn_4_14_using_eqn_5_3(r)
	return beta_for_g_eqn_4_14(g_for_r_eqn_5_3(r))
end

-- TODO eqn 4.14 .a=.b:
-- beta^2 = r d/dr normphi(r,0)

-- TODO eqn 4.15 
-- dbeta/dr = sqrt( f(r) (1 - beta^2) / r )
local function dr_beta_for_r_beta_eqn_4_15(r, beta, f)
	return math.sqrt(f * (1 - beta * beta) / r)
end

-- eqn 4.13:
-- (beta^2 + g(r)) r (dbeta/dr)^2 = 2 (1 - f(r)) * beta^4 + 2 (f(r) - g(r)) beta^2
-- dbeta/dr = beta sqrt( (2 (1 - f(r)) * beta^2 + 2 (f(r) - g(r))) / (r (beta^2 + g(r))) )
local function dr_beta_for_r_beta_f_g_eqn_4_13(r, beta, f, g)
	return math.sqrt( 
		(2 * (1 - f) * beta^4 + 2 * (f - g) * beta^2) 
		/ (r * (beta^2 + g)) 
	)
end

-- not used.  looks bad when used for plotting Fig 3.b
local function S_eqn_D_19(alpha) 
	return math.log(alpha / alphaeff) 
		/ math.log(2 / 5 * math.log(10) * (mu_for_alpha_eqn_D_11(alpha) - mu0))	
end


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

local function rho0_eqn_4_9_b()
	return (3 * M * lambda) / (4 * math.pi * R0_in_m^3)
end



--[[ 
text after eqn 4.5: 
beta = v(R,0) / c
R0 = 1 kpc
r = R / R0
z = Z / R0
normrho = rho / rho0
rho0 = rho(r=0,z=0)
normphi = phi / c^2
--]]




print[[


NGC 1560

]]

-- col 1 = dist from galaxy center, in arcseconds
-- col 2 = mu_B = luminosity
local _1992_Broeils_table_3 = table{
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

-- 1992 Broeis, table 4
-- col 1 is distance from galaxy center, in arcseconds
-- col 2 is velocity, in km/s
-- col 3 is velocity corrected for asymmetric drifts
local _1992_Broeils_table_4 = table{
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

-- SOOOOO
-- shitty paper
-- it never specifies to treat Y(r) = 1 for NGC 1560 ... 
-- ... and it doesn't even mention the Y(r) equations until after, in the NGC 3198 and NGC 3115 sections
-- (eqn 8.4a, 8.4b, 9.1b ...)
-- buuut some graphs only work correctly if you use this later function and this substitution
-- so I'll just set it here.
YOrder = -1

-- section 6, on spheroid fitting to NGC 1560 and then on solving the Abel equation for beta (= v/c?)
a = .373					-- kpc
b = .300					-- kpc
rs = 0.00000701				-- kpc
lbeta = 8.29				-- kpc ... this is the starting point of integration according to the text after eqn 5.3
v_at_lbeta = 80				-- km/s
beta_at_lbeta = 0.000267	-- unitless <-> v / c
M = 7.3e+10 * Msun			-- kg

--[[ take lambda from section 7:
lambda = 0.134
-- inferred from a,b,lambda
-- different from the rho0 in section 7, which uses a different a,b (and same lambda)
-- idk why I am even calculating this, I don't need it.
rho0 = rho0_eqn_4_9_b()		-- rho0 = 1.5810442712768e-19	-- kg/m^3
--]]
-- [[ alternatively, take rho0 from section 7, and calculate lambda
rho0 = 3.31e-20				-- kg/m^3
lambda = lambda_eqn_4_9_b()	-- lambda = 0.028053610392694
--]]

-- defined in Appendix D text after eqn D.11 ... but used for NGC 1560 graphs in section 6 graphs ... smh 
d = 3.0 * 1000 -- kpc


local alpha_in_arcsec_1992_Broeils_table_4, v_in_km_per_s_1992_Broeils, v_corr_in_km_per_s_1992_Broeils = matrix(_1992_Broeils_table_4):T():unpack()
local r_in_kpc_1992_Broeils_table_4 = alpha_in_arcsec_1992_Broeils_table_4:map(r_for_alpha)
local beta_corr_1992_Broeils = v_corr_in_km_per_s_1992_Broeils * 1000 / c_in_m_per_s

if calcRotCurveGraphs then
	local rvec, betavec = makeRotationCurve(function(r, beta)
			return dr_beta_for_r_beta_f_g_eqn_5_1(r, beta, f_for_r_eqn_5_3(r), g_for_r_eqn_5_3(r))
		end, 0, lbeta)
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig__1a_NGC_1560_normalized_rotation_curve_eqn_5.3.svg",
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
			beta_corr_1992_Broeils
		},
		key = 'bottom right',
		{using='1:2', title='2021 Ludwig'},
		{using='3:4', title='1992 Broeils', with='points pointtype 7'},
	}
	-- CHECK

	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig__2a_NGC_1560_rotation_velocity_of_spheriod_eqn_4.14.svg",
		style = 'data lines',
		xlabel = 'r (kpc)',
		xrange = {rvec[1], rvec[#rvec]},
		ylabel = 'v/c',
		yrange = {0, .0015},
		title = 'Rotation velocity of a spheroid',
		data = {
			rvec,
			rvec:map(beta_circ_for_r_eqn_4_14_using_eqn_5_3),
			betavec,
			r_in_kpc_1992_Broeils_table_4, 
			beta_corr_1992_Broeils
			--v_in_km_per_s_1992_Broeils * 1000 / c_in_m_per_s,
		},
		{using='1:2', title='2021 Ludwig, β-circular'},
		{using='1:3', title='2021 Ludwig, β'},
		{using='4:5', title='1992 Broeils, β-corrected', with='points pointtype 7'},	-- is in the 2021 Ludwig paper
		--{using='4:6', title='beta 1992 Broeils', with='points pointtype 7'},			-- not in the 2021 Ludwig paper
	}
	-- CHECK
end

local rvec = xvec * lbeta
-- if I have to normalize it ... then why do I call it "normrho_for_r_z..."?
local normphivec = rvec:map(normphi_for_r_z_eq_0_eqn_5_2_a) / math.abs(normphi_for_r_z_eq_0_eqn_5_2_a(rvec[1]))
local normrhovec = rvec:map(normrho_for_r_z_eq_0_eqn_5_2_b) / math.abs(normrho_for_r_z_eq_0_eqn_5_2_b(rvec[1]))
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__1b_NGC_1560_mass_density_and_potential_eqn_5.2.svg",
	style = 'data lines',
	xlabel = "r (kpc)",
	xrange = {rvec[1], rvec[#rvec]},
	yrange = {-1, 1},
	title = 'Mass density and potential',
	data = {rvec, normphivec, normrhovec},
	{using='1:2', title='ϕ / |ϕ_0|'},
	{using='1:3', title='ρ / |ρ_0|'},
	{'0', title='', lc='rgb "grey"'},
}
-- CHECK for phi
-- CHECK for rho

--[[
from Appendix D, after eqn D11
"As an example, the surface brightness μ B of the dwarf galaxy NGC 1560 analyzed in Sect. 7 can be adjusted to the observed values listed by Broeils [32] by taking..."
so I guess that means these values *aren't* supposed to match figures 3 and 4 in Section 7 ...
but then why do they match? and why do Section 7's numbers not match?
--]]

-- Appendix D, after eqn D.11:
mu0 = 22.28
alpha0 = 61.46	-- arcsec
s1 = 0.435
alpha1 = 99.05	-- arcsec
s2 = 1.144
Ms = -16.5
Ls = 3.06e+8 * Lsun	-- W
d = 3.0 * 1000 -- kpc
md = 10.9
-- Appendix D, after eqn D12:
r0 = r_for_alpha(alpha0)	-- r0 = 0.89389946522976
r1 = r_for_alpha(alpha1)	-- r1 = 1.440623853417
-- the Appendix D r0 is different from Section 7 r0 ... but so is the Appendix D vs Section 7 alpha0, so that makes sense.

-- Figure 3 subtext:
-- "The profiles extend to the last measurement taken at lrho = 5.13 kpc"
lrho = 5.13	-- kpc
alphae_check = alpha_for_r(lrho)	-- alphae_check = 352.71281868253 ... where the fig 3 xrange comes from
alphae = alphae_check

local alpha_in_arcsec_1992_Broeils_table_3, luminosity_1992_Broeils = matrix(_1992_Broeils_table_3):T():unpack()
local r_in_kpc_1992_Broeils_table_3 = alpha_in_arcsec_1992_Broeils_table_3:map(r_for_alpha)

-- CHECK
-- looks like figure 3 lhs
-- NOTICE this uses the alphae from Section 7, even though it uses the coefficients from Appendix D
local alphavec = xvec * alphae
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__3a_NGC_1560_luminosity_eqn_D11_using_appendix_D_numbers.svg",
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
	{using='3:4', title='1992 Broeils', with='points pointtype 7'}, 
}

-- -- section 7 1st paragraph:
mu0 = 22.27
alphaeff = 116.5	-- arcsec
reff = 1.69 		-- kpc, derived
alpha0 = 140.3		-- arcsec
r0 = 2.04 			-- kpc, derived.  different from the Appendix D r0

d_check = d_for_r_alpha(reff, alphaeff)		-- d_check = 2992.1675756017 --  kpc ~ 3.0 Mpc ... close to the 'd = 3.0 Mpc' for NGC 1560 in Appendix D (but not in Section 7 anywhere)
r0_check = r_for_alpha(alpha0)				-- r0_check = 2.04058078379 ... which actually matches

-- these are 4th order poly coeffs that are lin reg fitted to Sersec index profile
s0 = 0.360			-- s0 = b0
b2 = 0.0000344
b3 = -1.41e-7
b4 = -4.05e-10
d3 = -1.11e-9
d4 = 5.73e-11
bPolyOrder = 4
dPolyOrder = 4

-- parenthesis: (derivable)
alphae = 353.0		-- derived from alpha of d and lrho (lrho for NGC 1560)
se = 0.874			-- derived from alphae using eqn D20
b1 = 0.00245		-- derived from eqn D18
d2 = -3.22e-6		-- derived from eqn D18

Ms = -15.3
Ls = 1.02e+8 * Lsun	-- W
md = 12.1
rmax = 12.2	-- kpc

lrho_check = r_for_alpha(alphae)		-- lrho_check = 5.13417688295 vs 5.13 ... close
alphae_check = alpha_for_r(lrho)		-- alphae_check = 352.71281868253 - close
se_check = se_eqn_D_20()				-- se_check = 0.874 vs 0.874 - check
b1_check = b1_eqn_D_18_a()		-- b1_check = 0.0024569204655008 vs 0.00245 ... close
d2_check = d2_eqn_D_18_b()		-- d2_check = -3.2067929248811e-06 vs -3.22e-6 ... close

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


--[[
fig 2b subtext: "maximum radial distance 5.13 kpc" which happens to match lrho for NGC 1560 ...
but the graphs extend beyond r=8 ...
FAIL - what's the rmax of this?  shape doesn't look right either ...
why is it μ_B? what does the "B" mean? "Broeils"?  No, that's just the sampled data, 
 while this is the y axis label, applying also to the 2021 Ludwig graph ... 

is there ever any relation between mu(r) and rho(r) in general?
or do we always only ever have to use predefined mu(r) functions? how are these even generated?

[mu] = mag/arcsec^2
[arcsec] = [r] / [d] = [1] = unitless
mag ... is magnitude? absolute or apparent? is also unitless
so [mu] is really unitless ...

so rho is units of kg/m^3
Upsilon is units of kg/W
so rho/Upsilon = W/m^3

I guess we have ...
D.12.a: normrho(r) = 10^(-2/5 (mu(r) - mu0))
D.7: mu(R) = mu0 + 5 / (2 log (10)) * (R / (b_s^-s Reff))^(1/s)

should we say r = R / R0, reff = Reff / R0, so R/Ref = r/reff?
combined: normrho(r) = exp( -(r / (b_s^-s reff))^(1/s) )
and this fits with eqn 8.4: normrho(r) = Y(r) * exp(-(r / reff) ^ (1 / s(r)))
I guess that means R / (b_s^-s Reff) = r / reff ... does that mean "b_s^-s Reff" is a lone variable, and not two, or three?  I can't tell.

how do we invert that, to solve for mu using rho?
mu(r) = mu0 - 5/2 log10(normrho)

--]]
local rvec = xvec * lbeta
local normrhovec = rvec:map(normrho_for_r_z_eq_0_eqn_5_2_b) / math.abs(normrho_for_r_z_eq_0_eqn_5_2_b(rvec[1]))
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__2b_NGC_1560_luminosity_eqn_D11_using_section_7_numbers.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ_B (mag arcsec^-2)",	
	style = 'data lines',
	title = "Luminosity profile of NGC 1560",
	--xrange = {rvec[1], rvec[#rvec]},	-- graph goes to 6.4, xrange goes to 8 ... hmm ...
	xrange = {0, lbeta},
	yrange = {[3] = 'reverse', 31, 21},
	data = {
		rvec,
		rvec:size():lambda(function(i)
			local r = rvec[i]
			-- ok here we have normrho from eqn 5.2b, used in figure 1.b
			-- and that normrho inserted into eqn 8.4a solved for mu(r) in terms of normrho
			local normrho = normrhovec[i]
			return mu_for_r_eqn_8_4_a(r, normrho)
		end),
		r_in_kpc_1992_Broeils_table_3,
		luminosity_1992_Broeils,
	},
	{using='1:2', title='2021 Ludwig'}, 
	{using='3:4', title='1992 Broeils', with='points pointtype 7'}, 
}
-- CHECK

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


-- section 7 text:
lambda = 0.134
M = 1.52e+10 * Msun -- kg
rho0 = 3.31e-20	-- kg/m^3
l = 3
rs = 1.46e-6	-- kpc ... or, by eqn C.8, unitless because it is divided by R0 = 1 kpc
a = 7.19		-- kpc
b = 0.567		-- kpc
--rmax = 12.2	-- kpc repeated ... I think he repeats each of the 3 galaxies' rmax's def twice
alphamax = alpha_for_r(rmax)
-- alphamax = 838.81021207153
-- hmm, not the NGC 1560 graph alphamax for sure ... which is around 350
-- instead it turns out the alpha (arcsec) associated with the distance of lrho = 5.13 kpc , 
-- which is the furthest distance specified in the figure 3 subtext.

rmax_check = rmax_eqn_C_15()	-- rmax_check = 12.2202643935 vs 12.2 ... close

if calcGravPotGraphs then
	local rvec = xvec * rmax
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig__5b_NGC_1560_gravitational_potential.svg",
		style = 'data lines',
		xlabel = "r (kpc)",
		xrange = {rvec[1], rvec[#rvec]},
		format = {y = '%.2e'},
		title = 'Gravitatioanl potential of NGC 1560',
		data = {
			rvec,
			--rvec:map(normphi_for_r_z_eq_0_eqn_5_2_a)
			rvec:map(normphi_for_r_z_eq_0_eqn_C_17),	-- CLOSE
			rvec:map(dr_normphi_for_r_z_eq_0_eqn_C_18),	--
		},
		{using='1:2', title='normphi'},
		{using='1:3', title='d/dr normphi'},
		{'0', title='', lc='rgb "grey"'},	-- zero line
	}
	-- fail
	-- seems I need a different normphi(r,z) ...
	-- if I use the old a, b, rs then this has too big of an amplitude
	-- if I use the new a, b, rs then this has too small of an amplitude
end

-- eqn C5: lambda = 4 pi R0^3 rho0 / (3 M)
-- checking R0 based on NGC 1560's M, lambda, rho0:
R0_in_m_check = (3 * M * lambda / (4 * math.pi * rho0)) ^ (1/3) 	-- m	-- R0_in_m_check = 3.08008569838469e+19
R0_in_kpc_check = R0_in_m_check / m_in_pc / 1000	-- kpc	-- R0_in_kpc_check = 9.98187794103888e-01 kpc vs 1 kpc ... close
-- looks close compared to the paper's stated R0 = 1 kpc

rs_check = rs_eqn_4_9_a()			-- rs_check = 1.4551710314407e-06 vs 1.46e-6 ... close 
lambda_check = lambda_eqn_4_9_b()	-- lambda_check = 0.13473115517544 vs 0.134 ... close

-- Appendix D, after eqn D12
-- for NGC 1560, there is a d = 3 Mpc ... defined only in the Appendix D ...
-- for NGC 3198, there is a d = 9.2 Mpc
-- for NGC 3115, there is a d = 10 Mpc
-- ... so what do we use for NGC 1560?
r0_check = r_for_alpha(alpha0)		-- r0_check = 2.04058078379 vs 2.04 ... close
r1 = r_for_alpha(alpha1)


local alphavec = xvec * alphae

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__3b_NGC_1560_Sersic_index_eqn_D13.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = 'Sersic index of NGC 1560',
	xrange = {0, alphae},
	data = {alphavec, alphavec:map(s_for_alpha_eqn_D_13)},
	{using='1:2', title=''},
}
-- CHECK 
-- has a bump at alphaeff that can be fixed by re-deriving the b1, d1, d2 coeffs, just like NGC 3198


--[[
eqn D.12.a:
normrho(r) = 10^(-2/5 (mu(r) - mu0))
this looks relaly really good
but it doesn't fit the curve
maybe cuz the curve was really wrong?
but the curve is using d.12.a
does that mean the curve uses another eqn for rho?
--]]
-- [[
local normrho_for_lum_1992_Broeils = luminosity_1992_Broeils:map(normrho_for_mu_eqn_D_12_a)
normrho_for_lum_1992_Broeils = normrho_for_lum_1992_Broeils / normrho_for_lum_1992_Broeils[1] 
--]]

local rvec = makePow10Range(.01, rmax)
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__4a_NGC_1560_normalized_mass_density_eqn_5.2b.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 1560',
	log = 'x',
	data = {
		rvec,
		--rvec:map(normrho_for_r_NGC_1560),			-- works, but isn't provided in the paper
		rvec:map(normrho_for_r_z_eq_0_eqn_8_4_b),	-- works if YOrder = -1
		r_in_kpc_1992_Broeils_table_3,
		normrho_for_lum_1992_Broeils,
	},
	{using='1:2', title='2021 Ludwig'},
	{using='3:4', title='1992 Broeils', with='points pointtype 7'},
}
-- CHECK


-- TODO CAN'T DO THIS UNTIL YOU FIX NORMRHO!
--[[
Section 7 text beneath Fig 4:
	
	The left-hand graphic in Fig. 4 shows the normalized mass
	density profile ϱ(r, 0) obtained for NGC 1560 (the infor-
	mation in this figure is the same as displayed in Fig. 3).
	Now, ϱ(r, 0) can be substituted in the equation (C18) for
	∂φ(r, 0)/∂r and used to determine both the functions f(r)
	and g(r) in Abel’s equation (5.1).

and then you look at C.18, and just give up all hope. 
--]]
-- [[
if calcRotCurveGraphs then
	-- [[ this worked for figure 1 ... but no other figures
	local rvec, betavec = makeRotationCurve(
		function(r, beta)
			return dr_beta_for_r_beta_f_g_eqn_5_1(		-- has a big spike (due to f(r)?).  negative looks better but peaked in the wrong place.
			--return dr_beta_for_r_beta_f_g_eqn_4_13(	-- getting nans
			--return dr_beta_for_r_beta_eqn_4_15(r, beta, 	-- doesn't integrate across the full domain.  negative just flips it, still fails. 
				r, beta, 
				--f_for_r_eqn_5_3(r), g_for_r_eqn_5_3(r)
				f_for_r_eqn_4_12_a(r), g_for_r_eqn_4_12_b(r)
			)
		end,
		0, rmax)
	--]]
	--[[
	local rvec = xvec * rmax
	local betavec = 
	--]]
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig__4b_NGC_1560_normalized_rotation_curve.svg",
		xlabel = "r (kpc)",
		ylabel = "v / c",
		style = 'data lines',
		title = 'Normalized rotation curve of NGC 1560',
		xrange = {0, rmax},
		data = {
			rvec,
			betavec,
			r_in_kpc_1992_Broeils_table_4,
			beta_corr_1992_Broeils,
		},
		{using='1:2', title='2021 Ludwig'},
		{using='3:4', title='1992 Broeils', with='points pointtype 7'},
	}
end
--]]


makeGalaxyWidthGraphs('Fig__5a', '1560')


-- TODO gravitational potential



for _,k in ipairs(table.keys(Gstorage)) do Gstorage[k] = nil end
print[[


NGC 3198

]]

-- column 1 is in arcsec, 
-- column 2 is major mu in mag/arcsec^2
-- column 2 is minor mu in mag/arcsec^2
local _1987_Kent_table_2 = table{
	{0.00,		19.49,	19.49},
	{2.95,		19.58,	19.82},
	{5.91,		19.74,	20.24},
	{8.86,		19.91,	20.56},
	{11.82,		20.09,	20.79},
	{14.77,		20.24,	20.95},
	{17.72,		20.38,	21.03},
	{20.68,		20.51,	21.10},
	{23.63,		20.62,	21.19},
	{26.59,		20.71,	21.30},
	{29.54,		20.79,	21.40},
	{33.15,		20.88,	21.56},
	{37.05,		20.95,	21.83},
	{41.31,		21.00,	22.12},
	{45.99,		21.05,	22.44},
	{51.12,		21.10,	22.75},
	{56.77,		21.16,	22.97},
	{60*1.05,	21.25,	23.18},
	{60*1.16,	21.35,	23.46},
	{60*1.29,	21.45,	24.04},
	{60*1.43,	21.63,	24.58},
	{60*1.58,	21.89,	24.82},
	{60*1.75,	22.17,	math.nan},
	{60*1.94,	22.47,	math.nan},
	{60*2.14,	22.76,	math.nan},
	{60*2.37,	22.98,	math.nan},
	{60*2.62,	23.17,	math.nan},
	{60*2.90,	23.43,	math.nan},
	{60*3.20,	24.00,	math.nan},
	{60*3.54,	24.55,	math.nan},
	{60*3.91,	24.90,	math.nan},
	{60*4.32,	25.15,	math.nan},
	{60*4.78,	25.41,	math.nan},
	{60*5.28,	25.90,	math.nan},
}

-- col 1 is in arcminutes
-- col 2 is velocity in km/s
local _1989_Begeman_table_2 = table{
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
-- if alphae is the alpha of lrho, why don't they call it l_e? or alpha_rho? or anything else to connect the two? gah
alphae = 316.8		
se = 1.49
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

alphaeff_check = alpha_for_r(reff)		-- alphaeff_check = 22.420087635554 vs 22.4 ... check
reff_check = r_for_alpha(alphaeff)		-- reff_check = 0.99910403403053 vs 1.00 ... check
alpha0_check = alpha_for_r(r0)			-- alpha0_check = 154.02600205626 vs 154.0 ... check
r0_check = r_for_alpha(alpha0)			-- r0_check = 6.8688402339599 vs 6.87 ... check
se_check = se_eqn_D_20()				-- se_check = 1.49 vs 1.49 ... check

arctan_kspiral_in_deg = math.deg(math.atan(kspiral))
-- arctan_kspiral_in_deg = 5.7105931374996
-- theta = math.atan(ksprial) = math.rad(5.71)


rs = 1.20e-5		-- kpc ... or, by eqn C.8, unitless because it is divided by R0 = 1 kpc

-- used for normrho_for_r_z_eqn_5_2_b
a = 9.10			-- kpc - major radius bulge
b = 2.64			-- kpc - minor radius bulge
l = 3.0				-- range parameter
lambda = 0.00323	-- coefficient

M = 1.25e+11 * Msun	-- kg - total mass of galaxy
rho0 = 6.54e-21		-- kg/m^3 - central density
Upsilon = 6.50 * UpsilonSun	-- kg/W ? ... total mass-to-light ratio

alphamax = alpha_for_r(rmax)	
-- alphamax = 688.29669041151
-- but these graphs end at about alpha=315 ...

-- Figure 6 subtext:
-- "The profiles extend to the last measurement taken at lrho = 14.1 kpc"
lrho = 14.1	-- kpc
alphae_check = alpha_for_r(lrho)

-- alphae_check = 316.12323566131	-- tada - there's where the alpha rhs range of the graph comes from.  which is also alphae.  fuckin paper.
local alphavec = xvec * alphae

-- this is the last dot on the normalized velocity graph. its xmax is rmax, which is 30.7 though.
lbeta = 29.4	-- kpc


-- not used, unless we do evaluate s as a poly of 'r', but that takes remapping the coeffs, right?
-- wait, if alphae is derived from alpha of lrho then ... re = lrho ... so ... why are they even separate variables?
-- I'm pretty sure the answer is because the other two galaxies have fully different variable labels, which smh why....
lrho_check = r_for_alpha(alphae)	-- lrho_check = 14.130185624146 vs 14.1 ... close

rs_check = rs_eqn_4_9_a()			-- rs_check = 1.1966867034874e-05 vs 1.20e-5 ... check 
lambda_check = lambda_eqn_4_9_b()	-- lambda_check = 0.0032370645736991 vs 0.00323 ... check
rmax_check = rmax_eqn_C_15()		-- rmax_check = 30.72128372004 vs 30.7 ... close

local alpha_in_arcsec_1987_Kent_table_2, luminosity_1987_Kent_table_2, mu_minor_1987_Kent_table_2  = matrix(_1987_Kent_table_2):T():unpack()

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__6a_NGC_3198_luminosity_eqn_D11.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ (mag arcsec^-2)",
	style = 'data lines',
	title = "Luminosity profile of NGC 3198",
	xrange = {0, alphae},
	yrange = {[3] = 'reverse'},
	data = {
		alphavec,
		alphavec:map(mu_for_alpha_eqn_D_8),
		alpha_in_arcsec_1987_Kent_table_2,
		luminosity_1987_Kent_table_2,
	},
	{using='1:2', title='2021 Ludwig'}, 
	{using='3:4', title='1987 Kent', with='points pointtype 7'}, 
}
-- CHECK.  good.
-- though there is a hiccup in the sd->sb transition at alpha0
-- solution? re-derive the parenthesis values like b1

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__6b_NGC_3198_Sersic_index_eqn_D13.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = "Sersic index of NGC 3198",
	xrange = {0, alphae},
	data = {alphavec, alphavec:map(s_for_alpha_eqn_D_13)},
	{using='1:2', title=''}, 
}
-- CHECK.  good.

local r_in_kpc_1987_Kent_table_2 = alpha_in_arcsec_1987_Kent_table_2:map(r_for_alpha)
-- [[ this is my go-to for making normrho from lum for sampled data ...
local rho_1987_Kent_table_2 = luminosity_1987_Kent_table_2:map(normrho_for_mu_eqn_D_12_a)
local normrho_1987_Kent_table_2 = rho_1987_Kent_table_2 / rho_1987_Kent_table_2[1]
--]]

-- "The mass density profile extends to the last measurement taken at lrho = 1.41 kpc"
-- so looks like this mass density graph uses lrho as the rmax instead of rmax (like I think figure 3 used)
-- TODO what does this look like with normrho 8.4a?
local rvec = makePow10Range(.1, lrho)
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__7a_NGC_3198_normalized_mass_density_eqn_5.2b.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3198',
	log = 'xy',
	xrange = {.1, lrho},
	yrange = {.002, 2},
	data = {
		rvec,
		rvec:map(normrho_for_r_z_eq_0_eqn_8_4_b),
		r_in_kpc_1987_Kent_table_2,
		normrho_1987_Kent_table_2,
	},
	{using='1:2', title='2021 Ludwig'},
	{using='3:4', title='1987 Kent', with='points pointtype 7'},
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
--]]
local alpha_in_arcmin_1989_Begeman_table_2, v_in_km_per_s_1989_Begeman_table_2 = matrix(_1989_Begeman_table_2):T():unpack()
-- convert from arcminutes into kpc
local r_in_kpc_1989_Begeman = (alpha_in_arcmin_1989_Begeman_table_2 * 60):map(r_for_alpha)
-- convert from km/s to m/s
local beta_1989_Begeman = v_in_km_per_s_1989_Begeman_table_2 * 1000 / c_in_m_per_s

-- verify that the last sample point in kpc is about what lbeta in the paper is
lbeta_check = table.last(r_in_kpc_1989_Begeman)		-- lbeta_check = 29.437886716971 vs 29.4 ... close

-- now use beta_at_lbeta from the paper ... SINCE THE 2021 LUDWIG PAPER DIDN'T PROVIDE ITTTTTT>@I$#@JK$@#J:LK$@
beta_at_lbeta = table.last(beta_1989_Begeman)		-- beta_at_lbeta = 0.00049701050184525


--[[
Alright, we have the normalized mass density.
Now maybe we can use the "circular approximation"
of eqn 3.17 and 4.14: d/dr normphi(r) = sqrt(g(r))/r

end of Section 7:
	f(r) is directly proportional to ϱ(r) and
	g(r) depends on the integral of ϱ(r), so that these two contributions 
	are not in phase and impact differently the rotation velocity.

So can I use ϱ(r) to get f(r)?
Looks like it, from the text under eqn 5.2, nice and easy:
f(r) = (3/2) lambda rs normrho(r,0) r^2
--]]
local function f_for_r_eqn_5_3_based_on_normrho_for_r_z_eq_0_eqn_8_4_b(r)
	return 3/2 * lambda * rs * normrho_for_r_z_eq_0_eqn_8_4_b(r) * r * r
end

--[[
... and now we have to integrate to get g(r) = r d/dr normphi(r,0)
... ... ... using eqn C.18 ... smh ... we see d/dr normphi(r) is based on normrho(r) ... and ... delta(r)
and I still haven't managed to get the galaxy width root-finder to work, since the two functions don't intersect / IT HAS NO ROOTS!!!

... but wait ...
From section 8, bottom of page 186, first column:
	
	The circular velocity approximation, if applied to the gravitational
	potential shown in Fig. 8, gives a strongly oscillating velocity ...

but for circular velocity I still need a definition of normphi ...

--]]
local function g_for_r_eqn_5_3_based_on_normrho_for_r_z_eq_0_eqn_8_4_b(r)
	error'TODO'
end


if calcRotCurveGraphs then
	-- Looks like beta(lbeta) does match with the sample data, fig 7b says came from Begeman source
	-- but ... the value is not provided in this paper.
	-- well at least the 1989 Begeman paper looks correct
	local rvec, betavec = makeRotationCurve(function(r, beta)
			return dr_beta_for_r_beta_f_g_eqn_5_1(
				r, beta, 
				--f_for_r_eqn_5_3(r), g_for_r_eqn_5_3(r)
				f_for_r_eqn_4_12_a(r), g_for_r_eqn_4_12_b(r) 
			)
		end,
		0, rmax)
	--local rvec = xvec * rmax	-- equivalent, since rmax = 30.7, and lbeta = 29.4
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig__7b_NGC_3198_normalized_rotation_curve.svg",
		xlabel = "r (kpc)",
		ylabel = "v / c",
		style = 'data lines',
		title = 'Normalized rotation curve of NGC 3198',	-- why is this called "normalized" rotation curve?
		xrange = {0, rmax},
		data = {
			rvec,
			betavec,
			r_in_kpc_1989_Begeman,
			beta_1989_Begeman,
		},
		{using='1:2', title='2021 Ludwig'},
		{using='3:4', title='1989 Begeman', with='points pointtype 7'},
	}
end
-- FAIL
-- this will remain a FAIL until 8.b grav pot and d/dr are fixed


--[[
hmm, where's the equation for delta(r)?
-- all I seem to find are equalities using delta(r) ... so then you solve the nonlinear equation?
-- 6.9, C14 ...

I'm trying to generalize the rotation curve generation too quickly in this function:
--]]
makeGalaxyWidthGraphs('Fig__8a', '3198')


-- eqn 4.14
local function dr_normphi_for_r_g_eqn_4_14(r, g)
	return math.sqrt(g) / r
end

local rvec = xvec * rmax

--[[ From section 8, bottom of page 186, first column:
The circular velocity approximation, if applied to the gravitational
potential shown in Fig. 8, gives a strongly oscillating velocity ...

so I think that means I should use the circular velocity profile,
from eqn 3.17 and 4.14, which says beta(r) = sqrt(g(r)) = sqrt(r d/dr normphi(r))
which would mean that d/dr normphi(r) = sqrt(g(r))/r
... for ... i'm guessing ... g(r) being the circular velocity profile, in eqn 5.3 (as eqn 4.14 and eqn 3.17 use)?
--]]	
local dr_normphi_vec = rvec:map(function(r)
	return dr_normphi_for_r_g_eqn_4_14(r, g_for_r_eqn_5_3(r))
end)

if calcGravPotGraphs then	-- because it's so slow
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig__8b_NGC_3198_gravitational_potential_eqn_C7.svg",
		style = 'data lines',
		xlabel = "r (kpc)",
		format = {y = '%.2e'},
		title = 'Gravitational potential of NGC 3198',
		data = {
			rvec,
			rvec:map(normphi_for_r_z_eq_0_eqn_C_17),	-- CLOSE
			rvec:map(dr_normphi_for_r_z_eq_0_eqn_C_18),	--
		},
		{using='1:2', title='normphi'},
		{using='1:3', title='d/dr normphi'},
		{'0', title='', lc='rgb "grey"'},	-- zero line
	}
	-- if I could match phi(r), then I bet I could match d/dr phi(r)
	-- and then I bet I could match the rotation curve v/c = beta(r) = sqrt(r d/dr phi) ... ??? or is it the Abel function solution?
end


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
	ylabel = "v / c",
	style = 'data lines',
	title = 'Normalized rotation curve of NGC 3198',
	xrange = {0, rmax},
	data = {
		rvec, 
		rvec:map(v_for_r) / c_in_m_per_s,
		r_in_kpc_1989_Begeman,
		beta_1989_Begeman,
	},
	{using='1:2', title='2006 Cooperstock'},
	{using='3:4', title='1989 Begeman', with='points pointtype 7'},
}
--[[
CHECK
it looks just like the 2006 Cooperstock et al paper
but this graph doesnt' look like the 2021 Ludwig paper
--]]




for _,k in ipairs(table.keys(Gstorage)) do Gstorage[k] = nil end
print[[


NGC 3115

]]

-- 1st col is alpha^(1/4), in arcseconds, so ^4 it to get the arcseconds
-- 2nd col is the final processing major (axis?) luminosity in mag/arcsec^2 ... this is the one that the 2021 Ludwig paper uses
-- 3rd col is the final processing minor (axis?), which I don't think the paper uses ... and isn't complete anyways ...
local _1987_Capaccioli_et_al_table_9 = table{
	{0.65,	15.34,	15.35},
	{0.70,	15.36,	15.38},
	{0.75,	15.37,	15.42},
	{0.80,	15.39,	15.48},
	{0.85,	15.42,	15.57},
	{0.90,	15.49,	15.70},
	{0.95,	15.61,	15.89},
	{1.00,	15.75,	16.13},
	{1.05,	15.92,	16.40},
	{1.10,	16.09,	16.67},
	{1.15,	16.26,	16.93},
	{1.20,	16.43,	17.17},
	{1.25,	16.59,	17.41},
	{1.30,	16.76,	17.62},
	{1.35,	16.93,	17.82},
	{1.40,	17.11,	18.01},
	{1.45,	17.27,	18.20},
	{1.50,	17.42,	18.40},
	{1.55,	17.58,	18.61},
	{1.60,	17.71,	18.82},
	{1.65,	17.88,	19.01},
	{1.70,	17.97,	19.21},
	{1.75,	18.06,	19.43},
	{1.80,	18.18,	19.65},
	{1.85,	18.29,	19.88},
	{1.90,	18.36,	20.12},
	{1.95,	18.42,	20.36},
	{2.00,	18.50,	20.57},
	{2.05,	18.63,	20.79},
	{2.10,	18.82,	21.04},
	{2.15,	19.01,	21.23},
	{2.20,	19.19,	21.40},
	{2.25,	19.35,	21.56},
	{2.30,	19.51,	21.72},
	{2.35,	19.65,	21.91},
	{2.40,	19.80,	22.07},
	{2.45,	19.95,	22.21},
	{2.50,	20.08,	22.39},
	{2.55,	20.20,	22.59},
	{2.60,	20.34,	22.71},
	{2.65,	20.51,	22.84},
	{2.70,	20.67,	23.00},
	{2.75,	20.75,	math.nan},
	{2.80,	20.85,	23.27},
	{2.85,	21.01,	math.nan},
	{2.90,	21.15,	23.57},
	{2.95,	21.26,	math.nan},
	{3.00,	21.46,	23.87},
	{3.10,	21.83,	24.13},
	{3.20,	22.21,	24.40},
	{3.30,	22.58,	24.73},
	{3.40,	22.94,	25.02},
	{3.50,	23.30,	25.31},
	{3.60,	23.65,	25.60},
	{3.70,	23.98,	25.84},
	{3.80,	24.31,	26.10},
	{3.90,	24.63,	26.37},
	{4.00,	24.96,	26.63},
	{4.10,	25.29,	26.97},
	{4.20,	25.65,	27.29},
	{4.30,	26.01,	27.54},
	{4.40,	26.35,	27.81},
	{4.50,	26.63,	28.06},
	{4.60,	26.90,	28.31},
	{4.70,	27.21,	28.54},
	{4.80,	27.50,	28.76},
	{4.90,	27.78,	28.97},
	{5.00,	28.03,	29.15},
	{5.10,	28.27,	math.nan},
	{5.20,	28.49,	math.nan},
	{5.30,	28.69,	math.nan},
	{5.40,	28.90,	math.nan},
	{5.50,	29.14,	math.nan},
	{5.60,	29.42,	math.nan},
}

-- col 1 is radial dist in arcsec
-- col 2 is the number of samples
-- col 3 is velocity in km/s relative to the sun
-- col 4 is velocity error range in km/s 
-- col 5 is the velocity relative to NGC 3114
local _1980_Rubin_et_al_table_2 = table{
	-- NE data:
	{-95.6,	3,	910,	26,	248},
	{-90.8,	6,	912,	14,	250},
	{-85.7,	6,	894,	8,	232},
	{-80.6,	8,	967,	23,	305},
	{-75.5,	8,	971,	21,	309},
	{-70.3,	11,	938,	15,	276},
	{-65.0,	15,	938,	11,	276},
	{-59.9,	18,	931,	11,	269},
	{-54.9,	21,	937,	15,	275},
	{-50.0,	27,	921,	10,	259},
	{-44.9,	28,	914,	9,	252},
	{-39.8,	28,	915,	10,	253},
	{-34.6,	28,	923,	10,	261},
	{-29.5,	28,	923,	10,	261},
	{-24.4,	27,	935,	9,	273},
	{-19.3,	21,	931,	9,	269},
	{-14.1,	19,	916,	9,	254},
	{-9.1,	9,	866,	12,	204},
	{-3.8,	4,	766,	28,	104},
	-- center
	{0,		2,	685,	6,	23},
	-- SW data
	{3.8,	4,	583,	21,	79},
	{9.2,	10,	463,	26,	199},
	{14.1,	17,	415,	16,	247},
	{19.4,	23,	396,	12,	266},
	{24.3,	26,	394,	11,	268},
	{29.5,	27,	389,	9,	273},
	{34.7,	28,	390,	7,	274},
	{39.8,	28,	390,	7,	272},
	{44.9,	28,	398,	9,	264},
	{50.0,	27,	394,	8,	268},
	{55.0,	25,	394,	8,	268},
	{60.0,	20,	392,	11,	270},
	{65.1,	20,	393,	13,	269},
	{70.3,	18,	402,	17,	260},
	{75.4,	18,	400,	24,	262},
	{80.4,	14,	394,	23,	265},
	{85.5,	10,	376,	33,	286},
	{90.9,	4,	391,	17,	271},
	{96.0,	4,	401,	21,	261},
}



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
-- eqn 9.1 or 8.4 depends on Y(r), which depends on gammai ... which isn't defined for NGC 3115 
-- section 8 says "gamma0 = gammai = ...", section 9 says "gamma0 = 0.4" 
-- ... did section 9 mean to imply that "gammai = gamm0" still / always / for NGC 3115 as well?
gammai = gamma0

-- the paper says:
--y0 = 0.4
-- but the ridges on Fig 10.a don't rise nearly high enough (there's supposed to be concavity but with y0=0.4 there isn't)
-- but with y0 = 4.0 it looks right: ... typo anyone?
y0 = 4.0
-- I'm guessing this can be validated with some least-squares calculations but meh, I'm not going to do it.  I'm doing enough work as it is.

v = 0.3	-- v? nu? upsilon? 

-- Section 9, just before eqn 9.1: "These current rings can be represented taking two terms (i = 0, 1) in Eq. (8.1)." 
-- Likewise if you set this to 4 then you get extra peaks that shouldn't be there.
YOrder = 1

rs = 6.79e-6	-- kpc
a = 5.83		-- kpc
b = 0.282		-- kpc
l = 4.3
lambda = 0.432

M = 7.28e+10 * Msun	-- kg
rho0 = 5.09e-19		-- kg / m^3
Upsilon = 1.57 * UpsilonSun	-- kg / W

lrho = 47.7

lbeta = 4.64	-- kpc .. fig 10.b subtext


alphaeff_check = alpha_for_r(reff)		-- alphaeff_check = 1.879072384911 vs 1.88 ... close
reff_check = r_for_alpha(alphaeff)		-- reff_check = 0.091144972048593 vs 0.0911 ... close
alpha0_check = alpha_for_r(r0)			-- alpha0_check = 18.233808872243 vs 18.24 ... close
r0_check = r_for_alpha(alpha0)			-- r0_check = 0.88430015434379 vs 0.884 ... close
alphae_check = alpha_for_r(lrho)		-- alphae_check = 983.88312579865 vs 983.45 ... close
se_check = se_eqn_D_20()				-- se_check = 2.43 vs 2.43 ... check
b1_check = b1_eqn_D_18_a()				-- b1_check = 0.11562347077137 vs 0.116 ... close
d2_check = d2_eqn_D_18_b()				-- d2_check = -2.2101537775354e-06 vs -2.21e-6 ... close
d_check = d_for_r_alpha(reff, alphaeff)	-- d_check = 9995.0658771864 vs 10000 ... close

-- we can do this to fix the tiny gap in the piecewise graph, which for NGC 3115 only appears in the Sersic index graph
b1 = b1_check
d2 = d2_check

-- why is lrho the name of the dist of alphae, when r0 is the dist of alpha0 and reff is the dist of alphaeff?
-- why not call it 're' instead of 'lrho' ?
lrho_check = r_for_alpha(alphae)		-- lrho_check = 47.679001468717 vs 47.7 ... close
rmax_check = rmax_eqn_C_15()			-- rmax_check = 54.555763851156 vs 54.5 ... close


local alpha_in_arcsec_qtrt_1987_Capaccioli_et_al_table_9, luminosity_1987_Capaccioli_et_al_table_9 = matrix(_1987_Capaccioli_et_al_table_9):T():unpack()
local alpha_in_arcsec_1987_Capaccioli_et_al_table_9 = alpha_in_arcsec_qtrt_1987_Capaccioli_et_al_table_9:map(function(x) return x*x*x*x end)


local alphavec = xvec * alphae
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__9a_NGC_3115_luminosity_eqn_D11.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ (mag arcsec^-2)",
	style = 'data lines',
	title = "Luminosity profile of NGC 3115",
	xrange = {0, alphae},
	yrange = {[3] = 'reverse'},
	data = {
		alphavec,
		alphavec:map(mu_for_alpha_eqn_D_8),
		alpha_in_arcsec_1987_Capaccioli_et_al_table_9,
		luminosity_1987_Capaccioli_et_al_table_9,
	},
	{using='1:2', title='2021 Ludwig'}, 
	{using='3:4', title='1987 Capaccioli et al', with='points pointtype 7'}, 
}
-- CHECK

local alphavec = makePow10Range(0.1, alphae)	-- where does the rmax range come from?
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig__9b_NGC_3115_Sersic_index_eqn_D13.svg",
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

local r_in_kpc_1987_Capaccioli_et_al_table_9 = alpha_in_arcsec_1987_Capaccioli_et_al_table_9:map(r_for_alpha)
-- [[ I guess this is the way to convert from mu to rho ?  seems to work in fig 4 and here
local rho_1987_Capaccioli_et_al_table_9 = luminosity_1987_Capaccioli_et_al_table_9:map(normrho_for_mu_eqn_D_12_a)
local normrho_1987_Capaccioli_et_al_table_9 = rho_1987_Capaccioli_et_al_table_9 / rho_1987_Capaccioli_et_al_table_9[1]
--]]

-- two graph contestants for the normalized mass density, each with their own set of problems:

local rvec = makePow10Range(.01, lrho)

--[[ CHECK
with paper's "y0 = 0.4" this is a FAIL ...
xrange and yrange are good, peaks are in proper place, but peaks are wrong amplitude.
maybe cuz nowhere is 'gammai' specified?
but with my guess of a typo "y0 = 4.0" this is a CHECK -- looks correct.
--]]
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_10a_NGC_3115_normalized_mass_density_eqn_8.4b.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3115 (eqn 8.4b)',
	log = 'xy',
	xrange = {rvec[1], rvec[#rvec]},
	--yrange = {1e-6, 1},
	format = {y = '%.2e'},
	data = {
		rvec,
		rvec:map(normrho_for_r_z_eq_0_eqn_8_4_b),
		r_in_kpc_1987_Capaccioli_et_al_table_9,
		normrho_1987_Capaccioli_et_al_table_9,
	},
	{using='1:2', title='2021 Ludwig'},
	{using='3:4', title='1987 Capaccioli et al', with='points pointtype 7'},
}

--[=[ these two should match above.
-- CHECK
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_10a_NGC_3115_normalized_mass_density_eqn_8.4a.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3115 (eqn 8.4a)',
	log = 'xy',
	xrange = {rvec[1], rvec[#rvec]},
	--yrange = {1e-6, 1},
	format = {y = '%.2e'},
	data = {
		rvec,
		rvec:map(normrho_for_r_z_eq_0_eqn_8_4_a),
		r_in_kpc_1987_Capaccioli_et_al_table_9,
		normrho_1987_Capaccioli_et_al_table_9,
	},
	{using='1:2', title='2021 Ludwig'},
	{using='3:4', title='1987 Capaccioli et al', with='points pointtype 7'},
}


-- [[ CHECK
-- THIS FAILS WITH THE INCORRECT EQUATION IN 9.1b
-- unless of course I assumed wrong about how to fix 9.1a -- THEN IT WORKS
-- either way, now it looks right, now that I replaced the 9.1b denominator of r0 with reff  (and added a final piecewise boundary for s(r) = se)
-- wrong yrange, but right peaks 
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_10a_NGC_3115_normalized_mass_density_eqn_9.1b.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3115 (eqn 9.1b)',
	log = 'xy',
	xrange = {rvec[1], rvec[#rvec]},
	--yrange = {1e-6, 1},
	format = {y = '%.2e'},
	-- eqn 9.1 or 8.4 depends on Y(r), which depends on gammai ... which isn't defined for NGC 3115 ... so I'm guessing
	data = {
		rvec,
		rvec:map(normrho_for_r_z_eq_0_eqn_9_1_b),
		r_in_kpc_1987_Capaccioli_et_al_table_9,
		normrho_1987_Capaccioli_et_al_table_9,
	},
	{using='1:2', title='2021 Ludwig'},
	{using='3:4', title='1987 Capaccioli et al', with='points pointtype 7'},
}
--]]
--]=]

makeGalaxyWidthGraphs('Fig_11a', '3115')

local alpha_in_arcsec_1980_Rubin_et_al_table_2, _, _, _, v_in_km_per_s_1980_Rubin_et_al_table_2 = matrix(_1980_Rubin_et_al_table_2):T():unpack()
local r_in_kpc_1980_Rubin_et_al_table_2 = alpha_in_arcsec_1980_Rubin_et_al_table_2:map(r_for_alpha)
local beta_1980_Rubin_et_al_table_2 = v_in_km_per_s_1980_Rubin_et_al_table_2 * 1000 / c_in_m_per_s

lbeta_check = table.last(r_in_kpc_1980_Rubin_et_al_table_2)
beta_at_lbeta = table.last(beta_1980_Rubin_et_al_table_2)

if calcGravPotGraphs then	-- because it's so slow
	local rvec = xvec * rmax
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig_11b_NGC_3115_gravitational_potential_eqn_C7.svg",
		style = 'data lines',
		xlabel = "r (kpc)",
		format = {y = '%.2e'},
		title = 'Gravitational potential of NGC 3115',
		data = {
			rvec,
			rvec:map(normphi_for_r_z_eq_0_eqn_C_17),
			rvec:map(dr_normphi_for_r_z_eq_0_eqn_C_18),
		},
		{using='1:2', title='normphi'},
		{using='1:3', title='d/dr normphi'},
		{'0', title='', lc='rgb "grey"'},	-- zero line
	}
end

if calcRotCurveGraphs then
	local rvec = xvec * lbeta
	local rvec, betavec = makeRotationCurve(function(r, beta)
			return dr_beta_for_r_beta_f_g_eqn_5_1(
				r, beta, 
				--f_for_r_eqn_5_3(r), g_for_r_eqn_5_3(r)
				f_for_r_eqn_4_12_a(r), g_for_r_eqn_4_12_b(r)
			)
		end,
		0, lbeta)
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig_10b_NGC_3115_normalized_rotation_curve.svg",
		xlabel = "r (kpc)",
		ylabel = "v / c",
		style = 'data lines',
		title = 'Normalized rotation curve of NGC 3115',
		xrange = {0, lbeta},
		data = {
			rvec,
			betavec,
			r_in_kpc_1980_Rubin_et_al_table_2,
			beta_1980_Rubin_et_al_table_2,
		},
		{using='1:2', title='2021 Ludwig'},
		{using='(abs($3)):4', title='1980 Rubin et al', with='points pointtype 7'},
	}
end



-- fig 12 has adjusted schwarzschild radius
rs = 7.00e-6

-- normalized mass density with small increase in rs
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "Fig_12a_NGC_3115_normalized_mass_density_eqn_8.4a_rs_adj.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3115 (eqn 8.4a)',
	log = 'xy',
	xrange = {rvec[1], rvec[#rvec]},
	--yrange = {1e-6, 1},
	format = {y = '%.2e'},
	data = {
		rvec,
		rvec:map(normrho_for_r_z_eq_0_eqn_8_4_a),
		r_in_kpc_1987_Capaccioli_et_al_table_9,
		normrho_1987_Capaccioli_et_al_table_9,
	},
	{using='1:2', title='2021 Ludwig'},
	{using='3:4', title='1987 Capaccioli et al', with='points pointtype 7'},
}

-- normalized rotation curve with small increase in rs
if calcRotCurveGraphs then	
	local rvec = xvec * lbeta
	local rvec, betavec = makeRotationCurve(function(r, beta)
			return dr_beta_for_r_beta_f_g_eqn_5_1(
				r, beta,
				--f_for_r_eqn_5_3(r), g_for_r_eqn_5_3(r)
				f_for_r_eqn_4_12_a(r), g_for_r_eqn_4_12_b(r)
			)
		end,
		0, lbeta)
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig_12b_NGC_3115_normalized_rotation_curve_rs_adj.svg",
		xlabel = "r (kpc)",
		ylabel = "v / c",
		style = 'data lines',
		title = 'Normalized rotation curve of NGC 3115',
		xrange = {0, lbeta},
		data = {
			rvec,
			betavec,
			r_in_kpc_1980_Rubin_et_al_table_2,
			beta_1980_Rubin_et_al_table_2,
		},
		{using='1:2', title='2021 Ludwig'},
		{using='(abs($3)):4', title='1980 Rubin et al', with='points pointtype 7'},
	}
end

