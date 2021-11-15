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

--[[
equation 5.2.a
also equation C.7
φ = ϕ/c^2 ... WHY DO YOU HAVE TO USE TWO DIFFERENT PHIS!?!?!?!?
r, rs, a, b, z are in kpc ... or equivalently normalized by R0 = 1 kpc ...
--]]
local function normphi_for_r_z_eqn_5_2_a(r,z)
	return -rs / (2 * math.sqrt(r*r + (a + math.sqrt(b*b + z*z))^2))
end
local function normphi_for_r_z_eq_0_eqn_5_2_a(r)
	return normphi_for_r_z_eqn_5_2_a(r,0)
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
local eqn_C_9_root_expr_based_on_5_2_b, dz_eqn_C_9_root_expr_based_on_5_2_b	-- internal to the newton root-finder
local galacticWidth_for_r_eqn_C_9_based_on_eqn_5_2_b
local eqn_C_9_rhs_over_lambda		-- purely for debugging
local dr_beta_for_r_eqn_5_3
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
		local normphi_eqn_5_2_a_expr = normphi:eq(- rs / (2 * sqrt(r^2 + (a + tmp)^2)))
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
		local fdef_eqn_5_3 = f:eq( frac(3,2) * lambda * rs * normrho * r^2 )
		print('f for eqn 5.3:', fdef_eqn_5_3)

		fdef_eqn_5_3 = fdef_eqn_5_3:subst(normrho_eqn_5_2_b_expr, z:eq(0))()
		print('f for eqn 5.3:', fdef_eqn_5_3)

		local gdef_eqn_5_3 = g:eq(r * normphi:diff(r))
		print('g for eqn 5.3:', gdef_eqn_5_3)
		
		gdef_eqn_5_3 = gdef_eqn_5_3:subst(normphi_eqn_5_2_a_expr, z:eq(0))()
		print('g for eqn 5.3:', gdef_eqn_5_3)

		-- deriving eqn 5.3 ...
		local eqn_5_3 = eqn_5_1:clone()
		print('eqn 5.1 with substitutions:', eqn_5_3) 
		
		eqn_5_3 = eqn_5_3()
		print('eqn 5.1 with substitutions:', eqn_5_3) 
		
		eqn_5_3 = (eqn_5_3 - beta * g + r * beta^2 * beta:diff(r))()
		print('eqn 5.1 with substitutions:', eqn_5_3) 
		
		eqn_5_3 = (eqn_5_3 / (r * (g + beta^2)))()
		print('eqn 5.1 with substitutions:', eqn_5_3) 
		
		eqn_5_3 = eqn_5_3:subst(fdef_eqn_5_3, gdef_eqn_5_3)
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
		
		eqn_C_9_rhs_over_lambda = symmath.export.Lua:toFunc{output={eqn_C_9_rhs_over_lambda_expr}, input=argvars}
		eqn_C_9_root_expr_based_on_5_2_b = symmath.export.Lua:toFunc{output={eqn_C_9_root_expr_based_onn_5_2_b_expr}, input=argvars}
		dz_eqn_C_9_root_expr_based_on_5_2_b = symmath.export.Lua:toFunc{output={dz_eqn_C_9_root_expr}, input=argvars}
	
		dr_beta_for_r_eqn_5_3 = symmath.export.Lua:toFunc{output={eqn_5_3:rhs()}, input={r, beta}}
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
--[[
section 5 after eqn 5.3:
reference point: beta(lbeta) = beta_l ~ 0.000267 <-> v ~ 80 km/s
sure enough, on figure 4, looks like the sample point and the graph overlap at r = lbeta.
however the graph extents go further in each direction.
i'm guessing they explicit integrated it forwards and backwards.
--]]
do
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

--[[ only figure 4 rhs has graph data from lbeta=8.29 to rmax=12
-- but this graph appears to be closest to fig 1 lhs
	while r < 10 do
		local dbeta_dr = dr_beta_for_r_eqn_5_3(r, beta)
		beta = beta + dbeta_dr * dr
		r = r + dr
		rAndBetaVec:insert{r, beta}
	end
--]]	
	rAndBetaVec:sort(function(a,b) return a[1] < b[1] end) 

	local rvec, betavec = matrix(rAndBetaVec):T():unpack()

	-- now try plotting it ...
	-- should be figure 1 lhs, figure 2 lhs, or figure 4 rhs ... yeah they all have the same label
	-- CHECK
	-- this is figure 1 lhs
	gnuplot{
		terminal = 'svg size 1024,768 background rgb "white"',
		output = "Fig_1a_NGC_1560_normalized_rotation_curve_eqn_5.3.svg",
		style = 'data lines',
		xlabel = "r (kpc)",
		ylabel = 'v/c',
		title = "Normalized rotation curve of NGC 1560",
		xrange = {rvec[1], rvec[#rvec]},
		yrange = {(table.inf(betavec)), (table.sup(betavec))},
		data = {rvec, betavec},
		{using='1:2', title=''},
	}
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
	data = {alphavec, alphavec:map(mu_for_alpha_eqn_D_11)},
	{using='1:2', title=''}, 
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
local rvec = xvec * 6.4	-- xvec * lrho	

-- FAIL - what's the rmax of this?  shape doesn't look right either ...
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_luminosity_eqn_D11_using_section_7_numbers.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ_B (mag arcsec^-2)",	-- why is it μ_B? what does the "B" mean?
	style = 'data lines',
	title = "Luminosity profile of NGC 1560",
	--xrange = {rvec[1], rvec[#rvec]},	-- graph goes to 6.4, xrange goes to 8 ... hmm ...
	xrange = {0, 8},
	yrange = {[3] = 'reverse'},
	data = {rvec, rvec:map(function(r)
		local alpha = alpha_for_r_d(r, d)
		return mu_for_alpha_eqn_D_11(alpha)
	end)},
	{using='1:2', title=''}, 
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


M = 1.52e+10 * Msun -- kg
lambda = 0.134
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

local function normrho_for_r_z_eq_0_eqn_5_2_b(r)
	return normrho_for_r_z_eqn_5_2_b(r,0)
end

-- Appendix D, after eqn D12
-- for NGC 1560, there is no 'd'
-- for NGC 3198, there is a d = 9.2
-- for NGC 3115, there is a d = 10
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

-- [[ FAIL - inflection is too far to the right (though not as bad as D.12.a) ... until I changed something, and now it's 100% wrong.
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

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_gravitational_potential_eqn_C7.svg",
	style = 'data lines',
	xlabel = "r (kpc)",
	format = {y = '%.2e'},
	title = 'Gravitational potential of NGC 1560',
	data = {rvec, rvec:map(normphi_for_r_z_eq_0_eqn_5_2_a)},
	{using='1:2', title=''},
}
-- looks wrong.


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
TODOOOO THIS IS NEVER GIVEN!!!!!!
HOW CAN I REPRODUCE THIS GRAPH?!
something about Abel f & g function solution...
--]]
local function v_for_r(r)
	
end

local rvec = xvec * rmax
--[=[
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_3198_normalized_rotation_curve.svg",
	xlabel = "r (kpc)",
	ylabel = "v / c",
	style = 'data lines',
	title = 'Normalized rotation curve of NGC 3198',
	xrange = {0, rmax},
	data = {rvec, rvec:map(v_for_r)},
	{using='1:2', title=''},
}
--]=]



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
