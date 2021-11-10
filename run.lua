#!/usr/bin/env lua128

local range = require 'ext.range'
local gnuplot = require 'gnuplot'
local matrix = require 'matrix'
local symmath = require 'symmath'


local Gstorage = {}
setmetatable(_G, {
	__newindex = function(...)
		local t, k, v = ...
		local oldv = Gstorage[k]
		print(tostring(k)..' = '..tostring(v)..(oldv ~= nil and oldv ~= v and (' (from '..oldv..')') or ''))
		Gstorage[k] = v
	end,
	__index = Gstorage,
})


-- Section 5, after eqn 5.3 (sloppiest research paper ever)
Msun = 1.9891e+30	-- kg 
-- wikipedia:
--Msun = 1.98847e+30	-- kg

-- Appendix D, after eqn D11
Lsun = 3.828e+26 	-- W
-- wikipedia:
--Lsun = 3.828e+26	-- W

-- mass-to-light ratio
UpsilonSun = Msun / Lsun	-- kg / W ~ wikipedia: 5133 kg / W

-- me
m_in_pc = 648000 / math.pi * 149597870700


local n = 1000
local xvec = matrix{n}:lambda(function(i) return (i-1)/(n-1) end)	-- cell-centered


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

--[[
text from section 7:
"Using the variable Sérsic index profile, 
the calculated value of the absolute magnitude is Ms = −15.3,
the total luminosity is Ls = 1.02 × 108 Lsun
and the apparent magnitude is md = 12.1.
The total luminosity is calculated up to the maximum galactic radius rmax = 12.2 kpc 
that is estimated by the rotation velocity model described in the next paragraph."
--]]

-- from Appendix D, after eqn D11
-- "As an example, the surface brightness μ B of the dwarf galaxy NGC 1560 analyzed in Sect. 7 can be adjusted to the observed values listed by Broeils [32] by taking..."
-- so I guess that means these values *aren't* supposed to match figures 3 and 4 in Section 7 ...
-- but then why do they match? and why do Section 7's numbers not match?
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
local function mu_eqn_D_11(alpha) 
	return mu0 + (
		alpha < alpha0
		and 5 / (2 * math.log(10)) * (alpha / alpha1) ^ (1 / s1)
		or 5 / (2 * math.log(10)) * (alpha0 / alpha1) ^ (1 / s1) * (1 - s2 / s1 + s2 / s1 * (alpha / alpha0) ^ (1 / s2))
	)
end


-- Figure 3 subtext:
-- "The profiles extend tot he last measurement taken at lp = 5.13 kpc"
lp = 5.13	-- kpc
alpha_lp = alpha_for_r_d(lp, d)
-- alpha_lp = 352.71281868253 ... where the fig 3 xrange comes from

local alphavec = xvec * alpha_lp 

-- looks like figure 3 lhs.  check.  good.
-- .. except that figure 3 is supposed to be Section 7's numbers, and this is Appendix D's numbers ... smh
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_luminosity_eqn_D11_using_appendix_D_numbers.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ (mag arcsec^-2)",
	style = 'data lines',
	title = "Luminosity profile of NGC 1560",
	xrange = {0, alpha_lp},
	yrange = {[3] = 'reverse'},
	data = {alphavec, alphavec:map(mu_eqn_D_11)},
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
alphae = 353.0
se = 0.874
b1 = 0.00245
d2 = -3.22e-6	-- why can't this paper at least provide the d and b coeffs sequentially?

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
--]]

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_luminosity_eqn_D11_using_section_7_numbers.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ (mag arcsec^-2)",
	style = 'data lines',
	title = "Luminosity profile of NGC 1560",
	xrange = {0, alpha_lp},
	yrange = {[3] = 'reverse'},
	data = {alphavec, alphavec:map(mu_eqn_D_11)},
	{using='1:2', title=''}, 
}
-- does eqn D20 match up with se's definition in section 7?
-- eqn D20
se_check = math.log(alphae / alphaeff) / math.log(2 * math.log(10) / 5 * (mu_eqn_D_11(alphae) - mu0))	
-- using the section 7 alpha0: se_check = 0.49270654394305
-- using the Appendix D alpha0: se_check = 0.88439918653437 ... which matches section 7
-- ... WHY DOESN'T THE SECTION 7 ALPHAE MATCH THE SECTION 7 SE?!?!?!?

-- ... annnnd no.  eqn D20's se is half of what section 7's variable listing's se is.  but substutitnig the eqn D20 se makes things worse.
-- however, if we use the alpha0 from Appendix D, then we get se = 0.88439918653437 which is close to correct. 
se = se_check

-- how about eqn D18 vs b1 and d2 in section 7?
b1_check = 1 / (alphae - alpha0) * (2 * (se - s0)
	- range(2,4):mapi(function(i)
		return (i * (alphae - alpha0) + 2 * alpha0) * alpha0 ^ (i - 1) * _G['b'..i]
	end):sum()
	- range(3,4):mapi(function(j)
		return (j - 2) * (alphae - alpha0) ^ j * _G['d'..j]
	end):sum()
)
-- using Section 7 alpha0 we get b1 = 0.00569816109840886 ... double what it should be
-- ... annd once again, b1 is half what the paragraph says it is
b1 = b1_check
-- using the Appendix D alpha0 we get b1 = -0.00202970710142 ... which is the right order, but negative

d2_check = -1 / (alphae * alphae - alpha0 * alpha0) * ( se - s0
	+ (2 - 1) * alpha0 ^ 2 * b2
	+ (3 - 1) * alpha0 ^ 3 * b3
	+ (4 - 1) * alpha0 ^ 4 * b4
	+ (alphae + (3 - 1) * alpha0) * (alphae - alpha0) ^ (3 - 1) * d3
	+ (alphae + (4 - 1) * alpha0) * (alphae - alpha0) ^ (4 - 1) * d4
)
-- d2_check = -3.20679292488114e-06

-- ... and this is at least close 
d2 = d2_check
-- using the Appendix D alpha0 we get d2 = -1.06730771532024e-05  ... which is 3x larger ....


-- -- for the sake of the luminosity profile, section 7 still has no mention of s1, s2, or alpha1 ...
-- -- now the mu0 matches, but if I use this alpha0 then the graphs no longer match.
-- -- and using D8:
-- mu(alpha) = mu0 + (5 / (2 * math.log(10)) * (alpha / alphaeff) ** (1 / s)
-- -- ... but what is 's' ?  or would we use eqn D3 bs ~ 2s - 1/3 <=> s = (bs + 1/3) / 2
-- -- and even with eqn D3 we still aren't given a value of bs for NGC 1560 ... but we do have a way to calculate it ... what a rabbit hole
-- set output "luminosity_NGC_1560_eqn_D8.svg"
-- plot [0:350] mu(x) notitle
-- bleh .. I'm just going to assume Appendix D is used for NGC 1560.  I can only hope Appendix D also has hidden in the text somewhere the equivalent information for NGC 3198 and NGC 3115 ... because I doubt their section info does.

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

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_Sersic_index_eqn_D13.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = 'Sersic index of NGC 1560',
	xrange = {0, alpha_lp},
	data = {alphavec, alphavec:map(s_for_alpha_eqn_D_13)},
	{using='1:2', title=''},
}
-- FAIL. this looks incorrect.


-- how about using D19
local function S_eqn_D_19(alpha) 
	-- eqn D19 ... looks bad
	return math.log(alpha / alphaeff) / math.log(2 * math.log(10) / 5 * (mu_eqn_D_11(alpha) - mu0))	
end
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_Sersic_index_eqn_D19.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = 'Sersic index of NGC 1560',
	xrange = {0, alpha_lp},
	data = {alphavec, alphavec:map(S_eqn_D_19)},
	{using='1:2', title=''},
}



M = 1.52e+10 * Msun -- kg
lambda = 0.134
rho0 = 3.31e-20	-- kg/m^3
l = 3
rs = 1.46e-6	-- kpc
a = 7.19		-- kpc
b = 0.567		-- kpc
rmax = 12.2		-- kpc
alphamax = alpha_for_r_d(rmax, d)
-- alphamax = 838.81021207153
-- hmm, not the NGC 1560 graph alphamax for sure ... which is around 350
-- instead it turns out the alpha (arcsec) associated with the distance of lp = 5.13 kpc , 
-- which is the furthest distance specified in the figure 3 subtext.



-- eqn C5: lambda = 4 pi R0^3 rho0 / (3 M)
R0 = (3 * M * lambda / (4 * math.pi * rho0)) ^ (1/3) / m_in_pc / 1000	-- kpc
-- R0 = 3.08008569838469e+19 m
-- R0 = 9.98187794103888e-01 kpc

-- eqn C4
local function normrho_eqn_C_4(r,z) 
	return 1 / lambda * (
		(b*b * (
			a * r*r
			+ (a + 3 * math.sqrt(b*b + z*z))
			* (a + math.sqrt(b*b + z*z))
			* (a + math.sqrt(b*b + z*z))
		))
		/ (
			3
			* (r*r + (a + math.sqrt(b*b + z*z)) * (a + math.sqrt(b*b + z*z))) ^ (5/2)
			* (b*b + z*z) ^ (3/2)
		)
	)
end
local function normrho_z_eq_0_eqn_C_4(r)
	return normrho_eqn_C_4(r,0)
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
local function normrho_z_eq_0_eqn_D_12(r) 
	return 10 ^ ((-2/5) * (mu_eqn_D_11(r) - mu0)) * (
		r <= r0
		and math.exp(-(r / r1) ^ (1 / s1))
		or math.exp(-(r0 / r1) ^ (1 / s1)) * (1 - s2 / s1 + s2 / s1 * (r / r0) ^ (1 / s2))
	)
end

local rvec = makePow10Range(.01, rmax)
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_1560_normalized_mass_density_eqn_D12.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 1560',
	log = 'x',
	data = {rvec, rvec:map(normrho_z_eq_0_eqn_D_12)},
	{using='1:2', title=''},
}
-- FAIL. this looks incorrect.



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
alphae = 316.8
se = 1.49
b1 = 0.0499
d2 = 0.0000181

rmax = 30.7	-- kpc
Ms = -21
Ls = 1.93e+10 * Lsun
d = 9.2 * 1000 -- kpc
md = 8.85 -- apparent magnitude
rspiral = 4.0	-- kpc
kspiral = 0.1
-- theta = math.atan(ksprial) = math.rad(5.71)
gammai = 0.96 -- kpc
gamma0 = gammai	-- is this a ... function or a constant? is ri and yi functions or constants?
y0 = 8.0
nu = 1.4

-- eqn 8.2
local function r_for_i(i)
	-- 'i' is the index, right? the subscript of 'i' right? and not the imaginary unit *and* the denotation of r ...
	return rspiral * math.exp(2 * math.pi * i * kspiral)
end

-- eqn 8.3
local function y_for_i(i)
	-- to the nu?  or is that a reference? what is nu?
	-- why does the eqn say exp(...)^nu and not exp(nu * ...) ?
	return y0 * math.exp(2 * math.pi * i * kspiral * nu)	
end

-- eqn 8.1
local function Y_eqn_8_1(r)
	local sum = 1
	for i=0,4 do
-- I guess this is the constant they provided?
--		local gammai = gamma_for_i(i)
		sum = sum + (y_for_i(i) * gammai / math.pi) / ((r - r_for_i(i))^2 - gammai ^2)
	end
	return sum
end

-- eqn 8.5
-- basically same as S(alpha) from eqn D13, but with (r, r0, re) swapped with (alpha, alpha0, alphae)
local function s_for_r_eqn_8_5(r)
	--[[ since r and alpha are linear, can we just do this? no.  the coefficient of difference will scale with the polynomials ...
	if r <= r0 then
		return sb(r)
	end
	if r <= re then
		return sd(r)
	end
	return se
	--]]
	-- [[ convert to alpha and evaluate accordingly
	local alpha = alpha_for_r_d(r, d)
	return s_for_alpha_eqn_D_13(alpha)
	--]]
end
-- eqn 8.4
local function normrho_z_eq_0_eqn_8_4(r)
	--[[ 
	return Y_eqn_8_1(r) * 10 ^ (-2/5 * (mu_eqn_D_11(r) - mu0))
	--]]
	-- [[
	return Y_eqn_8_1(r) * math.exp(-(r / reff) ^ (1 / s_for_r_eqn_8_5(r)))
	--]]
	-- "where s(r) is the Sérsic profile adjusted to the luminosity:"
end

rs = 1.20e-5		-- kpc
a = 9.10			-- kpc - major radius bulge
b = 2.64			-- kpc - minor radius bulge
l = 3.0				-- range parameter
lambda = 0.00323	-- coefficient
rmax = 30.7			-- kpc - maximum radius
M = 1.25e+11 * Msun	-- kg - total mass of galaxy
rho0 = 6.54e-21		-- kg/m^3 - central density
Upsilon = 6.50 * UpsilonSun	-- total mass-to-light ratio

alphamax = alpha_for_r_d(rmax, d)	
-- alphamax = 688.29669041151
-- but these graphs end at about alpha=315 ...

-- Figure 6 subtext:
-- "The profiles extend to the last measurement taken at lp = 14.1 kpc"
lp = 14.1	-- kpc
alpha_lp = alpha_for_r_d(lp, d)
-- alpha_lp = 316.12323566131	-- tada - there's where the alpha rhs range of the graph comes from.  fuckin paper.
local alphavec = xvec * alpha_lp

-- eqn D8, doesn't require alpha1 (derived from r1) or alpha0 (derived from r0) or s2, like eqn D11 requires
-- but what is "s" ?
-- after eqn D1: "where s is the Sersic index"
-- but this 's' varies with alpha ... correct?
local function mu_for_eqn_D_8(alpha)
	local s = s_for_alpha_eqn_D_13(alpha)
	return mu0 + 5 / (2 * math.log(10)) * (alpha / alphaeff) ^ (1 / s)
end

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_3198_luminosity_eqn_D11.svg",
	xlabel = "α (arcsec)",
	ylabel = "μ (mag arcsec^-2)",
	style = 'data lines',
	title = "Luminosity profile of NGC 3198",
	xrange = {0, alpha_lp},
	yrange = {[3] = 'reverse'},
	data = {alphavec, alphavec:map(mu_for_eqn_D_8)},
	{using='1:2', title=''}, 
}
-- CHECK.  good.
-- though there is a hiccup in the sd->sb transition at alpha0
-- solution? re-derive the parenthesis values like b1

gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_3198_Sersic_index_eqn_D13.svg",
	xlabel = "α (arcsec)",
	ylabel = "s",
	style = 'data lines',
	title = "Sersic index of NGC 3198",
	xrange = {0, alpha_lp},
	data = {alphavec, alphavec:map(s_for_alpha_eqn_D_13)},
	{using='1:2', title=''}, 
}
-- CHECK.  good.
-- same as prev graph: fix the hiccup.

-- "The mass density profile extends to the last measurement taken at lp = 1.41 kpc"
-- so looks like this mass density graph uses lp as the rmax instead of rmax (like I think figure 3 used)
local rvec = makePow10Range(.1, lp)
gnuplot{
	terminal = 'svg size 1024,768 background rgb "white"',
	output = "NGC_3198_normalized_mass_density_eqn_C4.svg",
	xlabel = "r (kpc)",
	ylabel = "ρ / ρ0",
	style = 'data lines',
	title = 'Normalized mass density of NGC 3198',
	log = 'xy',
	xrange = {.1, lp},
	--yrange = {0, 1},
	data = {rvec, rvec:map(normrho_z_eq_0_eqn_8_4)},
	{using='1:2', title=''},
}
