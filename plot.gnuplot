#!/usr/bin/env gnuplot
set terminal svg size 1024,768 background rgb "white"


# Section 5, after eqn 5.3 (sloppiest research paper ever)
Msun = 1.9891e+30	# kg

# Appendix D, after eqn D11
Lsun = 3.828e+26 	# W

set output "NGC_1560_luminosity_eqn_D11.svg"
set xlabel "alpha (arcsec)"
set ylabel "mu (mag arcsec^-2)"
set title "Luminosity profile of NGC 1560"
set yrange reverse

# from section 7:
# Using the variable Sérsic index profile, 
# the calculated value of the absolute magnitude is Ms = −15.3,
# the total luminosity is Ls = 1.02 × 108 Lsun
# and the apparent magnitude is md = 12.1.
# The total luminosity is calculated up to the maximum galactic radius rmax = 12.2 kpc 
# that is estimated by the rotation velocity model described in the next paragraph.

# from Appendix D, after eqn D11
mu0 = 22.28
alpha0 = 61.46	# arcsec
s1 = 0.435
alpha1 = 99.05	# arcsec
s2 = 1.144
# r = d * (pi / (180. * 3600.)) * alpha
Ms = -16.5
Ls = 3.06e+8 * Lsun	# W
d = 3.0 # Mpc
md = 10.9
# after eqn D12:
r0 = d * (pi / (180. * 3600.)) * alpha0
r1 = d * (pi / (180. * 3600.)) * alpha1

# eqn D11:
mu(alpha) = mu0 + ( \
		alpha < alpha0 \
		? 5. / (2. * log(10.)) * (alpha / alpha1) ** (1. / s1) \
		: 5. / (2. * log(10.)) * (alpha0 / alpha1) ** (1. / s1) * (1. - s2 / s1 + s2 / s1 * (alpha / alpha0) ** (1. / s2)) \
	)

plot [0:350] mu(x) notitle
set yrange noreverse
# looks like figure 3 lhs.  good.


# # section 7 1st paragraph:
mu0 = 22.27
alphaeff = 116.5
reff = 1.69			# kpc
#alpha0 = 140.3		# but didn't we just define an alpha0? and doesn't it work better?
r0 = 2.04 			# kpc
s0 = 0.360
b2 = 0.0000344
b3 = -1.41e-7
b4 = -4.05e-10
d3 = -1.11e-9
d4 = 5.73e-11
alphae = 353.0
se = 0.874
b1 = 0.00245
d2 = -3.22e-6	# why can't this paper at least provide the d and b coeffs sequentially?

# does eqn D20 match up with se's definition in section 7?
se_check = log(alphae / alphaeff) / log(2. * log(10.) / 5. * (mu(alphae) - mu0))	# eqn D20
print(sprintf("se = %.14f", se_check))	# se = 0.49270654394305
# ... annnnd no.  eqn D20's se is half of what section 7's variable listing's se is.  but substutitnig the eqn D20 se makes things worse.
# however, if we use the alpha0 from Appendix D, then we get se = 0.88439918653437 which is close to correct. 
se = se_check

# how about eqn D18 vs b1 and d2 in section 7?
b1_check = 1. / (alphae - alpha0) * (2. * (se - s0) \
	- (2. * (alphae - alpha0) + 2. * alpha0) * alpha0 ** (2. - 1.) * b2 \
	- (3. * (alphae - alpha0) + 2. * alpha0) * alpha0 ** (3. - 1.) * b3 \
	- (4. * (alphae - alpha0) + 2. * alpha0) * alpha0 ** (4. - 1.) * b4 \
	- (3. - 2.) * (alphae - alpha0) ** 3. * d3 \
	- (4. - 2.) * (alphae - alpha0) ** 4. * d4 \
)
print(sprintf("b1 = %.14f", b1_check))	# b1 = 0.00569816109840886	
# ... annd once again, b1 is half what the paragraph says it is
b1 = b1_check
# using the Appendix D alpha0 we get b1 = -0.00202970710142 ... which is the right order, but negative

d2_check = -1. / (alphae * alphae - alpha0 * alpha0) * ( se - s0 \
	+ (2. - 1.) * alpha0 ** 2. * b2 \
	+ (3. - 1.) * alpha0 ** 3. * b3 \
	+ (4. - 1.) * alpha0 ** 4. * b4 \
	+ (alphae + (3. - 1.) * alpha0) * (alphae - alpha0) ** (3. - 1.) * d3 \
	+ (alphae + (4. - 1.) * alpha0) * (alphae - alpha0) ** (4. - 1.) * d4 \
)
print(sprintf("d2 = %.14e", d2_check))	# d2 = -3.20679292488114e-06
# ... and this is at least close 
d2 = d2_check
# using the Appendix D alpha0 we get d2 = -1.06730771532024e-05  ... which is 3x larger ....


# # for the sake of the luminosity profile, section 7 still has no mention of s1, s2, or alpha1 ...
# # now the mu0 matches, but if I use this alpha0 then the graphs no longer match.
# # and using D8:
# mu(alpha) = mu0 + (5. / (2. * log(10.)) * (alpha / alphaeff) ** (1. / s)
# # ... but what is 's' ?  or would we use eqn D3 bs ~ 2s - 1/3 <=> s = (bs + 1/3) / 2
# # and even with eqn D3 we still aren't given a value of bs for NGC 1560 ... but we do have a way to calculate it ... what a rabbit hole
# set output "luminosity_NGC_1560_eqn_D8.svg"
# plot [0:350] mu(x) notitle
# bleh .. I'm just going to assume Appendix D is used for NGC 1560.  I can only hope Appendix D also has hidden in the text somewhere the equivalent information for NGC 3198 and NGC 3115 ... because I doubt their section info does.

set output "NGC_1560_Sersic_index_eqn_D13.svg"
set xlabel "S"
set ylabel "alpha (arcsec)"

# Sb(alpha) = b0 + alpha * (b1 + alpha * (b2 + alpha * b3))	# eqn D14 def of Sb ... requires b0 which doesn't exist
# Sd(alpha) = d0 + alpha * (d1 + alpha * (d2 + alpha * d3))	# eqn D14 def of Sd ... requires d0, d1 which doesn't exist
Sb(alpha) = s0 + alpha * (b1 + alpha * (b2 + alpha * b3))	# eqn D17 def of Sb
Sd(alpha) = (dalpha = alphae - alpha, se + dalpha * dalpha * (d2 + dalpha * (d3 + dalpha * d4))) # eqn D17 def of Sd
S(alpha) = alpha <= alpha0 ? Sb(alpha) : (alpha <= alphae ? Sd(alpha) : se)	# eqn D13 ... looks bad

#S(alpha) = log(alpha / alphaeff) / log(2. * log(10.) / 5. * (mu(alpha) - mu0))	# eqn D19 ... looks bad

plot [0:350] S(x) notitle
# this looks incorrect.


set output "NGC_1560_mass_density.svg"
set xlabel "r (kpc)"
set ylabel "rho / rho0"
set title "Normalized mass density of NGC 1560"
set log x

M = 1.52e+10 * Msun # kg
lambda = 0.134
rho0 = 3.31e-20	# kg/m^3
_l = 3.
rs = 1.46e-6	# kpc
a = 7.19		# kpc
b = 0.567		# kpc
rmax = 12.2		# kpc

m_in_pc = 648000. / pi * 149597870700.

# eqn C5: lambda = 4 pi R0^3 rho0 / (3 M)
R0 = (3. * M * lambda / (4. * pi * rho0)) ** (1./3.) / m_in_pc / 1000.	# kpc
print(sprintf("R0 = %.14e", R0))
# R0 = 3.08008569838469e+19 m
# R0 = 9.98187794103888e-01 kpc

# eqn C4
# normrho(r,z) = 1. / lambda * ( \
# 	(b*b * ( \
# 		a * r*r \
# 		+ (a + 3. * sqrt(b*b + z*z)) \
# 		* (a + sqrt(b*b + z*z)) \
# 		* (a + sqrt(b*b + z*z)) \
# 	)) \
# 	/ ( \
# 		3. \
# 		* (r*r + (a + sqrt(b*b + z*z)) * (a + sqrt(b*b + z*z))) ** (5./2.) \
# 		* (b*b + z*z) ** (3./2.) \
# 	) \
# )
# normrho_zeq0(r) = normrho(r,0)

# Appendix D, after eqn D12
# for NGC 1560, there is no 'd'
# for NGC 3198, there is a d = 9.2
# for NGC 3115, there is a d = 10
# ... so what do we use for NGC 1560?
r0_check = d * (pi / (180. * 3600.)) * alpha0
print(sprintf("r1 = %.14f", r0_check))

r1 = d * (pi / (180. * 3600.)) * alpha1
# eqn D12
normrho_zeq0(r) = 10. ** ((-2./5.) * (mu(r) - mu0)) * ( \
	r <= r0 \
	? exp(-(r / r1) ** (1. / s1)) \
	: exp(-(r / r0) ** (1. / s1)) * (1. - s2 / s1 + s2 / s1 * (r / r0) ** (1. / s2)) \
)

plot [0.01:20] normrho_zeq0(x)
# this looks incorrect.
unset log x



# Alright, after enough nonsense and contradictions, and being split between two sections (Section 7 and Appendix D) that have contradicting information, 
# lets try starting fresh with Section 8, NGC 3198 ...
#
# "The analysis of the spiral galaxy NGC 3198 [11,34,35] follows the same procedure established for NGC 1560"
# Already we're off to a bad start.

mu0 = 19.49
alphaeff = 22.4
reff = 1.00	# kpc
alpha0 = 154.0
r0 = 6.87	# kpc
s0 = 0.586
b2 = -0.000433
b3 = 5.86e-8
b4 = 3.29e-9
b5 = -1.06e-11
b6 = 1.52e13
b7 = 2.90e-15
b8 = -1.75e-17
d3 = -4.96e-7
d4 = 1.85e-9
d5 = 1.07e-11
d6 = 2.04e-14
d7 = -1.75e-16
d8 = -1.20e-18
alphae = 316.8
se = 1.49
b1 = 0.0499
d2 = 0.0000181
rmax = 30.7	# kpc
Ms = -21
Ls = 1.93e+10 * Lsun
d = 9.2 # Mpc
md = 8.85 # apparent magnitude
