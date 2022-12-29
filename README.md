I'm trying to reproduce the graphs in 2021 Ludwig. 

Seems the paper isn't very clear where it gets its graphs from (no explicit description of the rotation curve graph whatsoever) or on what equation the parameters that it lists are used in.

Errors I think I've found:
* NGC 3115 y0 = 4.0 instead of 0.4 makes the graph look correct.  Especially since NGC 3198 has a value of 8.0, in the same range.
* Eqn 9.1b says `r / r0` but should say `r / reff` ... and maybe have a finally piecewise segment for `se`.

### NGC 1560 ###

For NGC 1560 there are two sets of functions/parameters, one in section 7 and one in appendix D.  seems the parameters and figures in one section are sometimes used in the other section, and there's no clear indication of which is used where.

CHECK:	Fig 1a: NGC 1560 - rotation velocity of a spheroid

CHECK:	Fig 1a: NGC 1560 - sampled rotation velocity of spheroid from 1992 Broeils 

![](Fig__1a_NGC_1560_normalized_rotation_curve_eqn_5.3.svg)

CHECK:	Fig 1b: NGC 1560 mass density

CHECK:	Fig 1b: NGC 1560 gravitational potential

![](Fig__1b_NGC_1560_mass_density_and_potential_eqn_5.2.svg)

CHECK:	Fig 2a: NGC 1560 circular velocity of a spheroid vs fig1a velocity of spheroid

CHECK:	Fig 2a: NGC 1560 - sampled rotation velocity of spheroid from 1992 Broeils 

![](Fig__2a_NGC_1560_rotation_velocity_of_spheriod_eqn_4.14.svg)

CHECK:	Fig 2b: NGC 1560 luminosity profile

CHECK:	Fig 2b: NGC 1560 - sampled luminosity from 1992 Broeils table 3

![](Fig__2b_NGC_1560_luminosity_eqn_D11_using_section_7_numbers.svg)

CHECK:	Fig 3a: NGC 1560 luminosity profile, adjusted, using the Appendix D variables .. this graph's equation is a mess of pieces from the paper.

CHECK:	Fig 3a: NGC 1560 - sampled luminosity from 1992 Broeils table 3

![](Fig__3a_NGC_1560_luminosity_eqn_D11_using_appendix_D_numbers.svg)

CHECK:	Fig 3b: NGC 1560 Sersic index ... I bet my code or the paper has a typo somewhere ... 

![](Fig__3b_NGC_1560_Sersic_index_eqn_D13.svg)

CHECK:	Fig 4a: NGC 1560 normalized mass density (using eqn D 13)

CHECK:	Fig 4a: NGC 1560 - sampled luminosity from 1992 Broeils

![](Fig__4a_derived_NGC_1560_density.svg)

![](Fig__4a_NGC_1560_normalized_mass_density_eqn_5.2b.svg)

TODO:	Fig 4b: NGC 1560 normalized rotation curve (well, it's not the eqn 5.3 ... so where do they get it from?)

This should be reproducable once the norm mass density and grav pot are done -- then just substitute both into Abel's equation (might not be the same as those provided in eqn 5.2) and ... viola?

![](Fig__4b_NGC_1560_normalized_rotation_curve.svg)

CHECK:	Fig 5a: NGC 1560 galactic width.

![](Fig__5a_NGC_1560_galactic_width_eqn_C9.svg)

FIXME:	Fig 5b: NGC 1560 gravitational potential 

TODO:	Fig 5b: NGC 1560 d/dr of gravitational potential

![](Fig__5b_NGC_1560_gravitational_potential.svg)

### NGC 3198 ###

Luckily for this galaxy the author didn't divide his information between two separate sections.

CHECK:	Fig 6a:	NGC 3198 luminosity profile function.

CHECK:	Fig 6a: NGC 3198 - sampled luminosity from 1987 Kent ... I have *a* 1987 Kent paper, but it doesn't have the same data, maybe it's not *the* Kent paper? 

![](Fig__6a_NGC_3198_luminosity_eqn_D11.svg)

CHECK:	Fig 6b: NGC 3198 Sersic index. Using the variables he provides does produce a discontinuity at the boundary of the piecewise segments of the function, and recalculating the coefficients using the equations provided fixed this.

![](Fig__6b_NGC_3198_Sersic_index_eqn_D13.svg)

CHECK: 	Fig 7a: NGC 3198 normalized mass density corrected for high mass-to-light ratio population.

TODO:	Fig 7a:	NGC 3198 normalized mass density matching measured values.

CHECK:	Fig 7a: NGC 3198 - sampled normalized mass density from 1987 Kent .. can't find the values in the paper

![](Fig__7a_derived_NGC_3198_density.svg)

![](Fig__7a_NGC_3198_normalized_mass_density_eqn_5.2b.svg)

FIXME:	Fig 7b: NGC 3198 normalized rotation curve

CHECK:	Fig 7b: NGC 3198 - sampled rotation curve from 1989 Begeman

![](Fig__7b_NGC_3198_normalized_rotation_curve.svg)

I also reconstructed the cited 2006 Cooperstock et al paper's rotation curve, however the fitted function looks like it is a different model than the 2021 Ludwig rotation curve graph.

Probably because I need to solve the Abel equations with NGC 3198's normrho and normph, and normphi is based on the galactic width ...

![](NGC_3198_derived_ram_pressure.svg)

![](NGC_3198_normalized_rotation_curve_2006_Cooperstock_Fig_3.svg)

CHECK:	Fig 8a: NGC 3198 galactic width.

![](Fig__8a_NGC_3198_galactic_width_eqn_C9.svg)

FIXME:	Fig 8b:	NGC 3198 gravitational potential 

FIXME:	Fig 8b:	NGC 3198 gravitational potential derivative.

![](Fig__8b_NGC_3198_gravitational_potential_eqn_C7.svg)

### NGC 3115 ###

CHECK:	Fig 9a:	NGC 3115 luminosity profile.

CHECK:	Fig 9a: NGC 3115 - sampled luminosity from 1987 Capaccioli et al.

![](Fig__9a_NGC_3115_luminosity_eqn_D11.svg)

CHECK:	Fig 9b:	NGC 3115 Sersic index.

![](Fig__9b_NGC_3115_Sersic_index_eqn_D13.svg)

CHECK:	Fig 10a: NGC 3115 normalized mass density corrected for high mass-to-light ratio population 

CHECK:	Fig 10a: NGC 3115 - sampled normalized mass density from 1987 Capaccioli et al.  1987 Capaccioli et al doesn't have these values.  How were they computed for 2021 Ludwig?

TODO	Fig 10a: NGC 3115 normalized mass density matching the samples 

![](Fig_10a_derived_NGC_3115_density.svg)

![](Fig_10a_NGC_3115_normalized_mass_density_eqn_8.4b.svg)

FIXME:	Fig 10b: NGC 3115 normalized rotation curve. 

CHECK:	Fig 10b: NGC 3115 - sampled rotation curve from 1980 Rubin.

![](Fig_10b_NGC_3115_normalized_rotation_curve.svg)

CHECK:	Fig 11a: NGC 3115 galactic width. 

![](Fig_11a_NGC_3115_galactic_width_eqn_C9.svg)

FIXME:	Fig 11b: NGC 3115 gravitational potential

FIXME:	Fig 11b: NGC 3115 gravitational potential derivative

![](Fig_11b_NGC_3115_gravitational_potential_eqn_C7.svg)

FIXME:	Fig 12a: NGC 3115 normalized mass density with modified schwarzschild radius. 

CHECK:	Fig 12a: NGC 3115 - sampled normalized mass density from 1987 Capaccioli et al.

![](Fig_12a_NGC_3115_normalized_mass_density_eqn_8.4a_rs_adj.svg)

FIXME:	Fig 12b: NGC 3115 normalized rotation curve with modified schwarzschild radius. 

CHECK:	Fig 12b: NGC 3115 - sampled rotation curve from 1980 Rubin.

![](Fig_12b_NGC_3115_normalized_rotation_curve_rs_adj.svg)
