I'm trying to reproduce the graphs in 2021 Ludwig. 

Seems the paper isn't very clear where it gets its graphs from (no explicit description of the rotation curve graph whatsoever) or on what equation the parameters that it lists are used in.

### NGC 1560 ###

For NGC 1560 there are two sets of functions/parameters, one in section 7 and one in appendix D.  seems the parameters and figures in one section are sometimes used in the other section, and there's no clear indication of which is used where.

I was able to reproduce the Section 7 Luminosity profile picture using the Appendix D variables.  

The Section 7 description of the Luminosity profile doesn't seem to match the Section 7 graph.

I haven't been able to reproduce the Sersic index. 

I haven't been able to reproduce the normalized mass density. 

And nowhere in the whole fucking paper is there a formula for the rotation curve, so I don't know where he came up with that graph. 

Haven't tried at the galactic width.  Can't find the function in the paper, only hints at numerically solving a PDE description of the width.

Haven't tried at the gravitational potential.


### NGC 3198 ###

Luckily for this galaxy the author didn't divide his information between two separate sections.

I was able to reproduce the luminosity profile function.

I was able to reproduce the Sersic index, though using the variables he provides does produce a discontinuity at the boundary of the piecewise segments of the function, and recalculating the coefficients using the equations provided fixed this.

I was able to reproduce the normalized mass density.

Once again I couldn't find the description of the rotation curve.  I did go back to the cited 2006 Cooperstock et al paper which includes samples as well as its own rotation curve function and was able to reproduce that papers' graph, however the fitted function looks like it is a different model than the 2021 Ludwig rotation curve graph.

Haven't found a galactic width function in the paper.  I didn't see any galactic width function in the 2006 Cooperstock paper as well.

Haven't tried the gravitational potential.
