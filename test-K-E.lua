#!/usr/bin/env lua128

local gnuplot = require 'gnuplot'
local matrix = require 'matrix'

local K = require 'K'
local E = require 'E'

local n = 1000
local u = matrix{n}:lambda(function(i) return (i-1)/(n-1) end)
local x = u * 2 - 1

-- verify K works
gnuplot{
	persist = true,
	style = 'data lines',
	data = {x, x:map(K)},
	{using='1:2', title='K'},
}

-- verify E works
gnuplot{
	persist = true,
	style = 'data lines',
	data = {x, x:map(E)},
	{using='1:2', title='E'},
}
