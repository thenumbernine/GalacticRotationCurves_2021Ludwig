#!/usr/bin/env lua128

-- integration format: integrate(f, xL, xR, n)
local integrateRectangular = require 'integraterectangular'
local integrateTrapezoid = require 'integratetrapezoid'
local integrateSimpson = require 'integratesimpson'
local integrateSimpson3_8 = require 'integratesimpson2'
local matrix = require 'matrix'
local gnuplot = require 'gnuplot'

local function f(x) return math.sqrt(1 - x*x) end
local xL = 0
local xR = 1
local rows = matrix()
for n=1,1000 do
	table.insert(rows, matrix{
		n,
--[[ TODO my gnuplot lib has no good way to print with higher-than-default precision
		('%.20f'):format(math.pi-4*integrateRectangular(f, xL, xR, n)),
		('%.20f'):format(math.pi-4*integrateTrapezoid(f, xL, xR, n)),
		('%.20f'):format(math.pi-4*integrateSimpson(f, xL, xR, n)),
		('%.20f'):format(math.pi-4*integrateSimpson3_8(f, xL, xR, n)),
--]]
-- [[
		(math.pi-4*integrateRectangular(f, xL, xR, n)),
		(math.pi-4*integrateTrapezoid(f, xL, xR, n)),
		(math.pi-4*integrateSimpson(f, xL, xR, n)),
		(math.pi-4*integrateSimpson3_8(f, xL, xR, n)),
--]]
	})
end

gnuplot{
	persist = true,
	style = 'data lines',
	data = rows:T(),
	log = 'xy',
	{using='1:2', title='rectangular'},
	{using='1:3', title='trapezoid'},
	{using='1:4', title='simpson'},
	{using='1:5', title='simpson2'},
}
