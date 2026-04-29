#import "@preview/cetz:0.3.4"

#set page(width: 180mm, height: 118mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))

#let garnet = rgb("#73000A")
#let charcoal = rgb("#363636")
#let black70 = rgb("#5C5C5C")
#let black50 = rgb("#A2A2A2")
#let black10 = rgb("#ECECEC")
#let horseshoe = rgb("#65780B")
#let atlantic = rgb("#466A9F")
#let white = rgb("#FFFFFF")

#align(center)[#text(size: 11pt, weight: "bold")[Stage 01 Beam Bending Error Gate]]
#v(2mm)

#cetz.canvas(length: 1cm, {
  import cetz.draw: *

  rect((1.250, 0.850), (14.500, 6.650), fill: white, stroke: charcoal + 0.65pt)
  line((1.250, 0.850), (14.500, 0.850), stroke: black10 + 0.45pt)
  content((1.090, 0.850), [0.0\%], anchor: "east")
  line((1.250, 2.010), (14.500, 2.010), stroke: black10 + 0.45pt)
  content((1.090, 2.010), [0.5\%], anchor: "east")
  line((1.250, 3.170), (14.500, 3.170), stroke: horseshoe + 0.9pt)
  content((1.090, 3.170), [1.0\%], anchor: "east")
  line((1.250, 4.330), (14.500, 4.330), stroke: black10 + 0.45pt)
  content((1.090, 4.330), [1.5\%], anchor: "east")
  line((1.250, 5.490), (14.500, 5.490), stroke: black10 + 0.45pt)
  content((1.090, 5.490), [2.0\%], anchor: "east")
  line((1.250, 6.650), (14.500, 6.650), stroke: black10 + 0.45pt)
  content((1.090, 6.650), [2.5\%], anchor: "east")
  rect((2.566, 0.850), (3.246, 5.236), fill: garnet, stroke: none)
  line((2.566, 5.236), (3.246, 5.236), stroke: charcoal + 0.6pt)
  content((2.906, 5.476), [1.891\%], anchor: "south")
  content((2.906, 0.530), [3pt M0], anchor: "north")
  rect((5.879, 0.850), (6.559, 5.254), fill: garnet, stroke: none)
  line((5.879, 5.254), (6.559, 5.254), stroke: charcoal + 0.6pt)
  content((6.219, 5.494), [1.898\%], anchor: "south")
  content((6.219, 0.530), [4pt M0], anchor: "north")
  rect((9.191, 0.850), (9.871, 3.096), fill: horseshoe, stroke: none)
  line((9.191, 3.096), (9.871, 3.096), stroke: charcoal + 0.6pt)
  content((9.531, 3.336), [0.968\%], anchor: "south")
  content((9.531, 0.530), [3pt M1], anchor: "north")
  rect((12.504, 0.850), (13.184, 3.096), fill: horseshoe, stroke: none)
  line((12.504, 3.096), (13.184, 3.096), stroke: charcoal + 0.6pt)
  content((12.844, 3.336), [0.968\%], anchor: "south")
  content((12.844, 0.530), [4pt M1], anchor: "north")
  content((7.875, 0.070), [Mesh / load case], anchor: "north")
  content((0.300, 3.750), [Euler-Bernoulli relative error], angle: 90deg)
  content((1.430, 3.410), [1\% pass gate], anchor: "west")
  rect((9.750, 6.930), (10.200, 7.230), fill: garnet, stroke: none)
  content((10.320, 7.080), [Coarse mesh], anchor: "west")
  rect((12.300, 6.930), (12.750, 7.230), fill: horseshoe, stroke: none)
  content((12.870, 7.080), [Gated mesh], anchor: "west")
})

#text(size: 8pt, fill: black70)[Source CSV: #raw("../results/timeseries.csv")]

#figure(
  table(
    columns: 4,
    [Load case], [Mesh], [Euler-Bernoulli error], [Verdict],
  [3pt], [M0], [1.891\%], [FAIL],
  [4pt], [M0], [1.898\%], [FAIL],
  [3pt], [M1], [0.968\%], [PASS],
  [4pt], [M1], [0.968\%], [PASS],
  ),
  caption: [Stage 01 midspan displacement error by load case and mesh.]
)
