#import "@preview/cetz:0.3.4"

#set page(width: 180mm, height: 118mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))
#let garnet = rgb("#73000A")
#let charcoal = rgb("#363636")
#let black70 = rgb("#5C5C5C")
#let black10 = rgb("#ECECEC")
#let atlantic = rgb("#466A9F")
#let horseshoe = rgb("#65780B")
#let white = rgb("#FFFFFF")

#align(center)[#text(size: 11pt, weight: "bold")[Stage 02 Cantilever Elastica Check]]
#v(2mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  rect((1.150, 0.850), (14.550, 6.850), fill: white, stroke: charcoal + 0.65pt)
  line((1.150, 0.850), (14.550, 0.850), stroke: black10 + 0.45pt)
  content((0.990, 0.850), [0.0], anchor: "east")
  line((1.150, 4.466), (14.550, 4.466), stroke: black10 + 0.45pt)
  content((0.990, 4.466), [0.2], anchor: "east")
  circle((14.550, 6.305), radius: 0.055, fill: atlantic, stroke: none)
  circle((14.550, 0.851), radius: 0.055, fill: garnet, stroke: none)
  content((14.550, 0.530), [1], anchor: "north")
  content((7.850, 0.070), [Dimensionless load alpha], anchor: "north")
  content((0.270, 3.850), [Tip sag, dy/L], angle: 90deg)
  line((10.050, 7.250), (10.800, 7.250), stroke: atlantic + 0.95pt)
  content((10.900, 7.250), [Elastica], anchor: "west")
  line((12.350, 7.250), (13.100, 7.250), stroke: garnet + 1.15pt)
  content((13.200, 7.250), [OpenRadioss], anchor: "west")
})

#text(size: 8pt, fill: black70)[Source CSV: #raw("../results/timeseries.csv")]

#figure(
  table(
    columns: 5,
    [alpha], [mesh], [dy error], [dx error], [verdict],
  [1], [coarse], [99.99\%], [100.00\%], [FAIL],
  ),
  caption: [Stage 02 tip displacement error by load and mesh.]
)
