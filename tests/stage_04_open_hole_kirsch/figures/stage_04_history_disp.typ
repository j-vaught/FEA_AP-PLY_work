#import "@preview/cetz:0.3.4"

#set page(width: 185mm, height: 118mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 8.8pt, fill: rgb("#363636"))
#let garnet = rgb("#73000A")
#let atlantic = rgb("#466A9F")
#let horseshoe = rgb("#65780B")
#let congaree = rgb("#1F414D")
#let charcoal = rgb("#363636")
#let black10 = rgb("#ECECEC")
#let white = rgb("#FFFFFF")

#align(center)[#text(size: 11pt, weight: "bold")[Stage 04 Open-Hole Kirsch History]]
#v(2mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  rect((1.450, 1.150), (14.950, 8.050), fill: white, stroke: charcoal + 0.65pt)
  line((1.450, 1.150), (14.950, 1.150), stroke: black10 + 0.45pt)
  content((1.270, 1.150), [0], anchor: "east")
  line((1.450, 2.875), (14.950, 2.875), stroke: black10 + 0.45pt)
  content((1.270, 2.875), [0.8175], anchor: "east")
  line((1.450, 4.600), (14.950, 4.600), stroke: black10 + 0.45pt)
  content((1.270, 4.600), [1.635], anchor: "east")
  line((1.450, 6.325), (14.950, 6.325), stroke: black10 + 0.45pt)
  content((1.270, 6.325), [2.453], anchor: "east")
  line((1.450, 8.050), (14.950, 8.050), stroke: black10 + 0.45pt)
  content((1.270, 8.050), [3.27], anchor: "east")
  line((1.450, 7.577), (5.950, 7.682), (14.950, 7.721), stroke: garnet + 1.05pt)
  line((1.450, 7.554), (5.950, 7.554), (14.950, 7.554), stroke: atlantic + 1.05pt)
  line((1.450, 7.480), (5.950, 7.480), (14.950, 7.480), stroke: horseshoe + 1.05pt)
  line((1.600, 8.650), (2.050, 8.650), stroke: garnet + 1.1pt)
  content((2.150, 8.650), [OpenRadioss recovered Kt], anchor: "west")
  line((5.850, 8.650), (6.300, 8.650), stroke: atlantic + 1.1pt)
  content((6.400, 8.650), [Howland target], anchor: "west")
  line((10.100, 8.650), (10.550, 8.650), stroke: horseshoe + 1.1pt)
  content((10.650, 8.650), [Kirsch infinite plate], anchor: "west")
  content((1.450, 0.730), [32], anchor: "north")
  content((14.950, 0.730), [128], anchor: "north")
  content((8.200, 0.430), [circumferential divisions], anchor: "north")
  content((0.400, 4.600), [stress concentration Kt], angle: 90deg)
})

#text(size: 8pt)[Verdict context: Stage 04 verdict is PASS; the primary gate is the Ntheta 64 Howland Kt comparison. The primary gate is Ntheta 64 against the Howland finite-width target. Source CSV: #raw("results/timeseries.csv").]
