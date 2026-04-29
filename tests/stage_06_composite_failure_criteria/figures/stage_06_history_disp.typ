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

#align(center)[#text(size: 11pt, weight: "bold")[Stage 06 Composite Failure Criteria History]]
#v(2mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  rect((1.450, 1.150), (14.950, 8.050), fill: white, stroke: charcoal + 0.65pt)
  line((1.450, 1.150), (14.950, 1.150), stroke: black10 + 0.45pt)
  content((1.270, 1.150), [0], anchor: "east")
  line((1.450, 2.875), (14.950, 2.875), stroke: black10 + 0.45pt)
  content((1.270, 2.875), [672], anchor: "east")
  line((1.450, 4.600), (14.950, 4.600), stroke: black10 + 0.45pt)
  content((1.270, 4.600), [1.34e+03], anchor: "east")
  line((1.450, 6.325), (14.950, 6.325), stroke: black10 + 0.45pt)
  content((1.270, 6.325), [2.02e+03], anchor: "east")
  line((1.450, 8.050), (14.950, 8.050), stroke: black10 + 0.45pt)
  content((1.270, 8.050), [2.69e+03], anchor: "east")
  line((1.450, 7.721), (1.885, 1.951), (2.321, 1.544), (2.756, 1.423), (3.192, 1.370), (3.627, 1.345), (4.063, 1.337), (4.498, 1.343), (4.934, 1.363), (5.369, 1.407), (5.805, 1.506), (6.240, 1.799), (6.676, 5.231), (7.111, 3.144), (7.547, 2.144), (7.982, 1.840), (8.418, 1.707), (8.853, 1.645), (9.289, 1.625), (9.724, 1.638), (10.160, 1.689), (10.595, 1.802), (11.031, 2.050), (11.466, 2.774), (11.902, 1.337), (12.337, 1.339), (12.773, 1.350), (13.208, 1.381), (13.644, 1.446), (14.079, 1.553), (14.515, 1.625), (14.950, 1.553), stroke: garnet + 1.05pt)
  line((1.450, 7.721), (1.885, 1.874), (2.321, 1.525), (2.756, 1.415), (3.192, 1.366), (3.627, 1.344), (4.063, 1.337), (4.498, 1.344), (4.934, 1.366), (5.369, 1.415), (5.805, 1.525), (6.240, 1.874), (6.676, 5.231), (7.111, 2.985), (7.547, 2.100), (7.982, 1.822), (8.418, 1.698), (8.853, 1.642), (9.289, 1.625), (9.724, 1.642), (10.160, 1.698), (10.595, 1.822), (11.031, 2.100), (11.466, 2.985), (11.902, 1.337), (12.337, 1.346), (12.773, 1.367), (13.208, 1.381), (13.644, 1.408), (14.079, 1.506), (14.515, 1.625), (14.950, 1.506), stroke: atlantic + 1.05pt)
  line((1.450, 7.721), (1.885, 1.874), (2.321, 1.525), (2.756, 1.415), (3.192, 1.366), (3.627, 1.344), (4.063, 1.337), (4.498, 1.344), (4.934, 1.366), (5.369, 1.415), (5.805, 1.525), (6.240, 1.874), (6.676, 5.231), (7.111, 2.985), (7.547, 2.100), (7.982, 1.822), (8.418, 1.698), (8.853, 1.642), (9.289, 1.625), (9.724, 1.642), (10.160, 1.698), (10.595, 1.822), (11.031, 2.100), (11.466, 2.985), (11.902, 1.337), (12.337, 1.366), (12.773, 1.417), (13.208, 1.381), (13.644, 1.417), (14.079, 1.612), (14.515, 1.625), (14.950, 1.612), stroke: congaree + 1.05pt)
  line((1.600, 8.650), (2.050, 8.650), stroke: garnet + 1.1pt)
  content((2.150, 8.650), [TSAIWU], anchor: "west")
  line((5.850, 8.650), (6.300, 8.650), stroke: atlantic + 1.1pt)
  content((6.400, 8.650), [HASHIN], anchor: "west")
  line((10.100, 8.650), (10.550, 8.650), stroke: congaree + 1.1pt)
  content((10.650, 8.650), [PUCK], anchor: "west")
  content((1.450, 0.730), [1], anchor: "north")
  content((14.950, 0.730), [32], anchor: "north")
  content((8.200, 0.430), [analytic path id], anchor: "north")
  content((0.400, 4.600), [reported strength (MPa)], angle: 90deg)
})

#text(size: 8pt)[Verdict context: Stage 06 verdict is INCONCLUSIVE due to LAW25-on-TYPE14 starter incompatibility; the figures here show /PROP/TYPE6 proxy runs, not the canonical solid path. Canonical TYPE14 solver samples were not produced; this plot shows the analytic envelope rows recorded in timeseries.csv. Source CSV: #raw("results/timeseries.csv").]
