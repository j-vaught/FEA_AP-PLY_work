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

#align(center)[#text(size: 11pt, weight: "bold")[Stage 03 Isotropic Dogbone E8 History]]
#v(2mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  rect((1.450, 1.150), (14.950, 8.050), fill: white, stroke: charcoal + 0.65pt)
  line((1.450, 1.150), (14.950, 1.150), stroke: black10 + 0.45pt)
  content((1.270, 1.150), [0], anchor: "east")
  line((1.450, 2.875), (14.950, 2.875), stroke: black10 + 0.45pt)
  content((1.270, 2.875), [79.09], anchor: "east")
  line((1.450, 4.600), (14.950, 4.600), stroke: black10 + 0.45pt)
  content((1.270, 4.600), [158.2], anchor: "east")
  line((1.450, 6.325), (14.950, 6.325), stroke: black10 + 0.45pt)
  content((1.270, 6.325), [237.3], anchor: "east")
  line((1.450, 8.050), (14.950, 8.050), stroke: black10 + 0.45pt)
  content((1.270, 8.050), [316.4], anchor: "east")
  line((1.450, 3.875), (1.937, 5.511), (2.267, 6.604), (14.950, 7.721), stroke: garnet + 1.05pt)
  line((1.450, 3.876), (1.937, 5.512), (2.267, 6.603), (14.950, 7.693), stroke: atlantic + 1.05pt)
  line((1.600, 8.650), (2.050, 8.650), stroke: garnet + 1.1pt)
  content((2.150, 8.650), [OpenRadioss mean stress], anchor: "west")
  line((5.850, 8.650), (6.300, 8.650), stroke: atlantic + 1.1pt)
  content((6.400, 8.650), [closed-form load/area], anchor: "west")
  content((1.450, 0.730), [0.06247], anchor: "north")
  content((14.950, 0.730), [1.103], anchor: "north")
  content((8.200, 0.430), [gauge strain (%)], anchor: "north")
  content((0.400, 4.600), [gauge stress (MPa)], angle: 90deg)
})

#text(size: 8pt)[Verdict context: Stage 03 verdict is PASS; these figures use the medium mesh that produced the reported metrics. Stress-strain points are extracted from final VTK frames for each load factor. Source CSV: #raw("results/timeseries.csv").]
