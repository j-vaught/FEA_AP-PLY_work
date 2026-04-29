#import "@preview/cetz:0.3.4"

#set page(width: 180mm, height: 112mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))
#let garnet = rgb("#73000A")
#let charcoal = rgb("#363636")
#let black10 = rgb("#ECECEC")
#let atlantic = rgb("#466A9F")
#let white = rgb("#FFFFFF")

#align(center)[#text(size: 11pt, weight: "bold")[Stage 04 Open-Hole Stress Concentration]]
#v(2mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  rect((1.0, 0.8), (14.8, 7.0), fill: white, stroke: charcoal + 0.65pt)
  line((1.0, 0.800), (14.8, 0.800), stroke: black10 + 0.45pt)
  content((0.84, 0.800), [2.90], anchor: "east")
  line((1.0, 2.350), (14.8, 2.350), stroke: black10 + 0.45pt)
  content((0.84, 2.350), [2.95], anchor: "east")
  line((1.0, 3.900), (14.8, 3.900), stroke: black10 + 0.45pt)
  content((0.84, 3.900), [3.00], anchor: "east")
  line((1.0, 5.450), (14.8, 5.450), stroke: black10 + 0.45pt)
  content((0.84, 5.450), [3.05], anchor: "east")
  line((1.0, 7.000), (14.8, 7.000), stroke: black10 + 0.45pt)
  content((0.84, 7.000), [3.10], anchor: "east")
  line((1.000, 5.315) (7.900, 6.863) (14.800, 7.444), stroke: garnet + 1.15pt)
  circle((1.000, 5.315), radius: 0.06, fill: atlantic, stroke: none)
  circle((7.900, 6.863), radius: 0.06, fill: atlantic, stroke: none)
  circle((14.800, 7.444), radius: 0.06, fill: atlantic, stroke: none)
  line((1.0, 4.985), (14.8, 4.985), stroke: charcoal + 0.8pt)
  content((14.7, 0.18), [Circumferential divisions], anchor: "north-east")
  content((0.22, 3.9), [$K_t$], angle: 90deg)
  content((14.6, 7.35), [Verdict: PASS], anchor: "east")
})

#figure(
  table(
    columns: 4,
    [$N_theta$], [$K_t$], [Howland error], [Ligament L2],
  [32], [3.0456], [0.351\%], [3.322\%],
  [64], [3.0956], [1.996\%], [3.768\%],
  [128], [3.1143], [2.614\%], [3.945\%],
  ),
  caption: [OpenRadioss recovered rim stress concentration compared with Howland's finite-width target.]
)

// raw_points (32, 3.045649) (64, 3.095581) (128, 3.114333)
