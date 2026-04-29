#import "@preview/cetz:0.3.4"

#set page(width: 180mm, height: 118mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))
#let garnet = rgb("#73000A")
#let charcoal = rgb("#363636")
#let black10 = rgb("#ECECEC")
#let atlantic = rgb("#466A9F")
#let white = rgb("#FFFFFF")

#align(center)[#text(size: 11pt, weight: "bold")[Stage 03 A36 Dogbone Stress-Strain]]
#v(2mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  rect((1.0, 0.8), (14.6, 7.0), fill: white, stroke: charcoal + 0.65pt)
  line((1.0, 0.800), (14.6, 0.800), stroke: black10 + 0.45pt)
  content((0.84, 0.800), [0], anchor: "east")
  line((1.0, 2.350), (14.6, 2.350), stroke: black10 + 0.45pt)
  content((0.84, 2.350), [100], anchor: "east")
  line((1.0, 3.900), (14.6, 3.900), stroke: black10 + 0.45pt)
  content((0.84, 3.900), [200], anchor: "east")
  line((1.0, 5.450), (14.6, 5.450), stroke: black10 + 0.45pt)
  content((0.84, 5.450), [300], anchor: "east")
  line((1.0, 7.000), (14.6, 7.000), stroke: black10 + 0.45pt)
  content((0.84, 7.000), [400], anchor: "east")
  line((1.700, 2.737) (2.121, 3.899) (2.406, 4.676) (13.364, 5.470), stroke: garnet + 1.15pt)
  circle((1.700, 2.737), radius: 0.06, fill: atlantic, stroke: none)
  circle((2.121, 3.899), radius: 0.06, fill: atlantic, stroke: none)
  circle((2.406, 4.676), radius: 0.06, fill: atlantic, stroke: none)
  circle((13.364, 5.470), radius: 0.06, fill: atlantic, stroke: none)
  content((7.8, 0.18), [Gauge strain], anchor: "north")
  content((0.22, 3.9), [Gauge stress (MPa)], angle: 90deg)
  content((14.5, 7.35), [Verdict: PASS], anchor: "east")
})

#figure(
  table(
    columns: 4,
    [Load factor], [Gauge strain], [Stress (MPa)], [Stress error],
  [0.50], [0.000625], [124.94], [0.044\%],
  [0.80], [0.001000], [199.96], [0.021\%],
  [1.00], [0.001254], [250.04], [0.014\%],
  [1.20], [0.011027], [301.29], [0.430\%],
  ),
  caption: [Stage 03 extracted gauge response from OpenRadioss final frames.]
)
