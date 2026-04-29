#import "@preview/cetz:0.4.1"

#set page(width: 180mm, height: 124mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))

#let garnet = rgb("#73000A")
#let charcoal = rgb("#363636")
#let black70 = rgb("#5C5C5C")
#let black10 = rgb("#ECECEC")
#let horseshoe = rgb("#65780B")
#let atlantic = rgb("#466A9F")

#let rows = csv("../results/timeseries.csv", row-type: dictionary)
#let row = rows.at(0)
#let alpha = float(row.alpha)
#let dy-fem = float(row.dy_fem_over_L)
#let dy-ref = float(row.dy_ref_over_L)
#let gate = 0.02
#let pct(v) = str(calc.round(v * 100, digits: 2)) + "%"

#align(center)[
  #text(size: 11pt, weight: "bold", fill: charcoal)[Stage 02 Cantilever Elastica Check]
]

#v(2mm)

#cetz.canvas(length: 1cm, {
  import cetz.draw: *

  let x0 = 1.35
  let y0 = 0.85
  let w = 13.0
  let h = 6.2
  let xmin = 0.0
  let xmax = 1.1
  let ymin = 0.0
  let ymax = 0.34
  let xm(v) = x0 + (v - xmin) / (xmax - xmin) * w
  let ym(v) = y0 + (v - ymin) / (ymax - ymin) * h
  let x = xm(alpha)
  let y-ref = ym(dy-ref)
  let y-fem = ym(dy-fem)
  let y-low = ym(dy-ref * (1.0 - gate))
  let y-high = ym(dy-ref * (1.0 + gate))
  let s = 0.08

  rect((x0, y0), (x0 + w, y0 + h), fill: white, stroke: charcoal + 0.65pt)
  for tick in (0.0, 0.1, 0.2, 0.3) {
    let y = ym(tick)
    line((x0, y), (x0 + w, y), stroke: black10 + 0.45pt)
    line((x0 - 0.08, y), (x0, y), stroke: charcoal + 0.55pt)
    content((x0 - 0.25, y), [#str(tick)], anchor: "east")
  }
  for tick in (0.0, 0.5, 1.0) {
    let xt = xm(tick)
    line((xt, y0), (xt, y0 - 0.08), stroke: charcoal + 0.55pt)
    content((xt, y0 - 0.35), [#str(tick)], anchor: "north")
  }

  // Two-percent acceptance gate around the Bisshopp-Drucker reference at the sampled alpha.
  line((x - 0.50, y-low), (x + 0.50, y-low), stroke: horseshoe + 0.85pt)
  line((x - 0.50, y-high), (x + 0.50, y-high), stroke: horseshoe + 0.85pt)
  line((x + 0.50, y-low), (x + 0.50, y-high), stroke: horseshoe + 0.85pt)
  content((x + 0.58, y-high), [2% gate], anchor: "west")

  line((x, y-fem), (x, y-ref), stroke: garnet + 1.0pt)
  rect((x - s, y-ref - s), (x + s, y-ref + s), fill: atlantic, stroke: none)
  rect((x - s, y-fem - s), (x + s, y-fem + s), fill: garnet, stroke: none)
  content((x - 0.15, y-ref + 0.34), [reference #str(calc.round(dy-ref, digits: 3))], anchor: "east")
  content((x - 0.15, y-fem + 0.18), [OpenRadioss #str(calc.round(dy-fem, digits: 5))], anchor: "east")
  content((x, y0 - 0.74), [#row.mesh, #row.solver_mode], anchor: "north")

  content((x0 + w / 2, y0 - 1.10), [Dimensionless load alpha], anchor: "north")
  content((x0 - 0.90, y0 + h / 2), [Tip sag dy/L], angle: 90deg)
  content((x0 + 0.10, y0 + h + 0.36), [#row.verdict: #pct(float(row.err_y)) #raw("dy") error], anchor: "west")

  rect((x0 + 8.65, y0 + h + 0.22), (x0 + 8.85, y0 + h + 0.42), fill: atlantic, stroke: none)
  content((x0 + 8.98, y0 + h + 0.32), [Bisshopp-Drucker], anchor: "west")
  rect((x0 + 11.15, y0 + h + 0.22), (x0 + 11.35, y0 + h + 0.42), fill: garnet, stroke: none)
  content((x0 + 11.48, y0 + h + 0.32), [OpenRadioss], anchor: "west")
})

#v(-2mm)

#text(size: 8pt, fill: black70)[Source: #raw("../results/timeseries.csv"); stage verdict INCONCLUSIVE in #raw("../results/results.json").]
