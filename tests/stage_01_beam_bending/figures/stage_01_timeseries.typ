#import "@preview/cetz:0.4.1"

#set page(width: 180mm, height: 118mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))

#let garnet = rgb("#73000A")
#let charcoal = rgb("#363636")
#let black70 = rgb("#5C5C5C")
#let black50 = rgb("#A2A2A2")
#let black10 = rgb("#ECECEC")
#let horseshoe = rgb("#65780B")
#let atlantic = rgb("#466A9F")

#let rows = csv("../results/timeseries.csv", row-type: dictionary)
#let row = rows.at(0)
#let err = float(row.rel_error_eb)
#let gate = 0.01
#let pct(v) = str(calc.round(v * 100, digits: 2)) + "%"
#let log10(v) = calc.ln(v) / calc.ln(10)

#align(center)[
  #text(size: 11pt, weight: "bold", fill: charcoal)[Stage 01 Beam Bending Error Gate]
]

#v(2mm)

#cetz.canvas(length: 1cm, {
  import cetz.draw: *

  let x0 = 1.35
  let y0 = 0.85
  let w = 13.0
  let h = 6.2
  let ymin = 0.001
  let ymax = 1.0
  let ym(v) = y0 + (log10(v) - log10(ymin)) / (log10(ymax) - log10(ymin)) * h
  let x = x0 + w * 0.52
  let barw = 1.15
  let yerr = ym(err)
  let ygate = ym(gate)

  // Plot frame and grid.
  rect((x0, y0), (x0 + w, y0 + h), fill: white, stroke: charcoal + 0.65pt)
  for tick in (0.001, 0.01, 0.1, 1.0) {
    let y = ym(tick)
    line((x0, y), (x0 + w, y), stroke: if tick == gate { horseshoe + 0.95pt } else { black10 + 0.45pt })
    line((x0 - 0.08, y), (x0, y), stroke: charcoal + 0.55pt)
    content((x0 - 0.32, y), [#pct(tick)], anchor: "east")
  }

  // Euler-Bernoulli closed-form reference is zero error; on a log axis it lies below the plotted floor.
  line((x0, y0 + 0.18), (x0 + w, y0 + 0.18), stroke: atlantic + 0.7pt)
  content((x0 + w - 0.1, y0 + 0.43), [Euler-Bernoulli reference: 0% error], anchor: "east")

  // Gated solver result.
  rect((x - barw / 2, y0), (x + barw / 2, yerr), fill: garnet, stroke: none)
  line((x - barw / 2, yerr), (x + barw / 2, yerr), stroke: charcoal + 0.65pt)
  content((x, yerr + 0.38), [#pct(err)], anchor: "south")
  content((x, y0 - 0.42), [#row.load_case #row.mesh], anchor: "north")

  // Axis labels and gate annotation.
  content((x0 + w / 2, y0 - 0.85), [Mesh / load case], anchor: "north")
  content((x0 - 0.95, y0 + h / 2), [Euler-Bernoulli relative error (%)], angle: 90deg)
  content((x0 + 0.15, ygate + 0.24), [1% pass gate], anchor: "west")

  // Compact legend.
  line((x0 + 8.4, y0 + h + 0.47), (x0 + 9.15, y0 + h + 0.47), stroke: horseshoe + 0.95pt)
  content((x0 + 9.25, y0 + h + 0.47), [Spec gate], anchor: "west")
  rect((x0 + 10.8, y0 + h + 0.31), (x0 + 11.25, y0 + h + 0.63), fill: garnet, stroke: none)
  content((x0 + 11.35, y0 + h + 0.47), [Fail mismatch], anchor: "west")
  line((x0 + 12.75, y0 + h + 0.47), (x0 + 13.5, y0 + h + 0.47), stroke: atlantic + 0.7pt)
  content((x0 + 13.6, y0 + h + 0.47), [Closed form], anchor: "west")
})

#v(1mm)

#text(size: 8pt, fill: black70)[Source: #raw("../results/timeseries.csv"). Verdict from #row.load_case #row.mesh is #row.verdict; mesh contains #row.elements elements and #row.nodes nodes.]
