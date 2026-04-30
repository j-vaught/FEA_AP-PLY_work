#import "@preview/cetz:0.3.4"
#set page(width: 202mm, height: 118mm, margin: 7mm)
#set text(font: "Libertinus Serif", size: 8pt, fill: rgb("#363636"))
#let garnet = rgb("#73000A")
#let horseshoe = rgb("#65780B")
#let honeycomb = rgb("#A49137")
#let atlantic = rgb("#466A9F")
#let congaree = rgb("#1F414D")
#let black10 = rgb("#ECECEC")
#let black50 = rgb("#A2A2A2")
#let charcoal = rgb("#363636")
#let white = rgb("#FFFFFF")

#let stages = (
  ("01", "PASS", "beam"),
  ("02", "INCONCLUSIVE", "alpha=1"),
  ("03", "PASS", "E8"),
  ("04", "PASS", "Kt"),
  ("05", "PASS", "damage"),
  ("06", "PASS", "cards"),
  ("07", "PASS", "UD"),
  ("08", "PASS", "angle"),
  ("09", "PASS", "CLT"),
  ("10", "INCONCLUSIVE", "KUBC"),
  ("11", "PASS", "Kok"),
  ("12", "INCONCLUSIVE", "cohesive"),
  ("13", "INCONCLUSIVE", "mesh"),
  ("14", "INCONCLUSIVE", "state"),
  ("15", "INCONCLUSIVE", "syntax"),
  ("16", "INCONCLUSIVE", "3.53h"),
)

#let chips = (
  ("MUMPS rebuild", atlantic, white, 1.65),
  ("LAW matrix", garnet, white, 1.65),
  ("orientation recipe", congaree, white, 2.16),
  ("Kok M1", horseshoe, white, 1.65),
  ("Kok M2", horseshoe, white, 1.65),
  ("Kok M3", horseshoe, white, 1.65),
  ("Kok M4", honeycomb, charcoal, 1.65),
  ("Kok M5", horseshoe, white, 1.65),
)

#align(center)[#text(size: 12pt, weight: "bold")[FEA_AP-PLY Final Stage Tally]]
#v(2mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  content((0.2, 6.42), [result], anchor: "east")
  content((0.2, 5.62), [criterion], anchor: "east")
  let cell-w = 0.92
  for (idx, item) in stages.enumerate() {
    let stage = item.at(0)
    let verdict = item.at(1)
    let criterion = item.at(2)
    let x = 0.55 + idx * 0.96
    let fill = if verdict == "PASS" { horseshoe } else if verdict == "FAIL" { garnet } else { honeycomb }
    let short = if verdict == "PASS" { "PASS" } else if verdict == "FAIL" { "FAIL" } else { "INC" }
    let text-fill = if verdict == "INCONCLUSIVE" { charcoal } else { white }
    content((x + 0.46, 7.05), [#stage], anchor: "center")
    rect((x, 6.05), (x + cell-w, 6.63), fill: fill, stroke: charcoal + 0.35pt)
    content((x + 0.46, 6.34), [#text(size: 6pt, fill: text-fill)[#short]], anchor: "center")
    rect((x, 5.26), (x + cell-w, 5.84), fill: white, stroke: black50 + 0.35pt)
    content((x + 0.46, 5.55), [#text(size: 5.5pt)[#criterion]], anchor: "center")
  }

  content((0.55, 4.45), [methodology discoveries and Kok port], anchor: "west")
  let x = 0.55
  for chip in chips {
    let label = chip.at(0)
    let fill = chip.at(1)
    let text-fill = chip.at(2)
    let width = chip.at(3)
    rect((x, 3.72), (x + width, 4.22), fill: fill, stroke: charcoal + 0.35pt)
    content((x + 0.08, 3.97), [#text(size: 7pt, fill: text-fill)[#label]], anchor: "west")
    x = x + width + 0.22
  }
})
