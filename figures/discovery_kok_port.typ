#import "@preview/cetz:0.3.4"
#set page(width: 166mm, height: 90mm, margin: 7mm)
#set text(font: "Libertinus Serif", size: 8pt, fill: rgb("#363636"))
#let garnet = rgb("#73000A")
#let horseshoe = rgb("#65780B")
#let honeycomb = rgb("#A49137")
#let atlantic = rgb("#466A9F")
#let congaree = rgb("#1F414D")
#let black10 = rgb("#ECECEC")
#let charcoal = rgb("#363636")
#let white = rgb("#FFFFFF")

#let labels = (
  ("M1 tow", horseshoe, white),
  ("M2 ply", horseshoe, white),
  ("M3 laminate", horseshoe, white),
  ("M4 CLI/export", honeycomb, charcoal),
  ("M5 Stage 11 PASS", horseshoe, white),
)

#align(center)[#text(size: 11pt, weight: "bold")[Discovery: Kok Geometry Port]]
#v(1mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  for (i, item) in labels.enumerate() {
    let label = item.at(0)
    let fill = item.at(1)
    let text-fill = item.at(2)
    let x = 0.6 + i * 2.95
    rect((x, 4.65), (x + 2.25, 5.50), fill: fill, stroke: charcoal + 0.4pt)
    content((x + 1.125, 5.08), [#text(size: 6.8pt, fill: text-fill)[#label]], anchor: "center")
    if i < labels.len() - 1 {
      line((x + 2.25, 5.08), (x + 2.95, 5.08), stroke: charcoal + 0.55pt)
    }
  }

  rect((0.8, 2.40), (13.80, 3.28), fill: black10, stroke: charcoal + 0.4pt)
  line((1.2, 2.84), (4.1, 2.84), stroke: garnet + 1.2pt)
  content((4.3, 2.84), [Stage 11 PASS: Ex 49.04 / Ey 53.29 / Gxy 22.53 GPa], anchor: "west")

  rect((1.8, 1.15), (12.8, 1.95), fill: white, stroke: charcoal + 0.35pt)
  content((7.3, 1.55), [The clean-room Kok port now runs through OpenRadioss and lands inside the 10% Kok 2022 modulus gate.], anchor: "center")
})
