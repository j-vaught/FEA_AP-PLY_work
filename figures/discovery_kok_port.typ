#import "@preview/cetz:0.3.4"
#set page(width: 166mm, height: 90mm, margin: 7mm)
#set text(font: "Libertinus Serif", size: 8pt, fill: rgb("#363636"))
#let garnet = rgb("#73000A")
#let horseshoe = rgb("#65780B")
#let honeycomb = rgb("#A49137")
#let atlantic = rgb("#466A9F")
#let congaree = rgb("#1F414D")
#let black10 = rgb("#ECECEC")
#let black30 = rgb("#C7C7C7")
#let black50 = rgb("#A2A2A2")
#let charcoal = rgb("#363636")
#let white = rgb("#FFFFFF")

#align(center)[#text(size: 11pt, weight: "bold")[Discovery: Kok Geometry Port]]
#v(1mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  let labels = ("M1 tow", "M2 ply", "M3 laminate", "M4 CLI/export")
  for i in range(4) {
    let x = 1.2 + i * 3.1
    rect((x, 4.6), (x + 2.25, 5.45), fill: if i < 3 { horseshoe } else { honeycomb }, stroke: charcoal + 0.4pt)
    content((x + 1.125, 5.02), [#text(fill: if i < 3 { white } else { charcoal })[#labels.at(i)]], anchor: "center")
    if i < 3 { line((x + 2.25, 5.02), (x + 3.1, 5.02), stroke: charcoal + 0.55pt) }
  }
  rect((1.2, 2.35), (10.85, 3.25), fill: black10, stroke: charcoal + 0.4pt)
  line((1.6, 2.8), (4.2, 2.8), stroke: garnet + 1.4pt)
  content((4.35, 2.8), [Stage 11 Ex = 16.7 GPa vs target 53.3 GPa], anchor: "west")
  content((6.2, 1.6), [Port is runnable through OpenRadioss; current clean-room geometry lacks the published stiffness.], anchor: "center")
})
