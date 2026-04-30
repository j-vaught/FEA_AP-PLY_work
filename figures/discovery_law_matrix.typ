#import "@preview/cetz:0.3.4"
#set page(width: 152mm, height: 90mm, margin: 7mm)
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

#align(center)[#text(size: 11pt, weight: "bold")[Discovery: LAW x PROP Compatibility]]
#v(1mm)
#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  let laws = ("LAW12", "LAW14", "LAW25", "LAW28", "LAW53", "LAW128")
  let props = ("TYPE14", "TYPE6", "HASHIN", "PUCK", "TSAIWU")
  for i in range(6) { content((1.75 + i * 1.45, 5.9), [#laws.at(i)], anchor: "center") }
  for j in range(5) { content((0.20, 5.15 - j * 0.78), [#props.at(j)], anchor: "west") }
  for i in range(6) {
    rect((1.30 + i * 1.45, 4.85), (2.20 + i * 1.45, 5.45), fill: garnet, stroke: charcoal + 0.3pt)
    content((1.75 + i * 1.45, 5.15), [#text(fill: white)[B3047]], anchor: "center")
    for j in range(4) {
      rect((1.30 + i * 1.45, 4.07 - j * 0.78), (2.20 + i * 1.45, 4.67 - j * 0.78), fill: horseshoe, stroke: charcoal + 0.3pt)
      content((1.75 + i * 1.45, 4.37 - j * 0.78), [#text(fill: white)[OK]], anchor: "center")
    }
  }
  content((0.55, 0.75), [Result: use LAW12 + TYPE6/SOL_ORTH for solid composite decks; TYPE14 remains blocked.], anchor: "west")
})
