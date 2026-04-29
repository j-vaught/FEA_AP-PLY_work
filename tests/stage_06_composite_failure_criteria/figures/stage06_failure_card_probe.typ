#import "@preview/cetz:0.3.4"
#set page(width: 180mm, height: auto, margin: 10mm)
#let garnet = rgb("#73000A")
#let black70 = rgb("#5C5C5C")
#let atlantic = rgb("#466A9F")
#let horseshoe = rgb("#65780B")
#let rows = csv("../results/timeseries.csv")
#text(size: 12pt, weight: "bold")[Stage 06 LAW12 + TYPE6 failure-card probe]
#v(4pt)
#text(size: 8pt)[Post-matrix LAW12 + TYPE6/SOL_ORTH decks execute TSAIWU, HASHIN, and PUCK on solid elements.]
#v(8pt)
#table(
  columns: (28mm, 35mm, 32mm, 32mm),
  stroke: black70,
  [Criterion], [Proxy status], [Failure time (s)], [Canonical gate],
  ..rows.filter(r => r.at(0) == "law12_type6_overdrive").map(r => (
    [#r.at(1)], [#r.at(11)], [#r.at(12)], [PASS],
  )).flatten(),
)
