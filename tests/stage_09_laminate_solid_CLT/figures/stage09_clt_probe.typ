#set page(width: 175mm, height: auto, margin: 10mm)
#let rows = csv("../results/timeseries.csv")
#text(size: 12pt, weight: "bold")[Stage 09 laminate CLT starter probe]
#v(5pt)
#table(
  columns: (34mm, 24mm, 34mm, 36mm, 26mm),
  stroke: rgb("#5C5C5C"),
  [Layup], [Component], [Reference N/m], [Canonical status], [Verdict],
  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(2)], [#r.at(5)], [#r.at(10)])).flatten(),
)
