#set page(width: 175mm, height: auto, margin: 10mm)
#let rows = csv("../results/timeseries.csv")
#text(size: 12pt, weight: "bold")[Stage 09 laminate CLT sweep]
#v(5pt)
#table(
  columns: (34mm, 24mm, 34mm, 34mm, 24mm),
  stroke: rgb("#5C5C5C"),
  [Layup], [Component], [CLT N/m], [FEM N/m], [Verdict],
  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(2)], [#r.at(3)], [#r.at(5)])).flatten(),
)
