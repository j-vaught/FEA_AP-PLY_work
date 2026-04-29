#set page(width: 175mm, height: auto, margin: 10mm)
#let rows = csv("../results/timeseries.csv")
#text(size: 12pt, weight: "bold")[Stage 08 ply-rotation starter probe]
#v(5pt)
#table(
  columns: (20mm, 32mm, 38mm, 28mm, 28mm),
  stroke: rgb("#5C5C5C"),
  [Theta], [Analytic Ex (GPa)], [Canonical status], [Proxy rc], [Verdict],
  ..rows.map(r => ([#r.at(0)], [#str(float(r.at(1)) / 1e9)], [#r.at(4)], [#r.at(6)], [#r.at(9)])).flatten(),
)
