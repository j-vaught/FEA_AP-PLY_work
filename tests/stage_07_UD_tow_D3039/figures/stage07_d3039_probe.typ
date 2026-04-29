#set page(width: 170mm, height: auto, margin: 10mm)
#let rows = csv("../results/timeseries.csv")
#text(size: 12pt, weight: "bold")[Stage 07 D3039 starter probe]
#v(5pt)
#table(
  columns: (22mm, 22mm, 42mm, 32mm, 30mm),
  stroke: rgb("#5C5C5C"),
  [Run], [Theta], [Canonical status], [Proxy rc], [Verdict],
  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(5)], [#r.at(7)], [#r.at(10)])).flatten(),
)
