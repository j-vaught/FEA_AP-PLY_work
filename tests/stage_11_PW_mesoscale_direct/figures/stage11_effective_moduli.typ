#set page(width: 170mm, height: auto, margin: 10mm)
#let rows = csv("../results/effective_moduli.csv")
#text(size: 12pt, weight: "bold")[Stage 11 AP-PLY effective moduli]
#v(5pt)
#table(
  columns: (32mm, 24mm, 28mm, 28mm, 26mm),
  stroke: rgb("#5C5C5C"),
  [Case], [Metric], [FEM GPa], [Target GPa], [Error %],
  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(2)], [#r.at(3)], [#r.at(4)])).flatten(),
)
