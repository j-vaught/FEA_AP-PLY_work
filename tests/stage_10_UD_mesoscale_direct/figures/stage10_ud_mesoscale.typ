#set page(width: 170mm, height: auto, margin: 10mm)
#let rows = csv("../results/effective_moduli.csv")
#text(size: 12pt, weight: "bold")[Stage 10 UD mesoscale moduli]
#v(5pt)
#table(
  columns: (34mm, 34mm, 34mm, 28mm, 24mm),
  stroke: rgb("#5C5C5C"),
  [Case], [FEM Pa], [Target Pa], [Error %], [Verdict],
  ..rows.map(r => ([#r.at(0)], [#r.at(2)], [#r.at(3)], [#r.at(4)], [#r.at(8)])).flatten(),
)
