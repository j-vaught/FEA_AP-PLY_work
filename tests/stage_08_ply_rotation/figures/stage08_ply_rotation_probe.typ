#set page(width: 175mm, height: auto, margin: 10mm)
#let rows = csv("../results/timeseries.csv")
#text(size: 12pt, weight: "bold")[Stage 08 ply-rotation sweep]
#v(5pt)
#table(
  columns: (18mm, 34mm, 34mm, 28mm, 26mm),
  stroke: rgb("#5C5C5C"),
  [Theta], [Analytic Ex (Pa)], [Solver Ex (Pa)], [Error %], [Verdict],
  ..rows.map(r => ([#r.at(0)], [#r.at(1)], [#r.at(2)], [#r.at(3)], [#r.at(13)])).flatten(),
)
