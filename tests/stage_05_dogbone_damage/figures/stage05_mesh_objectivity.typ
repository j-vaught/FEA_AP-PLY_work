#import "@preview/cetz:0.3.4"

#set page(width: 180mm, height: 116mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 9pt, fill: rgb("#363636"))
#let garnet = rgb("#73000A")
#let charcoal = rgb("#363636")
#let black10 = rgb("#ECECEC")
#let atlantic = rgb("#466A9F")
#let white = rgb("#FFFFFF")

#align(center)[#text(size: 11pt, weight: "bold")[Stage 05 LAW22 Mesh Objectivity]]
#v(2mm)
#figure(
  table(
    columns: 5,
    [Mesh], [Elements], [Peak force (N)], [Final EPSP max], [Wall-clock (s)],
    [coarse], [9504], [18637.07], [0.1228], [325.90],
    [medium], [27136], [18518.30], [0.1229], [709.65],
    [fine], [56000], [18575.06], [0.1228], [1340.45],
  ),
  caption: [OpenRadioss LAW22 notched-dogbone mesh sweep.]
)

#figure(
  table(
    columns: 4,
    [Check], [Value], [Tolerance], [Gating],
    [Coarse-medium pre RMSE], [4.843\%], [5.000\%], [yes],
    [Medium-fine pre RMSE], [3.736\%], [5.000\%], [yes],
    [Coarse-medium post RMSE], [4.957\%], [reported], [no],
    [Medium-fine post RMSE], [3.276\%], [reported], [no],
  ),
  caption: [Verdict: PASS. Post-onset divergence is reported because LAW22 is local CDM.]
)
