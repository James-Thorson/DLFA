
library(dsem)

nt = 40

data = data.frame(
  area1 = rep(NA, nt),
  area2 = rep(NA, nt),
  area3 = rep(NA, nt),
  F1 = rnorm(nt)
)

sem = "
  F1 -> area1, 1, NA, 0.707
  F1 -> area2, 0, NA, 0
  F1 -> area3, 0, NA, 0

  F1 -> area1, 1, NA, 0
  F1 -> area2, 1, NA, 0.707
  F1 -> area3, 1, NA, 0

  F1 -> area1, 2, NA, 0
  F1 -> area2, 2, NA, 0
  F1 -> area3, 2, NA, 0.707

  F1 <-> F1, 0, NA, 1
  area1 <-> area1, 0, NA, 0.707
  area2 <-> area2, 0, NA, 0.707
  area3 <-> area3, 0, NA, 0.707
"

simulator = dsem(
  sem = sem,
  tsdata = ts(data),
  control = dsem_control(

  )
)
