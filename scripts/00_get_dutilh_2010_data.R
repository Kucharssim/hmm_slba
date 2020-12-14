library(here)
library(readr)
load(url("http://www.gillesdutilh.com/materials/PTM/PTMdat.Rdata"))

table(experiment = td$experiment, subject = td$subject)

for(s in unique(td$subject)){
  sub <- subset(td, subject == s)
  sub$rt <- sub$RT / 1000
  sub$accumulator <- ifelse(sub$correct == 1, 1, 2)
  sub <- subset(sub, select = c("pacc", "rt", "accumulator"))
  readr::write_csv(x = sub, path = here("data", sprintf("dutilh_2010_subject_%s.csv", s)))
}
