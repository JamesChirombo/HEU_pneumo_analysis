## function to calculate diagnostic tests for the HUU group

calculate_diagnostics <- function(outcome,exposure){
  xtab <- table(factor(epiHuu[[exposure]],levels = c("1","0")),factor(epiHuu[[outcome]],levels = c("1","0")))
  xtabTest <- epiR::epi.tests(xtab)
  res <- data.frame(Est=c(xtabTest$elements$ppv,
                          xtabTest$elements$npv,
                          xtabTest$elements$lrpos,
                          youden=xtabTest$elements$youden$est),
                    Lower=c(xtabTest$elements$ppv.low,
                            xtabTest$elements$npv.low,
                            xtabTest$elements$lrpos.low,
                            xtabTest$elements$youden$lower),
                    Upper=c(xtabTest$elements$ppv.up,
                            xtabTest$elements$npv.up,
                            xtabTest$elements$lrpos.up,
                            xtabTest$elements$youden$upper))
  rownames(res) <- c("ppv","npv","lr","youden")
  return(res)
}