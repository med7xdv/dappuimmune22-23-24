# two tailed test, more conservative; within immune
interactionbranchtest <- data.frame(
  "interaction_yes" = c(13, 82),
  "interaction_no" = c(11, 69),
  row.names = c("Significant","non-Significant"),
  stringsAsFactors = FALSE
)

colnames(interactionbranchtest) <- c("Interaction", "Non-Interaction")

fisher.test(t(interactionbranchtest))

# number of trees under selection within immune system
interactiongenetest <- data.frame(
  "interact" = c(7, 4),
  "nointeract" = c(7, 4),
  row.names = c("Significant","non-Significant"),
  stringsAsFactors = FALSE
)
colnames(interactiongenetest) <- c("immune", "nonimmune")

fisher.test(t(interactiongenetest))

# between non immune and immune
immunenonimmunebranches <- data.frame(
  "immune" = c(24, 151),
  "nonimmune" = c(13, 113),
  row.names = c("Significant","non-Significant"),
  stringsAsFactors = FALSE
)
colnames(immunenonimmunebranches) <- c("immune", "nonimmune")

fisher.test(t(immunenonimmunebranches))

# number of trees under selection between immune and non immune
immunenonimmunegenes <- data.frame(
  "immune" = c(14, 8),
  "nonimmune" = c(9, 8),
  row.names = c("Significant","non-Significant"),
  stringsAsFactors = FALSE
)
colnames(immunenonimmunegenes) <- c("immune", "nonimmune")

fisher.test(t(immunenonimmunegenes))

# difference between toll and imd
tollvimd <- data.frame(
  "toll" = c(14, 106),
  "imd" = c(10, 45),
  row.names = c("Significant","non-Significant"),
  stringsAsFactors = FALSE
)
colnames(tollvimd) <- c("toll", "imd")

fisher.test(t(tollvimd))
 

