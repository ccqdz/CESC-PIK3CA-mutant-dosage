required <- c("data.table", "ggplot2", "patchwork", "scales", "grid")
to_install <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
message("Required packages are available.")
