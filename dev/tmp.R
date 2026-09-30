p <- c(0.15, 0.10, 0.10, 0.10, 0.10, 0.45)

n_precision_cont_po(p, beta = log(1.5), sd_x = 0.5, ratio_UL = 14, method = "whitehead")
n_precision_cont_po(p, beta = log(1.5), sd_x = 1, ratio_UL = 3.7, method = "whitehead")
n_precision_cont_po(p, beta = log(1.5), sd_x = 2, ratio_UL = 1.94, method = "whitehead")

n_precision_cont_po(p, beta = log(1.7), sd_x = 0.5, ratio_UL = 14, method = "whitehead")
n_precision_cont_po(p, beta = log(1.7), sd_x = 1, ratio_UL = 3.7, method = "whitehead")
n_precision_cont_po(p, beta = log(1.7), sd_x = 2, ratio_UL = 1.94, method = "whitehead")
