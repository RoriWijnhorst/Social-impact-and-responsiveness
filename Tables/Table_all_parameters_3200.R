## ===== 1) Truths (Sim. values) as in your new LaTeX table =====
sim_values <- data.frame(
  Parameter = c("B_0", "psi",
                "Sigma2_intercept", "Sigma2_psi", "Sigma2_epsilon",
                "Sigma2_x",
                "cor_1", "cor_2", "cor_3", "cor_4", "cor_5", "cor_6",
                "Sigma2_e", "Sigma2_ex"),
  Value = c(1.00, 0.30,
            0.20, 0.10, 0.01,
            1.00,
            -0.60, 0.00, -0.60, 0.60, -0.60, 0.00,
            0.60, 0.10),
  stringsAsFactors = FALSE
)

## ===== 2) Mapping from internal names -> LaTeX symbol + description =====
## Assumes your correlation parameters cor_1..cor_6 map as below.
## If your CSVs use different names, just tweak the left column ('Parameter').
param_map <- data.frame(
  Parameter = c("B_0", "psi",
                "Sigma2_intercept", "Sigma2_psi", "Sigma2_epsilon",
                "Sigma2_x",
                "cor_1", "cor_2", "cor_3", "cor_4", "cor_5", "cor_6",
                "Sigma2_e", "Sigma2_ex"),
  Symbol = c("$\\beta_0$", "$\\bar{\\psi}$",
             "$V_\\alpha$", "$V_\\psi$", "$V_\\epsilon$",
             "$V_\\chi$",
             "$r_{\\alpha\\epsilon}$", "$r_{\\alpha\\psi}$", "$r_{\\alpha\\chi}$",
             "$r_{\\epsilon\\psi}$", "$r_{\\chi\\psi}$", "$r_{\\chi\\epsilon}$",
             "$V_e$", "$V_{e_\\chi}$"),
  Description = c("Population mean", "Population response",
                  "Mean behaviour variance", "Responsiveness variance", "Residual impact variance",
                  "Impact trait variance",
                  "Corr: mean × res. impact", "Corr: mean × response", "Corr: mean × impact trait",
                  "Corr: res. impact × response", "Corr: impact trait × response", "Corr: impact trait × res. impact",
                  "Residual variance", "Measurement error"),
  stringsAsFactors = FALSE
)

## ===== 3) Your files and labels (as given) =====
files <- c("Output/res_3200_1600_2_1x_M1.csv", 
           "Output/res_3200_800_4_1x_M1.csv", 
           "Output/res_3200_400_8_1x_M1.csv", 
           "Output/res_3200_200_16_1x_M1.csv",
           "Output/res_3200_100_32_1x_M1.csv",
           "Output/res_3200_400_4_2x_M1.csv",
           "Output/res_3200_200_4_4x_M1.csv",
           "Output/res_3200_100_4_8x_M1.csv",
           "Output/res_3200_200_8_2x_M1.csv",
           "Output/res_3200_200_2_8x_M1.csv")

labels <-  c("1600_2_1x","800_4_1x", "400_8_1x", "200_16_1x", "100_32_1x", 
             "400_4_2x", "200_4_4x", "100_4_8x", 
             "200_8_2x", "200_2_8x")

## ===== 4) Build table of mean posterior medians per file =====
## Expects each CSV to have columns: Parameter, X50.  (X50. = posterior median)
## If your median column has a different name, change 'X50.' below.
results <- sim_values

for (i in seq_along(files)) {
  df <- read.csv(files[i], stringsAsFactors = FALSE)
  means <- aggregate(X50. ~ Parameter, data = df, FUN = mean)   # mean posterior median
  results <- merge(results, means, by = "Parameter", all.x = TRUE)
  colnames(results)[ncol(results)] <- labels[i]
}

## ===== 5) Attach symbols and descriptions; order rows like the LaTeX table =====
## Desired row order matching the LaTeX parameters above:
row_order <- c("B_0", "psi",
               "Sigma2_intercept", "Sigma2_psi", "Sigma2_epsilon",
               "Sigma2_x",
               "cor_1", "cor_2", "cor_3", "cor_4", "cor_5", "cor_6",
               "Sigma2_e", "Sigma2_ex")

results <- results[match(row_order, results$Parameter), ]

## Merge in symbol/description and arrange columns
results <- merge(param_map, results, by = "Parameter", all.y = TRUE, sort = FALSE)
results <- results[, c("Symbol", "Description", "Value", labels)]

## Optional: rename columns exactly like your LaTeX header later with kableExtra.
## For now, 'labels' stay as-is. You can subset/reorder to the 7 columns you want:
## Example mapping (edit to your exact 7):
## desired_cols <- c("400_2_1x","200_4_1x","100_8_1x","50_16_1x","100_4_2x","50_4_4x","100_2_4x")
## results_out <- results[, c("Symbol","Description","Value", desired_cols), drop = FALSE]

## View the final wide table
print(results)

desired_cols <- c(
  "1600_2_1x","800_4_1x","400_8_1x",
  "200_16_1x","100_32_1x",
  "400_4_2x","200_4_4x","100_4_8x",
  "200_8_2x"
)
stopifnot(length(desired_cols) == 9)

## ---- subset + format numbers ----
tab <- results[, c("Symbol","Description","Value", desired_cols)]
fmt3 <- function(x) if (is.numeric(x)) formatC(x, format = "f", digits = 3) else x
tab_fmt <- as.data.frame(lapply(tab, fmt3), stringsAsFactors = FALSE)

## ---- build the body rows in LaTeX ----
rows <- apply(tab_fmt, 1, function(row) paste(row, collapse = " & "))
rows_block <- paste(rows, "\\\\", collapse = "\n")

## ---- parse headers from labels like '200_16_1x' ----
parse_lab <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  ind <- parts[1]
  partners <- parts[2]
  reps <- parts[3]
  list(ind = ind, partners = partners, reps = reps)
}
hdr <- lapply(desired_cols, parse_lab)
hdr_ind <- vapply(hdr, `[[`, character(1), "ind")
hdr_partners <- vapply(hdr, `[[`, character(1), "partners")
hdr_reps <- vapply(hdr, `[[`, character(1), "reps")

## ---- make the three stacked header lines ----
## note: this matches your style (Parameter col blank; 'Individuals' etc. right over the next two cols)
header1 <- paste0(" & \\multicolumn{2}{r}{\\textit{Individuals}} & ",
                  paste(hdr_ind, collapse = " & "), " \\\\")
header2 <- paste0(" & \\multicolumn{2}{r}{\\textit{Social partners}} & ",
                  paste(hdr_partners, collapse = " & "), " \\\\")
header3 <- paste0(" & \\multicolumn{2}{r}{\\textit{Repeats}} & ",
                  paste(hdr_reps, collapse = " & "), " \\\\")

## ---- full LaTeX table text (12 columns total: ll + 10 c) ----
latex_output <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Mean model estimates (posterior medians) of 1000 simulated datasets under different sampling partitions.}\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{llcccccccccc}\n",  # ll + 10 'c' (Sim.value + 9 outcomes)
  "\\toprule\n",
  header1, "\n",
  header2, "\n",
  header3, "\n",
  "  \\midrule\n",
  "  Parameter & Description & Sim. value &  \\multicolumn{9}{c}{Model outcome} \\\\\n",
  "\\midrule\n",
  rows_block, "\n",
  "\\bottomrule\n",
  "\\end{tabular}%\n",
  "}\n",
  " \\label{fullpar}\n",
  "\\end{table}\n"
)

## ---- print so you can copy-paste into Overleaf ----
cat(latex_output)
