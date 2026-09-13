# GitHub upload checklist

Upload the full repository contents, including the hidden `.gitignore` file.

Required:
- `.gitignore`
- `README.md`
- `run_analysis.R`
- `analysis/`
- `data/`
- `docs/`
- `outputs/`
- `reproducibility/`
- `CITATION.cff.template`

Before public release:
1. Rename `CITATION.cff.template` to `CITATION.cff`.
2. Replace `repository-code: "TO COMPLETE"` with the GitHub repository URL.
3. Add `date-released`.
4. Choose a software license and replace `LICENSE_TO_CHOOSE.md` with a real `LICENSE`.
5. Do not upload the five derived `.rds` files unless you explicitly decide to redistribute them.
6. Run the locked analysis once locally and confirm:
   `FINAL ANALYSIS COMPLETE — ALL LOCKED QC ASSERTIONS PASSED`.

7. Confirm Dongmei Chen's CRediT roles before final manuscript submission.
