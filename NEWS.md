# survcompare 0.3.1.9000

* Cross-validation wrappers now forward `inner_cv` to their model trainers.
  This includes the `survcompare()` entry point, Cox-Lasso, survival random
  forests, and both ensemble models.
  Ensemble and stacking Cox fits also receive the requested tuning folds.
  Previously, Cox-Lasso silently used five tuning folds and forest-based
  trainers used three regardless of the wrapper's requested value. The default
  `inner_cv = 3` therefore changes Cox-Lasso results; explicitly request five
  folds to preserve its previous effective setting. Custom forest fold counts
  now take effect as requested. Outer folds and ensemble out-of-fold prediction
  partitions are unchanged.
  Both learners within an ensemble now use its requested `inner_cv`. The
  historical mixed five-fold Cox / three-fold forest setting cannot be
  reproduced with that single argument; use the earlier package version to
  reproduce that specific ensemble configuration.
* Cox-Lasso uses the full `cv.glmnet(nfolds = ...)` argument name. On glmnet
  versions supporting it, the 5000-iteration limit is passed through `control`
  and Cox ties are explicitly Breslow. Older glmnet retains its supported
  direct `maxit` argument and Breslow-only behavior.
* Regression tests are included in source-package builds.
