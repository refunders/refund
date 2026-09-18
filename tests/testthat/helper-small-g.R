#--------------------------------------
# Helper: muffle the S-C small-cluster-count guard (pffr_small_G_warning) for
# fixtures across the suite that deliberately use few clusters to keep fits
# fast while testing mechanics unrelated to sandwich calibration (AR(1)
# wiring, cache invalidation, autopolicy resolution, argument plumbing,
# etc.). Any OTHER warning from the wrapped call still propagates normally.
# The guard's own behaviour is exercised directly in test-pffr-small-g.R,
# which installs its own (more specific, innermost) handler and is
# unaffected by this helper.
#--------------------------------------

quiet_pffr <- function(...) {
  withCallingHandlers(
    pffr(...),
    pffr_small_G_warning = function(w) invokeRestart("muffleWarning")
  )
}
