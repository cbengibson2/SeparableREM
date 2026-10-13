#' Create a prompt configuration object
#'
#' Bundles everything that affects prompt construction into a single reusable
#' object. One `rem_cfg` can be passed to multiple [make_query_roleplay_sender()],
#' [make_query_predict()], and [make_query_roleplay_random()] calls. Swap
#' `provider` to replicate a study with a different model.
#'
#' @param provider An `llm_provider` object created by one of the
#'   `llm_provider_*()` functions.
#' @param context One or two sentences describing the network domain in second
#'   person (used as the system message).
#' @param action_verb Verb describing a directed event, e.g. `"attacked"`,
#'   `"emailed"`, `"called"`.
#' @param actor_noun Capitalised singular noun for nodes, e.g. `"Country"`,
#'   `"Person"`, `"Actor"`.
#' @param node_labels Optional named character vector mapping integer node IDs
#'   (as strings) to readable labels, e.g. `c("1" = "USA", "2" = "Russia")`.
#'   `NULL` generates generic `"<actor_noun> <id>"` labels.
#' @param history_window Maximum number of past events included in each prompt.
#'   `NULL` includes the full history.
#'
#' @return A `rem_cfg` object.
#' @export
rem_cfg <- function(provider,
                    context,
                    action_verb,
                    actor_noun,
                    node_labels    = NULL,
                    history_window = 30L) {
  if (!inherits(provider, "llm_provider"))
    stop("`provider` must be an <llm_provider> object.", call. = FALSE)
  structure(
    list(provider = provider, context = context, action_verb = action_verb,
         actor_noun = actor_noun, node_labels = node_labels,
         history_window = history_window),
    class = "rem_cfg"
  )
}

#' @export
print.rem_cfg <- function(x, ...) {
  cat("<rem_cfg>\n")
  cat("  Provider :", x$provider$name, "/", x$provider$model, "\n")
  cat("  Action   : <actor>", x$action_verb, "<actor>\n")
  cat("  Noun     :", x$actor_noun, "\n")
  cat("  History  :",
      if (is.null(x$history_window)) "full" else x$history_window,
      "events\n")
  if (!is.null(x$node_labels))
    cat("  Labels   :", length(x$node_labels), "named nodes\n")
  invisible(x)
}
