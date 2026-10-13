# Internal constructor ---------------------------------------------------------

new_provider <- function(name, model, call_fn) {
  structure(list(name = name, model = model, call = call_fn),
            class = "llm_provider")
}

# Parse the Retry-After header value (seconds) from a response, or return NULL.
retry_after_secs <- function(resp) {
  h <- httr2::resp_header(resp, "Retry-After")
  if (!is.null(h)) {
    n <- suppressWarnings(as.numeric(h))
    if (!is.na(n) && n > 0) {
      message("[llrem] Retry-After: ", n, "s")
      return(n)
    }
  }
  NULL  # fall through to backoff
}

# Execute an httr2 request, catching only network-level errors.
# Rate-limit and server errors are handled by req_retry inside the request.
execute_req <- function(req, extract_fn, max_retries = 3L) {
  for (attempt in seq_len(max_retries)) {
    resp <- tryCatch(
      httr2::req_perform(req),
      error = function(e) {
        message("[llrem] Request error (attempt ", attempt, "): ",
                conditionMessage(e))
        NULL
      }
    )
    if (!is.null(resp) && httr2::resp_status(resp) == 200L)
      return(extract_fn(resp))
    if (attempt < max_retries) Sys.sleep(min(2^(attempt + 2L), 120L))
  }
  warning("[llrem] LLM request failed after ", max_retries,
          " attempts; returning NA.", call. = FALSE)
  NA_character_
}

# Exported constructors -------------------------------------------------------

#' @export
print.llm_provider <- function(x, ...) {
  cat("<llm_provider>", x$name, "/", x$model, "\n")
  invisible(x)
}

#' Create a generic OpenAI-compatible LLM provider
#'
#' Suitable for any API that mirrors the OpenAI chat completions wire format.
#' Used internally by [llm_provider_openai()], [llm_provider_grok()],
#' [llm_provider_gemini()], and [llm_provider_ollama()].
#'
#' @param api_key API key string.
#' @param model Model identifier.
#' @param base_url Full chat-completions endpoint URL.
#' @param name Display name used in `print()` output and plot labels.
#' @param max_retries Maximum request attempts before returning `NA`.
#'
#' @return An `llm_provider` object.
#' @export
llm_provider_openai_compat <- function(api_key,
                                        model,
                                        base_url,
                                        name        = "OpenAI-compatible",
                                        max_retries = 5L) {
  force(api_key); force(model); force(base_url); force(name); force(max_retries)
  call_fn <- function(messages, max_tokens = 15L, temperature = 0.7) {
    req <- httr2::request(base_url) |>
      httr2::req_auth_bearer_token(api_key) |>
      httr2::req_headers("Content-Type" = "application/json") |>
      httr2::req_body_json(list(
        model       = model,
        messages    = messages,
        max_tokens  = max_tokens,
        temperature = temperature
      )) |>
      httr2::req_retry(
        max_tries    = max_retries,
        is_transient = \(r) httr2::resp_status(r) %in% c(429L, 500L, 503L),
        backoff      = \(i) min(2^i, 60L),
        after        = retry_after_secs
      )
    execute_req(req, \(resp) {
      trimws(httr2::resp_body_json(resp)$choices[[1L]]$message$content)
    }, max_retries = 3L)
  }
  new_provider(name, model, call_fn)
}

#' Create an OpenAI LLM provider
#'
#' @param api_key Defaults to `Sys.getenv("OPENAI_API_KEY")`.
#' @param model Default `"gpt-4o"`.
#' @param max_retries Maximum request attempts.
#'
#' @return An `llm_provider` object.
#' @export
llm_provider_openai <- function(api_key     = Sys.getenv("OPENAI_API_KEY"),
                                 model       = "gpt-4o",
                                 max_retries = 5L) {
  llm_provider_openai_compat(api_key, model,
                               "https://api.openai.com/v1/chat/completions",
                               name = "OpenAI", max_retries = max_retries)
}

#' Create an xAI Grok LLM provider
#'
#' @param api_key Defaults to `Sys.getenv("XAI_API_KEY")`.
#' @param model Default `"grok-3"`.
#' @param max_retries Maximum request attempts.
#'
#' @return An `llm_provider` object.
#' @export
llm_provider_grok <- function(api_key     = Sys.getenv("XAI_API_KEY"),
                               model       = "grok-3",
                               max_retries = 5L) {
  llm_provider_openai_compat(api_key, model,
                               "https://api.x.ai/v1/chat/completions",
                               name = "xAI/Grok", max_retries = max_retries)
}

#' Create a Google Gemini LLM provider
#'
#' Uses Google's OpenAI-compatible endpoint.
#'
#' @param api_key Defaults to `Sys.getenv("GEMINI_API_KEY")`.
#' @param model Default `"gemini-2.0-flash"`.
#' @param max_retries Maximum request attempts.
#'
#' @return An `llm_provider` object.
#' @export
llm_provider_gemini <- function(api_key     = Sys.getenv("GEMINI_API_KEY"),
                                 model       = "gemini-2.0-flash",
                                 max_retries = 5L) {
  llm_provider_openai_compat(
    api_key, model,
    "https://generativelanguage.googleapis.com/v1beta/openai/chat/completions",
    name = "Google/Gemini", max_retries = max_retries
  )
}

#' Create a local or remote Ollama LLM provider
#'
#' @param model       Model tag, e.g. `"llama3.2"` or `"mistral"`.
#' @param host        Hostname or IP of the Ollama server.  Defaults to
#'   `"localhost"`.  Override to target a remote machine, e.g.
#'   `host = "192.168.1.65"`.
#' @param port        Port the Ollama server is listening on (default 11434).
#' @param base_url    Full endpoint URL.  When supplied, `host` and `port` are
#'   ignored — useful as a manual escape hatch.
#' @param max_retries Maximum request attempts.
#'
#' @return An `llm_provider` object.
#' @export
llm_provider_ollama <- function(model       = "llama3.2",
                                 host        = "localhost",
                                 port        = 11434L,
                                 base_url    = NULL,
                                 max_retries = 3L) {
  url <- if (!is.null(base_url)) base_url else
    sprintf("http://%s:%d/v1/chat/completions", host, as.integer(port))
  llm_provider_openai_compat("ollama", model, url,
                               name = "Ollama", max_retries = max_retries)
}

#' Create an Anthropic Claude LLM provider
#'
#' Anthropic's API separates the system prompt as a top-level field and uses a
#' different response shape than the OpenAI format. This constructor handles
#' both translations internally.
#'
#' @param api_key Defaults to `Sys.getenv("ANTHROPIC_API_KEY")`.
#' @param model Default `"claude-3-5-sonnet-20241022"`.
#' @param anthropic_version API version header string.
#' @param max_retries Maximum request attempts.
#'
#' @return An `llm_provider` object.
#' @export
llm_provider_anthropic <- function(api_key           = Sys.getenv("ANTHROPIC_API_KEY"),
                                    model             = "claude-sonnet-4-6",
                                    anthropic_version = "2023-06-01",
                                    max_retries       = 5L) {
  force(api_key); force(model); force(anthropic_version); force(max_retries)
  call_fn <- function(messages, max_tokens = 15L, temperature = 0.7) {
    is_sys   <- vapply(messages, \(m) identical(m$role, "system"), logical(1L))
    sys_text <- paste(vapply(messages[is_sys], `[[`, character(1L), "content"),
                      collapse = "\n")
    body <- list(model = model, max_tokens = max_tokens,
                 temperature = temperature, messages = messages[!is_sys])
    if (nchar(sys_text) > 0L) body$system <- sys_text
    req <- httr2::request("https://api.anthropic.com/v1/messages") |>
      httr2::req_headers("x-api-key"         = api_key,
                         "anthropic-version" = anthropic_version,
                         "Content-Type"      = "application/json") |>
      httr2::req_body_json(body) |>
      httr2::req_retry(
        max_tries    = max_retries,
        is_transient = \(r) httr2::resp_status(r) %in% c(429L, 500L, 529L),
        backoff      = \(i) min(2^i, 60L),
        after        = retry_after_secs
      )
    execute_req(req, \(resp) {
      trimws(httr2::resp_body_json(resp)$content[[1L]]$text)
    }, max_retries = 3L)
  }
  new_provider("Anthropic", model, call_fn)
}

#' Create an AWS Bedrock LLM provider
#'
#' Uses Bedrock's Converse API with AWS SigV4 signing. Supports any model
#' available in your Bedrock account, including Anthropic Claude and Meta
#' Llama families.
#'
#' Credentials are resolved in order: explicit arguments → environment variables
#' (`AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_SESSION_TOKEN`,
#' `AWS_DEFAULT_REGION`).
#'
#' Common model IDs:
#' * `"anthropic.claude-haiku-4-5-20251001-v1:0"`
#' * `"anthropic.claude-3-5-sonnet-20241022-v2:0"`
#' * `"meta.llama3-70b-instruct-v1:0"`
#' * `"meta.llama3-8b-instruct-v1:0"`
#' * `"mistral.mistral-7b-instruct-v0:2"`
#'
#' @param model       Bedrock model ID string (see above).
#' @param region      AWS region. Defaults to `AWS_DEFAULT_REGION` env var,
#'   then `"us-east-1"`.
#' @param access_key  AWS access key ID. Defaults to `AWS_ACCESS_KEY_ID`.
#' @param secret_key  AWS secret access key. Defaults to `AWS_SECRET_ACCESS_KEY`.
#' @param session_token Optional STS session token. Defaults to
#'   `AWS_SESSION_TOKEN` (ignored if unset).
#' @param name        Display name for print output and plot labels.
#' @param max_retries Maximum request attempts.
#'
#' @return An `llm_provider` object.
#' @export
llm_provider_bedrock <- function(model,
                                  region        = Sys.getenv("AWS_DEFAULT_REGION",
                                                             unset = "us-east-1"),
                                  access_key    = Sys.getenv("AWS_ACCESS_KEY_ID"),
                                  secret_key    = Sys.getenv("AWS_SECRET_ACCESS_KEY"),
                                  session_token = Sys.getenv("AWS_SESSION_TOKEN"),
                                  name          = "AWS Bedrock",
                                  max_retries   = 5L) {
  force(model); force(region); force(access_key); force(secret_key)
  force(session_token); force(name); force(max_retries)

  tok <- if (nchar(session_token) > 0L) session_token else NULL

  call_fn <- function(messages, max_tokens = 15L, temperature = 0.7) {
    # Converse API: separate system messages from turns
    is_sys   <- vapply(messages, \(m) identical(m$role, "system"), logical(1L))
    sys_text <- paste(vapply(messages[is_sys], `[[`, character(1L), "content"),
                      collapse = "\n")
    turns <- lapply(messages[!is_sys], function(m) {
      list(role    = m$role,
           content = list(list(type = "text", text = m$content)))
    })

    body <- list(
      messages        = turns,
      inferenceConfig = list(maxTokens   = as.integer(max_tokens),
                             temperature = temperature)
    )
    if (nchar(sys_text) > 0L)
      body$system <- list(list(type = "text", text = sys_text))

    url <- sprintf(
      "https://bedrock-runtime.%s.amazonaws.com/model/%s/converse",
      region, utils::URLencode(model, reserved = TRUE)
    )

    req <- httr2::request(url) |>
      httr2::req_headers("Content-Type" = "application/json") |>
      httr2::req_body_json(body) |>
      httr2::req_auth_aws_v4(
        aws_access_key_id     = access_key,
        aws_secret_access_key = secret_key,
        aws_session_token     = tok,
        aws_region            = region,
        aws_service           = "bedrock"
      ) |>
      httr2::req_retry(
        max_tries    = max_retries,
        is_transient = \(r) httr2::resp_status(r) %in% c(429L, 500L, 503L),
        backoff      = \(i) min(2^i, 60L),
        after        = retry_after_secs
      )

    execute_req(req, \(resp) {
      trimws(httr2::resp_body_json(resp)$output$message$content[[1L]]$text)
    }, max_retries = 3L)
  }
  new_provider(name, model, call_fn)
}

#' Call an LLM provider
#'
#' Generic dispatch function. Passes a list of role/content message pairs to
#' the provider's internal call function and returns the text response.
#'
#' @param provider An `llm_provider` object.
#' @param messages List of role/content pairs in OpenAI format, e.g.
#'   `list(list(role = "system", content = "..."), list(role = "user", content = "..."))`.
#' @param max_tokens Token budget for the response.
#' @param temperature Sampling temperature.
#'
#' @return A character string, or `NA_character_` on failure.
#' @export
llm_call <- function(provider, messages, max_tokens = 15L, temperature = 0.7) {
  if (!inherits(provider, "llm_provider"))
    stop("`provider` must be an <llm_provider> object.", call. = FALSE)
  provider$call(messages = messages, max_tokens = max_tokens,
                temperature = temperature)
}
