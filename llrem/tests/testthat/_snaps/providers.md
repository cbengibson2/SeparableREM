# llm_call errors on non-provider input

    Code
      llm_call("not_a_provider", list())
    Condition
      Error:
      ! `provider` must be an <llm_provider> object.

