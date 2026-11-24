#!/usr/bin/env bash
set -euo pipefail

export HF_HOME=.huggingface
export HF_TOKEN_PATH=.huggingface/token
export HF_HUB_CACHE=.huggingface/hub
export HF_XET_CACHE=.huggingface/xet
export HF_HUB_DISABLE_PROGRESS_BARS=1
export HF_HUB_DISABLE_IMPLICIT_TOKEN=0
export HF_DEBUG=0

hf download \
    --repo-type '${repo_type}' \
    --revision '${revision}' \
    --local-dir download \
    --quiet \
    -- '${repo_id}' '${filename}'

printf '{"%s":{"huggingface_hub":"%s"}}\n' '${task.process}' "\$(hf --version)" > versions.yml
