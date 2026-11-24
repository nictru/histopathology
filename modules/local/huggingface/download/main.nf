process HUGGINGFACE_DOWNLOAD {
    tag "${repo_id}/${filename}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4d/4d7393e98d9456e1e55d2d19eed365afa849a09fa75b297a58eb275f956957cc/data'
        : 'community.wave.seqera.io/library/huggingface_hub_python:837d1249362fe322'}"

    input:
    tuple val(meta), val(repo_type), val(repo_id), val(revision), val(filename)

    output:
    tuple val(meta), path("download/${filename}", arity: '1'), emit: resource
    path "versions.yml", emit: versions_huggingface_hub, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!(repo_type in ['model', 'dataset', 'space'])) {
        error("Unsupported Hugging Face repository type: ${repo_type}")
    }
    if (!(repo_id ==~ /[A-Za-z0-9_.-]+(\/[A-Za-z0-9_.-]+)?/)) {
        error("Invalid Hugging Face repository ID: ${repo_id}")
    }
    if (!(revision ==~ /[0-9a-f]{40}/)) {
        error('Specify a full Hugging Face commit hash as the revision.')
    }
    validateHubFilename(filename)
    """
    set -euo pipefail

    export HF_HOME=.huggingface
    export HF_TOKEN_PATH=.huggingface/token
    export HF_HUB_CACHE=.huggingface/hub
    export HF_HUB_DISABLE_XET=1
    export HF_HUB_DISABLE_PROGRESS_BARS=1
    export HF_HUB_DISABLE_IMPLICIT_TOKEN=0
    export HF_DEBUG=0

    hf download \\
        --repo-type '${repo_type}' \\
        --revision '${revision}' \\
        --local-dir download \\
        --quiet \\
        -- '${repo_id}' '${filename}'

    printf '{"%s":{"huggingface_hub":"%s"}}\\n' '${task.process}' "\$(hf --version)" > versions.yml
    """

    stub:
    validateHubFilename(filename)
    """
    mkdir -p -- "\$(dirname -- 'download/${filename}')"
    touch -- 'download/${filename}'
    echo '{}' > versions.yml
    """
}

def validateHubFilename(filename) {
    if (!(filename ==~ /[A-Za-z0-9._+ -]+(\/[A-Za-z0-9._+ -]+)*/) || filename.tokenize('/').any { part -> part in ['.', '..'] }) {
        error("Invalid Hugging Face file path: ${filename}")
    }
}
