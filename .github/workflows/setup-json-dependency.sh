# Source after loading Spack.
setup_json_dependency() {
    local json_prefix json_config
    json_prefix="$(
        json_directory="$(mktemp -d)" || exit 1
        json_spec="$json_directory/spec.json"
        trap 'rm -f "$json_spec"; rmdir "$json_directory"' EXIT
        spack spec --json nlohmann-json@3.12.0 > "$json_spec" || exit 1
        spack install "$json_spec" >&2 || exit 1
        json_hash="$(spack python -c \
            'import sys; from spack.spec import Spec; print(Spec.from_json(sys.stdin.read()).dag_hash())' \
            < "$json_spec")" || exit 1
        spack location -i "/$json_hash"
    )" || return 1
    json_config="$(find "$json_prefix" -name nlohmann_jsonConfig.cmake -print -quit)"
    if [ -z "$json_config" ]; then
        echo "Missing nlohmann_jsonConfig.cmake in $json_prefix" >&2
        return 1
    fi
    export nlohmann_json_DIR="${json_config%/*}"
    printf 'nlohmann_json_DIR=%s\n' "$nlohmann_json_DIR" >> "$GITHUB_ENV"
}
