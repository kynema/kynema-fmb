# Source after loading Spack.
setup_json_dependency() {
    local json_prefix json_config
    spack install nlohmann-json@3.12.0
    json_prefix="$(spack location -i nlohmann-json@3.12.0)"
    json_config="$(find "$json_prefix" -name nlohmann_jsonConfig.cmake -print -quit)"
    if [ -z "$json_config" ]; then
        echo "Missing nlohmann_jsonConfig.cmake in $json_prefix" >&2
        return 1
    fi
    export nlohmann_json_DIR="${json_config%/*}"
    printf 'nlohmann_json_DIR=%s\n' "$nlohmann_json_DIR" >> "$GITHUB_ENV"
}
