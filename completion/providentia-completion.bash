_providentia_completion_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

_providentia_complete() {
    local cur confdir word

    cur="${COMP_WORDS[COMP_CWORD]}"
    confdir="${_providentia_completion_dir}/../configurations"

    # if --config_dir=... already appears earlier on the command line, use it
    # instead of the default. Resolved relative to the current working
    # directory unless absolute.
    for word in "${COMP_WORDS[@]}"; do
        if [[ "$word" == --config_dir=* ]]; then
            local override="${word#--config_dir=}"
            if [[ "$override" == /* ]]; then
                confdir="$override"
            else
                confdir="${PWD}/${override}"
            fi
        fi
    done

    if [[ "$cur" == --conf=* || "$cur" == --config=* ]]; then
        local prefix="${cur%%=*}="
        local value="${cur#*=}"
        local matches=()
        while IFS= read -r -d '' f; do
            matches+=("${prefix}$(basename "$f")")
        done < <(find "$confdir" -maxdepth 1 -name "${value}*.conf" -print0 2>/dev/null)
        COMPREPLY=("${matches[@]}")
        return
    fi

    COMPREPLY=()
}

# Strip '=' from word-breaks locally so --conf=<TAB> / --config=<TAB>
# are seen as one token instead of being split on '='.
_providentia_wrapper() {
    local old_wordbreaks="$COMP_WORDBREAKS"
    COMP_WORDBREAKS="${COMP_WORDBREAKS//=/}"
    _providentia_complete
    COMP_WORDBREAKS="$old_wordbreaks"
}

complete -F _providentia_wrapper -o nospace ./bin/providentia
complete -F _providentia_wrapper -o nospace providentia