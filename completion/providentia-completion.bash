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
        # Two separate -name tests (ANDed by find) rather than one
        # "${value}*.conf" pattern: once $value has already reached into the
        # ".conf" extension (e.g. typed/completed up to "station.c"), a single
        # concatenated pattern would require a *second* ".conf" after that,
        # which no real filename can satisfy - completion would silently stop
        # matching right as the user gets close to finishing the name.
        while IFS= read -r -d '' f; do
            matches+=("${prefix}$(basename "$f")")
        done < <(find "$confdir" -maxdepth 1 -name "${value}*" -name "*.conf" -print0 2>/dev/null)
        COMPREPLY=("${matches[@]}")
        return
    fi

    COMPREPLY=()
}

# Strip '=' from word-breaks globally (once) so --conf=<TAB> / --config=<TAB>
# are seen as one token instead of being split on '='. This must NOT be
# toggled per-call from inside the completion function: bash inserts the
# completed match into the line *after* the function returns, using
# whatever COMP_WORDBREAKS holds at that point. Restoring the old value
# before returning (as a previous version of this script did) makes that
# insertion split on '=' again, which replaces only the fragment after the
# last '=' and duplicates the "--conf=" prefix in the result. Default
# COMP_WORDBREAKS on Linux includes '=', so this bit silently mangled every
# completion there; it happened to be masked on Mac because the reporter's
# terminal runs zsh, which never executes this file at all.
COMP_WORDBREAKS="${COMP_WORDBREAKS//=/}"

complete -F _providentia_complete -o nospace ./bin/providentia
complete -F _providentia_complete -o nospace providentia