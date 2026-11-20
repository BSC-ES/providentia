_providentia_completion_dir="$(cd "$(dirname "${(%):-%N}")" && pwd)"

_providentia_zsh_complete() {
    local cur confdir word prefix value
    local -a matches

    cur="${words[CURRENT]}"
    confdir="${_providentia_completion_dir}/../configurations"

    # if --config_dir=... already appears earlier on the command line, use it
    # instead of the default. Resolved relative to the current working
    # directory unless absolute.
    for word in "${words[@]}"; do
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
        local -a all
        prefix="${cur%%=*}="
        value="${cur#*=}"
        # (N) = nullglob (no error if nothing matches), (:t) = basename only.
        # List by prefix first, then filter for the ".conf" suffix separately
        # (rather than a single "${value}*.conf" pattern): once $value has
        # already reached into the extension itself (e.g. "station.c"), a
        # concatenated pattern would require a *second* ".conf" after that,
        # which no real filename can satisfy - completion would silently
        # stop matching right as the user gets close to finishing the name.
        all=("${confdir}"/"${value}"*(N:t))
        matches=(${(M)all:#*.conf})
        compadd -P "$prefix" -- "${matches[@]}"
    fi
}

# only initialise completion if nothing has already done so (avoids clashing
# with oh-my-zsh, Prezto, etc. which typically call compinit themselves)
if ! (( ${+functions[compdef]} )); then
    autoload -Uz compinit
    compinit -C
fi

compdef _providentia_zsh_complete providentia
compdef _providentia_zsh_complete ./bin/providentia
compdef _providentia_zsh_complete bin/providentia