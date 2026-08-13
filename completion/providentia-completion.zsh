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
        prefix="${cur%%=*}="
        value="${cur#*=}"
        # (N) = nullglob (no error if nothing matches), (:t) = basename only
        matches=("${confdir}"/"${value}"*.conf(N:t))
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