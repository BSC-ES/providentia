# Dispatcher only - kept free of any shell-specific syntax that the other
# shell's parser would choke on, since both bash and zsh fully parse this
# file's syntax even for branches they don't execute.
if [[ -n "${ZSH_VERSION:-}" ]]; then
    _providentia_dispatch_dir="$(cd "$(dirname "${(%):-%N}")" && pwd)"
    source "${_providentia_dispatch_dir}/providentia-completion.zsh"
else
    _providentia_dispatch_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${_providentia_dispatch_dir}/providentia-completion.bash"
fi
unset _providentia_dispatch_dir