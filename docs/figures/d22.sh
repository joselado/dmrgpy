#!/usr/bin/env bash
# Run the figure calculations on D22 from this machine.
#
#   bash docs/figures/d22.sh push            # copy the code to D22
#   bash docs/figures/d22.sh run <id>...     # start each job, detached
#   bash docs/figures/d22.sh status          # RUNNING / DONE / FAILED per figure
#   bash docs/figures/d22.sh pull <id>...    # fetch data.npz, compute.log, PNGs
#   bash docs/figures/d22.sh plot <id>...    # run plot.py here
#
# push copies the git-tracked files and the (untracked) docs/figures/ to
# the remote checkout, without --delete. It never sends compiled files
# (*.so, *.a, mpscpp*/ITensor/lib), since the remote v3 extension is
# built there, and never sends results (data.npz, logs, DONE, PNGs) or the
# slot locks, so a push cannot overwrite a remote result. Results only
# flow back with pull. Several `run` jobs queue on the three slots of
# run.sh, so it is safe to start more than three at once.
set -eu
HOST="${D22_HOST:-ladovj1@D22-0161.org.aalto.fi}"
RDIR="${D22_DIR:-/l/ladovj1/dmrgpy-guide}"
RPY="${D22_PYTHON:-/l/ladovj1/venv/bin/python}"
# one SMT sibling of three P-cores (0, 1, 4); CPUs 4-7 (P-cores 2, 3) are
# taken by another project's jobs on D22, which do not share our locks
RCPUS="${D22_CPUS:-0 2 8}"
SSH=(ssh -o BatchMode=yes)
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
NOSEND=(--exclude='*.so' --exclude='*.a' --exclude='*.o'
        --exclude='mpscpp*/ITensor/lib/' --exclude='__pycache__/')
RESULTS=(--exclude='.locks/' --exclude='DONE' --exclude='FAILED'
         --exclude='RUNNING' --exclude='compute.log' --exclude='run.out'
         --exclude='data.npz' --exclude='*.png')

cmd="${1:-}"; shift || true
case "$cmd" in
push)
    cd "$repo"
    { git rev-parse HEAD | tr -d '\n'
      [ -n "$(git status --porcelain -- src)" ] && printf -- '-dirty'
      echo; } > "$here/COMMIT"
    # docs/figures/ goes in the second call, which never sends results
    git ls-files -z | grep -zv '^docs/figures/' | rsync -a -0 --files-from=- "${NOSEND[@]}" \
        -e "${SSH[*]}" ./ "$HOST:$RDIR/"
    rsync -a "${NOSEND[@]}" "${RESULTS[@]}" -e "${SSH[*]}" \
        docs/figures/ "$HOST:$RDIR/docs/figures/"
    echo "pushed $(cat "$here/COMMIT") to $HOST:$RDIR"
    ;;
run)
    [ $# -gt 0 ] || { echo "usage: d22.sh run <id>..." >&2; exit 2; }
    for id in "$@"; do
        [ -f "$here/$id/compute.py" ] || { echo "no $here/$id/compute.py" >&2; exit 2; }
        "${SSH[@]}" -f "$HOST" "cd $RDIR && test -f docs/figures/$id/compute.py && \
            rm -f docs/figures/$id/DONE docs/figures/$id/FAILED && \
            echo queued > docs/figures/$id/RUNNING && \
            PYTHON=$RPY FIGURE_CPUS=\"$RCPUS\" nohup bash docs/figures/run.sh $id \
            > docs/figures/$id/run.out 2>&1 &" < /dev/null
        echo "started $id"
    done
    ;;
status)
    "${SSH[@]}" "$HOST" "cd $RDIR/docs/figures && for d in */; do d=\${d%/}; \
        [ -f \$d/compute.py ] || continue; \
        if [ -f \$d/DONE ]; then s=\"DONE    \$(cat \$d/DONE)\"; \
        elif [ -f \$d/FAILED ]; then s=\"FAILED  \$(cat \$d/FAILED)\"; \
        elif [ -f \$d/RUNNING ]; then s=\"RUNNING \$(cat \$d/RUNNING)\"; \
        else s='-'; fi; printf '%-45s %s\n' \$d \"\$s\"; done" < /dev/null
    ;;
pull)
    [ $# -gt 0 ] || { echo "usage: d22.sh pull <id>..." >&2; exit 2; }
    for id in "$@"; do
        rsync -a -e "${SSH[*]}" --include='data.npz' --include='compute.log' \
            --include='run.out' --include='*.png' --include='DONE' \
            --include='FAILED' --exclude='*' \
            "$HOST:$RDIR/docs/figures/$id/" "$here/$id/"
        ls -l "$here/$id"
    done
    ;;
plot)
    [ $# -gt 0 ] || { echo "usage: d22.sh plot <id>..." >&2; exit 2; }
    for id in "$@"; do (cd "$here/$id" && python3 plot.py); done
    ;;
*)
    sed -n '2,16p' "$0"; exit 2
    ;;
esac
