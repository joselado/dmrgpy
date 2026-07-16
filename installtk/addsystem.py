import os
import platform


def _shellrc_path():
    shell = os.environ.get("SHELL", "")
    if "bash" in shell: shellrc = ".bashrc"
    elif "zsh" in shell: shellrc = ".zshrc"
    else:
        print("Unrecognised shell '"+shell+"', not editing any shell rc "
              "file. Set DMRGROOT manually if you need it (see "
              "examples/transverse_2d_ising_model).")
        return None

    if platform.system() == "Linux":
        return os.path.join(os.environ["HOME"], shellrc)

    # Mac: prefer whichever profile file already exists
    for profile in ("/.bash_profile", "/.bash_login", "/.profile"):
        path = os.environ["HOME"]+profile
        if os.path.isfile(path):
            return path
    return os.environ["HOME"]+"/.bash_profile" # none exist yet, create this one


def addrc():
    """Add (idempotently) a DMRGROOT export to the user's shell rc file.
    A handful of examples (e.g. transverse_2d_ising_model) read this
    environment variable directly."""
    shellrc = _shellrc_path()
    if shellrc is None:
        return

    pwd = os.path.dirname(os.path.realpath(__file__))+"/.."
    dmrgroot = os.path.realpath(pwd+"/src")
    marker = "export DMRGROOT=\""+dmrgroot+"\""

    try:
        current = open(shellrc).read()
    except FileNotFoundError:
        current = ""

    if marker in current:
        print("DMRGROOT already set correctly in "+shellrc)
        return

    route = "\n\n###############################\n"
    route += "  "+marker+"\n"
    route += "###############################\n"

    with open(shellrc, "a") as f:
        f.write(route)

    print("Added DMRGROOT="+dmrgroot+" to "+shellrc)
    print("Run 'source "+shellrc+"' or restart your shell to pick it up.")


if __name__ == "__main__":
    addrc()
