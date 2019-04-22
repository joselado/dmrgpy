import os

def cppversion():
    """
    Return the C++ version
    """
    os.system("g++ --version > /tmp/cpp.txt")
    out = open("/tmp/cpp.txt").read()
    out = out.split("\n")[0].split() # get the first line
    out = out[-1] # last one
    out = out.split(".")
    out = [int(o) for o in out]
    return out


def correct_version():
    try: out = cppversion() # get the version
    except: return False
    if out[0]>=6: return True
    else: return False


if __name__=="__main__":
    print(cppversion())
    print(correct_version())
