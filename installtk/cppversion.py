import os

def cppversion(gpp="g++"):
    """
    Return the C++ version
    """
    os.system(gpp+" --version > /tmp/cpp.txt")
    out = open("/tmp/cpp.txt").read()
    out = out.split("\n")[0].split() # get the first line
    out = out[-1] # last one
    out = out.split(".")
    out = [int(o) for o in out]
    return out


def correct_version(**kwargs):
    try: out = cppversion(**kwargs) # get the version
    except: return False
    if out[0]>=6: return True
    else: return False


if __name__=="__main__":
    print(cppversion())
    print(correct_version())
