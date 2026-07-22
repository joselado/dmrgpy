

def sites_text(self):
    text = []
    text = [str(self.ns)] # number of sites
    for si in self.sites: 
        text += [str(si)] # each site
    return text


def initialize(self):
    """Initialize the sites for the Julia calculation"""
    ls = sites_text(self)
    from .juliasession import Main, to_julia_strvec
    jlsites = Main.get_sites(to_julia_strvec(ls))
    self.jlsites = jlsites
