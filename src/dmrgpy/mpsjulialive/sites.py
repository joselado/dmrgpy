

def sites_text(self):
    text = []
    text = [str(self.ns)] # number of sites
    for si in self.sites: 
        text += [str(si)] # each site
    return text


def initialize(self):
    """Initialize the sites for the Julia calculation"""
    ls = sites_text(self)
    from .juliasession import Main
    Main.sites_string = ls
#    jlsites = Main.get_sites(Main.sites_string)
    jlsites = Main.eval('get_sites(sites_string)')
    self.jlsites = jlsites
