

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
    # the dense site set is the chain's reference; a conserved sector puts
    # a QN set built from the same codes on self.jlsites (sector.py), so a
    # rebuilt site set (setup_julia(), a fresh initialize()) is told about
    # the chain's sector again, as sites.py does for the session backends
    self.jlsites_dense = jlsites
    self._jl_sector_plan = None
    if getattr(self,"conserved_sector",None):
        self._apply_conserved_sector()
