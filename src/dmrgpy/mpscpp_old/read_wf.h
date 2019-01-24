// this does not work yet

auto read_wf(std::string name="psi_GS.mps") {
  auto sites = get_sites();
//  readFromFile("sites_file",sites);
  sites = get_sites() ;
  readFromFile("sites.sites",sites);
  auto psi = MPS(sites);
  readFromFile(name,psi);
  return psi ;
}
