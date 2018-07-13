// this does not work yet

auto read_wf() {
  auto sites = get_sites();
//  readFromFile("sites_file",sites);
  sites = get_sites() ;
  readFromFile("sites.dmrg",sites);
  auto psi = MPS(sites);
  readFromFile("psi_GS.dmrg",psi);
  return psi ;
}
