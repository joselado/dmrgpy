
// return the hamiltonian of the system
#include"get_field.h" // get local magnetic fields
#include"get_exchange.h" // get the exchange coupling
#include"get_hopping.h" // get the hoppings (in case there are)
#include"get_spinful_hopping.h" // get the hoppings (in case there are)
#include"get_pairing.h" // get the hoppings (in case there are)
#include"get_hubbard.h" // get hubbard
#include"get_vijkl.h" // get generalized interaction


static auto get_ampo=[](auto sites) {
    auto ampo = AutoMPO(sites); // create the MPO for the Hamiltonian
    if (get_bool("use_ampo_hamiltonian")) {
	    ampo = get_ampo_operator(ampo,"hamiltonian.in");
	    return ampo;
	    };
    cout << "Adding exchange couplings" << endl ;
    ampo = get_exchange(ampo); // add exchange couplings to the Hamiltonian
    cout << "Adding fermionic hoppings" << endl ;
    ampo = get_hopping(ampo); // add hopping to the Hamiltonian
    cout << "Adding fermionic spinful hoppings" << endl ;
    ampo = get_spinful_hopping(ampo); // add hopping to the Hamiltonian
    cout << "Adding Hubbard interaction" << endl ;
    ampo = get_hubbard(ampo); // add hubbard to the Hamiltonian
    cout << "Adding generalized interaction" << endl ;
    ampo = get_vijkl(ampo); // add generalized interaction
    cout << "Adding magnetic field" << endl ;
    ampo = get_field(sites,ampo); // add magnetic field to the Hamiltonian
    cout << "Adding superconducting pairing" << endl ;
    ampo = get_pairing(ampo); // add pairing to the Hamiltonian
    return ampo ;
}
;




static auto get_hamiltonian=[](auto sites) {
    auto ampo = get_ampo(sites) ; // get the ampo
    auto H = MPO(ampo);  // create the full Hamiltonian
      // workaround for generic interactions not supported by iTensor
      if (get_bool("use_multioperator_hamiltonian")) {
              cout << "WARNING, using multioperator hamiltonian" << endl ;
	      auto hm = get_multioperator("hamiltonian_multioperator",sites);
	      auto out = sum_mpo(H,hm); // sum the contributions
	      return out; // return the sum
      };
    return H ; // return the Hamiltonian
}
;










