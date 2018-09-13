
// return the hamiltonian of the system
#include"get_exchange.h" // get the exchange coupling
#include"get_hopping.h" // get the hoppings (in case there are)
#include"get_pairing.h" // get the hoppings (in case there are)
#include"get_hubbard.h" // get the hoppings (in case there are)


auto get_hamiltonian(auto sites) {
    auto ampo = AutoMPO(sites); // create the MPO for the Hamiltonian
    cout << "Adding exchange couplings" << endl ;
    ampo = get_exchange(ampo); // add exchange to the Hamiltonian
    cout << "Adding fermionic hoppings" << endl ;
    ampo = get_hopping(ampo); // add hopping to the Hamiltonian
    cout << "Adding Hubbard interaction" << endl ;
    ampo = get_hubbard(ampo); // add hubbard to the Hamiltonian
    cout << "Adding magnetic field" << endl ;
    ampo = get_field(sites,ampo); // add magnetic field to the Hamiltonian
    cout << "Adding superconducting pairing" << endl ;
    ampo = get_pairing(ampo); // add magnetic field to the Hamiltonian
    auto H = MPO(ampo);  // create the full Hamiltonian
    return H ;
}

