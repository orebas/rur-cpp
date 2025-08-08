#include "rur_solver.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>

void create_test_systems() {
    // Create simple system: x^2 - 1, y^2 - 1
    nlohmann::json simple_system;
    simple_system["variables"] = {"x", "y"};
    simple_system["prime"] = 100003;
    
    nlohmann::json polynomials = nlohmann::json::array();
    
    // x^2 - 1
    nlohmann::json poly1 = nlohmann::json::array();
    poly1.push_back({{"coefficient", 1}, {"exponents", {2, 0}}});   // x^2
    poly1.push_back({{"coefficient", -1}, {"exponents", {0, 0}}});  // -1
    
    // y^2 - 1
    nlohmann::json poly2 = nlohmann::json::array();
    poly2.push_back({{"coefficient", 1}, {"exponents", {0, 2}}});   // y^2
    poly2.push_back({{"coefficient", -1}, {"exponents", {0, 0}}});  // -1
    
    polynomials.push_back(poly1);
    polynomials.push_back(poly2);
    simple_system["polynomials"] = polynomials;
    
    std::ofstream file("simple_system.json");
    file << simple_system.dump(2) << std::endl;
    
    // Create cyclic-3 system
    nlohmann::json cyclic3;
    cyclic3["variables"] = {"x", "y", "z"};  
    cyclic3["prime"] = 100003;
    
    nlohmann::json polys3 = nlohmann::json::array();
    
    // x + y + z
    nlohmann::json p1 = nlohmann::json::array();
    p1.push_back({{"coefficient", 1}, {"exponents", {1, 0, 0}}});  // x
    p1.push_back({{"coefficient", 1}, {"exponents", {0, 1, 0}}});  // y
    p1.push_back({{"coefficient", 1}, {"exponents", {0, 0, 1}}});  // z
    
    // xy + yz + zx
    nlohmann::json p2 = nlohmann::json::array();
    p2.push_back({{"coefficient", 1}, {"exponents", {1, 1, 0}}});  // xy
    p2.push_back({{"coefficient", 1}, {"exponents", {0, 1, 1}}});  // yz
    p2.push_back({{"coefficient", 1}, {"exponents", {1, 0, 1}}});  // zx
    
    // xyz - 1
    nlohmann::json p3 = nlohmann::json::array();
    p3.push_back({{"coefficient", 1}, {"exponents", {1, 1, 1}}});   // xyz
    p3.push_back({{"coefficient", -1}, {"exponents", {0, 0, 0}}});  // -1
    
    polys3.push_back(p1);
    polys3.push_back(p2);
    polys3.push_back(p3);
    cyclic3["polynomials"] = polys3;
    
    std::ofstream file2("cyclic3.json");
    file2 << cyclic3.dump(2) << std::endl;
}

void test_json_system(const std::string& filename) {
    std::cout << "Testing system from " << filename << std::endl;
    
    std::ifstream file(filename);
    nlohmann::json j;
    file >> j;
    
    auto vars = j["variables"].get<std::vector<std::string>>();
    int prime = j["prime"].get<int>();
    
    std::cout << "Variables: ";
    for (const auto& var : vars) std::cout << var << " ";
    std::cout << std::endl;
    std::cout << "Prime: " << prime << std::endl;
    
    try {
        RURSolver<int> solver(prime, vars);
        
        // Load polynomials from JSON
        for (const auto& poly_json : j["polynomials"]) {
            auto poly = MultivariatePolynomial<int>::from_json(poly_json);
            std::cout << "Adding polynomial: " << poly.to_string(vars) << std::endl;
            solver.add_polynomial(poly);
        }
        
        std::cout << "Computing Gröbner basis..." << std::endl;
        solver.compute_groebner_basis();
        
        auto gb = solver.get_groebner_basis();
        std::cout << "Gröbner basis (" << gb.size() << " polynomials):" << std::endl;
        
        for (size_t i = 0; i < gb.size(); ++i) {
            std::cout << "  GB[" << i << "]: " << gb[i].to_string(vars) << std::endl;
        }
        
        std::cout << "Is zero-dimensional: " << (solver.is_zero_dimensional() ? "yes" : "no") << std::endl;
        
        if (solver.is_zero_dimensional()) {
            std::cout << "Computing RUR..." << std::endl;
            solver.compute_rur();
            std::cout << "Solution count: " << solver.get_solution_count() << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    std::cout << std::endl;
}

int main() {
    std::cout << "Creating test systems..." << std::endl;
    create_test_systems();
    
    std::cout << "\n=== Testing JSON-based polynomial input ===" << std::endl;
    test_json_system("simple_system.json");
    test_json_system("cyclic3.json");
    
    return 0;
}