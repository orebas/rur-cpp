#include <iostream>
#include <regex>
#include <string>

void test_regex(const std::string& input) {
    std::cout << "Testing input: '" << input << "'" << std::endl;
    
    // Current regex
    std::regex term_regex(R"([+-]?\d+(?:\*[a-zA-Z_][a-zA-Z0-9_]*(?:\^?\d+)?)*?)");
    std::sregex_iterator iter(input.begin(), input.end(), term_regex);
    std::sregex_iterator end;
    
    int count = 0;
    for (; iter != end; ++iter) {
        std::string term = iter->str();
        std::cout << "  Term " << count++ << ": '" << term << "'" << std::endl;
    }
    
    // Better regex that splits on +/- at the beginning of terms
    std::cout << "Using better approach:" << std::endl;
    
    // Split manually by finding +/- not at the start
    std::string modified = input;
    if (modified[0] != '+' && modified[0] != '-') {
        modified = "+" + modified;
    }
    
    // Find all positions of +/- that are not at the beginning
    std::vector<size_t> positions;
    positions.push_back(0);
    
    for (size_t i = 1; i < modified.length(); ++i) {
        if (modified[i] == '+' || modified[i] == '-') {
            positions.push_back(i);
        }
    }
    positions.push_back(modified.length());
    
    for (size_t i = 0; i < positions.size() - 1; ++i) {
        std::string term = modified.substr(positions[i], positions[i+1] - positions[i]);
        std::cout << "  Manual term " << i << ": '" << term << "'" << std::endl;
    }
}

int main() {
    test_regex("+1*y^2+100002");
    std::cout << std::endl;
    test_regex("+1*x^2+100002");
    std::cout << std::endl;
    test_regex("3*x^2+5*y*z-7");
    
    return 0;
}