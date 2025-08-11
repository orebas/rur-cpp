#!/bin/bash

# Convenience script for running RUR C++ tests with common options

set -e

if [ ! -f "build/rur_tests" ]; then
    echo "Error: rur_tests executable not found."
    echo "Please run './setup_testing.sh' first to build the test framework."
    exit 1
fi

cd build

case "${1:-all}" in
    "all")
        echo "Running all tests..."
        ./rur_tests
        ;;
    "basic")
        echo "Running basic operations tests..."
        ./rur_tests --gtest_filter="BasicOperations.*"
        ;;
    "univariate")
        echo "Running univariate system tests..."
        ./rur_tests --gtest_filter="UnivariateSystem.*"
        ;;
    "linear")
        echo "Running linear system tests..."
        ./rur_tests --gtest_filter="LinearSystem.*"
        ;;
    "nonlinear")
        echo "Running nonlinear system tests..."
        ./rur_tests --gtest_filter="NonlinearSystem.*"
        ;;
    "roots")
        echo "Running numerical root tests..."
        ./rur_tests --gtest_filter="NumericalRoots.*"
        ;;
    "integration")
        echo "Running integration tests..."
        ./rur_tests --gtest_filter="Integration.*"
        ;;
    "quick")
        echo "Running quick test suite..."
        ./rur_tests --gtest_filter="BasicOperations.MonomialCreation:LinearSystem.Simple2x2System:UnivariateSystem.QuadraticRoots"
        ;;
    "list")
        echo "Available tests:"
        ./rur_tests --gtest_list_tests
        ;;
    "verbose")
        echo "Running all tests with verbose output..."
        ./rur_tests --verbose
        ;;
    "xml")
        echo "Running all tests with XML output..."
        ./rur_tests --gtest_output=xml:test_results.xml
        echo "Results saved to test_results.xml"
        ;;
    "help")
        echo "Usage: ./run_tests.sh [category]"
        echo ""
        echo "Categories:"
        echo "  all         - Run all tests (default)"
        echo "  basic       - Basic polynomial operations"
        echo "  univariate  - Single-variable polynomial systems"  
        echo "  linear      - Linear equation systems"
        echo "  nonlinear   - Nonlinear polynomial systems"
        echo "  roots       - Numerical root-finding tests"
        echo "  integration - End-to-end integration tests"
        echo "  quick       - Run a few representative tests"
        echo "  list        - List all available tests"
        echo "  verbose     - Run all tests with verbose output"
        echo "  xml         - Run all tests with XML output"
        echo "  help        - Show this help message"
        echo ""
        echo "Custom filters:"
        echo "  ./run_tests.sh \"*Circle*\"     - Run tests matching pattern"
        echo ""
        echo "Direct gtest options:"
        echo "  cd build && ./rur_tests --gtest_filter=\"pattern\""
        echo "  cd build && ./rur_tests --gtest_help"
        ;;
    *)
        echo "Running tests with custom filter: $1"
        ./rur_tests --gtest_filter="$1"
        ;;
esac