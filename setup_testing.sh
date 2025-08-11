#!/bin/bash

# Setup script for RUR C++ Modern Testing Framework
# This script installs Google Test and configures the testing environment

set -e

echo "=== RUR C++ Testing Framework Setup ==="

# Check if vcpkg exists
if [ ! -d "vcpkg" ]; then
    echo "Error: vcpkg directory not found. Please ensure vcpkg is properly set up."
    exit 1
fi

echo "Installing Google Test via vcpkg..."
./vcpkg/vcpkg install gtest

echo "Checking Google Test installation..."
GTEST_LIB="./vcpkg/installed/x64-linux/lib/libgtest.a"
if [ -f "$GTEST_LIB" ]; then
    echo "✓ Google Test library found: $GTEST_LIB"
else
    echo "✗ Google Test library not found. Installation may have failed."
    echo "   Expected location: $GTEST_LIB"
    exit 1
fi

GTEST_INCLUDE="./vcpkg/installed/x64-linux/include/gtest/gtest.h"
if [ -f "$GTEST_INCLUDE" ]; then
    echo "✓ Google Test headers found"
else
    echo "✗ Google Test headers not found. Installation may have failed."
    echo "   Expected location: $GTEST_INCLUDE"
    exit 1
fi

echo ""
echo "Configuring build system..."
if [ -d "build" ]; then
    echo "Cleaning existing build directory..."
    rm -rf build
fi

mkdir build
cd build

echo "Running CMake configuration..."
cmake ..

if [ $? -eq 0 ]; then
    echo "✓ CMake configuration successful"
else
    echo "✗ CMake configuration failed"
    exit 1
fi

echo ""
echo "Building test framework..."
make rur_tests -j$(nproc)

if [ $? -eq 0 ]; then
    echo "✓ Build successful"
else
    echo "✗ Build failed"
    exit 1
fi

echo ""
echo "=== Setup Complete ==="
echo ""
echo "Run tests with:"
echo "  cd build"
echo "  ./rur_tests                           # All tests"
echo "  ./rur_tests --gtest_list_tests        # List available tests"  
echo "  ./rur_tests --gtest_filter=\"*Linear*\" # Run only linear tests"
echo "  ./rur_tests --gtest_help              # Show all options"
echo ""
echo "Test categories available:"
echo "  - BasicOperations.*    : Basic polynomial operations"
echo "  - UnivariateSystem.*   : Single-variable systems"
echo "  - LinearSystem.*       : Linear equation systems" 
echo "  - NonlinearSystem.*    : Nonlinear polynomial systems"
echo "  - NumericalRoots.*     : Root-finding algorithms"
echo "  - Integration.*        : End-to-end system tests"
echo ""
echo "For adding new tests, see: test/gtest/TESTING_TEMPLATE.md"