#!/bin/bash
# Setup script for RUR-CPP dependencies
# This script downloads and builds FLINT 3.3.1 with ARB support

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}RUR-CPP Dependency Setup${NC}"
echo -e "${GREEN}========================================${NC}"

# Check for required system dependencies
echo -e "\n${YELLOW}Checking system dependencies...${NC}"

check_command() {
    if ! command -v $1 &> /dev/null; then
        echo -e "${RED}Error: $1 is not installed${NC}"
        return 1
    else
        echo -e "${GREEN}✓${NC} $1 found"
        return 0
    fi
}

# Check for required build tools
MISSING_DEPS=0
check_command gcc || MISSING_DEPS=1
check_command g++ || MISSING_DEPS=1
check_command make || MISSING_DEPS=1
check_command cmake || MISSING_DEPS=1
check_command git || MISSING_DEPS=1
check_command autoconf || MISSING_DEPS=1
check_command libtool || MISSING_DEPS=1

# Check for required libraries using pkg-config or direct detection
echo -e "\n${YELLOW}Checking for required libraries...${NC}"

# Check for GMP
if pkg-config --exists gmp 2>/dev/null || [ -f /usr/include/gmp.h ]; then
    echo -e "${GREEN}✓${NC} GMP found"
else
    echo -e "${RED}Error: GMP not found${NC}"
    MISSING_DEPS=1
fi

# Check for MPFR
if pkg-config --exists mpfr 2>/dev/null || [ -f /usr/include/mpfr.h ]; then
    echo -e "${GREEN}✓${NC} MPFR found"
else
    echo -e "${RED}Error: MPFR not found${NC}"
    MISSING_DEPS=1
fi

# Check for Eigen3
if pkg-config --exists eigen3 2>/dev/null || [ -d /usr/include/eigen3 ]; then
    echo -e "${GREEN}✓${NC} Eigen3 found"
else
    echo -e "${RED}Error: Eigen3 not found${NC}"
    MISSING_DEPS=1
fi

if [ $MISSING_DEPS -eq 1 ]; then
    echo -e "\n${RED}Missing dependencies detected!${NC}"
    echo -e "Please install the required packages first:"
    echo -e "\n${YELLOW}Ubuntu/Debian:${NC}"
    echo -e "  sudo apt-get update"
    echo -e "  sudo apt-get install build-essential cmake git autoconf libtool"
    echo -e "  sudo apt-get install libgmp-dev libmpfr-dev libeigen3-dev libgtest-dev"
    echo -e "\n${YELLOW}Fedora/RHEL:${NC}"
    echo -e "  sudo dnf install gcc gcc-c++ make cmake git autoconf libtool"
    echo -e "  sudo dnf install gmp-devel mpfr-devel eigen3-devel gtest-devel"
    echo -e "\n${YELLOW}macOS (with Homebrew):${NC}"
    echo -e "  brew install cmake git autoconf automake libtool gmp mpfr eigen googletest"
    exit 1
fi

# Check if FLINT is already built
if [ -d "external/flint-install/lib" ] && [ -f "external/flint-install/lib/libflint.a" ]; then
    echo -e "\n${GREEN}FLINT 3.3.1 already installed, skipping build...${NC}"
    echo -e "${GREEN}Setup complete! You can now build the project.${NC}"
    exit 0
fi

# Download and build FLINT
echo -e "\n${YELLOW}Setting up FLINT...${NC}"

# Create external directory
mkdir -p external
cd external

# Clean up any partial builds
if [ -d "flint" ]; then
    echo "Cleaning up previous FLINT source..."
    rm -rf flint
fi
if [ -d "flint-install" ]; then
    echo "Cleaning up previous FLINT installation..."
    rm -rf flint-install
fi

# Clone FLINT from GitHub
echo "Cloning FLINT from GitHub..."
git clone --depth 1 --branch v3.1.3-p1 https://github.com/flintlib/flint.git
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to clone FLINT from GitHub${NC}"
    echo -e "${YELLOW}Please check your internet connection and try again${NC}"
    exit 1
fi

# Build FLINT
cd flint

# Bootstrap if needed (for git version)
if [ -f "bootstrap.sh" ]; then
    echo "Bootstrapping FLINT..."
    ./bootstrap.sh
    if [ $? -ne 0 ]; then
        echo -e "${RED}Error: FLINT bootstrap failed${NC}"
        exit 1
    fi
fi
echo "Configuring FLINT (this may take a minute)..."
./configure --prefix=$(pwd)/../flint-install --with-gmp --with-mpfr --disable-shared
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: FLINT configuration failed${NC}"
    exit 1
fi

echo "Building FLINT (this will take several minutes)..."
echo "Using $(nproc) parallel jobs..."
make -j$(nproc)
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: FLINT build failed${NC}"
    exit 1
fi

echo "Installing FLINT..."
make install
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: FLINT installation failed${NC}"
    exit 1
fi

cd ../..

# Verify installation
if [ -f "external/flint-install/lib/libflint.a" ]; then
    echo -e "\n${GREEN}========================================${NC}"
    echo -e "${GREEN}✓ Setup complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo -e "\nFLINT 3.3.1 has been successfully installed to external/flint-install/"
    echo -e "\nYou can now build the project with:"
    echo -e "  ${YELLOW}mkdir -p build && cd build${NC}"
    echo -e "  ${YELLOW}cmake -DCMAKE_BUILD_TYPE=Release ..${NC}"
    echo -e "  ${YELLOW}make -j\$(nproc)${NC}"
else
    echo -e "${RED}Error: FLINT installation verification failed${NC}"
    exit 1
fi