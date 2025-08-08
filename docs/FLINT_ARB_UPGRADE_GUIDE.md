# FLINT/ARB Support and Upgrade Guide

## Current Status

The rur-cpp project has been upgraded to support both FLINT with ARB functionality and a fallback mode using Eigen for polynomial root finding. By default, ARB support is **disabled** to avoid compatibility issues.

## Quick Start

### Building with Eigen Fallback (Default - Recommended)

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
ninja
```

This uses Eigen for polynomial root finding and should work reliably.

### Testing the Configuration

```bash
./test_flint_arb_support
```

This will show whether ARB support is enabled and test the polynomial root finding functionality.

## ARB Support Details

### Why ARB Support is Disabled by Default

The current vcpkg configuration provides:
- FLINT 2.9.0 (without integrated ARB)
- ARB 2.21.1 (separate library)

These versions have compatibility issues due to symbol mismatches between ARB expectations and FLINT exports. ARB was compiled expecting different FLINT symbol signatures than what FLINT 2.9.0 provides.

### Enabling ARB Support (Advanced)

If you want to enable ARB support and are willing to handle compatibility issues:

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_FLINT_ARB=ON
ninja
```

**Note**: This will likely fail with linking errors like:
```
undefined reference to `fmpz_add'
undefined reference to `fmpz_sub'
undefined reference to `arb_zero'
```

## Recommended Upgrade Paths

### Option 1: Use FLINT 3.x with Integrated ARB (Recommended)

FLINT 3.0+ integrates ARB functionality directly, eliminating compatibility issues.

#### Manual Installation of FLINT 3.x

1. **Remove vcpkg FLINT/ARB:**
   ```bash
   ./vcpkg/vcpkg remove flint arb
   ```

2. **Install FLINT 3.x from source:**
   ```bash
   # Install dependencies
   sudo apt-get install libgmp-dev libmpfr-dev
   
   # Clone and build FLINT 3.x
   git clone https://github.com/flintlib/flint.git
   cd flint
   git checkout v3.1.3  # or latest 3.x release
   
   ./configure --prefix=/usr/local --with-gmp --with-mpfr
   make -j$(nproc)
   sudo make install
   ```

3. **Update CMakeLists.txt to use system FLINT:**
   ```cmake
   # Replace vcpkg-specific paths with system paths
   find_library(FLINT_LIBRARY NAMES flint PATHS /usr/local/lib)
   find_path(FLINT_INCLUDE_DIR flint/flint.h PATHS /usr/local/include)
   ```

4. **Enable ARB support:**
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_FLINT_ARB=ON
   ```

#### Using Alternative Package Managers

**Conda/Mamba:**
```bash
mamba install flint  # This typically provides FLINT 3.x
```

**Homebrew (macOS):**
```bash
brew install flint   # Usually provides FLINT 3.x
```

### Option 2: Custom vcpkg Port for FLINT 3.x

Create a custom vcpkg port for FLINT 3.x:

1. Create `ports/flint3/portfile.cmake` in your vcpkg directory
2. Configure it to build FLINT 3.x with ARB support
3. Install via `./vcpkg install flint3`

### Option 3: Use Docker/Container with FLINT 3.x

Create a containerized environment with proper FLINT 3.x installation:

```dockerfile
FROM ubuntu:24.04
RUN apt-get update && apt-get install -y \
    build-essential cmake ninja-build \
    libgmp-dev libmpfr-dev libflint-dev
# FLINT package in Ubuntu 24.04+ typically includes ARB support
```

## Benefits of Full ARB Support

When properly configured with compatible FLINT/ARB versions, you get:

1. **Certified Root Finding**: Rigorous error bounds on polynomial roots
2. **Arbitrary Precision**: Work with polynomials requiring high precision
3. **Complex Arithmetic**: Native support for complex interval arithmetic
4. **Performance**: Optimized algorithms for polynomial root isolation

## Current Fallback Behavior

Without ARB support, the system:
- Uses Eigen for polynomial root finding
- Provides reasonable accuracy for most use cases
- Returns estimated error bounds (conservative estimates)
- Works reliably with the current build system

## Testing Your Configuration

The test suite `test_flint_arb_support` will:
- Show which backend is active (ARB or Eigen)
- Test basic polynomial root finding
- Verify error bound computation
- Validate configuration options

Example output:
```
=== FLINT/ARB Support Test Suite ===
ARB support: DISABLED (using Eigen fallback)
Testing FLINT/ARB polynomial root finding...
Found 2 roots for x^2 - 1 = 0:
  Root 0: 1 + 0i
  Root 1: -1 + 0i
✓ All tests passed!
```

## Troubleshooting

### Build Fails with ARB Enabled

This is expected with the current vcpkg setup. The ARB and FLINT versions are incompatible. Use the fallback mode or upgrade to FLINT 3.x.

### Performance Concerns

The Eigen fallback provides good performance for typical polynomial sizes. If you need maximum performance or certified error bounds, pursue the FLINT 3.x upgrade.

### Precision Requirements

If your application requires arbitrary precision polynomial arithmetic, upgrading to FLINT 3.x is recommended. The Eigen fallback uses double precision.

## Future Plans

- Monitor vcpkg for FLINT 3.x availability
- Consider providing custom vcpkg ports
- Evaluate alternative arbitrary precision libraries
- Add more configuration options for precision control

## Support

If you encounter issues:
1. First try the Eigen fallback mode (default)
2. Check the test output with `./test_flint_arb_support`
3. For FLINT 3.x upgrades, refer to FLINT documentation
4. Consider the containerized approach for complex setups