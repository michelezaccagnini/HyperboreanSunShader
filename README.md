# Hyperborean Sun
Code for the Revision product.

# License
TODO

# Build Instructions
## Prequisites
* download and install Visual Studio Community. Be sure to enable the option for "Tools for CMake".
* download and install CMake; add it to your system PATH.
* cd to the source root.
* run `git submodule update --init --recursive`.
* create an out-of-source build folder; for example `[SOURCE_ROOT]/build` and cd there.
* run `cmake [SOURCE_ROOT]`.
* run `cmake --build . --config Release -- -m`.
* everything will be built; find the executable in `[BUILD_FOLDER]/Release/hyperborean-sun.exe`.
* enjoy!
