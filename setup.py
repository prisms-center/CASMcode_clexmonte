from skbuild import setup

setup(
    name="libcasm-clexmonte",
    version="2.0a1",
    packages=[
        "libcasm",
        "libcasm.clexmonte",
        "libcasm.clexmonte.auto_configuration",
        "libcasm.clexmonte.semigrand_canonical",
    ],
    package_dir={"": "python"},
    cmake_install_dir="python/libcasm",
    include_package_data=False,
)
