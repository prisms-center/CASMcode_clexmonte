from skbuild import setup

setup(
    name="libcasm-clexmonte",
    version="2.0a2",
    packages=[
        "libcasm",
        "libcasm.clexmonte",
    ],
    package_dir={"": "python"},
    cmake_install_dir="python/libcasm",
    include_package_data=False,
)
