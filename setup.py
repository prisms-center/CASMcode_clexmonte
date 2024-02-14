from skbuild import setup

setup(
    name="libcasm-configuration",
    version="2.0a2",
    packages=[
        "libcasm",
        "libcasm.clusterography",
        "libcasm.configuration",
        "libcasm.configuration.io",
        "libcasm.enumerate",
        "libcasm.irreps",
        "libcasm.occ_events",
        "libcasm.sym_info",
    ],
    package_dir={"": "python"},
    cmake_install_dir="python/libcasm",
    include_package_data=False,
)
