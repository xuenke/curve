cc_library(
    name = "polynomial",
    hdrs = ["polynomial_fit.h"],
)
cc_binary(
    name = "fit",
    srcs = ["main.cc"],
    deps = [
        #":read_write",
        ":polynomial",
    ],
)
