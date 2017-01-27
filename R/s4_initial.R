# Alphabetr <- setClass("Alphabetr",
#                       slots = c(results = "data.frame",
#                                 pairs = "matrix",
#                                 data = "list",
#                                 pData = "list"))
#
# setMethod("show",
#           signature = "Alphabetr",
#           definition = function(object) {
#             cat("An object of class ", class(object), "\n", sep = "")
#             cat(" ", nrow(object@data$alpha), " wells containing",
#                 ncol(object@data$alpha), " unique alpha chains and ",
#                 ncol(object@data$beta), "unique beta chains", "\n", sep = "")
#             if (nrow(pairs) > 0) {
#               cat("No alpha-beta pairs have been determined.")
#             } else {
#               cat(nrow(pairs), " alpha-beta pairs have been determined.", sep = "")
#             }
#           })
#
