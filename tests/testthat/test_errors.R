context("Error handling")

expect_error(fic(wide=wide.glm, inds=c(inds1,0), inds0=inds0, focus=focus, X=X),
             "Length of `inds` must match number of parameters")
expect_error(fic(wide=wide.glm, inds=inds1, inds0=inds0[-1], focus=focus, X=X),
             "Length of `inds0` must match number of parameters")

expect_error(fic(wide=wide.glm, inds=inds, inds0=inds0, focus="nonexistent_function", X=vals.first),
             "not found")
